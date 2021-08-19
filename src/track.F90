!------------------------------------------------------------------------------
! This is a module for tracking important quantities during a simulation.
! It's purpose is to provide an easy way to dump certain quantities to 
! data files during the course of a simulation.  This is intended to only
! supplement the post-code.  If possible, please implement any time-intensive
! analysis that can be done from the configuration files in the post-code.
!
! The control line format for using this module during a simulation is:
!   [TYPE_IDENTIFIER]  [update_freq]  [write_freq]  [output_filename]  opt_keys
! Example:
!   TOTAL_NRG  2 100  data_totnrg  GNUPLOT
! 
! Each line thus identifies one quantity to monitor as a function of simulation
! step number.  This quantity is calculated and dumped to the specified output
! file.
!
! Other keywords:
!   PERMOLEC -- causes energies to be reported on a per-molecule basis
!
! Needed Improvements:
! 1) add gnuplot file functionality
!------------------------------------------------------------------------------

Module track

  Use defaults, Only: RDbl, strLen, lstrLen, TOTAL_INDEX, calToJ, &
      STATS_BLOCKSIZE, xlstrlen, zero
  Use file, Only: file_open
  Use utils, Only: filesrchstr, stripcmnt, split, toint, &
      toupper, allocerrdisplay, int2str, real2str, str2seq, &
      checkandstop, deallocerrdisplay, findint, findstr
  Use general, Only: genparams
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag, vector_getunitvec, vector_display, &
      vector_getnorm
  Use stats, Only: stats_getvalue, Statistics, stats_getblock, stats_getcavg, &
      stats_getstd, stats_init, stats_update, stats_display
  Use molecules, Only: molecules_getnsorbs, molecules_name
  Use config, Only: AtMolCoords, config_getnmoles, config_getnatoms, &
      config_isfixed, config_getatype, config_getnmoleslist
  Use simcell, Only: SimCell_Params
  Use interact, Only: Interaction_Model, interact_extract, interact_usampling
  Use storebase, Only: EnergyPlus, storebase_initcopy, storebase_zero, &
      storebase_clean, storebase_display, storebase_sumintra, storebase_totnrg
  Use storetop, Only: Forcefield_Results, storetop_totnrg, storetop_fastsum, &
      storetop_display, storetop_extract
  Use dmap, Only: DensityMap, dmap_idstring, dmap_init, dmap_dump, &
      dmap_update, dmap_clean, dmap_display
      
  Implicit None
  Save

  Private
  Public :: Tracking_Info, Tracking_Item, track_init, track_inititem, &
      track_process, track_clean, track_display

  Type Tracking_Info
    Logical                             :: screendump,usampling
    Integer                             :: nitems,screendump_freq
    Integer, Dimension(:), Pointer      :: lastupdate
    Type(Tracking_Item), Dimension(:), Pointer      :: item
  End Type Tracking_Info

  Type Tracking_Item
    Logical                   :: gnuplot_file,init,usampling,permolec
    Integer                   :: update_freq,write_freq
    Integer                   :: blocksize,itemno
    Integer, Dimension(3)     :: subset
    Character(len=strLen)     :: itemtype,filename
    Type(Statistics)          :: value
    Type(DensityMap), Pointer :: dmap
  End Type Tracking_Item

  Character(len=lstrLen), Parameter :: track_tag = 'Tracking Information'

  !** A list of accepted tracking key strings and descriptions
  Integer, Parameter       :: track_ntypes = 7
  Character(len=12), Dimension(track_ntypes), Parameter :: track_itemtypes = &
      (/'TOTAL_NRG   ',  &    !* total system energy
        'INTRA_NRG   ',  &    !* total intramolecular energy
        'COUL_NRG    ',  &    !* total Coulombic energy
        'NONCOUL_NRG ',  &    !* total Non-Coulombic energy
        'NONINTRA_NRG',  &    !* total Coulombic + Non-Coul. energy (no intra)
        'MOLEFRACTION',  &    !* mole fraction of a specified species
        dmap_idstring(1:12)/) !* density map analysis

Contains

  !----------------------------------------------------------------------------
  ! Initializes the parameters from the control file
  ! Requires:  params -- Tracking parameter data structure
  !            imodel -- interaction model information
  !            ctrl_filename -- control filename to find init info
  !----------------------------------------------------------------------------
  Logical Function track_init(params,imodel,ctrl_filename)
    Type(Tracking_Info), Intent(InOut)  :: params
    Type(Interaction_Model), Intent(In) :: imodel
    Character(*), Intent(In)            :: ctrl_filename

    Integer                               :: i,j,unit,error,ios
    Integer                               :: nfields,nlines,item
    Character(len=lstrLen)                :: line
    Character(len=255)                    :: text
    Character(len=strLen), Dimension(20)  :: fields

    !** Set default
    track_init = .False.
    params%usampling = interact_usampling(imodel)

    !** Open the control file if it is not already open
    unit = file_open(ctrl_filename)

    !** Find the Tracking information section or exit if not there
    If (filesrchstr(unit,track_tag,text,.True.) == 0) Return

    !** Determine the number of lines before the next spacer or EOF
    params%nitems = 0
    nlines = 0
    Do
      Read(unit,'(a)',IOSTAT=ios) line
      If (ios /= 0) Exit
      If (Index(line,'----') /= 0) Exit
      If (Index(Trim(line),'SCREENDUMP') == 0) params%nitems = params%nitems + 1
      nlines = nlines + 1
    End Do

    If (params%nitems > 0) track_init = .True.

    !** Allocate space for the items
    Allocate(params%item(params%nitems), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'tracking info')

    !** Allocate space for the last update information
    Allocate(params%lastupdate(params%nitems), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'tracking update info')
    params%lastupdate = 0

    !** Find the beginning of the section again
    If (filesrchstr(unit,track_tag,text,.True.) == 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' Could not find tracking section the second time'
      Stop
    End If

    !** Set screen-dump defaults
    params%screendump = .False.
    params%screendump_freq = genparams%iprint

    !** Read the parameters and process them
    item = 0
    Do i = 1,nlines
      Read(unit,'(a)',IOSTAT=ios) line
      If (ios /= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            ' Unexpected error while reading tracking lines'
        Stop
      End If

      !** Check for the screen-dump command line
      If (Index(Trim(line),'SCREENDUMP') /= 0) Then
        line = stripcmnt(line)
        nfields = split(line,fields)
        If (nfields > 2) Then
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              ' screen dump command line should read: SCREENDUMP [frequency]'
          Stop
        End If

        params%screendump = .True.
        If (nfields == 2) Then
          params%screendump_freq = toint(fields(2), &
              'could not convert tracking screen dump frequency')
        End If
        Cycle
      End If

      !** Process line
      item = item + 1
      line = stripcmnt(line)
      Call track_inititem(params%item(item),params%usampling,line)
    End Do

  End Function track_init

  !----------------------------------------------------------------------------
  ! Initialize a single tracking item
  ! Requires:  params -- Tracking item data structure
  !            usampling -- flag indicating Umbrella sampling 
  !            line -- input line for one item
  !----------------------------------------------------------------------------
  Subroutine track_inititem(params,usampling,line)
    Type(Tracking_Item), Intent(InOut)    :: params
    Logical, Intent(In)                   :: usampling
    Character(*), Intent(In)              :: line

    Integer                               :: i,j,nfields,nchunks,error
    Character(len=255)                    :: text
    Character(len=strLen), Dimension(20)  :: fields,chunks

    !** Set defaults
    params%usampling = usampling
    params%permolec = .False.
    params%gnuplot_file = .False.
    params%subset = 0  !** whole system
    params%blocksize = STATS_BLOCKSIZE
    params%init = .False.  !** means that it hasn't been started yet

    !** Split the line and interpret the contents
    text = stripcmnt(line)
    nfields = split(text,fields)    

    !** Give error if there aren't enough parameters
    If (nfields < 4) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' Need more than four parameters for each tracking line'
      Write(0,'(1x,2a)') 'read line: ',Trim(text)
      Stop
    End If

    !** Make sure the itemtype matches one of the supported types
    params%itemno = findstr(track_itemtypes,Trim(fields(1)))
    If (params%itemno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' Could not interpret tracking keyword: "',Trim(fields(1)),'"'
      Write(0,'(2x,20(1x,a,1x))') 'Supported keywords are: ', &
          (Trim(track_itemtypes(i)),i=1,track_ntypes)
      Stop
    End If

    !** Store the parameters
    params%itemtype = Trim(fields(1))
    params%update_freq = toint(fields(2))
    params%write_freq = toint(fields(3))
    params%filename = fields(4)

    !** Interpret the additional options
    Do i = 5,nfields
      nchunks = split(fields(i),chunks,'@')
      Select Case (ToUpper(chunks(1)))
      Case ('ATM')
        params%subset(3) = toint(chunks(2))
        If (nchunks /= 2) Then
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              ' example of ATM option: ATM@2'
          Stop
        End If
      Case ('SPC')
        params%subset(1) = toint(chunks(2))
        If (nchunks /= 2) Then
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              ' example of SPC option: SPC@2'
          Stop
        End If
      Case ('GNUPLOT')
        params%gnuplot_file = .True.
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            ' The GNUPLOT option is not yet supported'
        Stop
      Case ('BLOCKSIZE')
        params%blocksize = toint(chunks(2))
        If (nchunks /= 2) Then
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              ' example of BLOCKSIZE option: BLOCKSIZE@100'
          Stop
        End If
      Case ('PERMOLEC')
        params%permolec = .True.
      Case Default
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Unable to interpret tracking information option "', &
            Trim(fields(i)),'"'
        Stop
      End Select
    End Do

    !** Handle the density map case
    If (params%itemno == 7) Then  
      Allocate(params%dmap,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'density map')      
      Call dmap_init(params%dmap,params%write_freq,params%subset,params%filename)
      params%init = .True.
    End If

  End Subroutine track_inititem

  !----------------------------------------------------------------------------
  ! Initialize the item statistics.  Used the first time it is accessed
  ! Requires:  params -- Tracking item data structure
  !----------------------------------------------------------------------------
  Subroutine track_startitem(params)
    Type(Tracking_Item), Intent(InOut)       :: params

    Integer         :: unit

    params%init = .True.
    Call stats_init(params%value,params%itemtype,params%blocksize, &
        params%usampling,'f10.3')

    !** Don't open the tracking file if this is a density map type analysis
    If (params%itemno == 7) Return
    
    unit = file_open(params%filename)
    Write(unit,'(2a,3x,a,i5)') '# Tracking information for quantity: ', &
        Trim(params%itemtype),'Block size = ',params%blocksize
    Write(unit,'(2a)') '# Step number, instantaneous value, block average, ', &
        'cumulative average and standard deviation'

  End Subroutine track_startitem

  !----------------------------------------------------------------------------
  ! Process the simulation data and perform specified tracking options
  ! Requires:  params -- Tracking item data structure
  !            species -- coordinate data structure 
  !            simcell -- simulation cell information
  !            imodel -- interaction model information
  !            iter -- iteration number
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------------
  Subroutine track_process(params,species,simcell,imodel,iter,indent,unit)
    Type(Tracking_Info), Intent(InOut)          :: params
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(SimCell_Params), Intent(In)            :: simcell  
    Type(Interaction_Model), Intent(InOut)      :: imodel
    Integer, Intent(In)                         :: iter,indent,unit

    Integer                   :: i,last
    Logical                   :: summed,updated,writetofile
    Character(len=indent)     :: blank
    Character(len=strLen)     :: string1,string2
    Character(len=xlstrLen)   :: line

    blank = Repeat(' ',indent)

    !** Process the simulation information
    summed = .False.
    Do i = 1,params%nitems
      Select Case(params%item(i)%itemno)
      Case(1:5)
        updated = track_updatestats(params%item(i),species,imodel,iter, &
            writetofile,summed)
      Case(6)
        updated = track_updatecmplx(params%item(i),species,simcell,imodel, &
            iter,writetofile)
      Case(7)
        If (Mod(iter,params%item(i)%update_freq) == 0) Then
          updated = dmap_update(params%item(i)%dmap,species,simcell,imodel, &
              iter,writetofile)
        End If
      Case Default
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Unable to interpret tracking identifier "', &
            Trim(params%item(i)%itemtype),'"'
        Stop
      End Select

      !** Finish the updating 
      If (updated) params%lastupdate(i) = iter
      If ((writetofile).And.(params%lastupdate(i) > 0)) &
          Call track_dumptofile(params%item(i),iter)
    End Do

    !** Do screen dump if necessary
    If (params%screendump) Then
      If (Mod(iter,params%screendump_freq) == 0) Then

        Write(unit,'(3a)') blank,'Current tracking information ', &
            '(blocksize, inst., block avg., cum. avg., std dev.)'
        
        !** Output for each of the tracking items
        Do i = 1,params%nitems
          If ((params%lastupdate(i) == 0).Or.(.Not. params%item(i)%init)) Then
            string1 = int2str(i)
            Write(unit,'(2x,4a)') blank,'Item ',Trim(string1), &
                ' Not yet initialized'
            Cycle
          End If

          last = iter - params%lastupdate(i)
          string2 = 'freshly updated'
          If (last /= 0) Then
            string1 = int2str(last)
            Write(string2,'(3a)') 'updated ',Trim(string1),' steps ago'
          End If

          line = track_infoline(params%item(i))              
          Write(unit,'(2x,2a,2x,a)') blank,Trim(line),Trim(string2)
        End Do
        Write(unit,*)
      End If
    End If

  End Subroutine track_process

  !----------------------------------------------------------------------------
  ! Handle generic updating of energy statistics.  Returns True if the 
  ! statistics were updated.
  ! Requires:  params -- Tracking item data structure
  !            species -- coordinate data structure 
  !            imodel -- interaction model information
  !            iter -- iteration number
  !            writetofile -- flag indicating if line should be written to file
  !            summed -- flag indicating if the storage structure is summed
  !----------------------------------------------------------------------------
  Logical Function track_updatestats(params,species,imodel,iter, &
      writetofile,summed)
    Type(Tracking_Item), Intent(InOut)          :: params
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(Interaction_Model), Intent(InOut)      :: imodel
    Integer, Intent(In)                         :: iter
    Logical, Intent(Out)                        :: writetofile
    Logical, Intent(InOut)                      :: summed

    Integer                 :: i,unit,nmoles
    Logical, SAVE           :: firsttime = .True.
    Logical                 :: success
    Real(kind=RDbl)         :: quantity,addup,normfactor
    Type(EnergyPlus), SAVE  :: components

    !** Check if we need to write to file
    writetofile = .False.
    If (Mod(iter,params%write_freq) == 0) writetofile = .True.

    !** Return if iteration doesn't match update frequency
    track_updatestats = .False.
    If (Mod(iter,params%update_freq) /= 0) Return
    track_updatestats = .True.

    !** Initialize the temporary energy structure for the first time
    If (firsttime) Then
      firsttime = .False.
      Call storebase_initcopy(components,imodel%results(1)%total)
    End If

    !** Make sure there are molecules of this species if its species-specific
    If (params%subset(1) /= 0) Then
      nmoles = config_getnmoles(species,params%subset(1))
      If (nmoles == 0) Return
    End If

    !** Initialize the statistics and open file if it's the first time for item
    If (.Not. params%init) Call track_startitem(params)

    !** Make sure the storage structure is properly summed
    If (.Not. summed) Then
      Call storetop_fastsum(imodel%results(1))
      summed = .True.
    End If

    !** Get the energy quantity to update
    Select Case(params%itemno)
    Case(1)  !** 'TOTAL_NRG'
      !** Get the total energy (including intra)
      success = interact_extract(imodel,iter,.True.,.True.,.True., &
          components,normfactor,params%subset)
      quantity = storebase_totnrg(components,.True.)*calToJ
      Call checkandstop(success,__FILE__,__LINE__, &
          ' Unable to extract total energy for tracking')

    Case(2)  !** 'INTRA_NRG'
      !** Get the total intramolecular energy 
      success = interact_extract(imodel,iter,.True.,.False.,.False., &
          components,normfactor,params%subset)
      quantity = storebase_totnrg(components,.True.)*calToJ
      Call checkandstop(success,__FILE__,__LINE__, &
          ' Unable to extract total intramolecular energy for tracking')   

    Case(3)  !** 'COUL_NRG'
      !** Get the total coulombic energy 
      success = interact_extract(imodel,iter,.False.,.False.,.True., &
          components,normfactor,params%subset)
      quantity = storebase_totnrg(components,.True.)*calToJ
      Call checkandstop(success,__FILE__,__LINE__, &
          ' Unable to extract total coulombic energy for tracking')   
    
    Case(4)  !** 'NONCOUL_NRG'
      !** Get the total non-coulombic energy 
      success = interact_extract(imodel,iter,.False.,.True.,.False., &
          components,normfactor,params%subset)
      quantity = storebase_totnrg(components,.True.)*calToJ
      Call checkandstop(success,__FILE__,__LINE__, &
          ' Unable to extract total non-coulombic energy for tracking')   
    
    Case(5)  !** 'NONINTRA_NRG'
      !** Get the total non-intramolecular energy (coul + noncoul)
      success = interact_extract(imodel,iter,.False.,.True.,.True., &
          components,normfactor,params%subset)
      quantity = storebase_totnrg(components,.True.)*calToJ
      Call checkandstop(success,__FILE__,__LINE__, &
          ' Unable to extract total non-intra energy for tracking')   
    
    Case Default
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' Unable to interpret tracking identifier "', &
          Trim(params%itemtype),'"'
      Stop
    End Select

    !** Adjust to per-molecule basis if specified
    If (params%permolec) Then
      nmoles = config_getnmoles(species,params%subset(1))
      If (nmoles == 0) Then
        quantity = zero
      Else
        quantity = quantity/nmoles
      End If
    End If

    !** Update the statisitcs
    Call stats_update(params%value,quantity)

  End Function track_updatestats

  !----------------------------------------------------------------------------
  ! Handle generic updating of complex statistics (ie those that use a lot of
  ! configuration information).  Returns True if the statistics were updated.
  ! Requires:  params -- Tracking item data structure
  !            species -- coordinate data structure 
  !            simcell -- simulation cell information
  !            imodel -- interaction model information
  !            iter -- iteration number
  !            writetofile -- flag indicating if line should be written to file
  !----------------------------------------------------------------------------
  Logical Function track_updatecmplx(params,species,simcell,imodel, &
      iter,writetofile)
    Type(Tracking_Item), Intent(InOut)          :: params
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(SimCell_Params), Intent(In)            :: simcell  
    Type(Interaction_Model), Intent(InOut)      :: imodel
    Integer, Intent(In)                         :: iter
    Logical, Intent(Out)                        :: writetofile

    Integer                 :: i,unit,nspc,spc
    Logical                 :: firsttime = .True., success
    Real(kind=RDbl)         :: quantity
    Integer, Dimension(20)  :: molcount

    !** Return if iteration doesn't match update frequency
    writetofile = .False.
    If (Mod(iter,params%write_freq) == 0) writetofile = .True.

    !** Return if iteration doesn't match update frequency
    track_updatecmplx = .False.
    If (Mod(iter,params%update_freq) /= 0) Return
    track_updatecmplx = .True.

    !** Initialize the statistics and open file if it's the first time for item
    If (.Not. params%init) Call track_startitem(params)

    !** Get the energy quantity to update
    Select Case(params%itemno)
    Case(1:5)  
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' Wrong update routine called for "', &
          Trim(params%itemtype),'"'
      Stop
    Case(6)  !** 'MOLFRAC'
      !** Get the mole fraction of a specified species
      nspc = molecules_getnsorbs()
      molcount = 0
      Call config_getnmoleslist(species,molcount)
      Do spc = 1,nspc
        If (config_isfixed(species(spc))) molcount(spc) = 0
      End Do
      quantity = molcount(params%subset(1))/(1.0_RDbl*Sum(molcount(1:nspc)))
    
    Case Default
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' Unable to interpret tracking identifier "', &
          Trim(params%itemtype),'"'
      Stop
    End Select

    !** Update the statisitcs
    Call stats_update(params%value,quantity)

  End Function track_updatecmplx

  !----------------------------------------------------------------------------
  ! Dump a tracking line to the file specified for the item
  ! Requires:  params -- Tracking item data structure
  !            iter -- iteration number
  !----------------------------------------------------------------------------
  Subroutine track_dumptofile(params,iter)
    Type(Tracking_Item), Intent(InOut)       :: params
    Integer, Intent(In)                      :: iter

    Integer                 :: i,unit

    !** Open the output file if it's not already open
    unit = file_open(params%filename)

    !** Write the instanteous, block, cum. avg. and std dev. to file
    Write(unit,'(i8,4f12.5)') iter,stats_getvalue(params%value), &
        stats_getblock(params%value), stats_getcavg(params%value), &
        stats_getstd(params%value)

  End Subroutine track_dumptofile

  !----------------------------------------------------------------------------
  ! Return one-line feedback giving the statistics information
  ! Requires:  params -- Tracking item data structure
  !----------------------------------------------------------------------------
  Function track_infoline(params)
    Character(len=xlstrLen)                  :: track_infoline
    Type(Tracking_Item), Intent(InOut)       :: params

    Character(len=lstrLen)    :: xtra
    Character(len=strLen)     :: string1,string2,string3,string4,string5

    !** Skip out if this is a density map
    If (params%itemno == 7) Then
      Write(track_infoline,'(a,t14,a)') &
          Trim(params%itemtype),'no information available'
      Return
    End If

    !** Extra information
    xtra = track_xtrainfoline(params)

    !** Write the instanteous, block, cum. avg. and std dev. to string
    string1 = int2str(params%blocksize)
    string2 = real2str(stats_getvalue(params%value),8)
    string3 = real2str(stats_getblock(params%value),8)
    string4 = real2str(stats_getcavg(params%value),8)
    string5 = real2str(stats_getstd(params%value),8)
    Write(track_infoline,'(a,t14,a,4(2x,a8),2x,a)') &
        Trim(params%itemtype),Trim(string1),Trim(string2), &
        Trim(string3),Trim(string4),Trim(string5),Trim(xtra)

    Return

    !** Longer version
    Write(track_infoline,'(a,2x,2a,2x,4(a,f12.5))') &
        Trim(params%itemtype),'Freq: ',Trim(string1), &
        'Inst: ',stats_getvalue(params%value), &
        'Block avg: ',stats_getblock(params%value), &
        'Cum. avg: ',stats_getcavg(params%value), &
        'Std. dev.: ',stats_getstd(params%value)

  End Function track_infoline

  !----------------------------------------------------------------------------
  ! Return one-line feedback giving extra information for one tracking item
  ! Requires:  params -- Tracking item data structure
  !----------------------------------------------------------------------------
  Function track_xtrainfoline(params)
    Character(len=xlstrLen)               :: track_xtrainfoline
    Type(Tracking_Item), Intent(In)       :: params

    Character(len=lstrLen)              :: string1

    !** Extra information
    track_xtrainfoline = ''
    Select Case(params%itemno)
    Case(1)  !** 'TOTAL_NRG'
    Case(2)  !** 'INTRA_NRG'
    Case(3)  !** 'COUL_NRG'
    Case(4)  !** 'NONCOUL_NRG'
    Case(5)  !** 'NONINTRA_NRG'
    Case(6)  !** 'MOLEFRACTION'
    End Select
    If (params%subset(1) /= 0) Then
      string1 = int2str(params%subset(1))
      Write(track_xtrainfoline,'(3a)') Trim(track_xtrainfoline), &
          'SPC=',Trim(string1)
    End If

    If (params%subset(3) /= 0) Then
      string1 = int2str(params%subset(3))
      Write(track_xtrainfoline,'(3a)') Trim(track_xtrainfoline), &
          ' ATM=',Trim(string1)
    End If

    If (params%permolec) Then
      Write(track_xtrainfoline,'(a,1x,a)') Trim(track_xtrainfoline),'PerMolec'
    End If

  End Function track_xtrainfoline

  !-----------------------------------------------------------------------
  ! Displays the parameters in the Tracking data structure
  ! Requires:  params -- Tracking parameter data structure
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !-----------------------------------------------------------------------
  Subroutine track_display(params,indent,unit)
    Type(Tracking_Info), Intent(In)    :: params
    Integer, Intent(In)                :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=lstrLen)      :: xtra
    Character(len=strLen)       :: string1,string2

    blank = Repeat(' ',indent)

    Write(unit,'(2a)') blank,'Tracking Information:'

    string1 = int2str(params%nitems)
    Write(unit,'(3a)') blank,' Number of quantities to track: ',Trim(string1)
    Do i = 1,params%nitems
      !** Get extra information for this item
      xtra = track_xtrainfoline(params%item(i))

      string1 = int2str(params%item(i)%update_freq)
      string2 = int2str(params%item(i)%write_freq)

      Write(unit,'(a,2x,i3,2x,a,2(1x,a),2x,a,2x,a)') blank,i, &
          Trim(params%item(i)%itemtype), &
          Trim(string1),Trim(string2), &
          Trim(params%item(i)%filename),Trim(xtra)
    End Do

  End Subroutine track_display 

  !-----------------------------------------------------------------------
  ! Cleans the Tracking data structure
  ! Requires:  params -- Tracking parameter data structure
  !-----------------------------------------------------------------------
  Subroutine track_clean(params)
    Type(Tracking_Info), Intent(InOut)    :: params

    Integer               :: error

    Deallocate(params%item, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)

  End Subroutine track_clean
  
End Module track



