!------------------------------------------------------------------------------
! This module handles the general species-species interactions.
! It contains a data type that specifies the interaction 
! parameters for a specific species-species pair.  Essentially, this 
! simple data type includes options for either turning the interaction 
! off, using a pre-tabulated map, using an explicit calculation with
! cutoffs or using a more sophisticated cutoff-free method.
! The matrix of all the possible pairs is also declared in this module. 
!
! In keeping with the different needs of simulation types it includes 
! different routines for evaluating various types of multi-atom interactions:
! ssdriver_ssint -- calculates species-species interactions
! ssdriver_msint -- calculates molecule-species interactions
! ssdriver_asint -- calculates atom-species interactions
!
! all routines give Level_Pair as output, together these level_pairs make
! up a complete storage structure.
!------------------------------------------------------------------------------

Module ssdriver

  Use defaults, Only: RDbl,strLen,xlstrlen,dashedline,dashedline2,d_ss_file, &
      MAX_SORBS, MAX_ATOMS, zero,dbgflag,lstrlen
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use utils, Only: split,fileSrchStr,allocErrDisplay,deallocErrDisplay,combine, &
      checkandstop,swap
  Use file, Only: file_gettype,file_open
  Use molecules, Only: molecules_getnsorbs, molecules_name, &
      molecules_getnatomtypes, molecules_getnatoms, molecules_gettype
  Use config, Only: AtMolCoords,config_isfixed,config_getq
  Use simcell, Only: SimCell_Params, simcell_getmolectype
  Use ssbasic, Only: SSBasicParams, ssbasic_init, ssbasic_ssint, &
      ssbasic_display, ssbasic_clean, ssbasic_idstring, ssbasic_msint, &
      ssbasic_asint, ssbasic_getpotparameters
  Use ssmap, Only: SSMapParams, ssmap_init, ssmap_idstring, &
      ssmap_eval, ssmap_ssint, ssmap_msint, ssmap_asint, ssmap_Oldidstring, &
      ssmap_display, ssmap_clean, ssmap_boxinfo
  Use sssum, Only: SSSumParams, sssum_init, sssum_initCalc, sssum_idstring, &
      sssum_display, sssum_clean, sssum_ssint, sssum_asintHOT, Assignment(=), sssum_intraon
  Use xternal, Only: ExternalInfo, EXTERNAL_KEY, xternal_init, xternal_display, &
      xternal_int
  Use store, Only: Store_Level_Pair,store_init,store_terminate
  Use storebase, Only: EnergyPlus, storebase_disp
  Use storetop, Only: storetop_level
  Use storesym, Only: Symmetric_Results, storesym_ptrs, storesym_initsspair, &
      storesym_terminate, storesym_display

  Implicit None
  Save

  Private
  Public :: SpcSpc_Params, ssdriver_init, ssdriver_initstore, ssdriver_link, &
      ssdriver_ssint, ssdriver_display, ssdriver_clean, ssdriver_ison, &
      ssdriver_msint, ssdriver_asint, ssdriver_listparams, ssdriver_nullify, &
      ssdriver_getpotparameters, ssdriver_boxinfo

  Type SpcSpc_Params
    Logical                          :: off,coul
    Type(SSBasicParams), Pointer     :: basic 
    Type(SSMapParams), Pointer       :: map 
    Type(SSSumParams), Pointer       :: sum
    Type(ExternalInfo), Pointer      :: xtrnal
  End Type SpcSpc_Params

  Character(len=strLen), Parameter    ::  nocoulomb_string = "NCOUL"
  Character(len=strLen), Parameter    ::  coulomb_string = "COUL"

Contains

  !------------------------------------------------------------------------
  ! Initializes Spc-spc interaction for a single species-species  
  ! pair.  Reads the spc-spc file to get information
  ! Requires:  params -- spc-spc forcefield parameters
  !            coulomb_flag -- True => initialize coulomb interactions
  !            spc1 -- number for species 1
  !            spc2 -- number for species 2
  !            filename -- filename to read from
  !            simcell -- simulation cell information
  !------------------------------------------------------------------------
  Subroutine ssdriver_init(params,coulomb_flag,spc1,spc2,filename,simcell)
    Type(SpcSpc_Params), Intent(Out)      :: params
    Logical, Intent(In)                   :: coulomb_flag
    Integer, Intent(In)                   :: spc1,spc2
    Character(*), Intent(In)              :: filename
    Type(SimCell_Params), Intent(In)      :: simcell

    Integer                               :: unitno,lineno,error
    Integer                               :: nfields,nchunks,nsets
    Integer                               :: simcellspc,mapspc,s1,s2
    Logical                               :: foundit
    Character(len=strLen)                 :: name1,name2
    Character(len=lstrLen)                :: subline,name
    Character(len=MAX_ATOMS*lstrLen)      :: line
    Character(len=strLen), Dimension(20)  :: fields,chunks
    Integer, Dimension(:,:), Allocatable  :: atomlist
    Character(len=lstrLen), Dimension(:), Allocatable :: paramslist

    !** get the names
    name1 = molecules_name(spc1)
    name2 = molecules_name(spc2)

    !** nullify all the pointers
    Call ssdriver_nullify(params)
    params%coul = coulomb_flag
    params%off = .True.

    !** Read the species-species pair interaction file
    unitno = file_open(filename,110)

    !** Search the file for the correct line
    foundit = .False.
    lineno = 1
    Do While((lineno /= 0).And.(.Not. foundit))
      If (coulomb_flag) Then
        lineno = fileSrchStr(unitno,(/name1,name2,coulomb_string/),line,.True.)
      Else
        lineno = fileSrchStr(unitno,(/name1,name2,nocoulomb_string/), &
            line,.True.)
      End If

      !** Split the line into its components and interpret parameters
      nfields = split(line,fields)
      nchunks = split(fields(4),chunks,'@')
      Select Case (chunks(1))  !** should be the identifier
      Case (ssbasic_idstring)
        foundit = .True.
        params%off = .False.
        Allocate(params%basic, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ffresults%basic')
        Call ssbasic_init(params%basic,spc1,spc2,simcell,line)

      Case (ssmap_idstring)
        foundit = .True.
        params%off = .False.

        !** Second chunk should be name of species interacting via the map
        If (nchunks == 1) Then
          simcellspc = simcell_getmolectype(simcell)
          mapspc = simcellspc
          If ((simcellspc /= spc1).And.(simcellspc /= spc2)) mapspc = spc2
          name = molecules_name(mapspc)
          Write(0,'(2a,i5,2a)') __FILE__,":",__LINE__,' WARNING: ', &
              'map species name should be specified like: MAP@species_name'
          Write(0,'(2x,3a)') 'Will assume that map species is: "',Trim(name),'"'
        Else
          mapspc = molecules_gettype(chunks(2))
        End If

        !** Make sure the map species is the last species in molecules
        If (mapspc /= molecules_getnsorbs()) Then
          Write(0,'(2a,i5,2a)') __FILE__,":",__LINE__, &
              ' Due to current storage limitations, '
          Write(0,'(2x,a)') 'map species must be last species in control file'
          Stop
        End If

        !** Make sure the map species is the second spc argument to ssmap_init
        s1 = spc1
        s2 = mapspc 
        If (spc1 == mapspc) s1 = spc2

        Allocate(params%map, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ffresults%basic')
        Call ssmap_init(params%map,s1,s2,simcell,line)

      Case (ssmap_Oldidstring)
        Write(0,'(2a,i5,4a)') __FILE__,":",__LINE__,'  ', &
            Trim(ssmap_Oldidstring),' option no longer used, switch to ', &
            Trim(ssmap_idstring)
        Stop

      Case (sssum_idstring)
        foundit = .True.
        params%off = .False.
        Allocate(params%sum, stat=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ffresults%sum')
        Call sssum_init(params%sum,spc1,spc2,combine(fields(5:nfields)))

      Case ('OFF')
        foundit = .True.

      Case ("")
        Write(0,'(a,i4,a,2(1x,a))') __FILE__,__LINE__, &
            ': Could not find COUL or NCOUL interaction specification for'
        Write(0,'(2x,a,2(1x,a))') 'species pair: ',Trim(name1),Trim(name2)
        Write(0,'(2x,2a)') 'Please check species-species file: ',Trim(filename)
        Stop

      Case(EXTERNAL_KEY)
        foundit = .True.
        params%off = .False.
        Allocate(params%xtrnal, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ffresults%xtrnal')

        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' External spc-spc evaluations not yet ready, under construction '
!        Stop

        subline = combine(fields(4:))
        Call xternal_init(params%xtrnal,subline,(/spc1,spc2/),simcell)
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Stop

      Case Default
        Write(0,'(a,i4,3a)') __FILE__,__LINE__, &
            ': Unable to identify SS interaction string ',Trim(fields(4))
        Stop
      End Select

      !** Check for the possibility that the format is BASIC@EXTERNAL
      !** in this case, we'll get the BASIC parameters and use them
      !** instead to call an external program.
      If (nchunks > 1) Then
        If (chunks(2) == EXTERNAL_KEY) Then
          If (.Not. Associated(params%basic)) Then
            Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
                ' External evaluations not yet available with ',Trim(fields(4))
            Stop          
          End If

          !** Get the forcefield parameters
          Allocate(atomlist(100,4), STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'atomlist')
          Allocate(paramslist(100), STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'paramslist')
          Call ssdriver_listparams(params,(/spc1/),(/spc2/),nsets, &
              atomlist,paramslist)

          Allocate(params%xtrnal, STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'xtrnal')
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          subline = combine(fields(4:))
          Call xternal_init(params%xtrnal,subline,(/spc1,spc2/),simcell, &
              atomlist(1:nsets,:),paramslist(1:nsets))
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Stop          
        End If
      End If

    End Do

    !** Give error feedback if necessary
    If (.Not. foundit) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' Could not find interaction information for spc-spc pair:'
      Write(0,'(4x,5a)') Trim(name1),'--',Trim(name2),' in file: ',Trim(filename)
      Stop
    End If
    
    Close(unit=unitno)

  End Subroutine ssdriver_init

  !------------------------------------------------------------------------
  ! Nullifies the Spc-spc interaction pointer set
  ! Requires:  params -- spc-spc forcefield parameters
  !------------------------------------------------------------------------
  Subroutine ssdriver_nullify(params)
    Type(SpcSpc_Params), Intent(InOut)      :: params

    Nullify(params%basic)
    Nullify(params%map)
    Nullify(params%sum)
    Nullify(params%xtrnal)

  End Subroutine ssdriver_nullify

  !------------------------------------------------------------------------
  ! Initialize storage for one species-species entry in the symmetric 
  ! results storage structure.
  ! Requires:  params -- spc-spc forcefield parameters
  !            storage -- symmetric results data structure to initialize
  !            spc1 -- 1st index, number for species 1
  !            nmoles1 -- number of molecules to initialize 
  !            spc2 -- 2nd index, number for species 2
  !            nmoles2 -- number of molecules to initialize 
  !            nderivs -- number of derivatives in storage
  !            store_level -- string indicating desired storage level
  !------------------------------------------------------------------------
  Subroutine ssdriver_initstore(params,storage,spc1,nmoles1,spc2,nmoles2, &
      nderivs,storelevel)
    Type(SpcSpc_Params), Intent(In)        :: params
    Type(Symmetric_Results), Intent(InOut) :: storage
    Integer, Intent(In)                    :: spc1,nmoles1,spc2,nmoles2,nderivs
    Character(len=3), Intent(In)           :: storelevel

    Integer       :: level
    Logical       :: intra_also

    !** Get the level of detail from the passed string
    level = storetop_level(storelevel)

    !** initialize the storage
    If (Associated(params%basic)) Then
      Call storesym_initsspair(storage,spc1,nmoles1,spc2,nmoles2,1,level, &
          nderivs,.False.)

    Else If (Associated(params%map)) Then
      Call storesym_initsspair(storage,spc1,nmoles1,spc2,nmoles2,0,level, &
          nderivs,.False.)

    Else If (Associated(params%sum)) Then
      intra_also = sssum_intraon(params%sum)
!      Call storesym_initsspair(storage,spc1,nmoles1,spc2,nmoles2,1,level, &
!          nderivs,intra_also)
      Call storesym_initsspair(storage,spc1,nmoles1,spc2,nmoles2,0,level, &
          nderivs,.False.)

    Else If (Associated(params%xtrnal)) Then
      If (level == 4) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' external evaluations not available with ATM detail'
        Stop
      End If 
      If ((level == 3).And.((nmoles1 > 1).Or.(nmoles2 > 1))) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' external evaluations with more than 1 molecule and MOL detail'
      End If 
      Call storesym_initsspair(storage,spc1,nmoles1,spc2,nmoles2,1,level, &
          nderivs,.False.)

    Else If (params%off) Then
      !** No storage needed, terminate at the species level
      Call storesym_terminate(storage,spc1,spc2,0)

    Else 
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' Not equipped to handle anything but map and basic'
      Stop
    End If

  End Subroutine ssdriver_initstore

  !------------------------------------------------------------------------
  ! Links the first of the SpcSpc_Params structures to the second.  
  ! Useful for making the spc-spc array symmetric while saving memory.
  ! Requires:  params1 -- spc-spc forcefield parameters number 1
  !            params2 -- spc-spc forcefield parameters number 2
  !------------------------------------------------------------------------
  Subroutine ssdriver_link(params1,params2)
    Type(SpcSpc_Params), Intent(Out)  :: params1
    Type(SpcSpc_Params), Intent(In)   :: params2

    Call ssdriver_nullify(params1)

    params1%off = params2%off
    If (Associated(params2%basic)) Then
      params1%basic => params2%basic
    End If
    If (Associated(params2%map)) Then
      params1%map => params2%map
    End If
    If (Associated(params2%sum)) Then
      params1%sum => params2%sum
    End If

  End Subroutine ssdriver_link

  !----------------------------------------------------------------------------
  ! Evaluate a single SPECIES-SPECIES interaction.  It calculates the 
  ! interactions of all molecules of "spc1" with all molecules of "spc2".
  ! By default it will avoid double counting, this can be turned off by
  ! using the optional "allint" flag.
  ! Requires:  params -- spc-spc non-Coulombic forcefield parameters
  !            ffout -- symmetric results forcefield output
  !            spc1 -- first species number
  !            spc2 -- second species number
  !            species -- full coordinate storage structure
  !            simcell -- simulation cell information
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !            pcalcflag -- allow brute force map calculations if True
  !            allint -- force calculation of all atomic pairs if True
  !----------------------------------------------------------------------------
  Logical Function ssdriver_ssint(params,ffout,spc1,spc2,species, &
      simcell,fast,pcalcflag,allint)
    Type(SpcSpc_Params), Intent(InOut)            :: params
    Type(Symmetric_Results), Intent(InOut)        :: ffout
    Integer, Intent(In)                           :: spc1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In)   :: species
    Type(SimCell_Params), Intent(In)              :: simcell  
    Logical, Intent(In)                           :: fast,pcalcflag
    Logical, Intent(In), Optional                 :: allint

    Integer          :: s1,s2
    Logical          :: mapflag, ljflag, allflag

    !** Set the default flag values
    mapflag = .True.
    ljflag = .True.
    ssdriver_ssint = .True.  
    If (Present(allint)) Then
      allflag = allint
    Else
      allflag = .True.
    End If

    !** reorder species numbers such that spc1<=spc2
    s1 = spc1
    s2 = spc2
    If (s1 > s2) Call swap(s1,s2)

    !** Get the species-species interaction from pair interactions
    If (Associated(params%basic)) Then
      ssdriver_ssint = ssbasic_ssint(params%basic,ffout, &
          s1,s2,species,simcell,fast,allflag)
      If (.Not. ssdriver_ssint) Return
    End If

    !** Get the species - fixed species interactions from map interpolation
    If (ssmap_eval(params%map,fast)) Then
      ssdriver_ssint = ssmap_ssint(params%map,ffout%ab(s1,s2),s1,s2, &
          species,simcell,pcalcflag,params%coul)
      If (.Not. ssdriver_ssint) Return
    End If

    !** Get the species-species interaction from sophisticated summation
    If (Associated(params%sum)) Then
      ssdriver_ssint = sssum_ssint(params%sum,species,simcell, &
          s1,s2,fast,ffout)
      If (.Not. ssdriver_ssint) Return
    End If

    !** Get the species-species interactions from external evaluations
    If (Associated(params%xtrnal)) Then
      ssdriver_ssint = xternal_int(params%xtrnal,ffout,species,simcell, &
          (/s1,0,0/),(/s2,0,0/))
      If (.Not. ssdriver_ssint) Return
    End If

  End Function ssdriver_ssint

  !----------------------------------------------------------------------------
  ! Evaluate a single MOLECULE-SPECIES interaction.  It calculates the 
  ! interactions of "molec1" of "spc1" with all molecules of "spc2".
  ! By default it will avoid double counting, this can be turned off by
  ! using the optional "allint" flag.
  ! Requires:  params -- spc-spc non-Coulombic forcefield parameters
  !            ffout -- forcefield output 
  !            spc1 -- first species number
  !            mol1 -- first molecule number
  !            spc2 -- second species number
  !            species -- full coordinate storage structure
  !            simcell -- simulation cell information
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !            pcalcflag -- allow brute force map calculations if True
  !----------------------------------------------------------------------------
  Logical Function ssdriver_msint(params,ffout,spc1,mol1,spc2,species, &
      simcell,fast,pcalcflag)
    Type(SpcSpc_Params), Intent(InOut)          :: params
    Type(Symmetric_Results), Intent(InOut)      :: ffout
    Integer, Intent(In)                         :: spc1,mol1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(SimCell_Params), Intent(In)            :: simcell  
    Logical, Intent(In)                         :: fast,pcalcflag

    Logical                          :: mapflag, ljflag
    Integer                          :: branch
    Integer, Dimension(2)            :: indx
    Type(EnergyPlus), Pointer        :: nrgptr
    Type(Store_Level_Pair), Dimension(:), Pointer  :: abptr

    !** Set the default flag values
    mapflag = .True.
    ljflag = .True.
    ssdriver_msint = .True.

    !** Get the molecule-species interaction from pair interactions
    If (Associated(params%basic)) Then
      ssdriver_msint = ssbasic_msint(params%basic,ffout, &
          spc1,mol1,spc2,species,simcell,fast,.False.)
      If (.Not. ssdriver_msint) Return
    End If

    !** Get the molecule - fixed species interactions from map interpolation
    If (ssmap_eval(params%map,fast)) Then
      Call storesym_ptrs(ffout,(/spc1,mol1,0/),(/spc2,0,0/),branch, &
          indx,nrgptr,abptr)
      If (branch == -2) Then
        ssdriver_msint = ssmap_msint(params%map,ffout%ab(spc1,spc2),spc1, &
            mol1,spc2,species,simcell,pcalcflag,params%coul)        
      Else
        ssdriver_msint = ssmap_msint(params%map,abptr(indx(1)),spc1,mol1,spc2, &
            species,simcell,pcalcflag,params%coul)
      End If
      If (.Not. ssdriver_msint) Return
    End If

    !** Get the molecule-species interaction from sophisticated summation
    If (Associated(params%sum)) Then
      ssdriver_msint = sssum_ssint(params%sum,species,simcell, &
          spc1,spc2,fast,ffout)
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'ssdriver_msint is not implemented'
       Write(*,*) 'on the fly ewald sums are implemented for md type simulations only'
         Write(*,*) 'for which sssum_ssint call from ssdriver_ssint is made'
      Stop
    End If

    !** Get the molecule - species interactions from external evaluations
    If (Associated(params%xtrnal)) Then
      ssdriver_msint = xternal_int(params%xtrnal,ffout,species,simcell, &
          (/spc1,mol1,0/),(/spc2,0,0/))
      Call checkandstop(ssdriver_msint,__FILE__,__LINE__, &
          ' problem with EXTERNAL forcefield calculation for spc1,mol1,spc2:', &
          (/spc1,mol1,spc2/))
      If (.Not. ssdriver_msint) Return
    End If

  End Function ssdriver_msint

  !----------------------------------------------------------------------------
  ! Evaluate a single ATOM-SPECIES interaction.  It calculates the 
  ! interactions of "atom1", "molec1", "spc1" with all molecules of "spc2".
  ! By default it will avoid double counting, this can be turned off by
  ! using the optional "allint" flag.
  ! NOTE: for now I've hacked this to only return results in the EnergyPlus
  !       structure.  Perhaps later we'll switch to using the full structure.
  ! Requires:  params -- spc-spc non-Coulombic forcefield parameters
  !            ffout -- forcefield output 
  !            spc1 -- first species number
  !            mol1 -- first molecule number
  !            atm1 -- first atom number
  !            spc2 -- second species number
  !            species -- full coordinate storage structure
  !            simcell -- simulation cell information
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !            pcalcflag -- allow brute force map calculations if True
  !----------------------------------------------------------------------------
  Logical Function ssdriver_asint(params,ffout,spc1,mol1,atm1,spc2, &
      species,simcell,fast,pcalcflag,hot)
    Type(SpcSpc_Params), Intent(In)               :: params
    Type(EnergyPlus), Intent(InOut)               :: ffout
    Integer, Intent(In)                           :: spc1,mol1,atm1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In)   :: species
    Type(SimCell_Params), Intent(In)              :: simcell  
    Logical, Intent(In)                           :: fast,pcalcflag
    Real(Kind=RDbl), Dimension(:), Optional :: hot

    Logical                    :: mapflag, ljflag
    Real(Kind=RDbl)            :: q

    !** Set the default flag values
    mapflag = .True.
    ljflag = .True.
    ssdriver_asint = .True.

    !** Get the molecule-species interaction from pair interactions
    If (Associated(params%basic)) Then
             If (.Not.Present(hot)) Then
             ssdriver_asint = ssbasic_asint(params%basic,ffout, &
             spc1,mol1,atm1,spc2,species,simcell,fast)
             Else
             ssdriver_asint = ssbasic_asint(params%basic,ffout, &
            spc1,mol1,atm1,spc2,species,simcell,fast,hot)
             End If
      If (.Not. ssdriver_asint) Return
    End If

    !** Get the atom - fixed species interactions from map interpolation
    If (ssmap_eval(params%map,fast)) Then
      If (params%coul) Then
        Call config_getq(species,spc1,mol1,atm1,q)
        ssdriver_asint = ssmap_asint(params%map,ffout,spc1,mol1,atm1,spc2, &
            species,simcell,pcalcflag,q)
        If (.Not. ssdriver_asint) Return
      Else
        ssdriver_asint = ssmap_asint(params%map,ffout,spc1,mol1,atm1,spc2, &
            species,simcell,pcalcflag)
      End If
    End If

    !** Get the molecule-species interaction from sophisticated summation
    If (Associated(params%sum)) Then
             If (.Not.Present(hot)) Then
             ssdriver_asint = sssum_asintHOT(params%sum,spc1,mol1,atm1,spc2,species,simcell, &
             fast,ffout)
             Else
             ssdriver_asint = sssum_asintHOT(params%sum,spc1,mol1,atm1,spc2,species,simcell, &
             fast,ffout,hot)
             End If
      If (.Not. ssdriver_asint) Return
    End If


!!$L      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!!$L      Write(*,*) 'sssum_asint non-existant!!'
!!$L      Stop
  

    !** Get the molecule - species interactions from external evaluations
    If (Associated(params%xtrnal)) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'correct routine in xternal does not exist, must stop'
      Stop
      Call checkandstop(ssdriver_asint,__FILE__,__LINE__, &
          ' problem with EXTERNAL forcefield calculation for '// &
          'spc1,mol1,atm1,spc2:', (/spc1,mol1,atm1,spc2/))
      If (.Not. ssdriver_asint) Return
    End If

  End Function ssdriver_asint

  !----------------------------------------------------------------------------
  ! Check if a particular spc-spc interaction is turned on
  ! Requires:  params -- spc-spc non-Coulombic forcefield parameters
  !----------------------------------------------------------------------------
  Logical Function ssdriver_ison(params)
    Type(SpcSpc_Params), Intent(In)               :: params

    ssdriver_ison = (.Not. params%off)

  End Function ssdriver_ison

  !----------------------------------------------------------------------------
  ! Returns information about the forcefield parameters governing the 
  ! interactions between two subsets of the system.  This information can
  ! then be used elsewhere.
  ! Requires:  ssparams -- spc-spc forcefield parameters structure
  !            subset1 -- 1st subset 
  !            subset2 -- 2nd subset
  !            nsets -- number of potential sets returned
  !            list -- list of atom numbers for interaction (set,1:Natoms)
  !            params -- array of strings containing type and parameters
  !----------------------------------------------------------------------------
  Subroutine ssdriver_listparams(ssparams,subset1,subset2,nsets,list,params)
    Type(SpcSpc_Params), Intent(In)                   :: ssparams
    Integer, Dimension(:), Intent(In)                 :: subset1,subset2
    Integer, Intent(Out)                              :: nsets
    Integer, Dimension(:,:), Intent(Out)              :: list
    Character(len=lstrLen), Dimension(:), Intent(Out) :: params

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Not ready yet'
    Stop

  End Subroutine ssdriver_listparams

  !----------------------------------------------------------------------------
  ! Returns the box increments for a potential map if they are available.
  ! Requires:  ssparams -- spc-spc forcefield parameters structure
  !            boxincrements -- x,y,z increments for potential map
  !            boxsteps -- number of boxes in x,y,z directions
  !----------------------------------------------------------------------------
  Logical Function ssdriver_boxinfo(ssparams,boxincrements,boxsteps)
    Type(SpcSpc_Params), Intent(In)                 :: ssparams
    Real(Kind=RDbl), Dimension(3), Intent(Out)      :: boxincrements
    Integer, Dimension(3), Intent(Out)              :: boxsteps

    ssdriver_boxinfo = .False.

    If (Associated(ssparams%map)) Then
      ssdriver_boxinfo = ssmap_boxinfo(ssparams%map,boxincrements,boxsteps)
    End If

  End Function ssdriver_boxinfo

  !----------------------------------------------------------------------------
  ! Displays the ssdriverparams structure for one species pair
  ! Requires:  params -- noncoulombic spc-spc parameters
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine ssdriver_display(params,spc1,spc2,indent,unitno)
    Type(SpcSpc_Params), Intent(In)   :: params
    Integer, Intent(In)               :: spc1,spc2,indent
    Integer, Optional, Intent(In)     :: unitno

    Integer                           :: nassoc,unit
    Character(len=strLen)             :: name1,name2
    Character(len=indent)             :: blank

    blank = Repeat(' ',indent)

    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    nassoc = 0
    name1 = molecules_name(spc1)
    name2 = molecules_name(spc2)
    Write(unit,'(a,2(a,i1),a,4a,t35)', Advance='No') blank, &
        '(',spc1,',',spc2,')  ', Trim(name1),', ',Trim(name2), ': '

    If (Associated(params%basic)) Then
      nassoc = nassoc + 1
      Write(unit, '(2a)') blank,Trim(ssbasic_idstring)
      Call ssbasic_display(params%basic,indent+2,unit)
    End If

    If (Associated(params%map)) Then
      nassoc = nassoc + 1
      Write(unit, '(2a)') blank,Trim(ssmap_idstring)
      Call ssmap_display(params%map,indent+2,unit)
    End If

    If (Associated(params%sum)) Then
      nassoc = nassoc + 1
      Write(unit, '(2a)') blank,Trim(sssum_idstring)
      Call sssum_display(params%sum,indent+2,unit)
    End If

    If (Associated(params%xtrnal)) Then
      nassoc = nassoc + 1
      Write(unit, '(2a)') blank,Trim(EXTERNAL_KEY)
      Call xternal_display(params%xtrnal,indent+2,unit)
    End If

    If (params%off) Then
      Write(unit, '(2a)') blank,'OFF'
    End If
    
!LC    Write(unit,'(2a,i3)') blank,'Number of types associated: ',nassoc

  End Subroutine ssdriver_display

  !--------------------------------------------------------------------
  ! Cleans the species-species non-Coulombic parameters structure
  ! Requires:  params -- non-Coulombic parameters structure
  !--------------------------------------------------------------------
  Subroutine ssdriver_clean(params)
    Type(SpcSpc_Params), Intent(InOut)     :: params

    Integer        :: error

    If (Associated(params%basic)) Then
      Call ssbasic_clean(params%basic)
      Deallocate(params%basic, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'ssparams%basic')
    End If

    If (Associated(params%map)) Then
      Call ssmap_clean(params%map)
      Deallocate(params%map, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'ssparams%map')
    End If

    If (Associated(params%map)) Then
      Call sssum_clean(params%sum)
      Deallocate(params%sum, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'ssparams%sum')
    End If

  End Subroutine ssdriver_clean


!------------------------------------------------------------------------------
! Tina added
! Subroutine to get potential parameters, the array pot_params contains
! A, B, C, D, hicut, locut and returns zero if the parameters are not
! defined for the potential 
!------------------------------------------------------------------------------
 Subroutine ssdriver_getpotparameters(params, a1, a2, pot_params)
   Type(spcspc_params), Intent(in)               :: params
   Integer, Intent(IN)                          :: a1, a2
   Real(kind = Rdbl), Dimension(6),INTENT(OUT)  :: pot_params
 
   IF(ASSOCIATED(params%basic)) THEN
      Call ssbasic_getpotparameters(params%basic,a1,a2,pot_params)
   ELSE IF(ASSOCIATED(params%sum)) THEN
      pot_params = 0._RDbl
   ELSE
      Write(*,*) __FILE__,__LINE__
      Write(*,*) 'This routine should not be called!!!'
   END IF
 End Subroutine ssdriver_getpotparameters 


End Module ssdriver






