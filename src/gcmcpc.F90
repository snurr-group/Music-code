!-------------------------------------------------------------
! This module contains the different data structures and 
!routines for doing the analysis of results of a GCMC Run,
! ie. The post code stuff.
!-------------------------------------------------------------
Module gcmcpc

  Use defaults,    Only : strLen, RDbl, RSgl, MAX_SORBS, MAX_ATOMS, &
      MAX_SITES, MAX_DAT_FILES, d_ctrl_file, zero, one, INITIAL_SIZE
  Use file,        Only :  file_open, file_getunit, file_settag
  Use histogram,   Only :  Histogram_Params, histogram_init, histogram_display
  Use multihist,   Only :  MultiHist_Params, multihist_init, &
      multihist_normalize, multihist_getnbins, multihist_update
  Use stats,       Only :  Statistics,stats_init,stats_update, &
      stats_getblock, stats_getcavg, stats_getstd
  Use general,     Only :  general_getcontenttag, genparams, &
      general_getnoofconfigs, general_init, general_display
  Use config,      Only :  AtMolCoords, config_defaultinit, config_display, &
      config_getnatoms, config_initdefaultsize, config_allocfields, &
      config_checkandincr, config_init, config_isfixed
  Use datafile,    Only :  CONFILE, datafile_initin, datafile_gencontrolfile, &
      datafile_readconfig
  Use forcefield,  Only :  forcefield_init, forcefield_display  
  Use smap,        Only :  Smap_Params, smap_init, smap_getPosSiteType
  Use molecules,   Only : molecules_getnsorbs, molecules_gettype, &
      molecules_getnatoms, molecules_name, molecules_init, molecules_display, &
      molecs 
  Use molecule,    Only : MolecularParams, molecule_getname 
  Use atom,        Only : atom_init, atom_display
  Use simcell,     Only :  SimCell_Params, simcell_init, simcell_display, &
      simcell_maptouc
  Use utils,       Only : filesrchstr, genfilename, stripcmnt, split, toint, &
      toreal, maxn, toupper, getcom, isfileopen, allocErrDisplay
  Use vector,      Only : VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag

  Implicit None
  Save

  Private
  Public :: GcMCPCParams, GCMCPCVars, GCMCPCdat, postgcmc, gcmcpc_init, &
      gcmcpc_GetNoOfFiles, gcmcpc_initfilenames, gcmcpc_openIsoThermFiles, &
      gcmcpc_reinitstats, gcmcpc_getaverages
  
  !** The tag marking the beginning of the gcmcpc section
  Character(len=strLen), Parameter    :: default_gcmcpc_tag = &
      "GCMC Post Code Information"
  ! The tag marking the beginning of the GCMC section
  Character(len=strLen), Parameter    :: default_gcmc_tag = &
      "GCMC Information"  
  
  ! stuff to be read from postcode-control file
  Type GCMCPCVars
    Integer                       :: N0 
    Integer                       :: NN
    Character(len=strLen)         :: new_ctrlfile
    Character(len=strLen)         :: basename
    Integer                       :: nblocks
    Real(kind=RDbl)       :: startskip !** % of initial data to be skipped
    Real(kind=RDbl)       :: endskip   !** % of final data to be skipped
  End Type GCMCPCVars

  ! variables which pertain to averaging of one datafile
  Type GCMCPCdat
    !variables for averaging among all the data in the datafile
    Type(Statistics),Dimension(MAX_SORBS)::nmole_list
    Type(Statistics) ::energy

    !variables for finding standard devn and averages among the blocks
    Type(Statistics),Dimension(MAX_SORBS)::nmoleblock
    Type(Statistics) ::nrgblock
  End Type GCMCPCdat

  Type GCMCPcParams
    Type(GCMCPCVars) :: Var
    Type(GCMCPCdat)  :: dat
    
    Type(AtMolCoords), Dimension(:), Pointer :: sorbates
    Type(SimCell_Params)            :: scell 
    Character(len=strLen)         :: gcmc_tag
    Integer :: nsims,content_tag
    Integer,Dimension(MAX_SORBS):: molec_natoms
    Character(len=strlen),Dimension(MAX_SORBS):: molecnames
    Integer:: nads,nsorbates
    Character(len=strlen):: commentString
    
    !** List of fugacities (pressures) for each sorbate
    Real(kind=RDbl),Dimension(:,:),Pointer :: fuglists

    !list of datafiles to be processed
    Type(CONFILE),Dimension(MAX_DAT_FILES) :: dflist

    !* contains name of all output files for loading averages
    Character(len=strlen),Dimension(MAX_DAT_FILES,MAX_SORBS)   :: outfile 

    !* contains name of all output files for energy averages
    Character(len=strlen),Dimension(MAX_DAT_FILES)             :: en_file
    Character(len=2*strlen),Dimension(MAX_SORBS)             :: IsoThermLoad
    Character(len=2*strlen)             :: IsoThermNrg

  End Type GCMCPcParams
  
  Type(GCMCPcParams):: postgcmc
  
Contains
  !----------------------------------------------------------------
  !This subroutine does the necessary initialisations
  !for reading what is written in the file - config_filename
  !It reads the header, generates the original control file that 
  !was used for making the given  config file,and initialises all structures
  !---------------------------------------------------------------
  Subroutine gcmcpc_init(postfile)
    Character(len=strlen),Intent(in)::postfile
    Type(CONFILE)                 :: firstconfig
    Integer :: unitno
    postgcmc%gcmc_tag= "GCMC Post Code Information"

    !** read info from post code-control file
    Call gcmcpc_readPostFile(postfile,postgcmc%gcmc_tag)
    
    !*** generates the New Control File from the given config file
    firstconfig%name=genfilename(postgcmc%Var%basename, postgcmc%Var%N0)
    Write(*,*) 
    Write(*,*) "Reading the first config file("//Trim(firstconfig%name)//")"
    Write(*,*) " -and re-generating the control-file used for the simulation"
    Call datafile_initin(firstconfig,0)
    Call datafile_gencontrolfile(postgcmc%Var%new_ctrlfile,&
        firstconfig)
    Write(*,*)
    Call file_settag(postgcmc%Var%new_ctrlfile, d_ctrl_file)
    Close(firstconfig%unitno)
    
    !** initialises general,atom,molecs using new control file
    Write(*,*) 
    Write(*,*) "Reading from the re-generated control file"
    Write(*,*) 
    Call general_init(postgcmc%Var%new_ctrlfile)
    Write(*,*) "##################################################"
    Write(*,*) "THE GENERAL PARAMS USED FOR THE SIMULATION"
    Call general_display(6)
    write(*,*)
    Call atom_init()
    Write(*,*) "##################################################"
    Write(*,*) "THE ATOMS USED FOR THE SIMULATION"
    Call atom_display(6)

    Write(*,*) 
    Call molecules_init()
    Write(*,*) "##################################################"
    Write(*,*) "THE MOLECULES USED FOR THE SIMULATION"
    Call molecules_display(6)

    Write(*,*) 
    Call simcell_init(postgcmc%scell, postgcmc%Var%new_ctrlfile)
    Write(*,*) "##################################################"
    Write(*,*) "THE SIMULATION-CELL VARIABLES USED FOR THE SIMULATION"
    Call simcell_display(postgcmc%scell,2,6)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) 
    Call config_init(postgcmc%sorbates, postgcmc%scell, &
        postgcmc%Var%new_ctrlfile)
    Write(*,*) "##################################################"
    Write(*,*) "THE INITIAL CONFIGURATION USED FOR THE SIMULATION"
    Call config_display(postgcmc%sorbates)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Call forcefield_init(postgcmc%Var%new_ctrlfile, postgcmc%scell)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Call forcefield_display(postgcmc%sorbates, postgcmc%scell)
    Write(*,*) 
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__    

    Call gcmcpc_getMolecDetails(postgcmc%nsorbates,postgcmc%nads &
        ,postgcmc%molecnames,postgcmc%molec_natoms)

    !** Read the ctrl file section where gcmc was initialised
    Call gcmcpc_readgcmcsection(postgcmc%Var%new_ctrlfile)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
!    Call gcmc_initdisplay(postgcmc%gcmcparams,6)
!    Call gcmc_displaysimparams(postgcmc%gcmcparams,postgcmc%Var%N0,3)
!    postgcmc%nsims = gcmcpc_getnosims(postgcmc%gcmcparams)
    postgcmc%content_tag=general_getContentTag()



  End Subroutine gcmcpc_init
  
  !-----------------------------------------------------------------
  ! Reads the section containing specifics of the gcmc Post code,
  !and initialises the "var" part of postgcmc
  !Also sets the array dflist, which contains the name of files to be analysed
  !-----------------------------------------------------------------
  Subroutine gcmcpc_readPostFile(postfile,opt_tag)
    Character(*), Intent(in)  :: postfile
    Character(*), Optional, Intent(in)  :: opt_tag
    
    Character(len=strlen)::basename
    Character(len=strLen)   :: tag, line
    Integer                 :: unitno, lineno,i,nfiles
    
    !** Open the postfile if it is not opened
    unitno = isfileopen(postfile)
    If (unitno < 0) Then
      unitno = file_getunit(postfile)
      Open(file=postfile, unit=unitno)
    Endif
    
    If (Present(opt_tag)) Then
      tag = opt_tag
    Else
      tag = default_gcmcpc_tag
    End If
    !** Find the GCMC post code tag
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      !      Stop
    Endif
    
    !** Read the necessary stuff, and generate the list of config files    
    Write(*,*)
    Write(*,*) "Reading the post code control file"
    Read(unitno,*) postgcmc%Var%N0
    Read(unitno,*) postgcmc%Var%NN
    Read(unitno,*) postgcmc%Var%new_ctrlfile
    Read(unitno,*) postgcmc%Var%basename
    Read(unitno,*) postgcmc%Var%nblocks
    Read(unitno,*) postgcmc%Var%startskip
    If ((postgcmc%Var%startskip>100).Or.(postgcmc%Var%startskip<0)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,&
          "      Wrong value of startskip specified in control file"
      Stop
    Endif

    Read(unitno,*) postgcmc%Var%endskip
    If ((postgcmc%Var%endskip>100).Or.(postgcmc%Var%endskip<0).Or. &
        &((postgcmc%Var%endskip+postgcmc%Var%startskip)>100)) &
        Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,&
          "      Wrong value of endskip specified in control file"
      Stop
    Endif

    !** commentString will be the prefix for all comments in the output files
    Read(unitno,*) postgcmc%commentString
    postgcmc%commentString=Trim(stripcmnt(postgcmc%commentString))//"#"
    basename=postgcmc%Var%basename
    nfiles=gcmcpc_GetNoOfFiles()

    !*** Fixes the filename of config files to be processed
    Do i=1,nfiles
      If (nfiles > MAX_DAT_FILES) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Maximum no. of config files exceeded "
        Stop
      End If
      postgcmc%dflist(i)%name=genfilename(basename, postgcmc%Var%N0+i-1)
    End Do
    Close(unitno)

  End Subroutine gcmcpc_readPostFile
  
  
  !------------------------------------------
  ! Initialises the arrays and filenames for averaging varying quantities
  ! of the run, outfile contains the name of output files for each 
  ! config file and each molecule
  !------------------------------------------
  Subroutine gcmcpc_initfilenames()
    Character(len=strlen)                     ::tempstrg
    Integer :: i,j,k,nfiles
    
!    Call gcmcpc_getMolecDetails(postgcmc%nsorbates,postgcmc%nads &
!        ,postgcmc%molecnames,postgcmc%molec_natoms)
    nfiles=gcmcpc_GetNoOfFiles()
    Do i=1,nfiles
      Do j=1,postgcmc%nsorbates
        !only molecules with more than 0 atoms are counted. ie things like
        ! zeolites are excluded
        If (config_isfixed(postgcmc%sorbates(j))) Cycle
        If (postgcmc%molec_natoms(j)>0) Then
          tempstrg=postgcmc%dflist(i)%name
          postgcmc%outfile(i,j)= &
              Trim(tempstrg)//"."//Trim(postgcmc%molecnames(j))
        Endif
      End Do

    !*** names of the energy file for the i th config file
    tempstrg=genfilename(postgcmc%Var%basename,postgcmc%Var%N0+i-1)
    tempstrg=Trim(tempstrg)//".energy"
    postgcmc%en_file(i)=tempstrg

    End Do

    !*** names of the  file for the Overall Isotherms
      Do j=1,postgcmc%nsorbates
        !only molecules with more than 0 atoms are counted. ie things like
        ! zeolites are excluded
        If (postgcmc%molec_natoms(j)>0) Then
          tempstrg=Trim(postgcmc%Var%basename)//".IsoThermLoad"
          postgcmc%IsoThermLoad(j)= &
              Trim(tempstrg)//"."//Trim(postgcmc%molecnames(j))
        Endif
      End Do
    postgcmc%IsoThermNrg=Trim(postgcmc%Var%basename)//".IsoThermNrg"
    

  End Subroutine gcmcpc_initfilenames
  
  !------------------------------------------------------------------
  ! Reinitialises or Initialises all the stats variable, before
  ! starting to read a new config file
  !------------------------------------------------------------------  
  Subroutine gcmcpc_reInitStats()
    Integer::j,blocksize
    
    Call gcmcpc_getblocksize(blocksize)
    Do j=1,postgcmc%nsorbates
      Call stats_init(postgcmc%dat%nmole_list(j),' ',blocksize,'f15.3',.False.)
      Call stats_init(postgcmc%dat%nmoleblock(j),' ',postgcmc%Var%nblocks, &
          .False.,'f15.3')

    End Do
    Call stats_init(postgcmc%dat%energy,' ',blocksize,.False.,'e15.4')
    Call stats_init(postgcmc%dat%nrgblock,' ',postgcmc%Var%nblocks,&
        .False.,'e15.4')
  End Subroutine gcmcpc_reInitStats
  
  
  !------------------------------------------------------------------
  ! Gets the names of molecules, and number of atoms in each of them
  !------------------------------------------------------------------  
  Subroutine gcmcpc_getMolecDetails(nsorbates,nads,molecnames,molec_natoms)
    Integer,Intent(out):: nads
    Integer,Intent(out):: nsorbates
    Character(len=strlen),Dimension(:),Intent(inout)  :: molecnames
    Integer,Dimension(:),Intent(inout)  :: molec_natoms
    Integer :: i
    
    nsorbates = molecules_getnsorbs()
    nads=nsorbates
    Do i=1,nsorbates
      molecnames(i)=molecules_name(i)
      molec_natoms(i)=config_getnatoms(postgcmc%sorbates,i)
      If (molec_natoms(i)==0) nads=nads-1
    End Do
  End Subroutine gcmcpc_getMolecDetails
  
  !------------------------------------------
  ! Gives the number of files to be processed
  !------------------------------------------
  Integer Function gcmcpc_GetNoOfFiles()
    gcmcpc_GetNoOfFiles=(postgcmc%Var%NN-postgcmc%Var%N0)+1
  End Function gcmcpc_GetNoOfFiles

  !---------------------------------------------------------
  !gets the name of configfile to be opened
  !------------------------------------------------------------
  Subroutine gcmcpc_getconfigfile(sim_n,configfile)
    Type(confile), Intent(inout)       :: configfile
    Integer,Intent(in)::sim_n
    configfile%name=postgcmc%dflist(sim_n)%name
  End Subroutine gcmcpc_getconfigfile

  !---------------------------------------------------------
  !gets the fugacity used for generating the sim_n th config file
  ! for the given sorbate : sorb
  !------------------------------------------------------------
  Real(kind=RDbl) Function gcmcpc_getFugacity(sim_n,sorb)
    Integer,Intent(in)::sim_n,sorb
    Character(len=strLen) :: molecname
    gcmcpc_getFugacity=postgcmc%fuglists(sorb,sim_n)
  End Function gcmcpc_getFugacity

  !---------------------------------------------------------
  !gets the blocksize
  !------------------------------------------------------------
  Subroutine gcmcpc_getblocksize(blocksize)
    Integer,Intent(out)::blocksize
    Integer            ::ndats,nstart,nend,nblocks

    ndats=general_getnoofconfigs()
    nstart=Int(ndats*(postgcmc%Var%startskip/100))
    nend=ndats-Int(ndats*(postgcmc%Var%endskip/100))
    nblocks=postgcmc%Var%nblocks
    blocksize=Int((nend-nstart)/nblocks)
  End Subroutine gcmcpc_getblocksize


  !--------------------------------------------------------------------
  !Opens the output files and returns the unit numbers
  ! each sim_n denotes a config file
  !--------------------------------------------------------------------
  Subroutine gcmcpc_openoutfiles(sim_n,newunits,en_unit,prefx)  
    Integer,Intent(in) :: sim_n
    Integer,Dimension(:),Intent(out) :: newunits
    Integer,Intent(out) :: en_unit
    Character(*),Intent(in),Optional:: prefx
    Real (kind=RDbl)    :: pressure

    Character(len=strlen):: comnt1
    Character(len=strlen):: molecname
    Character(len=5):: temp_str
    Integer:: i
    
    If (Present(prefx)) Then 
      comnt1=Trim(prefx)
    Else 
      comnt1="#"
    Endif

    !** opens files for molecules with more than 0 atoms.
    Do i=1,postgcmc%nsorbates
      pressure = gcmcpc_getFugacity(sim_n,i)
      If (config_isfixed(postgcmc%sorbates(i))) Cycle
      If (postgcmc%molec_natoms(i)/=0) Then
        newunits(i)=file_getunit(postgcmc%outfile(sim_n,i))
        Open(file=postgcmc%outfile(sim_n,i),unit=newunits(i))
        Write(0,*) __FILE__,__LINE__, "Opening output file ", &
            postgcmc%outfile(sim_n,i)
        Write(newunits(i),'(a)') Trim(comnt1)//"###############################"
        Write(newunits(i),'(2a)') Trim(comnt1)//"GCMC Loading: &
            &Number of Molecules vs iteration number For ",&
            Trim(postgcmc%molecnames(i))
        Write(newunits(i),'(a,f15.5)') Trim(comnt1)//"Pressure  : ",pressure
        Write(newunits(i),'(a)') Trim(comnt1)//&
            "   Iter.No          Inst.  Block.Average  Cumul.Average         ST.Dev"
      End If
    End Do
    
    !*** opens the energy file for the sim_n th config file
    en_unit=file_getunit(postgcmc%en_file(sim_n))
    Open(file=postgcmc%en_file(sim_n),unit=en_unit)
    Write(en_unit,'(a)') Trim(comnt1)//"##################################"
    Write(en_unit,'(a)') Trim(comnt1)//"GCMC data :Energy vs &
        &Iteration Number "
    Write(en_unit,'(a)') Trim(comnt1)//&
    "   Iter.No          Inst.  Block.Average  Cumul.Average         ST.Dev"

    
  End Subroutine gcmcpc_openoutfiles
  
  !-----------------------------------------
  ! Writes the averaged data( no of molecs) into the files, so that 
  ! they can be plotted and viewed
  !-----------------------------------------
  Subroutine  gcmcpc_writeoutfiles(units,en_unit,itn)
    Integer,Dimension(MAX_SORBS),Intent(in):: units
    Integer,Intent(in):: en_unit
    Integer,Intent(in):: itn

    Integer :: i

    Do i=1,postgcmc%nsorbates
      If (postgcmc%molec_natoms(i)>0) Then
        Write(units(i),'(i10,4f15.3)') itn,postgcmc%dat%nmole_list(i)%inst, &
            postgcmc%dat%nmole_list(i)%blockavg, &
            postgcmc%dat%nmole_list(i)%cumavg, &
            Sqrt(postgcmc%dat%nmole_list(i)%cumdevsq)
        
      Endif
    End Do
    Write(en_unit,'(i10,4f15.3)') itn,postgcmc%dat%energy%inst, &
        postgcmc%dat%energy%blockavg, &
        postgcmc%dat%energy%cumavg, &
        Sqrt(postgcmc%dat%energy%cumdevsq)
    
  End Subroutine gcmcpc_writeoutfiles
  
  !-----------------------------------------
  !Closes the output files corresponding to sim_n
  !-----------------------------------------
  Subroutine gcmcpc_closeoutfiles(newunits,en_unit)  
    Integer,Intent(in) :: en_unit
    Integer,Dimension(:),Intent(in) :: newunits
    
    Integer :: i
    Do i=1,postgcmc%nsorbates
      If (postgcmc%molec_natoms(i)>0) Close(unit=newunits(i))
    End Do
    Close(en_unit)
  End Subroutine gcmcpc_closeoutfiles

  
  !-----------------------------------------
  !Closes the output files corresponding to sim_n
  !-----------------------------------------
  Subroutine gcmcpc_statsupdate(nmole_list,energy)  
    Integer,Dimension(:),Intent(in) :: nmole_list
    Real(kind=RDbl),Intent(in) :: energy
    
    Integer :: i
    
    Do i=1,postgcmc%nsorbates
      If (postgcmc%molec_natoms(i)>0) Then
        Call stats_update(postgcmc%dat%nmole_list(i),nmole_list(i)*1.0_RDbl)
      Endif
  End Do
  
  Call stats_update(postgcmc%dat%energy,energy)
    
  End Subroutine gcmcpc_statsupdate
  
  !---------------------------------------------------------
  !Opens the Isotherm files, writes their headers
  !------------------------------------------------------------
  Subroutine gcmcpc_openIsothermFiles()  
    Character(len=strlen):: comnt1
    Character(len=5):: temp_str
    Integer :: i,en_unit
    Integer,Dimension(MAX_SORBS)::newunits

      comnt1=postgcmc%commentString

    !** opens files for molecules with more than 0 atoms.
    Do i=1,postgcmc%nsorbates
      !only molecules with more than 0 atoms are counted. ie things like
      ! zeolites are excluded
      If (postgcmc%molec_natoms(i)>0) Then
        newunits(i)=file_getunit(postgcmc%IsoThermLoad(i))
        Open(file=postgcmc%IsoThermLoad(i),unit=newunits(i))
        Write(newunits(i),'(a)') Trim(comnt1)//"##############################"
        Write(newunits(i),'(2a)') Trim(comnt1)//"Isotherm From GCMC: &
            &Number of Molecules vs Fugacity ",Trim(postgcmc%molecnames(i))
        Write(newunits(i),'(a)') Trim(comnt1)//&
            "    Sim. No  Pressure        Loading         Loading/uc     STD(Loading) "
      Endif
    End Do
    !** opens file for writing energies into
    en_unit=file_getunit(postgcmc%IsoThermNrg)
    Open(file=postgcmc%IsoThermNrg,unit=en_unit)
    Write(en_unit,'(a)') Trim(comnt1)//"###############################"
    Write(en_unit,'(a)') Trim(comnt1)//"Energies From GCMC: &
        &Energy vs Fugacity "
    Write(en_unit,'(a)') Trim(comnt1)//&
            "    Sim.No  Total Energy    Energy/Molec    STD(Total) "

  End Subroutine gcmcpc_openIsothermFiles
  
  !---------------------------------------------------------
  !Calculates the averages in a datafile
  !sim_n= serial number of the run among the runs to be processed
  !------------------------------------------------------------
  Subroutine gcmcpc_getaverages(sim_n)
    Integer,Intent(in)           :: sim_n

    Type(confile)                :: configfile
    Integer                      :: blocksize,en_unit,ndats,nstart,nend,j,k!
    Integer                      :: unitno,itnScale
    Integer,Dimension(MAX_SORBS) :: newunits
    Integer,Dimension(MAX_SORBS) :: nmole_list
    Real(kind=RDbl)              :: energy,blockavg,cumnrg,cumnrg_permolec
    Real(kind=RDbl)              :: cavg,cavg_perunitcell,std_amongblocks
    Real(kind=RDbl)              :: cavgsum,pressure
    Real(kind = RDbl),Dimension(MAX_SORBS,MAX_SORBS,2)   :: potlist
    Real(kind = RDbl),Dimension(MAX_ATOMS)  :: intra,ke
    
    Call gcmcpc_getblocksize(blocksize)
    ndats=general_getnoofconfigs()
    !** itn.numbers from which averaging of data should start and end
    !** the values outside this range are neglected
    nstart=ndats*(postgcmc%Var%startskip/100)
    nend=ndats-Int(ndats*(postgcmc%Var%endskip/100))
    Call gcmcpc_getconfigfile(sim_n,configfile)
    
    Write(*,*)
    Write(*,*)"Reading and doing the average for data in configfile - "//&
        Trim(configfile%name)
    Write(*,*)
    
    Call gcmcpc_openoutfiles(sim_n,newunits,en_unit,postgcmc%commentString)
    Call datafile_initin(configfile,1)
    
    !** valuse were written to config file at evry "iconfig" step.
    itnScale=genparams%iconfig
    
    Do j=1,ndats

      Call datafile_readconfig(postgcmc%sorbates,configfile,nmole_list,&
          energy,potlist,intra,ke)
      If ((j>nstart).And.(j<(nend+1))) Then

        !**updates the averages stored in postgcmc, the gcmcparams structure
        Call gcmcpc_statsupdate(nmole_list,energy)
        
        !** write energy/loading to output files
        If (Mod(j-nstart, blocksize) == 0) Then
          Call gcmcpc_writeoutfiles(newunits,en_unit,j*itnScale)
          
          !** Update the deviation stats between the block avg values
          Do k=1,postgcmc%nsorbates
            If (postgcmc%molec_natoms(k)/=0) Then
              blockavg=stats_getblock(postgcmc%dat%nmole_list(k))
              Call stats_update(postgcmc%dat%nmoleblock(k),blockavg)
            Endif
          End Do
          
          blockavg=stats_getblock(postgcmc%dat%energy)
          Call stats_update(postgcmc%dat%nrgblock,blockavg)
        Endif
      Endif
      
    End Do
    
    Call gcmcpc_closeoutfiles(newunits,en_unit)
    
    !** writes the loading values for current fugacity into IsothermFile 
    cavgsum=0
    Do k=1,postgcmc%nsorbates
      !** avoid zeolite type of molecules
      If (postgcmc%molec_natoms(k)/=0) Then
        unitno=file_getunit(postgcmc%IsoThermLoad(k))
        cavg=stats_getcavg(postgcmc%dat%nmoleblock(k))
        cavg_perunitcell=cavg/postgcmc%scell%ncells
        
        ! std of loading average between the block values
        std_amongblocks=    stats_getstd(postgcmc%dat%nmoleblock(k)) 
        
        !total number of molecules/unit cell for all sorbates
        cavgsum=cavgsum+cavg
        pressure= gcmcpc_getFugacity(sim_n,k)
        Write(unitno,'(a,i4,a,f12.5,a,f12.5,a,f12.5,a,f12.5)')  &
            "   ",sim_n,"    ",pressure,"    ",cavg,"    ",&
            cavg_perunitcell,"    ",std_amongblocks
      Endif
    End Do
    
    !** enrgy values to Isotherm file
    unitno=file_getunit(postgcmc%IsoThermNrg)
    cumnrg=stats_getcavg(postgcmc%dat%nrgblock)
    Write(*,*) "cavg sum",cavgsum
    cumnrg_permolec=cumnrg/cavgsum
    
    std_amongblocks= stats_getstd(postgcmc%dat%nrgblock)
    Write(unitno,'(a,i4,a,f12.5,a,f12.5,a,f12.5)')  &
        "   ",sim_n,"    ",cumnrg,"    ",cumnrg_permolec,"    ",&
        std_amongblocks 
  End Subroutine gcmcpc_getaverages
  
  !---------------------------------------------------------
  !Reads the section of the control file that contains
  !GCMC Information
  !------------------------------------------------------------
  Subroutine gcmcpc_readgcmcsection(ctrlfile)
    Character(len=strLen),Intent(in) :: ctrlfile
    Character(len=strLen)::defaults_gcmcpc_teg,line,eostag
    Character(len=strLen)::pressureline,exclusionsites,smapname,sorbname,tag
    Integer :: gcmclineno,no_of_movetypes,i,j,unitno,error,sorbtype
    Real(kind=RDbl)::tk
    Integer :: nsims,niterations,blocksize,nsorbs

    Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    !** Open the ctrl_file if it is not opened
    unitno = file_open(ctrlfile)
    Rewind(unitno)

    !** Find the GCMC section
    gcmclineno = filesrchstr(unitno, default_gcmc_tag, line)
    If (gcmclineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif

    Read(unitno, *) niterations
    Read(unitno, *) tk
    Read(unitno, '(a)') eostag
    eostag = stripcmnt(eostag)
    Read(unitno, *) postgcmc%nsims
    Read(unitno, *) blocksize
    Read(unitno, *) nsorbs
    !** read the blankline at the end
    Read(unitno, *)
    
    Allocate(postgcmc%fuglists(postgcmc%nsorbates,postgcmc%nsims),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"fuglists-in-gcmcpc")
    
    !**I am trying to Read all the gcmcmoves params here, will it work??

    !** reading the first sorbate name
    Read(unitno, *) sorbname
    Do i=1, nsorbs

      !** sorbname for i>1 obtained from reading at the end of the loop 
      sorbtype=molecules_gettype(sorbname)

      !** different pressure values in the isotherm
      Read(unitno,'(a)') pressureline
      Read(unitno, *) smapname
      Write(*,*) smapname

      Call gcmcpc_getpressures(sorbtype,pressureline)      
      Read(unitno,'(a)') exclusionsites

      Write(*,*) "Excl list-Display to be Coded here"
      Read(unitno,*) no_of_movetypes

      Write(*,*) "Move weights, moves_init to be coded here "
      !** Hacking here, Assuming GCMC section will contain only 100 lines MAX

      Do j=1,100
        Read(unitno,*) sorbname
        sorbtype=molecules_gettype(sorbname)

        If(sorbtype>0) Exit

        If (j>99) Then
          Write(*,*) "Something screwed up here, Contact the guy who &
              &wrote this piece of code - Anonymous"
          Stop
        Endif
      End Do

    End Do

  End Subroutine gcmcpc_readgcmcsection
  
  !---------------------------------------------------------
  ! Parses the fugacity line to set the fugacities of the
  ! sorbate type.  It can either do the fugacities at equal
    ! intervals or read them from a file.
  !---------------------------------------------------------
  Subroutine gcmcpc_getpressures(sorbtype, pressureline)
    Integer, Intent(in) :: sorbtype 
    Character(*), Intent(in)  :: pressureline
    Real(kind=RDbl) :: startp, endp
    Integer         :: i, nfields
    Character(len=strLen)     :: strfields(10), filename, str, pressurefile
    
    ! Split the fields by commas
    nfields = split(pressureline, strfields, ",") 
    If (nfields == 1) Then
      Write(*,*)"The Pressure values are being read from the same &
          & file which was usec for the run, Make sure it was not changed"
      ! We need to read the points from a file
      ! Get the filename
      nfields = split(pressureline, strfields, " ")
      pressurefile = strfields(1)
      Write(*,*) "Reading Pressure from :"//Trim(pressurefile)
      Call gcmcpc_readpressures(postgcmc%fuglists(sorbtype,1:postgcmc%nsims), &
          pressurefile,sorbtype)
    Else

      ! We need to generate the range
      startp = toreal(strfields(1))

      ! Get the last field after separating it from the rest of the line
      str = strfields(2)
      nfields = split(str, strfields, " ")
      endp = toreal(strfields(1))

      ! Fill in the gcmcsorbs pressure data
      Call gcmcpc_genpressures(postgcmc%fuglists(sorbtype,1:postgcmc%nsims), &
          startp, endp)
    End If
  End Subroutine gcmcpc_getpressures

  !-----------------------------------------------------------
  ! Generates the actual fugacity list from the given intervals
  ! contains the start fugacity, end fugacity and no. of points
  !--------------------------------------------------------
  Subroutine gcmcpc_genpressures(fuglist, startp, endp)
    Real(kind=RDbl), Dimension(:),Intent(inout)  :: fuglist
    Real(kind=RDbl), Intent(in)                  :: startp, endp
    Real(kind=RDbl)     :: pressincr 
    Integer             :: i, npts

    ! Check to see what kind of a scale we are going to use
    ! A log scale is used if the ending fugacity is greater
    ! than the starting fugacity by 2 orders of magnitude
    npts = postgcmc%nsims
    If (Abs(Log10(endp/startp)) > 2) Then
      !use a log scale
      If (npts /= 1) Then
        pressincr = Exp(Log(endp/startp)/(npts-1.0_RDbl))
      Endif
      fuglist(1) = startp
      Do i=1, npts-1
        fuglist(i+1) = startp*(pressincr**i)
      End Do
    Else
      If (npts /= 1) Then
        pressincr = (endp - startp)/(npts - 1.0_RDbl)
      End If
      fuglist(1) = startp
      Do i=1, npts-1
        fuglist(i+1) = startp + pressincr*i
      End Do      
    Endif
    Return
  End Subroutine gcmcpc_genpressures

  !---------------------------------------------------------
  ! Read the fugacities from the file "filename"
  !---------------------------------------------------------
  Subroutine gcmcpc_readpressures(fuglist,pressurefile,sorbtype)
    Real(kind=RDbl), Dimension(:),Intent(inout)  :: fuglist
    Character(*), Intent(in) :: pressurefile
    Integer,Intent(in) :: sorbtype
    Integer     :: unitno, lineno, npts, i, error
    Character(len=strLen)    :: line, sorbname

    !** Open the file if not already open
    unitno = file_open(pressurefile)

    !** Find the sorbate name in the file
    sorbname = molecules_name(sorbtype)
    lineno = filesrchstr(unitno, Trim(sorbname), line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4, 4a)') __FILE__," : ",__LINE__, &
          " Could not find the sorbate ", Trim(sorbname), " in the file ",&
          Trim(pressurefile )
      Stop
    End If

    !** Read the no. of points
    Read(unitno, *) npts
    ! Check that the no. of points is the same as that specfied
    ! in the control file
    If (npts /= postgcmc%nsims) Then
      Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
          " The no. of points in the pressure file ", Trim(pressurefile), &
          " does not match that in the control file"
      Stop
    End If
    Read(unitno, *) (fuglist(i), i=1, npts)
    Close(unit=unitno)

  End Subroutine gcmcpc_readpressures


  
End Module gcmcpc






