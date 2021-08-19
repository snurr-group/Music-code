!---------------------------------------------------------------------------
! This contains a main program for doing hybridgcmc.  
! hybridgcmc needs a library of configurations to be maintained 
! during the simulation.
! This drivers runs two simulatiosn in parallel. One the main simulation 
! and the other for updating the library
! The library update simulation basically has an Ideal gas NVT md simulation 
! of one molecule of each of the sorbates
!---------------------------------------------------------------------------

Program music

  Use atom, Only: atom_display, atom_init, atom_sampleCF
  Use commandline, Only : commandline_init
  Use config, Only: config_writerestartfile, config_init, config_sampleCF, &
      AtMolCoords, config_display, config_dump, config_displaycoords, &
      config_makeconsistent
  Use datafile, Only: CONFILE, datafile_writeconfig, datafile_initout, &
      datafile_setExtraInfo
  Use defaults, Only: strLen, RDbl, d_ctrl_file, d_res_file, pcalls, &
      dashedline, MAX_ROUTINES, ROUT_MAIN, zero, dbgflag, one
  Use file, Only: file_getunit, file_settag, file_open
  Use general, Only: genparams, general_init, general_sampleCF
  Use hybridmc, Only: hybridmc_init,hybridmc_initdisplay,HYBRIDMC_Params, &
      hybridmc_dosim, hybridmc_libUpdtFreq, hybridmc_chknrgs, &
      hybridmc_getnosims,hybridmc_displaystats,hybridmc_displaysimparams, &
      hybridmc_NoOfLibsorbs, hybridmc_beginsim, hybridmc_temperature, &
      hybridmc_pressure
  Use interact, Only: Interaction_Model, interact_init, interact_display, &
      interact_chkrecalc
  Use md, Only: MDInfo, md_initEquil, md_simsummary, md_simdisplay, md_dosim, &
      md_display, md_sampleCF
  Use molecules, Only: molecules_display, molecules_init, molecules_sampleCF
  Use movie, Only: MovieInfo, movie_makemovie, movie_init, movie_display
  Use simcell, Only: SimCell_Params, simcell_display, simcell_init, &
      simcell_sampleCF
  Use storetop, Only: storetop_display
  Use utils, Only: genfilename, int2str, real2str, allocerrdisplay, &
      filesrchstr, cleanstring
  Use vector, Only: VecType
  Use wtime, Only: wtime_init, wtime_display, wtime_starttime, wtime_stoptime
  Use track, Only: Tracking_Info,track_init,track_process,track_display

  Implicit None

  !-----------------------------------------------------------------------------
  ! The variables common to both main simulations and library-update simulations
  !-----------------------------------------------------------------------------
  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** Tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action

  !** Simulation Cell Variable(s)
  Character(len=strLen),Parameter :: simcell_tag="Simulation Cell Information"
  Type(SimCell_Params)            :: scell

  !** Movie Parameters
  Type(MovieInfo) :: MovieParams
  Character(len=strLen), Parameter :: Movie_tag = "Movie Information"

  !** Tracking Parameters
  Logical              :: tracking
  Type(Tracking_Info)  :: trackinfo

  !** set the default display unit
  Integer, Parameter         :: dUnit = 6

  Logical                    :: lib_update_yes, rewindFlag
  Integer                    :: error, nmols, lineno, ctrlunit, nlibsorbs
  Integer                    :: simno, nsims, iter
  Character(len=strLen)      :: string1,string2
  Real(kind=RDbl)            :: simT, simP



  !-----------------------------------------------------------------------
  ! Varibales used for main simulation only
  !-----------------------------------------------------------------------
  !** Configuration Variable(s)
  Type(AtMolCoords), Dimension(:), Pointer :: species
  Character(len=strLen), Parameter         :: config_tag = & 
      "Configuration Initialization"

  !** Interaction/Forcefield Parameters and Statistics
  Type(Interaction_Model)                  :: imodel

  !** Configuration file
  Type(CONFILE) :: configfile
  Character(len=strLen) :: restartfile

  Integer :: lib_update_frequency

  !** hybridgcmc object
  Character(len=strLen),Parameter :: hybridmc_tag=&
      "HYBRID GCMC Information"
  Type(HYBRIDMC_Params) :: hybridmcparams
  Real(kind=RDbl) :: lib_mtime




  !-----------------------------------------------------------------------
  ! Varibales used for library update simulation simulation only
  !-----------------------------------------------------------------------
  Type(AtMolCoords), Dimension(:), Pointer :: lib_Species
  Type(MDInfo) :: lib_mdparams

  !** Interaction/Forcefield Parameters and Statistics
  Type(Interaction_Model)                  :: lib_imodel
  Character(len=strLen) :: lib_update_tag, lib_config_tag, lib_imodel_tag 
  Character(len=strLen) :: lib_moves_tag, lib_datafile_tag, lib_datafile_name
  Character(len=strLen) :: text 
  Integer :: lib_datafile_niters, lib_datafile_iconfig

  !** Configuration file
  Type(CONFILE) :: lib_configfile

  !-------------------------------------------------------------
  ! Do basic initialization, based on command line input
  !-------------------------------------------------------------
  Call commandline_init(ctrl_filename,action)

  If (action=="WriteSampleCtrlfile") Then

!!$    Write(0,'(a)') dashedline
!!$    Call general_sampleCF(0)
!!$    Write(0,'(a)') dashedline
!!$    Call atom_sampleCF(0)
!!$    Write(0,'(a)') dashedline
!!$    Call molecules_sampleCF(0)
!!$    Write(0,'(a)') dashedline
!!$    Call simcell_sampleCF(0)
!!$    Write(0,'(a)') dashedline
!!$    Call config_sampleCF(0)
!!$    Write(0,'(a)') dashedline
!!$    Call hybridmc_sampleCF(0)
!!$    Write(0,'(a)') dashedline
!!$    Stop
  Else If (action /= "DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  End If

  Call file_settag(ctrl_filename, d_ctrl_file)

  !-----------------------------------------------------------------
  ! Welcome the user and then initialize the main simulation details
  !-----------------------------------------------------------------
  Write(dUnit,'(a)') "Welcome to Music"
  Write(dUnit,'(2a)') "Reading from control file ",Trim(ctrl_filename)
  Call general_init(ctrl_filename)

  !----------------------------------------------------------------------------
  ! Initialize atoms, molecules, simulation cell, configuration, forcefield
  !----------------------------------------------------------------------------
  Call atom_init()
  Call atom_display(dUnit)

  nmols = molecules_init()
  Call molecules_display(dUnit)

  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,2,dUnit)

  !** Allocate the species structure
  Allocate(species(nmols), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'species structure')    

  Call config_init(species, scell, ctrl_filename, config_tag)
  Call config_display(species,dUnit)
  !  Call config_displaycoords(species,1,dUnit)

  Call interact_init(imodel,ctrl_filename,scell,species,'HMC')
!  Call interact_init(imodel,ctrl_filename,scell,species,'MC')  !** HACK for test
  Call interact_display(imodel,2,dUnit)

  !** If DOF is changed during interact, and molecule becomes internally_flexible
  !** Then we need to adjust the gccords too.  This is important for HGCMC
  Call config_makeConsistent(species)

  !** Initialize the tracking
  tracking = track_init(trackinfo,imodel,ctrl_filename)
  If (tracking) Call track_display(trackinfo,2,dunit)

  !** Initialize the movie
  Call movie_init(MovieParams,ctrl_filename,Movie_tag)
  Call movie_display(MovieParams,dunit,6)

  !---------------------------------------------------------------------
  ! Initialize the move type information. Also get back the section tag
  ! for library-update. 
  !---------------------------------------------------------------------
  Call hybridmc_init(hybridmcparams, imodel, species, scell, ctrl_filename, &
      hybridmc_tag, lib_update_tag)
  Call hybridmc_initdisplay(hybridmcparams,2,dUnit)
  Write(dUnit,'(a)') " Finished initializing main simulation details"

  !---------------------------------------------------------------------
  ! Initialize Library Updates if desired
  !---------------------------------------------------------------------
  !** Find out what sort of Library Update, if any, will be used
  If (cleanstring(lib_update_tag) == "NO_UPDATE" ) Then

    !** O.K then this is deliberately left out
    Write(dUnit,'(a)') " ************************ "
    Write(dUnit,*) "-- NO_UPDATE -- specified as the library updatetag "
    Write(dUnit,*) "The libraries used for hybridmcmc will not updated ", &
        "during the simulation"
    Write(dUnit,*) "Please make sure that the libraries are large enough.", &
        "  This will work only"
    Write(dUnit,*) "when the sorbate-molecules are small and ",&
        "do not have many degrees of freedom" 
    lib_update_yes = .False.

  Else
    Write(dUnit,'(2a)') "Looking for tag -"// Trim(lib_update_tag)//&
        " for initializing library files update" 

    !** Open the control file and search for tag
    ctrlunit = file_open(ctrl_filename,110)
    rewindflag = .True.
    lineno = filesrchstr(ctrlunit, lib_update_tag,text,rewindflag)

    !** Read the Library Update tag
    If (lineno/=0) Then
      lib_update_yes=.True.
    Else
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) " Somthing wrong with the given tag for libray update : ", &
          Trim(lib_update_tag)
      !** This was not intended, stop the simulation
      Write(0,*) " Stopping the simulation. You can modify the lib_update_tag &
          & in hybridmcmc section to NO_UPDATE and then re-run this"
      Stop
    End If
  End If

  !** Initialize the library-update simulations if desired
  If (lib_update_yes) Then
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Read(ctrlunit,'(a)') lib_config_tag
    Read(ctrlunit,'(a)') lib_imodel_tag
    Read(ctrlunit,'(a)') lib_moves_tag
    Read(ctrlunit,'(a)') lib_datafile_tag
    Read(ctrlunit,'(a)') lib_datafile_name
    Read(ctrlunit,*) lib_datafile_niters, lib_datafile_iconfig
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__

    !** intialize config  for this sub-simulation
    nlibsorbs=hybridmc_NoOfLibSorbs(hybridmcparams)
    Allocate(lib_species(nlibsorbs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'species structure')

    lib_config_tag=cleanstring(lib_config_tag)
    Call config_init(lib_species, scell, ctrl_filename, lib_config_tag)
    Call config_display(lib_species,dUnit)
!   Call config_displaycoords(species,1,dUnit)

    !** initialze interactions for sub-simulation
    Call interact_init(lib_imodel,ctrl_filename,scell,lib_species,&
        'MD', lib_imodel_tag)
    Call interact_display(lib_imodel,2,dUnit)

    !** Initialize the datafile used for sub-simulation
    lib_configfile%name=cleanstring(lib_datafile_name)
    Call datafile_initout(lib_configfile,ctrl_filename,lib_species,0, &
        lib_datafile_niters, lib_datafile_iconfig, "MD", lib_datafile_tag)

    !** Initialize the moves section (md) for sub-simulation 
    Call md_initEquil(lib_mdParams, lib_imodel, scell,lib_species, &
        ctrl_filename, lib_moves_tag, lib_configfile)
  End If

  nsims = hybridmc_getnosims(hybridmcparams)

  !----------------------------------------------
  ! Initialize the timer
  !----------------------------------------------
  Call wtime_init()
  Call wtime_starttime(ROUT_MAIN)

  !---------------------------------------------------------------
  ! Conduct the simulations, loop over simulations and iterations
  !---------------------------------------------------------------
  Do simno = 1,nsims
    !** Open the configuration file for the Main simulation and close previous
    configfile%name=&
        genfilename(genparams%configfile,simno+genparams%simstart-1)

    simT=hybridmc_temperature(hybridmcparams, simno)
    ! ** If multi-component , pressure of first compound
    simP=hybridmc_Pressure(hybridmcparams, simno) 
    Call datafile_setExtraInfo(configfile,"temperature",simT)
    Call datafile_setExtraInfo(configfile,"pressure",simP)
    Write(*,*) simT, simP
    Call datafile_initout(configfile,ctrl_filename,species, simno, &
        genparams%niterations, genparams%iconfig,"GCMC")

    !** Establish name of restartfile
    restartfile = genfilename(genparams%restartfile,simno+genparams%simstart-1)
    Call file_settag(restartfile, d_res_file)

    !** reset stat counters for this sim
    Call hybridmc_beginsim(hybridmcparams, simno, species, scell, 2, dUnit)

    !** Display the parameters for this simulation
    Call hybridmc_displaysimparams(hybridmcparams,simno,2,dUnit)

    !** Find out after how many iterations library files should be updated
    If (lib_update_yes) Then
      lib_update_frequency=hybridmc_libUpdtFreq(hybridmcparams)
    Else
      !** assign a large value so that it is never called, !** HACK
      lib_update_frequency=genparams%niterations+1
    End If

!!$    If (Associated(hybridmcparams%gcmcsorbs)) Then
!!$      lib_update_frequency=hybridmc_libUpdtFreq(hybridmcparams)

    Do iter = 1, genparams%niterations

      !-------------------------------------------------------
      ! Main Simulation Part
      !-------------------------------------------------------
      genparams%currentiteration = iter

      !** Do the NVTMC and GCMC simulations
      Call hybridmc_dosim(hybridmcparams, species, scell, configfile, simno)

      !** Dump to the standard output
      If (Mod(iter, genparams%iprint) == 0) Then
        string1 = int2str(iter)
        string2 = int2str(simno)
        Write(dUnit,'(1x,4a)') 'Iteration number ',Trim(string1), &
            ' of simulation ',Trim(string2)
        Call hybridmc_displaystats(hybridmcparams, species, imodel, &
            scell, simno, 1,6)
        Write(*,*) ! Put a space between outputs
      End If

      !** Handle the tracking if it's on
      If (tracking) Call track_process(trackinfo,species,scell,imodel, &
          iter,2,dunit)

      !** Dump to the config file
      If (Mod(iter, genparams%iconfig) == 0) Then
        Call datafile_writeconfig(species, imodel%spcstats, configfile)
      End If

      !** Dump to the crash file
      If (Mod(iter, genparams%icrash) == 0) Then
        Call config_writerestartfile(species, restartfile)
      End If

      !** Dump to movie file if need be
      Call movie_makemovie(movieParams,species,scell,iter,iter*one)

#ifdef DEBUG
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Call interact_chkrecalc(imodel,species,scell,1.0e-6_RDbl,2,dUnit)
!      Call hybridmc_chknrgs(hybridmcparams,imodel,species,scell,2,dUnit)
#endif

      !-------------------------------------------------------
      ! Library Update Simulation Part
      !-------------------------------------------------------
      If (lib_update_yes) Then
        If (Mod(iter, lib_update_frequency) == 0) Then
          lib_mtime=zero

          !** Note that this dosim is called only once. The mdmove should
          !** be initialized in such a way that it has enough number of
          !** iterations, and writes into a datafile enough number of times
          Call md_dosim(lib_mdParams, 1, scell, lib_species, &
              lib_configfile, lib_mtime)

          !** Display the final information FOR-DEBUG
          Call md_simsummary(lib_imodel,lib_mdParams,lib_species,iter,&
              lib_mtime,dUnit)

          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
!!$          Call hybridmc_libupdate()
        End If
      End If

    End Do  !** move iteration
  End Do  !** simulation number

  !---------------------------------------
  ! Get the elapsed time
  !---------------------------------------
  Call wtime_stoptime(ROUT_MAIN)
  Call wtime_display("Main Program", ROUT_MAIN)

End Program music

