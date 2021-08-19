!--------------------------------------------------------
! This contains the main program for GCMC simulation
!--------------------------------------------------------

Program music
  Use datafile, Only: CONFILE,datafile_initout,datafile_writeconfig, &
      datafile_close, datafile_setExtraInfo
  Use defaults, Only: strlen, RDbl, dashedline, d_res_file, d_ctrl_file, &
      ROUT_MAIN, one, dbgflag
  Use utils, Only: genfilename,int2str,allocerrdisplay
  Use file, Only: file_settag
  Use commandline, Only : commandline_init
  Use atom, Only: atom_init, atom_display, atom_sampleCF
  Use interact, Only: Interaction_Model, interact_init, interact_display
  Use config, Only: AtMolCoords, config_init, config_writerestartfile, &
      config_display, config_config2xyz, config_sampleCF, config_getnmoles
  Use gcmc, Only: GCMC_Params, gcmc_init, gcmc_getnosims, &
      gcmc_displaystats, gcmc_dosim, gcmc_displaysimparams, &
      gcmc_initdisplay, gcmc_beginsim, gcmc_endsim, gcmc_temperature, &
      gcmc_pressure, gcmc_sampleCF
  Use general, Only: genparams, general_getContentTag, general_init, &
      general_sampleCF,general_setCurrentSimno
  Use molecules, Only: molecules_init, molecules_sampleCF, &
      molecules_display
  Use storestats, Only: Species_Stats, storestats_displaynrg
  Use simcell, Only: SimCell_Params, simcell_init, simcell_display, &
      simcell_sampleCF
  Use storetop, Only: storetop_display
  Use vector, Only: VecType
  Use wtime, Only: wtime_display, wtime_stoptime, wtime_starttime, wtime_init
  Use movie, Only: MovieInfo, movie_makemovie, movie_init, movie_display
  Use track, Only: Tracking_Info, Tracking_Item, track_init, track_clean, &
      track_process, track_display
  Use stopcriteria, Only: STOP_CHECK, stopcriteria_init, stopcriteria_check

  Implicit None
  !** The  config files parameters
  Type (CONFILE) :: conf_file

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action

  !** Simulation Cell Variable(s)
  Character(len=strLen),Parameter :: simcell_tag="Simulation Cell Information"
  Type(SimCell_Params)            :: scell

  !** Configuration Variable(s)
  Type(AtMolCoords), Dimension(:), Pointer :: species
  Character(len=strLen), Parameter         :: config_tag = & 
      "Configuration Initialization"

  !** Interaction/Forcefield Parameters and Statistics
  Type(Interaction_Model)                  :: imodel

  !** GCMC Simulation Parameters
  Type(GCMC_Params) :: gcmcparams
  Character(len=strLen), Parameter         :: gcmc_tag = &
      "GCMC Information"

  !** Movie Parameters
  Type(MovieInfo) :: MovieParams
  Character(len=strLen), Parameter :: Movie_tag = "Movie Information"

  !** Tracking Parameters
  Logical                          :: tracking
  Type(Tracking_Info)              :: trackinfo

  Integer                :: nsims,simno,iter,content_tag,nmols
  Integer                :: error,dunit
  Character(len=strLen)  :: string1,string2,string3
  Character(len=strlen)  :: restartfile,datfilename
  Real(kind=RDbl)        :: simT, simP
  !-------------------------------------------------------------
  ! Do basic initialization, based on command line input
  !-------------------------------------------------------------
  Call commandline_init(ctrl_filename,action)

  If (action=="WriteSampleCtrlfile") Then
    Write(0,'(a)') dashedline
    Call general_sampleCF(0)
    Write(0,'(a)') dashedline
    Call atom_sampleCF(0)
    Write(0,'(a)') dashedline
    Call molecules_sampleCF(0)
    Write(0,'(a)') dashedline
    Call simcell_sampleCF(0)
    Write(0,'(a)') dashedline
    Call config_sampleCF(0)
    Write(0,'(a)') dashedline
    Call gcmc_sampleCF(0)
    Write(0,'(a)') dashedline
    Stop
  Else If (action /= "DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  End If
      
  Call file_settag(ctrl_filename, d_ctrl_file)

  dunit = 6

  !----------------------------------------------
  ! Welcome the user
  !----------------------------------------------
  Write(dUnit,'(a)') "Welcome to MuSiC (GCMC implementation)"
  Write(dUnit,'(2a)') "Reading from control file ",Trim(ctrl_filename)

  !----------------------------------------------
  ! Initialize the general parameters
  !----------------------------------------------
  Call general_init(ctrl_filename)

  !----------------------------------------------------------------------------
  ! Initialize atoms, molecules, simulation cell, configuration, forcefield
  !----------------------------------------------------------------------------
  Call atom_init()
  Call atom_display(dunit)

  nmols = molecules_init()
  Call molecules_display(dunit)

  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,6,dunit)

  !** Allocate the sorbates structure
  Allocate(species(nmols), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'species structure')    

  Call config_init(species, scell, ctrl_filename, config_tag)
  Call config_display(species,dUnit)

  Call config_config2xyz(species,'initconfig.xyz',.True.,'Initial configuration')

  Call interact_init(imodel,ctrl_filename,scell,species,'MC')
  Call interact_display(imodel,2,dUnit)

  tracking = track_init(trackinfo,imodel,ctrl_filename)
  If (tracking) Call track_display(trackinfo,2,dUnit)

  !----------------------------------------------------------
  ! Initialize the move type information and auxiliaries
  !----------------------------------------------------------

  Call gcmc_init(gcmcparams, imodel, species, scell, ctrl_filename, gcmc_tag)
  Call gcmc_initdisplay(gcmcparams, 0)

  !** Initialize the movie
  Write(dunit,*)
  Call movie_init(MovieParams,ctrl_filename,Movie_tag)
  Call movie_display(MovieParams,dunit,2)
  
  !** Initialize the timer
  Write(dunit,*)
  Call wtime_init()
  Call wtime_starttime(ROUT_MAIN)

  Call stopcriteria_init(ctrl_filename)

  !----------------------------------------------------------
  ! Run the simulation(s)
  !----------------------------------------------------------

  nsims = gcmc_getnosims(gcmcparams)
 
  Do simno = 1,nsims
    !** Open the configuration file for the new simulation 
    conf_file%name = genfilename(genparams%configfile,&
        simno+genparams%simstart-1)
    content_tag = general_getContenttag()
    Call general_setCurrentSimno(simno)


    simT=gcmc_temperature(gcmcparams, simno)
    ! ** If multi-component , pressure of first compound
    simP=gcmc_Pressure(gcmcparams, simno) 
    Call datafile_setExtraInfo(conf_file,"temperature",simT)
    Call datafile_setExtraInfo(conf_file,"pressure",simP)
    Call datafile_initout(conf_file,ctrl_filename,species,simno,&
        genparams%niterations, genparams%iconfig,"GCMC")

    restartfile = genfilename(genparams%restartfile,simno+genparams%simstart-1)
    Call file_settag(restartfile, d_res_file)

    !** Begin the simulation
    Call gcmc_beginsim(gcmcparams,simno,species,scell,0,dunit)

    !** Display the parameters for this simulation
    Write(dunit,'(a)') dashedline
    string1 = int2str(simno)
    string2 = int2str(gcmc_getnosims(gcmcparams))
    string3 = int2str(genparams%niterations)
    Write(dunit,'(4a,3x,3a)') 'Beginning GCMC simulation number ', &
        Trim(string1),' of ',Trim(string2),'(',Trim(string3),' steps)'
    Write(dunit,'(a)') dashedline
    Call gcmc_displaysimparams(gcmcparams,simno,0,dunit)
    Write(*,'(/)')

    !** run the simulation
    Do iter = 1,genparams%niterations
      Call stopcriteria_check(species,imodel%spcstats,restartfile)
      genparams%currentiteration = iter

      !** Do the GCMC simulation
      Call gcmc_dosim(gcmcparams,species,scell,simno)
      
      !** Dump to the standard output
      If (Mod(iter, genparams%iprint) == 0) Then
        string1 = int2str(iter)
        string2 = int2str(iter*gcmcparams%niterations)
        string3 = int2str(simno)
        Write(*,'(2x,7a)') 'Macro-iteration number ',Trim(string1), &
        ' (',Trim(string2),' micro)',' of simulation ',Trim(string3)
        Call gcmc_displaystats(gcmcparams,imodel,species,scell,simno,2,dunit)
        Write(*,*) 
      End If

      !** Handle tracking
      If (tracking) Call track_process(trackinfo,species,scell, &
          imodel,iter,2,dunit)

      !** Dump to the config file
      If (Mod(iter, genparams%iconfig) == 0) Then
        Call datafile_writeconfig(species,imodel%spcstats, conf_file)
      End If
      
      !** Dump to the crash file
      If (Mod(iter, genparams%icrash) == 0) Then
        Call config_writerestartfile(species, restartfile)
      End If

      !** Dump to movie file if need be
      Call movie_makemovie(movieParams,species,scell,iter,iter*one)
    End Do

    !** End the simulation
    Call gcmc_endsim(gcmcparams,simno,imodel,species,scell,'',0,dunit)

    !** Close the configuration file
    Call datafile_close(conf_file)

  End Do

  Call config_config2xyz(species,'finalconfig.xyz',.True.,'Final configuration')

!!  Call storetop_display(imodel%ff%results,.True.,2,6)
  
  !---------------------------------------
  ! Get the elapsed time
  !---------------------------------------
  Call wtime_stoptime(ROUT_MAIN)
  Call wtime_display("Main Program", ROUT_MAIN)

End Program music
