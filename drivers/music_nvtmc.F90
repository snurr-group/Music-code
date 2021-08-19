!--------------------------------------------------------
! This contains the main program for NVTMC simulation
!--------------------------------------------------------

Program music
  Use file, Only: file_settag
  Use defaults, Only: strlen, RDbl, dashedline, d_res_file, d_ctrl_file, &
      ROUT_MAIN,one
  Use utils, Only: genfilename,int2str,allocerrdisplay
  Use commandline, Only: commandline_init
  Use wtime, Only: wtime_display, wtime_stoptime, wtime_starttime, wtime_init
  Use general, Only: genparams, general_getContentTag, general_init,  &
      general_samplecf
  Use atom, Only: atom_init, atom_display, atom_samplecf
  Use molecules, Only: molecules_init, molecules_display, molecules_samplecf
  Use simcell, Only: SimCell_Params, simcell_init, simcell_display, &
      simcell_samplecf
  Use config, Only: AtMolCoords, config_init, config_writerestartfile, &
      config_display, config_samplecf, config_config2xyz
  Use interact, Only: Interaction_Model, interact_init, interact_display
  Use nvtmc, Only: NVTMC_Params, nvtmc_init, nvtmc_getnosims, &
      nvtmc_displaystats, nvtmc_dosim, nvtmc_displaysimparams, &
      nvtmc_initdisplay, nvtmc_beginsim, nvtmc_endsim, nvtmc_sampleCF
  Use datafile, Only: CONFILE,datafile_initout,datafile_writeconfig, &
      datafile_close
  Use vector, Only: VecType
  Use movie, Only: MovieInfo, movie_init, movie_display, movie_makemovie
  Use storestats, Only: storestats_getnoncoul,storestats_getcoul
  Use track, Only: Tracking_Info,track_init,track_process,track_display

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

  !** NVTMC Simulation Parameters
  Type(NVTMC_Params) :: nvtmcparams
  Character(len=strLen), Parameter         :: nvtmc_tag = &
      "NVTMC Information"

  !** Movie Parameters
  Type(MovieInfo) :: MovieParams
  Character(len=strLen), Parameter :: Movie_tag = "Movie Information"

  !** Tracking Parameters
  Logical              :: tracking
  Type(Tracking_Info)  :: trackinfo

  Integer                :: nsims,simno,iter,content_tag,dunit,nspc,error
!  Real(kind=RDbl)        :: noncoul,coul,sum
  Real(kind=RDbl)        :: accep_ratio
  Character(len=strLen)  :: restartfile
  Character(len=strlen)  :: datfilename
  Character(len=strlen)  :: string1,string2,string3

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
    Call nvtmc_sampleCF(0)
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
  Write(dUnit,'(a)') "Welcome to MuSiC (NVTMC implementation)"
  Write(dUnit,'(2a)') "Reading from control file ",Trim(ctrl_filename)

  !----------------------------------------------
  ! Initialize the general parameters
  !----------------------------------------------

  Call general_init(ctrl_filename)

  !----------------------------------------------------------------------------
  ! Initialize atoms, molecules, simulation cell, configuration, forcefield
  !----------------------------------------------------------------------------

  Call atom_init()
  Call atom_display(6)

  nspc = molecules_init()
  Call molecules_display(6)

  !** Allocate the species structure
  Allocate(species(nspc), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'species structure')    

  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,2,6)

  Call config_init(species, scell, ctrl_filename, config_tag)
  Call config_display(species)

  Call config_config2xyz(species,'initconfig.xyz',.True.,'Initial configuration')

  Call interact_init(imodel,ctrl_filename,scell,species,'MC')
  Call interact_display(imodel,2,dUnit)

  !----------------------------------------------------------
  ! Initialize the move type information and auxiliaries
  !----------------------------------------------------------

  Call nvtmc_init(nvtmcparams, imodel, species, scell, ctrl_filename, nvtmc_tag)
  Call nvtmc_initdisplay(nvtmcparams,0)

  !** Initialize the tracking 
  tracking = track_init(trackinfo,imodel,ctrl_filename)
  If (tracking) Call track_display(trackinfo,2,dunit)

  !** Initialize the movie
  Call movie_init(MovieParams,ctrl_filename,Movie_tag)
  Call movie_display(MovieParams,dunit,2)

  !** Initialize the timer
  Call wtime_init()
  Call wtime_starttime(ROUT_MAIN)

  !----------------------------------------------------------
  ! Run the simulation(s)
  !----------------------------------------------------------

  nsims = nvtmc_getnosims(nvtmcparams)
 
  Do simno = 1,nsims
    !** Open the configuration file for the new simulation 
    conf_file%name = genfilename(genparams%configfile,simno+genparams%simstart-1)
    content_tag = general_getContenttag()
!    Call datafile_setExtraInfo(conf_file,"temperature",simT)
    Call datafile_initout(conf_file,ctrl_filename,species,simno,&
        genparams%niterations, genparams%iconfig,"NVTMC")

    restartfile = genfilename(genparams%restartfile,simno+genparams%simstart-1)
    Call file_settag(restartfile, d_res_file)

    !** Begin the simulation
    Call nvtmc_beginsim(nvtmcparams,simno,species,scell,0,dunit)

    !** Display the parameters for this simulation
    Write(dunit,'(a)') dashedline
    string1 = int2str(simno)
    string2 = int2str(nvtmc_getnosims(nvtmcparams))
    string3 = int2str(genparams%niterations)
    Write(dunit,'(4a,3x,3a)') 'Beginning NVTMC simulation number ', &
        Trim(string1),' of ',Trim(string2),'(',Trim(string3),' steps)'
    Write(dunit,'(a)') dashedline
    Call nvtmc_displaysimparams(nvtmcparams,simno,0,dunit)
    Write(*,'(/)')

    !** Dump initial configuration to movie file if need be
    Call movie_makemovie(movieParams,species,scell,0,0.0_RDbl) 

    !** run the simulation
    Do iter = 1,genparams%niterations

      genparams%currentiteration = iter

      !** Do the NVTMC simulation
      accep_ratio = nvtmc_dosim(nvtmcparams,species,scell,simno)

      !** Dump to movie file if need be
      Call movie_makemovie(movieParams,species,scell,iter,iter*one) 

      !** Handle the tracking if it's on
      If (tracking) Call track_process(trackinfo,species,scell,imodel, &
          iter,2,dunit)

      !** Dump to the standard output
      If (Mod(iter, genparams%iprint) == 0) Then
        string1 = int2str(iter)
        string2 = int2str(iter*nvtmcparams%niterations)
        string3 = int2str(simno)
        Write(*,'(2x,7a)') 'Macro-iteration number ',Trim(string1), &
        ' (',Trim(string2),' micro)',' of simulation ',Trim(string3)
        Call nvtmc_displaystats(nvtmcparams,imodel,species,scell,simno,2,dunit)

#ifdef DEBUG
        noncoul = storestats_getnoncoul(imodel%spcstats,1,2,'inst')
        coul = storestats_getcoul(imodel%spcstats,1,2,'INST')        
        Write(*,'(a,i8,a,3f10.4)') 'Iter: ',iter, &
            ' Species 1,2 INST: noncoul,coul,sum: ',noncoul,coul,(noncoul+coul)

        noncoul = storestats_getnoncoul(imodel%spcstats,1,2,'CAVG')
        coul = storestats_getcoul(imodel%spcstats,1,2,'CAVG')        
        Write(*,'(a,i8,a,3f10.4)') 'Iter: ',iter, &
            ' Species 1,2 CAVG: noncoul,coul,sum: ',noncoul,coul,(noncoul+coul)
#endif        

        Write(*,*) 
      End If

      !** Dump to the config file
      If (Mod(iter, genparams%iconfig) == 0) Then
        Call datafile_writeconfig(species,imodel%spcstats,conf_file)
      End If
      
      !** Dump to the crash file
      If (Mod(iter, genparams%icrash) == 0) Then
        Call config_writerestartfile(species, restartfile)
      End If
    End Do

    !** End the simulation
    Call nvtmc_endsim(nvtmcparams,imodel,simno,species,scell,'',0,dunit)

    !** Close the configuration file
    Call datafile_close(conf_file)

  End Do
  
  !---------------------------------------
  ! Get the elapsed time
  !---------------------------------------
  Call wtime_stoptime(ROUT_MAIN)
  Call wtime_display("Main Program", ROUT_MAIN)

End Program music
