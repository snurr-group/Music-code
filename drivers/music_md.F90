!-------------------------------------------------------------------
! This contains the main program for Molecular Dynamics simulations
!-------------------------------------------------------------------

Program music

  Use defaults, Only: strLen, RDbl, d_ctrl_file, d_res_file, pcalls, &
      dashedline, MAX_ROUTINES, one
  Use utils, Only: genfilename, int2str, real2str, allocerrdisplay
  Use file, Only: file_getunit, file_settag
  Use vector, Only: VecType
  Use commandline, Only : commandline_init
  Use general, Only: genparams, general_init, general_sampleCF
  Use wtime, Only: wtime_init, wtime_display, wtime_starttime, wtime_stoptime
  Use atom, Only: atom_display, atom_init, atom_sampleCF
  Use molecules, Only: molecules_display, molecules_init, molecules_sampleCF
  Use datafile, Only: CONFILE, datafile_writeconfig, datafile_initout
  Use simcell, Only: SimCell_Params, simcell_display, simcell_init, &
      simcell_sampleCF
  Use config, Only: config_writerestartfile, config_init, config_sampleCF, &
      AtMolCoords, config_display, config_dump, config_displaycoords, &
      config_changerp
  Use interact, Only: Interaction_Model, interact_init, interact_display
  Use md, Only: MDInfo, md_initEquil, md_simsummary, md_simdisplay, md_dosim, &
      md_display, md_sampleCF
  Use movie, Only: MovieInfo, movie_init, movie_display, movie_makemovie
  Use stopcriteria, Only: STOP_CHECK, stopcriteria_init, stopcriteria_check
  Use track, Only: Tracking_Info,track_init,track_process,track_display

  Implicit None

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** Tells what to do. Ex:  run simulation?, help?, write sample file?
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

  !** MD Simulation Parameters
  Type(MDInfo) :: mdParams
  Character(len=strLen), Parameter :: md_tag = "MD Information"

  !** Movie Parameters
  Type(MovieInfo)                  :: mainMovieParams
  Character(len=strLen), Parameter :: mainMovie_tag = "Movie Information"

  !** Configuration file
  Type(CONFILE) :: configfile

  !** MD Simulation time
  Real(Kind=RDbl)            :: mtime

  !** Tracking Parameters
  Logical              :: tracking
  Type(Tracking_Info)  :: trackinfo

  !** set the default display unit
  Integer, Parameter         :: dUnit = 6

  Integer                    :: nsims, simno, iter, error, nmols
  Integer                    :: configunitno, restartunitno   
  Character(len=strLen)      :: comment, filename, restartfile

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
    Call md_sampleCF(0)
    Write(0,'(a)') dashedline
    Stop
  Else If (action /= "DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  End If

  Call file_settag(ctrl_filename, d_ctrl_file)

  !----------------------------------------------
  ! Welcome the user
  !----------------------------------------------
  Write(dUnit,'(a)') "Welcome to Music"
  Write(dUnit,'(2a)') "Reading from control file ",Trim(ctrl_filename)

  !----------------------------------------------
  ! Initialize the general parameters
  !----------------------------------------------
  Call general_init(ctrl_filename)

  !----------------------------------------------------------------------------
  ! Initialize atoms, molecules, simulation cell, configuration, forcefield
  !----------------------------------------------------------------------------
  Call atom_init()
  Call atom_display(dUnit)

  nmols = molecules_init()

  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,2,dUnit)
  Call molecules_display(dUnit)

  !** Allocate the species structure
  Allocate(species(nmols), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'species structure')    

  Call config_init(species, scell, ctrl_filename, config_tag)
  Call config_display(species,dUnit)
  Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
  !  Call config_displaycoords(species,1,dUnit)
  !  Call config_displaycoords(species,2,dUnit)
  !  stop

  Call interact_init(imodel,ctrl_filename,scell,species,'MD')
  Call interact_display(imodel,2,dUnit)

  !----------------------------------------------------------
  ! Initialize the move type information
  !----------------------------------------------------------
  Call wtime_init()
  Call wtime_starttime(MAX_ROUTINES)
  Call wtime_init(1)
  Call wtime_init(2)
  Call wtime_init(3)

  Call stopcriteria_init(ctrl_filename)

  simno = 1

  !** Initialize the file tags and names and 
  !** Open the configuration file
  configfile%name =  genfilename(genparams%configfile, &
      simno+genparams%simstart-1)
!!$
!!$  !** need to write a routine for extracting md-temeperature
!!$  Call datafile_setExtraInfo(configfile,mdParams%simT,"temperature")
  Call datafile_initout(configfile,ctrl_filename,species, simno, &
      genparams%niterations, genparams%iconfig,"MD")

  !** Initialize the MD stuff
  Call md_initEquil(mdParams, imodel, scell,species,ctrl_filename, &
      md_tag,configfile)

  !** Write to the restart file
  restartfile = genfilename(genparams%restartfile,simno+genparams%simstart-1)
  Call file_settag(restartfile, d_res_file)
  restartunitno  = file_getunit(restartfile)
  Call config_writerestartfile(species,restartfile)

  !** Display the MD information
  Call md_display(mdParams,dUnit)

  !** Initialize the tracking
  tracking = track_init(trackinfo,imodel,ctrl_filename)
  If (tracking) Call track_display(trackinfo,2,dunit)

  !** Initialize the movie
  Call movie_init(mainMovieParams,ctrl_filename,mainMovie_tag)
  Call movie_display(mainMovieParams,dUnit,6)

  !MDEBUG
  Write(dUnit,'(1x,a)') "Begining Simulation"
  Write(dUnit,'(1x,a,i10)') "Number of Iterations to Perform :", &
      genparams%niterations

  nsims = 1

  Do simno=1, nsims
    Write(*,'(1x,2a,i4,a,i4)') __FILE__," : ",__LINE__," Start of sim :", simno
    mtime = 0.0_RDbl

    Do iter = 1, genparams%niterations

      Call stopcriteria_check(species,imodel%spcstats,restartfile)

      Call wtime_starttime(1)
      Call md_dosim(mdParams,1,scell,species,configfile,mtime)
      Call wtime_stoptime(1)

      !** Dump to std io
      If (Mod(iter, genparams%iprint) == 0) Then
        Call md_simdisplay(imodel,mdParams,species,scell,iter,mtime,dUnit)
        Write(dunit,*)
      End If

      !** Handle the tracking if it's on
      If (tracking) Call track_process(trackinfo,species,scell,imodel, &
          iter,2,dunit)

      !** Dump to the config file
      If (Mod(iter, genparams%iconfig) == 0) Then
        Call datafile_writeconfig(species,imodel%spcstats,configfile,mtime)
      End If

      !** Dump to the crash file
      If (Mod(iter+1, genparams%icrash) == 0) Then
        Call config_writerestartfile(species,restartfile)
      End If

#ifdef IDEAL_GAS
      !If you are doing IG simulations molecules tend to fly away from origin
      ! to fix this periodically you can shift them close to origin
      ! this improves intra nrg calculations
      !** Dump to the crash file
      If (Mod(iter, 1000) == 0) Then
        Call config_changerp(species,scell)
      End If
#endif

      !** Dump to movie file if need be
      Call movie_makemovie(mainmovieParams,species,scell,iter,iter*one)

    End Do
  End Do

  !** Display the final information
  Call md_simsummary(imodel,mdParams,species,iter,mtime,dUnit)
  Write(*,*) "pcalls : ",pcalls

  !---------------------------------------
  ! Get the elapsed time
  !---------------------------------------
  Call wtime_stoptime(MAX_ROUTINES)
  Call wtime_display("Main Program", MAX_ROUTINES)

!!$  Call wtime_display('md',1)
!!$  Call wtime_display('integrate',2)
!!$  Call wtime_display('integrator',3)
!!$  Call wtime_display('gear6',4)
  !  Call wtime_display('pairwise',2,6)
End Program music
