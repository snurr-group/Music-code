!-----------------------------------------------
! This contains the main program.  
!-----------------------------------------------

Program music
  Use file, Only: file_settag
  Use defaults, Only: strLen, RDbl, d_ctrl_file, dashedline, d_res_file, &
      ROUT_MAIN
  Use utils, Only: genfilename
  Use wtime, Only: wtime_init, wtime_starttime, wtime_display, wtime_stoptime
  Use general, Only: general_init, general_getfiledesc, genparams
  Use atom, Only: atom_init, atom_display
  Use molecules, Only: molecules_init, molecules_display
  Use simcell, Only: Simcell_Params, simcell_init, simcell_display
  Use forcefield, Only: forcefield_init, forcefield_display
  Use config, Only: config_init, config_display
  Use visxyz, Only:
  Use datafile, Only:
  Use mepbaker, Only: mepbaker_dosim, MEPBAKER_Params, mepbaker_init, &
      mepbaker_display

  Implicit None

  !** Program call variables
  Integer               :: nargs, iargc
  Character(len=strLen) :: progname, ctrl_filename

  !** Simulation Cell Variable(s)
  Character(len=strLen),Parameter :: simcell_tag="Simulation Cell Information"
  Type(SimCell_Params)            :: scell

  !** Configuration Variable(s)
  Type(AtMolCoords), Dimension(:), Pointer :: sorbates
  Character(len=strLen), Parameter         :: config_tag = & 
      "Configuration Initialization"

  !** Minimum Energy Path Parameters
  Type(MEPBAKER_Params) :: mepparams
  Character(len=strLen), Parameter         :: mepbaker_tag = &
      "Baker Minimum Energy Path"

  Integer        :: nplanes, planeno, iter, restartunitno
  Character(len=strLen)  :: restartfile

  !----------------------------------------------------------------------------
  ! Get the name of the control file from the command line
  !----------------------------------------------------------------------------
  nargs = iargc()    
                     
  If (nargs < 1) Then
    Call getarg(0, progname)
    Write(0,*) 'Usage: ',Trim(progname),' run_filename'
    Stop
  End If 
  Call getarg(1, ctrl_filename)
  Call file_settag(ctrl_filename,d_ctrl_file)

  !----------------------------------------------
  ! Initialize the timer
  !----------------------------------------------
  Call wtime_init()
  Call wtime_starttime(ROUT_MAIN)

  !----------------------------------------------
  ! Initialize the general parameters
  !----------------------------------------------
  Call general_init(ctrl_filename)

  !----------------------------------------------------------------------------
  ! Initialize atoms, molecules, zeolite
  !----------------------------------------------------------------------------
  Call atom_init()
  Call atom_display(6)

  Call molecules_init()
  Call molecules_display(6)

  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,6)

  !----------------------------------------------------------
  ! Initialize the configuration
  !----------------------------------------------------------
  Call config_init(sorbates, scell, ctrl_filename, config_tag)
  Call config_display(sorbates)

  !----------------------------------------------------------
  ! Initialize the forcefield parameters
  !----------------------------------------------------------  
  Call forcefield_init(ctrl_filename, scell)
  Call forcefield_display(sorbates, scell)

  !----------------------------------------------------------
  ! Initialize the move type information
  !----------------------------------------------------------
  Call mepbaker_init(mepparams, sorbates, scell, ctrl_filename)
  Call mepbaker_display(mepparams, 0)

  !----------------------------------------------------------
  ! Start the simulation
  !----------------------------------------------------------
  Write(*,'(///,a)') dashedline
  Write(*,'(a,/,a//)') &
      "Beginning Simulation of :", Trim(general_getfiledesc())

  restartfile = &
      genfilename(genparams%restartfile, genparams%simstart)
  Call file_settag(restartfile, d_res_file)
  
!!$  Call mepbaker_FiniteDiff(mepparams, sorbates, scell, 1.0e-7_RDbl)
  Call mepbaker_dosim(mepparams, sorbates, scell)


  !---------------------------------------
  ! Get the elapsed time
  !---------------------------------------
  Call wtime_stoptime(ROUT_MAIN)
  Call wtime_display("Main", ROUT_MAIN)

End Program music
