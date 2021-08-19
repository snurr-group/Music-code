!-----------------------------------------------
! This contains the main program for calculating the configuration
! Integral of a chain molecules in Ideal gas phase at a set of temperatures  
!-----------------------------------------------

Program music


#ifdef USEONLY
  Use utils,only:
  Use defaults, Only: RDbl, strLen, dashedline, ROUT_MAIN, d_ctrl_file
  Use file, Only: file_settag
  Use simcell, Only: SimCell_Params, simcell_init, simcell_display
  Use atom, Only: atom_init, atom_display
  Use molecules, Only: molecules_init, molecules_display,&
      molecules_displaynrg
  Use wtime,Only:wtime_init,wtime_starttime,wtime_stoptime,wtime_display
  Use config,Only:AtMolCoords,config_init,config_display,&
      config_writerestartfile
  Use confinteg, Only : CONFINTEG_Params,confinteg_init, confinteg_getInteg, &
      confinteg_initdisplay
  Use forcefield, Only : forcefield_init,forcefield_display
#else
  Use utils
  Use defaults
  Use file
  Use simcell
  Use atom
  Use molecules
  Use wtime
  Use config
#endif

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

  !** GCMC Simulation Parameters
  Type(CONFINTEG_Params)                   :: integparams
  Character(len=strLen), Parameter         :: confinteg_tag = &
      "Configuration Integration Info"

  Integer        :: nsims, simno, iter, restartunitno,content_tag
  Character(len=strlen)  :: datfilename

  !** Temporary variables for testing the energies
  Integer                :: sorb1,molec1,sorb2,i,j,k
  Logical                :: fast, mapflag
  Real(kind=RDbl)        :: pot
  Character(len=strLen)           :: comment,filename

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
  Call file_settag(ctrl_filename, d_ctrl_file)

  !----------------------------------------------------------------------------
  ! Initialize atoms, molecules, zeolite
  !----------------------------------------------------------------------------
  Call atom_init()
  Call atom_display(6)
  Write(*,*) " Finished atoms , Going to Initialise molecules "

  Call molecules_init()
  Call molecules_display(6)
  Write(*,*) " Finished molecules , Going to Initialise zeolites  "

  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,6)
  Write(*,*) " Finished simcell, Going to Initialise config  "

  !----------------------------------------------------------
  ! Initialize the configuration
  !----------------------------------------------------------
  Call config_init(sorbates, scell, ctrl_filename, config_tag)
  Call config_display(sorbates)
  Write(*,*) " Finished config, Going to Initialise forcefield  "

  !----------------------------------------------------------
  ! Initialize the forcefield information
  !----------------------------------------------------------
  Call forcefield_init(sorbates, scell)
  Call forcefield_display(sorbates, scell)
  Write(*,*) " Finished forcefield, Going to Initialise confinteg section  "

  !----------------------------------------------------------
  ! Initialize the move type information
  !----------------------------------------------------------
  Call confinteg_init(integparams, sorbates, scell, ctrl_filename, &
      confinteg_tag)
  Call confinteg_initdisplay(integparams, 6)
  Write(*,*) " Finished confinteg , Going to do simulation  "

  !----------------------------------------------
  ! Initialize the timer
  !----------------------------------------------
  Call wtime_init()
  Call wtime_starttime(ROUT_MAIN)
 
  Call confinteg_getInteg(sorbates(1),integparams)

  Call wtime_stoptime(ROUT_MAIN)
  Call wtime_display("Main Program", ROUT_MAIN)
End Program music
