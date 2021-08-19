!-----------------------------------------------
! This contains the main program that will generate a restart file 
!that can be used by an Md Run. This basically does gcmc.
!This starts from an initial configuration and generates a restart file 
!containing the required number of molecules ( "Nmolec") of type "NAMEmolec"  
!A restart file is written to "stopfile" which can be used by further programs
!-----------------------------------------------

Program music
  Use file
  Use defaults
  Use utils
  Use wtime
  Use general
  Use atom
  Use molecules
  Use simcell
  Use forcefield
  Use config
  Use gcmc
  Use nvtmc
  Use visxyz
  Use ssinteraction
  Use datafile

  Implicit None
  !** The  config files parameters
   Type (CONFILE) :: conf_file

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
  Type(GCMC_Params) :: gcmcparams
  Character(len=strLen), Parameter         :: gcmc_tag = &
      "GCMC Information"

  !** NVTMC Simulation Parameters
  Type(NVTMC_Params) :: nvtmcparams
  Character(len=strLen), Parameter         :: nvtmc_tag = &
      "NVTMC Information"
  
  
  !*** Parameters for stopping the program at the rquired point
  Integer,              Parameter      :: Nmolec=32
  Character(len=strLen),Parameter :: NAMEmolec="Ethane"
  Character(len=strLen), Parameter         :: stopfile = "stopfile"
  
  Integer        :: nsims, simno, iter, restartunitno,content_tag
  Character(len=strLen)  :: restartfile

  !** Temporary variables for testing the energies
  Integer                :: i,j,k,typemol,currentNmolec
  Logical                :: mapflag
  Real(kind=RDbl)        :: pot

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
  Call file_settag(ctrl_filename,_ctrl_file)

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
  Call zeolite_init()

  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,6)

  !----------------------------------------------------------
  ! Initialize the configuration
  !----------------------------------------------------------
  Call config_init(sorbates, scell, ctrl_filename, config_tag)
  Call display(sorbates)

  !----------------------------------------------------------
  ! Initialize the forcefield information
  !----------------------------------------------------------
  Call forcefield_init(ctrl_filename, scell)
  Call forcefield_display(sorbates, scell)

  !----------------------------------------------------------
  ! Initialize the move type information
  !----------------------------------------------------------
  Call gcmc_init(gcmcparams, sorbates, scell, ctrl_filename, gcmc_tag)
  Call display(gcmcparams, 0)
  
  nsims = gcmc_getnosims(gcmcparams)
  Call nvtmc_init(nvtmcparams, ctrl_filename, nvtmc_tag)
  Call display(nvtmcparams, 0)
 
 Do simno=1, nsims

    !** Open the configuration file for the new simulation and close
    !** the previous one if any
    conf_file%name=genfilename(genparams%configfile,&
        simno+genparams%simstart-1)
    content_tag=general_getContenTtag()
    Call datafile_initout(conf_file,ctrl_filename,simno,content_tag)
    restartfile = genfilename(genparams%restartfile,simno+genparams%simstart-1)

    Call file_settag(restartfile, _res_file)

    !** Display the parameters for this simulation
    Call gcmc_displaysimparams(gcmcparams, simno, 0)

    Call random_init(400913420)
    Call random_display()
    
    Do iter = 1, genparams%niterations

      genparams%currentiteration = iter

      !** Do the NVTMC and GCMC simulations
      Call nvtmc_dosim(nvtmcparams, sorbates, scell, 1)
      Call gcmc_dosim(gcmcparams, sorbates, scell, simno)
      
      !** Dump to the standard output
      If (Mod(iter, genparams%iprint) == 0) Then
        Write(*,'(a, i4)') "Simulation Number :", simno
        Write(*,'(a, i9)') "Iteration Number  :", iter
        Call gcmc_displaystats(gcmcparams, 0)
        Call nvtmc_displaystats(nvtmcparams, 1, 0)
        Write(*,'(//)') ! Put some spaces between outputs
      End If

      !** Dump to the config file
      If (Mod(iter, genparams%iconfig) == 0) Then
        Call datafile_writeconfig(sorbates, conf_file)
      End If

      !** Dump to the crash file
      If (Mod(iter, genparams%icrash) == 0) Then
        Call config_writerestartfile(sorbates, restartfile)
      End If
      
      !*** checks for number of molecs of "NAMEmolec"
      typemol= molecules_gettype(Trim(NAMEmolec))
      currentNmolec=config_getnmoles_array(sorbates,typemol)

      
      If (currentNmolec==Nmolec) Then
        Write(*,*) " stoping the program at nmoles = ",Nmolec," of ", &
            Trim(NAMEmolec)
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Call   gcmc_stop(sorbates,Trim(stopfile),"YES")
      Endif
      
  End Do
  
    Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Call molecules_displaynrg()

    Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Call forcefield_getssint(sorbates, scell, .False., pot, mapflag)
    Call molecules_displaynrg()
    
  
  End Do
End Program music


