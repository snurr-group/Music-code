! This is the driver file for reading config files and doing averages
! and any other post run calculations. It can also be used while the run 
! is still in progress
Program music
  Use atom, Only:
  Use commandline, Only: commandline_init
  Use config, Only: AtMolCoords, config_display
  Use datafile, Only: CONFILE, datafile_initin, datafile_readconfig
  Use defaults, Only: MAX_SORBS, strlen, RDbl, dashedline
  Use file, Only: file_open, file_close
  Use general, Only: general_getnoofconfigs, genparams
  Use molecules, Only: molecules_display
  Use post, Only: PostInfo, post_init, post_getnfiles, post_display, &
      post_getConfigFilename, post_reinitStats, &
      post_writeAvgstoFile, post_readandavg, post_samplecf, post_calcMSD, &
      post_initconfileReading, post_calcRadProf
  Use simcell, Only: SimCell_Params
  Use stats, Only:
  Use utils, Only: genfilename


  Implicit None

  Type (CONFILE)::configfile
  Integer :: nfile, simno, tempunitnum

  !** file which contains all informations about what to do 
  Character(len=strLen) :: postfile
  Type(PostInfo) :: postParams

  !** action based on command line input .Ex:  run ?, help?, write sample file?
  Character(len=strLen) :: action

  !-------------------------------------------------------------
  ! Do basic initialization, based on command line input
  !-------------------------------------------------------------
  Call commandline_init(postfile,action)

  If (action=="WriteSampleCtrlfile") Then
    Write(0,'(a)') dashedline
    Call post_sampleCF(0)
    Stop
  Elseif(action/="DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Endif

  !** Initialize the general post stuff, which includes initializing
  !** the atoms, molecules, simcell, configuration, and forcefield
  !** Also initialize the analysis-types.
  Call post_init(postParams, postfile )
  Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
      "Finished Reading the post code control file"

  Call post_display(postParams, 6)


  !** Read data given in each file and do the analysis
  nfile = post_getnfiles(postParams )
  
  !** Display some information
  Call molecules_display(6)
  Call config_display(postParams%species,6)

  Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
      "Going to read each datafile and average"

  !** For each file read through all the available data, do ensemble averages
  Do simno=1,nfile

    !** (Re)Initialize all statistics variables in the postParams 
    Call post_reInitStats(postParams)


    ! The datafiles have been already intialized in post.F90
    ! we just need to make sure that the datafile corresponding to this simno
    ! is ready to be read from
    Call post_initconfileReading(postParams, simno)

    Call post_readAndAvg( postParams, simno )
    Write(*,*) " Finished reading thru the config file "

    !** Write whatever averages were calculated for this config file
    Write(*,*) " Writing the results "
    Call post_writeAvgstoFile( postParams, simno )
    
  End Do

  !** Calculate the diffusivities from and MD simulation if asked for
  If (Associated(postParams%MSD)) Then
    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__," Calculating MSD"  
    Call post_calcMSD(postParams)
    ! Calculate the radial profiles from an MD simulation if asked for
    If (Associated(postParams%MSD%mdpp%radial)) Then
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,&
          " Calculating Radial Profiles"  
      Call post_calcRadProf(postParams)
    Endif
  End If


  Write(*,*) 
  Write(*,*) 'Finished post-processing'

End Program music











