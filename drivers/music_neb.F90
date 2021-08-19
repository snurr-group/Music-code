!-------------------------------------------------------------------
! This contains the main program for the Nudged Elastic Band Method
!-------------------------------------------------------------------

Program music

  Use defaults, Only: strLen, RDbl, d_ctrl_file, d_res_file, pcalls, &
      dashedline, MAX_ROUTINES
  Use utils, Only: genfilename, int2str, real2str, allocerrdisplay, ToUpper
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
  Use config, Only: AtMolCoords, config_init, config_sampleCF, &
      config_display, config_dump, config_displaycoords, &
      config_writerestartfile, config_config2xyz, config_config2struc
  Use neb, Only: NEB_Info,neb_initbasic,neb_initrest,neb_nimages,neb_configtag, &
      neb_movechain,neb_images2xyz,neb_sampleCF,neb_display,neb_clean, &
      neb_displaystats, neb_modimages2xyz, neb_modimage2xyz
  Use interact, Only: Interaction_Model, interact_init, interact_display
  Use movie, Only: MovieInfo, movie_init, movie_display
  Use readstruc, Only: Structure,readstruc_visxyz,readstruc_natoms, &
      readstruc_recenter

  Implicit None

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** Tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action

  !** Simulation Cell Variable(s)
  Character(len=strLen),Parameter :: simcell_tag="Simulation Cell Information"
  Type(SimCell_Params)            :: scell

  !** NEB parameters
  Integer                         :: nimages
  Type(NEB_Info)                  :: nebinfo

  !** Configuration Variable(s)
  Type(AtMolCoords), Dimension(:,:), Pointer :: image
  Character(len=strLen)                      :: config_tag

  !** Interaction/Forcefield Parameters and Statistics
  Type(Interaction_Model)                  :: imodel

  !** Movie Parameters
  Type(MovieInfo)                  :: mainMovieParams
  Character(len=strLen), Parameter :: mainMovie_tag = "Movie Information"

  !** Configuration file
  Type(CONFILE) :: configfile

  !** set the default display unit
  Integer, Parameter         :: dUnit = 6

  Integer                    :: i, nsims, simno, iter, error, nmols
  Integer                    :: length
  Character(len=strLen)      :: comment, filename
  Character(len=strLen)      :: string, string1, string2
  Character(len=strLen), Dimension(:), Allocatable  :: restartfile

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
    Call neb_sampleCF(0)
    Write(0,'(a)') dashedline
    Stop
  Else If (action /= "DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  End If
      
  !----------------------------------------------
  ! Welcome the user 
  !----------------------------------------------
  Write(dUnit,'(a)') "Welcome to Music (NEB implementation)"
  Write(dUnit,'(2a)') "Reading from control file ",Trim(ctrl_filename)
  
  !----------------------------------------------
  ! Basic Initialization
  !----------------------------------------------
  Call general_init(ctrl_filename)

  Call file_settag(ctrl_filename, d_ctrl_file)

  !** Initialize atoms 
  Call atom_init()
  Call atom_display(dUnit)

  !** Initialize molecules
  nmols = molecules_init()
  Call molecules_display(dUnit)

  !** Initialize simulation cell
  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,2,dUnit)

  !** Trim the name of the restart file if required
  length = Len(Trim(genparams%restartfile))
  If (genparams%restartfile(length-3:length) == '.res') Then
    genparams%restartfile = genparams%restartfile(1:length-4)
  End If

  !----------------------------------------------
  ! NEB-Specific Initialization
  !----------------------------------------------

  !** Read the NEB section of the control file
  Call neb_initbasic(nebinfo,ctrl_filename)

  !** Allocate the space for the images
  nimages = neb_nimages(nebinfo)
  Allocate(image(nimages,nmols), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'images structure')

  !** Initialize the start configuration
  config_tag = neb_configtag(nebinfo,'START')
  Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
  Write(*,*) config_tag
  Call config_init(image(1,:), scell, ctrl_filename, config_tag)
  Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
  Call config_display(image(1,:),dUnit)

  !** Initialize the end configuration
  config_tag = neb_configtag(nebinfo,'END')
  Call config_init(image(nimages,:), scell, ctrl_filename, config_tag)
  Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
  Call config_display(image(nimages,:),dUnit)

  !** Turn the intermediate configuration init flags off (pseudo-hack)
  Do i = 2,nimages-1
    image(i,:)%isinit = .False.
  End Do

  !** Initialize the intra-image interactions
  Call interact_init(imodel,ctrl_filename,scell,image(1,:),'MD')
  Call interact_display(imodel,2,dUnit)

  !** Initialize the remainder of the NEB parameters
  Call neb_initrest(nebinfo,image,imodel,scell)
  Call neb_display(nebinfo,2,dUnit)

  !** Dump the initial chain to a movie file for visualization
  Call neb_images2xyz(nebinfo,image,'chain0.xyz')

  !** Create the restart filenames and write the first files
  Allocate(restartfile(nimages), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'restart file names')
  Do i = 1,nimages
    string = int2str(i)
    Write(restartfile(i),'(3a)') Trim(genparams%restartfile),Trim(string),'.res'
    Call config_writerestartfile(image(i,:),restartfile(i))
  End Do 

  !---------------------------------------
  ! Run the simulation
  !---------------------------------------

  Do iter = 1,genparams%niterations
    !** Move the chain
    Call neb_movechain(nebinfo,image,imodel,scell,2,dunit)

    !** Dump to the crash file
    If (Mod(iter, genparams%icrash) == 0) Then
      Do i = 2,nimages-1
        Call config_writerestartfile(image(i,:),restartfile(i))
      End Do
    End If

    !** Display to screen if desired
    If (Mod(iter,genparams%iprint) == 0) Then
      string1 = int2str(iter)
      string2 = int2str(genparams%niterations)
      Write(dunit,'(2x,4a)') 'Iteration number ',Trim(string1), &
          ' of ',Trim(string2)
      Call neb_displaystats(nebinfo,image,2,dunit)
      Write(dunit,*) 
    End If

    !** Make a movie of each image
    Do i = 2,nimages-1
      string = int2str(i)
      Write(filename,'(3a)') 'image',Trim(string),'.xyz'
      Write(comment,'(a,i4)') 'Iteration number: ',iter
      Call config_config2xyz(image(i,:),filename,.False.,comment)
!      Call neb_modimage2xyz(nebinfo,image,i,filename,comment)  !** with tangents
    End Do
  End Do

  !** Make a movie of the final chain
  string = int2str(genparams%niterations)
  Write(filename,'(3a)') 'chain',Trim(string),'.xyz'
  Call neb_images2xyz(nebinfo,image,filename)

  !** Make a movie of the final chain WITH tangents included
  Write(filename,'(3a)') 'endchain_wtangents.xyz'
  Call neb_modimages2xyz(nebinfo,image,filename)

  !---------------------------------------
  ! Get the elapsed time
  !---------------------------------------
  Call wtime_stoptime(MAX_ROUTINES)
  Call wtime_display("Main Program", MAX_ROUTINES)

End Program music
