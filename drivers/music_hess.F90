!-------------------------------------------------------------------
! This is the driver reading a molecule file, and 
! calculate the hessian and writing hessian to hess.out file
! and writing their nrg into a data file
! hessian is calculated by numerically differentiating the forces
! input  : a control file, initial config in music.xyz
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
  Use molecules, Only: molecules_display, molecules_init, &
      molecules_sampleCF, molecules_getnatoms
  Use datafile, Only: CONFILE, datafile_writeconfig, datafile_initout
  Use simcell, Only: SimCell_Params, simcell_display, simcell_init, &
      simcell_sampleCF
  Use config, Only: config_writerestartfile, config_init, config_sampleCF, &
      AtMolCoords, config_display, config_dump, config_displaycoords, &
      config_changerp, config_xyz2config
  Use interact, Only: Interaction_Model, interact_init, interact_display
  !  Use movie, Only: MovieInfo, movie_init, movie_display, movie_makemovie
  Use stopcriteria, Only: STOP_CHECK, stopcriteria_init, stopcriteria_check
  Use ffcall, Only: ffcall_nrg, ffcall_force
  Implicit None

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  Character(len=strLen) :: xyzfile="music.xyz"

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
  !  Type(MDInfo) :: mdParams
  !  Character(len=strLen), Parameter :: md_tag = "MD Information"

  !** Movie Parameters
  !  Type(MovieInfo)                  :: mainMovieParams
  !  Character(len=strLen), Parameter :: mainMovie_tag = "Movie Information"

  !** Configuration file
  Type(CONFILE) :: configfile

  !** MD Simulation time
  Real(Kind=RDbl)            :: mtime, newx, oldx, dx

  Real(Kind=RDbl)            :: intranrg, conv_fac

  !** set the default display unit
  Integer, Parameter         :: dUnit = 6

  Integer                    :: nsims, simno, iter, error, nmols, nconfs, k
  Integer                    :: configunitno, restartunitno
  Integer :: natoms, a, i, j, s, m, ncols, index
  Character(len=strLen)      :: comment, filename, restartfile

  Logical :: success=.False.

  Integer , Parameter :: AMAX=50
  Real(kind=RDbl), Parameter :: DXDEFAULT=0.01_RDbl
  Real(kind=RDbl)::H(AMAX*3,AMAX*3)
  Real(kind=RDbl)::F0(AMAX*3), F1(AMAX*3), F2(AMAX*3)

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
  Call config_xyz2config(species,scell,xyzfile,success)

  Call interact_init(imodel,ctrl_filename,scell,species,'MD')
  Call interact_display(imodel,2,dUnit)

  !----------------------------------------------------------
  ! Initialize the move type information
  !----------------------------------------------------------
  Call stopcriteria_init(ctrl_filename)

  simno = 1

  natoms=molecules_getnatoms(1)
  index=0
  dx=DXDEFAULT
  Do a=1, natoms
    Do j=1,3

      index=index+1
      oldx=species(1)%coords(a,1)%rp%comp(j)
      newx=oldx+dx
      species(1)%coords(a,1)%r%comp(j)=newx
      species(1)%coords(a,1)%rp%comp(j)=newx
      Call ffcall_force(imodel,scell,species,F2)
      Call ffcall_nrg(imodel,scell,species,"INTRA",intranrg)
      
      newx=oldx-dx
      species(1)%coords(a,1)%rp%comp(j)=newx
      species(1)%coords(a,1)%r%comp(j)=newx
      Call ffcall_force(imodel,scell,species,F0)


      species(1)%coords(a,1)%r%comp(j)=oldx
      species(1)%coords(a,1)%rp%comp(j)=oldx    
      Call ffcall_force(imodel,scell,species,F1)
      Call ffcall_nrg(imodel,scell,species,"INTRA",intranrg)


      Do k=1,natoms*3
        ! F is -ve of gradient (force)
        ! d2V/dxdxy = d(-fx)/dy
        H(k,index)= (F0(k)-F2(k))/(2*dx)
      End Do

    End Do
  End Do

  ! convert to hartree/bohr/bohr from kcal/ang/ang
  conv_fac=(one/627.51)*(0.529*0.529)
  Do i=1,natoms*3
    Do j=1,natoms*3
      H(i,j)=H(i,j)*conv_fac 
    End Do
  End Do

  Write(*,*) "Force constants in Cartesian coordinates: "
  ! write to out put
  ! write 5 per column
  !
  Do s=1,(natoms*3/5)+1
  ! loop over group of columns


    index=(s-1)*5+1
    Write(*,'(3x,5i14)') index, index+1, index+2, index+3, index+4
    If (s==(natoms*3/5+1)) Then
      ! decide how many columns in last group
      ncols=3*natoms-index+1
    Else
      ncols=5
    Endif
!    Do m=1,natoms*3 ! full matrix
    Do m=index,natoms*3 ! gaussian style,
    Select Case(ncols)
    Case(1)
      Write(*,'(i3,e16.5)') m, H(m,index)
    Case(2)
      Write(*,'(i3,2e16.5)') m, H(m,index), H(m,index+1)
    Case(3)
      Write(*,'(i3,3e16.5)') m, H(m,index), H(m,index+1), H(m,index+2)
    Case(4)
      Write(*,'(i3,4e16.5)') m, H(m,index), H(m,index+1), H(m,index+2), &
          H(m,index+3)
    Case(5)
      Write(*,'(i3,5e14.5)') m, H(m,index), H(m,index+1), H(m,index+2), &
          H(m,index+3), H(m,index+4)
    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select


    End Do
  End Do
  
  !---------------------------------------
  ! Get the elapsed time
  !---------------------------------------
  !  Call wtime_stoptime(MAX_ROUTINES)
  !  Call wtime_display("Main Program", MAX_ROUTINES)

End Program music
