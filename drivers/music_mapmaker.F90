!-----------------------------------------------
! This contains the main program for generating maps
! 
!-----------------------------------------------

Program music
  Use atom, Only: atom_init, atom_display, atom_samplecf
  Use commandline, Only: commandline_init
  Use config, Only: config_init, config_display, AtMolCoords, config_samplecf
  Use defaults, Only: strLen, RDbl, d_ctrl_file, dashedline
  Use file, Only: file_settag
  Use general, Only: general_init, general_samplecf
  Use interact, Only: Interaction_Model, interact_init, interact_display
  Use mapmaker, Only: MapInfo, mapmaker_idstring, mapmaker_display, &
      mapmaker_init, mapmaker_makemap
  Use molecules, Only: molecules_init, molecules_display, &
      molecules_dump2xyz, molecules_samplecf
  Use simcell, Only: simcell_init, simcell_display, simcell_getell, &
      SimCell_Params, simcell_samplecf
  Use utils, Only: allocErrDisplay
  Use vector, Only: VecType
  Use visxyz, Only:
  Use wtime, Only:

  Implicit None

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** Tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action


  !** Simulation Cell Variable(s)
  Character(len=strLen),Parameter :: simcell_tag="Simulation Cell Information"
  Type(SimCell_Params)            :: scell

  !** Configuration Variable(s)
  Type(AtMolCoords), Dimension(:), Pointer :: sorbates
  Character(len=strLen), Parameter         :: config_tag = & 
      "Configuration Initialization"

  !** Map making parameters
  Type(MapInfo) :: ptmaps  ! Pretabulated maps variable
  Character(len=strLen) :: mm_tag = mapmaker_idstring
  Integer        :: nsims, simno, iter, restartunitno,content_tag, error
  Character(len=strLen)  :: restartfile
  Character(len=strlen)  :: datfilename

  !** Temporary variables for testing the energies
  Integer                :: sorb1,molec1,sorb2,i,j,k,nspecies
  Logical                :: fast, mapflag
  Real(kind=RDbl)        :: pot
  Type(VecType), Dimension(100)   :: amolec1
  Type(VecType), Dimension(1,100) :: asorb2
  Character(len=strLen)           :: comment,filename

  !** Checking the maps vs old ones
  Real(Kind=RDbl), Dimension(3) :: ell


  !** Interaction/Forcefield Parameters and Statistics
  Type(Interaction_Model)                  :: imodel

  !-------------------------------------------------------------
  ! Do basic initialization, based on command line input
  !-------------------------------------------------------------
  Call commandline_init(ctrl_filename,action)

  If (action=="WriteSampleCtrlfile") Then
    Write(0,*) "Incomplete control file"
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
  ! Initialize the general parameters
  !----------------------------------------------
  Call general_init(ctrl_filename)

  !----------------------------------------------------------------------------
  ! Initialize atoms, molecules, zeolite
  !----------------------------------------------------------------------------
  Call atom_init()
  Call atom_display(6)

  nspecies=molecules_init()

  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,2,6)

  !** Allocate the species structure
  Allocate(sorbates(nspecies), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'species structure')    

  Call molecules_display(6)
  Call molecules_dump2xyz(2,"sili.xyz")


  !----------------------------------------------------------
  ! Initialize the configuration
  !----------------------------------------------------------
  Call config_init(sorbates, scell, ctrl_filename, config_tag)

  !----------------------------------------------------------
  ! Initialize the forcefield information
  !----------------------------------------------------------
  Call interact_init(imodel,ctrl_filename,scell,sorbates,'MC')
  Call interact_display(imodel,2,6)
!!$
!!$  Call forcefield_init(ctrl_filename, scell)
!!$  Call forcefield_display(sorbates, scell)

  !----------------------------------------------------------
  ! Initialize the map making parameters
  !----------------------------------------------------------
  Call mapmaker_init(ptmaps,imodel, scell,sorbates,ctrl_filename,mm_tag)
  Call mapmaker_display(ptmaps,6)

  Call config_display(sorbates,6)

  !-----------------------------------------------------------
  ! Be brave, try to make a map
  !-----------------------------------------------------------
  Call mapmaker_makemap(ptmaps%maps(1),scell,sorbates)

!!$  !-----------------------------------------------------------
!!$  ! Compare this one to an old one
!!$  !-----------------------------------------------------------
!!$  ell = simcell_getell(scell,.True.)
!!$  Call maps_pmapinit(map2,ell,'sili.Methane.pmap')
!!$  Call maps_pmapinit(map1,ell,'sili.Methane.LJ.map')
!!$  Call maps_compare(map1,map2,0.10,6)

End Program music




