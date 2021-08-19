!----------------------------------------------------------------------
! This module implements the different kinds of moves that can be done
! on united atom alkane molecules. As of now, the abillity to deal
! with branched molecules is very limited. Most of the routines work
! only with linear alkanes - Shaji 20/Dec/2002
!  
! NOTES on various issues added by shaji, these were added when 
! they were noticed and will be continuously updated
! - energy return policy : any move routine should finally return the energy
!   my making a subinteract_int call
! - use of %max_nrg : the final subinteract_int call from a move routine 
!   may be done with a higher %max_nrg to avoid unnecessary 
!   rejection/Stop especially While dealing With large alkanes 
!----------------------------------------------------------------------
Module brmoves
  Use angle, Only: Angle_Library_Params, angle_getbiasindex, angle_CountLib, &
      angle_indexFromVal,angle_init
  Use auxmoveparams, Only: AuxMoveObjects 
  Use branchedcoords,Only : NodeType, branchedcoords_toxyz, &
      branchedcoords_place, branchedcoords_placenode, branchedcoords_unplace, &
      branchedcoords_update_ta, branchedcoords_findnode, &
      branchedcoords_markhigheratoms, branchedcoords_copyvalues
  Use bmap,Only     : BiasMap_Params,bmap_init,bmap_getbiasindex,&
      bmap_getncubelets,bmap_getbiaswt,bmap_getbiaspt,bmap_getcellwt
  Use bbmodel, Only: BENDING_KEY
  Use cavitylist, Only: Cavity_Params, cavitylist_checkcavity, &
      cavitylist_updatemolecinfo, cavitylist_checkandIncrMemory, &
      cavitylist_getcubeindices
  Use config, Only : AtMolCoords, config_getnatoms, config_allocfields, &
      config_copymolec, config_delmol, config_setnmoles, config_checkandincr, &
      config_getnmoles
  Use defaults, Only: RDbl, strLen, lstrLen, pi,twopi, zero, one, &
      scalepe, kcalmole_kb, degTorad, NO_OF_INTRA_POTS,default_MAX_EXP,&
      TOTAL_INDEX,MAX_SORBS, MAX_ATOMS, d_ctrl_file, dbgflag
  Use file,Only     : file_getunit,file_open,file_gettype
  Use frame, Only: Frame_Library_Params, frame_getbiasfactor, & 
      frame_CountLib, frame_init 
  Use general, only : general_getCurrentSimno
  Use ipmodel, Only: INTRAPAIR_KEY
  Use interact, Only: Interaction_Model, interact_listparams, &
      interact_simpleint
  Use matrix,Only   : MatrixType,matrix_display,Operator(*)
  Use molecule, Only: MolecularParams 
  Use molecules,Only: molecules_getnatoms, molecules_getcomdefn, & 
      molecules_getpointer, molecules_getnsorbs 
  Use random,Only   : rranf, random_gettrial 
  Use simcell, Only : SimCell_Params, simcell_getzeoell, simcell_geteff, &
      simcell_pbc, simcell_getmolectype, &
      simcell_maptouc, simcell_getcellorigin, simcell_uniformRandPt
  Use storebase, Only: EnergyPlus, storebase_display, storebase_init 
  Use subinteract, Only: Subset_Interactions,subinteract_int, &
      subinteract_changenmoles
  Use tormodel, Only: TORSION_KEY
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, swap , &
      toint, tolower, toupper, allocErrDisplay, toreal, cleanstring, &
      checkandstop
  Use vector, Only : VecType, IntVecType, Assignment(=), Operator(+), &
      Operator(-), Operator(*), Operator(/), vector_getdistsq, &
      vector_getnormsq, vector_filedisplay

  Implicit None
  Save

  Private
  Public :: &  ! params
      BranchedFTrans_Params,BranchedRTRans_Params, &
      BranchedRRot_Params,  BranchedRIns_Params, BranchedCBIns_Params, &
      BranchedCBDel_Params , BranchedRDel_Params, BranchedReGrow_Params,  &
      BranchedIDFlip_Params, &
      !--------- inits
  brmoves_FTransinit, brmoves_RTransinit, brmoves_RRotInit, &
      brmoves_RDelInit, brmoves_reGrowInit,       brmoves_CBInsinit, &
      brmoves_CBDelInit, brmoves_RInsinit,  brmoves_IDFlipinit,  &
      ! init-auxiliary
  brmoves_checkOppMoves, brmoves_initTrialArrays,  &
      !--------- moves
  brmoves_randomtranslate, brmoves_fixedtranslate, &
      brmoves_randomrotate, brmoves_rdelete, brmoves_cbdelete , &
      brmoves_reGrow, brmoves_cbinsert, brmoves_rinsert, brmoves_idflip,  &
      brmoves_adjustdeltatrans, brmoves_adjustdeltarot, brmoves_postadjust,  &
      !---------- display
  brmoves_displayparams ! , &


  !** Branched gc random length Rotation params
  Type BranchedRRot_Params
    Logical           :: scale_rot
    Real(kind=RDbl)   :: deltarot
    Type(VecType), Dimension(:), Pointer   :: com_defn
  End Type BranchedRRot_Params

  !** Branched gc fixed length translation params
  Type BranchedFTrans_Params
    Logical           :: constrained
    Type(MatrixType)  :: constrtrans ! Transforming to the consstr system.
    Type(MatrixType)  :: constrinvtrans
    Type(VecType)     :: normal
    Logical           :: scale_trans
    Real(kind=RDbl)   :: deltatrans
    Type(VecType), Dimension(:), Pointer   :: com_defn
  End Type BranchedFTrans_Params

  !** Branched gc random length translation params
  Type BranchedRTrans_Params
    Logical           :: constrained
    Type(MatrixType)  :: constrtrans ! Transforming to the consstr system.
    Type(MatrixType)  :: constrinvtrans
    Type(VecType)     :: normal
    Logical           :: scale_trans
    Real(kind=RDbl)   :: deltatrans
    Real(kind=RDbl)   :: max_nrg
    Type(SimCell_Params), Pointer :: scell
    Type(VecType), Dimension(:), Pointer   :: com_defn
  End Type BranchedRTrans_Params

  !** Branched gc random Insertion params
  Type BranchedRIns_Params
    Type(SimCell_Params), Pointer :: scell
    Type(VecType), Dimension(:), Pointer   :: com_defn
  End Type BranchedRIns_Params

  !** Branched gc biased Insertion params
  Type BranchedCBIns_Params
    Type(BiasMap_Params), Pointer :: bmap
    Character(len=strLen)         :: biasfilename
    Type(SimCell_Params), Pointer :: scell

    Real(kind=RDbl)   :: rti, tk  

    Real(kind=RDbl) :: max_nrg
    !** max_nrg to be set from main progfram for each run at the beginning 
    !** of the run

  End Type BranchedCBIns_Params


  !** Branched gc biased Insertion params
  Type BranchedIDFlip_Params
    Type(BiasMap_Params), Pointer :: bmap
    Character(len=strLen)         :: biasfilename
    Type(SimCell_Params), Pointer :: scell

    Real(kind=RDbl)   :: rti, tk  

    Real(kind=RDbl) :: max_nrg
    !** max_nrg to be set from main progfram for each run at the beginning 
    !** of the run

  End Type BranchedIDFlip_Params


  !** Branched gc cut and regrow params
  Type BranchedReGrow_Params
    Type(BiasMap_Params), Pointer :: bmap
    Character(len=strLen)         :: biasfilename
    Type(SimCell_Params), Pointer :: scell

    ! for storing the molecule during deletion 
    Type(AtMolCoords)         :: oldcoords,suggested

    Real(kind=RDbl)   :: rti , tk 

    Real(kind=RDbl) :: max_nrg

  End Type BranchedReGrow_Params

  !** Branched gc biased deletion params
  ! opposite move for BIns_Params
  Type BranchedCBDel_Params
    Type(BiasMap_Params), Pointer :: bmap
    Character(len=strLen)         :: biasfilename
    Type(SimCell_Params), Pointer :: scell
    Real(kind=RDbl)   :: rti , tk 
    Real(kind=RDbl) :: max_nrg
  End Type BranchedCBDel_Params

  Integer, Parameter   :: COS_EXP_MAX_PARAMS = 6   
  Real(kind=RDbl), Parameter   :: CB_IP_HICUT = 13.0_RDbl 
  Real(kind=RDbl), Parameter   :: CB_IP_LOCUT = 0.01_RDbl
  Real(kind=RDbl), Parameter   :: CB_IP_HC2 =CB_IP_HICUT *CB_IP_HICUT 
  Real(kind=RDbl), Parameter   :: CB_IP_LC2 =CB_IP_LOCUT *CB_IP_LOCUT

  ! these are the libraries used during insertion/deletion/regrow
  Type AuxLibraries
    Logical :: initialized = .False.
    Integer                              ::  sorbtype
    Integer                              ::  ntors, natoms, nangles, nframes
    Real(kind=RDbl) , Dimension(:,:),Pointer  ::  tor
    Type(Frame_Library_Params), Dimension(:), Pointer ::  frame_lib
    Type(Angle_Library_Params), Dimension(:), Pointer ::  angle_lib
    Integer, Dimension(:), Pointer :: angle_lib_list, tor_list
    ! these lists return the angle_library for a given atom

    Real(kind=RDbl) , Dimension(:,:),Pointer  ::  ip_ABvals
    ! first dimesnion is num of ip potentials, second is 2 : A and B from LJ
    Integer, Dimension(:), Pointer :: n_ip_pots
    Integer, Dimension(:,:), Pointer :: ip_pots_atomlist, ip_pots_index
    ! n_ip_pots = tells how many ippots at each node
    ! ip_pots_atomlist = tells which are those atoms at each node
    ! ip_pots_index = corrspondong index in ip_ABvals array

  End Type AuxLibraries

  ! These arrays are allocated during the first move allocation
  ! They are initialized when the move for that particular sorbate is called
  Logical :: auxlibs_sized = .False.
  Logical, Dimension(:), Pointer :: auxlibs_initialized
  Type(AuxLibraries), Dimension(:), Pointer :: auxlibs_array

  !** Branched gc biased Insertion params
  Type BranchedRDel_Params
    Type(SimCell_Params), Pointer :: scell
  End Type BranchedRDel_Params

  Interface brmoves_displayparams
    Module Procedure brmoves_displayFTrans
    Module Procedure brmoves_displayRTrans
    Module Procedure brmoves_displayRRot
    Module Procedure brmoves_displayCBIns
    Module Procedure brmoves_displayCBDel
    Module Procedure brmoves_displayRegrow
    Module Procedure  brmoves_displayIDFlip
  End Interface

  Interface brmoves_postadjust
    Module Procedure brmoves_postadjustRegrow
    Module Procedure brmoves_postadjustRot
    Module Procedure brmoves_postadjustTrans
    Module Procedure brmoves_postadjustCBIns
    Module Procedure brmoves_postadjustCBDel
  End Interface


  Interface brmoves_checkOppMoves
    Module procedure brmoves_checkInsDel
  End Interface

  !** default value of trials for each of these angles, they are later set 
  !** to the required value in brmoves_inittrialarrays
  Integer :: NPsiTrialsSS   = 8
  Integer :: NThetaTrialsSS = 8
  Integer :: NPhiTrialsDS   = 8
  Integer :: NHLTrials      = 10 

  Integer,Parameter :: DEF_indent=10 ! default display indentation 
  Integer,Parameter :: MAX_TRIAL_NO = 1000, MAX_HLFAS_TRIALS = 1000

  ! to be read from ctrlfile, MAX_STOP_NRG
  Real(kind=RDbl)       :: MAX_STOP_NRG =80.00 ! kcal/mol
  Real(kind=RDbl),Parameter      :: VERYLRGMAXNRG =1.00e5


  !** SSTrialArr(i,:) i=1->costheta, i=2->theta, i=3-> psi
  Real(kind=RDbl),Dimension(3,MAX_TRIAL_NO) :: SSTrialArr
  Real(kind=RDbl),Dimension(MAX_TRIAL_NO)   :: DSTrialArr  
  Real(kind=RDbl)::SSdcos,SSdpsi,DSdphi
  Logical :: trial_arrays_inited = .False.
  Logical :: MAP_EXISTS = .False.
  Integer ::  MAP_SPC=-1

  Type(EnergyPlus) :: ffout


  !** temperory arrays that will be used by the insert/delete routines
  Real(kind=RDbl),Dimension(MAX_TRIAL_NO) :: trial_e,trial_exp_e,trial_intra
  Real(kind=RDbl),Dimension(MAX_TRIAL_NO) :: theta1,phi1,cum_prob
  Real(kind=RDbl),Dimension(MAX_TRIAL_NO,MAX_SORBS) :: trial_ncoul_ar
  Real(kind=RDbl),Dimension(MAX_HLFAS_TRIALS) :: bmap_bias_array
  Real(kind=RDbl),Dimension(MAX_HLFAS_TRIALS,3) :: startposarray

  ! for initializing some specific info
  Character(len=strLen),Parameter :: CBGCMC_STRING= "CBGCMC Specific Info"
  Logical :: firsttime_in_brmoves = .True.


  ! common cavity bias object and some associated variables
  Logical :: cavity_initialized=.False. , avoid_cubes
  Type(Cavity_Params), Pointer :: cavity
  Integer, Dimension(MAX_ATOMS) :: cubes2avoid=0


  ! HLFAS or LLFAS
  Logical :: hlfas=.False., llfas=.True., hlfasON=.False.
  Integer :: hlfas_simno
  Character(len=5),Parameter :: hlfas_string="HLFAS", llfas_string="LLFAS"

Contains

  !-------------------------------------------------------------------------
  ! update cavity information after bdelete
  !-------------------------------------------------------------------------
  Subroutine brmoves_postadjustCBDel(delpar, species, spc, mol)
    Type(BranchedCBDel_Params), Pointer                   :: delpar
    Type(AtMolCoords), Dimension(:), Intent(in)    :: species
    Integer, Intent(in)                               :: spc, mol

    If (cavity_initialized) Then
      Call cavitylist_updateMolecInfo(cavity, species(spc), spc, mol, 110 )
    Endif

  End Subroutine brmoves_postadjustCBDel

  !-------------------------------------------------------------------------
  ! update cavity information after binsert 
  !-------------------------------------------------------------------------
  Subroutine brmoves_postadjustCBIns(inspar, species, spc, mol)
    Type(BranchedCBIns_Params),Pointer                   :: inspar
    Type(AtMolCoords), Dimension(:), Intent(in)      :: species
    Integer, Intent(in)                                :: spc, mol
    If (cavity_initialized) Then
      Call brmoves_postadjustCavity(species, spc, mol)
    Endif
  End Subroutine brmoves_postadjustCBIns


  !-------------------------------------------------------------------------
  ! update cavity information after binsert or linsert.
  ! takes RigidInsBias as input
  !-------------------------------------------------------------------------
  Subroutine brmoves_postadjustCavity(species, spc, mol)
    Type(AtMolCoords), Dimension(:), Intent(in)         :: species
    Integer, Intent(in)                                 :: spc, mol

    !** Also check whether memory in cubelist is alright
    !** has to be done before updatemolecinfo
    Call cavitylist_checkandIncrMemory(cavity, species(spc), spc )

    Call cavitylist_updateMolecInfo(cavity,species(spc), spc, &
        mol, 101 )

  End Subroutine brmoves_postadjustCavity


  !-------------------------------------------------------------------------
  ! update cavity information after rrotate
  !-------------------------------------------------------------------------
  Subroutine brmoves_postadjustRot(rotpar, species, spc, mol)
    Type(BranchedRRot_Params),Pointer                :: rotpar
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: species
    Integer, Intent(in)                              :: spc, mol

    If (cavity_initialized) Then
      Call cavitylist_updateMolecInfo(cavity, species(spc), spc, mol, 111 )
    Endif
  End Subroutine brmoves_postadjustRot

  !-------------------------------------------------------------------------
  ! update cavity information after rtranslate 
  !-------------------------------------------------------------------------
  Subroutine brmoves_postadjustTrans(transpar, species, spc, mol)
    Type(BranchedRTrans_Params),Pointer                :: transpar
    Type(AtMolCoords), Dimension(:), Intent(in)        :: species
    Integer, Intent(in)                                :: spc, mol

    If (cavity_initialized) Then
      Call cavitylist_updateMolecInfo(cavity, species(spc), spc, mol, 111 )
    Endif

  End Subroutine brmoves_postadjustTrans

  !-------------------------------------------------------------------------
  ! update cavity information after rtranslate 
  !-------------------------------------------------------------------------
  Subroutine brmoves_postadjustRegrow(regrowpar, species, spc, mol)
    Type(BranchedRegrow_Params),Pointer                :: regrowpar
    Type(AtMolCoords), Dimension(:), Intent(in)        :: species
    Integer, Intent(in)                                :: spc, mol

    If (cavity_initialized) Then
      Call cavitylist_updateMolecInfo(cavity, species(spc), spc, mol, 111 )
    Endif

  End Subroutine brmoves_postadjustRegrow



  !------------------------------------------------------------------------
  ! Initializes temperory array used in this module
  !------------------------------------------------------------------------
  Subroutine brmoves_basicInit( scell, indent)
    Type(SimCell_Params), Intent(IN)             :: scell
    Integer, Intent(in) :: indent
    Integer                   :: unitno, lineno , NPsiSS, NthetaSS, NPhiDS

    Integer                   :: disk_samp_num, sphere_samp_num, nfields
    Character(len=strLen)     :: line
    Character(len=strLen),dimension(strLen)     :: fields
    Character(len=indent)     :: blank

    Real(kind=RDbl) :: new_max_nrg
    Logical :: nocbgcmap, problem

    firsttime_in_brmoves=.False.
    blank=Repeat(" ",indent)
    Write(*,'(3a,i4)') blank,__FILE__," : ",__LINE__
    Write(*,'(2a)') blank, "CBGCMC needs some arrays to be initialized during "
    Write(*,'(2a)') blank, "the first time a move is called, &
        &That is being done now"

    ! initialize map_spc
    !** find the index of the map species, these interactions should &
    !** be evaluates first. This can be done neater if we could calll 
    !** forceield to find out map species
    MAP_SPC = simcell_getmolectype(scell)
    If (MAP_SPC>0)   Then
      MAP_EXISTS =.True.
      If (MAP_SPC>molecules_getnsorbs()) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Else
      MAP_EXISTS = .False.
      nocbgcmap=.False.
#ifdef NOCBGCPMAP
      nocbgcmap=.True.
      MAP_SPC=0
#endif
      If (.not.nocbgcmap) Then
        Write(*,'(2a)')blank, "----- no map exists for cbgcmc ------- "
        Write(*,'(2a)')blank, " This is a check: if you want to do cbgcmc without a map"
        Write(*,'(2a)')blank, " then compile with -DNOCBGCPMAP flag in Makefile "
        Write(*,'(2a)')blank, " OR, if not, make sure simcell section &
            & has the sorb name, not filename"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Endif


    Call file_gettype(d_ctrl_file, line, unitno)
    !** Find the CBGCMC section
    Write(0,'(3a)')blank, "Looking for CBGMC section in ctrlfile : ", &
        Trim(line)
    rewind(unitno)
    lineno = filesrchstr(unitno, CBGCMC_STRING, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ",  CBGCMC_STRING, " in the control file"
      Stop
    Endif

    ! number of trials in disk sampling, and sphere sampling 
    Read(unitno,*) disk_samp_num, sphere_samp_num
    NpsiSS=Int(Sqrt(sphere_samp_num*one))+1
    NthetaSS=NpsiSS
    Write(*,'(2a,i6)') blank, &
        "Actual number of sphere sampling points used : ", NpsiSS**2
    NPhiDS=disk_samp_num

    If (.Not.trial_arrays_inited) &
        Call brmoves_initTrialArrays(NPsiSS, NthetaSS, NPhiDS,DEF_indent)
    Read(unitno,*)  new_max_nrg
    If (new_max_nrg>MAX_STOP_NRG) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) " cbgcmc section should specify a max_nrg less than: ",&
          MAX_STOP_NRG
      Write(*,*) " other wise you could run into problems"
      stop
    Endif
    MAX_STOP_NRG= new_max_nrg


    ! decide whether HLFAS, LLFAS, or both
    ! if hlfas = .True. and hlfasON=.true. then HLFAS will be used
    problem=.False.
    Read(unitno,'(a)') line
    line=cleanstring(line)
    nfields=split(line,fields,",")
    Select Case (Trim(toupper(fields(1))))
    Case(hlfas_string)
      hlfas=.True.
      If (nfields==1) Then
        ! all sims hlfas
        llfas=.False.
        hlfasON=.True.
        hlfas_simno=0
      Elseif(nfields==2) Then
        ! decide when to turn on HLFAS
        hlfas_simno=toint(fields(2))
        llfas=.True.
        hlfasON=.False.
      Else
        problem=.True.
      Endif
      Read(unitno,*) NHLtrials

    Case(llfas_string)
      if (nfields>1)  Then
        problem=.True.
      else
        llfas=.true.
        hlfas=.false.
        hlfasON=.False.
      endif
    Case default
      problem=.True.
    End Select

    If (problem) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "problem with reading cbgcmc specific section (Line-3) "
      Write(*,*) "string used : ", line
      stop
    endif

    
  End Subroutine brmoves_basicInit

  !-------------------------------------------------------------------------
  ! Check the binsert and bdel moves are properly initialized 
  !-------------------------------------------------------------------------
  Subroutine brmoves_CheckInsDel(inspar, delpar, ctrlfilename)
    Type(BranchedCBIns_Params),Pointer :: inspar
    Type(BranchedCBDel_Params),Pointer :: delpar
    Character(len=strLen), Intent(in) :: ctrlfilename

    Integer :: error

    delpar%max_nrg=inspar%max_nrg

    If (Trim(inspar%biasfilename)=="NULL") Then
      delpar%biasfilename=inspar%biasfilename
    Else
      delpar%biasfilename=inspar%biasfilename
      delpar%bmap=>inspar%bmap
      delpar%rti=inspar%rti
      delpar%tk=inspar%tk
    Endif

  End Subroutine brmoves_CheckInsDel

  !------------------------------------------------------
  ! Read the information that will be required to do branched fixed length
  ! translation type of moves.
  !------------------------------------------------------
  Subroutine brmoves_FTransInit(params, sorbtype, filename)
    Type(BranchedFTrans_Params), Intent(inout) :: params
    Integer, Intent(in)                     :: sorbtype
    Character(*), Intent(in)                :: filename

    Integer    :: unitno, natoms, error, scale_trans

    !** Get the delta displacements and "whether to scale or not".
    unitno = file_getunit(filename)
    Read(unitno,*) params%deltatrans, scale_trans

    !** Set the flag to decide whether we scale the jump lengths or not
    params%scale_trans = .False.
    If (scale_trans == 1) Then
      params%scale_trans = .True.
    End If

    !** Initialize the molecular definition
    natoms = molecules_getnatoms(sorbtype)
    Allocate(params%com_defn(natoms), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Could not allocate memory for 'com_defn'"
      Stop
    End If
    Call molecules_getcomdefn(sorbtype, params%com_defn)
    !    params%com_defn = molecules_getcomdefn(sorbtype)

  End Subroutine brmoves_FTransInit

  !------------------------------------------------------
  ! Read the information that will be required to do branched random length
  ! translation type of moves.
  !           auxmv -- contains info about cavity list
  !------------------------------------------------------
  Subroutine brmoves_RTransInit(params, simcell, auxmv, sorbtype, filename)
    Type(BranchedRTrans_Params), Intent(inout) :: params
    Type(SimCell_Params), Target, Intent(in) :: simcell
    Integer, Intent(in)                     :: sorbtype
    Character(*), Intent(in)                :: filename
    Type(AuxMoveObjects),Pointer           :: auxmv 
    Integer    :: unitno, natoms, error, scale_trans

    !** point the cavity bias object if needed, 
    !** cavity is a defined as a unversal varaible in this Module
    If (auxmv%cavBiasON) Then
      If (.Not.cavity_initialized) Then
        Nullify(cavity)
        cavity=>auxmv%cavity
        cavity_initialized=.True.
      Endif
    Endif

    !** Set the constrained flag to .False.
    params%constrained = .False.

    !** Set the simcell pointer
    params%scell => simcell

    params%max_nrg=VERYLRGMAXNRG

    !** Get the delta displacements and "whether to scale or not".
    unitno = file_getunit(filename)
    Read(unitno,*) params%deltatrans, scale_trans

    !** Set the flag to decide whether we scale the jump lengths or not
    params%scale_trans = .False.
    If (scale_trans == 1) Then
      params%scale_trans = .True.
    End If

    !** Initialize the molecular definition
    natoms = molecules_getnatoms(sorbtype)
    Allocate(params%com_defn(natoms), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Could not allocate memory for 'com_defn'"
      Stop
    End If
    Call molecules_getcomdefn(sorbtype, params%com_defn)
    !    params%com_defn = molecules_getcomdefn(sorbtype)

  End Subroutine brmoves_RTransInit

  !------------------------------------------------------
  ! Read the information that will be required to do fixed angle 
  ! rotation type of moves.
  !  auxmv -- contains info about cavity list
  !------------------------------------------------------
  Subroutine brmoves_RRotInit(params, sorbtype, auxmv, filename)
    Type(BranchedRRot_Params), Intent(inout) :: params
    Integer, Intent(in)        :: sorbtype
    Character(*), Intent(in)   :: filename
    Type(AuxMoveObjects),Pointer           :: auxmv 

    Integer    :: unitno, natoms, error, scale_rot

    !** point the cavity bias object if needed
    !** cavity is a defined as a unversal varaible in this Module
    If (auxmv%cavBiasON) Then
      If (.Not.cavity_initialized) Then
        Nullify(cavity)
        cavity=>auxmv%cavity
        cavity_initialized=.True.
      Endif
    Endif

    !** Get the delta displacements
    unitno = file_getunit(filename)
    Read(unitno,*) params%deltarot, scale_rot

    !** Set the flag to decide whether we scale the jump lengths or not
    params%scale_rot = .False.
    If (scale_rot == 1) Then
      params%scale_rot = .True.
    End If

    !** Initialize the molecular definition
    natoms = molecules_getnatoms(sorbtype)
    Allocate(params%com_defn(natoms), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Could not allocate memory for 'com_defn'"
      Stop
    End If
    Call molecules_getcomdefn(sorbtype, params%com_defn)
    !    params%com_defn = molecules_getcomdefn(sorbtype)

  End Subroutine brmoves_RRotInit

  !------------------------------------------------------
  ! Initialize all info related to brmoves, that need imodel
  ! e.g., angle and frame libraries, torsion params 
  ! tk - temperature, K.
  !------------------------------------------------------
  Subroutine brmoves_initAuxLibs(auxlib, sorbate, sorbtype, imodel, tk, indent)
    Type(AuxLibraries),Intent(inout) :: auxlib
    Type(AtMolCoords), Intent(in) :: sorbate
    Integer, Intent(in) :: sorbtype,indent
    Type(Interaction_Model), Pointer :: imodel
    Real(kind=RDbl),  Intent(in) :: tk

    Integer :: natoms, max_params, error, nsets, nlibs
    Integer :: a, l, c, p, nfields, other_atom
    Integer, Dimension(3)                 :: subset1, subset2
    Integer, Dimension(:,:), Pointer                :: list
    Integer, Dimension(:), Pointer                :: atom_num
    Character(len=lstrLen), Dimension(:), Pointer   :: params
    Character(len=lstrLen), Dimension(10) :: fields   
    Character(len=lstrLen) :: modelname     
    Character(len=indent) :: blank

    Real(kind=RDbl) :: k_theta, theta_eq, sigma, eps, A_lj, B_lj 
    Real(kind=RDbl), Dimension(2) :: param_array
    Type(NodeType), Pointer   :: node 
    Logical :: found
    Logical,Dimension(:),Pointer :: is_higheratom
    blank=Repeat(" ",indent)
    natoms = molecules_getnatoms(sorbtype)
    If (natoms<3) Then
      Write(*,*) "angles library not initiaized for species:", sorbtype
      auxlibs_initialized(sorbtype)=.True.
      Return
    Endif

    ! assumes maximum possible number of angles/torsions/intrpair
    max_params = natoms*natoms
    Allocate(list(max_params, 5), params(max_params), atom_num(max_params), &
        STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Allocate(auxlib%angle_lib_list(natoms), auxlib%tor_list(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    ! assumes that coords are initialized for at least one molecule
    If(Associated(sorbate%gcoords(1)%branchedcoords%root%C)) Then 

      ! assumes linear alkne, assumes that there is only one angle 
      ! potential type throughout the chain
      subset1=(/ sorbtype, 0, 0 /)  
      subset2=(/ 0, 0, 0 /)

      ! subset1, subset2 needed for calling interact.
      Call interact_listparams(imodel, subset1, subset2, BENDING_KEY, &
          nsets, list, params )

      ! find number of libraries, atom_num returns the atoms on which these
      ! angle/libraries exist, works only for linear alkane.
      ! vague routine, works only if initially nlibs=0
      nlibs =0
      Call angle_CountLib(sorbate%gcoords(1)%branchedcoords%root%C, nlibs, &
          atom_num)

      If (nsets/=nlibs) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      auxlib%nangles=nlibs
    Endif
    If (nlibs > 0) Then
      Allocate(auxlib%angle_lib(nlibs),STAT=error) 
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      ! initialize each of the libraries, assumes linear alkane.
      Do l=1, nlibs
        Call branchedcoords_findnode(sorbate%gcoords(1)%branchedcoords%root &
            , atom_num(l), node, found)
        If (.Not.found) Then
          Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
              " Could not find the node corresponding to this atom"
          Stop
        End If
        nfields=split(params(l),fields)
        If (nfields<3) then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) " not enough parameters for bond angle"
          stop
        else
          modelname=cleanstring(fields(1))
          k_theta=toreal(fields(2))
          ! convert from degrees to radians 
          theta_eq=toreal(fields(3))*degToRad 
          param_array=(/ k_theta, theta_eq /)
        endif

        Call angle_init(node, auxlib%angle_lib(l), tk, param_array, modelname) 
        auxlib%angle_lib_list(atom_num(l))=l ! for easy reverse mapping

      End Do
    End If

    !** this is the stuff for inserting a library of 2 bracnhing atoms
    ! not tested as of now, 11/20/2002 shaji. 
!!$      lib_num = 0 
!!$      Call frame_CountLib(sorbates(sorbtype)%gcoords(1)&
!!$          %branchedcoords%root%C,lib_num,atom_num)
!!$      params%framelib_num = lib_num  
!!$
!!$      If (lib_num > 0) Then
!!$        Allocate(Params%frame_lib(lib_num),STAT=error) 
!!$        If (error /= 0) Then
!!$          Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
!!$              " Could not allocate 'frame_lib' "
!!$          Stop
!!$        End If
!!$
!!$        Do j=1, lib_num
!!$          Call find_node(sorbates(sorbtype)%gcoords(1)%branchedcoords &
!!$              %root,atom_num(j),node,found)
!!$          If (found /= 1) Then
!!$            Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
!!$                " Could not find the branched point "
!!$            Stop
!!$          End If
!!$          Call frame_init(node, params%frame_lib(j), &
!!$              1,tk) 
!!$        End Do
!!$      End If
!!$    End If
!!$


    If (natoms<4) Then
      Write(*,'(2a,i3)') blank, "Finished initalizing libraries for cbgcmc &
          &for sorb: ",sorbtype
      auxlibs_initialized(sorbtype)=.True.
      Return
    Endif


    !** initialise the torsion parameters that will be stored in the 
    !** binsert-params for faster reference

    ! subset1, subset2 needed for calling interact.
    Call interact_listparams(imodel, subset1, subset2, TORSION_KEY, &
        nsets, list, params )
    auxlib%ntors=nsets

    Allocate(auxlib%tor(nsets,0:COS_EXP_MAX_PARAMS-1), &
        auxlib%tor_list(natoms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"torsion-params")
    Do l=1,nsets
      ! given atom number of the 4th atom in a torsion list we 
      ! should be able to tell which is the tor set right for us
      ! i.e., given atom num= list(l,4) we return tor(l)
      auxlib%tor_list(list(l,4))=l

      ! extra check, may be not required, make sure -shaji
      If (list(l,4)<list(l,1)) Then
        Write(*,*) "specify atom torsion in molecule file in proper fashion"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      nfields=split(params(l),fields)
      If (nfields/=7) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) " not enough parameters for bond torsion"
        Write(*,*) " modelname and 5 coefficients expected"
        Write(*,*) " but got only : "//Trim(params(l))
        stop
      else

        modelname=toupper(cleanstring(fields(1)))
        If (trim(modelname)/="COSEXPANSION") then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "wrong modelname", Trim(modelname)
          stop
        endif
        ! should be in units of kcal
        Do c=0,5
          auxlib%tor(l,c)=toreal(fields(c+2))
        end do
      endif
    End Do



    If (natoms<5) Then
      auxlibs_initialized(sorbtype)=.True.
      Write(*,'(2a,i3)') blank, "Finished initalizing libraries for cbgcmc &
          &for sorb: ",sorbtype
      Return
    Endif

    !** initialise the ipmodel parameters that will be stored in the 
    !** binsert-params for faster reference

    ! subset1, subset2 needed for calling interact.
    Call interact_listparams(imodel, subset1, subset2, INTRAPAIR_KEY, &
        nsets, list, params )

    ! here we will do it bit differently compared to torsion/angle
    ! we will allocate arrays for each atom
    ! last index is two because we have to store A and B
    Allocate(auxlib%ip_ABvals(nsets,2), auxlib%n_ip_pots(natoms),&
        auxlib%ip_pots_atomlist(natoms,natoms), &
        auxlib%ip_pots_index(natoms,natoms),STAT=error)
    auxlib%ip_ABvals(:,:) = zero
    auxlib%n_ip_pots(:)= 0
    auxlib%ip_pots_atomlist(:,:) = 0
    auxlib%ip_pots_index(:,:) = 0

    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"ip_list-params")

    Do l=1,nsets
      nfields=split(params(l),fields)
      If (nfields/=3) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) " not enough parameters for intra pair"
        Write(*,*) " modelname and 2 coefficients expected"
        Write(*,*) " but got only : "//Trim(params(l))
        stop
      Else

        modelname=toupper(cleanstring(fields(1)))
        If (Trim(modelname)/="LJ") Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "wrong modelname for ip in cbgcmc ", Trim(modelname)
          stop
        Endif
        ! store atom-numbers in proper fashion 
        If (list(l,2)<list(l,1)) Call swap( list(l,2), list(l,1) )
        ! should be in units of Ang. and K (sigma and epsilon)
        sigma=toreal(fields(2))
        eps=toreal(fields(3))
        B_lj = 4.0000_Rdbl * eps * (sigma**6) * kcalmole_kb
        A_lj = B_lj * (sigma**6)
        auxlib%ip_ABvals(l,1)=A_lj
        auxlib%ip_ABvals(l,2)=B_lj
      Endif
#ifdef DEBUG
      Write(*,*) "List : ",list(l,1),list(l,2)
      Write(*,*) " A, B", auxlib%ip_ABvals(l,1),  auxlib%ip_ABvals(l,2)
#endif
    End Do


    ! now go through atom lsit and initialize rest
    Allocate(is_higheratom(natoms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Do a=1,natoms
      found=.False.
      is_higheratom=.False.
      Call  branchedcoords_markhigheratoms(&
          sorbate%gcoords(1)%branchedcoords%root, is_higheratom,a,found)
      If (.Not.found) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__ 
        Stop
      Endif
      Do l=1,nsets
        If (list(l,1)==a) Then
          other_atom=list(l,2)
        Elseif(list(l,2)==a) Then
          other_atom=list(l,1)
        Else
          Cycle
        Endif
        If (is_higheratom(other_atom)) Then

          ! increase num of pots at that atom
          auxlib%n_ip_pots(a)= auxlib%n_ip_pots(a)+1

          ! mark that atom
          auxlib%ip_pots_atomlist(a, auxlib%n_ip_pots(a) ) = other_atom

          ! mark the potential
          auxlib%ip_pots_index(a, auxlib%n_ip_pots(a)) = l

        Endif

      End Do

    End Do

    Deallocate(list, params, atom_num, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Write(*,'(2a,i3)') blank, "Finished initalizing libraries for cbgcmc &
        &for sorb: ",sorbtype
    auxlibs_initialized(sorbtype)=.True.

  End Subroutine brmoves_initAuxLibs

  !------------------------------------------------------------------------
  ! Auxlibrary libraries are sized here. They are initialized later
  ! when the corresponding moves are called
  !------------------------------------------------------------------------
  Subroutine brmoves_sizeAuxLibs()
    Integer :: nspc, spc, error
    auxlibs_sized=.True.
    nspc=molecules_getnsorbs()
    Allocate(auxlibs_array(nspc),  auxlibs_initialized(nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do spc=1,nspc
      auxlibs_initialized(spc)=.False.
      auxlibs_array(spc)%initialized=.False.
      auxlibs_array(spc)%sorbtype=spc
      auxlibs_array(spc)%natoms=molecules_getnatoms(spc)
    End Do


  End Subroutine brmoves_sizeAuxLibs

  !------------------------------------------------------------------------
  ! Reads the information required for doing Insertion type of move and
  ! Initilises the relecvant variables
  !  auxmv -- contains info about cavity list
  !------------------------------------------------------------------------
  Subroutine brmoves_CBInsInit(params,simcell, auxmv, sorbtype, &
      sorbates, filename)
    Type(BranchedCBIns_Params), Intent(inout) :: params
    Type(SimCell_Params), Target, Intent(in)  :: simcell
    Integer, Intent(in)         :: sorbtype
    Character(*), Intent(in)    :: filename
    Type(AtMolCoords), Dimension(:), Intent(in):: sorbates 
    Type(AuxMoveObjects),Pointer           :: auxmv 

    Character(len=strLen)   :: biasfilename
    Real(kind=RDbl), Dimension(3) :: eff
    Real(kind=RDbl)               :: tk
    Integer    :: j, unitno, natoms, error, lib_num
    Integer, Dimension(sorbates(sorbtype)%natoms) :: atom_num
    Logical :: found

    Type(NodeType), Pointer   :: node 

    If (.Not.auxlibs_sized) Call brmoves_sizeAuxLibs()


    unitno = file_getunit(filename)
    !** point the cavity bias object if needed
    !** cavity is a defined as a unversal varaible in this Module
    If (auxmv%cavBiasON) Then
      If (.Not.cavity_initialized) Then
        Nullify(cavity)
        cavity=>auxmv%cavity
        cavity_initialized=.True.
      Endif
    Endif

    !** Set the simcell pointer
    params%scell => simcell

    !** Read the bias filename
    Read(unitno,*) biasfilename
    biasfilename=Trim(stripcmnt(biasfilename))
    params%biasfilename = biasfilename
    Read(unitno,*) tk
    params%tk=tk
    !** Set the temperature parameter
    params%rti     = (1.0_RDbl)/(tk*kcalmole_kb)
    params%max_nrg=MAX_STOP_NRG  ! dubious, should be fixed later, shaji
    !** Initialize the bias map
    eff = simcell_geteff(simcell, .True.)
    Call bmap_init(params%bmap, eff, tk, biasfilename, simcell)

    !** Initialize the molecular definition
    natoms = molecules_getnatoms(sorbtype)

  End Subroutine brmoves_CBInsInit


  !------------------------------------------------------------------------
  ! Reads the information required for doing identity change 
  ! Initilises the relecvant variables
  ! Works only for linear alkanes
  !  auxmv -- contains info about cavity list
  !------------------------------------------------------------------------
  Subroutine brmoves_idflipInit(params,scell, auxmv, spc, species, ctrlfile)
    Type(BranchedIdFlip_Params), Intent(inout) :: params
    Type(SimCell_Params), Target, Intent(in)  :: scell
    Integer, Intent(in)         :: spc
    Character(*), Intent(in)    :: ctrlfile
    Type(AtMolCoords), Dimension(:), Intent(in):: species 
    Type(AuxMoveObjects),Pointer           :: auxmv 

    Character(len=strLen)   :: biasfilename
    Real(kind=RDbl), Dimension(3) :: eff
    Real(kind=RDbl)               :: tk
    Integer    :: unitno

    If (.Not.auxlibs_sized) Call brmoves_sizeAuxLibs()
    unitno = file_getunit(ctrlfile)

    !** point the cavity bias object if needed
    !** cavity is a defined as a unversal varaible in this Module
    If (auxmv%cavBiasON) Then
      If (.Not.cavity_initialized) Then
        Nullify(cavity)
        cavity=>auxmv%cavity
        cavity_initialized=.True.
      Endif
    Endif

    !** Set the simcell pointer
    params%scell => scell

    !** Read the bias filename
    Read(unitno,*) biasfilename
    biasfilename=Trim(stripcmnt(biasfilename))
    params%biasfilename = biasfilename
    Read(unitno,*) tk
    params%tk=tk

    !** Set the temperature parameter
    params%rti     = (1.0_RDbl)/(tk*kcalmole_kb)
    params%max_nrg=MAX_STOP_NRG  ! dubious, should be fixed later, shaji

    !** Initialize the bias map
    eff = simcell_geteff(scell, .True.)
    Call bmap_init(params%bmap, eff, tk, biasfilename,scell)

  End Subroutine brmoves_idflipInit


  !------------------------------------------------------------------------
  ! Reads the information required for doing Insertion type of move and
  ! Initilises the relecvant variables
  ! auxmv -- contains info about cavity list
  !------------------------------------------------------------------------
  Subroutine brmoves_reGrowInit(params,simcell, auxmv, sorbtype, &
      sorbates, filename)
    Type(BranchedRegrow_Params), Intent(inout) :: params
    Type(SimCell_Params), Target, Intent(in)  :: simcell
    Integer, Intent(in)         :: sorbtype
    Character(*), Intent(in)    :: filename
    Type(AtMolCoords), Dimension(:), Intent(in):: sorbates 
    Type(AuxMoveObjects),Pointer           :: auxmv 

    Character(len=strLen)   :: biasfilename
    Real(kind=RDbl), Dimension(3) :: eff
    Integer    :: j, unitno, natoms, error, lib_num, found
    Integer, Dimension(sorbates(sorbtype)%natoms) :: atom_num 

    Type(NodeType), Pointer   :: node 

    If (.Not.auxlibs_sized) Call brmoves_sizeAuxLibs()

    !** point the cavity bias object if needed
    !** cavity is a defined as a unversal varaible in this Module
    If (auxmv%cavBiasON) Then
      If (.Not.cavity_initialized) Then
        Nullify(cavity)
        cavity=>auxmv%cavity
        cavity_initialized=.True.
      Endif
    Endif

    !** Set the simcell pointer
    params%scell => simcell

    !** Read the bias filename
    unitno = file_getunit(filename)
    Read(unitno,*) biasfilename
    biasfilename=Trim(stripcmnt(biasfilename))
    params%biasfilename = biasfilename
    Read(unitno,*) params%tk

    !** Set the temperature parameter
    params%rti     = (1.0_RDbl)/(params%tk*kcalmole_kb)

    !** Initialize the bias map
    eff = simcell_geteff(simcell, .True.)
    Call bmap_init(params%bmap, eff, params%tk, biasfilename, simcell)

    !** Initialize the molecular definition
    natoms = molecules_getnatoms(sorbtype)

    ! Initialize the variable used to store the old coordinates and
    ! suggested coordinates during the cut regrow move
    Call config_allocfields(params%oldcoords, sorbtype, 1)
    Call config_allocfields(params%suggested, sorbtype, 1)

    params%max_nrg=MAX_STOP_NRG  ! dubious, should be fixed later, shaji

  End Subroutine brmoves_reGrowInit

  !------------------------------------------------------------------------
  ! Initializes the information required for doing Insertion type of move 
  ! with no bias
  !------------------------------------------------------------------------
  Subroutine brmoves_RInsInit(params, simcell, sorbtype)
    Type(BranchedRIns_Params), Intent(inout) :: params
    Type(SimCell_Params), Target, Intent(in)  :: simcell
    Integer, Intent(in)         :: sorbtype

    Integer    :: natoms

    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    stop

  End Subroutine brmoves_RInsInit

  !------------------------------------------------------------------------
  ! Reads the information required for doing Deletion type of move and
  ! Initilises the relevant variables- for biased case
  !  auxmv -- contains info about cavity list
  !------------------------------------------------------------------------
  Subroutine brmoves_CBDelInit(params, &
      simcell, auxmv, sorbtype, sorbates, filename)
    Type(BranchedCBDel_Params), Intent(inout) :: params
    Type(SimCell_Params), Target, Intent(in)  :: simcell
    Integer, Intent(in)         :: sorbtype
    Real(kind=RDbl)             :: tk
    Character(*), Intent(in)    :: filename
    Type(AtMolCoords), Dimension(:), Intent(in):: sorbates  
    Type(AuxMoveObjects),Pointer           :: auxmv 

    Character(len=strLen)   :: biasfilename
    Real(kind=RDbl), Dimension(3) :: eff
    Integer    :: j, unitno, natoms, error, lib_num, found
    Integer, Dimension(sorbates(sorbtype)%natoms) :: atom_num  

    Type(NodeType), Pointer   :: node 

    !** point the cavity bias object if needed
    !** cavity is a defined as a unversal varaible in this Module
    If (auxmv%cavBiasON) Then
      If (.Not.cavity_initialized) Then
        Nullify(cavity)
        cavity=>auxmv%cavity
        cavity_initialized=.True.
      Endif
    Endif

    !** Set the simcell pointer
    params%scell => simcell

    If (.Not.auxlibs_sized) Call brmoves_sizeAuxLibs()

  End Subroutine brmoves_CBDelInit

  !------------------------------------------------------------------------
  ! Reads the information required for doing Deletion type of move and
  ! Initilises the relevant variables - for unbiased case.
  !------------------------------------------------------------------------
  Subroutine brmoves_RDelInit(params, simcell)
    Type(BranchedRDel_Params), Intent(inout) :: params
    Type(SimCell_Params), Target, Intent(in)  :: simcell
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  End Subroutine brmoves_RDelInit


  !-----------------------------------------------------------
  ! Converts to xyz coords and do periodic bc of the coordinates
  ! natoms=total no of atoms in molecule
  ! atomno = specific atom for which pbs is required (Optional)
  ! NOTE : at present "atomno" is not used from many possible places ,fix this
  !-----------------------------------------------------------
  Subroutine brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms,atomno)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(IN)             :: scell
    Integer , Intent(in) :: molec,sorbtype,natoms
    Integer , Intent(in),Optional :: atomno

    Call branchedcoords_toxyz  & 
        (sorbates(sorbtype)%gcoords(molec)%branchedcoords, & 
        sorbates(sorbtype)%coords(1:natoms,molec)%rp)

    If (Present(atomno)) Then
      Call simcell_pbc(scell, sorbates(sorbtype)%coords(atomno,molec)%rp, &
          sorbates(sorbtype)%coords(atomno, molec)%r, &
          sorbates(sorbtype)%coords(atomno, molec)%cr)
    Else
      Call simcell_pbc(scell, sorbates(sorbtype)%coords(1:natoms,molec)%rp, &
          sorbates(sorbtype)%coords(1:natoms, molec)%r, &
          sorbates(sorbtype)%coords(1:natoms, molec)%cr)
    Endif

  End Subroutine brmoves_xyz_pbc

  !-----------------------------------------------------------
  ! Calculates the nrg of the atom with all other sorbs
  ! nrg_array( i,j)  -> i=1, 2          1=NONCOUL 2=COUL
  !                  -> j=1, nsorbs
  ! As of now COUL is not calculated and zero is returned
  !-----------------------------------------------------------
  Subroutine brmoves_get_atom_sorbNrg(species,scell,imodel,spc,&
      mol,atm,fast,noncoulnrg,nrg_calc_flag,nrg_arr,max_nrg)
    Type(Interaction_Model), Pointer         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: species
    Type(SimCell_Params), Intent(IN)             :: scell
    Integer, Intent(in) :: spc,mol,atm
    Logical ,Intent(in) :: fast
    Real(kind=RDbl),Intent(out)     :: noncoulnrg
    Logical ,Intent(out) :: nrg_calc_flag
    Real(kind=RDbl),Dimension(:,:),Intent(inout),Optional     :: nrg_arr
    Real(kind=RDbl),Intent(in),Optional     :: max_nrg

    Integer, Dimension(3):: subset1, subset2
    Logical :: successflag, incavity
    Logical, Save  :: first_time = .True.
    Integer :: spc2,nspc, loop_index, natoms
    Real(kind=RDbl) :: pot,tempmaxnrg

    If (first_time) Then
      first_time = .False.
      Call storebase_init(ffout,0)
    Endif

    If (.Not.fast) Then
      Write(*,*) "You might be making a mistake, fast=.true. by default"
      Write(*,*) "Make sure everything else is OK, recompile and run.. "
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    nrg_calc_flag=.True.
    noncoulnrg=zero
    nspc=molecules_getnsorbs()  

    tempmaxnrg=MAX_STOP_NRG
    If (Present(max_nrg)) tempmaxnrg=max_nrg

    !** Check whether this is in a cavity
    If (cavity_initialized) Then
      If (avoid_cubes) Then
        natoms=molecules_getnatoms(spc)
        incavity= cavitylist_checkcavity(cavity, &
            species(spc)%coords(atm,mol)%r, cubes2avoid(1:natoms))
      Else
        incavity= cavitylist_checkcavity(cavity, &
            species(spc)%coords(atm,mol)%r)
      Endif
      If (.Not. (incavity) )  Then 
        nrg_calc_flag=.False.     
        Return
      End If
    Endif


    Do loop_index=1, nspc
      spc2=loop_index

      ! first calculate map, then ssbasic, so interchange first spc and map_spc
      If (MAP_EXISTS) Then
        If (loop_index==1) Then
          spc2=MAP_SPC
        Elseif (loop_index==MAP_SPC) Then
          spc2=1
        Endif
      Endif

      pot = 0.0_RDbl
      subset1=(/ spc, mol, atm /)
      subset2=(/ spc2, 0, 0 /)

      successflag = interact_simpleint(imodel,ffout,subset1,subset2,species, &
          scell,fast)
      If (.Not. (successflag) )  Then 
        nrg_calc_flag=.False.     
        Return
      End If

      pot= ffout%nrg
      noncoulnrg = noncoulnrg + pot
      If (Present(nrg_arr)) Then
        nrg_arr(2,spc2) = zero   ! now coul is zero by default
        nrg_arr(1,spc2) = pot
      Endif

      If (noncoulnrg>(tempmaxnrg)) Then

        If (.Not.Present(max_nrg)) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
        nrg_calc_flag=.False.     
        Return
      Endif
      tempmaxnrg=tempmaxnrg-pot

    End Do
  End Subroutine brmoves_get_atom_sorbNrg


  !-----------------------------------------------------------
  ! Grows a linear alkane in the simulation cell uses a bias map 
  ! The angles are generated using various biasing techniques 
  ! It returns the biasfactor also
  ! **biasfactor** is the factor to be used in the acceptance criteria 
  ! to remove the bias of this particular insertion. 
  ! It is essentially the inverse of our bias for this configuration
  ! nrg_calc_flag=.FALSE. if nrg is too high or could not be calculated
  !           subint -- subset interactions (molecule-system)
  !-----------------------------------------------------------
  Subroutine brmoves_cbinsert(sorbates,sorbtype,molec,params, & 
      biasfactor, nrg_calc_flag,subint)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)                              :: sorbtype,molec
    Type(BranchedCBIns_Params), Intent(in)        :: params
    Real(kind=RDbl), Intent(out)                     :: biasfactor
    Logical, Intent(inout)                           :: nrg_calc_flag  
    Type(Subset_Interactions), Intent(InOut)         :: subint    

    Real(kind=RDbl)                      :: ufinal, intra 
    Real(kind=RDbl) , Dimension(2,Size(sorbates))     :: nrg_arr

    Integer, Dimension(MAX_ATOMS)                     :: atom_num 
    Real(kind=RDbl)      :: bias, rti, u, temp_max_nrg
    Real(kind=RDbl)      :: first_atom_nrg, second_atom_nrg, third_atom_nrg
    Real(kind=RDbl)      :: fourth_atom_nrg, current_total_nrg,tmpIntra
    Real(kind=RDbl) ,Dimension(2,Size(sorbates,1)) :: tmpNrgArr 
    Logical              :: fast, no_intra
    Integer              :: natoms,nsorbs,k
    Type(NodeType), Pointer       :: root_ptr
    natoms=molecules_getnatoms(sorbtype)

    If  (firsttime_in_brmoves) &
        Call brmoves_basicInit(params%scell,DEF_indent)
    If (.Not.auxlibs_initialized(sorbtype)) Then
      Call brmoves_initAuxLibs(&
          auxlibs_array(sorbtype), sorbates(sorbtype), sorbtype, &
          subint%imodel, params%tk,DEF_indent)
    Endif

    !** For easy refernce to the co-ordinates of molecule
    !** make sure everything is correct with the f90-pointer business here
    root_ptr=> sorbates(sorbtype)%gcoords(molec)%branchedcoords%root

    nsorbs=molecules_getnsorbs()
    nrg_calc_flag =  .False.
    fast = .True.
    rti = params%rti 
    biasfactor=one
    ufinal=zero
    nrg_arr(1:2,1:nsorbs) = zero
    temp_max_nrg=params%max_nrg

    If (cavity_initialized) avoid_cubes=.False.

    !**
    !**       INSERT FIRST ATOM
    !** 
    atom_num(1) = root_ptr%atom_num  
    Call brmoves_insert_atom1(params%bmap, params%scell, subint%imodel, &
        sorbates, sorbtype, molec, atom_num(1), rti, bias, u, nrg_calc_flag, &
        temp_max_nrg,tmpNrgArr )
    If (.Not.nrg_calc_flag) Return

    first_atom_nrg=u
    current_total_nrg=u
    Do k=1,nsorbs
      nrg_arr(1,k)=tmpNrgArr(1,k) !only noncoul is used here
    End Do

    biasfactor=biasfactor*bias

    If(natoms>1) Then
      !**
      !**       INSERT SECOND ATOM
      !** 
      temp_max_nrg=params%max_nrg-first_atom_nrg
      atom_num(2)=root_ptr%C%atom_num  
      Call  brmoves_insert_atom2( params%scell, sorbates, subint%imodel, sorbtype,&
          molec, atom_num(2), rti, bias, second_atom_nrg ,nrg_calc_flag, &
          temp_max_nrg, tmpNrgArr)
      If (.Not.nrg_calc_flag) Return

      biasfactor=biasfactor*bias
      ufinal=first_atom_nrg+second_atom_nrg
      Do k=1,nsorbs
        nrg_arr(1,k)=nrg_arr(1,k) + tmpNrgArr(1,k) !only noncoul is used here
      End Do
    Endif

    If(natoms>2) Then
      !**
      !**       INSERT ATOM-3 
      !** 
      atom_num(3) = root_ptr%C%C%atom_num  
      temp_max_nrg=params%max_nrg-ufinal
      Call  brmoves_insert_atom3( params%scell, sorbates, subint%imodel, sorbtype,&
          molec, atom_num(3), rti, bias, third_atom_nrg ,nrg_calc_flag, &
          temp_max_nrg, intra, tmpNrgArr)

      If (.Not.nrg_calc_flag) Return 

      biasfactor=biasfactor*bias
      ufinal = ufinal + third_atom_nrg
      Do k=1,nsorbs
        nrg_arr(1,k)=nrg_arr(1,k) + tmpNrgArr(1,k) !only noncoul is used here
      End Do
      If (dbgflag) Write(*,*) atom_num(3), ufinal, intra, biasfactor
    Endif

    !**
    !**       INSERT ATOMS 4 Onwards 
    !** 
    If (natoms>3) Then
      If(Associated(root_ptr%C%C%C)) Then
        temp_max_nrg=params%max_nrg-ufinal
        Call brmoves_insertnode(sorbates,  params%scell, subint%imodel, sorbtype, &
            molec, root_ptr%C%C%C, bias, fourth_atom_nrg, nrg_calc_flag, &
            temp_max_nrg, params%rti, tmpIntra, tmpNrgArr)

        If (nrg_calc_flag) Then
          ufinal=ufinal+fourth_atom_nrg
          biasfactor=biasfactor*bias
          intra=intra+tmpIntra
          Do k=1,nsorbs
            nrg_arr(1,k)=nrg_arr(1,k) + tmpNrgArr(1,k) !only noncoul is used here
          End Do
        Endif

      End If
    Endif
    If (ufinal>params%max_nrg) Then
!!$      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$      Write(*,*) "legitimate insertion exceeds max_nrg, &
!!$          & consider increasing the max_nrg"
      nrg_calc_flag=.False.
      Return
    Endif

    If (nrg_calc_flag) Then  
      ! recalculate the energy to stick to general music policy of 
      ! returning the energies in subint, recalc_flag=.True.
      no_intra=.False.
      Call brmoves_xyz_pbc(sorbates, params%scell, molec, &
          sorbtype, natoms)        
      
      ! see notes in the header to see the reason for increasing max_nrg
      nrg_calc_flag = subinteract_int(subint, sorbates, params%scell, &
          fast, .True., no_intra, (/50*Abs(params%max_nrg)/), &
          (/sorbtype,molec,0/))

      If (.Not.nrg_calc_flag) Then
        ! we have serious problems ??, 
        Write(*,*) "ufinal, fourth_atom_nrg, intra",ufinal, &
            fourth_atom_nrg, intra
        Write(*,*) "biasfactor, tmpintra,params%max_nrg",biasfactor, &
            tmpintra,params%max_nrg
        Write(*,*) sorbtype, molec
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    endif

    Return

  End Subroutine brmoves_cbinsert




  !-----------------------------------------------------------
  ! removes spc1,mol1 and adds spc2,mol2 at he same place
  ! If spc1 longer than spc2 then the extra atoms are deleted
  ! If spc2 longer than spc1 then the extra atoms are added 
  ! Works only for linear alkanes
  ! nrg_calc_flag=.FALSE. if nrg is too high or could not be calculated
  !           subint -- subset interactions (molecule-system)
  ! subint -- nrg of spc2,mol2 with system
  ! mol2 -- molecule to be added at the end of the list
  !-----------------------------------------------------------
  Subroutine brmoves_idflip(params, species,spc1, mol1, spc2, mol2, & 
      biasfactor, nrg_calc_flag,subint)
    Type(BranchedIDFlip_Params), Intent(in)        :: params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: species
    Integer, Intent(in)                              :: spc1,spc2, mol1,mol2 

    Real(kind=RDbl), Intent(out)                     :: biasfactor
    Logical, Intent(inout)                           :: nrg_calc_flag  
    Type(Subset_Interactions), Intent(InOut)         :: subint    

    Logical              :: fast, no_intra, succ_flag
    Integer              :: na1, na2, nsorbs,k, natoms
    Type(NodeType), Pointer       :: root1_ptr,  root2_ptr
    Type(NodeType), Pointer       :: last1_ptr, last2_ptr

    If  (firsttime_in_brmoves) &
        Call brmoves_basicInit(params%scell,DEF_indent)

    If (.Not.auxlibs_initialized(spc2)) Then
      Call brmoves_initAuxLibs(&
          auxlibs_array(spc2), species(spc2), spc2, &
          subint%imodel, params%tk,DEF_indent)
    Endif
    If (.Not.auxlibs_initialized(spc1)) Then
      Call brmoves_initAuxLibs(&
          auxlibs_array(spc1), species(spc1), spc1, &
          subint%imodel, params%tk,DEF_indent)
    Endif

    fast=.True.
    na1=molecules_getnatoms(spc1)
    na2=molecules_getnatoms(spc2)

    ! some cavity bias things may be used during delete, the atoms 
    ! occupied be the existing molecule should not be included in the 
    ! cavity definition
    If (cavity_initialized) Then
      Call cavitylist_getCubeIndices(cavity, &
          species(spc1), mol1, na1, cubes2avoid)
      avoid_cubes=.True.
    Endif

    ! config check is done against current size
    Call config_checkandincr(species,spc2,mol2-1)

    !** copy from 1 to 2
    species(spc2)%gcoords(mol2)%branchedcoords%startpos= &
        species(spc1)%gcoords(mol1)%branchedcoords%startpos
    root1_ptr=> species(spc1)%gcoords(mol1)%branchedcoords%root
    root2_ptr=> species(spc2)%gcoords(mol2)%branchedcoords%root

    ! copy coordsvalues, placed data etc.
    Call branchedcoords_copyvalues(root1_ptr, root2_ptr, last1_ptr,last2_ptr)
    Call brmoves_xyz_pbc(species, params%scell, mol2, spc2, na2)        

    If (na2>na1) Then

      !** Remove the specified molecule and its space
      Call config_delmol(species,spc1,mol1)

      !** Create space for a new molecules of spc2 
      Call config_setnmoles(species,spc2,mol2)
      succ_flag = subinteract_changenmoles(subint,spc2,mol2)
      Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

      ! adds more atoms in spc2,mol2
      Call brmoves_idflipADD(species, params%scell, subint%imodel, &
          spc2, mol2, last2_ptr, biasfactor,nrg_calc_flag, params%max_nrg, &
          params%rti )

    Elseif (na2<na1) Then

      ! find the bias for deleting atoms from spc1, mol1
      Call brmoves_idflipDEL(species, params%scell, subint%imodel, &
          spc1, mol1, last1_ptr, biasfactor,nrg_calc_flag, params%max_nrg, &
          params%rti )

      ! do this to be consistent with idchange routine
      !** Remove the specified molecule and its space
      Call config_delmol(species,spc1,mol1)

      !** Create space for a new molecules of spc2 
      Call config_setnmoles(species,spc2,mol2)
      succ_flag = subinteract_changenmoles(subint,spc2,mol2)
      Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

      If (.not.nrg_calc_flag) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    If (.Not.nrg_calc_flag) Return 
    no_intra=.False.
    Call brmoves_xyz_pbc(species, params%scell, mol2, &
        spc2, na2)        

    ! see notes in the header to see the reason for increasing max_nrg
    nrg_calc_flag = subinteract_int(subint, species, params%scell, &
        fast, .True., no_intra, (/50*Abs(params%max_nrg)/), (/spc2,mol2,0/))

    If (.Not.nrg_calc_flag) Then
      ! we have serious problems ??, 
      Write(*,*) "biasfactor,params%max_nrg",biasfactor, params%max_nrg
      Write(*,*) spc2, mol2
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    Return

  End Subroutine brmoves_idflip

!!$
!!$  !------------------------------------------------------------------
!!$  ! deletes atoms from ptr onwards, for the purpose of idflip
!!$  ! ptr warning : ptr should point to memory in sorbates(spc2,mol2), 
!!$  ! not some other random place
!!$  !------------------------------------------------------------------
!!$  Subroutine brmoves_idflipDEL(sorbates, scell, imodel, spc,mol,&
!!$      ptr, biasfactor, nrg_calc_flag, max_nrg, rti)
!!$
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Integer, Intent(in)   :: spc,mol
!!$    Real(kind=RDbl), Intent(inout)   :: biasfactor
!!$    Logical, Intent(out)   :: nrg_calc_flag
!!$    Type(NodeType), Pointer       :: ptr
!!$    Type(Simcell_Params), Intent(in) :: scell
!!$    Type(Interaction_Model), Pointer         :: imodel       
!!$    Real(kind=RDbl) , Intent(in) :: max_nrg,rti
!!$    Integer :: nsorbs,k
!!$    Real(kind=RDbl)      :: nxtIntra, intra, next_atom_nrg, bias, ufinal
!!$    Real(kind=RDbl) ,Dimension(2,Size(sorbates,1)) :: tmpNrgArr, nrg_arr 
!!$
!!$    nsorbs=molecules_getnsorbs()
!!$    ufinal=zero
!!$    intra=zero
!!$    nrg_arr(1:2,1:nsorbs)=zero
!!$    biasfactor=one
!!$
!!$    Call  brmoves_deletenode(sorbates, scell, imodel, spc,mol, &
!!$        ptr%C, bias, next_atom_nrg, nrg_calc_flag, max_nrg, rti, nxtIntra, &
!!$        tmpNrgArr)
!!$    If (nrg_calc_flag) Then
!!$      biasfactor =  biasfactor * bias 
!!$      intra=intra+nxtIntra
!!$      Do k=1,nsorbs
!!$        nrg_arr(1,k)=tmpNrgArr(1,k)
!!$        ufinal=ufinal+nrg_arr(1,k)
!!$      End Do
!!$      ufinal=ufinal+nxtIntra
!!$
!!$    Else
!!$      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$      Stop
!!$    Endif
!!$
!!$  End Subroutine brmoves_idflipDEL
!!$
!!$
!!$  !------------------------------------------------------------------
!!$  ! adds atoms from ptr onwards, for the purpose of idflip
!!$  ! ptr warning : ptr should point to memory in sorbates(spc,mol), 
!!$  ! not some other place
!!$  !------------------------------------------------------------------
!!$  Subroutine brmoves_idflipADD(sorbates, scell, imodel, spc,mol,&
!!$      ptr, biasfactor, nrg_calc_flag, max_nrg, rti)
!!$
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Integer, Intent(in)   :: spc,mol
!!$    Real(kind=RDbl), Intent(inout)   :: biasfactor
!!$    Logical, Intent(out)   :: nrg_calc_flag
!!$    Type(NodeType), Pointer       :: ptr
!!$    Type(Simcell_Params), Intent(in) :: scell
!!$    Type(Interaction_Model), Pointer         :: imodel       
!!$    Real(kind=RDbl) , Intent(in) :: max_nrg,rti
!!$    Integer :: nsorbs, k
!!$    Real(kind=RDbl)      :: nxtIntra, intra, next_atom_nrg, bias, ufinal
!!$    Real(kind=RDbl) ,Dimension(2,Size(sorbates,1)) :: tmpNrgArr, nrg_arr 
!!$
!!$    nsorbs=molecules_getnsorbs()
!!$    ufinal=zero
!!$    intra=zero
!!$    nrg_arr(1:2,1:nsorbs)=zero
!!$    biasfactor=one
!!$
!!$    Call  brmoves_insertnode(sorbates, scell, imodel, spc,mol, &
!!$        ptr%C, bias, next_atom_nrg, nrg_calc_flag, max_nrg, rti, nxtIntra, &
!!$        tmpNrgArr)
!!$
!!$    If (nrg_calc_flag) Then
!!$      biasfactor =  biasfactor * bias 
!!$      intra=intra+nxtIntra 
!!$      Do k=1,nsorbs
!!$        nrg_arr(1,k)=tmpNrgArr(1,k)
!!$        ufinal=ufinal+nrg_arr(1,k)
!!$      End Do
!!$      ufinal=ufinal+nxtIntra 
!!$
!!$    Else
!!$      Return
!!$    Endif
!!$
!!$  End Subroutine brmoves_idflipADD

  !------------------------------------------------------------------
  ! deletes atoms from ptr onwards, for the purpose of idflip
  ! ptr warning : ptr should point to memory in sorbates(spc2,mol2), 
  ! not some other random place
  !------------------------------------------------------------------
  Subroutine brmoves_idflipDEL(sorbates, scell, imodel, spc,mol,&
      ptr, biasfactor, nrg_calc_flag, max_nrg, rti)

    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)   :: spc,mol
    Real(kind=RDbl), Intent(inout)   :: biasfactor
    Logical, Intent(out)   :: nrg_calc_flag
    Type(NodeType), Pointer       :: ptr
    Type(Simcell_Params), Intent(in) :: scell
    Type(Interaction_Model), Pointer         :: imodel       
    Real(kind=RDbl) , Intent(in) :: max_nrg,rti
    Integer :: nsorbs,k, natoms, atom_num
    Real(kind=RDbl)      :: nxtIntra, intra, next_atom_nrg, bias, ufinal
    Real(kind=RDbl)      :: second_atom_nrg, third_atom_nrg
    Real(kind=RDbl) ,Dimension(2,Size(sorbates,1)) :: tmpNrgArr, nrg_arr 
    Type(NodeType), Pointer       :: nextptr, insptr
    nsorbs=molecules_getnsorbs()
    natoms=molecules_getnatoms(spc)
    ufinal=zero
    intra=zero
    nrg_arr(1:2,1:nsorbs)=zero
    biasfactor=one

    ! idflip for natoms<4 is not well tested -shaji 29 sep 2003
    ! works only for linear alkanes numbered from 1-n
    nextptr=>ptr%C
    atom_num=nextptr%atom_num

    If (atom_num==2) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
      insptr=>nextptr
      Call brmoves_delete_atom2( scell, sorbates, imodel, spc, mol, &
          atom_num, rti, bias, second_atom_nrg , &
          nrg_calc_flag, tmpNrgArr,  max_nrg)
      If (.not.nrg_calc_flag) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      endif
      biasfactor=biasfactor*bias
      nextptr=>insptr%C
      If (Associated(insptr%C)) atom_num=insptr%C%atom_num
    Endif

    If (atom_num==3) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
      insptr=>nextptr
      Call brmoves_delete_atom3( scell, sorbates, imodel, spc, mol, &
          atom_num, rti, bias, third_atom_nrg , &
          nrg_calc_flag, nxtIntra, tmpNrgArr)
      If (.not.nrg_calc_flag) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      endif
      biasfactor=biasfactor*bias
      ! this might be pointing to a null pointer, that is fine
      nextptr=>insptr%C

      If (Associated(insptr%C)) atom_num=insptr%C%atom_num

    Endif

    If (Associated(nextptr)) Then
      Call  brmoves_deletenode(sorbates, scell, imodel, spc,mol, &
          nextptr, bias, next_atom_nrg, nrg_calc_flag, max_nrg, rti, &
          nxtIntra, tmpNrgArr)
      biasfactor =  biasfactor * bias 
    Endif

    If (nrg_calc_flag) Then
      intra=intra+nxtIntra
      Do k=1,nsorbs
        nrg_arr(1,k)=tmpNrgArr(1,k)
        ufinal=ufinal+nrg_arr(1,k)
      End Do
      ufinal=ufinal+nxtIntra

    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

  End Subroutine brmoves_idflipDEL


  !------------------------------------------------------------------
  ! adds atoms from ptr onwards, for the purpose of idflip
  ! ptr warning : ptr should point to memory in sorbates(spc,mol), 
  ! not some other place
  !------------------------------------------------------------------
  Subroutine brmoves_idflipADD(sorbates, scell, imodel, spc,mol,&
      ptr, biasfactor, nrg_calc_flag, max_nrg, rti)

    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)   :: spc,mol
    Real(kind=RDbl), Intent(inout)   :: biasfactor
    Logical, Intent(out)   :: nrg_calc_flag
    Type(NodeType), Pointer       :: ptr
    Type(Simcell_Params), Intent(in) :: scell
    Type(Interaction_Model), Pointer         :: imodel       
    Real(kind=RDbl) , Intent(in) :: max_nrg,rti
    Integer :: nsorbs, k, natoms, atom_num
    Real(kind=RDbl)      :: nxtIntra, intra, next_atom_nrg, bias, ufinal
    Real(kind=RDbl)      :: second_atom_nrg, third_atom_nrg
    Real(kind=RDbl) ,Dimension(2,Size(sorbates,1)) :: tmpNrgArr, nrg_arr 
    Type(NodeType), Pointer       :: nextptr, insptr

    nsorbs=molecules_getnsorbs()
    ufinal=zero
    intra=zero
    nrg_arr(1:2,1:nsorbs)=zero
    biasfactor=one

    ! idflip for natoms<4 is not well tested -shaji 29 sep 2003
    ! works only for linear alkanes numbered from 1-n
    nextptr=>ptr%C
    atom_num=nextptr%atom_num
    If (atom_num==2) Then
      insptr=>nextptr
      Call  brmoves_insert_atom2( scell, sorbates, imodel, spc, mol, &
          atom_num, rti, bias, second_atom_nrg ,nrg_calc_flag, &
          max_nrg, tmpNrgArr)
      If (.not.nrg_calc_flag) Return

      biasfactor=biasfactor*bias
      nextptr=>insptr%C
      If (Associated(insptr%C)) Then
        atom_num=insptr%C%atom_num
      Endif
    Endif

    If (atom_num==3) Then
      insptr=>nextptr
      Call  brmoves_insert_atom3(scell, sorbates, imodel, spc, mol, &
          atom_num, rti, bias, third_atom_nrg ,nrg_calc_flag, &
          max_nrg, intra, tmpNrgArr)
      If (.not.nrg_calc_flag) Return
      biasfactor=biasfactor*bias
      nextptr=>insptr%C
      If (Associated(insptr%C)) Then
        atom_num=insptr%C%atom_num
      Endif
    Endif


    If (Associated(nextptr)) Then
      Call  brmoves_insertnode(sorbates, scell, imodel, spc,mol, &
          nextptr, bias, next_atom_nrg, nrg_calc_flag, max_nrg, rti, &
          nxtIntra,tmpNrgArr)
      biasfactor =  biasfactor * bias 
    Endif

    If (nrg_calc_flag) Then
      intra=intra+nxtIntra 
      Do k=1,nsorbs
        nrg_arr(1,k)=tmpNrgArr(1,k)
        ufinal=ufinal+nrg_arr(1,k)
      End Do
      ufinal=ufinal+nxtIntra 

    Else
      Return
    Endif
  End Subroutine brmoves_idflipADD


  !------------------------------------------------------------------
  ! Inserts one node of the molecule and returns the bias of
  ! and energy of the insertion
  ! **biasfactor** is the factor to be used in the acceptance criteria 
  ! to remove the bias of this particular node-insertion. 
  ! It is essentially the inverse of our bias for this node-configuration
  ! Note: both the bond angle and dihedral potentials have to exist
  !------------------------------------------------------------------
  Recursive Subroutine brmoves_insertnode(sorbates, scell, imodel, sorbtype,&
      molec, ptr, biasfactor, ufinal, nrg_calc_flag, max_nrg, rti, intra, &
      nrg_arr)

    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)   :: sorbtype,molec
    Real(kind=RDbl), Intent(out)   :: biasfactor, ufinal
    Logical, Intent(out)   :: nrg_calc_flag
    Type(NodeType), Pointer       :: ptr
    Type(Simcell_Params), Intent(in) :: scell
    Type(Interaction_Model), Pointer         :: imodel       
    Real(kind=RDbl) , Intent(in) :: max_nrg,rti
    Real(kind=RDbl) ,Intent(out) :: intra !!!,coul, noncoul
    Real(kind=RDbl) , Intent(out), Dimension(:,:) :: nrg_arr

    Integer              :: k, nsorbs, index, atomno, parent_num
    Real(kind=RDbl)      :: lib_theta,theta_wt,temp_max_nrg
    Real(kind=RDbl)      :: ubend, intrapot, nxtIntra, next_atom_nrg, bias 
    Real(kind=RDbl) ,Dimension(2,Size(sorbates,1)) :: tmpNrgArr 

    Type(Angle_Library_Params), Pointer :: angle_lib   

    atomno=ptr%atom_num
    nsorbs = molecules_getnsorbs()
    ufinal=zero
    intra=zero
    nrg_arr(1:2,1:nsorbs)=zero
    biasfactor=one

    Call branchedcoords_placenode(ptr)
    parent_num=ptr%P%atom_num
    !** Obtains the angle library for this node, the mid atom for the
    !** angle is at ptr%P%atom_num
    Call brmoves_find_Anglelib(auxlibs_array(sorbtype), &
        parent_num, angle_lib) 

    !**  Retrieves a bond angle according to cumulative probability
    index = angle_getbiasindex(angle_lib,lib_theta,ubend,theta_wt)
    ptr%P%baPC = lib_theta
    ufinal = ufinal + ubend
    intra=intra+ubend
    If (ufinal>max_nrg) Then
      nrg_calc_flag=.False.
      return
    Endif
    temp_max_nrg=max_nrg-ufinal
    biasfactor =  biasfactor / (Real(angle_lib%num_ba)*theta_wt)
    nrg_calc_flag=.False.

    Call brmoves_dsgetbias(scell, sorbates, imodel, sorbtype, molec, &
        rti, ptr, bias, intrapot, tmpNrgArr, .True. , nrg_calc_flag, &
        temp_max_nrg)

    If (nrg_calc_flag) Then
      biasfactor =  biasfactor * bias 

      intra=intra+intrapot 
      Do k=1,nsorbs
        nrg_arr(1,k)=tmpNrgArr(1,k)
        ufinal=ufinal+nrg_arr(1,k)
      End Do
      ufinal=ufinal+intrapot

    Else
      Return
    Endif

    !** beauty of recursion
    If (Associated(ptr%C)) Then
      Call  brmoves_insertnode(sorbates, scell, imodel, sorbtype, molec, &
          ptr%C, bias, next_atom_nrg, nrg_calc_flag, max_nrg, rti, nxtIntra, &
          tmpNrgArr)
      If (nrg_calc_flag) Then
        ufinal=ufinal+next_atom_nrg
        intra=intra+nxtIntra

        Do k=1,nsorbs
          nrg_arr(2,k) =zero !coul
          nrg_arr(1,k) =nrg_arr(1,k)+tmpNrgArr(1,k) ! noncoul
        End Do
        biasfactor=biasfactor*bias
      Else
        Return
      Endif

    End If

  End Subroutine brmoves_insertnode

  !-----------------------------------------------------------
  ! Generates a center-of-mass of the molecule 
  ! uniformly from the simulation volume, and then using random eulerian 
  ! angles generates the  coordinates.  It also returns the bias of the
  ! insertion in "biasfactor"
  !-----------------------------------------------------------
  Subroutine brmoves_rinsert(sorbates,sorb,molec,params, biasfactor)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)                              :: sorb,molec
    Real(kind=RDbl), Intent(out)                     :: biasfactor
    Type(BranchedRIns_Params), Intent(in)            :: params

    Integer              :: natoms
    Real(kind=RDbl)      :: costheta
    Type(VecType)        :: startpos 

    ! ** bias factor = N_c *  bias  = N_c *  (1/ N_c) =one.
    biasfactor=one

    !** Get the position of the center of mass unifromly from the 
    !** sim.volume
    startpos = simcell_uniformRandPt(params%scell)

    !** 3.) Translate the pt. to the appropriate unit cell
    sorbates(sorb)%gcoords(molec)%branchedcoords%startpos = startpos 
    sorbates(sorb)%gcoords(molec)%branchedcoords%root%baPC = 0.0_RDbl 
    sorbates(sorb)%gcoords(molec)%branchedcoords%root%taC  = 0.0_RDbl 
    sorbates(sorb)%gcoords(molec)%branchedcoords%root%C%taC = 0.0_RDbl 

    !** If natoms is =1, we are done
    natoms=molecules_getnatoms(sorb)
    If (natoms==1) Then
      sorbates(sorb)%coords(1,molec)%rp=startpos 
      Return
    Endif

    !** Generate the eulerian angles at random
    ! Sample theta from a cosine distribution.
    costheta = 1 - rranf() * 2.0
    sorbates(sorb)%gcoords(molec)%branchedcoords%root%baPC = Acos(costheta)
    sorbates(sorb)%gcoords(molec)%branchedcoords%root%taC  = pi - rranf()*twopi
    sorbates(sorb)%gcoords(molec)%branchedcoords%root%C%taC = pi - rranf()*twopi

    !** Generate the xyz coordinates
    Call branchedcoords_place(sorbates(sorb)% &
        gcoords(molec)%branchedcoords%root) 
    Call branchedcoords_toxyz  &
        (sorbates(sorb)%gcoords(molec)%branchedcoords, &
        sorbates(sorb)%coords(1:natoms,molec)%rp)
  End Subroutine brmoves_rinsert

  !-----------------------------------------------------------------
  ! This routine gets the biasfactor for deleting a molecule
  ! This is the opposite move for "brmoves_cbinsert"
  ! **biasfactor** is the factor to be used in the deletion-acceptance 
  ! criteria to account for the bias used during the insertion of this 
  ! molecule. It is essentially equal to our bias for this configuration, 
  ! see note at end of Module
  !           subint -- subset interactions (molecule-system)
  !-----------------------------------------------------------------
  Subroutine brmoves_cbdelete(sorbates, sorbtype, molec, &
      params, biasfactor, subint,nrg_calc_flag)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)   ::   sorbtype,molec
    Type(BranchedCBDel_Params), Intent(in)   :: params
    Real(kind=RDbl), Intent(out)   :: biasfactor
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Logical, Intent(out)           :: nrg_calc_flag

    Real(kind=RDbl)                   :: ufinal, intra
    Real(kind=RDbl), Dimension(2,Size(sorbates))   :: nrg_arr 
    Integer              :: natoms, nsorbs, k
    Integer, Dimension(sorbates(sorbtype)%natoms)   :: atom_num 
    Real(kind=RDbl)   :: first_atom_nrg, second_atom_nrg, third_atom_nrg
    Real(kind=Rdbl)   :: fourth_atom_nrg, bias,rti,tmpIntra
    Real(kind=RDbl)   :: temp_max_nrg, ustored, ufresh
    Logical           :: fast, no_intra
    Real(kind=RDbl) ,Dimension(2,Size(sorbates,1)) :: tmpNrgArr 
    Type(NodeType), Pointer       :: root_ptr

    If  (firsttime_in_brmoves) &
        Call brmoves_basicInit(params%scell,DEF_indent)

    If (.Not.auxlibs_initialized(sorbtype)) Then
      Call brmoves_initAuxLibs(&
          auxlibs_array(sorbtype), sorbates(sorbtype), sorbtype, &
          subint%imodel, params%tk,DEF_indent)
    Endif


    !** For easy refernce to the co-ordinates of molecule
    !** make sure everything is correct with the f90-pointer business here
    root_ptr=> sorbates(sorbtype)%gcoords(molec)%branchedcoords%root

    !** Set the temperature parameter and other initializations
    rti = params%rti
    natoms = molecules_getnatoms(sorbtype)
    nsorbs = molecules_getnsorbs()
    ufinal=zero
    biasfactor=one
    fast=.true.
    intra=zero
    nrg_arr(1:2,1:nsorbs)=zero
    nrg_calc_flag=.False.

    temp_max_nrg= MAX_STOP_NRG

    ! some cavity bias things may be used during delete, the atoms 
    ! occupied be the existing molecule should not be included in the 
    ! cavity definition
    If (cavity_initialized) Then
      Call cavitylist_getCubeIndices(cavity, &
          sorbates(sorbtype), molec, natoms, cubes2avoid)
      avoid_cubes=.True.
    Endif

    !**
    !** Get the Bias of the ATOM-1
    !**
    atom_num(1)=root_ptr%atom_num
    Call  brmoves_delete_atom1(params%bmap, params%scell, subint%imodel, &
        sorbates, sorbtype, molec, atom_num(1), rti, bias, first_atom_nrg , &
        nrg_calc_flag ,tmpNrgArr)
    biasfactor = biasfactor * bias
    ufinal=ufinal+first_atom_nrg
    Do k=1,nsorbs
      nrg_arr(1,k)=tmpNrgArr(1,k)
    End Do

    If(natoms>1) Then
      !**
      !** Get the Bias of the ATOM-2
      !**
      atom_num(2)= root_ptr%C%atom_num  
      Call brmoves_delete_atom2( params%scell, sorbates, subint%imodel, &
          sorbtype, molec, atom_num(2), rti, bias, second_atom_nrg , &
          nrg_calc_flag, tmpNrgArr, params%max_nrg)
      biasfactor=biasfactor*bias
      ufinal=ufinal+second_atom_nrg
      Do k=1,nsorbs
        nrg_arr(1,k)=nrg_arr(1,k)+tmpNrgArr(1,k)
      End Do

      If (Associated(sorbates(sorbtype)%gcoords(molec) &
          %branchedcoords%root%C%L)) Then 
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__ 
        Stop
      End If

    Endif


    If(natoms>2) Then
      !**
      !** Get the Bias of the ATOM-3
      !**
      atom_num(3) = root_ptr%C%C%atom_num  
      Call   brmoves_delete_atom3( params%scell, sorbates, subint%imodel, sorbtype,&
          molec, atom_num(3), rti, bias, third_atom_nrg ,&
          nrg_calc_flag, tmpIntra,tmpNrgArr)
      ufinal=ufinal+third_atom_nrg
      biasfactor=biasfactor*bias
      intra=tmpIntra
      Do k=1,nsorbs
        nrg_arr(1,k)=nrg_arr(1,k)+tmpNrgArr(1,k)
      End Do
    endif

    !**
    !** delete atoms 4 onwards
    !**
    If(natoms>3) Then
      If(Associated(root_ptr%C%C%C)) Then
        Call brmoves_deletenode(sorbates, params%scell, subint%imodel, &
            sorbtype, molec, root_ptr%C%C%C, bias, fourth_atom_nrg, &
            nrg_calc_flag, params%max_nrg, params%rti, tmpIntra, tmpNrgArr )
        If (nrg_calc_flag) Then
          ufinal=ufinal+fourth_atom_nrg
          biasfactor=biasfactor*bias
          intra=intra+tmpIntra
          Do k=1,nsorbs
            nrg_arr(1,k)=nrg_arr(1,k)+tmpNrgArr(1,k)
          End Do
        Else
          Write(*,*) "DELETE_WARNING-something wrong"
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          stop
        Endif
      End If
    endif

    If (nrg_calc_flag) Then  
      ! recalculate the energy to stick to general music policy of 
      ! returning the energies in subint, recalc_flag=.False.
      no_intra=.False.
      Call brmoves_xyz_pbc(sorbates, params%scell, molec, &
          sorbtype, natoms)    

      ! see notes in the header to see the reason for increasing max_nrg
      nrg_calc_flag = subinteract_int(subint, sorbates, params%scell, &
          fast, .False., no_intra, (/50*Abs(params%max_nrg)/), &
          (/sorbtype,molec,0/))

      If (.Not.nrg_calc_flag) Then
        ! we have serious problems ??
        ! maybe due to MAX_STOM_NRG, max_nrg conflict
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    endif

    Return


  End Subroutine brmoves_cbdelete

  !------------------------------------------------------------------
  ! Deletes one node of the molecule and returns the bias of
  ! and energy of the deletion 
  ! Note: both the bond angel and dihedral potentials have to exist
  !------------------------------------------------------------------
  Recursive Subroutine brmoves_deletenode(sorbates, scell, imodel, sorbtype,&
      molec, ptr, biasfactor,ufinal,nrg_calc_flag,  max_nrg, rti, intra, &
      nrg_arr)
    !    Type(Subset_Interactions), Intent(InOut)       :: subint  
    Type(Interaction_Model), Pointer         :: imodel    
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)   :: sorbtype,molec
    Real(kind=RDbl), Intent(out)   :: biasfactor, ufinal
    Logical, Intent(out)   :: nrg_calc_flag
    Type(NodeType), Pointer       :: ptr
    Type(Simcell_Params), Intent(in) :: scell
    Real(kind=RDbl) , Intent(in) :: rti
    Real(kind=RDbl), Intent(out)       :: intra !!,coul,noncoul
    Real(kind=RDbl), Intent(out), Dimension(:,:) :: nrg_arr
    Real(kind=RDbl) , Intent(in) :: max_nrg

    Integer              :: k, nsorbs, index, atomno, parent_num
    Real(kind=RDbl)      :: temp_max_nrg,ubend,theta_wt,bias
    Real(kind=RDbl)      :: next_atom_nrg,tmpIntra
    Real(kind=RDbl),Dimension(2,Size(sorbates,1))    :: tmpNrgArr
    Type(Angle_Library_Params), Pointer :: angle_lib   


    nsorbs = molecules_getnsorbs()
    Call branchedcoords_placenode(ptr)
    biasfactor=one
    ufinal=zero
    intra=zero
    nrg_arr(1:2,1:nsorbs)=zero

    !If nrg more than this then exponential is essentially zero
    !such values dont add anything to parttion function
    !SDEBUG
    temp_max_nrg= MAX_STOP_NRG
!!$    temp_max_nrg= max_nrg
    !SDEBUG
    atomno = ptr%atom_num 
    parent_num=ptr%P%atom_num
    Call brmoves_find_Anglelib(auxlibs_array(sorbtype), &
        parent_num, angle_lib)

    !    index = ptr%P%angle_index  
    index=angle_indexFromVal(angle_lib,ptr%P%baPC)
    ubend = angle_lib%angle(index)%u
    theta_wt=angle_lib%angle(index)%wt

    ufinal = ufinal + ubend
    intra=intra+ubend
    biasfactor =  biasfactor*Real(angle_lib%num_ba)*theta_wt 

    Call brmoves_dsgetbias(scell, sorbates, imodel, sorbtype, molec, &
        rti, ptr, bias,  tmpIntra, tmpNrgArr, .False. , nrg_calc_flag, &
        temp_max_nrg)

    If (nrg_calc_flag) Then
      biasfactor =  biasfactor * bias 
      intra=intra+tmpIntra
      Do k=1,nsorbs
        nrg_arr(1,k)=tmpNrgArr(1,k)
        ufinal=ufinal+nrg_arr(1,k)
      End Do
      ufinal=ufinal+tmpIntra
    Else
      Write(*,*) "DELETE_WARNING-something wrong"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      !!      Stop
    Endif

    If (Associated(ptr%C)) Then
      Call brmoves_deletenode(sorbates, scell, imodel, sorbtype, &
          molec, ptr%C, bias, next_atom_nrg, nrg_calc_flag, &
          max_nrg, rti, tmpIntra, tmpNrgArr )

      intra=intra+tmpIntra
      ufinal=ufinal+next_atom_nrg
      Do k=1,nsorbs
        nrg_arr(1,k)=nrg_arr(1,k)+tmpNrgArr(1,k)
      End Do
      biasfactor=biasfactor*bias

    End If

  End Subroutine brmoves_deletenode

  !-----------------------------------------------------------------
  ! This routine gets the bias-weight for deleting a molecule
  ! with no bias , this value is actually one.
  !-----------------------------------------------------------------
  Subroutine brmoves_rdelete(sorbates,params, biasfactor)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Real(kind=RDbl), Intent(out)                     :: biasfactor
    Type(BranchedRDel_Params), Intent(in)                    :: params

    !** unbiased 
    biasfactor = one
    Return
  End Subroutine brmoves_rdelete

  !-------------------------------------------------------------
  ! Translates the linear molecule and generates
  ! first atom is translated randomly, keeping other gencoords same
  ! xyzcoords.
  !           subint -- subset interactions (molecule-system)
  !-------------------------------------------------------------
  Subroutine brmoves_randomtranslate(params, sorbates,sorb,molec,succ_flag, &
      subint)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)                              :: sorb,molec
    Type(BranchedRTrans_Params), Intent(in) :: params
    Logical, Intent(out)                             :: succ_flag
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Type(VecType)       :: deltadisp, startpos 
    Real(kind=RDbl)     :: deltatrans, xdelta, ydelta, zdelta
    Integer :: natoms
    Logical :: fast, skip_intra_flag

    natoms=molecules_getnatoms(sorb)
    deltatrans = params%deltatrans

    !** Calculate the energy of molecule at its current position
    !** calculate intra nrg also
    fast = .True.
    skip_intra_flag = .False.

    ! see notes in the header to see the reason for increasing max_nrg
    succ_flag = subinteract_int(subint,sorbates,params%scell,fast, &
        .False.,skip_intra_flag,(/50*Abs(params%max_nrg)/),(/sorb,molec,0/))
    If (.Not.succ_flag) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Problem with nrg calculation of an existing structure"
      Stop
    Endif

    If (params%constrained) Then
      !** Do a coordinate transformation to the constrained plane
      !** system where the normal is along the x-axis
      !***DEBUG STUFF
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) " is this part initialised?"
      Stop

      startpos = sorbates(sorb)%gcoords(molec)%branchedcoords%startpos 
      startpos = params%constrtrans*startpos 

      !** Get the displacement in the y-z plane
      ydelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      zdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      deltadisp = (/0.0_RDbl, ydelta, zdelta/)
      startpos = startpos + deltadisp 

      !** Tranform back to the original coordinate system
      sorbates(sorb)%gcoords(molec)%branchedcoords%startpos = &
          params%constrinvtrans*startpos 
    Else
      xdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      ydelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      zdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      deltadisp = (/xdelta, ydelta, zdelta/)
      sorbates(sorb)%gcoords(molec)%branchedcoords%startpos =&
          sorbates(sorb)%gcoords(molec)%branchedcoords%startpos + deltadisp
    End If

    !** Generate the xyz coordinates
    Call branchedcoords_place(sorbates(sorb)% &
        gcoords(molec)%branchedcoords%root) 
    Call branchedcoords_toxyz  &
        (sorbates(sorb)%gcoords(molec)%branchedcoords, &
        sorbates(sorb)%coords(1:natoms,molec)%rp)

    !** do periodic b. condns before energy calcs
    !** Update the other arrays in species using PBCs
    Call simcell_pbc(params%scell, &
        sorbates(sorb)%coords(1:natoms,molec)%rp, &
        sorbates(sorb)%coords(1:natoms,molec)%r, &
        sorbates(sorb)%coords(1:natoms,molec)%cr)

    succ_flag = subinteract_int(subint,sorbates,params%scell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/sorb,molec,0/))

  End Subroutine brmoves_randomtranslate


  !-----------------------------------------------------------------
  ! This routine shifts the center of mass of a molecule
  ! with generalized coordinates "gcoords" by amount "deltacom" and 
  ! generates the "xyzcoords"
  !----------------------------------------------------------------
  Subroutine brmoves_fixedtranslate(sorbates,sorb,molec,params,deltacom)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)                              :: sorb,molec
    Type(BranchedFTrans_Params), Intent(in) :: params
    Type(VecType), Intent(in) :: deltacom

    Integer :: natoms

    natoms=molecules_getnatoms(sorb)
    !** Displace the molecule
    sorbates(sorb)%gcoords(molec)%branchedcoords%startpos =&
        sorbates(sorb)%gcoords(molec)%branchedcoords%startpos + deltacom

    !** Generate the new coordinates
    Call branchedcoords_place(sorbates(sorb)% &
        gcoords(molec)%branchedcoords%root) 
    Call branchedcoords_toxyz  &
        (sorbates(sorb)%gcoords(molec)%branchedcoords, &
        sorbates(sorb)%coords(1:natoms,molec)%rp)

  End Subroutine brmoves_fixedtranslate

  !-------------------------------------------------------------
  ! Rotate the molecule by generating new eulerian angles and 
  ! generates xyzcoords.
  !-------------------------------------------------------------
  Subroutine brmoves_randomrotate(sorbates,sorb,molec,params)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)                              :: sorb,molec
    Type(BranchedRRot_Params), Intent(in) :: params

    Integer             :: natoms
    Real(kind=RDbl)     :: deltarot, phi, psi
    Real(kind=RDbl)     :: costheta

    !** Check the no. of atoms and return if no. of atoms is 1
    natoms = molecules_getnatoms(sorb)
    If (natoms == 1) Return

    deltarot = params%deltarot

    !** Sample phi and psi uniformly between -pi,pi but theta
    !** needs to be sampled from a cosine distribution
    !** See Allen and Tildesley, Pg. 132-133, "Molecular Liquids"
    phi = sorbates(sorb)%gcoords(molec)%branchedcoords%root%taC +  &
        ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
    sorbates(sorb)%gcoords(molec)%branchedcoords%root%taC = phi - &
        Anint(phi/twopi)*twopi
    psi = sorbates(sorb)%gcoords(molec)%branchedcoords%root%C%taC + &
        ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
    sorbates(sorb)%gcoords(molec)%branchedcoords%root%C%taC = psi - &
        Anint(psi/twopi)*twopi

    costheta = Cos(sorbates(sorb)%gcoords(molec)%branchedcoords%root%baPC) + &
        ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
    costheta = costheta - Anint(costheta/2.0_RDbl)*2.0_RDbl
    sorbates(sorb)%gcoords(molec)%branchedcoords%root%baPC = Acos(costheta)

    !** Generate the new coordinates
    Call branchedcoords_place(sorbates(sorb)% &
        gcoords(molec)%branchedcoords%root) 
    Call branchedcoords_toxyz  &
        (sorbates(sorb)%gcoords(molec)%branchedcoords, &
        sorbates(sorb)%coords(1:natoms,molec)%rp)
  End Subroutine brmoves_randomrotate

  !-------------------------------------------------------------
  ! Returns a pointer to the angle_lib for node at  atom_num.
  ! the mid atom of the ablge is atom_num
  !-------------------------------------------------------------
  Subroutine brmoves_find_Anglelib(auxlib, atom_num, angle_lib)
    Type(AuxLibraries),Intent(inout) :: auxlib
    Type(Angle_Library_Params), Pointer :: angle_lib
    Integer, Intent(In)        :: atom_num
    Logical :: found

    Integer  :: angle_num
    angle_num= auxlib%angle_lib_list(atom_num)
    angle_lib => auxlib%angle_lib(angle_num)

  End Subroutine brmoves_find_Anglelib




  !-------------------------------------------------------------
  ! Returns a pointer (frame_lib) to the frame_lib with atom_num.
  !-------------------------------------------------------------
  Subroutine brmoves_find_InsFramelib(params,atom_num,frame_lib,found)
    Type(BranchedCBIns_Params), Intent(in) :: params
    Type(Frame_Library_Params), Pointer :: frame_lib
    Integer, Intent(In)        :: atom_num
    Integer, Intent(Out)       :: found  ! 1 if successful, 0 if not

    Integer  :: i

    found = 0
!!$
!!$    Do i=1, params%framelib_num
!!$      If(params%frame_lib(i)%atom_num == atom_num) Then
!!$        frame_lib => params%frame_lib(i)
!!$        found = 1
!!$      End If
!!$    End Do
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    stop
  End Subroutine brmoves_find_Insframelib


  !-------------------------------------------------------------
  ! Returns a pointer (angle_lib) to the frame_lib with atom_num.
  !-------------------------------------------------------------
  Subroutine brmoves_find_DelAnglelib(params,atom_num,angle_lib)
    Type(BranchedCBDel_Params), Intent(in) :: params
    Type(Angle_Library_Params), Pointer :: angle_lib
    Integer, Intent(In)        :: atom_num
    Logical       :: found 

    Integer  :: i

    found = .false.
!!$
!!$    Do i=1, params%anglelib_num
!!$      If(params%angle_lib(i)%atom_num == atom_num) Then
!!$        angle_lib => params%angle_lib(i)
!!$        found = .true.
!!$        exit
!!$      End If
!!$    End Do
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    stop
    If (.Not.found) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Could not find the angle_lib "
      Stop
    End If


  End Subroutine brmoves_find_DelAnglelib

  !-------------------------------------------------------------
  ! Returns a pointer (frame_lib) to the frame_lib with atom_num.
  !-------------------------------------------------------------
  Subroutine brmoves_find_DelFramelib(params,atom_num,frame_lib,found)
    Type(BranchedCBDel_Params), Intent(in) :: params
    Type(Frame_Library_Params), Pointer :: frame_lib
    Integer, Intent(In)        :: atom_num
    Integer, Intent(Out)       :: found  ! 1 if successful, 0 if not

    Integer  :: i

    found = 0
!!$
!!$    Do i=1, params%framelib_num
!!$      If(params%frame_lib(i)%atom_num == atom_num) Then
!!$        frame_lib => params%frame_lib(i)
!!$        found = 1
!!$      End If
!!$    End Do
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    stop
  End Subroutine brmoves_find_DelFramelib


  !-------------------------------------------------------------
  ! display params for fixed-translation moves
  !-------------------------------------------------------------
  Subroutine brmoves_displayFTrans(params,unitno)
    Type(BranchedFTrans_Params),Intent(in) :: params
    Integer,Intent(in)                                :: unitno
    Write(unitno,*) "      Branched Molecule Fixed-Translation Parameters :"
    If (params%constrained) Then
      Write(unitno,*) "      CONSTRAINED by matrix : "
      Call matrix_display(params%constrtrans,'f8.3',unitno)
      Write(unitno,*) "      The vector normal to the constarined plain is : "
      Call vector_filedisplay(params%normal)
    Endif
    If (params%scale_trans) Then
      Write(unitno,'(a)') "      Scaling of Translation length allowed."
      Write(unitno,'(a,f12.5)')"      Initial translation length : ",&
          params%deltatrans
    Else
      Write(unitno,'(a)') "      Scaling of Translation length not allowed."
      Write(unitno,'(a,f12.5)')"      Initial translation length : ",&
          params%deltatrans
    Endif
    Write(unitno,*)
  End Subroutine brmoves_displayFTrans



  !-------------------------------------------------------------
  ! display params for fixed-translation moves
  !-------------------------------------------------------------
  Subroutine brmoves_displayIDFlip(params,indent,unitno)
    Type(BranchedIDFlip_Params),Intent(in) :: params
    Integer,Intent(In)                          :: indent,unitno

    Character(len=indent)       :: blank
    Character(len=strLen)       :: string
    blank = Repeat(' ',indent)

    Write(unitno,'(2a)') blank, "Branched Molecule IDFlip Parameters :"
    Write(unitno,'(2a)') blank, "Bias map file    :"//Trim(params%biasfilename)
  End Subroutine brmoves_displayIDFlip

  !-------------------------------------------------------------
  ! display params for random-translation moves
  !-------------------------------------------------------------
  Subroutine brmoves_displayRTrans(params,indent,unitno)
    Type(BranchedRTrans_Params),Intent(in) :: params
    Integer,Intent(in)                                :: indent, unitno
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string

    blank = Repeat(' ',indent)
    Write(unitno,'(2a)') blank, "Branched Molecule Rondom-Translation Parameters :"
    If (params%constrained) Then
      Write(unitno,*) "      CONSTRAINED by matrix : "
      Call matrix_display(params%constrtrans,'f8.3',unitno)
      Write(unitno,*) "      The vector normal to the constarined plain is : "
      Call vector_filedisplay(params%normal)
    Endif
    If (params%scale_trans) Then
      Write(unitno,'(2a)') blank," Scaling of Translation length allowed."
      Write(unitno,'(2a,f12.5)')blank, " Initial translation length : ",&
          params%deltatrans
    Else
      Write(unitno,'(2a)') blank, " Scaling of Translation length not allowed."
      Write(unitno,'(2a,f12.5)') blank, " Initial translation length : ",&
          params%deltatrans
    Endif
    Write(unitno,*)
  End Subroutine brmoves_displayRTrans


  !-------------------------------------------------------------
  ! display params for Random-Rotation moves
  !-------------------------------------------------------------
  Subroutine brmoves_displayRRot(params,unitno)
    Type(BranchedRRot_Params),Intent(in) :: params
    Integer,Intent(in)                                :: unitno
    Write(unitno,*) "      Branched Molecule Rondom-Rotation Parameters :"
    If (params%scale_rot) Then
      Write(unitno,'(a)') "      Scaling of Rotation increment allowed."
      Write(unitno,'(a,f12.5)')"      Initial rotation-delta : ",&
          params%deltarot
    Else
      Write(unitno,'(a)') "      Scaling of Rotation increment not allowed."
      Write(unitno,'(a,f12.5)')"      Initial rotation-delta : ",&
          params%deltarot
    Endif
    Write(unitno,*) 
  End Subroutine brmoves_displayRRot

  !-------------------------------------------------------------
  ! display params for biased insertion moves
  !-------------------------------------------------------------
  Subroutine brmoves_displayCBIns(params,unitno)
    Type(BranchedCBIns_Params),Intent(in) :: params
    Integer ,Intent(in)                                :: unitno
    Write(unitno,*) "      Branched Molecule Biased-Insertion Parameters :"
    Write(unitno,*) "      Bias map file    :"//Trim(params%biasfilename)
    Write(unitno,*)
  End Subroutine brmoves_displayCBIns

  !-------------------------------------------------------------
  ! display params for biased insertion moves
  !-------------------------------------------------------------
  Subroutine brmoves_displayReGrow(params,unitno)
    Type(BranchedReGrow_Params),Intent(in) :: params
    Integer ,Intent(in)                                :: unitno
    Write(unitno,*) "      Branched Moxlecule Biased-Insertion Parameters : "
    Write(unitno,*) "      Bias map file    : "//Trim(params%biasfilename)
    Write(unitno,*)
  End Subroutine brmoves_displayReGrow

  !-------------------------------------------------------------
  ! display params for deletion moves
  !-------------------------------------------------------------
  Subroutine brmoves_displayCBDel(params,unitno)
    Type(BranchedCBDel_Params),Intent(in) :: params
    Integer,Intent(in)                :: unitno
    Write(unitno,*) "      Branched Molecule Deletion Parameters :"
    Write(unitno,*) "      Bias map file    :"//Trim(params%biasfilename)
    Write(unitno,*)
  End Subroutine brmoves_displayCBDel

  !------------------------------------------------------------
  ! Adjust the jump length to try and keep the acceptance ratio
  ! close to 50%
  !------------------------------------------------------------
  Subroutine brmoves_adjustdeltatrans(params, ratio)
    Type(BranchedRTrans_Params),Intent(inout) :: params
    Real(kind=RDbl), Intent(in)       :: ratio

    If (.Not. params%scale_trans) Return

    !** these values should not be har coded here, should go to control file
    If (ratio < 0.49) Then
      params%deltatrans = Max(params%deltatrans*0.95_RDbl, 0.01_RDbl)
      !      Write(delta_unit,*) params%deltatrans 
    Else If (ratio > 0.51) Then
      params%deltatrans = Min(params%deltatrans*1.05_RDbl, 1.0_RDbl)
      !      Write(delta_unit,*) params%deltatrans 
    Endif
  End Subroutine brmoves_adjustdeltatrans


  !------------------------------------------------------------
  ! adjusts the rotational angle increment to maintain 50% acceptance ratio
  !------------------------------------------------------------
  Subroutine brmoves_adjustdeltarot(params, ratio)
    Type(BranchedRRot_Params),Intent(inout) :: params
    Real(kind=RDbl), Intent(in)       :: ratio

    If (.Not. params%scale_rot) Return

    !** these values should not be hard coded here, should go to control file
    If (ratio < 0.49) Then
      params%deltarot = Max(params%deltarot*0.95_RDbl, 0.01_RDbl)
    Else If (ratio > 0.51) Then
      params%deltarot = Min(params%deltarot*1.05_RDbl, 1.0_RDbl)
    Endif
  End Subroutine brmoves_adjustdeltarot


  !-----------------------------------------------------------
  ! all atoms above atomn number="atomno" are removed and regrown.
  ! The "biasfactor" will remove all the biases used
  ! nrg_calc_flag=.FALSE. if nrg is too high or could not be calculated
  ! nrg_arr(i,j,k) -> i=1 old, i=2 new, j=1 ncoul, j=2 coul, k=1,nsorbs
  ! intra(1) old intra, intra(2) new intra
  !           subint -- subset interactions (molecule-system)
  !-----------------------------------------------------------
  Subroutine brmoves_regrow(sorbates,sorbtype,molec,atomno,params, & 
      biasfactor,nrg_calc_flag, subint)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)                              :: sorbtype,molec,atomno
    Type(BranchedReGrow_Params), Pointer             :: params
    Real(kind=RDbl), Intent(out)                     :: biasfactor
    Logical, Intent(inout)                           :: nrg_calc_flag   
    Type(Subset_Interactions), Intent(InOut)         :: subint    

    Integer :: k,natoms,nsorbs,ins_atom,del_atom
    Real(kind=RDbl) :: forward_bias, reverse_bias, bias,temp_max_nrg
    Real(kind=RDbl) :: intranrg,usecond,uthird,unode, ufinal
    Real(kind=RDbl),Dimension(2,Size(sorbates,1)) :: tmpNrgArr
    Real(kind=RDbl), Dimension(2,2,Size(sorbates,1)) :: nrg_arr
    Real(kind=RDbl), Dimension(2,NO_OF_INTRA_POTS) :: intra
    Type(NodeType), Pointer       :: root_ptr,node_ptr,delete_ptr,insert_ptr
    Logical :: found, fast, skip_intra_flag

    Real(kind=RDbl) :: new_nrgT, new_nrgI, old_nrgT, old_nrgI

    If  (firsttime_in_brmoves) &
        Call brmoves_basicInit(params%scell,DEF_indent)


    If (.Not.auxlibs_initialized(sorbtype)) Then
      Call brmoves_initAuxLibs(&
          auxlibs_array(sorbtype), sorbates(sorbtype), sorbtype, &
          subint%imodel, params%tk,DEF_indent)
    Endif

    !** default return values
    biasfactor=one
    ufinal=zero

    !** Calculate the energy of molecule at its current position, recalc=false
    fast = .True.
    skip_intra_flag = .False.

    ! see notes in the header to see the reason for increasing max_nrg
    nrg_calc_flag = subinteract_int(subint,sorbates,params%scell,fast, &
        .False.,skip_intra_flag,(/50*Abs(params%max_nrg)/), &
        (/sorbtype,molec,0/))

    If (.not.nrg_calc_flag) then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    endif

    forward_bias=one ! following the typical monte carlo notation
    reverse_bias=one ! see notes at end of module

    natoms=molecules_getnatoms(sorbtype)
    nsorbs=molecules_getnsorbs()
    intra(1:2,NO_OF_INTRA_POTS)=zero
    nrg_arr(1:2,1:2,1:nsorbs)=zero


    ! some cavity bias things may be used during delete, the atoms 
    ! occupied be the existing molecule should not be included in the 
    ! cavity definition
    If (cavity_initialized) Then
      Call cavitylist_getCubeIndices(cavity, &
          sorbates(sorbtype), molec, natoms, cubes2avoid)
      avoid_cubes=.True.
    Endif

    temp_max_nrg=MAX_STOP_NRG
!!$    temp_max_nrg=params%max_nrg

    !** Save the old coordinates
    Call config_copymolec(sorbates(sorbtype), molec, params%oldcoords, 1)

    !** get the root_ptr and node_ptr at atom=atomno
    root_ptr=> sorbates(sorbtype)%gcoords(molec)%branchedcoords%root
    Call branchedcoords_findnode(root_ptr,atomno,node_ptr,found)

    If (.Not.found) Then
      Write(*,*) "Could not find the node for atom number: ",atomno
      Stop
    Endif


    !** Suggest a new configuration, get its energy and bias 
    !** Copy current config to oldcoords 
    !** ATOM-2 
    !**
    If ( (atomno==1) .And. (natoms>=2) ) Then
      ins_atom=2
      If (root_ptr%C%atom_num/=ins_atom) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      Call  brmoves_insert_atom2( params%scell, sorbates, subint%imodel, &
          sorbtype, molec, ins_atom, params%rti, bias, usecond , &
          nrg_calc_flag, temp_max_nrg, tmpNrgArr)
      If (nrg_calc_flag) Then
        forward_bias=forward_bias*bias
        Do k=1,nsorbs
          nrg_arr(2,1,k)=nrg_arr(2,1,k)+tmpNrgArr(1,k)
        End Do
      Else
        Return
      Endif
    Endif

    !**
    !** ATOM-3 
    !**
    If ( (atomno<=2) .And. (natoms>=3) ) Then
      ins_atom=3
      If (root_ptr%C%C%atom_num/=ins_atom) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      Call  brmoves_insert_atom3( params%scell, sorbates, subint%imodel, &
          sorbtype, molec, ins_atom, params%rti, bias, uthird ,nrg_calc_flag, &
          temp_max_nrg, intranrg, tmpNrgArr)

      If (nrg_calc_flag) Then
        forward_bias=forward_bias*bias
        Do k=1,nsorbs
          nrg_arr(2,1,k)=nrg_arr(2,1,k)+tmpNrgArr(1,k)
        End Do
        intra(2,TOTAL_INDEX)=intra(2,TOTAL_INDEX)+intranrg
      Else
        Return
      Endif
    Endif

    !**
    !** ATOM-4 Onwards recursive insertion
    !**
    If (natoms>=4) Then
      If (atomno<=3) Then
        ins_atom=4
        If (root_ptr%C%C%C%atom_num/=ins_atom) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
        insert_ptr=>root_ptr%C%C%C
      Else
        ins_atom=node_ptr%C%atom_num
        insert_ptr=>node_ptr%C
      Endif

      Call brmoves_insertnode(sorbates,  params%scell, subint%imodel, &
          sorbtype, molec,  insert_ptr, bias, unode, nrg_calc_flag, &
          temp_max_nrg, params%rti, intranrg, tmpNrgArr)

      If (nrg_calc_flag) Then
        forward_bias=forward_bias*bias
        ufinal=zero
        Do k=1,nsorbs
          nrg_arr(2,1,k)=nrg_arr(2,1,k)+tmpNrgArr(1,k)
          ufinal=ufinal+ nrg_arr(2,1,k)
        End Do
        intra(2,TOTAL_INDEX)=intra(2,TOTAL_INDEX)+intranrg
        ufinal=ufinal+ intra(2,TOTAL_INDEX)
        If (ufinal>params%max_nrg) Then
!!$          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$          Write(*,*) "legitimate insertion exceeds max_nrg, &
!!$              & consider increasing the max_nrg"
          nrg_calc_flag=.False.
          Return
        Endif
        !** Calculate the energy of molecule at its current position,
        ! see notes in the header to see the reason for increasing max_nrg
        nrg_calc_flag = subinteract_int(subint,sorbates,params%scell,fast, &
            .True.,skip_intra_flag,(/50*Abs(params%max_nrg)/),&
            (/sorbtype,molec,0/))
        If (.Not.nrg_calc_flag) Then
          ! we have serious problems ??
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif

      Else
        Return
      Endif

    Else
      !** Calculate the energy of molecule at its current position,
      ! see notes in the header to see the reason for increasing max_nrg
      nrg_calc_flag = subinteract_int(subint,sorbates,params%scell,fast, &
          .True.,skip_intra_flag,(/50*Abs(params%max_nrg)/),&
          (/sorbtype,molec,0/))
      If (.Not.nrg_calc_flag) Then
        ! we have serious problems ??
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Endif

    !** Save the suggested coordinates, so that they wont be erased 
    !** during the deletion trials
    Call config_copymolec(sorbates(sorbtype),molec, params%suggested,1)

    !** Copy the oldcoords back to sorbates so that the deletion bias of 
    !**existing structure can be calculated    
    Call config_copymolec(params%oldcoords,1,sorbates(sorbtype),molec)

    !**
    !** ATOM-2
    !**
    If ( (atomno==1) .And. (natoms>=2) ) Then
      del_atom=2
      If (root_ptr%C%atom_num/=del_atom) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      Call  brmoves_delete_atom2( params%scell, sorbates, subint%imodel, &
          sorbtype,  molec, del_atom, params%rti, bias, usecond , &
          nrg_calc_flag, tmpNrgArr, params%max_nrg)

      If (nrg_calc_flag) Then
        reverse_bias=reverse_bias*bias
        Do k=1,nsorbs
          nrg_arr(1,1,k)=nrg_arr(1,1,k)+tmpNrgArr(1,k)
        End Do
      Else
        Write(*,*) "DELETE_WARNING-something wrong"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Endif


    !**
    !** ATOM-3
    !**
    If ( (atomno<=2) .And. (natoms>=3) ) Then
      del_atom=3
      If (root_ptr%C%C%atom_num/=del_atom) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      Call brmoves_delete_atom3( params%scell, sorbates, subint%imodel, &
          sorbtype, molec, del_atom, params%rti, bias, uthird,&
          nrg_calc_flag,  intranrg, tmpNrgArr)

      If (nrg_calc_flag) Then
        reverse_bias=reverse_bias*bias
        Do k=1,nsorbs
          nrg_arr(1,1,k)=nrg_arr(1,1,k)+tmpNrgArr(1,k)
        End Do
        intra(1,TOTAL_INDEX)=intra(1,TOTAL_INDEX)+intranrg
      Else
        Write(*,*) "DELETE_WARNING-something wrong"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Endif


    !**
    !** ATOM-4 onwards recursive deletion
    !**
    If (natoms>=4) Then
      If (atomno<=3) Then

        del_atom=4
        If (root_ptr%C%C%C%atom_num/=del_atom) Then
          Write(*,*) del_atom
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
        delete_ptr=>root_ptr%C%C%C

      Else

        del_atom=node_ptr%C%atom_num
        delete_ptr=>node_ptr%C

      Endif
      Call brmoves_deletenode(sorbates, params%scell, subint%imodel, &
          sorbtype, molec, delete_ptr, bias, unode, nrg_calc_flag, &
          params%max_nrg, params%rti, intranrg, tmpNrgArr )

      If (nrg_calc_flag) Then
        reverse_bias=reverse_bias*bias
        Do k=1,nsorbs
          nrg_arr(1,1,k)=nrg_arr(1,1,k)+tmpNrgArr(1,k)
        End Do
        intra(1,TOTAL_INDEX)=intra(1,TOTAL_INDEX)+intranrg
      Else
        Write(*,*) "DELETE_WARNING-something wrong"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Endif

    !** final bias
    biasfactor=reverse_bias*forward_bias

    !   !** Finally copy the suggested coords back to sorbates
    Call config_copymolec(params%suggested,1,sorbates(sorbtype),molec)

  End Subroutine brmoves_regrow


  !------------------------------------------------------------
  ! Inserts the first atom of chain, returns energy and bias factor 
  ! biasfactor = 1/(bias for that particular position [from bmap])
  !------------------------------------------------------------
  Subroutine brmoves_insert_atom1(bmap, scell,imodel, sorbates,sorbtype,&
      molec, atomno, rti, biasfactor, ufirst ,nrg_calc_flag, max_nrg ,nrg_arr)
    Type(BiasMap_Params), Intent(in) :: bmap
    Type(Simcell_Params), Intent(in) :: scell
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer , Intent(in)             :: sorbtype, molec, atomno
    Real(kind=RDbl),Intent(in)       :: rti
    Real(kind=RDbl),Intent(out)      :: biasfactor,ufirst
    Logical,Intent(out)              :: nrg_calc_flag
    Real(kind=RDbl),Intent(in)       :: max_nrg
    Real(kind=RDbl),Dimension(:,:),Intent(out) :: nrg_arr
    !    Type(Subset_Interactions), Intent(InOut)         :: subint    
    Type(Interaction_Model), Pointer         :: imodel    
    Integer         :: index, ncubelets, icellx, icelly, icellz
    Integer         :: natoms, first_atom_num
    Type(VecType)   :: startpos
    Logical         :: fast, insertflag
    Type(NodeType), Pointer       :: root_ptr

    ! decide HLFAS or LLFAS
    If (hlfasON) Then
      ! do nothing
    Else
      ! check hlfas is asked for
      If (hlfas)Then
        If (general_getCurrentSimno()>=hlfas_simno) Then
          hlfasON=.True.
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "Shifting from LLFAS to HLFAS "
          Write(*,*) " ----- "
        Endif
      Endif
    Endif
    !** for easy refernce to the coord of the molecule
    root_ptr => sorbates(sorbtype)%gcoords(molec)%branchedcoords%root

    nrg_calc_flag=.False.
    biasfactor=one
    ufirst=zero
    fast=.True.

    ! Place the first atom  
    Call branchedcoords_unplace( root_ptr)    
    Call branchedcoords_placenode( root_ptr)    
    first_atom_num = root_ptr%atom_num  
    If (first_atom_num/=atomno) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    insertflag=.True.
    ! generate a new position and find its bias
    If (hlfasON) Then
      Call brmoves_HLFASbias(scell, sorbates, bmap, imodel, sorbtype, &
      molec, rti, root_ptr, biasfactor, insertflag, nrg_calc_flag, max_nrg)
    Else
      Call brmoves_LLFASbias(scell, bmap, sorbates, sorbtype, molec, &
          root_ptr, biasfactor, insertflag)
    Endif

    !** Now calculate the energy of its interaction with all other sorbs
    Call brmoves_get_atom_sorbNrg(sorbates, scell, imodel, sorbtype,&
        molec, atomno, fast, ufirst , nrg_calc_flag,nrg_arr, max_nrg)
!    Write(*,*) nrg_calc_flag, nrg_arr, ufirst

  End Subroutine brmoves_insert_atom1


  !-------------------------------------------------------------------------
  ! given a node pointer and bmap, does an HLFAS sampling 
  ! biasfactor will have two factors f1 anf f2
  ! f1 comes from map bias of the point
  ! f2 comes from rosenbluth wt (based on full sorb-sorb energy) 
  ! If move is a delete move, returns the bias of existing posn
  ! If move is insert, returns the bias energies and New positions 
  !--------------------------------------------------------------------------
  Subroutine brmoves_HLFASbias(scell, sorbates, bmap, imodel, sorbtype, &
      molec, rti, rootptr, biasfactor, insertflag, nrg_calc_flag, opt_max_nrg)

    Type(Simcell_Params), Intent(in) :: scell
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(BiasMap_Params), Intent(in) :: bmap
    Type(Interaction_Model), Pointer         :: imodel     
    Integer, Intent(in)                              :: sorbtype,molec
    Real(kind=RDbl), Intent(in)                      :: rti
    Type(NodeType), Pointer                          :: rootPtr
    Real(kind=RDbl) , Intent(out)                    :: biasfactor
    Logical,Intent(in)                               :: insertflag 
    Logical,Intent(out)                              :: nrg_calc_flag
    Real(kind=RDbl) , Intent(in),Optional            :: opt_max_nrg

    Logical :: yes_ins, yes_del, fast, temp_nrg_flag
    Integer :: i, k, atomno, trial_index
    Integer :: success_num, natoms, nsorbs, ncubelets
    Integer :: index, icellx, icelly, icellz
    Real(kind=RDbl) :: partition, max_nrg, bfac, f1, f2
    Real(kind=RDbl) :: u, ufirst, exponent_temp, wt
    Real(kind=RDbl),Dimension(2,Size(sorbates,1)) :: tmpNrgArr
    Real(kind=RDbl),Dimension(3) :: startPos
    Type(VecType) :: startposVec
    Type(VecType) :: orig_pos

    If (Present(opt_max_nrg)) Then
      max_nrg=opt_max_nrg
    Else
      max_nrg=VERYLRGMAXNRG
    Endif

    If (insertflag) Then
      yes_ins=.True.
    Else
      yes_ins=.False.
    Endif
    yes_del=.Not.(yes_ins)

    !** While deleting calculate the index corresponding to current theta 
    !** and psi in the array SSTrialArr
    If (yes_del) Then
      orig_pos= sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos
    Endif

    natoms=molecules_getnatoms(sorbtype)
    nsorbs=molecules_getnsorbs()
    fast=.True.
    partition = 0.0_RDbl
    success_num = 0       
    atomno=rootptr%atom_num

    ! ufirst is the total energy of the new atom that will be added
    biasfactor=one
    ufirst =zero
    nrg_calc_flag=.false.
    ncubelets  = bmap_getncubelets(bmap)

    !**build up the node ; use sphere sampling
    Do i=1,NHLtrials

      ! if deleting last one is current configuration
      If ((yes_del).and.(i==NHLtrials))then
        ! Map the position to the unit cell; convert the vector type to an array
        startpos = simcell_maptouc(scell, orig_pos)

        ! Get the weigth of cubelet where the COM lies
        ! note this bfac is inverse of insert bfac
        wt = bmap_getcellwt(bmap, startpos)
        bfac = (ncubelets*wt)

        sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos=orig_pos

      else
        ! get a trial position using bmap
        index = bmap_getbiasindex(bmap)
        bfac = 1.0_RDBl/bmap_getbiaswt(bmap,index)/ncubelets
        If (yes_del) bfac = 1/bfac
        startposVec%comp(1:3)   = bmap_getbiaspt(bmap, index)

        !** Select one of the unit cells in the simulation cell
        icellx = Int(rranf()*scell%nx)
        icelly = Int(rranf()*scell%ny)
        icellz = Int(rranf()*scell%nz)
        !** Translate the pt. to the appropriate unit cell
        startposVec = startposVec + &
            simcell_getcellorigin(scell, icellx, icelly, icellz)

        sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos  = &
            startposVec%comp
      endif


      ! convert to xyz xoords, do periodic bc
      Call brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms,atomno)
      temp_nrg_flag=.true.
      u=zero

      ! calculate the nrg of the second atom  , note max_nrg is the maximum 
      ! allowed at this point
      Call brmoves_get_atom_sorbNrg(sorbates,scell,imodel, sorbtype,&
          molec,atomno,fast,u,temp_nrg_flag,tmpNrgArr,max_nrg)

      ! make sure only good places are sampled, over-cautious here
      If ( (u>max_nrg) ) temp_nrg_flag=.False.

      If (temp_nrg_flag) Then
        success_num = success_num + 1     
        startposArray(success_num,1:3) = &
            sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos%comp(1:3)
        exponent_temp= -(u*rti)
        bmap_bias_array(success_num)=bfac
        !** make sure energy is not too low and exponent is not too big
        If (exponent_temp>default_MAX_EXP) Then
          exponent_temp=default_MAX_EXP
        Elseif (exponent_temp<-default_MAX_EXP) Then
          exponent_temp=-default_MAX_EXP
        Endif
        trial_e(success_num)=u

        Do k=1,nsorbs
          trial_ncoul_ar(success_num,k)=tmpNrgArr(1,k)
        End Do

        trial_exp_e(success_num) = Exp(exponent_temp)
        partition = trial_exp_e(success_num) + partition

      Else
        ! this should not have happened if the grid was fine enough
        ! so we will use the current value instead of the mid point of 
        ! the intervals. delCurrent will be set later..
        If ((yes_del).And.(i==NHLtrials)) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif

      Endif

    End Do !** end of Sphere sampling loop 

    If (success_num == 0) Then
      If (yes_ins) Then
        !** If all trials were bad, then return... Dead End!!
        nrg_calc_flag=.False.
        Return
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Else
      !** 'cumulate' and normalise the trial_exp_e array
      wt = 0.0_RDbl
      Do i=1,success_num
        wt = wt + trial_exp_e(i)/partition
        cum_prob(i) = wt 
      End Do
    Endif

    If (yes_ins) Then
      !** select the interval for the position of second atom
      trial_index = random_gettrial(1,success_num,cum_prob)             

      ! use the new position
      sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos%comp(1:3)&
          =startposArray(trial_index,1:3)

    Else
      trial_index=success_num ! the last one
      sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos =orig_pos
    Endif

    ! fix everything
    Call brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms)
    nrg_calc_flag=.true.
    If (yes_ins) then

      ! bmap bias
      f1=bmap_bias_array(trial_index)

      ! rosenbluth bias
      f2= partition/(Real(NHLTrials)* &
          trial_exp_e(trial_index))

      ! note bias factor anyway need values from midpoint only
      biasfactor =  f1*f2
    Else

      ! bmap bias
      f1=bmap_bias_array(success_num)

      ! rosenbluth bias
      f2= (   (Real(NHLTrials)) *   &
          trial_exp_e(success_num)  ) / partition

      ! note bias factor anyway need values from midpoint only
      biasfactor =  f1*f2

    Endif

  End Subroutine brmoves_HLFASbias


  !-------------------------------------------------------------------------
  ! given the rootptr, picks a position based on the biasmap
  ! usses only sorb-zeo interactions for biasing
  ! If move is a delete move, returns the bias of exisiting coord
  ! If move is insert, returns the bias and New positions for the next node
  ! , Note that the coordinates 
  ! are returned as part of "sorbates" structure
  !--------------------------------------------------------------------------
  Subroutine brmoves_LLFASbias(scell, bmap, sorbates, sorbtype, molec, &
      rootPtr, biasfactor, insertflag)

    Type(Simcell_Params), Intent(in) :: scell
    Type(BiasMap_Params), Intent(in) :: bmap
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)          :: sorbtype,molec
    Type(NodeType), Pointer                          :: rootPtr
    Real(kind=RDbl) , Intent(out)                    :: biasfactor
    Logical,Intent(in)                               :: insertflag 

    Logical :: yes_ins, yes_del
    Integer :: natoms, atomno, ncubelets, index, icellx, icelly, icellz
    Real(kind=RDbl) :: wt
    Real(kind=RDbl),Dimension(3)  :: startpos
    Type(VecType)  :: startposVec
    yes_ins=.False.
    If (insertflag) yes_ins=.True.
    yes_del=.Not.(yes_ins)

    natoms=molecules_getnatoms(sorbtype)
    atomno=rootptr%atom_num
    biasfactor=one

    ncubelets  = bmap_getncubelets(bmap)
    If (yes_del) Then

      ! Map the position to the unit cell; convert the vector type to an array
      startpos = simcell_maptouc(scell, &
          sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos)

      ! Get the weigth of cubelet where the COM lies
      wt = bmap_getcellwt(bmap, startpos)
      biasfactor = biasfactor* (ncubelets*wt)

    Else

      ! get a trial position using bmap
      index = bmap_getbiasindex(bmap)
      biasfactor = 1.0_RDBl/bmap_getbiaswt(bmap,index)/ncubelets
      startposVec%comp(1:3)   = bmap_getbiaspt(bmap, index)

      !** Select one of the unit cells in the simulation cell
      icellx = Int(rranf()*scell%nx)
      icelly = Int(rranf()*scell%ny)
      icellz = Int(rranf()*scell%nz)
      !** Translate the pt. to the appropriate unit cell
      startposVec = startposVec + &
          simcell_getcellorigin(scell, icellx, icelly, icellz)

      sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos  = &
          startposVec%comp
    Endif

    ! convert to xyz xoords, do periodic bc
    Call brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms,atomno)

  End Subroutine brmoves_LLFASbias


  !------------------------------------------------------------
  ! Inserts the second atom of chain, returns energy and bias factor 
  ! biasfactor = 1/(bias for that particular position over the other 
  !                           positions sampled using sphere sampling)
  !------------------------------------------------------------
  Subroutine brmoves_insert_atom2( scell, sorbates, imodel,sorbtype,&
      molec, atomno, rti, biasfactor, usecond ,nrg_calc_flag, &
      max_nrg, nrg_arr)
    Type(Simcell_Params), Intent(in) :: scell
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(Interaction_Model), Pointer         :: imodel    
    Integer , Intent(in)             :: sorbtype, molec, atomno
    Real(kind=RDbl),Intent(in)       :: rti 
    Real(kind=RDbl),Intent(out)      :: biasfactor,usecond
    Logical,Intent(out)              :: nrg_calc_flag
    Real(kind=RDbl),Intent(in)       :: max_nrg
    Real(kind=RDbl) , Dimension(:,:), Intent(out)    :: nrg_arr

    Integer                       :: k,second_atom_num
    Type(NodeType), Pointer       :: root_ptr

    !** for easy refernce to the coord of the molecule
    root_ptr => sorbates(sorbtype)%gcoords(molec)%branchedcoords%root

    Call branchedcoords_placenode(root_ptr%C)    
    second_atom_num = root_ptr%C%atom_num  
    usecond =zero

    If( second_atom_num/=atomno) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    ! Do disk sampling find a trial value of theta and psi for root_ptr%C
    ! If no success then nrg_calc_flag will be .False.
    Call brmoves_ssgetbias(scell, sorbates, imodel, sorbtype, molec, rti, &
        root_ptr%C, biasfactor, nrg_arr, .True. , nrg_calc_flag, max_nrg )

    If (nrg_calc_flag) Then
      Do k=1,Size(sorbates,1)
        usecond=usecond+nrg_arr(1,k)
      End Do
    Endif

  End Subroutine brmoves_insert_atom2



  !------------------------------------------------------------
  ! Inserts the third atom of chain, returns energy and bias factor 
  ! biasfactor = 1/(bias for that particular position over the other 
  !                           positions sampled using disk sampling)
  !------------------------------------------------------------
  Subroutine brmoves_insert_atom3( scell, sorbates, imodel, sorbtype,&
      molec, atomno, rti, biasfactor, uthird ,nrg_calc_flag, max_nrg, &
      intra, nrg_arr)
    Type(Simcell_Params), Intent(in) :: scell
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    !    Type(Subset_Interactions), Intent(InOut)         :: subint
    Type(Interaction_Model), Pointer         :: imodel    
    Integer , Intent(in)             :: sorbtype, molec, atomno
    Real(kind=RDbl),Intent(in)       :: rti 
    Real(kind=RDbl),Intent(out)      :: biasfactor,uthird
    Logical,Intent(out)              :: nrg_calc_flag
    Real(kind=RDbl),Intent(in)       :: max_nrg
    Real(kind=RDbl) , Intent(out)    :: intra
    Real(kind=RDbl),Dimension(:,:),Intent(out) :: nrg_arr    

    Integer         :: index, k, third_atom_num, parent_num, nsorbs
    Real(kind=RDbl) :: ubend, temp_max_nrg, bias, tempIntra, lib_theta
    Real(kind=RDbl) :: theta_wt
    Type(NodeType), Pointer       :: root_ptr
    Type(Angle_Library_Params), Pointer :: angle_lib 


    !** for easy refernce to the coord of the molecule
    root_ptr => sorbates(sorbtype)%gcoords(molec)%branchedcoords%root

    nsorbs=molecules_getnsorbs()
    biasfactor=one
    uthird=zero
    intra=zero
    nrg_arr(1:2,1:nsorbs)=zero

    !This should be ideally in _cbinsert
    If (Associated(root_ptr%C%L)) Then 
      !yet to be debugged
      !SDEBUG
      !  Call branchedmoves_insertframe(sorbates, nsorbs, sorbtype, molec, &
      !  root_ptr%C, & 
      !            params, biasfactor,ufinal,in_flag)            
      !SDEBUG

      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

    Call branchedcoords_placenode (root_ptr%C%C)

    ! parent is usually atom-2 of the chain
    parent_num=root_ptr%C%atom_num

    !** Obtains the angle library for this node, the angle is at atom-2
    Call brmoves_find_Anglelib(auxlibs_array(sorbtype), &
        parent_num, angle_lib) 

    !**  Retrieves a bond angle according to its cumulative probability   
    index = angle_getbiasindex(angle_lib,lib_theta,ubend,theta_wt)          
    root_ptr%C%baPC = lib_theta

    biasfactor =  biasfactor /(Real(angle_lib%num_ba)*theta_wt)

    uthird = uthird + ubend
    intra=ubend
    If (uthird>max_nrg) Then
      nrg_calc_flag=.False.
      Return
    Endif

    third_atom_num = root_ptr%C%C%atom_num  
    If (third_atom_num/=atomno) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif


    !** disk samling
    temp_max_nrg= max_nrg - uthird 
    Call brmoves_dsgetbias(scell, sorbates, imodel, sorbtype, molec, &
        rti, root_ptr%C%C, bias, tempIntra, nrg_arr, .True. , nrg_calc_flag, &
        temp_max_nrg)

    If (nrg_calc_flag) Then
      biasfactor=biasfactor*bias
      Do k=1,nsorbs
        uthird=uthird+nrg_arr(1,k)
      End Do
    Endif

  End Subroutine brmoves_insert_atom3

  !------------------------------------------------------------
  ! Gives the bias and energy of the first atom of chain
  ! biasfactor = (bias for that particular position [from bmap] ) 
  !------------------------------------------------------------
  Subroutine brmoves_delete_atom1(bmap, scell, imodel, sorbates, sorbtype,&
      molec, atomno, rti, biasfactor, ufirst ,nrg_calc_flag,nrg_arr )
    Type(BiasMap_Params), Intent(in) :: bmap
    Type(Simcell_Params), Intent(in) :: scell
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer , Intent(in)             :: sorbtype, molec, atomno
    Real(kind=RDbl),Intent(in)       :: rti
    Real(kind=RDbl),Intent(out)      :: biasfactor,ufirst
    Logical,Intent(out)              :: nrg_calc_flag
    Real(kind=RDbl) , Intent(out), Dimension(:,:) :: nrg_arr
    Type(Interaction_Model), Pointer         :: imodel    
    Integer         :: first_atom_num
    Logical         :: fast, insertflag
    Type(NodeType), Pointer       :: root_ptr

    !** for easy refernce to the coord of the molecule
    root_ptr => sorbates(sorbtype)%gcoords(molec)%branchedcoords%root

    nrg_calc_flag=.True.
    biasfactor=one
    ufirst=zero
    fast=.True.

    Call branchedcoords_unplace(root_ptr)    
    Call branchedcoords_placenode(root_ptr)

    first_atom_num = root_ptr%atom_num  
    If (first_atom_num/=atomno) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    insertflag=.False.

    ! find the bias of current position
    If (hlfasON) Then
      Call brmoves_HLFASbias(scell, sorbates, bmap, imodel, sorbtype, &
      molec, rti, root_ptr, biasfactor, insertflag, nrg_calc_flag)
     Else
      Call brmoves_LLFASbias(scell, bmap, sorbates, sorbtype, molec, &
          root_ptr, biasfactor, insertflag)
    Endif

    ! calculate energy
    Call brmoves_get_atom_sorbNrg(sorbates,scell,imodel, sorbtype,&
        molec,atomno,fast,ufirst, nrg_calc_flag, nrg_arr)

    If (.Not.nrg_calc_flag) Then
      Write(*,*) "DELETE_WARNING-something wrong"
      Write(*,*) "Problem with first atom deletion in cbgcmc, &
          & something terribly wrong.. "
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      !!      Stop
    Endif

  End Subroutine brmoves_delete_atom1


  !------------------------------------------------------------
  ! Gives the bias and energy of the 2nd atom of chain
  ! biasfactor = (bias for this particular position over the other 
  ! positions sampled using sphere sampling during its insertion )
  !------------------------------------------------------------
  Subroutine brmoves_delete_atom2( scell, sorbates, imodel, sorbtype,&
      molec, atomno, rti, biasfactor, usecond ,nrg_calc_flag,nrg_arr,max_nrg)
    Type(Simcell_Params), Intent(in) :: scell
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    !    Type(Subset_Interactions), Intent(InOut)         :: subint    
    Type(Interaction_Model), Pointer         :: imodel    
    Integer , Intent(in)             :: sorbtype, molec, atomno
    Real(kind=RDbl),Intent(in)       :: rti, max_nrg
    Real(kind=RDbl),Intent(out)      :: biasfactor,usecond
    Logical,Intent(out)              :: nrg_calc_flag
    Real(kind=RDbl) , Intent(out), Dimension(:,:) :: nrg_arr

    Integer         :: k,nsorbs, second_atom_num
    Real(kind=RDbl) :: temp_max_nrg

    Type(NodeType), Pointer       :: root_ptr

    !** for easy refernce to the coord of the molecule
    root_ptr => sorbates(sorbtype)%gcoords(molec)%branchedcoords%root

    nsorbs=molecules_getnsorbs()
    biasfactor=one
    usecond=zero
    nrg_arr(1:2,1:nsorbs)=zero
    second_atom_num = root_ptr%C%atom_num  
    If (second_atom_num/=atomno) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Call branchedcoords_placenode(root_ptr%C)    

    ! If nrg more than this then exponential is essentially zero
    ! such values dont add anything to parttion function
    !SDEBUG
    temp_max_nrg= MAX_STOP_NRG

    !SDEBUG
    ! Do disk sampling find a trial value of theta and psi for root_ptr%C
    ! If no success then nrg_calc_flag will be .False.
    Call brmoves_ssgetbias(scell, sorbates, imodel, sorbtype, molec, rti, &
        root_ptr%C, biasfactor, nrg_arr, .False. , nrg_calc_flag, &
        temp_max_nrg )

    If (nrg_calc_flag) Then
      Do k=1,Size(sorbates,1)
        usecond=usecond+nrg_arr(1,k)
      End Do
    Else
      Write(*,*) "DELETE_WARNING-something wrong"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

  End Subroutine brmoves_delete_atom2

  !------------------------------------------------------------
  ! Gives the bias and energy of the 3rd atom of chain
  ! biasfactor = (bias for this particular position over the other 
  ! positions sampled using disk sampling during its insertion )
  !------------------------------------------------------------
  Subroutine brmoves_delete_atom3( scell, sorbates, imodel, sorbtype,&
      molec, atomno, rti, biasfactor, uthird ,nrg_calc_flag, &
      intra, nrg_arr)
    Type(Simcell_Params), Intent(in) :: scell
    !    Type(Subset_Interactions), Intent(InOut)         :: subint
    Type(Interaction_Model), Pointer         :: imodel    
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer , Intent(in)             :: sorbtype, molec, atomno
    Real(kind=RDbl),Intent(in)       :: rti 
    Real(kind=RDbl),Intent(out)      :: biasfactor,uthird
    Logical,Intent(out)              :: nrg_calc_flag
    Real(kind=RDbl),Intent(out)     :: intra
    Real(kind=RDbl) , Intent(out), Dimension(:,:) :: nrg_arr
    Integer         :: k,nsorbs, third_atom_num, index, parent_num
    Real(kind=RDbl) :: theta_wt, temp_max_nrg, ubend, bias, tmpIntra

    Type(Angle_Library_Params), Pointer :: angle_lib 
    Type(NodeType), Pointer       :: root_ptr

    !** for easy refernce to the coord of the molecule
    root_ptr => sorbates(sorbtype)%gcoords(molec)%branchedcoords%root

    nsorbs=molecules_getnsorbs()
    nrg_calc_flag=.True.
    biasfactor=one
    uthird=zero
    intra=zero
    nrg_arr(1:2,1:nsorbs)=zero

    ! If nrg more than this then exponential is essentially zero
    ! such values dont add anything to parttion function
    !SDEBUG
    temp_max_nrg= MAX_STOP_NRG
!!$    temp_max_nrg= max_nrg
    !SDEBUG

    Call branchedcoords_placenode(root_ptr%C%C)

    third_atom_num = root_ptr%C%C%atom_num  
    If (third_atom_num/=atomno) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    ! parent atom in the branched coords, thats were the angle is
    !** Obtains the angle library for this node, the angle is at atom-2
    parent_num=root_ptr%C%atom_num
    Call brmoves_find_Anglelib(auxlibs_array(sorbtype), &
        parent_num, angle_lib) 

    !**  find the probability/nrg of the current angle at atom-2
    index = angle_indexFromVal(angle_lib,root_ptr%C%baPC)
    theta_wt=angle_lib%angle(index)%wt
    ubend=angle_lib%angle(index)%u
    uthird = uthird + ubend
    intra=intra + ubend
    biasfactor =  biasfactor*real(angle_lib%num_ba)*theta_wt

    Call brmoves_dsgetbias(scell, sorbates, imodel, sorbtype, molec, &
        rti, root_ptr%C%C, bias,  tmpIntra, nrg_arr,  .False., nrg_calc_flag, &
        temp_max_nrg)

    If (nrg_calc_flag) Then
      biasfactor=biasfactor * bias
      Do k=1,nsorbs
        uthird=uthird+nrg_arr(1,k)
      End Do
      Return
    Else
      Write(*,*) "DELETE_WARNING-something wrong"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

  End Subroutine brmoves_delete_atom3

  !---------------------------------------------
  ! Initialise the arrays that contain the trial values for dihedral angles
  ! Psi and Phi vary from -pi to pi
  ! costheta varies from -1 to 1
  ! SS refers to sphere sampling, DS  refers to disk sampling
  !---------------------------------------------
  Subroutine brmoves_initTrialArrays(NPsiSS, NThetaSS, NPhiDS,indent)
    Integer , Intent(in) :: NPsiSS, NThetaSS, NPhiDS,indent
    Integer :: i,j,count
    Real(kind=RDbl) :: theta,costheta,dcostheta,psi,dpsi,phi,dphi
    Character(len=indent)     :: blank
    blank=Repeat(" ",indent)
    trial_arrays_inited=.True.
    !** Setting these variable which is global to this module
    NPsiTrialsSS   = NPsiSS
    NThetaTrialsSS = NThetaSS
    NPhiTrialsDS   = NPhiDS

    Write(*,'(2a)') blank,"Initialising The SSTrialArr "

    count=0
    dcostheta=(2*one)/NThetaTrialsSS
    dpsi=twopi/NPsiTrialsSS
    SSdcos=dcostheta
    SSdpsi=dpsi
#ifdef DEBUG
    Write(*,*) "dcostheta", dcostheta
    Write(*,*) "dpsi", dpsi
    Write(*,*) "NThetaTrialsSS,NPsiTrialsSS" ,NThetaTrialsSS,NPsiTrialsSS
#endif
    Write(*,'(2a)') blank, "Initializing The angle trial values "
    Do i=1,NPsiTrialsSS
      psi=(-pi)+(dpsi/2)+(i-1)*dpsi
      Do j=1,NThetaTrialsSS
        costheta=(-one)+(dcostheta/2)+(j-1)*dcostheta
        theta=Acos(costheta)
        count=count+1
        SSTrialArr(1,count)=costheta
        SSTrialArr(2,count)=theta
        SStrialArr(3,count)=psi
      End do
    End do
#ifdef DEBUG
    Write(*,*) "End Points cos   : ", SSTrialArr(1,1),SSTrialArr(1,count)
    Write(*,*) "End Points theta : ", SSTrialArr(2,1),SSTrialArr(2,count)
    Write(*,*) "End Points psi   : ", SSTrialArr(3,1),SSTrialArr(3,count)
#endif
    count=0
    dphi=twopi/NPhiTrialsDS
    DSdphi=dphi
#ifdef DEBUG
    Write(*,*) "dphi", dphi
    Write(*,*) "NPhiTrialsDS",NPhiTrialsDS
#endif
    Do i=1,NPhiTrialsDS
      phi=(-pi)+(dphi/2)+(i-1)*dphi
      count=count+1
      DSTrialArr(count)=phi
    End Do
#ifdef DEBUG
    Write(*,*) "End Points ", DSTrialArr(1),DSTrialArr(count)
#endif
  End Subroutine brmoves_initTrialArrays



  !-------------------------------------------------------------------------
  ! given a node pointer does a sphere sampling for placing the atom, 
  ! usses only sorb-sorb interactions for biasing
  ! If move is a delete move, returns the bias and energies of the 
  ! existing coordinate . If move is insert, returns the bias energies 
  ! and New positions for the next node, Note that the coordinates 
  ! are returned as part of "sorbates" structure
  ! This is a very tough routine to understand, Should be rewritten!!
  !--------------------------------------------------------------------------
  Subroutine brmoves_ssgetbias(scell, sorbates, imodel, sorbtype, molec, rti, &
      nodePtr, biasfactor, nrg_arr, insertflag,nrg_calc_flag,&
      opt_max_nrg)

    Type(Simcell_Params), Intent(in) :: scell
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    !    Type(Subset_Interactions), Intent(InOut)         :: subint
    Type(Interaction_Model), Pointer         :: imodel     
    Integer, Intent(in)                              :: sorbtype,molec
    Real(kind=RDbl), Intent(in)                      :: rti
    Type(NodeType), Pointer                          :: nodePtr
    Real(kind=RDbl) , Intent(out)                    :: biasfactor
    Real(kind=RDbl) , Intent(out), Dimension(:,:)    :: nrg_arr
    Logical,Intent(in)                               :: insertflag 
    Logical,Intent(out)                              :: nrg_calc_flag
    Real(kind=RDbl) , Intent(in),Optional            :: opt_max_nrg

    Logical :: yes_ins, yes_del, fast, delIntValIsBad, temp_nrg_flag
    Integer :: i,k,sstrials,current_index,delCurrent,atomno, trial_index
    Integer :: n_theta, n_psi, success_num, natoms, nsorbs
    Real(kind=RDbl) :: orig_theta, orig_psi, partition, costheta, max_nrg
    Real(kind=RDbl) :: u, usecond, exponent_temp, wt
    Real(kind=RDbl),Dimension(2,Size(sorbates,1)) :: tmpNrgArr

    If (Present(opt_max_nrg)) Then
      max_nrg=opt_max_nrg
    Else
      max_nrg=VERYLRGMAXNRG
    Endif

    If (insertflag) Then
      yes_ins=.True.
    Else
      yes_ins=.False.
    Endif
    yes_del=.Not.(yes_ins)

    !** While deleting calculate the index corresponding to current theta 
    !** and psi in the array SSTrialArr
    If (yes_del) Then
      orig_theta=nodePtr%P%baPC
      orig_psi=nodePtr%P%taC

      n_theta=Int((1+Cos(orig_theta))/SSdcos)+1
      n_psi=Int((orig_psi+pi)/SSdpsi)+1
      current_index=(n_psi - 1)*(NThetaTrialsSS)+n_theta

      ! this flag tells whether the midpoint of the region in which current 
      ! theta and phi lies succeds during nrg calculation. It should.., 
      ! If we Use a fine grid 
      delIntValIsBad=.False.

    Endif

    sstrials= NPsiTrialsSS * NThetaTrialsSS 
    natoms=molecules_getnatoms(sorbtype)
    nsorbs=molecules_getnsorbs()
    fast=.True.
    partition = 0.0_RDbl
    success_num = 0       
    atomno=nodePtr%atom_num

    ! usecond is the total energy of the new atom that will be added
    biasfactor=one
    usecond =zero
    nrg_calc_flag=.false.

    !**build up the node ; use sphere sampling
    Do i=1,sstrials
      !** The array SSTrialArr contains precalculated values of costheta
      !** theta and Psi
      costheta = SSTrialArr(1,i) 

      nodePtr%P%baPC = SSTrialArr(2,i) 
      nodePtr%P%taC =  SSTrialArr(3,i) 

      ! convert to xyz xoords, do periodic bc
      Call brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms,atomno)
      temp_nrg_flag=.true.
      u=zero

      ! calculate the nrg of the second atom  , note max_nrg is the maximum 
      ! allowed at this point
      Call brmoves_get_atom_sorbNrg(sorbates,scell,imodel, sorbtype,&
          molec,atomno,fast,u,temp_nrg_flag,tmpNrgArr,max_nrg)

      ! make sure only good places are sampled, over-cautious here
      If ( (u>max_nrg) ) temp_nrg_flag=.False.

      If (temp_nrg_flag) Then
        success_num = success_num + 1     
        theta1(success_num) = nodePtr%P%baPC  
        phi1(success_num) = nodePtr%P%taC 

        exponent_temp= -(u*rti)
        !** make sure energy is not too low and exponent is not too big
        If (exponent_temp>default_MAX_EXP) Then
          exponent_temp=default_MAX_EXP
        Endif
        trial_e(success_num)=u

        ! does 1 mean ncoul ?
        Do k=1,nsorbs
          trial_ncoul_ar(success_num,k)=tmpNrgArr(1,k)
        End Do

        trial_exp_e(success_num) = Exp(exponent_temp)
        partition = trial_exp_e(success_num) + partition

        ! delCurrent tells which is the interval in which our current 
        ! theta and psi are situated, this is the index in trial_e() array,
        ! current_index was the index in SSTrialArr array
        If ((yes_del).And.(i==current_index)) delCurrent=success_num

      Else

        ! this should not have happened if the grid was fine enough
        ! so we will use the current value instead of the mid point of 
        ! the intervals. delCurrent will be set later..
        If ((yes_del).And.(i==current_index)) Then
          delIntValIsBad=.True.
        Endif

      Endif

    End Do !** end of Sphere sampling loop 

    If (success_num == 0) Then
      If (yes_ins) Then
        !** If all trials were bad, then return... Dead End!!
        nrg_calc_flag=.False.
        Return
      Endif
    Else
      !** 'cumulate' and normalise the trial_exp_e array
      wt = 0.0_RDbl
      Do i=1,success_num
        wt = wt + trial_exp_e(i)/partition
        cum_prob(i) = wt 
      End Do
    Endif


    If (yes_ins) Then
      !** select the interval for the position of second atom
      trial_index = random_gettrial(1,success_num,cum_prob)             

      ! find a theta and psi in the vicinity using uniform sampling
      costheta=Cos(theta1(trial_index))+( 2*( rranf() ) - 1 )*SSdcos/2
      nodePtr%P%baPC = Acos(costheta)
      nodePtr%P%taC = phi1(trial_index)+( 2*( rranf() ) - 1 )*SSdpsi/2

    Else
      nodePtr%P%baPc = orig_theta
      nodePtr%P%taC = orig_psi
    Endif

    Call brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms)
    temp_nrg_flag=.true.
    Call brmoves_get_atom_sorbNrg(sorbates,scell,imodel,sorbtype,&
        molec,atomno,fast,u,temp_nrg_flag,tmpNrgArr,max_nrg)

    ! make sure only good places are sampled, over-cautious here
    If ( (u>max_nrg) ) temp_nrg_flag=.False.

    If (yes_ins) Then
      ! If the new point is bad then
      ! use the mid point as the trial point, and use the old energies
      If (temp_nrg_flag) Then
        usecond=u
        Do k=1,nsorbs
          nrg_arr(1,k)=tmpNrgArr(1,k)
        End Do
      Else
        nodePtr%P%baPC = theta1(trial_index)
        nodePtr%P%taC  = phi1(trial_index)
        usecond=trial_e(trial_index) 
        Call brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms)
        Do k=1,nsorbs
          nrg_arr(1,k)=trial_ncoul_ar(trial_index,k)
        End Do
      Endif

      ! note bias factor anyway need values from midpoint only
      biasfactor =  biasfactor * partition/(Real(sstrials)* &
          trial_exp_e(trial_index))
      nrg_calc_flag=.True.

    Else
      ! I guess we are doing deletions..
      If (temp_nrg_flag) Then
        ! if the interval in which current theta and phi are there 
        ! happen to be bad Then..
        If (delIntValIsBad) Then
          delCurrent=success_num+1
          exponent_temp= -(u*rti)
          ! make sure energy is not too low and exponent is not too big
          If (exponent_temp>default_MAX_EXP) Then
            exponent_temp=default_MAX_EXP
          Endif
          trial_e(delCurrent)=u
          trial_exp_e(delCurrent) =Exp(exponent_temp)
          partition=partition+trial_exp_e(delCurrent)
        Endif

        usecond=u
        Do k=1,nsorbs
          nrg_arr(1,k)=tmpNrgArr(1,k)
        end do
        nrg_calc_flag=.true.  
        biasfactor =  biasfactor * (   (Real(sstrials)) *   &
            trial_exp_e(delCurrent)  ) / partition

      Else
        Write(*,*) "DELETE_WARNING-something wrong"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "deletion not possible, something wrong!! "
        !!        Stop
      Endif

    Endif

  End Subroutine brmoves_ssgetbias

  !-------------------------------------------------------------------------
  ! intra pair potential of atomno with all previously placed atoms 
  !--------------------------------------------------------------------------
  Subroutine  brmoves_getintrapair(atoms_vec, spc, auxlib, atomno, nrg, &
      ip_calc_flag)
    Type(VecType), Dimension(:), Intent(in) :: atoms_vec
    Type(AuxLibraries), Intent(in)          :: auxlib
    Integer, Intent(in)                     :: spc,atomno
    Real(kind = RDbl), Intent(OUT)            :: nrg
    Logical,  Intent(OUT)            ::  ip_calc_flag

    Integer :: npots, p, other_atom, pot_index
    Real(kind = RDbl) :: A_lj ,B_lj ,r2i,r6i,r12i,r2,pot
    Type(VecType) :: rvec

    nrg=zero
    ip_calc_flag = .True.

    ! no of pots at that atom
    npots=auxlib%n_ip_pots(atomno)
    If (npots==0) Return 

    Do p=1,npots

      ! the other atom in the intra pair
      other_atom=auxlib%ip_pots_atomlist(atomno,p)

      rvec=atoms_vec(other_atom) - atoms_vec(atomno)
      ! index in the array where params are stored 
      pot_index=auxlib%ip_pots_index(atomno,p)

      A_lj=auxlib%ip_ABvals(pot_index,1)
      B_lj=auxlib%ip_ABvals(pot_index,2)

      pot = 0.0_RDbl

      !** Get the square of the distance
      r2 = vector_getnormsq(rvec)

      !**Check if it is within the cut-off radius
      If (r2 > CB_IP_HC2) Cycle ! dont return... carry on calculations
      If (r2 < CB_IP_LC2) Then  ! let's return...  energy is too high
        ip_calc_flag = .False.
        Return
      End If

      r2i = 1.0_RDbl/r2
      r6i = r2i*r2i*r2i
      r12i= r6i*r6i
      pot = A_lj*r12i - B_lj*r6i
      nrg=nrg+pot


    End Do

  End Subroutine brmoves_getintrapair

  !-------------------------------------------------------------------------
  ! given a node pointer does a disk-sampling for placing the atom. 
  ! If atomno==3 Takes only sorb-sorb interactions for biasing
  ! If atomno==4 Takes sorb-sorb interactions + Torsion nrg for biasing
  ! If atomno >4 Takes sorb-sorb interactions + Torsion nrg + Intrpair
  ! If move is a delete move returns the bias and energies of the 
  ! existing coordinate 
  ! If move is insert returns the bias energies and new phi-values for 
  ! the next node
  ! This is a very tough routine to understand, Should be rewritten!!
  !SDEBUG
  ! Break this into smaller movetypes... e.g. scell, sorbtype, molec, 
  ! natoms, nsorbs cosmodel, anglibarr should be part of this smaller 
  ! move Type SHAJI
  !SDEBUG
  !--------------------------------------------------------------------------
  Subroutine brmoves_dsgetbias(scell, sorbates, imodel, sorbtype, molec, &
      rti, nodePtr, biasfactor, intra, nrg_arr,insertflag,nrg_calc_flag,&
      opt_max_nrg)

    Type(Simcell_Params), Intent(in) :: scell
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(Interaction_Model), Pointer         :: imodel
    Integer, Intent(in)                              :: sorbtype,molec
    Real(kind=RDbl), Intent(in)                      :: rti
    Type(NodeType), Pointer                          :: nodePtr
    Real(kind=RDbl) , Intent(out)                    :: biasfactor,intra
    Real(kind=RDbl) , Intent(out), Dimension(:,:)    :: nrg_arr
    Logical,Intent(in)                               :: insertflag 
    Logical,Intent(out)                              :: nrg_calc_flag
    Real(kind=RDbl) , Intent(in),Optional            :: opt_max_nrg

    Logical :: yes_ins, yes_del, fast, delIntValIsBad, temp_nrg_flag, ip_calc_flag
    Integer :: i,j,k,dstrials,current_index,delCurrent,atomno, trial_index
    Integer :: success_num, natoms, nsorbs, MAX_P, tor_index
    Real(kind=RDbl) :: orig_phi, partition, max_nrg, temp_max_nrg
    Real(kind=RDbl) :: utors, uip, usscoul, utorsplusip, intrapot, ncoulnrg
    Real(kind=RDbl) :: ucurrent, totalu, exponent_temp, wt
    Real(kind=RDbl) :: costac, costerm
    Real(kind=RDbl),Dimension(2,Size(sorbates,1)) :: tmpNrgArr

    If (Present(opt_max_nrg)) Then
      max_nrg=opt_max_nrg
    Else
      max_nrg=VERYLRGMAXNRG
    Endif

    If (insertflag) Then
      yes_ins=.True.
    Else
      yes_ins=.False.
    Endif
    yes_del=.Not.(yes_ins)

    ! If  deleting,  calculate the index corresponding to current phi value  
    ! in the array DSTrialArr
    If (yes_del) Then
      orig_phi=nodePtr%P%taC
      current_index=Int((orig_phi+pi)/DSdphi)+1

      ! this flag tells whether the midpoint of the region in which current 
      ! phi lies succeds during nrg calculation. It should.., 
      ! If we Use a fine grid 
      delIntValIsBad=.False.

    Endif

    dstrials= NPhiTrialsDS
    natoms=molecules_getnatoms(sorbtype)
    nsorbs=molecules_getnsorbs()
    fast=.True.
    partition = 0.0_RDbl
    success_num = 0       
    atomno=nodePtr%atom_num
    MAX_P=COS_EXP_MAX_PARAMS

    ! usecond is the total energy of the new atom that will be added
    biasfactor=one
    ucurrent =zero
    intra=zero
    nrg_calc_flag=.False.

    !**build up the node ; use disk sampling
    Do i=1,dstrials
      temp_max_nrg=max_nrg
      !** The array DSTrialArr contains precalculated values of costheta
      !** theta and Psi
      nodePtr%P%taC =  DSTrialArr(i) 
      ! convert to xyz xoords, do periodic bc
      Call brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms,atomno)

      temp_nrg_flag=.true.
      utors=zero
      uip=zero
      usscoul=zero

      ! Torsion energy, HACK- assuming linear alkane
      If (atomno>3) Then
        tor_index= auxlibs_array(sorbtype)%tor_list(atomno) 
        costerm=one
        costac=Cos(nodePtr%P%taC)
        Do j=0,MAX_P-1
          utors= utors + auxlibs_array(sorbtype)%tor(tor_index,j)*costerm
          costerm=costerm*costac
        End Do
      Endif

      If (atomno>4) Then
        ip_calc_flag = .True.
        intrapot = zero 
        ! get intrapair potential with all atoms placed before.
        Call brmoves_getintrapair(&
            sorbates(sorbtype)%coords(1:natoms,molec)%rp, sorbtype, &
            auxlibs_array(sorbtype), atomno, intrapot, ip_calc_flag)
        If ( .Not. ip_calc_flag) temp_nrg_flag=.False.
        uip = intrapot
      End If

      utorsplusip=utors+uip
      If (utorsplusip>temp_max_nrg) Then
        temp_nrg_flag=.false.
      Else
        temp_max_nrg=temp_max_nrg-utorsplusip
      Endif

      If (.Not.temp_nrg_flag) Then
        If ((yes_del).And.(i==current_index)) Then
          delIntValIsBad=.True.
        Endif
        Cycle
      Endif

      Call brmoves_get_atom_sorbNrg(sorbates, scell, imodel, sorbtype,&
          molec, atomno, fast, ncoulnrg, temp_nrg_flag, tmpNrgArr, &
          temp_max_nrg)

      usscoul=ncoulnrg
      ! make sure only good places are sampled, over-cautious here
      If ( (usscoul>temp_max_nrg) ) temp_nrg_flag=.False.
      totalu=usscoul+utorsplusip
      If (temp_nrg_flag) Then
        success_num = success_num + 1     
        phi1(success_num) = nodePtr%P%taC 

        exponent_temp= -(totalu*rti)
        !** make sure energy is not too low and exponent is not too big
        If (exponent_temp>default_MAX_EXP) Then
          exponent_temp=default_MAX_EXP
        Elseif(exponent_temp<(-default_MAX_EXP)) Then
!          exponent_temp=(-default_MAX_EXP)
        Endif
        trial_e(success_num)=totalu
        trial_intra(success_num) = utorsplusip
        !SDEBUG ?????????????????? CORRECT ?
        Do k=1,nsorbs
          trial_ncoul_ar(success_num,k)=tmpNrgArr(1,k)
        End Do
        !SDEBUG
        trial_exp_e(success_num) = Exp(exponent_temp)
        partition = trial_exp_e(success_num) + partition

        ! delCurrent tells which is the interval in which our current 
        ! phi is situated, this is the index in trial_e() array,
        ! current_index was the index in DSTrialArr array
        If ((yes_del).And.(i==current_index)) delCurrent=success_num

      Else

        ! this should not have happened if the grid was fine enough
        ! so we will use the current value instead of the mid point of 
        ! the intervals. delCurrent will be set later..
        If ((yes_del).And.(i==current_index)) Then
          delIntValIsBad=.True.
        Endif

      Endif

    End Do !** end of disk sampling loop 

    If (success_num == 0) Then
      If (yes_ins) Then
        !** If all trials were bad, then return... Dead End!!
        nrg_calc_flag=.False.
        Return
      Endif
    Else
      !** 'cumulate' and normalise the trial_exp_e array
      wt = 0.0_RDbl
      Do i=1,success_num
        wt = wt + trial_exp_e(i)/partition
        cum_prob(i) = wt 
      End Do
    Endif
    If (yes_ins) Then
      !** select the interval for the position of second atom
      trial_index = random_gettrial(1,success_num,cum_prob)             
      ! find phi in the vicinity using uniform sampling
      nodePtr%P%taC = phi1(trial_index)+( 2*( rranf() ) - 1 )*(DSdphi/2)
    Else
      nodePtr%P%taC = orig_phi
    Endif

    Call brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms)
    temp_nrg_flag=.true.
    utors=zero
    uip=zero
    usscoul=zero

    If (atomno>3) Then
      tor_index= auxlibs_array(sorbtype)%tor_list(atomno) 
      costerm=one
      costac=Cos(nodePtr%P%taC)
      Do j=0,MAX_P-1
        utors= utors + auxlibs_array(sorbtype)%tor(tor_index,j)*costerm
        costerm=costerm*costac
      End Do
    Endif

    If (atomno>4) Then
      ip_calc_flag = .True.
      intrapot = zero 

      ! get intrapair potential with all atoms placed before.
      Call brmoves_getintrapair(&
          sorbates(sorbtype)%coords(1:natoms,molec)%rp, sorbtype, &
          auxlibs_array(sorbtype), atomno, intrapot, ip_calc_flag)
      If ( .Not. ip_calc_flag) Then
        temp_nrg_flag=.False.
        If (yes_del) Then
          Write(*,*) "DELETE_WARNING-something wrong"
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          !!          Stop
        Endif
      Endif
      uip = intrapot
    End If

    utorsplusip=utors+uip
    If ((utorsplusip>max_nrg).and.(yes_ins)) Then
      temp_nrg_flag=.False.
    Else
      temp_max_nrg= max_nrg-utorsplusip
    Endif

    If (temp_nrg_flag) Then
      ! calculate the nrg of the second atom  , note max_nrg is the maximum 
      ! allowed at this point
      Call brmoves_get_atom_sorbNrg(sorbates, scell, imodel, sorbtype,&
          molec, atomno, fast, ncoulnrg, temp_nrg_flag, tmpNrgArr, &
          temp_max_nrg)
    Endif

    usscoul=ncoulnrg
    ! make sure only good places are sampled, over-cautious here
    totalu=usscoul+utorsplusip

    If (yes_ins) Then
      ! If the new point is bad then
      ! use the mid point as the trial point, and use the old energies
      If (temp_nrg_flag) Then
        ucurrent=totalu
        intra=utorsplusip
        Do k=1,nsorbs
          nrg_arr(1,k)=tmpNrgArr(1,k)
        End Do
      Else
        nodePtr%P%taC  = phi1(trial_index)
        ucurrent=trial_e(trial_index) 
        intra=trial_intra(trial_index)
        Call brmoves_xyz_pbc(sorbates,scell,molec,sorbtype,natoms)
        Do k=1,nsorbs
          nrg_arr(1,k)=trial_ncoul_ar(trial_index,k)
        End Do
      Endif

      ! note bias factor anyway need values from midpoint only
      biasfactor =  biasfactor * partition/(Real(dstrials)* &
          trial_exp_e(trial_index))
      nrg_calc_flag=.True.

    Else
      ! I guess we are doing deletions..
      If (temp_nrg_flag) Then
        ! if the interval in which curren phi is there happen to be bad Then..
        If (delIntValIsBad) Then
          delCurrent=success_num+1
          exponent_temp= -(totalu*rti)
          ! make sure energy is not too low and exponent is not too big
          If (exponent_temp>default_MAX_EXP) Then
            exponent_temp=default_MAX_EXP
          Endif
          trial_e(delCurrent)=totalu
          trial_exp_e(delCurrent) =Exp(exponent_temp)
          partition=partition+trial_exp_e(delCurrent)
        Endif

        ucurrent=totalu
        intra=utorsplusip
        Do k=1,nsorbs
          nrg_arr(1,k)=tmpNrgArr(1,k)
        end do
        nrg_calc_flag=.true.  
        biasfactor =  biasfactor * (   (Real(dstrials)) *   &
            trial_exp_e(delCurrent)  ) / partition

      Else
        Write(*,*) "DELETE_WARNING-something wrong"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "deletion not possible, something wrong!! "
        Write(*,*) temp_max_nrg, max_nrg, utors, uip, ncoulnrg, ip_calc_flag
        !!        Stop
      Endif

    Endif
    !SDEBUG
    If (yes_del) Then
!!!      Write(*,'(a,i3,a,i3)') "DS Bias, molec : ", molec, "atom : ", atomno
!!!      Write(*,*) "partition ", partition
    Endif
    !SDEBUG
  End Subroutine brmoves_dsgetbias



End Module brmoves

!** bias factors used
!** state with N   molecules : m
!** state with N+1 molecules : n
!** alpha_mn       = probability of attemting a move from m to n
!** alpha_nm       = probability of attemting a move from n to m
!** rho_m, rho_n   = probability of the particular config dictated by 
!**                  statistical mechanics
!** Insertion Acceptance
!** P_Acc(m->n) = min[1, (rho_n X alpha_nm)/ ( rho_m X alpha_mn )
!** Deletion Acceptance
!** P_Acc(n->m) = min[1, (rho_m X alpha_mn)/ ( rho_n X alpha_nm )
!** biasfactor in insertions is part of alpha_mn, so it is in the denominator
!** biasfactor = 1/(all the bias towards that configuration)
!** biasfactor in deletions is part of alpha_mn, so it is in the nominator
!** biasfactor = (all the bias towards the configuration being deleted )



!---------------------------------------------------------------
! BELOW ARE ROUTINES WRITTEN BY YI FOR DEALING WITH BRANCHES IN 
! ALKANES, THEY ARE NOT TESTED, AND WONT COMPILE WITH THE CURRENT CODE.
! STORED HERE FOR FUTURE REFERNCE. DONT DELETE -Shaji
!---------------------------------------------------------------
!!$  !-----------------------------------------------------------
!!$  ! Inserts first atom of the molecule according to a
!!$  ! bias map and randomly chooses angular degrees of freedom.
!!$  ! It also returns the bias of the insertion
!!$  ! in "biasfactor"
!!$  !-----------------------------------------------------------
!!$  Subroutine brmoves_integral(sorbates, sorbtype, tk, ufinal, biasfactor) 
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates 
!!$    Integer, Intent(in)   :: sorbtype   
!!$    Real(kind=RDbl), Intent(in)    :: tk                  
!!$    Real(kind=RDbl), Intent(out)   :: ufinal, biasfactor  
!!$
!!$    Integer, Parameter   :: MAX_PARAMS = 6   
!!$    Integer              :: natoms, i, j, k,index, nx, ny, nz
!!$    Integer              :: icellx,icelly,icellz, ncubelets
!!$    Integer              :: num,low,high,mid,flag 
!!$    Integer              :: atomt1, atomt2, atomt3, atomt4   
!!$    Integer, Dimension(sorbates(sorbtype)%natoms)   :: atom_num  
!!$    Real(kind=RDbl),Dimension(0:MAX_PARAMS-1) :: cn    
!!$    Real(kind=RDbl)      :: costheta,ktheta,thetaeq,theta 
!!$    Real(kind=RDbl)      :: costac,costerm   
!!$    Real(kind=RDbl)      :: partition,key,wt,prob,rti 
!!$    Real(kind=RDbl)      :: u,energy_1
!!$    Real(kind=RDbl),Dimension(ba_ntrials*2) :: trial_e,ta_e,theta1,phi1,cum_prob
!!$    Type(VecType)        :: startpos
!!$    Logical              :: fast, mapflag
!!$
!!$    Type(MolecularParams), Pointer :: molecule
!!$
!!$    !** builds up the molecule:
!!$    Call branchedcoords_place(sorbates(sorbtype)%gcoords(1)%branchedcoords%root)
!!$
!!$    !** Set the temperature parameter
!!$    rti = 1.0_RDbl/(kcalmole_kb*tk)  
!!$  
!!$    If(Associated(sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%L)) Then 
!!$     Call brmoves_BrNdInteg(sorbates, sorbtype, &
!!$          sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C, &
!!$          rti, ufinal, biasfactor)
!!$     Return 
!!$    End If  
!!$
!!$    !** Retrieves the molecule pointer
!!$    Call molecules_getpointer(sorbtype,molecule)
!!$    !** Obtains the potential parameters
!!$    atomt1 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%atom_type
!!$    atomt2 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%atom_type
!!$    atomt3 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%atom_type
!!$!LC    ktheta = molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$!LC    thetaeq = molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$    
!!$    partition = 0.0_RDbl
!!$    Do i=1,ba_ntrials
!!$      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/real(ba_ntrials)
!!$      theta = Acos(costheta)
!!$      trial_e(i) = 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$      partition = exp(-trial_e(i)*rti) + partition
!!$    End Do
!!$
!!$    wt = 0.0_RDbl
!!$    Do i=1,ba_ntrials
!!$        wt = wt + exp(-trial_e(i)*rti)/partition
!!$        cum_prob(i) = wt
!!$    End Do
!!$
!!$    low = random_gettrial(1,ba_ntrials,cum_prob)
!!$    costheta = ( rranf() * 2.0_RDbl/real(ba_ntrials) ) - &
!!$         1.0_RDbl + (2.0_RDbl*(low-1))/real(ba_ntrials)
!!$    theta = Acos(costheta)
!!$    ufinal = 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$    biasfactor =  real(ba_ntrials)*exp(-trial_e(low)*rti)/partition 
!!$    sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%baPC = theta 
!!$
!!$    If (.NOT. Associated(sorbates(sorbtype)%gcoords(1)%branchedcoords% &
!!$        root%C%C%C)) Then 
!!$       Return 
!!$    End If 
!!$
!!$    !** Locate the fourth atom  
!!$    !** Chooses the bond angle
!!$    !** Retrieves the molecule pointer
!!$    Call molecules_getpointer(sorbtype,molecule)
!!$    !** Obtains the potential parameters
!!$    atomt1 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%atom_type
!!$    atomt2 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%atom_type
!!$    atomt3 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%C%atom_type
!!$!LC    ktheta = molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$!LC    thetaeq = molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$
!!$   !** Obtains the dihedral angle potential
!!$    atomt1 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%atom_type
!!$    atomt2 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%atom_type
!!$    atomt3 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%atom_type
!!$    atomt4 = sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%C%atom_type
!!$
!!$    Do i=0,MAX_PARAMS-1
!!$!LC     cn(i)= molecule%torsion%torparams(atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)
!!$    End Do
!!$
!!$    partition = 0.0_RDbl
!!$    Do i=1,ba_ntrials
!!$      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/real(ba_ntrials)
!!$      theta = Acos(costheta)
!!$      trial_e(i) = 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$
!!$     !Sample phi  from a cosine distribution.
!!$      phi1(i) = pi - rranf()*twopi
!!$      costac  = cos(phi1(i))  
!!$      costerm= 1.0_RDbl
!!$      ta_e(i) = 0.0_RDbl 
!!$        Do j=0,MAX_PARAMS-1
!!$          ta_e(i)= ta_e(i) + cn(j)*costerm
!!$          costerm=costerm*costac
!!$        End Do
!!$      trial_e(i) = trial_e(i)  + ta_e(i)      
!!$      partition = exp(-(trial_e(i))*rti) + partition
!!$    End Do
!!$
!!$    wt = 0.0_RDbl
!!$    Do i=1,ba_ntrials  
!!$      wt = wt + exp(-trial_e(i)*rti)/partition
!!$      cum_prob(i) = wt
!!$    End Do
!!$
!!$    low = random_gettrial(1,ba_ntrials,cum_prob)
!!$    costheta = ( rranf() * 2.0_RDbl/real(ba_ntrials) ) -  &
!!$               1.0_RDbl + (2.0_RDbl*(low-1))/real(ba_ntrials)
!!$    theta = Acos(costheta)
!!$    ufinal = ufinal + 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)  & 
!!$             + ta_e(low)  
!!$    biasfactor = biasfactor*real(ba_ntrials)*exp(-trial_e(low)*rti)/partition 
!!$    sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%baPC = theta
!!$    sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%taC  = phi1(low) 
!!$
!!$    If(Associated(sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%C%L)) Then 
!!$     Call brmoves_BrNdInteg(sorbates, sorbtype, &
!!$          sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%C, &
!!$          rti, ufinal, biasfactor)
!!$    Else If(Associated(sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%C%C)) Then 
!!$     Call brmoves_NodeIntegral(sorbates, sorbtype, &
!!$          sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%C%C, &
!!$          rti, ufinal, biasfactor)
!!$    End If
!!$
!!$   Return 
!!$  End Subroutine brmoves_integral 

!!$  !-----------------------------------------------------------
!!$  ! Calculate Z/omelga ratio from each individual node 
!!$  ! ufinal and biasfactor are returned: only deals with
!!$  ! two 
!!$  !-----------------------------------------------------------
!!$  Recursive Subroutine brmoves_BrNdInteg(sorbates, &
!!$            sorbtype,ptr,rti,ufinal,biasfactor)
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Type(NodeType), Pointer       :: ptr                   
!!$    Integer, Intent(in)   :: sorbtype   
!!$    Real(kind=RDbl), Intent(in)    :: rti                   
!!$    Real(kind=RDbl), Intent(inout)   :: ufinal, biasfactor  
!!$
!!$    Integer, Parameter   :: MAX_PARAMS = 6   
!!$    Integer              :: natoms, i, j, k,index, nx, ny, nz
!!$    Integer              :: icellx,icelly,icellz, ncubelets
!!$    Integer              :: num,low,high,mid,flag 
!!$    Integer              :: atomt1, atomt2, atomt3, atomt4   
!!$    Integer, Dimension(ba_ntrials) :: index1, index2      
!!$    Real(kind=RDbl),Dimension(0:MAX_PARAMS-1) :: cn    
!!$    Real(kind=RDbl)      :: costheta,theta,phi 
!!$    Real(kind=RDbl)      :: ktheta1,thetaeq1,ktheta2,thetaeq2, &
!!$                            ktheta3,thetaeq3 
!!$    Real(kind=RDbl)      :: sqrdis,costac,costerm   
!!$    Real(kind=RDbl)      :: partition,key,wt,prob  
!!$    Real(kind=RDbl)      :: u,pot,energy_1
!!$    Real(kind=RDbl),Dimension(ba_ntrials*TA_NO_OF_TRIALS) :: &
!!$                    trial_e,ta_e,phi1,cum_prob,ba_num,ta_num 
!!$    Real(kind=RDbl),Dimension(ba_ntrials*TA_NO_OF_TRIALS) :: &
!!$                    theta1,ba_e1,theta2,ba_e2,theta3,ba_e3 
!!$    Type(VecType)        :: startpos
!!$    Logical              :: fast, ljflag
!!$    Type(MolecularParams), Pointer :: molecule
!!$
!!$    natoms = config_getnatoms(sorbates, sorbtype) 
!!$    !** Chooses the bond angle
!!$    !** Retrieves the molecule pointer
!!$    Call molecules_getpointer(sorbtype,molecule)
!!$    !** Obtains the potential parameters
!!$    atomt1 = ptr%P%atom_type
!!$    atomt2 = ptr%atom_type
!!$    atomt3 = ptr%C%atom_type 
!!$!LC    ktheta1= molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$!LC    thetaeq1= molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$    atomt1 = ptr%P%atom_type
!!$    atomt2 = ptr%atom_type
!!$    atomt3 = ptr%L%atom_type 
!!$!LC    ktheta2= molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$!LC    thetaeq2= molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$    atomt1 = ptr%C%atom_type
!!$    atomt2 = ptr%atom_type
!!$    atomt3 = ptr%L%atom_type 
!!$!LC    ktheta3= molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$!LC    thetaeq3= molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$
!!$    Do i=1,ba_ntrials
!!$      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/real(ba_ntrials)
!!$      theta1(i) = Acos(costheta)
!!$      ba_e1(i) = 0.5_RDbl*ktheta1*(theta1(i)-thetaeq1)*(theta1(i)-thetaeq1)
!!$    End Do 
!!$
!!$    Do i=1,ba_ntrials
!!$      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/real(ba_ntrials)
!!$      theta2(i) = Acos(costheta)
!!$      ba_e2(i) = 0.5_RDbl*ktheta2*(theta2(i)-thetaeq2)*(theta2(i)-thetaeq2)
!!$    End Do 
!!$
!!$    Do i=1,ba_ntrials
!!$      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/real(ba_ntrials)
!!$      theta3(i) = Acos(costheta)
!!$      ba_e3(i) = 0.5_RDbl*ktheta3*(theta3(i)-thetaeq3)*(theta3(i)-thetaeq3)
!!$    End Do 
!!$
!!$   !** Deals with the case where no dihedral potential 
!!$    If (.NOT. Associated(ptr%P%P)) Then 
!!$    partition = 0.0_RDbl 
!!$    num = 0 
!!$   !** locate both ptr%C and ptr%L  
!!$    Do i=1,ba_ntrials
!!$    Do j=1,ba_ntrials
!!$      u =  ba_e1(i) + ba_e2(j) 
!!$      ptr%baPC = theta1(i) 
!!$      ptr%baPL = theta2(j) 
!!$      ptr%taC  = 0.0       
!!$      ptr%taL  = pi - rranf()*twopi 
!!$   !** computes  ptr%baCL  
!!$      Call branchedcoords_toxyz  &
!!$           (sorbates(sorbtype)%gcoords(1)%branchedcoords, &
!!$           sorbates(sorbtype)%coords(1:natoms,1)%rp)
!!$    !** Since the node is a branch point, then additional bond angle to be  
!!$    !** sampled(baCL) 
!!$    !** Since the node is a branch point, then additional bond angle to be
!!$    !** sampled(baCL): dis is the distance between ptr%C and ptr%L
!!$    sqrdis = vector_getdistsq(sorbates(sorbtype)&
!!$               %coords(ptr%C%atom_num,1)%rp, &
!!$               sorbates(sorbtype)%coords&
!!$               (ptr%L%atom_num,1)%rp)
!!$    costheta = (ptr%C%bond_length*ptr%C%bond_length + &
!!$                 ptr%L%bond_length * ptr%L%bond_length - &
!!$                 sqrdis)/(2.0_RDbl*ptr%C%bond_length* ptr%L%bond_length)
!!$    theta = Acos(costheta)
!!$
!!$     If ( theta<= (2.0_RDbl*pi- ptr%baPC - ptr%baPL) &
!!$          .AND. theta   >= ABS(ptr%baPC - ptr%baPL) ) Then 
!!$     Else 
!!$        Write(0,'(2a,i4,a,f16.2,6f16.10)') __FILE__,": ",__LINE__, &
!!$        " angle out of range " 
!!$      stop
!!$     End If      
!!$
!!$    u = u + 0.5_RDbl*ktheta3*(theta-thetaeq3)* &
!!$        (theta-thetaeq3)  
!!$     If(u*rti < 30.0_RDbl) Then 
!!$      num = num + 1  
!!$      index1(num) = i 
!!$      index2(num) = j 
!!$      phi1(num) = ptr%taC 
!!$      trial_e(num) = u 
!!$      partition = exp(-(trial_e(num))*rti) + partition
!!$     End If  
!!$    End Do
!!$    End Do
!!$
!!$    wt = 0.0_RDbl
!!$    Do i=1,num  
!!$      wt = wt + exp(-trial_e(i)*rti)/partition
!!$      cum_prob(i) = wt
!!$    End Do
!!$
!!$    low = random_gettrial(1,num,cum_prob)
!!$    index = index1(low) 
!!$    costheta = ( rranf() * 2.0_RDbl/real(ba_ntrials) ) -  &
!!$        1.0_RDbl + (2.0_RDbl*(index-1))/Real(ba_ntrials)
!!$    theta = Acos(costheta)
!!$    u =  0.5_RDbl*ktheta1*(theta-thetaeq1)*(theta-thetaeq1)
!!$    ptr%P%baPC = theta 
!!$    index = index2(low) 
!!$    costheta = ( rranf() * 2.0_RDbl/Real(ba_ntrials) ) - &
!!$        1.0_RDbl + (2.0_RDbl*(index-1))/Real(ba_ntrials)
!!$    theta = Acos(costheta)
!!$    u = u + 0.5_RDbl*ktheta2*(theta-thetaeq2)*(theta-thetaeq2)
!!$    ptr%P%baPL = theta 
!!$    ptr%P%taC = phi1(low)  
!!$    Call branchedcoords_toxyz  &
!!$         (sorbates(sorbtype)%gcoords(1)%branchedcoords, &
!!$          sorbates(sorbtype)%coords(1:natoms,1)%rp)
!!$    sqrdis = vector_getdistsq(sorbates(sorbtype)& 
!!$             %coords(ptr%atom_num,1)%rp, &
!!$             sorbates(sorbtype)%coords& 
!!$             (ptr%P%C%atom_num,1)%rp)   
!!$    costheta = (ptr%bond_length*ptr%bond_length + & 
!!$               ptr%P%C%bond_length * ptr%P%C%bond_length - & 
!!$               sqrdis)/(2.0_RDbl*ptr%bond_length* ptr%P%C%bond_length)  
!!$    theta = Acos(costheta)
!!$    u = u + 0.5_RDbl*ktheta3*(theta-thetaeq3)*(theta-thetaeq3)
!!$    ufinal = ufinal + u  
!!$    biasfactor = biasfactor*real(ba_ntrials)*exp(-trial_e(low)*rti)& 
!!$                 /partition 
!!$    End If
!!$
!!$    If(Associated(ptr%C)) Then
!!$     Call brmoves_NodeIntegral(sorbates, sorbtype, ptr%P%C%C, &
!!$          rti, ufinal, biasfactor)
!!$    End If
!!$
!!$  End Subroutine brmoves_BrNdInteg  
!!$
!!$
!!$  !-----------------------------------------------------------
!!$  ! Calculate Z/omelga ratio from each individual note 
!!$  ! ufinal and biasfactor are returned 
!!$  !-----------------------------------------------------------
!!$  Recursive Subroutine brmoves_NodeIntegral(sorbates, sorbtype, ptr, &
!!$            rti, ufinal, biasfactor)
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Type(NodeType), Pointer       :: ptr                   
!!$    Integer, Intent(in)   :: sorbtype   
!!$    Real(kind=RDbl), Intent(in)    :: rti                   
!!$    Real(kind=RDbl), Intent(inout)   :: ufinal, biasfactor  
!!$
!!$    Integer, Parameter   :: MAX_PARAMS = 6   
!!$    Integer              :: natoms, i, j, k,index, nx, ny, nz
!!$    Integer              :: icellx,icelly,icellz, ncubelets
!!$    Integer              :: num,low,high,mid,flag 
!!$    Integer              :: atomt1, atomt2, atomt3, atomt4   
!!$    Real(kind=RDbl),Dimension(0:MAX_PARAMS-1) :: cn    
!!$    Real(kind=RDbl)      :: costheta,ktheta,thetaeq,theta,phi 
!!$    Real(kind=RDbl)      :: costac,costerm   
!!$    Real(kind=RDbl)      :: partition,key,wt,prob  
!!$    Real(kind=RDbl)      :: u,pot,energy_1
!!$    Real(kind=RDbl),Dimension(ba_ntrials*TA_NO_OF_TRIALS) :: &
!!$                    trial_e,ta_e,ba_e,theta1,phi1,cum_prob  
!!$    Integer,Dimension(ba_ntrials*TA_NO_OF_TRIALS):: &
!!$                    ba_num, ta_num                 
!!$    Type(VecType)        :: startpos
!!$    Logical              :: fast, ljflag
!!$    Type(MolecularParams), Pointer :: molecule
!!$    Real(kind=RDbl), Dimension(NO_OF_INTRA_POTS) :: potlist
!!$
!!$    natoms = config_getnatoms(sorbates, sorbtype) 
!!$    !** Chooses the bond angle
!!$    !** Retrieves the molecule pointer
!!$    Call molecules_getpointer(sorbtype,molecule)
!!$    !** Obtains the potential parameters
!!$    atomt1 = ptr%P%P%atom_type
!!$    atomt2 = ptr%P%atom_type
!!$    atomt3 = ptr%atom_type 
!!$!LC    ktheta = molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$!LC    thetaeq = &
!!$!LC        molecule%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$
!!$   !** Obtains the dihedral angle potential
!!$    atomt1 = ptr%P%P%P%atom_type
!!$    atomt2 = ptr%P%P%atom_type
!!$    atomt3 = ptr%P%atom_type
!!$    atomt4 = ptr%atom_type     
!!$    Do i=0,MAX_PARAMS-1
!!$!LC      cn(i) = &
!!$!LC          molecule%torsion%torparams(atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)
!!$    End Do
!!$
!!$    partition = 0.0_RDbl
!!$    Do i=1,ba_ntrials
!!$      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/real(ba_ntrials)
!!$      theta1(i) = Acos(costheta)
!!$      ba_e(i) = 0.5_RDbl*ktheta*(theta1(i)-thetaeq)*(theta1(i)-thetaeq)
!!$      phi1(i) = pi - rranf()*twopi       
!!$      costac  = cos(phi1(i))
!!$      costerm= 1.0_RDbl
!!$      ta_e(i) = 0.0_RDbl 
!!$        Do j=0,MAX_PARAMS-1
!!$          ta_e(i)= ta_e(i) + cn(j)*costerm
!!$          costerm=costerm*costac
!!$        End Do
!!$    trial_e(i) = ba_e(i) + ta_e(i) 
!!$    partition = exp(-(trial_e(i))*rti) + partition
!!$    End Do
!!$
!!$    wt = 0.0_RDbl
!!$    Do i=1,ba_ntrials   
!!$      wt = wt + exp(-trial_e(i)*rti)/partition
!!$      cum_prob(i) = wt
!!$    End Do
!!$
!!$    low = random_gettrial(1,ba_ntrials,cum_prob)
!!$    costheta = ( rranf() * 2.0_RDbl/real(ba_ntrials) ) - &
!!$               1.0_RDbl + (2.0_RDbl*(low-1))/real(ba_ntrials)
!!$    theta = Acos(costheta)
!!$    phi = phi1(low)   
!!$    ufinal = ufinal + 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$    ufinal = ufinal + ta_e(low)  
!!$
!!$    ptr%P%baPC = theta     
!!$    ptr%P%taC  = phi    
!!$    Call branchedcoords_toxyz  &
!!$         (sorbates(sorbtype)%gcoords(1)%branchedcoords, &
!!$         sorbates(sorbtype)%coords(1:natoms,1)%rp)
!!$    ljflag = .True.   
!!$    Call ssbasic_IntraGetinteraction(sorbates(sorbtype),&
!!$         ssnoncoulparams(sorbtype,sorbtype)%basic,molecule,1,&
!!$         ptr%atom_num,pot,ljflag)
!!$    ufinal = ufinal + pot
!!$    !** Fast flag is hacked
!!$    fast = .True.
!!$    Call intramolecular_getinteraction(sorbtype,sorbates(sorbtype),fast, &
!!$        potlist)
!!$    ufinal = ufinal + potlist(INTRAPAIR_INDEX) 
!!$    biasfactor = biasfactor*real(ba_ntrials)*exp(-trial_e(low)*rti)& 
!!$                 /partition 
!!$
!!$    If(Associated(ptr%P%L)) Then
!!$     Call brmoves_BrNdInteg(sorbates, sorbtype, &
!!$          ptr%P%L, rti, ufinal, biasfactor)
!!$    Else If(Associated(ptr%C)) Then
!!$     Call brmoves_NodeIntegral(sorbates, sorbtype, ptr%C, &
!!$          rti, ufinal, biasfactor)
!!$    End If
!!$
!!$  End Subroutine brmoves_NodeIntegral  


!!$!*** I Have not checked it so comment it out
!!$  !------------------------------------------------------------------
!!$  ! Inserts one fragment of the molecule and returns the bias of
!!$  ! and energy of the insertion
!!$  ! Note: both the bond angel and dihedral potentials have to exist
!!$  ! Note: ptr is already fixed, ptr%C and ptr%L are to be fixed    
!!$  !------------------------------------------------------------------
!!$  Recursive Subroutine brmoves_insertframe(sorbates,nsorbs,sorbtype,&
!!$            molec,ptr,params, biasfactor,ufinal,in_flag) 
!!$
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Integer, Intent(in)   :: nsorbs,sorbtype,molec
!!$    Type(BranchedBIns_Params), Intent(in) :: params
!!$    Real(kind=RDbl), Intent(inout)   :: biasfactor, ufinal
!!$    Logical, Intent(out)   :: in_flag
!!$    Type(NodeType), Pointer       :: ptr
!!$
!!$    Integer, Parameter   :: MAX_PARAMS = 6
!!$    Integer              :: natoms, num, low, index, flag, i, j
!!$    Integer              :: atomt1, atomt2, atomt3, atomt4
!!$    Integer              :: found  
!!$    Integer,Dimension(sorbates(sorbtype)%natoms) :: atom_num   
!!$    Real(kind=RDbl)      :: costheta,ktheta,theta,thetaeq 
!!$    Real(kind=RDbl)      :: costac,costerm
!!$    Real(kind=RDbl)      :: partition,wt,rti,bias 
!!$    Real(kind=RDbl)      :: u,noncoulpot,coulpot
!!$    Real(kind=RDbl)      :: dis1,dis2,dis3 
!!$    Real(kind=RDbl),Dimension(TA_NO_OF_TRIALS*2) :: trial_e,theta1,cum_prob
!!$    Real(kind=RDbl),Dimension(TA_NO_OF_TRIALS*2) :: taC,taL      
!!$    Real(kind=RDbl),Dimension(0:MAX_PARAMS-1) :: cn1, cn2 
!!$    Logical              :: mapflag, fast        
!!$
!!$    Type(MolecularParams), Pointer :: molecule
!!$    Type(Frame_Library_Params), Pointer :: frame_lib 
!!$
!!$    in_flag = .True. 
!!$    Call branchedcoords_placenode(ptr)
!!$    !** Set the temperature parameter
!!$    rti = params%rti
!!$    !** Obtains the frame library for this node  
!!$    Call brmoves_find_Insframelib(params,ptr%atom_num,frame_lib,found) 
!!$    If (found == 0) Then 
!!$      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
!!$       " Could not find the frame_lib "        
!!$      Stop
!!$    End If
!!$
!!$    !** Chooses the two dihedral angles   
!!$    Call molecules_getpointer(sorbtype,molecule)
!!$    !** Obtains the dihedral angle potential
!!$!    If (Associated(ptr%P%P)) Then 
!!$!    atomt1 = ptr%P%P%atom_type
!!$!    atomt2 = ptr%P%atom_type
!!$!    atomt3 = ptr%atom_type
!!$!    atomt4 = ptr%C%atom_type
!!$!    Do i=0,MAX_PARAMS-1
!!$!LC       cn1(i)= molecule%torsion%torparams & 
!!$!LC             (atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)  
!!$!    End Do
!!$!    atomt1 = ptr%P%P%atom_type
!!$!    atomt2 = ptr%P%atom_type
!!$!    atomt3 = ptr%atom_type
!!$!    atomt4 = ptr%L%atom_type
!!$!    Do i=0,MAX_PARAMS-1
!!$!LC       cn2(i)= molecule%torsion%torparams & 
!!$!LC             (atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)  
!!$!    End Do
!!$!    End If 
!!$
!!$    !** Two-step sampling:                               
!!$    !** Step one: selects the frame so one dihedral angle 
!!$    !** is determined, 
!!$    !** Get the values of bond angles using the frame library
!!$
!!$    index = INT(rranf()*frame_lib%num_frame) + 1 
!!$    biasfactor = biasfactor * frame_lib%frame(index)%biasfactor
!!$    ptr%baPC= frame_lib%frame(index)%baPC
!!$    ptr%baPL= frame_lib%frame(index)%baPL
!!$    ptr%baCL= frame_lib%frame(index)%baCL
!!$     
!!$    ufinal = ufinal + frame_lib%frame(index)%u  
!!$    ptr%num  = index  
!!$
!!$    !** Step two: the other is then chosen randomly 
!!$    !** Note: taL cannot be decided without taC  
!!$    partition = 0.0_RDbl
!!$          num = 0
!!$    Do i=1,TA_NO_OF_TRIALS
!!$     !**build up the node
!!$     !Samples taC  
!!$      ptr%taC =  pi - rranf()*twopi
!!$     ! computes taL
!!$      Call branchedcoords_update_ta(ptr)
!!$
!!$      u  = 0.0_RDbl
!!$      If (Associated(ptr%P%P)) Then 
!!$        costac   = cos(ptr%taC)
!!$        costerm= 1.0_RDbl
!!$        Do j=0,MAX_PARAMS-1
!!$          u= u + cn1(j)*costerm
!!$          costerm=costerm*costac
!!$        End Do
!!$        costac   = cos(ptr%taL)
!!$        costerm= 1.0_RDbl
!!$        Do j=0,MAX_PARAMS-1
!!$          u= u + cn2(j)*costerm
!!$          costerm=costerm*costac
!!$        End Do
!!$      End If 
!!$
!!$      natoms = ptr%L%atom_num
!!$      Call branchedcoords_toxyz  &
!!$      (sorbates(sorbtype)%gcoords(molec)%branchedcoords, &
!!$        sorbates(sorbtype)%coords(1:natoms,molec)%rp)
!!$      Call simcell_pbc(params%scell, sorbates(sorbtype)%coords(1:natoms,molec)%rp, &
!!$         sorbates(sorbtype)%coords(1:natoms, molec)%r, &
!!$         sorbates(sorbtype)%coords(1:natoms, molec)%cr)
!!$
!!$    ! Calculate the energy
!!$      fast = .true.
!!$      flag = 0
!!$
!!$      atom_num(1)=ptr%C%atom_num 
!!$      atom_num(2)=ptr%L%atom_num 
!!$    
!!$        Do j=1, nsorbs
!!$         mapflag = .True.
!!$         noncoulpot = 0.0_RDbl
!!$         coulpot = 0.0_RDbl 
!!$         mapflag = ssnoncoul_asdriver(sorbates, params%scell, sorbtype, & 
!!$                  molec, atom_num(1:2), j, fast, noncoulpot, & 
!!$                  sorbates(sorbtype)%afast(1:natoms, molec), &
!!$                  sorbates(j)%afast)
!!$
!!$         u = u + noncoulpot
!!$
!!$         If (.NOT. (mapflag) )  Then
!!$              flag = 1
!!$              Exit
!!$         End If
!!$
!!$         mapflag = sscoul_msdriver(sorbates, params%scell, sorbtype, & 
!!$                  molec, j, fast, coulpot, & 
!!$                  sorbates(sorbtype)%afast(1:natoms, molec), &
!!$                  sorbates(j)%afast)
!!$
!!$         u = u + coulpot
!!$         If (.NOT. (mapflag) )  Then
!!$              flag = 1
!!$              Exit
!!$         End If
!!$
!!$        End Do
!!$
!!$        If (flag == 0 .AND. u*rti < 30.0_RDbl ) Then
!!$!       If (flag == 0  ) Then
!!$         num = num + 1
!!$         taC(num) =  ptr%taC          
!!$         taL(num) =  ptr%taL         
!!$         trial_e(num) = u
!!$         If (trial_e(num) < -30)  trial_e(num)= -30
!!$         partition = exp(-trial_e(num)*rti) + partition
!!$        End If
!!$
!!$      End Do
!!$
!!$   !nflag measures the success of the energy calculation
!!$     If (num == 0) Then
!!$         in_flag =  .False. 
!!$         Return
!!$     End If
!!$
!!$     wt = 0.0_RDbl
!!$     Do i=1,num
!!$       wt = wt + exp(-trial_e(i)*rti)/partition
!!$       cum_prob(i) = wt
!!$     End Do
!!$
!!$    low = random_gettrial(1,num,cum_prob)
!!$    biasfactor = biasfactor * &
!!$                 partition/(real(TA_NO_OF_TRIALS)*exp(-trial_e(low)*rti))
!!$    ufinal = ufinal + trial_e(low)  
!!$    ptr%taC = taC(low)
!!$    ptr%taL = taL(low)
!!$
!!$    Call branchedcoords_toxyz  &
!!$       (sorbates(sorbtype)%gcoords(molec)%branchedcoords, &
!!$        sorbates(sorbtype)%coords(1:natoms,molec)%rp)
!!$    Call simcell_pbc(params%scell, sorbates(sorbtype)%coords(1:natoms,molec)%rp, &
!!$        sorbates(sorbtype)%coords(1:natoms, molec)%r, &
!!$        sorbates(sorbtype)%coords(1:natoms, molec)%cr)
!!$
!!$    If (Associated(ptr%C)) Then   
!!$     If (Associated(ptr%C%L)) Then   
!!$       Call  brmoves_insertframe(sorbates,nsorbs,sorbtype,molec, &
!!$             ptr%C,params,biasfactor,ufinal,in_flag)
!!$     End If  
!!$    Else 
!!$     Call  brmoves_insertnode(sorbates,nsorbs,sorbtype,molec, &
!!$           ptr%C,params,biasfactor,ufinal,in_flag) 
!!$    End If
!!$    If (Associated(ptr%L)) Then   
!!$     If (Associated(ptr%L%L)) Then   
!!$       Call  brmoves_insertframe(sorbates,nsorbs,sorbtype,molec, &
!!$             ptr%C,params,biasfactor,ufinal,in_flag)
!!$     End If  
!!$    Else 
!!$     Call  brmoves_insertnode(sorbates,nsorbs,sorbtype,molec, &
!!$           ptr%L,params,biasfactor,ufinal,in_flag) 
!!$    End If
!!$
!!$    in_flag = .True.
!!$ 
!!$    Return
!!$  End Subroutine brmoves_insertframe 
!!$
!!$  !------------------------------------------------------------------
!!$  ! Deletes one fragment of the molecule and returns the bias of
!!$  ! and energy of the insertion
!!$  ! Note: both the bond angel and dihedral potentials have to exist
!!$  !------------------------------------------------------------------
!!$  Recursive Subroutine brmoves_deleteframe(sorbates,nsorbs,sorbtype,&
!!$            molec,ptr,params,biasfactor,ufinal)
!!$
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Integer, Intent(in)   :: nsorbs,sorbtype,molec
!!$    Type(BranchedBDel_Params), Intent(in)   :: params
!!$    Real(kind=RDbl), Intent(inout)   :: biasfactor, ufinal
!!$    Type(NodeType), Pointer       :: ptr
!!$
!!$    Integer, Parameter   :: MAX_PARAMS = 6
!!$    Integer              :: natoms, num, index, low, flag, i, j
!!$    Integer              :: atomt1, atomt2, atomt3, atomt4
!!$    Integer              :: found  
!!$    Integer, Dimension(sorbates(sorbtype)%natoms)   :: atom_num  
!!$    Real(kind=RDbl)      :: orig_e,orig_baPC,orig_baPL,orig_taC,orig_taL  
!!$    Real(kind=RDbl)      :: costheta,ktheta,theta,thetaeq 
!!$    Real(kind=RDbl)      :: costac,costerm
!!$    Real(kind=RDbl)      :: partition,wt,rti,bias 
!!$    Real(kind=RDbl)      :: u,noncoulpot,coulpot
!!$    Real(kind=RDbl),Dimension(TA_NO_OF_TRIALS*2) :: trial_e,theta1,cum_prob
!!$    Real(kind=RDbl),Dimension(TA_NO_OF_TRIALS*2) :: baPC,baPL,taC,taL  
!!$    Real(kind=RDbl),Dimension(0:MAX_PARAMS-1) :: cn1, cn2 
!!$    Logical              :: mapflag, fast        
!!$
!!$    Type(MolecularParams), Pointer :: molecule
!!$    Type(Frame_Library_Params), Pointer :: frame_lib  
!!$
!!$    Call branchedcoords_placenode(ptr)
!!$    !** Set the temperature parameter
!!$    rti = params%rti
!!$    !** Obtains the frame library for this node
!!$    Call brmoves_find_DelFramelib(params,ptr%atom_num,frame_lib,found)
!!$    If (found == 0) Then
!!$      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
!!$       " Could not find the frame_lib "
!!$      Stop
!!$    End If
!!$
!!$    natoms = config_getnatoms(sorbates, sorbtype) 
!!$    !** Chooses the two dihedral angles   
!!$    Call molecules_getpointer(sorbtype,molecule)
!!$    !** Obtains the dihedral angle potential
!!$    If (Associated(ptr%P%P)) Then 
!!$    atomt1 = ptr%P%P%atom_type
!!$    atomt2 = ptr%P%atom_type
!!$    atomt3 = ptr%atom_type
!!$    atomt4 = ptr%C%atom_type
!!$    Do i=0,MAX_PARAMS-1
!!$!LC       cn1(i)= molecule%torsion%torparams & 
!!$!LC             (atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)  
!!$    End Do
!!$
!!$    atomt1 = ptr%P%P%atom_type
!!$    atomt2 = ptr%P%atom_type
!!$    atomt3 = ptr%atom_type
!!$    atomt4 = ptr%L%atom_type
!!$    Do i=0,MAX_PARAMS-1
!!$!LC       cn2(i)= molecule%torsion%torparams & 
!!$!LC             (atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)  
!!$    End Do
!!$    End If 
!!$
!!$    !** Adds the hard degrees of potential   
!!$    ufinal = ufinal + frame_lib%frame(ptr%num)%u    
!!$    
!!$    !** Chooses one dihedral angle: no need to choose a frame 
!!$    !** from one dihedral angle and three known bond angles, 
!!$    !** the other dihedral is determined 
!!$    partition = 0.0_RDbl
!!$    Do i=1,TA_NO_OF_TRIALS
!!$
!!$     If ( i == 1 ) Then 
!!$      orig_taC =  ptr%taC          
!!$      orig_taL =  ptr%taL          
!!$     Else 
!!$     !**build up the node
!!$        ptr%taC =  pi - rranf()*twopi
!!$       ! computes taL 
!!$        Call branchedcoords_update_ta(ptr) 
!!$        Call branchedcoords_toxyz  &
!!$            (sorbates(sorbtype)%gcoords(molec)%branchedcoords, &
!!$             sorbates(sorbtype)%coords(1:natoms,molec)%rp)
!!$        Call simcell_pbc(params%scell, & 
!!$             sorbates(sorbtype)%coords(1:natoms,molec)%rp, &
!!$             sorbates(sorbtype)%coords(1:natoms, molec)%r, &
!!$             sorbates(sorbtype)%coords(1:natoms, molec)%cr)
!!$     End If
!!$
!!$      !**u does not include the bonded potential  
!!$      !**since the bonded energy does not contribute to 
!!$      !**the biasfactor of the soft degrees of freedom  
!!$      u = 0.0_RDbl        
!!$      If (Associated(ptr%P%P)) Then 
!!$        costac   = cos(ptr%taC)
!!$        costerm= 1.0_RDbl
!!$        Do j=0,MAX_PARAMS-1
!!$          u= u + cn1(j)*costerm
!!$          costerm=costerm*costac
!!$        End Do
!!$        costac   = cos(ptr%taL)
!!$        costerm= 1.0_RDbl
!!$        Do j=0,MAX_PARAMS-1
!!$            u= u + cn2(j)*costerm
!!$          costerm=costerm*costac
!!$        End Do
!!$      End If 
!!$
!!$    ! Calculate the energy
!!$      fast = .true.
!!$      flag = 0
!!$
!!$      atom_num(1)=ptr%C%atom_num 
!!$      atom_num(2)=ptr%L%atom_num 
!!$
!!$        Do j=1, nsorbs
!!$         mapflag = .True.
!!$         noncoulpot = 0.0_RDbl
!!$         coulpot = 0.0_RDbl 
!!$        mapflag = ssnoncoul_asdriver(sorbates, params%scell, sorbtype, & 
!!$                  molec, atom_num(1:2), j, fast, noncoulpot, & 
!!$                  sorbates(sorbtype)%afast(1:natoms, molec), &
!!$                  sorbates(j)%afast)
!!$         u = u + noncoulpot
!!$         If (.NOT. (mapflag) )  Then
!!$              flag = 1
!!$              Exit
!!$         End If
!!$
!!$        mapflag = sscoul_msdriver(sorbates, params%scell, sorbtype, & 
!!$                  molec, j, fast, coulpot, & 
!!$                  sorbates(sorbtype)%afast(1:natoms, molec), &
!!$                  sorbates(j)%afast)
!!$         u = u + coulpot
!!$         If (.NOT. (mapflag) )  Then
!!$              flag = 1
!!$              Exit
!!$         End If
!!$
!!$        End Do
!!$
!!$        If ( i == 1 ) Then 
!!$            orig_e = u
!!$            ufinal = ufinal + u
!!$            partition = exp(-u*rti) 
!!$        Else If (flag == 0 .AND. u*rti < 30.0_RDbl ) Then
!!$         partition = exp(-u*rti) + partition
!!$        End If
!!$
!!$      End Do
!!$
!!$    biasfactor = (biasfactor/frame_getbiasfactor  &
!!$                    (frame_lib,ptr%num ) ) *  &  
!!$                 partition/(real(TA_NO_OF_TRIALS)*exp(-orig_e*rti))   
!!$    ptr%taC  =  orig_taC 
!!$    ptr%taL  =  orig_taL 
!!$
!!$    Call branchedcoords_toxyz  &
!!$       (sorbates(sorbtype)%gcoords(molec)%branchedcoords, &
!!$        sorbates(sorbtype)%coords(1:natoms,molec)%rp)
!!$    Call simcell_pbc(params%scell, &
!!$        sorbates(sorbtype)%coords(1:natoms,molec)%rp, &
!!$        sorbates(sorbtype)%coords(1:natoms, molec)%r, &
!!$        sorbates(sorbtype)%coords(1:natoms, molec)%cr)
!!$
!!$    If (Associated(ptr%C)) Then
!!$     If (Associated(ptr%C%L)) Then
!!$       Call  brmoves_deleteframe(sorbates,nsorbs,sorbtype,molec,&
!!$             ptr%C,params,biasfactor,ufinal)
!!$     End If
!!$    Else
!!$     Call  brmoves_deletenode(sorbates,nsorbs,sorbtype,molec,&
!!$           ptr%L,params,biasfactor,ufinal)
!!$    End If
!!$    If (Associated(ptr%L)) Then
!!$     If (Associated(ptr%L%L)) Then
!!$       Call  brmoves_deleteframe(sorbates,nsorbs,sorbtype,molec, &
!!$             ptr%C,params,biasfactor,ufinal)
!!$     End If
!!$    Else
!!$     Call  brmoves_deletenode(sorbates,nsorbs,sorbtype,molec, &
!!$           ptr%L,params,biasfactor,ufinal)
!!$    End If
!!$
!!$    Return
!!$  End Subroutine brmoves_deleteframe 
!---------------------------------------------------------------
! ABOVE ARE ROUTINES WRITTEN BY YI FOR DEALING WITH BRANCHES IN 
! ALKANES, THEY ARE NOT TESTED, AND WONT COMPILE WITH THE CURRENT CODE.
! STORED HERE FOR FUTURE REFERNCE. DONT DELETE -Shaji
!---------------------------------------------------------------



