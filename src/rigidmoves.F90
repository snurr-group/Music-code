!-------------------------------------------------------------------------
! This module implements the different kinds of moves that can be done
! on rigid molecules.  The move routines for the move types that are used 
! in the Monte Carlo routines, should all return the biasfactor.
!
! Current rigid move types are:
!    Random translation - move the molecule a bounded random distance
!    Random rotation - rotate the molecule by bounded random increment
!    Biased rotation/translation - Rosenbluth-style multi-config biasing
!    Reinsertion - reinsert the molecule with new rotation 
!    random insertion - insertion of randomly generated molecule
!    biased insertion - energy-biased insertion of random molecule
!    library insertion - insertion of molecule from a library
!    random deletion - deletion of molecule
!    biased deletion - energy-biased deletion of molecule
!
! Needed Improvements: 
! 1) finish "requires" comments for the routines
! 2) remove simcell pointer from the move data structures
! There needs to be sweeping changes in this module.  An examination of
! the move type data structures indicates that they are very similar
! and could be consolidated.  I would also suggest a separation of 
! biasing from the perturbations themselves.  
!-------------------------------------------------------------------------

Module rigidmoves
  Use auxmoveparams, Only: AuxMoveObjects
  Use bmap, Only: BiasMap_Params,bmap_init,bmap_getbiasindex,&
      bmap_getncubelets,bmap_getbiaswt,bmap_getbiaspt,bmap_getcellwt, &
      bmap_getcellCOM,bmap_display
  Use cavitylist, Only: Cavity_Params, cavitylist_checkcavity, &
      cavitylist_updatemolecinfo, cavitylist_checkandIncrMemory
  Use config, Only: AtMolCoords, config_getTotMoles, config_dumpmol, &
      config_config2xyz, config_rotatexyz, config_getMolecCOM, &
      config_rotateAroundCom, config_copymolec, config_allocfields, &
      config_subset2xyz
  Use conflib, Only: ConfLibrary, conflib_basicInit, conflib_dataInit, &
      conflib_getrandomconfig
  Use defaults, Only: RDbl, strLen, lstrLen, pi,twopi, zero, one, &
      scalepe, kcalmole_kb, STATS_BLOCKSIZE, DEFAULT_MAX_EXP
  Use file, Only: file_open
  Use general, Only: genparams
  Use subinteract, Only: Subset_Interactions,subinteract_int, &
      subinteract_oldnrg,subinteract_newnrg
  Use matrix, Only: MatrixType,matrix_display,Operator(*),Assignment(=)
  Use molecules, Only: molecules_getnatoms, molecules_getcomdefn, &
      molecules_isflexible, molecules_name
  Use random, Only: rranf,random_getnewseed
  Use rigidcoords, Only: rigidcoords_toxyz, rigidcoords_changeref, &
      rigidcoords_setfrom_rp, RigidMolec, rigidcoords_getgencoords, &
      rigidcoords_setgencoords
  Use simcell, Only: SimCell_Params,simcell_getzeoell,simcell_geteff,simcell_maptouc,&
      simcell_getcellorigin, simcell_uniformRandPt, simcell_pbc
  Use stats, Only: Statistics, stats_init, stats_update, &
      stats_getcavg, Assignment(=)
  Use storebase, Only: storebase_display
  Use storetop, Only: storetop_display
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, real2str, &
      toint, tolower, allocErrDisplay, checkandstop, int2str, cleanstring
  Use vector, Only: VecType, Assignment(=), Operator(+), &
      Operator(-), Operator(*), Operator(/), vector_filedisplay

  Implicit None
  Save

  Private
  Public :: RigidRTRans_Params,RigidRRot_Params,RigidReinsert_Params, &
      RigidRIns_Params,RigidLIns_Params,RigidBIns_Params,&
      RigidBDel_Params,RigidRDel_Params,RigidBiasRT_Params, &
      rigidmoves_RTransinit,rigidmoves_RRotInit, &
      rigidmoves_ReinsertInit,rigidmoves_RInsinit,rigidmoves_BInsinit, &
      rigidmoves_RDelInit, rigidmoves_BDelInit,&
      rigidmoves_LInsInit, rigidmoves_biasrtInit, rigidmoves_randomtranslate, &
      rigidmoves_randomrotate,rigidmoves_reinsert, &
      rigidmoves_rinsert,rigidmoves_binsert,rigidmoves_linsert, &
      rigidmoves_rdelete, rigidmoves_bdelete, rigidmoves_biasrt, &
      rigidmoves_rtretry, &
      rigidmoves_adjustdeltatrans,rigidmoves_adjustdeltarot, &
      rigidmoves_displayparams, rigidmoves_adjustBiasRT, &
      rigidmoves_dyninfo, rigidmoves_checkOppMoves, &
      rigidmoves_postadjust, rigidInsBias

  !** Rigid gc random length translation params
  Type RigidRTrans_Params
    Logical           :: constrained
    Type(MatrixType)  :: constrtrans ! Transforming to the consstr system.
    Type(MatrixType)  :: constrinvtrans
    Type(VecType)     :: normal
    Logical           :: scale_trans
    Real(kind=RDbl)   :: deltatrans
    Type(SimCell_Params), Pointer :: simcell
    Type(Cavity_Params), Pointer :: cavity
    Real(kind=Rdbl)   :: max_nrg
  End Type RigidRTrans_Params

  !** Rigid gc random length Rotation params
  Type RigidRRot_Params
    Logical           :: scale_rot
    Real(kind=RDbl)   :: deltarot
    Type(SimCell_Params), Pointer :: simcell
    Type(Cavity_Params), Pointer :: cavity
    Real(kind=Rdbl)   :: max_nrg
  End Type RigidRRot_Params

  !** Rigid gc reinsertion parameters
  Type RigidReinsert_Params
    Real(kind=Rdbl)   :: max_nrg
  End Type RigidReinsert_Params

  !** Rigid Biased Rotation/Translation parameters
  Type RigidBiasRT_Params
    Integer           :: ntrials
    Logical           :: scale_rot,scale_trans
    Real(kind=RDbl)   :: deltarot,deltatrans
    Real(kind=Rdbl)   :: max_nrg,biastemp,beta
  End Type RigidBiasRT_Params

  !** Rigid gc random Insertion params
  Type RigidRIns_Params
    Real(kind=Rdbl)   :: max_nrg
    Type(SimCell_Params), Pointer :: simcell
  End Type RigidRIns_Params

  !** Rigid gc biased Insertion params
  Type RigidBIns_Params
    Character(len=strLen)         :: type_of_bias ! ** MAP_BIAS or CAVITY_BIAS
    !       or MAP_CAVITY_BIAS
    Type(RigidInsBias),Pointer    :: insbias
    Real(kind=RDbl)               :: max_nrg      
    Type(SimCell_Params), Pointer :: simcell
  End Type RigidBIns_Params

  !** For biased insertion from a library of configurations 
  !** there are 3   biases : volume bias from bmap
  !**                      : config bias from the conflib-library
  !**                      : bias from cavity bias
  Type RigidLIns_Params
    Character(len=strLen)         :: type_of_bias !** MAP_BIAS or CAVITY_BIAS
                                                  !** or MAP_CAVITY_BIAS
    Type(RigidInsBias),Pointer    :: insbias
    Real(kind=RDbl)               :: max_nrg      
    Type(SimCell_Params), Pointer :: simcell
    Type(ConfLibrary), Pointer    :: lib
    Character(len=strLen)         :: confFileName
    Type(AtMolCoords),Pointer     :: temp_ref_struc 
  End Type RigidLIns_Params

  !** Rigid gc random deletion params
  Type RigidRDel_Params
    Real(kind=RDbl)               :: max_nrg      
    Type(SimCell_Params), Pointer :: simcell
  End Type RigidRDel_Params

  !** Rigid gc biased deletion params
  !** opposite move for BIns_Params, OR LIns_Params
  Type RigidBDel_Params
    Character(len=strLen)         :: type_of_bias !** MAP_BIAS or CAVITY_BIAS
                                                  !** or MAP_CAVITY_BIAS
    Type(RigidInsBias),Pointer    :: insbias
    Real(kind=RDbl)               :: max_nrg
    Type(SimCell_Params), Pointer :: simcell
  End Type RigidBDel_Params

  !** This is the bias for insertion used by Binsert, and Linsert 
  !** this Type will be used by _binsert, _bdelete, _linsert movetypes
  Type RigidInsBias
    Character(len=strLen)         :: type_of_bias !** MAP_BIAS or CAVITY_BIAS
                                                  !** or MAP_CAVITY_BIAS
    Real(kind=Rdbl)               :: max_nrg
    Type(SimCell_Params), Pointer :: simcell
    Type(BiasMap_Params), Pointer :: bmap
    Character(len=strLen)         :: biasfilename
    Real(kind=Rdbl)               :: bias_temp

    !** cav_prob = running average of number of cavities found in cavNTrials
    Type(Cavity_Params), Pointer  :: cavity
    Type(Statistics) , Dimension(:), Pointer :: cav_prob
    Integer                       :: cavNTrials
  End Type RigidInsBias

  Interface rigidmoves_displayparams
    Module Procedure rigidmoves_displayRTrans
    Module Procedure rigidmoves_displayRRot
    Module Procedure rigidmoves_displayReinsert
    Module Procedure rigidmoves_displayBiasRT
    Module Procedure rigidmoves_displayBIns
    Module Procedure rigidmoves_displayBDel
    Module Procedure rigidmoves_displayLins
  End Interface

  Interface rigidmoves_checkOppMoves
    Module procedure rigidmoves_checkBInsDel
    Module procedure rigidmoves_checkLInsDel
  End Interface

  Interface rigidmoves_postadjust
    Module Procedure rigidmoves_postadjustRot
    Module Procedure rigidmoves_postadjustTrans
    Module Procedure rigidmoves_postadjustBIns
    Module Procedure rigidmoves_postadjustLIns
    Module Procedure rigidmoves_postadjustDel
  End Interface

  Interface rigidmoves_dyninfo
    Module Procedure rigidmoves_dyninfoRTrans
    Module Procedure rigidmoves_dyninfoRRot
    Module Procedure rigidmoves_dyninfoBiasRT
  End Interface

  !** In cavitybias we have an array which keeps track of p(n)=probability
  !** of finding a cavity at n molecules. This is the defaults dimension 
  !** of that (used only for initialization)
  Integer, Parameter :: CAVBIAS_MAXMOLEC = 200

  !** These probably should be moved to somewhere else
  Real(kind=RDbl), Parameter        :: MIN_TRANS_DELTA = 0.01_RDbl
  Real(kind=RDbl), Parameter        :: MAX_TRANS_DELTA = 2.0_RDbl

  !** This too, for setting the maximum during moves
  Real(kind=RDbl), Parameter        :: VERY_LRG_NRG = 1.0e10_RDbl


#ifdef BENT_GCMC
! For molecules with a  bent geometry. For now it works with the hack 
! that atom-1  is the anchor-atom, instead of COM. This part should be 
! integrated to regular code after more testing. (Tested for benzene)

  Integer:: anchor_atom=1

#endif

Contains    

  !----------------------------------------------------------------------
  ! Read the information that will be required to do rigid RANDOM length
  ! TRANSLATION type of moves.
  ! Requires: params -- parameters to be initialized
  !           simcell -- simulation cell information
  !           auxmv -- contains info about cavity list
  !           spc -- species number
  !           filename -- file where information is located
  !----------------------------------------------------------------------
  Subroutine rigidmoves_RTransInit(params, simcell, auxmv, spc, filename)
    Type(RigidRTrans_Params), Intent(inout) :: params
    Integer, Intent(in)                     :: spc
    Character(*), Intent(in)                :: filename
    Type(SimCell_Params), Target, Intent(in)  :: simcell
    Type(AuxMoveObjects),Pointer           :: auxmv 

    Integer    :: unitno, scale_trans

    !** Initialize the cavity bais pointer
    Nullify(params%cavity)
    If (auxmv%cavBiasON)     params%cavity=>auxmv%cavity

    !** Set the constrained flag to .False.
    params%constrained = .False.

    !** Set the simcell pointer
    params%simcell => simcell

    params%max_nrg = VERY_LRG_NRG

    !** Get the delta displacements and "whether to scale or not".
    unitno = file_open(filename)
    Read(unitno,*) params%deltatrans, scale_trans

    !** Set the flag to decide whether we scale the jump lengths or not
    params%scale_trans = .False.
    If (scale_trans == 1) Then
      params%scale_trans = .True.
    End If

  End Subroutine rigidmoves_RTransInit

  !-----------------------------------------------------------------------
  ! Read the information that will be required to do RANDOM angle ROTATION
  ! type of moves.
  ! Requires: params -- parameters to be initialized
  !           simcell -- simulation cell information
  !           spc -- species number
  !           filename -- file where information is located
  !-----------------------------------------------------------------------
  Subroutine rigidmoves_RRotInit(params, simcell ,auxmv, spc, filename)
    Type(RigidRRot_Params), Intent(inout) :: params
    Type(SimCell_Params), Target, Intent(In)  :: simcell
    Integer, Intent(in)        :: spc
    Character(*), Intent(in)   :: filename
    Type(AuxMoveObjects),Pointer           :: auxmv 

    Integer    :: unitno, scale_rot

    !** Nullify the cavity bais pointer
    Nullify(params%cavity)
    If (auxmv%cavBiasON)     params%cavity=>auxmv%cavity

    !** Get the delta displacements
    unitno = file_open(filename)
    Read(unitno,*) params%deltarot, scale_rot

    !** Set the simcell pointer
    params%simcell => simcell

    params%max_nrg = VERY_LRG_NRG

    !** Set the flag to decide whether we scale the jump lengths or not
    params%scale_rot = .False.
    If (scale_rot == 1) Then
      params%scale_rot = .True.
    End If

  End Subroutine rigidmoves_RRotInit

  !-----------------------------------------------------------------------
  ! Read the information that will be required to do reinsertion move.
  ! This move is essentially a combined large translation and rotation.
  ! Requires: params -- parameters to be initialized
  !           spc -- species number
  !           filename -- file where information is located
  !-----------------------------------------------------------------------
  Subroutine rigidmoves_ReinsertInit(params,spc,filename)
    Type(RigidReinsert_Params), Intent(InOut) :: params
    Integer, Intent(In)                       :: spc
    Character(*), Intent(In)                  :: filename

    params%max_nrg = VERY_LRG_NRG

  End Subroutine rigidmoves_ReinsertInit

  !-----------------------------------------------------------------------
  ! Read the information that will be required to do a Biased combination
  ! rotation and translation move.
  ! Requires:  params -- parameters to be initialized
  !            spc -- species number
  !            filename -- file where information is located
  !            nums -- alternative supply of initialization parameters
  !-----------------------------------------------------------------------
  Subroutine rigidmoves_biasrtInit(params,filename,nums)
    Type(RigidBiasRT_Params), Intent(InOut)              :: params
    Character(*), Intent(In)                             :: filename
    Real(kind=RDbl), Dimension(:), Intent(In), Optional  :: nums

    Integer    :: unit,tempint

    !** Set defaults
    params%max_nrg = VERY_LRG_NRG
    params%scale_rot = .False.
    params%scale_trans = .False.

    !** Convert the input into parameters for the move type
    If (Present(nums)) Then
      params%ntrials = Int(nums(1))
      params%biastemp = Int(nums(2))
      params%deltatrans = nums(3)
      If (Abs(nums(4) - 1.0_RDbl) <= 1.0e-5_RDbl) params%scale_trans = .True.
      params%deltarot = nums(5)
      If (Abs(nums(6) - 1.0_RDbl) <= 1.0e-5_RDbl) params%scale_rot = .True.
    Else
      !** Read the translation and rotation perturbation sizes
      unit = file_open(filename)
      Read(unit,*) params%ntrials
      Read(unit,*) params%biastemp
      Read(unit,*) params%deltatrans, tempint
      If (tempint == 1) params%scale_rot = .True.
      Read(unit,*) params%deltarot, tempint
      If (tempint == 1) params%scale_trans = .True.
    End If

    !** Calculate beta
    params%beta = 1.0_RDbl/(kcalmole_kb*params%biastemp)  !** 1.0/(Rgas*tk)

  End Subroutine rigidmoves_biasrtInit

  !------------------------------------------------------------------------
  ! Initializes the information required for doing Insertion type of move 
  ! with NO BIAS
  ! Requires: params -- Rigid biased insertion parameters
  !           simcell -- simulation cell information
  !           spc -- species number
  !------------------------------------------------------------------------
  Subroutine rigidmoves_RInsInit(params, simcell, auxmv, spc)
    Type(RigidRIns_Params), Intent(InOut)     :: params
    Type(SimCell_Params), Target, Intent(In)  :: simcell
    Integer, Intent(In)                       :: spc
    Type(AuxMoveObjects),Pointer           :: auxmv 

    params%max_nrg = VERY_LRG_NRG

    !** Set the simcell pointer
    params%simcell => simcell

    !** Check for internal flexibility
    If (molecules_isflexible(spc)) Then
      Write(0,*) " I assume you want to have rigid insertion?"
      Write(0,*) " Add the following in Molecule file for spc ", spc
      Write(0,*) "Molecule_DOF: 6"
      Write(0,*) "degree of freedom should be 6 only, &
          & (5 for 2 atom, 3 for 1 atom)"
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " ERROR: DOF is greater than 6, insertions not correctly biased"
      Stop
    End If

  End Subroutine rigidmoves_RInsInit

  !------------------------------------------------------------------------
  ! Reads the information required for doing Insertion type of move and
  ! Initilises the relevant variables
  ! Requires: params -- Rigid biased insertion parameters
  !           scell -- simulation cell information
  !           spc -- species number
  !           filename -- file to read intialization params from
  !------------------------------------------------------------------------
  Subroutine rigidmoves_BInsInit(params, scell, auxmv, spc, filename)
    Type(RigidBIns_Params), Intent(InOut)     :: params
    Type(SimCell_Params), Target, Intent(In)  :: scell
    Integer, Intent(In)                       :: spc
    Character(*), Intent(In)                  :: filename
    Type(AuxMoveObjects),Pointer           :: auxmv 

    Integer    :: unitno

    params%max_nrg = VERY_LRG_NRG  !** maximum energy HACK for now

    !** Set the simcell pointer
    params%simcell => scell

    !** Read the relevant info from ctrlfile
    unitno = file_open(filename)
    Write(*,*) "Initializing biasing details : cavity, bmap etc "
    !** Initializing biasing details : cavity, bmap etc
    Call rigidmoves_initInsBias(params%insbias,scell,auxmv,params%max_nrg,&
        params%type_of_bias,unitno)

    !** Check for internal flexibility
    If (molecules_isflexible(spc)) Then
      Write(0,*) " I assume you want to have rigid insertion?"
      Write(0,*) " Add the following in Molecule file for spc ", spc
      Write(0,*) "Molecule_DOF: 6"
      Write(0,*) "degree of freedom should be 6 only, &
          & (5 for 2 atom, 3 for 1 atom)"
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " ERROR: DOF is greater than 6, insertions not correctly biased"
      Stop
    End If

  End Subroutine rigidmoves_BInsInit

  !------------------------------------------------------------------------
  ! Reads the information required for doing Insertion type of move and
  ! Initilises the relevant variables for biased insertions 
  ! of molecules, whose configurations are taken from a previously made 
  ! library
  !------------------------------------------------------------------------
  Subroutine rigidmoves_LInsInit(params,scell,auxmv,spc,ctrlfile)
    Type(RigidLIns_Params), Intent(InOut)     :: params
    Type(SimCell_Params), Target, Intent(In)  :: scell
    Type(AuxMoveObjects),Pointer           :: auxmv 
    Integer, Intent(In)                       :: spc
    Character(*), Intent(In)                  :: ctrlfile

    Integer                       :: unitno, error, nfields    
    Character(len=strLen)         :: biasfilename,confFileName, molecname
    Character(len=strLen)         :: tempstring, spc_origname
    Character(len=strLen),Dimension(strLen)  :: fields 
    Logical                       :: withVels

    !** a very useful check
    If (.Not.molecules_isflexible(spc)) Then
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Some Routine may not work for rigid molecules"
      Stop    
    End If

    !** Set the simcell pointer
    params%simcell => scell

    params%max_nrg = VERY_LRG_NRG

    !** Read the relevant info from ctrlfile
    unitno = file_open(ctrlfile)

    Call rigidmoves_initInsBias(params%insbias, scell, auxmv, params%max_nrg,&
        params%type_of_bias,unitno)

    Allocate(params%lib,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Read datafile from which library will be initialized
    !** it is possible that a different name was used for the species 
    !** during library generations, it is also possible that we 
    !** need the chiral image of original library
    !** The format to be read: filename, spc_orig_name, chiral
    Read(unitno,'(a)') tempstring
    tempstring = cleanstring(tempstring)
    nfields = split(tempstring, fields, ",")
    confFileName = cleanstring(fields(1))
    params%confFileName = confFileName
    If (nfields > 1) Then
      spc_origname=cleanstring(fields(2))
    Else
      spc_origname=molecules_name(spc)
    End If

    !** Initialize the library with a single species type
    withVels = .True.
    molecname = molecules_name(spc)
    Call conflib_basicInit(params%lib,(/molecname/),withVels, &
        .True.,(/spc_origname/))
    Call conflib_dataInit(params%lib,molecname,confFileName)

    !** Set the 'chiral image' flag.  What's this?!?
    params%lib%chiral_image = .False.
    If (nfields>2) Then
      If (tolower(Trim(cleanstring(fields(3))))=="chiral") Then
        !** HACK, reaches into library
        params%lib%chiral_image = .True.
      Else
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,*) "wrong string here, expecting 'chiral': ", fields(3)
        Stop
      End If
    End If

    !** Initialize the temp_ref_struc to be used for copying from library
    Allocate(params%temp_ref_struc,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Call config_allocfields(params%temp_ref_struc,spc,1)

  End Subroutine rigidmoves_LInsInit

  !------------------------------------------------------------------------
  ! Reads the information required for doing Deletion type of move and
  ! Initializes the relevant variables- for BIASED CASE
  ! Requires: params -- Rigid biased insertion parameters
  !           simcell -- simulation cell information
  !           spc -- species number
  !           filename -- file to read initialization parameters from
  !------------------------------------------------------------------------
  Subroutine rigidmoves_BDelInit(params,scell,spc,filename)
    Type(RigidBDel_Params), Intent(InOut)     :: params
    Type(SimCell_Params), Target, Intent(In)  :: scell
    Integer, Intent(In)                       :: spc

    Character(*), Intent(in)    :: filename

    Integer    :: unitno

    !** Set the simcell pointer
    params%simcell => scell
    params%max_nrg= VERY_LRG_NRG

    unitno = file_open(filename)

    ! usually this is not required for a regular run
    ! if we are debuggin with only delete and noinsert then, it might be required
#ifdef DEBUG
    !** intialise the biasing info : this is required ?
    !** this should be done while equatins insertions deletions
    !** but then we should be able to have independent deletion moves
    !** for debugging purposes
    Call rigidmoves_initInsBias(params%insbias,scell,params%max_nrg,&
        params%type_of_bias,unitno)
#endif

  End Subroutine rigidmoves_BDelInit

  !------------------------------------------------------------------------
  ! Reads the information required for doing Deletion type of move and
  ! Initilises the relevant variables - for unbiased case.
  ! Requires: params -- Rigid biased insertion parameters
  !           simcell -- simulation cell information
  !           spc -- species number
  !------------------------------------------------------------------------
  Subroutine rigidmoves_RDelInit(params,simcell,spc)
    Type(RigidRDel_Params), Intent(InOut)     :: params
    Type(SimCell_Params), Target, Intent(In)  :: simcell
    Integer, Intent(In)                       :: spc

    params%max_nrg = VERY_LRG_NRG

    !** Set the simcell pointer
    params%simcell => simcell


  End Subroutine rigidmoves_RDelInit

  !------------------------------------------------------------------------
  ! Initialises the specifics required for insertion biases : energy 
  ! map bias, cavity bias. This is intended for use by _LInsinit, _Binsinit
  ! ans _BDelInit
  !------------------------------------------------------------------------
  Subroutine rigidmoves_initInsBias(params, scell, auxmv, max_nrg, bias_type, &
      unitno)
    Type(RigidInsBias), Pointer  :: params
    Type(SimCell_Params), Target, Intent(In)  :: scell
    Type(AuxMoveObjects),Pointer           :: auxmv 
    Real(kind=RDbl), Intent(in) :: max_nrg
    Character(*), Intent(out)                 :: bias_type 
    Integer, Intent(in) :: unitno

    Character(len=lstrLen)   :: text
    Integer    :: error, nfields
    Character(len=strLen), Dimension(10) :: fields
    Real(kind=RDbl), Dimension(3) :: eff
    Real(kind=RDbl)               :: tk

    Allocate(params,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Set the simcell pointer
    params%simcell => scell

    params%max_nrg=max_nrg
    Nullify(params%bmap)
    Nullify(params%cavity)

    !** if no bmap then no bias, -negative temp implies no bias...
    params%bias_temp= -1.00

    !** Now decide what exactly is the type of bias
    Read(unitno,'(a)') text
    text=Trim(stripcmnt(text))
    nfields=split(text,fields,",")
    If (nfields==1) Then
      If(Trim(fields(1))=="CAVITY") Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        ! cavity bias but no bmap
        bias_type="CAVITY_BIAS"
        params%biasfilename="NULL"
        Allocate(params%cavity,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        Read(unitno,*) params%cavNTrials
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Else
        ! initialise simple bmap bias
        Allocate(params%bmap,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        bias_type="MAP_BIAS"
        params%biasfilename=fields(1)
        Read(unitno,*) tk
        params%bias_temp=tk
        eff = simcell_geteff(scell, .True.)
        Call bmap_init(params%bmap, eff, tk, params%biasfilename, scell)

      Endif
    Elseif (nfields==2) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      If(Adjustl(Trim(fields(2)))=="CAVITY") Then
        ! initialise bmap+cavity
        Allocate(params%bmap,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        Allocate(params%cavity,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        bias_type="MAP_CAVITY_BIAS"
        params%biasfilename=fields(1)
        Read(unitno,*) tk
        params%bias_temp=tk
        Read(unitno,*) params%cavNTrials
        eff = simcell_geteff(scell, .True.)
        Call bmap_init(params%bmap, eff, tk, params%biasfilename, scell)
      Else
        Write(*,*) " In the control file put bmapname and CAVITY &
            & keyword in proper order, CAVITY Key word second"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End If
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    params%type_of_bias = bias_type

  End Subroutine rigidmoves_initInsBias


  !---------------- Move routines -----------------------------------------

  !------------------------------------------------------------------------
  ! Pick a RANDOM displacement based on the move parameters, TRANSLATE 
  ! the center-of-mass of the molecule and generate xyzcoords.
  ! Requires:  species -- species data structure
  !            spc -- species number
  !            mol -- molecule number
  !            simcell -- simulation cell information
  !            params -- move parameters
  !            succ_flag -- success flag, returns TRUE if successful
  !            subint -- subset interactions (molecule-system)
  !------------------------------------------------------------------------
  Subroutine rigidmoves_randomtranslate(params,species,spc,mol,simcell, &
      succ_flag,subint)
    Type(RigidRTrans_Params), Intent(In)           :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Type(SimCell_Params), Intent(In)               :: simcell
    Logical, Intent(Out)                           :: succ_flag
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Integer             :: natoms    
    Logical             :: fast, skip_intra_flag
    Type(VecType)       :: dvec, com
    Real(kind=RDbl)     :: deltatrans, xdelta, ydelta, zdelta

    natoms = molecules_getnatoms(spc)
    deltatrans = params%deltatrans

    !** Calculate the energy of molecule at its current position
    fast = .True.
    skip_intra_flag = .True.
    succ_flag = subinteract_int(subint,species,simcell,fast, &
        .False.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))
    Call checkandstop(succ_flag,__FILE__,__LINE__, &
        'Problem with nrg calculation of an existing structure')

    !** It is possible that the molecule is internally flexible 
    !** Then we need to re-initialize the gcoords
    If (species(spc)%gcoords(mol)%rigidcoords%internally_flexible) Then
      com = config_getMolecCOM(species, spc, mol)
      Call rigidcoords_setfrom_rp(species(spc)%gcoords(mol)%rigidcoords, spc, &
          species(spc)%coords(1:natoms,mol)%rp,com)
    End If

    If (params%constrained) Then
      !** Do a coordinate transformation to the constrained plane
      !** system where the normal is along the x-axis
      !***DEBUG STUFF
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) " is this part initialised?"
      Stop
      com = species(spc)%gcoords(mol)%rigidcoords%com      
      com = params%constrtrans*com

      !** Get the displacement in the y-z plane
      ydelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      zdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      dvec = (/0.0_RDbl, ydelta, zdelta/)
      com = com + dvec

      !** Transform back to the original coordinate system
      species(spc)%gcoords(mol)%rigidcoords%com = &
          params%constrinvtrans*com

      !** Generate the xyz coordinates
      Call rigidcoords_toxyz(species(spc)%gcoords(mol)%rigidcoords,&
          species(spc)%coords(1:natoms,mol)%rp)
  
      !** do periodic b. condns before energy calcs
      !** Update the other arrays in species using PBCs
      Call simcell_pbc(simcell, &
          species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms,mol)%r, &
          species(spc)%coords(1:natoms,mol)%cr)

    Else  !** normal method
      xdelta = ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      ydelta = ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      zdelta = ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      dvec = (/xdelta, ydelta, zdelta/)
      Call rigidmoves_transperturb(species,spc,mol,simcell,dvec,.True.)

    End If

    !** Evaluate new energy
    succ_flag = subinteract_int(subint,species,simcell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))

  End Subroutine rigidmoves_randomtranslate

  !-----------------------------------------------------------------------
  ! ROTATE the molecule by generating new RANDOM eulerian angles and 
  ! regenerates xyzcoords.
  ! Requires:  params -- move parameters
  !            species -- species data structure
  !            spc -- species number
  !            mol -- molecule number
  !            simcell -- simulation cell information
  !            succ_flag -- success flag, returns TRUE if successful
  !            subint -- subset interactions (molecule-system)
  !-----------------------------------------------------------------------
  Subroutine rigidmoves_randomrotate(params,species,spc,mol,simcell, &
      succ_flag, subint)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: species
    Integer, Intent(in)                            :: spc,mol
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(RigidRRot_Params), Intent(in)             :: params
    Logical, Intent(out)                           :: succ_flag
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Integer             :: natoms
    Logical             :: fast, skip_intra_flag, rotateVels
    Real(kind=RDbl)     :: deltarot, dphi, dpsi, dcostheta, phi, psi, theta
    Type(RigidMolec), Pointer :: rigidcoord

    !** Calculate the energy of molecule at its current position
    fast = .True.
    skip_intra_flag = .True.
    succ_flag = subinteract_int(subint,species,simcell,fast, &
        .False.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))

    !** Check the no. of atoms and warn if no. of atoms is 1
    !** The energy update structure works on the assumption that moves returns
    !** both new and old energies, so you can not return if natoms=1 
    !** before calculating energies
    natoms = molecules_getnatoms(spc)
    If (natoms == 1) Then
      Write(0,*) "make sure this is some debug run, and not a real one"
      Write(0,*) "No need to rotate a 1-atom molecule"
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    End If

    If (.Not.succ_flag) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "Problem with nrg calculation of an existing structure"
      Write(0,*) 'Maximum energy: ',params%max_nrg
      Write(0,*) 'Returned energy: ',subinteract_oldnrg(subint,.False.)
      Stop
    End If

    !** for clarity
    rigidcoord => species(spc)%gcoords(mol)%rigidcoords

    deltarot = params%deltarot

    If (rigidcoord%internally_flexible) Then
      !** It is possible that the molecule is internally flexible 
      ! the gencoords sytem is pretty much useless here
      ! the rotations are done only in %rp, %v

      !** choose a unit vector direction randomly
      phi = rranf()*pi
      theta = Acos( (rranf()*2.0_RDbl) - one )

      !** find the psi angle by which you want to rotate the molecule
      psi = (rranf()*2.0_RDbl - 1.0_RDbl)*deltarot

      rotateVels = .True.
      !** rotate the values in %rp, %v by psi around unitvector along phi, theta
      Call config_rotateAroundCOM(species, spc, mol, phi, theta, psi, &
          rotateVels)

      !** Apply PBCs to the new xyz coordinates
      Call simcell_pbc(simcell, &
          species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms,mol)%r, &
          species(spc)%coords(1:natoms,mol)%cr)

    Else
      !** Sample phi and psi uniformly between -pi,pi but theta
      !** needs to be sampled from a cosine distribution
      !** See Allen and Tildesley, Pg. 132-133, "Molecular Liquids"
      dphi = ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
      dpsi = ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
      dcostheta = ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
      Call rigidmoves_rotperturb(species,spc,mol,simcell, &
          dphi,dpsi,dcostheta,.True.)

    End If

    !** Evaluate new energy
    succ_flag = subinteract_int(subint,species,simcell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))

  End Subroutine rigidmoves_randomrotate

  !------------------------------------------------------------------------
  ! Performs a REINSERTION move.  This is essentially a reinsertion an 
  ! existing molecules in a random orientation in a random location in the
  ! simulation cell.  biasfactor = 1.0 in this case.  Essentially identical
  ! to the random_insertion move except that it operates on an existing 
  ! molecule.  This move is used in simulations where N_molecs is fixed.
  ! Requires:  params -- insertion parameters
  !            species -- species data structure for whole system
  !            spc -- species number
  !            mol -- molecule number
  !            simcell -- simulation cell information
  !            biasfactor -- returned bias factor
  !            succ_flag -- success flag, returns TRUE if successful
  !            subint -- subset interactions (molecule-system)
  !------------------------------------------------------------------------
  Subroutine rigidmoves_reinsert(params,species,spc,mol,simcell, &
      biasfactor,succ_flag,subint)
    Type(RigidReinsert_Params), Intent(In)         :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Type(SimCell_Params), Intent(In)               :: simcell
    Real(kind=RDbl), Intent(Out)                   :: biasfactor
    Logical, Intent(Out)                           :: succ_flag
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Integer              :: natoms
    Logical              :: fast,skip_intra_flag
    Real(kind=RDbl)      :: costheta
    Type(VecType)        :: com

    !** Calculate the energy of molecule at its current position
    fast = .True.
    skip_intra_flag = .True.
    succ_flag = subinteract_int(subint,species,simcell,fast, &
        .False.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))

    If (.Not.succ_flag) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "Problem with nrg calculation of an existing structure"
      Write(0,*) 'Maximum energy: ',params%max_nrg
      Write(0,*) 'Returned energy: ',subinteract_oldnrg(subint,.False.)
      Stop
    End If

    !** bias factor = N_c * bias  =  N_c * (1/ N_c) = 1.0
    biasfactor = one

    !** Get the position of the center of mass unifromly from the 
    !** simulation volume
    com = simcell_uniformRandPt(simcell)

    !** 3.) Translate the pt. to the appropriate unit cell
    species(spc)%gcoords(mol)%rigidcoords%com = com
    species(spc)%gcoords(mol)%rigidcoords%theta = 0.0_RDbl
    species(spc)%gcoords(mol)%rigidcoords%phi   = 0.0_RDbl
    species(spc)%gcoords(mol)%rigidcoords%psi   = 0.0_RDbl

    !** If natoms is = 1, we are done
    natoms = molecules_getnatoms(spc)
    If (natoms /= 1) Then
      !** Generate the eulerian angles at random
      !** Sample theta from a cosine distribution.
      costheta = 1 - rranf() * 2.0
      species(spc)%gcoords(mol)%rigidcoords%theta = Acos(costheta)
      species(spc)%gcoords(mol)%rigidcoords%phi = pi - rranf()*twopi
      species(spc)%gcoords(mol)%rigidcoords%psi = pi - rranf()*twopi
    End If

    !** Generate the xyz coordinates
    Call rigidcoords_toxyz(species(spc)%gcoords(mol)%rigidcoords,&
        species(spc)%coords(1:natoms,mol)%rp)

    !** apply periodic boundary conditions before forcefield evaluation
    Call simcell_pbc(simcell, &
        species(spc)%coords(1:natoms,mol)%rp, &
        species(spc)%coords(1:natoms,mol)%r, &
        species(spc)%coords(1:natoms,mol)%cr)

    !** get new energies
    succ_flag = subinteract_int(subint,species,simcell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))

  End Subroutine rigidmoves_reinsert

  !------------------------------------------------------------------------
  ! Performs random translations and rotations from a given starting 
  ! configuration until either the maximum number of moves is reached
  ! or the energy is calculable.  This is useful for taking a post-idchange
  ! configuration that does not have a calculable energy and converting it
  ! to a more reasonable configuration.
  ! NOTE: for now, it uses the biased rot/trans parameters to specify
  ! the size of the random translations and rotations, but it could be
  ! made independent.
  ! Requires:  params -- biased rotation/translation parameters
  !            maxmoves -- maximum number of 'retrys'
  !            species -- species data structure for whole system
  !            spc -- species number
  !            mol -- molecule number
  !            simcell -- simulation cell information
  !            succ_flag -- success flag, returns TRUE if successful
  !            subint -- subset interactions (molecule-system)
  !------------------------------------------------------------------------
  Subroutine rigidmoves_rtretry(params,maxmoves,species,spc,mol,simcell, &
      succ_flag,subint)
    Type(RigidBiasRT_Params), Intent(In)           :: params
    Integer, Intent(In)                            :: maxmoves
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Type(SimCell_Params), Intent(In)               :: simcell
    Logical, Intent(Out)                           :: succ_flag
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Integer              :: natoms,i,n
    Logical              :: fast,no_intra,success
    Real(kind=RDbl)      :: dphi,dpsi,dcostheta,dx,dy,dz,randnum,cum,nrg
    Type(VecType)        :: dvec
    Type(RigidMolec), Pointer         :: gcoords
    Real(kind=RDbl), Dimension(6)     :: genarray

    !** Set defaults
    no_intra = .True.   !** ignore intramolecular energies
    success = .False.   !** False => no trial moves were under max_nrg
    gcoords => species(spc)%gcoords(mol)%rigidcoords
    natoms = molecules_getnatoms(spc)

    !** Save the initial generalized coordinates
    Call rigidcoords_getgencoords(gcoords,genarray)

    !** Generate random moves from the initial coordinates
    success = .False.
    Do i = 1,maxmoves
      !** Reset configuration to the original
      If (i /= 1) Then
        Call rigidcoords_setgencoords(gcoords,genarray)
      End If

      !** Perturb rotational generalized coordinates
      dphi = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltarot
      dpsi = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltarot
      dcostheta = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltarot
      Call rigidmoves_rotperturb(species,spc,mol,simcell, &
          dphi,dpsi,dcostheta,.False.)

      !** Perturb COM generalized coordinates and generate xyz configuration
      dx = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltatrans
      dy = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltatrans
      dz = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltatrans
      dvec = (/dx,dy,dz/)
      Call rigidmoves_transperturb(species,spc,mol,simcell,dvec,.True.)

      !** Generate full xyz coordinate, apply PBCs and get energy
      Call rigidcoords_setgencoords(gcoords,genarray)
      Call rigidcoords_toxyz(gcoords,species(spc)%coords(1:natoms,mol)%rp)
      Call simcell_pbc(simcell,species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms,mol)%r, &
          species(spc)%coords(1:natoms,mol)%cr)
      success = subinteract_int(subint,species,simcell,fast, &
          .True.,no_intra,(/params%max_nrg/),(/spc,mol,0/))
      If (success) Exit
    End Do 

    !** Make sure that at least one of the trials was successful else undo
    If (.Not. success) Then
      Call rigidcoords_setgencoords(gcoords,genarray)
      Call rigidcoords_toxyz(gcoords,species(spc)%coords(1:natoms,mol)%rp)

      !** Apply periodic boundary conditions before forcefield evaluation
      Call simcell_pbc(simcell,species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms,mol)%r, &
          species(spc)%coords(1:natoms,mol)%cr)

      !** Make sure energy has been evaluated 
      succ_flag = subinteract_int(subint,species,simcell,fast, &
          .True.,no_intra,(/params%max_nrg/),(/spc,mol,0/))
    End If

  End Subroutine rigidmoves_rtretry

  !------------------------------------------------------------------------
  ! Performs a biased combination rotation and perturbation using the 
  ! Rosenbluth weightings.  (old config = o; new config = n)
  !   prob. generating o->n        W_(o->n)    f(o)
  !   ------------------------- = --------- x ------- = biasfactor
  !   prob. generating n->o        W_(n->o)    f(n)
  ! 
  ! Where:
  !   W_(o->n) = Rosenbluth factor for moving from state o to state n
  !            = Sum(k=1,nconfigs) f(k)
  !   W_(n->o) = Rosenbluth factor for moving from state n to state o
  !            = f(o) + Sum(l=1,nconfigs-1) f(l)
  !       f(k) = Exp(-beta*V(k)), Boltzmann factor for state k
  ! See Frenkel&Smit pg.274 for a confusing explanation
  !
  ! Requires:  params -- biased rotation/translation parameters
  !            species -- species data structure for whole system
  !            spc -- species number
  !            mol -- molecule number
  !            simcell -- simulation cell information
  !            biasfactor -- returned bias factor
  !            succ_flag -- success flag, returns TRUE if successful
  !            subint -- subset interactions (molecule-system)
  !------------------------------------------------------------------------
  Subroutine rigidmoves_biasrt(params,species,spc,mol,simcell,biasfactor, &
      succ_flag,subint)
    Type(RigidBiasRT_Params), Intent(In)           :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Type(SimCell_Params), Intent(In)               :: simcell
    Real(kind=RDbl), Intent(Out)                   :: biasfactor
    Logical, Intent(Out)                           :: succ_flag
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Integer              :: natoms,i,n
    Logical              :: fast,no_intra,success
    Real(kind=RDbl)      :: Wo2n,Wn2o
    Real(kind=RDbl)      :: dphi,dpsi,dcostheta,dx,dy,dz,randnum,cum,nrg
    Type(VecType)        :: dvec
    Type(RigidMolec), Pointer                      :: gcoords
    Real(kind=RDbl), Dimension(params%ntrials)     :: fo2n,fn2o
    Real(kind=RDbl), Dimension(params%ntrials,6)   :: o2nset,n2oset

    !** Set defaults
    fast = .True.  !** HaCK
    Wo2n = 0.0_RDbl
    biasfactor = 1.0_RDbl
    no_intra = .True.   !** ignore intramolecular energies
    success = .False.   !** False => no trial moves were under max_nrg
    gcoords => species(spc)%gcoords(mol)%rigidcoords
    natoms = molecules_getnatoms(spc)

    !** Save the initial generalized coordinates and get weighting factor
    fn2o(1) = rigidmoves_getwt(params,species,spc,mol,simcell, &
        succ_flag,subint,no_intra,nrg)
    Wn2o = fn2o(1)
    Call rigidcoords_getgencoords(gcoords,n2oset(1,:))
!    Write(*,*) 'incoming energy, beta: ',nrg,params%beta

    !** Return now if the weighting factor is zero (ie above max_nrg)
    !** In this case, we can't do this biasing because biasfactor is undefined
    If (Abs(fn2o(1)) < 1.0e-12_RDbl) Return

    !** Generate ntrials configurations going from o->n
    Do i = 1,params%ntrials
      !** Reset configuration to the original
      If (i /= 1) Then
        Call rigidcoords_setgencoords(gcoords,n2oset(1,:))
      End If

      !** Perturb rotational generalized coordinates
      dphi = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltarot
      dpsi = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltarot
      dcostheta = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltarot
      Call rigidmoves_rotperturb(species,spc,mol,simcell, &
          dphi,dpsi,dcostheta,.False.)

      !** Perturb COM generalized coordinates and generate xyz configuration
      dx = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltatrans
      dy = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltatrans
      dz = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltatrans
      dvec = (/dx,dy,dz/)
      Call rigidmoves_transperturb(species,spc,mol,simcell,dvec,.True.)

      !** Save the generalized coordinates of the configuration
      Call rigidcoords_getgencoords(gcoords,o2nset(i,:))

      !** Get weighting factor 
      fo2n(i) = rigidmoves_getwt(params,species,spc,mol,simcell, &
          succ_flag,subint,no_intra,nrg)
      If (succ_flag) success = .True.  
      Wo2n = Wo2n + fo2n(i)
!      Write(*,*) 'forward: ',i,nrg,succ_flag,fo2n(i),Wo2n
    End Do 

    !** Make sure that at least one of the trials was successful
    If (.Not. success) Then
      !** Reset configuration and return without finishing biasing
      Call rigidcoords_setgencoords(gcoords,n2oset(1,:))
      Call rigidcoords_toxyz(gcoords,species(spc)%coords(1:natoms,mol)%rp)

      !** Apply periodic boundary conditions before forcefield evaluation
      Call simcell_pbc(simcell,species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms,mol)%r, &
          species(spc)%coords(1:natoms,mol)%cr)

      !** Make sure energy has been evaluated 
      succ_flag = subinteract_int(subint,species,simcell,fast, &
          .True.,no_intra,(/params%max_nrg/),(/spc,mol,0/))
      Return

      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " WARNING: none of the o->n trial moves had energies under max_nrg"
      Write(0,'(2x,a)') " Consider increasing ntrials"
    End If

    !** Pick a new configuration (n) from the list
    !** multiplication by the Wo2n factor avoids normalization of fo2n
    randnum = rranf()*Wo2n
    cum = 0.0_RDbl
    Do n = 1,params%ntrials
      cum = cum + fo2n(n)
      If (cum >= randnum) Exit
    End Do
!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,*) 'picked configuration: ',n

    !** Generate ntrials-1 configurations going from n->o
    success = .False.
    Do i = 2,params%ntrials
      !** Reset configuration to the new configuration
      Call rigidcoords_setgencoords(gcoords,o2nset(n,:))

      !** Perturb rotational generalized coordinates
      dphi = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltarot
      dpsi = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltarot
      dcostheta = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltarot
      Call rigidmoves_rotperturb(species,spc,mol,simcell, &
          dphi,dpsi,dcostheta,.False.)

      !** Perturb COM generalized coordinates and generate xyz configuration
      dx = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltatrans
      dy = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltatrans
      dz = ((rranf()*2.0_RDbl) - 1.0_RDbl)*params%deltatrans
      dvec = (/dx,dy,dz/)
      Call rigidmoves_transperturb(species,spc,mol,simcell,dvec,.True.)

      !** Get weighting factor 
      fn2o(i) = rigidmoves_getwt(params,species,spc,mol,simcell, &
          succ_flag,subint,no_intra,nrg)
      If (succ_flag) success = .True.  
      Wn2o = Wn2o + fn2o(i)
!      Write(*,*) 'reverse: ',i,nrg,succ_flag,fn2o(i),Wn2o
    End Do

    !** Make sure that at least one of the reverse trials was successful
    If (.Not. success) Then
      !** Reset configuration and return without finishing biasing
      Call rigidcoords_setgencoords(gcoords,n2oset(1,:))
      Call rigidcoords_toxyz(gcoords,species(spc)%coords(1:natoms,mol)%rp)

      !** Apply periodic boundary conditions before forcefield evaluation
      Call simcell_pbc(simcell,species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms,mol)%r, &
          species(spc)%coords(1:natoms,mol)%cr)

      !** Make sure energy has been evaluated 
      succ_flag = subinteract_int(subint,species,simcell,fast, &
          .True.,no_intra,(/params%max_nrg/),(/spc,mol,0/))
      Return

      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " WARNING: none of the o->n trial moves had energies under max_nrg"
      Write(0,'(2x,a)') " Consider increasing ntrials"
    End If 

    !** Calculate the final biasing factor (see notes above)
    biasfactor = (Wo2n/Wn2o)*(fn2o(1)/fo2n(n))

    !** Reset configuration to the new configuration
    Call rigidcoords_setgencoords(gcoords,o2nset(n,:))
    Call rigidcoords_toxyz(gcoords,species(spc)%coords(1:natoms,mol)%rp)

    !** Apply periodic boundary conditions before forcefield evaluation
    Call simcell_pbc(simcell,species(spc)%coords(1:natoms,mol)%rp, &
        species(spc)%coords(1:natoms,mol)%r, &
        species(spc)%coords(1:natoms,mol)%cr)

    !** Make sure energy has been evaluated 
    succ_flag = subinteract_int(subint,species,simcell,fast, &
        .True.,no_intra,(/params%max_nrg/),(/spc,mol,0/))

#ifdef DEBUG
    nrg = subinteract_newnrg(subint,no_intra)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Wo2n, Wn2o: ',Wo2n,Wn2o
    Write(*,*) 'f(o),f(n): ',fn2o(1),fo2n(n)
    Write(*,*) 'Final biasing factor: ',biasfactor
    Write(*,*) 'Final energy of biased configuration: ',nrg
#endif

  End Subroutine rigidmoves_biasrt

  !------------------------------------------------------------------------
  ! Gets the weighting function for a specified molecule.  Used by _biasrt
  ! Requires:  params -- biased rotation/translation parameters
  !            species -- species data structure for whole system
  !            spc -- species number
  !            mol -- molecule number
  !            simcell -- simulation cell information
  !            succ_flag -- success flag, returns TRUE if successful
  !            subint -- subset interactions (molecule-system)
  !            no_intra -- flags evaluation of intramolecular energies
  !            energy -- optional returned energy
  !------------------------------------------------------------------------
  Real(Kind=RDbl) Function rigidmoves_getwt(params,species,spc,mol,simcell, &
      succ_flag,subint,no_intra,energy)
    Type(RigidBiasRT_Params), Intent(In)           :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Type(SimCell_Params), Intent(In)               :: simcell
    Logical, Intent(Out)                           :: succ_flag
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Logical, Intent(In)                            :: no_intra
    Real(kind=RDbl), Intent(Out), Optional         :: energy

    Logical              :: fast = .True.  !** HaCK
    Real(kind=RDbl)      :: nrg,toexp

    !** Default
    rigidmoves_getwt = 0.0_RDbl
    nrg = params%max_nrg

    !** Evaluate energy
    succ_flag = subinteract_int(subint,species,simcell,fast, &
        .True.,no_intra,(/params%max_nrg/),(/spc,mol,0/))
    
    !** Calculate the weighting function for this configuration
    If (succ_flag) Then
      nrg = subinteract_newnrg(subint,no_intra)
    
      toexp = -1.0_RDbl*params%beta*nrg
      If (toexp > default_MAX_EXP) toexp = default_MAX_EXP   
      If (toexp < -default_MAX_EXP) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
            " WARNING: using underflow damping in rigidmoves_getwt"
        Write(0,'(2x,a)') " Consider increasing the bias temperature"
        toexp = -default_MAX_EXP
      End If
      rigidmoves_getwt = Exp(toexp)
    End If

    If (Present(energy)) energy = nrg

  End Function rigidmoves_getwt

  !------------------------------------------------------------------------
  ! Generates a center-of-mass for the molecule uniformly from the 
  ! simulation volume, and with random eulerian angles generates the xyz 
  ! coordinates.  biasfactor = 1.0 in this case
  ! Requires: species -- species data structure for whole system
  !           spc -- species number
  !           mol -- molecule number
  !           params -- insertion parameters
  !           biasfactor -- returned bias factor
  !           succ_flag -- success flag, false if energies too high
  !           subint -- subset interactions (molecule-system)
  ! NOTE: need to rewrite this routine so it returns energies etc.
  !------------------------------------------------------------------------
  Subroutine rigidmoves_rinsert(params, species, spc, mol, biasfactor, &
      succ_flag,subint)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Type(RigidRIns_Params), Intent(In)             :: params
    Real(kind=RDbl), Intent(Out)                   :: biasfactor
    Logical, Intent(Out)                           :: succ_flag    
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Integer              :: natoms
    Logical              :: fast, no_intra
    Real(kind=RDbl)      :: costheta
    Type(VecType)        :: com

    !** bias factor = N_c *  bias  = N_c *  (1/ N_c) =one.
    biasfactor = one

    !** Get the position of the center of mass unifromly from the 
    !** simulation volume
    com = simcell_uniformRandPt(params%simcell)

    !** Translate the pt. to the appropriate unit cell
    species(spc)%gcoords(mol)%rigidcoords%com = com
    species(spc)%gcoords(mol)%rigidcoords%theta = 0.0_RDbl
    species(spc)%gcoords(mol)%rigidcoords%phi   = 0.0_RDbl
    species(spc)%gcoords(mol)%rigidcoords%psi   = 0.0_RDbl

    !** If natoms is = 1, we are done
    natoms = molecules_getnatoms(spc)
    If (natoms /= 1) Then
      !** Generate the eulerian angles at random
      !** Sample theta from a cosine distribution.
      costheta = 1 - rranf() * 2.0
      species(spc)%gcoords(mol)%rigidcoords%theta = Acos(costheta)
      species(spc)%gcoords(mol)%rigidcoords%phi = pi - rranf()*twopi
      species(spc)%gcoords(mol)%rigidcoords%psi = pi - rranf()*twopi
    End If

    !** Generate the xyz coordinates
    Call rigidcoords_toxyz(species(spc)%gcoords(mol)%rigidcoords,&
        species(spc)%coords(1:natoms,mol)%rp)

    !** Apply periodic boundary conditions before forcefield evaluation
    Call simcell_pbc(params%simcell, &
        species(spc)%coords(1:natoms,mol)%rp, &
        species(spc)%coords(1:natoms,mol)%r, &
        species(spc)%coords(1:natoms,mol)%cr)

    !** Get the new energies
    fast = .True.
    no_intra = .False.  !** could be some intra interactions
    succ_flag = subinteract_int(subint,species,params%simcell,fast, &
        .True.,no_intra,(/params%max_nrg/),(/spc,mol,0/))

  End Subroutine rigidmoves_rinsert

  !----------------------------------------------------------------------
  ! Generates a center-of-mass of the molecule according to a
  ! bias map , cavity bias OR both, and then using random eulerian angles 
  ! generates the xyz coordinates.  
  ! It also returns :
  ! 1) the bias of the insertion in "biasfactor", 
  ! 2) energies of the inserted molecule (Only pairwise energies, assumes 
  !    that intra energies are zero )
  ! Requires: species -- the species coordinates
  !           spc -- the specific species number
  !           mol -- the molecule number
  !           params -- insertion parameters
  !           biasfactor -- the returned bias factor
  !           succ_flag -- success flag, false if energies too high
  !           subint -- subset interactions (molecule-system)
  ! NOTE: "fast" is hardcoded to be .true.
  !----------------------------------------------------------------------
  Subroutine rigidmoves_binsert(params, species, spc, mol, biasfactor, &
      succ_flag,subint)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Real(kind=RDbl), Intent(Out)                   :: biasfactor
    Type(RigidBIns_Params), Intent(InOut)          :: params
    Logical, Intent(Out)                           :: succ_flag    
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Integer         :: totmol, natoms
    Logical         :: fast, skip_intra_flag, newxyz

    natoms=molecules_getnatoms(spc)

    !** Get the COM in the simcell using one of the biasing methods
    !** also get the biasfactor.
    Select Case (Trim(params%type_of_bias))
    Case("MAP_BIAS")

      species(spc)%gcoords(mol)%rigidcoords%com = &
          rigidmoves_comFromBmap(params%insbias%bmap, params%simcell , &
          biasfactor )

      !** assign random eularian angles and get new values of %rp
      newxyz=.True.
      Call rigidmoves_reorient(species(spc),spc,mol,newxyz)

    Case("CAVITY_BIAS","MAP_CAVITY_BIAS")
      !** O.K, we have cavity bias, we need to make many trial attempts 
      !** and one of the trials will be accepted. The accepted trial config
      !** is returned in species
      totmol=config_getTotMoles(species)
      Call rigidmoves_getConfigFromCavBias(params%insbias,species, &
          spc,mol,totmol, biasfactor, succ_flag )
      If (.Not.succ_flag) Return

    Case Default
      Write(*,*) "This type of rigid bias insertion not yet coded fullly"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select

    !** apply periodic boundary conditions before forcefield evaluation
    Call simcell_pbc(params%simcell, &
        species(spc)%coords(1:natoms,mol)%rp, &
        species(spc)%coords(1:natoms,mol)%r, &
        species(spc)%coords(1:natoms,mol)%cr)

    !** Get the new energies
    fast = .True.
    skip_intra_flag = .False.  !** could be some intra interactions

    succ_flag = subinteract_int(subint,species,params%simcell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))
    ! Call checkandstop(succ_flag,__FILE__,__LINE__,'insert calc problem')

  End Subroutine rigidmoves_binsert

  !----------------------------------------------------------------------
  ! Generates a center-of-mass of the molecule according to a
  ! bias map , cavity bias OR both, and then uses random eulerian angles  
  ! and inserts a configuration taken from a library
  ! It also returns :
  ! 1) the bias of the insertion in "biasfactor", 
  ! 2) energies of the inserted molecule (Only pairwise energies, assumes 
  !    that intra energies are zero )
  ! Requires: species -- the species coordinates
  !           spc -- the specific species number
  !           mol -- the molecule number
  !           params -- insertion parameters
  !           biasfactor -- the returned bias factor
  !           succ_flag -- success flag, false if energies too high
  !           subint -- subset interactions (molecule-system)
  ! NOTE: "fast" is hardcoded to be .true.
  !----------------------------------------------------------------------
  Subroutine rigidmoves_linsert(params, species, spc, mol, biasfactor, &
      succ_flag,subint)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Real(kind=RDbl), Intent(Out)                   :: biasfactor
    Type(RigidLIns_Params), Intent(InOut)          :: params
    Logical, Intent(Out)                           :: succ_flag    
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Integer         :: totmol, natoms, a
    Real(kind=Rdbl) :: intra 
    Logical         :: fast, skip_intra_flag, rotatevels, copy_accels
    Type(VecType)   :: v0, r0, com  

    natoms=molecules_getnatoms(spc)

    !** Get the COM in the simcell using one of the biasing methods
    !** also get the biasfactor.
    Select Case (Trim(params%type_of_bias))
    Case("MAP_BIAS")

      com = rigidmoves_comFromBmap(params%insbias%bmap, params%simcell , &
          biasfactor )

      !** conflib gives back a molecule with zero com vel and position 
      !** at origin
      Call conflib_getrandomconfig(params%lib, &
          params%temp_ref_struc%coords(1:natoms,1)%rp,&
          params%temp_ref_struc%coords(1:natoms,1)%v, r0, v0, intra )

 !SDEBUG
      !** take image about z axis, assumes molecule is centered at origin
      If (params%lib%chiral_image) Then
          params%temp_ref_struc%coords(1:natoms,1)%rp%comp(3)= (-one)* &
              params%temp_ref_struc%coords(1:natoms,1)%rp%comp(3)
          params%temp_ref_struc%coords(1:natoms,1)%v%comp(3)= (-one)* &
              params%temp_ref_struc%coords(1:natoms,1)%v%comp(3)
      Endif
 !SDEBUG
      !** rotate the vels and positions, add com posn and velocity
      Call rigidmoves_getNewConfigOri(params%temp_ref_struc, spc, 1, com, v0)

      copy_accels=.False.
      Call config_copymolec(params%temp_ref_struc,1 ,species(spc), mol, &
          copy_accels)

    Case("CAVITY_BIAS","MAP_CAVITY_BIAS")

      totmol=config_getTotMoles(species)
      Call rigidmoves_ConfFromCavAndLib(params%insbias,species, &
          spc,mol,totmol, biasfactor, succ_flag, params%lib, &
          params%temp_ref_struc )

      If (.Not.succ_flag) Return

    Case Default
      Write(*,*) "This type of rigid bias insertion not yet coded fullly"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select

    !** apply periodic boundary conditions before forcefield evaluation
    Call simcell_pbc(params%simcell, &
        species(spc)%coords(1:natoms,mol)%rp, &
        species(spc)%coords(1:natoms,mol)%r, &
        species(spc)%coords(1:natoms,mol)%cr)

    !** Get the new energies
    fast = .True.
    skip_intra_flag = .False.  !** could be some intra interactions
    succ_flag = subinteract_int(subint,species,params%simcell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))
    ! Call checkandstop(succ_flag,__FILE__,__LINE__,'insert energy calc problem')

  End Subroutine rigidmoves_linsert

  !-----------------------------------------------------------------
  ! This routine handles the deletion of a molecule and gets the 
  ! bias-weight for deleting a molecule with no bias , this value is 
  ! actually one.
  ! Requires:  species -- the species coordinates
  !            spc -- the specific species number
  !            mol -- the molecule number
  !            params -- move parameters
  !            biasfactor -- the returned bias factor
  !            succ_flag -- success flag, false if energies too high
  !            subint -- subset interactions (molecule-system)
  !-----------------------------------------------------------------
  Subroutine rigidmoves_rdelete(params, species, spc, mol, biasfactor, &
      nrg_calc_flag,subint)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Type(RigidRDel_Params), Intent(In)             :: params
    Real(kind=RDbl), Intent(Out)                   :: biasfactor
    Logical, Intent(InOut)                         :: nrg_calc_flag    
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Logical                        :: fast,skip_intra_flag

    !** unbiased 
    biasfactor = one

    !** HACK
    fast = .True.
    !!    skip_intra_flag = .True.       ! intra nrgs will be zero  ????
    skip_intra_flag = .False.
    nrg_calc_flag = subinteract_int(subint,species,params%simcell,fast, &
        .False.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))

  End Subroutine rigidmoves_rdelete

  !-----------------------------------------------------------------
  ! This routine gets the bias-weight for deleting a molecule
  ! This is the opposite move for "rigidmoves_binsert"
  ! returns energy of existing molecule, biasfactor
  !           subint -- subset interactions (molecule-system)
  !-----------------------------------------------------------------
  Subroutine rigidmoves_bdelete(params, species, spc, mol, biasfactor, &
      nrg_calc_flag, subint)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: species
    Integer, Intent(in)                            :: spc,mol
    Type(RigidBDel_Params), Intent(in)             :: params
    Real(kind=RDbl), Intent(out)                   :: biasfactor
    Logical, Intent(inout)                         :: nrg_calc_flag    
    Type(Subset_Interactions), Intent(InOut)       :: subint    

    Integer                        :: ncubelets, totmol
    Real(kind=RDbl)                :: wt, pc_n
    Real(kind=RDbl), Dimension(3)  :: uccom
    Logical                        :: fast,skip_intra_flag

    biasfactor=one

    !** get map bias
    Select Case (Trim(params%type_of_bias))
    Case("CAVITY_BIAS")
      ! nothing to be done for now 
    Case ("MAP_BIAS","MAP_CAVITY_BIAS")

    !** Map the position to the unit cell. Also convert the 
    !** vector type to an array
! regular case
#ifndef BENT_GCMC
    uccom = simcell_maptouc(params%simcell,  &
        species(spc)%gcoords(mol)%rigidcoords%com)
#endif

! for molecule with bent shape
#ifdef BENT_GCMC
! For molecules with a  bent geometry. For now it works with the hack 
! that atom-1  is the anchor-atom, instead of COM. This part should be 
! integrated to regular code after more testing. (Tested for benzene)
    uccom = simcell_maptouc(params%simcell,  &
        species(spc)%coords(anchor_atom,mol)%rp)
#endif

    !** Get the weight of cubelet where the COM lies
    wt = bmap_getcellwt(params%insbias%bmap, uccom)

    !** Calculate the biasfactor
    ncubelets  = bmap_getncubelets(params%insbias%bmap)
    biasfactor = ncubelets*wt

    Case Default
      Write(*,*) "This type of rigid bias deletion not yet coded fullly"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select



    !** get cavity bias
    Select Case (Trim(params%type_of_bias))
    Case("MAP_BIAS")
      ! nothing to be done for now 
    Case ("CAVITY_BIAS","MAP_CAVITY_BIAS")
      ! modify biasfactor for cavity
      totmol=config_getTotMoles(species)
      pc_n= stats_getcavg(params%insbias%cav_prob(totmol)) / &
          params%insbias%cavNTrials
      biasfactor=biasfactor /pc_n

    Case Default
      ! error message already above
    End Select

    !** HACK
    fast = .True.
    !!    skip_intra_flag = .True.       ! intra nrgs will be zero  ????
    skip_intra_flag = .False.
    nrg_calc_flag = subinteract_int(subint,species,params%simcell,fast, &
        .False.,skip_intra_flag,(/params%max_nrg/),(/spc,mol,0/))

  End Subroutine rigidmoves_bdelete

  !--------------------------------------------------------------------------
  !** Some generic routines used by several rigid moves
  !--------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Perturb the rigid body rotational coordinates of the molecule by 
  ! specified increments.
  ! Requires:  species -- species data structure
  !            spc -- species number
  !            mol -- molecule number
  !            simcell -- simulation cell information
  !            dphi -- phi angle increment
  !            dpsi -- psi angle increment
  !            dcostheta -- cos(theta) angle increment
  !            newxyz -- flags creation of xyz coordinate from gen coords
  ! NOTE: This routine should be generalized to also handle rotations of
  ! the molecule's atomic velocities.  Note that no rigidmoves data type
  ! is used here, perhaps this belongs in config instead
  !-----------------------------------------------------------------------
  Subroutine rigidmoves_rotperturb(species,spc,mol,simcell, &
      dphi,dpsi,dcostheta,newxyz)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Type(SimCell_Params), Intent(In)               :: simcell
    Real(kind=RDbl), Intent(In)                    :: dphi,dpsi,dcostheta
    Logical, Intent(In)                            :: newxyz

    Integer                   :: natoms
    Real(kind=RDbl)           :: phi, psi, costheta
    Type(RigidMolec), Pointer :: rigidcoord

    natoms = molecules_getnatoms(spc)
    If (natoms == 1) Return

    !** Use a pointer to the generalized coordinates for clarity
    rigidcoord => species(spc)%gcoords(mol)%rigidcoords

    !** Sample phi and psi uniformly between -pi,pi but theta
    !** needs to be sampled from a cosine distribution
    !** See Allen and Tildesley, Pg. 132-133, "Molecular Liquids"
    phi = rigidcoord%phi + dphi
    psi = rigidcoord%psi + dpsi
    costheta = Cos(rigidcoord%theta) + dcostheta

    !** Map into a single range    
    rigidcoord%phi = phi - Anint(phi/twopi)*twopi
    rigidcoord%psi = psi - Anint(psi/twopi)*twopi
    
    !** This is dangerous near theta=zero, could lead to 
    !** 180 Deg rotation. - Shaji.
    costheta = costheta - Anint(costheta/2.0_RDbl)*2.0_RDbl
    rigidcoord%theta = Acos(costheta)

    !** Generate the new xyz coordinates if requested
    If (newxyz) Then
      Call rigidcoords_toxyz(rigidcoord,species(spc)%coords(1:natoms,mol)%rp)
  
      !** Apply PBCs to the generated xyz coordinates
      Call simcell_pbc(simcell, &
          species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms,mol)%r, &
          species(spc)%coords(1:natoms,mol)%cr)
    End If

  End Subroutine rigidmoves_rotperturb

  !-----------------------------------------------------------------------
  ! Perturb the rigid body COM coordinates of the molecule by given vector
  ! Requires:  species -- species data structure
  !            spc -- species number
  !            mol -- molecule number
  !            simcell -- simulation cell information
  !            dvec -- displacement vector
  !            newxyz -- flags creation of xyz coordinate from gen coords
  ! NOTE: Note that no rigidmoves data type
  ! is used here, perhaps this belongs in config instead
  !-----------------------------------------------------------------------
  Subroutine rigidmoves_transperturb(species,spc,mol,simcell,dvec,newxyz)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(VecType), Intent(In)                      :: dvec
    Logical, Intent(In)                            :: newxyz

    Integer         :: natoms

    !** Displace the COM coordinate
    species(spc)%gcoords(mol)%rigidcoords%com = &
        species(spc)%gcoords(mol)%rigidcoords%com + dvec

    !** Generate the new xyz coordinates if requested
    If (newxyz) Then
      natoms = molecules_getnatoms(spc)
      Call rigidcoords_toxyz(species(spc)%gcoords(mol)%rigidcoords, &
          species(spc)%coords(1:natoms,mol)%rp)
  
      !** Apply PBCs to the generated xyz coordinates
      Call simcell_pbc(simcell, &
          species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms,mol)%r, &
          species(spc)%coords(1:natoms,mol)%cr)
    End If

  End Subroutine rigidmoves_transperturb

  !-----------------------------------------------------------------------
  ! Generates new random orientations for a molecule, works on rigidcoords
  ! should be used only with rigidcoords having no internal flexibility
  ! Requires:  species -- species data structure for whole system
  !            spc -- species number
  !            mol -- molecule number
  !            newxyz -- flags creation of xyz coordinate from gen coords
  ! NOTE: calls to this routine should be replaced by generation of 
  ! phi,psi,costheta angles through the whole range, then calling 
  ! rigidmoves_rotperturb.
  !-----------------------------------------------------------------------
  Subroutine rigidmoves_reorient(species, spc, mol, newxyz )
    Type(AtMolCoords), Intent(InOut) :: species
    Integer, Intent(In)              :: spc, mol
    Logical, Intent(In)              :: newxyz

    Integer         :: natoms,i
    Real(kind=RDbl) :: costheta
#ifdef BENT_GCMC
! For molecules with a  bent geometry. For now it works with the hack 
! that atom-1  is the anchor-atom, instead of COM. This part should be 
! integrated to regular code after more testing. (Tested for benzene)
    Type(VecType) :: disp
#endif
    !** If natoms is = 1, we are done
    natoms = molecules_getnatoms(spc)

    If (natoms == 1) Then
      !** fill in rest of the rigidcoords
      species%gcoords(mol)%rigidcoords%theta = zero
      species%gcoords(mol)%rigidcoords%phi   = zero
      species%gcoords(mol)%rigidcoords%psi   = zero

      !** Generate the xyz coordinates
      If (newxyz) species%coords(1,mol)%rp = &
          species%gcoords(mol)%rigidcoords%com 

    Else 
      !** Generate the eulerian angles at random
      !** Sample theta from a cosine distribution.
      costheta = 1 - rranf() * 2.0
      species%gcoords(mol)%rigidcoords%theta = Acos(costheta)
      species%gcoords(mol)%rigidcoords%phi = pi - rranf()*twopi
      species%gcoords(mol)%rigidcoords%psi = pi - rranf()*twopi

      ! regular case
#ifndef BENT_GCMC
      !** Generate the xyz coordinates
      If (newxyz) Call rigidcoords_toxyz(species%gcoords(mol)%rigidcoords,&
          species%coords(1:natoms,mol)%rp)
#endif

      ! for molecule with bent shape
#ifdef BENT_GCMC
! For molecules with a  bent geometry. For now it works with the hack 
! that atom-1  is the anchor-atom, instead of COM. This part should be 
! integrated to regular code after more testing. (Tested for benzene)
      !** genreate xyz coords
      If (newxyz) Call rigidcoords_toxyz(species%gcoords(mol)%rigidcoords,&
          species%coords(1:natoms,mol)%rp)

      ! vector by which the whole molecules should be displaced so that 
      ! anchor_atom reachs the com position
      disp=species%gcoords(mol)%rigidcoords%com-&
          species%coords(anchor_atom,mol)%rp

      ! change com, rest of gcoords remain same
      species%gcoords(mol)%rigidcoords%com= &
          species%gcoords(mol)%rigidcoords%com+disp

      ! change rp, rest (like r) will be changed by calling program
      Do i=1,natoms
        species%coords(i,mol)%rp= species%coords(i,mol)%rp+disp
      End Do
#endif
    End If

  End Subroutine rigidmoves_reorient

  !----------------------------------------------------------------------
  ! Generates new random orientations for a molecule, works on %rp
  ! This should be used while working with rigidcoords with 
  ! internal flexibility. It works only when tempspecies contains the molecule
  ! whose com coincides with origin. Bit of a **HACK***. 
  ! Requires:  tempspecies -- a atmolcoords structure holding one molecule 
  !            spc -- the specific species number
  !            mol -- index of molecule in tempspecies, (will be 1 always)
  !            com -- com of the new position
  !            v0 -- new velocity of center of mass
  !----------------------------------------------------------------------
  Subroutine rigidmoves_getNewConfigOri(tempspecies, spc, mol, com, v0)
    Type(AtMolCoords), Intent(InOut) :: tempspecies
    Type(VecType), Intent(in)        :: com,v0
    Integer, Intent(in)              :: spc, mol

    Real(kind=RDbl) :: phi, theta, psi
    Integer         :: natoms,a
    Logical         :: rotateVels
#ifdef BENT_GCMC
    ! For molecules with a  bent geometry. For now it works with the hack 
    ! that atom-1  is the anchor-atom, instead of COM. This part should be 
    ! integrated to regular code after more testing. (Tested for dm-cyclobutane)
    Type(VecType) :: disp
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) "bent gcmc and linsert move is a deadly combination -shaji"
    Write(*,*) "comment this out to continue"
    Stop
#endif
    natoms=molecules_getnatoms(spc)

    !** Generate the eulerian angles at random
    !** Sample theta from a cosine distribution.
    theta = Acos( 1 - (rranf() * 2.0) )
    phi = pi - rranf()*twopi
    psi = pi - rranf()*twopi

    !** not useful, still storing anyway
    Call rigidcoords_changeref(tempspecies%gcoords(mol)%rigidcoords, &
        tempspecies%coords(1:natoms,mol)%rp )

    tempspecies%gcoords(mol)%rigidcoords%theta = theta
    tempspecies%gcoords(mol)%rigidcoords%phi = phi 
    tempspecies%gcoords(mol)%rigidcoords%psi = psi
    tempspecies%gcoords(mol)%rigidcoords%com = com

    !** If natoms is = 1, we are done
    natoms = molecules_getnatoms(spc)
    If (natoms == 1) Then
      !      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      !      Stop
      !** Generate the xyz coordinates
      tempspecies%coords(1,mol)%rp = zero 
    Else
      ! ** rotate %rp, %v around origin
      ! ** Note : this molecule has its com coinciding with origin
      rotateVels=.True.
      Call config_rotatexyz(tempspecies, mol, theta,phi,psi,rotateVels)
    End If
#ifndef BENT_GCMC
    ! add the com_vel obtained from library and com_posn from bmap
    Do a=1,natoms
      tempspecies%coords(a,mol)%v = tempspecies%coords(a,mol)%v + v0
      tempspecies%coords(a,mol)%rp = tempspecies%coords(a,mol)%rp + &
          tempspecies%gcoords(mol)%rigidcoords%com 
    End Do
#endif

#ifdef BENT_GCMC
    ! For molecules with a  bent geometry. For now it works with the hack 
    ! that atom-1  is the anchor-atom, instead of COM. This part should be 
    ! integrated to regular code after more testing. (Tested for benzene)

    ! vector by which the whole molecules should be displaced so that 
    ! anchor_atom reachs the com position
    disp=(-one)*tempspecies%coords(anchor_atom,mol)%rp

 ! change com, rest of gcoords remain same
      tempspecies%gcoords(mol)%rigidcoords%com= &
          tempspecies%gcoords(mol)%rigidcoords%com+disp

    Do a=1,natoms
      tempspecies%coords(a,mol)%v = tempspecies%coords(a,mol)%v + v0
      tempspecies%coords(a,mol)%rp = tempspecies%coords(a,mol)%rp + &
          tempspecies%gcoords(mol)%rigidcoords%com 
    End Do

#endif

  End Subroutine rigidmoves_getNewConfigOri

  !----------------------------------------------------------------------
  ! Generates a center-of-mass of the molecule according to a
  ! bias map , return the biasfactor
  ! Requires: b_map -- bias map 
  !           scell - simcell details
  !           biasfactor -- the returned bias factor
  !----------------------------------------------------------------------
  Type(VecType) Function rigidmoves_comFromBmap(b_map, scell ,biasfactor)
    Real(kind=RDbl), Intent(Out)                     :: biasfactor
    Type(BiasMap_Params), Pointer                    :: b_map
    Type(SimCell_Params), Pointer                    :: scell

    Integer              :: icellx,icelly,icellz, ncubelets, index
    Type(VecType)        :: com
    Real(kind=RDbl)      :: wt

    !** Get the position of the center of mass using the bias map
    !** 1.) Pick a point in the unit cell according to the biasmap
    index = bmap_getbiasindex(b_map)
    If ((index<1).Or.(index>b_map%numtab)) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "some error",index
    End If

    ncubelets  = bmap_getncubelets(b_map)
    wt = bmap_getbiaswt(b_map,index)
    biasfactor = 1.0_RDBl/wt/ncubelets
    com = bmap_getbiaspt(b_map, index)

    !** 2.) Select one of the unit cells in the simulation cell
    icellx = Int( rranf() * scell%nx )
    icelly = Int( rranf() * scell%ny )
    icellz = Int( rranf() * scell%nz )

    !** 3.) Translate the pt. to the appropriate unit cell
    com = com + simcell_getcellorigin(scell, icellx, icelly, icellz)

    rigidmoves_comFromBmap=com
  End Function rigidmoves_comFromBmap


  !----------------------------------------------------------------------
  ! Generates the coordinates of a molecule from cavity bias 
  ! , returns the biasfactor
  ! Requires: params -- insertion parameters
  !           species -- the coordinates
  !           spc, mol -- on which the insert is to performed
  !           totmol -- current total occupancy ( sum of nmoles of non-fixed)
  !           biasfactor -- the returned bias factor
  !           cavSuccFlag -- tells whether at least one cavity was found or not
  ! Note : 
  !      1.) usually a check for lying in a cavity is done Only for the 
  ! center of mass. But we generate a random configuration and check 
  ! whether all atoms lie in a cavity
  !      2.) works ONLY for rigidcoords with no internal_flexibility
  !----------------------------------------------------------------------
  Subroutine rigidmoves_getConfigFromCavBias(insbias, species, spc, &
      mol, totmol, biasfactor,  cavSuccFlag)
    Type(RigidInsBias), Intent(InOut)                     :: insbias
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: species
    Integer, Intent(in)                              :: spc, mol, totmol
    Real(kind=RDbl), Intent(out)                     :: biasfactor
    Logical, Intent(out)                            :: cavSuccFlag

    Integer :: n_success, n_trials, i, a, ix, iy, iz, trial_index, ncubelets
    Integer :: index, error, natoms
    Real(kind=RDbl) :: pc_n
    Logical :: is_cavity, newxyz

    Type(VecType),Dimension(:,:), Pointer, Save      :: temp_vec_array
    Type(VecType),Dimension(:), Pointer, Save        :: com_array
    Real(kind=RDbl), Dimension(:), Pointer, Save :: map_bias, phi, theta, psi 
    Logical :: only_cavity
    Integer, Save :: max_natoms
    Logical, Save :: first_time=.True.

    n_trials=insbias%cavNTrials
    natoms=molecules_getnatoms(spc)

    !** Allocate these storage arrays
    If (first_time) Then
      first_time=.False.
      max_natoms=natoms
      Allocate(map_bias(n_trials), phi(n_trials), theta(n_trials), & 
          psi(n_trials), com_array(n_trials), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Allocate(temp_vec_array(n_trials,max_natoms), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Endif

    If (natoms>max_natoms) Then
      max_natoms=natoms
      DeAllocate(temp_vec_array, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Allocate(temp_vec_array(n_trials,max_natoms), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Endif


    ! decide type of bias for com
    If (Associated(insbias%cavity)) Then
      If (Associated(insbias%bmap)) Then
        only_cavity=.False.
      Else
        only_cavity=.True.
      Endif
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    n_success=0
    Do i=1,n_trials
      
      If (only_cavity) Then
        !** Get the position of the center of mass unifromly from the 
        !** simulation volume
        species(spc)%gcoords(mol)%rigidcoords%com= &
            simcell_uniformRandPt(insbias%simcell)
        biasfactor=one
      Else
        species(spc)%gcoords(mol)%rigidcoords%com = rigidmoves_comFromBmap(&
            insbias%bmap, insbias%simcell , biasfactor )
      Endif

      newxyz = .True.

      Call rigidmoves_reorient(species(spc),spc,mol,newxyz)

      !** make sure all atoms are within one unit cell, we will be 
      !** using %r While checking for cavities
      Call simcell_pbc(insbias%simcell, &
          species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms,mol)%r, &
          species(spc)%coords(1:natoms,mol)%cr)

      is_cavity=.True.
      Do a=1,natoms
        !** 5. Check whether com lies in  a cavity or not
        is_cavity = cavitylist_checkcavity(insbias%cavity, &
            species(spc)%coords(a,mol)%r)
        If (.Not.is_cavity) Exit
      End Do

      If (is_cavity) Then
        ! store indices for this com
        n_success=n_success+1
        map_bias(n_success) = biasfactor
        temp_vec_array(n_success,1:natoms) = &
            species(spc)%coords(1:natoms,mol)%rp
        phi(n_success) = species(spc)%gcoords(mol)%rigidcoords%phi
        theta(n_success) = species(spc)%gcoords(mol)%rigidcoords%theta
        psi(n_success) = species(spc)%gcoords(mol)%rigidcoords%psi
        com_array(n_success) = species(spc)%gcoords(mol)%rigidcoords%com
      Endif
    End Do

    !** Update the running average of number of cavities found at 
    !** occupancy=totmol-1 ( i.e. while trying for inserting totmol'th mol)
    Call stats_update(insbias%cav_prob(totmol), n_success*one)

    !** Now we have recorded all the successes
    !** If there is no success, then Mezei (1980) says do usual insertion
    !** But I will just an error flag telling that no low energy sites, 
    !** so insertion failed
    If (n_success==0) Then
      cavSuccFlag=.False.
      Return
    Else
      cavSuccFlag=.True.

      !** Choose one of the successes randomly
      trial_index=Int(rranf()*n_success)+1
      If (trial_index<=n_success) Then
        biasfactor=map_bias(trial_index)
        species(spc)%coords(1:natoms,mol)%rp= &
            temp_vec_array(trial_index,1:natoms)         
        species(spc)%gcoords(mol)%rigidcoords%phi = phi(trial_index)
        species(spc)%gcoords(mol)%rigidcoords%theta = theta(trial_index)
        species(spc)%gcoords(mol)%rigidcoords%psi = psi(trial_index)
        species(spc)%gcoords(mol)%rigidcoords%com = com_array(trial_index)
        pc_n=stats_getcavg(insbias%cav_prob(totmol)) / n_trials

        !** bias factor = dV_c / V  * 1/wt * Pc(N)
        !** biasfactor returned from _comfrombmap = dV_c / V  * 1/wt 
        biasfactor = biasfactor * pc_n

      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

    Endif

  End Subroutine rigidmoves_getConfigFromCavBias

  !----------------------------------------------------------------------
  ! Generates the coordinates of a molecule from cavity bias and
  ! returns the biasfactor, uses a library
  ! Requires: params -- insertion parameters
  !           species -- the coordinates
  !           spc, mol -- on which the insert is to performed
  !           totmol -- current total occupancy ( sum of nmoles of non-fixed)
  !           biasfactor -- the returned bias factor
  !           cavSuccFlag -- tells whether at least one cavity was found or not
  !           lib -- the library object from which the trial config is chosen 
  !           tempspc -- a temperory memory used for holding trial molecule.
  ! Note : 
  !      1.) usually a check for lying in a cavity is done Only for the 
  ! center of mass. But we generate a random configuration and check 
  ! whether all atoms lie in a cavity
  !      2.) works ONLY for rigidcoords WITH internal_flexibility
  !----------------------------------------------------------------------
  Subroutine rigidmoves_ConfFromCavAndLib(insbias, species, spc, &
      mol, totmol, biasfactor,  cavSuccFlag, lib, tempspc)
    Type(RigidInsBias), Intent(InOut)                     :: insbias
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: species
    Integer, Intent(in)                              :: spc, mol, totmol
    Real(kind=RDbl), Intent(out)                     :: biasfactor
    Logical, Intent(out)                             :: cavSuccFlag
    Type(ConfLibrary), Pointer                       :: lib
    Type(AtMolCoords),  Intent(inout)                :: tempspc


    Integer :: n_success, n_trials, i, a, ix, iy, iz, trial_index, ncubelets
    Integer :: index, error, natoms
    Real(kind=RDbl) :: pc_n
    Real(kind=Rdbl) :: intra 
    Type(VecType)   :: v0, r0, com  
    Logical :: is_cavity, newxyz, rotatevels, copy_accels

    Type(VecType),Dimension(:,:), Pointer, Save      :: temp_pos_array
    Type(VecType),Dimension(:,:), Pointer, Save      :: temp_vel_array
    Type(VecType),Dimension(:), Pointer, Save        :: com_array
    Real(kind=RDbl), Dimension(:), Pointer, Save :: map_bias, phi, theta, psi 
    Integer, Save :: max_natoms
    Logical, Save :: first_time=.True.
    Logical       :: only_cavity 

    n_trials=insbias%cavNTrials
    natoms=molecules_getnatoms(spc)

    !** Allocate these storage arrays
    If (first_time) Then
      first_time=.False.
      max_natoms=natoms
      Allocate(map_bias(n_trials), phi(n_trials), theta(n_trials), & 
          psi(n_trials), com_array(n_trials), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Allocate(temp_pos_array(n_trials,max_natoms), &
          temp_vel_array(n_trials,max_natoms),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Endif

    !** this routine might be visited by more than one specie-type
    If (natoms>max_natoms) Then
      max_natoms=natoms
      DeAllocate(temp_pos_array, temp_vel_array, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Allocate(temp_pos_array(n_trials,max_natoms), &
          temp_vel_array(n_trials,max_natoms),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Endif

    ! decide type of bias for com
    If (Associated(insbias%cavity)) Then
      If (Associated(insbias%bmap)) Then
        only_cavity=.False.
      Else
        only_cavity=.True.
      Endif
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    n_success=0

    Do i=1,n_trials

      If (only_cavity) Then
        !** Get the position of the center of mass unifromly from the 
        !** simulation volume
        com = simcell_uniformRandPt(insbias%simcell)
        species(spc)%gcoords(mol)%rigidcoords%com=com
        biasfactor=one
      Else
        !** only cavity bias
        com = & 
            rigidmoves_comFromBmap(insbias%bmap, insbias%simcell , biasfactor)
        species(spc)%gcoords(mol)%rigidcoords%com=com
      Endif

      newxyz=.True.

      !** conflib gives back a molecule with zero com vel and position 
      !** at origin
      Call conflib_getrandomconfig(lib, tempspc%coords(1:natoms,1)%rp, &
          tempspc%coords(1:natoms,1)%v, r0, v0, intra )

 !SDEBUG
      !** take image about z axis, assumes molecule is centered at origin
      If (lib%chiral_image) Then
          tempspc%coords(1:natoms,1)%rp%comp(3)= (-one)* &
              tempspc%coords(1:natoms,1)%rp%comp(3)
          tempspc%coords(1:natoms,1)%v%comp(3)= (-one)* &
              tempspc%coords(1:natoms,1)%v%comp(3)
      Endif
 !SDEBUG

      !** rotate the vels and positions, add com posn and velocity
      ! if BENT_GCMC compiler option is used com will actually be shifted
      Call rigidmoves_getNewConfigOri(tempspc, spc, 1, com, v0)

      !** make sure all atoms are within one unit cell, we will be 
      !** using %r While checking for cavities
      Call simcell_pbc(insbias%simcell, &
          tempspc%coords(1:natoms,1)%rp, tempspc%coords(1:natoms,1)%r, &
          tempspc%coords(1:natoms,1)%cr)

      is_cavity=.True.
      Do a=1,natoms
        !** 5. Check whether com lies in  a cavity or not
        is_cavity = cavitylist_checkcavity(insbias%cavity, &
            tempspc%coords(a,1)%r)
        If (.Not.is_cavity) Exit
      End Do

      If (is_cavity) Then
        ! store indices for this com
        n_success=n_success+1
        map_bias(n_success) = biasfactor
        temp_vel_array(n_success,1:natoms) = tempspc%coords(1:natoms,1)%v
        temp_pos_array(n_success,1:natoms) = tempspc%coords(1:natoms,1)%rp
        phi(n_success) = tempspc%gcoords(1)%rigidcoords%phi
        theta(n_success) = tempspc%gcoords(1)%rigidcoords%theta
        psi(n_success) = tempspc%gcoords(1)%rigidcoords%psi
        com_array(n_success) = com
      Endif
    End Do

    !** Update the running average of number of cavities found at 
    !** occupancy=totmol-1 ( i.e. while trying for inserting totmol'th mol)
    Call stats_update(insbias%cav_prob(totmol), n_success*one)

    !** Now we have recorded all the successes
    !** If there is no success, then Mezei (1980) says do usual insertion
    !** But I will just an error flag telling that no low energy sites, 
    !** so insertion failed
    If (n_success==0) Then
      cavSuccFlag=.False.
      Return
    Else
      cavSuccFlag=.True.

      !** Choose one of the successes randomly
      trial_index=Int(rranf()*n_success)+1
      If (trial_index<=n_success) Then
        biasfactor=map_bias(trial_index)
        species(spc)%coords(1:natoms,mol)%rp= &
            temp_pos_array(trial_index,1:natoms)         
        species(spc)%coords(1:natoms,mol)%v= &
            temp_vel_array(trial_index,1:natoms)         
        species(spc)%gcoords(mol)%rigidcoords%phi = phi(trial_index)
        species(spc)%gcoords(mol)%rigidcoords%theta = theta(trial_index)
        species(spc)%gcoords(mol)%rigidcoords%psi = psi(trial_index)
        species(spc)%gcoords(mol)%rigidcoords%com = com_array(trial_index)
        pc_n=stats_getcavg(insbias%cav_prob(totmol)) / n_trials

        !** bias factor = dV_c / V  * 1/wt * Pc(N)
        !** biasfactor returned from _comfrombmap = dV_c / V  * 1/wt 
        biasfactor = biasfactor * pc_n

      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

    Endif

  End Subroutine rigidmoves_ConfFromCavAndLib


  !---------------- Move parameter display routines ---------------------------

  !-------------------------------------------------------------
  ! display params for random-translation moves
  !-------------------------------------------------------------
  Subroutine rigidmoves_displayRTrans(params,indent,unitno)
    Type(RigidRTrans_Params), Intent(In) :: params
    Integer, Intent(In)                  :: indent,unitno

    Character(len=indent)       :: blank
    Character(len=strLen)       :: string

    blank = Repeat(' ',indent)

    Write(unitno,'(2a)') blank,"Rigid Molecule Random-Translation Parameters :"
    If (params%constrained) Then
      Write(unitno,'(2a)') blank,"CONSTRAINED by matrix : "
      Call matrix_display(params%constrtrans,'f8.3',unitno)
      Write(unitno,'(2a)') blank,"The vector normal to constrained plain is: "
      Call vector_filedisplay(params%normal,unitno)
    End If

    If (params%scale_trans) Then
      Write(unitno,'(2a)') blank,"Scaling of Translation length allowed."
    Else
      Write(unitno,'(2a)') blank,"Scaling of Translation length not allowed."
    End If
    Write(unitno,'(2a,f12.5)') blank,"Initial translation length: ",&
        params%deltatrans

    string = real2str(params%max_nrg,8)
    Write(unitno,'(4a)') blank,"Maximum allowed energy: ",&
        Trim(string)," kcal/mol"

  End Subroutine rigidmoves_displayRTrans

  !-------------------------------------------------------------
  ! display params for Random-Rotation moves
  !-------------------------------------------------------------
  Subroutine rigidmoves_displayRRot(params,indent,unitno)
    Type(RigidRRot_Params), Intent(In) :: params
    Integer, Intent(In)                :: indent,unitno

    Character(len=indent)       :: blank
    Character(len=strLen)       :: string

    blank = Repeat(' ',indent)

    Write(unitno,'(2a)') blank,"Rigid Molecule Random-Rotation Parameters :"
    If (params%scale_rot) Then
      Write(unitno,'(2a)') blank,"Scaling of Rotation increment allowed."
    Else
      Write(unitno,'(2a)') blank,"Scaling of Rotation increment not allowed."
    End If
    Write(unitno,'(2a,f12.5)') blank,"Initial rotation-delta: ",&
        params%deltarot

    string = real2str(params%max_nrg,8)
    Write(unitno,'(4a)') blank,"Maximum allowed energy: ",&
        Trim(string)," kcal/mol"

  End Subroutine rigidmoves_displayRRot

  !-------------------------------------------------------------
  ! Display params for Reinsertion moves
  ! Requires:  params -- move parameters
  !            indent -- indent
  !            unitno -- unit number to dump into
  !-------------------------------------------------------------
  Subroutine rigidmoves_displayReinsert(params,indent,unitno)
    Type(RigidReinsert_Params), Intent(In)       :: params
    Integer, Intent(In)                          :: indent,unitno

    Character(len=indent)       :: blank
    Character(len=strLen)       :: string

    blank = Repeat(' ',indent)

    Write(unitno,'(2a)') blank,"Rigid Molecule Reinsertion Parameters :"

    string = real2str(params%max_nrg,8)
    Write(unitno,'(4a)') blank,"Maximum allowed energy: ",&
        Trim(string)," kcal/mol"
    Write(unitno,'(a)') blank,"No parameters for this move"

  End Subroutine rigidmoves_displayReinsert

  !-------------------------------------------------------------
  ! Display params for Biased Rotation/Translation move type
  ! Requires:  params -- move parameters
  !            indent -- indent
  !            unitno -- unit number to dump into
  !-------------------------------------------------------------
  Subroutine rigidmoves_displayBiasRT(params,indent,unitno)
    Type(RigidbiasRT_Params), Intent(In)        :: params
    Integer,Intent(In)                          :: indent,unitno

    Character(len=indent)       :: blank
    Character(len=strLen)       :: string

    blank = Repeat(' ',indent)

    Write(unitno,'(2a)') blank,"Rigid Molecule Biased Rotation/Translation:"
    string = int2str(params%ntrials)
    Write(unitno,'(1x,3a)') blank,'Number of trial moves: ',Trim(string)

    If (params%scale_rot) Then
      Write(unitno,'(1x,2a)') blank,"Scaling of Rotation increment allowed."
    Else
      Write(unitno,'(1x,2a)') blank,"Scaling of Rotation increment not allowed."
    End If
    string = real2str(params%deltarot,5)
    Write(unitno,'(1x,3a)') blank,"Initial rotation-delta: ",Trim(string)

    If (params%scale_trans) Then
      Write(unitno,'(1x,2a)') blank,"Scaling of Translation length allowed."
    Else
      Write(unitno,'(1x,2a)') blank,"Scaling of Translation length not allowed."
    End If
    string = real2str(params%deltatrans,5)
    Write(unitno,'(1x,3a)') blank,"Initial translation length : ",Trim(string)

    string = real2str(params%max_nrg,8)
    Write(unitno,'(4a)') blank,"Maximum allowed energy: ",&
        Trim(string)," kcal/mol"

  End Subroutine rigidmoves_displayBiasRT

  !-------------------------------------------------------------
  ! display params for biased insertion moves from a library of configs
  !-------------------------------------------------------------
  Subroutine rigidmoves_displayLIns(params,indent,unitno)
    Type(RigidLIns_Params), Intent(In)   :: params
    Integer, Intent(In)                  :: indent,unitno

    Character(len=strLen)       :: string
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    Write(unitno,'(2a)') blank,"Rigid  Biased-Library-Insertion Parameters :"
    Call rigidmoves_displayInsBias(params%insbias, indent+2, unitno)
    Write(unitno,'(3a)') blank,"Config_library File: ", &
        Trim(params%confFileName)

    string = real2str(params%max_nrg,8)
    Write(unitno,'(4a)') blank,"Maximum allowed energy: ",&
        Trim(string)," kcal/mol"

  End Subroutine rigidmoves_displayLIns

  !-------------------------------------------------------------
  ! display params for biased insertion moves
  !-------------------------------------------------------------
  Subroutine rigidmoves_displayBIns(params,indent,unitno)
    Type(RigidBIns_Params), Intent(In)   :: params
    Integer, Intent(In)                  :: indent,unitno

    Character(len=strLen)       :: string
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    Write(unitno,'(2a)') blank,"Rigid Molecule Biased-Insertion Parameters :"
    Call rigidmoves_displayInsBias(params%insbias, indent+2, unitno)

  End Subroutine rigidmoves_displayBIns

  !-------------------------------------------------------------
  ! Display params for deletion moves
  !-------------------------------------------------------------
  Subroutine rigidmoves_displayBDel(params,indent,unitno)
    Type(RigidBDel_Params), Intent(In) :: params
    Integer, Intent(In)                :: indent,unitno

    Character(len=strLen)       :: string
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    Write(unitno,'(2a)') blank,"Rigid Molecule Deletion Parameters :"
    Call rigidmoves_displayInsBias(params%insbias, indent+2, unitno)

  End Subroutine rigidmoves_displayBDel

  !-------------------------------------------------------------
  ! Display params insertion bias; used by lins, bins, bdel 
  !-------------------------------------------------------------
  Subroutine rigidmoves_displayInsBias(insbias,indent,unitno)
    Type(RigidInsBias), Intent(In)  :: insbias
    Integer, Intent(In)             :: indent,unitno

    Character(len=strLen)         :: string
    Character(len=indent)         :: blank

    blank = Repeat(' ',indent)

    Write(unitno,'(3a)') blank,"Type of bias: ",Trim(insbias%type_of_bias)
    Write(unitno,'(2a)') blank,"Details of Bias for insertion: "

    If (Associated(insbias%bmap)) Then
      Write(unitno,'(3a)') blank,"Bias map file: "//&
          Trim(insbias%biasfilename)
      Write(unitno,'(2a,f8.3)') blank,"Bias Temperature, K: ",insbias%bias_temp
      Write(unitno,'(2a)') blank,"Bias map details: "
      Call bmap_display(insbias%bmap,.False.,indent+2,unitno)
    End If

    If (Associated(insbias%cavity)) Then
      Write(unitno,'(2a,i5)') blank,"Number of trials during cavity: ", &
          insbias%cavNtrials
    End If

    string = real2str(insbias%max_nrg,8)
    Write(unitno,'(4a)') blank,"Maximum allowed energy: ", &
        Trim(string)," kcal/mol"
    
  End Subroutine rigidmoves_displayInsBias


  !---------------- Dynamic display routines ------------------------------

  !----------------------------------------------------------------------
  ! Single line display information for dynamic Translation parameters
  ! Requires: params -- random translation parameters
  !----------------------------------------------------------------------
  Function rigidmoves_dyninfoRTrans(params)
    Character(len=strLen)               :: rigidmoves_dyninfoRTrans
    Type(RigidRTrans_Params),Intent(In) :: params

    Write(rigidmoves_dyninfoRTrans,'(a,f6.3)') 'Max disp: ',&
        params%deltatrans

  End Function rigidmoves_dyninfoRTrans

  !----------------------------------------------------------------------
  ! Single line display information for dynamic Rotation parameters
  ! Requires: params -- random rotation parameters
  !----------------------------------------------------------------------
  Function rigidmoves_dyninfoRRot(params)
    Character(len=strLen)               :: rigidmoves_dyninfoRRot
    Type(RigidRRot_Params),Intent(In)   :: params

    Write(rigidmoves_dyninfoRRot,'(a,f6.3)') 'Max disp: ',&
        params%deltarot

  End Function rigidmoves_dyninfoRRot

  !-------------------------------------------------------------------------
  ! Single line display information for dynamic Biased rotation/translation 
  ! parameters.
  ! Requires: params -- random rotation parameters
  !----------------------------------------------------------------------
  Function rigidmoves_dyninfoBiasRT(params)
    Character(len=strLen)                :: rigidmoves_dyninfoBiasRT
    Type(RigidBiasRT_Params), Intent(In) :: params

    Character(len=strLen)               :: string1,string2

    string1 = real2str(params%deltatrans,4)
    string2 = real2str(params%deltarot,4)
    Write(rigidmoves_dyninfoBiasRT,'(2a,1x,a)') 'Max dx,rot: ',Trim(string1), &
        Trim(string2)

  End Function rigidmoves_dyninfoBiasRT

  !---------------- Parameter adjustment routines ----------------------------

  !-------------------------------------------------------------------------
  ! Adjust the jump length to try and keep the acceptance ratio close to 50%
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_adjustdeltatrans(params, ratio)
    Type(RigidRTrans_Params),Intent(inout) :: params
    Real(kind=RDbl), Intent(in)       :: ratio

    If (.Not. params%scale_trans) Return

    !** these values should not be hard coded here, should go into control file
    If (ratio < 0.49) Then
      params%deltatrans = Max(params%deltatrans*0.95_RDbl, MIN_TRANS_DELTA)
    Else If (ratio > 0.51) Then
      params%deltatrans = Min(params%deltatrans*1.05_RDbl, MAX_TRANS_DELTA)
    End If

  End Subroutine rigidmoves_adjustdeltatrans

  !-------------------------------------------------------------------------
  ! Adjust the jump length to try and keep the acceptance ratio close to 50%
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_adjustdeltarot(params, ratio)
    Type(RigidRRot_Params),Intent(inout) :: params
    Real(kind=RDbl), Intent(in)       :: ratio

    If (.Not. params%scale_rot) Return

    !** these values should not be hard coded here, should go to control file
    If (ratio < 0.49) Then
      params%deltarot = Max(params%deltarot*0.95_RDbl, MIN_TRANS_DELTA)
    Else If (ratio > 0.51) Then                                        
      params%deltarot = Min(params%deltarot*1.05_RDbl, MAX_TRANS_DELTA)
    End If

  End Subroutine rigidmoves_adjustdeltarot

  !-------------------------------------------------------------------------
  ! Adjust the maximum rotation and translation to try to keep the 
  ! acceptance ratio close to 50%
  ! Requires:  params -- move parameters
  !            ratio -- current acceptance ratio
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_adjustBiasRT(params,ratio)
    Type(RigidBiasRT_Params), Intent(InOut) :: params
    Real(kind=RDbl), Intent(In)             :: ratio

    If (.Not. params%scale_rot) Return

    !** Adjust rotation maximum if necessary
    If (params%scale_rot) Then
      If (ratio < 0.49) Then
        params%deltarot = Max(params%deltarot*0.95_RDbl, MIN_TRANS_DELTA)
      Else If (ratio > 0.51) Then                                        
        params%deltarot = Min(params%deltarot*1.05_RDbl, MAX_TRANS_DELTA)
      End If
    End If

    !** Adjust translation maximum if necessary
    If (params%scale_trans) Then
      If (ratio < 0.49) Then
        params%deltatrans = Max(params%deltatrans*0.95_RDbl, MIN_TRANS_DELTA)
      Else If (ratio > 0.51) Then                                        
        params%deltatrans = Min(params%deltatrans*1.05_RDbl, MAX_TRANS_DELTA)
      End If
    End If

  End Subroutine rigidmoves_adjustBiasRT

  !-------------------------------------------------------------------------
  ! Check the binsert and bdel moves are properly initialized 
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_CheckBInsDel(insert_par, delete_par, species, &
      auxmv, ctrlfilename)
    Type(RigidBIns_Params),Pointer                     :: insert_par
    Type(RigidBDel_Params),Pointer                     :: delete_par
    Type(AtMolCoords), Dimension(:), Intent(inout)     :: species 
    Type(AuxMoveObjects),Pointer           :: auxmv
    Character(len=strLen), Intent(in)                  :: ctrlfilename

    Type(RigidInsBias),Pointer :: inspar, delpar
    Integer :: error

    delete_par%type_of_bias=insert_par%type_of_bias

    inspar=>insert_par%insbias
    Allocate(delete_par%insbias,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    delpar=>delete_par%insbias

    Call rigidmoves_equateInsBias(inspar, delpar,species, auxmv, ctrlfilename)

  End Subroutine rigidmoves_CheckBInsDel

  !-------------------------------------------------------------------------
  ! Check the linsert and bdel moves are properly initialized 
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_CheckLInsDel(insert_par, delete_par, species, &
      auxmv, ctrlfilename)
    Type(RigidLIns_Params),Pointer                     :: insert_par
    Type(RigidBDel_Params),Pointer                     :: delete_par
    Type(AtMolCoords), Dimension(:), Intent(inout)     :: species 
    Type(AuxMoveObjects),Pointer           :: auxmv 
    Character(len=strLen), Intent(in)                  :: ctrlfilename

    Type(RigidInsBias),Pointer :: inspar, delpar
    Integer :: error

    delete_par%type_of_bias=insert_par%type_of_bias

    inspar=>insert_par%insbias
    Allocate(delete_par%insbias,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    delpar=>delete_par%insbias

    Call rigidmoves_equateInsBias(inspar, delpar, species, auxmv, ctrlfilename)

  End Subroutine rigidmoves_CheckLInsDel

  !----------------------------------------------------------
  ! equates both inspar and delpar, initializes cavity-params
  !----------------------------------------------------------
  Subroutine rigidmoves_equateInsBias(inspar, delpar, species, auxmv, ctrlfile)
    Type(RigidInsBias),Pointer     :: inspar, delpar
    Type(AtMolCoords), Dimension(:), Intent(inout)     :: species 
    Type(AuxMoveObjects), Pointer           :: auxmv 
    Character(len=strLen), Intent(in)                  :: ctrlfile

    Logical :: init_cavity
    Integer :: i,error, newsize

    ! ** nullify all deletion parameters
    Nullify(delpar%bmap,delpar%cavity)
    delpar%max_nrg=inspar%max_nrg

    If (Trim(inspar%biasfilename)=="NULL") Then
      delpar%biasfilename=inspar%biasfilename
    Else
      delpar%biasfilename=inspar%biasfilename
      delpar%bmap=>inspar%bmap
      delpar%bias_temp=inspar%bias_temp
    Endif
    
    init_cavity=.False.
    If ((Trim(inspar%type_of_bias)=="MAP_CAVITY_BIAS").Or.&
        (Trim(inspar%type_of_bias)=="CAVITY_BIAS"))  init_cavity=.True.

    If (init_cavity) Then
      Write(*,*) "Initializing cavity"
      delpar%type_of_bias=inspar%type_of_bias
      delpar%cavNTrials =inspar%cavNTrials

      Nullify(delpar%cavity)

!!$
!!$      Call cavitylist_init(inspar%cavity, species, inspar%simcell, &
!!$          ctrlfile )
!!$
!!$      delpar%cavity=>inspar%cavity

      delpar%cavity=>auxmv%cavity
      inspar%cavity=>auxmv%cavity
      

      !** Allocate the array that stores the probability of finding 
      !** a cavity at occupation N
      newsize=Max(CAVBIAS_MAXMOLEC, config_getTotMoles(species))
      Allocate(inspar%cav_prob(newsize),STAT=error)
      Do i=1,Size(inspar%cav_prob)

        !** - STATS_BLOCKSIZE defined in defaults
        Call stats_init(inspar%cav_prob(i),"Prob Cav"//Trim(int2str(i)), &
            STATS_BLOCKSIZE,.False.)

        !** - We assume that we will be able to find 50% cavities
        !**  at the beginning of the simulation
        !** This is just to avoid divison by zero errors later
        Call stats_update(inspar%cav_prob(i),inspar%cavNTrials/2.00_Rdbl)

      End Do

      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      delpar%cav_prob=> inspar%cav_prob

    Else

      delpar%type_of_bias=inspar%type_of_bias

    Endif

  End Subroutine rigidmoves_equateInsBias

  !-------------------------------------------------------------------------
  ! update cavity information after rtranslate 
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_postadjustTrans(transpar, species, spc, mol)
    Type(RigidRTrans_Params),Pointer                   :: transpar
    Type(AtMolCoords), Dimension(:), Intent(in)        :: species
    Integer, Intent(in)                                :: spc, mol

    If (Associated(transpar%cavity)) Then
      Call cavitylist_updateMolecInfo(transpar%cavity, species(spc), &
          spc, mol, 111 )
    Endif

  End Subroutine rigidmoves_postadjustTrans

  !-------------------------------------------------------------------------
  ! update cavity information after rrotate
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_postadjustRot(rotpar, species, spc, mol)
    Type(RigidRRot_Params),Pointer                   :: rotpar
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: species
    Integer, Intent(in)                              :: spc, mol

    If (Associated(rotpar%cavity)) Then
      Call cavitylist_updateMolecInfo(rotpar%cavity, species(spc), &
          spc, mol, 111 )
    Endif
  End Subroutine rigidmoves_postadjustRot


  !-------------------------------------------------------------------------
  ! update cavity information after bdelete
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_postadjustDel(delpar, species, spc, mol)
    Type(RigidBDel_Params), Pointer                   :: delpar
    Type(AtMolCoords), Dimension(:), Intent(in)    :: species
    Integer, Intent(in)                               :: spc, mol

    If (Associated(delpar%insbias%cavity)) Then
      Call cavitylist_updateMolecInfo(delpar%insbias%cavity, species(spc), &
          spc, mol, 110 )
    Endif

  End Subroutine rigidmoves_postadjustDel

  !-------------------------------------------------------------------------
  ! update cavity information after binsert 
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_postadjustBIns(inspar, species, spc, mol)
    Type(RigidBIns_Params),Pointer                   :: inspar
    Type(AtMolCoords), Dimension(:), Intent(in)      :: species
    Integer, Intent(in)                                :: spc, mol

    If (Associated(inspar%insbias)) Then
      If (Associated(inspar%insbias%cavity)) Then
        Call rigidmoves_postadjustCavity(inspar%insbias, species, spc, mol)
      Endif
    Endif
  End Subroutine rigidmoves_postadjustBIns

  !-------------------------------------------------------------------------
  ! update cavity information after binsert 
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_postadjustLIns(inspar, species, spc, mol)
    Type(RigidLIns_Params),Pointer                   :: inspar
    Type(AtMolCoords), Dimension(:), Intent(in)      :: species
    Integer, Intent(in)                                :: spc, mol

    If (Associated(inspar%insbias)) Then
      If (Associated(inspar%insbias%cavity)) Then
        Call rigidmoves_postadjustCavity(inspar%insbias, species, spc, mol)
      Endif
    Endif
  End Subroutine rigidmoves_postadjustLIns


  !-------------------------------------------------------------------------
  ! update cavity information after binsert or linsert.
  ! takes RigidInsBias as input
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_postadjustCavity(insbias, species, spc, mol)
    Type(RigidInsBias),Pointer                   :: insbias
    Type(AtMolCoords), Dimension(:), Intent(in)         :: species
    Integer, Intent(in)                                 :: spc, mol

    Integer :: totmol, i
    Integer, Save :: outputcount=0, debugUnit=0


    !** Also check whether memory in cubelist is alright
    !** has to be done before updatemolecinfo
    Call cavitylist_checkandIncrMemory(insbias%cavity, species(spc), spc )
    totmol=config_getTotMoles(species)

    Call cavitylist_updateMolecInfo(insbias%cavity,species(spc), spc, &
        mol, 101 )

    !** Check whether cav_prob array is in good shape
    If (totmol>=Size(insbias%cav_prob)) Call rigidmoves_reSetCavProb(insbias)

    !SDEBUG This whole thing should finally go to a display routine
    !** Write cav_prob array to screen only totmol+/-5
    If (Trim(genparams%displaymode)=="VERBOSE") Then
      If (debugUnit==0) debugUnit=file_open("Prob_Cavity.dat ")
      outputcount=outputcount+1
      ! writing max 2000 times of debugInfo
      If (outputcount<2000) Then
        Write(debugUnit,*) '  ------- Cavity Bias information ---------- '
        Write(debugUnit,'(a,t30,a, i2)') "Occupancy", &
            "Cavity Probability For spc ",       spc
        Do i=Max(1,totmol-5), Min(totmol+5, Size(insbias%cav_prob))
          Write(debugUnit,'(a,i6,t30,f16.10)') "MOL#", i, &
              stats_getcavg(insbias%cav_prob(i))/insbias%cavNTrials
        End Do
        Write(debugUnit,*)
      Endif
    Endif

  End Subroutine rigidmoves_postadjustCavity


  !-------------------------------------------------------------------------
  ! Incerease the size of cav_prob array in RInsert_Params by 100 
  !-------------------------------------------------------------------------
  Subroutine rigidmoves_resetCavProb(insbias)
    Type(RigidInsBias),Pointer                     :: insbias

    Type(Statistics) , Dimension(:), Pointer        :: temp_cav_prob
    Integer :: currentsize, error, i

    ! point to original memory
    temp_cav_prob=> insbias%cav_prob

    currentsize=size(insbias%cav_prob)

    ! release old pointer, and reallocate some new memory
    Nullify(insbias%cav_prob)
    Allocate(insbias%cav_prob(currentsize+100), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Do i=1, currentsize
      ! this should copy all the fields of stats objects
      ! Using stats_copy overloading on "=" 
      insbias%cav_prob(i)=temp_cav_prob(i)
    End Do

    Do i=currentsize+1,currentsize+100

      !** - STATS_BLOCKSIZE defined in defaults
      Call stats_init(insbias%cav_prob(i),"Prob Cav"//&
          Trim(int2str(i)), STATS_BLOCKSIZE,.False.)

      !** - We assume that we will be able to find 50% cavities
      !**  at the beginning of the simulation
      !** This is just to avoid divison by zero errors later
      Call stats_update(insbias%cav_prob(i),&
          insbias%cavNTrials/2.00_Rdbl)

    End Do


    !** Release unused space
    Deallocate(temp_cav_prob, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

  End Subroutine rigidmoves_resetCavProb

End Module rigidmoves











