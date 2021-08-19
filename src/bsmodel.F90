!------------------------------------------------------------------------------
! This module contains the objects and methods for dealing with bond 
! stretching.  Like other potential model data structures, it contains a 
! list of pointers to specific bond stretching potential types.  The pointer
! that is initialized directs the code to do that type of calculation.
!
! Needed Improvements:
! 1) remove usage of molecule pointer, this is hack
!------------------------------------------------------------------------------

Module bsmodel

  Use defaults, Only: RDbl,strLen,lstrLen,STRETCH_INDEX
  Use utils, Only: toupper,combine,deallocerrdisplay,allocerrdisplay, &
      fileSrchwoComment,split,toint,stripcmnt,findint
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/),vector_getnorm
  Use morse, Only: MorseParams, morse_idstring, morse_init, &
      morse_getinteraction, morse_display
  Use harmonic, Only: HarmonicParams, harmonic_idstring, harmonic_init, &
      harmonic_getinteraction, harmonic_display
  Use atom, Only: atom_gettypename,atom_getntypes
  Use molecule, Only: MolecularParams, &  !!HACK
      molecule_getnatomtypes,molecule_getatomcoords, molecule_getnthatom
  Use molecules, Only: molecules_name
  Use connects, Only: connects_makeTypeList
  Use storebase, Only: storebase_inc,storebase_nderivs
  Use store, Only: Store_Level_Pair,store_nderivs,store_idbranch
  Use simcell, Only: SimCell_Params, simcell_minimage

  Implicit None
  Save

  Private
  Public :: BondStretchModel, StretchInfo, bsmodel_initparams, bsmodel_eval, &
      bsmodel_getinteraction, bsmodel_bsinit, bsmodel_bsint, isfast, &
      bsmodel_display, bsmodel_bsdisplay,bsmodel_cleanup, bsmodel_bscleanup, &
      bsmodel_isinit, bsmodel_bsinfo, STRETCH_KEY, bsmodel_setHarK

  Character(len=strLen), Parameter  :: STRETCH_KEY = 'STRETCH'

  !** Bond stretching data type, points to potential types
  Type BondStretchModel
    Integer                       :: i
    Type(MorseParams), Pointer    :: mo
    Type(HarmonicParams), Pointer :: har
  End Type BondStretchModel

  !** handles all 2-body interactions, excluding intrapair interactions
  Type StretchInfo
    Integer                                :: nstretch
    Logical                                :: fast
    Character(40)                          :: model
    Integer, Dimension(:,:), Pointer       :: bslist
    Logical, Dimension(:,:), Pointer       :: bsflag
    Type(BondStretchModel), Dimension(:), Pointer :: bsparams
  End Type StretchInfo

  Interface isfast
    Module procedure bsmodel_isfastbs
  End Interface

Contains

  !-------------------------------------------------------------------
  ! Initialize the bond-bending parameters, assumes the input is 
  ! already allocated.
  ! Requires: bspntr -- bond stretching type, to be initialized
  !           modelName -- character string defining model type
  !           bond length -- bond length
  !           params -- array of strings containing model parameters
  !-------------------------------------------------------------------
  Subroutine bsmodel_initparams(bspntr,modelName,blength,params)
    Type(BondStretchModel), Intent(InOut)   :: bspntr
    Character(*), Intent(In)                :: modelName
    Real(Kind=RDbl), Intent(In)             :: blength
    Character(*), Dimension(:), Intent(In)  :: params

    Integer                      :: error

    !** Nullify all the pointers first
    Nullify(bspntr%har)
    Nullify(bspntr%mo)

    Select Case (Toupper(trim(modelName)))

    Case (harmonic_idstring)
      Allocate(bspntr%har,STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'har')

      If (Index(toupper(combine(params)),'CALC') /= 0) Then
        Call harmonic_init(bspntr%har,params,blength)
      Else
        Call harmonic_init(bspntr%har,params)
      End If

    Case (morse_idstring)
      Allocate(bspntr%mo,STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'mo')

      If (Index(toupper(combine(params)),'CALC') /= 0) Then
        Call morse_init(bspntr%mo,params,blength)
      Else
        Call morse_init(bspntr%mo,params)
      End If

    Case Default
      Write(0,*) 
      Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
          " : No such bond stretch model exists"

    End Select

  End Subroutine bsmodel_initparams

  !----------------------------------------------------------------------------
  ! Initialize the bond stretch interaction
  ! Requires: ptr -- pointer to bond stretching information 
  !           molecule -- pointer to molecule data structure
  !           filename -- file where bond bending information can be found
  !----------------------------------------------------------------------------
  Recursive Subroutine bsmodel_bsinit(ptr,molecule,filename)
    Type(StretchInfo), Pointer       :: ptr
    Type(MolecularParams), Pointer   :: molecule
    Character(*), Intent(In)         :: filename
 
    Integer                                  :: i,j,error,unitno,pos
    Integer                                  :: nstretch,usebs,nmatch
    Integer                                  :: atom1,atom2
    Integer                                  :: natypes,ntypes,atype1,atype2
    Real(kind=RDbl)                          :: blength
    Type(VecType)                            :: vec1,vec2
    Character(len=strLen)                    :: text, key
    Character(len=strLen), Dimension(strLen) :: params
    Integer, Dimension(molecule%natoms)      :: atypes
    Integer, Dimension(:,:), Allocatable     :: conlist
    Character(len=lstrLen), Dimension(:,:), Allocatable :: pairparams

    If(.Not.Associated(ptr)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          " : passed pointer must associated"
      Stop
    End If

    !** Open the file containing the intramolecular interaction information
    unitno = file_open(filename,110)
 
    key = "Bond_Stretch"

    !** Read the type of interation for the molecule
    usebs = fileSrchwoComment(unitno,key,text,.True.)

    If (usebs == 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "No bond stretching for this molecule- Something is wrong"
      Stop
    End If

    !** Get the number of atom types
    natypes = atom_getntypes()

    !** initialize the atype flag pointer
    Allocate(ptr%bsflag(natypes,natypes),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bsflag')
    ptr%bsflag = .FALSE.
    
    !** Initialize the number of bends
    ptr%nstretch = 0

    j = split(text,params)
    ptr%model = trim(params(3))
    !** Fast or slow?
    If (Trim(toupper(params(4))) == "FAST") Then
      ptr%fast = .True.
    Else If (trim(toupper(params(4))) == "SLOW") Then
      ptr%fast = .False.
    Else
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          " : Must indicate 'Fast' or 'Slow' for ", key
      Stop
    End If

    !** Look for 'LISTED', 'GENERATED' or another filename
    If (toupper(Trim(params(2))) == 'LISTED') Then

      !** Bond stretching parameters are listed, get number
      Read(unitno,'(a)') text
      j = split(text,params)
      nstretch = toint(params(1))
      ptr%nstretch = nstretch

      !** initialize the model information pointer
      Allocate(ptr%bsparams(ptr%nstretch),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bsparams')
    
      !** initialize the bslist
      Allocate(ptr%bslist(nstretch,2),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bslist')

      Do i = 1,nstretch
        Read(unitno,'(a)') text
        j = split(text,params)
        atom1 = toint(params(1))
        atom2 = toint(params(2))
        ptr%bslist(i,:) = (/atom1, atom2/)

        !** Indicate that this atype doublet has interactions.  Note, 
        !** this is only used in a hack in intramolecular_reconnect
        atype1 = molecule%atoms(atom1)%atype
        atype2 = molecule%atoms(atom2)%atype
        ptr%bsflag(atype1,atype2) = .TRUE.
        ptr%bsflag(atype2,atype1) = .TRUE.
    
        !** get bond length in case we need it
        vec1 = molecule_getatomcoords(molecule,atom1)
        vec2 = molecule_getatomcoords(molecule,atom2)
        blength = vector_getnorm(vec1-vec2)
    
        Call bsmodel_initparams(ptr%bsparams(i),ptr%model,blength,params(3:j))
      End Do
    
    Else If (toupper(Trim(params(2))) == 'GENERATE') Then

      !** We need to generate the bond stretching list
      !** A list of atom pairs is needed for the list
    
      !** Read in the number of types list
      Read(unitno,'(a)') text
      j = split(text,params)
      ntypes = toint(params(1))

      !** Allocate temporary storage for atomtype-based parameters
      Allocate(pairparams(natypes,natypes),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      !** Get the number of atom types
      natypes = molecule_getnatomtypes(molecule,atypes)

      !** Allocate temporary storage for the connectivity list.
      Allocate(conlist(molecule%natoms*Int(0.5*molecule%natoms+1),2),STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i5,a,i12)') __FILE__,":",__LINE__, &
            " Could not allocate conlist. Required size : ", &
            molecule%natoms*Int(0.5*molecule%natoms+1)
        Stop
      End If

      !** Get all the individual atom types
      natypes = molecule_getnatomtypes(molecule,atypes,.True.)
    
      !** Read the interaction types and process
      Do i = 1, ntypes

        !** Read in the type and store it
        !** The type should consist of 2 atom type names and the
        !** parameters corresponding to the model that was specified
        Read(unitno,'(a)') text
        text = stripcmnt(text)
        j = split(text,params)
        atype1 = atom_gettypename(params(1))
        atype2 = atom_gettypename(params(2))
        ptr%bsflag(atype1,atype2) = .TRUE.
        ptr%bsflag(atype2,atype1) = .TRUE.
        pos = findint((/atype1,atype2/),0)
        If (pos /= 0) Then
          Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
              " Could not recognized atom type ",params(pos)
          Stop
        End If

        !** save the parameters
        pairparams(atype1,atype2) = combine(params(3:j))
        pairparams(atype2,atype1) = combine(params(3:j))

        !** With the specific type initialized, we need to build 
        !** the list of atom pairs that correspond to that type
        nmatch = connects_makeTypeList(molecule%connect, atypes, &
            (/atype1,atype2/), conlist(ptr%nstretch+1:Size(conlist,1),:))
    
        !** Update the number of stretch pairs
        ptr%nstretch = ptr%nstretch + nmatch
      End Do

      !** Check to see if this makes sense.
      If (ptr%nstretch == 0) Then
        Write(0,'(a,i6,3a)') __FILE__,__LINE__, &
            " : ERROR, found no bond stretch pairs for ",Trim(molecule%name), &
            " though STRETCH was specified"
        Stop
      End If

      !** initialize the model information pointer
      Allocate(ptr%bsparams(ptr%nstretch),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bsparams')
    
      !** Allocate the bslist. nstretch, each with 2 atoms
      Allocate(ptr%bslist(ptr%nstretch,2),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bslist')

      !** Initialize the bond stretching parameters
      Do i = 1,ptr%nstretch
        atom1 = conlist(i,1)
        atom2 = conlist(i,2)
        ptr%bslist(i,:) = (/atom1, atom2/)
        vec1 = molecule_getatomcoords(molecule,atom1)
        vec2 = molecule_getatomcoords(molecule,atom2)
        blength = vector_getnorm(vec1-vec2)
        atype1 = molecule_getnthatom(molecule, atom1)
        atype2 = molecule_getnthatom(molecule, atom2)
        j = split(pairparams(atype1,atype2),params)
        Call bsmodel_initparams(ptr%bsparams(i),ptr%model,blength,params(1:j))
      End Do

      !** Clean up temporary storage
      Deallocate(conlist,STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'conlist')
      Deallocate(pairparams,STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'pairparams')

    Else
      !** we assume this must be a file. Open it and read it.
      Close(unit=unitno)
      Call bsmodel_bsinit(ptr,molecule,Trim(params(2)))
    End If

  End Subroutine bsmodel_bsinit

  !----------------------------------------------------------------------------
  ! Returns True if the given pointer has a bond stretch interaction
  ! Requires: bspntr -- bond stretching pointer set
  !----------------------------------------------------------------------------
  Logical Function bsmodel_isinit(bspntr)
    Type(BondStretchModel), Intent(In) :: bspntr

    !** Default value
    bsmodel_isinit = .False.

    !** Check each interaction type
    If (Associated(bspntr%har)) bsmodel_isinit = .True.
    If (Associated(bspntr%mo)) bsmodel_isinit = .True.

  End Function bsmodel_isinit

  !----------------------------------------------------------------------------
  ! Interaction energy, potential AND forces (optional).  Either two atomic 
  ! coordinate vectors or a separation vector can be passed in.  
  ! Requires: bspntr -- bond stretching pointer set
  !           coords -- atomic coordinates (2), or separation vector
  !           sepvec -- flag indicating that separation vectors are passed
  !           u -- potential 
  !           vf -- optional forces on the two atoms
  !----------------------------------------------------------------------------
  Subroutine bsmodel_getinteraction(bspntr,coords,sepvec,u,vf)
    Type(BondStretchModel), Intent(In)                 :: bspntr
    Type(VecType), Dimension(:), Intent(In)            :: coords
    Logical, Intent(In)                                :: sepvec
    Real(kind=RDbl), Intent(Out)                       :: u       
    Type(VecType), Optional, Intent(Out), Dimension(2) :: vf

    If (Present(vf)) Then
      If (Associated(bspntr%har)) Then
        Call harmonic_getinteraction(bspntr%har,coords,sepvec,u,vf)
      Else If (Associated(bspntr%mo)) Then
        Call morse_getinteraction(bspntr%mo,coords,sepvec,u,vf)
      End If
    Else
      If (Associated(bspntr%har)) Then
        Call harmonic_getinteraction(bspntr%har,coords,sepvec,u)
      Else If (Associated(bspntr%mo)) Then
        Call morse_getinteraction(bspntr%mo,coords,sepvec,u)
      End If
    End If
    
  End Subroutine bsmodel_getinteraction

  !----------------------------------------------------------------------------
  ! Calculate the bond stretch interaction using actual coordinates of atoms.
  ! Periodic boundary conditions will be applied to the atomic coordinates if
  ! the simcell is passed.
  ! Requires: params -- bond stretching information structure
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Subroutine bsmodel_bsint(params,ffout,coords,simcell)
    Type(StretchInfo), Intent(In)              :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell

    Integer                      :: i,nderivs,level
    Integer                      :: atom1, atom2  !atom numbers in molecule
    Real(kind=RDbl)              :: u,upart
    Type(VecType), Dimension(1)  :: sepvec
    Type(VecType), Dimension(2)  :: force

    !** Get the type of the passed storage
    level = store_idbranch(ffout)

    !** Loop over all the stored stretching pairs in the molecule
    nderivs = store_nderivs(ffout)
    Do i = 1,params%nstretch
      atom1 = params%bslist(i,1)
      atom2 = params%bslist(i,2)

      !** get the separation vector
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,coords(atom2),coords(atom1),sepvec(1))
      Else
        sepvec(1) = coords(atom2) - coords(atom1)
      End If

      !** Do the evaluation and store the results
      If (nderivs == 0) Then
        Call bsmodel_getinteraction(params%bsparams(i),sepvec,.True.,u)
      Else If (nderivs == 1) Then
        Call bsmodel_getinteraction(params%bsparams(i),sepvec,.True.,u,force)

        !** Store the results according to the storage form
        Select Case(level)
        Case(0)
          Call storebase_inc(ffout%binfo(atom1),force(1))
          Call storebase_inc(ffout%binfo(atom2),force(2))
        Case(10)
          Call storebase_inc(ffout%mi(atom1)%total,force(1))
          Call storebase_inc(ffout%mi(atom2)%total,force(2))
        Case Default
          Write(0,'(2a,i5,a,i3)') __FILE__,":",__LINE__, &
              ' Unable to handle storage type ',level
          Stop        
        End Select
      Else
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            ' Second order terms not available'
        Stop        
      End If

      !** Store the potential energy results according to the storage form
      Select Case(level)
      Case(-1)
        Call storebase_inc(ffout%total,u,STRETCH_INDEX)
      Case(0)
        upart = 0.5_RDbl*u
        Call storebase_inc(ffout%binfo(atom1),upart,STRETCH_INDEX)
        Call storebase_inc(ffout%binfo(atom2),upart,STRETCH_INDEX)
      Case(10)
        upart = 0.5_RDbl*u
        Call storebase_inc(ffout%mi(atom1)%total,upart,STRETCH_INDEX)
        Call storebase_inc(ffout%mi(atom2)%total,upart,STRETCH_INDEX)
      Case Default
        Write(0,'(2a,i5,a,i3)') __FILE__,":",__LINE__, &
            ' Unable to handle storage type ',level
        Stop        
      End Select

    End Do

  End Subroutine bsmodel_bsint

#ifdef PREVIOUS
  !----------------------------------------------------------------------------
  ! Calculate the bond stretch interaction using actual coordinates of atoms.
  ! Periodic boundary conditions will be applied to the atomic coordinates if
  ! the simcell is passed.
  ! Requires: params -- bond stretching information structure
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Subroutine bsmodel_bsint(params,ffout,coords,simcell)
    Type(StretchInfo), Intent(In)              :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell

    Integer                      :: i,nderivs
    Integer                      :: atom1, atom2  !atom numbers in molecule
    Real(kind=RDbl)              :: u,upart
    Type(VecType), Dimension(1)  :: sepvec
    Type(VecType), Dimension(2)  :: grad

    !** Loop over all the stored stretching pairs in the molecule
    nderivs = store_nderivs(ffout)
    Do i = 1,params%nstretch
      atom1 = params%bslist(i,1)
      atom2 = params%bslist(i,2)

      !** get the separation vector
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,coords(atom2),coords(atom1),sepvec(1))
      Else
        sepvec(1) = coords(atom2) - coords(atom1)
      End If

      !** Do the evaluation and store the results
      If (nderivs == 0) Then
        Call bsmodel_getinteraction(params%bsparams(i),sepvec,.True.,u)
      Else If (nderivs == 1) Then
        Call bsmodel_getinteraction(params%bsparams(i),sepvec,.True.,u,grad)
        Call storebase_inc(ffout%mi(atom1)%total,grad(1))
        Call storebase_inc(ffout%mi(atom2)%total,grad(2))
      Else
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            ' Second order terms not available'
        Stop        
      End If

      upart = 0.5_RDbl*u
      Call storebase_inc(ffout%mi(atom1)%total,upart,STRETCH_INDEX)
      Call storebase_inc(ffout%mi(atom2)%total,upart,STRETCH_INDEX)
    End Do

  End Subroutine bsmodel_bsint
#endif

  !----------------------------------------------
  ! Checks if the BS interaction is fast or slow
  !----------------------------------------------
  Logical Function bsmodel_isfastbs(stretch)
    Type(StretchInfo) :: stretch

    bsmodel_isfastbs = stretch%fast
  End Function bsmodel_isfastbs

  !--------------------------------------------------------------------
  ! Determines if the interaction should be evaluated
  ! Requires:  ptr -- pointer to interaction parameters structure
  !            fast -- proposed evaluation speed
  !--------------------------------------------------------------------
  Logical Function bsmodel_eval(ptr,fast)
    Type(StretchInfo), Pointer    :: ptr
    Logical, Intent(In)           :: fast

    bsmodel_eval = .False.
    If (.Not. Associated(ptr)) Return
    bsmodel_eval = (ptr%fast .And. fast)

  End Function bsmodel_eval

  !----------------------------------------------------------------------------
  ! Returns all tabulated bond-stretching information
  ! Requires:  info -- stretching data structure for a single species
  !            nsets -- number of potential sets returned
  !            list -- list of atom number pairs in molecule  (pair_index,1:2)
  !            params -- array of strings containing type and parameters
  !----------------------------------------------------------------------------
  Subroutine bsmodel_bsinfo(info,nsets,list,params)
    Type(StretchInfo), Intent(In)                       :: info
    Integer, Intent(Out)                                :: nsets
    Integer, Dimension(:,:), Intent(Out)                :: list
    Character(len=lstrLen), Dimension(:), Intent(Out)   :: params

    Integer                 :: i
    Character(len=lstrLen)  :: paramset

    !** copy the list
    nsets = info%nstretch 

    list(1:nsets,1:2) = info%bslist(1:nsets,1:2)

    !** get the parameters
    Do i = 1,nsets
      paramset = bsmodel_display(info%bsparams(i))
      Write(params(i),'(a,1x,a)') Trim(STRETCH_KEY),Trim(paramset)
    End Do

  End Subroutine bsmodel_bsinfo

  !-------------------------------------------------------------------
  ! Display the bsparams index
  ! Requires: bspntr -- bond stretching pointer set
  !-------------------------------------------------------------------
!!$  Character(len=strLen) Function bsmodel_display(bspntr)
  Function bsmodel_display(bspntr)
    Type(BondStretchModel), Intent(In) :: bspntr

    Character(len=strLen) :: bsmodel_display
    Character(len=strLen) :: vals
    Character(len=strLen) :: mname

    vals = ""

    If (Associated(bspntr%har)) Then
      mname="Harmonic"
      vals = trim(harmonic_display(bspntr%har))
    Else If (Associated(bspntr%mo)) Then
      mname="Morse"
      vals = trim(morse_display(bspntr%mo))
    Else If (.Not.Associated(bspntr%mo).And. &
        .Not.Associated(bspntr%har)) Then
      mname="unassociated"
    Else
      mname="undetermined"
    End If

    Write(bsmodel_display,'(a,2x,a)') trim(mname),trim(vals)

  End Function bsmodel_display

  !-----------------------------------------------
  ! Display the information about bond stretch
  !-----------------------------------------------
  Subroutine bsmodel_bsdisplay(ptr,spc,unit,indent)
    Type(StretchInfo), Pointer :: ptr
    Integer, Intent(In)        :: spc,unit,indent
                             
    Integer                  :: i, atom1, atom2, atype1,atype2
    Character(len=indent)    :: blank
    Character(len=strLen)    :: name
    Character(len=lstrLen)   :: junk

    blank = Repeat(' ',indent)

    name = molecules_name(spc)

    If (Associated(ptr)) Then
      Write(unit,'(3a)') blank, &
          'Bond Stretch Information for ',Trim(name)
      If (ptr%fast) Then
        Write(unit,'(2x,2a,l1)') blank,'Fast Interaction'
      Else
        Write(unit,'(2x,2a,l1)') blank,'Slow Interaction'
      End If
      Do i = 1,ptr%nstretch
        atom1 = ptr%bslist(i,1)
        atom2 = ptr%bslist(i,2)
        Write(unit,'(2a,i3,1x,i3,2a)') blank,'Bond stretch pair :',atom1,atom2
        junk = bsmodel_display(ptr%bsparams(i))
        Write(unit,'(4x,3a)') blank,'Model and Params: ',Trim(junk)
      End Do
    Else
      Write(unit,'(2a)') 'No bond stretch information for ',Trim(name)
    End If

  End Subroutine bsmodel_bsdisplay

  !----------------------------------------------------------------
  ! Deallocate pointers and the like
  ! Requires: bspntr -- bond stretching pointer set
  !----------------------------------------------------------------
  Subroutine bsmodel_cleanup(bspntr)
    Type(BondStretchModel), Dimension(:), Pointer :: bspntr

    Integer        :: i,error,nlist

    nlist = Size(bspntr)

    !** Deallocate all allocated pointers
    Do i = 1,nlist
      If (Associated(bspntr(i)%har)) Then
        Deallocate(bspntr(i)%har,STAT=error)
        If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'har')
      End If

      If (Associated(bspntr(i)%mo)) Then
        Deallocate(bspntr(i)%mo,STAT=error)
        If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'mo')
      End If
    End Do

  End Subroutine bsmodel_cleanup

  !--------------------------------------------------------
  ! Cleanup bs allocations
  !--------------------------------------------------------
  Subroutine bsmodel_bscleanup(stretch)

    Type(StretchInfo) :: stretch
    Integer :: error

    Call bsmodel_cleanup(stretch%bsparams)
    Deallocate(stretch%bsparams,STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'bsparams')

    Deallocate(stretch%bslist,STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'bslist')

    Deallocate(stretch%bsflag,STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'bsflag')

  End Subroutine bsmodel_bscleanup



  !----------------------------------------------------------------------------
  ! set harmonic params to a set valu
  !----------------------------------------------------------------------------
  Subroutine bsmodel_setHarK(bspntr, alist,harK)
    Type(StretchInfo), Pointer                 :: bspntr
    Integer, Dimension(2), Intent(In)            :: alist
    Real(kind=RDbl), Intent(in)                       :: harK
    Integer :: index,atom1, atom2 ,i
    index=0
    Do i = 1,bspntr%nstretch
      index=i
      atom1 = bspntr%bslist(i,1)
      atom2 = bspntr%bslist(i,2)
      If ( (atom1==alist(1)).And.(atom2==alist(2)) ) Exit
    End Do
    If (index==0) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    If (Associated(bspntr%bsparams(index)%har)) Then
      bspntr%bsparams(index)%har%k=harK
    Else 
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If
  End Subroutine bsmodel_setHarK
  

End Module bsmodel
