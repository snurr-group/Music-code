!------------------------------------------------------------------------------
! This module contains the objects and methods for dealing with 4-body
! interactions.  Like other potential model data structures, it contains a 
! list of pointers to specific 4-body potential types.  The pointer
! that is initialized directs the code to do that type of calculation.
!
! Needed Improvements
! 1) Write code for the GENERATE option
!------------------------------------------------------------------------------

Module tormodel

  Use defaults, Only: RDbl,strLen,lstrLen,TORSION_INDEX,degTorad,xlstrLen
  Use utils, Only: toupper,combine,deallocerrdisplay,allocerrdisplay, &
      filesrchwocomment,filesrchstrAll,split,toint,stripcmnt,findint,chklogical
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/),vector_getnorm
  Use cosexpansion, Only: CosExpansionModel, cosexpansion_init, &
      cosexpansion_getinteraction, cosexpansion_display, cosexpansion_idstring
  Use cosexpansiondr, Only: CosExpansionModelDr, cosexpansiondr_init, &
      cosexpansiondr_getinteraction, cosexpansiondr_display
  Use torjorg, Only: JorgTorsionModel,torjorg_idstring, torjorg_init, &
      torjorg_getinteraction,torjorg_display,torjorg_clean
  Use mondello, Only : MondelloModel,mondello_init,mondello_display,&
      mondello_getinteraction, mondello_idstring
  Use outofplane, Only : outofplaneModel, outofplane_init, outofplane_display, &
      outofplane_getinteraction,outofplane_cleanup, outofplane_idstring
  Use atom, Only: atom_gettypename,atom_getntypes,atom_nbonds
  Use molecules, Only: molecules_name, molecules_getatype, &
      molecules_getatomcoords
  Use connects, Only: connects_makeTypeList
  Use storebase, Only: storebase_inc
  Use store, Only: Store_Level_Pair,store_nderivs,store_idbranch
  Use simcell, Only: SimCell_Params, simcell_minimage

  Implicit None
  Save

  Private
  Public :: TorsionModel, TorsionInfo, tormodel_torinit, tormodel_torint, &
      tormodel_tordisplay, tormodel_torinfo, tormodel_torcleanup, isfast, &
      tormodel_eval, TORSION_KEY, tormodel_setA1A2

  Character(len=strLen), Parameter  :: TORSION_KEY = 'TORSION'
  Character(len=strLen), Parameter  :: tormodel_idstring = 'Bond_Torsion'

  Type TorsionModel
    Type(CosExpansionModel), Pointer   :: cosexp
    Type(CosExpansionModelDr), Pointer :: cosexpdr
    Type(JorgTorsionModel), Pointer    :: jorg
    Type(MondelloModel), Pointer       :: mond
    Type(OutofplaneModel), Pointer     :: outofplane
  End Type TorsionModel

  !** handles all 4-body interactions  
  Type TorsionInfo         
    Integer                                   :: ntorsion
    Logical                                   :: fast
    Character(strLen)                         :: model
    Integer, Dimension(:,:), Pointer          :: torlist
    Logical, Dimension(:,:,:,:), Pointer      :: torflag
    Type(TorsionModel), Dimension(:), Pointer :: torparams
  End Type TorsionInfo

  Interface isfast
    Module procedure tormodel_isfasttor
  End Interface

Contains

  !----------------------------------------------------------------------
  ! Initialize the 4-body parameters, assumes the input is allocated
  ! Requires:  torparams -- 4-body pointer set, to be initialized
  !            modelName -- character string defining model type
  !            params -- array of strings containing model parameters
  !            acoords -- coordinates for the four atoms, might be needed
  !----------------------------------------------------------------------
  Subroutine tormodel_init(torparams,modelname,params,acoords)
    Type(TorsionModel), Intent(Out)         :: torparams
    Character(*), Intent(In)                :: modelname
    Character(*), Dimension(:), Intent(In)  :: params
    Type(VecType), Dimension(4), Intent(In) :: acoords

    Integer              :: error

    !** Nullify all the pointers first
    Nullify(torparams%cosexp)
    Nullify(torparams%cosexpdr)
    Nullify(torparams%jorg)
    Nullify(torparams%mond)
    Nullify(torparams%outofplane)

    !** Select appropriate model to initialize and do it
    Select Case (toupper(Trim(modelname)))
    Case(cosexpansion_idstring)
      Allocate(torparams%cosexp,Stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call cosexpansion_init(torparams%cosexp,params)

    Case('COSEXPANSIONDR')
      Allocate(torparams%cosexpdr,Stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call cosexpansiondr_init(torparams%cosexpdr,params)

    Case(torjorg_idstring)
      Allocate(torparams%jorg,Stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call torjorg_init(torparams%jorg,params)

    Case(mondello_idstring)
      Allocate(torparams%mond,Stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call mondello_init(torparams%mond,params,acoords)

    Case(outofplane_idstring)
      Allocate(torparams%outofplane,Stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call outofplane_init(torparams%outofplane,params)

    Case Default
      Write(0,'(2a,i4,2a)') __FILE__, " : ", __LINE__, &
          "No such 4-body potential model recognized: ",Trim(modelname)
      Stop
      
    End Select

  End Subroutine tormodel_init

  !----------------------------------------------------------------------------
  ! Initialize ALL the torsional angle info
  ! Requires:  ptr -- pointer to torsion information 
  !            spc -- species number
  !            filename -- file where torsion information can be found
  !----------------------------------------------------------------------------
  Subroutine tormodel_torinit(ptr,spc,filename)
    Type(TorsionInfo), Pointer              :: ptr
    Integer, Intent(In)                     :: spc
    Character(*), Intent(In)                :: filename

    Integer                    :: i,j,error,unitno,section
    Integer                    :: usetor,idx,natypes,nsections
    Logical                    :: fast,last_fast
    Character(len=strLen)      :: model
    Character(len=255)         :: text
    Character(len=strLen), Dimension(strLen) :: params

    !** Make sure passed pointer is associated
    If (.Not. Associated(ptr)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          " : passed pointer must associated"
      Stop
    End If

    !** Open the file containing the intramolecular interaction information
    unitno = file_open(filename,110)

    !** Get the number of occurances of the identifier key
    nsections = filesrchstrAll(unitno,tormodel_idstring,.True.)
    If (nsections == 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "No torsion information for this molecule- Something is wrong"
      Stop
    End If

    !** Get the number of atom types and initialize the atype flag pointer
    natypes = atom_getntypes()
    Allocate(ptr%torflag(natypes,natypes,natypes,natypes),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    ptr%torflag = .FALSE.

    !** Loop through sections and estimate number of torsion sets
    ptr%ntorsion = 0  
    Rewind(unit=unitno)
    Do section = 1,nsections
      usetor = fileSrchwoComment(unitno,tormodel_idstring,text,.False.)
      If (usetor == 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            '  Unexpected failure to find torsion model section'
        Stop
      End If

      If (Index(ToUpper(text),'LISTED') /= 0) Then
        Read(unitno,'(a)') text
        j = split(text,params)
        ptr%ntorsion = ptr%ntorsion + toint(params(1))

      Else If (Index(ToUpper(text),'GENERATE') /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            '  Not prepared to estimate number of torsion for GENERATE option'
        Stop
      Else 
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            '  Could not find recognized keyword'
        Stop
      End If
    End Do

    !** Allocate torsion parameter list, 4 atoms in each set
    Allocate(ptr%torlist(ptr%ntorsion,4),STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a,i4,a)') __FILE__,": ",__LINE__, &
          " : Error allocating array torsion%torlist of size (", &
          ptr%ntorsion," x 4 )"
      Stop
    End If

    !** Initialize the model information pointer array
    Allocate(ptr%torparams(ptr%ntorsion),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'torparams')

    !** Handle initialization for each torsion section
    Rewind(unit=unitno)
    idx = 1
    Do section = 1,nsections
      !** Call the section initialization routine 
      Call tormodel_initsection(ptr,spc,idx,filename,fast)
      
      !** Check 'fast' flags for consistency (must all be the same)
      If (section > 1) Then
        If (.Not. chklogical(last_fast,fast)) Then
          Write(*,*) 'here ',last_fast,fast
          Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              '  All torsion interactions must have same "speed"'
          Stop
        End If
      End If
      last_fast = fast

    End Do

  End Subroutine tormodel_torinit

  !----------------------------------------------------------------------------
  ! Initialize the torsional angle info from a single section
  ! Requires:  ptr -- pointer to torsion information 
  !            spc -- species number
  !            idx -- starting index in the parameters list
  !            filename -- file where torsion information can be found
  !            fast -- returns interaction 'speed' for external checking
  !----------------------------------------------------------------------------
  Recursive Subroutine tormodel_initsection(ptr,spc,idx,filename,fast)
    Type(TorsionInfo), Pointer       :: ptr
    Integer, Intent(In)              :: spc
    Integer, Intent(InOut)           :: idx
    Character(*), Intent(In)         :: filename
    Logical, Intent(Out)             :: fast

    Integer                    :: error, unitno, i, j, n
    Integer                    :: atom1, atom2, atom3, atom4, atm
    Integer                    :: atype1, atype2, atype3, atype4
    Integer                    :: usetor, nlisted
    Character(len=strLen)      :: model
    Character(len=255)         :: text
    Type(VecType), Dimension(4)              :: acoords
    Character(len=strLen), Dimension(strLen) :: params

    !** Make sure passed pointer is associated
    If (.Not. Associated(ptr)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          " : passed pointer must associated"
      Stop
    End If

    !** Open the file containing the intramolecular interaction information
    unitno = file_open(filename,110)

    !** Find the next idstring in this file
    usetor = fileSrchwoComment(unitno,tormodel_idstring,text,.False.)
    If (usetor == 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          '  Unexpected failure to find torsion model section'
      Stop
    End If
    
    !** Split the line and get the model name
    j = split(text,params)
    model = trim(params(3))
    
    !** Fast or slow?
    fast = .True.
    If (Trim(toupper(params(4))) == "FAST") Then
      ptr%fast = .True.
    Else If (trim(toupper(params(4))) == "SLOW") Then
      ptr%fast = .False.
      fast = .False.
    Else
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          " : Must indicate 'Fast' or 'Slow' with ",Trim(tormodel_idstring)
      Stop
    End If
    
    !** Look for 'LISTED', 'GENERATED' or another filename    
    If (toupper(Trim(params(2))) == 'LISTED') Then
    
      Read(unitno,'(a)') text
      j = split(text,params)
      nlisted = toint(params(1))
      
      !** Read the listed torsion set and associated parameters
      Do i = idx,(idx + nlisted - 1)
        Read(unitno,'(a)') text
        j = split(text,params)
        atom1 = toint(params(1))
        ptr%torlist(i,1) = atom1
        atom2 = toint(params(2))
        ptr%torlist(i,2) = atom2
        atom3 = toint(params(3))
        ptr%torlist(i,3) = atom3
        atom4 = toint(params(4))
        ptr%torlist(i,4) = atom4

        !** Flag interactions for this atype quartet
        atype1 = molecules_getatype(spc,atom1)
        atype2 = molecules_getatype(spc,atom2)
        atype3 = molecules_getatype(spc,atom3)
        atype4 = molecules_getatype(spc,atom4)
        ptr%torflag(atype1,atype2,atype3,atype4) = .True.

        !** Get the atom coordinates in case _init wants to do calculation
        Do atm = 1,4
          acoords(atm) = molecules_getatomcoords(spc,ptr%torlist(i,atm))
        End Do

        !** Account for some special interaction type (?)
        If (Trim(toupper(model)) == "SPECIAL") Then
          ptr%model = params(5)
          Call tormodel_init(ptr%torparams(i),ptr%model,params(6:j),acoords)
        Else
          ptr%model = model
          Call tormodel_init(ptr%torparams(i),ptr%model,params(5:j),acoords) 
        End If
       End Do

       !** increment the starting index
       idx = idx + nlisted
    
    Else If (toupper(Trim(params(2))) == 'GENERATE') Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' : Sorry, GENERATE option is not available for torsions yet'
      Stop      

    Else
      !** assume it's a file we are reading in.
      Close(unitno)
      Call tormodel_initsection(ptr,spc,idx,filename,fast)
    End If

  End Subroutine tormodel_initsection

  !--------------------------------------------------------------------
  ! Get interations by passing the parameters and the coordinates 
  ! to the appropriate model
  ! Requires: torparams -- 4-body pointer set
  !           coords -- atomic coordinates (4)
  !           u -- potential 
  !           vf -- optional forces on the two atoms
  !--------------------------------------------------------------------
  Subroutine tormodel_getinteraction(torparams,coords,ifvec,u,vf)
    Type(TorsionModel), Intent(In)                     :: torparams
    Type(VecType), Dimension(:), Intent(In)            :: coords
    Logical, Intent(In)                                :: ifvec
    Real(kind=RDbl), Intent(OUT)                       :: u            
    Type(VecType), Dimension(:), Intent(OUT), Optional :: vf

    If (Present(vf)) Then
      If (Associated(torparams%cosexp)) Then
        Call cosexpansion_getinteraction(torparams%cosexp,coords,ifvec,u,vf)
      End If

      If (Associated(torparams%cosexpdr)) Then
        Call cosexpansiondr_getinteraction(torparams%cosexpdr,coords,ifvec,u,vf)
      End If

      If (Associated(torparams%jorg)) Then
        Call torjorg_getinteraction(torparams%jorg,coords,ifvec,u,vf)
      End If

      If (Associated(torparams%mond)) Then
        Call mondello_getinteraction(torparams%mond,coords,ifvec,u,vf)
      End If
      
      If (Associated(torparams%outofplane)) Then
        Call outofplane_getinteraction(torparams%outofplane,coords,ifvec,u,vf)
      End If

    Else
      If (Associated(torparams%cosexp)) Then
        Call cosexpansion_getinteraction(torparams%cosexp,coords,ifvec,u)
      End If

      If (Associated(torparams%cosexpdr)) Then
        Call cosexpansiondr_getinteraction(torparams%cosexpdr,coords,ifvec,u)
      End If

      If (Associated(torparams%jorg)) Then
        Call torjorg_getinteraction(torparams%jorg,coords,ifvec,u)
      End If

      If (Associated(torparams%mond)) Then
        Call mondello_getinteraction(torparams%mond,coords,ifvec,u)
      End If

      If (Associated(torparams%outofplane)) Then
        Call outofplane_getinteraction(torparams%outofplane,coords,ifvec,u)
      End If
    End If
    
  End Subroutine tormodel_getinteraction

  !----------------------------------------------------------------------
  ! Get torsion interactions using atomic coordinates as input
  ! Requires: params -- torsion forcefield parameters
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  ! NOTE: simcell usage not yet implemented, see bsmodel for example
  !----------------------------------------------------------------------------
  Subroutine tormodel_torint(params,ffout,coords,simcell)
    Type(TorsionInfo), Intent(In)              :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell

    Integer                     :: i, level
    Integer                     :: nderivs, atom1, atom2, atom3, atom4
    Real(kind=RDbl)             :: u,upart
    Type(VecType), Dimension(4) :: force,tempcoords
     Type(VecType), Dimension(3)       :: sepvec

    !** Get the type of the passed storage
    level = store_idbranch(ffout)

    !** loop through the bending triplets and get interactions 
    nderivs = store_nderivs(ffout)
    Do i = 1,params%ntorsion
      atom1 = params%torlist(i,1)
      atom2 = params%torlist(i,2)
      atom3 = params%torlist(i,3)
      atom4 = params%torlist(i,4)
      tempcoords(1) = coords(atom1)
      tempcoords(2) = coords(atom2)
      tempcoords(3) = coords(atom3)
      tempcoords(4) = coords(atom4)

      !** Get the separation vectors
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,coords(atom2),coords(atom1),sepvec(1))
        Call simcell_minimage(simcell,coords(atom3),coords(atom2),sepvec(2))
        Call simcell_minimage(simcell,coords(atom4),coords(atom3),sepvec(3))
      Else
        sepvec(1) = coords(atom2) - coords(atom1)
        sepvec(2) = coords(atom3) - coords(atom2)
        sepvec(3) = coords(atom4) - coords(atom3)
      End If


      If (nderivs == 0) Then
        Call tormodel_getinteraction(params%torparams(i),sepvec,.True.,u)
      Else If (nderivs == 1) Then
        Call tormodel_getinteraction(params%torparams(i),sepvec,&
            .True.,u,force)

        !** Store the results according to the storage form
        Select Case(level)
        Case(0)
          Call storebase_inc(ffout%binfo(atom1),force(1))
          Call storebase_inc(ffout%binfo(atom2),force(2))
          Call storebase_inc(ffout%binfo(atom3),force(3))
          Call storebase_inc(ffout%binfo(atom4),force(4))
        Case(10)
          Call storebase_inc(ffout%mi(atom1)%total,force(1))
          Call storebase_inc(ffout%mi(atom2)%total,force(2))
          Call storebase_inc(ffout%mi(atom3)%total,force(3))
          Call storebase_inc(ffout%mi(atom4)%total,force(4))
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
        Call storebase_inc(ffout%total,u,TORSION_INDEX)
      Case(0)
        upart = u/4.0_RDbl
        Call storebase_inc(ffout%binfo(atom1),upart,TORSION_INDEX)
        Call storebase_inc(ffout%binfo(atom2),upart,TORSION_INDEX)
        Call storebase_inc(ffout%binfo(atom3),upart,TORSION_INDEX)
        Call storebase_inc(ffout%binfo(atom4),upart,TORSION_INDEX)
      Case(10)
        upart = u/4.0_RDbl
        Call storebase_inc(ffout%mi(atom1)%total,upart,TORSION_INDEX)
        Call storebase_inc(ffout%mi(atom2)%total,upart,TORSION_INDEX)
        Call storebase_inc(ffout%mi(atom3)%total,upart,TORSION_INDEX)
        Call storebase_inc(ffout%mi(atom4)%total,upart,TORSION_INDEX)
      Case Default
        Write(0,'(2a,i5,a,i3)') __FILE__,":",__LINE__, &
            ' Unable to handle storage type ',level
        Stop        
      End Select

    End Do

  End Subroutine tormodel_torint

#ifdef PREVIOUS
  !----------------------------------------------------------------------
  ! Get torsion interactions using atomic coordinates as input
  ! Requires: params -- torsion forcefield parameters
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  ! NOTE: simcell usage not yet implemented, see bsmodel for example
  !----------------------------------------------------------------------------
  Subroutine tormodel_torint(params,ffout,coords,simcell)
    Type(TorsionInfo), Intent(In)              :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell

    Integer                     :: i ! Everybody's favorite counter
    Integer                     :: nderivs, atom1, atom2, atom3, atom4
    Real(kind=RDbl)             :: u,upart
    Type(VecType), Dimension(4) :: force,tempcoords
    Type(VecType), Dimension(3)       :: sepvec

    !** loop through the bending triplets and get interactions 
    nderivs = store_nderivs(ffout)
    Do i = 1,params%ntorsion
      atom1 = params%torlist(i,1)
      atom2 = params%torlist(i,2)
      atom3 = params%torlist(i,3)
      atom4 = params%torlist(i,4)
      tempcoords(1) = coords(atom1)
      tempcoords(2) = coords(atom2)
      tempcoords(3) = coords(atom3)
      tempcoords(4) = coords(atom4)

      !** Get the separation vectors
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,coords(atom2),coords(atom1),sepvec(1))
        Call simcell_minimage(simcell,coords(atom3),coords(atom2),sepvec(2))
        Call simcell_minimage(simcell,coords(atom4),coords(atom3),sepvec(3))
      Else
        sepvec(1) = coords(atom2) - coords(atom1)
        sepvec(2) = coords(atom3) - coords(atom2)
        sepvec(3) = coords(atom4) - coords(atom3)
      End If

      If (nderivs == 0) Then
        Call tormodel_getinteraction(params%torparams(i),sepvec,.True.,u)
      Else If (nderivs == 1) Then
        Call tormodel_getinteraction(params%torparams(i),sepvec, &
            .True.,u,force)
        Call storebase_inc(ffout%mi(atom1)%total,force(1))
        Call storebase_inc(ffout%mi(atom2)%total,force(2))
        Call storebase_inc(ffout%mi(atom3)%total,force(3))
        Call storebase_inc(ffout%mi(atom4)%total,force(4))
      Else
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            ' Second order terms not available'
        Stop        
      End If

      upart = u/4.0_RDbl
      Call storebase_inc(ffout%mi(atom1)%total,upart,TORSION_INDEX)
      Call storebase_inc(ffout%mi(atom2)%total,upart,TORSION_INDEX)
      Call storebase_inc(ffout%mi(atom3)%total,upart,TORSION_INDEX)
      Call storebase_inc(ffout%mi(atom4)%total,upart,TORSION_INDEX)
    End Do

  End Subroutine tormodel_torint
#endif

  !----------------------------------------------------------------------
  ! Get torsion angle interactions
  ! Requires: ptr -- torsion information structure
  !           coords -- atomic coordinates
  !           u -- potential
  !           vf -- optional forces
  !----------------------------------------------------------------------
  Subroutine tormodel_gettorint(ptr,coords,u,vf)
    Type(TorsionInfo), Pointer               :: ptr
    Type(VecType), Intent(IN), Dimension(:)  :: coords
    Real(kind=RDbl), Intent(OUT)             :: u
    Type(VecType), Intent(Out), Optional,Dimension(:) :: vf

    Integer                     :: i ! Everybody's favorite counter
    Integer                     :: atom1, atom2, atom3, atom4
    Real(kind=RDbl)             :: utemp
    Type(VecType), Dimension(4) :: ftemp
    Type(VecType), Dimension(Size(coords)) :: tempcoords
    

    !** initialize energy and forces
    u = 0.0_RDbl
    If (Present(vf)) Then
      Do i = 1, Size(vf)
        vf(i) = 0.0_RDbl
      End Do
    End If

    Do i = 1, ptr%ntorsion
      atom1 = ptr%torlist(i,1)
      atom2 = ptr%torlist(i,2)
      atom3 = ptr%torlist(i,3)
      atom4 = ptr%torlist(i,4)
      tempcoords(1) = coords(atom1)
      tempcoords(2) = coords(atom2)
      tempcoords(3) = coords(atom3)
      tempcoords(4) = coords(atom4)

      If (Present(vf)) Then
        Call tormodel_getinteraction(ptr%torparams(i),tempcoords(1:4), .False., &
            utemp,ftemp)
        vf(atom1) = ftemp(1) + vf(atom1)
        vf(atom2) = ftemp(2) + vf(atom2)
        vf(atom3) = ftemp(3) + vf(atom3)
        vf(atom4) = ftemp(4) + vf(atom4)
      Else
        Call tormodel_getinteraction(ptr%torparams(i),tempcoords(1:4),.False.,utemp)
      End If
      u = utemp + u
    End Do

  End Subroutine tormodel_gettorint

  !----------------------------------------------
  ! Checks if the TOR interaction is fast or slow
  !----------------------------------------------
  Logical Function tormodel_isfasttor(torsion)
    Type(TorsionInfo) :: torsion

    tormodel_isfasttor = torsion%fast
  End Function tormodel_isfasttor

  !--------------------------------------------------------------------
  ! Determines if the interaction should be evaluated
  ! Requires:  ptr -- pointer to interaction parameters structure
  !            fast -- proposed evaluation speed
  !--------------------------------------------------------------------
  Logical Function tormodel_eval(ptr,fast)
    Type(TorsionInfo), Pointer    :: ptr
    Logical, Intent(In)           :: fast

    tormodel_eval = .False.
    If (.Not. Associated(ptr)) Return
    tormodel_eval = (ptr%fast .And. fast)

  End Function tormodel_eval

  !----------------------------------------------------------------------------
  ! Returns all tabulated torsion information
  ! Requires:  info -- bending data structure for a single species
  !            nsets -- number of potential sets returned
  !            list -- list of atom number pairs in molecule  (pair_index,1:2)
  !            params -- array of strings containing type and parameters
  !----------------------------------------------------------------------------
  Subroutine tormodel_torinfo(info,nsets,list,params)
    Type(TorsionInfo), Intent(In)                       :: info
    Integer, Intent(Out)                                :: nsets
    Integer, Dimension(:,:), Intent(Out)                :: list
    Character(len=lstrLen), Dimension(:), Intent(Out)   :: params

    Integer          :: i

    !** copy the list
    nsets = info%ntorsion
    list(1:nsets,1:4) = info%torlist(1:nsets,1:4)

    !** get the parameters
    Do i = 1,nsets
      params(i) = tormodel_display(info%torparams(i))
    End Do

  End Subroutine tormodel_torinfo

  !--------------------------------------------------------------------
  ! Display the model types and parameters
  ! Requires: torparams -- 4-body pointer set
  !--------------------------------------------------------------------
  Function tormodel_display(torparams)
    Character(len=lstrLen)           :: tormodel_display
    Type(TorsionModel)               :: torparams

    Character(len=lstrLen) :: val,mname

    If (Associated(torparams%cosexp)) Then
      mname = "Cosexpansion"
      val = Trim(cosexpansion_display(torparams%cosexp))
      Write(tormodel_display,'(2a)') Trim(mname), Trim(val)
    End If

    If (Associated(torparams%cosexpdr)) Then
      mname = "CosexpansionDR"
      val = Trim(cosexpansiondr_display(torparams%cosexpdr))
      Write(tormodel_display,'(2a)') Trim(mname), Trim(val)
    End If

    If (Associated(torparams%jorg)) Then
      mname = "Jorgensen Torsion Model"
      val = Trim(torjorg_display(torparams%jorg))
      Write(tormodel_display,'(2a)') Trim(mname), Trim(val)
    End If

    If (Associated(torparams%mond)) Then
      mname = "MondelloModel"
      val = Trim(mondello_display(torparams%mond))
      Write(tormodel_display,'(2a)') Trim(mname), Trim(val)
    End If

    If (Associated(torparams%outofplane)) Then
      mname = "OutofPlane"
      val = Trim(outofplane_display(torparams%outofplane))
      Write(tormodel_display,'(2a)') Trim(mname), Trim(val)
    End If

  End Function tormodel_display

  !-----------------------------------------------------------------
  ! Display the information about torsion
  ! Requires: torinfo -- intra-pair information structure
  !           spc -- species number
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !-----------------------------------------------------------------
  Subroutine tormodel_tordisplay(torinfo,spc,indent,unit)
    Type(TorsionInfo)              :: torinfo
    Integer, Intent(In)            :: spc,indent,unit

    Integer                        :: i,atom1,atom2,atom3,atom4
    Character(len=indent)          :: blank
    Character(len=xlstrLen)        :: string

    blank = Repeat(' ',indent)

    string = molecules_name(spc)
    Write(unit,'(4a)') blank,'Four-Body Forcefield Information for "', &
        Trim(string),'"'
    If (torinfo%fast) Then
      Write(unit,'(2x,2a,l1)') blank,'Fast Interaction'
    Else
      Write(unit,'(2x,2a,l1)') blank,'Slow Interaction'
    End If
    
    Do i = 1,torinfo%ntorsion
      atom1 = torinfo%torlist(i,1)
      atom2 = torinfo%torlist(i,2)
      atom3 = torinfo%torlist(i,3)
      atom4 = torinfo%torlist(i,4)

      Write(unit,'(2x,2a,4i3)') blank, 'Atom set:',atom1,atom2,atom3,atom4
      string = tormodel_display(torinfo%torparams(i))
      Write(unit,'(4x,3a)') blank,'Model and Params: ',Trim(string)
    End Do

  End Subroutine tormodel_tordisplay


  !--------------------------------------------------------
  ! Clean up your mess
  ! Requires: torparams -- 4-body pointer set
  !--------------------------------------------------------
  Subroutine tormodel_clean(torparams)
    Type(TorsionModel)       :: torparams

    Integer           :: error

    If (Associated(torparams%cosexp)) Then
      Deallocate(torparams%cosexp,STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'cosexp')
    End If

    If (Associated(torparams%cosexpdr)) Then
      Deallocate(torparams%cosexpdr,STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'cosexpdr')
    End If
    
    If (Associated(torparams%mond)) Then
      Deallocate(torparams%mond,STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'mond')
    End If

    If (Associated(torparams%outofplane)) Then
      Deallocate(torparams%outofplane,STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'outofplane')
    End If
      
  End Subroutine tormodel_clean

  !------------------------------------------
  ! Cleans up torsion stuff
  !------------------------------------------
  Subroutine tormodel_torcleanup(torsion)
    Type(TorsionInfo) :: torsion

    Integer :: i,error

    Do i = 1,torsion%ntorsion
      Call tormodel_clean(torsion%torparams(i))
    End Do

    Deallocate(torsion%torflag,stat=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'torflag')

    Deallocate(torsion%torlist,stat=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'torlist')

  End Subroutine tormodel_torcleanup


  !----------------------------------------------------------------------------
  ! set cos exp a1, a2 to value passed here
  ! a1=val, a2=-val a3--a6=same as old
  !----------------------------------------------------------------------------
  Subroutine tormodel_setA1A2(torpntr, alist,val)
    Type(TorsionInfo), Pointer                 :: torpntr
    Integer, Dimension(4), Intent(In)            :: alist
    Real(kind=RDbl), Intent(in)                       :: val
    Integer :: index,atom1, atom2 ,atom3, atom4, i
    index=0
    Do i = 1,torpntr%ntorsion
      index=i
      atom1 = torpntr%torlist(i,1)
      atom2 = torpntr%torlist(i,2)
      atom3 = torpntr%torlist(i,3)
      atom4 = torpntr%torlist(i,4)
      If ( (atom1==alist(1)).And.(atom2==alist(2)).And.(atom3==alist(3)).And.(atom4==alist(4)) ) Exit
    End Do
    If (index==0) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    If (Associated(torpntr%torparams(index)%cosexp)) Then
      torpntr%torparams(index)%cosexp%cn(0)=val
      torpntr%torparams(index)%cosexp%cn(1)=-val
    Else 
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

  End Subroutine tormodel_setA1A2


End Module tormodel












