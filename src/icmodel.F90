!------------------------------------------------------------------------------
! This module contains the objects and methods for dealing with intra-pair
! coulombic interaction. It considers partial charges with in a molecule 
! and calculates all pairwise interactions for pairs spacified in molecule file
!------------------------------------------------------------------------------

Module icmodel

  Use defaults, Only: RDbl,strLen,lstrLen,INTRACOUL_INDEX,degTorad,MAX_ATOMS
  Use utils, Only: toupper,combine,deallocerrdisplay,allocerrdisplay, &
      filesrchstr,split,toint,stripcmnt,findint
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/),vector_getnorm,mag
  Use coul, Only: Coul_Params, coul_displayparams, coul_init, & 
      coul_calc_interaction, coul_idstring, coul_isinit, &
      coul_displayCutoffs, coul_getcutoff ! , Assignment(=)
  Use atom, Only: atom_gettypename,atom_getntypes,atom_getname
  Use molecule, Only: MolecularParams, & !HACK!
      molecule_getnatomtypes,molecule_getatomcoords
  Use molecules, Only: molecules_getnatomtypes
  Use connects, Only: connects_makeTypeList
  Use storebase, Only: storebase_inc
  Use store, Only: Store_Level_Pair,store_nderivs,store_idbranch
  Use simcell, Only: SimCell_Params, simcell_minimage

  Implicit None

  Private
  Public :: IntracoulModel, IntracoulInfo, icmodel_icinit, icmodel_icint, &
      icmodel_iccleanup, icmodel_isinit, icmodel_getcutoff, icmodel_copy, &
      icmodel_icdisplay, icmodel_icinfo, icmodel_displaycutoffs, isfast, &
      icmodel_eval, INTRACOUL_KEY

  Character(len=strLen), Parameter  :: INTRACOUL_KEY = 'INTRACOUL'

  Type IntracoulModel
    Type(Coul_Params), Pointer         :: coul
  End Type IntracoulModel

  Type IntracoulInfo
    Integer                                   :: nipairs, naway
    Logical                                   :: fast
    Character(len=strLen)                     :: model
    Integer, Dimension(:,:), Pointer          :: iclist
    Type(IntracoulModel), Dimension(:), Pointer   :: icparams
    Type(IntracoulModel), Dimension(:,:), Pointer :: atype_list
  End Type IntracoulInfo

  Interface icmodel_getinteraction
    Module Procedure icmodel_getinteractionSngl
    Module Procedure icmodel_getinteractionMult
  End Interface

  Interface isfast
    Module procedure icmodel_isfastip
  End Interface

Contains
  !----------------------------------------------------------------------
  ! Initialize the intra pair parameters, assumes the input is 
  ! already allocated.
  ! Requires: params -- intra-pair pointer set, to be initialized
  !           modelName -- character string defining model type
  !           specs_line -- string containing model parameters
  ! If only the IntracoulModel variable is passed, it only 
  ! nullifies the two pointers and then exits.
  !----------------------------------------------------------------------
  Subroutine icmodel_initparams(params,modelName,specs_line)
    Type(IntracoulModel), Intent(InOut)     :: params
    Character(*), Intent(In), Optional      :: modelName,specs_line

    Integer                      :: error

    !** Nullify all the pointers first
    Nullify(params%coul)


    !** Check if the modelName is present. If not, return.
    If (.Not.Present(modelName)) Return

    !** Figure out which model to initialize
    Select Case (Toupper(trim(modelName)))
    Case (coul_idstring)
      Allocate(params%coul,STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'Coul')

      Call coul_init(params%coul,specs_line)

    Case Default
      Write(0,*) 
      Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
          " : No such thing exists"

    End Select

  End Subroutine icmodel_initparams

  !----------------------------------------------------------------------------
  ! Initialize intra-coulombic
  ! Requires: ptr -- pointer to intra-coul information 
  !           molecule -- pointer to molecule data structure
  !           filename -- file where intra-coul information can be found
  !----------------------------------------------------------------------------
  Recursive Subroutine icmodel_icinit(ptr,molec,filename,simcell)
    Type(IntracoulInfo), Pointer               :: ptr
    Type(MolecularParams), Pointer             :: molec
    Character(*), Intent(In)                   :: filename
    Type(Simcell_Params), Intent(In), Optional :: simcell

    Integer                          :: error,i,j,k,m
    Integer                          :: atom1,atom2,ninteract,unitno
    Integer                          :: atype1,atype2,nipairs
    Integer                          :: useIP,naway,natypes,n
    Real(Kind=RDbl)                  :: cut, dist, dist2
    Type(VecType)                    :: sepvec
    Character(len=strLen)            :: modelname,atom_name1,atom_name2
    Character(len=lstrLen)           :: key, specs_line
    Character(len=255)               :: text
    Character(len=strLen), Dimension(strLen) :: params
    Integer, Dimension(molec%natoms)         :: alist, atypes
    Integer, Dimension(Int(0.5*molec%natoms*(molec%natoms-1)),2) :: conlist

    If(.Not.Associated(ptr)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          " : passed pointer must associated"
      Stop
    End If

    !** nullify the atom-based parameters list, we may not use it
    Nullify(ptr%atype_list)

    !** Open the file containing the intramolecular interaction information
    unitno = file_open(filename,110)

    key = "Intracoul_Info"

    !** Read the type of interation for the molecule
    useip = filesrchstr(unitno,key,text,.True.)

    If (useip == 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "No intracoul potentials for this molecule- Something is wrong"
      Stop
    End If

    j = split(text,params)
    ptr%model = trim(params(3))
    modelname = trim(params(3))
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

      !** intra-coul parameters are listed, get number
      Read(unitno,'(a)') text
      j = split(text,params)
      nipairs = toint(params(1))
      ptr%nipairs = nipairs

      !** initialize the model information pointer
      Allocate(ptr%icparams(ptr%nipairs),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'icparams')
    
      !** initialize the intra-coul list
      Allocate(ptr%iclist(nipairs,2),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'iclist')

      Do i = 1,nipairs
        Read(unitno,'(a)') text
        j = split(text,params)
        atom1 = toint(params(1))
        atom2 = toint(params(2))
        ptr%iclist(i,:) = (/atom1, atom2/)
        atom_name1 = atom_getname(molec%atoms(atom1)%atype)
        atom_name2 = atom_getname(molec%atoms(atom2)%atype)
!	write(*,'(2(a,2x),a)') Trim(atom_name1),Trim(atom_name2), &
!            Trim(combine(params(3:j)))
        Write(specs_line,'(2(a,2x),a)') Trim(atom_name1),Trim(atom_name2), &
            Trim(combine(params(3:j)))
        Call icmodel_initparams(ptr%icparams(i),ptr%model,specs_line)
                ptr%icparams(i)%coul%charge1=molec%atoms(atom1)%q0
                ptr%icparams(i)%coul%charge2=molec%atoms(atom2)%q0
                ptr%icparams(i)%coul%atom_number1=atom1
                ptr%icparams(i)%coul%atom_number2=atom2
      End Do

    Else If (toupper(Trim(params(2))) == 'GENERATE') Then
              Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          " : Generate option not implemented for IntraCoulombic interactions"
      Stop
    
    Else
             Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          " : This option not implemented for IntraCoulombic interactions"
      Stop
    End If

  End Subroutine icmodel_icinit

  !-------------------------------------------------------------------
  ! Initialize pointer set by copying from an existing set
  ! Requires: newset -- new intra-pair pointer set, to be initialized
  !           oldset -- old intra-pair pointer set
  !-------------------------------------------------------------------
  Subroutine icmodel_copy(newset,oldset)
    Type(IntracoulModel), Intent(Out)     :: newset
    Type(IntracoulModel), Intent(In)      :: oldset

    Integer           :: error
    
    Nullify(newset%coul)


    If (Associated(oldset%coul)) Then
      Allocate(newset%coul,STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'coul')
      newset%coul = oldset%coul
    End If    

  End Subroutine icmodel_copy

  !----------------------------------------------------------------------------
  ! Returns True if the interaction is initialized, else False
  ! Requires: params -- intra-coul pointer set
  !----------------------------------------------------------------------------
  Logical Function icmodel_isinit(params)
    Type(IntracoulModel), Intent(In) :: params

    !** Set the flag to false
    icmodel_isinit = .False.

    If (Associated(params%coul)) Then
      icmodel_isinit = coul_isinit(params%coul)
    End If

  End Function icmodel_isinit

  !----------------------------------------------------------------------------
  ! Returns the value of the cutoff. May return the neighbor cutoff if
  ! NEIGHBOR is present and True.
  ! Requires: params -- intra-coul pointer set
  !           neighbor -- flag to get neighbor cutoff
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function icmodel_getcutoff(params,neighbor)
    Type(IntracoulModel), Intent(In) :: params
    Logical, Optional, Intent(In)    :: neighbor

    !** Pick the correct interaction

    If (Associated(params%coul)) Then
      If (Present(neighbor)) Then
        icmodel_getcutoff = coul_getcutoff(params%coul,neighbor)
      Else
        icmodel_getcutoff = coul_getcutoff(params%coul)
      End If
    End If

  End Function icmodel_getcutoff

  !----------------------------------------------------------------------------
  ! Calculate the intracoul potential for atom pair atype1, atype2 and
  ! return the potential u and (optionally) the force vectors vf
  ! Requires: params -- intra-coul pointer set
  !           sepvec -- separation vector for a single pair
  !           u -- calculated potential energy
  !           vf -- calculated force, optional
  !----------------------------------------------------------------------------
  Subroutine icmodel_getinteractionSngl(params,sepvec,u,vf)
    Type(IntracoulModel), Intent(In)      :: params
    Type(VecType), Intent(In)             :: sepvec
    Real(kind=RDbl), Intent(OUT)          :: u
    Type(VecType), Intent(Out), Optional  :: vf

    Logical       :: loFlag   

    !** check to see which model is associated for this given pair
    If (Associated(params%coul)) Then
      If (Present(vf)) Then
        Call coul_calc_interaction(params%coul,sepvec,u,loFlag,vf)
      Else
        Call coul_calc_interaction(params%coul,sepvec,u,loFlag)
      Endif
      
!      If (.Not. loFlag) Then 
!        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
!            ' The distance is too close, LJ calculation failed '
!        Stop
!      End If
    End If


  End Subroutine icmodel_getinteractionSngl

  !----------------------------------------------------------------------------
  ! Choose the model and get the interaction. Works on a dimensioned
  ! list of separation vectors
  ! Requires: params -- array of intra-coul pointer set, one for each pair of charges
  !           sepvec -- separation vectors
  !           u -- calculated potential energy
  !           vf -- calculated force, optional
  !----------------------------------------------------------------------------
  Subroutine icmodel_getinteractionMult(params,sepvec,u,vf)
    Type(IntracoulModel), Dimension(:), Intent(In)     :: params
    Type(VecType), Dimension(:), Intent(In)            :: sepvec
    Real(kind=RDbl), Intent(Out)                       :: u
    Type(VecType), Dimension(:), Intent(Out), Optional :: vf

    Logical          :: loFlag   
    Integer          :: ic
    Real(kind=RDbl)  :: pot

    u = 0.0_RDbl

    !** Loop through the pairs and call the appropriate model
    Do ic = 1,Size(params)

      If (Associated(params(ic)%coul)) Then
        Call coul_calc_interaction(params(ic)%coul,sepvec(ic),pot,loFlag,vf(ic))
!        If( .Not. loFlag) Then 
!          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
!              ' The distance is too close '
!          Stop
!        End If
      End If
      u = u + pot

    End Do

  End Subroutine icmodel_getinteractionMult

  !----------------------------------------------------------------------------
  ! Calculate the intra-coul interaction using actual coordinates of atoms.
  ! Periodic boundary conditions will be applied to the atomic coordinates if
  ! the simcell is passed.
  ! Requires: params -- intra-coul forcefield parameters
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Subroutine icmodel_icint(params,ffout,coords,simcell)
    Type(IntracoulInfo), Intent(In)            :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell

    Integer                      :: i,nderivs,level
    Integer                      :: atom1, atom2  !atom numbers in molecule
    Real(kind=RDbl)              :: u,upart
    Type(VecType)                :: sepvec,force


    !** Get the type of the passed storage
    level = store_idbranch(ffout)

    !** Loop through the pairs
    nderivs = store_nderivs(ffout)
    Do i = 1, params%nipairs
      atom1 = params%iclist(i,1)
      atom2 = params%iclist(i,2)

      !** Get the separation vector
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,coords(atom1),coords(atom2),sepvec)
      Else
        sepvec = coords(atom1) - coords(atom2)
      End If

!	print*, "Coul interaction",  atom1, atom2

      !** Do the evaluation and store the results
      If (nderivs == 0) Then
        Call icmodel_getinteraction(params%icparams(i),sepvec,u)
      Else If (nderivs == 1) Then
        Call icmodel_getinteraction(params%icparams(i),sepvec,u,force)


        !** Store the results according to the storage form
        Select Case(level)
        Case(0)
          Call storebase_inc(ffout%binfo(atom1),force)
          Call storebase_inc(ffout%binfo(atom2),(force*(-1.0_RDbl)))
        Case(10)
          Call storebase_inc(ffout%mi(atom1)%total,force)
          Call storebase_inc(ffout%mi(atom2)%total,(force*(-1.0_RDbl)))
        Case Default
          Write(0,'(2a,i5,a,i3)') __FILE__,":",__LINE__, &
              ' Unable to handle storage type ',level
          Stop        
        End Select
            
      Else
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            ' Second order terms not presently available'
        Stop        
      End If

      !** Store the potential energy results according to the storage form
      Select Case(level)
      Case(-1)
        Call storebase_inc(ffout%total,u,INTRACOUL_INDEX)
      Case(0)
        upart = 0.5_RDbl*u
        Call storebase_inc(ffout%binfo(atom1),upart,INTRACOUL_INDEX)
        Call storebase_inc(ffout%binfo(atom2),upart,INTRACOUL_INDEX)
      Case(10)
        upart = 0.5_RDbl*u
        Call storebase_inc(ffout%mi(atom1)%total,upart,INTRACOUL_INDEX)
        Call storebase_inc(ffout%mi(atom2)%total,upart,INTRACOUL_INDEX)
      Case Default
        Write(0,'(2a,i5,a,i3)') __FILE__,":",__LINE__, &
            ' Unable to handle storage type ',level
        Stop        
      End Select
    End Do
 
  End Subroutine icmodel_icint

#ifdef PREVIOUS
  !----------------------------------------------------------------------------
  ! Calculate the intra-coul interaction using actual coordinates of atoms.
  ! Periodic boundary conditions will be applied to the atomic coordinates if
  ! the simcell is passed.
  ! Requires: params -- intra-coul forcefield parameters
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Subroutine icmodel_icint(params,ffout,coords,simcell)
    Type(IntracoulInfo), Intent(In)            :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell

    Integer                      :: i,nderivs
    Integer                      :: atom1, atom2  !atom numbers in molecule
    Real(kind=RDbl)              :: u,upart
    Type(VecType)                :: sepvec,force

    !** Loop through the pairs
    nderivs = store_nderivs(ffout)
    Do i = 1, params%nipairs
      atom1 = params%iclist(i,1)
      atom2 = params%iclist(i,2)

      !** Get the separation vector
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,coords(atom1),coords(atom2),sepvec)
      Else
        sepvec = coords(atom1) - coords(atom2)
      End If

      !** Do the evaluation and store the results
      If (nderivs == 0) Then
        Call icmodel_getinteraction(params%icparams(i),sepvec,u)
      Else If (nderivs == 1) Then
        Call icmodel_getinteraction(params%icparams(i),sepvec,u,force)
            
        Call storebase_inc(ffout%mi(atom1)%total,force)
        Call storebase_inc(ffout%mi(atom2)%total,(force*(-1.0_RDbl)))
      Else
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            ' Second order terms not presently available'
        Stop        
      End If

      upart = 0.5_RDbl*u
      Call storebase_inc(ffout%mi(atom1)%total,upart,INTRACOUL_INDEX)
      Call storebase_inc(ffout%mi(atom2)%total,upart,INTRACOUL_INDEX)
    End Do

  End Subroutine icmodel_icint
#endif

  !----------------------------------------------
  ! Checks if the IC interaction is fast or slow
  !----------------------------------------------
  Logical Function icmodel_isfastip(intracoul)
    Type(IntracoulInfo) :: intracoul

    icmodel_isfastip = intracoul%fast
  End Function icmodel_isfastip

  !--------------------------------------------------------------------
  ! Determines if the interaction should be evaluated
  ! Requires:  ptr -- pointer to interaction parameters structure
  !            fast -- proposed evaluation speed
  !--------------------------------------------------------------------
  Logical Function icmodel_eval(ptr,fast)
    Type(IntracoulInfo), Pointer    :: ptr
    Logical, Intent(In)             :: fast

    icmodel_eval = .False.
    If (.Not. Associated(ptr)) Return
    icmodel_eval = (ptr%fast .And. fast)

  End Function icmodel_eval

  !--------------------------------------------------------
  ! Display the model types and parameters
  ! Requires: params -- intra-coul pointer set
  !--------------------------------------------------------
!!$  Character(len=lstrLen) Function icmodel_display(params)
  Function icmodel_display(params)
    Character(len=lstrLen)              :: icmodel_display
    Type(IntracoulModel), Intent(In)    :: params

    !** COUL values
    If (Associated(params%coul)) Then
      Write(icmodel_display,'(a,2x,a)') "COUL", &
          Trim(coul_displayparams(params%coul))
    End If


  End Function icmodel_display

  !--------------------------------------------------------
  ! Display the model cutoffs
  ! Requires: params -- intra-coul pointer set
  !--------------------------------------------------------
!!$  Character(len=lstrLen) Function icmodel_displayCutoffs(params)
  Function icmodel_displayCutoffs(params)
    Character(len=lstrLen)              :: icmodel_displayCutoffs
    Type(IntracoulModel), Intent(In)    :: params

    !** COUL values
    If (Associated(params%coul)) Then
      Write(icmodel_displayCutoffs,'(a)') &
          Trim(coul_displayCutoffs(params%coul))
    End If


  End Function icmodel_displayCutoffs

  !----------------------------------------------------------------------------
  ! Returns all tabulated bond-stretching information
  ! Requires:  info -- stretching data structure for a single species
  !            nsets -- number of potential sets returned
  !            list -- list of atom number pairs in molecule  (pair_index,1:2)
  !            params -- array of strings containing type and parameters
  !----------------------------------------------------------------------------
  Subroutine icmodel_icinfo(info,nsets,list,params)
    Type(IntracoulInfo), Intent(In)                     :: info
    Integer, Intent(Out)                                :: nsets
    Integer, Dimension(:,:), Intent(Out)                :: list
    Character(len=lstrLen), Dimension(:), Intent(Out)   :: params

    Integer          :: i

    !** copy the list
    nsets = info%nipairs 
!LC    Write(*,*) nsets, Size(list,1), Size(list,2)

    list(1:nsets,1:2) = info%iclist(1:nsets,1:2)

    !** get the parameters
    Do i = 1,nsets
      params(i) = icmodel_display(info%icparams(i))
    End Do

  End Subroutine icmodel_icinfo

  !------------------------------------------------------------------------
  ! Display the information about intra-molecular pairwise interactions
  ! Requires: icinfo -- intra-coul information structure
  !           spc -- species number
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !------------------------------------------------------------------------
  Subroutine icmodel_icdisplay(icinfo,spc,indent,unit)
    Type(IntracoulInfo)            :: icinfo
    Integer, Intent(In)            :: spc,indent,unit

    Integer                        :: i,j,natypes
    Character(len=indent)          :: blank
    Character(len=strLen)          :: name1,name2
    Integer, Dimension(MAX_ATOMS)  :: atyp

    blank = Repeat(' ',indent)

    !** atype-based list if it exists    
    If (Associated(icinfo%atype_list)) Then
      natypes = molecules_getnatomtypes(spc,atyp)
      Do i = 1, natypes
        Do j = i, natypes
          name1 = atom_getname(atyp(i))
          name2 = atom_getname(atyp(j))
    
          If (icmodel_isinit(icinfo%atype_list(atyp(i),atyp(j)))) Then
            Write(unit,'(a,3(a,2x))') blank,Trim(name1),Trim(name2), &
                Trim(icmodel_display(icinfo%atype_list(atyp(i),atyp(j))))
            Write(unit,'(4x,3a)') blank,"Cutoffs : ",&
                Trim(icmodel_displayCutoffs(icinfo%atype_list(atyp(i),atyp(j))))
          End If
    
        End Do
      End Do
    End If

    !** individual pairs and their parameters
    If (icinfo%nipairs < 20) Then
      Write(unit,'(2a)') blank,'List of individual pairs (atom1,atom2,info)'
      Do i = 1,icinfo%nipairs
        Write(unit,'(2x,a,i2,a,2i4,4x,a)') blank,i,':',icinfo%iclist(i,1), &
            icinfo%iclist(i,2),Trim(icmodel_display(icinfo%icparams(i)))
      End Do
    Else
      Write(unit,'(2x,3a)') blank,'number of intra-coul > 20,', &
          ' not showing intra-coul doublets'
    End If

  End Subroutine icmodel_icdisplay

  !--------------------------------------------------------
  ! Cleanup the intracoul pointer set
  ! Requires: params -- intra-coul pointer set
  !--------------------------------------------------------
  Subroutine icmodel_cleanup(params)
    Type(IntracoulModel), Intent(InOut) :: params

    Integer          :: error
    
    If (Associated(params%coul)) Then
      Deallocate(params%coul,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'coul')
    End If

  End Subroutine icmodel_cleanup

  !------------------------------------------
  ! Cleanup coul info
  !------------------------------------------
  Subroutine icmodel_iccleanup(intracoul)
    Type(IntracoulInfo)            :: intracoul

    Integer       :: i,j,error

    Do i = 1,Size(intracoul%icparams)
      Call icmodel_cleanup(intracoul%icparams(i))
    End Do

    Do i = 1,Size(intracoul%atype_list,1)
      Do j = 1,Size(intracoul%atype_list,2)
        Call icmodel_cleanup(intracoul%atype_list(i,j))
      End Do
    End Do

    Deallocate(intracoul%iclist,stat=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'charges')

  End Subroutine icmodel_iccleanup

End Module icmodel



















