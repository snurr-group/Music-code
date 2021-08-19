!---------------------------------------------------------------------------
! This module handles choosing the correct bond length constraints
! model, including the Ciccotti rigid unit constraints which need the
! bond constraint information.  The bond length constraint information
! is read from the control file and processed in this module
!
! Needed Improvements:
! 1) remove usage of molecule pointer, this is hack
! 2) storage not fully incorporate, test and fix
!---------------------------------------------------------------------------

Module conmodel

  Use defaults, Only: RDbl,strLen,lstrLen,TORSION_INDEX,degTorad, &
      MAX_ATOMTYPES
  Use utils, Only: toupper,combine,deallocerrdisplay,allocerrdisplay, &
      filesrchstr,split,toint,stripcmnt,findint,toreal
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/),vector_getnorm,mag
  Use evansmorriss, Only: EvansMorrissModel, evansmorriss_getnconstr, &
      evansmorriss_init, evansmorriss_getconstraints, &
      evansmorriss_getpenalty, evansmorriss_display, evansmorriss_cleanup
  Use ciccotti, Only: ciccotti_mapforces,CiccottiRigidSetInfo,ciccotti_init, &
      ciccotti_display, ciccotti_regenSet, ciccotti_getdofred
  Use atom, Only: atom_gettypename, atom_getname, atom_invmass
  Use molecule, Only: MolecularParams, &  !!HACK!!
      molecule_getnatomtypes,molecule_getatomcoords
  Use molecules, Only : molecules_getnthatom,molecules_name
  Use connects, Only: connects_makeXList
  Use storebase, Only: storebase_inc
  Use store, Only: Store_Level_Pair,store_nderivs
  Use simcell, Only: SimCell_Params, simcell_minimage

  Implicit None
  Save

  Private
  Public :: ConstraintModel, ConstraintInfo, conmodel_getnconstr, isfast, &
      conmodel_coninit, conmodel_getconint, conmodel_getpenalty, conmodel_eval, &
      conmodel_concleanup, conmodel_display, conmodel_postadjust,  &
      conmodel_getdofred, conmodel_condisplay, conmodel_initarray, CONSTRAINT_KEY

  Character(len=strLen), Parameter  :: CONSTRAINT_KEY = 'CONSTRAINT'

  Type ConstraintModel
    Type(EvansMorrissModel), Pointer    :: em
    Type(CiccottiRigidSetInfo), Pointer :: ru
  End Type ConstraintModel

  Type ConstraintInfo
    Logical               :: fast
    Character(len=strLen) :: model
    Type(ConstraintModel) :: conparams
  End Type ConstraintInfo

  Interface isfast
    Module procedure conmodel_isfastcon
  End Interface

Contains

  !----------------------------------------------------------------------------
  ! Get the number of constraints from the model
  !----------------------------------------------------------------------------
  Integer Function conmodel_getnconstr(conparams)
    Type(ConstraintModel) :: conparams
    If (Associated(conparams%em)) Then
      conmodel_getnconstr = evansmorriss_getnconstr(conparams%em)
    End If
  End Function conmodel_getnconstr

  !----------------------------------------------------------------------------
  ! Get the reduction in number of degrees of freedom from the model
  !----------------------------------------------------------------------------
  Integer Function conmodel_getdofred(conparams)
    Type(ConstraintModel) :: conparams

    conmodel_getdofred = 0

    If (Associated(conparams%em)) Then
      conmodel_getdofred = conmodel_getdofred +  &
          evansmorriss_getnconstr(conparams%em)
    End If

    If (Associated(conparams%ru)) Then
      conmodel_getdofred = conmodel_getdofred +  &
          ciccotti_getdofred(conparams%ru)
    End If

  End Function conmodel_getdofred

  !----------------------------------------------------------------------------
  ! Initialize the constraints model and nullify the pointers.
  !----------------------------------------------------------------------------
  Subroutine conmodel_initarray(conparams,model)
    Type(ConstraintModel),Intent(inout) :: conparams ! Not pointer
    Integer                             :: error
    Character(*)                        :: model

    Nullify(conparams%em)
    Nullify(conparams%ru)

    If (Index(toupper(Trim(model)),'EVANSMORRISS') /= 0) Then
      If (.Not.Associated(conparams%em)) Then
        Allocate(conparams%em,stat=error)
        If (error /= 0) Then
          Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
              " : Could not allocate memory for EvansMorriss 'conparams'"
          Stop
        End If
      End If
    End If

    If (Index(toupper(Trim(model)),'CICCOTTI') /= 0) Then
      If (.Not.Associated(conparams%ru)) Then
        Allocate(conparams%ru,stat=error)
        If (error /= 0) Then
          Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
              " : Could not allocate memory for Ciccotti 'conparams'"
          Stop
        End If
        Nullify(conparams%ru%unit)  
      End If
    End If

   If (Index(toupper(Trim(model)),'EVANSMORRISS') == 0) Then
      If (Index(toupper(Trim(model)),'CICCOTTI') == 0) Then
        Write(0,'(2a,i4,2a)') __FILE__, " : ", __LINE__, &
              " : Unknown constraint model " , toupper(Trim(model))
                  Stop
       End If
   End If

  End Subroutine Conmodel_initarray

  !----------------------------------------------------------------------------
  ! Find the bond length constraint pairs and initialize the constraint
  ! models and parameters.
  ! Requires: conparams -- the constraint data structure
  !           genflag -- determines if constraints are to be generated
  !           natoms -- number of atoms in molecule
  !           atomcoords -- coordinates of all atoms in molecule 
  !           setlist -- a list of "set" numbers, one for each atom
  !           typelist -- a list of "type" numbers, one for each atom
  !           atomlist -- a list of atom types, one for each atom number
  !           invmass -- a list of inverse masses, one for each atom number
  !           pairlist -- a list of connected pairs 
  !           unitno -- that unit number where the constraint info is
  !----------------------------------------------------------------------------
  Subroutine conmodel_initparams(conparams,genflag,natoms,atomcoords, &
      setlist,typelist,atomlist,invmass,pairlist,unitno)
    Type(ConstraintModel),Intent(inout)      :: conparams !Not pointer again...
    Logical                                  :: genflag
    Integer, Intent(In)                      :: natoms
    Type(VecType), Dimension(natoms), Intent(In)  :: atomcoords
    Integer, Dimension(natoms), Intent(In)        :: setlist,typelist,atomlist
    Real(kind=RDbl), Dimension(natoms), Intent(In):: invmass
    Integer, Dimension(:,:), Intent(In)           :: pairlist
    Integer, Intent(In)                           :: unitno
                                            
    Integer                                  :: i,j,atype1,atype2,ncon
    Integer                                  :: atom1,atom2
    Integer, Dimension(natoms*3)             :: atom1list,atom2list
    Real(Kind=RDbl)                          :: locut,hicut,dist
    Real(Kind=RDbl), Dimension(natoms*3)     :: lengthlist
    Character(len=strLen), Dimension(strLen) :: params
    Character(len=strLen)                    :: text
    Real(kind=RDbl), Dimension(MAX_ATOMTYPES,MAX_ATOMTYPES) :: baselength
    Real(kind=RDbl), Dimension(MAX_ATOMTYPES,MAX_ATOMTYPES) :: tolerance

    !** Read the number of lines in constraint specification
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    j = split(text,params)
    ncon = toint(params(1))

    If (ncon > (natoms*3)) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Number of constraints greater than 3N (why?)"
        Stop
    End If

    !** Define the length constraints here
    If (.Not.genflag) Then

      !** read the listed constraints
      Do i = 1, ncon
        Read(unitno,'(a)') text
        text = stripcmnt(text)
        j = split(text,params)
        atom1list(i) = toint(params(1))
        atom2list(i) = toint(params(2))
        If ((atom1list(i) > natoms).Or.(atom2list(i) > natoms)) Then
          Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
              " : Atom numbers in constraints must not exceed Natoms"
          Stop
        End If
        If (Index(toupper(Trim(params(3))),'CALC') == 0) Then
          lengthlist(i) = toreal(params(3))
        Else
!LC          Write(*,*) i,atom1list(i),atom2list(i)
!LC          Write(*,'(2i2,3f8.3)') i,atom1list(i),atomcoords(atom1list(i))
!LC          Write(*,'(2i2,3f8.3)') i,atom2list(i),atomcoords(atom2list(i))
          lengthlist(i) = mag(atomcoords(atom1list(i)) - &
              atomcoords(atom2list(i)))
        End If
      End Do

    Else
      tolerance = 0.0_RDbl
      baselength = 0.0_RDbl

      !** read the generic atom type based constraint types
      Do i = 1,ncon
        Read(unitno,'(a)') text
        text = stripcmnt(text)
        j = split(text,params)
        atype1 = atom_gettypename(trim(params(1)))
        atype2 = atom_gettypename(trim(params(2)))
        If (Index(toupper(Trim(params(3))),'CALCULATE') == 0) Then
          baselength(atype1,atype2) = toreal(params(3))
          baselength(atype2,atype1) = baselength(atype1,atype2)
        Else
          baselength(atype1,atype2) = mag(atomcoords(atom1list(i) - &
              atom2list(i)))
          baselength(atype2,atype1) = baselength(atype1,atype2)
        End If
        If (j == 4) Then
          tolerance(atype1,atype2) = toreal(params(4))
          tolerance(atype2,atype1) = tolerance(atype1,atype2)
        End If
        If (j > 4) Then
          Write(0,'(2a,i4,2a)') __FILE__, " : ", __LINE__, &
              " : Too many parameters on line: ",Trim(text)
          Stop
        End If
      End Do

      !** translate the generic types to a full list
      !** note that only pairs that are in the connects structure are checked
      ncon = 0
      Do i = 1,Size(pairlist,1)
        atom1 = pairlist(i,1)
        atom2 = pairlist(i,2)
        dist = mag(atomcoords(atom1) - atomcoords(atom2))
        locut = baselength(atomlist(atom1),atomlist(atom2)) - &
            tolerance(atomlist(atom1),atomlist(atom2))
        hicut = baselength(atomlist(atom1),atomlist(atom2)) + &
            tolerance(atomlist(atom1),atomlist(atom2))

        If ((dist >= locut).And.(dist <= hicut)) Then
          ncon = ncon + 1
          atom1list(ncon) = atom1
          atom2list(ncon) = atom2
          lengthlist(ncon) = dist
        End If
      End Do

      If (ncon > (natoms*3)) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Number of constraints greater than 3N (why?)"
        Stop
      End If
    End If    

!LC    Write(*,'(2a,i4)') __FILE__," : ",__LINE__                        
!LC    Write(*,*) 'number of bonding connections: ',ncon                 
!LC    Do i = 1,ncon                                                     
!LC      Write(*,'(2i4,f8.3)') atom1list(i),atom2list(i),lengthlist(i)   
!LC    End Do                                                            

    If (Associated(conparams%em)) Then
      Call evansmorriss_init(conparams%em,natoms,ncon,invmass, &
          atom1list(1:ncon),atom2list(1:ncon),lengthlist(1:ncon))
    End If

    If (Associated(conparams%ru)) Then
      Call ciccotti_init(conparams%ru,natoms,atomcoords,invmass,setlist, &
          typelist,atom1list(1:ncon),atom2list(1:ncon),lengthlist(1:ncon))
    End If

  End Subroutine conmodel_initparams

  !----------------------------------------------------------------------------
  ! Initialize constraints
  !----------------------------------------------------------------------------
  Recursive Subroutine conmodel_coninit(ptr,molecule,filename)

    Type(ConstraintInfo), Pointer    :: ptr
    Type(MolecularParams), Pointer   :: molecule
    Character(*), Intent(In)         :: filename

    Integer                                     :: error,j,nparams,a
    Logical                                     :: found_speed
    Character(len=strLen)                       :: key
    Character(len=255)                          :: text
    Character(len=strLen), Dimension(strLen)    :: params
    Integer                                     :: usecon, unitno
    Integer, Dimension(molecule%natoms)         :: alist,setlist,typelist

    !** number of pairs can be larger than natoms, especially if 
    !** you have a larger molecule, with many connections. 
    !** The previous definition will work strictly only if it is a 
    !** linear moleule
    ! PREVIOUS :   Integer, Dimension(molecule%natoms,2)       :: pairlist
    ! NEW HACK : ->
    Integer, Dimension(molecule%natoms*3,2)       :: pairlist

    Type(VecType), Dimension(molecule%natoms)   :: atomcoords
    Real(kind=RDbl), Dimension(molecule%natoms) :: invmass

    If(.Not.Associated(ptr)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          " : passed pointer must associated"
      Stop
    End If

    !** Open the file containing the intramolecular interaction information
    unitno = file_open(filename,110)

    key = "Bond_Constraints"

    !** Read the type of interation for the molecule
    usecon = filesrchstr(unitno,key,text,.True.)
    
    If (usecon /= 0) Then

      !MDEBUG
      Write(*,*) "Initializing constraints for ",trim(molecule%name)

      Allocate(ptr,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            " : Error allocating pointer intrapair"
        Stop
      End If

      text = stripcmnt(text)
      nparams = split(text,params)
      ptr%model = ''
      found_speed = .False.
      !** search through the string and extract model(s) and
      !** integration "speed".  Initialize models as necessary.
      !** The first two entries will be important elsewhere
      Do j = 3,nparams  
        Select Case (toupper(params(j)))
        Case ('FAST')
          ptr%fast = .True.
          found_speed = .True.
        Case ('SLOW')
          ptr%fast = .False.
          found_speed = .True.

        Case Default
          Write(ptr%model,'(a,1x,a)') &
              Trim(ptr%model),Trim(params(j))
        End Select
      End Do

      !** Fast or slow?  One or the other must be declared
      If (.Not.(found_speed)) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            " : Must indicate 'Fast' or 'Slow' for ", key
        Stop
      End If

      !** Initialize model type
      Call conmodel_initarray(ptr%conparams, &
          ptr%model)

      !** get the pair list
      j = connects_makeXList(molecule%connect,1,pairlist)

      !** get various arrays from the molecule structure
      alist = molecule%atoms(1:molecule%natoms)%atype
      setlist = molecule%atoms(1:molecule%natoms)%set
      typelist = molecule%atoms(1:molecule%natoms)%type
      atomcoords = molecule%atoms(1:molecule%natoms)%r

      !** get the atomic inverse mass list
      Do a = 1,molecule%natoms
        invmass(a) = atom_invmass(molecule%atoms(a)%atype)
      End Do

      !** generate or listed? Second position must be the indicator
      If (toupper(Trim(params(2))) == "GENERATE") Then
        Call conmodel_initparams(ptr%conparams,.True., &
            molecule%natoms,atomcoords,setlist,typelist,alist,invmass, &
            pairlist(1:j,1:2),unitno)
      Else If (toupper(Trim(params(2))) == "LISTED") then
        Call conmodel_initparams(ptr%conparams,.False., &
            molecule%natoms,atomcoords,setlist,typelist,alist,invmass, &
            pairlist(1:j,1:2),unitno)
      Else
        !** assume it's a file we are reading in.
        Close(unitno)
        Call conmodel_coninit(ptr,molecule,trim(params(2)))
      End If
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "No constraints for this molecule- Something wrong ?"
    End If

  End Subroutine conmodel_coninit

  !----------------------------------------------------------------------------
  ! Calculate the constraint interactions
  ! Requires: ptr -- the constraint model data structure
  !           coords -- coordinates of one molecule
  !           v -- velocities of one molecule
  !           u -- potential energies of one molecule
  !           f -- forces on one molecule
  !           all_flag -- if TRUE, include possible Ciccotti calls
  !----------------------------------------------------------------------------
  Subroutine conmodel_getconint(ptr,coords,v,u,f,all_flag)
    Type(ConstraintInfo), Pointer              :: ptr
    Type(VecType), Dimension(:), Intent(In)    :: coords,v
    Type(VecType), Dimension(:), Intent(InOut) :: f
    Real(kind=RDbl), Intent(Out)               :: u
    Logical, Intent(In)                        :: all_flag

    !** maps forces, also does evansmorris constraints, which is built into 
    !** this vesrion of ciccotti integration.
    If ((all_flag).And.(Associated(ptr%conparams%ru))) Then
      Call ciccotti_mapforces(ptr%conparams%ru,coords,v,f)
    End If

    !** This is a redundant call in case of ciccotti integration if all of
    !** the bond length constraints are inside the Ciccotti rigid units.
    !** It increased the CPU time by 14% on the case that I checked
    If (Associated(ptr%conparams%em)) Then
      Call evansmorriss_getconstraints(ptr%conparams%em,coords,v,u,f)
    End If

  End Subroutine conmodel_getconint

  !----------------------------------------------------------------------------
  ! Calculates the penalty for the constraint model and adjusts parameters
  !----------------------------------------------------------------------------
  Subroutine conmodel_getpenalty(ptr,coords,v,sumobjf)
    Type(ConstraintInfo), Pointer              :: ptr
    Type(VecType), Dimension(:), Intent(InOut) :: coords, v
    Real(Kind=RDbl), Intent(Out)               :: sumobjf

    If (Associated(ptr%conparams%em)) Then
      Call evansmorriss_getpenalty(ptr%conparams%em,coords,v,sumobjf)
    End If
    
  End Subroutine conmodel_getpenalty


  !----------------------------------------------------------------------------
  ! Do necessary adjustments AFTER the integration for a WHOLE molecule
  ! Requires:  conparams -- constraint model information pointer
  !            coords -- principal (non-simcell) coordinate vectors
  !            v -- corresponding veclocity vectors
  !----------------------------------------------------------------------------
  Subroutine conmodel_postadjust(conparams,coords,v)
    Type(ConstraintModel), Intent(In)          :: conparams
    Type(VecType), Dimension(:), Intent(InOut) :: coords, v

    If (Associated(conparams%ru)) Then
      Call ciccotti_regenSet(conparams%ru,coords,v)
    End If
    
  End Subroutine conmodel_postadjust

  !--------------------------------------------------------------
  ! Checks if the CON interaction is fast or slow
  !--------------------------------------------------------------
  Logical Function conmodel_isfastcon(constraint)
    Type(ConstraintInfo) :: constraint

    conmodel_isfastcon = constraint%fast

  End Function conmodel_isfastcon

  !--------------------------------------------------------------------
  ! Determines if the interaction should be evaluated
  ! Requires:  ptr -- pointer to interaction parameters structure
  !            fast -- proposed evaluation speed
  !--------------------------------------------------------------------
  Logical Function conmodel_eval(ptr,fast)
    Type(ConstraintInfo), Pointer    :: ptr
    Logical, Intent(In)              :: fast

    conmodel_eval = .False.
    If (.Not. Associated(ptr)) Return
    conmodel_eval = (ptr%fast .And. fast)

  End Function conmodel_eval

  !----------------------------------------------------------------------------
  ! Displays the information about the model
  !----------------------------------------------------------------------------
  Subroutine conmodel_display(conparams,unit,indent)
    Type(ConstraintModel), Intent(In) :: conparams
    Integer, Intent(In) :: unit, indent
    Character(len=indent) :: blank
    Integer :: i

    blank = ""
    Do i = 1, indent
      blank = blank//" "
    End Do

    If (Associated(conparams%em)) Then
      Call evansmorriss_display(conparams%em,unit,indent+2)
    End If

    If (Associated(conparams%ru)) Then
      Call ciccotti_display(conparams%ru,unit,indent+2)
    End If

  End Subroutine conmodel_display

  !----------------------------------------------------------------------------
  ! Display information about the constraints
  !----------------------------------------------------------------------------
  Subroutine conmodel_condisplay(ptr,spc,unit,indent)
    Type(ConstraintInfo), Pointer  :: ptr
    Integer, Intent(In)            :: spc, unit, indent 

    Integer               :: i
    Character(len=indent) :: blank
    Character(len=strLen) :: name

    blank = Repeat(' ',indent)

    name = molecules_name(spc)

    Write(unit,'(3a)') blank,"Constraint information for ",Trim(name)
    Write(unit,'(2x,3a)') blank,"Model name : ",ptr%model
    Call conmodel_display(ptr%conparams,unit,indent+2)

  End Subroutine conmodel_condisplay

  !----------------------------------------------------------------------------
  ! Cleanup the constraints model
  !----------------------------------------------------------------------------
  Subroutine conmodel_cleanup(conparams)
    Type(ConstraintModel),Intent(inout) :: conparams
    Integer :: error

    If (Associated(conparams%em)) Then
      Call evansmorriss_cleanup(conparams%em)
      Deallocate(conparams%em,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'conparams%em'"
        Stop
      End If
    End If

  End Subroutine conmodel_cleanup

  !------------------------------------------
  ! Cleanup constraint info
  !------------------------------------------
  Subroutine conmodel_concleanup(constraint)
    Type(ConstraintInfo) :: constraint

    Call conmodel_cleanup(constraint%conparams)
    
  End Subroutine conmodel_concleanup



End Module conmodel






