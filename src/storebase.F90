!------------------------------------------------------------------------------
! This module provides the basic energy plus derivatives storage structure.
! From this structure is the storage of the forcefield built.  The structure
! is also used in many low-level forcefield routines to return results.  This 
! module also supplies routines for manipulating this basic structure.
!
! Description of contents:
!   nrg --      scalar non-intramolecular energy 
!   force --    vector force
!   hess --     tensor (3x3) 2nd derivatives
!   intranrg -- set of scalar intramolecular energies
! NOTE: that the intramolecular and "normal" energies are stored separately,
! in other words, %nrg does NOT include the intramolecular energy, however,
! %force DOES include the intramolecular forces.  Maybe this is ideologically
! inconsistent?
!
! Important routines:
!   storebase_init -- initilizes and zeros the basic structure 
!   storebase_zero -- zeros the basic structure 
!   storebase_nderivs -- returns the number of derivatives associated
!   storebase_inc -- performs various forms of incrementing based on the call
!   storebase_sum -- sums over an array of basic structures
!   storebase_disp -- returns a string with basic structure values
!------------------------------------------------------------------------------

Module storebase

  Use defaults, Only: RDbl,strLen,xlstrLen,NO_OF_INTRA_POTS,TOTAL_INDEX
  Use utils, Only: toreal,allocErrDisplay,deallocErrDisplay,real2str,int2str
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getdistsq 

  Implicit None

  Private
  Public :: EnergyPlus, storebase_scalarmult, storebase_init, storebase_inc, &
      storebase_nderivs, storebase_zero, storebase_sum, storebase_sumarray, &
      storebase_initcopy, storebase_copy, storebase_nrg, storebase_null, &
      storebase_forces, storebase_hess, storebase_incAllArray, storebase_sumincarray, &
      storebase_display, storebase_disp, storebase_clean, &
      storebase_sumintra, storebase_incinv, storebase_chkintra, &
      storebase_iszero, storebase_scalenrg, storebase_chkequal, &
      storebase_projectout, storebase_project, storebase_totnrg, &
      storebase_subtract, storebase_add

  !** Might be wise to use VecType and MatrixType for force and hess?
  Type EnergyPlus
    Real(kind=RDbl)                          :: nrg     !* energy
    Real(kind=RDbl), Dimension(:),   Pointer :: force   !* forces
    Real(kind=RDbl), Dimension(:), Pointer :: hess    !* hessian (2nd deriv)
    Real(kind=RDbl), Dimension(:), Pointer   :: intranrg
  End Type EnergyPlus

  Interface storebase_inc
    Module Procedure storebase_incAll
    Module Procedure storebase_incAllArray
    Module Procedure storebase_incNrg
    Module Procedure storebase_incIntra
    Module Procedure storebase_incForce
   Module Procedure storebase_incHess
  End Interface

  Interface storebase_incinv
    Module Procedure storebase_incinvAll
  End Interface

  Interface storebase_clean
    Module Procedure storebase_cleanSingle
    Module Procedure storebase_cleanArray
  End Interface 

Contains

  !----------------------------------------------------------------------------
  ! Initializes the basic Energy and Derivatives sub-structure
  ! Requires:  base -- Energy and Derivatives structure
  !            nderivs -- number of derivatives to store (0,1,2)
  !            intra -- flags initialization of intra potentials, default False
  !----------------------------------------------------------------------------
  Subroutine storebase_init(base,nderivs,intra)
    Type(EnergyPlus), Intent(Out)          :: base
    Integer, Intent(In)                    :: nderivs
    Logical, Intent(In), Optional          :: intra

    Integer     :: error

    !** nullify the pointers in the structure
    Call storebase_null(base)

    !** initialize gradients if necessary
    If (nderivs > 0) Then
      Allocate(base%force(3),stat=error)    
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'base%force')
    End If

    !** initialize 2nd derivatives if necessary
    If (nderivs > 1) Then
      Allocate(base%hess(4),stat=error)    
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'base%hess')
    End If

    !** initialize intra-molecular interaction storage if necessary
    If (Present(intra)) Then
      If (intra) Then
        Allocate(base%intranrg(NO_OF_INTRA_POTS),stat=error)    
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'base%intranrg')
      End If
    End If

    !** Zero everything
    Call storebase_zero(base)

  End Subroutine storebase_init

  !----------------------------------------------------------------------------
  ! Initialize and copy a basic Energy and Derivatives sub-structure
  ! Requires:  base -- Energy and Derivatives structure
  !            old -- Energy and Derivatives structure
  !----------------------------------------------------------------------------
  Subroutine storebase_initcopy(base,old)
    Type(EnergyPlus), Intent(Out)          :: base
    Type(EnergyPlus), Intent(In)           :: old

    Integer     :: i,error

    base%nrg = old%nrg

    If (Associated(old%force)) Then
      Allocate(base%force(3),stat=error)    
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'base%force')
      base%force = old%force
    Else
      Nullify(base%force)
    End If

    If (Associated(old%hess)) Then
      Allocate(base%hess(4),stat=error)    
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'base%hess')
      base%hess = old%hess
    Else
      Nullify(base%hess)
    End If

    If (Associated(old%intranrg)) Then
      Allocate(base%intranrg(NO_OF_INTRA_POTS),stat=error)    
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'base%intranrg')
      Do i = 1,NO_OF_INTRA_POTS
        base%intranrg(i) = old%intranrg(i)
      End Do
    Else
      Nullify(base%intranrg)
    End If

  End Subroutine storebase_initcopy

  !----------------------------------------------------------------------------
  ! Nullify the pointers in the storage base
  ! Requires:  base -- base structure to nullify
  !----------------------------------------------------------------------------
  Subroutine storebase_null(base)
    Type(EnergyPlus), Intent(InOut)          :: base

    Nullify(base%force)
    Nullify(base%hess)
    Nullify(base%intranrg)

  End Subroutine storebase_null

  !----------------------------------------------------------------------------
  ! Copy a basic Energy and Derivatives sub-structure.  Assumes that all the
  ! fields in the original structure are initialized in the destination.
  ! Requires:  dest -- Energy and Derivatives structure (must be initialized)
  !            orig -- Energy and Derivatives structure 
  !----------------------------------------------------------------------------
  Subroutine storebase_copy(dest,orig)
    Type(EnergyPlus), Intent(InOut)       :: dest
    Type(EnergyPlus), Intent(In)          :: orig

    dest%nrg = orig%nrg

    If (Associated(orig%force)) Then
      dest%force = orig%force
    End If

    If (Associated(orig%hess)) Then
      dest%hess = orig%hess
    End If

    If (Associated(orig%intranrg)) Then
      dest%intranrg = orig%intranrg
    End If

  End Subroutine storebase_copy

  !----------------------------------------------------------------------------
  ! Zeros the basic Energy and Derivatives structure
  ! Requires:  base -- Energy and Derivatives structure
  !----------------------------------------------------------------------------
  Subroutine storebase_zero(base)
    Type(EnergyPlus), Intent(InOut)          :: base

    base%nrg = 0.0_RDbl
    If (Associated(base%force)) base%force = 0.0_RDbl
    If (Associated(base%hess)) base%hess = 0.0_RDbl
    If (Associated(base%intranrg)) base%intranrg = 0.0_RDbl

  End Subroutine storebase_zero

  !----------------------------------------------------------------------------
  ! Scales the entire structure by a scalar
  ! Requires:  base -- Energy and Derivatives structure
  !            factor -- factor to scale energy by
  !----------------------------------------------------------------------------
  Subroutine storebase_scalarmult(base,factor)
    Type(EnergyPlus), Intent(InOut)          :: base
    Real(kind = RDbl), Intent(In)            :: factor

    base%nrg = factor*base%nrg
    If (Associated(base%force)) base%force = factor*base%force
    If (Associated(base%hess)) base%hess = factor*base%hess
    If (Associated(base%intranrg)) base%intranrg = factor*base%intranrg

  End Subroutine storebase_scalarmult

  !----------------------------------------------------------------------------
  ! Projects out a specified component of the force vector. 
  ! Requires:  base -- Energy and Derivatives structure
  !            vec -- unit vector along which to remove force component
  !----------------------------------------------------------------------------
  Subroutine storebase_projectout(base,vec)
    Type(EnergyPlus), Intent(InOut)          :: base
    Type(VecType), Intent(In)                :: vec

    If (.Not. Associated(base%force)) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' force component not initialized, unexpected'
      Stop
    End If

    base%force = base%force - base%force*vec

  End Subroutine storebase_projectout

  !----------------------------------------------------------------------------
  ! Projects the force vector onto a specified unit vector.  Keep only the
  ! component of the force vector in the direction of the given unit vector.
  ! Requires:  base -- Energy and Derivatives structure
  !            vec -- unit vector to project force vector with
  !----------------------------------------------------------------------------
  Subroutine storebase_project(base,vec)
    Type(EnergyPlus), Intent(InOut)          :: base
    Type(VecType), Intent(In)                :: vec

    If (.Not. Associated(base%force)) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' force component not initialized, unexpected'
      Stop
    End If

    base%force = (base%force * vec) * vec

  End Subroutine storebase_project

  !----------------------------------------------------------------------------
  ! Add two base structures together and return the result.  This initializes
  ! a new structure at the same level as "base1"
  ! Requires:  base1 -- 1st Energy and Derivatives structure
  !            base2 -- 2nd Energy and Derivatives structure
  !----------------------------------------------------------------------------
  Type(EnergyPlus) Function storebase_sum(base1,base2)
    Type(EnergyPlus), Intent(In)    :: base1,base2

    Call storebase_init(storebase_sum,storebase_nderivs(base1))
    storebase_sum%nrg = base1%nrg + base2%nrg
    If (Associated(base1%force)) storebase_sum%force = base1%force + base2%force
    If (Associated(base1%hess)) storebase_sum%hess = base1%hess + base2%hess
    If (Associated(base1%intranrg)) Then
      storebase_sum%intranrg = base1%intranrg + base2%intranrg
    End If

  End Function storebase_sum

  !----------------------------------------------------------------------------
  ! Subtract the second base structure from the first.  Assumes that all
  ! associated entries in the 1st are also associated in the second
  ! Requires:  base1 -- 1st Energy and Derivatives structure
  !            base2 -- 2nd Energy and Derivatives structure
  !----------------------------------------------------------------------------
  Subroutine storebase_subtract(base1,base2)
    Type(EnergyPlus), Intent(InOut) :: base1
    Type(EnergyPlus), Intent(In)    :: base2

    base1%nrg = base1%nrg - base2%nrg
    If (Associated(base1%force)) base1%force = base1%force - base2%force
    If (Associated(base1%hess)) base1%hess = base1%hess - base2%hess
    If (Associated(base1%intranrg)) Then
      base1%intranrg = base1%intranrg - base2%intranrg
    End If

  End Subroutine storebase_subtract

  !----------------------------------------------------------------------------
  ! Add the second base structure from the first.  Assumes that all
  ! associated entries in the 1st are also associated in the second
  ! Requires:  base1 -- 1st Energy and Derivatives structure
  !            base2 -- 2nd Energy and Derivatives structure
  !----------------------------------------------------------------------------
  Subroutine storebase_add(base1,base2)
    Type(EnergyPlus), Intent(InOut) :: base1
    Type(EnergyPlus), Intent(In)    :: base2

    base1%nrg = base1%nrg + base2%nrg
    If (Associated(base1%force)) base1%force = base1%force + base2%force
    If (Associated(base1%hess)) base1%hess = base1%hess + base2%hess
    If (Associated(base1%intranrg)) Then
      base1%intranrg = base1%intranrg + base2%intranrg
    End If

  End Subroutine storebase_add

  !----------------------------------------------------------------------------
  ! Sum a collection of base structures
  ! Requires:  base -- Energy and Derivatives structure
  !            total -- output summed structure, must be pre-initialized
  !----------------------------------------------------------------------------
  Subroutine storebase_sumarray(base,total)
    Type(EnergyPlus), Dimension(:), Intent(In)  :: base
    Type(EnergyPlus), Intent(InOut)             :: total

    Integer          :: i

    Call storebase_zero(total)
    Do i = 1,Size(base)
      Call storebase_incAll(total,base(i))
    End Do

  End Subroutine storebase_sumarray

  !----------------------------------------------------------------------------
  ! Sum a collection of base structures and use to increment another
  ! Requires:  base -- Energy and Derivatives structure
  !            total -- output summed structure, must be pre-initialized
  !----------------------------------------------------------------------------
  Subroutine storebase_sumincarray(base,total)
    Type(EnergyPlus), Dimension(:), Intent(In)  :: base
    Type(EnergyPlus), Intent(InOut)             :: total

    Integer          :: i

    Do i = 1,Size(base)
      Call storebase_incAll(total,base(i))
    End Do

  End Subroutine storebase_sumincarray

  !----------------------------------------------------------------------------
  ! Returns the highest associated derivative pointer in the storage structure
  ! Requires:  base -- Energy and Derivatives structure
  !----------------------------------------------------------------------------
  Integer Function storebase_nderivs(base)
    Type(EnergyPlus), Intent(In)           :: base

    If (Associated(base%hess)) Then
      storebase_nderivs = 2
    Else If (Associated(base%force)) Then
      storebase_nderivs = 1
    Else
      storebase_nderivs = 0
    End If

  End Function storebase_nderivs

  !----------------------------------------------------------------------------
  ! Returns True if the intra component is initialized
  ! Requires:  base -- Energy and Derivatives structure
  !----------------------------------------------------------------------------
  Logical Function storebase_chkintra(base)
    Type(EnergyPlus), Intent(In)           :: base

    If (Associated(base%intranrg)) Then
      storebase_chkintra = .True.
    Else
      storebase_chkintra = .False.
    End If

  End Function storebase_chkintra

  !----------------------------------------------------------------------------
  ! Returns the energy
  ! Requires:  base -- Energy and Derivatives structure
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function storebase_nrg(base)
    Type(EnergyPlus), Intent(In)           :: base

    storebase_nrg = base%nrg

  End Function storebase_nrg

  !----------------------------------------------------------------------------
  ! Return the total energy (%nrg + %intranrg(TOTAL))
  ! Requires:  base -- Energy and Derivatives structure
  !            resum -- if True, resum the total intramolecular energy
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function storebase_totnrg(base,resum)
    Type(EnergyPlus), Intent(In)           :: base
    Logical, Intent(In)                    :: resum

    Integer           :: i
    Real(kind=RDbl)   :: intra

    If (resum) Then
      intra = 0.0_RDbl
      Do i = 1,NO_OF_INTRA_POTS
        If (i == TOTAL_INDEX) Cycle
        intra = intra + base%intranrg(i)
      End Do
    Else
      intra = base%intranrg(TOTAL_INDEX)
    End If

    storebase_totnrg = base%nrg + intra

  End Function storebase_totnrg

  !----------------------------------------------------------------------------
  ! Return an array of forces from a passed array of basic structures
  ! Requires:  base -- array of Energy and Derivatives structures
  !----------------------------------------------------------------------------
  Function storebase_forces(base)
    Type(EnergyPlus), Dimension(:), Intent(In)    :: base
    Type(VecType), Dimension(Size(base))          :: storebase_forces

    Integer              :: i

    Do i = 1,Size(base)
      storebase_forces(i) = base(i)%force
    End Do

  End Function storebase_forces
  !----------------------------------------------------------------------------
  ! Return an array of 2nd derivatives from a passed array of basic structures
  ! Requires:  base -- array of Energy and Derivatives structures
  !----------------------------------------------------------------------------
 Logical  Function storebase_hess(base,hhot)
    Type(EnergyPlus), Dimension(:), Intent(In)    :: base
     Real(kind=RDbl), Dimension(:), Intent(Out) :: hhot
   
 
    Integer              :: i, j
    storebase_hess=.False.

!!$L Something has to be written here

!    Do i = 1,size(hhot)
!     hhot = base%hess(i)
!    End Do
    
   storebase_hess=.True.

  End Function storebase_hess



  !----------------------------------------------------------------------------
  ! INCREMENT a basic storage type with another of the same type.  Assumes 
  ! that all the fields in the 'add' structure are initialized in the 
  ! 'base' structure
  ! Requires:  base -- Energy and Derivatives structure
  !            add -- base structure to add
  !----------------------------------------------------------------------------
  Subroutine storebase_incAll(base,add)
    Type(EnergyPlus), Intent(InOut)           :: base
    Type(EnergyPlus), Intent(In)              :: add

!LC    Write(*,*) base%nrg,add%nrg
    base%nrg = base%nrg + add%nrg

    If (Associated(add%force)) Then
      base%force(1) = base%force(1) + add%force(1)
      base%force(2) = base%force(2) + add%force(2)
      base%force(3) = base%force(3) + add%force(3)
    End If

    If (Associated(add%hess)) Then
      base%hess = base%hess + add%hess
    End If

    If (Associated(add%intranrg)) Then
      base%intranrg = base%intranrg + add%intranrg
    End If

  End Subroutine storebase_incAll

  !----------------------------------------------------------------------------
  ! Add the "inverse" of a basic storage type to another of the same type. 
  ! This means that scalars add and the forces subtract.
  ! Requires:  base -- Energy and Derivatives structure
  !            minus -- base structure to subtract
  !----------------------------------------------------------------------------
  Subroutine storebase_incinvAll(base,minus)
    Type(EnergyPlus), Intent(InOut)           :: base
    Type(EnergyPlus), Intent(In)              :: minus

    base%nrg = base%nrg + minus%nrg

    If (Associated(base%force)) Then
      base%force = base%force - minus%force
    End If

    !** not sure about this (subtract?)
    If (Associated(base%hess)) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Stop
      base%hess = base%hess - minus%hess
    End If

    If (Associated(base%intranrg)) Then
      base%intranrg = base%intranrg + minus%intranrg
    End If

  End Subroutine storebase_incinvAll

  !----------------------------------------------------------------------------
  ! INCREMENT an array of basic storage types with another of the same type.
  ! The arrays MUST be the same size.
  ! Requires:  base -- Energy and Derivatives structure
  !            add -- base structure to add
  !----------------------------------------------------------------------------
  Subroutine storebase_incAllArray(base,add)
    Type(EnergyPlus), Dimension(:), Intent(InOut)   :: base
    Type(EnergyPlus), Dimension(:), Intent(In)      :: add

    Integer               :: i

    Do i = 1,Size(base)
      Call storebase_incAll(base(i),add(i))
    End Do

  End Subroutine storebase_incAllArray

  !----------------------------------------------------------------------------
  ! INCREMENT an ENERGY value in the basic store type
  ! Requires:  base -- Energy and Derivatives structure
  !            nrg -- energy value to add
  !----------------------------------------------------------------------------
  Subroutine storebase_incNrg(base,nrg)
    Type(EnergyPlus), Intent(InOut)           :: base
    Real(kind=RDbl), Intent(In)               :: nrg

    base%nrg = base%nrg + nrg

  End Subroutine storebase_incNrg

  !----------------------------------------------------------------------------
  ! INCREMENT an INTRAMOLECULAR ENERGY value in the basic store type
  ! Requires:  base -- Energy and Derivatives structure
  !            nrg -- energy value to add
  !            index -- intramolecular energy index
  !----------------------------------------------------------------------------
  Subroutine storebase_incIntra(base,nrg,index)
    Type(EnergyPlus), Intent(InOut)           :: base
    Real(kind=RDbl), Intent(In)               :: nrg
    Integer, Intent(In)                       :: index

    base%intranrg(index) = base%intranrg(index) + nrg

  End Subroutine storebase_incIntra

  !----------------------------------------------------------------------------
  ! INCREMENT a FORCE value in the basic store type
  ! Requires:  base -- Energy and Derivatives structure
  !            force -- force vector to add
  !----------------------------------------------------------------------------
  Subroutine storebase_incForce(base,force)
    Type(EnergyPlus), Intent(InOut)           :: base
    Type(VecType), Intent(In)                 :: force

    base%force(1) = base%force(1) + force%comp(1)
    base%force(2) = base%force(2) + force%comp(2)
    base%force(3) = base%force(3) + force%comp(3)

  End Subroutine storebase_incForce

  !----------------------------------------------------------------------------
  ! INCREMENT HESS value in the basic store type
  ! Requires:  base -- Energy and Derivatives structure
  !            Hess -- 2nd derivatives array
  !----------------------------------------------------------------------------
  Subroutine storebase_incHess(base,hess)
    Type(EnergyPlus), Intent(InOut)           :: base
    Real(kind=RDbl), Dimension(Size(base%hess)), Intent(In)     :: hess

    base%hess(1) = base%hess(1) + hess(1)
    base%hess(2) = base%hess(2) + hess(2)
    base%hess(3) = base%hess(3) + hess(3)
    base%hess(4) = base%hess(4) + hess(4)

  End Subroutine storebase_incHess

  !----------------------------------------------------------------------------
  ! Sum the energies in the intramolecular set and put in total index
  ! Requires:  base -- Energy and Derivatives structure
  !            inctotal -- if True, increment the total energy also
  !----------------------------------------------------------------------------
  Subroutine storebase_sumintra(base,inctotal)
    Type(EnergyPlus), Intent(InOut)           :: base
    Logical, Intent(In)                       :: inctotal

    Integer           :: i

    base%intranrg(TOTAL_INDEX) = 0.0_RDbl
    Do i = 1,NO_OF_INTRA_POTS
      If (i == TOTAL_INDEX) Cycle
      base%intranrg(TOTAL_INDEX) = base%intranrg(TOTAL_INDEX) + base%intranrg(i)
    End Do

    If (inctotal) Then
      base%nrg = base%nrg + base%intranrg(TOTAL_INDEX)
    End If

  End Subroutine storebase_sumintra

  !----------------------------------------------------------------------------
  ! Scale just the energy
  ! Requires:  base -- Energy and Derivatives structure
  !            factor -- factor to scale energy by
  !----------------------------------------------------------------------------
  Subroutine storebase_scalenrg(base,factor)
    Type(EnergyPlus), Intent(InOut)          :: base
    Real(kind = RDbl), Intent(In)            :: factor

    base%nrg = base%nrg*factor

  End Subroutine storebase_scalenrg

  !----------------------------------------------------------------------------
  ! Returns a flag indicating if the storage contains non-zero numbers within
  ! a given tolerance.  Used for displaying only non-zero entries
  ! Requires:  base -- basic storage structure  
  !            tolerance -- tolerance for calling something zero
  !----------------------------------------------------------------------------
  Logical Function storebase_iszero(base,tolerance)
    Type(EnergyPlus), Intent(In)      :: base
    Real(kind=RDbl), Intent(In)       :: tolerance

    Integer               :: i,nderivs
    Real(kind=RDbl)       :: total

    total = Abs(base%nrg)
    nderivs = storebase_nderivs(base)

    !** Get the sum of the intra part 
    If (storebase_chkintra(base)) Then
      Do i = 1,NO_OF_INTRA_POTS
        total = total + Abs(base%intranrg(i))
      End Do
    End If
    
    !** Get the force part
    If (nderivs > 0) Then
      Do i = 1,3
        total = total + Abs(base%force(i))
      End Do
    End If

    !** Get the hessian part
    If (nderivs > 1) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Stop
    End If
    
    !** Compare versus tolerance and return without writing if too small
    storebase_iszero = .True.
    If (total > tolerance) storebase_iszero = .False.

  End Function storebase_iszero

  !----------------------------------------------------------------------------
  ! Compares all the initialized entries in two structures.  For an entry that
  ! doesn't match its analog in the other structure, it adds a note to the 
  ! returned string.
  ! Requires:  base1 -- 1st basic storage structure  
  !            base2 -- 2nd basic storage structure  
  !            tolerance -- tolerance for calling equal
  !            output -- string containing notes when entries are not equal
  !----------------------------------------------------------------------------
  Logical Function storebase_chkequal(base1,base2,tolerance,output)
    Type(EnergyPlus), Intent(In)      :: base1,base2
    Real(kind=RDbl), Intent(In)       :: tolerance
    Character(*), Intent(InOut)       :: output

    Integer               :: i,j,nderivs
    Real(kind=RDbl)       :: diff
    Character(len=strLen) :: string1,string2

    storebase_chkequal = .True.

    diff = base1%nrg - base2%nrg
    If (Abs(diff) > tolerance) Then
      string1 = real2str(diff)
      Write(output,'(3a)') Trim(output),' NRG: ',Trim(string1)
      storebase_chkequal = .False.
    End If

    If (Associated(base1%force)) Then
      If (.Not. Associated(base2%force)) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            " Force component of second storage base structure not associated"
        Stop        
      End If
      diff = 0.0_RDbl
      Do i = 1,3
        diff = diff + (base1%force(i) - base2%force(i))**2
      End Do
      If (diff > tolerance) Then
        string1 = real2str(diff)
        Write(output,'(3a)') Trim(output),' FRC: ',Trim(string1)
        storebase_chkequal = .False.
      End If
    End If 

    If (Associated(base1%hess)) Then
      If (.Not. Associated(base2%hess)) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            " Hessian component of second storage base structure not associated"
        Stop        
      End If
      diff = 0.0_RDbl
              Do j = 1,4
          diff = diff + (base1%hess(j) - base2%hess(j))**2
        End Do

      If (diff > tolerance) Then
        string1 = real2str(diff)
        Write(output,'(3a)') Trim(output),' HESS: ',Trim(string1)
        storebase_chkequal = .False.
      End If 
    End If 

    If (Associated(base1%intranrg)) Then
      If (.Not. Associated(base2%intranrg)) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            " intra component of second storage base structure not associated"
        Stop        
      End If
      Do i = 1,NO_OF_INTRA_POTS
        diff = base1%intranrg(i) - base2%intranrg(i)
        If (diff > tolerance) Then
          string1 = real2str(diff)
          string2 = int2str(i)
          Write(output,'(5a)') Trim(output),' INTRA(',Trim(string2), &
              '): ',Trim(string1)
          storebase_chkequal = .False.
        End If
      End Do
    End If 

  End Function storebase_chkequal

  !----------------------------------------------------------------------------
  ! Display the basic storage structure 
  ! Requires:  base -- basic storage structure
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine storebase_display(base,indent,unitno)
    Type(EnergyPlus), Intent(In)      :: base
    Integer, Intent(In)               :: indent
    Integer, Optional, Intent(In)     :: unitno

    Integer                           :: i,j,unit
    Character(len=indent)             :: blank
    Character(len=strLen)             :: string

    blank = Repeat(' ',indent)
    
    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    Write(unit,'(2a,f8.3)') blank,'Energy: ',base%nrg

    If (Associated(base%force)) Then
      Write(unit,'(2a,3f8.3)') blank,'Force: ',(base%force(i),i=1,3)
    End If 

    If (Associated(base%hess)) Then
      Do i = 1,4
        Write(unit,'(2a,3f8.3)') blank,'Hessian: ',(base%hess(i))
      End Do
    End If 

    If (Associated(base%intranrg)) Then
      Write(unit,'(a,a16,7(a5,4x))') blank,'Intramolecular: ', &
          'bs','bb','tor','ip','con','total'
      Write(unit,'(a,16x)',Advance='No') blank
      Do i = 1,NO_OF_INTRA_POTS
        string = real2str(base%intranrg(i),8)
        Write(unit,'(a9)',Advance='No') Trim(string)
      End Do
      Write(unit,*) 
    End If 
    
  End Subroutine storebase_display

  !----------------------------------------------------------------------------
  ! Returns the basic storage structure in a one-line format
  ! Requires:  base -- basic storage structure
  !----------------------------------------------------------------------------
  Function storebase_disp(base)
    Type(EnergyPlus), Intent(In)      :: base
    Character(len=xlstrLen)           :: storebase_disp

    Integer                           :: i,j
    Character(len=strLen)             :: string

    string = real2str(base%nrg,11)
    Write(storebase_disp,'(a,i1,2a)') 'nd=',storebase_nderivs(base), &
        ' NRG=',Trim(string)

    If (Associated(base%force)) Then
      Write(storebase_disp,'(2a,3e11.3)') Trim(storebase_disp), &
!          ' HACKFORCE=',(base%force(i)*100,i=1,3)
          ' FRC=',(base%force(i),i=1,3)
    End If 

    If (Associated(base%hess)) Then
      Write(storebase_disp,'(2a)') Trim(storebase_disp), &
          ' HESSIAN'
    End If 

    If (Associated(base%intranrg)) Then
!      Write(storebase_disp,'(2a)') Trim(storebase_disp), &
!          ' INTRA'
      Write(storebase_disp,'(2a,10f8.4)') Trim(storebase_disp), &
          ' INTRA=',(base%intranrg(i),i=1,NO_OF_INTRA_POTS)
!      Write(storebase_disp,'(2a,10e12.5)') Trim(storebase_disp), &
!          ' INTRA=',(base%intranrg(i),i=1,NO_OF_INTRA_POTS)
    End If 
    
  End Function storebase_disp

  !----------------------------------------------------------------------------
  ! Clean the basic storage structure 
  ! Requires:  base -- basic storage structure
  !----------------------------------------------------------------------------
  Subroutine storebase_cleanSingle(base)
    Type(EnergyPlus), Intent(InOut)  :: base

    Integer            :: error

    If (Associated(base%force)) Then
      Deallocate(base%force, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'force')
    End If 

    If (Associated(base%hess)) Then
      Deallocate(base%hess, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'hess')
    End If 

    If (Associated(base%intranrg)) Then
      Deallocate(base%intranrg, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'intranrg')
    End If 

  End Subroutine storebase_cleanSingle

  !----------------------------------------------------------------------------
  ! Clean an array of basic storage structures
  ! Requires:  base -- basic storage structures
  !----------------------------------------------------------------------------
  Subroutine storebase_cleanArray(base)
    Type(EnergyPlus), Dimension(:), Intent(InOut)  :: base

    Integer            :: i

    Do i = 1,Size(base)
      Call storebase_cleanSingle(base(i))
    End Do
    
  End Subroutine storebase_cleanArray


End module storebase

