!---------------------------------------------------------------
! This module defines the data structure of a 3-D coordinate and 
! the operations that accompany it. All operations are defined on 
! vectors with components of type real.  Vector of type integer
! is also provided including assignment operations that convert
! back and forth between the real and integer types
!----------------------------------------------------------------
Module genvector

  Use defaults, Only: RDbl, strLen

  Implicit None
  Save

  Private
  Public :: GenVecType, IntGenVecType, Assignment(=), Operator(+), &
      Operator(-), Operator(*), Operator(/), genvector_display

  ! When allocatable arrays are allowed in
  ! structures we can make this module more
  ! general.
  Integer, Parameter  :: MAX_GENVEC_SIZE = 30

  Type IntGenVecType
    Integer, Dimension(MAX_GENVEC_SIZE)   :: comp
    Integer    :: nsize
  End Type IntGenVecType

  Type GenVecType
    Real(kind=RDbl), Dimension(MAX_GENVEC_SIZE) :: comp
    Integer          :: nsize
  End Type GenVecType

  Interface Assignment(=)
    Module Procedure genvector_vectype_eq_vectype
    Module Procedure genvector_vectype_eq_intvectype
!!$    Module Procedure genvector_intvectype_eq_vectype
    Module Procedure genvector_vec_eq_scalar
    Module Procedure genvector_vectype_eq_rarray
    Module Procedure genvector_vectype_eq_intarray
    Module Procedure genvector_rarray_eq_vectype
  End Interface

  Interface Operator(+)
    Module Procedure genvector_add
  End Interface

  Interface Operator(-)
    Module Procedure genvector_subtract
  End Interface

  Interface Operator(*)
    Module Procedure genvector_dotprod
    Module Procedure genvector_vecxscalar
  End Interface

  Interface Operator(/)
    Module Procedure genvector_intdivide
    Module Procedure genvector_realdivide
  End Interface

  Interface genvector_display
!    Module Procedure genvector_filedisplay
    Module Procedure genvector_strdisplay
  End Interface

  Interface unitvec
    Module Procedure genvector_getunitvec
  End Interface

  Interface mag
    Module Procedure genvector_getnorm
  End Interface

  Interface isinplane
    Module Procedure genvector_isingenplane
  End Interface

  Interface getplanedist
    Module Procedure genvector_getgenplanedist
  End Interface

Contains
  !--------------------------------------------------------
  ! Initialize the vector "vec1".  This routine just sets
  ! the size of "vec1" to "nsize"
  !--------------------------------------------------------
  Subroutine genvector_init(vec1, nsize)
    Type(GenVecType), Intent(inout) :: vec1
    Integer, Intent(in)             :: nsize

    Integer    :: i

    vec1%nsize = nsize
    Do i=1, nsize
      vec1%comp(i) = 0.0_RDbl
    End Do
  End Subroutine genvector_init

  !---------------------------------------------------
  ! Assign one vector type "vec2" to "vec1"
  !---------------------------------------------------
  Subroutine genvector_vectype_eq_vectype(vec1, vec2)
    Type(GenVecType), Intent(out) :: vec1
    Type(GenVecType), Intent(in) :: vec2

    Integer     :: i

    vec1%nsize = vec2%nsize
    Do i=1, vec2%nsize
      vec1%comp(i) = vec2%comp(i)
    End Do
  End Subroutine genvector_vectype_eq_vectype
  
  !------------------------------------------------------------
  ! Assign an array to a real vector type
  !------------------------------------------------------------
  Subroutine genvector_vectype_eq_rarray(vec1, array1)
    Type(GenVecType), Intent(out) :: vec1
    Real(kind=RDbl), Dimension(:), Intent(in) :: array1
    Integer   :: arraysize, i
    
    arraysize = Size(array1, 1)
    If (arraysize > MAX_GENVEC_SIZE) Then
      Write(0,'(a,i4)') 'Array is too big.  Vector size is fixed at :', &
          MAX_GENVEC_SIZE
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If
    vec1%nsize = arraysize
    Do i=1, arraysize
      vec1%comp(i) = array1(i)
    Enddo
    Return
  End Subroutine genvector_vectype_eq_rarray

  !------------------------------------------------------------
  ! Assign an array to a integer vector type
  !------------------------------------------------------------
  Subroutine genvector_vectype_eq_intarray(vec1, array1)
    Type(IntGenVecType), Intent(out) :: vec1
    Integer, Dimension(:), Intent(in) :: array1
    Integer   :: arraysize, i
    
    arraysize = Size(array1, 1)
    If (arraysize > MAX_GENVEC_SIZE) Then
      Write(0,'(a,i4)') 'Array is too big.  Vector size is fixed at :', &
          MAX_GENVEC_SIZE
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If
    vec1%nsize = arraysize
    Do i=1, arraysize
      vec1%comp(i) = array1(i)
    Enddo
    Return
  End Subroutine genvector_vectype_eq_intarray

  !----------------------------------------------
  ! Assign a real vector type to an array 
  !----------------------------------------------
  Subroutine genvector_rarray_eq_vectype(array1, vec1)
    Real(kind=RDbl), Dimension(:), Intent(out) :: array1
    Type(GenVecType), Intent(in) :: vec1

    Integer   :: arraysize, i, nsize
    
    arraysize = Size(array1, 1)
    If (arraysize < vec1%nsize) Then
      Write(0,'(a,i4)') 'Array is too small.  Vector size is :', vec1%nsize
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If
    Do i=1, vec1%nsize
      array1(i) = vec1%comp(i)
    Enddo
    Return
  End Subroutine genvector_rarray_eq_vectype
  
  !----------------------------------------------------------
  ! Assign an integer vector type to a real vector type
  !----------------------------------------------------------
  Subroutine genvector_vectype_eq_intvectype(vec1, intvec2)
    Type(GenVecType), Intent(out)   :: vec1
    Type(IntGenVecType), Intent(in) :: intvec2
    Integer           :: i

    Do i = 1, intvec2%nsize
      vec1%comp(i) = intvec2%comp(i)
    End Do
    vec1%nsize = intvec2%nsize
    Return
  End Subroutine genvector_vectype_eq_intvectype

  !----------------------------------------------------------
  ! Assign an real vector type to an integer vector type
  !----------------------------------------------------------
  Subroutine genvector_intvectype_eq_vectype(intvec1, vec2)
    Type(IntGenVecType), Intent(out)    :: intvec1
    Type(GenVecType), Intent(in)  :: vec2
    Integer           :: i

    Do i = 1, vec2%nsize
      intvec1%comp(i) = vec2%comp(i)
    End Do
    intvec1%nsize = vec2%nsize
    Return
  End Subroutine genvector_intvectype_eq_vectype

  !-------------------------------------------
  ! Initialize a vector's components to "num1"
  !-------------------------------------------
  Subroutine genvector_vec_eq_scalar(vec1, num1)
    Type(GenVecType), Intent(inout) :: vec1
    Real(kind=RDbl), Intent(in) :: num1
    Integer   :: i
    
    Do i=1, vec1%nsize
      vec1%comp(i) = num1
    End Do
    Return
  End Subroutine genvector_vec_eq_scalar

  !--------------------------------------------------
  ! Add two vector structures
  !--------------------------------------------------
  Type(GenVecType) Function genvector_add(vec1, vec2)
    Type(GenVecType), Intent(in)  :: vec1, vec2
    Integer               :: i
    
    Do i=1, vec1%nsize
      genvector_add%comp(i) = vec1%comp(i) + vec2%comp(i)
    Enddo
    genvector_add%nsize = vec1%nsize
    Return
  End Function genvector_add

  !--------------------------------------------------
  ! Subtract two vector structures
  !--------------------------------------------------
  Type(GenVecType) Function genvector_subtract(vec1, vec2)
    Type(GenVecType), Intent(in)  :: vec1, vec2
    Integer               :: i
    
    Do i=1, vec2%nsize
      genvector_subtract%comp(i) = vec1%comp(i) - vec2%comp(i)
    Enddo
    genvector_subtract%nsize = vec1%nsize
    Return
  End Function genvector_subtract

  !--------------------------------------------------
  ! Get the product of a scalar with a vector
  !--------------------------------------------------
  Type(GenVecType) Function genvector_vecxscalar(vec1, num1)
    Type(GenVecType), Intent(in)  :: vec1
    Real(kind=RDbl), Intent(in)  :: num1
    Integer      :: i
    
    Do i=1, vec1%nsize
      genvector_vecxscalar%comp(i) = vec1%comp(i)*num1
    End Do
    genvector_vecxscalar%nsize = vec1%nsize
    Return
  End Function genvector_vecxscalar

  !--------------------------------------------------
  ! Get the dot product of two vectors
  !--------------------------------------------------
  Real(kind=RDbl) Function genvector_dotprod(vec1, vec2)
    Type(GenVecType), Intent(in)  :: vec1, vec2
    Real(kind=RDbl)  :: dotprod
    Integer               :: i

    dotprod = 0.0_RDbl
    Do i=1, vec2%nsize
      dotprod = dotprod + vec1%comp(i)*vec2%comp(i)
    Enddo
    genvector_dotprod = dotprod
    Return
  End Function genvector_dotprod
 
  !--------------------------------------------------
  ! Divide a vector by an int
  !--------------------------------------------------
  Type(GenVecType) Function genvector_intdivide(vec1, int1)
    Type(GenVecType), Intent(in)  :: vec1
    Integer, Intent(in)        :: int1
    Integer   :: i
    
    Do i=1, vec1%nsize
      genvector_intdivide%comp(i) = vec1%comp(i)/int1
    End Do
    genvector_intdivide%nsize = vec1%nsize
    Return
  End Function genvector_intdivide

  !--------------------------------------------------
  ! Divide a vector by a real
  !--------------------------------------------------
  Type(GenVecType) Function genvector_realdivide(vec1, real1)
    Type(GenVecType), Intent(in)  :: vec1
    Real(kind=RDbl), Intent(in) :: real1
    Integer      :: i
    
    Do i=1, vec1%nsize
      genvector_realdivide%comp(i) = vec1%comp(i)/real1
    Enddo
    genvector_realdivide%nsize = vec1%nsize
    Return
  End Function genvector_realdivide

  !----------------------------------------------------------
  ! Gets the unit vector in the direction of "vec"
  !----------------------------------------------------------
  Type(GenVecType) Function genvector_getunitvec(vec)
    Type(GenVecType), Intent(in) :: vec
    Real(kind=RDbl) :: magnitude

    magnitude = genvector_getnorm(vec)
    genvector_getunitvec = genvector_realdivide(vec,magnitude)
    genvector_getunitvec%nsize = vec%nsize

    Return
  End Function genvector_getunitvec

  !-----------------------------------------------------
  ! Get the norm square of a vector
  !-----------------------------------------------------    
  Real(kind=RDbl) Function genvector_getnormsq(vec1)
    Type(GenVecType), Intent(in)  :: vec1
    Real(kind=RDbl)   :: normsq
    Integer             :: i

    normsq = 0.0_RDbl
    
    Do i=1, vec1%nsize
      normsq = normsq + vec1%comp(i)**2
    End Do

    genvector_getnormsq = normsq
    Return
  End Function genvector_getnormsq

  !-----------------------------------------------------
  ! Get the norm of a vector "vec1"
  !-----------------------------------------------------
  Real(kind=RDbl) Function genvector_getnorm(vec1)
    Type(GenVecType), Intent(in) :: vec1
    Real(kind=RDbl)    :: norm

    genvector_getnorm =  Sqrt(genvector_getnormsq(vec1))
    Return
  End Function genvector_getnorm

  !-----------------------------------------------------
  ! Gets the nth component  of vec 
  !-----------------------------------------------------
  Real(kind=RDbl) Function genvector_getcomp_n(vec1,n)
    Type(GenVecType), Intent(in) :: vec1
    Integer ,Intent(in):: n
    genvector_getcomp_n =  vec1%comp(n)
  End Function genvector_getcomp_n

  !--------------------------------------------------------------
  ! Get the distance between two coordinates "vec1" and "vec2"
  !--------------------------------------------------------------
  Real(kind=RDbl) Function genvector_getdist(vec1, vec2)
    Type(GenVecType), Intent(IN) :: vec1, vec2
    Type(GenVecType)   :: temp
    
    temp = vec1 - vec2
    genvector_getdist = genvector_getnorm(temp)

    Return
  End Function genvector_getdist

  !-----------------------------------------------------
  ! Get the squared distance between two coordinates
  !-----------------------------------------------------
  Real(kind=RDbl) Function genvector_getdistsq(vec1, vec2)
    Type(GenVecType), Intent(IN) :: vec1, vec2
    Type(GenVecType)             :: temp
    
    temp = vec1 - vec2
    genvector_getdistsq = genvector_getnormsq(temp)

    Return
  End Function genvector_getdistsq

  !---------------------------------------------------------------------
  ! Performs Gram-Schmidt orthogonalization on the array of 
  ! vectors "vec1" to generate "orthovecs".  See 
  ! "Linear Algebra and its Applications" by Strang, Pg. 172, 3rd Ed.
  !---------------------------------------------------------------------
  Subroutine genvector_gsortho(vecs, orthovecs)
    Type(GenVecType), Dimension(:), Intent(in) :: vecs
    Type(GenVecType), Dimension(:), Intent(out):: orthovecs

    Integer    :: nvecs, i, j
    Real(kind=RDbl) :: dotprod

    !** The no. of vectors should be the same as the dimension of each vector
    nvecs = vecs(1)%nsize
    Do i=1, nvecs
      orthovecs(i) = vecs(i)
      Do j=1, i-1
        dotprod = orthovecs(i)*orthovecs(j)
        orthovecs(i) = orthovecs(i) - orthovecs(j)*dotprod
      End Do
      If (genvector_getnormsq(orthovecs(i)) < 1.0e-8_RDbl) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Initial Vectors are not linearly independent "
        Stop
      Else
        orthovecs(i) = genvector_getunitvec(orthovecs(i))        
      End If
    End Do
  End Subroutine genvector_gsortho

  !--------------------------------------------------------------
  ! This function checks if the points in the array "xyzcoords"
  ! lie in the plane normal to "normal" going throught the point
  ! "pt"
  !--------------------------------------------------------------
  Logical Function genvector_isingenplane(xyzcoords, pt, normal, opttol)
    Type(GenVecType), Intent(in):: xyzcoords
    Type(GenVecType), Intent(in)      :: pt, normal
    Real(kind=RDbl), Optional, Intent(in)  :: opttol

    Type(GenVecType)       :: inplanevec
    Real(kind=Rdbl)        :: dotprod, tol

    If (Present(opttol)) Then
      tol = opttol
    Else
      tol = 1.0e-5_RDbl
    End If

    inplanevec = xyzcoords - pt
      
    ! If inplanevec is actually in the plane then its
    ! dot product with the normal should be zero
    dotprod = inplanevec*normal
    If (Abs(dotprod) > tol) Then
      genvector_isingenplane = .False.
      Return
    End If
    genvector_isingenplane = .True.
    Return
  End Function genvector_isingenplane

  !-----------------------------------------------------------------
  ! This function gets the distance of the point "pt1" to the
  ! plane defined by the point "planept" and normal "planenormal"
  ! along the normal "normal1". 
  !-----------------------------------------------------------------
  Real(kind=RDbl) Function genvector_getgenplanedist(pt1, normal1, &
      planept, planenormal)
    Type(GenVecType), Intent(in) :: pt1, normal1, planept, planenormal
    
    genvector_getgenplanedist = planenormal*(planept - pt1)/(planenormal*normal1)
  End Function genvector_getgenplanedist

  !----------------------------------------------------------
  ! Write out the components of the vector to unit "unitno"
  !----------------------------------------------------------
  Subroutine genvector_filedisplay(vec1, unitno)
    Type(GenVecType), Intent(in) :: vec1
    Integer, Optional, Intent(in) :: unitno
    Integer        :: i, unit, vec_size
    Character(len=strLen) :: strformat

    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    vec_size = vec1%nsize
    If (vec_size < 10) Then
      Write(strformat, '(a,i1,2a)') "(1x, ",vec_size, "f8.3)"
    Else
      Write(strformat, '(a,i2,a)') "(1x, ",vec_size, "f8.3)"
    End If

    Write(unit,strformat) (vec1%comp(i), i=1, vec_size)
  End Subroutine genvector_filedisplay

  !--------------------------------------------------------
  ! This function returns a string with the components of
  ! the vector formatted according to "fmt".
  !--------------------------------------------------------
!!$  Character(len=strLen) Function genvector_strdisplay(vec1, fmt)
  Function genvector_strdisplay(vec1, fmt)
    Character(len=strLen) :: genvector_strdisplay
    Type(GenVecType), Intent(in)  :: vec1
    Character(*), Intent(in)   :: fmt
    Character(len=strLen)      :: strformat
    Integer                    :: i, vec_size

    strformat = fmt
    
    vec_size = vec1%nsize
    If (vec_size < 10) Then
      Write(strformat, '(a,i1,2a)') "(1x, ",vec_size, Trim(fmt),")"
    Else
      Write(strformat, '(a,i2,2a)') "(1x, ",vec_size, Trim(fmt),")"
    End If
    Write(genvector_strdisplay, strformat) (vec1%comp(i), i=1, vec_size)
    genvector_strdisplay = Adjustl(genvector_strdisplay)

  End Function genvector_strdisplay

End Module genvector
