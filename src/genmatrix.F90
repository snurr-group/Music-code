!----------------------------------------------------------
! This module defines various matrix operations.  Although
! I have tried to keep the matrix dimensions general a lot
! of routines work on only 3x3 matrices.
!----------------------------------------------------------
Module genmatrix

  Use genvector, Only: GenVecType
  Use defaults, Only: RDbl, strLen
  Implicit None
  Save

  Private
  Public :: Assignment(=), Operator(*), GenMatrixType, genmatrix_display

  Integer, Parameter   :: MAX_GENROWS = 30
  Integer, Parameter   :: MAX_GENCOLS = 30
  
  Type GenMatrixType
    Real(kind=RDbl), Dimension(MAX_GENROWS, MAX_GENCOLS)  :: comp
    Integer       :: ncols, nrows
  End Type GenMatrixType

  Interface Assignment(=)
    Module Procedure genmatrix_m_eq_real
    Module Procedure genmatrix_m_eq_2darray
  End Interface

  Interface Operator(*)
    Module Procedure genmatrix_m_mult_v
    Module Procedure genmatrix_v_mult_m
    Module Procedure genmatrix_m_mult_m
  End Interface

  Interface genmatrix_display
    Module Procedure genmatrix_display
  End Interface

Contains
  !-------------------------------------------------------
  ! This function initializes the matrix from a 2-D array
  !-------------------------------------------------------
  Subroutine genmatrix_m_eq_2darray(m1, array)
    Type(GenMatrixType), Intent(out) :: m1
    Real(kind=RDbl), Dimension(:,:), Intent(in)   :: array
    Integer  :: r, c, dim1, dim2

    dim1 = Size(array,1)
    dim2 = Size(array,2)
    Do c=1, dim2
      Do r=1, dim1
        m1%comp(r, c) = array(r, c)
      End Do
    End Do
    Return
  End Subroutine genmatrix_m_eq_2darray

  !-------------------------------------------------------
  ! This function assigns the scalar "r1" to each element 
  ! of the matrix "m1" 
  !-------------------------------------------------------
  Subroutine genmatrix_m_eq_real(m1, r1)
    Type(GenMatrixType), Intent(out) :: m1
    Real(kind=RDbl), Intent(in):: r1
    Integer     :: r, c
    
    Do c=1, m1%ncols
      Do r=1, m1%nrows
        m1%comp(r, c) = r1
      End Do
    End Do
  End Subroutine genmatrix_m_eq_real

  !--------------------------------------------------
  ! This function multiplies the matrix "m1" by the
  ! vector "vec1".  Could have written in terms of
  ! dot products but vectors and matrices are so
  ! intimate that they have access to each others
  ! private members.
  !--------------------------------------------------
  Type(GenVecType) Function genmatrix_m_mult_v(m1, vec1)
    Type(GenMatrixType), Intent(in) :: m1
    Type(GenVecType), Intent(in)    :: vec1
    Type(GenVecType)      :: temp
    Integer            :: r, c
    
    Call genvector_init(temp, vec1%nsize)
    Do c=1, m1%ncols
      Do r=1, m1%nrows
        temp%comp(r) = temp%comp(r) + m1%comp(r,c)*vec1%comp(c)
      End Do
    End Do
    genmatrix_m_mult_v = temp
    genmatrix_m_mult_v%nsize = m1%nrows
    Return
  End Function genmatrix_m_mult_v
  
  !---------------------------------------------------
  ! This function multiplies the vector "vec1" by the
  ! matrix "m1"
  !---------------------------------------------------
  Type(GenVecType) Function genmatrix_v_mult_m(vec1, m1)
    Type(GenMatrixType), Intent(in) :: m1
    Type(GenVecType), Intent(in)    :: vec1
    Type(GenVecType)      :: temp
    Integer            :: r, c
    
    Call genvector_init(temp, vec1%nsize)
    Do c=1, m1%ncols
      Do r=1, m1%nrows
        temp%comp(c) = temp%comp(c) + vec1%comp(r)*m1%comp(r,c)
      End Do
    End Do
    genmatrix_v_mult_m = temp
    genmatrix_v_mult_m%nsize = m1%ncols
    Return
  End Function genmatrix_v_mult_m

  !-------------------------------------------------
  ! Multiplies matrices "m1" and "m2"
  !-------------------------------------------------
  Type(GenMatrixType) Function genmatrix_m_mult_m(m1, m2)
    Type(GenMatrixType), Intent(in) :: m1, m2
    Type(GenMatrixType)   :: temp
    Integer            :: r, c, i
    
    temp = 0.0_RDbl
    Do c=1, m1%ncols
      Do r=1, m1%nrows
        Do i=1, m1%ncols
          temp%comp(r,c) = temp%comp(r,c) + m1%comp(r,i)*m2%comp(i,c)
        End Do
      End Do
    End Do
    genmatrix_m_mult_m = temp
    genmatrix_m_mult_m%nrows = m1%nrows
    genmatrix_m_mult_m%ncols = m2%ncols
    Return
  End Function genmatrix_m_mult_m

  !----------------------------------------------------------
  ! Gets the transpose of a matrix "m1"
  !----------------------------------------------------------
  Type(GenMatrixType) Function genmatrix_transpose(m1)
    Type(GenMatrixType), Intent(in)  :: m1
    Integer     :: r, c
    
    Do r=1, m1%nrows
      Do c=1, m1%ncols
        genmatrix_transpose%comp(c,r) = m1%comp(r,c)
      Enddo
    End Do
    genmatrix_transpose%nrows = m1%ncols
    genmatrix_transpose%ncols = m1%nrows
    Return
  End Function genmatrix_transpose

  !----------------------------------------------------------
  ! Gets the determinant of a matrix of dimension "size X size"
  ! The default size is 3.
  !----------------------------------------------------------
  Recursive Function genmatrix_getdet(m1, optsize) Result(det)
    Type(GenMatrixType), Intent(in)  :: m1
    Integer, Optional, Intent(in) :: optsize
    Real(kind=RDbl)      :: det   
    Type(GenMatrixType)     :: minor
    Integer              :: dim, c, sign
    
    If (Present(optsize)) Then
      dim = optsize
    Else
      dim = 3
    End If

    If (dim == 1) Then
      det = m1%comp(1,1)
      Return
    Endif

    det = 0.0_RDbl
    Do c=1, dim
      minor = genmatrix_getminor(m1, 1, c, dim)
      sign = (-1)**(1+c)
      det = det + m1%comp(1,c)*sign*genmatrix_getdet(minor, dim-1)
    End Do
    Return
  End Function genmatrix_getdet

  !-----------------------------------------------------
  ! Takes a matrix "m1" of dimension "dim X dim" and returns 
  ! the minor of element (row,col).
  ! The minor is matrix resulting from removing row "row" and
  ! and column "col"
  !-----------------------------------------------------
  Type(GenMatrixType) Function genmatrix_getminor(m1, row, col, dim)
    Type(GenMatrixType), Intent(in)  :: m1
    Integer, Intent(in)           :: row, col
    Integer, Intent(in)           :: dim
    Integer    :: r, c, rowno, colno

    rowno = 0
    genmatrix_getminor = 0.0_RDbl
    Do r=1, dim
      If (r == row) Cycle
      rowno = rowno + 1
      colno = 0
      Do c=1, dim
        If (c == col) Cycle
        colno = colno + 1
        genmatrix_getminor%comp(rowno, colno) = m1%comp(r, c)
      End Do
    End Do
    genmatrix_getminor%nrows = rowno
    genmatrix_getminor%ncols = colno
    Return
  End Function genmatrix_getminor

  !----------------------------------------------------------
  ! Gets the inverse of a matrix "m1"
  !----------------------------------------------------------
  Type(GenMatrixType) Function genmatrix_getinv(m1)
    Type(GenMatrixType), Intent(in)  :: m1
    Type(GenMatrixType) :: minor
    Real(kind=RDbl)  :: det, minordet
    Integer     :: r, c, dim, sign

    genmatrix_getinv = 0.0_RDbl
    dim = m1%nrows
    det = genmatrix_getdet(m1, dim)
    Do r=1, dim
      Do c=1, dim
        minor = genmatrix_getminor(m1, r, c, dim)
        minordet = genmatrix_getdet(minor, dim-1)
        sign = (-1)**(r+c)
        genmatrix_getinv%comp(c,r) = sign*minordet/det
      End Do
    End Do
    genmatrix_getinv%nrows = dim
    genmatrix_getinv%ncols = dim
    Return
  End Function genmatrix_getinv

  !-----------------------------------------------------------------------
  ! Generates a matrix to transform coodinates from basis "basis1" to
  ! "basis2".  The inverse transformation is just the transpose.
  !-----------------------------------------------------------------------
  Type(GenMatrixType) Function genmatrix_genbasistransmat(basis1, basis2)
    Type(GenVecType), Dimension(:), Intent(in) :: basis1, basis2

    Integer         :: nvecs, i, j
    Real(kind=RDbl) :: dotprod

    If (basis1(1)%nsize /= basis2(1)%nsize) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " The dimensions of the two bases are not the same"
      Stop
    End If
    nvecs = basis1(1)%nsize

    Do i=1, nvecs
      Do j=1, nvecs
        dotprod = basis2(i)*basis1(j)
        genmatrix_genbasistransmat%comp(i,j) = dotprod
      End Do
    End Do
    genmatrix_genbasistransmat%ncols = nvecs
    genmatrix_genbasistransmat%nrows = nvecs
  End Function genmatrix_genbasistransmat

  !------------------------------------------------------
  ! A quick and dirty routine to print the elements
  ! of the matrix "m1" out to unit "unitno".  The default
  ! unit is 6 (standard output)
  !------------------------------------------------------
  Subroutine genmatrix_display(m1, unitno)
    Type(GenMatrixType), Intent(in) :: m1
    Integer, Optional, Intent(in) :: unitno
    Integer  :: r, c, unit
    Character(len=strLen) :: strformat

    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    Endif

    If (m1%ncols < 10) Then
      Write(strformat, '(a,i1,2a)') "(1x, ",m1%ncols, "f8.3)"
    Else
      Write(strformat, '(a,i2,a)') "(1x, ",m1%ncols, "f8.3)"
    End If

    Do r=1, m1%nrows
      Write(unit, strformat) (m1%comp(r,c), c=1, m1%ncols)
    End Do
  End Subroutine genmatrix_display
  
End Module genmatrix

