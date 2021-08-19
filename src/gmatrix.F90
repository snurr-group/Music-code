!-------------------------------------------------------------------------
! This module handles operations on a generic matrix type.  It uses
! a type definition that contains a pointer to a 2D array.  This allows
! flexible matrix manipulations.
!
! Currently included manipulations include:
!   matrix multiplication, subtraction and addition, assignment from array
!   identity, transpose, normalization, eigeninfo generation
!   simularity transform, projector-based removal of vector components
!
! NOTE: since the type uses pointer, the matrices must be allocated
! and deallocated.  Hopefully we can improve this when the Fortran2000
! standard is fully adopted and allocatables can be used in data types.
!-------------------------------------------------------------------------
Module gmatrix

  Use defaults, Only: strLen,RDbl,zeroTolerance
  Use utils, Only: int2str,real2str,allocErrdisplay,deallocErrdisplay

  Implicit None
  Save

  Private
  Public :: GenericMatrix,gmatrix_init,gmatrix_size,gmatrix_free, &
      gmatrix_checkinit, gmatrix_setmtx,gmatrix_iden,gmatrix_eigeninfo, &
      gmatrix_display,gmatrix_simxform,gmatrix_einfo,gmatrix_copy, &
      Operator(*),Operator(+),Operator(-),gmatrix_normcols,gmatrix_clean, &
      gmatrix_projectout,gmatrix_isorthonorm,gmatrix_makeorthonorm, &
      gmatrix_issymm

  Type GenericMatrix
    Logical               :: symmetric,orthonorm
    Integer               :: nrows,ncols
    Real(kind=RDbl), Dimension(:,:), Pointer  :: mtx
  End Type GenericMatrix

!** Cannot use a equality assignment here because it conflicts with
!** the presence and usage of functions that return a gmatrix
!  Interface Assignment(=)
!    Module Procedure gmatrix_equateinit
!  End Interface

  Interface Operator(*)
    Module Procedure gmatrix_m_mult_m
  End Interface

  Interface Operator(+)
    Module Procedure gmatrix_m_plus_m
  End Interface

  Interface Operator(-)
    Module Procedure gmatrix_m_minus_m
  End Interface

Contains

  !-----------------------------------------------------------------------
  ! Initialize a matrix to a specified size, allocates only if necessary
  !-----------------------------------------------------------------------
  Subroutine gmatrix_size(matrix,nrows,ncols)
    Type(GenericMatrix), Intent(InOut) :: matrix
    Integer, Intent(In)                :: nrows,ncols

    Integer               :: error

    If (.Not.gmatrix_checkinit(matrix)) Then
      matrix%nrows = nrows
      matrix%ncols = ncols
      Allocate(matrix%mtx(nrows,ncols), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Else If ((matrix%nrows /= nrows).Or.(matrix%ncols /= ncols)) Then
      Call gmatrix_free(matrix)
      matrix%nrows = nrows
      matrix%ncols = ncols
      Allocate(matrix%mtx(nrows,ncols), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    
    End If

    matrix%mtx = 0.0_RDbl
    
  End Subroutine gmatrix_size

  !-----------------------------------------------------------------------
  ! Initialize a matrix to a specified size
  !-----------------------------------------------------------------------
  Type(GenericMatrix) Function gmatrix_init(nrows,ncols)
    Integer, Intent(In)   :: nrows,ncols

    Integer               :: error

    gmatrix_init%nrows = nrows
    gmatrix_init%ncols = ncols
    gmatrix_init%symmetric = .False.
    gmatrix_init%orthonorm = .False.

    Allocate(gmatrix_init%mtx(nrows,ncols), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    gmatrix_init%mtx = 0.0_RDbl
    
  End Function gmatrix_init

  !-----------------------------------------------------------------------
  ! Deallocate a matrix 
  !-----------------------------------------------------------------------
  Subroutine gmatrix_free(matrix)
    Type(GenericMatrix), Intent(InOut) :: matrix

    Integer               :: error

    Deallocate(matrix%mtx, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)    
    
  End Subroutine gmatrix_free

  !-----------------------------------------------------------------------
  ! Copy a matrix
  !-----------------------------------------------------------------------
  Type(GenericMatrix) Function gmatrix_copy(matrix)
    Type(GenericMatrix), Intent(In) :: matrix

    gmatrix_copy = gmatrix_init(matrix%nrows,matrix%ncols)
    gmatrix_copy%symmetric = matrix%symmetric
    gmatrix_copy%orthonorm = matrix%orthonorm
    gmatrix_copy%mtx = matrix%mtx
    
  End Function gmatrix_copy

#ifdef THISTOO   !conflicts with use of functions
  !-----------------------------------------------------------------------
  ! Copy a matrix
  !-----------------------------------------------------------------------
  Subroutine gmatrix_equateinit(new,old)
    Type(GenericMatrix), Intent(Out) :: new
    Type(GenericMatrix), Intent(In)  :: old

    new = gmatrix_init(old%nrows,old%ncols)
    new%symmetric = old%symmetric
    new%orthonorm = old%orthonorm
    
  End Subroutine gmatrix_equateinit
#endif

  !-----------------------------------------------------------------------
  ! Check the initialization of a matrix
  !-----------------------------------------------------------------------
  Logical Function gmatrix_checkinit(matrix)
    Type(GenericMatrix), Intent(In) :: matrix

    gmatrix_checkinit = .FALSE.

    If (Associated(matrix%mtx)) gmatrix_checkinit = .TRUE.
    
  End Function gmatrix_checkinit

  !-----------------------------------------------------------------------
  ! Produce an identity matrix
  ! Requires:  dim -- dimension of matrix to create
  !-----------------------------------------------------------------------
  Type(GenericMatrix) Function gmatrix_iden(dim)
    Integer, Intent(In)   :: dim

    Integer               :: i,error

    gmatrix_iden%nrows = dim
    gmatrix_iden%ncols = dim
    gmatrix_iden%orthonorm = .TRUE.
    gmatrix_iden%symmetric = .TRUE.

    Allocate(gmatrix_iden%mtx(dim,dim), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    gmatrix_iden%mtx = 0.0_RDbl
    Do i = 1,dim
      gmatrix_iden%mtx(i,i) = 1.0_RDbl
    End Do
    
  End Function gmatrix_iden

  !-------------------------------------------------------------------------
  ! Check a passed matrix for orthonormality.  This entails looping through
  ! the possible combinations of two columns, taking their dot products and
  ! comparing with some tolerance (optional).  It also checks that the 
  ! column vectors are normalized.
  ! Requires: matrix -- generalized matrix structure
  !           tolerance -- expected tolerance on orthogonality
  !-------------------------------------------------------------------------
  Logical Function gmatrix_isorthonorm(matrix,tolerance)
    Type(GenericMatrix), Intent(In)        :: matrix
    Real(kind=RDbl), Intent(In), Optional  :: tolerance

    Integer               :: i,j,k,n
    Real(kind=RDbl)       :: dotprod,tol,norm

    If (Present(tolerance)) Then
      tol = tolerance
    Else
      tol = zeroTolerance
    End If

    n = 0
    gmatrix_isorthonorm = .True.

    Do i = 1,matrix%ncols
      norm = 0.0_RDbl
      Do j = 1,matrix%nrows
        norm = norm + matrix%mtx(i,j)**2
      End Do
      If (Abs(1.0_RDbl - Sqrt(norm)) > tol) Then
        Write(*,*) __FILE__," : ",__LINE__, &
            'non-normalized column: ',i,j,norm
        gmatrix_isorthonorm = .False.
      End If
    End Do

    Do i = 1,matrix%ncols
      Do j = i+1,matrix%ncols
        dotprod = 0.0_RDbl
        Do k = 1,matrix%nrows
          dotprod = dotprod + matrix%mtx(k,i)*matrix%mtx(k,j) 
        End Do
        If (dotprod > tol) Then
          Write(*,*) __FILE__," : ",__LINE__, &
              'non-orthogonal column pair: ',i,j,dotprod
          n = n + 1
          gmatrix_isorthonorm = .False.
        End If
      End Do
    End Do

  End Function gmatrix_isorthonorm

  !-------------------------------------------------------------------------
  ! Check a passed matrix for symmetry.
  ! Requires: matrix -- generalized matrix structure
  !           tolerance -- expected tolerance on element difference
  !-------------------------------------------------------------------------
  Logical Function gmatrix_issymm(matrix,tolerance)
    Type(GenericMatrix), Intent(In)        :: matrix
    Real(kind=RDbl), Intent(In), Optional  :: tolerance

    Integer               :: i,j
    Real(kind=RDbl)       :: diff,tol

    gmatrix_issymm = .True.

    If (Present(tolerance)) Then
      tol = tolerance
    Else
      tol = zeroTolerance
    End If

    Do i = 1,matrix%nrows
      Do j = i+1,matrix%ncols
        diff = matrix%mtx(i,j) - matrix%mtx(j,i)
        If (Abs(diff) > tol) Then
          Write(*,*) __FILE__," : ",__LINE__, &
              'non-matching pair: ',i,j,diff
          gmatrix_issymm = .False.
        End If
      End Do
    End Do

  End Function gmatrix_issymm

  !-----------------------------------------------------------------------
  ! Create a generic matrix using an input 2D array
  ! Requires: array -- Real 2D array of any size
  !           orthonorm -- logical indication of orthonormality (optional)
  !-----------------------------------------------------------------------
  Type(GenericMatrix) Function gmatrix_setmtx(array,orthonorm)
    Real(kind=RDbl), Dimension(:,:), Intent(In) :: array
    Logical, Intent(In), Optional               :: orthonorm

    Integer               :: error

    gmatrix_setmtx%symmetric = .False.
    If (Present(orthonorm)) Then
      gmatrix_setmtx%orthonorm = orthonorm
    Else
      gmatrix_setmtx%orthonorm = .False.
    End If

    gmatrix_setmtx%nrows = Size(array,1)
    gmatrix_setmtx%ncols = Size(array,2)

    Allocate(gmatrix_setmtx%mtx(gmatrix_setmtx%nrows, &
        gmatrix_setmtx%ncols), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    gmatrix_setmtx%mtx = array
    
  End Function gmatrix_setmtx

  !--------------------------------------------------------------------------
  ! Find the transpose of a matrix (not adjoint).  matrix(m,n) -> matrix(n,m)
  ! Simply copies the matrix if it is flagged as being symmetric.  
  ! Requires: matrix -- matrix to be transposed
  !--------------------------------------------------------------------------
  Type(GenericMatrix) Function gmatrix_transpose(matrix)
    Type(GenericMatrix), Intent(In)      :: matrix

    Integer             :: i,j

    If (matrix%symmetric) Then
      gmatrix_transpose = gmatrix_copy(matrix)
    Else
      gmatrix_transpose = gmatrix_init(matrix%ncols,matrix%nrows)
      Do i = 1,matrix%nrows
        Do j = 1,matrix%ncols
          gmatrix_transpose%mtx(j,i) = matrix%mtx(i,j)
        End Do
      End Do
    End If
    
  End Function gmatrix_transpose

  !-----------------------------------------------------------------------
  ! Invert a matrix
  ! Requires: matrix -- matrix to be inverted
  !-----------------------------------------------------------------------
  Type(GenericMatrix) Function gmatrix_invert(matrix)
    Type(GenericMatrix), Intent(In)      :: matrix

    Type(GenericMatrix)    :: temp

    If (matrix%orthonorm) Then
      gmatrix_invert = gmatrix_transpose(matrix)
    Else
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " No provision for non-orthogonal matrix inversion"
      Stop
    End If

#if DEBUG
    !** just a check
    temp = gmatrix_invert*matrix
    Write(*,'(2a,i4)') __FILE__," This is U^(-1)*U (inversion check)",__LINE__
    Call gmatrix_display(temp,'f8.3',2,6)
    Call gmatrix_free(temp)
#endif
    
  End Function gmatrix_invert

  !------------------------------------------------------------------------
  ! Performs matrix multiplication: resultant(l,n) = m1(l,m)*m2(m,n).
  ! Note that it creates and allocates a NEW matrix
  ! Requires:  m1 -- matrix with size (l,m)
  !            m2 -- matrix with size (m,n)
  !------------------------------------------------------------------------
  Type(GenericMatrix) Function gmatrix_m_mult_m(m1, m2)
    Type(GenericMatrix), Intent(In) :: m1, m2

    Integer               :: l,m,n

    !** make sure the multiplication makes sense 
    If (m1%ncols /= m2%nrows) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Matrix sizes do not permit multiplication"
      Write(0,'(3x,4(a,i3),a)') '(',m1%nrows,',',m1%ncols,')   (', &
          m2%nrows,',',m2%ncols,')'
      Stop      
    End If
    
    gmatrix_m_mult_m = gmatrix_init(m1%nrows,m2%ncols)

    Do l = 1,m1%nrows
      Do n = 1,m2%ncols
        Do m = 1,m1%ncols
          gmatrix_m_mult_m%mtx(l,n) = gmatrix_m_mult_m%mtx(l,n) + &
              m1%mtx(l,m) * m2%mtx(m,n)
        End Do
      End Do
    End Do

    Return

  End Function gmatrix_m_mult_m

  !------------------------------------------------------------------------
  ! Performs matrix ADDITION
  ! Note that it creates and allocates a NEW matrix
  ! Requires:  m1 -- matrix with size (m,n)
  !            m2 -- matrix with size (m,n)
  !------------------------------------------------------------------------
  Type(GenericMatrix) Function gmatrix_m_plus_m(m1, m2)
    Type(GenericMatrix), Intent(In) :: m1, m2

    Integer               :: m,n

    !** make sure the multiplication makes sense 
    If ((m1%nrows /= m2%nrows).Or.(m1%ncols /= m2%ncols)) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Matrix sizes do not permit subtraction"
      Write(0,'(3x,4(a,i3),a)') '(',m1%nrows,',',m1%ncols,')   (', &
          m2%nrows,',',m2%ncols,')'
      Stop      
    End If
    
    gmatrix_m_plus_m = gmatrix_init(m1%nrows,m1%ncols)

    Do m = 1,m1%nrows
      Do n = 1,m1%ncols
        gmatrix_m_plus_m%mtx(m,n) = m1%mtx(m,n) + m2%mtx(m,n)
      End Do
    End Do

  End Function gmatrix_m_plus_m

  !------------------------------------------------------------------------
  ! Performs matrix SUBTRACTION
  ! Note that it creates and allocates a NEW matrix
  ! Requires:  m1 -- matrix with size (m,n)
  !            m2 -- matrix with size (m,n)
  !------------------------------------------------------------------------
  Type(GenericMatrix) Function gmatrix_m_minus_m(m1, m2)
    Type(GenericMatrix), Intent(In) :: m1, m2

    Integer               :: m,n

    !** make sure the multiplication makes sense 
    If ((m1%nrows /= m2%nrows).Or.(m1%ncols /= m2%ncols)) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Matrix sizes do not permit subtraction"
      Write(0,'(3x,4(a,i3),a)') '(',m1%nrows,',',m1%ncols,')   (', &
          m2%nrows,',',m2%ncols,')'
      Stop      
    End If
    
    gmatrix_m_minus_m = gmatrix_init(m1%nrows,m1%ncols)

    Do m = 1,m1%nrows
      Do n = 1,m1%ncols
        gmatrix_m_minus_m%mtx(m,n) = m1%mtx(m,n) - m2%mtx(m,n)
      End Do
    End Do

  End Function gmatrix_m_minus_m

  !------------------------------------------------------------------------
  ! Performs matrix multiplication: resultant(l,n) = m1(l,m)*m2(m,n).
  ! Note that it uses an already allocated matrix
  ! Requires:  m1 -- matrix with size (l,m)
  !            m2 -- matrix with size (m,n)
  !------------------------------------------------------------------------
  Subroutine gmatrix_mtxmult(m1,m2,resultant)
    Type(GenericMatrix), Intent(In)    :: m1, m2
    Type(GenericMatrix), Intent(InOut) :: resultant

    Integer               :: l,m,n

    !** make sure the multiplication makes sense 
    If (m1%ncols /= m2%nrows) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Matrix sizes do not permit multiplication"
      Write(0,'(3x,4(a,i3),a)') '(',m1%nrows,',',m1%ncols,')   (', &
          m2%nrows,',',m2%ncols,')'
      Stop      
    End If

    !** make sure the passed resultant has the correct size
    If ((resultant%nrows /= m1%nrows).Or.(resultant%ncols /= m2%ncols)) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Incoming resultant matrix does not have correct size"
      Write(0,'(3x,4(a,i3),a)') '(',resultant%nrows,',',resultant%ncols, &
          ') should be  (', m1%nrows,',',m2%ncols,')'
      Stop      
    End If

    resultant%mtx = 0.0_RDbl

    Do l = 1,m1%nrows
      Do n = 1,m2%nrows
        Do m = 1,m1%ncols
          resultant%mtx(l,n) = resultant%mtx(l,n) + &
              m1%mtx(l,m) * m2%mtx(m,n)
        End Do
      End Do
    End Do

  End Subroutine gmatrix_mtxmult

  !------------------------------------------------------------------------
  ! Normalizes the individual column vectors that form a matrix
  ! Requires:  matrix -- matrix to be operated on
  !------------------------------------------------------------------------
  Subroutine gmatrix_normcols(matrix)
    Type(GenericMatrix), Intent(InOut) :: matrix

    Integer            :: i,j
    Real(kind=RDbl)    :: norm

    Do j = 1,matrix%ncols
      norm = 0.0_RDbl
      Do i = 1,matrix%nrows
        norm = norm + matrix%mtx(i,j)**2
      End Do
      norm = 1.0_RDbl/Sqrt(norm)
      Do i = 1,matrix%nrows
        matrix%mtx(i,j) = matrix%mtx(i,j)*norm
      End Do
    End Do

  End Subroutine gmatrix_normcols

  !------------------------------------------------------------------------
  ! Makes the column vectors of a matrix orthonormal using a simple
  ! procedure.  1st: normalizes them  2nd: apply (Gram-Schmidt?) technique,
  ! v'_i = v_i - (v_k .dot. v_i)*v_k for k<i,  This may not work on some
  ! systems, see Numerical Recipes for better methods.
  !------------------------------------------------------------------------
  Subroutine gmatrix_gramschmidt(matrix)
    Type(GenericMatrix), Intent(InOut) :: matrix

    Integer            :: i,j,k
    Real(kind=RDbl)    :: dotprod

    Call gmatrix_normcols(matrix)

    Do j = 2,matrix%ncols

      Do k = 1,j-1
        dotprod = 0.0_RDbl
        Do i = 1,matrix%nrows
          dotprod = dotprod + matrix%mtx(i,k)*matrix%mtx(i,j)
        End Do

        Do i = 1,matrix%nrows
          matrix%mtx(i,j) = matrix%mtx(i,j) - dotprod*matrix%mtx(i,k)
        End Do
      End Do

      dotprod = 0.0_RDbl
      Do i = 1,matrix%nrows
        dotprod = dotprod + matrix%mtx(i,j-1)*matrix%mtx(i,j)
      End Do
      Write(*,*) j,dotprod
    End Do

  End Subroutine gmatrix_gramschmidt

  !------------------------------------------------------------------------
  ! Makes the column vectors of a matrix orthonormal using SVD algorithm
  !------------------------------------------------------------------------
  Subroutine gmatrix_makeorthonorm(matrix)
    Type(GenericMatrix), Intent(InOut) :: matrix

    Integer            :: m,n
    Real(kind=RDbl), Dimension(matrix%ncols)              :: singularvalues
    Real(kind=RDbl), Dimension(matrix%ncols,matrix%ncols) :: vmatrix

    m = matrix%nrows
    n = matrix%ncols
!    Call svdcmp(matrix%mtx,m,n,m,n,singularvalues,vmatrix)
    Write (*,*) "This requires an algorithm from numerical recipes"
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  End Subroutine gmatrix_makeorthonorm

  !---------------------------------------------------------------------------
  ! Get the eigenvectors and eigenvalues of a matrix
  ! Requires:  matrix -- the matrix to process
  !            evectors -- eigenvectors in columns of matrix (pre-associated)
  !            evalues -- eigenvalues of matrix
  !---------------------------------------------------------------------------
  Subroutine gmatrix_eigeninfo(matrix,evectors,evalues)
    Type(GenericMatrix), Intent(In)            :: matrix
    Type(GenericMatrix), Intent(InOut)         :: evectors
    Real(kind=RDbl), Dimension(:), Intent(Out) :: evalues

    Integer            :: error,dim
    Real(kind=RDbl), Dimension(matrix%nrows)   :: subdiagonal

    If (matrix%nrows /= matrix%ncols) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Passed matrix is not square"
      Stop
    End If    

    dim = matrix%nrows

    !** Convert to a tridiagonal form
    Call tred2(dim,dim,matrix%mtx,evalues,subdiagonal,evectors%mtx)

    !** solve the tridiagonal form for the eigenvalues and eigenvectors
    Call tql2(dim,dim,evalues,subdiagonal,evectors%mtx,error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not diagonalize the matrix"
      Stop
    End If

    evectors%orthonorm = .True.

  End Subroutine gmatrix_eigeninfo

  !---------------------------------------------------------------------------
  ! Get the eigenvectors and eigenvalues of an array
  ! Requires:  matrix -- the matrix to process (2D array)
  !            evectors -- eigenvectors in columns of matrix (pre-associated)
  !            evalues -- eigenvalues of matrix
  ! NOTE: this routine should NOT be here, just a test
  !---------------------------------------------------------------------------
  Subroutine gmatrix_einfo(matrix,evectors,evalues)
    Real(kind=RDbl), Dimension(:,:), Intent(In)  :: matrix
    Real(kind=RDbl), Dimension(:,:), Intent(Out) :: evectors
    Real(kind=RDbl), Dimension(:), Intent(Out)   :: evalues

    Integer            :: error,dim
    Real(kind=RDbl), Dimension(Size(matrix,1))   :: subdiagonal

    If (Size(matrix,1) /= Size(matrix,2)) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Passed matrix is not square"
      Stop
    End If    

    dim = Size(matrix,1)

    !** Convert to a tridiagonal form
    Call tred2(dim,dim,matrix,evalues,subdiagonal,evectors)

    !** solve the tridiagonal form for the eigenvalues and eigenvectors
    Call tql2(dim,dim,evalues,subdiagonal,evectors,error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not diagonalize the matrix"
      Stop
    End If

  End Subroutine gmatrix_einfo

  !-----------------------------------------------------------------------
  ! Perform a similarity transform on a matrix: resultant = U^(-1)*M*U
  ! Requires: matrix -- the matrix to be operated on
  !           umatrix -- the operating matrix
  ! NOTE: may wish to rewrite this routine so it does the whole 
  !       thing at once, this would probably be faster
  !-----------------------------------------------------------------------
  Subroutine gmatrix_simxform(matrix,umatrix)
    Type(GenericMatrix), Intent(InOut)       :: matrix
    Type(GenericMatrix), Intent(In)          :: umatrix

    Type(GenericMatrix)        :: temp,uinv

    temp = matrix*umatrix
    uinv = gmatrix_invert(umatrix)

    Call gmatrix_mtxmult(uinv,temp,matrix)

    Call gmatrix_free(uinv)    
    Call gmatrix_free(temp)    

  End Subroutine gmatrix_simxform

  !--------------------------------------------------------------------------
  ! Project out a set of vectors from a matrix 
  !   projection_matrix = P = I - F*F^(T)
  !   resultant = P*M*P, where M = object matrix
  ! Requires: object -- the matrix to be operated on
  !           modifier -- matrix containing projection vectors in columns (F)
  ! NOTE: the modifier MUST be orthonormal
  !--------------------------------------------------------------------------
  Subroutine gmatrix_projectout(object,modifier)
    Type(GenericMatrix), Intent(InOut)       :: object
    Type(GenericMatrix), Intent(InOut)       :: modifier

    Type(GenericMatrix)        :: projector,identity,transpose,temp,dyad

    identity = gmatrix_iden(object%nrows) 
    transpose = gmatrix_transpose(modifier)
    dyad = modifier*transpose
    projector = identity - dyad

    temp = object*projector

    Call gmatrix_mtxmult(projector,temp,object)

    Call gmatrix_free(identity)
    Call gmatrix_free(transpose)
    Call gmatrix_free(dyad)
    Call gmatrix_free(projector)
    Call gmatrix_free(temp)
    
  End Subroutine gmatrix_projectout

  !-----------------------------------------------------------------------
  ! Displays the data in the generic matrix structure
  ! Requires: matrix -- generalized matrix to display
  !           fmt -- format of entries, if 'a', will pick the shortest 
  !           indent -- number of spaces from left margin
  !           unit -- unit to dump into
  !-----------------------------------------------------------------------
  Subroutine gmatrix_display(matrix,fmt,indent,unit)
    Type(GenericMatrix), Intent(In)  :: matrix
    Character(*), Intent(In)         :: fmt
    Integer, Intent(In)              :: indent,unit

    Integer                     :: i,j,k,lo,hi
    Integer                     :: max_line_length
    Character(len=strLen)       :: strformat,string
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    max_line_length = 10   
    If (Trim(fmt) == 'a') Then
      string = int2str(max_line_length)
      Write(strformat,'(a,3a)') '(2x,a,i5,',Trim(string),Trim(fmt),')'
    Else
      string = int2str(max_line_length)
      Write(strformat,'(a,3a)') '(2x,a,i5,',Trim(string),Trim(fmt),')'
    End If

    Write(unit,'(2a,2i4)') blank," matrix dimensions (row,col): ", &
        matrix%nrows,matrix%ncols

    If (Associated(matrix%mtx)) Then
      Do i = 1,matrix%nrows
        Do j = 1,matrix%ncols,max_line_length
          lo = j
          hi = Min((j+max_line_length-1),matrix%ncols)
          If (Trim(fmt) == 'a') Then
            string = real2str(matrix%mtx(i,k),1)
            Write(unit,strformat) blank,i,(Trim(string)//' ',k=lo,hi)
          Else
            Write(unit,strformat) blank,i,(matrix%mtx(i,k),k=lo,hi)
          End If
        End Do
      End Do
    End If

  End Subroutine gmatrix_display

  !-----------------------------------------------------------------------
  ! Cleans the generic matrix structure
  ! Requires: matrix -- generalized matrix structure
  !-----------------------------------------------------------------------
  Subroutine gmatrix_clean(matrix)
    Type(GenericMatrix), Intent(InOut)  :: matrix

    Integer               :: error

    If (Associated(matrix%mtx)) Then
      Deallocate(matrix%mtx, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'mtx')
    End If

  End Subroutine gmatrix_clean

End Module gmatrix
