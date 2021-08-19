!-------------------------------------------------------------------------
! This module handles the definition and utilization of a local cartesian
! coordinate system.  The functionality was previously in molecule.
!-------------------------------------------------------------------------

Module localcoords

  Use defaults, Only: RDbl, strLen, zeroTolerance, lstrLen
  Use utils, Only: split, stripcmnt, findint, tolower, getpath, filesrchstr, &
      toint, isfileopen, toupper, toreal, filesrchwocomment, real2str, int2str
  Use file, Only: file_getunit
  Use vector, Only: VecType, mag, vector_getnormsq, vector_display, &
      vector_getcomp, vector_iscollinear, vector_crossprod,&
      Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use matrix, Only: MatrixType, matrix_getinv, Assignment(=), Operator(*), &
      matrix_transpose, matrix_display, matrix_identity

  Implicit None
  Save

  Private
  Public :: LocalCoordInfo, localcoords_init, localcoords_getaxes, &
      localcoords_display, localcoords_clean

  !** Stores the information for determining the local coordinate
  !** system in terms of the center-of-mass coordinates of the atoms
  Type LocalCoordInfo
    Integer               :: nvects   ! Maximum of 3
    Integer, Dimension(3) :: atomnos  ! Atom Nos. which form the basis
    Real(kind=RDbl), Dimension(3) :: xcoeff ! Coeffs to get (1,0,0) from atom
                                            ! vectors
    Real(kind=RDbl), Dimension(3) :: ycoeff ! Coeffs to get (0,1,0) from atom
                                            ! vectors
    Real(kind=RDbl), Dimension(3) :: zcoeff ! Coeffs to get (0,0,1) from atom
                                            ! vectors
  End Type LocalCoordInfo

Contains
  !----------------------------------------------------------------------------
  ! Generates the local coordinate information by getting coefficients
  ! which when multiplied by certain atom vectors give the x, y and z axes.
  ! It stores the vectors and their respective coefficients for
  ! the three different coordinate axes.  These coefficients and vectors
  ! are later used for the purpose of calculating generalized eulerian angles
  ! from cartesian coordinates.  
  ! Requires:  lcoords -- local coordinate data structure to initialize
  !            natoms -- number of atoms in molecule
  !            avecs -- atomic position vectors for each atom
  !            com -- center of mass
  !            geometry -- returned string specifying geometry of molecule
  !            atomnos -- optional atom numbers to use for construction
  !----------------------------------------------------------------------------
  Subroutine localcoords_init(lcoords,natoms,avecs,com,geometry,atomnos)
    Type(LocalCoordInfo), Intent(InOut)           :: lcoords
    Integer, Intent(In)                           :: natoms
    Type(VecType), Dimension(:), Intent(In)       :: avecs
    Type(VecType), Intent(In)                     :: com
    Character(*), Intent(Out)                     :: geometry
    Integer, Dimension(3), Intent(In), Optional   :: atomnos

    Integer                         :: i, natom2
    Logical                         :: found
    Type(VecType)                   :: vec1, vec2, vec3
    Type(MatrixType)                :: vecmatrix, invvecmatrix
    Type(VecType), Dimension(3)     :: e, coeff
    Real(kind=Rdbl), Dimension(3,3) :: vecarray

    !** Set Defaults
    lcoords%nvects = 0
    lcoords%atomnos = 0
    lcoords%xcoeff = 0.0_RDbl
    lcoords%ycoeff = 0.0_RDbl
    lcoords%zcoeff = 0.0_RDbl

    !** We need at least 3 atoms to generate a local coordinate system
    If (natoms < 3) Then
      If (natoms == 1) Then
        geometry = "Spherical"
      Else
        geometry = "Linear"
      End If
      Return
    End If
    
    !** Find the atoms that will form the basis for the local coordinate
    !** system.  We need to find at least 3 non-collinear atoms
    ! Pick an atom which is not at the origin
    Do i = 1,natoms
      vec1 = avecs(i) - com
      If (mag(vec1) > zeroTolerance) Then
        lcoords%atomnos(1) = i
        Exit
      End If
    End Do

    !** Find another atom whose vectors are not collinear
    found = .False.
    Do i = 1,natoms
      vec2 = avecs(i) - com
      natom2 = i
!LC      Write(*,'(2i4,3f8.3,5x,3f8.3)') molecule%bodyaxes%atomnos(1),i,vec1,vec2
      If (.Not. vector_iscollinear(vec1, vec2)) Then
        found = .True.
        Exit
      End If
    End Do
    If (found) Then
      lcoords%atomnos(2) = natom2
    Else
      geometry = "Linear"
      Return
    End If

    !** Vector 3 is just the cross product of vec1 and vec2
    vec3 = vector_crossprod(vec1, vec2)

    !** Now we need to get the coefficients c(i)1, c(i)2, c(i)3 such that
    !** c(i)1*vec1 + c(i)2*vec2 + c(i)3*vec3 = ei where i=1, 3
    e(1) = (/1.0_Rdbl, 0.0_RDbl, 0.0_Rdbl/)
    e(2) = (/0.0_Rdbl, 1.0_RDbl, 0.0_Rdbl/)
    e(3) = (/0.0_Rdbl, 0.0_RDbl, 1.0_Rdbl/)
    vecarray(1:3,1) = vec1
    vecarray(1:3,2) = vec2
    vecarray(1:3,3) = vec3

    !** Construct the matrix from the 2-D array
    vecmatrix = vecarray

    !** Solve for the coefficients
    invvecmatrix = matrix_getinv(vecmatrix)
    Do i = 1,3
      coeff(i) = invvecmatrix*e(i)
    End Do

    lcoords%xcoeff = coeff(1)
    lcoords%ycoeff = coeff(2)
    lcoords%zcoeff = coeff(3)

  End Subroutine localcoords_init

  !----------------------------------------------------------------------------
  ! Generates the vectors in the non-molecule frame that form the body axes.
  ! These can be used to calculate the eulerian angles.  Note that the order
  ! of passed atomic vectors must be the same as it was during _init.
  ! Requires:  lcoords -- local coordinate data structure to initialize
  !            avecs -- atomic position vectors for each atom in molecule
  !            com -- center of mass
  !            xvec,yvec,zvec -- vectors forming body axes
  !----------------------------------------------------------------------------
  Subroutine localcoords_getaxes(lcoords,avecs,com,xvec,yvec,zvec)
    Type(LocalCoordInfo), Intent(In)              :: lcoords
    Type(VecType), Dimension(:), Intent(In)       :: avecs
    Type(VecType), Intent(In)                     :: com
    Type(VecType), Intent(Out)                    :: xvec,yvec,zvec

    Integer                      :: i,atom1,atom2
    Type(VecType), Dimension(3)  :: atvec

    !** Initialization
    xvec = 0.0_RDbl
    yvec = 0.0_RDbl
    zvec = 0.0_RDbl

    !** Get the construction materials
    atom1  = lcoords%atomnos(1)
    atom2  = lcoords%atomnos(2)
    atvec(1) = avecs(atom1) - com
    atvec(2) = avecs(atom2) - com
    atvec(3) = vector_crossprod(atvec(1),atvec(2))

    !** Generate in the fixed coordinate system what the position of 
    !** the body axes would be.  The coordinates of the body axes
    !** as expressed in the fixed coordinates are given by "xf, yf, zf"
    Do i = 1,3
      xvec  = xvec + atvec(i)*lcoords%xcoeff(i)
      yvec  = yvec + atvec(i)*lcoords%ycoeff(i)
      zvec  = zvec + atvec(i)*lcoords%zcoeff(i)
    End Do

  End Subroutine localcoords_getaxes

  !----------------------------------------------------------------------------
  ! Write information about the local coordinates to the specified unit
  ! Requires:  lcoords -- local coordinate data structure
  !            unit -- unit number to dump into
  !            indent -- indentation from left margin
  !----------------------------------------------------------------------------
  Subroutine localcoords_display(lcoords,unit,indent)
    Type(LocalCoordInfo), Intent(InOut) :: lcoords
    Integer, Intent(In)                 :: unit,indent

    Integer                             :: i
    Character(len=indent)               :: blank
    Character(len=strLen)               :: string

    blank = Repeat(' ',indent)
    
    Write(unit,'(2a,3i4)') blank,'body axes atom nos  : ', &
        (lcoords%atomnos(i),i=1,lcoords%nvects)

  End Subroutine localcoords_display

  !----------------------------------------------------------------------------
  ! Clean the local coordinates
  ! Requires:  lcoords -- local coordinate data structure
  !----------------------------------------------------------------------------
  Subroutine localcoords_clean(lcoords)
    Type(LocalCoordInfo), Intent(InOut) :: lcoords
    !** nothing to do
  End Subroutine localcoords_clean

End Module localcoords
