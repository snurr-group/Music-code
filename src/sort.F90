!--------------------------------------------------------
! This module contains routines to sort an array of type
! "SortArray_Type"
!--------------------------------------------------------
Module sort
  Use defaults, Only: RDbl

  Implicit None
  Save

  Private
  Public :: SortArray_Type, sort_init, sort_heapsort

  Integer, Parameter  :: MAX_KEYS = 5
  
  Type SortArray_Type
    Integer         :: index
    Integer         :: nfields
    Real(kind=RDbl), Dimension(MAX_KEYS) :: keylist
  End Type SortArray_Type

  Interface Assignment(=)
    Module Procedure sort_s1_eq_s2
  End Interface

  Interface sort_init
    Module Procedure sort_init1
    Module Procedure sort_init2
    Module Procedure sort_init3
    Module Procedure sort_init4
    Module Procedure sort_init5
  End Interface

Contains

  !------------------------------------------------------------
  ! Generates an array of "SortArray_Type" given an array of
  ! keys that needs to be sorted
  !------------------------------------------------------------
  Subroutine sort_init1(keyarray1, sortarray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: keyarray1
    Type(SortArray_Type), Dimension(:), Intent(out) :: sortarray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      sortarray(i)%index = i
      sortarray(i)%keylist(1) = keyarray1(i)
      sortarray(i)%nfields = 1
    End Do
  End Subroutine sort_init1

  !------------------------------------------------------------
  ! Generates an array of "SortArray_Type" given an array of
  ! keys that needs to be sorted
  !------------------------------------------------------------
  Subroutine sort_init2(keyarray1, keyarray2, sortarray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: keyarray1, keyarray2
    Type(SortArray_Type), Dimension(:), Intent(out) :: sortarray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      sortarray(i)%index = i
      sortarray(i)%keylist(1) = keyarray1(i)
      sortarray(i)%keylist(2) = keyarray2(i)
      sortarray(i)%nfields = 2
    End Do
  End Subroutine sort_init2

  !------------------------------------------------------------
  ! Generates an array of "SortArray_Type" given an array of
  ! keys that needs to be sorted
  !------------------------------------------------------------
  Subroutine sort_init3(keyarray1, keyarray2, keyarray3, sortarray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: &
        keyarray1, keyarray2, keyarray3
    Type(SortArray_Type), Dimension(:), Intent(out) :: sortarray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      sortarray(i)%index = i
      sortarray(i)%keylist(1) = keyarray1(i)
      sortarray(i)%keylist(2) = keyarray2(i)
      sortarray(i)%keylist(3) = keyarray3(i)
      sortarray(i)%nfields = 3
    End Do
  End Subroutine sort_init3

  !------------------------------------------------------------
  ! Generates an array of "SortArray_Type" given an array of
  ! keys that needs to be sorted
  !------------------------------------------------------------
  Subroutine sort_init4(keyarray1, keyarray2, keyarray3, keyarray4, &
      sortarray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: &
        keyarray1, keyarray2, keyarray3, keyarray4
    Type(SortArray_Type), Dimension(:), Intent(out) :: sortarray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      sortarray(i)%index = i
      sortarray(i)%keylist(1) = keyarray1(i)
      sortarray(i)%keylist(2) = keyarray2(i)
      sortarray(i)%keylist(3) = keyarray3(i)
      sortarray(i)%keylist(4) = keyarray4(i)
      sortarray(i)%nfields = 4
    End Do
  End Subroutine sort_init4

  !------------------------------------------------------------
  ! Generates an array of "SortArray_Type" given an array of
  ! keys that needs to be sorted
  !------------------------------------------------------------
  Subroutine sort_init5(keyarray1, keyarray2, keyarray3, keyarray4, &
      keyarray5, sortarray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: &
        keyarray1, keyarray2, keyarray3, keyarray4, keyarray5
    Type(SortArray_Type), Dimension(:), Intent(out) :: sortarray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      sortarray(i)%index = i
      sortarray(i)%keylist(1) = keyarray1(i)
      sortarray(i)%keylist(2) = keyarray2(i)
      sortarray(i)%keylist(3) = keyarray3(i)
      sortarray(i)%keylist(4) = keyarray4(i)
      sortarray(i)%keylist(5) = keyarray5(i)
      sortarray(i)%nfields = 5
    End Do
  End Subroutine sort_init5


  !----------------------------------------------------------------------
  ! Compares two elements and returns true if "s1" is ordinally
  ! less than the corresponding key fields of "s2".
  !----------------------------------------------------------------------
  Logical Function sort_isLessThan(s1, s2)
    Type(SortArray_Type), Intent(in) :: s1, s2
    
    Integer    :: nfields, i

    nfields = s1%nfields
    Do i=1, nfields
      If ((s1%keylist(i) > s2%keylist(i))) Then
        sort_isLessThan = .False.
        Return
      Else If ((s1%keylist(i) < s2%keylist(i))) Then
        sort_isLessThan = .True.
        Return
      End If
    End Do
    ! All the fields are equal so we return false
    sort_isLessThan = .False.
  End Function sort_isLessThan

  !-------------------------------------------------------------
  ! Assign one element of the array "s1" to "s2"
  !-------------------------------------------------------------
  Subroutine sort_s1_eq_s2(s1, s2)
    Type(SortArray_Type), Intent(out) :: s1
    Type(SortArray_Type), Intent(in)  :: s2

    s1%index = s2%index
    s1%nfields = s2%nfields
    s1%keylist(1:s2%nfields) = s2%keylist(1:s2%nfields)
  End Subroutine sort_s1_eq_s2

  !---------------------------------------------------------
  ! Use "heapsort" to sort the array "arr".  
  ! Ref. "Numerical Recipes in Fortran 90".  On completion
  ! "arr" is an array sorted by "arr(i)%key", "arr(i)%index"
  ! is the value of the initial index of the value "arr(i)%key"
  ! before sorting.  The size of the array is "n"
  !---------------------------------------------------------
  Subroutine sort_heapsort(arr, n)
    Type(SortArray_Type), Dimension(:), Intent(inout) :: arr
    Integer, Intent(in)      :: n

    Type(SortArray_Type)    :: temp
    Integer :: i
    
    !** Build the initial heap. We only need to look at the
    !** elements in the first half
    Do i=n/2,1,-1
      Call sift_down(i,n)
    End Do
    
    !** Construct the sorted list
    Do i=n,2,-1
      ! Swap the first element which is the largest element
      ! with the last element of the array and decrease the
      ! array length by one
      temp = arr(1)
      arr(1) = arr(i)
      arr(i) = temp
      
      ! Now reconstruct the heap
      Call sift_down(1,i-1)
    End Do
    
  Contains
    !-------------------------------------------------------
    ! This routine inserts the element at position "l"
    ! at the correct position in the heap of size 1 to "r"
    !-------------------------------------------------------
    Subroutine sift_down(l,r)
      Integer, Intent(IN) :: l,r
      
      Integer    :: j,jold
      Type(SortArray_Type)   :: a
      
      a=arr(l)
      jold=l      ! The index of the empty root
      j=l+l       ! The left child of the empty root
      Do
        If (j > r) Exit
        If (j < r) Then
          If (sort_isLessThan(arr(j), arr(j+1))) j=j+1
        End If
        If (.Not. (sort_isLessThan(a, arr(j)))) Exit
        arr(jold)=arr(j)
        jold=j
        j=j+j
      End Do
      arr(jold) = a
    End Subroutine sift_down
  End Subroutine sort_heapsort

End Module sort

