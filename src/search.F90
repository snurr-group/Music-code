!---------------------------------------------------------------
! This module contains routines to search an array of type
! "SearchArray_Type" by the criteria stored in "keylist"
!---------------------------------------------------------------
Module search
  Use defaults, Only: RDbl

  Implicit None
  Save

  Private
  Public :: SearchArray_Type, MAX_SRCHKEYS, search_init, search_binarysearch

  Integer, Parameter  :: MAX_SRCHKEYS = 5
  
  Type SearchArray_Type
    Integer         :: nfields
    Real(kind=RDbl), Dimension(MAX_SRCHKEYS) :: elemlist
  End Type SearchArray_Type

  Interface Assignment(=)
    Module Procedure search_s1_eq_s2
  End Interface

  Interface search_init
    Module Procedure search_init1
    Module Procedure search_init2
    Module Procedure search_init3
    Module Procedure search_init4
    Module Procedure search_init5
  End Interface

Contains
  !------------------------------------------------------------
  ! Generates an array of "SearchArray_Type" given an array of
  ! keys that needs to be searched
  !------------------------------------------------------------
  Subroutine search_init1(keyarray1, searcharray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: keyarray1
    Type(SearchArray_Type), Dimension(:), Intent(out) :: searcharray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      searcharray(i)%elemlist(1) = keyarray1(i)
      searcharray(i)%nfields = 1
    End Do
  End Subroutine search_init1

  !------------------------------------------------------------
  ! Generates an array of "SearchArray_Type" given an array of
  ! keys that needs to be searched
  !------------------------------------------------------------
  Subroutine search_init2(keyarray1, keyarray2, searcharray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: keyarray1, keyarray2
    Type(SearchArray_Type), Dimension(:), Intent(out) :: searcharray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      searcharray(i)%elemlist(1) = keyarray1(i)
      searcharray(i)%elemlist(2) = keyarray2(i)
      searcharray(i)%nfields = 2
    End Do
  End Subroutine search_init2

  !------------------------------------------------------------
  ! Generates an array of "SearchArray_Type" given an array of
  ! keys that needs to be searched
  !------------------------------------------------------------
  Subroutine search_init3(keyarray1, keyarray2, keyarray3, searcharray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: &
        keyarray1, keyarray2, keyarray3
    Type(SearchArray_Type), Dimension(:), Intent(out) :: searcharray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      searcharray(i)%elemlist(1) = keyarray1(i)
      searcharray(i)%elemlist(2) = keyarray2(i)
      searcharray(i)%elemlist(3) = keyarray3(i)
      searcharray(i)%nfields = 3
    End Do
  End Subroutine search_init3

  !------------------------------------------------------------
  ! Generates an array of "SearchArray_Type" given an array of
  ! keys that needs to be searched
  !------------------------------------------------------------
  Subroutine search_init4(keyarray1, keyarray2, keyarray3, keyarray4, &
      searcharray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: &
        keyarray1, keyarray2, keyarray3, keyarray4
    Type(SearchArray_Type), Dimension(:), Intent(out) :: searcharray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      searcharray(i)%elemlist(1) = keyarray1(i)
      searcharray(i)%elemlist(2) = keyarray2(i)
      searcharray(i)%elemlist(3) = keyarray3(i)
      searcharray(i)%elemlist(4) = keyarray4(i)
      searcharray(i)%nfields = 4
    End Do
  End Subroutine search_init4

  !------------------------------------------------------------
  ! Generates an array of "SearchArray_Type" given an array of
  ! keys that needs to be searched
  !------------------------------------------------------------
  Subroutine search_init5(keyarray1, keyarray2, keyarray3, keyarray4, &
      keyarray5, searcharray, n)
    Real(kind=Rdbl), Dimension(:), Intent(in)       :: &
        keyarray1, keyarray2, keyarray3, keyarray4, keyarray5
    Type(SearchArray_Type), Dimension(:), Intent(out) :: searcharray
    Integer, Intent(in)     :: n

    Integer     :: i
    
    Do i=1, n
      searcharray(i)%elemlist(1) = keyarray1(i)
      searcharray(i)%elemlist(2) = keyarray2(i)
      searcharray(i)%elemlist(3) = keyarray3(i)
      searcharray(i)%elemlist(4) = keyarray4(i)
      searcharray(i)%elemlist(5) = keyarray5(i)
      searcharray(i)%nfields = 5
    End Do
  End Subroutine search_init5


  !----------------------------------------------------------------------
  ! Compares two elements and returns true if each key field of "s1" is
  ! is less than the corresponding key fields of "s2"
  !----------------------------------------------------------------------
  Logical Function search_isLessThan(s1, s2)
    Type(SearchArray_Type), Intent(in) :: s1, s2
    
    Integer    :: nfields, i

    nfields = s1%nfields
    Do i=1, nfields
      If ((s1%elemlist(i) > s2%elemlist(i))) Then
        search_isLessThan = .False.
        Return
      Else If ((s1%elemlist(i) < s2%elemlist(i))) Then
        search_isLessThan = .True.
        Return
      End If
    End Do
    ! All the fields are equal so we return false
    search_isLessThan = .False.
  End Function search_isLessThan

  !----------------------------------------------------------------------
  ! Compares two elements and returns true if "s1" comes AFTER "s2" in
  ! an ordinally
  !----------------------------------------------------------------------
  Logical Function search_isMoreThan(s1, s2)
    Type(SearchArray_Type), Intent(in) :: s1, s2
    
    Integer    :: nfields, i

    nfields = s1%nfields
    Do i=1, nfields
      If ((s1%elemlist(i) > s2%elemlist(i))) Then
        search_isMoreThan = .True.
        Return
      Else If ((s1%elemlist(i) < s2%elemlist(i))) Then
        search_isMoreThan = .False.
        Return
      End If
    End Do
    ! All the fields are equal so we return false
    search_isMoreThan = .False.
  End Function search_isMoreThan

  !-------------------------------------------------------------
  ! Assign one element of the array "s1" to "s2"
  !-------------------------------------------------------------
  Subroutine search_s1_eq_s2(s1, s2)
    Type(SearchArray_Type), Intent(out) :: s1
    Type(SearchArray_Type), Intent(in)  :: s2

    s1%nfields = s2%nfields
    s1%elemlist(1:s2%nfields) = s2%elemlist(1:s2%nfields)
  End Subroutine search_s1_eq_s2

  !---------------------------------------------------------
  ! Use "binarysearch" to search the array "arr" for the keys  
  ! stored in "elemlist". The size of the array is "n"
  !---------------------------------------------------------
  Subroutine search_binarysearch(arr, srchkeylist, n, index)
    Type(SearchArray_Type), Dimension(:), Intent(in) :: arr
    Real(kind=RDbl), Dimension(:), Intent(in) :: srchkeylist
    Integer, Intent(in)      :: n
    Integer, Intent(out)     :: index

    Type(SearchArray_Type)   :: srchkey
    Integer :: i, nkeys, low, high, mid

    !** Create an element of SearchArray_Type
    nkeys = Size(srchkeylist,1)
    If (nkeys > MAX_SRCHKEYS) Then
      Write(0,'(1x,2a,i4,a,i4)') __FILE__," : ",__LINE__, &
          " Can't search on more than ", MAX_SRCHKEYS
      Stop
    End If
    Do i=1, nkeys
      srchkey%elemlist(i) = srchkeylist(i)
    End Do
    srchkey%nfields = nkeys

    !** Now do the binary search and get the index where "srchkey" would
    !** appear ordinally
    low = 1
    high = n
    mid = (low + high)/2
    Do
      If (low > high) Exit
      If (search_isLessThan(srchkey, arr(mid))) Then
        high = mid - 1
      Else If (search_isMoreThan(srchkey, arr(mid))) Then
        low = mid + 1
      Else
        high = mid
        low  = mid+1
      End If
      mid = (low + high)/2
    End Do
    index = low
  End Subroutine search_binarysearch
  
  !---------------------------------------------------------------
  ! This routine does a linear search for the keys in arr
  !---------------------------------------------------------------
  Subroutine search_linearsearch(arr, keylist, n)
    Type(SearchArray_Type), Dimension(:), Intent(inout) :: arr
    Real(kind=RDbl), Dimension(:), Intent(in) :: keylist
    Integer, Intent(in)      :: n

    Type(SearchArray_Type)    :: temp
    Integer :: i
  End Subroutine search_linearsearch

End Module search

