!------------------------------------------------------------------------------
! This module provides a recursive hierarchy of data structures that store
! information between one "group" and another "group".  The separate data 
! structure, and therefore the definition of the stored information is 
! provided by another module called "storebase".  By making the structure
! here recursive, these groups can be infinitely divided into sub-groups.  
! For example, application of this storage structure for forcefield 
! information storage yields the following group definition hierarchy:
!           full_system <- species <- molecule <- atom 
! Each group (e.g. molecule) is made up of sub-groups (e.g. atom).  Since 
! interactions are defined between groups, the structure must be defined 
! between two groups.  Note however, that this module has no knowledge of
! species, molecules or atoms.  The data structure is called a "level-pair" 
! (Store_Level_Pair).  Each data structure contains the total interaction 
! information plus a collection of possible lower-level data structures that 
! contain more detailed information.  Only one of the lower-level data 
! structures can be associated for each level pair.  The recursive structure 
! can be terminated by not associating another type of itself.
!
! The form of the associated lower-level structure is dictated by the need 
! to provide various levels of immediate detail.  The simpliest form is one
! containing detail for only all of the basic elements in the first group.
! This form compresses or integrates over all elements of second group to
! level only detail for the first group.  This structure is initialized when
! detail for the lowest level of the hierarchy is required, but when one
! wants to save memory by not storing all the detail.  For example, MD 
! simulations require atomic forces on each atomic, but the size of the
! system may prohibit storage of all specific atom-atom interactions.
!
! The lower-level level-pair structures are derived from the containing 
! structure by operating on each of the level pairs with either a "minus" or 
! an "identity" operator.  For example, the containing structure is 'spc-spc', 
! then the 3 possible more detailed structures would be: 'mol-mol', 'mol-spc'
! and 'spc-mol'.  Note that the first takes the form of a 2D array and the 
! last two are 1D arrays.  Since only one of the 3 structures is allocated,
! the interaction tree can take a large number of different forms.  This
! provides flexibility.  Note that the number of levels could also be 
! adjusted, for example a group (grp) level could be inserted between the
! atom and the molecules levels (this would be difficult).
!
! When dealing with storage in 2D arrays, data is only stored in the upper
! triangle (m x n, n>=m) to avoid duplication. -- (not always true?, note %mm)
!
! for now, the whole structure is assumed to have the same nderivatives
! initialization level.  Also, it is assumed that only one of the branch 
! structures is initialized.  
!
! Important routines:
!   store_init -- initilizes a single level pair structure
!   store_terminate -- terminate a branch by initializing only total quantity
!   store_zero -- recursively zero level pair structures
!   store_idbranch -- returns an integer indicating which branch is allocated
!   store_copy -- smart copy between two pre-initialized level pair structures
!   store_sum -- fill in the total information for a level pair (opt recursive)
!   store_extract -- extract a single group1-group2 interaction
!   store_compress -- turns group-group storage into group storage 
!   store_resize -- resizes a recursive level pair structure
!------------------------------------------------------------------------------

Module store

  Use defaults, Only: RDbl,strLen,xlstrLen,dbgflag,NO_OF_INTRA_POTS
  Use utils, Only: allocErrDisplay,deallocErrDisplay,int2str
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_display
  Use storebase, Only: EnergyPlus, storebase_init, storebase_scalarmult, &
      storebase_nderivs, storebase_forces, storebase_zero, storebase_sum, &
      storebase_display, storebase_clean, storebase_sumarray, &
      storebase_incAllArray, storebase_disp, storebase_incinv, &
      storebase_initcopy, storebase_copy, storebase_chkintra, &
      storebase_sumintra, storebase_iszero, storebase_scalenrg, &
      storebase_inc, storebase_chkequal, storebase_null, &
      storebase_projectout, storebase_subtract, storebase_add

  Implicit None

  Private
  Public :: Store_Level_Pair, store_init, store_nderivs, &
      store_extract, store_copy, store_sum, store_pack2D, store_unpack2Dinc, &
      store_getforces, store_clean, store_terminate, &
      store_display, store_idbranch, store_zero, store_flipcompress, &
      store_unpack2Dfast, store_initcopy, store_resize, &
      store_scalenrgs, store_scaleall, store_disp, store_chkequal, &
      store_projectout, store_basicptr, store_initcopyidx, store_subtract, &
      store_add

  !** This structure holds interaction between group of level1 with 
  !** group of level2 (eg mol-spc).  Designators mm,mi,im hold different
  !** descents into more detail information.  "binfo" holds basic element
  !** information for the FIRST group when the storage structure must 
  !** be truncated at this level and before the basic element level.
  Type Store_Level_Pair
    Integer, Dimension(2)                           :: sizes
    Type(EnergyPlus), Pointer                       :: total !* sum
    Type(EnergyPlus), Dimension(:), Pointer         :: binfo !* basic info
    Type(Store_Level_Pair), Dimension(:,:), Pointer :: mm    !* minus,minus
    Type(Store_Level_Pair), Dimension(:), Pointer   :: mi    !* minus,identity
    Type(Store_Level_Pair), Dimension(:), Pointer   :: im    !* identity,minus
  End Type Store_Level_Pair

  !** NOTE, the "sizes" array is there because I was having problems
  !** using the Size intrinsic on Store_Level_Pair, especially when
  !** try to use this intrinsic from a higher level. see _sum and _compress

  Interface store_terminate
    Module Procedure store_terminateSingle
    Module Procedure store_terminateArray
    Module Procedure store_terminate2DArray
  End Interface

Contains
  !----------------------------------------------------------------------------
  ! Initialize a single storage level
  ! Requires:  storage -- storage level pair structure
  !            sizes -- array sizes (size1,size2)
  !            initlevel -- the storage level to be initialized
  !            nderivs -- number of derivatives to store (0,1,2)
  !            intra -- flags initialization of intra potentials, default False
  ! Either the "minus-minus", "minus,identity", "identity,minus" or the
  ! "binfo" levels will initialized.
  !   initlevel = 11  => initialize "minus,minus" sub-storage
  !   initlevel = 10  => initialize "minus,identity" sub-storage
  !   initlevel = 01  => initialize "identity,minus" sub-storage
  !   initlevel = 00  => initialize "binfo" sub-storage
  !   initlevel = -1  => don't initialize any of the branches
  !----------------------------------------------------------------------------
  Subroutine store_init(storage,sizes,initlevel,nderivs,intra)
    Type(Store_Level_Pair), Intent(Out)  :: storage
    Integer, Dimension(:), Intent(In)    :: sizes
    Integer, Intent(In)                  :: initlevel,nderivs
    Logical, Intent(In), Optional        :: intra

    Integer            :: i,j,error

    storage%sizes = 0

    !** nullify the pointers first
    Call store_null(storage)

    !** initialize the summed total quantity
    Allocate(storage%total,stat=error)
    If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%total')
    If (Present(intra)) Then
      Call storebase_init(storage%total,nderivs,intra)
    Else
      Call storebase_init(storage%total,nderivs)
    End If

    !** initialize the detailed storage
    Select Case (initlevel)
    Case (11) 
      storage%sizes = sizes
      Allocate(storage%mm(sizes(1),sizes(2)),stat=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%mm')        
      Do i = 1,sizes(1)
        Do j = 1,sizes(2)
          Call store_terminate(storage%mm(i,j),nderivs,intra)
        End Do
      End Do

    Case (10) 
      storage%sizes(1) = sizes(1)
      Allocate(storage%mi(sizes(1)),stat=error)    
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%mi')        
      Do i = 1,sizes(1)
        Call store_terminate(storage%mi(i),nderivs,intra)
      End Do

    Case (01) 
      storage%sizes(1) = sizes(1)
      Allocate(storage%im(sizes(1)),stat=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%im')        
      Do i = 1,sizes(1)
        Call store_terminate(storage%im(i),nderivs,intra)
      End Do

    Case (00) 
      storage%sizes(1) = sizes(1)
      Allocate(storage%binfo(sizes(1)),stat=error)    
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%binfo')        
      If (Present(intra)) Then
        Do i = 1,sizes(1)
          Call storebase_init(storage%binfo(i),nderivs,intra)
        End Do
      Else
        Do i = 1,sizes(1)
          Call storebase_init(storage%binfo(i),nderivs)
        End Do
      End If

    Case (-1) 
      !** Do nothing, branchless level pair

    Case Default
      Write(0,'(2a,i6,a,i9)') __FILE__,":",__LINE__, &
          " could not interpret initialization level", initlevel
      Stop
    End Select

    If (error /= 0) Then
      Write(0,'(2a,i6,a,i2,a)') __FILE__,":",__LINE__, &
          ' Could not allocate memory for requested storage (', initlevel,')'
      If (initlevel == 11) Then
        Write(0,'(3x,5a)') 'Requested array size was: (', &
            Trim(int2str(sizes(1))),',',Trim(int2str(sizes(2))),')'
      Else
        Write(0,'(3x,3a)') 'Requested array size was: (', &
            Trim(int2str(sizes(1))),')'
      End If
      Write(0,'(3x,a)') 'Please consider reducing storage level in control file'
      Stop
    End If

  End Subroutine store_init

  !----------------------------------------------------------------------------
  ! Recursively initialize and copy a storage structure
  ! Requires:  old -- forcefield results structure to copy
  !            image -- new forcefield results structure (must be deallocated)
  !            recurse -- flag indicating desired recursive operation
  !            eraseold -- flag prompting cleaning of old structure
  !----------------------------------------------------------------------------
  Recursive Subroutine store_initcopy(image,old,recurse,eraseold)
    Type(Store_Level_Pair), Intent(InOut)  :: image
    Type(Store_Level_Pair), Intent(In)     :: old
    Logical, Intent(In)                    :: recurse
    Logical, Intent(In), Optional          :: eraseold

    Integer            :: i,j,level,nderivs,error
    Logical            :: intra

    !** Nullify the pointers or clean the structure first if requested
    If (Present(eraseold)) Then
      If (eraseold) Then
        Call storebase_clean(image%total)
        Call store_clean(image)
      Else
        Call store_null(image)
      End If
    Else
      Call store_null(image)
    End If

    !** Copy the basic information
    image%sizes = old%sizes

    !** Allocate the total quantity pointer
    Allocate(image%total,stat=error)
    If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%total')

    !** Get the initialization level of the old structure and allocate
    level = store_idbranch(old)
    If (level == -1) Then
      !** Just copy the total quantity
      Call storebase_initcopy(image%total,old%total)
      Return
    End If
    nderivs = storebase_nderivs(old%total)
    intra = Associated(old%total%intranrg)
    Call store_init(image,old%sizes,level,nderivs,intra)

    !** Copy the total quantity
    Call storebase_copy(image%total,old%total)

    If (.Not. recurse) Return

    !** Copy the contained information recursively
    Select Case(level)
    Case (00)
      Do i = 1,image%sizes(1)
        Call storebase_initcopy(image%binfo(i),old%binfo(i))
      End Do

    Case (01)
      Do i = 1,image%sizes(1)
        Call store_initcopy(image%im(i),old%im(i),recurse)
      End Do

    Case (10)
      Do i = 1,image%sizes(1)
        Call store_initcopy(image%mi(i),old%mi(i),recurse)
      End Do

    Case (11)
      Do i = 1,image%sizes(1)
        Do j = 1,image%sizes(2)
          Call store_initcopy(image%mm(i,j),old%mm(i,j),recurse)
        End Do
      End Do

    End Select    

  End Subroutine store_initcopy

  !----------------------------------------------------------------------------
  ! Initialize and copy just a portion of a storage structure.  It acts only
  ! on a single piece of the next sub-levelpair.  For now, it only works with
  ! the binfo structure.
  ! Requires:  old -- forcefield results structure to copy
  !            image -- new forcefield results structure (must be deallocated)
  !            idx -- integers giving start and end indices to copy
  !            resum -- if True, resum the level pair, otherwise copy total
  ! Is this the best way to write this routine?  
  !----------------------------------------------------------------------------
  Subroutine store_initcopyidx(image,old,idx,resum)
    Type(Store_Level_Pair), Intent(InOut)  :: image
    Type(Store_Level_Pair), Intent(In)     :: old
    Integer, Dimension(2), Intent(In)      :: idx
    Logical, Intent(In)                    :: resum

    Integer            :: i,j,level,nderivs,error,nindices
    Logical            :: intra

    !** Nullify the new structure
    Call store_null(image)

    !** Copy the basic information
    image%sizes = old%sizes

    !** Allocate the total quantity pointer
    Allocate(image%total,stat=error)
    If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%total')

    !** Get the initialization level of the old structure and allocate
    level = store_idbranch(old)
    nderivs = storebase_nderivs(old%total)
    intra = Associated(old%total%intranrg)

    Select Case (level)
    Case(-1)   
      !** Copy the total quantity
      Call storebase_initcopy(image%total,old%total)
      Return

    Case(0)
      nindices = idx(2) - idx(1) + 1
      Call store_init(image,(/nindices/),level,nderivs,intra)
      j = 0
      Do i = idx(1),idx(2)
        j = j + 1
        Call storebase_copy(image%binfo(j),old%binfo(i))  
      End Do

    Case(10)
      nindices = idx(2) - idx(1) + 1
      Call store_init(image,(/nindices/),level,nderivs,intra)
      j = 0
      Do i = idx(1),idx(2)
        j = j + 1
        Call store_initcopy(image%mi(j),old%mi(i),.True.)
      End Do

    Case(01)
      nindices = idx(2) - idx(1) + 1
      Call store_init(image,(/nindices/),level,nderivs,intra)
      j = 0
      Do i = idx(1),idx(2)
        j = j + 1
        Call store_initcopy(image%im(j),old%im(i),.True.)
      End Do

    Case(11)
      Write(0,'(2a,i6,a,i3)') __FILE__,":",__LINE__, &
          ' store_initcopyidx not working for: ',level
      Stop        
    End Select
    
    !** Either resum the structure or just copy the total
    If (resum) Then
      Call store_sum(image,.True.)
    Else 
      Call storebase_copy(image%total,old%total)
    End If

  End Subroutine store_initcopyidx

  !----------------------------------------------------------------------------
  ! Recursively zero a level pair structure's contents
  ! Requires:  storage -- storage level pair structure
  !            recurse -- flag indicating desired recursive operation
  !----------------------------------------------------------------------------
  Recursive Subroutine store_zero(storage,recurse)
    Type(Store_Level_Pair), Intent(InOut)  :: storage
    Logical, Intent(In)                    :: recurse

    Integer            :: i,j

    Call storebase_zero(storage%total)

    Select Case(store_idbranch(storage))
    Case (00)
      If (recurse) Then
        Do i = 1,storage%sizes(1)
          Call storebase_zero(storage%binfo(i))
        End Do
      End If

    Case (01)
      If (recurse) Then
        Do i = 1,storage%sizes(1)
          Call store_zero(storage%im(i),recurse)
        End Do
      End If

    Case (10)
      If (recurse) Then
        Do i = 1,storage%sizes(1)
          Call store_zero(storage%mi(i),recurse)
        End Do
      End If

    Case (11)
      If (recurse) Then
        Do i = 1,storage%sizes(1)
          Do j = 1,storage%sizes(2)
            Call store_zero(storage%mm(i,j),recurse)
          End Do
        End Do
      End If

    End Select    

  End Subroutine store_zero

  !----------------------------------------------------------------------------
  ! Nullify the pointers of a storage level
  ! Requires:  storage -- storage level pair structure
  !----------------------------------------------------------------------------
  Subroutine store_null(storage)
    Type(Store_Level_Pair), Intent(InOut)  :: storage

    !** nullify the pointers in the structure
    Nullify(storage%total)    
    Nullify(storage%binfo)    
    Nullify(storage%mm)    
    Nullify(storage%mi)
    Nullify(storage%im)

  End Subroutine store_null

  !----------------------------------------------------------------------------
  ! Terminate a single storage level structure by only initializing the
  ! 'total' component of the structure.
  ! Requires:  storage -- array of storage level pair structures
  !            nderivs -- number of derivatives to store (0,1,2), default 0
  !            intra -- flags initialization of intra potentials, default False
  !----------------------------------------------------------------------------
  Subroutine store_terminateSingle(storage,nderivs,intra)
    Type(Store_Level_Pair), Intent(InOut)  :: storage
    Integer, Intent(In)                    :: nderivs
    Logical, Intent(In), Optional          :: intra

    Integer       :: error

    storage%sizes = 0

    !** terminate the tree by nullifying the pointers
    Call store_null(storage)

    !** Allocate the total quantity pointer
    Allocate(storage%total,stat=error)
    If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%total')

    !** initialize the summed total quantity
    If (Present(intra)) Then
      Call storebase_init(storage%total,nderivs,intra)
    Else
      Call storebase_init(storage%total,nderivs)
    End If

  End Subroutine store_terminateSingle

  !----------------------------------------------------------------------------
  ! Terminate an array of storage level structures but only initializing the
  ! 'total' component of each structure.
  ! Requires:  storage -- array of storage level pair structures
  !            nderivs -- number of derivatives to store (0,1,2)
  !            intra -- flags initialization of intra potentials, default False
  !----------------------------------------------------------------------------
  Subroutine store_terminateArray(storage,nderivs,intra)
    Type(Store_Level_Pair), Dimension(:), Intent(InOut)  :: storage
    Integer, Intent(In)                                  :: nderivs
    Logical, Intent(In), Optional                        :: intra

    Integer            :: i,error

    Do i = 1,Size(storage)
      storage(i)%sizes = 0

      !** Allocate the total quantity pointer
      Allocate(storage(i)%total,stat=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%total')

      !** initialize the summed total quantity
      If (Present(intra)) Then
        Call storebase_init(storage(i)%total,nderivs,intra)
      Else
        Call storebase_init(storage(i)%total,nderivs)
      End If

      !** terminate the tree by nullifying the pointers
      Call store_null(storage(i))
    End Do

  End Subroutine store_terminateArray

  !----------------------------------------------------------------------------
  ! Terminate an array of storage level structures but only initializing the
  ! 'total' component of each structure.
  ! Requires:  storage -- array of storage level pair structures
  !            nderivs -- number of derivatives to store (0,1,2)
  !            intra -- flags initialization of intra potentials, default False
  !----------------------------------------------------------------------------
  Subroutine store_terminate2DArray(storage,nderivs,intra)
    Type(Store_Level_Pair), Dimension(:,:), Intent(InOut)  :: storage
    Integer, Intent(In)                                    :: nderivs
    Logical, Intent(In), Optional                          :: intra

    Integer            :: i,j,error

    Do i = 1,Size(storage,1)
      Do j = 1,Size(storage,2)
        storage(i,j)%sizes = 0

        !** Allocate the total quantity pointer
        Allocate(storage(i,j)%total,stat=error)
        If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'%total')

        !** initialize the summed total quantity
        If (Present(intra)) Then
          Call storebase_init(storage(i,j)%total,nderivs,intra)
        Else
          Call storebase_init(storage(i,j)%total,nderivs)
        End If

        !** terminate the tree by nullifying the pointers
        Call store_null(storage(i,j))
      End Do
    End Do

  End Subroutine store_terminate2DArray

  !----------------------------------------------------------------------------
  ! Copy one level pair structure into another.  The new structure must first
  ! be initialized.  If the initializations of the two structure do not match,
  ! then the routine will attempt to restructure the information to perform
  ! the copy.  If information is missing, then it will return False.  This 
  ! routine assumes that the "original" structure is fully sum-updated, i.e.
  ! that it has alread been run through store_sum.
  ! Requires:  dest -- destination level pair structure (must be initialized)
  !            orig -- origin level pair structure
  !            recurse -- flag indicating desired recursive operation
  !----------------------------------------------------------------------------
  Logical Recursive Function store_copy(dest,orig) Result(res)
    Type(Store_Level_Pair), Intent(InOut)  :: dest
    Type(Store_Level_Pair), Intent(In)     :: orig

    Integer                      :: i,j,dinfo,oinfo,n
    Integer, Dimension(2)        :: dlevels,olevels

    res = .False.

    !** Copy the total quantity
    Call storebase_copy(dest%total,orig%total)

    !** Get structure of the two level-pairs
    dinfo = store_idbranch(dest,dlevels)
    oinfo = store_idbranch(orig,olevels)

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,*) 'levels (orig,dest): ',oinfo,dinfo

    !** Do the actual copying if possible
    If (dinfo == oinfo) Then
      !** structures are identical at this level, try to copy
      If (dinfo <= 0) Then
        Call storebase_copy(dest%total,orig%total)
        If (dinfo == 0) Then
          Do i = 1,Size(orig%binfo)
            Call storebase_copy(dest%binfo(i),orig%binfo(i))
          End Do
        End If
        res = .True.

      Else
        Select Case(oinfo)
        Case (01)
          Do i = 1,orig%sizes(1)
            res = store_copy(dest%im(i),orig%im(i))
          End Do
    
        Case (10)
          Do i = 1,orig%sizes(1)
            res = store_copy(dest%mi(i),orig%mi(i))
          End Do
    
        Case (11)
          Do i = 1,orig%sizes(1)
            Do j = 1,orig%sizes(2)
              res = store_copy(dest%mm(i,j),orig%mm(i,j))
            End Do
          End Do
    
        End Select
      End If

    Else 
      !** Compress into a single base structure, always possible
      If (dinfo == -1) Then
        Call storebase_copy(dest%total,orig%total)
        res = .True.
        Return
      End If

      !** Insufficient detail in the origin structure, return False
      If ((oinfo == 0).Or.(oinfo == -1)) Return

      !** Compress into an array of base structures
      If (dinfo == 0) Then
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) 'not ready yet'
        Stop
        n = Size(dest%binfo)  !change later?
!LC        Call store_compress(orig,n,dest%binfo)
        res = .True.
        Return
      End If

      !** may be possible, but haven't written code for it yet
      Write(0,'(2a,i6,a,i3)') __FILE__,":",__LINE__, &
          ' Need to write some code for this functionality'
      Stop        
    End If

  End Function store_copy

  !----------------------------------------------------------------------------
  ! Get a pointer to the EnergyPlus structure that stores a given interaction.
  ! This routine is controlled by a pair of paths.  
  ! Requires:  lp -- storage level pair structure
  !            path1 -- level 1 path
  !            path2 -- level 2 path
  !            ptr -- returned pointer
  !----------------------------------------------------------------------------
  Logical Recursive Function store_basicptr(lp,path1,path2,ptr) Result(res)
    Type(Store_Level_Pair), Intent(In)   :: lp
    Integer, Dimension(:), Intent(In)    :: path1,path2
    Type(EnergyPlus), Pointer            :: ptr

    Integer                     :: level

    !** Set default return value
    res = .False.

    !** Check to see if we've reached the desired level pair
    If ((path1(1) == 0).And.(path2(1) == 0)) Then
      ptr => lp%total
      res = .True.
      Return
    End If

    !** Otherwise, descend further into the level-pair structures
    level = store_idbranch(lp)
    Select Case(level)
    Case (11)   !** descend into both levels
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'NOT TESTED'
      Stop

      res = store_basicptr(lp%mm(path1(1),path2(1)),path1(2:),path2(2:),ptr)
      If (.Not. res) Return

    Case (10)   !** descend into only FIRST level
      If (path1(1) /= 0) Then
        res = store_basicptr(lp%mi(path1(1)),path1(2:),path2,ptr)
        If (.Not. res) Return
      End If

    Case (01)   !** descend into only SECOND level
      If (path2(1) /= 0) Then
        res = store_basicptr(lp%im(path2(1)),path1,path2(2:),ptr)
        If (.Not. res) Return
      End If

    Case (0)    !** descend into condensed storage
      If (path1(1) /= 0) Then
        res = .True.
        ptr => lp%binfo(path1(1))
      End If

    End Select

  End Function store_basicptr

  !----------------------------------------------------------------------------
  ! Sum over the detail provided in the detail of the level-pair.  Can operate
  ! recursively if desired (recommended for safety).
  ! Requires:  storage -- storage level pair structure
  !            recurse -- flag indicating desired recursive operation
  !            total -- optional summed result (pre-initialized)
  !----------------------------------------------------------------------------
  Recursive Subroutine store_sum(storage,recurse,total)
    Type(Store_Level_Pair), Intent(InOut)   :: storage
    Logical, Intent(In)                     :: recurse
    Type(EnergyPlus), Intent(Out), Optional :: total

    Integer            :: i,j

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,*) 'entering store_sum ',store_idbranch(storage)
!    Write(*,*) ' ',store_idbranch(storage)

    Select Case(store_idbranch(storage))
    Case (00)
      Call storebase_sumarray(storage%binfo,storage%total)

    Case (01)
      Call storebase_zero(storage%total)
      Do i = 1,storage%sizes(1)
        If (recurse) Call store_sum(storage%im(i),recurse)
        Call storebase_inc(storage%total,storage%im(i)%total)
      End Do

    Case (10)
      Call storebase_zero(storage%total)
      Do i = 1,storage%sizes(1)
        If (recurse) Call store_sum(storage%mi(i),recurse)
        Call storebase_inc(storage%total,storage%mi(i)%total)
      End Do

    Case (11)
      Call storebase_zero(storage%total)
      Do i = 1,storage%sizes(1)
        Do j = 1,storage%sizes(2)
          If (recurse) Call store_sum(storage%mm(i,j),recurse)
          Call storebase_inc(storage%total,storage%mm(i,j)%total)
        End Do
      End Do

    End Select

    If (Present(total)) total = storage%total
    
  End Subroutine store_sum

  !----------------------------------------------------------------------------
  ! Subtract the second base structure from the first.  Assumes that anything
  ! initialized in the 2nd structure is initialized in the 1st.
  ! Requires:  storage1 -- 1st level pair structure
  !            storage2 -- 2nd level pair structure
  !            recurse -- flag indicating desired recursive operation
  !----------------------------------------------------------------------------
  Recursive Subroutine store_subtract(storage1,storage2,recurse)
    Type(Store_Level_Pair), Intent(InOut)   :: storage1
    Type(Store_Level_Pair), Intent(In)      :: storage2
    Logical, Intent(In)                     :: recurse

    Integer            :: i,j

    !** Handle the total quantity
    Call storebase_subtract(storage1%total,storage2%total)

    Select Case(store_idbranch(storage1))
    Case (00)
      Do i = 1,storage2%sizes(1)
        Call storebase_subtract(storage1%binfo(i),storage2%binfo(i))
      End Do

    Case (01)
      Do i = 1,storage2%sizes(1)
        If (recurse) Call store_subtract(storage1%im(i),storage2%im(i),recurse)
      End Do

    Case (10)
      Do i = 1,storage2%sizes(1)
        If (recurse) Call store_subtract(storage1%mi(i),storage2%mi(i),recurse)
      End Do

    Case (11)
      Do i = 1,storage2%sizes(1)
        Do j = 1,storage2%sizes(2)
          If (recurse) Call store_subtract(storage1%mm(i,j), &
              storage2%mm(i,j),recurse)
        End Do
      End Do

    End Select

  End Subroutine store_subtract

  !----------------------------------------------------------------------------
  ! Add the second base structure to the first.  Assumes that the
  ! initialization structures are identical
  ! Requires:  storage1 -- 1st level pair structure
  !            storage2 -- 2nd level pair structure
  !            recurse -- flag indicating desired recursive operation
  !----------------------------------------------------------------------------
  Recursive Subroutine store_add(storage1,storage2,recurse)
    Type(Store_Level_Pair), Intent(InOut)   :: storage1
    Type(Store_Level_Pair), Intent(In)      :: storage2
    Logical, Intent(In)                     :: recurse

    Integer            :: i,j

    !** Handle the total quantity
    Call storebase_add(storage1%total,storage2%total)

    Select Case(store_idbranch(storage1))
    Case (00)
      Do i = 1,storage1%sizes(1)
        Call storebase_add(storage1%binfo(i),storage2%binfo(i))
      End Do

    Case (01)
      Do i = 1,storage1%sizes(1)
        If (recurse) Call store_add(storage1%im(i),storage2%im(i),recurse)
      End Do

    Case (10)
      Do i = 1,storage1%sizes(1)
        If (recurse) Call store_add(storage1%mi(i),storage2%mi(i),recurse)
      End Do

    Case (11)
      Do i = 1,storage1%sizes(1)
        Do j = 1,storage1%sizes(2)
          If (recurse) Call store_add(storage1%mm(i,j), &
              storage2%mm(i,j),recurse)
        End Do
      End Do

    End Select

  End Subroutine store_add

  !----------------------------------------------------------------------------
  ! The recursive portion of storetop_extract to operate on the level pair
  ! data structure.  Returns the info structure when it can, otherwise
  ! returns False.  If the first index pair of the path array is not zero 
  ! then it indicates the next descent level(s) for that side of the level 
  ! pair.  If it is zero, this means the desired level pair structure has 
  ! been reached for that side of the level pair.
  ! Requires:  lp -- Level pair data structure corresponding to path1,path2
  !            path1 -- level 1 path
  !            path2 -- level 2 path
  !            info -- returned energy plus possible derivatives structure
  ! Path specification:
  !   the 'paths' take the form of a 2D integer array with the second
  !   dimension have a fixed size of two.  Example:
  !     path1 =  1 1 0   
  !              3 8 0
  !   will instruct the routine to descend into indices 1->3 of the first
  !   level.  At each of those sub-levels, it will descend into indices
  !   1->8.  The zeros at the end of the path indicate that the descending
  !   is to be terminated.  Before each descent, the head (first column) of
  !   the path is cropped-off.
  ! Needed Improvements
  ! 1) could add option for resumming to add security
  !----------------------------------------------------------------------------
  Logical Recursive Function store_extract(lp,path1,path2,info) Result(res)
    Type(Store_Level_Pair), Intent(In)     :: lp
    Integer, Dimension(:,:), Intent(In)    :: path1,path2
    Type(EnergyPlus), Intent(InOut)        :: info

    Integer                     :: i,j,level

    !** Set default return value
    res = .False.

#if DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'before next descent, level = ',store_idbranch(lp), &
        '   sizes = ',lp%sizes
    Write(*,*) 'path1: ',path1(:,1),'   path2: ',path2(:,1)
    Write(*,*) '       ',path1(:,2),'          ',path2(:,2)
#endif

    !** Check to see if we've reached the desired level pair
    If ((path1(1,1) == 0).And.(path2(1,1) == 0)) Then
      Call storebase_inc(info,lp%total)
      res = .True.
      Return
      Write(*,*) 'successful - lp%nrg,info%nrg: ',lp%total%nrg,info%nrg
    End If

    !** Otherwise, descend further into the level-pair structures
    level = store_idbranch(lp)
    Select Case(level)
    Case (11)   !** descend into both levels
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'NOT TESTED'
      Stop

      If (path1(1,1) == 0) Then
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Stop
      Else If (path2(1,1) == 0) Then
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Stop
      Else
        Do i = path1(1,1),path1(1,2)
          Do j = path2(1,1),path2(1,2)
            res = store_extract(lp%mm(i,j),path1(2:,:),path2(2:,:),info)
            If (.Not. res) Return
          End Do
        End Do
      End If

    Case (10)   !** descend into only FIRST level
      If (path1(1,1) /= 0) Then
        Do i = path1(1,1),path1(1,2)
          res = store_extract(lp%mi(i),path1(2:,:),path2,info)
          If (.Not. res) Return
        End Do
      End If

    Case (01)   !** descend into only SECOND level
      If (path2(1,1) /= 0) Then
        Do j = path2(1,1),path2(1,2)
          res = store_extract(lp%im(j),path1,path2(2:,:),info)
          If (.Not. res) Return
        End Do
      End If

    End Select

  End Function store_extract

  !-------------------------------------------------------------------------     
  ! Projects out a specified component of a force vector.  Uses integer 
  ! array "paths".  If the first index pair of the path array is not zero 
  ! then it indicates the next descent level(s) for that side of the level 
  ! pair.  If it is zero, this means the desired level pair structure has 
  ! been reached for that side of the level pair.  See the description of
  ! _extract for double array path instructions.
  ! Requires:  lp -- Level pair data structure
  !            path1 -- level 1 path
  !            path2 -- level 2 path
  !            vec -- unit vector along which to remove force component
  !-------------------------------------------------------------------------     
  Logical Recursive Function store_projectout(lp,path1,path2,vec) Result(res)
    Type(Store_Level_Pair), Intent(InOut)    :: lp
    Integer, Dimension(:,:), Intent(In)      :: path1,path2
    Type(VecType), Intent(In)                :: vec

    Integer                     :: i,j,level

    !** Set default return value
    res = .False.

    !** Check to see if we've reached the desired level pair
    If ((path1(1,1) == 0).And.(path2(1,1) == 0)) Then
      res = .True.
      Call storebase_projectout(lp%total,vec)
      Return
    End If

    !** Otherwise, descend further into the level-pair structures
    level = store_idbranch(lp)
    Select Case(level)
    Case (11)   !** descend into both levels
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'NOT TESTED'
      Stop

      If (path1(1,1) == 0) Then
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Stop
      Else If (path2(1,1) == 0) Then
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Stop
      Else
        Do i = path1(1,1),path1(1,2)
          Do j = path2(1,1),path2(1,2)
            res = store_projectout(lp%mm(i,j),path1(2:,:),path2(2:,:),vec)
            If (.Not. res) Return
          End Do
        End Do
      End If

    Case (10)   !** descend into only FIRST level
      If (path1(1,1) /= 0) Then
        Do i = path1(1,1),path1(1,2)
          res = store_projectout(lp%mi(i),path1(2:,:),path2,vec)
          If (.Not. res) Return
        End Do
      End If

    Case (01)   !** descend into only SECOND level
      If (path2(1,1) /= 0) Then
        Do j = path2(1,1),path2(1,2)
          res = store_projectout(lp%im(j),path1,path2(2:,:),vec)
          If (.Not. res) Return
        End Do
      End If

    Case (00)   !** descend into condensed atom info structure
      Do i = path1(1,1),path1(1,2)
        res = .True.
        Call storebase_projectout(lp%binfo(i),vec)
      End Do

    End Select

  End Function store_projectout

  !----------------------------------------------------------------------------
  ! Extract a 1D array of just the forces on each of the base elements in
  ! the first group of the level pair.  The order of the base elements in
  ! the 1D array corresponds to that provided by _compress.
  ! Requires:  lp -- Level pair data structure corresponding to path1,path2
  !            n -- number of base elements to get gradients for (MUST MATCH!)
  !            flip -- if True, instead extract for base elements in 2nd level 
  !            opt_nmoles -- I dont know where else this routine is used, 
  !                          so this will just assume that lp Contains 
  !                          forces of opt_nmoles molecules
  ! Possible Improvements:
  ! 1) bypassing _compress and putting grads directly in would be faster
  !----------------------------------------------------------------------------
  Function store_getforces(lp,n,flip,opt_nmoles)
    Type(Store_Level_Pair), Intent(In) :: lp
    Integer, Intent(In)                :: n
    Logical, Intent(In)                :: flip
    Type(VecType), Dimension(n)        :: store_getforces
    Integer, Intent(In), Optional      :: opt_nmoles

    Integer                            :: i,nbases
    Type(EnergyPlus), Dimension(n)     :: results

    !** Allocate results array so we can use _compress to get gradients
    Do i = 1,n
      Call storebase_init(results(i),1,.True.)
    End Do
    
#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    nbases = 0
    Call store_count(lp,nbases)
    Write(*,*) 'real nbases, passed final size ',nbases,n
    Write(*,*) flip
#endif
    
    !** Compress information from below into 1D results array
    If (flip) Then
      Call store_flipcompress(lp,1,n,results)
    Else
      If (Present(opt_nmoles)) Then
        Call store_compress(lp,1,n,results,opt_nmoles)
      Else
        Call store_compress(lp,1,n,results)
      Endif
    End If

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) n

    Do i = 1,n
      Write(*,*) i,'  ',Trim(storebase_disp(results(i)))
    End Do
    !    stop
#endif

    !** Extract gradients
    store_getforces = storebase_forces(results)    

    !** Deallocate results array 
    Call storebase_clean(results)

  End Function store_getforces

  !----------------------------------------------------------------------------
  ! Contracts information from a Level Pair data structure into just 
  ! interactions for the fundamental unit of the FIRST level of the pair
  ! (atom-detail).  This means that the details of the atoms in the first
  ! level interacting with individual atoms in the second level is lost, all
  ! interactions are summed onto those atoms of the first level.  This is
  ! useful for getting net forces on the atoms of the first level.  The output
  ! is a 1D array of base storage structures.  The resulting array
  ! has elements ordered in the most fundamental sequence.  For example:
  ! {(mol1,atm1),(mol1,atm2),...(mol1,atmN),...(mol2,atm1),...(molN,atmN)}.
  ! Definitions: "compressing" means to sum multiple group contributions. 
  !              This is implemented by repassing same results section
  !              "compiling" means record results for new group
  !              This is implemented by repassing only part of incoming results
  ! NOTE that this routine only INCREMENTS the passed 1D storage array.
  ! Also, I originally wrote this just to pass chunks of the results array, 
  ! but I found it easier to debug if I passed the wrong array plus indices.
  ! Requires:  lp -- Level pair data structure corresponding to path1,path2
  !            start -- starting index to increment in results array
  !            end -- ending index to increment in results array
  !            results -- the resultant 1D array of base storage types
  !            nmoles -- says how many molecules storage should be taken from
  !                   lp and put it into results. If absent the whole storage
  !                   from lp will be compressed and put into results
  !----------------------------------------------------------------------------
  Recursive Subroutine store_compress(lp,start,last,results,opt_nmoles)
    Type(Store_Level_Pair), Intent(In)            :: lp
    Integer, Intent(In)                           :: start,last
    Type(EnergyPlus), Dimension(:), Intent(InOut) :: results
    Integer, Intent(In), Optional                 :: opt_nmoles 

    Integer            :: i,j,s1,s2
    Integer            :: hi,lo,length,nbases
    Type(VecType)      :: sumvec

    !** Number of incoming array entries to increment at this level
    nbases = last - start + 1

    !** Look at initialized branch and operate accordingly
    Select Case(store_idbranch(lp))
    Case (00)   !** just copy (storage truncated)
      Call storebase_incAllArray(results(start:last),lp%binfo)

    Case (01)   !** compress through 2nd level
      Do i = 1,Size(lp%im)
!        Write(*,*) 'im-associated, compressing along 2nd level: ',start,last
        Call store_compress(lp%im(i),start,last,results)
      End Do

    Case (10)   !** compile along 1st level
      !** it should not be size, it should be only the number of molecules
      If (Present(opt_nmoles)) Then
        s1 = opt_nmoles
      Else
        s1 = Size(lp%mi)
      End If
      length = nbases/s1

      Do i = 1,s1
        lo = start + length*(i-1)
        hi = lo + length - 1
!        Write(*,*) 'mi-associated, compressing along 1st level: ',lo,hi
        Call store_compress(lp%mi(i),lo,hi,results)
      End Do

    Case (11)   !** compile 1st level, compress 2nd level
      s1 = Size(lp%mm,1)
      s2 = Size(lp%mm,2)
      length = nbases/s1

      Do i = 1,s1
        lo = start + length*(i-1)
        hi = lo + length - 1
        Do j = 1,s2
!          Write(*,'(a,2i3,a,2i3)') 'mm-associated, ',i,j, &
!         ' compiling along 1st level: ',lo,hi
          Call store_compress(lp%mm(i,j),lo,hi,results)
        End Do
      End Do

#ifdef DEBUG
      sumvec = VecType(0.0_RDbl)
      Do i = 1,8
        Write(*,*) Trim(storebase_disp(results(i)))
        sumvec = sumvec + results(i)%force
      End Do
      Write(*,*) Trim(vector_display(sumvec,'e14.4'))
#endif

    Case Default   !** no branches, increment the results
!      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Call storebase_inc(results(start),lp%total)
      Return

      !** Can use this for debugging IF the structure is allocated symmetrically
      If (start == last) Then
        Call storebase_inc(results(start),lp%total)
!        Write(*,*) 'base ',start,Trim(storebase_disp(lp%total))
      Else
        Write(0,'(2a,i6,a,2i8)') __FILE__,":",__LINE__, &
            " Error: range passed to _compress was wrong ",start,last
        Stop
      End If

    End Select

  End Subroutine store_compress

  !----------------------------------------------------------------------------
  ! Contracts information from a Level Pair data structure into just 
  ! interactions for the fundamental unit of the SECOND level of the pair
  ! (atom-detail).  Thus, it is similar, but opposite to the function of 
  ! _compress.  The output is a 1D array of base storage structures. with
  ! elements ordered in the most fundamental sequence.  For example:
  ! {(mol1,atm1),(mol1,atm2),...(mol1,atmN),...(mol2,atm1),...(molN,atmN)}.
  ! Definitions: "compressing" means to sum multiple group contributions. 
  !              This is implemented by repassing same results section
  !              "compiling" means record results for new group
  !              This is implemented by repassing only part of incoming results
  ! NOTE that this routine only INCREMENTS the passed 1D storage array.
  ! Requires:  lp -- Level pair data structure corresponding to path1,path2
  !            start -- starting index to increment in results array
  !            end -- ending index to increment in results array
  !            results -- the resultant 1D array of base storage types
  !----------------------------------------------------------------------------
  Recursive Subroutine store_flipcompress(lp,start,last,results)
    Type(Store_Level_Pair), Intent(In)            :: lp
    Integer, Intent(In)                           :: start,last
    Type(EnergyPlus), Dimension(:), Intent(InOut) :: results

    Integer            :: i,j,s1,s2
    Integer            :: hi,lo,length,nbases
    Type(VecType)      :: sumvec

    !** number of incoming array entries to increment at this level
    nbases = last - start + 1

    !** Look at initialized branch and operate accordingly
    Select Case(store_idbranch(lp))
    Case (00)   !** just copy (storage truncated)
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' Not tested, please look at code'
      Stop
      !** I think that all the binfo array items need to be compressed
      !** onto a single item in results, think about this

    Case (10)   !** compress through 2nd level
      Do i = 1,Size(lp%mi)
!        Write(*,*) 'mi-associated, compressing along 2nd level: ',start,last
        Call store_flipcompress(lp%mi(i),start,last,results)
      End Do

    Case (01)   !** compile along 1st level
      s1 = Size(lp%im)
      length = nbases/s1

      Do i = 1,s1
        lo = start + length*(i-1)
        hi = lo + length - 1
!        Write(*,*) 'im-associated, compressing along 1st level: ',lo,hi
        Call store_flipcompress(lp%im(i),lo,hi,results)
      End Do

    Case (11)   !** compress 1st level, compile 2nd level
      s1 = Size(lp%mm,1)
      s2 = Size(lp%mm,2)
      length = nbases/s2

      Do j = 1,s2
        lo = start + length*(j-1)
        hi = lo + length - 1
        Do i = 1,s1
!          Write(*,'(a,2i3,a,2i3)') 'mm-associated, ',i,j, &
!              ' compressing along 1st level: ',lo,hi
          Call store_flipcompress(lp%mm(i,j),lo,hi,results)
        End Do
      End Do

#ifdef DEBUG
      sumvec = VecType(0.0_RDbl)
      Do i = 1,8
        Write(*,*) Trim(storebase_disp(results(i)))
        sumvec = sumvec + results(i)%grad
      End Do
      Write(*,*) Trim(vector_display(sumvec,'e14.4'))
#endif

    Case Default   !** no branches

      Call storebase_inc(results(start),lp%total)
      Return

      !** Can use this for debugging IF the structure is allocated symmetrically
      If (start == last) Then
        Call storebase_inc(results(start),lp%total)
!        Write(*,*) 'base ',start,Trim(storebase_disp(lp%total))

      Else
        Write(0,'(2a,i6,a,2i8)') __FILE__,":",__LINE__, &
            " Error: range passed to _flipcompress was wrong ",start,last
        Stop
      End If

    End Select

  End Subroutine store_flipcompress

  !-------------------------------------------------------------------------  
  ! Resizes level-pair structures given new size(s)
  ! Requires:  lp -- Level pair data structure
  !            sizes -- desired new size(s)
  !-------------------------------------------------------------------------  
  Subroutine store_resize(lp,sizes)
    Type(Store_Level_Pair), Intent(InOut)  :: lp
    Integer, Dimension(2), Intent(In)      :: sizes
    
    Integer                         :: i,j,initlevel,nderivs
    Integer                         :: max1,max2,error
    Logical                         :: intra,decrease
    Type(Store_Level_Pair), Pointer :: lptmp
    Character(len=xlstrLen)         :: junk

#ifdef DEBUG
    junk = ''
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'BEFORE resizing'
    Call store_display(lp,junk,.False.,2,6)
#endif

    !** Make sure we're increasing the size
    decrease = .False.
    If (sizes(1) < lp%sizes(1)) decrease = .True.
    If (Size(sizes) > 1) Then
      If (sizes(2) < lp%sizes(2)) decrease = .True.
    End If
    If (decrease) Then
      Write(0,'(2a,i6,a,i9)') __FILE__,":",__LINE__, &
          ' Only prepared to increase storage size'
      Stop
    End If

    !** Make a copy of the existing structure
    Allocate(lptmp, STAT=error)
    If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'lptmp')        
    Call store_initcopy(lptmp,lp,.True.)

    !** Erase the old structure
    Call store_clean(lp)   

    !** Create the top of a new structure
    initlevel = store_idbranch(lptmp)
    nderivs = storebase_nderivs(lptmp%total)
    intra = storebase_chkintra(lptmp%total)
    Call store_init(lp,sizes,initlevel,nderivs,intra)

    !** Copy the total quantity, no reason for it to change
    Call storebase_copy(lp%total,lptmp%total)

    !** Fill in the lower parts using the temporary structure
    max1 = lptmp%sizes(1)
    Select Case(initlevel)
    Case (00)   
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'NOT TESTED'
      Stop

      !** Copy those branches that are preserved
      Do i = 1,max1
        Call storebase_initcopy(lp%binfo(i),lptmp%binfo(i))
      End Do

      !** Make sure the rest are zeroed
      Do i = max1,lp%sizes(1)
        Call storebase_zero(lp%binfo(i))
      End Do

    Case (10)
      !** Copy those branches that are preserved
      Do i = 1,max1
        Call store_initcopy(lp%mi(i),lptmp%mi(i),.True.)
      End Do

      !** Crude way to preserve the shape, but makes sure new entries are zero
      Do i = max1+1,lp%sizes(1)
        Call store_initcopy(lp%mi(i),lptmp%mi(max1),.True.)
        Call store_zero(lp%mi(i),.True.)
      End Do

    Case (01)   
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'NOT TESTED'
      Stop

      !** Copy those branches that are preserved
      Do i = 1,max1
        Call store_initcopy(lp%im(i),lptmp%im(i),.True.)
      End Do

      !** Crude way to preserve the shape, but makes sure new entries are zero
      Do i = max1+1,lp%sizes(1)    
        Call store_initcopy(lp%im(i),lptmp%im(max1),.True.)
        Call store_zero(lp%im(i),.True.)
      End Do

    Case (11)   
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'NOT TESTED'
      Stop

      max2 = lptmp%sizes(2)
      !** Copy those branches that are preserved
      Do i = 1,max1
        Do j = 1,max2
          Call store_initcopy(lp%mm(i,j),lptmp%mm(i,j),.True.)
        End Do
      End Do

      !** Crude way to preserve the shape, but makes sure new entries are zero
      Do i = max1+1,lp%sizes(1)    
        Do j = max2+1,lp%sizes(2)
          Call store_initcopy(lp%mm(i,j),lptmp%mm(max1,max2),.True.)
          Call store_zero(lp%mm(i,j),.True.)
        End Do
      End Do

    Case (-1)
      !nothing to do here

    Case Default
      Write(0,'(2a,i6,a,i9)') __FILE__,":",__LINE__, &
          ' Should not be here, could not understand lp init level'
      Stop

    End Select

    !** Eliminate the temporary storage
    Call store_clean(lptmp)
    Deallocate(lptmp, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'lptmp')        

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'AFTER resizing'
    Call store_display(lp,junk,.False.,2,6)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Stop
#endif

  End Subroutine store_resize

  !----------------------------------------------------------------------------
  ! Returns the highest associated derivative pointer in the storage structure
  ! Requires:  lp -- level pair storage structure
  !----------------------------------------------------------------------------
  Integer Function store_nderivs(lp)
    Type(Store_Level_Pair), Intent(In)          :: lp

    store_nderivs = storebase_nderivs(lp%total)  

  End Function store_nderivs

  !----------------------------------------------------------------------------
  ! Counts the number of initialized "bases" of the storage tree. NOTE, this
  ! routine only incrments the 'count' variables, it must be zeroed outside.
  ! Requires:  lp -- level pair storage structure
  !            count -- returned number of initialized bases
  !----------------------------------------------------------------------------
  Recursive Subroutine store_count(lp,count)
    Type(Store_Level_Pair), Intent(In)   :: lp
    Integer, Intent(InOut)               :: count

    Integer          :: i,j

    If (Associated(lp%im)) Then   
      Do i = 1,Size(lp%im)
        Call store_count(lp%im(i),count)
      End Do
    Else If (Associated(lp%mi)) Then   
      Do i = 1,Size(lp%mi)
        Call store_count(lp%mi(i),count)
      End Do
    Else If (Associated(lp%mm)) Then  
      Do i = 1,Size(lp%mm,1)
        Do j = 1,Size(lp%mm,2)
          Call store_count(lp%mm(i,j),count)
        End Do
      End Do
    Else If (Associated(lp%binfo)) Then
      count = count + Size(lp%binfo)
    Else
      count = count + 1
    End If    

  End Subroutine store_count

  !----------------------------------------------------------------------------
  ! Returns integer notation for the associated storage level pair
  ! Requires:  lp -- level pair storage structure
  !            levels -- optional two integer output
  !----------------------------------------------------------------------------
  Integer Function store_idbranch(lp,levels)
    Type(Store_Level_Pair), Intent(In)   :: lp
    Integer, Dimension(2), Optional      :: levels

    If (Associated(lp%im)) Then   
      store_idbranch = 1
      If (Present(levels)) Then
        levels(1) = 0
        levels(2) = 1
      End If
    Else If (Associated(lp%mi)) Then   
      store_idbranch = 10
      If (Present(levels)) Then
        levels(1) = 1
        levels(2) = 0
      End If
    Else If (Associated(lp%mm)) Then  
      store_idbranch = 11
      If (Present(levels)) Then
        levels(1) = 1
        levels(2) = 1
      End If
    Else If (Associated(lp%binfo)) Then
      store_idbranch = 0
      If (Present(levels)) Then
        levels(1) = 0
        levels(2) = 0
      End If
    Else
      store_idbranch = -1
      If (Present(levels)) Then
        levels(1) = -1
        levels(2) = -1     
      End If
    End If

  End Function store_idbranch

  !----------------------------------------------------------------------------
  ! This routine packs the basic EnergyPlus structures from an mm-associated
  ! level pair storage structure into a 1D array based on a mask.  Note that 
  ! the inside loop runs over the row index, this is consistent with the 
  ! intrinsic Pack function.
  ! Requires:  lp -- Level pair data structure, mm must be associated
  !            mask -- a 2D array, turning packing on/off for each index pair
  !            array -- output array of basic EnergyPlus structures
  !----------------------------------------------------------------------------
  Subroutine store_pack2D(lp,mask,array)
    Type(Store_Level_Pair), Intent(In)          :: lp
    Logical, Dimension(:,:), Intent(In)         :: mask                
    Type(EnergyPlus), Dimension(:), Intent(Out) :: array

    Integer         :: i,j,idx

    idx = 0
    Do j = 1,Size(lp%mm,2)
      Do i = 1,Size(lp%mm,1)
        If (mask(i,j)) Then
          idx = idx + 1
          Call storebase_initcopy(array(idx),lp%mm(i,j)%total)
!          array(idx) = lp%mm(i,j)%total  !this actually copies pointer locations
        End If
      End Do
    End Do

  End Subroutine store_pack2D

  !----------------------------------------------------------------------------
  ! This routine unpacks an array of basic EnergyPlus structures into two
  ! mm-associated level pair storage structures using a mask.  The 1D array
  ! of EnergyPlus structures is assumed to be ordered like this:
  ! {(mol1,atm1),(mol1,atm2),...(mol1,atmN),...(mol2,atm1),...(molN,atmN)}.
  ! Additionally, the quantities in this array are assumed to apply to the
  ! 'lp1' structure.  The mirror-image of the quantities (equal and opposite
  ! reaction) is added to the 'lp2' structure.  This means the same energies
  ! are added, but that the negative of the force/gradients are added.
  !
  ! The unpacking only INCREMENTS the storage structure.  Note that the inside 
  ! loop runs over the row index, this is consistent with the intrinsic Pack 
  ! function.
  ! Requires:  array -- input array of basic EnergyPlus structures
  !            mask -- a 2D array, turning unpacking on/off for each index pair
  !            lp1 -- Level pair data structure, mm must be associated
  !            lp2 -- Level pair data structure, mm must be associated
  !----------------------------------------------------------------------------
  Subroutine store_unpack2Dinc(array,mask,lp1,lp2)
    Type(EnergyPlus), Dimension(:), Intent(In)  :: array
    Logical, Dimension(:,:), Intent(In)         :: mask                
    Type(Store_Level_Pair), Intent(InOut)       :: lp1,lp2

    Integer         :: i,j,idx

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Do i = 1,16
      Write(*,*) i,Trim(storebase_disp(array(i)))
    End Do
    Write(*,*) 
#endif

    idx = 0
    Do j = 1,Size(mask,2)
      Do i = 1,Size(mask,1)
        If (mask(i,j)) Then
          idx = idx + 1
!          Write(*,*) i,j,Trim(storebase_disp(lp1%mm(i,j)%total))
!          Write(*,*) i,j,Trim(storebase_disp(array(idx)))
          Call storebase_inc(lp1%mm(i,j)%total,array(idx))

          Call storebase_incinv(lp2%mm(i,j)%total,array(idx))
!          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!          Write(*,*) i,j,Trim(storebase_disp(array(idx)))
!          Stop
        End If
      End Do
    End Do

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) idx

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Do i = 1,idx
      Write(*,*) i,Trim(storebase_disp(array(i)))
    End Do
    Write(*,*) 
#endif

  End Subroutine store_unpack2Dinc

  !----------------------------------------------------------------------------
  ! attempt to speed by doing increments internally
  !----------------------------------------------------------------------------
  Subroutine store_unpack2Dfast(array,mask,lp1,lp2,nderivs)
    Type(EnergyPlus), Dimension(:), Intent(In)  :: array
    Logical, Dimension(:,:), Intent(In)         :: mask                
    Type(Store_Level_Pair), Intent(InOut)       :: lp1,lp2
    Integer, Intent(In)                         :: nderivs

    Integer         :: i,j,idx

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Do i = 1,16
      Write(*,*) i,Trim(storebase_disp(array(i)))
    End Do
    Write(*,*) 
#endif

    idx = 0
    Do j = 1,Size(mask,2)
      Do i = 1,Size(mask,1)
        If (mask(i,j)) Then
          idx = idx + 1
!          Write(*,*) i,j,Trim(storebase_disp(lp1%mm(i,j)%total))
!          Write(*,*) i,j,Trim(storebase_disp(array(idx)))
!          Call storebase_inc(lp1%mm(i,j)%total,array(idx))

          lp1%mm(i,j)%total%nrg = lp1%mm(i,j)%total%nrg + array(idx)%nrg
          lp1%mm(i,j)%total%force = lp1%mm(i,j)%total%force + array(idx)%force

          lp2%mm(i,j)%total%nrg = lp2%mm(i,j)%total%nrg + array(idx)%nrg
          lp2%mm(i,j)%total%force = lp2%mm(i,j)%total%force - array(idx)%force

          !** increment the total also
          lp1%total%nrg =  lp1%total%nrg + array(idx)%nrg
          lp1%total%force = lp1%total%force + array(idx)%force
          lp2%total%nrg =  lp2%total%nrg + array(idx)%nrg
          lp2%total%force = lp2%total%force - array(idx)%force

!          Call storebase_incinv(lp2%mm(i,j)%total,array(idx))
!          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!          Write(*,*) i,j,Trim(storebase_disp(array(idx)))
!          Stop
        End If
      End Do
    End Do

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) idx

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Do i = 1,idx
      Write(*,*) i,Trim(storebase_disp(array(i)))
    End Do
    Write(*,*) 
#endif

  End Subroutine store_unpack2Dfast

  !-------------------------------------------------------------------------     
  ! Scales all energies in the level pair
  ! Requires:  lp -- Level pair data structure
  !            factor -- factor to scale energies by
  !-------------------------------------------------------------------------     
  Recursive Subroutine store_scalenrgs(lp,factor)
    Type(Store_Level_Pair), Intent(InOut)    :: lp
    Real(kind = RDbl), Intent(In)            :: factor

    Integer       :: i,j,level

    !** scale the summed quantity (or the base, if no branches associated)
    Call storebase_scalenrg(lp%total,factor)

    level = store_idbranch(lp)
    Select Case(level)
    Case(01)  !** Descend into the im-Associated level-pair
      Do i = 1,lp%sizes(1)
        Call store_scalenrgs(lp%im(i),factor)
      End Do

    Case(10)  !** Descend into the mi-Associated level-pair
      Do i = 1,lp%sizes(1)
        Call store_scalenrgs(lp%mi(i),factor)
      End Do

    Case(11)  !** Descend into the mm-Associated level-pair
      Do i = 1,lp%sizes(1)
        Do j = 1,lp%sizes(2)
          Call store_scalenrgs(lp%mm(i,j),factor)
        End Do
      End Do

    Case(0)   !** Descend into the simple array
      Do i = 1,Size(lp%binfo)
        Call storebase_scalenrg(lp%binfo(i),factor)
      End Do

    End Select

  End Subroutine store_scalenrgs

  !-------------------------------------------------------------------------     
  ! Scales all interactions in the level pair
  ! Requires:  lp -- Level pair data structure
  !            factor -- factor to scale energies by
  !-------------------------------------------------------------------------     
  Recursive Subroutine store_scaleall(lp,factor)
    Type(Store_Level_Pair), Intent(InOut)    :: lp
    Real(kind = RDbl), Intent(In)            :: factor

    Integer       :: i,j,level

    !** scale the summed quantity (or the base, if no branches associated)
    Call storebase_scalarmult(lp%total,factor)

    level = store_idbranch(lp)
    Select Case(level)
    Case(01)  !** Descend into the im-Associated level-pair
      Do i = 1,lp%sizes(1)
        Call store_scaleall(lp%im(i),factor)
      End Do

    Case(10)  !** Descend into the mi-Associated level-pair
      Do i = 1,lp%sizes(1)
        Call store_scaleall(lp%mi(i),factor)
      End Do

    Case(11)  !** Descend into the mm-Associated level-pair
      Do i = 1,lp%sizes(1)
        Do j = 1,lp%sizes(2)
          Call store_scaleall(lp%mm(i,j),factor)
        End Do
      End Do

    Case(0)   !** Descend into the simple array
      Do i = 1,Size(lp%binfo)
        Call storebase_scalarmult(lp%binfo(i),factor)
      End Do

    End Select

  End Subroutine store_scaleall

  !----------------------------------------------------------------------------
  ! Compares all the initialized entries in two structures.  If a unit
  ! number and indent specific are included, the routine will dump output. 
  ! Otherwise, only the True/False answer will be returned.
  ! Requires:  lp1 -- 1st Level Pair structure  
  !            lp2 -- 2nd Level Pair structure  
  !            tol -- tolerance for calling equal
  !            indent -- no. of spaces from the left margin
  !            unit -- optional display unit number
  !----------------------------------------------------------------------------
  Logical Recursive Function store_chkequal(lp1,lp2,outstring,tol, &
      indent,unit) Result(res)
    Type(Store_Level_Pair), Intent(In)      :: lp1,lp2
    Character(len=xlstrLen), Intent(InOut)  :: outstring
    Real(kind=RDbl), Intent(In)             :: tol
    Integer, Intent(In), Optional           :: indent,unit

    Integer                      :: i,j,level,nd,nstrings
    Logical                      :: equal,dump
    Character(len=strLen)        :: string1,string2
    Character(len=xlstrLen)      :: string
    Character(len=xlstrlen), Dimension(100)  :: stringset

    dump = .False.
    If ((Present(indent)).And.(Present(unit))) dump = .True.

    !** Set default return value
    res = .True.

    !** Get level pair characteristics
    nd = store_nderivs(lp1)
    level = store_idbranch(lp1)

    !** Check for obvious errors
    If (level /= store_idbranch(lp2)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__,' : ',__LINE__, &
          ' Passed level pairs do not have same structure shape'
      Stop            
    End If
    If (lp1%sizes(1) /= lp2%sizes(1)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__,' : ',__LINE__, &
          ' Passed level pairs do not have same structure size'
      Stop            
    End If

    !** Do comparison based on form of level pair
    Select Case(level)
    Case(01)  !** Descend into the im-Associated level-pair
      Do i = 1,lp1%sizes(1)
        string = outstring
        string1 = int2str(nd)
        Write(string,'(4a)') Trim(string),'(',Trim(string1),')'
        string1 = int2str(i)
        Write(string,'(a,1x,2a)') Trim(string),'im:',Trim(string1)
        If (dump) Then
          equal = store_chkequal(lp1%im(i),lp2%im(i),string,tol,indent,unit)
        Else
          equal = store_chkequal(lp1%im(i),lp2%im(i),string,tol)
        End If
        If (.Not. equal) res = .False.
      End Do

    Case(10)  !** Descend into the mi-Associated level-pair
      Do i = 1,lp1%sizes(1)
        string = outstring
        string1 = int2str(nd)
        Write(string,'(4a)') Trim(string),'(',Trim(string1),')'
        string1 = int2str(i)
        Write(string,'(a,1x,2a)') Trim(string),'mi:',Trim(string1)
        If (dump) Then
          equal = store_chkequal(lp1%mi(i),lp2%mi(i),string,tol,indent,unit)
        Else
          equal = store_chkequal(lp1%mi(i),lp2%mi(i),string,tol)
        End If
        If (.Not. equal) res = .False.
      End Do

    Case(11)  !** Descend into the mm-Associated level-pair
      Do i = 1,lp1%sizes(1)
        Do j = 1,lp1%sizes(2)
          string = outstring
          string1 = int2str(nd)
          Write(string,'(4a)') Trim(string),'(',Trim(string1),')'
          string1 = int2str(i)
          string2 = int2str(j)
          Write(string,'(a,1x,4a)') Trim(string),'mm:',Trim(string1), &
              ',',Trim(string2)
          If (dump) Then
            equal = store_chkequal(lp1%mm(i,j),lp2%mm(i,j),string, &
                tol,indent,unit)
          Else
            equal = store_chkequal(lp1%mm(i,j),lp2%mm(i,j),string,tol)
          End If
          If (.Not. equal) res = .False.
        End Do
      End Do

    Case(0)   !** Descend into the simple array
      If (Size(stringset) < Size(lp1%binfo)) Then
        Write(0,'(2a,i6,a,i9)') __FILE__,":",__LINE__, &
            ' Temporary string storage not large enough' 
        Stop
      End If

      nstrings = 0
      Do i = 1,Size(lp1%binfo)
        string = outstring
        string1 = int2str(i)
        If (.Not. storebase_chkequal(lp1%binfo(i),lp2%binfo(i),tol,string)) Then
          res = .False.
          If (dump) Then 
            nstrings = nstrings + 1
            string1 = int2str(i)
            Write(stringset(nstrings),'(a3,1x,a)') Trim(string1),Trim(string)
          End If
        End If
      End Do

      !** Dump the strings if any were recorded
      If (nstrings > 0) Then
        string1 = int2str(Size(lp1%binfo))
        string2 = int2str(indent)
        Write(string,'(3a)') '(',Trim(string2),'x,3a)'
        Write(unit,string) Trim(outstring), &
          ' Condensed atom info, size ',Trim(string1)
        Write(string,'(3a)') '(',Trim(string2),'x,a)'
        Do i = 1,nstrings
          Write(unit,string) Trim(stringset(i))
        End Do
      End If

    Case(-1)  !** Tree is terminated here
      string = outstring
      If (.Not. storebase_chkequal(lp1%total,lp2%total,tol,string)) Then
        res = .False.
        If (dump) Write(unit,'(20a)') (' ',j=1,indent),Trim(string)
      End If

    End Select

  End Function store_chkequal

  !----------------------------------------------------------------------------
  ! Display the level pair structure recursively.  This is just a front-end
  ! for store_display that doesn't require a string input.
  ! Requires:  lp -- Level pair data structure
  !            skip -- True => will skip display if nrg,grad(1) < tolerance
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine store_disp(lp,skip,indent,unitno)
    Type(Store_Level_Pair), Intent(In)     :: lp
    Logical, Intent(In)                    :: skip
    Integer, Intent(In)                    :: indent
    Integer, Optional, Intent(In)          :: unitno

    Character(len=xlstrLen)           :: outstring

    outstring = 'debug: '
    Call store_display(lp,outstring,skip,indent,unitno)

  End Subroutine store_disp

  !----------------------------------------------------------------------------
  ! Display the level pair structure recursively.  Operates such that it 
  ! doesn't display anything until it reaches a terminated level pair.  Also
  ! has an option to not display results if they are close to zero.
  ! Requires:  lp -- Level pair data structure
  !            string -- beginning of output line to display (will be added to)
  !            skip -- True => will skip display if nrg,grad(1) < tolerance
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  ! Gives output that looks like this (for example):
  !   s:1(1) mi:30(1) mi:4 nd=1 NRG=0.000 GRAD=  3.2592 -6.2798 -5.0163 INTRA
  !  
  !  "s:1(1)" means that it's species 1 with 1 derivative in its 'total' 
  !  "mi:30" means the next level is a descent into index 30 of the 1st level
  !  "nd=1" gives the number of derivatives at this base of the tree
  !  "NRG=" and "GRAD=" give values for the energy and gradient
  !  "INTRA" and "HESSIAN" flag presence of intramolecular nrg and hessian 
  !----------------------------------------------------------------------------
  Recursive Subroutine store_display(lp,string,skip,indent,unitno)
    Type(Store_Level_Pair), Intent(In)     :: lp
    Character(len=xlstrLen), Intent(InOut) :: string
    Logical, Intent(In)                    :: skip
    Integer, Intent(In)                    :: indent
    Integer, Optional, Intent(In)          :: unitno

    Integer                           :: i,j,unit,nd,level
    Character(len=indent)             :: blank
    Character(len=strLen)             :: string1,string2
    Character(len=xlstrLen)           :: outstring,shortstr

    blank = Repeat(' ',indent)
    
    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    nd = storebase_nderivs(lp%total)

#ifdef DEBUG
    If (store_idbranch(lp) /= -1) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) store_idbranch(lp),lp%sizes
    End If

    !** Use this to display the total quantity and check associations of it
    shortstr = storebase_disp(lp%total)
    Write(unit,'(3a)') blank,'Top: ',Trim(shortstr)
#endif

    level = store_idbranch(lp)
    Select Case(level)
    Case(01)  !** Descend into the im-Associated level-pair
      Do i = 1,lp%sizes(1)
        string1 = int2str(nd)
        Write(outstring,'(4a)') Trim(string),'(',Trim(string1),')'
        string1 = int2str(i)
        Write(outstring,'(a,1x,2a)') Trim(outstring),'im:',Trim(string1)
        Call store_display(lp%im(i),outstring,skip,indent,unit)
      End Do

    Case(10)  !** Descend into the mi-Associated level-pair
      Do i = 1,lp%sizes(1)
        string1 = int2str(nd)
        Write(outstring,'(4a)') Trim(string),'(',Trim(string1),')'
        string1 = int2str(i)
        Write(outstring,'(a,1x,2a)') Trim(outstring),'mi:',Trim(string1)
        Call store_display(lp%mi(i),outstring,skip,indent,unit)
      End Do

    Case(11)  !** Descend into the mm-Associated level-pair
      Do i = 1,lp%sizes(1)
        Do j = 1,lp%sizes(2)
          string1 = int2str(nd)
          Write(outstring,'(4a)') Trim(string),'(',Trim(string1),')'
          string1 = int2str(i)
          string2 = int2str(j)
          Write(outstring,'(a,1x,4a)') Trim(outstring),'mm:',Trim(string1), &
              ',',Trim(string2)
          Call store_display(lp%mm(i,j),outstring,skip,indent,unit)
        End Do
      End Do

    Case(0)   !** Descend into the simple array
      string1 = int2str(Size(lp%binfo))
      Write(unit,'(2a,1x,2a)') blank,Trim(string), &
          'Condensed atom info, size ',Trim(string1)
      Do i = 1,Size(lp%binfo)
        If (skip) Then
          If (storebase_iszero(lp%binfo(i),1.0e-10_RDbl)) Cycle
        End If
        string1 = int2str(i)
        outstring = storebase_disp(lp%binfo(i))
        Write(unit,'(a,3x,a4,1x,a)') blank,Trim(string1),Trim(outstring)
      End Do

    Case(-1)  !** If the tree is terminated here, then display the storage
      !** if desired, check base structure for zeros before displaying
      If (skip) Then
        If (storebase_iszero(lp%total,1.0e-10_RDbl)) Return
      End if

      !** actually display the storage results
      outstring = storebase_disp(lp%total)
      Write(unit,'(2a,1x,a)') blank,Trim(string),Trim(outstring)

    End Select
    
  End Subroutine store_display

  !----------------------------------------------------------------------------
  ! Clean a storage structure recursively
  ! Requires:  storage -- storage level pair structure
  !----------------------------------------------------------------------------
  Recursive Subroutine store_clean(storage)
    Type(Store_Level_Pair), Intent(InOut)  :: storage

    Integer            :: i,j,error

    Call storebase_clean(storage%total)

    If (Associated(storage%binfo)) Then
      Do i = 1,Size(storage%binfo)
        Call storebase_clean(storage%binfo(i))
      End Do
      Deallocate(storage%binfo, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'binfo')
    End If 

    If (Associated(storage%mm)) Then
      Do i = 1,Size(storage%mm,1)
        Do j = 1,Size(storage%mm,2)
          Call store_clean(storage%mm(i,j))
        End Do
      End Do
      Deallocate(storage%mm, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'mm')
    End If 

    If (Associated(storage%mi)) Then
      Do i = 1,Size(storage%mi)
        Call store_clean(storage%mi(i))
      End Do
      Deallocate(storage%mi, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'mi')
    End If 

    If (Associated(storage%im)) Then
      Do i = 1,Size(storage%im)
        Call store_clean(storage%im(i))
      End Do
      Deallocate(storage%im, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'im')
    End If 

  End Subroutine store_clean

End Module store
