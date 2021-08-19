!------------------------------------------------------------------------------
! This module handles the storage and extraction of forcefield quantities.  
! It provides a flexible data structure for storing the interactions of the
! entire system at various levels of detail.  It may also be used to store
! iteractions between two subsets of the entire system.  The interaction
! information is dictated by the basic storage data type (EnergyPlus) which
! is defined in the storebase module.  This abstraction of the interaction
! information makes it possible to store just interaction energies or also 
! their associated 1st (forces), 2nd and 3rd derivatives.
!
! The hierarchy of scale is: 
!   atom (atm) -> molecule (mol) -> species (spc) -> full_system (ful)
!
! More specifically, the data structures contain interaction information for
! one of these scale levels (e.g. molecule) interacting with another, perhaps
! similar scale level (e.g. molecule or species).  This data structure is
! called a "level pair" (Store_Level_Pair).  It and its associated routines
! are given in the store module.  This module also uses the symmetric 
! collection of level pairs and their associated routines from storesym.
!
! When dealing with storage in 2D arrays, data is only stored in the upper
! triangle (m x n, n>=m) to avoid duplication.
!
! Important routines:
!   storetop_init -- initilizes the overall Forcefield_Results structure, but
!                    just the top of it, not the tree structure.
!   storetop_zero -- zero the whole storage structure
!   storetop_initcopy -- initialize and copy another structure
!   storetop_copy -- just copy values in another structure
!   storetop_extract -- extract the interaction between two system subsets 
!   storetop_display -- very useful, displays all the stored interactions in
!                       a condensed, readable format.  Good for debugging.
!   storetop_fastforces -- extract forces from storage structure
!   storetop_changenmoles -- checks the number of molecules of a certain
!                            species, increments/decrements the size if needed
!   storetop_delmol -- deletes the storage information of one molecule
!   storetop_update -- use a partial storage structure to update the 
!                      corresponding elements in a full system structure
!
! Needed Improvements:
! 1) consider storing full initialization information in here.  This would 
!    allow _init to initialize the full storage structure.  It's rather
!    disturbing that the this data type doesn't know its initilization 
!    characteristics.
! 2) work on speeding the routines, especially for MC-related operations
! 3) work on getting the only spc-spc storage level completely functional
!------------------------------------------------------------------------------

Module storetop

  Use defaults, Only: RDbl,strLen,xlstrLen,NO_OF_INTRA_POTS,INCR_STORAGE_SIZE, &
      dbgflag,lstrLen
  Use utils, Only: toreal,allocErrDisplay,deallocErrDisplay,ToUpper,int2str, &
      checkandstop,getdepth
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getdistsq 
  Use storebase, Only: EnergyPlus, storebase_scalarmult, storebase_init, &
      storebase_nderivs, storebase_display, storebase_clean, storebase_null, &
      storebase_initcopy, storebase_copy, storebase_zero, storebase_disp, &
      storebase_inc, storebase_chkequal, storebase_sumintra, storebase_totnrg, &
      storebase_subtract, storebase_add
  Use storesym, Only: Symmetric_Results, storesym_init, storesym_fastzero, &
      storesym_getforces, storesym_fastsum, storesym_initcopy, storesym_copy, &
      storesym_display, storesym_setmap, storesym_ptrs, storesym_zero, &
      storesym_sum, storesym_update, store_levels, storesym_extract, &
      storesym_allnrgs, storesym_delmol, storesym_resize, storesym_clean, &
      storesym_null, storesym_dumpmol, storesym_scalenrgs, storesym_scaleall, &
      storesym_depth, storesym_chkequal, storesym_fillsub, storesym_projectout, &
      storesym_basicptrs, storesym_subtract, storesym_add
  Use store, Only: Store_Level_Pair, store_init, store_extract, store_disp, &
      store_getforces, store_sum, store_resize, store_scalenrgs, store_scaleall,&
      store_zero, store_idbranch, store_initcopy, store_copy, store_terminate, &
      store_display, store_clean, store_chkequal, store_projectout, &
      store_basicptr, store_subtract, store_add

  Implicit None

  Private
  Public :: Forcefield_Results, storetop_init, storetop_getforces, &
      storetop_extract, storetop_clean, storetop_level, storetop_display, &
      storetop_sum, storetop_zero, storetop_changenmoles, storetop_fastzero, &
      storetop_fastforces, storetop_fastsum, storetop_initcopy, storetop_copy,&
      storetop_setmap, storetop_ptrs, storetop_update, storetop_allnrgs, &
      storetop_chksizes, storetop_delmol, storetop_chklink, storetop_totnrg, &
      storetop_dumpmol, storetop_scalenrgs, storetop_scaleall, &
      storetop_chkequal, storetop_fillsub, storetop_projectout, &
      storetop_basicptrs, storetop_subtract, storetop_add, &
      storetop_nmolesDecr, storetop_delpossible

  !** Contains all the forcefield results for a single forcefield
  Type Forcefield_Results
    Integer                      :: nspc,nderivs,storelevel
    Logical                      :: ncoulon,coulon  !* not used yet
    Logical                      :: subsetonly         
    Integer, Dimension(3)        :: subset
    Type(EnergyPlus)             :: total  !* NOTE: total%nrg wo intra term
    Integer, Dimension(:), Pointer                :: nmoles_sizes
    Integer, Dimension(:), Pointer                :: nmoles_system
    Type(Symmetric_Results)                       :: ncoul,coul !* spc-spc nrgs
    Logical, Dimension(:), Pointer                :: intraon    !* on/off
    Type(Store_Level_Pair), Dimension(:), Pointer :: intra      !* intra nrgs
  End Type Forcefield_Results
  !** 'nmoles_sizes' contains the number of molecules in each species, 
  !** corresponding to the dimensioned sizes of the storage structure
  !** 'nmoles_system' contains the number of molecules of each species
  !** actually in the system whose interaction are stored.

  !** The 'subsetonly' flag indicates that the storage structure contains
  !** only the interaction of a subset of the system with the whole system.
  !** In this case, the molecule number (and perhaps later atom numbers)
  !** of the first level in the level pairs may be mapped into the first
  !** index.  This keeps the storage size to a minimum, but requires that
  !** we keep track of which subset the storage represents.  The subset
  !** that is interation with the whole system is stored in 'subset'

Contains
  !----------------------------------------------------------------------------
  ! Initializes the top of a forcefield results storage structure.  The 
  ! remainder of the structure is initilized by call store_init at various
  ! levels.
  ! Requires:  ffresults -- forcefield results structure
  !            nspc -- number of species
  !            storelevel -- string indicating desired storage level
  !            nderivs -- opt. number of derivatives to store at top (0,1,2)
  !----------------------------------------------------------------------------
  Subroutine storetop_init(ffresults,nspc,nmoleslist,storelevel,nderivs)
    Type(Forcefield_Results), Intent(Out)  :: ffresults
    Integer, Intent(In)                    :: nspc
    Integer, Dimension(:), Intent(In)      :: nmoleslist
    Character(len=3), Intent(In)           :: storelevel
    Integer, Intent(In), Optional          :: nderivs

    Integer     :: spc,error

    ffresults%nspc = nspc
    ffresults%nderivs = nderivs
    ffresults%storelevel = storetop_level(storelevel)

    !** Set the characteristics flags
    ffresults%ncoulon = .False.
    ffresults%coulon = .False.
    ffresults%subsetonly = .False.
    ffresults%subset = 0

    !** Allocate the total number of molecules of each species SIZE list
    Allocate(ffresults%nmoles_sizes(nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__, &
        'ffresults%nmoles_sizes')
    Do spc = 1,nspc
      ffresults%nmoles_sizes(spc) = nmoleslist(spc)
    End Do

    !** Allocate the total number of molecules of each species IN SYSTEM list
    Allocate(ffresults%nmoles_system(nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__, &
        'ffresults%nmoles_system')
    Do spc = 1,nspc
      ffresults%nmoles_system(spc) = nmoleslist(spc)
    End Do

    !** Allocate the sum total system energy storage
    If (Present(nderivs)) Then
      Call storebase_init(ffresults%total,nderivs,.True.)
    Else
      Call storebase_init(ffresults%total,0,.True.)
    End If

    !** Deal with the intramolecular storage first
    !** Allocate the intramolecular on/off flags
    Allocate(ffresults%intraon(nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ffresults%intraon')
    ffresults%intraon = .False.

    !** Allocate the intramolecular spc-based storage
    Allocate(ffresults%intra(nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ffresults%intra')

    !** Initialize the symmetric coulombic and non-coulombic storage
    Call storesym_init(ffresults%ncoul,ffresults%nspc,ffresults%nderivs)
    Call storesym_init(ffresults%coul,ffresults%nspc,ffresults%nderivs)

  End Subroutine storetop_init

  !----------------------------------------------------------------------------
  ! Nullify the pointers in the structure
  ! Requires:  storage -- forcefield results structure
  !----------------------------------------------------------------------------
  Subroutine storetop_null(storage)
    Type(Forcefield_Results), Intent(InOut) :: storage

    Call storebase_null(storage%total)
    Nullify(storage%nmoles_sizes)
    Nullify(storage%nmoles_system)

!    Call storesym_null(storage%ncoul)
!    Call storesym_null(storage%coul)

    Nullify(storage%intraon)
    Nullify(storage%intra)

  End Subroutine storetop_null

  !----------------------------------------------------------------------------
  ! Copies a forcefield structure, or a part of it, into a new structure.  This
  ! can be used to copy results or prepare a structure for holding subset
  ! interactions.
  ! Requires:  image -- new forcefield results structure (must be deallocated)
  !            old -- forcefield results structure to copy
  !            subset -- integer array identifying subset to copy
  !----------------------------------------------------------------------------
  Subroutine storetop_initcopy(image,old,subset)
    Type(Forcefield_Results), Intent(Out) :: image
    Type(Forcefield_Results), Intent(In)  :: old
    Integer, Dimension(3), Intent(In)     :: subset

    Integer                  :: i,spc,error

    !** Nullify the parts of the copy
    Call storetop_null(image)

    spc = subset(1)
    image%subsetonly = old%subsetonly
    image%subset = old%subset
    image%nspc = old%nspc
    image%nderivs = old%nderivs

    image%storelevel =  old%storelevel

    !** Allocate the total number of molecules of each species SIZE list
    Allocate(image%nmoles_sizes(image%nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'image%nmoles_sizes')
    Do i = 1,image%nspc
      image%nmoles_sizes(i) = old%nmoles_sizes(i)
    End Do

    !** Allocate the total number of molecules of each species SYSTEM list
    Allocate(image%nmoles_system(image%nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'image%nmoles_system')
    Do i = 1,image%nspc
      image%nmoles_system(i) = old%nmoles_system(i)
    End Do

    !** Initialize and copy the total interactions
    Call storebase_initcopy(image%total,old%total)

    !** Allocate intramolecular storage
    Allocate(image%intraon(image%nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'image%intraon')
    image%intraon = .False.
    Allocate(image%intra(image%nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'image%intra')

    !** Copy intramolecular storage
    If (subset(1) == 0) Then
      Do i = 1,image%nspc
        image%intraon(i) = old%intraon(i)
        Call store_initcopy(image%intra(i),old%intra(i),.True.)
      End Do
    Else
      image%intraon(spc) = old%intraon(spc)
      Do i = 1,image%nspc
        If (i == spc) Then
          Call store_initcopy(image%intra(i),old%intra(spc),.True.)
        Else 
          Call store_terminate(image%intra(i),0)
        End If
      End Do
    End If

    !** Initialize and copy the symmetric coulombic and non-coulombic storage
    Call storesym_initcopy(image%ncoul,old%ncoul,subset)
    Call storesym_initcopy(image%coul,old%coul,subset)

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Call storetop_display(image,.False.,2,6)

  End Subroutine storetop_initcopy

  !----------------------------------------------------------------------------
  ! Copies a forcefield structure, or a part of it, into another structure.
  ! The structure may be initialized to different shapes.
  ! Requires:  dest -- receiving forcefield structure (already init)
  !            orig -- forcefield results structure to copy
  !----------------------------------------------------------------------------
  Logical Function storetop_copy(dest,orig)
    Type(Forcefield_Results), Intent(InOut) :: dest
    Type(Forcefield_Results), Intent(In)    :: orig

    Integer                  :: i,spc
    Integer                  :: error1,error2

    storetop_copy = .True.

    !** check the number of species
    If (orig%nspc /= dest%nspc) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' forcefield results structure do not having matching nspc'
      Stop      
    End If

    !** check the number of derivatives
    If (orig%nderivs /= dest%nderivs) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' forcefield results structure do not having matching nderivs'
      Stop      
    End If

    Call storebase_copy(dest%total,orig%total)

    !** Copy intramolecular interactions
    Do spc = 1,orig%nspc
      If (.Not. orig%intraon(spc)) Cycle
      storetop_copy = store_copy(dest%intra(spc),orig%intra(spc))
!      Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storetop_copy  !**DEBUG
      If (.Not. storetop_copy) Return
    End Do

    !** Copy the non-coulombic and coulombic interactions
    storetop_copy = storesym_copy(dest%coul,orig%coul)
!    Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storetop_copy  !**DEBUG_INSERT
    If (.Not. storetop_copy) Return
    storetop_copy = storesym_copy(dest%ncoul,orig%ncoul)
!    Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storetop_copy  !**DEBUG_INSERT
    If (.Not. storetop_copy) Return

  End Function storetop_copy

  !----------------------------------------------------------------------------
  ! "Updates" one forcefield structure using another or a part of another.
  ! Used during acceptance of MC moves to update the full forcefield storage
  ! Requires:  dest -- receiving forcefield structure (already init)
  !            orig -- forcefield results structure to copy
  !            subset -- specification of which subset interaction is stored
  !            skip_intra -- skips the update of the intramolecular part
  !----------------------------------------------------------------------------
  Logical Function storetop_update(dest,orig,subset,skip_intra)
    Type(Forcefield_Results), Intent(InOut) :: dest
    Type(Forcefield_Results), Intent(In)    :: orig
    Integer, Dimension(:), Intent(In)       :: subset
    Logical, Intent(In)                     :: skip_intra

    Integer                  :: i,spc
    Integer                  :: error1,error2

    !** Set default
    storetop_update = .True.

    !** Check the number of species
    If (orig%nspc /= dest%nspc) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' forcefield results structure do not having matching nspc'
      Stop      
    End If

    !** Check the number of derivatives
    If (orig%nderivs /= dest%nderivs) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' forcefield results structure do not having matching nderivs'
      Stop      
    End If

    !** Update the destination structure using the origin structure
    If (.Not. orig%subsetonly) Then
      storetop_update = storetop_copy(dest,orig)

    Else  !** update using a mol-system subset
      !** Copy intramolecular interactions
      If (.Not. skip_intra) Then
        Do spc = 1,orig%nspc
          If (.Not. orig%intraon(spc)) Cycle
          storetop_update = store_copy(dest%intra(spc)%mi(subset(2)), &
              orig%intra(spc)%mi(1))
          If (.Not. storetop_update) Return
        End Do
      End If
  
      !** Copy the non-coulombic and coulombic interactions
      storetop_update = storesym_update(dest%ncoul,orig%ncoul, &
          orig%nmoles_system,subset)
      If (.Not. storetop_update) Return
  
      storetop_update = storesym_update(dest%coul,orig%coul, &
          orig%nmoles_system,subset)
      If (.Not. storetop_update) Return
    End If

  End Function storetop_update

  !----------------------------------------------------------------------------
  ! Set the subset mapping.  When the storage structure does not contain the 
  ! entire set of interactions within the system, the 'subsetonly' flag should
  ! be True.  This means that the interactions stored are only those of the
  ! specified subset with the entire system.
  ! Requires:  storage -- Forcefield results storage structure
  !            mapflag -- new setting 
  !            subset -- new subset specification
  !----------------------------------------------------------------------------
  Subroutine storetop_setmap(storage,mapflag,subset)
    Type(Forcefield_Results), Intent(InOut) :: storage
    Logical, Intent(In)                     :: mapflag
    Integer, Dimension(3), Intent(In)       :: subset

    !** Return immediately without setting if the subset is the whole system
    If (subset(1) == 0) Return

    storage%subsetonly = mapflag
    storage%subset = subset
    Call storesym_setmap(storage%coul,mapflag,subset)
    Call storesym_setmap(storage%ncoul,mapflag,subset)

  End Subroutine storetop_setmap

  !----------------------------------------------------------------------------
  ! Returns necessary pointer(s) to EnergyPlus structures at the desired 
  ! species-molecule-atom levels.  
  ! Requires:  storage -- structure to extract pointer(s) from
  !            intra -- flag indicating if intramolecular pointer(s) sought
  !            coul -- flag indicating if coulombic/noncoul pointer(s) sought
  !            subset1 -- 1st species-molecule-atom subset 
  !            subset2 -- 2nd species-molecule-atom subset 
  !            abptr -- pointer containing interactions on subset1 from subset2
  !            baptr -- pointer containing interactions on subset1 from subset2
  ! a "subset" is a 1D array (size > 0) with indices: species,molecule,atom
  !----------------------------------------------------------------------------
  Subroutine storetop_basicptrs(storage,intra,coul,subset1,subset2,abptr,baptr)
    Type(Forcefield_Results), Intent(InOut)  :: storage
    Logical, Intent(In)                      :: intra,coul
    Integer, Dimension(3), Intent(In)        :: subset1,subset2
    Type(EnergyPlus), Pointer                :: abptr
    Type(EnergyPlus), Pointer, Optional      :: baptr

    Logical                 :: success
    Integer, Dimension(3)   :: path

    If (intra) Then   !** Get Intramolecular storage access
      path = (/subset1(2:3),0/)
      success = store_basicptr(storage%intra(subset1(1)),path,path,abptr)
      Call checkandstop(success,__FILE__,__LINE__, &
          ' Could not find specified intramolecular storage pointer')

    Else  !** Get Non-Coulombic or Coulombic access from symmetric part
      If (Present(baptr)) Then
        If (coul) Then
          success = storesym_basicptrs(storage%coul,subset1,subset2,abptr,baptr)
        Else
          success = storesym_basicptrs(storage%ncoul,subset1,subset2,abptr,baptr)
        End If 

      Else
        If (coul) Then
          success = storesym_basicptrs(storage%coul,subset1,subset2,abptr)
        Else
          success = storesym_basicptrs(storage%ncoul,subset1,subset2,abptr)
        End If
      End If

    End If

  End Subroutine storetop_basicptrs

  !----------------------------------------------------------------------------
  ! Returns necessary pointer(s) to level-pair structures at the desired 
  ! species-molecule-atom levels.  Used in forcefield routines to access the
  ! correct part of the storage structure.
  ! Requires:  storage -- structure to extract pointer(s) from
  !            intra -- flag indicating if intramolecular pointer(s) sought
  !            coul -- flag indicating if coulombic/noncoul pointer(s) sought
  !            subset1 -- 1st species-molecule-atom subset 
  !            subset2 -- 2nd species-molecule-atom subset 
  !            indices -- correct array indices in abptr, baptr
  !            abptr -- pointer containing interactions on subset1 from subset2
  !            baptr -- pointer containing interactions on subset1 from subset2
  ! a "subset" is a 1D array (size > 0) with indices: species,molecule,atom
  !----------------------------------------------------------------------------
  Subroutine storetop_ptrs(storage,intra,coul,subset1,subset2, &
      indices,abptr,baptr)
    Type(Forcefield_Results), Intent(InOut)                 :: storage
    Logical, Intent(In)                                     :: intra,coul
    Integer, Dimension(3), Intent(In)                       :: subset1,subset2
    Integer, Dimension(2), Intent(Out)                      :: indices    
    Type(Store_Level_Pair), Dimension(:), Pointer           :: abptr
    Type(Store_Level_Pair), Dimension(:), Pointer, Optional :: baptr

    Integer                      :: depth1,branch
    Type(EnergyPlus), Pointer    :: nrgptr
    
    If (intra) Then   !** Get Intramolecular storage access
      If (Present(baptr)) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
            ' only one passed pointer need for intramolecular storage access'
        Stop      
      End If

      depth1 = storesym_depth(subset1)

      Select Case(depth1)
      Case (1)
        indices(1) = 0
        abptr => storage%intra(subset1(1))%mi
  
      Case (2)
        If (storage%subsetonly) Then
          indices(1) = 1
          abptr => storage%intra(subset1(1))%mi
        Else
          indices(1) = subset1(2)
          abptr => storage%intra(subset1(1))%mi
        End If

      Case (3)
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) 'NOT TESTED'
        Write(*,*) "should it be %im OR %mi ?"
        Stop
        If (storage%subsetonly) Then
          indices(1) = subset1(3)
          abptr => storage%intra(subset1(1))%im(1)%im
        Else
          indices(1) = subset1(3)
          abptr => storage%intra(subset1(1))%im(subset1(2))%im
        End If

      Case Default
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Stop
      End Select

    Else  !** Get Non-Coulombic or Coulombic access from symmetric part
      If (Present(baptr)) Then
        If (coul) Then
          Call storesym_ptrs(storage%coul,subset1,subset2,branch, &
              indices,nrgptr,abptr,baptr)
        Else
          Call storesym_ptrs(storage%ncoul,subset1,subset2,branch, &
              indices,nrgptr,abptr,baptr)
        End If 

      Else
        If (coul) Then
          Call storesym_ptrs(storage%coul,subset1,subset2,branch,indices, &
              nrgptr,abptr)
        Else
          Call storesym_ptrs(storage%ncoul,subset1,subset2,branch,indices, &
              nrgptr,abptr)
        End If
      End If

    End If

  End Subroutine storetop_ptrs

  !----------------------------------------------------------------------------
  ! Zeros the storage structure using calls to other store modules
  ! Requires:  ffresults -- forcefield results structure
  !            skip_intra -- True => skip intramolecular interactions
  !----------------------------------------------------------------------------
  Subroutine storetop_zero(ffresults,skip_intra)
    Type(Forcefield_Results), Intent(InOut)  :: ffresults
    Logical, Intent(In)                      :: skip_intra

    Integer              :: i,j

    Call storebase_zero(ffresults%total)

    If (.Not. skip_intra) Then
      Do i = 1,ffresults%nspc
        Call store_zero(ffresults%intra(i),.True.)
      End Do
    End If

    Call storesym_zero(ffresults%ncoul)
    Call storesym_zero(ffresults%coul)

  End Subroutine storetop_zero

  !----------------------------------------------------------------------------
  ! Zeros the storage structure using calls to other store modules.  This is 
  ! a faster, but less general version of storetop_zero
  ! Requires:  ffresults -- forcefield results structure
  !            skip_intra -- True => skip intramolecular interactions
  !----------------------------------------------------------------------------
  Subroutine storetop_fastzero(ffresults,skip_intra)
    Type(Forcefield_Results), Intent(InOut)  :: ffresults
    Logical, Intent(In)                      :: skip_intra

    Integer              :: i,j

    !** zero the top
    Call storebase_zero(ffresults%total)

    !** zero the intramolecular part
    If (.Not. skip_intra) Then
      Do i = 1,ffresults%nspc
        Call store_zero(ffresults%intra(i),.True.)
      End Do
    End If

    !** zero the non-coulombic and coulombic parts
    Call storesym_fastzero(ffresults%ncoul)
    Call storesym_fastzero(ffresults%coul)

  End Subroutine storetop_fastzero

  !----------------------------------------------------------------------------
  ! Returns the total energy.  Note that this quantity DOES NOT contain the
  ! intramolecular energies.
  ! Requires:  ffresults -- forcefield results structure
  !            intra -- if True, include intramolecular interactions in total
  !            nrgs -- energy components, intra included, must be initialized
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function storetop_totnrg(ffresults,intra,nrgs)
    Type(Forcefield_Results), Intent(In)      :: ffresults
    Logical, Intent(In)                       :: intra
    Type(EnergyPlus), Intent(InOut), Optional :: nrgs

    !** Add total intramolecular energy if requested
    storetop_totnrg = storebase_totnrg(ffresults%total,intra)

    !** Return the energy components if necessary
    If (Present(nrgs)) Then
      Call storebase_copy(nrgs,ffresults%total)
    End If

  End Function storetop_totnrg

  !----------------------------------------------------------------------------
  ! Extracts a specified energy from the forcefield results storage structure.
  ! Returns False if the quantity is not available.
  ! Requires:  ffresults -- forcefield results structure
  !            intra -- flag indicating if intramolecular interactions sought
  !            ncoul -- flag indicating if Non-Coulombic interactions sought
  !            coul -- flag indicating if coulombic interactions sought
  !            info -- returned energy plus possible derivatives structure
  !            subset1 -- 1st spc,mol,atm specifier
  !            subset2 -- 2nd spc,mol,atm specifier
  ! a "subset" is a 1D array (size > 0) with indices: species,molecule,atom
  ! Usage: the spc,molec or atm designators can also be set to zero.  This
  !        will tell the routine not descend to that level.  For example,
  !        spc1=1, molec1=0, spc2=2, molec2=0 => get overall interaction of 
  !        species 1 with species 2.
  !----------------------------------------------------------------------------
  Logical Function storetop_extract(ffresults,intra,ncoul,coul,info, &
      subset1,subset2)
    Type(Forcefield_Results), Intent(In)         :: ffresults
    Logical, Intent(In)                          :: intra,ncoul,coul
    Type(EnergyPlus), Intent(InOut)              :: info
    Integer, Dimension(3), Intent(In)            :: subset1 
    Integer, Dimension(3), Intent(In), Optional  :: subset2

    Integer                          :: spc,spc1,spc2,depth
    Integer, Dimension(4,2)          :: path1,path2,zeropath
    Character(len=xlstrLen)          :: string

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'BEFORE EXTRACTION'
    Call storetop_display(ffresults,.True.,2,6)
!    Call storetop_display(ffresults,.False.,2,6)
#endif

    !** Check depth desired, return false if not available
    depth = getdepth(subset1)

    If (Present(subset2)) depth = Max(depth,getdepth(subset2))
    If ((ffresults%storelevel - 1) < depth) Then
      storetop_extract = .False.
      Return
    End If

    !** Set default
    storetop_extract = .True.

    !** Zero the output structure
    Call storebase_zero(info)

    !** Relabel the subsets as paths, indices will be peeled off later
    zeropath = 0  !** full system
    path1(:,1) = (/subset1,0/)
    path1(:,2) = (/subset1,0/)
    path2 = 0
    If (Present(subset2)) Then
      path2(:,1) = (/subset2,0/)
      path2(:,2) = (/subset2,0/)
    End If
    spc1 = path1(1,1)
    spc2 = path2(1,1)

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Seeking interaction energy for subsets: ',subset1,'   ', &
        subset2,'      ',intra,ncoul,coul
    Write(*,'(a,4i2,a,4i2)') 'path1: ',path1(:,1),'   path2: ',path2(:,1)
    Write(*,'(a,4i2,a,4i2)') '       ',path1(:,2),'          ',path2(:,2)
#endif

    !** Handle return for simple system-system case
    If ((spc1 == 0).And.(spc2 == 0)) Then
      If ((intra) .And. (coul) .And. (ncoul)) Then
        Call storebase_copy(info,ffresults%total)
        Return
      End If
    End If

    !** Get the intramolecular part of the interactions
    If (intra) Then
      If (spc1 == 0) Then
        Do spc = 1,ffresults%nspc
          If (.Not. ffresults%intraon(spc)) Cycle
          storetop_extract = store_extract(ffresults%intra(spc),path1, &
              zeropath,info)
          If (.Not. storetop_extract) Return
        End Do
      Else

        If (ffresults%intraon(spc1)) Then
          storetop_extract = store_extract(ffresults%intra(spc1), &
              path1(2:,:),path2(2:,:),info)
          If (.Not. storetop_extract) Return
        End If
      End If
    End If

    !** Get the Non-Coulombic part of the interations
    If (ncoul) Then
      storetop_extract = storesym_extract(ffresults%ncoul,info,path1,path2)
!      Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storetop_extract  !**DEBUG
      If (.Not. storetop_extract) Return
    End If

    !** Get the Coulombic part of the interations
    If (coul) Then
      storetop_extract = storesym_extract(ffresults%coul,info,path1,path2)
!      Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storetop_extract  !**DEBUG
      If (.Not. storetop_extract) Return
    End If

#ifdef DEBUG
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'AFTER EXTRACTION'    
      Call storebase_display(info,2,6)
      Call storetop_display(ffresults,.True.,2,6)
      Stop
#endif

  End Function storetop_extract

  !-------------------------------------------------------------------------  
  ! Projects out a specified component of a force vector
  ! Requires:  storage -- Forcefield results storage structure
  !            subset -- subset identifier (spc,mol,atm)
  !            vec -- unit vector along which to remove force component
  !-------------------------------------------------------------------------    
  Subroutine storetop_projectout(storage,subset,vec)
    Type(Forcefield_Results), Intent(InOut) :: storage
    Integer, Dimension(:), Intent(In)       :: subset
    Type(VecType), Intent(In)               :: vec

    Integer                  :: depth
    Logical                  :: success
    Integer, Dimension(3,2)  :: path1,zeropath
    
    !** Check depth
    If (getdepth(subset) /= 3) Then
      Write(0,'(2a,i6,a,i3)') __FILE__,":",__LINE__, &
          ' currently can only handle non-atm-depth subsets, think about this?'
      Stop        
    End If

    !** Project intramolecular interactions
    If (storage%intraon(subset(1))) Then
      zeropath = 0
      path1(:,1) = (/subset(2:3),0/)
      path1(:,2) = (/subset(2:3),0/)
      success = store_projectout(storage%intra(subset(1)),path1,zeropath,vec)  
      If (.Not. success) Then
        Write(0,'(2a,i6,a,i3)') __FILE__,":",__LINE__, &
            ' unexpected difficulty during intramolecular force projection'
        Stop        
      End If
    End If

    !** Project Non-Coulombic and Coulombic interactions
    Call storesym_projectout(storage%ncoul,subset,vec)  
    Call storesym_projectout(storage%coul,subset,vec)  

  End Subroutine storetop_projectout

  !----------------------------------------------------------------------------
  ! Extract the gradients on each atom of a whole species from the storage 
  ! structure.  This is needed for MD simulations.  This routine requests
  ! 1D arrays of gradients from store_getforces.  Each array contains a gradient
  ! for each atom of the species.  The gradients are summed for each atom
  ! and reshaped into a 2D array (atom,molec) that is returned.
  ! Requires:  storage -- full storage structure
  !            spc -- species number
  !            nmoles -- number of molecules of species
  !            natoms -- number of atoms in a molecule of species
  !            flags -- flags for non-coulomb, coulomb and intra grads
  ! NOTE: nmoles, natoms MUST match storage structure, however, 
  ! the sizes of grads don't need to match (can be larger)
  !----------------------------------------------------------------------------
  Subroutine storetop_getforces(storage,spc,nmoles,natoms,grads)
    Type(Forcefield_Results), Intent(In)       :: storage
    Integer, Intent(In)                        :: spc,nmoles,natoms
    Type(VecType), Dimension(:,:), Intent(Out) :: grads

    Integer                                 :: a,m,n,s,spc1,spc2,idx
    Type(VecType), Dimension(natoms*nmoles) :: newgrads

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'WARNING: I am not 100% certain storetop_getforces works'
    Stop

    !** Calculate total number of atoms in spc
    n = nmoles*natoms
!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,*) 'nmoles,natoms ',nmoles,natoms,n

    !** reform 1D gradients from each of the spc-spc boxes into 2D grads array
    Do s = 1,storage%nspc

      !** pick upper triangle pair of the 2D spc-spc storage arrays
      If (s >= spc) Then
        spc1 = spc
        spc2 = s
      Else
        spc1 = spc
        spc2 = s
      End If
      
!      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!      Write(*,*) spc1,spc2
!      Write(*,*) storage%ncoul%on(spc1,spc2)

      !** non-coulombic 
      If (storage%ncoul%on(spc1,spc2)) Then
        idx = 0
        newgrads = store_getforces(storage%ncoul%ab(spc1,spc2),n,.False.)
        Do m = 1,nmoles
          Do a = 1,natoms
            idx = idx + 1
            grads(a,m) = grads(a,m) + newgrads(idx)
          End Do
        End Do

        idx = 0
        If (store_idbranch(storage%ncoul%ba(spc1,spc2)) >= 0) Then
          newgrads = store_getforces(storage%ncoul%ba(spc1,spc2),n,.True.)
          Do m = 1,nmoles
            Do a = 1,natoms
              idx = idx + 1
              grads(a,m) = grads(a,m) + newgrads(idx)
            End Do
          End Do
        End If
      End If

      !** coulombic 
      If (storage%coul%on(spc1,spc2)) Then
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) 'Shouldnt be here, not ready yet'
        stop
        idx = 0
        newgrads = store_getforces(storage%coul%ab(spc1,spc2),n,.False.)
        Do m = 1,nmoles
          Do a = 1,natoms
            idx = idx + 1
            grads(a,m) = grads(a,m) + newgrads(idx)
          End Do
        End Do

        idx = 0
        newgrads = store_getforces(storage%coul%ba(spc1,spc2),n,.False.)
        Do m = 1,nmoles
          Do a = 1,natoms
            idx = idx + 1
            grads(a,m) = grads(a,m) + newgrads(idx)
          End Do
        End Do
      End If

    End Do

    !** intramolecular
    If (storage%intraon(spc)) Then
      idx = 0
      newgrads = store_getforces(storage%intra(spc),n,.False.)
      Do m = 1,nmoles
        Do a = 1,natoms
          idx = idx + 1
          grads(a,m) = grads(a,m) + newgrads(idx)
        End Do
      End Do
    End If

  End Subroutine storetop_getforces

  !----------------------------------------------------------------------------
  ! An attempt to speed _getforces by directly accessing the storage.  Extracts
  ! the forces from the storage structure.
  ! Requires:  storage -- full storage structure
  !            spc -- species number
  !            nmoles -- number of molecules of species
  !            natoms -- number of atoms in a molecule of species
  !            forces -- the output forces on each (atom,molecule)
  !----------------------------------------------------------------------------
  Subroutine storetop_fastforces(storage,spc,nmoles,natoms,forces)
    Type(Forcefield_Results), Intent(In)         :: storage
    Integer, Intent(In)                          :: spc,nmoles,natoms
    Type(VecType), Dimension(:,:), Intent(InOut) :: forces

    Integer                                 :: a,m,n,s,spc1,spc2,idx
    Type(VecType), Dimension(natoms*nmoles) :: newforces

    !** Calculate total number of atoms in spc
    n = nmoles*natoms

    !** non-coulombic 
    Call storesym_getforces(storage%ncoul,spc,forces)

    !** coulombic 
    Call storesym_getforces(storage%coul,spc,forces)

    !** intramolecular
    If (storage%intraon(spc)) Then
      idx = 0
      !** size of intra might be more than required; values from only first
      !** nmoles are required 
      newforces = store_getforces(storage%intra(spc),n,.False.,nmoles)
      Do m = 1,nmoles
        Do a = 1,natoms
          idx = idx + 1
          forces(a,m) = forces(a,m) + newforces(idx)
        End Do
      End Do
    End If

  End Subroutine storetop_fastforces

  !----------------------------------------------------------------------------
  ! Returns the integer storage level given the string
  ! Requires:  string -- string indicating storage level
  !----------------------------------------------------------------------------
  Integer Function storetop_level(string)
    Character(len=3), Intent(In)        :: string

    Integer               :: i

    storetop_level = 0

    Do i = 1,Size(store_levels)
      If (ToUpper(string) == ToUpper(store_levels(i))) Then
        storetop_level = i
        Exit
      End If
    End Do

  End Function storetop_level

  !----------------------------------------------------------------------------
  ! Update the storage structure so that all the sums at each level are 
  ! consistent with the information at lower levels.  This is not done 
  ! automatically during the evaluations because it would lead to duplication
  ! of effort.
  ! Requires:  storage -- Forcefield results storage structure
  !----------------------------------------------------------------------------
  Subroutine storetop_sum(storage)
    Type(Forcefield_Results), Intent(InOut)  :: storage

    Integer            :: i

    Call storebase_zero(storage%total)

    !** intramolecular summing
    Do i = 1,storage%nspc
      Call store_sum(storage%intra(i),.True.)
      Call storebase_inc(storage%total,storage%intra(i)%total)
    End Do

    !** non-coulombic and coulombic summing
    Call storesym_sum(storage%ncoul)
    Call storebase_inc(storage%total,storage%ncoul%total)
    Call storesym_sum(storage%coul)
    Call storebase_inc(storage%total,storage%coul%total)

  End Subroutine storetop_sum

  !----------------------------------------------------------------------------
  ! Subtract the second base structure from the first.  Assumes that the
  ! initialization structures are identical
  ! Requires:  storage1 -- 1st Forcefield results storage structure
  !            storage2 -- 2nd Forcefield results storage structure
  !----------------------------------------------------------------------------
  Subroutine storetop_subtract(storage1,storage2)
    Type(Forcefield_Results), Intent(InOut)  :: storage1
    Type(Forcefield_Results), Intent(In)     :: storage2

    Integer          :: i

    Call storebase_subtract(storage1%total,storage2%total)

    !** Intramolecular subtraction
    Do i = 1,storage1%nspc
      Call store_subtract(storage1%intra(i),storage2%intra(i),.True.)
    End Do

    !** Non-coulombic and coulombic subtraction
    Call storesym_subtract(storage1%ncoul,storage2%ncoul)
    Call storesym_subtract(storage1%coul,storage2%coul)

  End Subroutine storetop_subtract

  !----------------------------------------------------------------------------
  ! Add the second base structure to the first.  Assumes that the
  ! initialization structures are identical
  ! Requires:  storage1 -- 1st Forcefield results storage structure
  !            storage2 -- 2nd Forcefield results storage structure
  !----------------------------------------------------------------------------
  Subroutine storetop_add(storage1,storage2)
    Type(Forcefield_Results), Intent(InOut)  :: storage1
    Type(Forcefield_Results), Intent(In)     :: storage2

    Integer          :: i

    Call storebase_add(storage1%total,storage2%total)

    !** Intramolecular addition
    Do i = 1,storage1%nspc
      Call store_add(storage1%intra(i),storage2%intra(i),.True.)
    End Do

    !** Non-coulombic and coulombic addition
    Call storesym_add(storage1%ncoul,storage2%ncoul)
    Call storesym_add(storage1%coul,storage2%coul)

  End Subroutine storetop_add

  !----------------------------------------------------------------------------
  ! This is a faster version of storetop_fast.  The speed increase comes from
  ! using the storesym_fastsum routine which assumes that mm contributions 
  ! are already summed and does all the work in a single routine.
  ! Requires:  storage -- Forcefield results storage structure
  !----------------------------------------------------------------------------
  Subroutine storetop_fastsum(storage)
    Type(Forcefield_Results), Intent(InOut)  :: storage

    Integer            :: i

    Call storebase_zero(storage%total)

    !** intramolecular summing 
    !** (will need to change, if the intramolecular can ever truly be SPC-detail)
    Do i = 1,storage%nspc
      Call storebase_zero(storage%intra(i)%total)
      Call store_sum(storage%intra(i),.True.)
      Call storebase_inc(storage%total,storage%intra(i)%total)
    End Do

    !** non-coulombic and coulombic summing 
    Call storesym_fastsum(storage%ncoul)
    Call storesym_fastsum(storage%coul)

    !** update total with non-coulombic and coulombic
    Call storebase_inc(storage%total,storage%ncoul%total)
    Call storebase_inc(storage%total,storage%coul%total)

  End Subroutine storetop_fastsum

  !-------------------------------------------------------------------------
  ! Scales all the species-species based energies
  ! Requires:  storage  -- Forcefield results storage structure
  !            factors  -- 2D array of spc-spc energy scaling factors
  ! Optional:  whichnrg -- Logical, dimension(3), with entries corresponding
  !                        to which energies to scale. By default, all are
  !                        scaled. The format is (/ Noncoulombic, Coulombic,
  !                        Intramolecular /). Passing (/ T,T,F /) would 
  !                        scale noncoulombic and coulombic forces but
  !                        not intramolecular
  !-------------------------------------------------------------------------
  Subroutine storetop_scalenrgs(storage,factors,whichnrg)
    Type(Forcefield_Results), Intent(InOut)         :: storage
    Real(kind = RDbl), Dimension(:,:), Intent(In)   :: factors
    Logical, Dimension(3), Intent(In), Optional     :: whichnrg

    Integer         :: i,j
    Logical         :: coul, noncoul, intra

    !** Check for optional arguments
    coul = .True.
    noncoul = .True.
    intra = .True.
    If (Present(whichnrg)) noncoul = whichnrg(1)
    If (Present(whichnrg)) coul = whichnrg(2)
    If (Present(whichnrg)) intra = whichnrg(3)

    !** Size checking
    If (Size(factors,1) /= storage%nspc) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' Wrongly sized scaling factors array passed, index 1'
      Stop            
    End If
    If (Size(factors,2) /= storage%nspc) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' Wrongly sized scaling factors array passed, index 2'
      Stop            
    End If
    
    If (intra) Then
      Do i = 1,storage%nspc
        !** scale the intramolecular interactions
        If (storage%intraon(i)) Then
          Call store_scalenrgs(storage%intra(i),factors(i,1))  !** HACK!!
        End If
      End Do
    End If

    !** scale the Non-Coulombic and Coulombic interactions
    If (coul) Call storesym_scalenrgs(storage%ncoul,factors)
    If (noncoul) Call storesym_scalenrgs(storage%coul,factors)

  End Subroutine storetop_scalenrgs

  !-------------------------------------------------------------------------
  ! Scales all the species-species based interactions
  ! Requires:  storage -- Forcefield results storage structure
  !            factor -- scaling factor
  !-------------------------------------------------------------------------
  Subroutine storetop_scaleall(storage,factor)
    Type(Forcefield_Results), Intent(InOut)   :: storage
    Real(kind = RDbl), Intent(In)             :: factor

    Integer         :: i,j

    Do i = 1,storage%nspc
      !** scale the intramolecular interactions
      If (storage%intraon(i)) Then
        Call store_scaleall(storage%intra(i),factor)  
      End If
    End Do

    !** scale the Non-Coulombic and Coulombic interactions
    Call storesym_scaleall(storage%ncoul,factor)
    Call storesym_scaleall(storage%coul,factor)

  End Subroutine storetop_scaleall

  !-------------------------------------------------------------------------     
  ! Returns all the species-based energies
  ! Requires:  storage -- Forcefield results storage structure
  !            noncoul -- 2D array of spc-spc Non-Coulombic energies
  !            coul -- 2D array of spc-spc Coulombic energies
  !            intra -- 1D array of spc intramolecular energies
  !-------------------------------------------------------------------------     
  Subroutine storetop_allnrgs(storage,noncoul,coul,intra)
    Type(Forcefield_Results), Intent(In)            :: storage
    Real(kind = RDbl), Dimension(:,:), Intent(Out)  :: noncoul
    Real(kind = RDbl), Dimension(:,:), Intent(Out)  :: coul
    Type(EnergyPlus), Dimension(:), Intent(InOut)   :: intra

    Integer           :: i

    !** Copy the intramolecular energies
    Do i = 1,storage%nspc
      Call storebase_copy(intra(i),storage%intra(i)%total)
    End Do

    !** Get the Non-Coulombic energies
    Call storesym_allnrgs(storage%ncoul,noncoul)

    !** Get the Coulombic energies
    Call storesym_allnrgs(storage%coul,coul)

  End Subroutine storetop_allnrgs

  !------------------------------------------------------------------------- 
  ! This routine compares two storage structures to make sure they think
  ! they are storing the same number of molecules.  If the 'image' structure
  ! is not consistent with the reference structure, it is resized.  
  ! Requires:  image -- forcefield results storage structure to check
  !            ref -- reference forcefield results storage structure
  !            subset -- the subset for which the image structure is valid
  !------------------------------------------------------------------------- 
  Subroutine storetop_chksizes(image,ref,subset)
    Type(Forcefield_Results), Intent(InOut) :: image
    Type(Forcefield_Results), Intent(In)    :: ref
    Integer, Dimension(3), Intent(In)       :: subset

    Integer        :: i,outunit
    Logical        :: disagree, resizeonly
    Character(len=lstrLen)   :: filename
    
    !** Update number of molecules in image system
    Do i = 1,ref%nspc
      image%nmoles_system(i) = ref%nmoles_system(i)
    End Do

    !** Check for disagreement in the sizes of the storage
    disagree = .False.
    Do i = 1,ref%nspc
      If (image%nmoles_sizes(i) /= ref%nmoles_sizes(i)) Then
        disagree = .True.
        Exit
      End if
    End Do
    If (.Not. disagree) Return

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'RESIZING image '
    Write(*,*) 'ref sizes    ',ref%nmoles_sizes
    Write(*,*) 'image sizes  ',image%nmoles_sizes

    filename = 'reference.txt'
    Write(*,*) 'Writing reference storage to: ',Trim(filename)
    outunit = file_open(filename,0)
    Call storetop_display(ref,.False.,2,outunit)
    Close(unit=outunit)
#endif

    !** Resize the image storage
    Call storetop_clean(image)  
    Call storetop_initcopy(image,ref,subset)  

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'After RESIZING the image '
    Write(*,*) 'ref sizes    ',ref%nmoles_sizes
    Write(*,*) 'image sizes  ',image%nmoles_sizes
    Call storetop_display(image,.False.,2,6)

    filename = 'image.txt'
    Write(*,*) 'Writing image storage to: ',Trim(filename)
    outunit = file_open(filename,0)
    Call storetop_display(image,.False.,2,outunit)
    Close(unit=outunit)
#endif

    !** Return if nothing else is needed, ie the image is the full system
    If (subset(1) == 0) Return

    !** Fill in the "gaps" in the subset storage.  The subset storage needs
    !** to be generic for any subset of the specified depth, but it was 
    !** copied from a specific subset.  This means that it has uninitialized
    !** self-interactions that need to be initialized as others are, "filled".
    Call storetop_fillsub(image,subset)  

    !** Label the image storage structure as a mol-map subset
    Call storetop_setmap(image,.True.,subset)

  End Subroutine storetop_chksizes

  !------------------------------------------------------------------------- 
  ! Fill in the "gaps" in the subset storage.  The subset storage needs
  ! to be generic for any subset of the specified depth, but it was 
  ! copied from a specific subset.  This means that it has uninitialized
  ! self-interactions that need to be initialized as others are, "filled".
  ! Requires:  image -- forcefield results storage structure to check
  !            subset -- the subset for which the image structure is valid
  !------------------------------------------------------------------------- 
  Subroutine storetop_fillsub(image,subset)
    Type(Forcefield_Results), Intent(InOut) :: image
    Integer, Dimension(3), Intent(In)       :: subset

    Call storesym_fillsub(image%ncoul,subset)
    Call storesym_fillsub(image%coul,subset)

  End Subroutine storetop_fillsub

  !---------------------------------------------------------------------------   
  ! This routine checks the storage space allocated for the number of 
  ! molecules and increments the space by INCR_SIZE (defaults) if needed.
  ! It also keeps track of how many molecules are actually in the storage
  ! structure.  This routine is allowed to change the counter for the number 
  ! of molecules of that species in the system by as much as one.
  !
  ! If no molecule is specified and the new number of molecules is less,
  ! then it is assumed that (1) the molecule to be deleted was the last of 
  ! that species and (2) its interactions were never entered into the full 
  ! system storage.  In this case, simply changing the number of molecules 
  ! is sufficient to erase it
  !
  ! WARNING: For situations where the storage lacks sufficient detail to
  ! simply delete a molecule's interactions, it will do as much as it can
  ! (rearrange mol-info if possible) and return no error message.  It is
  ! assumed that the appropriate interactions have already been subtracted
  ! out of the total.  This is a potential problem for symmetric portion
  ! of the storage with SPC-detail.
  ! Requires:  storage -- Forcefield results storage structure
  !            nmoles -- new number of molecules of species spc in system
  !            spc -- the species number
  !            delmol -- optional molecule to delete
  !-----------------------------------------------------------------------
  Subroutine storetop_changenmoles(storage,spc,nmoles,delmol)
    Type(Forcefield_Results), Intent(InOut) :: storage
    Integer, Intent(In)                     :: spc,nmoles
    Integer, Intent(In), Optional           :: delmol

    Integer       :: difference,sizediff
    Logical       :: increment

    !** Change the number of molecules of species in SYSTEM counter
    increment = .False.
    difference = nmoles - storage%nmoles_system(spc)

    If (difference == 0) Then
      !** Do nothing, counters and sizes are ok
      Return

    Else If (difference == -1) Then
      !** Remove molecule from storage structure if the number was specified
      If (Present(delmol)) Then
        !** Delete the molecule's interactions from the structure
        Call storetop_delmol(storage,spc,delmol)
      End If

      storage%nmoles_system(spc) = nmoles

    Else If (difference == 1) Then
      storage%nmoles_system(spc) = nmoles

      sizediff = storage%nmoles_system(spc) - storage%nmoles_sizes(spc)
      If (sizediff == 1) Then
        increment = .True.
      Else If (sizediff > 1) Then
        Write(0,'(2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(a)') 'number of molecules incremented ', &
            'by more than one, unexpected'
      End If

    Else
      Write(0,'(2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(a)') 'number of molecules incremented or decremented ', &
          'by more than one, unexpected'
      Write(0,'(a,3i3)') 'spc, new nmoles, nmoles change: ',spc,nmoles,difference
      Stop      
    End If

    If (storage%subsetonly) Then
      Write(0,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(0,'(a)') 'storetop_changenmoles should not be used on subset storage'
      Write(0,'(a)') 'to change the shape of the structure, copy instead'
      Stop
    End If

    !** Actually resize the storage structure SIZE
    If (increment) Then
      !** new storage size
      storage%nmoles_sizes(spc) = storage%nmoles_sizes(spc) + INCR_STORAGE_SIZE

#ifdef DEBUG
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'resizing storage for species number: ',spc
      Write(*,*) 'nmoles list (for each species): ',storage%nmoles_sizes
      Write(*,*) 'BEFORE non-coul/coul resizing'
!      Call storetop_display(storage,.False.,2,6)
      Call storetop_display(storage,.True.,2,6)
#endif

      !** Resize intramolecular storage first
      If (storage%intraon(spc)) Then 
        Call store_resize(storage%intra(spc),(/storage%nmoles_sizes(spc),0/))
      End If

      !** Resize symmetric storage structures
      Call storesym_resize(storage%ncoul,storage%nmoles_sizes(spc),spc)
      Call storesym_resize(storage%coul,storage%nmoles_sizes(spc),spc)

#ifdef DEBUG
      Write(*,*) 'AFTER resizing'
      Call storetop_display(storage,.False.,2,6)
!      Stop
#endif
    End If

  End Subroutine storetop_changenmoles

  !-----------------------------------------------------------------------
  ! decrease the number of molecules of a species that will be considered
  ! It does not do any checks, use storetop_changenmoles for anything more 
  ! sophisticated It is assumed that there is more space than nmoles
  ! in storage
  ! Requires:  storage -- Forcefield results storage structure
  !            spc -- the species number
  !            nmoles -- new number of molecules of species spc in system
  !-----------------------------------------------------------------------
  Subroutine storetop_nmolesDecr(storage,spc,nmoles)
    Type(Forcefield_Results), Intent(InOut) :: storage
    Integer, Intent(In)                     :: spc,nmoles

    Integer       :: difference

    difference = nmoles - storage%nmoles_system(spc)
    storage%nmoles_system(spc) = nmoles


  End Subroutine storetop_nmolesDecr

  !-------------------------------------------------------------------------  
  ! This function returns a flag saying if it is possible to delete just a
  ! molecule from the storage structure.  Sometimes the structure of the
  ! storage does not allow this.  SPC-detail with forces is an example.
  ! In that situation, attempts to delete a molecule's storage by the
  ! move-last-to-empty-spot + zero-last method do not work because of the
  ! way that ssbasic enters molecule-molecule interactions into the 
  ! symmetric portion of the storage.  If it's not possible all the 
  ! system interactions must be reevaluated.
  ! Requires:  storage -- Forcefield results storage structure
  !-------------------------------------------------------------------------    
  Logical Function storetop_delpossible(storage)
    Type(Forcefield_Results), Intent(InOut) :: storage

    storetop_delpossible = .True.

    If ((storage%storelevel < 3).And.(storage%nderivs > 0)) Then
      storetop_delpossible = .False.
    End If

  End Function storetop_delpossible

  !-------------------------------------------------------------------------  
  ! This routine removes a molecule from the storage and compresses the
  ! remaining storage to eliminate the gap.  It DOES NOT change the number
  ! of molecules in system counters, this is done in _changenmoles
  ! Requires:  storage -- Forcefield results storage structure
  !            spc -- the species number
  !            mol -- the molecule number
  !-------------------------------------------------------------------------    
  Subroutine storetop_delmol(storage,spc,mol)
    Type(Forcefield_Results), Intent(InOut) :: storage
    Integer, Intent(In)                     :: spc,mol

    Integer        :: nmoles
    Logical        :: okflag

    !** Delete one molecule's intramolecular interactions first
    If (storage%intraon(spc)) Then
      nmoles = storage%nmoles_system(spc)

      !** Copy last molecule if 'mol' is not already the last molecule
      If (mol /= nmoles) Then
        okflag = store_copy(storage%intra(spc)%mi(mol), &
            storage%intra(spc)%mi(nmoles))
        Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
      End If
      Call store_zero(storage%intra(spc)%mi(nmoles),.True.)
    End If

    !** Delete one molecule's intramolecular symmetric storage
    Call storesym_delmol(storage%ncoul,storage%nmoles_system,spc,mol)
    Call storesym_delmol(storage%coul,storage%nmoles_system,spc,mol)

  End Subroutine storetop_delmol

  !----------------------------------------------------------------------------
  ! Compares all the initialized entries in two structures.  If a unit
  ! number and indent specifier are included, the routine will dump output. 
  ! Otherwise, only the True/False answer will be returned.  Differences are
  ! reported as (storage1 - storage2).
  ! Requires:  storage1 -- 1st basic storage structure  
  !            storage2 -- 2nd basic storage structure  
  !            tolerance -- tolerance for calling equal
  !            indent -- no. of spaces from the left margin
  !            unit -- optional unit number
  !----------------------------------------------------------------------------
  Logical Function storetop_chkequal(storage1,storage2,tolerance,indent,unit)
    Type(Forcefield_Results), Intent(In) :: storage1,storage2
    Real(kind=RDbl), Intent(In)          :: tolerance
    Integer, Intent(In), Optional        :: indent,unit

    Integer                      :: i,spc
    Logical                      :: equal,dump
    Character(len=strLen)        :: string
    Character(len=xlstrLen)      :: startstring,output

    storetop_chkequal = .True.

    dump = .False.
    If ((Present(indent)).And.(Present(unit))) dump = .True.

    If (dump) Then
      Write(unit,'(30a)') (' ',i=1,indent), &
          'Checking equality of two storage structures, will report differences'
    End If

    !** Compare the total energy
    output = ''
    If (.Not. storebase_chkequal(storage1%total,storage2%total, &
        tolerance,output)) Then
      Write(unit,'(30a)') (' ',i=1,indent+2), &
          'Total interactions not equal: ',Trim(output),' (may be unsummed)'
    End If

    !** Compare the intramolecular parts 
    If (dump) Write(unit,'(30a)') (' ',i=1,indent+2),'Intramolecular problems:'
    Do spc = 1,storage1%nspc
      If (.Not. storage1%intraon(spc)) Cycle
      If (.Not. storage2%intraon(spc)) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            ' The two passed storage structures do not have same intra shape'
        Stop
      End If
      string = int2str(spc)
      Write(startstring,'(2a))') 's:',Trim(string)

      If (dump) Then
        equal = store_chkequal(storage1%intra(spc),storage2%intra(spc), &
            startstring,tolerance,indent+4,unit)
      Else
        equal = store_chkequal(storage1%intra(spc),storage2%intra(spc), &
            startstring,tolerance)
      End If
      If (.Not. equal) storetop_chkequal = .False.
    End Do

    !** Compare the Non-Coulombic interactions
    If (dump) Write(unit,'(30a)') (' ',i=1,indent+2),'Non-Coulombic problems:'
    If (dump) Then
      equal = storesym_chkequal(storage1%ncoul,storage2%ncoul,tolerance, &
          indent+4,unit)
    Else
      equal = storesym_chkequal(storage1%ncoul,storage2%ncoul,tolerance)
    End If
    If (.Not. equal) storetop_chkequal = .False.

    !** Compare the Coulombic interactions
    If (dump) Write(unit,'(30a)') (' ',i=1,indent+2),'Coulombic problems:'
    If (dump) Then
      equal = storesym_chkequal(storage1%coul,storage2%coul,tolerance, &
          indent+4,unit)
    Else
      equal = storesym_chkequal(storage1%coul,storage2%coul,tolerance)
    End If
    If (.Not. equal) storetop_chkequal = .False.

  End Function storetop_chkequal

  !-------------------------------------------------------------------------  
  ! This routine checks to see if two storage structures are illegally
  ! linked together.  It's for debugging and only checks part of the 
  ! intramolecular storage at the moment.
  ! Requires:  storage1 -- 1st Forcefield results storage structure
  !            storage2 -- 2nd Forcefield results storage structure
  !-------------------------------------------------------------------------    
  Subroutine storetop_chklink(storage1,storage2)
    Type(Forcefield_Results), Intent(InOut) :: storage1
    Type(Forcefield_Results), Intent(In)    :: storage2

    Logical             :: linked
    Real(kind=RDbl)     :: Save,value

!    storage1%ncoul%ab(1,1)%mi(1)%im(1)%total%nrg = 1.111e-4_RDbl

    value = 6.666e-1_RDbl
    save = storage1%intra(1)%mi(1)%total%intranrg(1) 
    storage1%intra(1)%mi(1)%total%intranrg(1) = value

    linked = .False.
    If (storage2%intra(1)%mi(1)%total%intranrg(1) == value) linked = .True.

    If (linked) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'Found linked storage'
      Write(*,*) 'storage1:'
      Call storetop_display(storage1,.True.,2,6)
      Write(*,*) 'storage2:'
      Call storetop_display(storage2,.True.,2,6)
      Stop
    Else
      storage1%intra(1)%mi(1)%total%intranrg(1) = save
    End If

  End Subroutine storetop_chklink

  !----------------------------------------------------------------------------
  ! Pseudo-Hack for dumping information about a molecule
  ! Requires:  storage -- Forcefield results storage structure
  !            spc -- species number
  !            mol -- molecule number
  !            indent -- no. of spaces from the left margin
  !            unitno -- unit number
  !----------------------------------------------------------------------------
  Subroutine storetop_dumpmol(storage,spc,mol,indent,unitno)
    Type(Forcefield_Results), Intent(In)  :: storage
    Integer, Intent(In)                   :: spc,mol,indent,unitno

    Call storesym_dumpmol(storage%ncoul,spc,mol,indent,unitno)

  End Subroutine storetop_dumpmol

  !----------------------------------------------------------------------------
  ! Display the level pair structure recursively.  This uses a compressed form
  ! for displaying the results.
  ! Requires:  storage -- Forcefield results storage structure
  !            skip -- True => will skip display if nrg,grad(1) < tolerance
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine storetop_display(storage,skip,indent,unitno)
    Type(Forcefield_Results), Intent(In)  :: storage
    Logical, Intent(In)                   :: skip
    Integer, Intent(In)                   :: indent
    Integer, Optional, Intent(In)         :: unitno

    Integer                           :: i,j,unit,s1,s2,spc
    Character(len=indent)             :: blank
    Character(len=strLen)             :: string1,string2
    Character(len=xlstrLen)           :: startstring,string

    blank = Repeat(' ',indent)
    
    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If
    
    If (skip) Then
      Write(unit,'(2a)') blank,'Displaying only non-zero entries in storage'
    End If

    !** Display number of molecules of each species in SYSTEM
    Write(unit,'(2a)',Advance='No') blank,'molecules of each species: '
    Do spc = 1,storage%nspc
      string1 = int2str(storage%nmoles_system(spc))
      string2 = int2str(spc)
      Write(unit,'(4a)',Advance='No') Trim(string1),'(',Trim(string2),')'
      If (spc /= storage%nspc) Write(unit,'(a)',Advance='No') ','
    End Do
    Write(unit,*)

    !** Display space for number of molecules of each species in storage
    Write(unit,'(2a)',Advance='No') blank,'space for molecules of each species: '
    Do spc = 1,storage%nspc
      string1 = int2str(storage%nmoles_sizes(spc))
      string2 = int2str(spc)
      Write(unit,'(4a)',Advance='No') Trim(string1),'(',Trim(string2),')'
      If (spc /= storage%nspc) Write(unit,'(a)',Advance='No') ','
    End Do
    Write(unit,*)

    If (storage%subsetonly) Then
      string1 = int2str(storage%subset)
      Write(unit,'(4a)') blank,'This is a subset-only structure ',& 
          'storing subset: ',Trim(string1)
    End If

    string = storebase_disp(storage%total)
    Write(unit,'(3a)') blank,'Top: ',Trim(string)

    Write(unit,'(2a)') blank,'Non-Coulombic storage structure'
    Call storesym_display(storage%ncoul,skip,indent+2,unit)

!    Write(unit,'(2a)') blank,'Coulombic storage structure'
!    Call storesym_display(storage%coul,skip,indent+2,unit)

    Write(unit,'(2a)') blank,'Intramolecular storage structure'
    Do i = 1,storage%nspc
      string = 'OFF'
      If (storage%intraon(i)) string = 'ON'
      Write(unit,'(2x,2a,i3,3x,a)') blank,'Species: ',i,Trim(string)
      If (storage%intraon(i)) Then
        Write(startstring,'(a,i1))') 's:',i
        Call store_display(storage%intra(i),startstring,skip,indent+2,unit)
      End If
    End Do

    Write(unit,*) 

  End Subroutine storetop_display

  !----------------------------------------------------------------------------
  ! Clean the forcefield results storage structure 
  ! Requires:  storage -- storage level pair structure
  !----------------------------------------------------------------------------
  Subroutine storetop_clean(storage)
    Type(Forcefield_Results), Intent(InOut)  :: storage

    Integer            :: i,j,error

    !** Clean the total energy
    Call storebase_clean(storage%total)

    !** Clean the number of molecules SIZES list
    Deallocate(storage%nmoles_sizes, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'nmoles_sizes')    

    !** Clean the number of molecules in SYSTEM list
    Deallocate(storage%nmoles_system, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'nmoles_system')    

    !** Clean the Non-Coulombic and Coulombic interactions
    Call storesym_clean(storage%ncoul)
    Call storesym_clean(storage%coul)

    !** Clean the intramolecular energy
    If (Associated(storage%intra)) Then
      Do i = 1,Size(storage%intra)
        Call store_clean(storage%intra(i))
      End Do
      Deallocate(storage%intra, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'intra')
    End If 

    Deallocate(storage%intraon, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'intraon')
    
  End Subroutine storetop_clean

End Module storetop
