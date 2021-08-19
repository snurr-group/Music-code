!------------------------------------------------------------------------------
! This module provides the data type that bundles the level pairs from the 
! store module into a symmetric form that is applicable to species, molecules
! and atoms.  The "symmetric" name comes from the storage of both "A with B"
! and "B with A" interactions for each species-species pair.  The separate
! "ab" and "ba" storage is needed because force interactions between atoms are
! equal and opposite.  We don't always store full atom-atom detail so
! it's not possible to simply store one interaction and invert it when needed.
!
! It is this routine that contains information about how interactions are
! stored.  For each species-species pair the "level", "form" and number of
! atoms and molecules of each of the two species is stored.  The "level"
! indicates how much detail of the interactions are stored.  Interactions
! can be stored on an atom-atom, molecule-molecule or species-species level.  
! More detail requires more memory and more processing time but allows 
! direct access to more fundamental forcefield interaction results.  
!
! The symmetric nature of each of the species-species storage sets is 
! indicated in the "form" field.  Currently there is an option to use one-sided
! asymmetric storage or the full two-sided (ab,ba) storage.  The asymmetric
! storage is used for situations where the interactions of the second species
! with the first are not explicitly stored.  This is the case for map-type
! interactions where there is no need to store forces on the atoms of
! fixed species because they will not be moved.  The symmetric storage form
! is used for interactions between two mobile species.
!
! More details about the symmetric (form 1) structure for two mobile species:
! 'ab' stores forces on molecule A due to interaction with molecule B
! 'ba' stores forces on molecule B due to interaction with molecule A 
! For a given species-species pair, the extent of the storage usage varies.  
! If the two species are identical, only the "upper triangle" (m1<m2) of 
! the AB and the "upper triangle" (m1<m2) of the BA storage are used for 
! energies.  Note that the use of "upper triangle" here does not imply a
! 2D array, but rather a level pair nested in another level pair.  For 
! dissimilar species pairs, the full structure must be used because all 
! molecule pairs are unique.  In some cases (Ewald), interactions
! of a molecule with itself are included.  
!
! Only the upper triangles (m x n, n>=m) of 2D arrays are used.  Also note
! that energies are stored in the 'ab' portion of the symmetric form storage
! structures.  Note that the storage of the energies in only the 'ab' side
! means that the 'ba' side could be eliminated if only the energies are being
! calculated.  This should be implemented in the future.
!
! Here's a few examples of how interactions are stored:
!  variables: spc1=species1, spc2=species2, m1=molecule1, 
!             m2=molecule2, a1=atom1, a2=atom2
!
!  1) Interaction of m1 with m2 at ATM init level (atom-atom):
!     force or energy on a1 of m1 due to a2 of m2: 
!       storage%ab(spc1,spc2)%mi(m1)%im(m2)%mm(a1,a2)
!     force or energy on a2 of m2 due to a1 of m1: 
!       storage%ab(spc1,spc2)%im(m2)%mi(m1)%mm(a2,a1)
!
!  2) Interaction of m1 with m2 at MOL init level (molecule-molecule):
!     force or energy on a1 of m1 due to m2: 
!       storage%ab(spc1,spc2)%mi(m1)%im(m2)%binfo(a1)
!     force or energy on a2 of m2 due to m1: 
!       storage%ba(spc1,spc2)%im(m2)%mi(m1)%binfo(a2)
!     NOTE: in contrast to example 1, atom-atom interactions are not stored
!
!  3) Interaction of spc1 with spc2 at SPC init level (species-species):
!     force or energy on a1 of m1 due to spc2: 
!       storage%ab(spc1,spc2)%mi(m1)%binfo(a1)
!     force or energy on a2 of m2 due to spc1: 
!       storage%ba(spc1,spc2)%im(m2)%binfo(a2)
!     if no derivatives are stored, then it has a simplier form:
!       storage%ab(spc1,spc2)%total, storage%ba(spc1,spc2)%total
!       store the spc-spc energies.
!     NOTE: no molecule-molecule or atom-atom interactions available
!
! Needed Improvements:
! 1) For the map-type form, it is assumed that 1st species is map species.
!    to fix this it be necessary to define a map form starting with %im
! 2) Implement a 3rd form type for two mobile species with just 'ab' storage.
!    this will speed calculations when only energies are used.
! 3) Should intramolecular interactions really be in the symmetric structure?
!    they are in the Ewald case, this seems inconsistent.
!------------------------------------------------------------------------------

Module storesym

  Use defaults, Only: RDbl,xlstrLen,NO_OF_INTRA_POTS,INCR_STORAGE_SIZE,dbgflag, &
      strLen
  Use utils, Only: toreal,allocErrDisplay,deallocErrDisplay,ToUpper,int2str, &
      checkandstop,swap,getdepth
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getdistsq 
  Use molecules, Only: molecules_getnatoms
  Use storebase, Only: EnergyPlus, storebase_scalarmult, storebase_init, &
      storebase_nderivs, storebase_display, storebase_clean, storebase_inc, &
      storebase_copy, storebase_zero, storebase_disp, storebase_null, &
      storebase_chkequal, storebase_subtract, storebase_add
  Use store, Only: Store_Level_Pair, store_init, store_terminate, &
      store_clean, store_getforces, store_scalenrgs, store_scaleall, &
      store_idbranch, store_initcopy, store_copy, store_zero, store_sum, &
      store_resize, store_extract, store_display, store_disp, store_chkequal, &
      store_projectout, store_basicptr, store_initcopyidx, store_subtract, &
      store_add

  Implicit None

  Private
  Public :: Symmetric_Results, Symmetric_Shape, storesym_init, &
      storesym_initsspair, storesym_ptrs, storesym_terminate, storesym_sum, &
      storesym_fastzero, storesym_getforces, storesym_fastsum, storesym_zero, &
      storesym_initcopy, storesym_display, storesym_copy, storesym_setmap, &
      storesym_update, storesym_extract, store_levels, storesym_basicptrs, &
      storesym_delmol, storesym_allnrgs, storesym_resize, storesym_nderivs, &
      storesym_clean, storesym_inc, storesym_null, storesym_dumpmol, &
      storesym_scalenrgs, storesym_scaleall, storesym_depth, storesym_chkequal, &
      storesym_fillsub, storesym_projectout, storesym_subtract, &
      storesym_add

  !** These are identifier strings for the storage levels (DON'T CHANGE)
  Character(len=3), Dimension(5), Save  :: store_levels = &
      (/'FUL','SPC','MOL','ATM','END'/)

  Type Symmetric_Shape
    Integer      :: level,form
    Logical      :: intra_also
    Integer      :: natoms1,natoms2
    Integer      :: nmoles1,nmoles2
  End Type Symmetric_Shape
  !** "level" is: 4=ATM, 3=MOL, 2=SPC
  !** if intra_also is True, then storage includes intramolecular interactions
  !** "form" takes one of the following values
  !** 0 => one-sided storage of type used in map routines
  !** 1 => two-sided, symmetric storage of type used in ssbasic routines
  !** IMPORTANT: nmoles refers to the size of the structure, this may not
  !** match how many molecules are actually in the system!!

  !** Contains storage for symmetric interactions on a species-species level
  !** 'ab' contains interactions for species A interacting with species B.
  !** 'ba' is naturally the opposite.  They are different because, at the
  !** atomic level, two interacting atoms experience opposite forces due to
  !** the other in the pair.  
  Type Symmetric_Results
    Integer                      :: nspc,nderivs
    Logical                      :: subsetonly  
    Integer, Dimension(3)        :: subset
    Type(EnergyPlus)             :: total
    Logical, Dimension(:,:), Pointer                :: on   !* on/off
    Type(Symmetric_Shape), Dimension(:,:), Pointer  :: init !* init structure
    Type(Store_Level_Pair), Dimension(:,:), Pointer :: ab   !* spc-spc nrgs
    Type(Store_Level_Pair), Dimension(:,:), Pointer :: ba   !* spc-spc nrgs
  End Type Symmetric_Results

  !** The 'subsetonly' flag indicates that the storage structure contains
  !** only the interaction of a subset of the system with the whole system.
  !** In this case, the molecule number (and perhaps later atom numbers)
  !** of the first level in the level pairs may be mapped into the first
  !** index.  This keeps the storage size to a minimum, but requires that
  !** we keep track of which subset the storage represents.  The subset
  !** that is interation with the whole system is stored in 'subset'

Contains

  !----------------------------------------------------------------------------
  ! Initializes a symmetric results storage structure
  ! Requires:  storage -- structure to initialize
  !            nderivs -- number of derivatives in structure
  !            nspc -- number of species
  !----------------------------------------------------------------------------
  Subroutine storesym_init(storage,nspc,nderivs)
    Type(Symmetric_Results), Intent(Out)   :: storage
    Integer, Intent(In)                    :: nderivs,nspc

    Integer     :: i,j,error

    storage%nspc = nspc
    storage%nderivs = nderivs
    storage%subsetonly = .False.

    Call storebase_init(storage%total,nderivs)

    !** Allocate and initialize the on/off flags
    Allocate(storage%on(nspc,nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'storage%ncon')
    storage%on = .False.

    !** Allocate and initialize the on/off flags
    Allocate(storage%init(nspc,nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'storage%init')
    storage%init(:,:)%form = -1
    storage%init(:,:)%level = -1
    storage%init(:,:)%intra_also = .False.

    !** Allocate the storage arrays
    Allocate(storage%ab(nspc,nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'storage%ab')
    Allocate(storage%ba(nspc,nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'storage%ba')

    !** --Consider removing this to minimize memory leak possibilities--
    !** Terminate each allocated level pair for safety
    Do i = 1,nspc
      Do j = 1,nspc
        Call store_terminate(storage%ab(i,j),0)
        Call store_terminate(storage%ba(i,j),0)
      End Do
    End Do

  End Subroutine storesym_init

  !----------------------------------------------------------------------------
  ! Initializes and copies all or part of a symmetric results storage structure
  ! Requires:  old -- forcefield results structure to copy
  !            image -- new forcefield results structure (must be deallocated)
  !            subset -- integer array identifying subset to copy
  !----------------------------------------------------------------------------
  Subroutine storesym_initcopy(image,old,subset)
    Type(Symmetric_Results), Intent(Out)  :: image
    Type(Symmetric_Results), Intent(In)   :: old
    Integer, Dimension(3), Intent(In)     :: subset

    Integer     :: spc1,spc2,m,m1,m2,mol,index1,index2
    Integer     :: nspc,depth,nderivs,nmoles

    !** Check the depth of the copying 
    depth = storesym_depth(subset)
    If (depth > 2) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Not setup to discriminately copy on atom level'
      Stop
    End If

    nspc = old%nspc
    image%nspc = nspc

    !** Allocate the basics (allocates down to spc-spc level)
    Call storesym_init(image,nspc,old%nderivs)

    !** Copy in the easy stuff
    image%on = old%on
    image%nderivs = old%nderivs
    image%subsetonly = old%subsetonly
    image%subset = old%subset

    Call storebase_copy(image%total,old%total)

    !** Loop over the species pairs and copy desired sections
    Do spc1 = 1,image%nspc
      Do spc2 = spc1,image%nspc
        image%init(spc1,spc2) = old%init(spc1,spc2)
        If (.Not. image%on(spc1,spc2)) Cycle

        If (depth == 0) Then  !** simple, just copy the whole spc-spc 
          Call store_initcopy(image%ab(spc1,spc2),old%ab(spc1,spc2),.True.)
          Call store_initcopy(image%ba(spc1,spc2),old%ba(spc1,spc2),.True.)
          Cycle
        End If

        !** Skip if this spc-spc pair is not pertinent
        If ((spc1 /= subset(1)).And.(spc2 /= subset(1))) Then
          Call storesym_terminate(image,spc1,spc2,nderivs)
          image%on(spc1,spc2) = .False.
          Cycle
        End If

        nderivs = storebase_nderivs(old%ab(spc1,spc2)%total)

        !** Copy appropriate portion of AB and BA storage 
        If (subset(1) == spc1) Then

          !** Handle SPC detail storage simply and cycle
          If (old%init(spc1,spc2)%level == 2) Then
            If (nderivs > 0) Then
              Call store_initcopyidx(image%ab(spc1,spc2),old%ab(spc1,spc2), &
                  (/subset(2),subset(2)/),.True.)
              Call store_initcopy(image%ba(spc1,spc2),old%ba(spc1,spc2),.True.)
            Else
              Call store_initcopy(image%ab(spc1,spc2),old%ab(spc1,spc2),.True.)
              Call store_initcopy(image%ba(spc1,spc2),old%ba(spc1,spc2),.True.)
            End If
            Cycle
          End If

          !** partially initialize AB structure
          nmoles = old%init(spc1,spc2)%nmoles1
          Call store_initcopy(image%ab(spc1,spc2),old%ab(spc1,spc2),.False.)

          !** Do AB copying based on the similarity or dissimilarity of spc1,spc2
          If (spc1 == spc2) Then
            Call store_init(image%ab(spc1,spc2)%mi(1),(/nmoles/),01,nderivs)
            !** Copy relevant elements in the initialized upper triangle 
            Do m = 1,nmoles
              m1 = m
              m2 = subset(2)
              If (m1 > m2) Call swap(m1,m2)  !** copy from upper triangle only
              Call store_initcopy(image%ab(spc1,spc2)%mi(1)%im(m), &
                  old%ab(spc1,spc2)%mi(m1)%im(m2),.True.)
            End Do

          Else  
            !** Initialize and copy only specified molecule in 'ab'
            mol = subset(2)
            Call store_initcopy(image%ab(spc1,spc2)%mi(1), &
                old%ab(spc1,spc2)%mi(mol),.True.)
          End If
          
          !** skip BA form if not initialized
          If (old%init(spc1,spc2)%form == 0) Cycle

          !** partially initialize BA structure
          Call store_initcopy(image%ba(spc1,spc2),old%ba(spc1,spc2),.False.)

          !** Do BA copying based on the similarity or dissimilarity of spc1,spc2
          If (spc1 == spc2) Then

            !** Copy relevent elements in the initialized upper triangle 
            Do m = 1,nmoles
              m1 = m
              m2 = subset(2)
              If (m1 > m2) Call swap(m1,m2)  !** copy from upper triangle only
              Call store_init(image%ba(spc1,spc2)%im(m),(/1/),10,nderivs)  
              Call store_initcopy(image%ba(spc1,spc2)%im(m)%mi(1), &
                  old%ba(spc1,spc2)%im(m2)%mi(m1),.True.)
            End Do
          
          Else
            !** Initialize and copy only specified molecule in 'ba'
            If (store_idbranch(image%ba(spc1,spc2)) == 1) Then
              Do m = 1,old%ba(spc1,spc2)%sizes(1)
                Call store_init(image%ba(spc1,spc2)%im(m),(/1/),10,nderivs)  
                Call store_initcopy(image%ba(spc1,spc2)%im(m)%mi(1), &
                    old%ba(spc1,spc2)%im(m)%mi(mol),.True.)
              End Do
            End If
          End If

        Else  !** copy the whole spc-spc set
          Call store_initcopy(image%ab(spc1,spc2),old%ab(spc1,spc2),.True.)
          Call store_initcopy(image%ba(spc1,spc2),old%ba(spc1,spc2),.True.)
        End If
        
      End Do 
    End Do

  End Subroutine storesym_initcopy

  !------------------------------------------------------------------------- 
  ! Fill in the "gaps" in the subset storage.  The subset storage needs
  ! to be generic for any subset of the specified depth, but it was 
  ! copied from a specific subset.  This means that it has uninitialized
  ! self-interactions that need to be initialized as others are, "filled".
  ! Requires:  image -- forcefield results storage structure to check
  !            subset -- the subset for which the image structure is valid
  !------------------------------------------------------------------------- 
  Subroutine storesym_fillsub(storage,subset)
    Type(Symmetric_Results), Intent(InOut)  :: storage
    Integer, Dimension(3), Intent(In)       :: subset

    Integer      :: mol,spc,level,natoms,nderivs

    If ((subset(2) == 0).Or.(subset(3) /= 0)) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' have not thought about non-mol level subsets'
      Stop
    End If

    spc = subset(1)
    mol = subset(2)
    level = storage%init(spc,spc)%level 
    nderivs = storage%nderivs
    natoms = storage%init(spc,spc)%natoms1

    !** SPC detail storage is already ok, leave
    If (storage%init(spc,spc)%level == 2) Return

    !** Leave if the spc-spc pair is not inititialized
    If (storage%init(spc,spc)%form < 0) Return

    !** Clean up whatever is currently there
    Call store_clean(storage%ab(spc,spc)%mi(1)%im(mol))

    !** Initialize mol-mol entry at appropriate depth and nderivs
    !** remember that the mol-mol entry is now (1,mol)
    If (level > 3) Then
      Call store_init(storage%ab(spc,spc)%mi(1)%im(mol), &
          (/natoms,natoms/),11,nderivs)
    Else If ((level == 3).And.(nderivs > 0)) Then
      Call store_init(storage%ab(spc,spc)%mi(1)%im(mol), &
          (/natoms/),0,nderivs)
    Else
      Call store_terminate(storage%ab(spc,spc)%mi(1)%im(mol),nderivs)
    End If

    !** Do the same for the BA storage if initialized
    If (storage%init(spc,spc)%form > 0) Then

      !** Clean up whatever is currently there
      Call store_clean(storage%ba(spc,spc)%im(mol)%mi(1))

      !** Initialize mol-mol entry at appropriate depth and nderivs
      If (level > 3) Then
        Call store_init(storage%ba(spc,spc)%im(mol)%mi(1), &
            (/natoms,natoms/),11,nderivs)
      Else If ((level == 3).And.(nderivs > 0)) Then
        Call store_init(storage%ba(spc,spc)%im(mol)%mi(1), &
            (/natoms/),0,nderivs)
      Else
        Call store_terminate(storage%ba(spc,spc)%im(mol)%mi(1),nderivs)
      End If

    End If

  End Subroutine storesym_fillsub

  !----------------------------------------------------------------------------
  ! Copies one structure into another.  They need not have the same 
  ! initialization structure, but anything initialized in the "orig" structure
  ! will be assumed to be also initialized in the destination structure 
  ! Requires:  dest -- receiving symmetric storage structure (already init)
  !            orig -- forcefield results structure to copy
  !----------------------------------------------------------------------------
  Logical Function storesym_copy(dest,orig)
    Type(Symmetric_Results), Intent(InOut)  :: dest
    Type(Symmetric_Results), Intent(In)     :: orig

    Integer     :: i,spc1,spc2

    storesym_copy = .True.

    !** check the number of species
    If (orig%nspc /= dest%nspc) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' symmetric results structures do not having matching nspc'
      Stop      
    End If

    !** Loop over the species pairs and copy desired sections
    Do spc1 = 1,orig%nspc
      Do spc2 = spc1,orig%nspc
        If (.Not. orig%on(spc1,spc2)) Cycle
!        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!        Write(*,*) spc1,spc2

        !** Copy the summed top total
        Call storebase_copy(dest%ab(spc1,spc2)%total,orig%ab(spc1,spc2)%total)

        !** Copy the AB portion
        storesym_copy = store_copy(dest%ab(spc1,spc2),orig%ab(spc1,spc2))
#if DEBUG    
        If (.Not. storesym_copy) Then
          Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storesym_copy  !**DEBUG
          Write(*,*) 'DESTINATION'
          Call store_disp(dest%ab(spc1,spc2),.False.,2,6)
    
          Write(*,*) 'ORIGIN'
          Call store_disp(orig%ab(spc1,spc2),.False.,2,6)
        End If
#endif
        If (.Not. storesym_copy) Return

        !** skip the 'ba' part if it's not used
        If (orig%init(spc1,spc2)%form < 1) Cycle

        !** Copy the summed top total
        Call storebase_copy(dest%ba(spc1,spc2)%total,orig%ba(spc1,spc2)%total)

        !** Copy the BA portion    
        storesym_copy = store_copy(dest%ba(spc1,spc2),orig%ba(spc1,spc2))
!        Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storesym_copy  !**DEBUG
        If (.Not. storesym_copy) Return
        
      End Do 
    End Do

  End Function storesym_copy

  !----------------------------------------------------------------------------
  ! Copies one structure into another.  They need not have the same 
  ! initialization structure, but anything initialized in the "orig" structure
  ! will be assumed to be also initialized in the destination structure.  This
  ! is used mostly for copying interactions from temporary storage into the
  ! full forcefield storage structure.
  ! Requires:  dest -- receiving symmetric storage structure (already init)
  !            orig -- forcefield results structure to copy
  !            nmoleslist -- number of molecules of each species
  !            subset -- specification of which subset interaction is stored
  ! NOTE: This has only been written for temporary (orig) structures containing
  ! molecule-system interactions. It probably won't work for other cases.
  !----------------------------------------------------------------------------
  Logical Function storesym_update(dest,orig,nmoleslist,subset)
    Type(Symmetric_Results), Intent(InOut)  :: dest
    Type(Symmetric_Results), Intent(In)     :: orig
    Integer, Dimension(:), Intent(In)       :: nmoleslist
    Integer, Dimension(3), Intent(In)       :: subset

    Integer              :: i,mol1,mol2,spc1,spc2,depth
    Integer              :: om1,om2     !* origin mol1,mol2
    Integer              :: dm1,dm2     !* destination mol1,mol2
    Logical              :: samespc
    Real(kind=RDbl)      :: nrg1,nrg2

    !** Set default
    storesym_update = .True.

    !** Give warning if 'orig' is not a temporary structure
    If (.Not. orig%subsetonly) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' only designed for updating from a temporary structure'
      Stop      
    End If

    !** Check the number of species
    If (orig%nspc /= dest%nspc) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' symmetric results structures do not having matching nspc'
      Stop      
    End If

    !** Get the depth of the origin structure subset
    depth = storesym_depth(subset)

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Going to update full storage with subset: ',subset
    Write(*,*) 'BEFORE SUB-UPDATE: temp storesym storage:'
    Call storesym_display(orig,.False.,2,6)
!    Call storesym_display(orig,.True.,2,6)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'BEFORE SUB-UPDATE: Full storesym storage:'
!    Call storesym_display(dest,.False.,2,6)
    Call storesym_display(dest,.True.,2,6)
    Write(*,*) 'subset: ',subset
#endif

    !** Just copy the whole structure if origin has full system detail 
    If (depth == 0) Then 
      storesym_update = storesym_copy(dest,orig)
      Return
    End If

    !** Loop over the species pairs and copy desired sections
    Do spc1 = 1,orig%nspc
      Do spc2 = spc1,orig%nspc
        If (.Not. orig%on(spc1,spc2)) Cycle
        samespc = .False.
        If (spc1 == spc2) samespc = .True.

        !** Make sure the orig and dest structure have the same form
        If (dest%init(spc1,spc2)%form /= orig%init(spc1,spc2)%form) Then
          Write(0,'(2a,i6,a,i9)') __FILE__,":",__LINE__, &
              ' structures passed with different forms'
          Stop
        End If

        !** Skip spc-spc pairs not pertinent to the subset we're copying 
        If ((spc1 /= subset(1)) .And. (spc2 /= subset(1))) Cycle

#ifdef DEBUG
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) 'spc-spc pair: ',spc1,spc2
        Write(*,*) 'nmoles spc1: ',nmoleslist(spc1)
        Write(*,*) 'nmoles spc2: ',nmoleslist(spc2)
#endif

        !** Handle the easy possibility, subset only has SPC detail
        If (depth == 1) Then
          storesym_update = store_copy(dest%ab(spc1,spc2),orig%ab(spc1,spc2))
!          Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storesym_update  !**DEBUG
          If (.Not. storesym_update) Return
          If (orig%init(spc1,spc2)%form > 0) Then
            storesym_update = store_copy(dest%ba(spc1,spc2),orig%ba(spc1,spc2))
!          Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storesym_update  !**DEBUG
            If (.Not. storesym_update) Return
          End If
          Cycle
        End If

        !** Handle the copy/updating for MOL subsets in SPC-detail storage
        If (orig%init(spc1,spc2)%level == 2) Then
          If (orig%nderivs > 0) Then
            !** AB portion, only one molecule, copy it
            storesym_update = store_copy(dest%ab(spc1,spc2)%mi(subset(2)), &
                orig%ab(spc1,spc2)%mi(1))

            If (orig%init(spc1,spc2)%form > 0) Then
              !** BA portion, copy all molecules
              Do dm2 = 1,nmoleslist(spc2)
                storesym_update = store_copy(dest%ba(spc1,spc2)%im(dm2), &
                    orig%ba(spc1,spc2)%im(dm2))
                If (.Not. storesym_update) Return
              End Do
            End If
          Else
            Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
            Stop  !** HERE
            storesym_update = store_copy(dest%ab(spc1,spc2),orig%ab(spc1,spc2))
            If (orig%init(spc1,spc2)%form > 0) Then
              storesym_update = store_copy(dest%ba(spc1,spc2),orig%ba(spc1,spc2))
              If (.Not. storesym_update) Return
            End If
          End If
          Cycle
        End If

        !** Handle the copy/updating for MOL or ATM subsets depending on form
        Select Case(orig%init(spc1,spc2)%form)
        Case(0)  !** Deal with the easy case, map-type storage
          storesym_update = store_copy(dest%ab(spc1,spc2)%mi(subset(2)), &
              orig%ab(spc1,spc2)%mi(1))
!          Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storesym_update  !**DEBUG
          If (.Not. storesym_update) Return
          Cycle

        Case(1)  !** symmetric storage
          If (subset(1) == spc1) Then   !** subset species is spc1 or both are
            om1 = 1
            Do mol2 = 1,nmoleslist(spc2)
              om2 = mol2
              dm1 = subset(2)
              dm2 = mol2

              !** make sure it goes into upper triangle
              If ((samespc).And.(dm2 < dm1)) Call swap(dm1,dm2)
#ifdef DEBUG              
              nrg1 = orig%ab(spc1,spc2)%mi(om1)%im(om2)%total%nrg
              nrg2 = dest%ab(spc1,spc2)%mi(dm1)%im(dm2)%total%nrg
!              If ((Abs(nrg1) > 1.0e-10).Or.(Abs(nrg2) > 1.0e-10)) Then
                Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
                Write(*,*) spc1,om1,'  ',spc2,om2,' -> ',spc1,dm1,'  ',spc2,dm2
!                Write(*,*) spc1,spc2,'  ',om1,om2,' -> ',spc1,spc2,'  ',dm1,dm2
                Write(*,*) nrg1,' onto ',nrg2,'  Case 1'
!              End If
#endif
              storesym_update = store_copy(dest%ab(spc1,spc2)%mi(dm1)%im(dm2), &
                  orig%ab(spc1,spc2)%mi(om1)%im(om2))
              If (.Not. storesym_update) Return
              storesym_update = store_copy(dest%ba(spc1,spc2)%im(dm2)%mi(dm1), &
                  orig%ba(spc1,spc2)%im(om2)%mi(om1))
              If (.Not. storesym_update) Return
            End Do

          Else  !** second species matches subset species
            om2 = 1
            dm2 = subset(2)
            Do om1 = 1,nmoleslist(spc1)
              dm1 = om1

              If (samespc) Then  !** not expected
                Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
                Stop
              End If
#ifdef DEBUG              
              nrg1 = orig%ab(spc1,spc2)%mi(om1)%im(om2)%total%nrg
              nrg2 = dest%ab(spc1,spc2)%mi(dm1)%im(dm2)%total%nrg
!              If ((Abs(nrg1) > 1.0e-10).Or.(Abs(nrg2) > 1.0e-10)) Then
                Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
                Write(*,*) spc1,om1,'  ',spc2,om2,' -> ',spc1,dm1,'  ',spc2,dm2
!                Write(*,*) spc1,spc2,'  ',om1,om2,' -> ',spc1,spc2,'  ',dm1,dm2
                Write(*,*) nrg1,' onto ',nrg2,'  Case 2'
!              End If
#endif
              storesym_update = store_copy(dest%ab(spc1,spc2)%mi(dm1)%im(dm2), &
                  orig%ab(spc1,spc2)%mi(om1)%im(om2))
              If (.Not. storesym_update) Return
              storesym_update = store_copy(dest%ba(spc1,spc2)%im(dm2)%mi(dm1), &
                  orig%ba(spc1,spc2)%im(om2)%mi(om1))
              If (.Not. storesym_update) Return
            End Do

          End If

        Case Default
          Write(0,'(2a,i6,a,i9)') __FILE__,":",__LINE__, &
              ' Unable to recognize form number'
          Stop
        End Select
    
      End Do   !* spc2 loop
    End Do   !* spc1 loop

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'AFTER SUB-UPDATE: Full storesym storage:'
!    Call storesym_display(dest,.False.,2,6)
    Call storesym_display(dest,.True.,2,6)
#endif

  End Function storesym_update

  !----------------------------------------------------------------------------
  ! Initializes a species-species pair symmetric results storage structure
  ! Requires:  storage -- structure to initialize
  !            spc1 -- 1st index, number for species 1
  !            nmoles1 -- number of molecules to initialize 
  !            spc2 -- 2nd index, number for species 2
  !            nmoles2 -- number of molecules to initialize 
  !            init_type -- identifies the initialization method and structure
  !            level -- depth of initialization
  !            nderivs -- number of derivatives for the structure to hold
  !            intra_also -- also include intramolecular storage if spc1=spc2
  ! available initalization options:
  ! 0 => one-sided storage of type used in map routines
  ! 1 => two-sided, symmetric storage of type used in ssbasic routines
  !----------------------------------------------------------------------------
  Subroutine storesym_initsspair(storage,spc1,nmoles1,spc2,nmoles2, &
      init_type,level,nderivs,intra_also)
    Type(Symmetric_Results), Intent(InOut)   :: storage
    Integer, Intent(In)                      :: spc1,nmoles1,spc2,nmoles2
    Integer, Intent(In)                      :: init_type,level,nderivs
    Logical, Intent(In)                      :: intra_also

    Integer       :: a,m1,m2
    Integer       :: natoms1,natoms2

    !** Get number of atoms in the two species
    natoms1 = molecules_getnatoms(spc1)
    natoms2 = molecules_getnatoms(spc2)

    !** Store the initialization structure shape
    storage%init(spc1,spc2)%form = init_type
    storage%init(spc1,spc2)%level = level
    storage%init(spc1,spc2)%intra_also = intra_also
    storage%init(spc1,spc2)%natoms1 = natoms1
    storage%init(spc1,spc2)%natoms2 = natoms2
    storage%init(spc1,spc2)%nmoles1 = nmoles1
    storage%init(spc1,spc2)%nmoles2 = nmoles2
!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,*) spc1,spc2,init_type

    Select Case(init_type)
    Case (0)   !** map-type storage structure   (form 0)
      !** Note that the %im level is not initialized here because
      !** there is assumed to be only one molecule of species 2

      !** Don't need equal and opposite interactions, terminate 'ba'
      Call store_terminate(storage%ba(spc1,spc2),nderivs)

      !** Initialized based on detail level
      Select Case(level)
      Case(2,3)   !** SPC-detail and MOL-detail
        !** Initialize condensed atom storage if needed for forces
        If (nderivs > 0) Then
          Call store_init(storage%ab(spc1,spc2),(/nmoles1/),10,nderivs)
          Do m1 = 1,nmoles1
            Call store_init(storage%ab(spc1,spc2)%mi(m1),(/natoms1/),0,nderivs)
          End Do

        Else
          If (level == 3) Then
            Call store_init(storage%ab(spc1,spc2),(/nmoles1/),10,nderivs)
          Else
            Call store_terminate(storage%ab(spc1,spc2),nderivs)
          End If
        End If

      Case(4)   !** ATM-detail
        Call store_init(storage%ab(spc1,spc2),(/nmoles1/),10,nderivs)
        Do m1 = 1,nmoles1
          Call store_init(storage%ab(spc1,spc2)%mi(m1), &
              (/natoms1/),10,nderivs)
        End Do

      Case Default
        Write(0,'(2a,i4,a,i2)') __FILE__,": ",__LINE__, &
            ' Could not interpret storage level integer: ',level
        Stop
      End Select

    Case (1)   !** SYMMETRIC ssbasic-type storage structure (form 1)
      !** Initialized based on detail level
      Select Case(level)
      Case(2)   !** SPC-detail
        If (nderivs > 0) Then
          !** Initialize condensed atom storage for forces
          Call store_init(storage%ab(spc1,spc2),(/nmoles1/),10,nderivs)
          Do m1 = 1,nmoles1
            Call store_init(storage%ab(spc1,spc2)%mi(m1),(/natoms1/),0,nderivs)
          End Do
          Call store_init(storage%ba(spc1,spc2),(/nmoles2/),01,nderivs)
          Do m2 = 1,nmoles2
            Call store_init(storage%ba(spc1,spc2)%im(m2),(/natoms2/),0,nderivs)
          End Do
        Else
          Call store_terminate(storage%ab(spc1,spc2),nderivs)
          Call store_terminate(storage%ba(spc1,spc2),nderivs)
        End If

      Case(3)   !** MOL-detail
        !** 'ab' structure
        Call store_init(storage%ab(spc1,spc2),(/nmoles1/),10,nderivs)
        Do m1 = 1,nmoles1
          Call store_init(storage%ab(spc1,spc2)%mi(m1),(/nmoles2/),01,nderivs)
        End Do

        !** 'ba' structure
        Call store_init(storage%ba(spc1,spc2),(/nmoles2/),01,nderivs)
        Do m2 = 1,nmoles2
          Call store_init(storage%ba(spc1,spc2)%im(m2),(/nmoles1/),10,nderivs)
        End Do

        !** Initialize 'ab' atom-based storage
        If (nderivs > 0) Then
          Do m1 = 1,nmoles1
            Do m2 = 1,nmoles2
              If ((spc1 == spc2).And.(m2 <= m1)) Cycle
              If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle
              Call store_init(storage%ab(spc1,spc2)%mi(m1)%im(m2), &
                  (/natoms1/),0,nderivs)
            End Do
          End Do

          !** Initialize 'ba' atom-based storage
          Do m2 = 1,nmoles2
            Do m1 = 1,nmoles1
              If ((spc1 == spc2).And.(m2 <= m1)) Cycle
              If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle
              Call store_init(storage%ba(spc1,spc2)%im(m2)%mi(m1), &
                  (/natoms2/),0,nderivs)
            End Do
          End Do
        End If

      Case(4)   !** ATM-detail
        !** 'ab' structure
        Call store_init(storage%ab(spc1,spc2),(/nmoles1/),10,nderivs)
        Do m1 = 1,nmoles1
          Call store_init(storage%ab(spc1,spc2)%mi(m1),(/nmoles2/),01,nderivs)
        End Do

        !** 'ba' structure
        Call store_init(storage%ba(spc1,spc2),(/nmoles2/),01,nderivs)
        Do m2 = 1,nmoles2
          Call store_init(storage%ba(spc1,spc2)%im(m2),(/nmoles1/),10,nderivs)
        End Do

        !** Initialize 'ab' atom-based storage
        Do m1 = 1,nmoles1
          Do m2 = 1,nmoles2
            If ((spc1 == spc2).And.(m2 <= m1)) Cycle
            If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle
            Call store_init(storage%ab(spc1,spc2)%mi(m1)%im(m2), &
                (/natoms1,natoms2/),11,nderivs)
          End Do
        End Do

        !** Initialize 'ba' atom-based storage
        Do m2 = 1,nmoles2
          Do m1 = 1,nmoles1
            If ((spc1 == spc2).And.(m2 <= m1)) Cycle
            If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle
            Call store_init(storage%ba(spc1,spc2)%im(m2)%mi(m1), &
                (/natoms1,natoms2/),11,nderivs)
          End Do
        End Do

      Case Default
        Write(0,'(2a,i4,a,i2)') __FILE__,": ",__LINE__, &
            ' Could not interpret storage level integer: ',level
        Stop
      End Select

    End Select  !** (select based on form)

    Return

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call store_disp(storage%ab(spc1,spc2),.False.,2,6)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT    
    Call store_disp(storage%ba(spc1,spc2),.False.,2,6)
    Stop

  End Subroutine storesym_initsspair

  !----------------------------------------------------------------------------
  ! Get a pointer to the EnergyPlus structure that stores a given interaction.
  ! This routine is controlled by a pair of paths.  
  ! Requires:  storage -- storage level pair structure
  !            path1 -- level 1 path
  !            path2 -- level 2 path
  !            abptr -- returned pointer for A interactions from B
  !            baptr -- returned pointer for B interactions from A
  !----------------------------------------------------------------------------
  Logical Function storesym_basicptrs(storage,path1,path2, &
      abptr,baptr) Result(res)
    Type(Symmetric_Results), Intent(In)  :: storage
    Integer, Dimension(:), Intent(In)    :: path1,path2
    Type(EnergyPlus), Pointer            :: abptr
    Type(EnergyPlus), Pointer, Optional  :: baptr

    !** Set default return value
    res = .False.

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(0,'(a)') 'storesym_basicptrs is unfinished'
    Stop

  End Function storesym_basicptrs

  !----------------------------------------------------------------------------
  ! Zeros the storage structure using calls to store_zero.  This call path is
  ! more general than storesym_fastzero, but will be slower.
  ! Requires:  storage -- symmetric results structure
  !----------------------------------------------------------------------------
  Subroutine storesym_zero(storage)
    Type(Symmetric_Results), Intent(InOut)  :: storage

    Integer              :: i,j

    Call storebase_zero(storage%total)

    Do i = 1,storage%nspc
      Do j = i,storage%nspc
        If (storage%on(i,j)) Then
          Call store_zero(storage%ab(i,j),.True.)
          If (storage%init(i,j)%form < 1) Cycle
          Call store_zero(storage%ba(i,j),.True.)
        End If
      End Do
    End Do

  End Subroutine storesym_zero

  !----------------------------------------------------------------------------
  ! Zero a symmetric results storage structure.  Meant to be fast, that is 
  ! why it directly accesses data types from storebase and store.  
  ! ASSUMES that the forces are associated.
  ! Requires:  storage -- symmetric results structure
  !----------------------------------------------------------------------------
  Subroutine storesym_fastzero(storage)
    Type(Symmetric_Results), Intent(InOut)  :: storage

    Integer              :: i,j,a1,a2,m1,m2
    Integer              :: spc1,spc2
    Integer              :: natoms1,natoms2
    Integer              :: nmoles1,nmoles2
    Logical              :: intra_also

    Call storebase_zero(storage%total)

    !** zero the upper triangle of the storage matrix
    Do spc1 = 1,storage%nspc
      Do spc2 = spc1,storage%nspc

        !** skip if this storage not on
        If (.Not. storage%on(spc1,spc2)) Cycle

        natoms1 = storage%init(spc1,spc2)%natoms1
        natoms2 = storage%init(spc1,spc2)%natoms2
        nmoles1 = storage%init(spc1,spc2)%nmoles1
        nmoles2 = storage%init(spc1,spc2)%nmoles2
        intra_also = storage%init(spc1,spc2)%intra_also

#ifdef DEBUG
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) 'storesym_fastzero: ',storage%on(spc1,spc2),spc1,spc2
        Write(*,*) 'form:  ',storage%init(spc1,spc2)%form
        Write(*,*) 'level: ',storage%init(spc1,spc2)%level
!        stop
#endif

        If (storage%init(spc1,spc2)%level == 4) Then  !** ATM detail
          
          Select Case(storage%init(spc1,spc2)%form)
          Case (0)   !** map-type storage structure
            storage%ab(spc1,spc2)%total%nrg = 0.0_RDbl
            storage%ab(spc1,spc2)%total%force = 0.0_RDbl
  
            Do m1 = 1,nmoles1
              storage%ab(spc1,spc2)%mi(m1)%total%nrg = 0.0_RDbl
              storage%ab(spc1,spc2)%mi(m1)%total%force = 0.0_RDbl
              Do a1 = 1,natoms1
                storage%ab(spc1,spc2)%mi(m1)%mi(a1)%total%nrg = 0.0_RDbl
                storage%ab(spc1,spc2)%mi(m1)%mi(a1)%total%force = 0.0_RDbl
              End Do
            End Do
  
          Case (1)   !** symmetric ssbasic-type storage structure
            intra_also = storage%init(spc1,spc2)%intra_also
            storage%ab(spc1,spc2)%total%nrg = 0.0_RDbl
            storage%ab(spc1,spc2)%total%force = 0.0_RDbl
  
            storage%ba(spc1,spc2)%total%nrg = 0.0_RDbl
            storage%ba(spc1,spc2)%total%force = 0.0_RDbl
  
            Do m1 = 1,nmoles1
              Do m2 = 1,nmoles2

                If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle

                If ((spc1 /= spc2).Or.(m1 <= m2)) Then   !** zero ab storage
                  storage%ab(spc1,spc2)%mi(m1)%im(m2)%total%nrg = 0.0_RDbl
                  storage%ab(spc1,spc2)%mi(m1)%im(m2)%total%force = 0.0_RDbl
                          
                  Do a1 = 1,natoms1
                    Do a2 = 1,natoms2
                      storage%ab(spc1,spc2)%mi(m1)%im(m2)%mm(a1,a2)%total%nrg = &
                          0.0_RDbl
                      storage%ab(spc1,spc2)%mi(m1)%im(m2)%mm(a1,a2)%total%force =&
                          0.0_RDbl
                    End Do
                  End Do
                End If
  
                If ((spc1 /= spc2).Or.(m1 <= m2)) Then   !** zero ba storage
                  storage%ba(spc1,spc2)%im(m2)%mi(m1)%total%nrg = 0.0_RDbl
                  storage%ba(spc1,spc2)%im(m2)%mi(m1)%total%force = 0.0_RDbl
  
                  Do a1 = 1,natoms1
                    Do a2 = 1,natoms2
                      storage%ba(spc1,spc2)%im(m2)%mi(m1)%mm(a1,a2)%total%nrg = &
                          0.0_RDbl
                      storage%ba(spc1,spc2)%im(m2)%mi(m1)%mm(a1,a2)%total%force =&
                          0.0_RDbl
                    End Do
                  End Do
                End If
  
              End Do
            End Do
  
          End Select

        Else If (storage%init(spc1,spc2)%level == 3) Then   !** MOL detail

          Select Case(storage%init(spc1,spc2)%form)
          Case (0)   !** map-type storage structure   
            storage%ab(spc1,spc2)%total%nrg = 0.0_RDbl
            storage%ab(spc1,spc2)%total%force = 0.0_RDbl
  
            Do m1 = 1,nmoles1
              storage%ab(spc1,spc2)%mi(m1)%total%nrg = 0.0_RDbl
              storage%ab(spc1,spc2)%mi(m1)%total%force = 0.0_RDbl
              Do a1 = 1,natoms1
                storage%ab(spc1,spc2)%mi(m1)%binfo(a1)%nrg = 0.0_RDbl
                storage%ab(spc1,spc2)%mi(m1)%binfo(a1)%force = 0.0_RDbl
              End Do
            End Do
  
          Case (1)   !** symmetric ssbasic-type storage structure
            intra_also = storage%init(spc1,spc2)%intra_also
            storage%ab(spc1,spc2)%total%nrg = 0.0_RDbl
            storage%ab(spc1,spc2)%total%force = 0.0_RDbl
  
            storage%ba(spc1,spc2)%total%nrg = 0.0_RDbl
            storage%ba(spc1,spc2)%total%force = 0.0_RDbl
  
            Do m1 = 1,nmoles1
              Do m2 = 1,nmoles2
  
                If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle

                If ((spc1 /= spc2).Or.(m1 <= m2)) Then   !** zero ab storage
                  storage%ab(spc1,spc2)%mi(m1)%im(m2)%total%nrg = 0.0_RDbl
                  storage%ab(spc1,spc2)%mi(m1)%im(m2)%total%force = 0.0_RDbl
                        
                  Do a1 = 1,natoms1
                    storage%ab(spc1,spc2)%mi(m1)%im(m2)%binfo(a1)%nrg = 0.0_RDbl
                    storage%ab(spc1,spc2)%mi(m1)%im(m2)%binfo(a1)%force = 0.0_RDbl
                  End Do
                End If

                If ((spc1 /= spc2).Or.(m1 <= m2)) Then   !** zero ba storage
                  storage%ba(spc1,spc2)%im(m2)%mi(m1)%total%nrg = 0.0_RDbl
                  storage%ba(spc1,spc2)%im(m2)%mi(m1)%total%force = 0.0_RDbl
                
                  Do a2 = 1,natoms2
                    storage%ba(spc1,spc2)%im(m2)%mi(m1)%binfo(a2)%nrg = 0.0_RDbl
                    storage%ba(spc1,spc2)%im(m2)%mi(m1)%binfo(a2)%force = 0.0_RDbl
                  End Do
                End If
  
              End Do
            End Do
  
          End Select

        Else If (storage%init(spc1,spc2)%level == 2) Then   !** SPC detail
          !** Zero 'ab' portion
          storage%ab(spc1,spc2)%total%nrg = 0.0_RDbl
          storage%ab(spc1,spc2)%total%force = 0.0_RDbl
          Do m1 = 1,storage%ab(spc1,spc2)%sizes(1)
            storage%ab(spc1,spc2)%mi(m1)%total%nrg = 0.0_RDbl
            Do a1 = 1,storage%ab(spc1,spc2)%mi(m1)%sizes(1)
              storage%ab(spc1,spc2)%mi(m1)%binfo(a1)%nrg = 0.0_RDbl
              storage%ab(spc1,spc2)%mi(m1)%binfo(a1)%force = 0.0_RDbl
            End Do
          End Do

          !** Zero 'ba' portion if it exists
          If (storage%init(spc1,spc2)%form > 0) Then
            storage%ba(spc1,spc2)%total%nrg = 0.0_RDbl
            storage%ba(spc1,spc2)%total%force = 0.0_RDbl
            Do m2 = 1,storage%ba(spc1,spc2)%sizes(1)
              storage%ba(spc1,spc2)%im(m2)%total%nrg = 0.0_RDbl
              Do a2 = 1,storage%ba(spc1,spc2)%im(m2)%sizes(1)
                storage%ba(spc1,spc2)%im(m2)%binfo(a2)%nrg = 0.0_RDbl
                storage%ba(spc1,spc2)%im(m2)%binfo(a2)%force = 0.0_RDbl
              End Do
            End Do
          End If

        Else 
          Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              ' storesym_fastzero not prepared to handle this type of storage'
          Write(0,'(a,i2)') 'form:  ',storage%init(spc1,spc2)%form
          Write(0,'(a,i2)') 'level: ',storage%init(spc1,spc2)%level
          Stop

        End If

      End Do
    End Do

  End Subroutine storesym_fastzero

  !----------------------------------------------------------------------------
  ! Subtract the second base structure from the first.  Assumes that anything
  ! initialized in the 2nd structure is initialized in the 1st.
  ! Requires:  storage1 -- 1st symmetric results structure
  !            storage2 -- 2nd symmetric results structure
  !----------------------------------------------------------------------------
  Subroutine storesym_subtract(storage1,storage2)
    Type(Symmetric_Results), Intent(InOut)  :: storage1
    Type(Symmetric_Results), Intent(In)     :: storage2

    Integer              :: i,j

    !** Handle the total quantity
    Call storebase_subtract(storage1%total,storage2%total)

    Do i = 1,storage1%nspc
      Do j = i,storage1%nspc
        If (storage2%on(i,j)) Then
          Call store_subtract(storage1%ab(i,j),storage2%ab(i,j),.True.)
          If (storage2%init(i,j)%form < 1) Cycle
          Call store_subtract(storage1%ba(i,j),storage2%ba(i,j),.True.)
        End If
      End Do
    End Do

  End Subroutine storesym_subtract

  !----------------------------------------------------------------------------
  ! Add the second base structure to the first.  Assumes that the
  ! initialization structures are identical.
  ! Requires:  storage1 -- 1st symmetric results structure
  !            storage2 -- 2nd symmetric results structure
  !----------------------------------------------------------------------------
  Subroutine storesym_add(storage1,storage2)
    Type(Symmetric_Results), Intent(InOut)  :: storage1
    Type(Symmetric_Results), Intent(In)     :: storage2

    Integer              :: i,j

    !** Handle the total quantity
    Call storebase_add(storage1%total,storage2%total)

    Do i = 1,storage1%nspc
      Do j = i,storage1%nspc
        If (storage1%on(i,j)) Then
          Call store_add(storage1%ab(i,j),storage2%ab(i,j),.True.)
          If (storage1%init(i,j)%form < 1) Cycle
          Call store_add(storage1%ba(i,j),storage2%ba(i,j),.True.)
        End If
      End Do
    End Do

  End Subroutine storesym_add

  !----------------------------------------------------------------------------
  ! Sum the symmetric storage recursively
  ! Requires:  storage -- symmetric results structure
  !----------------------------------------------------------------------------
  Subroutine storesym_sum(storage)
    Type(Symmetric_Results), Intent(InOut)  :: storage

    Integer              :: i,j

    Call storebase_zero(storage%total)

    Do i = 1,storage%nspc
      Do j = i,storage%nspc
        If (storage%on(i,j)) Then
          Call store_sum(storage%ab(i,j),.True.)
          Call storebase_inc(storage%total,storage%ab(i,j)%total)
          If (storage%init(i,j)%form < 1) Cycle
          Call store_sum(storage%ba(i,j),.True.)
          Call storebase_inc(storage%total,storage%ba(i,j)%total)
        End If
      End Do
    End Do

  End Subroutine storesym_sum

  !----------------------------------------------------------------------------
  ! An attempt to sum faster.  Assumes that the mol-mol level sums are already 
  ! correct.  The intermediate summations are not updated here, only spc-spc.
  ! It also ONLY SUMS THE ENERGIES.
  ! Requires:  storage -- symmetric results structure
  !----------------------------------------------------------------------------
  Subroutine storesym_fastsum(storage)
    Type(Symmetric_Results), Intent(InOut)  :: storage

    Integer              :: i,j,a1,a2,m1,m2
    Integer              :: spc1,spc2
    Integer              :: natoms1,natoms2
    Integer              :: nmoles1,nmoles2
    Logical              :: intra_also

    !** Zero the total storage quantity
    Call storebase_zero(storage%total)

    !** Sum the upper triangle of the storage matrix
    Do spc1 = 1,storage%nspc
      Do spc2 = spc1,storage%nspc

        !** skip if this storage not on
        If (.Not. storage%on(spc1,spc2)) Cycle

        !** Zero the spc-spc level total or use it if it's the extent of detail
        If ((storage%init(spc1,spc2)%level == 2).And.(storage%nderivs == 0)) Then
          Call storebase_inc(storage%total,storage%ab(spc1,spc2)%total)
          If (storage%init(spc1,spc2)%form >= 1) Then
            Call storebase_inc(storage%total,storage%ba(spc1,spc2)%total)
          End If
          Cycle
        Else
          storage%ab(spc1,spc2)%total%nrg = 0.0_RDbl
        End If

        natoms1 = storage%init(spc1,spc2)%natoms1
        natoms2 = storage%init(spc1,spc2)%natoms2
        nmoles1 = storage%init(spc1,spc2)%nmoles1
        nmoles2 = storage%init(spc1,spc2)%nmoles2
        intra_also = storage%init(spc1,spc2)%intra_also

#ifdef DEBUG
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) 'on,spc1,spc2: ',storage%on(spc1,spc2),spc1,spc2
        Write(*,*) storage%init(spc1,spc2)%form
        Write(*,*) storage%init(spc1,spc2)%level
#endif

        If (storage%init(spc1,spc2)%level == 4) Then     !** ATM detail

          Select Case(storage%init(spc1,spc2)%form)
          Case (0)   !** map-type storage structure
            Do m1 = 1,nmoles1
              storage%ab(spc1,spc2)%total%nrg =  &
                  storage%ab(spc1,spc2)%total%nrg + &
                  storage%ab(spc1,spc2)%mi(m1)%total%nrg 
            End Do
  
          Case (1)   !** symmetric ssbasic-type storage structure
            Do m1 = 1,nmoles1
              Do m2 = 1,nmoles2
                If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle
  
                !** ab storage
                storage%ab(spc1,spc2)%total%nrg = &
                    storage%ab(spc1,spc2)%total%nrg + &
                    storage%ab(spc1,spc2)%mi(m1)%im(m2)%total%nrg

                If (storage%init(spc1,spc2)%form < 1) Cycle
                
                !** ba storage
                storage%ba(spc1,spc2)%total%nrg = &
                    storage%ba(spc1,spc2)%total%nrg + &
                    storage%ba(spc1,spc2)%im(m2)%mi(m1)%total%nrg 
              End Do
            End Do
  
          End Select

        Else If (storage%init(spc1,spc2)%level == 3) Then     !** MOL detail
          !** note, currently this is identical to ATM detail procedure
          !** because the summing is assumed to be accurate at the MOL level

          Select Case(storage%init(spc1,spc2)%form)
          Case (0)   !** map-type storage structure
!            Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!            Write(*,*) spc1,spc2,storage%ab(spc1,spc2)%total%nrg 
            Do m1 = 1,nmoles1
              storage%ab(spc1,spc2)%total%nrg =  &
                  storage%ab(spc1,spc2)%total%nrg + &
                  storage%ab(spc1,spc2)%mi(m1)%total%nrg 
            End Do
  
          Case (1)   !** symmetric ssbasic-type storage structure
            Do m1 = 1,nmoles1
              Do m2 = 1,nmoles2

                If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle
  
                !** AB storage
                storage%ab(spc1,spc2)%total%nrg = &
                    storage%ab(spc1,spc2)%total%nrg + &
                    storage%ab(spc1,spc2)%mi(m1)%im(m2)%total%nrg

                If (storage%init(spc1,spc2)%form < 1) Cycle
                
                !** BA storage
                storage%ba(spc1,spc2)%total%nrg = &
                    storage%ba(spc1,spc2)%total%nrg + &
                    storage%ba(spc1,spc2)%im(m2)%mi(m1)%total%nrg 
              End Do
            End Do
  
          End Select

        Else If (storage%init(spc1,spc2)%level == 2) Then     !** SPC detail
          !** Unless there are stored derivatives, there is nothing to do
          If (storage%nderivs > 0) Then
            Select Case(storage%init(spc1,spc2)%form)
            Case (0)   !** map-type storage structure
!              Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!              Write(*,*) spc1,spc2,storage%ab(spc1,spc2)%total%nrg 
              Do m1 = 1,nmoles1
                storage%ab(spc1,spc2)%total%nrg =  &
                    storage%ab(spc1,spc2)%total%nrg + &
                    storage%ab(spc1,spc2)%mi(m1)%total%nrg 
              End Do
  
            Case (1)   !** symmetric ssbasic-type storage structure
              Do m1 = 1,nmoles1
                storage%ab(spc1,spc2)%total%nrg =  &
                    storage%ab(spc1,spc2)%total%nrg + &
                    storage%ab(spc1,spc2)%mi(m1)%total%nrg 
              End Do
              Do m2 = 1,nmoles2
                storage%ba(spc1,spc2)%total%nrg =  &
                    storage%ba(spc1,spc2)%total%nrg + &
                    storage%ba(spc1,spc2)%im(m2)%total%nrg 
              End Do
            End Select
          End If
          
        Else
          Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              ' storesym_fastsum not prepared to handle this type of storage'
          Write(0,'(a,i2)') 'form:  ',storage%init(spc1,spc2)%form
          Write(0,'(a,i2)') 'level: ',storage%init(spc1,spc2)%level
          Stop

        End If

        !** Increment the total storage with the spc-spc value
        Call storebase_inc(storage%total,storage%ab(spc1,spc2)%total)
        If (storage%init(spc1,spc2)%form >= 1) Then
          Call storebase_inc(storage%total,storage%ba(spc1,spc2)%total)
        End If

      End Do
    End Do

  End Subroutine storesym_fastsum

  !----------------------------------------------------------------------------
  ! Terminate an entry pair in a symmetric storage structure. Does not touch
  ! the total field because it applies to all of the pairs.
  ! Requires:  storage -- array of storage level pair structures
  !            spc1 -- 1st index, number for species 1
  !            spc2 -- 2nd index, number for species 2
  !            nderivs -- number of derivatives to store (0,1,2)
  !            intra -- flags initialization of intra potentials, default False
  !----------------------------------------------------------------------------
  Subroutine storesym_terminate(storage,spc1,spc2,nderivs,intra)
    Type(Symmetric_Results), Intent(InOut)  :: storage
    Integer, Intent(In)                     :: spc1,spc2,nderivs
    Logical, Intent(In), Optional           :: intra

    !** initialize and zero total storage and sub-storage
    If (Present(intra)) Then
      Call store_terminate(storage%ab(spc1,spc2),nderivs,intra)
      Call store_terminate(storage%ba(spc1,spc2),nderivs,intra)
    Else
      Call store_terminate(storage%ab(spc1,spc2),nderivs)
      Call store_terminate(storage%ba(spc1,spc2),nderivs)
    End If

  End Subroutine storesym_terminate

  !----------------------------------------------------------------------------
  ! Returns a pair of pointers to level-pair structures at the desired 
  ! species-molecule-atom level.  The returned pointers are "symmetric", one 
  ! contains interactions of A on B and on contains interactions for B on A.  
  ! This is needed to store interactions such as forces on two moving groups.
  ! May also only return a SINGLE pointer to a level-pair structures at the 
  ! desired species-molecule-atom level.  Use this when the second group does 
  ! not need to have its own symmetric interactions stored, for example, when
  ! it will not be moved.
  ! Requires:  storage -- structure to extract pointer(s) from
  !            sub1 -- 1st species-molecule-atom subset 
  !            sub2 -- 2nd species-molecule-atom subset 
  !            branch -- indicates the initialization type of the pointers
  !            indices -- correct array indices in abptr, baptr
  !            nrgptr -- pointer to an EnergyPlus plus
  !            abptr -- pointer containing interactions on subset1 from subset2
  !            baptr -- pointer containing interactions on subset1 from subset2
  ! a "subset" is a 1D array with indices: species,molecule,atom number
  ! NOTE: unfortunately, Fortran does not allow pointers to be assigned to
  ! anything but other pointers and "targets", so I cannot point to single
  ! level pairs but must point to arrays of level pairs and return the desired
  ! index or indices in the array(s).  We should fix it if standards change.
  ! Where's real polymorphism when you need it??
  ! NOTE2: For spc1=spc2 cases, the convention is that mol1<=mol2, 
  ! ssbasic is currently setup so that this is not a
  ! problem.  But, if the setup of ssbasic were to change, it would be 
  ! necessary to check this criteria in this routine.
  !----------------------------------------------------------------------------
  Subroutine storesym_ptrs(storage,sub1,sub2,branch,indices,nrgptr,abptr,baptr)
    Type(Symmetric_Results), Intent(InOut)                  :: storage
    Integer, Dimension(3), Intent(In)                       :: sub1,sub2
    Integer, Intent(Out)                                    :: branch
    Integer, Dimension(2), Intent(Out)                      :: indices
    Type(EnergyPlus), Pointer                               :: nrgptr
    Type(Store_Level_Pair), Dimension(:), Pointer           :: abptr
    Type(Store_Level_Pair), Dimension(:), Pointer, Optional :: baptr

    Integer                 :: depth1,depth2,unit
    Integer                 :: spc1,spc2,mol1,mol2
    Logical                 :: found
    Integer, Dimension(3)   :: subset1,subset2

    !** Copy the subsets so we can change them
    subset1 = sub1
    subset2 = sub2

    !** If the storage structure is only subset-system, map the molecule 
    !** number of the appropriate subset to the first entry.
    If (storage%subsetonly) Then
      If (storage%subset(1) == subset1(1)) Then
        subset1(2) = 1  
      Else
        subset2(2) = 1  
      End If
    End If

    !** Set the internal subsets such that spc1 <= spc2
    If (subset1(1) > subset2(1)) Call swap(subset1,subset2)

    !** Get the depths of the two subsets
    depth1 = storesym_depth(subset1)
    depth2 = storesym_depth(subset2)

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,'(a,3i2,3x,i2)') 'subset1: ',subset1,depth1
!    Write(*,'(a,3i2,3x,i2)') 'subset2: ',subset2,depth2

    !** Set the species numbers and zero the returned indices
    spc1 = subset1(1)
    spc2 = subset2(1)
    indices = 0

    !** Find the correct level pair pointer(s)
    found = .False.
    branch = -2
    Select Case(depth1)
    Case (1)    !** SPC-?
      Select Case(depth2)
      Case (1)    !** SPC-SPC
      Case (2)    !** SPC-MOL
      End Select

    Case (2)    !** MOL-?
      mol1 = subset1(2)

      Select Case(depth2)
      Case (1)    !** MOL-SPC
        found = .True.
        If (storage%init(spc1,spc2)%level == 2) Then   !** SPC detail case
          If (storage%nderivs == 0) Then
            branch = -2
          Else
            abptr => storage%ab(spc1,spc2)%mi
            indices(1) = mol1
            branch = store_idbranch(abptr(indices(1)))
            If (Present(baptr)) Then
              Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
              Stop
            End If
          End If

        Else
!          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!          Write(*,*) 'assigned to storage%ab(spc1,spc2)%mi ',spc1,spc2
          abptr => storage%ab(spc1,spc2)%mi
          indices(1) = mol1
          branch = store_idbranch(abptr(indices(1)))
          If (Present(baptr)) Then
            Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
            Stop
          End If
        End If

      Case (2)    !** MOL-MOL
        found = .True.
        mol2 = subset2(2)
        If (storage%init(spc1,spc2)%level == 2) Then   !** SPC detail case
          If (storage%nderivs == 0) Then
            nrgptr => storage%ab(spc1,spc2)%total
            branch = -2
          Else
            !** If it's same species, apply m1<m2 convention 
            If ((spc1 == spc2).And.(mol1 > mol2)) Call swap(mol1,mol2)

!            Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!            Write(*,*) 'getting pointers for: ',spc1,mol1,'     ',spc2,mol2
            abptr => storage%ab(spc1,spc2)%mi
            indices(1) = mol1
            If (Present(baptr)) baptr => storage%ba(spc1,spc2)%im
            indices(2) = mol2
            branch = store_idbranch(abptr(indices(1)))
!            Write(*,*) 'branch is: ',branch
!            Write(*,*) 'indices are: ',indices
          End If 

        Else
!          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!          Write(*,*) 'getting pointers for: ',spc1,mol1,'     ',spc2,mol2
          abptr => storage%ab(spc1,spc2)%mi(mol1)%im
          indices(1) = mol2
          branch = store_idbranch(abptr(indices(1)))
  
          If (Present(baptr)) baptr => storage%ba(spc1,spc2)%im(mol2)%mi
          indices(2) = mol1
!          Write(*,*) 'indices are: ',indices
        End If

!        Call storesym_display(storage,.False.,2,6)
#ifdef DEBUG
        If (spc1 /= spc2) Then
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Write(*,*) 'dumping forcefield structure'
          unit = file_open('ffinside.txt')
          Call storesym_display(storage,.False.,2,unit)
        End If
        Close(unit=unit)
#endif
#ifdef DEBUG
        Write(*,*) spc1,spc2,mol1,mol2
        Write(*,*) 'ab ', &
            store_idbranch(storage%ab(spc1,spc2)%mi(mol1)%im(mol2)), &
            storage%ab(spc1,spc2)%mi(mol1)%im(mol2)%sizes(1)
        Write(*,*) spc1,spc2,mol2,mol1
        Write(*,*) 'ba ', &
            store_idbranch(storage%ba(spc1,spc2)%im(mol2)%mi(mol1)), &
            storage%ba(spc1,spc2)%im(mol2)%mi(mol1)%sizes(1)
!        Stop
#endif
      End Select
    End Select

    If (.Not. found) Then
      Write(0,'(2a,i4,a,2i2)') __FILE__," : ",__LINE__, &
          ' storesym_ptrs not equipped for this case: ',depth1,depth2
      Stop
    End If

  End Subroutine storesym_ptrs

  !----------------------------------------------------------------------------
  ! The symmetric counterpart to storetop_extract.  Gets a summed interaction
  ! quantity between two different system subsets.  Assumes that mol-mol
  ! summed quantities are accurate!
  ! Requires:  ffresults -- forcefield results structure
  !            info -- returned energy plus possible derivatives structure
  !            subset1 -- 1st path specifier: spc,mol,atm
  !            subset2 -- 2nd path specifier: spc,mol,atm
  ! Usage: the spc,molec or atm designators can also be set to zero.  This
  !        will tell the routine not descend to that level.  For example,
  !        spc1=1, molec1=0, spc2=2, molec2=0 => get overall interaction of 
  !        species 1 with species 2.
  ! Needed Improvements:
  ! 1) assumes that 2nd spc of map form is the map species, restrictive
  ! 2) only subset1 with full system extraction possible
  !----------------------------------------------------------------------------
  Logical Function storesym_extract(storage,info,subset1,subset2)
    Type(Symmetric_Results), Intent(In)      :: storage
    Type(EnergyPlus), Intent(InOut)          :: info
    Integer, Dimension(3), Intent(In)        :: subset1,subset2

    Integer                      :: spc1,spc2,s1,s2,mol,a,match,depth
    Character(len=strLen)        :: string1
    Integer, Dimension(4,2)      :: path1,path2

    storesym_extract = .True.

    !** Copy path information
    spc1 = subset1(1)
    mol = subset1(2)
    spc2 = subset2(1)
    path1(:,1) = (/subset1,0/)
    path1(:,2) = (/subset1,0/)
    path2(:,1) = (/subset2,0/)
    path2(:,2) = (/subset2,0/)

    !** Handle simple system-system case
    If ((spc1 == 0) .And. (spc2 == 0)) Then
      Call storebase_copy(info,storage%total)
      Return
    End If

#ifdef DEBUG
    If (dbgflag) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Call storesym_display(storage,.False.,2,6)
    End If
#endif

    !** Find the interaction info
    If (spc2 == 0) Then    !** subset1 with full system

      Do s1 = 1,storage%nspc
        Do s2 = s1,storage%nspc
          !** Skip out if this spc-spc pair is not pertinent

          If ((s1 /= spc1).And.(s2 /= spc1)) Cycle

          If (.Not. storage%on(s1,s2)) Cycle

          !** Identify matching species
          match = 0
          If (s1 == spc1) match = 1
          If (s2 == spc1) match = match + 2

          !** Restore paths that may have been changed
          path1(:,1) = (/subset1,0/)
          path1(:,2) = (/subset1,0/)
          path2(:,1) = (/subset2,0/)
          path2(:,2) = (/subset2,0/)

          !** Handle the simple, species-only case, not form dependent
          !** assuming that structure is accurately summed here
          If (mol == 0) Then

            Call storebase_inc(info,storage%ab(s1,s2)%total)
            Cycle
          End If

#ifdef DEBUG         
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Write(*,*) 'SPECIES PAIR: ',s1,s2
          Write(*,*) 'mol: ',mol
          Write(*,*) 'form, level: ',storage%init(s1,s2)%form, &
              storage%init(s1,s2)%level
          Write(*,*) 'nmoles1, nmoles2: ',storage%init(s1,s2)%nmoles1, &
              storage%init(s1,s2)%nmoles2
          Write(*,*) 'current energy (before extraction): ',info%nrg
#endif

          !** Extract based on the form of the structure
          Select Case(storage%init(s1,s2)%form)
          Case (0)  !** one-sided form used for maps

            !     storesym_extract = store_extract(storage%ab(s1,s2), &
            !     path1(2:,:),path2(2:,:),info)   !* should be enough, try it

            If (storage%init(s1,s2)%level == 3) Then  !** MOL detail
              !** is this NECESSARY? mol-level sums should be ok! try above
              storesym_extract = store_extract(storage%ab(s1,s2), &
                  path1(2:,:),path2(2:,:),info)
              ! Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storesym_extract
              !              Write(*,*) 'one side ',info%nrg
              If (.Not. storesym_extract) Return

            Else If (storage%init(s1,s2)%level == 4) Then  !** ATM detail
              If (path1(3,1) /= 0) Then  !** just one atom
                storesym_extract = store_extract(storage%ab(s1,s2), &
                    path1(2:,:),path2(2:,:),info)

              Else  !** all the atoms
                path1(3,1) = 1
                path1(3,2) = storage%init(s1,s2)%natoms1
                storesym_extract = store_extract(storage%ab(s1,s2), &
                    path1(2:,:),path2(2:,:),info)
                !                  Write(*,'(2a,i4,l2)') __FILE__," : ",__LINE__,storesym_extract
                !                  Write(*,*) 'one side ',a,info%nrg
                If (.Not. storesym_extract) Return
              End If

            Else 
              !** Check depth desired, return false if not available
              depth = Max(getdepth(path1(:,1)),getdepth(path2(:,1)))
              If ((storage%init(s1,s2)%level - 1) < depth) Then
                storesym_extract = .False.
                Return
              End If

              Write(0,'(2a,i6,a,i2)') __FILE__,":",__LINE__, &
                  " Unprepared for this detail level in storesym_extract ", &
                  storage%init(s1,s2)%level
              Stop              
            End If

          Case (1)  !** symmetric form 
#if DEBUG
            Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
            Write(*,*) 'match number: ',match
            Write(*,*) 'path1: ',path1(:,1),'   path2: ',path2(:,1)
            Write(*,*) '       ',path1(:,2),'          ',path2(:,2)
#endif

            If (match == 1) Then  !** spc we're interested in is 1st
              !** Extract the 'ab' portion
              path2(2,1) = 1
              path2(2,2) = storage%init(s1,s2)%nmoles2
              storesym_extract = store_extract(storage%ab(s1,s2), &
                  path1(2:,:),path2(2:,:),info)
              If (.Not. storesym_extract) Return

              !** Extract the 'ba' portion, need to specify multiple descents 
              !** to get to desired level pairs since %im is first
              path2(2,1) = 1
              path2(2,2) = storage%init(s1,s2)%nmoles2
              storesym_extract = store_extract(storage%ba(s1,s2), &
                  path1(2:,:),path2(2:,:),info)
              If (.Not. storesym_extract) Return

            Else If (match == 2) Then  !** spc we want is 2nd, reverse path
              !** Extract the 'ab' portion, reverse path order, specify 
              !** multiple descents to get to desired level pairs
              path2(2,1) = 1
              path2(2,2) = storage%init(s1,s2)%nmoles1
              storesym_extract = store_extract(storage%ab(s1,s2), &
                  path2(2:,:),path1(2:,:),info)
              If (.Not. storesym_extract) Return

              !** Extract the 'ba' portion, reverse path order
              path2(2,1) = 1
              path2(2,2) = storage%init(s1,s2)%nmoles1
              storesym_extract = store_extract(storage%ba(s1,s2), &
                  path2(2:,:),path1(2:,:),info)
              If (.Not. storesym_extract) Return

            Else If (match == 3) Then  !** spc we want is both  
              !** Must visit all entries where one of the molecules is mol,
              !** except for mol-mol entries.

              !** Extract the 'ab' portion in four parts to skip self-interaction
              !** picture a nmoles x nmoles array, extract all entries in a cross
              !** such that m1=mol or m2=mol, skipping mol,mol entry
              !** AB, indices 1 -> mol-1
              path2(2,1) = 1
              path2(2,2) = mol-1
              If (path2(2,1) <= path2(2,2)) Then
                storesym_extract = store_extract(storage%ab(s1,s2), &
                    path1(2:,:),path2(2:,:),info)
                If (.Not. storesym_extract) Return
                storesym_extract = store_extract(storage%ab(s1,s2), &
                    path2(2:,:),path1(2:,:),info)
                If (.Not. storesym_extract) Return
              End If

              !** AB, indices mol+1 -> nmoles
              path2(2,1) = mol+1
              path2(2,2) = storage%init(s1,s2)%nmoles2
              If (path2(2,1) <= path2(2,2)) Then
                storesym_extract = store_extract(storage%ab(s1,s2), &
                    path1(2:,:),path2(2:,:),info)
                If (.Not. storesym_extract) Return
                storesym_extract = store_extract(storage%ab(s1,s2), &
                    path2(2:,:),path1(2:,:),info)
                If (.Not. storesym_extract) Return
              End If

              !** Extract the 'ba' portion, same story
              !** BA, indices 1 -> mol-1
              path2(2,1) = 1
              path2(2,2) = mol-1
              If (path2(2,1) <= path2(2,2)) Then
                storesym_extract = store_extract(storage%ba(s1,s2), &
                    path1(2:,:),path2(2:,:),info)
                If (.Not. storesym_extract) Return
                storesym_extract = store_extract(storage%ba(s1,s2), &
                    path2(2:,:),path1(2:,:),info)
                If (.Not. storesym_extract) Return
              End If

              !** BA, indices mol+1 -> nmoles
              path2(2,1) = mol+1
              path2(2,2) = storage%init(s1,s2)%nmoles2
              If (path2(2,1) <= path2(2,2)) Then
                storesym_extract = store_extract(storage%ba(s1,s2), &
                    path1(2:,:),path2(2:,:),info)
                If (.Not. storesym_extract) Return
                storesym_extract = store_extract(storage%ba(s1,s2), &
                    path2(2:,:),path1(2:,:),info)
                If (.Not. storesym_extract) Return
              End If

            End If

          Case Default
            Write(0,'(2a,i4,a,i2)') __FILE__,": ",__LINE__, &
                ' Could not interpret form ',storage%init(s1,s2)%form
            Stop
          End Select

        End Do
      End Do

    Else If ((spc1 /= 0).And.(spc2 /= 0)) Then
      ! there was a stopsign here before, I know this works for the case of 
      ! spc level storage and if you want nrg of one spc with another -shaji

      If (storage%on(spc1,spc2)) Then
        Call storebase_inc(info,storage%ab(spc1,spc2)%total)
      Else
        !** Return successfully if the storage is not turned on
        Return
      Endif
      !#endif

    Else If (spc1 == 0) Then    !** full system with subset2
      Write(0,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(0,*) 'not ready for subset1 = full system'
      Stop

    End If

  End Function storesym_extract

  !-------------------------------------------------------------------------  
  ! Projects out a specified component of a force vector
  ! Requires:  storage -- Forcefield results storage structure
  !            subset -- subset identifier (spc,mol,atm)
  !            vec -- unit vector along which to remove force component
  !-------------------------------------------------------------------------    
  Subroutine storesym_projectout(storage,subset,vec)
    Type(Symmetric_Results), Intent(InOut)  :: storage
    Integer, Dimension(:), Intent(In)       :: subset
    Type(VecType), Intent(In)               :: vec

    Integer                  :: depth,s1,s2,nmols,natms
    Logical                  :: success
    Integer, Dimension(3,2)  :: path1,path2

    Do s1 = 1,storage%nspc
      Do s2 = s1,storage%nspc
        !** Skip out if this spc-spc pair is not pertinent
        If ((s1 /= subset(1).And.(s2 /= subset(1)))) Cycle
        If (.Not. storage%on(s1,s2)) Cycle

        !** Set up the paths to handle all forces on subset from other species
        path1(:,1) = (/subset(2:3),0/)
        path1(:,2) = (/subset(2:3),0/)
        nmols = storage%init(s1,s2)%nmoles2
        natms = storage%init(s1,s2)%natoms2
        path2(:,1) = (/1,1,0/)
        path2(:,2) = (/nmols,natms,0/)

        !** Handle the forces that are stored as forces on A from B
        success = store_projectout(storage%ab(s1,s2),path1,path2,vec) 
        Call checkandstop(success,__FILE__,__LINE__, &
            ' unexpected difficulty during force projection')

        !** Handle the forces that are stored as forces on B from A
        If (storage%init(s1,s2)%form > 0) Then
          success = store_projectout(storage%ba(s1,s2),path2,path1,vec) 
          Call checkandstop(success,__FILE__,__LINE__, &
              ' unexpected difficulty during force projection')
        Else
          Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              ' Storage usage not yet available for interact_allint'
          Stop
        End If
      End Do
    End Do

  End Subroutine storesym_projectout

  !----------------------------------------------------------------------------
  ! Extract the forces and put them into a 2-D array with indices given by
  ! (atom,molecule).  Useful for filling in species acceleration arrays.
  ! Requires:  storage -- structure to initialize
  !            spc -- species number for which to get forces
  !            forces -- 2D forces array to increment
  !----------------------------------------------------------------------------
  Subroutine storesym_getforces(storage,spc,forces)
    Type(Symmetric_Results), Intent(In)           :: storage
    Integer, Intent(In)                           :: spc
    Type(VecType), Dimension(:,:), Intent(InOut)  :: forces

    Integer            :: s,spc1,spc2
    Integer            :: a1,a2,m1,m2
    Integer            :: natoms1,natoms2
    Integer            :: nmoles1,nmoles2
    Logical            :: intra_also,reverse

    !** loop through the other species and get forces
    Do s = 1,storage%nspc

      !** pick upper triangle pair of the 2D spc-spc storage arrays
      !** the 'reverse' flag indicates that spc1 /= spc and we are
      !** looking for interactions stored in the BA portion
      If (s >= spc) Then
        spc1 = spc
        spc2 = s
        reverse = .False.
      Else
        spc1 = s
        spc2 = spc
        reverse = .True.
      End If

      If (.Not. storage%on(spc1,spc2)) Cycle

      natoms1 = storage%init(spc1,spc2)%natoms1
      natoms2 = storage%init(spc1,spc2)%natoms2
      nmoles1 = storage%init(spc1,spc2)%nmoles1
      nmoles2 = storage%init(spc1,spc2)%nmoles2
      intra_also = storage%init(spc1,spc2)%intra_also
      
      !** Handle the ATM-depth initialization case
      If (storage%init(spc1,spc2)%level == 4) Then  
      
        Select Case(storage%init(spc1,spc2)%form)
        Case (0)   !** map-type storage structure
          Do m1 = 1,nmoles1
            Do a1 = 1,natoms1
              forces(a1,m1) = forces(a1,m1) + &
                  storage%ab(spc1,spc2)%mi(m1)%mi(a1)%total%force
            End Do
          End Do
        
        Case (1)   !** symmetric ssbasic-type storage structure
          Do m1 = 1,nmoles1

            If ((.Not. reverse).Or.(spc1 == spc2)) Then
              Do m2 = 1,nmoles2  !** ab storage
                If ((spc1 == spc2).And.(m2 < m1)) Cycle
                If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle

                Do a1 = 1,natoms1
                  Do a2 = 1,natoms2
                    forces(a1,m1) = forces(a1,m1) + &
                        storage%ab(spc1,spc2)%mi(m1)%im(m2)%mm(a1,a2)%total%force
                  End Do
                End Do

              End Do
            End If

            If ((reverse).Or.(spc1 == spc2)) Then
              Do m2 = 1,nmoles2        !** ba storage
                If ((spc1 == spc2).And.(m2 < m1)) Cycle
                If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle
  
                Do a1 = 1,natoms1
                  Do a2 = 1,natoms2
                    forces(a2,m2) = forces(a2,m2) + &
                        storage%ba(spc1,spc2)%im(m2)%mi(m1)%mm(a1,a2)%total%force
                  End Do
                End Do
              End Do
            End If      

          End Do
        
        End Select

      !** Handle the MOL-depth initialization case
      Else If (storage%init(spc1,spc2)%level == 3) Then
      
        Select Case(storage%init(spc1,spc2)%form)  
        Case (0)   !** map-type storage structure
          Do m1 = 1,nmoles1
            Do a1 = 1,natoms1
              forces(a1,m1) = forces(a1,m1) + &
                  storage%ab(spc1,spc2)%mi(m1)%binfo(a1)%force
            End Do
          End Do
        
        Case (1)   !** symmetric ssbasic-type storage structure
!          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!          Write(*,*) spc1,spc2,intra_also,spc
          Do m1 = 1,nmoles1

            If ((.Not. reverse).Or.(spc1 == spc2)) Then
              Do m2 = 1,nmoles2  !** ab storage
                If ((spc1 == spc2).And.(m2 < m1)) Cycle
                If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle

                Do a1 = 1,natoms1
                  forces(a1,m1) = forces(a1,m1) + &
                      storage%ab(spc1,spc2)%mi(m1)%im(m2)%binfo(a1)%force
                End Do

              End Do
            End If

            If ((reverse).Or.(spc1 == spc2)) Then
              Do m2 = 1,nmoles2        !** ba storage
                If ((spc1 == spc2).And.(m2 < m1)) Cycle
                If ((spc1 == spc2).And.(m1 == m2).And.(.Not. intra_also)) Cycle

                Do a2 = 1,natoms2
                  forces(a2,m2) = forces(a2,m2) + &
                      storage%ba(spc1,spc2)%im(m2)%mi(m1)%binfo(a2)%force
                End Do

              End Do
            End If

          End Do
        
        End Select

      !** Handle the SPC-depth initialization case
      Else If (storage%init(spc1,spc2)%level == 2) Then
        Select Case(storage%init(spc1,spc2)%form)  
        Case (0)   !** map-type storage structure
          Do m1 = 1,nmoles1
            Do a1 = 1,natoms1
              forces(a1,m1) = forces(a1,m1) + &
                  storage%ab(spc1,spc2)%mi(m1)%binfo(a1)%force
            End Do
          End Do
        
        Case (1)   !** symmetric ssbasic-type storage structure
!          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!          Write(*,*) spc1,spc2,intra_also,spc
          If ((.Not. reverse).Or.(spc1 == spc2)) Then
            Do m1 = 1,nmoles1
              Do a1 = 1,natoms1
                forces(a1,m1) = forces(a1,m1) + &
                    storage%ab(spc1,spc2)%mi(m1)%binfo(a1)%force
              End Do
            End Do
          End If

          If ((reverse).Or.(spc1 == spc2)) Then
            Do m2 = 1,nmoles2        !** ba storage
              Do a2 = 1,natoms2
                forces(a2,m2) = forces(a2,m2) + &
                    storage%ba(spc1,spc2)%im(m2)%binfo(a2)%force
              End Do
            End Do
          End If

        End Select

      Else
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' storesym_getforces not prepared to handle this type of storage'
        Write(0,'(a,i2)') 'form:  ',storage%init(spc1,spc2)%form
        Write(0,'(a,i2)') 'level: ',storage%init(spc1,spc2)%level
        Stop

      End If

    End Do    

  End Subroutine storesym_getforces

  !----------------------------------------------------------------------------
  ! Set the subset mapping.  When the storage structure does not contain the 
  ! entire set of interactions within the system, the 'subsetonly' flag should
  ! be True.  This means that the interactions stored are only those of the
  ! specified subset with the entire system.
  ! Requires:  storage -- Symmetric results storage structure
  !            mapflag -- new setting 
  !            subset -- new subset specification
  !----------------------------------------------------------------------------
  Subroutine storesym_setmap(storage,mapflag,subset)
    Type(Symmetric_Results), Intent(InOut)  :: storage
    Logical, Intent(In)                     :: mapflag
    Integer, Dimension(3), Intent(In)       :: subset

    storage%subsetonly = mapflag
    storage%subset = subset

  End Subroutine storesym_setmap

  !-------------------------------------------------------------------------  
  ! Returns all the species-based energies
  ! Requires:  storage -- Symmetric results storage structure
  !            nrgs -- 2D array of spc-spc energies
  !-------------------------------------------------------------------------  
  Subroutine storesym_allnrgs(storage,nrgs)
    Type(Symmetric_Results), Intent(In)             :: storage
    Real(kind = RDbl), Dimension(:,:), Intent(Out)  :: nrgs

    Integer           :: i,j,spc1,spc2,size1,size2

    size1 = Size(nrgs,1)
    size2 = Size(nrgs,2)

    If ((size1 < storage%nspc).Or.(size2 < storage%nspc)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed 2D array is not large enough'
      Stop
    End If

    Do i = 1,storage%nspc
      Do j = 1,storage%nspc
        spc1 = i
        spc2 = j
        If (spc1 > spc2) Then  !** select upper triangle
          spc1 = j
          spc2 = i
        End If

        nrgs(i,j) = storage%ab(spc1,spc2)%total%nrg

        If (storage%init(spc1,spc2)%form > 0) Then
          nrgs(i,j) = nrgs(i,j) + storage%ba(spc1,spc2)%total%nrg
        End If

      End Do
    End Do

  End Subroutine storesym_allnrgs

  !-------------------------------------------------------------------------  
  ! This routine removes a molecule from the storage, copies the last 
  ! molecule in the list into the gap and zeros all the interactions
  ! in the former position of the last molecule.
  ! Requires:  storage -- Forcefield results storage structure
  !            nmoleslist -- A list of nmoles for each species in system
  !            spc -- the species number
  !            mol -- the molecule number
  ! Needed Improvements:
  ! 1) assumes that nothing is in the 'ba' parts
  !-------------------------------------------------------------------------    
  Subroutine storesym_delmol(storage,nmoleslist,spc,mol)
    Type(Symmetric_Results), Intent(InOut)  :: storage
    Integer, Dimension(:), Intent(In)       :: nmoleslist
    Integer, Intent(In)                     :: spc,mol

    Integer           :: spc1,spc2,m1,m2,match
    Integer           :: nmoles1,nmoles2,nmoles
    Logical           :: okflag

    Do spc1 = 1,storage%nspc
      Do spc2 = spc1,storage%nspc

        !** Skip out if no deletion needs to be done
        If (.Not. storage%on(spc1,spc2)) Cycle 
        If ((spc /= spc1).And.(spc /= spc2)) Cycle

        !** Get number of molecules in each species
        nmoles1 = nmoleslist(spc1)
        nmoles2 = nmoleslist(spc2)

#ifdef DEBUG
        Write(*,*) 'spc1,spc2: ',spc1,spc2,'   ',storage%init(spc1,spc2)%form
        Write(*,*) 'nmoles1, nmoles2: ',nmoles1,nmoles2
        Write(*,*) 'form: ',storage%init(spc1,spc2)%form
#endif

        !** Handle the simple SPC-detail case
        !** This doesn't seem to work, presumably because of the way that
        !** storesym_ptrs puts the new interactions into the storage.
        !** The code is setup not to use this, but instead to reevaluate
        !** the full system interactions after a deletion
        If (storage%init(spc1,spc2)%level == 2) Then
          If (spc1 /= spc) Cycle
          If (storage%nderivs > 0) Then

            !** AB portion: Move end molecule to vacancy if necessary
            If (mol /= nmoles1) Then
              okflag = store_copy(storage%ab(spc1,spc2)%mi(mol), &
                  storage%ab(spc1,spc2)%mi(nmoles1))
              Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
            End If
          
            !** Zero last molecule's storage, it was moved or should be deleted
            Call store_zero(storage%ab(spc1,spc2)%mi(nmoles1),.True.)

            If ((spc2 == spc).And.(storage%init(spc1,spc2)%form > 0)) Then
              !** BA portion: Move end molecule to vacancy if necessary
              If (mol /= nmoles1) Then
                okflag = store_copy(storage%ba(spc1,spc2)%im(mol), &
                    storage%ba(spc1,spc2)%im(nmoles1))
                Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
              End If
          
              !** Zero last molecule's storage, it was moved or should be deleted
              Call store_zero(storage%ba(spc1,spc2)%im(nmoles1),.True.)
            End If
          Else
            !** nothing to do, no molecular detail
          End If
          Cycle
        End If

        !** Do the deletion depending on the form of storage
        Select Case(storage%init(spc1,spc2)%form)
        Case (0)  !** map-type
          If (spc1 /= spc) Then
            Write(0,'(2a,i5,2a)') __FILE__,":",__LINE__, &
                ' Cannot currently delete fixed map species molecules'
            Stop
          End If

          !** Assuming spc == spc1
          nmoles = nmoles1 

          !** Move end molecule to vacancy if necessary
          If (mol /= nmoles) Then
            okflag = store_copy(storage%ab(spc1,spc2)%mi(mol), &
                storage%ab(spc1,spc2)%mi(nmoles))
            Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
          End If
          
          !** Zero the last molecule's storage, it was moved or should be deleted
          Call store_zero(storage%ab(spc1,spc2)%mi(nmoles1),.True.)
          
        Case (1)  !** fully symmetric
          !** Identify matching species, determine if copying is necessary
          match = 0
          If ((spc1 == spc).And.(mol < nmoles1)) match = 1
          If ((spc2 == spc).And.(mol < nmoles2)) match = match + 2

          !** Do deletion in 'AB' portion
          !** The relationship of 'spc' to spc1 and spc2 determines how
          !** the matrix of mol-mol interactions needs to be compressed.
          !** Molecules of spc1 and spc2 form the rows and columns, resp.
          If (match == 1) Then 
            !** Just need to delete a ROW of (mol1,mol2) matrix
            Do m2 = 1,nmoles2
              okflag = store_copy(storage%ab(spc1,spc2)%mi(mol)%im(m2), &
                  storage%ab(spc1,spc2)%mi(nmoles1)%im(m2))
              Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
            End Do
            
          Else If (match == 2) Then 
            !** Just need to delete a COLUMN of (mol1,mol2) matrix
            Do m1 = 1,nmoles1
              okflag = store_copy(storage%ab(spc1,spc2)%mi(m1)%im(mol), &
                  storage%ab(spc1,spc2)%mi(m1)%im(nmoles2))
              Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
            End Do

          Else If (match == 3) Then 
            !** Need to delete BOTH row and column of (mol1,mol2) matrix
            nmoles = nmoles1
            If (nmoles /= nmoles2) Then
              Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
              Stop
            End If

            !** Move row first, skipping 'mol' and (nmoles,nmoles) entry
            Do m2 = 1,nmoles-1
              If (m2 == mol) Cycle
              If (nmoles > m2) Cycle  !** confine to upper triangle
              If (mol > m2) Then   !** flip to maintain m1<m2 convention
!                Write(*,*) 'copying ',nmoles,m2,' -> ',m2,mol,'   flip'
                okflag = store_copy(storage%ab(spc1,spc2)%mi(m2)%im(mol), &
                    storage%ab(spc1,spc2)%mi(nmoles)%im(m2))
              Else
!                Write(*,*) 'copying ',nmoles,m2,' -> ',mol,m2
                okflag = store_copy(storage%ab(spc1,spc2)%mi(mol)%im(m2), &
                    storage%ab(spc1,spc2)%mi(nmoles)%im(m2))
              End If
              Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
            End Do

            !** Now, move column, skipping 'mol' and (nmoles,nmoles) entry   
            Do m1 = 1,nmoles-1
              If (m1 == mol) Cycle
              If (m1 > nmoles) Cycle  !** confine to upper triangle
              If (m1 > mol) Then   !** flip to maintain m1<m2 convention
!                Write(*,*) 'copying ',m1,nmoles,' -> ',mol,m1
                okflag = store_copy(storage%ab(spc1,spc2)%mi(mol)%im(m1), &
                    storage%ab(spc1,spc2)%mi(m1)%im(nmoles))
              Else
!                Write(*,*) 'copying ',m1,nmoles,' -> ',m1,mol
                okflag = store_copy(storage%ab(spc1,spc2)%mi(m1)%im(mol), &
                    storage%ab(spc1,spc2)%mi(m1)%im(nmoles))
              End If
              Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
            End Do

            !** Lastly, move (nmoles,nmoles) diagonally up to (mol,mol)
!            Write(*,*) 'copying ',nmoles,nmoles,' -> ',mol,mol
            okflag = store_copy(storage%ab(spc1,spc2)%mi(mol)%im(mol), &
                storage%ab(spc1,spc2)%mi(nmoles)%im(nmoles))
            Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
          End If 

          !** Do deletion in 'BA' portion if pertinent
          !** The relationship of 'spc' to spc1 and spc2 determines how
          !** the matrix of mol-mol interactions needs to be compressed.
          !** Molecules of spc1 and spc2 form the rows and columns, resp.
          If (storage%init(spc1,spc2)%form > 0) Then
            If (match == 1) Then 
              !** Just need to delete a COLUMN of (mol2,mol1) matrix
              Do m2 = 1,nmoles2
                okflag = store_copy(storage%ba(spc1,spc2)%im(m2)%mi(mol), &
                    storage%ba(spc1,spc2)%im(m2)%mi(nmoles1))
                Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
              End Do
              
            Else If (match == 2) Then 
              !** Just need to delete a ROW of (mol2,mol1) matrix
              Do m1 = 1,nmoles1
                okflag = store_copy(storage%ba(spc1,spc2)%im(mol)%mi(m1), &
                    storage%ba(spc1,spc2)%im(nmoles2)%mi(m1))
                Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
              End Do
  
            Else If (match == 3) Then 
              !** Need to delete BOTH row and column of (mol2,mol1) matrix
              nmoles = nmoles1
              If (nmoles /= nmoles2) Then
                Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
                Stop
              End If
  
              !** Move column first, skipping 'mol' and (nmoles,nmoles) entry
              Do m2 = 1,nmoles-1
                If (m2 == mol) Cycle
                If (nmoles > m2) Cycle  !** confine to upper triangle
                If (mol > m2) Then   !** flip to maintain m1<m2 convention
!                  Write(*,*) 'copying ',m2,nmoles,' -> ',mol,m2
                  okflag = store_copy(storage%ba(spc1,spc2)%im(mol)%mi(m2), &
                      storage%ba(spc1,spc2)%im(m2)%mi(nmoles))
                Else
!                  Write(*,*) 'copying ',m2,nmoles,' -> ',m2,mol
                  okflag = store_copy(storage%ba(spc1,spc2)%im(m2)%mi(mol), &
                      storage%ba(spc1,spc2)%im(m2)%mi(nmoles))
                End If
                Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
              End Do
  
              !** Now, move row, skipping 'mol' and (nmoles,nmoles) entry   
              Do m1 = 1,nmoles-1
                If (m1 == mol) Cycle
                If (m1 > nmoles) Cycle !** confine to upper triangle
                If (m1 > mol) Then   !** flip to maintain m1<m2 convention
!                  Write(*,*) 'copying ',nmoles,m1,' -> ',m1,mol
                  okflag = store_copy(storage%ba(spc1,spc2)%im(m1)%mi(mol), &
                      storage%ba(spc1,spc2)%im(nmoles)%mi(m1))
                Else
!                  Write(*,*) 'copying ',nmoles,m1,' -> ',mol,m1
                  okflag = store_copy(storage%ba(spc1,spc2)%im(mol)%mi(m1), &
                      storage%ba(spc1,spc2)%im(nmoles)%mi(m1))
                End If
                Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
              End Do
  
              !** Lastly, move (nmoles,nmoles) diagonally up to (mol,mol)
!              Write(*,*) 'copying ',nmoles,nmoles,' -> ',mol,mol
              okflag = store_copy(storage%ba(spc1,spc2)%im(mol)%mi(mol), &
                  storage%ba(spc1,spc2)%im(nmoles)%mi(nmoles))
              Call checkandstop(okflag,__FILE__,__LINE__,' level pair copying')
            End If 
          End If  !* if BA portion initialized
          
          !** Recalculate the 'match' if we want to just delete the last molecule
          If (match == 0) Then
            If (spc1 == spc) match = 1
            If (spc2 == spc) match = match + 2
          End If

          !** Zero the last molecule's AB storage, it was deleted or moved 
          If ((match == 1).Or.(match == 3)) Then   !**delete last row
            Do m2 = 1,nmoles2
              Call store_zero(storage%ab(spc1,spc2)%mi(nmoles1)%im(m2),.True.) 
            End Do
          End If
          If ((match == 2).Or.(match == 3)) Then   !**delete last column
            Do m1 = 1,nmoles1
              Call store_zero(storage%ab(spc1,spc2)%mi(m1)%im(nmoles2),.True.) 
            End Do
          End If

          !** Zero the last molecule's BA storage, it was deleted or moved 
          If (storage%init(spc1,spc2)%form > 0) Then
            If ((match == 1).Or.(match == 3)) Then   !**delete last column
              Do m2 = 1,nmoles2
!                Write(*,*) spc1,spc2,'deleting column',m2,nmoles1
                Call store_zero(storage%ba(spc1,spc2)%im(m2)%mi(nmoles1),.True.) 
              End Do
            End If
            If ((match == 2).Or.(match == 3)) Then   !**delete last row
              Do m1 = 1,nmoles1
!                Write(*,*) spc1,spc2,'deleting row ',nmoles2,m1
                Call store_zero(storage%ba(spc1,spc2)%im(nmoles2)%mi(m1),.True.) 
              End Do
            End If
          End If

        Case Default
          Write(0,'(2a,i4,a,i2)') __FILE__,": ",__LINE__, &
              ' Could not interpret form ',storage%init(spc1,spc2)%form
          Stop
        End Select

      End Do  !* spc2 loop
    End Do  !* spc1 loop

  End Subroutine storesym_delmol

  !-------------------------------------------------------------------------   
  ! This routine resizes the storage for a single species
  ! Requires:  storage -- Symmetric results storage structure
  !            nmoles -- NEW number of molecules of species spc
  !            spc -- the species number
  !-------------------------------------------------------------------------   
  Subroutine storesym_resize(storage,nmoles,spc)
    Type(Symmetric_Results), Intent(InOut)  :: storage
    Integer, Intent(In)                     :: nmoles,spc

    Integer                  :: spc1,spc2,size
    Integer                  :: nderivs,form,level,nmoles2,nmoles1
    Logical                  :: success,intra_also
    Type(Store_Level_Pair)   :: abtmp,batmp

    !** Get the number of derivatives
    nderivs = storage%nderivs

    !** Perform resizing for each species-species pair
    Do spc1 = 1,storage%nspc

      !** just a check, could comment out to speed
      Do spc2 = 1,spc1-1
        If (storage%on(spc1,spc2)) Then
          Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              ' Not setup to handle spc1 > spc2 resizing'
          Stop
        End If
      End Do

      Do spc2 = spc1,storage%nspc
        If (.Not. storage%on(spc1,spc2)) Cycle 
        intra_also = storage%init(spc1,spc2)%intra_also

        !** Skip unchanged pairs of species
        If ((spc1 /= spc).And.(spc2 /= spc)) Cycle

        !** Set the NEW number of molecules for each species
        nmoles1 = storage%init(spc1,spc2)%nmoles1
        nmoles2 = storage%init(spc1,spc2)%nmoles2
        If (spc1 == spc) nmoles1 = nmoles
        If (spc2 == spc) nmoles2 = nmoles

#ifdef DEBUG
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) 'resizing spc pair: ',spc1,spc2,nmoles
        Write(*,*) 'form:    ',storage%init(spc1,spc2)%form
        Write(*,*) 'intraon: ',storage%init(spc1,spc2)%intra_also
#endif 

        Select Case(storage%init(spc1,spc2)%form)
        Case (0)  !** map-type
          !** copy old stuff
          Call store_initcopy(abtmp,storage%ab(spc1,spc2),.True.)

          !** Destroy the old structure
          Call store_clean(storage%ab(spc1,spc2))

          form = storage%init(spc1,spc2)%form
          level = storage%init(spc1,spc2)%level
!          Write(*,*) 'spc,nmoles: ',spc1,nmoles1,spc2,nmoles2
!          Write(*,*) 'form,level,nderivs: ',form,level,nderivs
          
          !** Initialize the new, large structure
          Call storesym_initsspair(storage,spc1,nmoles1,spc2,nmoles2, &
              form,level,nderivs,intra_also)

          !** Copy old contents into new structure
          success = store_copy(storage%ab(spc1,spc2),abtmp)
          Call checkandstop(success,__FILE__,__LINE__,' Failure during copy')

          !** Destroy the temporary structure
          Call store_clean(abtmp)

        Case (1)  !** fully symmetric
          If (nmoles1 < storage%init(spc1,spc2)%nmoles1) Then
            Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
                ' new nmoles1 is less than existing, not allowed'
            Stop
          End If

          !** Copy old stuff to temporary storage
          Call store_initcopy(abtmp,storage%ab(spc1,spc2),.True.)
          Call store_initcopy(batmp,storage%ba(spc1,spc2),.True.)

          !** Destroy the old structure
          Call store_clean(storage%ab(spc1,spc2))
          Call store_clean(storage%ba(spc1,spc2))

          !** Setup form of new structure
          form = storage%init(spc1,spc2)%form
          level = storage%init(spc1,spc2)%level
          
          !** Initialize the new, larger structure
          Call storesym_initsspair(storage,spc1,nmoles1,spc2,nmoles2, &
              form,level,nderivs,intra_also)

          !** Copy old contents into new structure
          success = store_copy(storage%ab(spc1,spc2),abtmp)
          Call checkandstop(success,__FILE__,__LINE__,' Failure during copy')
          success = store_copy(storage%ba(spc1,spc2),batmp)
          Call checkandstop(success,__FILE__,__LINE__,' Failure during copy')

          !** Destroy the temporary structure
          Call store_clean(abtmp)
          Call store_clean(batmp)

        End Select

      End Do
    End Do
      
  End Subroutine storesym_resize

  !----------------------------------------------------------------------------
  ! Increment the the interaction storage of a particular molecule and atom.
  ! Only work on the 'ab' portion of the storage for now.  Assumes that the
  ! form of the storage is map-type.
  ! Requires:  storage -- Symmetric results storage structure
  !            subset1 -- first spc,mol,atm numbers
  !            subset2 -- second spc,mol,atm numbers
  !            nrg -- energy value to add to storage
  !            force -- optional force value to add to storage
  ! NOTE: much of this is not tested, don't trust it
  !----------------------------------------------------------------------------
  Subroutine storesym_inc(storage,subset1,subset2,nrg,force)
    Type(Symmetric_Results), Intent(InOut) :: storage
    Integer, Dimension(3), Intent(In)      :: subset1,subset2
    Real(kind=RDbl), Intent(In)            :: nrg         
    Type(VecType), Intent(In), Optional    :: force

    Integer            :: idx
    Integer            :: spc1,mol1,atm1,spc2,mol2,atm2

    !** Change the 'subset' input into a more understandable form, slow?
    spc1 = subset1(1)
    mol1 = subset1(2)
    atm1 = subset1(3)
    spc2 = subset2(1)
    mol2 = subset2(2)
    atm2 = subset2(3)

    Select Case (storage%init(spc1,spc2)%form)
    Case (0)    !** map-type

      Select Case (storage%init(spc1,spc2)%level)
      Case (2)    !** species-level
        If (storage%nderivs == 0) Then
          Call storebase_inc(storage%ab(spc1,spc2)%total,nrg)
        Else
          idx = storage%init(spc1,spc2)%natoms1*(mol1 - 1) + atm1
          Call storebase_inc(storage%ab(spc1,spc2)%binfo(idx),nrg)
          Call storebase_inc(storage%ab(spc1,spc2)%binfo(idx),force)
        End If
  
      Case (3)    !** molecule-level
        If (storage%nderivs == 0) Then
          Call storebase_inc(storage%ab(spc1,spc2)%mi(mol1)%total,nrg)
        Else
          Call storebase_inc(storage%ab(spc1,spc2)%mi(mol1)%binfo(atm1),nrg)
          Call storebase_inc(storage%ab(spc1,spc2)%mi(mol1)%binfo(atm1),force)
          !** make sure the MOL-level sum is ok
          Call storebase_inc(storage%ab(spc1,spc2)%mi(mol1)%total,nrg)
        End If
  
      Case (4)    !** atom-level
        If (storage%nderivs == 0) Then
          Call storebase_inc(storage%ab(spc1,spc2)%mi(mol1)%mi(atm1)%total,nrg)
        Else
          Call storebase_inc(storage%ab(spc1,spc2)%mi(mol1)%mi(atm1)%total,nrg)
          Call storebase_inc(storage%ab(spc1,spc2)%mi(mol1)%mi(atm1)%total,force)
          !** make sure the MOL-level sum is ok
          Call storebase_inc(storage%ab(spc1,spc2)%mi(mol1)%total,nrg)
        End If
  
      Case Default
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Could not interpret storage level'
        Stop
      End Select

    Case (1)    !** fully symmetric type
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'NOT TESTED'
      Stop

      Select Case (storage%init(spc1,spc2)%level)
      Case (2)    !** species-level
        If (storage%nderivs == 0) Then
          Call storebase_inc(storage%ab(spc1,spc2)%total,nrg)
        Else
          idx = storage%init(spc1,spc2)%natoms1*(mol1 - 1) + atm1
          Call storebase_inc(storage%ab(spc1,spc2)%binfo(idx),nrg)
          Call storebase_inc(storage%ab(spc1,spc2)%binfo(idx),force)
        End If
  
      Case (3)    !** molecule-level
        If (storage%nderivs == 0) Then
          Call storebase_inc(storage%ab(spc1,spc2)%mi(mol1)%im(mol2)%total,nrg)
        Else
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Write(*,*) spc1,mol1,atm1,'       ',spc2,mol2,atm2    
          Call storebase_inc( &
              storage%ab(spc1,spc2)%mi(mol1)%im(mol2)%binfo(atm1),nrg)
          Call storebase_inc( &
              storage%ab(spc1,spc2)%mi(mol1)%im(mol2)%binfo(atm1),force)
        End If
  
      Case (4)    !** atom-level
        If (storage%nderivs == 0) Then
          Call storebase_inc( &
              storage%ab(spc1,spc2)%mi(mol1)%im(mol2)%mm(atm1,atm2)%total,nrg)
        Else
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Write(*,*) spc1,mol1,atm1,'       ',spc2,mol2,atm2    
          Call storebase_inc( &
              storage%ab(spc1,spc2)%mi(mol1)%im(mol2)%mm(atm1,atm2)%total,nrg)
          Call storebase_inc( &
              storage%ab(spc1,spc2)%mi(mol1)%im(mol2)%mm(atm1,atm2)%total,force)
        End If
  
      Case Default
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Could not interpret storage level'
        Stop
      End Select

    Case Default
      Write(0,'(2a,i4,a,i2)') __FILE__,": ",__LINE__, &
          ' Cannot currently deal with form: ',storage%init(spc1,spc2)%form
      Stop
    End Select

  End Subroutine storesym_inc

  !----------------------------------------------------------------------------
  ! Get the number of derivatives in the structure
  ! Requires:  storage -- Symmetric results storage structure
  !----------------------------------------------------------------------------
  Integer Function storesym_nderivs(storage)
    Type(Symmetric_Results), Intent(In)   :: storage

    storesym_nderivs = storage%nderivs

  End Function storesym_nderivs

  !-------------------------------------------------------------------------     
  ! Scales all the species-species based energies
  ! Requires:  storage -- Forcefield results storage structure
  !            factors -- 2D array of spc-spc energy scaling factors
  !-------------------------------------------------------------------------     
  Subroutine storesym_scalenrgs(storage,factors)
    Type(Symmetric_Results), Intent(InOut)          :: storage
    Real(kind = RDbl), Dimension(:,:), Intent(In)   :: factors

    Integer         :: i,j

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

    !** scale the interactions
    Do i = 1,storage%nspc
      Do j = 1,storage%nspc
        If (.Not. storage%on(i,j)) Cycle

        Call store_scalenrgs(storage%ab(i,j),factors(i,j))
        
        If (storage%init(i,j)%form > 0) Then
          Call store_scalenrgs(storage%ba(i,j),factors(i,j))
        End If

      End Do
    End Do

  End Subroutine storesym_scalenrgs

  !-------------------------------------------------------------------------     
  ! Scales all the species-species based interactions 
  ! Requires:  storage -- Forcefield results storage structure
  !            factor -- scaling factor
  !-------------------------------------------------------------------------     
  Subroutine storesym_scaleall(storage,factor)
    Type(Symmetric_Results), Intent(InOut)   :: storage
    Real(kind = RDbl), Intent(In)            :: factor

    Integer         :: i,j

    !** scale the interactions
    Do i = 1,storage%nspc
      Do j = 1,storage%nspc
        If (.Not. storage%on(i,j)) Cycle

        Call store_scaleall(storage%ab(i,j),factor)
        
        If (storage%init(i,j)%form > 0) Then
          Call store_scaleall(storage%ba(i,j),factor)
        End If

      End Do
    End Do

  End Subroutine storesym_scaleall

  !----------------------------------------------------------------------------
  ! Pseudo-Hack for dumping information about a molecule
  ! Requires:  storage -- Symmetric results storage structure
  !            spc -- species number
  !            mol -- molecule number
  !            indent -- no. of spaces from the left margin
  !            optunit -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine storesym_dumpmol(storage,spc,mol,indent,optunit)
    Type(Symmetric_Results), Intent(In)   :: storage
    Integer, Intent(In)                   :: spc,mol,indent
    Integer, Optional, Intent(In)         :: optunit

    Integer                    :: a,natoms,spc2,unitno
    Character(len=indent)      :: blank

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Call storesym_display(storage,.False.,indent,unitno)

    Do spc2 = 1,storage%nspc
      If (.Not. storage%on(spc,spc2)) Cycle

      !** skip non-maptype pairs
      If (storage%init(spc,spc2)%form /= 0) Cycle

      Write(unitno,'(a,i3)') 'species 2 = ',spc2
      Do a = 1,storage%init(spc,spc2)%natoms1
        Write(unitno,'(a,i4,3x,3f11.5)') blank,a, &
            storage%ab(spc,spc2)%mi(mol)%binfo(a)%force
      End Do
    End Do

  End Subroutine storesym_dumpmol

  !----------------------------------------------------------------------------
  ! Returns the 'depth' of a system subset specification.  This is the level
  ! at which the specification ends.  Returned depths and examples of subsets:
  !   0 = full system      subset = (/0,0,0/)
  !   1 = species level    subset = (/a,0,0/)
  !   2 = molecule level   subset = (/a,b,0/)
  !   3 = atom level       subset = (/a,b,c/)
  ! Requires:  subset -- subset specification
  !            idstring -- optional returned string identifier
  !----------------------------------------------------------------------------
  Integer Function storesym_depth(subset,idstring)
    Integer, Dimension(3), Intent(In)        :: subset
    Character(len=3), Intent(Out), Optional  :: idstring

    Integer         :: i

    Do i = 3,0,-1
      storesym_depth = i
      If (i == 0) Exit
      If (subset(i) /= 0) Exit
    End Do    

    If (Present(idstring)) Then
      idstring = store_levels(storesym_depth)
    End If

  End Function storesym_depth

  !----------------------------------------------------------------------------
  ! Check two symmetric storage structure to see if they are equal.  If a unit
  ! number and indent specific are included, the routine will dump output. 
  ! Otherwise, only the True/False answer will be returned.
  ! Requires:  storage1 -- 1st Symmetric results storage structure
  !            storage2 -- 2nd Symmetric results storage structure
  !            tol -- tolerance for calling equal
  !            indent -- no. of spaces from the left margin
  !            unit -- optional display unit number
  !----------------------------------------------------------------------------
  Logical Function storesym_chkequal(storage1,storage2,tol,indent,unit)
    Type(Symmetric_Results), Intent(In)   :: storage1,storage2
    Real(kind=RDbl), Intent(In)           :: tol
    Integer, Intent(In), Optional         :: indent,unit

    Integer                       :: i,j,k,s1,s2
    Logical                       :: success,dump,equal
    Character(len=strLen)         :: string1,string2,string3
    Character(len=xlstrLen)       :: startstring,string

    storesym_chkequal = .True.

    dump = .False.
    If ((Present(indent)).And.(Present(unit))) dump = .True.

    !** Compare the total energy
    string = ''
    If (.Not. storebase_chkequal(storage1%total,storage1%total,tol,string)) Then
      If (dump) Then
        Write(unit,'(20a)') (' ',k=1,indent), &
            'Total interactions not equal: ',Trim(string)
        storesym_chkequal = .False.
      End If
    End If
    
    If (dump) Write(unit,'(20a)') (' ',k=1,indent),'AB storage structure'
    Do i = 1,storage1%nspc
      Do j = i,storage1%nspc
        If (.Not. storage1%on(i,j)) Cycle
        If (.Not. storage2%on(i,j)) Then
          Write(0,'(1x,2a,i4,a)') __FILE__,' : ',__LINE__, &
              ' Passed Symmetric Storage structures do not have same shape'
          Stop            
        End If
        string1 = int2str(i)
        string2 = int2str(j)    
        Write(startstring,'(4a))') 's:',Trim(string1),',',Trim(string2)

        If (dump) Then
          equal = store_chkequal(storage1%ab(i,j),storage2%ab(i,j), &
              startstring,tol,indent,unit)
        Else
          equal = store_chkequal(storage1%ab(i,j),storage2%ab(i,j), &
              startstring,tol)
        End If
        If (.Not. equal) storesym_chkequal = .False.

      End Do
    End Do

    If (dump) Write(unit,'(20a)') (' ',k=1,indent),'BA storage structure'
    Do i = 1,storage1%nspc
      Do j = i,storage1%nspc
        If (.Not. storage1%on(i,j)) Cycle
        If (.Not. storage2%on(i,j)) Then
          Write(0,'(1x,2a,i4,a)') __FILE__,' : ',__LINE__, &
              ' Passed Symmetric Storage structures do not have same shape'
          Stop            
        End If
        string1 = int2str(i)
        string2 = int2str(j)    
        Write(startstring,'(4a))') 's:',Trim(string1),',',Trim(string2)
        If (dump) Then
          equal = store_chkequal(storage1%ba(i,j),storage2%ba(i,j), &
              startstring,tol,indent,unit)
        Else
          equal = store_chkequal(storage1%ba(i,j),storage2%ba(i,j), &
              startstring,tol)
        End If
        If (.Not. equal) storesym_chkequal = .False.
      End Do
    End Do

  End Function storesym_chkequal

  !----------------------------------------------------------------------------
  ! Display the symmetric storage structure.
  ! Requires:  storage -- Symmetric results storage structure
  !            skip -- True => will skip display if nrg,grad(1) < tolerance
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine storesym_display(storage,skip,indent,unitno)
    Type(Symmetric_Results), Intent(In)   :: storage
    Logical, Intent(In)                   :: skip
    Integer, Intent(In)                   :: indent
    Integer, Optional, Intent(In)         :: unitno

    Integer                           :: i,j,unit,s1,s2
    Character(len=indent)             :: blank
    Character(len=strLen)             :: string1,string2,string3
    Character(len=xlstrLen)           :: startstring,string

    blank = Repeat(' ',indent)
    
    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    string = storebase_disp(storage%total)
    Write(unit,'(3a)') blank,'Top: ',Trim(string)
    
    Write(unit,'(2a)') blank,'AB storage structure'
    Do i = 1,storage%nspc
      Do j = i,storage%nspc
        string1 = int2str(i)
        string2 = int2str(j)    
        string3 = 'OFF'
        If (storage%on(i,j)) string3 = 'ON'
        If (storage%init(i,j)%intra_also) Then
          Write(string3,'(a,3x,a)') Trim(string3),'INTRA ALSO'
        End If
        Write(unit,'(2x,3a,1x,a,2x,2(2x,a,i2),3x,a)') blank,'Species Pair: ', &
            Trim(string1),Trim(string2),'Form: ',storage%init(i,j)%form,&
            'Level: ',storage%init(i,j)%level,Trim(string3)

        !** Display interaction information using recursive descent if ON
        If (storage%on(i,j)) Then
          startstring = storebase_disp(storage%ab(i,j)%total)
          Write(unit,'(2x,3a)') blank,'spc-spc top: ', Trim(startstring)
          Write(startstring,'(4a))') 's:',Trim(string1),',',Trim(string2)
          Call store_display(storage%ab(i,j),startstring,skip,indent+2,unit)
        End If

        !** Error checking
        If ((.Not.storage%on(i,j)).And. &
            (store_idbranch(storage%ab(i,j)) /= -1)) Then
          Write(unit,'(2x,2a)') blank,'WARNING: spc-spc off, but initialized'
        End If
      End Do
    End Do

    Return

    Write(unit,'(2a)') blank,'BA storage structure'
    Do i = 1,storage%nspc
      Do j = i,storage%nspc
        string1 = int2str(i)
        string2 = int2str(j)
        string3 = 'OFF'
        If (storage%on(i,j)) string3 = 'ON'
        If (storage%init(i,j)%intra_also) Then
          Write(string3,'(a,3x,a)') Trim(string3),'INTRA ALSO'
        End If
        Write(unit,'(2x,3a,1x,a,2x,2(2x,a,i2),3x,a)') blank,'Species Pair: ', &
            Trim(string1),Trim(string2),'Form: ',storage%init(i,j)%form,&
            'Depth: ',storage%init(i,j)%level,Trim(string3)

        !** Display interaction information using recursive descent if ON
        If (storage%on(i,j)) Then
          startstring = storebase_disp(storage%ba(i,j)%total)
          Write(unit,'(2x,3a)') blank,'spc-spc top: ', Trim(startstring)
          Write(startstring,'(4a))') 's:',Trim(string1),',',Trim(string2)
          Call store_display(storage%ba(i,j),startstring,skip,indent+2,unit)
        End If

        !** Error checking
        If ((.Not.storage%on(i,j)).And. &
            (store_idbranch(storage%ba(i,j)) /= -1)) Then
          Write(unit,'(2x,2a)') blank,'WARNING: spc-spc off, but initialized'
        End If
      End Do
    End Do

  End Subroutine storesym_display

  !----------------------------------------------------------------------------
  ! Nullify the symmetric storage structure
  ! Requires:  storage -- Symmetric results storage structure
  !----------------------------------------------------------------------------
  Subroutine storesym_null(storage)
    Type(Symmetric_Results), Intent(InOut)   :: storage

    Integer         :: spc1,spc2,error

    Call storebase_null(storage%total)
    Nullify(storage%on)
    Nullify(storage%init)
    Nullify(storage%ab)
    Nullify(storage%ba)

  End Subroutine storesym_null

  !----------------------------------------------------------------------------
  ! Clean the symmetric storage structure
  ! Requires:  storage -- Symmetric results storage structure
  !----------------------------------------------------------------------------
  Subroutine storesym_clean(storage)
    Type(Symmetric_Results), Intent(InOut)   :: storage

    Integer         :: spc1,spc2,error

    Call storebase_clean(storage%total)

    Do spc1 = 1,storage%nspc
      Do spc2 = 1,storage%nspc
        Call store_clean(storage%ab(spc1,spc2))
        Call store_clean(storage%ba(spc1,spc2))
      End Do
    End Do

    Deallocate(storage%on, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'on')    
    Deallocate(storage%init, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'init')    
    Deallocate(storage%ab, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'ab')    
    Deallocate(storage%ba, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'ba')    

  End Subroutine storesym_clean

End Module storesym
