!------------------------------------------------------------------------------
! This modules handles calculating interactions that require 
! atom-system summations, like Ewald.
! NOTE: Forces are scaled in this routine to accelerations
! before being returned
!-----------------------------------------------------------------------------
Module sssum

  Use defaults, Only: RDbl, strLen, lstrLen, MAX_SORBS, zero, scalef, &
      INCR_SIZE
  Use vector, Only: VecType, Operator(+), Assignment(=), Operator(*), &
      Operator(-),mag, Operator(/)
  Use ewald, Only: EwaldParams, ewald_init, ewald_initCalc, ewald_idstring, &
      ewald_getline, ewald_display, ewald_getinteractionHOT, &
      ewald_getinteractionSngl, ewald_getCoreShellCorr
  Use molecules, Only: molecules_name, molecules_getnsorbs, &
      molecules_getnatoms, molecules_ischarged, molecules_getcharges, &
      molecules_getnatomtypes, molecules_isCoreShell
  Use utils, Only: stripcmnt, split, tolower, toupper, combine, findint,&
      allocErrDisplay
  Use config, Only: AtMolCoords, config_getnmoles, config_totalvecs, &
      config_isfixed
  Use simcell, Only: SimCell_Params, simcell_getell, simcell_getnuc
  Use atom, Only: atom_invmass, atom_invmassMult
  Use storebase, Only: EnergyPlus, storebase_inc
  Use storesym, Only: Symmetric_Results, storesym_nderivs, storesym_inc

  Implicit None
  Save

  Private
  Public :: SSSumParams, sssum_init, sssum_initCalc, sssum_idstring, &
       sssum_getmsinteraction, sssum_display, Assignment(=), &
       sssum_getssinteraction, sssum_clean, sssum_ssint, sssum_asintHOT, sssum_intraon

  Interface sssum_getmsinteraction
     Module Procedure sssum_getmsint
     Module Procedure sssum_getmsintHOT
  End Interface

  Interface sssum_getssinteraction
     Module Procedure sssum_getssint
!!$     Module Procedure sssum_getssintHOT
  End Interface
  
  Interface Assignment(=)
    Module Procedure sssum_copyParams
  End Interface

  Type SSSumParams 
    Logical                :: fast       ! Fast or slow interaction
    Logical                :: intraOn    ! Allows intramolecular sum calculations 
    Logical                :: fixedOnly  ! Only include FIXED sorbates in the calc
    Character(len=strLen)  :: model
    Character(len=lstrLen) :: line
    Type(EwaldParams), Pointer :: ewald
!!$    Type(DSumParams)        :: dsum
  End Type SSSumParams

  Type(VecType), Dimension(:), Pointer   :: atvec
  Real(Kind=RDbl), Dimension(:), Pointer :: charge,clists

  !** Identifying string that we look for in the sorbate-sorbate
  !** interaction file
  Character(len=strLen), Parameter :: sssum_idstring = "SUM"

  !** Arrays that are used to pass the atom vectors and charges
  Real(Kind=RDbl), Dimension(:), Allocatable  :: charges
  Real(Kind=RDbl), Dimension(:), Allocatable  :: clist

  !** Information about the interaction types
  Integer, Dimension(:), Allocatable :: ewaldList
  Integer :: newald
  Logical :: firstKvec

Contains

  !----------------------------------------------------------------------------
  ! Initializes the data for this module
  ! Things for ewald are a little strange. We only call ewald_init once in
  ! order to allocate the kvec array. After the first call, we check to make
  ! sure the subsequent sorb-sorb pairs that use ewald have the same 
  ! parameters.
  !----------------------------------------------------------------------------
  Subroutine sssum_init(params,sorb1,sorb2,line)
    Type(SSSumParams), Intent(InOut) :: params
    Integer, Intent(In)              :: sorb1, sorb2
    Character(*), Intent(In)         :: line

    Integer                                  :: error, nfields, nsorbs, i
    Logical                                  :: endFlag
    Character(len=lstrLen)                   :: newline, checkLine
    Character(len=strLen), Dimension(strLen) :: fields

    !** Note that we haven't computed the kvecs yet
    firstKvec = .False.

    !** Store the line and strip out comments
    newline = stripcmnt(line)
    params%line = newline
    nfields = split(newline,fields)

    !** Initialize the fields that govern the interaction between two
    !** molecules
    params%fast = .FALSE.
    params%intraOn = .False.
    params%fixedOnly = .False.

    !** Initialize the do loop flag
    endFlag = .False.

    !** Loop through the fields to find the important information
    Do i = 1, nfields

      Select Case (toupper(fields(i)))

      Case ('FAST')
        params%fast = .True.

      Case ('INTRAON')
        params%intraOn = .True.

      Case ('FIXED')
        params%fixedOnly = .True.
        
      Case (ewald_idstring)
        !** Using Ewald
        params%model = fields(i)
        newline = combine(fields(i:nfields))

        !** Check to see if the ewaldList is allocated
        If (.Not.Allocated(ewaldList)) Then
          nsorbs = molecules_getnsorbs()
          Allocate(ewaldList(nsorbs),stat=error)
          If (error /= 0) Then
            Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
                " Could not allocate ewaldList"
            Stop
          End If
          ewaldList = 0
          newald = 0
        End If
        !** We won't use the sorb1 sorb2 info for now.

        !** Check the line stored in Ewald, 
        !** see If it matches the previous.
        checkLine = ewald_getline()
        If (checkLine == newline) Then
          !** Add the species pair to ewaldList
          If (findint(ewaldList,sorb1) == 0) Then
            newald = newald+1
            ewaldList(newald) = sorb1
          End If
          If (findint(ewaldList,sorb2) == 0) Then
            newald = newald+1
            ewaldList(newald) = sorb2
          End If
          !** Ok, same as before. Continue.
          Return
        Else If (checkLine(1:5) /= newline(1:5)) Then
          !** It hasn't been initialized yet.
          !** Calling ewald for the first time. Initialize it.
          Allocate(params%ewald,stat=error)
          If (error/=0) Then
            Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
                " Could not allocate space for params%ewald."
            Stop
          End If
          Call ewald_init(newline)
          !** Add the sorbate pair to ewaldList
          If (findint(ewaldList,sorb1) == 0) Then
            newald = newald+1
            ewaldList(newald) = sorb1
          End If
          If (findint(ewaldList,sorb2) == 0) Then
            newald = newald+1
            ewaldList(newald) = sorb2
          End If
        Else
          !** Uh-oh, they don't match like they should
          !** Report the error and stop
          Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
              " EWALD parameter lines MUST match. Conflicting lines: "
          Write(0,'(3x,a)') newline
          Write(0,'(3x,a)') checkLine
          Stop
        End If

        !** set the exit flag to end do loop
        endFlag = .True.
    
!!$    Case(tolower(dsum_idstring))
!!$      !** Using direct summation
!!$      Call dsum_init(newline)

      Case Default
        Write(0,'(2a,i5,6a)') __FILE__,":",__LINE__, &
            " Unrecognized summation parameter ",trim(fields(i)), &
            " for spc-spc pair ",Trim(molecules_name(sorb1)),"-", &
            Trim(molecules_name(sorb2))
        Write(0,'(2a)') "NOTE: All model parameters must come AFTER the ",&
            "model name"
        Stop
      End Select

      If (endFlag) exit

    End Do
    IF(.NOT. params%fast) THEN
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) 'Simulation will not run without the fast flag in the sorb_sorb_file!'
        STOP
    END IF

  End Subroutine sssum_init

  !----------------------------------------------------------------------------
  ! Copies the contents of one variable of type SSSumParams to another of
  ! the same type
  !----------------------------------------------------------------------------
  Subroutine sssum_copyParams(copy,orig)
    Type(SSSumParams), Intent(In) :: orig
    Type(SSSumParams), Intent(Out):: copy

    copy%fast = orig%fast
    copy%model = orig%model
    copy%line = orig%line
    copy%ewald => orig%ewald
    copy%fixedOnly = orig%fixedOnly
  End Subroutine sssum_copyParams


  !----------------------------------------------------------------------------
  ! Initializes anything that needs to be done immediately before the 
  ! sorbate-sorbate interactions are calculated.
  !----------------------------------------------------------------------------
  Subroutine sssum_initCalc(params,sorbPair,sorbates,sparams,fast)
    Type(SSSumParams), Dimension(:), Intent(InOut) :: params
    Integer, Dimension(:,:), Intent(In)            :: sorbPair
    Type(AtMolCoords), Dimension(:), Intent(In)    :: sorbates
    Type(SimCell_Params), Intent(In)               :: sparams
    Logical, Intent(In)                            :: fast

    Integer                              :: lastIndx, i, natoms, error, nmoles
    Integer                              :: j, currentSize, requiredSize
    Integer                              :: nwald, thisSpec
    Logical, Save                        :: firsttime = .True.
    Real(Kind=RDbl), Dimension(3)        :: ell
    Integer, Dimension(Size(sorbates,1)) :: specList

           
    !** Get the size of sorbates 
    requiredSize = config_totalvecs(sorbates)

    If (.Not. firsttime) Then 
      !** If the temporary arrays are allocated, then make sure
      !** they are the right size.
      currentSize = Size(clists,1)
      If (currentSize < requiredSize) Then
        
        !** Increment by more than the required amount so that
        !** during simulations where N is not fixed we do not
        !** need to reallocate each time
        currentSize = currentSize + INCR_SIZE
        
        !** Allocate clists
        Deallocate(clists,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_dealloc")
        Allocate(clists(currentSize),STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_alloc")
        
        !** Allocate charge
        Deallocate(charge,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"charge_dealloc")
        Allocate(charge(currentSize),STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"charge_alloc")
        
        !** Allocate atvec
        Deallocate(atvec,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"atvec_dealloc")
        Allocate(atvec(currentSize),STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"atvec_alloc")
      End If

    Else
      firsttime = .False.
      Allocate(clists(requiredSize),stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i5,a,i6)') __FILE__,":",__LINE__, &
            " Could not allocate clists of size ",requiredSize
      End If
      Allocate(atvec(requiredSize),stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i5,a,i6)') __FILE__,":",__LINE__, &
            " Could not allocate atvec of size ",requiredSize
      End If
      Allocate(charge(requiredSize),stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i5,a,i6)') __FILE__,":",__LINE__, &
            " Could not allocate charge of size ",requiredSize
      End If
    End If
     
    !** Initialize all the necessary summation interaction types

!!$      !MDEBUG
!!$      Write(0,*) "Initializing ewald"

    !** Get the edge lengths that ewald needs for calculating the sum.
    !** The length we need should correspond to the length of the
    !** simulation cell included in the ewald calculation.
    ell = simcell_getell(sparams)
    
    !** Reshape the structures into arrays. Here we loop through 
    !** each of the molecule types that are included in the Ewald
    !** calculation, get their coordinates and charges, and place
    !** them all in 2 large arrays, one for charges, the other for
    !** coords.
    lastIndx = 0
    nwald = 0

    Do i = 1, Size(params,1)
      
      !** Check to see if we are only using fixed molecule types in
      !** the summation. If we are, check the firstKvec logical. If
      !** it is true, we've already calculated the Kvec array and can
      !** return.
      If (params(i)%fixedOnly) Then
        If (firstKvec) Return
      End If

      Do j = 1, Size(sorbPair,2)

        !** Get the species number
        thisSpec = sorbPair(i,j)

        !** Check to see if we are looking for fixed sorbates
        !** only. If so, then cycle if this isn't fixed.
        If (params(i)%fixedOnly) Then
          If (.Not.config_isfixed(sorbates(thisSpec))) Cycle
        End If

        !** Add the species to the list so we don't repeat
        If (nwald == 0) Then
          nwald = nwald + 1
          specList(nwald) = thisSpec
        Else If (findint(specList(1:nwald),thisSpec) == 0) Then
          nwald = nwald + 1
          specList(nwald) = thisSpec
        Else
          Cycle
        End If
      
        !** Get the number of atoms and molecules for this type
        natoms = molecules_getnatoms(thisSpec)
        nmoles = config_getnmoles(sorbates(thisSpec))

        !** position vectors
        atvec(lastIndx+1:lastIndx+natoms*nmoles) = &
            Reshape(sorbates(thisSpec)%coords(1:natoms,1:nmoles)%r, &
            (/natoms*nmoles/))
        !** atom charges
        Call molecules_getcharges(thisSpec,clists(1:natoms))
        charge(lastIndx+1:lastIndx+natoms*nmoles) = &
            Reshape(clists(1:natoms),(/natoms*nmoles/),clists(1:natoms))
        lastIndx = lastIndx+natoms*nmoles
      End Do ! Do j = 1, Size(sorbPairs,2)
    End Do ! Do i = 1, Size(params,1)
      
    If (nwald == 0) Return

    !** Call the Ewald initCalc routine
    Call ewald_initCalc(atvec(1:lastIndx),charge(1:lastIndx),ell)

    firstKvec = .True.

  End Subroutine sssum_initCalc

  !----------------------------------------------------------------------------
  ! Calculates the potential and forces due to sorb1 and sorb2 interacting
  ! using a summation type of calculation. Returns the accelerations on 
  ! each of the sorbates.
  !
  ! Some notes: The way in which the code is currently set up, we are actually
  ! double counting potentials for the Ewald summation because we loop over
  ! all atoms in both sorbates, so we must multiply the FINAL potential by
  ! 1/2. This would be true of the forces as well IF we added atvecs2 back
  ! into asorb2. We are not doing that, however, so we don't need the 1/2.
  ! Requires:  params -- summation parameters
  !            sparams -- simulation cell information
  !            species -- species coordinate structure
  !            spc1 -- first species number
  !            spc2 -- second species number
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !            ffout -- symmetric results forcefield output
  !----------------------------------------------------------------------------
  Logical Function sssum_ssint(params,species,sparams,spc1,spc2,fast,ffout)
    Type(SSSumParams), Intent(In)               :: params
    Type(SimCell_Params), Intent(In)            :: sparams
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: spc1, spc2
    Logical, Intent(In)                         :: fast
    Type(Symmetric_Results), Intent(InOut)      :: ffout

    Integer             :: natoms1, natoms2, a, nmoles1, nmoles2, m, i, j, jskipatom
    Integer             :: natypes1, natypes2, csAtom
    Integer             :: low, med, high,currentSize,requiredSize,error
    Integer             :: natvecs1, natvecs2, aa, mm
    Integer             :: level, nderivs
    Logical             :: sameSpc, getf
    Logical             :: iscoreshell, selfTerm, isCore
    Real(Kind=Rdbl)     :: utemp, amass, ucorr
    Integer, Dimension(Size(species(spc1)%coords,1))         :: atypes1
    Integer, Dimension(Size(species(spc2)%coords,1))         :: atypes2
    Real(Kind=RDbl), Dimension(Size(species(spc1)%coords)) :: charges1
    Real(Kind=RDbl), Dimension(Size(species(spc2)%coords)) :: charges2
    Real(Kind=RDbl), Dimension(Size(species(spc1)%coords,1)) :: atInvMass1
    Real(Kind=RDbl), Dimension(Size(species(spc2)%coords,1)) :: atInvMass2
    Type(VecType)                                            :: fcorr,ftemp
    Type(VecType), Dimension(Size(species(spc2)%coords,1)* &
        Size(species(spc2)%coords,2))                        :: atvecs2
    Type(VecType), Dimension(Size(atvecs2,1))                :: ftemp2
    Integer, Dimension(Size(atvecs2,1)) :: skipMe !** atom numbers to skip
    Integer                             :: nskip  !** molecule number to skip

    Type(VecType)      :: dummy

    !** this is a hack that keeps the storesym_inc call the same
    !** NAG compiler isn't flexible enough to handle different calls
    dummy = VecType(0.0_RDbl)

    !** Set the atom and molecules to skip to zero
    skipMe = 0
    nskip = 0

    !** Check size of temporary array for holding charges
    currentSize = Size(clists,1)
    requiredSize= config_totalvecs(species)

    !** Size the array that holds the charges if necessary
    If (currentSize < requiredSize) Then
      Deallocate(clists,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_dealloc")
      Allocate(clists(requiredSize+INCR_SIZE),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_alloc")
    End If

    !** Check to see if we are calculating forces
    getf = .False.
    nderivs = storesym_nderivs(ffout)
    If (nderivs > 0) getf = .True.
!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,*) nderivs,getf

    !** zero the potential and the forces
    ftemp = VecType(zero)

    !** Check if the species are the same.
    !** We need to check the atvecs and make sure that for
    !** spc1 == spc2, molec1 == molec2, the interaction
    !** between atom1 == atom2 is not calculated.
    sameSpc = .False.
    If (spc1 == spc2) sameSpc = .True.

    !** We need the atom masses in order to scale the forces. Get the
    !** atom types so we can look up the atom masses later.
    natypes1 = molecules_getnatomtypes(spc1,atypes1,.True.)
    natypes2 = molecules_getnatomtypes(spc2,atypes2,.True.)
    
    !** Get the number of atoms and molecules for each type
    natoms1 = molecules_getnatoms(spc1)
    natoms2 = molecules_getnatoms(spc2)
    nmoles1 = config_getnmoles(species,spc1)
    nmoles2 = config_getnmoles(species,spc2)

    !** Calculate the total number of entries in the temporary arrays,
    !** e.g., natoms1*nmoles1
    natvecs1 = natoms1*nmoles1
    natvecs2 = natoms2*nmoles2

    !** Get the inverse mass for each of the species we are using
    Call atom_invmassMult(atypes1,atInvMass1)
    If (.Not.sameSpc) Call atom_invmassMult(atypes2,atInvMass2)

    !** Get the charges on molec1's atoms
    Call molecules_getcharges(spc1,charges1)

    !** Covert spc2 to a 1D array of atom vectors. The contents of atvecs2
    !** will look like (/ (atom1,mol1),(atom2,mol1),...,(atom1,mol2),... /)
    !** such that atvecs2(1) corresponds to atom1 of molecule1, atvecs2(2)
    !** corresponds to atom2 of molecule 1, and atvecs2((y-1)*natoms+x) 
    !** corresponds to atom x of molecule y, etc.
    atvecs2(1:natvecs2) = Reshape( &
        species(spc2)%coords(1:natoms2,1:nmoles2)%r, (/natvecs2/))
        
    !** Get the atom charges for spc2 and reshape into an array. Its 
    !** entries are in an order identical to that of atvecs2.
    Call molecules_getcharges(spc2,clists)
    charges2(1:natvecs2) = Reshape(clists(1:natoms2),(/natvecs2/), &
        clists(1:natoms2))

    !** Loop through all the molecules of spc1
    Do m = 1, nmoles1

      !** Loop through all the atoms of spc1
      Do a = 1, natoms1

        !** Check to see if this is the same spc
        If (sameSpc) Then

          !** Check to see if we should be calculating intramolecular
          !** potentials, too
          If (params%intraOn) Then

            !** Calculate the interaction for just molecule m
            !** Take atom a out of the atvecs2 list

            !** Check to see if we are dealing with a core-shell interaction
            !** If we are, then store the partner in csAtom. If isCore is
            !** .True., then atom a is a core. If isCore is .False., then
            !** atom a is a shell. If isCoreShell is .False., it is not
            !** a core/shell atom.
            isCoreShell = molecules_isCoreShell(spc1,a,isCore,csAtom)
            If (.Not.isCoreShell) csAtom = a

            !** We will be skipping at least one, possibly two atoms
            !** in our summation. Figure out the entry position in the
            !** atvecs2 array. Remember, the atvecs2 array is 1D.
            skipMe(1) = (m-1)*natoms2+a
            nskip = 1
            If (isCoreShell) Then
              skipMe(2) = (m-1)*natoms2+csAtom
              nskip = nskip + 1
            End If

            !** Call ewald to get the potential and the forces
            If (getf) Then

              !** Zero the potential and forces.
              utemp = 0.0_RDbl
              ftemp = VecType(0.0_RDbl)
              ftemp2 = VecType(0.0_RDbl)

              !** Call ewald to get the interaction between atom1 and 
              !** all of atvecs2. It returns the potential and forces.
              !** We want to include the self-term in the calculation,
              !** so we pass the selfTerm flag as .True.
              Call ewald_getinteractionSngl(sparams, &
                  species(spc1)%coords(a,m)%r, atvecs2, charges1(a), &
                  charges2, skipMe(1:nskip), .True., utemp, ftemp, ftemp2)

              !** We have two choices. Since we are calculating the potential
              !** and force on each atom, we have double counted potentials
              !** and must multiply by 1/2. For forces, we may either only
              !** update the acceleration for atom a, or we may update
              !** both a and all other atoms but multiply by 1/2.

              !** Update the acceleration on a only. We skip the other atoms
              !** since we will calculate their force later. No need for 
              !** 1/2 term. Scale the forces to accelerations.
!LC              aspc1(a,m) = aspc1(a,m) + ftemp*atInvMass1(a)*scalef

              !** Now get the core-shell correction, if necessary. Since we
              !** call it for both shells and cores, we need only to add
              !** the force to atom a.
              If (isCoreShell) Then
                Call ewald_getCoreShellCorr(sparams, &
                    species(spc1)%coords(a,m)%r, &
                    species(spc2)%coords(csAtom,m)%r, charges1(a), &
                    charges2(csAtom),ucorr,fcorr)
                utemp = utemp - ucorr
!LC                aspc1(a,m) = aspc1(a,m) - fcorr*scalef*atInvMass1(a)
                ftemp = ftemp - fcorr
              End If

!              Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!              Write(*,*) spc1,m,a,'       ',spc2
!              Write(*,*) utemp*0.5_RDbl
!              Write(*,*) ftemp

              !** Store it appropriately                 
              Call storesym_inc(ffout,(/spc1,m,a/),(/spc2,0,0/), &
                  utemp*0.5_RDbl,ftemp)

            Else  !** potential only
              !** Call ewald to get the total potential.
              Call ewald_getinteractionSngl(sparams, &
                  species(spc1)%coords(a,m)%r,atvecs2, &
                  charges1(a),charges2,skipMe(1:nskip),.True.,utemp)

              !** Get the core-shell correction if necessary
              If (isCoreShell) Then
                Call ewald_getCoreShellCorr(sparams, &
                    species(spc1)%coords(a,m)%r, &
                    species(spc2)%coords(csAtom,m)%r, charges1(a), &
                    clists(csAtom),ucorr)
                utemp = utemp - ucorr
              End If

              !** Store it appropriately                  
              Call storesym_inc(ffout,(/spc1,m,a/),(/spc2,0,0/), &
                  utemp*0.5_RDbl,dummy)
            End If

            !** Update the potential
!LC            pot = utemp + pot

          Else If (.Not.params%intraOn) Then
            !** We are not calculating intramolecular electrostatics.
            !** We must remove the molecule which contains the atom
            !** we are calculating the electrostatics for from the list
            nskip = natoms1

!            Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!            Stop
!FIXME           

              Do jskipatom=1, natoms1
              skipMe(jskipatom) = (m-1)*natoms1+ jskipatom  !
              End Do

            If (getf) Then
              !** Call ewald to get the potential and forces. Even though
              !** we will not be calculating intramolecular electrostatics,
              !** we still need to take into account the self term.
              Call ewald_getinteractionSngl(sparams, &
                  species(spc1)%coords(a,m)%r,atvecs2, &
                  charges1(a),charges2,skipMe(1:natoms1),.True.,utemp,ftemp,ftemp2)

              !** Update the forces
!LC              aspc1(a,m) = aspc1(a,m)+ftemp*scalef*atInvMass1(a)

              !** store it appropriately                  
              Call storesym_inc(ffout,(/spc1,m,a/),(/spc2,0,0/), &
                  utemp*0.5_RDbl,ftemp)

            Else

                  Call ewald_getinteractionSngl(sparams, &
                  species(spc1)%coords(a,m)%r,atvecs2(1:natoms2*(m-1)), &
                  charges1(a),charges2(1:natoms2*(m-1)),skipMe(1:natoms1), &
                  .True.,utemp)

              !** store it appropriately
              Call storesym_inc(ffout,(/spc1,m,a/),(/spc2,0,0/), &
                  utemp*0.5_RDbl,dummy)
            End If

            !** Update the potential
!LC            pot = pot + utemp

          End If   !** If (params%intraOn) Then

        Else    !** If (sameSpc) Then

          !** We are not calculating the interaction among a single
          !** spcate type. Therefore, we need to include both the
          !** species.

          !** WARNING: THIS DOES NOT WORK IF BOTH MOLECULE TYPES ARE 
          !** MOBILE

          !** Zero the potentials and forces
          ftemp = 0.0_RDbl
          utemp = 0.0_RDbl

          !** Check to see if we are calculating for fixed molecule
          !** types only in the ewald sum. If so, the non-fixed molecules
          !** have been left out of the recipricol space sum so we
          !** must not include the self term in the ewald calculation.
          selfTerm = .True.
          If (params%fixedOnly) Then
            If (.Not.config_isfixed(species(spc1))) selfTerm = .False.
          Else
            Write(0,'(a,i5,a)') __FILE__,__LINE__, &
                ": WARNING - Ewald will NOT work for two NONFIXED molecule &
                &types"
            Stop
          End If

          !** We aren't skipping anything, though nskip can't be 0
          skipMe = 0
          nskip = 1

          If (getf) Then
            Call ewald_getinteractionSngl(sparams, &
                species(spc1)%coords(a,m)%r,atvecs2(1:natvecs2), &
                charges1(a),charges2(1:natvecs2),skipMe(1:nskip), &
                selfTerm,utemp,ftemp,ftemp2)

            If (.Not. config_isfixed(species(spc2))) Then
              !** Update the forces
!LC              aspc1(a,m) = aspc1(a,m) + 0.5_RDbl*ftemp*scalef*atInvMass1(a)
              Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
              Write(*,*) 'not prepared to update spc2 forces'
              Stop
              Do aa = 1, natoms2
                Do mm = 1, nmoles2
!LC                  aspc2(aa,mm) = aspc2(aa,mm) - &
!LC                      0.5_RDbl*ftemp*scalef*atInvMass2(aa)
                End Do
              End Do

            Else
              !** We won't be double counting since one type is fixed, 
              !** so we don't need to multiply by 1/2.
!LC              aspc1(a,m) = aspc1(a,m) + ftemp*scalef*atInvMass1(a)

              !** store it appropriately
              Call storesym_inc(ffout,(/spc1,m,a/),(/spc2,0,0/),utemp,ftemp)
            End If

          Else   !** potential only
            Call ewald_getinteractionSngl(sparams, &
                species(spc1)%coords(a,m)%r,atvecs2(1:natvecs2), &
                charges1(a),charges2(1:natvecs2),skipMe(1:nskip), &
                selfTerm,utemp)

            !** store it appropriately
            Call storesym_inc(ffout,(/spc1,m,a/),(/spc2,0,0/),utemp,dummy)
          End If

          !** Double the potential if one is fixed only. This is to correct
          !** for the 1/2 factor at the end of this subroutine.
!LC          If (params%fixedOnly) utemp = 2.0_RDbl*utemp
!LC          pot = pot + utemp

        End If   !** If (sameSpc) Then
      End Do   !** Do a = 1, natoms1
    End Do   !** Do m = 1, nmoles1

    !** Finally, for the ewald sum, we must multiply by 1/2 to get the 
    !** correct potential. The 1/2 arises from the double-counting of
    !** interactions.
!LC    pot = pot*0.5_RDbl    !accounted for elsewhere, see _inc's

    !** if we've made it this far, then we're ok
    sssum_ssint = .True.

  End Function sssum_ssint


  !----------------------------------------------------------------------------
  ! Atom Species Ewald summation
  !----------------------------------------------------------------------------
  Logical Function sssum_asintHOT(params,spc1,mol1,atm1,spc2,species,sparams,fast,ffout,hot)
    Type(SSSumParams), Intent(In)               :: params
    Integer, Intent(In)                          :: spc1,mol1,atm1,spc2
    Type(SimCell_Params), Intent(In)            :: sparams
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Logical, Intent(In)                         :: fast
    Type(EnergyPlus), Intent(InOut)      :: ffout
    Real(Kind=RDbl), Dimension(:), Optional :: hot

    Integer             :: natoms1, natoms2, a, nmoles1, nmoles2, m, i, j
    Integer             :: natypes1, natypes2, csAtom
    Integer             :: low, med, high,currentSize,requiredSize,error
    Integer             :: natvecs1, natvecs2, aa, mm
    Integer             :: level, nderivs
    Logical             :: sameSpc, getf
    Logical             :: iscoreshell, selfTerm, isCore
    Real(Kind=Rdbl)     :: utemp, amass, ucorr
    Integer, Dimension(Size(species(spc1)%coords,1))         :: atypes1
    Integer, Dimension(Size(species(spc2)%coords,1))         :: atypes2
    Real(Kind=RDbl), Dimension(Size(species(spc1)%coords,1)) :: charges1
    Real(Kind=RDbl), Dimension(Size(species(spc2)%coords,1)) :: charges2
    Real(Kind=RDbl), Dimension(Size(species(spc1)%coords,1)) :: atInvMass1
    Real(Kind=RDbl), Dimension(Size(species(spc2)%coords,1)) :: atInvMass2
    Type(VecType)                                            :: fcorr,ftemp
    Type(VecType), Dimension(Size(species(spc2)%coords,1)* &
        Size(species(spc2)%coords,2))                        :: atvecs2
    Type(VecType), Dimension(Size(atvecs2,1))                :: ftemp2
    Integer, Dimension(Size(atvecs2,1)) :: skipMe !** atom numbers to skip
    Integer                             :: nskip  !** molecule number to skip

    Type(VecType)      :: dummy
    Real(Kind=Rdbl), Dimension(1:4)  :: hess
    Logical                    :: ewaldflag

      If (.Not.Present(hot)) Then
       Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'This routine serves one purpose: to Calculate Higher Order Derivatives for map generation'
       Write(*,*) 'However hots are not found. Are you using this call for something else?'
       Stop
      End If

    ewaldflag = .False.

    !** this is a hack that keeps the storesym_inc call the same
    !** NAG compiler isn't flexible enough to handle different calls
    dummy = VecType(0.0_RDbl)

    !** Check size of temporary array for holding charges
    currentSize = Size(clists,1)
    requiredSize= config_totalvecs(species)

    !** Size the array that holds the charges if necessary
    If (currentSize < requiredSize) Then
      Deallocate(clists,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_dealloc")
      Allocate(clists(requiredSize+INCR_SIZE),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_alloc")
    End If
     
    !** Check if the species are the same.
    !** We need to check the atvecs and make sure that for
    !** spc1 == spc2, molec1 == molec2, the interaction
    !** between atom1 == atom2 is not calculated.
    sameSpc = .False.
    If (spc1 == spc2) Then
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'This routine is aimed for map generation only'
       Write(*,*) 'The probe and zeolite can not be the same species:   ', spc1, spc2
       Stop
      End If

    !** We need the atom masses in order to scale the forces. Get the
    !** atom types so we can look up the atom masses later.
    natypes1 = molecules_getnatomtypes(spc1,atypes1,.True.)
    natypes2 = molecules_getnatomtypes(spc2,atypes2,.True.)
    
    !** Get the number of atoms and molecules for each type
    natoms1 = molecules_getnatoms(spc1)
    natoms2 = molecules_getnatoms(spc2)
    nmoles1 = config_getnmoles(species,spc1)
    nmoles2 = config_getnmoles(species,spc2)

    !** Calculate the total number of entries in the temporary arrays,
    !** e.g., natoms1*nmoles1
    natvecs1 = natoms1*nmoles1
    natvecs2 = natoms2*nmoles2

    !** Get the inverse mass for each of the species we are using
    Call atom_invmassMult(atypes1,atInvMass1)
    If (.Not.sameSpc) Call atom_invmassMult(atypes2,atInvMass2)

    !** Get the charges on molec1's atoms
    Call molecules_getcharges(spc1,charges1)

    !** Covert spc2 to a 1D array of atom vectors. The contents of atvecs2
    !** will look like (/ (atom1,mol1),(atom2,mol1),...,(atom1,mol2),... /)
    !** such that atvecs2(1) corresponds to atom1 of molecule1, atvecs2(2)
    !** corresponds to atom2 of molecule 1, and atvecs2((y-1)*natoms+x) 
    !** corresponds to atom x of molecule y, etc.
    atvecs2(1:natvecs2) = Reshape( &
        species(spc2)%coords(1:natoms2,1:nmoles2)%r, (/natvecs2/))

    !** Get the atom charges for spc2 and reshape into an array. Its 
    !** entries are in an order identical to that of atvecs2.
    Call molecules_getcharges(spc2,clists)
    charges2(1:natvecs2) = Reshape(clists(1:natoms2),(/natvecs2/), &
        clists(1:natoms2))

    m = mol1
    a = atm1

          !** Zero the potentials
          utemp = 0.0_RDbl
          ftemp = VecType(0.0_RDbl)

            Call ewald_getInteractionHOT(sparams,species(spc1)%coords(a,m)%r, &
            atvecs2(1:natvecs2), charges1(a),charges2(1:natvecs2),utemp,hot,ewaldflag)

           If (.Not.(ewaldflag)) Then 
             sssum_asintHOT = .False.
           Return
           End If


            !** Updating energy
            Call storebase_inc(ffout,utemp)
              
            !** Updating forces
            ftemp%comp(1)=-hot(2)
            ftemp%comp(2)=-hot(3)
            ftemp%comp(3)=-hot(4)
         
           Call storebase_inc(ffout,ftemp)

            hess(1) = hot(5)
            hess(2) = hot(6)
            hess(3) = hot(7)
            hess(4) = hot(8)

           !** Updating hessian in storage
            Call storebase_inc(ffout,hess)

    sssum_asintHOT = .True.

  End Function sssum_asintHOT

  !----------------------------------------------------------------------------
  ! Calculates the potential and forces due to sorb1 and sorb2 interacting
  ! using a summation type of calculation. Returns the accelerations on 
  ! each of the sorbates.
  !
  ! Some notes: The way in which the code is currently set up, we are actually
  ! double counting potentials for the Ewald summation because we loop over
  ! all atoms in both sorbates, so we must multiply the FINAL potential by
  ! 1/2. This would be true of the forces as well IF we added atvecs2 back
  ! into asorb2. We are not doing that, however, so we don't need the 1/2.
  ! ERASE THIS SOON
  !----------------------------------------------------------------------------
  Subroutine sssum_getssint(params,sorbates,sparams,sorb1,sorb2,&
      fast,pot,asorb1,asorb2)
    Type(SSSumParams), Intent(In)           :: params
    Type(SimCell_Params), Intent(In)        :: sparams
    Type(AtMolCoords), Dimension(:), Intent(In) :: sorbates
    Integer, Intent(In)                     :: sorb1, sorb2
    Logical, Intent(In)                     :: fast
    Real(Kind=RDbl), Intent(Out)            :: pot
    Type(VecType), Dimension(:,:), Intent(InOut), Optional :: asorb1, asorb2
    Integer :: natoms1, natoms2, a, nmoles1, nmoles2, m, i, j
    Real(Kind=RDbl), Dimension(Size(sorbates(sorb1)%coords,1)) :: charges1
    Real(Kind=RDbl), Dimension(Size(sorbates(sorb2)%coords,1)) :: charges2
    Integer, Dimension(Size(sorbates(sorb1)%coords,1)) :: atypes1
    Integer, Dimension(Size(sorbates(sorb2)%coords,1)) :: atypes2
    Real(Kind=RDbl), Dimension(Size(sorbates(sorb1)%coords,1)) :: atInvMass1
    Real(Kind=RDbl), Dimension(Size(sorbates(sorb2)%coords,1)) :: atInvMass2

    Integer :: natypes1, natypes2, csAtom
    Type(VecType), Dimension(Size(sorbates(sorb2)%coords,1)* &
        Size(sorbates(sorb2)%coords,2)) :: atvecs2
    Type(VecType) :: ftemp
    Type(VecType), Dimension(Size(atvecs2,1)) :: ftemp2
    Real(Kind=Rdbl) :: utemp, amass, ucorr
    Logical :: sameSorb, getf, iscoreshell, selfTerm, isCore
    Type(VecType) :: fcorr
    Integer :: low, med, high,currentSize,requiredSize,error
    Integer :: natvecs1, natvecs2
    Integer :: aa, mm

    !** Molecules and atoms to skip in our calculation
    Integer, Dimension(Size(atvecs2,1)) :: skipMe ! atom numbers to skip
    Integer :: nskip  ! molecule number to skip

    !** Set the atom and molecules to skip to zero
    skipMe = 0
    nskip = 0

    !** Check to see if the temporary array for holding charges
    !** is of correct size
    currentSize = Size(clists,1)
    requiredSize= config_totalvecs(sorbates)

    If (currentSize < requiredSize) Then
      Deallocate(clists,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_dealloc")
      Allocate(clists(requiredSize+INCR_SIZE),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_alloc")
    End If
    
    !** Check to see if we are calculating forces
    getf = .False.
    If (Present(asorb1)) getf = .True.

    !** zero the potential and the forces
    pot = zero
    ftemp = VecType(zero)

    !** Check if the sorbates are the same.
    !** We need to check the atvecs and make sure that for
    !** sorb1 == sorb2, molec1 == molec2, the interaction
    !** between atom1 == atom2 is not calculated.
    sameSorb = .False.
    If (sorb1 == sorb2) sameSorb = .True.

    !** We need the atom masses in order to scale the forces. Get the
    !** atom types so we can look up the atom masses later.
    natypes1 = molecules_getnatomtypes(sorb1,atypes1,.True.)
    natypes2 = molecules_getnatomtypes(sorb2,atypes2,.True.)
    
    !** Get the number of atoms and molecules for each type
    natoms1 = molecules_getnatoms(sorb1)
    natoms2 = molecules_getnatoms(sorb2)
    nmoles1 = config_getnmoles(sorbates,sorb1)
    nmoles2 = config_getnmoles(sorbates,sorb2)

    !** Calculate the total number of entries in the temporary arrays,
    !** e.g., natoms1*nmoles1
    natvecs1 = natoms1*nmoles1
    natvecs2 = natoms2*nmoles2

    !** Get the inverse mass for each of the species we are using
    Call atom_invmassMult(atypes1,atInvMass1)
    If (.Not.sameSorb) Call atom_invmassMult(atypes2,atInvMass2)

    !** Get the charges on molec1's atoms
    Call molecules_getcharges(sorb1,charges1)

    !** Covert sorb2 to a 1D array of atom vectors. The contents of atvecs2
    !** will look like (/ (atom1,mol1),(atom2,mol1),...,(atom1,mol2),... /)
    !** such that atvecs2(1) corresponds to atom1 of molecule1, atvecs2(2)
    !** corresponds to atom2 of molecule 1, and atvecs2((y-1)*natoms+x) 
    !** corresponds to atom x of molecule y, etc.
    atvecs2(1:natvecs2) = Reshape( &
        sorbates(sorb2)%coords(1:natoms2,1:nmoles2)%r, (/natvecs2/))

    !** Get the atom charges for sorb2 and reshape into an array. Its 
    !** entries are in an order identical to that of atvecs2.
    Call molecules_getcharges(sorb2,clists)
    charges2(1:natvecs2) = Reshape(clists(1:natoms2),(/natvecs2/), &
        clists(1:natoms2))

    !** Loop through all the molecules of sorb1
    Do m = 1, nmoles1

      !** Loop through all the atoms of sorb1
      Do a = 1, natoms1
        !** Check to see if this is the same sorb
        If (sameSorb) Then
          !** Check to see if we should be calculating intramolecular
          !** potentials, too
          If (params%intraOn) Then

            !** Calculate the interaction for just molecule m
            !** Take atom a out of the atvecs2 list

            !** Check to see if we are dealing with a core-shell interaction
            !** If we are, then store the partner in csAtom. If isCore is
            !** .True., then atom a is a core. If isCore is .False., then
            !** atom a is a shell. If isCoreShell is .False., it is not
            !** a core/shell atom.
            isCoreShell = molecules_isCoreShell(sorb1,a,isCore,csAtom)
            If (.Not.isCoreShell) csAtom = a

            !** We will be skipping at least one, possibly two atoms
            !** in our summation. Figure out the entry position in the
            !** atvecs2 array. Remember, the atvecs2 array is 1D.
            skipMe(1) = (m-1)*natoms2+a
            nskip = 1
            If (isCoreShell) Then
              skipMe(2) = (m-1)*natoms2+csAtom
              nskip = nskip + 1
            End If

            !** Call ewald to get the potential and the forces
            If (getf) Then

              !** Zero the potential and forces.
              utemp = 0.0_RDbl
              ftemp = VecType(0.0_RDbl)
              ftemp2 = VecType(0.0_RDbl)
              
              !** Call ewald to get the interaction between atom1 and 
              !** all of atvecs2. It returns the potential and forces.
              !** We want to include the self-term in the calculation,
              !** so we pass the selfTerm flag as .True.
              Call ewald_getinteractionSngl(sparams, &
                  sorbates(sorb1)%coords(a,m)%r, atvecs2, charges1(a), &
                  charges2, skipMe(1:nskip), .True., utemp, ftemp, ftemp2)

              !** We have two choices. Since we are calculating the potential
              !** and force on each atom, we have double counted potentials
              !** and must multiply by 1/2. For forces, we may either only
              !** update the acceleration for atom a, or we may update
              !** both a and all other atoms but multiply by 1/2.

              !** Update the acceleration on a only. We skip the other atoms
              !** since we will calculate their force later. No need for 
              !** 1/2 term. Scale the forces to accelerations.
              asorb1(a,m) = asorb1(a,m) + ftemp*atInvMass1(a)*scalef

!!$              !** Update the forces by adding ftemp and reshaping ftemp2
!!$              !** Remember, these are the same sorbates, so we must add
!!$              !** the forces to ONLY asorb1! We add the 1/2 term here.
!!$              asorb1(a,m) = asorb1(a,m) + 0.5_RDbl*ftemp*atInvMass1(a)*scalef
!!$              Do aa = 1, natoms1
!!$                Do mm = 1, nmoles1
!!$                  If (findint(skipMe(1:nskip),i) /= 0) Cycle
!!$                  asorb1(aa,mm) = asorb1(aa,mm) + &
!!$                      0.5_RDbl*ftemp2((mm-1)*natoms1+aa)*scalef*atInvMass1(aa)
!!$                End Do
!!$              End Do

              !** Now get the core-shell correction, if necessary. Since we
              !** call it for both shells and cores, we need only to add
              !** the force to atom a.
              If (isCoreShell) Then
                Call ewald_getCoreShellCorr(sparams, &
                    sorbates(sorb1)%coords(a,m)%r, &
                    sorbates(sorb2)%coords(csAtom,m)%r, charges1(a), &
                    charges2(csAtom),ucorr,fcorr)
                utemp = utemp - ucorr
                asorb1(a,m) = asorb1(a,m) - fcorr*scalef*atInvMass1(a)
!!$                asorb1(csAtom,m) = asorb1(csAtom,m) + &
!!$                    0.5_RDbl*fcorr*scalef*atInvMass1(csAtom)
              End If

            Else
              !** Call ewald to get the total potential.
              Call ewald_getinteractionSngl(sparams, &
                  sorbates(sorb1)%coords(a,m)%r,atvecs2, &
                  charges1(a),charges2,skipMe(1:nskip),.True.,utemp)
              !** Get the core-shell correction if necessary
              If (isCoreShell) Then
                Call ewald_getCoreShellCorr(sparams, &
                    sorbates(sorb1)%coords(a,m)%r, &
                    sorbates(sorb2)%coords(csAtom,m)%r, charges1(a), &
                    clists(csAtom),ucorr)
                utemp = utemp - ucorr
              End If
            End If

            !** Update the potential
            pot = utemp + pot

          Else If (.Not.params%intraOn) Then
            !** We are not calculating intramolecular electrostatics.
            !** We must remove the molecule which contains the atom
            !** we are calculating the electrostatics for from the
            !** list
            nskip = natoms1
            Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
            Stop
!FIXME            skipMe = (m-1)*natoms1+(/1:natoms1/)

            If (getf) Then
              !** Call ewald to get the potential and forces. Even though
              !** we will not be calculating intramolecular electrostatics,
              !** we still need to take into account the self term.
              Call ewald_getinteractionSngl(sparams, &
                  sorbates(sorb1)%coords(a,m)%r,atvecs2, &
                  charges1(a),charges2,skipMe,.True.,utemp,ftemp,ftemp2)

              !** Update the forces
              asorb1(a,m) = asorb1(a,m)+ftemp*scalef*atInvMass1(a)
!!$              Do aa = 1, natoms1
!!$                Do mm = 1, nmoles1
!!$                  asorb1(aa,mm) = asorb1(aa,mm) + &
!!$                      0.5_RDbl*ftemp2((mm-1)*natoms1+aa)*scalef*atInvMass1(aa)
!!$                End Do
!!$              End Do

            Else
              Call ewald_getinteractionSngl(sparams, &
                  sorbates(sorb1)%coords(a,m)%r,atvecs2(1:natoms2*(m-1)), &
                  charges1(a),charges2(1:natoms2*(m-1)),skipMe(1:nskip), &
                  .True.,utemp)
            End If

            !** Update the potential
            pot = pot + utemp

          End If ! If (params%intraOn) Then

        Else ! If (sameSorb) Then

          !** We are not calculating the interaction among a single
          !** sorbate type. Therefore, we need to include both the
          !** sorbates.

          !** WARNING: THIS DOES NOT WORK IF BOTH MOLECULE TYPES ARE 
          !** MOBILE

          !** Zero the potentials and forces
          ftemp = 0.0_RDbl
          utemp = 0.0_RDbl

          !** Check to see if we are calculating for fixed molecule
          !** types only in the ewald sum. If so, the non-fixed molecules
          !** have been left out of the recipricol space sum so we
          !** must not include the self term in the ewald calculation.
          selfTerm = .True.
          If (params%fixedOnly) Then
            If (.Not.config_isfixed(sorbates(sorb1))) selfTerm = .False.
          Else
            Write(0,'(a,i5,a)') __FILE__,__LINE__, &
                ": WARNING - Ewald will NOT work for two NONFIXED molecule &
                &types"
            Stop
          End If

          !** We aren't skipping anything, though nskip can't be 0
          skipMe = 0
          nskip = 1

          If (getf) Then
            Call ewald_getinteractionSngl(sparams, &
                sorbates(sorb1)%coords(a,m)%r,atvecs2(1:natvecs2), &
                charges1(a),charges2(1:natvecs2),skipMe(1:nskip), &
                selfTerm,utemp,ftemp,ftemp2)
            If (.Not.config_isfixed(sorbates(sorb2))) Then
              !** Update the forces
              asorb1(a,m) = asorb1(a,m) + 0.5_RDbl*ftemp*scalef*atInvMass1(a)
              Do aa = 1, natoms2
                Do mm = 1, nmoles2
                  asorb2(aa,mm) = asorb2(aa,mm) - &
                      0.5_RDbl*ftemp*scalef*atInvMass2(aa)
                End Do
              End Do
            Else
              !** We won't be double counting since one type is fixed, 
              !** so we don't need to multiply by 1/2.
              asorb1(a,m) = asorb1(a,m) + ftemp*scalef*atInvMass1(a)
            End If
          Else 
            Call ewald_getinteractionSngl(sparams, &
                sorbates(sorb1)%coords(a,m)%r,atvecs2(1:natvecs2), &
                charges1(a),charges2(1:natvecs2),skipMe(1:nskip), &
                selfTerm,utemp)
          End If
          !** Double the potential if one is fixed only. This is to correct
          !** for the 1/2 factor at the end of this subroutine.
          If (params%fixedOnly) utemp = 2.0_RDbl*utemp
          pot = pot + utemp
        End If ! If (sameSorb) Then

      End Do ! Do a = 1, natoms1

    End Do ! Do m = 1, nmoles1

    !** Finally, for the ewald sum, we must multiply by 1/2 to get the 
    !** correct potential. The 1/2 arises from the double-counting of
    !** interactions.
    pot = pot*0.5_RDbl

  End Subroutine sssum_getssint

  !----------------------------------------------------------------------------
  ! Calculates the potential and forces due to sorb1 and sorb2 interacting
  ! using a summation type of calculation. Returns ONLY the pot and forces
  ! on SORB1, MOLEC1 due to SORB2.
  !----------------------------------------------------------------------------
  Subroutine sssum_getmsint(params,sorbates,sparams,sorb1,molec1,sorb2,&
      fast,pot,asorb1,asorb2)
    Type(SSSumParams), Intent(In)        :: params
    Type(SimCell_Params), Intent(In) :: sparams
    Type(AtMolCoords), Dimension(:), Intent(In) :: sorbates
    Integer, Intent(In)                  :: sorb1, sorb2, molec1
    Logical, Intent(In)                  :: fast
    Real(Kind=RDbl), Intent(Out)         :: pot
    Type(VecType), Dimension(:,:), Intent(InOut), Optional :: asorb1,asorb2
    Integer :: natoms1, natoms2, a, nmoles2,error, currentVecSize
    Real(Kind=RDbl), Dimension(Size(sorbates(sorb1)%coords,1)) :: charges1
    Real(Kind=RDbl), Dimension(Size(sorbates(sorb2)%coords,1)) :: charges2
    Type(VecType), Dimension(Size(sorbates(sorb2)%coords,1)* &
        Size(sorbates(sorb2)%coords,2)) &
        :: atvecs2
    Type(VecType) :: fvec
    Type(VecType), Dimension(Size(atvecs2,1)) :: fvec2
    Real(Kind=Rdbl) :: utemp

    !** zero the potential
    pot = zero
    
    !** Check to see if the temeprory array clists is of correct size
    currentVecSize=config_totalvecs(sorbates)
    If (Size(clists,1)<currentVecSize) Then
      Deallocate(clists,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_dealloc")
      Allocate(clists(currentVecSize),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_alloc")
    End If
    
    !MDEBUG 
    If (sorb1==sorb2) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          "Sorry, can't calculate for sorb1=sorb2 yet"
    End If

    !** Get the number of atoms and molecules
    natoms1 = molecules_getnatoms(sorb1)
    natoms2 = molecules_getnatoms(sorb2)
    nmoles2 = config_getnmoles(sorbates,sorb2)

    !** Get the charges on molec1's atoms
    Call molecules_getcharges(sorb1,charges1)

    !** Covert sorb2 to an array of atom vectors
    atvecs2 = Reshape(sorbates(sorb2)%coords(1:natoms2,1:nmoles2)%r, &
        (/natoms2*nmoles2/))

    !** Get the atom charges for sorb2 and reshape into an array
    Call molecules_getcharges(sorb2,clists)
    charges2 = Reshape(clists(1:natoms2),(/natoms2*nmoles2/),clists(1:natoms2))

    !** Loop through all the atoms of sorb1
    Do a = 1, natoms1
      If (Associated(params%ewald)) Then 
        If (Present(asorb1)) Then
          Call ewald_getinteractionSngl(sparams, &
              sorbates(sorb1)%coords(a,molec1)%r, &
              atvecs2,charges1(a),charges2,(/0/),.True.,utemp,fvec)
        Else
          Call ewald_getinteractionSngl(sparams, &
              sorbates(sorb1)%coords(a,molec1)%r, &
              atvecs2,charges1(a),charges2,(/0/),.True.,utemp)
        End If
      End If
      pot = utemp + pot
    End Do

  End Subroutine sssum_getmsint

  !----------------------------------------------------------------------------
  ! Calculates the potential and forces due to sorb1 and sorb2 interacting
  ! using a summation type of calculation. Returns ONLY the pot and forces
  ! on SORB1, MOLEC1 due to SORB2.
  !----------------------------------------------------------------------------
  Subroutine sssum_getmsintHOT(params,sorbates,sparams,sorb1,molec1, &
       sorb2,fast,pot,hot)
    Type(SSSumParams), Intent(In) :: params
    Type(AtMolCoords), Dimension(:), Intent(In) :: sorbates
    Type(SimCell_Params), Intent(In) :: sparams
    Integer, Intent(In) :: sorb1, sorb2, molec1
    Logical, Intent(In) :: fast
    Real(Kind=RDbl), Intent(Out) :: pot
    Real(Kind=RDbl), Dimension(:), Intent(Out) :: hot
    Integer :: natoms1, natoms2, a, nmoles2, error, currentVecSize
    Real(Kind=RDbl), Dimension(Size(sorbates(sorb1)%coords,1)) :: charges1
    Real(Kind=RDbl), Dimension(Size(sorbates(sorb2)%coords,1)) :: charges2
    Type(VecType), Dimension(Size(sorbates(sorb2)%coords,1)* &
        Size(sorbates(sorb2)%coords,2)) &
        :: atvecs2
    Real(Kind=Rdbl) :: utemp

    
    !** Check to see if the temeprory array clists is of correct size
    currentVecSize=config_totalvecs(sorbates)
    If (Size(clists,1)<currentVecSize) Then
      Deallocate(clists,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_dealloc")
      Allocate(clists(currentVecSize),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"clists_alloc")
    Endif
    
!!$    !MDEBUG
!!$    Write(*,*) "Arrived in sssum_msgetintHOT"

    !** zero the potential
    pot = zero

    !MDEBUG 
    If (sorb1==sorb2) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          "Sorry, can't calculate for sorb1=sorb2 yet"
    End If

    !** Get the number of atoms and molecules
    natoms1 = molecules_getnatoms(sorb1)
    natoms2 = molecules_getnatoms(sorb2)
    nmoles2 = config_getnmoles(sorbates,sorb2)

    !** Get the charges on molec1's atoms
    Call molecules_getcharges(sorb1,charges1)

    !** Covert sorb2 to an array of atom vectors
    atvecs2 = Reshape(sorbates(sorb2)%coords(1:natoms2,1:nmoles2)%r, &
        (/natoms2*nmoles2/))

    !** Get the atom charges for sorb2 and reshape into an array
    Call molecules_getcharges(sorb2,clists)
    charges2 = Reshape(clists(1:natoms2),(/natoms2*nmoles2/),clists(1:natoms2))

    !** Loop through all the atoms of sorb1
    Do a = 1, natoms1
      If (Associated(params%ewald)) Then 
        Call ewald_getinteractionHOT(sparams,sorbates(sorb1)%coords(a,molec1)%r, &
            atvecs2,charges1(a),charges2,utemp,hot)
      End If
      pot = utemp + pot
    End Do

  End Subroutine sssum_getmsintHOT

  !----------------------------------------------------------------------------
  ! Returns the intraOn flag which indicates if the interaction will include
  ! interactions with the molecules for like species pairs.
  ! Requires:  params -- summation parameters
  !----------------------------------------------------------------------------
  Logical Function sssum_intraon(params)
    Type(SSSumParams), Intent(In)               :: params

    sssum_intraon = params%intraon

  End Function sssum_intraon

  !----------------------------------------------------------------------------
  ! Displays information about the sum parameters
  ! Requires: params -- SS sum spc-spc parameters
  !           indent -- no. of spaces from the left margin
  !           unit -- output unit number
  !----------------------------------------------------------------------------
  Subroutine sssum_display(params,indent,unit)
    Type(SSSumParams), Intent(In)        :: params
    Integer, Intent(In)                  :: indent
    Integer, Optional, Intent(In)        :: unit

    Integer                     :: i
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    !** Report the model name and the other relevant parameters
    Write(unit,'(3a)') blank,"Sum Model Name : ",params%model
    If (params%fast) Then
      Write(unit,'(2a)') blank,"Fast           : Yes"
    Else
      Write(unit,'(2a)') blank,"Fast           : No"
    End If
    If (params%fixedOnly) Then
      Write(unit,'(2a)') blank,"Fixed Only     : Yes"
    Else
      Write(unit,'(2a)') blank,"Fixed Only     : No"
    End If
    If (params%intraOn) Then
      Write(unit,'(2a)') blank,"Include Intramolecular Interactions : Yes"
    Else
      Write(unit,'(2a)') blank,"Include Intramolecular Interactions : No"
    End If

    !** Now call the individual models to get the info from them
    If (Associated(params%ewald)) Then
      Call ewald_display(unit,indent+2)
    End If

  End Subroutine sssum_display

  !----------------------------------------------------------------------------
  ! Cleans the structure
  ! Requires: params -- SS sum spc-spc parameters
  !----------------------------------------------------------------------------
  Subroutine sssum_clean(params)
    Type(SSSumParams), Intent(InOut)        :: params

  End Subroutine sssum_clean

End Module sssum


