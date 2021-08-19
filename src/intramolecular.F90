!-----------------------------------------------------------------------------
! Divies out responsibility to calculate all intramolecular potentials. In all
! cases, the data structures for each interaction type contain a list of 
! interactions.  This list actually takes the form of two arrays, a 1D array
! of pointer sets and another, 2D array containing the atom numbers for the 
! interaction.  The pointer sets contain a list of pointers to various potential
! types.  For example, the pointer set for the intra-pair interactions contains
! a pointer each for the Lennard-Jones and the Buckingham potentials.  Those
! pointers that are initialized allow their corresponding potential interactions
! to be calculated.
!
! This module has facilities to handle 2-body (stretching, bs), 3-body 
! (bond bending, bb), 4-body (torsion, tor), constraint (con) and intra-pair
! interactions.  All new potential terms should be inserted into the routines 
! appropriate to the number of atoms (bodies) involved in the potential.  Note
! that when assigning potential energies to groups of atoms, the X-body 
! interaction is divided equally into X parts, one for each atom in the set.
!
! There are two potential calculation routines in this module:
!   intramolecular_spcint or intramolecular_molint
!     used to calculate the intramolecular interactions using the principal
!     coordinates.  Pass simcell in for molecules that fill the simulation 
!     cell like zeolites.  
!
! Note that this routine depends intimately on molecule, and directly accesses
! objects from molecule. This means changes in molecule requires changes in
! this routine.
!
! Needed Improvements:
! 1) more comments and requires statements
! 2) rename stretching, bending, torsion to 2,3,4-body to generalize
! 3) directly access molecule's data type -- bad OOP (see _init)
! 4) add another intramolecular_molint with EnergyPlus input for spc-spc storage
!-------------------------------------------------------------------------------

Module intramolecular

  Use defaults, Only: RDbl, strLen, lstrLen, zeroTolerance, degTorad, &
      scalef, scalepe, zero, d_ctrl_file, MAX_ATOMS, xlstrlen, &
      STRETCH_INDEX, BENDING_INDEX, TORSION_INDEX, &
      CONSTRAINT_INDEX, INTRAPAIR_INDEX, INTRACOUL_INDEX, &
      TOTAL_INDEX, NO_OF_INTRA_POTS
  Use utils, Only: split, stripcmnt, findint, tolower, getpath, filesrchstr, &
      toint, isfileopen, toupper, toreal, combine, allocErrDisplay, &
      deallocErrDisplay
  Use file, Only: file_getunit, file_getname, file_open
  Use atom, Only: atom_gettypename, atom_getmass, atom_nbonds, &
      atom_getntypes, atom_getsymbol, atom_getname, atom_invmass
  Use molecule, Only: molecule_getnthatom, MolecularParams, &
      molecule_getnatomtypes,molecule_getatomcoords
  Use molecules, Only: molecules_getpointer, molecules_getnatoms, &
      molecules_name, molecules_getfilename, molecules_getmasses, &
      molecules_getdof, molecules_changedof, molecules_dofIsOK, &
      molecules_getnatomtypes
  Use connects, Only: MAX_CONNECTIONS,connects_makeTypeList, &
      connects_makeXList
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag, vector_display, vector_getnormsq, &
      vector_getcomp, vector_iscollinear, vector_crossprod, vector_angle, &
      vector_getcom, vector_getnorm
  Use matrix, Only: MatrixType, matrix_getinv, Assignment(=), Operator(*)
  Use stats, Only: Statistics, stats_init, stats_update
  Use config, Only: AtMolCoords, config_genfromparent, config_getnmoles
  Use simcell, Only: SimCell_Params, simcell_minimage, simcell_getmolectype
  Use bsmodel, Only: StretchInfo, bsmodel_bsinit, bsmodel_bsint, &
      bsmodel_bsdisplay, bsmodel_bscleanup,bsmodel_bsinfo, isfast, &
      bsmodel_eval, STRETCH_KEY, bsmodel_setHarK
  Use bbmodel, Only: BendingInfo, bbmodel_bbinit, bbmodel_getparams, isfast, &
      bbmodel_bbdisplay, bbmodel_bbint, bbmodel_bbcleanup, bbmodel_bbinfo, &
      bbmodel_eval, BENDING_KEY, bbmodel_setHarK
  Use tormodel, Only: TorsionInfo, tormodel_torinit, tormodel_tordisplay, &
      tormodel_torint, tormodel_torinfo, tormodel_torcleanup, isfast, &
      tormodel_eval, TORSION_KEY, tormodel_setA1A2
  Use ipmodel, Only: IntrapairModel, IntrapairInfo, ipmodel_ipinit, &
      ipmodel_ipdisplay, ipmodel_displayCutoffs, ipmodel_copy, isfast, &
      ipmodel_ipint, ipmodel_ipcleanup, ipmodel_isinit, ipmodel_eval, &
      ipmodel_getcutoff, ipmodel_ipinfo, INTRAPAIR_KEY
  Use icmodel, Only: IntracoulModel, IntracoulInfo, icmodel_icinit, &
      icmodel_icdisplay, icmodel_displayCutoffs, icmodel_copy, isfast, &
      icmodel_icint, icmodel_iccleanup, icmodel_isinit, icmodel_eval, &
      icmodel_getcutoff, icmodel_icinfo, INTRACOUL_KEY
  Use conmodel, Only: ConstraintModel, ConstraintInfo, conmodel_getnconstr, &
      conmodel_initarray, conmodel_coninit, conmodel_getconint, isfast, &
      conmodel_getpenalty, conmodel_concleanup, conmodel_condisplay, &
      conmodel_postadjust, conmodel_getdofred, conmodel_eval, CONSTRAINT_KEY
  Use xternal, Only: ExternalInfo, EXTERNAL_KEY, xternal_init, &
      xternal_intraint, xternal_display
  Use ciccotti, Only: CiccottiRigidSetInfo
  Use storetop, Only: storetop_level
  Use storebase, Only: storebase_disp,storebase_init
  Use store, Only: Store_Level_Pair,store_init,store_terminate,store_display

  Implicit None
  Save

  Private
  Public :: IntramolecularInfo, intramolecular_init, intramolecular_initstore, &
      intramolecular_sampleCF, intramolecular_clean, intramolecular_display, &
      intramolecular_molint, intramolecular_spcint, &
      intramolecular_getpenalty, intramolecular_mconint, &
      intramolecular_postadjust, intramolecular_changespeed, &
      intramolecular_intrainfo, intramolecular_hasint, &
      intramolecular_torIndex, STRETCH_KEY, BENDING_KEY, TORSION_KEY, &
      INTRAPAIR_KEY, INTRACOUL_KEY,CONSTRAINT_KEY, intramolecular_setparams

  Type IntramolecularInfo
    Integer                       :: nterms,spc
    Type(StretchInfo), Pointer    :: stretch    !* Stretch info (2-body)
    Type(BendingInfo), Pointer    :: bending    !* Bending info (3-body)
    Type(TorsionInfo), Pointer    :: torsion    !* Torsional info (4-body)
    Type(IntrapairInfo), Pointer  :: intrapair  !* 1-4 Interaction info
    Type(IntracoulInfo), Pointer  :: intracoul !* Internal Coulombic Info
    Type(ConstraintInfo), Pointer :: constraint !* bond-length constraints
    Type(CiccottiRigidSetInfo), Pointer :: rigidset  !* Rigid Body constraints
    Type(ExternalInfo), Pointer   :: xtrnal !* external interaction evaluations
  End Type IntramolecularInfo

  Character(len=strLen), Parameter  :: intramolecular_infotag = 'INTRA:'

Contains 

  !----------------------------------------------------------------------------
  ! Initialize the intramolecular interactions for a single species
  ! Requires:  ffparams -- intramolecular forcefield parameters
  !            spc -- the molecule type number
  !            intra_filename -- the filename to take info from
  !            simcell -- the simulation cell
  !----------------------------------------------------------------------------
  Subroutine intramolecular_init(ffparams,spc,intra_filename,simcell)
    Type(IntramolecularInfo), Intent(Out)       :: ffparams
    Integer, Intent(In)                         :: spc
    Character(*), Intent(In)                    :: intra_filename
    Type(SimCell_Params), Intent(In), Optional  :: simcell  

    Integer                                  :: i,error,m
    Integer                                  :: natoms,dof,nterms,nmoles
    Integer                                  :: unitno,lineno,nfields,nchunks
    Logical                                  :: skiperror
    Character(len=255)                       :: line
    Character(len=strLen)                    :: dof_origin
    Character(len=strLen)                    :: molec_name,molec_filename
    Character(len=lstrLen)                   :: subline
    Character(len=5*strLen)                  :: molec_filewPath
    Character(len=strLen), Dimension(2)      :: srchstr
    Character(len=strLen), Dimension(20)     :: fields,chunks
    Type(MolecularParams), Pointer           :: molecule

    !** Nullify the pointer set
    Call intramolecular_nullify(ffparams)

    skiperror = .False.
    ffparams%spc = spc
    natoms = molecules_getnatoms(spc)
    molec_name = molecules_name(spc)
    molec_filename = molecules_getfilename(spc)

    !** the music philosophy is that molecules are stored in MOLSDIR directory
    molec_filename = molecules_getfilename(spc)
    molec_filewPath = Trim(getpath('MOLSDIR'))//Trim(molec_filename)

    Write(*,'(2a)') 'Initializing intra molecular info from &
        &molecule file: '//Trim(molec_filewPath)
    Call molecules_getpointer(spc,molecule)

    !** Open the file containing the intramolecular interaction information
    unitno = file_open(intra_filename,110)

    !** find the molecule name and intrainfo tag in the file
    srchstr = (/molec_name,intramolecular_infotag/)
    lineno = filesrchstr(unitno,srchstr,line,.True.)
    If (lineno == 0) Then
      Write(0,'(2a,i4,6a)') __FILE__,": ",__LINE__, &
          ' Could not find search strings: ',Trim(srchstr(1)), ', ',&
          Trim(srchstr(2)), ' in file ', Trim(intra_filename)
      Stop
    End If
    Close(unit=unitno)

    !** warn if necessary
    nfields = split(line,fields)
    If ((Index(Toupper(line),STRETCH_KEY) /= 0).And. &
        (Index(Toupper(line),CONSTRAINT_KEY)/= 0)) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ': code is not currently equipped to handle both constraints', &
          ' and bond stretching.'
      Stop
    End If

    !** initialize the interactions
    nterms = 0
    nfields = split(line,fields)
    Do i = 3,nfields
      nchunks = split(fields(i),chunks,'@')
      Select Case(toupper(chunks(1)))

      Case(STRETCH_KEY)
        If (natoms > 1) Then
          !** Allocate the stretch pointer
          Allocate(ffparams%stretch,STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bs pointer')
          Call bsmodel_bsinit(ffparams%stretch,molecule,molec_filewPath)
          nterms = nterms + 1

          !MDEBUG
          !LC          Write(0,'(2a)') "Initialized Bond Stretch for ",Trim(molec_name)
        End If

      Case(BENDING_KEY)
        If (natoms > 2) Then
          !** Allocate the bending pointer
          Allocate(ffparams%bending,STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bb pointer')
          Call bbmodel_bbinit(ffparams%bending,spc,molec_filewPath)
          nterms = nterms + 1

          !MDEBUG
          !LC          Write(*,'(2a)') "Initialized Bond Bending for ",Trim(molec_name)
        End If

      Case(TORSION_KEY)
        If (natoms > 3) Then
          !** Allocate the torsion pointer
          Allocate(ffparams%torsion,STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'torsion pointer')
          Call tormodel_torinit(ffparams%torsion,spc,molec_filewPath)
          nterms = nterms + 1

          !MDEBUG
          !LC          Write(*,'(2a)') "Initialized Torsion for ",Trim(molec_name)
        End If

      Case(INTRAPAIR_KEY)
        If (natoms > 4) Then
          !** Allocate the intrapair pointer
          Allocate(ffparams%intrapair,STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ip pointer')
          nterms = nterms + 1

          !** Check to see if we need the simcell
          If (Present(simcell)) Then
            If (spc == simcell_getmolectype(simcell)) Then
              Call ipmodel_ipinit(ffparams%intrapair,molecule, &
                  molec_filewPath,simcell)
            Else
              Call ipmodel_ipinit(ffparams%intrapair,molecule,molec_filewPath)
            End If
          Else
            Call ipmodel_ipinit(ffparams%intrapair,molecule,molec_filewPath)
          End If
        End If

      Case(INTRACOUL_KEY)
        If (natoms > 1) Then
          !** Allocate the intrapair pointer
          Allocate(ffparams%intracoul,STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ic pointer')
          nterms = nterms + 1

          !** Check to see if we need the simcell
          If (Present(simcell)) Then
            If (spc == simcell_getmolectype(simcell)) Then
              Call icmodel_icinit(ffparams%intracoul,molecule, &
                  molec_filewPath,simcell)
            Else
            Call icmodel_icinit(ffparams%intracoul,molecule,molec_filewPath)
            End If
          Else
            Call icmodel_icinit(ffparams%intracoul,molecule,molec_filewPath)
          End If


          !MDEBUG
          !LC          Write(*,'(2a)') "Initialized IntraCoulombic for ",Trim(molec_name)
        End If


      Case(CONSTRAINT_KEY)
        If (natoms > 1) Then
          !** Allocate the constraint pointer
          Allocate(ffparams%constraint,STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'con pointer')
          Call conmodel_coninit(ffparams%constraint,molecule,molec_filewPath)
          nterms = nterms + 1

          !** Modify the DOF if appropriate 
          dof = molecules_getdof(spc,dof_origin)
          If (Trim(Toupper(dof_origin)) /= 'MOLECULE_FILE') Then
            dof = dof - conmodel_getdofred(ffparams%constraint%conparams)
            Call molecules_changedof(spc,dof,'Modifed_3N_Calculation')
            If (.Not.molecules_dofIsOK(spc)) Then
              Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
              stop
            Endif
            Write(*,'(3a,i4)') __FILE__," Changed DOF ",Trim(molec_name),dof
          End If

          !MDEBUG
!LC          Write(*,'(2a)') "Initialized contraints for ",Trim(molec_name)
        End If

      Case(EXTERNAL_KEY)
        skiperror = .True.
        !** Allocate pointer
        Allocate(ffparams%xtrnal,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'external pointer')

        !** could pass simcell in here and carry pointer in %xtrnal
        subline = combine(fields(i:))
        If (Present(simcell)) Then
          Call xternal_init(ffparams%xtrnal,subline,(/spc/),simcell)
        Else
          Call xternal_init(ffparams%xtrnal,subline,(/spc/))
        End If
        nterms = nterms + 1

      Case Default
        If (.Not. skiperror) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              ': could not interpret keyword ',Trim(fields(i))
          Write(0,'(a,5(1x,a))') 'Currently understood keys: ', &
              Trim(STRETCH_KEY),Trim(BENDING_KEY),Trim(TORSION_KEY), &
              Trim(INTRAPAIR_KEY),Trim(INTRACOUL_KEY),Trim(CONSTRAINT_KEY), &
              Trim(EXTERNAL_KEY)
          Stop
        End If

      End Select
    End Do

    ffparams%nterms = nterms

  End Subroutine intramolecular_init

  !----------------------------------------------------------------------------
  ! Nullify the intramolecular pointers
  ! Requires:  ffparams -- intramolecular forcefield parameters
  !----------------------------------------------------------------------------
  Subroutine intramolecular_nullify(ffparams)
    Type(IntramolecularInfo), Intent(InOut)  :: ffparams

    Nullify(ffparams%stretch)
    Nullify(ffparams%bending)
    Nullify(ffparams%torsion)
    Nullify(ffparams%intrapair)
    Nullify(ffparams%intracoul)
    Nullify(ffparams%constraint)
    Nullify(ffparams%rigidset)
    Nullify(ffparams%xtrnal)

  End Subroutine intramolecular_nullify

  !----------------------------------------------------------------------------
  ! Initialize the intramolecular storage for a single species
  ! Requires:  ffparams -- intramolecular forcefield parameters
  !            storage -- forcefield results storage structure
  !            spc -- species number
  !            nderivs -- number of derivatives in storage
  !            nmoles -- number of molecules to initialize 
  !            storelevel -- integer indicating desired storage level
  !----------------------------------------------------------------------------
  Subroutine intramolecular_initstore(ffparams,storage,spc,nderivs, &
      nmoles,storelevel)
    Type(IntramolecularInfo), Intent(In)   :: ffparams
    Type(Store_Level_Pair), Intent(Out)    :: storage 
    Integer, Intent(In)                    :: spc,nderivs,nmoles,storelevel

    Integer           :: m,natoms

    natoms = molecules_getnatoms(spc)

    !** initialize the storage structure
    If (ffparams%nterms > 0) Then

      !** initialize molecule-based storage or terminate
      If (storelevel > 2) Then
        Call store_init(storage,(/nmoles/),10,nderivs,.True.)
      Else
        Call store_terminate(storage,nderivs,.True.)        
      End If

      !** initialize atom-based storage if desired
      If (storelevel > 3) Then
        Do m = 1,nmoles
          Call store_init(storage%mi(m),(/natoms/),10,nderivs,.True.)
        End Do

      Else If (storelevel == 3) Then
        !** initialize condensed atom storage if needed for forces
        If (nderivs > 0) Then
          Do m = 1,nmoles
            Call store_init(storage%mi(m),(/natoms/),0,nderivs,.True.)
          End Do
        End If

      End If

    Else
      !** No storage needed, terminate at the species level
      Call store_terminate(storage,0,.True.)

    End If

  End Subroutine intramolecular_initstore

  !----------------------------------------------------------------------------
  ! Get the intramolecular interactions for a SPECIES.  If the simulation
  ! cell information is also passed, then it will be used to calculate minimum
  ! image distances.
  ! Requires:  ffparams -- intramolecular forcefield parameters
  !            ffout -- forcefield output 
  !            spc -- species number
  !            species -- the coordinates for all the molecules of type spc
  !            fast -- a logical flag that indicates if
  !                   fast or slow interactions are being calculated
  !            simcell -- optional simulation cell information
  !----------------------------------------------------------------------------
  Subroutine intramolecular_spcint(ffparams,ffout,spc,species,fast,simcell)
    Type(IntramolecularInfo), Intent(InOut)      :: ffparams
    Type(Store_Level_Pair), Intent(InOut)        :: ffout
    Integer, Intent(In)                          :: spc
    Type(AtMolCoords), Intent(In)                :: species
    Logical, Intent(In)                          :: fast
    Type(SimCell_Params), Intent(In), Optional   :: simcell

    Integer           :: m,nmoles,natoms
    Logical           :: success

    If (ffparams%nterms == 0) Return

    !** Initialize some important numbers
    nmoles = config_getnmoles(species)
    natoms = molecules_getnatoms(spc)

    Do m = 1,nmoles
      If (Present(simcell)) Then
        success = intramolecular_molint(ffparams,ffout%mi(m), &
            species%coords(1:natoms,m)%rp,fast,simcell)
      Else
        success = intramolecular_molint(ffparams,ffout%mi(m), &
            species%coords(1:natoms,m)%rp,fast)
      End If
      If (.Not. success) Then
        Write(0,'(2a,i5,a,2i3)') __FILE__,":",__LINE__, &
            " Intramolecular evaluation unsuccessful (spc,mol): ",spc,m 
        Stop
      End If
    End Do

  End Subroutine intramolecular_spcint

  !----------------------------------------------------------------------------
  ! Get the intramolecular interactions for a MOLECULE.  If the simulation
  ! cell information is also passed, then it will be used to calculate minimum
  ! image distances.  Note that this routine does not include constraint
  ! interactions.
  ! Requires:  ffparams -- intramolecular forcefield parameters
  !            ffout -- forcefield output 
  !            coords -- coordinates for whole molecule indexed by atm num
  !            fast -- a logical flag that indicates if
  !                   fast or slow interactions are being calculated
  !            simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Logical Function intramolecular_molint(ffparams,ffout,coords,fast,simcell)
    Type(IntramolecularInfo), Intent(InOut)      :: ffparams
    Type(Store_Level_Pair), Intent(InOut)        :: ffout
    Type(VecType), Dimension(:), Intent(In)      :: coords
    Logical, Intent(In)                          :: fast
    Type(SimCell_Params), Intent(In), Optional   :: simcell

    !** Set default to True since most cases don't have failure feedback
    intramolecular_molint = .True.

    !** Check if simcell is present then evaluate individual terms
    If (Present(simcell)) Then
      If (bsmodel_eval(ffparams%stretch,fast)) Then
        Call bsmodel_bsint(ffparams%stretch,ffout,coords,simcell)
      End If

      If (bbmodel_eval(ffparams%bending,fast)) Then
        Call bbmodel_bbint(ffparams%bending,ffout,coords,simcell)
      End If

      If (tormodel_eval(ffparams%torsion,fast)) Then
          Call tormodel_torint(ffparams%torsion,ffout,coords,simcell)
      End If

      If (ipmodel_eval(ffparams%intrapair,fast)) Then
        Call ipmodel_ipint(ffparams%intrapair,ffout,coords,simcell)
      End If

      If (icmodel_eval(ffparams%intracoul,fast)) Then
        Call icmodel_icint(ffparams%intracoul,ffout,coords,simcell)
      End If

      If (Associated(ffparams%xtrnal)) Then
        intramolecular_molint = xternal_intraint(ffparams%xtrnal, &
            ffout,coords,ffparams%spc,simcell)
      End If

    Else   !** simcell not present, don't use PBCs

      If (bsmodel_eval(ffparams%stretch,fast)) Then
        Call bsmodel_bsint(ffparams%stretch,ffout,coords)
!        Call store_disp(ffout,.False.,2,6)
!        stop
      End If

      If (bbmodel_eval(ffparams%bending,fast)) Then
        Call bbmodel_bbint(ffparams%bending,ffout,coords)
      End If

      If (tormodel_eval(ffparams%torsion,fast)) Then
        Call tormodel_torint(ffparams%torsion,ffout,coords)
      End If

      If (ipmodel_eval(ffparams%intrapair,fast)) Then
        Call ipmodel_ipint(ffparams%intrapair,ffout,coords)
      End If

      If (icmodel_eval(ffparams%intracoul,fast)) Then
        Call icmodel_icint(ffparams%intracoul,ffout,coords)
      End If

      If (Associated(ffparams%xtrnal)) Then
        intramolecular_molint = xternal_intraint(ffparams%xtrnal, &
            ffout,coords,ffparams%spc)
      End If

    End If

  End Function intramolecular_molint

  !----------------------------------------------------------------------------
  ! Get the intramolecular CONSTRAINT interactions for a MOLECULE.  
  ! Requires:  ffparams -- intramolecular forcefield parameters
  !            coords -- coordinates for whole molecule indexed by atm num
  !            v -- velocities for whole molecule indexed by atm num
  !            u -- potential energy from constraints
  !            accels -- accelerations on each molecule
  !            fast -- a logical flag that indicates if
  !                   fast or slow interactions are being calculated
  !            all_flag -- flag indicating if all constraints evaluated
  !----------------------------------------------------------------------------
  Subroutine intramolecular_mconint(ffparams,coords,v,u,accels,all_flag)
    Type(IntramolecularInfo), Intent(In)       :: ffparams
    Type(VecType), Dimension(:), Intent(In)    :: coords,v
    Real(kind=RDbl), Intent(Out)               :: u
    Type(VecType), Dimension(:), Intent(InOut) :: accels
    Logical, Intent(In)                        :: all_flag

    Call conmodel_getconint(ffparams%constraint,coords,v,u,accels,all_flag)

  End Subroutine intramolecular_mconint

  !----------------------------------------------------------------------------
  ! Calculate the penalty for the bond length constraints
  ! Requires: ffparams -- intramolecular forcefield information
  !           coords -- atomic coordinates
  !           v -- atomic velocities
  !           sumobjf -- sum of object function
  !----------------------------------------------------------------------------
  Subroutine intramolecular_getpenalty(ffparams,coords,v,sumobjf)
    Type(IntramolecularInfo)                   :: ffparams
    Type(VecType), Dimension(:), Intent(InOut) :: coords, v
    Real(Kind=RDbl), Intent(Out)               :: sumobjf

    If (Associated(ffparams%constraint)) Then 
      Call conmodel_getpenalty(ffparams%constraint,coords,v,sumobjf)
    End If

  End Subroutine intramolecular_getpenalty

  !----------------------------------------------------------------------------
  ! Finds the Torsion Index, given the atom-numbers
  ! Requires: ptr -- pointer to torsion information 
  !           molecule -- pointer to molecule data structure
  !           atlist -- the index of the atoms
  ! This could be expensive: dont call frequently
  !----------------------------------------------------------------------------
  Integer Function intramolecular_torIndex(ptr,molecule,atlist)
    Type(TorsionInfo), Pointer       ::  ptr
    Type(MolecularParams), Pointer   :: molecule
    Integer, Dimension(4), Intent(in):: atlist 
    Integer :: i

    intramolecular_torindex=0

    Do i=1,ptr%ntorsion
      !** Check for matching with the torlist, with atoms in ascending order
      If ((atlist(1)==ptr%torlist(i,1)) .And.  &
          (atlist(2)==ptr%torlist(i,2)) .And.  &
          (atlist(3)==ptr%torlist(i,3)) .And.  &
          (atlist(4)==ptr%torlist(i,4))) Then
        intramolecular_torIndex=i
        Return
      Endif

      !** Check for matching with the torlist, with atoms in descending order
      If ((atlist(4)==ptr%torlist(i,1)) .And.  &
          (atlist(3)==ptr%torlist(i,2)) .And.  &
          (atlist(2)==ptr%torlist(i,3)) .And.  &
          (atlist(1)==ptr%torlist(i,4))) Then
        intramolecular_torIndex=i
        Return
      Endif

    End do
  End Function Intramolecular_torIndex

  !----------------------------------------------------------------------------
  ! Checks to see if the intramolecular interaction is turned on.  If a 
  ! specific identifier string is passed it will only check that interaction
  ! type.  Otherwise, it will determine if the ANY intramolecular interactions
  ! are turned on.
  ! Requires:  ffparams -- intramolecular forcefield information
  !            whichOne -- string identifier for model type
  !----------------------------------------------------------------------------
  Logical Function intramolecular_hasint(ffparams,whichOne)
    Type(IntramolecularInfo), Intent(In)  :: ffparams
    Character(*), Intent(In), Optional    :: whichOne

    intramolecular_hasint = .False.
    If (Present(whichOne)) Then
      Select Case (Trim(ToUpper(whichOne)))
      Case ('STRETCH')
        If (Associated(ffparams%stretch)) &
            intramolecular_hasint = .True.
      Case ('BENDING')
        If (Associated(ffparams%bending)) &
            intramolecular_hasint = .True.
      Case ('CONSTRAINT')
        If (Associated(ffparams%constraint)) &
            intramolecular_hasint = .True.
      Case ('TORSION')
        If (Associated(ffparams%torsion)) & 
            intramolecular_hasint = .True.
      Case ('INTRAPAIR')
        If (Associated(ffparams%intrapair)) &
            intramolecular_hasint = .True.
      Case ('INTRACOUL')
        If (Associated(ffparams%intracoul)) &
            intramolecular_hasint = .True.
      Case ('EXTERNAL')
        If (Associated(ffparams%xtrnal)) &
            intramolecular_hasint = .True.
      End Select
    Else
      If (ffparams%nterms > 0) intramolecular_hasint = .True.
    End If

  End Function intramolecular_hasint

  !----------------------------------------------------------------
  ! Cleanup intramolecular allocations
  ! Requires:  ffparams -- intramolecular parameters structure  
  !----------------------------------------------------------------
  Subroutine intramolecular_clean(ffparams)
    Type(IntramolecularInfo), Intent(InOut)    :: ffparams

    Integer       :: error

    !** Cleanup bond-stretching stuff
    If (Associated(ffparams%stretch)) Then
      Call bsmodel_bscleanup(ffparams%stretch)
      Deallocate(ffparams%stretch, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'stretch')
    End If

    !** Cleanup bond-bending stuff
    If (Associated(ffparams%bending)) Then
      Call bbmodel_bbcleanup(ffparams%bending)
      Deallocate(ffparams%bending, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'bending')
    End If

    !** Cleanup torsion stuff
    If (Associated(ffparams%torsion)) Then
      Call tormodel_torcleanup(ffparams%torsion)
      Deallocate(ffparams%torsion, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'torsion')
    End If

    !** Cleanup intrapair stuff
    If (Associated(ffparams%intrapair)) Then
      Call ipmodel_ipcleanup(ffparams%intrapair)
      Deallocate(ffparams%intrapair, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'intrapair')
    End If

    !** Cleanup intracoul stuff
    If (Associated(ffparams%intracoul)) Then
      Call icmodel_iccleanup(ffparams%intracoul)
      Deallocate(ffparams%intracoul, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'intracoul')
    End If


    !** Cleanup constraint stuff
    If (Associated(ffparams%constraint)) Then
      Call conmodel_concleanup(ffparams%constraint)
      Deallocate(ffparams%constraint, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'constraint')
    End If

  End Subroutine intramolecular_clean

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine intramolecular_sampleCF(unitno)
    Integer, Intent(In) :: unitno
    
    Write(unitno,'(2a,1x,5(3a))') intramolecular_infotag, ' Molecule_name',&
        ' [',STRETCH_KEY,']',' [',BENDING_KEY,']', &
        ' [',TORSION_KEY,']',' [',CONSTRAINT_KEY,']', &
        ' [',INTRAPAIR_KEY,']', ' [',INTRACOUL_KEY,']'
  End Subroutine intramolecular_sampleCF

  !----------------------------------------------------------------------------
  ! Changes all the interaction types to "fast" or "slow"
  !----------------------------------------------------------------------------
  Subroutine intramolecular_changeSpeed(ptr,fast)
    Type(IntramolecularInfo), Intent(InOut) :: ptr
    Logical, Intent(In) :: fast

    If (Associated(ptr%stretch)) ptr%stretch%fast = fast
    If (Associated(ptr%bending)) ptr%bending%fast = fast
    If (Associated(ptr%torsion)) ptr%torsion%fast = fast
    If (Associated(ptr%intrapair)) ptr%intrapair%fast = fast
    If (Associated(ptr%intracoul)) ptr%intracoul%fast = fast
    If (Associated(ptr%constraint)) ptr%constraint%fast = fast

  End Subroutine intramolecular_changeSpeed

  !----------------------------------------------------------------------------
  ! Resets information for a given species's specific intramolecular
  ! interaction information.  
  ! Requires:  iparams -- intramolecular forcefield information
  !            itype -- interaction type via keywords defined at top
  !            nsets -- number of potential sets returned
  !            list -- list of atom numbers for interaction (set,1:Natoms)
  !            params -- array of strings containing type and parameters
  ! NOTE : STRETCH has different return format!! need to be fixed at some point
  !----------------------------------------------------------------------------
  Subroutine intramolecular_setparams(iparams,description,alist,val)
    Type(IntramolecularInfo), Intent(InOut)                :: iparams
    Character(len=strlen),Dimension(:), Intent(In)    :: description
    Integer, Dimension(:), Intent(in)                 :: alist
    Real(kind=RDbl), Intent(in) :: val

    !** default nothing found
    Select Case (Trim(description(1)))
    Case (STRETCH_KEY)

      If (Associated(iparams%stretch)) Then
        If (Trim(description(2))=="KR") Then
          Call bsmodel_setHarK(iparams%stretch,alist(1:2),val)
        Else
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Case (BENDING_KEY)
      If (Associated(iparams%bending)) Then
        If (Trim(description(2))=="KTHETA") Then
          Call bbmodel_setHarK(iparams%bending,alist(1:3),val)
        Else
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Case (TORSION_KEY)
      If (Associated(iparams%torsion)) Then
        If (Trim(description(2))=="TORA") Then
          Call tormodel_setA1A2(iparams%torsion,alist(1:4),val)
        Else
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Case Default
      Write(0,'(2a,i6,a,i9,a)') __FILE__,":",__LINE__, &
          " Could not interpret potential key: ",Trim(description(1))
      Stop
    End Select

  End Subroutine intramolecular_setparams


  !----------------------------------------------------------------------------
  ! Returns information about a given species's specific intramolecular
  ! interaction information.  
  ! Requires:  ffparams -- intramolecular forcefield information
  !            itype -- interaction type via keywords defined at top
  !            nsets -- number of potential sets returned
  !            list -- list of atom numbers for interaction (set,1:Natoms)
  !            params -- array of strings containing type and parameters
  ! NOTE : STRETCH has different return format!! need to be fixed at some point
  !----------------------------------------------------------------------------
  Subroutine intramolecular_intrainfo(ffparams,itype,nsets,list,params)
    Type(IntramolecularInfo), Intent(In)                :: ffparams
    Character(*), Intent(In)                            :: itype
    Integer, Intent(Out)                                :: nsets
    Integer, Dimension(:,:), Intent(Out)                :: list
    Character(len=lstrLen), Dimension(:), Intent(Out)   :: params

    !** default nothing found
    nsets=0

    Select Case (Trim(itype))
    Case (STRETCH_KEY)
      If (Associated(ffparams%stretch)) &
          Call bsmodel_bsinfo(ffparams%stretch,nsets,list,params)
    Case (BENDING_KEY)
      If (Associated(ffparams%bending)) &
          Call bbmodel_bbinfo(ffparams%bending,nsets,list,params)
    Case (TORSION_KEY)
      If (Associated(ffparams%torsion)) &
          Call tormodel_torinfo(ffparams%torsion,nsets,list,params)
    Case (INTRAPAIR_KEY)
      If (Associated(ffparams%intrapair))&
          Call ipmodel_ipinfo(ffparams%intrapair,nsets,list,params)
    Case (INTRACOUL_KEY)
      If (Associated(ffparams%intracoul)) &
          Call icmodel_icinfo(ffparams%intracoul,nsets,list,params)
    Case (CONSTRAINT_KEY)
      !      Call intramolcular_coninfo(ffparams%constraint,nsets,list,params)
      Write(0,'(2a,i6,a,i9,a)') __FILE__,":",__LINE__, &
          " Constraint info routine does not yet exist, sorry"
      Stop
    Case Default
      Write(0,'(2a,i6,a,i9,a)') __FILE__,":",__LINE__, &
          " Could not interpret potential key: ",Trim(itype)
      Stop
    End Select

  End Subroutine intramolecular_intrainfo

  !----------------------------------------------------------------------------
  ! Do any necessary post-integration adjustments for a single species
  !----------------------------------------------------------------------------
  Subroutine intramolecular_postadjust(intrainfo,simcell,species)
    Type(IntramolecularInfo), Intent(In)       :: intrainfo
    Type(SimCell_Params), Intent(In)           :: simcell  
    Type(AtMolCoords), Intent(InOut)           :: species

    Integer                    :: m

    If (Associated(intrainfo%constraint)) Then 
      Do m = 1,Size(species%coords,2)
        Call conmodel_postadjust(intrainfo%constraint%conparams, &
            species%coords(:,m)%rp, &
            species%coords(:,m)%v)
      End Do
      Call config_genfromparent(simcell,species)
    End If

  End Subroutine intramolecular_postadjust

  !----------------------------------------------------------------------------
  ! Dump the fields of the Intramolecular Info structure
  ! Requires: ffparams -- intramolecule info structure for a single species
  !           spc -- species number
  !           indent -- no. of spaces from the left margin
  !           unit -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine intramolecular_display(ffparams,spc,indent,unit)
    Type(IntramolecularInfo), Intent(In)   :: ffparams
    Integer, Intent(In)                    :: spc,indent
    Integer, Optional, Intent(In)          :: unit

    Integer                        :: i
    Character(len=indent)          :: blank
    Character(len=strLen)          :: line

    blank = Repeat(' ',indent)
    
    !** Bond stretching info
    If (Associated(ffparams%stretch)) Then
      line = ""//" On "//Trim(ffparams%stretch%model)
      If (ffparams%stretch%fast) Write(line,'(2a)') Trim(line)," Fast "
      Write(unit,'(4a,i5,a)') blank,'Bond Stretching   : ',Trim(line), &
          '  ',ffparams%stretch%nstretch,' stretch pairs'
      If (ffparams%stretch%nstretch < 100) Then
        Call bsmodel_bsdisplay(ffparams%stretch,spc,unit,indent+2)
      Else
        Write(unit,'(2x,2a)') blank,'nstretch > 100, not showing stretch pairs'
      End If
    Else
      Write(unit,'(2a)') blank,'Bond Stretching   :  OFF'
    End If

    !** Bond bending information
    If (Associated(ffparams%bending)) Then
      line = ""//" On "//trim(ffparams%bending%model)
      If (ffparams%bending%fast) Write(line,'(2a)') Trim(line)," Fast "
      Write(unit,'(3a,i5,a)') blank,'Bond Bending      : ',Trim(line), &
          ffparams%bending%nbends,' bending triplets'
      If (ffparams%bending%nbends < 100) Then
        Call bbmodel_bbdisplay(ffparams%bending,spc,unit,indent+2)
      Else
        Write(unit,'(2x,3a)') blank,'nbends > 100,', &
            ' not showing bond bending triplets'
      End If
    Else
      Write(unit,'(2a)') blank,'Bond Bending      :  OFF'
    End If

    !** Torsion angle information
    If (Associated(ffparams%torsion)) Then
      line = ""//" On "
      If (ffparams%torsion%fast) Write(line,'(2a)') Trim(line)," Fast "
      Write(unit,'(3a,i5,a)') blank,'Four-Body Interactions: ',Trim(line), &
          ffparams%torsion%ntorsion,' four-body atom sets'
      If (ffparams%torsion%ntorsion < 100) Then
        Call tormodel_tordisplay(ffparams%torsion,spc,indent+2,unit)
      Else
        Write(unit,'(2x,3a)') blank,'ntorsion > 100,', &
            ' not showing torsion quadruplets'
      End If
    Else
      Write(unit,'(2a)') blank,'Torsion Angles    :  OFF'
    End If

    !** Report the intrapair interaction details
    If (Associated(ffparams%intrapair)) Then
      !** basics
      line = ""//" On "//Trim(ffparams%intrapair%model)
      If (ffparams%intrapair%fast) line = Trim(line)//" Fast "
      Write(unit,'(3a,i7,a)') blank,'Intra Pairwise Pot: ',Trim(line), &
          ffparams%intrapair%nipairs,' intra pairs'
      !** more detail
      Call ipmodel_ipdisplay(ffparams%intrapair,spc,indent+2,unit)

    Else
      Write(unit,'(3a)') blank,'Intra Pairwise Pot:  OFF'
    End If

    !** Report the intracoul interaction details
    If (Associated(ffparams%intracoul)) Then
      !** basics
      line = ""//" On "//Trim(ffparams%intracoul%model)
      If (ffparams%intracoul%fast) line = Trim(line)//" Fast "
      Write(unit,'(3a,i7,a)') blank,'Intra Coulombic Pot: ',Trim(line), &
          ffparams%intracoul%nipairs,' intra pairs of charges'

      !** more detail
      Call icmodel_icdisplay(ffparams%intracoul,spc,indent+2,unit)
    Else
      Write(unit,'(3a)') blank,'Intra Coulombic Pot:  OFF'
    End If

    !** Report the constraint details
    If (Associated(ffparams%constraint)) Then
      Write(unit,'(2a,l1)') blank,'Bond Length Constr: ', &
          Associated(ffparams%constraint)
      Call conmodel_condisplay(ffparams%constraint, &
          spc,unit,indent+2)
    Else
      Write(unit,'(2a)')  blank,'Constraints       :  OFF'
    End If

    !** Report usage of external evaluation if they are present
    If (Associated(ffparams%xtrnal)) Then
      Call xternal_display(ffparams%xtrnal,indent+2,unit)
    End If

  End Subroutine intramolecular_display

End Module intramolecular

#ifdef MUST_BE_REWRITTEN

  !--------------------------------------------------------------------
  ! There are parts of the code that need to know what parameters are used 
  ! for bondbending, given the atomnumber this routine returns the parameters
  !--------------------------------------------------------------------
  Subroutine intramolecular_getbbparams(model,species,atomlist,params)
    Character(len=strlen),Intent(out) :: model
    Integer, Intent(in)                :: species
    Integer,Dimension(3), Intent(in)   :: atomlist
    Real(kind=RDbl), Dimension(:), Intent(out)  :: params
    Integer :: i, paramno
    paramno = -1

    ! Note : here atomlist contains the atom-numbers and not atomtypes
    Do i=1,intraint(species)%bending%nbends
      !check the middel atom first
      If (atomlist(2) == intraint(species)%bending%bblist(i,2)) Then
        
        If (atomlist(1) == intraint(species)%bending%bblist(i,1)) Then
          
          If (atomlist(3) == intraint(species)%bending%bblist(i,3)) Then
            paramno=i
            Exit
          Endif
          
          ! the atomlist could  match in the reverse order too
        Elseif (atomlist(3) == intraint(species)%bending%bblist(i,1)) Then
          
          If (atomlist(1) == intraint(species)%bending%bblist(i,3)) Then
            paramno=i
            Exit
          Endif
      
        Endif

      Endif

    End Do
    If (paramno==-1) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    ! gets the params value from bbmodel
    Call bbmodel_getparams(intraint(species)%bending%bbparams(paramno), &
        model,params)
  End Subroutine intramolecular_getbbparams

#endif








