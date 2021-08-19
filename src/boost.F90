!------------------------------------------------------------------------------
! This modules contains the routines and datatypes for boosting forces 
! and potentials, as applied to accelerated molecular dynamics. For 
! references, see:
!
! 1. Voter, A. "A method for accelerating the molecular dynamics simulation
!    of infrequent events," J. Chem. Phys., 106, 4665-4676 (1997).
! 2. Voter, A. "Hyperdynamics: Accelerated Molecular Dynamics of Infrequent
!    Events," Phys. Rev. Lett., 78, 3908-3911 (1997).
! 3. Wang, J., Somnath, P., and Fichthorn, K. "Accelerated molecular 
!    dynamics of rare events using the local boost method," Phys. Rev. B,
!    63, 085403 (2001).
! 4. Steiner, M. and Genilloud, P., "Simple bias potential for boosting
!    molecular dynamics with the hyperdynamics scheme," Phys. Rev. B, 57,
!    10236-10239 (1998).
! 5. Pal, S. and Fichthorn, K., "Accelerated molecular dynamics of rare 
!    events using the local boost method," Phys. Rev. B, 74, 77-83,
!    (1999).
!
! The idea is normal potentials and forces are passed to this module. The
! routines in this module will search for molecules to boost, apply the
! boost, and returned the modified forces, potentials, and time scales.
! For MD, that means that after the forces are obtained, they are passed
! to this routine before the integration is finished such that new
! positions and velocities are based on these modified forces.
!
! For more information on implementation, force derivatives, etc., see
! Sanborn, Martin, "On the Study of Diffusion in Zeolites Using Molecular
!    Dynamics Techniques," Ph.D. thesis, Northwestern University, 2002.
!
! Needed Improvements
! 1) boosted is declared in this module, should be external
!------------------------------------------------------------------------------
Module boost

  Use defaults, Only: RDbl, strLen, dashedline, scalef, STATS_BLOCKSIZE, &
      scalepe
  Use utils, Only: stripcmnt, split, toreal, filesrchstr
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(*), &
      Operator(/), Operator(-), mag
  Use molecules, Only: molecules_getnsorbs, molecules_gettype, &
      molecules_getnatomtypes, molecules_name, molecules_getnatoms
  Use config, Only: AtMolCoords,config_isfixed,config_getnmoles,config_dumpmol
  Use atom, Only: atom_getmass, atom_invmass
  Use stats, Only: Statistics, stats_init, stats_update, stats_getvalue, &
      stats_getcavg, stats_getblock, stats_getstd

  Implicit None
  Save

  Private
  Public :: BoostParams, BoostType, StBoostVmax, StBoostVmin, StBoost, &
      boost_init, boost_calcBoost, boost_sampleCF, boost_display, &
      boost_displayStats, boost_getBoostMolname, boost_calcBoostHack

  !** References indicate the type of smoothing used.
  Type BoostType
    Type(StBoostVmax), Pointer :: bVmax     ! Reference 3
    Type(StBoostVmin), Pointer :: bVmin     ! Reference 3
    Type(StBoostVmin), Pointer :: bVminMol  ! Reference 3
    Type(StBoost), Pointer     :: bV        ! Reference 4
  End Type BoostType

  Type BoostParams
    Logical, Dimension(:), Pointer :: isBoost ! indx by mtype, True if boosted
    Real(Kind=RDbl)       :: vb               ! boosting threshold potential
    Type(BoostType)       :: type             ! Type of boost potential
    Character(len=strLen) :: molecName        ! name of the molecule to boost
    Type(Statistics)      :: bstats           ! Pontential boost stats
    Type(Statistics)      :: tstats           ! Time step stats
    Type(Statistics)      :: dvstats          ! Potential difference stats
    Type(Statistics)      :: sstats           ! Smoothing function stats
    Type(Statistics)      :: drstats          ! derivative value stats
  End Type BoostParams

  !** This is based on reference (3) 
  Type StBoostVmax
    Real(Kind=RDbl) :: vb     ! boosting threshold potential
    Real(kind=RDbl) :: gamma    ! smoothing parameter, has energy units
    Real(Kind=RDbl) :: c        ! smoothing parameter, dimensionless
  End Type StBoostVmax

  !** This is based on reference (3)
  Type StBoostVmin
    Real(Kind=RDbl), Dimension(:), Pointer :: vb
    Integer         :: nvb
    Real(kind=RDbl) :: n        ! smoothing parameter, dimensionless
    Real(Kind=RDbl) :: c        ! smoothing parameter, weird units
  End Type StBoostVmin

  !** This is based on reference (4)
  Type StBoost
    Real(Kind=RDbl) :: vb
    Real(Kind=RDbl) :: alpha
    Real(Kind=RDbl) :: gamma
  End Type StBoost

  Character(len=strLen), Parameter :: vbUnits = "kJ/mol"

  Type(Statistics) :: boostStat

  !** Make an interface for the calcboost command.
  Interface boost_calcBoost
    Module Procedure boost_calcBoostSingle
    Module Procedure boost_calcBoostMult
  End Interface

Contains

  !----------------------------------------------------------------------------
  ! Reads in the required information from the control file. It assumes it
  ! must find the boosted MD section in the given file
  ! Requires:  params -- boost parameters
  !            species -- full system coordinate structure
  !            unitno -- unit to read input parameters from
  !            isBoost -- False => just nullify and leave, True => initialize
  !----------------------------------------------------------------------------
  Subroutine boost_init(params,species,unitno,isBoost)
    Type(BoostParams), Intent(InOut)            :: params
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: unitno
    Logical, Intent(In), Optional               :: isBoost

    Integer                                  :: error, mtype, i
    Integer                                  :: nfields, lineno
    Logical                                  :: boost
    Character(len=strLen)                    :: text
    Character(len=strLen), Dimension(strLen) :: fields

    !** Check the boost flag. If it's false, we just want to nullify
    !** stuff and leave
    boost = .True.
    If (Present(isBoost)) boost = isBoost

    !** Initialize and nullify the boost type pointers
    Call boost_nullify(params)

    !** Size the boost type array
    Allocate(params%isBoost(molecules_getnsorbs()), stat = error)
    If (error /= 0) Then
      Write(0,'(2a,i6,a,i6)') __FILE__,":",__LINE__, &
          " Could not allocate boosted%isBoost of size ",molecules_getnsorbs()
      Stop
    End If

    !** Set all entries to false
    params%isBoost = .False.

    !** Leave if we aren't boosting
    If (.Not. boost) Return

    !** Read in the molecule name to boost
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    nfields = split(text, fields)
    Do i = 1, nfields
      mtype = molecules_gettype(Trim(fields(i)))
      !** Check to make sure that molecule isn't fixed
      If (config_isfixed(species(mtype))) Then
        Write(0,'(2a,i6,2a)') __FILE__,":",__LINE__, &
            " Ignoring request to boost molecule with FIXED config ", &
            Trim(fields(i))
      Else 
        params%isBoost(mtype) = .True.
        params%molecName = Trim(fields(i))
      End If
    End Do

    !** Read in the type of boost to be used
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    nfields = split(text,fields)
    Select Case (fields(1))

    Case ('STATICVMAX')
      !** initialize the pointer
      Allocate(params%type%bVmax,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            " Could not allocate pointer params%type%bVmax"
        Stop
      End If

      !** Read in the boosting threshhold in kJ/mol, convert to kcal/mol
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bVmax%vb = toreal(fields(1))/scalepe

      !** Read in the gamma parameter
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bVmax%gamma = toreal(fields(1))

      !** Read in the c parameter
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bVmax%c = toreal(fields(1))

    Case ('STATICVMIN')
      !** Warn the users about the molecule size
      Write(0,*) __FILE__,__LINE__, &
          " Warning: this won't necessarily work for molecules with &
          & natoms > 1!"

      !** Initialize the pointer
      Allocate(params%type%bVmin,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            " Could not allocate pointer params%type%bVmin"
        Stop
      End If

      !** Read in all the boost potentials
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields,',')
      params%type%bVmin%nvb = nfields
      
      !** Allocate the array to hold the boost potentials in
      Allocate(params%type%bVmin%vb(nfields),stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i6,a,i3)') __FILE__,":",__LINE__, &
            " Could not allocate params%type%bVmin%vb of size ",nfields
        Stop
      End If

      !** Store the potentials in kcal/mol (read in as kJ/mol)
      Do i = 1, nfields
        params%type%bVmin%vb(i) = toreal(fields(i))/scalepe
      End Do

      !** Read in the smoothing parameter n
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bVmin%n = toreal(fields(1))

      !** Read in the smoothing parameter c
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bVmin%c = toreal(fields(1))

    Case ('STATICVMINMOL')

      !** Initialize the pointer
      Allocate(params%type%bVminMol,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            " Could not allocate pointer params%type%bVminMol"
        Stop
      End If

      !** Read in all the boost potentials
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields,',')
      params%type%bVminMol%nvb = nfields
      
      !** Allocate the array to hold the boost potentials in
      Allocate(params%type%bVminMol%vb(nfields),stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i6,a,i3)') __FILE__,":",__LINE__, &
            " Could not allocate params%type%bVminMol%vb of size ",nfields
        Stop
      End If

      !** Store the potentials in kcal/mol (read in as kJ/mol)
      Do i = 1, nfields
        params%type%bVminMol%vb(i) = toreal(fields(i))/scalepe
      End Do

      !** Read in the smoothing parameter n
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bVminMol%n = toreal(fields(1))

      !** Read in the smoothing parameter c
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bVminMol%c = toreal(fields(1))

    Case('STATIC')
      !** initialize the pointer
      Allocate(params%type%bV,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            " Could not allocate pointer params%type%bV"
        Stop
      End If
      !** Read in the boosting threshhold, convert from kJ/mol to kcal/mol
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bV%vb = toreal(fields(1))/scalepe
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bV%alpha = toreal(fields(1))
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      nfields = split(text,fields)
      params%type%bV%gamma = toreal(fields(1))

    Case Default
      Write(0,'(2a,i6,2a)') __FILE__,":",__LINE__, &
          " Could not recoginize boost type ",trim(fields(1))
      Stop
    End Select

    !** Initialize the boost statistics
    Call stats_init(params%bstats,'Boost Statistics',STATS_BLOCKSIZE,.False.)
    Call stats_init(params%tstats,'Timestep Statistics',STATS_BLOCKSIZE,.False.)
    Call stats_init(params%dvstats,'Delta V Statistics',STATS_BLOCKSIZE,.False.)
    Call stats_init(params%sstats,'Smoothing Statistics',STATS_BLOCKSIZE,.False.)
    Call stats_init(params%drstats,'Derivative Statistics', &
        STATS_BLOCKSIZE,.False.)

  End Subroutine boost_init

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file input to the given unitno
  !----------------------------------------------------------------------------
  Subroutine boost_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a)') dashedline
    Write(unitno,'(a,t30,a)') 'Character','# Molecule to apply boost to'
    Write(unitno,'(a,t30,a)') 'Real','# Boosting threshhold, kcal'
    Write(unitno,'(a,t30,a)') '[STATICVMAX,STATICVMIN]', &
        '# Type of boost and smoothing to use'
    Write(unitno,'(a)') &
        '# IF STATICVMAX, delete this line and include the following: '
    Write(unitno,'(a,t30,a)') 'Real','# Boosting threshhold, kcal'
    Write(unitno,'(a,t30,a)') 'Real','# Value of gamma, a smoothing parameter'
    Write(unitno,'(a,t30,a)') 'Real','# Value of c, a smoothing parameter'
    Write(unitno,'(a)') '# END IF STATICVMAX'
    Write(unitno,'(a)') &
        '# IF STATICVMIN, delete this line and include the following: '
    Write(unitno,'(a,t30,a)') 'Real,Real,Real,...', &
        '# List of boost threshholds, kcal'
    Write(unitno,'(a,t30,a)') 'Real','# Value of n, a smoothing parameter'
    Write(unitno,'(a,t30,a)') 'Real','# Value of c, a smoothing parameter'
    Write(unitno,'(a)') '# END IF STATICVMIN'
    
  End Subroutine boost_sampleCF

  !----------------------------------------------------------------------------
  ! Allocates and nullifies all the boost type pointers
  !----------------------------------------------------------------------------
  Subroutine boost_nullify(params)
    Type(BoostParams), Intent(InOut)            :: params

    Integer :: error

    !** Nullify each
    Nullify(params%type%bVmax)
    Nullify(params%type%bVmin)
    Nullify(params%type%bVminMol)
    Nullify(params%type%bV)

  End Subroutine boost_nullify

  !----------------------------------------------------------------------------
  ! Returns the name of the molecule to be boosted, as specified in the 
  ! control file.
  !----------------------------------------------------------------------------
  Function boost_getBoostMolname(params)
    Type(BoostParams), Intent(In) :: params
    Character(len=strLen) :: boost_getBoostMolname

    !** Get the molecule name
    boost_getBoostMolname = params%molecName

  End Function boost_getBoostMolname

  !----------------------------------------------------------------------------
  ! Calculates the boosted potential and forces if required. Currently
  ! this assumes that you have passed and array of potentials corresponding
  ! to the entries in sorbates. Furthermore, we assume only ONE!!!!!! molecule
  ! in the sorbates array. This is a hack, but I need to get this working fast.
  ! Requires:  species -- coordinate storage structure for whole system
  !            noncoul -- spc-spc pair noncoul energies
  !            coul -- spc-spc pair coul energies
  !            intra -- spc intramolecular energies (NOT USED!)
  !            beta -- beta in unit of kcals
  !            deltat -- 
  !            pscale -- scaling factor for spc-spc energies
  !----------------------------------------------------------------------------
  Subroutine boost_calcBoostSingle(params,species,noncoul,coul,beta, &
      deltat,pscale)
    Type(BoostParams), Intent(InOut)               :: params
    Type(AtMolCoords), Intent(InOut), Dimension(:) :: species
    Real(Kind=RDbl), Dimension(:,:), Intent(In)    :: noncoul,coul
    Real(Kind=RDbl), Intent(InOut)                 :: deltat
    Real(Kind=RDbl), Intent(In)                    :: beta
    Real(Kind=RDbl), Dimension(:,:), Intent(Out)   :: pscale

    Integer                     :: nmoles, nmtypes, natoms, i, natypes, m, j
    Integer, Dimension(100)     :: atypes
    Real(Kind=RDbl)             :: totPot, pot 
    Type(VecType), Dimension(1) :: accel, pos

    !** Get the number of molecule types
    nmtypes = molecules_getnsorbs()

    !** Loop through the different species
    Do i = 1, nmtypes
      If (params%isBoost(i)) Then
        !** Loop through molecules
        Do m = 1, config_getnmoles(species,i)

          !** Get the atom types (needed for calculating forces)
          natypes = molecules_getnatomtypes(i,atypes,.True.)
          
          !** Get the number of atoms
          natoms = molecules_getnatoms(i)

          !** Calculate the total potential of the molecule
          !** For zeolite-adsorbate systems with rigid (or mostly
          !** rigid) molecules, it is the noncoulombic and coulombic
          !** interactions that we want to boost since they are
          !** large compared to intramolecular interactions.
          Do j = 1, nmtypes
            totPot = noncoul(i,j) + coul(i,j) 
          End Do
          pot = totPot

          !** Just figure out which boost to call
          If (Associated(params%type%bVmax)) Then
            Call boost_boostVmax(params,species(i)%coords(:,m)%r, &
                species(i)%afast(:,m), atypes(1:natoms), pot, beta, deltat)
          Else If (Associated(params%type%bVmin)) Then
            Call boost_boostVmin(params,species(i)%coords(:,m)%r, &
                species(i)%afast(:,m), atypes(1:natoms), pot, beta, deltat)
          Else If (Associated(params%type%bVminMol)) Then
            Call boost_boostVminMol(params,species(i)%coords(:,m)%r, &
                species(i)%afast(:,m), atypes(1:natoms), pot, beta, deltat)
          Else If (Associated(params%type%bV)) Then
            Call boost_boostV(params,species(i)%coords(:,m)%r, &
                species(i)%afast(:,m), atypes(1:natoms), pot, beta, deltat)
          End If

          !** store the potential energy scaling factor for return
          Do j = 1,nmtypes
            pscale(i,j) = 1.0_RDbl/totPot*pot
          End Do

        End Do
      End If
    End Do

  End Subroutine boost_calcBoostSingle

  !----------------------------------------------------------------------------
  ! Calculates the boosted potential and forces if required. Currently
  ! this assumes that you have passed and array of potentials corresponding
  ! to the entries in sorbates. Furthermore, we assume only ONE!!!!!! molecule
  ! in the sorbates array. This is a hack, but I need to get this working fast.
  ! Requires:  species -- coordinate storage structure for whole system
  !            noncoul -- spc-spc pair noncoul energies
  !            coul -- spc-spc pair coul energies
  !            intra -- spc intramolecular energies (NOT USED!)
  !            beta -- beta in unit of kcals
  !            deltat -- 
  !            pscale -- scaling factor for spc-spc energies
  !----------------------------------------------------------------------------
  Subroutine boost_calcBoostHack(params,species,energies,beta, &
      deltat,pscale)
    Type(BoostParams), Intent(InOut)               :: params
    Type(AtMolCoords), Intent(InOut), Dimension(:) :: species
    Real(Kind=RDbl), Dimension(:), Intent(InOut)   :: energies
    Real(Kind=RDbl), Intent(InOut)                 :: deltat
    Real(Kind=RDbl), Intent(In)                    :: beta
    Real(Kind=RDbl), Dimension(:,:), Intent(Out)   :: pscale

    Integer                     :: nmoles, nmtypes, natoms, i, natypes, m, j
    Integer, Dimension(100)     :: atypes
    Real(Kind=RDbl)             :: totPot, pot
    Type(VecType), Dimension(1) :: accel, pos

    !** Get the number of molecule types
    nmtypes = molecules_getnsorbs()

    !** Loop through the different species
    Do i = 1, nmtypes
      If (params%isBoost(i)) Then
        !** Loop through molecules
        Do m = 1, config_getnmoles(species,i)

          !** Get the atom types (needed for calculating forces)
          natypes = molecules_getnatomtypes(i,atypes,.True.)
          
          !** Get the number of atoms
          natoms = molecules_getnatoms(i)

          !** Total potential
          totPot = 0.0_RDbl
          Do j = 1, natoms
            totPot = totPot + energies(j)
          End Do

          !** Just figure out which boost to call
          If (Associated(params%type%bVminMol)) Then
            Call boost_boostVminMolHack(params,species(i)%coords(:,m)%r, &
                species(i)%afast(:,m), atypes(1:natoms), energies,beta,deltat)
          End If

          !** Total potential after boosting
          pot = 0.0_RDbl
          Do j = 1,natoms
            pot = pot + energies(j)
          End Do

          !** store the potential energy scaling factor for return
          Do j = 1,nmtypes
            pscale(i,j) = 1.0_RDbl/totPot*pot
          End Do

        End Do
      End If
    End Do

  End Subroutine boost_calcBoostHack

  !----------------------------------------------------------------------------
  ! Calculates the boosted potential and forces if required. This version
  ! works with multiple molecules to boost. Only ONE of the molecules is 
  ! boosted at any given time, and the molecule that is boosted will change
  ! depending upon the potential energy of the individual molecule.
  ! Requires:  species -- coordinate storage structure for whole system
  !            pot -- total potential for molecule
  !            beta -- beta in unit of kcals
  !            deltat -- time step after boosting
  !            pscale -- scaling factor for spc-spc energies
  !----------------------------------------------------------------------------
  Subroutine boost_calcBoostMult(params,species,pot,beta,deltat,pscale)
    Type(BoostParams), Intent(InOut)               :: params
    Type(AtMolCoords), Intent(InOut), Dimension(:) :: species
    Real(Kind=RDbl), Dimension(:), Intent(InOut)   :: pot
    Real(Kind=RDbl), Intent(InOut)                 :: deltat
    Real(Kind=RDbl), Intent(In)                    :: beta
    Real(Kind=RDbl), Dimension(:,:), Intent(Out)   :: pscale

    Integer                     :: nmoles, nmtypes, natoms, i, natypes, m, j
    Integer, Dimension(100)     :: atypes
    Real(Kind=RDbl)             :: totPot
    Type(VecType), Dimension(1) :: accel, pos

    !** Get the number of molecule types
    nmtypes = molecules_getnsorbs()

    !** Loop through the different species
    Do i = 1, nmtypes
      If (params%isBoost(i)) Then

        !** Get the total potential
        totPot = Sum(pot)

        !** Get the atom types (needed for calculating forces)
        natypes = molecules_getnatomtypes(i,atypes,.True.)
        
        !** Get the number of atoms
        natoms = molecules_getnatoms(i)
        
        !** Just figure out which boost to call
        If (Associated(params%type%bVmin)) Then
          Call boost_boostVminMult(params,species(i)%coords(:,:)%r, &
              species(i)%afast(:,:), atypes(1:natoms), pot, beta, deltat, m)
        End If

        !** store the potential energy scaling factor for return
        !** Hmm... This won't work too well here, since we are 
        !** scaling only one molecule. 
        Do j = 1,nmtypes
          pscale(i,j) = 1.0_RDbl/totPot*pot(m)
        End Do
        
      End If
    End Do

  End Subroutine boost_calcBoostMult

  !----------------------------------------------------------------------------
  ! This routine used the boost form given by reference 3
  !----------------------------------------------------------------------------
  Subroutine boost_boostVmax(params,coords,accels,atypes,pot,beta,deltat)
    Type(BoostParams), Intent(InOut)            :: params
    Type(VecType), Dimension(:), Intent(In) :: coords
    Type(VecType), Dimension(:), Intent(InOut) :: accels
    Integer, Dimension(:), Intent(In) :: atypes
    Real(Kind=RDbl), Intent(InOut) :: pot
    Real(Kind=RDbl), Intent(In)    :: beta
    Real(Kind=RDbl), Intent(InOut) :: deltat

    Real(Kind=RDbl) :: deltav
    Real(kind=RDbl) :: potBoost
    Real(Kind=RDbl) :: vmax
    Real(Kind=RDbl) :: s, f, dsdr
    Type(VecType)   :: force
    Real(kind=RDbl) :: gamma, c, scale
    Integer :: natoms, a

    !** Get the number of atoms (it's the size of the coords array)
    natoms = Size(coords,1)

    !** In the case of one molecule, vmax = pot
    vmax = pot

    !** Next we need deltaV, the difference between the threshhold
    !** potential
    deltav = params%vb - vmax

    !** Theta is a delta function that turns the boost on or off. It
    !** better be on at this point
    If (deltav < 0) Return

    !** Set the smoothing parameters to make this look neater
    c = params%type%bVmax%c
    gamma = params%type%bVmax%gamma

    !** Calculate the smoothing function f
    f = c*Exp(-gamma/deltav) / (1.0_RDbl+Exp(-gamma/deltav))

    !** Need to calculate s
    s = 1.0_RDbl + f

    !** Now we have the total boost
    potBoost = vmax*(1.0_RDbl-s)/s

    !** Calculate the new potential
    pot = potBoost + pot
    
    !** Need to back out the force in order to calculate the acceleration
    Do a = 1, natoms
      scale = atom_getmass(atypes(a))/scalef
      force = accels(a)*scale
    
      !** Calculate the boosted force. First we need the
      !** derivative of s
      dsdr = gamma*f*(f-1.0_RDbl)/(c*deltav*deltav)

      !** Calculate the boosted force
      force = force*(1.0_RDbl + (1.0_RDbl-s)/s - pot/(s*s)*dsdr)

      !** Get the new accelerations
      accels(a) = force*scalef*atom_invmass(atypes(a))

    End Do

    !** Finally, the time
    deltat = deltat*exp(beta*potBoost)

    !** Update the boost statistic with the boost value
    Call stats_update(params%bstats,potBoost*scalepe)
    Call stats_update(params%tstats,deltat)

  End Subroutine boost_boostVmax

  !----------------------------------------------------------------------------
  ! This routine uses the boost suggested by reference (3), specificaly
  ! the one based on minimum energy differences
  !----------------------------------------------------------------------------
  Subroutine boost_boostVmin(params,coords,accels,atypes,pot,beta,deltat)
    Type(BoostParams), Intent(InOut)           :: params
    Type(VecType), Dimension(:), Intent(In) :: coords
    Type(VecType), Dimension(:), Intent(InOut) :: accels
    Integer, Dimension(:), Intent(In) :: atypes
    Real(Kind=RDbl), Intent(InOut) :: pot
    Real(Kind=RDbl), Intent(In)    :: beta
    Real(Kind=RDbl), Intent(InOut) :: deltat
    Real(Kind=RDbl) :: g, deltaVmin, dgdr, n, c, dVminn, potBoost
    Real(kind=RDbl) :: scale, Vthresh, origPot, dVminnm1
    Type(VecType)   :: force
    Integer :: natoms, i, a, whichvb

    !** Get the number of atoms
    natoms = Size(coords,1)

    !** Zero the potential boost (for statistics)
    potBoost = 0.0_RDbl
    deltaVmin = 0.0_RDbl
    g = 0.0_RDbl

    !** Loop through the atoms
    Do a = 1, natoms

      !** Loop through the boosts and figure out which to use
      deltaVmin = params%type%bVmin%vb(1) - pot
      whichvb = 1
      Do i = 2, params%type%bVmin%nvb
        If ((params%type%bVmin%vb(i) - pot) < deltaVmin) Then
          deltaVmin = params%type%bVmin%vb(i) - pot
          whichvb = i
          Vthresh = params%type%bVmin%vb(i)
        End If
      End Do

      !** Check the delta function to make sure we should be doing this
      If (deltaVmin < 0) Then
        Cycle
      End If

      !** Define c and n to make things neater
      c = params%type%bVmin%c
      n = params%type%bVmin%n
      dVminn = deltaVmin**n
      dVminnm1 = dVminn/deltaVmin

      !** Calculate the value of g
      g = c*dVminn/(1.0_RDbl + c*(dVminn))

      !** Calculate the constant part of dgdr. 
      dgdr = c*n*dVminnm1/(1.0_RDbl+c*dVminn) &
          *(c*dVminn/(1.0_RDbl+c*dVminn) - 1.0_RDbl)
      
      !** Calculate the boost
      potBoost = deltaVmin*g

      !** Calculate the boosted potential
      origPot = pot
      pot = pot + potBoost

      !** Back out the force from the acceleration
      scale = atom_getmass(atypes(a))/scalef
      force = accels(a)*scale

      !** Calculate the boosted force
      force = force*(1.0_RDbl + deltaVmin*dgdr - g)

      !** Get the new acceleration
      accels(a) = force*scalef*atom_invmass(atypes(a))

      !** The new time step
      deltat = deltat*Exp(beta*potBoost)

    End Do

    !** Update the boost statistic with the boost value
    Call stats_update(params%bstats,potBoost*scalepe)
    Call stats_update(params%tstats,deltat)
    Call stats_update(params%dvstats,deltaVmin*scalepe)
    Call stats_update(params%sstats,g)
    Call stats_update(params%drstats,dgdr)

  End Subroutine boost_boostVmin

  !----------------------------------------------------------------------------
  ! This routine uses the boost suggested by reference (3), specificaly
  ! the one based on minimum energy differences. This assumes a single atom
  ! molecule, sorry.
  !----------------------------------------------------------------------------
  Subroutine boost_boostVminMult(params,coords,accels,atypes,pot,beta, &
      deltat,whichmol)
    Type(BoostParams), Intent(InOut)           :: params
    Type(VecType), Dimension(:,:), Intent(In) :: coords
    Type(VecType), Dimension(:,:), Intent(InOut) :: accels
    Integer, Dimension(:), Intent(In) :: atypes
    Real(Kind=RDbl), Dimension(:), Intent(InOut) :: pot
    Real(Kind=RDbl), Intent(In)    :: beta
    Integer, Intent(Out) :: whichmol  ! molecule that is boosted
    Real(Kind=RDbl), Intent(InOut) :: deltat
    Real(Kind=RDbl) :: g, deltaVmin, dgdr, n, c, dVminn, potBoost
    Real(kind=RDbl) :: scale, Vthresh, origPot, dVminnm1
    Type(VecType)   :: force
    Integer :: nmoles, i, m
    Integer :: whichvb   ! Which of the thresholds yields smallest delta Vb

    !** Get the number of atoms
    nmoles = Size(coords,2)

    !** Zero the potential boost (for statistics)
    potBoost = 0.0_RDbl
    deltaVmin = 0.0_RDbl
    g = 0.0_RDbl

    !** Calculate the first vmin
    !** Loop through the boosts and figure out which to use
    deltaVmin = params%type%bVmin%vb(1) - pot(1)
    whichvb = 1
    whichmol= 1
    Do i = 2, params%type%bVmin%nvb
      If ((params%type%bVmin%vb(i) - pot(1)) < deltaVmin) Then
        deltaVmin = params%type%bVmin%vb(i) - pot(1)
        whichvb = i
        whichmol= m
        Vthresh = params%type%bVmin%vb(i)
      End If
    End Do

    !** Loop through the remaining molecules
    Do m = 2, nmoles

      !** Loop through the boosts and figure out which to use
      Do i = 1, params%type%bVmin%nvb
        If ((params%type%bVmin%vb(i) - pot(m)) < deltaVmin) Then
          deltaVmin = params%type%bVmin%vb(i) - pot(m)
          whichvb = i
          whichmol= m
          Vthresh = params%type%bVmin%vb(i)
        End If
      End Do
    End Do

    !** Check the delta function to make sure we should be doing this
    If (deltaVmin < 0) Then
      Return
    End If

    !** Define c and n to make things neater
    c = params%type%bVmin%c
    n = params%type%bVmin%n
    dVminn = deltaVmin**n
    dVminnm1 = dVminn/deltaVmin
    
    !** Calculate the value of g
    g = c*dVminn/(1.0_RDbl + c*(dVminn))
    
    !** Calculate the constant part of dgdr. 
    dgdr = c*n*dVminnm1/(1.0_RDbl+c*dVminn) &
        *(c*dVminn/(1.0_RDbl+c*dVminn) - 1.0_RDbl)
    
    !** Calculate the boost
    potBoost = deltaVmin*g
    
    !** Calculate the boosted potential
    origPot = pot(m)
    pot(m) = pot(m) + potBoost
    
    !** Back out the force from the acceleration
    scale = atom_getmass(atypes(1))/scalef
    force = accels(1,m)*scale
    
    !** Calculate the boosted force
    force = force*(1.0_RDbl + deltaVmin*dgdr - g)
    
    !** Get the new acceleration
    accels(1,m) = force*scalef*atom_invmass(atypes(1))
    
    !** The new time step
    deltat = deltat*Exp(beta*potBoost)
    
    !** Update the boost statistic with the boost value
    Call stats_update(params%bstats,potBoost*scalepe)
    Call stats_update(params%tstats,deltat)
    Call stats_update(params%dvstats,deltaVmin*scalepe)
    Call stats_update(params%sstats,g)
    Call stats_update(params%drstats,dgdr)

  End Subroutine boost_boostVminMult

  !----------------------------------------------------------------------------
  ! This routine uses the boost suggested by reference (3), specificaly
  ! the one based on minimum energy differences. This version uses the
  ! TOTAL molecule energy to calculate the boost. If the total molecule
  ! energy is below the threshold, the entire molecule is boosted. The
  ! boost is divided equally among the atoms of the molecule in computing
  ! the modified forces. And imporvement would be dividing the boost
  ! potential according to the distribution of the (magnitude of the)
  ! forces on the atoms.
  !----------------------------------------------------------------------------
  Subroutine boost_boostVminMol(params,coords,accels,atypes,pot,beta,deltat)
    Type(BoostParams), Intent(InOut)           :: params
    Type(VecType), Dimension(:), Intent(In) :: coords
    Type(VecType), Dimension(:), Intent(InOut) :: accels
    Integer, Dimension(:), Intent(In) :: atypes
    Real(Kind=RDbl), Intent(InOut) :: pot  ! TOTAL potential on the molecule
    Real(Kind=RDbl), Intent(In)    :: beta
    Real(Kind=RDbl), Intent(InOut) :: deltat
    Real(Kind=RDbl) :: g, deltaVmin, dgdr, n, c, dVminn, potBoost, totalmass
    Real(kind=RDbl) :: scale, Vthresh, origPot, dVminnm1, dVmin, dtorig
    Type(VecType)   :: force, totalForce, unitvec
    Integer :: natoms, i, a, whichvb

    !** Keep the original dt
    dtorig = deltat

    !** Get the number of atoms
    natoms = Size(coords,1)

    !** Get the total mass
    totalmass = 0.0_RDbl
    Do a = 1, natoms
      totalmass = atom_getmass(atypes(a)) + totalmass
    End Do

    !** Zero the potential boost (for statistics)
    potBoost = 0.0_RDbl
    deltaVmin = 0.0_RDbl
    g = 0.0_RDbl

    !** Loop through the boosts and figure out which to use
    deltaVmin = params%type%bVminMol%vb(1) - pot
    whichvb = 1
    Do i = 2, params%type%bVminMol%nvb
      If ((params%type%bVminMol%vb(i) - pot) < deltaVmin) Then
        deltaVmin = params%type%bVminMol%vb(i) - pot
        whichvb = i
        Vthresh = params%type%bVminMol%vb(i)
      End If
    End Do

    !** Check the delta function to make sure we should be doing this
    If (deltaVmin < 0) Then
      Return
    End If

    !** Define c and n to make things neater
    c = params%type%bVminMol%c
    n = params%type%bVminMol%n
    dVminn = deltaVmin**n
    dVminnm1 = dVminn/deltaVmin
    
    !** Calculate the value of g
    g = c*dVminn/(1.0_RDbl + c*(dVminn))

    !MDEBUG
    !** Make it 1.
!    g = 1.0_RDbl

    !** Calculate the boost
    potBoost = deltaVmin*g
    
    !** Calculate the boosted potential
    origPot = pot
    pot = pot + potBoost
    
    !** The new time step
    deltat = deltat*Exp(beta*potBoost)

    
!!$! This calculates using the COM force to get the forces.
!!$
!!$    !** Get the total COM force
!!$    totalForce = 0.0_RDbl
!!$    Do a = 1, natoms
!!$
!!$      !** Back out the force from the acceleration
!!$      scale = atom_getmass(atypes(a))/scalef
!!$      force = accels(a)*scale
!!$      totalForce = totalForce + force*atom_getmass(atypes(a))/totalmass
!!$
!!$    End Do
!!$
!!$!    !MDEBUG
!!$!    !** Report the total force, mag and direction
!!$!    Write(0,'(a,i6,a,4f10.3)') __FILE__,__LINE__," force : ",mag(totalForce),&
!!$!        totalForce/mag(totalForce)
!!$
!!$    !** Get the constants for the derivatives based on the boost potential
!!$    dVmin = deltaVmin
!!$    dVminn = dVmin**n
!!$    dVminnm1 = dVminn/dVmin
!!$    
!!$    !** Calculate the constant part of dgdr. 
!!$    dgdr = c*n*dVminnm1/(1.0_RDbl+c*dVminn) &
!!$        *(c*dVminn/(1.0_RDbl+c*dVminn) - 1.0_RDbl)
!!$
!!$    !** Calculate the boosted force on the center of the molecule
!!$    totalforce = totalforce*(1.0_RDbl + dVmin*dgdr - g)
!!$    
!!$    !** Now modify the force using the total potential
!!$    Do a = 1, natoms
!!$      
!!$      !** Assign forces to each of the atoms
!!$      force = totalForce*atom_getmass(atypes(a))/totalmass
!!$
!!$      !** Get the new acceleration
!!$      accels(a) = force*scalef*atom_invmass(atypes(a))
!!$
!!$    End Do
!!$
!!$!    !MDEBUG
!!$!    !** Report the new boosted force on the COM
!!$!    Do a = 1, natoms
!!$!      !** Back out the force from the acceleration
!!$!      scale = atom_getmass(atypes(a))/scalef
!!$!      force = accels(a)*scale
!!$!      totalForce = totalForce + force*atom_getmass(atypes(a))/totalmass
!!$!    End Do
!!$!    Write(0,'(a,i6,a,4f10.3)') __FILE__,__LINE__," force : ",mag(totalForce),&
!!$!        totalForce/mag(totalForce)
!!$
!!$! This calculates using the boost potential divided by natoms

    !** Calculate the new deltaVmin, based on the number of atoms in 
    !** the molecule
    dVmin = deltaVmin/Real(natoms,Kind=RDbl)
    
    !** Get the constants for the derivatives based on the new 
    !** (distributed) boost potential
    dVminn = dVmin**n
    dVminnm1 = dVminn/dVmin
    
    !** Calculate the constant part of dgdr. 
    dgdr = c*n*dVminnm1/(1.0_RDbl+c*dVminn) &
        *(c*dVminn/(1.0_RDbl+c*dVminn) - 1.0_RDbl)

    !** Loop through the atoms
    Do a = 1, natoms

      !** Back out the force from the acceleration
      scale = atom_getmass(atypes(a))/scalef
      force = accels(a)*scale

!      Write(0,'(a,i6,i4,3e15.3)') __FILE__,__LINE__,a,force

      !** Calculate the boosted force on the center of the molecule
      force = force*(1.0_RDbl + dVmin*dgdr - g)

!      Write(0,'(a,i6,i4,3e15.3)') __FILE__,__LINE__,a,force

      !** Get the new acceleration
      accels(a) = force*scalef*atom_invmass(atypes(a))

!      accels(a) = 0.0_RDBL

    End Do

! This calculates using the magnitude of the new force and the unit vectors
!!$    !** Calculate the new deltaVmin
!!$    dVmin = deltaVmin
!!$
!!$    !** Get the constants for the derivatives based on the new 
!!$    !** (distributed) boost potential
!!$    dVminn = dVmin**n
!!$    dVminnm1 = dVminn/dVmin
!!$    
!!$    !** Calculate the constant part of dgdr. 
!!$    dgdr = c*n*dVminnm1/(1.0_RDbl+c*dVminn) &
!!$        *(c*dVminn/(1.0_RDbl+c*dVminn) - 1.0_RDbl)
!!$
!!$    !** Zero our force summation
!!$    totalForce = 0.0_RDbl
!!$
!!$    !** Loop through the atoms
!!$    Do a = 1, natoms
!!$
!!$      !** Back out the force from the acceleration
!!$      scale = atom_getmass(atypes(a))/scalef
!!$      force = accels(a)*scale
!!$
!!$      !** Update the sum
!!$      totalForce = force + totalForce
!!$
!!$    End Do
!!$
!!$    !** Calculate the boosted force on the center of the molecule
!!$    totalForce = totalForce*(1.0_RDbl + dVmin*dgdr - g)
!!$
!!$    !** Calculate the force magnitude
!!$    magForce = mag(totalForce)
!!$
!!$    !** Loop through and calculate the new forces based on the directions
!!$    !** of the old
!!$    Do a = 1, natoms
!!$
!!$      !** Back out the force from the acceleration
!!$      scale = atom_getmass(atypes(a))/scalef
!!$      force = accels(a)*scale
!!$
!!$      !** Unit vector of the force
!!$      unitvec = force/mag(force)
!!$
!!$      !** Get the new acceleration
!!$      accels(a) = magForce*unitvec*scalef*atom_invmass(atypes(a))
!!$
!!$      accels(a) = 0.0_RDbl
!!$
!!$    End Do

    !** Update the boost statistic with the boost value
    Call stats_update(params%bstats,potBoost*scalepe)
    Call stats_update(params%tstats,deltat)
    Call stats_update(params%dvstats,deltaVmin*scalepe)
    Call stats_update(params%sstats,g)
    Call stats_update(params%drstats,dgdr)

  End Subroutine boost_boostVminMol

  !----------------------------------------------------------------------------
  ! This routine uses the boost suggested by reference (3), specificaly
  ! the one based on minimum energy differences. This version uses the
  ! TOTAL molecule energy to calculate the boost. If the total molecule
  ! energy is below the threshold, the entire molecule is boosted. The
  ! boost is divided equally among the atoms of the molecule in computing
  ! the modified forces. And imporvement would be dividing the boost
  ! potential according to the distribution of the (magnitude of the)
  ! forces on the atoms.
  !----------------------------------------------------------------------------
  Subroutine boost_boostVminMolHack(params,coords,accels,atypes,pots,beta, &
      deltat)
    Type(BoostParams), Intent(InOut)           :: params
    Type(VecType), Dimension(:), Intent(In) :: coords
    Type(VecType), Dimension(:), Intent(InOut) :: accels
    Integer, Dimension(:), Intent(In) :: atypes
    Real(Kind=RDbl), Dimension(:), Intent(InOut) :: pots
    Real(Kind=RDbl), Intent(In)    :: beta
    Real(Kind=RDbl), Intent(InOut) :: deltat
    Real(Kind=RDbl) :: g, deltaVmin, dgdr, n, c, dVminn, potBoost, pot
    Real(kind=RDbl) :: scale, Vthresh, origPot, dVminnm1, dVmin, dtorig
    Type(VecType)   :: force, totalForce, unitvec
    Integer :: natoms, i, a, whichvb

    !** Keep the original dt
    dtorig = deltat

    !** Get the number of atoms
    natoms = Size(coords,1)

    !** Zero the potential boost (for statistics)
    potBoost = 0.0_RDbl
    deltaVmin = 0.0_RDbl
    g = 0.0_RDbl

    !** Get the total potential
    pot = 0.0_RDbl
    Do i = 1, natoms
      pot = pot + pots(i)
    End Do

    !** Loop through the boosts and figure out which to use
    deltaVmin = params%type%bVminMol%vb(1) - pot
    whichvb = 1
    Do i = 2, params%type%bVminMol%nvb
      If ((params%type%bVminMol%vb(i) - pot) < deltaVmin) Then
        deltaVmin = params%type%bVminMol%vb(i) - pot
        whichvb = i
        Vthresh = params%type%bVminMol%vb(i)
      End If
    End Do

    !** Check the delta function to make sure we should be doing this
    If (deltaVmin < 0) Then
      Return
    End If

    !** Define c and n to make things neater
    c = params%type%bVminMol%c
    n = params%type%bVminMol%n
    dVminn = deltaVmin**n
    dVminnm1 = dVminn/deltaVmin
    
    !** Calculate the value of g
    g = c*dVminn/(1.0_RDbl + c*(dVminn))

!!$    !MDEBUG
!!$    !** Make it 1.
!!$    g = 1.0_RDbl

    !** Calculate the boost
    potBoost = deltaVmin*g
    
    !** Calculate the boosted potential
    origPot = pot
    pot = pot + potBoost
    
    !** The new time step
    deltat = deltat*Exp(beta*potBoost)

    
! This calculates using the boost potential divided by natoms

    !** Loop through the atoms
    Do a = 1, natoms

      !MDEBUG
      Write(0,'(a,i6,a,i3,5f8.3)') __FILE__,__LINE__,"oldpot ",a,pots(a),pots(a)/origpot,potBoost/origpot*pots(a),potBoost,pot

      !** Adjust the potential on each atom
      pots(a) = pots(a) + potBoost/origpot*pots(a)

      !MDEBUG
      Write(0,'(a,i6,a,i3,3f8.3)') __FILE__,__LINE__,"newpot ",a,pots(a)

      !** Calculate the new deltaVmin, based on the relative potential
      !** of each atom
      dVmin = deltaVmin/pot*pots(a)
      
      !** Get the constants for the derivatives based on the new 
      !** (distributed) boost potential
      dVminn = dVmin**n
      dVminnm1 = dVminn/dVmin
      
      !** Calculate the constant part of dgdr. 
      dgdr = c*n*dVminnm1/(1.0_RDbl+c*dVminn) &
          *(c*dVminn/(1.0_RDbl+c*dVminn) - 1.0_RDbl)

      !** Back out the force from the acceleration
      scale = atom_getmass(atypes(a))/scalef
      force = accels(a)*scale

!!$      If (dtorig*1.10_RDbl >= deltat) Then
!!$        Write(0,'(a,i6,i4,3e15.3)') __FILE__,__LINE__,a,force
!!$      End If

      !** Calculate the boosted force on the center of the molecule
      force = force*(1.0_RDbl + dVmin*dgdr - g)

!!$      If (dtorig*1.10_RDbl >= deltat) Then
!!$        Write(0,'(a,i6,i4,3e15.3)') __FILE__,__LINE__,a,force
!!$      End If

      !** Get the new acceleration
      accels(a) = force*scalef*atom_invmass(atypes(a))
!!$
!!$      accels(a) = 0.0_RDBL

    End Do

    !** Update the boost statistic with the boost value
    Call stats_update(params%bstats,potBoost*scalepe)
    Call stats_update(params%tstats,deltat)
    Call stats_update(params%dvstats,deltaVmin*scalepe)
    Call stats_update(params%sstats,g)
    Call stats_update(params%drstats,dgdr)

  End Subroutine boost_boostVminMolHack

  !----------------------------------------------------------------------------
  ! This routine uses a simple static boost with minimal smoothing
  ! as described in reference 4. Gives a very flat potential surface.
  !----------------------------------------------------------------------------
  Subroutine boost_boostV(params,coords,accels,atypes,pot,beta,deltat)
    Type(BoostParams), Intent(InOut)           :: params
    Type(VecType), Dimension(:), Intent(In) :: coords
    Type(VecType), Dimension(:), Intent(InOut) :: accels
    Integer, Dimension(:), Intent(In) :: atypes
    Real(Kind=RDbl), Intent(InOut) :: pot  ! TOTAL potential on the molecule
    Real(Kind=RDbl), Intent(In)    :: beta
    Real(Kind=RDbl), Intent(InOut) :: deltat
    Real(Kind=RDbl) :: deltaV, potBoost, vmeb
    Real(kind=RDbl) :: scale, dtorig
    Real(Kind=RDbl) :: alpha, gamma, g, h
    Type(VecType)   :: force
    Integer :: natoms, i, a

    !** Keep the original dt
    dtorig = deltat

    !** Get the number of atoms
    natoms = Size(atypes,1)

    !** Zero the potential boost (for statistics)
    potBoost = 0.0_RDbl
    deltaV = 0.0_RDbl
    g = 0.0_RDbl
    h = 0.0_RDbl

    !** Loop through the boosts and figure out which to use
    deltaV = params%type%bV%vb - pot

    !** Check the delta function to make sure we should be doing this
    If (deltaV < 0) Then
      Return
    End If

    !** Calculate the smoothing. It is usually 1 unless we are close
    !** to the actual potential surface
    alpha = params%type%bV%alpha
    gamma = params%type%bV%gamma
    vmeb = -deltaV
    g = 1.0_RDbl/(1.0_RDbl + Exp(-alpha*vmeb))
    h = Exp(-gamma*gamma*alpha*alpha*vmeb*vmeb)/alpha

    !** Calculate the boost
    potBoost = vmeb*g + h + deltaV
    
    !** Calculate the boosted potential
    pot = pot + potBoost
    
    !** The new time step
    deltat = deltat*Exp(beta*potBoost)
    
    !** Need to zero the accelerations
    Do a = 1, natoms

      !** Back out the force from the acceleration
      scale = atom_getmass(atypes(a))/scalef
      force = accels(a)*scale
      
      !** Smoothed forces
      force = force*(deltaV*alpha*Exp(-alpha*vmeb)*g*g - g &
          + 2.0_RDbl*gamma*gamma*alpha*alpha*vmeb*h)

      !** Get the new acceleration
      accels(a) = force*scalef*atom_invmass(atypes(a))

    End Do

    !** Update the boost statistic with the boost value
    Call stats_update(params%bstats,potBoost*scalepe)
    Call stats_update(params%tstats,deltat)
    Call stats_update(params%dvstats,deltaV*scalepe)
    Call stats_update(params%sstats,g)
    Call stats_update(params%drstats,h)

  End Subroutine boost_boostV

  !----------------------------------------------------------------------------
  ! Displays information about the boost potential and its parameters
  !----------------------------------------------------------------------------
  Subroutine boost_display(params,unitno,nspaces)
    Type(BoostParams), Intent(In)           :: params
    Integer, Intent(in) :: unitno, nspaces
    Integer :: i
    Character(len=nspaces):: space

    !** Make the spacing
    space = ""
    Do i = 1, nspaces
      space = " "//space
    End Do

    !** Write the important information to the given unitno
    Write(unitno,'(2a)') space, dashedline
    Write(unitno,'(2a)') space,'Boosted MD Information'
    Do i = 1, Size(params%isBoost)
      If (.Not.params%isBoost(i)) Cycle
      Write(unitno,'(3a)') space,'Boosted Molecule : ',Trim(params%molecName)
    End Do
    If (Associated(params%type%bVmax)) Then
      Write(unitno,'(2a)') space,'Boost Type       : STATICVMAX'
      Write(unitno,'(2a,f8.3,1x,a)') space, 'Boost threshhold : ', &
          params%type%bVmax%vb*scalepe,Trim(vbunits)
      Write(unitno,'(2a,f8.3)') space,'Smoothing gamma  : ', &
          params%type%bVmax%gamma
      Write(unitno,'(2a,f8.3)') space,'Smoothing c      : ', &
          params%type%bVmax%c
    Else If (Associated(params%type%bVmin)) Then
      Write(unitno,'(2a)') space,'Boost Type       : STATICVMIN'
      Write(unitno,'(2a,i3,a)') space,'Boost threshholds (',&
          params%type%bVmin%nvb,' total)'
      Do i = 1, params%type%bVmin%nvb
        Write(unitno,'(2x,a,f8.3,1x,a)') space, &
            params%type%bVmin%vb(i)*scalepe, Trim(vbunits)
      End Do
      Write(unitno,'(2a,f8.3)') space,'Smoothing n      : ', &
          params%type%bVmin%n
      Write(unitno,'(2a,f8.3)') space,'Smoothing c      : ', &
          params%type%bVmin%c
    Else If (Associated(params%type%bVminMol)) Then
      Write(unitno,'(2a)') space,'Boost Type       : STATICVMINMOL'
      Write(unitno,'(2a,i3,a)') space,'Boost threshholds (',&
          params%type%bVminMol%nvb,' total)'
      Do i = 1, params%type%bVminMol%nvb
        Write(unitno,'(2x,a,f8.3,1x,a)') space, &
            params%type%bVminMol%vb(i)*scalepe, Trim(vbunits)
      End Do
      Write(unitno,'(2a,f8.3)') space,'Smoothing n      : ', &
          params%type%bVminMol%n
      Write(unitno,'(2a,f8.3)') space,'Smoothing c      : ', &
          params%type%bVminMol%c
    Else If (Associated(params%type%bV)) Then
      Write(unitno,'(2a)') space,'Boost Type       : STATIC'
      Write(unitno,'(2a)') space,'Boost threshholds ( 1 total)'
      Write(unitno,'(2x,a,f8.3,1x,a)') space, &
            params%type%bV%vb*scalepe, Trim(vbunits)
    End If

  End Subroutine boost_display

  !----------------------------------------------------------------------------
  ! Displays boost potential statistics
  !----------------------------------------------------------------------------
  Subroutine boost_displayStats(params,unitno,nspaces)
    Type(BoostParams), Intent(In)           :: params
    Integer, Intent(in) :: unitno, nspaces
    Integer :: i
    Character(len=nspaces):: space
    Logical :: isFirst

    !** Make the spacing
    space = ""
    Do i = 1, nspaces
      space = " "//space
    End Do

    !** This flag determines whether or not we display the boost information
    !** for the first time (i.e., the first boosted molecule).
    !** Unless we find a boosted molecule in the list, we will not write
    !** the line 'Boosting Energy' to the screen
    isFirst = .True.

    Do i = 1, Size(params%isBoost)
      If (.Not.params%isBoost(i)) Cycle
      If (isFirst) Then
        !** This is the first boosted molecule to be displayed. Write
        !** the header line first, and change isFirst to FALSE so we
        !** don't write it again.
        isFirst = .False.
        Write(unitno,'(2a)') space,'Boosting Energy'
      End If
      Write(unitno,'(a,2x,a12,a3,4f15.6)') space,Trim(molecules_name(i)), &
          " : ", &
          stats_getvalue(params%bstats), stats_getcavg(params%bstats), &
          stats_getblock(params%bstats), stats_getstd(params%bstats)
    End Do

    !** Write the time step statistics
    If (.Not.isFirst) Then 
      Write(unitno,'(2a)') space,'Boosted Statistics'
      Write(unitno,'(a,2x,a15,4f15.6)') space,'Timestep Size : ', &
          stats_getvalue(params%tstats), stats_getcavg(params%tstats), &
          stats_getblock(params%tstats), stats_getstd(params%tstats)
      Write(unitno,'(a,2x,a15,4f15.6)') space,'Delta V       : ', &
          stats_getvalue(params%dvstats), stats_getcavg(params%dvstats), &
          stats_getblock(params%dvstats), stats_getstd(params%dvstats)
      Write(unitno,'(a,2x,a15,4f15.6)') space,'Smoothing Func: ', &
          stats_getvalue(params%sstats), stats_getcavg(params%sstats), &
          stats_getblock(params%sstats), stats_getstd(params%sstats)
      Write(unitno,'(a,2x,a15,4f15.6)') space,'Derivative Val: ', &
          stats_getvalue(params%drstats), stats_getcavg(params%drstats), &
          stats_getblock(params%drstats), stats_getstd(params%drstats)
    End If

  End Subroutine boost_displayStats

End Module boost
