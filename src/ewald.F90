!------------------------------------------------------------------------------
! This modules contains routines for those nasty ewald summations.
! Good reference for ewald: Allen and Tildesley, Computer Simulation of 
! Liquids, p 155-162; Amit's notes for doing part of the calculation in the
! complex plane, making a double summation into a single summation
!------------------------------------------------------------------------------

Module ewald

  Use defaults, Only: Pi, RDbl, zero, strLen, lstrLen, one, twopi, &
      zeroTolerance, scalepe, scalef, lengthUnit, LOW_DIST_CUTOFF
  Use utils, Only: erfc, stripcmnt, tolower, split, toreal, toint, findint
  Use vector, Only: VecType, Operator(*), Assignment(=), Operator(+), &
      Operator(-), Operator(/), mag, vector_sumVecs
  Use fundcell, Only: Fundamental_Cell
  Use simcell, Only: Simcell_Params, simcell_minimage

  Implicit None

  Private
  Public :: EwaldParams, ewald_idstring, ewald_init, ewald_initCalc, &
      ewald_getline, ewald_display, ewald_getinteractionSngl, &
      ewald_getinteractionHOT, ewald_getCoreShellCorr, epots

  Interface ewald_rwald
     Module Procedure ewald_rwaldHOT
     Module Procedure ewald_rwaldFvec
  End Interface

  Interface ewald_skwald
     Module Procedure ewald_skwaldHOT
     Module Procedure ewald_skwaldFvec
  End Interface

  Interface ewald_selfterm
    Module Procedure ewald_selftermSngl
    Module Procedure ewald_selftermMult
  End Interface

  Interface ewald_getInteraction
    Module Procedure ewald_getInteractionSngl
    Module Procedure ewald_getInteractionMult
    Module Procedure ewald_getInteractionHOT
  End Interface

  !** These were the cutoffs in the original Snurr group Ewald code.
  !** They are here for historic preservation reasons. 
  Real(Kind=RDbl), Parameter :: FACTOR_LOCUT = 1.0E-10_RDbl
  Real(Kind=RDbl), Parameter :: SFACTOR_LOCUT = 1.0E-13_RDbl
  Real(Kind=RDbl), Parameter :: SELF_LOCUT = 1.0E-6_RDbl

  !** Interaction type identifier string
  Character(len=strLen), Parameter :: ewald_idstring = "EWALD"

  !** Converts from e^2/Ang to kcal/mol. Remember that the charges'
  !** values are in elementry charge units e = 1.602177E-19 C. 
  !** Usual derivations (this one included) drop the 4*Pi*eps0 term
  !** from the equations, where eps0 = 8.85419E-12 C^2/(J*m).
  !** This yields units of e^2/Ang. The conversion factor, therefore, is 
  Real(Kind=RDbl), Parameter :: tokcal = 332.064_RDbl

  Type EwaldParams
    Character(len=lstrLen) :: line
    Real(kind=RDbl)        :: kappa ! kappaParam/Min(SimBox Length) (AKA alpha)
    Real(kind=RDbl)        :: kappaParam  ! Gaussian width
    Real(kind=RDbl)        :: locut     ! Low cutoff for reciprocal space sum
    Real(Kind=RDbl)        :: rcut      ! Real space spherical cutoff
    Logical                :: useSF     ! If TRUE, use the structure factor
    Integer                :: kmax  ! No. of recip. space images in each dir
    Integer                :: nk, nkmax ! No. of accepted kvecs, max no. kvecs
    Real(Kind=RDbl)        :: length    ! Min(box edge lengths in x,y,z)
    Type(VecType), Dimension(:), Pointer   :: k
    Real(Kind=RDbl), Dimension(:), Pointer :: s
    Real(Kind=RDbl), Dimension(:), Pointer :: si
  End Type EwaldParams

  !** This is used in the older double sum implementation of the Ewald sum
  Type KwaldVecs
    Real(Kind=RDbl) :: k, kx, ky, kz, kxy, kxz, kyz, kxyz
  End Type KwaldVecs

  !** Units for the energy cutoff
  Character(len=strLen), Parameter :: cutUnits = "Ang"

  !** Ewald parameters variable
  Type(EwaldParams), Save :: eparams

  !** Ewald kvectors array for double-sum implementation
  Type(KwaldVecs), Dimension(:), Allocatable, Save :: kvec
  Real(Kind=RDbl), Save :: igamma, ibeta, ialpha
  Logical, Save :: firstKvec

  !** Debugging
  Real(Kind=RDbl), Dimension(4) :: epots

Contains
  !----------------------------------------------------------------------------
  ! Initializes - you guessed it - the ewald data type from the given unit no.
  !----------------------------------------------------------------------------
  Subroutine ewald_init(line)
    Character(*), Intent(InOut) :: line
    Character(len=strLen), Dimension(strLen) :: fields, chunks
    Integer :: nfields, nchunks, i, error

    !** No nasty comments, please. Remove everything after the comment char
    line = stripcmnt(line)

    !** Store the line from the interaction file. We may want it later.
    eparams%line = line

    !** Set a default cutoff of -1.0. This indicates there is no cutoff for
    !** the real space summation. Some people use one, so it may be useful
    !** for comparing results.
    eparams%rcut = -1.0_RDbl

    !** Set the default for the structure factor. Default is FALSE
    eparams%useSF = .FALSE.

    !** Split the line into fields and then extract the data
    nfields = split(line,fields)
    Do i = 1, nfields
      nchunks = split(fields(i),chunks,'@')

      Select Case(tolower(chunks(1)))
      Case ('kmax')
        !** Number of k "boxes" in each direction (x,y,z).
        eparams%kmax = toint(chunks(2))
        !MDEBUG
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            " WARNING: nkmax is hard coded to 50000"
        !** Limit the number of reciprocal space vectors. There is nothing
        !** special about this limit, it's just that the higher the number
        !** of k vectors, the more computationally intensive kwald is.
        eparams%nkmax = 50000  
      Case ('kappa')
       !** Width of the gaussian distribution.
        eparams%kappaParam = toreal(chunks(2))
      Case ('locut')
        !** Low end cutoff for the reciprocal space k vectors.
        eparams%locut = toreal(chunks(2))
      Case ('rcut')
        !** Real space cutoff. If it is less than zero, it is disabled.
        eparams%rcut = toreal(chunks(2))
      Case ('sfactor')
        !** We want to use the structure factor.
        eparams%useSF = .TRUE.
      End Select
    End Do

    !** Set the firstKvec flag to TRUE. We need to calculate the 
    !** double sum kvecs only once.
    firstKvec = .True.

    !** Allocate space for the k vector array
    Allocate(eparams%k(eparams%nkmax),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Cannot allocate space for eparams%k'
      Stop              
    End If

    !** ..and the s array. This holds the real values for the k-space sum
    Allocate(eparams%s(eparams%nkmax),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Cannot allocate space for eparams%s'
      Stop              
    End If
    
    !** ..and the si array. This hold imaginary values for the k-space sum
    Allocate(eparams%si(eparams%nkmax),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Cannot allocate space for eparams%si'
      Stop              
    End If

  End Subroutine ewald_init

  !----------------------------------------------------------------------------
  ! Calculates the structure factor for a given configuration. It requires
  ! the vectors of all ionic atoms (atvecs) and their respective charges 
  ! (charges). Also, the simulation cell box size must be passed (ell), 
  ! which is used to calculate the kappa parameter (sometimes called alpha).
  !----------------------------------------------------------------------------
  Subroutine ewald_initCalc(atvecs,charges,ell)
    Type(VecType), Dimension(:), Intent(In)   :: atvecs
    Real(Kind=RDbl), Dimension(:), Intent(In) :: charges
    Real(Kind=RDbl), Dimension(:), Intent(In) :: ell

    !** Get the minimum edge length of the simulation cell.
    eparams%length = Minval(ell)

    If (eparams%useSF) Then
      !** Call the structure factor routine. This calculates the k space
      !** vectors we will use later when computing the Ewald summation. If
      !** the coordinates of the ions are fixed, this needs to be calculated
      !** only once. If they change, it must be recalculated.
      Call ewald_sfactor(atvecs,charges,ell)
    Else
      !** Call the kvec routine. It calculates the k vectors for the
      !** reciprocal space summation. Note this does not use the structure
      !** factor, so the overall ewald calculation will be slower due to 
      !** the double sum.
      If (firstKvec) Then
        Call ewald_kvec(ell)
         firstKvec = .False.
      End If
    End If

!!$    !MDEBUG
!!$    !** Report the number of ks
    Write(0,*) "Total number of K vectors calculated: ",eparams%nk
    Write(0,*) "Total number of atoms used for K vec: ",Size(atvecs,1)

  End Subroutine ewald_initCalc

  !----------------------------------------------------------------------------
  ! Calculates the potential at atvec1 with q=charge1 due to the presence of
  ! of charges2 at atvecs2. Returns the potential and forces at atvec1 in
  ! pot and fvec1, respectively, and the forces at atvecs2. Potential is 
  ! returned in unit of kcal/mol, and forces in Ang/ps^2. 
  !
  ! If only the potential at a point (not an ion center) is needed, this
  ! routine may still be used IF selfTerm=.False., atvec1 is the point, and 
  ! charge1=1.0. This would be useful for creating pretabulated electrostatic
  ! maps.
  !----------------------------------------------------------------------------
  Subroutine ewald_getInteractionSngl(sparams,atvec1,atvecs2,charge1, &
      charges2,skip,selfTerm,pot,fvec1,fvec2)
    Type(SimCell_Params), Intent(In) :: sparams ! Simulation cell info
    Type(VecType), Intent(In) :: atvec1         ! Position of ion 1
    Type(VecType), Dimension(:), Intent(In) :: atvecs2 ! Coords of other ions
    Real(Kind=RDbl), Intent(In) :: charge1      ! Charge on ion 1
    Real(Kind=RDbl), Dimension(:), Intent(In) :: charges2 ! Charges on others
    Integer, Dimension(:), Intent(In) :: skip   ! Entries in atvecs2 to skip
    Real(Kind=RDbl), Intent(Out) :: pot         ! Total potential, kcal/mol
    Logical, Intent(In) :: selfTerm             ! If TRUE, calc self term
    Real(Kind=RDbl) :: rpot, kpot, spot         ! Real, recip, and self pot
    Type(VecType), Intent(Out), Optional :: fvec1      ! Force on ion 1
    Type(VecType), Dimension(:), Intent(Out), Optional :: fvec2 ! Force on ions
    Real(Kind=RDbl), Dimension(1,Size(atvecs2,1)) :: dist2  ! Ion separation^2
    Type(VecType), Dimension(1,Size(atvecs2,1)) :: sepvec   ! Ion sep. vectors
    Logical :: getf  ! If TRUE, calculate forces
    Integer :: natoms2  ! Number of entries in atvecs2

    !** Check for forces. If fvec1 is present, calculate the forces.
    getf = .False.
    If (Present(fvec1)) getf = .True.

    !** Zero the real, reciprocal, self term, and total potentials and forces
    pot = zero
    rpot = zero
    kpot = zero
    spot = zero
    If (Present(fvec1)) fvec1 = zero
    If (Present(fvec2)) fvec2 = VecType(zero)
    
    !** Get the number of entries in atvecs2
    natoms2 = Size(atvecs2,1)

    !** Get the separation distance between the two sets of coords.
    !** These distances should be the minimum image distances.
    Call simcell_minimage(sparams,(/atvec1/),atvecs2,sepvec,dist2)

    If (getf) Then
      !** Call real space sum. It will loop over all the separation vectors
      !** to calculate the total potential and forces. Energies are returned
      !** in units of e^2/Ang, and forces in kcal/Ang/mol.
      If (Present(fvec2)) Then
        !** Get the forces on fvec2 as well.
        Call ewald_rwald(sepvec(1,1:natoms2),dist2(1,1:natoms2),charge1, &
            charges2,skip,rpot,fvec1,fvec2)
      Else
        !** We only need the forces on fvec1
        Call ewald_rwald(sepvec(1,1:natoms2),dist2(1,1:natoms2),charge1, &
            charges2,skip,rpot,fvec1)
      End If

      !** Get recipricol space sum. This uses the k vectors that were
      !** calculated in the ewald_initCalc routine. The force calculated
      !** is only on atvec1 since this potential is only due to the 
      !** background charge.
      If (eparams%useSF) Then
        Call ewald_skwald(atvec1,charge1,kpot,fvec1)
      Else
        Call ewald_kwald(atvec1,atvecs2,charge1,charges2,kpot,fvec1)
      End If

    Else
      !** No forces are needed. Just calculate the total potential.
      Call ewald_rwald(sepvec(1,1:natoms2),dist2(1,1:natoms2),charge1, &
          charges2,skip,rpot)
      If (eparams%useSF) Then
        Call ewald_skwald(atvec1,charge1,kpot)
      Else
        Call ewald_kwald(atvec1,atvecs2,charge1,charges2,kpot)
      End If
    End If
    
    !** Sum the real and reciprocal space potentials.
    pot = rpot + kpot


    !** Get the self term if it is requested, and subtract it from the
    !** total potential. There is no force contribution from the self
    !** term UNLESS the charges on the ions are moving during the 
    !** simulation (which is almost never, so don't worry about it).
    If (selfTerm) Then
      spot = ewald_selfterm(charge1)
      pot = pot - spot
    End If

    !** The potential returned should be in terms of kcal/mol, but are in
    !** e^2/Ang. Multiply by the conversion factor. Forces are already in
    !** kcal/mol/Ang. See the parameter definition of tokcal at the top
    !** of this file for more information about the conversion factor.
    pot = pot*tokcal

    !** Debugging potentials
    epots(1) = epots(1) + rpot*tokcal
    epots(2) = epots(2) + kpot*tokcal
    epots(3) = epots(3) + spot*tokcal

  End Subroutine ewald_getInteractionSngl

  !----------------------------------------------------------------------------
  ! Calculates the potential at atvec1 with q=charge1 due to the presence of
  ! of charges2 at atvecs2. Returns the potential and forces ONLY at atvec1.
  ! We NEGLECT the self term here because we assume this is used for making
  ! an electrostatic map. If you use it for something else, you'll have to 
  ! add a flag to include the self term.
  !----------------------------------------------------------------------------
  Subroutine ewald_getInteractionHOT(sparams,atvec1,atvecs2, &
      charge1,charges2,pot,hot,success)
    Type(SimCell_Params), Intent(In) :: sparams
    Type(VecType), Intent(In) :: atvec1
    Type(VecType), Dimension(:), Intent(In) :: atvecs2
    Real(Kind=RDbl), Intent(In) :: charge1
    Real(Kind=RDbl), Dimension(:), Intent(In) :: charges2
    Real(Kind=RDbl), Intent(Out) :: pot
    Real(Kind=RDbl) :: rpot, kpot, spot, kpot1, kpot2
     Real(Kind=RDbl), Intent(Out), Dimension(:) :: hot
    Logical, Optional :: success

    Real(Kind=RDbl), Dimension(1,Size(atvecs2,1)) :: dist2
    Type(VecType), Dimension(1,Size(atvecs2,1)) :: sepvec
    Integer :: natoms2
    Integer :: i
    Real(Kind=RDbl) :: lowcut, lowcut2

    !** Zero the potential and forces
    pot = zero
    rpot = zero
    kpot1 = zero
    kpot2 = zero
    spot = zero
    hot(:) = zero
    success=.False.

    lowcut = LOW_DIST_CUTOFF
    lowcut2 = lowcut*lowcut

    !** Get the number of atoms in atvecs2
    natoms2 = Size(atvecs2,1)

    !** Get the separation distance between the two sets of coords.
    !** We use the minimum image convention.
    Call simcell_minimage(sparams,(/atvec1/),atvecs2,sepvec,dist2)

    !** Check if we meer low cut condition, so to avoid overflows

    Do i=1, Size(atvecs2,1)
       If(dist2(1,i)<lowcut2) Return
    End Do

    !** Call real space sum
    Call ewald_rwaldHOT(sepvec(1,1:natoms2),dist2(1,1:natoms2),charge1, &
        charges2,rpot,hot)
    !** Get recipricol space sum
    Call ewald_skwaldHOT(atvec1,charge1,kpot,hot)
    
    !** The molecule has not been included in the reciprical sum
    !** so we don't need to subtract the self term. 
    pot = pot + rpot + kpot

    !** The potential returned should be in terms of kcal/mol
    !** Multiply by the conversion factor
    pot = pot*tokcal

    !** So should the HOT. Multiply by the conversion factor.

   hot=hot*tokcal
   hot(1)=pot

    success=.True.

  End Subroutine ewald_getInteractionHOT

  !----------------------------------------------------------------------------
  ! Gets the total interaction between atvecs1 and atvecs2. It returns the
  ! potentials (pot) and optionally the forces (fvec) between atvecs1 and
  ! atvecs2. The entries in pot and fvec correspond to the entries in atvecs1.
  !----------------------------------------------------------------------------
  Subroutine ewald_getInteractionMult(sparams,atvecs1,atvecs2,&
      charges1,charges2,pot,fvec)
    Type(SimCell_Params), Intent(In)           :: sparams
    Type(VecType), Dimension(:), Intent(In)    :: atvecs1, atvecs2
    Real(Kind=RDbl), Dimension(:), Intent(In)  :: charges1, charges2
    Real(Kind=RDbl), Dimension(:), Intent(Out) :: pot
    Type(VecType), Intent(Out), Dimension(:), Optional :: fvec
    Real(Kind=RDbl) :: kpot, rpot
    Type(VecType), Dimension(Size(atvecs1,1),Size(atvecs2,1)) :: sepvec
    Real(Kind=RDbl), Dimension(Size(atvecs1,1),Size(atvecs2,1)) :: dist2
    Logical :: getf
    Integer :: natoms2, a
    Real(kind=RDbl), Dimension(Size(charges1,1)) :: spot

    !** MDEBUG
    Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
        "Sorry, this routine does not work yet."
    Stop

!!$    !** Size of atvecs2
!!$    natoms2 = Size(atvecs2,1)
!!$
!!$    !** Zero the potential and forces
!!$    pot = zero
!!$    rpot = zero
!!$    kpot = zero
!!$    spot = zero
!!$    getf = .False.
!!$    If (Present(fvec)) getf = .True.
!!$    If (getf) fvec = zero
!!$
!!$    !** Get the separation distance between the two sets of coords
!!$    Call simcell_minimage(sparams,atvecs1,atvecs2,sepvec,dist2)
!!$
!!$    !** Get the self terms
!!$    spot = ewald_selfterm(charges1)
!!$
!!$    !** Loop through the atoms
!!$    Do a = 1, Size(atvecs1,1)
!!$
!!$      If (getf) Then
!!$        !** Call real space sum
!!$        Call ewald_rwald(sepvec(a,1:natoms2),dist2(a,1:natoms2),charges1(a), &
!!$            charges2,rpot,fvec)
!!$        !** Get recipricol space sum
!!$        Call ewald_skwald(atvecs1(a),charges1(a),kpot,fvec)
!!$      Else
!!$        Call ewald_rwald(sepvec(a,1:natoms2),dist2(a,1:natoms2),charges1(a), &
!!$            charges2,rpot)
!!$        Call ewald_skwald(atvecs1(a),charges1(a),kpot)
!!$      End If
!!$
!!$      !** sum the potential
!!$      !pot = pot + rpot + kpot - spot(a)
!!$      !** In gcmc, there is no self term since the molecule being
!!$      !** inserted is not included in the original reciprical sum
!!$      pot = pot + rpot + kpot
!!$      !MDEBUG
!!$      Write(*,*) "rpot, kpot ",rpot,kpot
!!$
!!$    End Do
!!$
!!$    !** Convert the potential to kcal/mol
!!$    !** multiply by the conversion factor to change from
!!$    !** eV^2/Ang to kcal/mol
!!$    pot = pot*tokcal
!!$
!!$    !** Covert the forces too
!!$    fvec = fvec*tokcal

  End Subroutine ewald_getInteractionMult

  !----------------------------------------------------------------------------
  ! Returns the potential and force on atom1 due to a core/shell partner 
  ! that was included in the reciprical space sum. These values should be
  ! SUBTRACTED from the final potential and forces to obtain the correct
  ! value. atom1 is located at atvec1 and its corresponding core or shell
  ! partner is located at atvec2.
  !----------------------------------------------------------------------------
  Subroutine ewald_getCoreShellCorr(sparams,atvec1,atvec2,q1,q2,pot,fvec)
    Type(SimCell_Params), Intent(In)     :: sparams ! Simulation cell info
    Type(VecType), Intent(In)            :: atvec1, atvec2 ! Shell/core coords
    Real(Kind=RDbl), Intent(In)          :: q1, q2         ! Shell/core charges
    Real(Kind=RDbl), Intent(Out)         :: pot           ! Potential, kcal/mol
    Type(VecType), Intent(Out), Optional :: fvec          ! Force, Ang/ps^2
    Real(Kind=RDbl) :: dist2, r      ! Distance (r) and distance^2 between ions
    Real(Kind=RDbl) :: kr            ! kappa*r
    Real(Kind=RDbl) :: erfkr         ! Complementary error function
    Real(Kind=RDbl) :: derfkr        ! Derivative of the complementary erf
    Type(VecType)   :: sepvec        ! Separation vector between ions 
    Type(VecType)   :: unitvec       ! sepvec/r

    !** Zero the potential and force
    pot = 0.0_RDbl
    If (Present(fvec)) fvec = 0.0_RDbl

    !** Check to see if we need a correction. If the charge on either
    !** is zero, then no correction is needed because it will not have
    !** contributed to the sum.
    If ((sqrt(q1*q1) <=zeroTolerance).Or.(sqrt(q2*q2) <=zeroTolerance)) Return

    !** Calculate the distance between the shell/core pair. Use the
    !** minimum image convention to calculate the distance.
    Call simcell_minimage(sparams,atvec1,atvec2,sepvec,dist2)

    !** Calculate the potential. If sqrt(dist2) == 0, then we use the
    !** familiar self-term form of a gaussian centered at r = 0.
    If (Sqrt(dist2) <= zeroTolerance) Then
      !** Note that the term is almost identical to the self term for a charge
      !** q1, but multiply by q1*q2 rather than q1*q1.
      pot = ewald_selftermSngl(q1)/q1*q2
      !** There is no correction to the force if r = 0.
    Else
      !** We are a small, non-zero distance away so the extra 
      !** charge is a gaussian distribution not centered at r = 0.
      !** This is similar to the real space potential term, but we use
      !** the error function instead of the complementary error function.
      r = Sqrt(dist2)
      kr = eparams%kappa*r
      erfkr = (1.0_RDbl - erfc(kr))
      pot = q1*q2/r*erfkr
      !** Calculate the force contribution (see ewald_rwaldFvec for deriv.)
      If (Present(fvec)) Then
        unitvec = sepvec/r

        !** Calculate the value of the error function derivative,
        !** d/dr(erf(x)) = 2/sqrt(pi) * Exp(-x*x) * d(x)/dr 
        derfkr = 2.0_RDbl/Sqrt(Pi) * Exp(-kr*kr) * eparams%kappa
        
        !** Force = -grad(potential)
        fvec = (-1.0_RDbl)*unitvec*q1*q2*(-erfkr/(r*r) + derfkr/r)
      End If
    End If

    !** Change to kcal/mol from e^2/A
    pot = pot*tokcal
    If (Present(fvec)) fvec = fvec*tokcal

    !** Debugging potential
    epots(4) = epots(4) + pot*tokcal

  End Subroutine ewald_getCoreShellCorr

  !----------------------------------------------------------------------------
  ! Calculates the structure factor for the ewald sum. It requires the 
  ! position vectors atvecs of all the atoms in the system, their 
  ! corresponding charges (in units of e), and the edge lengths of the
  ! simulation box.
  !
  ! The structure factor is used because it converts the Ewald summation from
  ! a double sum over i and j to a single summation over i, so it scales with
  ! N instead of N^2.
  !----------------------------------------------------------------------------
  Subroutine ewald_sfactor(atvecs,charges,ell)
    Type(VecType), Intent(In), Dimension(:) :: atvecs
    Real(Kind=RDbl), Intent(In), Dimension(:) :: charges
    Real(Kind=RDbl), Dimension(:), Intent(In) :: ell

    Real(Kind=RDbl) :: ialpha, ibeta, igamma
    Real(Kind=RDbl) :: kalpha, kbeta, kgamma
    Real(Kind=RDbl) :: ikap, ikapPiSq, iVolPi
    Real(Kind=RDbl) :: ka2pi, kb2pi, kg2pi
    Real(Kind=RDbl) :: ksqx, ksqy, ksqz, ksq
    Real(Kind=RDbl) :: factor
    Integer :: kmax, x, y, z, j

    !** Counts the number of accepted k vectors.
    j = 0

    !** Set the number of k vectors to zero. 
    !** The actual number will vary with the coordinates passed and
    !** the kappa and cutoff values.
    eparams%nk = 0

    !** set kmax. This is the number of "boxes" to loop over in each
    !** of the three directions (x,y,z)
    kmax = eparams%kmax

    !** The values of the simulation cell edge lengths
    ialpha = one/ell(1)
    ibeta  = one/ell(2)
    igamma = one/ell(3)

    !** Set kappa. It has units of 1/Length.
    eparams%kappa = eparams%kappaParam/Min(ell(1),ell(2),ell(3))

    !** Inverse of kappa. Units of Length.
    ikap = one/eparams%kappa

    !** (Pi/kappa)^2. Units of Length^2
    ikapPiSq = (pi*ikap)**2

    !** 1/Volume * 1/Pi. Units of 1/Length^3
    iVolPi = ialpha*ibeta*igamma/pi

    !** Loop over the reciprocal space images
    Do x = -kmax,kmax
      
      kalpha = Real(x,RDbl)*ialpha     ! Units of 1/Length
      ka2pi = kalpha*twopi             ! Units of 1/Length
      ksqx = kalpha**2                 ! Units of 1/Length^2
 
      Do y = -kmax,kmax
        
        kbeta = Real(y,RDbl)*ibeta
        kb2pi = kbeta*twopi
        ksqy = kbeta**2
                
        Do z = -kmax,kmax

          !** We skip calculating k vectors for the central image
          If (x == 0 .And. y == 0 .And. z == 0) Cycle
          kgamma = Real(z,RDbl)*igamma
          kg2pi = kgamma*twopi
          ksqz = kgamma**2
          ksq = ksqx + ksqy + ksqz    ! Units of 1/Length^2
          factor = iVolPi*Exp(-ikapPiSq*ksq)/ksq
          !** Check the cutoff. If it's above, then add it to the kvecs
          If (factor > eparams%locut) Then
            j = j+1
            Call ewald_sfSum(atvecs,charges,(/ka2pi,kb2pi,kg2pi/), &
                factor)
          End If

        End Do
      End Do
    End Do

  End Subroutine ewald_sfactor

  !----------------------------------------------------------------------------
  ! Calculates the structure factor sum for the given reciprocal space k
  ! vectors. It loops over all the ions located at atvecs with charge values
  ! charges. The real and imaginary space sums are stored separately.
  !----------------------------------------------------------------------------
  Subroutine ewald_sfSum(atvecs,charges,kvals,factor)
    Type(VecType), Intent(In), Dimension(:) :: atvecs
    Real(kind=RDbl), Intent(In), Dimension(:) :: charges
    Real(Kind=RDbl), Dimension(:), Intent(In) :: kvals
    Real(Kind=RDbl) :: smag, sum, sumi, rdotk
    Real(Kind=Rdbl), Intent(InOut) :: factor
    Integer :: i, nk

    !** Zero the sums
    sum = zero
    sumi = zero

    !** Calculate the real and imaginary parts of the structure factor
    Do i = 1, Size(atvecs,1)
      rdotk = kvals*atvecs(i)  ! Dimensionless
      sum = sum + Cos(rdotk)*charges(i) ! Units of |e|
      sumi = sumi + Sin(rdotk)*charges(i) ! Units of |e|
    End Do

    !** Verify that the value is within our desired bounds. Note that
    !** 1E-3 is hard coded. It is possible to add a second cutoff
    !** parameter to be specified, but I'm not sure it's necessary.
    smag = factor**2 * (sum**2 + sumi**2) ! Units of |e|^2/Length^2
    If (smag <= eparams%locut*1.0e-3_RDbl) Return

    !** Increment the number of k vectors calculated by 1
    eparams%nk = eparams%nk + 1 
    nk = eparams%nk

    If (nk > eparams%nkmax) Then
      !** We've exceeded the limit of kvectors. Report the error and
      !** exit.
      Write(0,'(2a,i4,a,i4)') __FILE__, ":", __LINE__, &
          "Maximum number of k vectors exceed. NKMAX=",eparams%nkmax
      Stop
    Else
      !** Store the structure factor values
      eparams%s(nk)  = sum*factor  ! Units of |e|/Length
      eparams%si(nk) = sumi*factor ! Units of |e|/Length
      eparams%k(nk)  = kvals       ! Units of 1/Length
    End If
  End Subroutine ewald_sfSum

  !----------------------------------------------------------------------------
  ! Calculates the surface term for the ewald sum.
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function ewald_surfaceTerm(eparams)
    Type(EwaldParams), Intent(In) :: eparams

    !** Warn this isn't working yet
    Write(0,'(a,i5,a)') __FILE__,__LINE__, &
        " : WARNING, the surface term is not implemented yet"

    ewald_surfaceTerm = 0.0_RDbl
  End Function ewald_surfaceTerm

  !----------------------------------------------------------------------------
  ! Calculates the self term for the Ewald summation. The self term
  ! arises from the contribution of an ion's reciprocal space (background)
  ! charge with its real space self. This must be subtracted to get the 
  ! proper value. It is returned in e^2/Ang. 
  !
  ! Notes: The value of the self-term does not correspond exactly to any single
  ! term in the reciprocal space sum, i.e., you won't find this value if
  ! you look at all the values calculated in ewald_kwald. It's odd, I know.
  ! In fact, it may be orders of magnitude larger than the reciprocal space
  ! potential altogether.
  !
  ! Old jibber-jabber:
  !  The equation is taken from eqn(1.7) on p.28 of :
  !  Proc. R. Soc. Lond. A., 373, 27-56, 1980
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function ewald_selftermSngl(charge)
    Real(kind=RDbl), Intent(In) :: charge

    !** Calculate the self-term potential. This is a gaussian centered
    !** at r = 0.
    ewald_selftermSngl = 2.0_RDbl*eparams%kappa/Sqrt(Pi)*charge*charge

    !** This is the previous self term calculation. It's not used here,
    !** but the code may be of use in the future.
!!$    Real(Kind=RDbl) :: a, dxi
!!$
!!$    a = eparams%kappa * eparams%length
!!$    dxi = 2.0_RDbl*a/sqrt(pi) 
!!$    ! left out in A&T because it's small: 
!!$    !n = 1
!!$    !dxi = dxi + erfc(a*n)/(1.0e0*n) + 1.0/pi*Exp(-(pi*n/a)**2) 
!!$    ewald_selftermSngl  = dxi * charge*charge/eparams%length

  End Function ewald_selftermSngl

  !----------------------------------------------------------------------------
  ! Same as selftermSngl, but takes a dimensioned array of charges
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function ewald_selftermMult(charge)
    Real(kind=RDbl), Dimension(:), Intent(In) :: charge
    Real(kind=RDbl) :: dxi
    Integer         :: i

    !** zero the value
    ewald_selftermMult = zero

    dxi = 2.0_RDbl*eparams%kappa/sqrt(Pi)

    ! left out in A&T because it's small: 
    !a = eparams%kappa * eparams%length
    !dxi = - 2.0*a/sqrt(pi)
    !n = 1
    !dxi = dxi + erfc(a*n)/(1.0e0*n) + 1.0/pi*Exp(-(pi*n/a)**2) 

    !** loop through the charges and add them to the summation
    Do i = 1, Size(charge,1)
      ewald_selftermMult  = dxi * charge(i) * charge(i) + ewald_selftermMult
    End Do

  End Function ewald_selftermMult

  !----------------------------------------------------------------------------
  ! Calculates the real part of the ewald sum between a single atom and
  ! a group of other atoms. Units returned are e^2/Ang (pot) and 
  ! kcal/mol/Ang (fvec). The values contained in the integer array skip
  ! correspond to separation vectors that we should skip using in the 
  ! potential and force calculation.
  !
  ! Note: the force vectors are NOT zeroed in this routine. Any new forces
  ! calculated here are added to those already passed in fvec1 and fvec2.
  !
  ! The real space summation is
  !
  ! U = Sum(over i/=j to N) zi * zj / rij * erfc(kappa*rij)
  !----------------------------------------------------------------------------
  Subroutine ewald_rwaldFvec(sepvec,dist2,charge1,charges2,skip, &
      pot,fvec1,fvec2)
    Real(Kind=RDbl), Intent(In)               :: charge1
    Real(Kind=RDbl), Dimension(:), Intent(In) :: charges2
    Integer, Dimension(:), Intent(In)         :: skip
    Type(VecType), Dimension(:), Intent(In)   :: sepvec
    Real(Kind=RDbl), Dimension(:), Intent(In) :: dist2
    Real(Kind=RDbl), Intent(Out)      :: pot
    Type(VecType), Optional, Intent(InOut) :: fvec1
    Type(VecType), Optional, Dimension(:), Intent(InOut) :: fvec2

    Integer :: i
    Real(Kind=RDbl) :: ri, ri2, kr
    Real(Kind=RDbl) :: kappa   ! Gaussian distribution width
    Real(Kind=RDbl) :: erfckr  ! complementary error func, erfc(k*r)
    Real(Kind=RDbl) :: dist    ! mag(sepvec(i))
    Real(Kind=RDbl) :: derfkr  ! derivative of the error func
    Type(VecType)   :: unitvec ! Direction of the force, sepvec(i)/dist
    Type(VecType)   :: force   ! the force vector calculated

    !** zero the potential energy
    pot = zero
 
    !** Define kappa for neatness
    kappa = eparams%kappa

    !** Loop through r to calculate the potential and 
    !** (optionally) the force between atom 1 and atoms 2
    Do i = 1,Size(sepvec,1)

      !** Check to see if we should skip this entry. It may be a
      !** core/shell partner to ion 1, or may be a copy of ion 1, or 
      !** something else.
      If (findint(skip,i) /= 0) Cycle

      !** Get the distance
      dist = Sqrt(dist2(i))

      !** Check for a spherical cutoff, if it has been used. If the
      !** cutoff value is negative, it is ignored.
      If ((dist > eparams%rcut).And.(eparams%rcut > 0.0_RDbl)) Cycle

      !** Complementary error function value 
      kr = kappa*dist
      erfckr = erfc(kr)
      
      !** Define 1/r
      ri = 1.0_RDbl/dist

      !** Calculate the potential
      pot = pot + charge1*charges2(i)*erfckr*ri

      !** Calculate fvec if necessary. The derivative is
      !** d(pot)/dr = z1*z2*(-1/r^2 + erf(k*r)/r^2 - 1/r*d(erf(k*r))/dr)
      !**           = z1*z2*(-erfc(k*r)/r^2 - 1/r*d(erf(k*r))/dr)
      If (Present(fvec1)) Then
        !** Get the unit vector for the force direction
        unitvec = sepvec(i)/dist

        !** Define 1/r^2
        ri2 = ri*ri

        !** Calculate the value of the error function derivative,
        !** d/dr(erf(x)) = 2/sqrt(pi) * Exp(-x*x) * d(x)/dr 
        derfkr = 2.0_RDbl/Sqrt(Pi) * Exp(-kr*kr) * kappa

        !** Force = -grad(potential)
        force = (-1.0_RDbl)*unitvec*charge1*charges2(i)*(-erfckr*ri2-derfkr*ri)

        !** Convert the force to kcal/mol/Ang and add the force
        !** to the total force vector.
        fvec1 = fvec1 + force*tokcal
        fvec2(i) = fvec2(i) - force*tokcal
      
      End If
    End Do

  End Subroutine ewald_rwaldFvec

  !----------------------------------------------------------------------------
  ! Calculates the real part of the ewald sum between a single atom and
  ! a group of other atoms and the higher order terms (HOT). Remember, 
  ! the HOT are the derivatives, so Force = -HOT(1).
  !----------------------------------------------------------------------------
  Subroutine ewald_rwaldHOT(sepvec,dist2,charge1,charges2, &
      pot,hot)
    Real(Kind=RDbl), Intent(In)               :: charge1
    Real(Kind=RDbl), Dimension(:), Intent(In) :: charges2
    Type(VecType), Dimension(:), Intent(In)   :: sepvec
    Real(Kind=RDbl), Dimension(:), Intent(In) :: dist2
    Real(Kind=RDbl), Intent(Out)      :: pot
    Real(Kind=RDbl), Dimension(:), Intent(InOut) :: hot

    Integer :: i
    Real(kind=RDbl), Parameter :: risqrtpi = 0.5641896_RDbl
    Real(Kind=RDbl) :: ri1, ri2, ri3, ri4, ri5, ri6, kappa3, ri7
    Real(Kind=RDbl) :: erfkr, dist, kr, ek2r2, temp

    !** zero the sums, etc.
    pot = zero

    !** Loop through r to calculate the potential and 
    !** (optionally) the force between atom 1 and atoms 2
    Do i = 1,Size(sepvec,1)
      dist = Sqrt(dist2(i))
!!$      If (dist < SELF_LOCUT) Then
!!$        eparams%self = .True.
!!$        Cycle
!!$      End If
      kr = eparams%kappa*dist
      erfkr = erfc(kr)
      ri1 = one/dist
      pot = pot + charge1*charges2(i)*erfkr*ri1

      !** Calculate the HOT if necessary
      ri2 = ri1*ri1
      ri3 = ri2*ri1
      ek2r2 = Exp(-kr*kr)
      temp = -ri3*erfkr - 2.0_RDbl*eparams%kappa*ri2*risqrtpi*ek2r2 

!      hot(1) = pot
      hot(2) = hot(2) + sepvec(i)%comp(1)*(charges2(i)*temp)*charge1
      hot(3) = hot(3) + sepvec(i)%comp(2)*(charges2(i)*temp)*charge1
      hot(4) = hot(4) + sepvec(i)%comp(3)*(charges2(i)*temp)*charge1
      
      ri4 = ri3*ri1
      ri5 = ri4*ri1
      kappa3 = eparams%kappa*eparams%kappa*eparams%kappa
      temp = 3.0_RDbl*ri5*erfkr + 6.0_RDbl*eparams%kappa*risqrtpi*ri4*ek2r2 &
          + 4.0_RDbl*kappa3*risqrtpi*ri2*ek2r2
      hot(5) = hot(5) + &
          charge1*charges2(i)*sepvec(i)%comp(1)*sepvec(i)%comp(2)*temp
      hot(6) = hot(6) + &
          charge1*charges2(i)*sepvec(i)%comp(1)*sepvec(i)%comp(3)*temp
      hot(7) = hot(7) + &
          charge1*charges2(i)*sepvec(i)%comp(2)*sepvec(i)%comp(3)*temp
      
      ri6 = ri5*ri1
      ri7 = ri6*ri1
      temp = -15.0_RDbl*ri7*erfkr -30.0_RDbl*eparams%kappa*risqrtpi*ri6*ek2r2 &
          - 20.0_RDbl*kappa3*risqrtpi*ri4*ek2r2 &
          - 8.0_RDbl*kappa3*eparams%kappa*eparams%kappa*risqrtpi*ri2*ek2r2
      hot(8) = hot(8) + charge1*charges2(i)*sepvec(i)%comp(1) &
          *sepvec(i)%comp(2)*sepvec(i)%comp(3)*temp
    End Do

  End Subroutine ewald_rwaldHOT

  !----------------------------------------------------------------------------
  ! Calculates the kwald portion using the structure factor
  ! returning the potential and forces given the separation vector r
  ! in units of e^2/Ang (pot) and kcal/mol/Ang (force)
  !----------------------------------------------------------------------------
  Subroutine ewald_skwaldFvec(r,q,pot,fvec)
    Type(VecType), Intent(In) :: r
    Real(Kind=RDbl), Intent(In) :: q
    Real(Kind=RDbl), Intent(Out) :: pot
    Type(VecType), Intent(InOut), Optional :: fvec

    Integer :: i
    Type(VecType) :: k
    Real(Kind=RDbl) :: rdotk,crdotk,srdotk,realPart,imagPart
    
    !** Initialize the potential and forces
    pot = zero

    !** Loop over the k vecs
    Do i = 1, eparams%nk
      k = eparams%k(i)      ! Units of 1/Length
      rdotk = r*k           ! Dimensionless
      crdotk = cos(rdotk)
      srdotk = sin(rdotk)
      ! sum and sumi have already been doubled to account for inversion
      ! in reciprocal space? I don't think we've doubled any results - 
      ! we actually sum over all the kvecs in reciprocal space. 
      realPart = eparams%s(i)*crdotk - eparams%si(i)*srdotk  ! |e|/Length
      imagPart = eparams%s(i)*srdotk - eparams%si(i)*crdotk  ! |e|/Length

      !** Potential = realPart*charge
      pot = pot + realPart*q          ! Units of |e|^2/Length

      !** Add to the force vector. Note k*imagPart*q has units |e|^2/Length^2
      If (Present(fvec)) fvec = fvec + k*imagPart*q*tokcal

    End Do

  End Subroutine ewald_skwaldFvec

  !----------------------------------------------------------------------------
  ! Calculates the kwald portion using the structure factor
  ! returning the potential and forces given the separation vector r
  ! Returns the HOT (higher order terms), i.e., the derivatives. Remember,
  ! Force = -HOT(1)
  !----------------------------------------------------------------------------
  Subroutine ewald_skwaldHOT(r,q,pot,hot)
    Type(VecType), Intent(In) :: r
    Real(Kind=RDbl), Intent(In) :: q
    Real(Kind=RDbl), Intent(Out) :: pot
    Real(Kind=RDbl), Dimension(:), Intent(InOut) :: hot 

    Integer :: i
    Type(VecType) :: k
    Real(Kind=RDbl) :: rdotk,crdotk,srdotk,realPart,imagPart

    !** Zero the potential
    pot = zero
    
    !** Loop over the k vecs
    Do i = 1, eparams%nk
      k = eparams%k(i)
      rdotk = r*k
      crdotk = cos(rdotk)
      srdotk = sin(rdotk)
      ! sum and sumi have already been doubled to account for inversion
      ! in reciprocal space?
      realPart = eparams%s(i)*crdotk - eparams%si(i)*srdotk
      imagPart = eparams%s(i)*srdotk - eparams%si(i)*crdotk


      pot = pot + realPart*q

      !** Calculate the HOT if necessary, i.e. when making an emap
      !** clean this up too! It's messy accessing the components
!      hot(1) = pot
      hot(2) = hot(2) - k%comp(1)*imagPart*q
      hot(3) = hot(3) - k%comp(2)*imagPart*q
      hot(4) = hot(4) - k%comp(3)*imagPart*q
      hot(5) = hot(5) - k%comp(1)*k%comp(2)*realPart*q
      hot(6) = hot(6) - k%comp(1)*k%comp(3)*realPart*q
      hot(7) = hot(7) - k%comp(2)*k%comp(3)*realPart*q
      hot(8) = hot(8) -k%comp(1)*k%comp(2)*k%comp(3)*imagPart*q

    End Do
  End Subroutine ewald_skwaldHOT
  
  !----------------------------------------------------------------------------
  ! This is the reciprocal space summation of the Ewald sum. It requires the
  ! coordinates of all the atoms in the simulation cell. Also, ewald_kvec 
  ! MUST be run before this routine.
  !----------------------------------------------------------------------------
  Subroutine ewald_kwald(atvec,atvecs2,charge1,charges2,pot,fvec)
    Type(VecType), Intent(in) :: atvec
    Type(VecType), Dimension(:), Intent(In) :: atvecs2
    Real(Kind=RDbl), Intent(In) :: charge1
    Real(Kind=RDbl), Dimension(:), Intent(In) :: charges2
    Type(VecType), Intent(InOut), Optional :: fvec
    Real(Kind=RDbl), Intent(Out) :: pot
    Integer :: natoms, i, kmax, kx, ky, kz, ksq, ksqmax, totk
    Real(Kind=RDbl) :: sum, sumi, vd, vdx, vdy, vdz
    Complex, Dimension(Size(atvecs2,1),-eparams%kmax:eparams%kmax) :: eikx, &
        eiky, eikz 
    Complex :: eikr
    Complex, Dimension(Size(atvecs2,1)) :: etemp

    !** Zero the potential. Don't zero the forces - those are summed.
    pot = 0.0_RDbl

    !** Set a kmax variable for cleanliness
    kmax = eparams%kmax

    !** Set ksqmax
    ksqmax = kmax*kmax + 8

    !** Set natoms to the size of atvecs2
    natoms = Size(atvecs2,1)

    !** Construct exp(ik . r) for all ions and k-vectors 
    !** Calculate kx, ky, kz = 0 , -1 and 1 explicitly 
    Do i = 1, natoms
      
      !***e^0 = 1
      eikx(i,0) = (1.0, 0.0)
      eiky(i,0) = (1.0, 0.0)
      eikz(i,0) = (1.0, 0.0)
      
      !***Note:  no minimum image convention applied.
      eikx(i,1) = Cmplx (Cos(twopi*(atvec%comp(1)-atvecs2(i)%comp(1))*ialpha),&
          Sin(twopi*(atvec%comp(1)-atvecs2(i)%comp(1))*ialpha))
      eiky(i,1) = Cmplx (Cos(twopi*(atvec%comp(2)-atvecs2(i)%comp(2))*ibeta), &
          Sin(twopi*(atvec%comp(2)-atvecs2(i)%comp(2))*ibeta))
      eikz(i,1) = Cmplx (Cos(twopi*(atvec%comp(3)-atvecs2(i)%comp(3))*igamma),&
          Sin(twopi*(atvec%comp(3)-atvecs2(i)%comp(3))*igamma))
      
      eikx(i,-1) = Conjg (eikx(i,1))
      eiky(i,-1) = Conjg (eiky(i,1))
      eikz(i,-1) = Conjg (eikz(i,1))
      
    End Do
    
    !** Calculate remaining kx, ky and kz by recurrence 
    !** Message from THE MAN: "This looks cool to me (2-28-91)"

    Do kx = 2,KMAX
      Do i = 1, natoms
        eikx(i, kx) = eikx(i,kx-1) * eikx(i,1)
        eikx(i,-kx) = Conjg (eikx(i,kx))
      End do
    End do
    
    Do ky = 2,KMAX
      Do i = 1, natoms
        eiky(i, ky) = eiky(i,ky-1) * eiky(i,1)
        eiky(i,-ky) = Conjg (eiky(i,ky))
      End do
    End do
    
    Do kz = 2,KMAX
      Do i = 1, natoms
        eikz(i, kz) = eikz(i,kz-1) * eikz(i,1)
        eikz(i,-kz) = Conjg (eikz(i,kz))
      End do
    End do

    !** Initialize the summation variables
    totk = 0
    vd   = 0.0e0
    vdx  = 0.0e0
    vdy  = 0.0e0
    vdz  = 0.0e0
!!$    vdxy = 0.0e0
!!$    vdxz = 0.0e0
!!$    vdyz = 0.0e0
!!$    vdxyz = 0.0e0

    !** Sum over all vectors
    Do kx = -KMAX,KMAX
      
      Do ky = -KMAX,KMAX
        
        Do i = 1, natoms
          etemp(i) = eikx(i,kx) * eiky(i,ky)
        End do
        
        Do kz = -KMAX,KMAX
          
          ksq = kx*kx + ky*ky + kz*kz
          
          If ((ksq < KSQMAX) .And. (ksq > 0)) Then
            
            totk = totk + 1
            sum  = 0.0
            sumi = 0.0
            
            Do i = 1, natoms
              eikr = etemp(i) * eikz(i,kz)
              sum  = sum + charges2(i) * Real(eikr)
              sumi = sumi + charges2(i) * Aimag(eikr)
            End Do

            vd = vd + kvec(totk)%k*sum
	
            !** Note these are the DERIVATIVES
            vdx = vdx + kvec(totk)%kx*sumi
            vdy = vdy + kvec(totk)%ky*sumi
            vdz = vdz + kvec(totk)%kz*sumi
!!$            vdxy = vdxy + kvec(totk)%kxy*sum
!!$            vdxz = vdxz + kvec(totk)%kxz*sum
!!$            vdyz = vdyz + kvec(totk)%kyz*sum
!!$            vdxyz = vdxyz + kvec(totk)%kxyz*sumi
            
          End If
        End Do        !kz loop
      End Do        !ky loop
    End Do        !kx loop
        
    !** Transfer to potential and force arrays. Force = -grad(pot)
    pot = vd*charge1
    If (Present(fvec)) fvec = fvec + charge1*tokcal*(/ -vdx, -vdy, -vdz /)

  End Subroutine ewald_kwald

  !----------------------------------------------------------------------------
  ! Calculates the k vectors for the reciprocal space portion of the Ewald
  ! summation. This subroutine uses the original formulation of the Ewald 
  ! Sum (includes a double sum). It is slower than the using the structure
  ! factor. It must be run before ewald_getinteraction.
  !----------------------------------------------------------------------------
  Subroutine ewald_kvec(ell)
    Real(Kind=RDbl), Dimension(3), Intent(In) :: ell
    Real(Kind=RDbl) :: volume, b, kappa, temp, alpha, beta, gamma
    Integer :: kmax, kx, ky, kz, ksq, ksqmax, error, totk
    Real(Kind=RDbl) :: rkx, rky, rkz, rksq
    Real(Kind=RDbl), Parameter :: twopi = 6.2831853
    
    !** Check to see if kvec is allocated. If not, do so.
    If (.Not.Allocated(kvec)) Then
      Allocate(kvec(eparams%nkmax),stat=error)
      If (error /= 0) Then
        Write(0,'(a,i5,a,i6)') __FILE__,__LINE__, &
            ": Could not allocate kvecs array of size ",eparams%nkmax
        Stop
      End If
    End If

    !** Set some parameters for the sum
    !** Set kappa. It has units of 1/Length.
    eparams%kappa = eparams%kappaParam/Min(ell(1),ell(2),ell(3))
    kappa = eparams%kappa
    alpha = ell(1)
    beta = ell(2)
    gamma = ell(3)

    !** Set the ialpha, ibeta, and igamma parameters
    ialpha = 1.0_RDbl/alpha
    ibeta = 1.0_RDbl/beta
    igamma = 1.0_RDbl/gamma
    
    b = 1.0_RDbl/(4.0_RDbl*kappa*kappa)
    volume = alpha*beta*gamma

    !** Set kmax and ksqmax
    kmax = eparams%kmax
    ksqmax = kmax*kmax + 8

    !** Loop through all the k vectors
    totk = 0
    Do kx = -kmax, kmax
      rkx = twopi*Real(kx,kind=RDbl)*1.0d0/alpha
      
      Do ky = -kmax, kmax
        rky = twopi*Real(ky,kind=RDbl)*1.0d0/beta
        
        Do kz = -kmax, kmax
          rkz = twopi*Real(kz,kind=RDbl)*1.0d0/gamma
          
          ksq = kx**2 + ky**2 + kz**2
          
          If ((ksq < KSQMAX) .And. (ksq /= 0)) Then
            
            totk = totk + 1

            If (totk > eparams%nkmax) Then
              Write(*,*) 'KVEC IS TOO SMALL'
              Stop
            Endif
            
            rksq = rkx**2 + rky**2 + rkz**2
            
            temp = 1.0_RDbl/volume * 2.0_RDbl*twopi/rksq * Exp(-rksq*b)
            
            kvec(totk)%k = temp
            kvec(totk)%kx = -rkx * temp
            kvec(totk)%ky = -rky * temp
            kvec(totk)%kz = -rkz * temp
            kvec(totk)%kxy = -rkx * rky * temp
            kvec(totk)%kxz = -rkx * rkz * temp
            kvec(totk)%kyz = -rky * rkz * temp
            kvec(totk)%kxyz = rkx * rky * rkz * temp
            
          End If
          
        End Do   !kz loop
      End Do      !ky loop
    End Do         !kx loop

    eparams%nk = totk

  End Subroutine ewald_kvec

  !----------------------------------------------------------------------------
  ! Cleans up any allocated memory
  !----------------------------------------------------------------------------
  Subroutine ewald_cleanup()
    Integer :: error
    
    !** Deallocate the k array
    Deallocate(eparams%k,stat=error)
    If (error/=0) Then
      Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
          " Could not deallocate eparams%k"
      Stop
    End If

    !** ..and the s array
    Deallocate(eparams%s,stat=error)
    If (error/=0) Then
      Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
          " Could not deallocate eparams%s"
      Stop
    End If

    !** ..and the si array
    Deallocate(eparams%si,stat=error)
    If (error/=0) Then
      Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
          " Could not deallocate eparams%si"
      Stop
    End If

  End Subroutine ewald_cleanup

  !----------------------------------------------------------------------------
  ! Returns the stored line that was read from the interaction parameters
  ! files
  !----------------------------------------------------------------------------
!!$  Character(len=lstrLen) Function ewald_getline()
  Function ewald_getline()
    Character(len=lstrLen) :: ewald_getline
    ewald_getline = eparams%line
  End Function ewald_getline

  !----------------------------------------------------------------------------
  ! Displays information about the Ewald parameters being used.
  !----------------------------------------------------------------------------
  Subroutine ewald_display(unit,nspace)
    Integer, Intent(In) :: unit, nspace
    Integer :: i
    Character(len=nspace) :: spc

    spc = ''
    Do i=1, nspace
      spc = spc//" "
    End Do

    Write(unit,'(2a)') spc,"Ewald Parameters"
    Write(unit,'(2a,f10.4)') spc,"kappaParams   : ",eparams%kappaParam
    If (eparams%rcut >= 0.0_RDbl) Then
      Write(unit,'(2a,f10.4,1x,a)') spc,"Real space cutoff : ",eparams%rcut, &
          lengthUnit
    End If
    Write(unit,'(2a,i10)')  spc,"Kmax in x, y, z   : ",eparams%kmax
    Write(unit,'(2a,i10)') spc,"Total k vec limit : ",eparams%nkmax
    If (eparams%useSF) Then
      Write(unit,'(2a)') spc,"Using the structure factor for Ewald"
      Write(unit,'(2a,e10.4,1x,a)') spc,"Factor locut : ",eparams%locut, &
          cutUnits
      Write(unit,'(2a,e10.4,1x,a)') spc,"smag locut   : ", &
          eparams%locut*1.0e-3_RDbl, "(e/Ang)^2"
    Else
      Write(unit,'(2a)') spc,"Using double summation in k-space for Ewald"
      Write(unit,'(2a,i10)') spc,"Kmax**2      : ",eparams%kmax**2
    End If

  End Subroutine ewald_display

End Module ewald
