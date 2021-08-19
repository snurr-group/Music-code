!------------------------------------------------------------------------------
! Data types and routines for using the Berendsen thermostat. It is being
! used here as it was implemented in Mitchell, P. J. and Fincham, D. "Shell
! model simulations by adiabatic dynamics," J. Phys. Condens. Matter, 5
! (1993) 1031-1038.
!
! Velocities are scaled at each step by a factor of beta, whose form is
!
!   beta = [ 1 + (Tr/Ta - 1) ( delta/tau ) ]^1/2
!
! where Tr is the required temperature, Ta is the instantaneous temperature,
! delta is the time step, and tau is the relaxation time (usually 0.4 ps).
!------------------------------------------------------------------------------
Module berendsen

  Use defaults, Only: RDbl, strLen, scaleke, kjmole_kb
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use utils, Only: stripcmnt, split, toreal, toint
  Use atom, Only: atom_getmass
  Use molecules, Only: molecules_getnthatom, molecules_isCoreShell, &
      molecules_getdof

  Implicit None

  Private
  Public :: BerendsenParams, berendsen_init, berendsen_display, &
      berendsen_rescale

  Type BerendsenParams
    Integer :: scaleStep
    Real(Kind=RDbl) :: Tref
    Real(Kind=RDbl) :: dt
    Real(Kind=RDbl) :: tau
  End Type BerendsenParams

Contains

  !----------------------------------------------------------------------------
  ! Initializes the Berendsen thermostat by reading in values from
  ! the control file
  !----------------------------------------------------------------------------
  Subroutine berendsen_init(br,unitno)
    Type(BerendsenParams), Intent(InOut) :: br
    Integer, Intent(In) :: unitno
    Character(len=strLen) :: text
    Character(len=strLen), Dimension(strLen) :: fields
    Integer :: nfields

    !** Read in the relaxation time
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    nfields = split(text,fields)
    br%tau = toreal(fields(1))

    !** Read in the scale step. This is recommended to be 1, so
    !** rescaling is every step.
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    nfields = split(text,fields)
    br%scaleStep = toint(fields(1))

  End Subroutine berendsen_init

  !----------------------------------------------------------------------------
  ! Rescales the given velocity using the Berendsen thermostat.
  !----------------------------------------------------------------------------
  Subroutine berendsen_rescale(br,mType,v,Tr,dt)
    Type(BerendsenParams), Intent(In) :: br
    Type(VecType), Dimension(:,:), Intent(InOut) :: v
    Real(Kind=RDbl), Intent(In) :: dt, Tr
    Integer, Intent(In) :: mType

    Real(Kind=RDbl) :: scale, Tinst, cMass, sMass
    Integer :: a, m, atype, shell
    Type(VecType) :: vrel, vcom
    Logical :: isCore

    !** Calculate the instantaneous temperature
    Tinst = berendsen_comTemp(v,mType)

    !** Calculate the scaling factor
    scale = Sqrt(1.0_RDbl + (Tr/Tinst - 1.0_RDbl)*dt/br%tau)
    
    !** Scale the velocities
    Do a = 1, Size(v,1)
      Do m = 1, Size(v,2)
        !** Get the atom type and mass
        atype = molecules_getnthatom(mType,a)
        cMass = atom_getmass(atype)
        !** Check to see if this is a core-shell pair. If so, we
        !** must add the COM velocity, not the core and shell 
        !** velocities individually.
        If (molecules_isCoreShell(mType,a,isCore,shell)) Then
          If (.Not.isCore) Cycle
          !** We need the COM velocity
          sMass = atom_getmass(molecules_getnthatom(mType,shell))
          vcom = scale*(v(a,m)*cMass + v(shell,m)*sMass) / (cMass + sMass)
          !** Calculate the core-shell relative velocity, which must be
          !** preserved
          vrel = v(a,m) - v(shell,m)
          !** With the center of mass velocity, we can calculate the 
          !** new velocities for the cores and shells
          v(a,m) = ((cMass+sMass)/cMass*vcom + sMass/cMass*vrel)/ &
              (1.0_RDbl + sMass/cMass)
          v(shell,m) = v(a,m) - vrel
        Else
          v(a,m) = v(a,m)*scale
        End If
      End Do
    End Do

    !MDEBUG
    Write(0,'(a,i4,2(a,f10.3),a)') __FILE__,__LINE__,": Original T=", &
        Tinst," K; New T=",berendsen_comTemp(v,mType)," K"

  End Subroutine berendsen_rescale

  !----------------------------------------------------------------------------
  ! Calculates the com temperature
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function berendsen_comtemp(vel,mType)
    Type(VecType), Dimension(:,:), Intent(In) :: vel
    Integer, Intent(In) :: mType
    Integer :: a,m, atype,partner
    Real(Kind=RDbl) :: mass, pMass, ke
    Logical :: isCore
    Type(VecType) :: v

    !** Zero the kinetic energy sum
    ke = 0.0_RDbl

    !** Loop over all the atoms and molecules
    Do m = 1, Size(vel,2)
      Do a = 1, Size(vel,1)
        !** Get the atom type and mass
        atype = molecules_getnthatom(mType,a)
        mass = atom_getmass(atype)
        !** Check to see if this is a core-shell pair. If so, we
        !** must add the COM velocity, not the core and shell 
        !** velocities individually.
        If (molecules_isCoreShell(mType,a,isCore,partner)) Then
          If (.Not.isCore) Cycle
          !** We need the COM velocity
          pMass = atom_getmass(molecules_getnthatom(mType,partner))
          v = (vel(a,m)*mass + vel(partner,m)*pMass) / (mass + pMass)
          mass = mass + pMass
        Else
          v = vel(a,m)
        End If
        !** Sum the kinetic energy
        ke = v*v*mass + ke
      End Do
    End Do

    berendsen_comtemp = ke*scaleke / (kjmole_kb* &
        Real(molecules_getdof(mType)*Size(vel,2)))

  End Function berendsen_comtemp

  !----------------------------------------------------------------------------
  ! Displays information about the Berendsen thermostat.
  !----------------------------------------------------------------------------
  Subroutine berendsen_display(br,dunit,nspc)
    Type(BerendsenParams), Intent(In) :: br
    Integer, Intent(In) :: dunit, nspc
    Character(nspc) :: spaces
    Integer :: i

    spaces = ""

    !** Set the spaced
    Do i = 1, nspc
      spaces = spaces//" "
    End Do

    !** Write information to the screen
    Write(dunit,'(2a)') spaces,"Berendsen Parameters"
    Write(dunit,'(2a,f10.2,a)') spaces,"Relaxation Time : ",br%tau," ps"
    Write(dunit,'(2a,i6,a)') spaces,"Rescale every ",br%scaleStep," steps"

  End Subroutine berendsen_display

End Module berendsen
