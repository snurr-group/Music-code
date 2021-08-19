!------------------------------------------------------------------------------
! This module initializes and dishes out work to the correct integrator
!------------------------------------------------------------------------------
Module integrator

#ifdef USEONLY
  Use defaults, Only: RDbl
  Use velverlet, Only: VelVerletParams, velverlet_init, velverlet_integrate, &
      velverlet_extraenergy, velverlet_display
  Use gear6, Only: GearParams, gear6_init, gear6_integrate, & 
      gear6_extraenergy, gear6_display 
  Use config, Only: AtMolCoords
  Use simcell, Only: SimCell_Params
#else
  Use defaults
  Use velverlet
  Use gear6
  Use config
  Use simcell
#endif

  Implicit None
  Save

#ifdef PRIVATE
  Private
  Public :: IntegratorType, integrator_init, integrator_integrate, &
      integrator_display, integrator_extraenergy
#endif

  Type IntegratorType
    Type(GearParams), Pointer      :: gear6
!    Type(MTSParams), Pointer :: mts
    Type(VelVerletParams), Pointer :: vverlet
  End Type IntegratorType

Contains

  !----------------------------------------------------------------------------
  ! Initializes the correct pointer for the given integrator
  !----------------------------------------------------------------------------
  Subroutine integrator_init(int,intName,dt,ensemble,sorbates,unitno,T)
    Type(IntegratorType), Intent(InOut) :: int
    Character(*), Intent(In) :: intName, ensemble
    Type(AtMolCoords), Intent(In), Dimension(:) :: sorbates
    Real(kind=RDbl), Intent(In) :: dt
    Real(kind=RDbl), Intent(In), Optional :: T
    Integer, Intent(In) :: unitno
    Integer :: error
    
    Select Case (intName)
      
    Case ("VelocityVerlet")
      Allocate(int%vverlet,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            "Could not allocate memory for 'int%vverlet'"
        Stop
      End If
      If (Present(T)) Then 
        Call velverlet_init(int%vverlet,ensemble,unitno,T)
      Else
        Call velverlet_init(int%vverlet,ensemble,unitno)
      End If
      nullify(int%gear6)

    Case ("Gear6")
      Allocate(int%gear6,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            "Could not allocate memory for 'int%gear6'"
        Stop
      End If
      If (Present(T)) Then
        Call gear6_init(int%gear6,sorbates,dt,ensemble,unitno,T)
      Else
        Call gear6_init(int%gear6,sorbates,dt,ensemble,unitno)
      End If
      nullify(int%vverlet)

    Case ("MTS")
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          "Multiple time step method not yet implemented"
      Stop

    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          "Don't know how to handle integrator ",trim(intName)
      Write(0,'(a)') "Allowed Types (case sensitive): VelocityVerlet, Gear6"
      Stop

    End Select

  End Subroutine integrator_init

  !----------------------------------------------------------------------------
  ! Call the correct integrator
  !----------------------------------------------------------------------------
  Subroutine integrator_integrate(int,ensemble,simcell,sorbates,dt,fast,step)
    Type(IntegratorType), Intent(InOut) :: int
    Type(SimCell_Params), Intent(InOut) :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: sorbates
    Character(*), Intent(In) :: ensemble
    Real(kind=RDbl), Intent(In) :: dt
    Integer, Intent(In) :: step
    Logical, Intent(In) :: fast

    If (Associated(int%vverlet)) Then
      Call velverlet_integrate(int%vverlet,simcell,sorbates,ensemble,dt, &
          fast,step)
    Else If (Associated(int%gear6)) Then
      Call gear6_integrate(int%gear6,simcell,sorbates,ensemble,dt,fast,step)
    End If

  End Subroutine integrator_integrate

  !----------------------------------------------------------------------------
  ! Displays extra energy information for extended systems
  !----------------------------------------------------------------------------
  Subroutine integrator_extraEnergy(int,nrg,avgNrg,eType)
    Type(IntegratorType), Intent(In) :: int
    Real(Kind=RDbl), Intent(Out) :: nrg, avgNrg
    Character(*), Intent(In) :: eType
    
    If (Associated(int%vverlet)) Then
      Call velverlet_extraEnergy(int%vverlet,nrg,avgNrg,eType)
    Else If (Associated(int%gear6)) Then
      Call gear6_extraEnergy(int%gear6,nrg,avgNrg,eType)
    End If
  End Subroutine integrator_extraEnergy

  !----------------------------------------------------------------------------
  ! Displays information about the integrator
  !----------------------------------------------------------------------------
  Subroutine integrator_display(int,unitno)
    Type(IntegratorType), Intent(In) :: int
    Integer, Intent(In) :: unitno

    If (Associated(int%vverlet)) Then
      Write(unitno,'(6x,a)') "Velocity Verlet Integrator Information"
      Call velverlet_display(int%vverlet,unitno)
    Else If (Associated(int%gear6)) Then
      Write(unitno,'(6x,a)') "6th Order Gear Integrator Information"
      Call gear6_display(int%gear6,unitno)
    End If
    
  End Subroutine integrator_display

  Subroutine integrator_cleanup

  End Subroutine integrator_cleanup

End Module integrator
