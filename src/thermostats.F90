!------------------------------------------------------------------------------
! This module initializes the thermostats and interfaces with them
!------------------------------------------------------------------------------
Module thermostats

  Use defaults, Only: RDbl, strLen, STATS_BLOCKSIZE
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use utils, Only: stripcmnt, split, tolower, toupper, allocErrDisplay
  Use nosehoover, Only: nosehoover_init, NoseHooverParams, &
      nosehoover_display, nosehoover_sampleCF
  Use gaussts, Only: GaussTSParams, gaussts_init, gaussts_display
  Use velrescale, Only: VelRescaleParams, velrescale_init, &
      velrescale_display, velrescale_rescale, velrescale_sampleCF
  Use berendsen, Only: BerendsenParams, berendsen_init, berendsen_display, &
      berendsen_rescale
  Use stats, Only: Statistics, stats_getvalue, stats_getcavg, stats_getstd, &
      stats_init, stats_update
  Use interact, Only: Interaction_Model

  Implicit None

  Private
  Public :: ThermostatInfo, thermostats_T, thermostats_getke, & 
      thermostats_getpe, thermostats_isNH, thermostats_isGAuss, &
      thermostats_init, &
      thermostats_rescale, thermostats_display, thermostats_updateenergy, &
      thermostats_sampleCF, thermostats_clean

  Type ThermostatInfo
    Real(kind=RDbl)                 :: T
    Type(Statistics)                :: ke
    Type(Statistics)                :: pe
    Type(NoseHooverParams), Pointer :: nh
    Type(VelRescaleParams), Pointer :: vr
    Type(BerendsenParams), Pointer  :: br
    Type(GaussTSParams), Pointer  :: gauss
  End Type ThermostatInfo

Contains

  !----------------------------------------------------------------------------
  ! Initializes the correct pointer for the given thermostat
  ! Requires:  thermostat -- thermostat data structure to initialize
  !            simT -- simulation temperature
  !            unitno -- unit number to read parameters from
  !----------------------------------------------------------------------------
  Subroutine thermostats_init(thermostat,simT,unitno)
    Type(ThermostatInfo), Intent(InOut) :: thermostat
    Integer, Intent(In)                 :: unitno
    Real(kind=RDbl), Intent(In)         :: simT

    Integer               :: error, j
    Character(len=strLen) :: name, text
    Character(len=strLen), Dimension(strLen) :: params
    Logical :: species_wise
    thermostat%T = simT

    Read(unitno,'(a)') text
    text = stripcmnt(text)
    j = split(text,params)
    name = params(1)

    !** Nullify the pointer set
    Call thermostats_nullset(thermostat)

    Select Case (Toupper(name))
      ! one thermostat for each species
    Case ("SPC_NOSEHOOVER")
      Allocate(thermostat%nh,stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      species_wise=.True.
      Call nosehoover_init(thermostat%nh,unitno,species_wise)

    Case ("NOSEHOOVER")
      Allocate(thermostat%nh,stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call nosehoover_init(thermostat%nh,unitno)


    Case ("VELOCITYRESCALE")
      Allocate(thermostat%vr,stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call velrescale_init(thermostat%vr,unitno)

    Case ("BERENDSEN")
      Allocate(thermostat%br,stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call berendsen_init(thermostat%br,unitno)

    Case ("GAUSS")
      Allocate(thermostat%gauss,stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call gaussts_init(thermostat%gauss,unitno)
    Case Default
      Write(0,'(2a,i4,3a)') __FILE__,": ",__LINE__, &
          " Unrecognized thermostat ",Trim(name)," specified in control file"
      Write(0,'(2a)') "Accepted thermostats are: ", &
          "SPC_NOSEHOOVER, NOSEHOOVER, VELOCITYRESCALE, BERENDSEN, GAUSS"
      Stop

    End Select

    Call stats_init(thermostat%pe,"Thermostat PE",STATS_BLOCKSIZE, &
        .False.,'f10.3')
    Call stats_init(thermostat%ke,"Thermostat KE",STATS_BLOCKSIZE, &
        .False.,'f10.3')

  End Subroutine thermostats_init

  !--------------------------------------------------------------------
  ! Nullify the pointer set
  ! Requires:  thermostat -- thermostat data structure 
  !--------------------------------------------------------------------
  Subroutine thermostats_nullset(thermostat)
    Type(ThermostatInfo), Intent(InOut) :: thermostat

    Nullify(thermostat%nh)
    Nullify(thermostat%vr)
    Nullify(thermostat%br)
    Nullify(thermostat%gauss)

  End Subroutine thermostats_nullset

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine thermostats_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a,t30,a)') '[NoseHoover,VelocityRescale]', &
        '# Name of the thermostat to be used'
    Write(unitno,'(a)') &
        '# If NoseHoover is specified, remove this line and include:'
    Call nosehoover_sampleCF(unitno)
    Write(unitno,'(a)') '# End NoseHoover section (remove this line)'
    Write(unitno,'(a)') &
        '# If VelocityRescale, remove this line and include:'
    Call velrescale_sampleCF(unitno)
    Write(unitno,'(a)') '# End VelocityRescale section (remove this line)'
  End Subroutine thermostats_sampleCF

  !----------------------------------------------------------------------------
  ! Returns true if the thermostat type is Nose-Hoover
  !----------------------------------------------------------------------------
  Logical Function thermostats_isNH(thermostat)
    Type(ThermostatInfo), Intent(In) :: thermostat
    
    thermostats_isNH = .False.
    If (Associated(thermostat%nh)) thermostats_isNH = .True.
  End Function thermostats_isNH

  !----------------------------------------------------------------------------
  ! Returns true if the thermostat type is Nose-Hoover
  !----------------------------------------------------------------------------
  Logical Function thermostats_isGauss(thermostat)
    Type(ThermostatInfo), Intent(In) :: thermostat
    
    thermostats_isGauss = .False.
    If (Associated(thermostat%gauss)) thermostats_isGauss = .True.
  End Function thermostats_isGauss

  !----------------------------------------------------------------------------
  ! Interfaces with the velocity rescaling routine
  !----------------------------------------------------------------------------
  Subroutine thermostats_rescale(imodel,thermostat,sorb,v,dt)
    Type(Interaction_Model), Intent(In)          :: imodel
    Type(ThermostatInfo), Intent(InOut)          :: thermostat
    Integer, Intent(In)                          :: sorb
    Type(VecType), Dimension(:,:), Intent(InOut) :: v
    Real(Kind=RDbl), Intent(In)                  :: dt

    !** Choose the correct rescaling method
    If (Associated(thermostat%vr)) Then
      Call velrescale_rescale(imodel,thermostat%vr,sorb,v,thermostat%T)
    Else If (Associated(thermostat%br)) Then
      Call berendsen_rescale(thermostat%br,sorb,v,thermostat%T,dt)
    End If

  End Subroutine thermostats_rescale

  !----------------------------------------------------------------------------
  ! Returns the thermostat set-point temperature
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function thermostats_T(thermostat)
    Type(ThermostatInfo), Intent(In) :: thermostat

    thermostats_T = thermostat%T
  End Function thermostats_T

  !----------------------------------------------------------------------------
  ! Updates the thermostat energies
  !----------------------------------------------------------------------------
  Subroutine thermostats_updateEnergy(thermostat,eType,nrg)
    Type(ThermostatInfo), Intent(InOut) :: thermostat
    Character(*), Intent(In)            :: eType
    Real(Kind=RDbl), Intent(In)         :: nrg

    Select Case(tolower(eType))
      
    Case ('ke')
      Call stats_update(thermostat%ke,nrg)
    Case ('pe')
      Call stats_update(thermostat%pe,nrg)
    End Select

  End Subroutine thermostats_updateEnergy

  !----------------------------------------------------------------------------
  ! Returns the thermostat kinetic energy type requested by sType
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function thermostats_getke(thermostat,sType)
    Type(ThermostatInfo), Intent(In) :: thermostat
    Character(*), Intent(In)         :: sType

    Select Case(tolower(sType))
      
    Case ('inst')
      thermostats_getke = stats_getvalue(thermostat%ke)
    Case ('cavg')
      thermostats_getke = stats_getcavg(thermostat%ke)
    Case ('std')
      thermostats_getke = stats_getstd(thermostat%ke)
    End Select

  End Function thermostats_getke

  !----------------------------------------------------------------------------
  ! Returns the thermostat potential energy type requested by sType
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function thermostats_getpe(thermostat,sType)
    Type(ThermostatInfo), Intent(In) :: thermostat
    Character(*), Intent(In)         :: sType

    Select Case(tolower(sType))
      
    Case ('inst')
      thermostats_getpe = stats_getvalue(thermostat%pe)
    Case ('cavg')
      thermostats_getpe = stats_getcavg(thermostat%pe)
    Case ('std')
      thermostats_getpe = stats_getstd(thermostat%pe)
    End Select

  End Function thermostats_getpe

  !----------------------------------------------------------------------------
  ! Displays information about the thermostat, i.e., the name and calls
  ! the approriate module to display the rest of the information.
  !----------------------------------------------------------------------------
  Subroutine thermostats_display(thermostat,unitno,nskip)
    Type(ThermostatInfo), Intent(In) :: thermostat
    Integer, Intent(In), Optional    :: unitno
    Integer, Intent(In)              :: nskip

    Integer               :: dUnit
    Character(len=strLen) :: frmt

    !** Use the default unit (6) if none is passed
    If (.Not.Present(unitno)) Then
      dUnit = 6
    Else
      dUnit = unitno
    End If

    !** Define the format string
    Write(frmt,'(a,i2,a)') "(",nskip,"x,2a)"

    Write(dUnit,frmt) "Thermostat Information"
    If (Associated(thermostat%nh)) Then
      If(thermostat%nh%species_wise) Then
        Write(dUnit,frmt) "Name : ","NoseHoover (species wise)"
      else
        Write(dUnit,frmt) "Name : ","NoseHoover"
      endif
      Call nosehoover_display(thermostat%nh,dUnit)
    Else If (Associated(thermostat%vr)) Then
      Write(dUnit,frmt) "Name : ","VelocityRescale"
      Call velrescale_display(thermostat%vr,dUnit,8)
    Else If (Associated(thermostat%br)) Then
      Write(dUnit,frmt) "Name : ","Berendsen"
      Call berendsen_display(thermostat%br,dUnit,8)
    Else If (Associated(thermostat%gauss)) Then
      Write(dUnit,frmt) "Name : ","Gauss"
      Call gaussts_display(thermostat%gauss,dUnit)
    End If

  End Subroutine thermostats_display

  !----------------------------------------------------------------------------
  ! Cleans up the thermostat variable by deallocating anything allocated
  ! Requires:  thermostat -- thermostat data structure 
  !----------------------------------------------------------------------------
  Subroutine thermostats_clean(thermostat)
    Type(ThermostatInfo), Intent(InOut) :: thermostat

    If (Associated(thermostat%nh)) Then
      !** Nothing
    End If
    If (Associated(thermostat%vr)) Then
      !** Nothing
    End If
    If (Associated(thermostat%br)) Then
      !** Nothing
    End If
    If (Associated(thermostat%gauss)) Then
      !** Nothing
    End If

  End Subroutine thermostats_clean

End Module thermostats

