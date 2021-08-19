!------------------------------------------------------------------------------
! This module handles integration using the Leap Frog algorithm and
! all other functionality associated with this algorithm.  
! refer to - Allen Tildesely page-231 [eqns 7.56-7.58] and
!          - Brown and Clarke Mol. Phys. 1984 vol.51 1243-1252 eqns [3.1-3.11]
!
! The Leap Frog algorithm:              where:                       
! x(0+dt)   = x(0)      + dt* v(0+dt/2)                x = positions
! v(0+dt/2) = v(0-dt/2) + dt* a(x(0))                  v = velocity 
! v(0)      = [v(0-dt/2) + v(0+dt/2)]/2 
!                                             a = acceleration (force/mass)
! Notes : 
!     - curretnly implemented for only NVT (constraint based thermostat)
!     - after the integration is done we have x(dt) v(0+dt/2) stored 
!       v(0+dt) will be calculated only ion next step, 
!       this might lead to apparent enery-non conservation, 
!------------------------------------------------------------------------------

Module leapfrog

  Use config, Only: AtMolCoords, config_getnmoles, config_getconfig, &
      config_isfixed, config_display, config_dump, config_getSpcCOMVel
  Use defaults, Only: RDbl, one,zero, scaleke, scalef, kjmole_kb, caltoj
  Use molecules, Only: molecules_getnatoms, molecules_getnsorbs, &
      molecules_getdof, molecules_atommass
  Use nemd, Only: NEMD_IS_ON, NEMDParams, nemd_updateNrg, nemd_getpe, &
      nemd_updateflux
  Use simcell, Only: SimCell_Params
  Use subinteract, Only: Subset_Interactions, subinteract_int
  Use thermostats, Only: ThermostatInfo, thermostats_init, &
      thermostats_rescale, thermostats_display, thermostats_sampleCF, &
      thermostats_clean, thermostats_isGauss, thermostats_T
  Use utils, Only: toupper, deallocerrdisplay, allocerrdisplay, tolower
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use fixedsorbate, Only: fixedsorbate_check

  Implicit None
  Save

  Private
  Public :: LeapfrogParams, halfvels, leapfrog_init, leapfrog_integrate, &
      leapfrog_extraenergy, leapfrog_display, leapfrog_sampleCF
  Type halfvels
    Type(VecType), Dimension(:,:),Pointer :: v
  End Type halfvels

  Type LeapfrogParams
    Type(ThermostatInfo), Pointer :: thermostat
    Type(halfvels), Dimension(:),Pointer :: halfv
    Logical :: half_vels_placed
  End Type LeapfrogParams

Contains

  !--------------------------------------------------------------------------
  ! Initializes the leapfrog params
  ! Requires:  lf -- Velocity Verlet data structure
  !            ensemble -- string identifying ensemble to integrate in
  !            unitno -- unit to read information from
  !            T -- optional temperature, if needed
  !--------------------------------------------------------------------------
  Subroutine leapfrog_init(lf,species,ensemble,unitno,T)
    Type(LeapfrogParams), Intent(InOut)  :: lf
    Type(AtMolCoords), Intent(In), Dimension(:) :: species
    Character(*), Intent(In)              :: ensemble
    Integer, Intent(In)                   :: unitno
    Real(kind=RDbl), Intent(In), Optional :: T

    Integer         :: error, spc, nspc, natoms, nmoles

    Nullify(lf%thermostat)

    If (toupper(ensemble) == "NVT") Then

      Allocate(lf%thermostat, stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      If (.Not.Present(T)) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            "Must pass parameter T if using NVT ensemble"
      End If

      !** Call thermostat to read in the information
      Write(0,'(1x,2a,i4)') __FILE__," : Calling thermo stat init",__LINE__
      Call thermostats_init(lf%thermostat,T,unitno)
      If (.Not.thermostats_isGauss(lf%thermostat)) Then
        Write(*,*) " leap frog needs gauss thermostat"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        stop
      Endif
    Else
      ! nothing to do
    Endif

    ! Check for nemd
    If (NEMD_IS_ON) Then
      If (associated(lf%thermostat)) Then
        If (.not.(thermostats_isgauss(lf%thermostat))) Then
          Write(*,*) "NEMD requires gauss thermostat"
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          stop
        endif
      else
        Write(*,*) "NEMD requires gauss thermostat and NVT ensemble"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        stop
      endif
    endif

    nspc=molecules_getnsorbs()
    ! they will be placed later in _inithalfvels
    lf%half_vels_placed=.False.
    Allocate(lf%halfv(nspc),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle
      If (fixedsorbate_check(spc)) Cycle
      nmoles=config_getnmoles(species,spc)
      natoms=molecules_getnatoms(spc)
      Allocate(lf%halfv(spc)%v(natoms,nmoles),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End Do

  End Subroutine leapfrog_init

  !----------------------------------------------------------------------------
  ! Writes a sample of the require control file information to unit unitno
  ! Requires:  unitno -- unit to write sample to
  !----------------------------------------------------------------------------
  Subroutine leapfrog_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(2a)') '# If NVT was specified, remove this line and', &
        'include the following:'
    Call thermostats_sampleCF(unitno)

  End Subroutine leapfrog_sampleCF

  !--------------------------------------------------------------------
  ! initializes half-vels by assuimg they are equal to current velocities
  !--------------------------------------------------------------------
  Subroutine leapfrog_inithalfvels(lf,species)
    Type(LeapfrogParams), Intent(InOut)           :: lf
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer :: nspc, spc,nmoles, natoms, a , m

    nspc=molecules_getnsorbs()
    lf%half_vels_placed=.True.
    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle
      If (fixedsorbate_check(spc)) Cycle
      nmoles=config_getnmoles(species,spc)
      natoms=molecules_getnatoms(spc)
      Do a=1,natoms
        Do m=1,nmoles
          lf%halfv(spc)%v(a,m)=species(spc)%coords(a,m)%v
        end do
      end do
    End Do

  End Subroutine leapfrog_inithalfvels

  !----------------------------------------------------------------------------
  ! Based on the given ensemble, the correct leapfrog routine is called
  ! Requires:  subint -- subset interactions
  !            lf -- leapfrog parameters
  !            simcell -- simulation cell information
  !            species -- configuration to move
  !            ensemble -- a string defining the ensemble, why?!
  !            dt -- timestep in picoseconds
  !            fast -- flag indicating interaction type to be calculated
  !            success -- indicates whether integrationw as successful
  !----------------------------------------------------------------------------
  Subroutine leapfrog_integrate(subint,lf,simcell,species,ensemble,&
      dt, fast, success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(LeapfrogParams), Intent(InOut)           :: lf
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Character(*), Intent(In)                       :: ensemble
    Logical, Intent(In)                            :: fast
    Real(kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(out)                           :: success

    If (.Not.lf%half_vels_placed) Then
      ! put some reasonable velocities in there
      Call leapfrog_inithalfvels(lf,species)
    Endif

    Select Case (toupper(ensemble))
    Case ("NVT")
      !** Figure out which thermostat to use and call the proper routine
      If (Associated(lf%thermostat%gauss)) Then
        If(NEMD_IS_ON) Then
          Call leapfrog_gaussNemdInteg(subint,lf,simcell,species,dt,success)
        Else
          Call leapfrog_gaussIntegrate(subint,lf,simcell,species,dt,success)
        Endif
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        stop
      Endif
    Case ("NVE")
      !** Call the NVE integrator
      Call leapfrog_nveIntegrate(subint,lf,simcell,species,dt,success)

    End Select

  End Subroutine leapfrog_integrate


  !----------------------------------------------------------------------------
  ! This conducts a step of Leap Frog integration in NVT ensemble with
  ! the Gaussian thermostat
  ! Requires:  subint -- subset interactions
  !            lf -- Velocity Verlet data structure
  !            simcell -- simulation cell information
  !            species -- configuration to move
  !            dt -- timestep in picoseconds
  !            success -- indicates whether integrationw as successful
  ! Note : at the end of this algorithm both KE, and PE are the ones 
  !        corresponding to previous time step 
  ! Note : species contains r(t+dt), v(t), a(t), subint contains PE(t)
  ! Note : uses only "fast" interactions
  !----------------------------------------------------------------------------
  Subroutine leapfrog_gaussIntegrate(subint, lf, simcell, species, dt, &
      success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(LeapfrogParams), Intent(InOut)           :: lf
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(out)                           :: success

    Integer                    :: nspc, nmoles, natoms,spc,molec
    Real(kind=RDbl)            :: Tk, B
    Type(AtMolCoords), Pointer :: single
    Type(VEcType) :: v0_dash
    Logical :: nrg_success, fast

    nspc = Size(species)
    success = .True.
    fast = .True.

    !** Calculate the forces at current positions
    Tk = thermostats_T(lf%thermostat)
    nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,Tk/),(/0,0,0/))
    If (.Not. nrg_success) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Interaction evaluation unsuccessful'
      Stop    
    End If

    Do spc = 1, nspc

      !** Assign the pointer
      Call config_getconfig(species,spc,single)
      nmoles = config_getnmoles(species,spc)
      natoms = molecules_getnatoms(spc)

      !** Skip if it's fixed
      If (config_isfixed(single)) Cycle
      If (fixedsorbate_check(spc)) Cycle

      !** find v_dash
      Do molec = 1, nmoles

        !** Calc the velocity at current position buy integrating only for 
        !** 1/2 t With regular accelartions Only
        Call leapfrog_halfIntegVel(lf%halfv(spc)%v(1:natoms,molec),&
            single%afast(1:natoms,molec), &
            single%coords(1:natoms,molec)%v,dt/2)
      End Do

      !** calulate streaming velocity of v_dash
      v0_dash=config_getSpcCOMVel(species,spc)

      Call leapfrog_getscale(lf,species,spc,v0_dash,B)

      Do molec = 1, nmoles

        !** Calc the velocity at 1/2 time step and current vels
        Call leapfrog_gaussIntegVel(lf%halfv(spc)%v(1:natoms,molec),&
            single%afast(1:natoms,molec), &
            single%coords(1:natoms,molec)%v,dt,B)

        !** Calc the new positions for the principle coords
        Call leapfrog_integratePos(single%coords(1:natoms,molec)%rp,&
            lf%halfv(spc)%v(1:natoms,molec),dt)

      End Do
    End Do

  End Subroutine leapfrog_gaussIntegrate


  !----------------------------------------------------------------------------
  ! This conducts a step of Leap Frog integration in NVT ensemble with
  ! the Gaussian thermostat and NEMD field
  ! Requires:  subint -- subset interactions
  !            lf -- Velocity Verlet data structure
  !            simcell -- simulation cell information
  !            species -- configuration to move
  !            dt -- timestep in picoseconds
  !            success -- indicates whether integrationw as successful
  ! Note : at the end of this algorithm both KE, and PE are the ones 
  !        corresponding to previous time step 
  ! Note : species contains r(t+dt), v(t), a(t), subint contains PE(t)
  ! Note : uses only "fast" interactions
  !----------------------------------------------------------------------------
  Subroutine leapfrog_gaussNemdInteg(subint, lf, simcell, species, dt, &
      success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(LeapfrogParams), Intent(InOut)           :: lf
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(out)                           :: success

    Integer                    :: a, m, nspc, nmoles, natoms,spc
    Real(kind=RDbl)            :: Tk, B, max_v, mass, mult
    Type(AtMolCoords), Pointer :: single
    Type(VEcType) :: v0_dash
    Logical :: nrg_success, fast

    nspc = Size(species)
    success = .True.
    fast = .True.

    !** Calculate the forces at current positions
    Tk = thermostats_T(lf%thermostat)
    nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,Tk/),(/0,0,0/))
    If (.Not. nrg_success) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Interaction evaluation unsuccessful'
      Stop    
    End If

    ! potential nrg of color field, not really needed, but helpful in checking 
    ! hamiltonian conservation
    Call nemd_updateNrg(species)

    ! update flux values for reporting on the screen
    ! n gear-6 max_ is evaluated properly
    ! here we are using approximate value, anyway
    ! this way of updating flux is Only an extra check
    max_v = 10.00 ! ang/ps
    Call nemd_updateFlux(species,simcell,max_v*dt)


    Do spc = 1, nspc


      !** Assign the pointer
      Call config_getconfig(species,spc,single)

      !** Skip if it's fixed
      If (config_isfixed(single)) Cycle

      nmoles = config_getnmoles(species,spc)
      natoms = molecules_getnatoms(spc)


      If (NEMDParams%sorb==spc) Then

        ! add nemd force
        Do a=1,natoms
          !** Update accels with NEMD force also, so we add nemd accels 
          !** to regular accels, so this is needed only once when c==1
          ! units : mass = amu
          !         scalef = (amu. ang^2/ ps^2) / (Kcal/mol) 
          !         mult = (ang^2/ps^2)/(Kcal/mol)
          ! nemd_force*mult = ang/ps^2
          mass=molecules_AtomMass(spc,a)  
          mult=(scalef*one/mass)          ! 
          Do m=1,nmoles
            ! NEMDParams%force(a) = force from color field in Kcal/mol/ang.
            single%afast(a,m) = single%afast(a,m) + (NEMDParams%force(a)*mult)
          End Do
        End Do
      Endif

      !** find v_dash
      Do m = 1, nmoles

        !** Calc the velocity at current position buy integrating only for 
        !** 1/2 t With regular accelartions Only
        Call leapfrog_halfIntegVel(lf%halfv(spc)%v(1:natoms,m),&
            single%afast(1:natoms,m), &
            single%coords(1:natoms,m)%v,dt/2)
      End Do

      !** calulate streaming velocity of v_dash
      v0_dash=config_getSpcCOMVel(species,spc)

      !** get the scaling factor B
      Call leapfrog_getscale(lf,species,spc,v0_dash,B)

      Do m = 1, nmoles

        !** Calc the velocity at 1/2 time step and current vels
        Call leapfrog_gaussIntegVel(lf%halfv(spc)%v(1:natoms,m),&
            single%afast(1:natoms,m), &
            single%coords(1:natoms,m)%v,dt,B)

        !** Calc the new positions for the principle coords
        Call leapfrog_integratePos(single%coords(1:natoms,m)%rp,&
            lf%halfv(spc)%v(1:natoms,m),dt)

      End Do
    End Do

  End Subroutine leapfrog_gaussNemdInteg


  !----------------------------------------------------------------------------
  ! calcualte the scaling factor beta (B)
  ! Requires:  lf -- Velocity Verlet data structure
  !            species -- configuration to move
  !----------------------------------------------------------------------------
  Subroutine leapfrog_getScale(lf,species,spc,v0_dash,B)
    Type(LeapfrogParams), Intent(InOut)           :: lf
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(VecType), Intent(inout) :: v0_dash
    Integer, Intent(In)                 :: spc
    Real(kind=RDbl),Intent(out) :: B

    Integer :: nmoles, natoms,a , m
    Real(kind=RDbl) :: dof,totaldof, totalKE, v2_sum, amass, ke_exp, Tk
    Type(VEcType) :: v

    dof=molecules_getdof(spc)
    nmoles=config_getnmoles(species,spc)
    natoms=molecules_getnatoms(spc)

    If (lf%thermostat%gauss%absolute) Then
      ! For thermostating wrt lab-frame,otherwise thermostat is 
      ! wrt streaming velocities
      totaldof=dof*nmoles
      v0_dash=zero
    Else
      totaldof=dof*nmoles-3 ! without com vel
    Endif

    Tk=thermostats_T(lf%thermostat)
    ! calculate kin energy wrt streaming vels
    totalKE=zero
    Do a=1,natoms
      amass=molecules_atommass(spc,a)
      v2_sum=zero
      Do m=1,nmoles
        v=species(spc)%coords(a,m)%v-v0_dash
        v2_sum=v2_sum+(v*v)
      End Do
      totalKE=totalKE+(amass)*v2_sum
    End Do
    totalKE=totalKE*scaleKe/2 ! 1/2 mv^2, also scales to KJ/mol

    ! expected kinetic energy
    ke_exp=totaldof*kjmole_kb*Tk/2

    ! beta
    B=Sqrt(ke_exp/totalKe)

  End Subroutine leapfrog_getScale





  !----------------------------------------------------------------------------
  ! This conducts a step of Leap Frog integration in NVE ensemble
  ! Requires:  subint -- subset interactions
  !            lf -- Velocity Verlet data structure
  !            simcell -- simulation cell information
  !            species -- configuration to move
  !            dt -- timestep in picoseconds
  !            success -- indicates whether integrationw as successful
  ! Note : at the end of this algorithm both Ke, and PE are the ones 
  !        corresponding to previous time step 
  ! Note : species contains r(t+dt), v(t), a(t), subint contains PE(t)
  ! Note : uses only "fast" interactions
  !----------------------------------------------------------------------------
  Subroutine leapfrog_nveIntegrate(subint, lf, simcell, species, dt, &
      success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(LeapfrogParams), Intent(InOut)           :: lf
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(out)                           :: success

    Integer                    :: nspc, nmoles, natoms,spc,molec
    Real(kind=RDbl)            :: TKelvin
    Type(AtMolCoords), Pointer :: single

    Logical :: nrg_success, fast

    nspc = Size(species)
    success = .True.
    fast = .True.

    !** Calculate the forces at current positions
    TKelvin = 300.0_Rdbl
    nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
    If (.Not. nrg_success) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Interaction evaluation unsuccessful'
      Stop    
    End If

    Do spc = 1, nspc

      !** Assign the pointer
      Call config_getconfig(species,spc,single)
      nmoles = config_getnmoles(species,spc)
      natoms = molecules_getnatoms(spc)

      !** Skip if it's fixed
      If (config_isfixed(single)) Cycle
      If (fixedsorbate_check(spc)) Cycle

      Do molec = 1, nmoles
        !** Calc the velocity at 1/2 time step and current vels
        Call leapfrog_integrateVel(lf%halfv(spc)%v(1:natoms,molec),&
            single%afast(1:natoms,molec), &
            single%coords(1:natoms,molec)%v,dt)

        !** Calc the new positions for the principle coords
        Call leapfrog_integratePos(single%coords(1:natoms,molec)%rp,&
            lf%halfv(spc)%v(1:natoms,molec),dt)

        ! maybe not required ? integrateonestep in integrate.F90 
        ! takes care of this
        !        !** Calc the new positions for the simcell coords
        !        Call leapfrog_integratePos(single%coords(1:natoms,molec)%r,&
        !            lf%halfv(spc)%v(1:natoms,molec),dt)

      End Do
    End Do

  End Subroutine leapfrog_nveIntegrate

  !----------------------------------------------------------------------------
  ! Performs the position integration of leapfrog
  ! x(dt) = x(0) + dt * v(0+dt/2)     
  ! Requires:  pos -- array of position vectors, will be changed
  !            vel -- array of velocity vectors
  !            dt -- timestep (in ps)
  !----------------------------------------------------------------------------
  Subroutine leapfrog_integratePos(pos,vel,dt)
    Type(VecType), Intent(InOut), Dimension(:)   :: pos
    Type(VecType), Intent(In), Dimension(:)      :: vel
    Real(kind=RDbl), Intent(InOut)               :: dt
    Integer         :: i
    Do i = 1,Size(pos,1)
      pos(i) = pos(i) + vel(i)*dt
    End Do
  End Subroutine leapfrog_integratePos

  !----------------------------------------------------------------------------
  ! Performs the velocity integration of velocity verlet for 1/2 time step
  ! v(dt/2) = v(-dt/2) + dt * a(x(0))
  ! v(0)  = [v(dt/2) +v(-dt/2)]/2
  ! Requires:  hvel -- array of 1/2 velocity vectors 
  !            acc -- array of acceleration vectors (force/mass) 
  !            vel -- array of velocity vectors
  !            dt -- timestep (in ps)
  !----------------------------------------------------------------------------
  Subroutine leapfrog_integrateVel(hvel,acc,vel,dt)
    Type(VecType), Intent(InOut), Dimension(:)   :: hvel
    Type(VecType), Intent(In), Dimension(:)      :: acc
    Type(VecType), Intent(InOut), Dimension(:)   :: vel
    Real(kind=RDbl), Intent(In)                  :: dt

    Type(VecType) :: tempvel ! temperorily stores v(dt/2)
    Integer :: i

    Do i = 1,Size(vel,1)
      tempvel = hvel(i) + acc(i)*dt
      vel(i) = (tempvel+hvel(i))/2
      hvel(i)=tempvel
    End Do
  End Subroutine leapfrog_integrateVel


  !----------------------------------------------------------------------------
  ! Performs the first 1/2 velocity integration of gaussian thermostat 
  ! leapfrogvelocity for 1/2 time step
  ! v_dash(0) = v(-dt/2) + dt * a(x(0))
  ! Requires:  hvel -- array of 1/2 velocity vectors 
  !            acc -- array of acceleration vectors (force/mass) 
  !            vel -- array of velocity vectors
  !            dt -- timestep (in ps)
  ! Note : dt passed here is half-step
  !----------------------------------------------------------------------------
  Subroutine leapfrog_halfIntegVel(hvel,acc,vel,dt)
    Type(VecType), Intent(InOut), Dimension(:)   :: hvel
    Type(VecType), Intent(In), Dimension(:)      :: acc
    Type(VecType), Intent(InOut), Dimension(:)   :: vel
    Real(kind=RDbl), Intent(In)                  :: dt
    Integer :: i
    Do i = 1,Size(vel,1)
      vel(i) = hvel(i) + acc(i)*dt
    End Do
  End Subroutine leapfrog_halfIntegVel


  !----------------------------------------------------------------------------
  ! Performs the full velocity integration of gaussian leapfrog
  ! v(dt/2) = v(-dt/2)*(2B-1) + dt * a(x(0)) *B
  ! v(0)  = [v(dt/2) +v(-dt/2)]/2
  ! Requires:  hvel -- array of 1/2 velocity vectors 
  !            acc -- array of acceleration vectors (force/mass) 
  !            vel -- array of velocity vectors
  !            dt -- timestep (in ps)
  !----------------------------------------------------------------------------
  Subroutine leapfrog_GaussIntegVel(hvel,acc,vel,dt,B)
    Type(VecType), Intent(InOut), Dimension(:)   :: hvel
    Type(VecType), Intent(In), Dimension(:)      :: acc
    Type(VecType), Intent(InOut), Dimension(:)   :: vel
    Real(kind=RDbl), Intent(In)                  :: dt,B

    Type(VecType) :: tempvel ! temperorily stores v(dt/2)
    Integer :: i
    Real(kind=RDbl) :: a1, a2
    a1= (2*B-1)
    a2=B*dt
    Do i = 1,Size(vel,1)
      tempvel = a1*hvel(i) + a2*acc(i)
      vel(i) = (tempvel+hvel(i))/2
      hvel(i)=tempvel
    End Do
  End Subroutine leapfrog_GaussIntegVel



  !----------------------------------------------------------------------------
  ! Add any additional energy information to pe, ke, and their averages
  ! Requires:  lf -- Velocity Verlet data structure
  !----------------------------------------------------------------------------
  Subroutine leapfrog_extraEnergy(lf,nrg,avgNrg,eType)
    Type(LeapfrogParams), Intent(In) :: lf
    Real(Kind=RDbl), Intent(Out)      :: nrg, avgNrg
    Character(*), Intent(In)          :: eType

    !** Nothing to do right now
    nrg = zero
    avgNrg = zero

    !** Add energies from color-field MD if any
    !** In color field energies are stored in kcal/mol; 
    If (NEMD_IS_ON) Then
      Select Case(tolower(eType))
      Case ('pe')
        nrg = nrg+ nemd_getpe('inst')*caltoj
        avgNrg = avgNrg + nemd_getpe('cavg')*caltoj
      Case ('ke')
        ! nothing here
      Case Default
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End Select
    Endif

  End Subroutine leapfrog_extraEnergy

  !----------------------------------------------------------------------------
  ! Displays information about this integration routine
  ! Requires:  lf -- Velocity Verlet data structure
  !----------------------------------------------------------------------------
  Subroutine leapfrog_display(lf,unitno)
    Type(LeapfrogParams), Intent(In) :: lf
    Integer, Intent(In) :: unitno

    If (Associated(lf%thermostat)) Then
      Call thermostats_display(lf%thermostat,unitno,8)
    Else
      Write(unitno,'(6x,a)') 'No Extra Parameters'
    End If
  End Subroutine leapfrog_display

  !----------------------------------------------------------------------------
  ! Clean the data structure
  ! Requires:  lf -- Velocity Verlet data structure
  !----------------------------------------------------------------------------
  Subroutine leapfrog_clean(lf)
    Type(LeapfrogParams), Intent(InOut)  :: lf

    Integer       :: error

    If (Associated(lf%thermostat)) Then
      Call thermostats_clean(lf%thermostat)
      Deallocate(lf%thermostat,STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)    
    End If

  End Subroutine leapfrog_clean

End Module Leapfrog


