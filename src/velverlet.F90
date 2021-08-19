!------------------------------------------------------------------------------
! This module handles integration using the Velocity Verlet algorithm and
! all other functionality associated with this algorithm.  
!
! The Velocity Verlet algorithm:              where:                       
! x(dt) = x(0) + dt*v(0) + (1/2)dt^2*a(0)     x = positions                
! v(dt) = v(0) + dt/2*[a(x(0)) + a(x(dt))]    v = velocity                 
!                                             a = acceleration (force/mass)
! Note that the routines here split the velocity integration
! into two steps, one at x(0) and one at x(dt), separated by the interactions
! evaluation.
!
! Needed Improvements:
! 1) uses If/Else for fast/slow, consolidate by using config_getaccel
! 2) still lacking good comments
! 3) consider removing use of config_getconfig
! 4) why not store the ensemble type in the data structure?
! 5) why are the simcell coordinates explicitly moved rather than be generated?
!------------------------------------------------------------------------------

Module velverlet

  Use defaults, Only: RDbl
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use utils, Only: toupper, deallocerrdisplay
  Use config, Only: AtMolCoords, config_getnmoles, config_getconfig, &
      config_isfixed, config_display, config_dump
  Use simcell, Only: SimCell_Params
  Use subinteract, Only: Subset_Interactions, subinteract_int
  Use thermostats, Only: ThermostatInfo, thermostats_init, &
      thermostats_rescale, thermostats_display, thermostats_sampleCF, &
      thermostats_clean
  Use molecules, Only: molecules_getnatoms
  Use fixedsorbate, Only: fixedsorbate_check

  Implicit None
  Save

  Private
  Public :: VelVerletParams, velverlet_init, velverlet_integrate, &
      velverlet_integratePos, velverlet_halfIntegrateVel, &
      velverlet_integrateOnePos, velverlet_extraenergy, velverlet_display, &
      velverlet_sampleCF

  Type VelVerletParams
    Type(ThermostatInfo), Pointer :: thermostat
  End Type VelVerletParams

Contains

  !--------------------------------------------------------------------------
  ! Initializes the verlet params
  ! Requires:  vv -- Velocity Verlet data structure
  !            ensemble -- string identifying ensemble to integrate in
  !            unitno -- unit to read information from
  !            T -- optional temperature, if needed
  !--------------------------------------------------------------------------
  Subroutine velverlet_init(vv,ensemble,unitno,T)
    Type(VelVerletParams), Intent(InOut)  :: vv
    Character(*), Intent(In)              :: ensemble
    Integer, Intent(In)                   :: unitno
    Real(kind=RDbl), Intent(In), Optional :: T

    Integer         :: error

    Nullify(vv%thermostat)

    If (toupper(ensemble) == "NVT") Then
      Allocate(vv%thermostat, stat=error)
      If (error /=0) Then
        Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
            "Could not allocate memory for 'vv%thermostat'"
      End If
      If (.Not.Present(T)) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            "Must pass parameter T if using NVT ensemble"
      End If

      !** Call thermostat to read in the information
      Write(0,'(1x,2a,i4)') __FILE__," : Callinf thermo ",__LINE__
      Call thermostats_init(vv%thermostat,T,unitno)
    End If

  End Subroutine velverlet_init

  !----------------------------------------------------------------------------
  ! Writes a sample of the require control file information to unit unitno
  ! Requires:  unitno -- unit to write sample to
  !----------------------------------------------------------------------------
  Subroutine velverlet_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(2a)') '# If NVT was specified, remove this line and', &
        'include the following:'
    Call thermostats_sampleCF(unitno)

  End Subroutine velverlet_sampleCF

  !----------------------------------------------------------------------------
  ! Based on the given ensemble, the correct velocity verlet routine is called
  ! Requires:  subint -- subset interactions
  !            vv -- velocity verlet parameters
  !            simcell -- simulation cell information
  !            species -- configuration to move
  !            ensemble -- a string defining the ensemble, why?!
  !            dt -- timestep in picoseconds
  !            fast -- flag indicating interaction type to be calculated
  !            success -- indicates whether integrationw as successful
  !----------------------------------------------------------------------------
  Subroutine velverlet_integrate(subint,vv,simcell,species,ensemble,&
      dt, fast, success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(VelVerletParams), Intent(InOut)           :: vv
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Character(*), Intent(In)                       :: ensemble
    Logical, Intent(In)                            :: fast
    Real(kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(out)                           :: success
    
    Select Case (toupper(ensemble))
    Case ("NVT")
      !** Figure out which thermostat to use and call the proper routine
      If (Associated(vv%thermostat%vr) .Or. Associated(vv%thermostat%br)) Then
        Call velverlet_vrIntegrate(subint,vv,simcell,species,dt,fast,success)
      Else If (Associated(vv%thermostat%nh)) Then
        Call velverlet_nhIntegrate(subint,vv,simcell,species,dt,fast,success)
      End If

    Case ("NVE")
      !** Call the NVE integrator
      Call velverlet_nveIntegrate(subint,simcell,species,dt,fast,success)
          
    End Select

  End Subroutine velverlet_integrate

  !------------------------------------------------------------------------
  ! This conducts a step of velocity verlet integration in NVE ensemble
  ! Requires:  subint -- subset interactions
  !            simcell -- simulation cell information
  !            species -- configuration to move
  !            dt -- timestep in picoseconds
  !            fast -- flag indicating interaction type to be calculated
  !            success -- indicates whether integration was successful
  !------------------------------------------------------------------------
  Subroutine velverlet_nveIntegrate(subint,simcell,species,dt,fast,success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Logical, Intent(In)                            :: fast
    Real(kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(out)                           :: success

    Integer                    :: nspc, nmolecules, natoms,i,molec,a
    Real(kind=RDbl)            :: TKelvin
    Type(AtMolCoords), Pointer :: snglspc

    Logical :: nrg_success

    nspc = Size(species)
    success = .True.

    If (fast) Then
      Do i = 1, nspc

        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)

        !** Skip this species if it is fixed in position
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        !** Get the number of atoms and molecules
        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)
        
        Do molec = 1, nmolecules
          !** Calc the new positions for the principle coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%rp, &
              snglspc%coords(1:natoms,molec)%v,snglspc%afast(1:natoms,molec),dt)

          !** Calc the new positions for the simcell coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%r, &
              snglspc%coords(1:natoms,molec)%v,snglspc%afast(1:natoms,molec),dt)

          !** Calc the velocity at 1/2 time step
          Call velverlet_halfIntegrateVel(snglspc%coords(1:natoms,molec)%v, &
              snglspc%afast(1:natoms,molec),dt)

        End Do
      End Do

      !** Calculate the intermediate interactions, store as accelerations
      TKelvin = 300.0_RDbl
      nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
      If (.Not. nrg_success) Then
        success = nrg_success
        Return
      End If

!      Call config_dump(species,1,2,6)

      Do i = 1, nspc

        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        !** Skip if it's fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        Do molec = 1, nmolecules
          !** Calc the velocity at final time step
          Call velverlet_halfIntegrateVel(snglspc%coords(1:natoms,molec)%v, &
              snglspc%afast(1:natoms,molec),dt)
        End Do
      End Do

    Else   ! i.e., if not fast

      Do i = 1, nspc

        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        !** Skip if it's fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        Do molec = 1, nmolecules

          !** Calc the new positions for the principle coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%rp, &
              snglspc%coords(1:natoms,molec)%v,snglspc%aslow(1:natoms,molec),dt)

          !** Calc the new positions for the simcell coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%r, &
              snglspc%coords(1:natoms,molec)%v,snglspc%aslow(1:natoms,molec),dt)

          !** Calc the velocity at 1/2 time step
          Call velverlet_halfIntegrateVel(snglspc%coords(1:natoms,molec)%v, &
              snglspc%aslow(1:natoms,molec),dt)

          !** Zero the forces
          Do a = 1, natoms
            snglspc%aslow(a,molec) = 0.0_RDbl
          End Do
        End Do
      End Do

      !** Calculate the intermediate interactions, store as accelerations
      TKelvin = 300.0_Rdbl
      nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
      If (.Not. nrg_success) Then
        success = nrg_success
        Return
      End If

      Do i = 1, nspc
        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        
        !** Skip if it's fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        Do molec = 1, nmolecules
          !** Calc the velocity at final time step
          Call velverlet_halfIntegrateVel( &
              species(i)%coords(1:natoms,molec)%v, &
              species(i)%aslow(1:natoms,molec),dt)
        End Do
      End Do

    End If

  End Subroutine velverlet_nveIntegrate

  !----------------------------------------------------------------------------
  ! This conducts a step of velocity verlet integration in NVT ensemble with
  ! the Nose-Hoover thermostat
  ! Requires:  subint -- subset interactions
  !            vv -- Velocity Verlet data structure
  !            simcell -- simulation cell information
  !            species -- configuration to move
  !            dt -- timestep in picoseconds
  !            fast -- flag indicating interaction type to be calculated
  !            success -- indicates whether integrationw as successful
  !----------------------------------------------------------------------------
  Subroutine velverlet_nhIntegrate(subint,vv,simcell,species,dt,fast, success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(VelVerletParams), Intent(InOut)           :: vv
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Logical, Intent(In)                            :: fast
    Real(kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(out)                           :: success

    Integer                    :: nspc, nmolecules, natoms,i,molec
    Real(kind=RDbl)            :: TKelvin
    Type(AtMolCoords), Pointer :: snglspc

    Logical :: nrg_success

    Write (*,*) " does NH thermostat work with velverlet ?"
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop

    nspc = Size(species)
    success = .True.

    If (fast) Then
      Do i = 1, nspc

        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        !** Skip if it's fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        Do molec = 1, nmolecules
          !** Calc the new positions for the principle coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%rp, &
              snglspc%coords(1:natoms,molec)%v,snglspc%afast(1:natoms,molec),dt)

          !** Calc the new positions for the simcell coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%r, &
              snglspc%coords(1:natoms,molec)%v,snglspc%afast(1:natoms,molec),dt)

          !** Calc the velocity at 1/2 time step
          Call velverlet_halfIntegrateVel(snglspc%coords(1:natoms,molec)%v, &
              snglspc%afast(1:natoms,molec),dt)

        End Do
      End Do

      !** Calculate the intermediate forces
      TKelvin = 300.0_Rdbl
      nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
      If (.Not. nrg_success) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Interaction evaluation unsuccessful'
        Stop    
      End If

      Do i = 1, nspc
        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        !** Skip if it's fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        Do molec = 1, nmolecules
          !** Calc the velocity at final time step
          Call velverlet_halfIntegrateVel(snglspc%coords(1:natoms,molec)%v, &
              snglspc%afast(1:natoms,molec),dt)
        End Do
      End Do

    Else

      Do i = 1, nspc
        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        !** Skip if it's fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        Do molec = 1, nmolecules

          !** Calc the new positions for the principle coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%rp, &
              snglspc%coords(1:natoms,molec)%v,snglspc%aslow(1:natoms,molec),dt)

          !** Calc the new positions for the simcell coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%r, &
              snglspc%coords(1:natoms,molec)%v,snglspc%aslow(1:natoms,molec),dt)

          !** Calc the velocity at 1/2 time step
          Call velverlet_halfIntegrateVel(snglspc%coords(1:natoms,molec)%v, &
              snglspc%aslow(1:natoms,molec),dt)

        End Do
      End Do

      !** Calculate the intermediate interactions, store as accelerations
      TKelvin = 300.0_RDbl
      nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
      If (.Not. nrg_success) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Interaction evaluation unsuccessful'
        Stop    
      End If

      Do i = 1, nspc

        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        
        !** Skip if it's fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        Do molec = 1, nmolecules
          !** Calc the velocity at final time step
          Call velverlet_halfIntegrateVel( &
              species(i)%coords(1:natoms,molec)%v, &
              species(i)%aslow(1:natoms,molec),dt)
        End Do
      End Do

    End If
 
  End Subroutine velverlet_nhIntegrate

  !-------------------------------------------------------------------------
  ! This conducts a step of velocity verlet integration in NVT ensemble with
  ! the velocity rescaling thermostat.
  ! Requires:  subint -- subset interactions
  !            vv -- Velocity Verlet data structure
  !            simcell -- simulation cell information
  !            species -- configuration to move
  !            dt -- timestep in picoseconds
  !            fast -- flag indicating interaction type to be calculated
  !            success -- indicates whether integrationw as successful
  !-------------------------------------------------------------------------
  Subroutine velverlet_vrIntegrate(subint,vv,simcell,species,dt,fast,success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(VelVerletParams), Intent(InOut)           :: vv
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Logical, Intent(In)                            :: fast
    Real(kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(out)                           :: success

    Integer                    :: nspc, nmolecules, natoms,i,molec,a
    Real(kind=RDbl)            :: TKelvin
    Type(AtMolCoords), Pointer :: snglspc

    Logical :: nrg_success

    nspc = Size(species)
    success = .True.

    If (fast) Then
      Do i = 1, nspc

        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        !** Skip if it's fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        !** Call the velocity rescaling thermostat first
        Call thermostats_rescale(subint%imodel,vv%thermostat,i, &
            species(i)%coords(1:natoms,1:nmolecules)%v,dt)

        Do molec = 1, nmolecules
          !** Calc the new positions for the principle coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%rp, &
              snglspc%coords(1:natoms,molec)%v,snglspc%afast(1:natoms,molec),dt)

          !** Calc the new positions for the simcell coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%r, &
              snglspc%coords(1:natoms,molec)%v,snglspc%afast(1:natoms,molec),dt)

          !** Calc the velocity at 1/2 time step
          Call velverlet_halfIntegrateVel(snglspc%coords(1:natoms,molec)%v, &
              snglspc%afast(1:natoms,molec),dt)

        End Do
      End Do

      !** Calculate the intermediate interactions, store as accelerations
      TKelvin = 300.0_RDbl
      nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
      If (.Not. nrg_success) Then
        success = nrg_success
        Return
      End If

      Do i = 1, nspc

        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        !** Skip if it's fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        Do molec = 1, nmolecules
          !** Calc the velocity at final time step
          Call velverlet_halfIntegrateVel(snglspc%coords(1:natoms,molec)%v, &
              snglspc%afast(1:natoms,molec),dt)
        End Do
      End Do

    Else  ! i.e., not fast

      Do i = 1, nspc

        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        !** Skip if it is fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        !** Call the velocity rescaling thermostat first
        Call thermostats_rescale(subint%imodel,vv%thermostat,i, &
            snglspc%coords(1:natoms,1:nmolecules)%v,dt)

        Do molec = 1, nmolecules

          !** Calc the new positions for the principle coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%rp, &
              snglspc%coords(1:natoms,molec)%v,snglspc%aslow(1:natoms,molec),dt)

          !** Calc the new positions for the simcell coords
          Call velverlet_integratePos(snglspc%coords(1:natoms,molec)%r, &
              snglspc%coords(1:natoms,molec)%v,snglspc%aslow(1:natoms,molec),dt)

          !** Calc the velocity at 1/2 time step
          Call velverlet_halfIntegrateVel(snglspc%coords(1:natoms,molec)%v, &
              snglspc%aslow(1:natoms,molec),dt)

          !** Zero the forces
          Do a = 1, natoms
            snglspc%aslow(a,molec) = 0.0_RDbl
          End Do
        End Do
      End Do

      !** Calculate the intermediate interactions, store as accelerations
      TKelvin = 300.0_RDbl
      nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
      If (.Not. nrg_success) Then
        success = nrg_success
        Return
      End If

      Do i = 1, nspc

        !** Assign the pointer
        Call config_getconfig(species,i,snglspc)
        
        !** Skip if it is fixed
        If (config_isfixed(snglspc)) Cycle
        If (fixedsorbate_check(i)) Cycle

        nmolecules = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)

        Do molec = 1, nmolecules
          !** Calc the velocity at final time step
          Call velverlet_halfIntegrateVel( &
              species(i)%coords(1:natoms,molec)%v, &
              species(i)%aslow(1:natoms,molec),dt)
        End Do
      End Do

    End If
 
  End Subroutine velverlet_vrIntegrate

  !----------------------------------------------------------------------------
  ! Performs the full position integration of velocity verlet.
  ! x(dt) = x(0) + dt*v(0) + (1/2)dt^2*a(0)     
  ! Requires:  pos -- array of position vectors, will be changed
  !            vel -- array of velocity vectors
  !            acc -- array of acceleration vectors (force/mass) 
  !            dt -- timestep (in ps)
  !----------------------------------------------------------------------------
  Subroutine velverlet_integratePos(pos,vel,acc,dt)
    Type(VecType), Intent(InOut), Dimension(:)   :: pos
    Type(VecType), Intent(In), Dimension(:)      :: vel
    Type(VecType), Intent(In), Dimension(:)      :: acc
    Real(kind=RDbl), Intent(InOut)               :: dt

    Integer         :: i
    Real(kind=RDbl) :: dt2

    dt2 = dt*dt*0.5_RDbl
    
    Do i = 1,Size(pos,1)
      pos(i) = pos(i) + vel(i)*dt + acc(i)*dt2
    End Do

  End Subroutine velverlet_integratePos

  !----------------------------------------------------------------------------
  ! Performs the full position integration of velocity verlet on a SINGLE coord
  ! x(dt) = x(0) + dt*v(0) + (1/2)dt^2*a(0)     
  ! Requires:  pos -- position vector, will be changed
  !            vel -- velocity vector
  !            acc -- acceleration vector (force/mass) 
  !            dt -- timestep (in ps)
  !----------------------------------------------------------------------------
  Subroutine velverlet_integrateOnePos(pos,vel,acc,dt)
    Type(VecType), Intent(InOut)    :: pos
    Type(VecType), Intent(In)       :: vel
    Type(VecType), Intent(In)       :: acc
    Real(kind=RDbl), Intent(InOut)  :: dt

    Integer         :: i
    Real(kind=RDbl) :: dt2

    dt2 = dt*dt*0.5_RDbl
    pos = pos + vel*dt + acc*dt2

  End Subroutine velverlet_integrateOnePos

  !----------------------------------------------------------------------------
  ! Performs the velocity integration of velocity verlet for 1/2 time step
  ! v(dt/2) = v(0) + dt/2 * a(x(0))
  ! Requires:  vel -- array of velocity vectors
  !            acc -- array of acceleration vectors (force/mass) 
  !            dt -- timestep (in ps)
  !----------------------------------------------------------------------------
  Subroutine velverlet_halfIntegrateVel(vel,acc,dt)
    Type(VecType), Intent(InOut), Dimension(:)   :: vel
    Type(VecType), Intent(In), Dimension(:)      :: acc
    Real(kind=RDbl), Intent(In)                  :: dt

    Integer :: i

    Do i = 1,Size(vel,1)
      vel(i) = vel(i) + acc(i)*0.5_RDbl*dt
    End Do

  End Subroutine velverlet_halfIntegrateVel

  !----------------------------------------------------------------------------
  ! Performs the velocity integration of velocity verlet for the full step
  ! v(dt) = v(0) + dt/2*[a(x(0)) - a(x(dt))]   ** NOTE the minus sign!
  ! Requires:  vel -- array of velocity vectors
  !            acc -- array of acceleration vectors (force/mass) 
  !            oldAcc -- array of OLD acceleration vectors (force/mass) 
  !            dt -- timestep (in ps)
  ! This routine is not currently in use
  !----------------------------------------------------------------------------
  Subroutine velverlet_integrateVel(pos,vel,acc,oldAcc,dt)
    Type(VecType), Intent(InOut), Dimension(:)   :: pos
    Type(VecType), Intent(InOut), Dimension(:)   :: vel
    Type(VecType), Intent(InOut), Dimension(:)   :: acc, oldAcc
    Real(kind=RDbl), Intent(InOut)               :: dt

    Integer           :: i

    Do i = 1,Size(vel)
      vel(i) = vel(i) + (oldAcc(i) - acc(i))*0.5_RDbl*dt
    End Do

  End Subroutine velverlet_integrateVel

  !----------------------------------------------------------------------------
  ! Add any additional energy information to pe, ke, and their averages
  ! Requires:  vv -- Velocity Verlet data structure
  !----------------------------------------------------------------------------
  Subroutine velverlet_extraEnergy(vv,nrg,avgNrg,eType)
    Type(VelVerletParams), Intent(In) :: vv
    Real(Kind=RDbl), Intent(Out)      :: nrg, avgNrg
    Character(*), Intent(In)          :: eType
    
    !** Nothing to do right now
    nrg = 0.0_Rdbl
    avgNrg = 0.0_RDbl

  End Subroutine velverlet_extraEnergy

  !----------------------------------------------------------------------------
  ! Displays information about this integration routine
  ! Requires:  vv -- Velocity Verlet data structure
  !----------------------------------------------------------------------------
  Subroutine velverlet_display(vv,unitno)
    Type(VelVerletParams), Intent(In) :: vv
    Integer, Intent(In) :: unitno
    
    If (Associated(vv%thermostat)) Then
      Call thermostats_display(vv%thermostat,unitno,8)
    Else
      Write(unitno,'(6x,a)') 'No Extra Parameters'
    End If
  End Subroutine velverlet_display

  !----------------------------------------------------------------------------
  ! Clean the data structure
  ! Requires:  vv -- Velocity Verlet data structure
  !----------------------------------------------------------------------------
  Subroutine velverlet_clean(vv)
    Type(VelVerletParams), Intent(InOut)  :: vv

    Integer       :: error

    If (Associated(vv%thermostat)) Then
      Call thermostats_clean(vv%thermostat)
      Deallocate(vv%thermostat,STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)    
    End If

  End Subroutine velverlet_clean

End Module Velverlet


