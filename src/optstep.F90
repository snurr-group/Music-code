!------------------------------------------------------------------------------
! This module handles optimizations.  It's purpose is to perform the selected 
! optimization step(s) on a single configuration.
! 
! At the moment, it is very limited.  Currently only does steepest-descent
! and a modified velocity-verlet using the forces in storage.
!
! Need Improvements:
! 1) do forcefield evaluations in this module, consistent with other moves ??
! 2) add proper initialization from strings/control file
! 3) incorporate this module into 'moves'?
! 4) split module into different step types when there are enough
! 5) add quasi-Newton steps
! 6) add conjugate gradient steps
!------------------------------------------------------------------------------

Module optstep

  Use defaults, Only: RDbl, strLen, lstrLen, kjmole_kb, scaleke, scalepe, &
      zero, shortdashedline, TUnit, velUnit, nrgUnit
  Use file, Only: file_getunit, file_open, file_close
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, &
      toupper, allocerrdisplay, deallocerrdisplay, int2str, real2str, str2seq
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnorm, mag, vector_getunitvec
  Use config, Only: AtMolCoords, config_getnmoles, &
      config_getnatoms, config_isfixed, config_kineticEnergy, config_getaccel, &
      config_applypbcs
  Use interact, Only: Interaction_Model
  Use forcefield, Only: Forcefield_Info, forcefield_putaccels2
  Use storetop, Only: Forcefield_Results, storetop_fastforces
  Use simcell, Only: SimCell_Params, simcell_pbc, simcell_getfillsorb
  Use velverlet, Only: velverlet_integratePos, velverlet_halfIntegrateVel, &
      velverlet_integrateOnePos

  Implicit None
  Save

  Private
  Public :: Optimization_Step, optstep_init, optstep_simpleinit, &
      optstep_move, optstep_postmove, optstep_setstep, optstep_stepsize, &
      optstep_stpdescmove, optstep_clean, optstep_display, &
      Steepest_Descent, Mod_Steepest_Descent, Mod_Velocity_Verlet, &
      stpdesc_idstring, modstpdesc_idstring, modvv_idstring

  Type Optimization_Step
    Type(Steepest_Descent), Pointer       :: stpdesc
    Type(Mod_Steepest_Descent), Pointer   :: modstpdesc
    Type(Mod_Velocity_Verlet), Pointer    :: modvv
  End Type Optimization_Step

  Type Steepest_Descent
    Real(kind=RDbl)      :: max_stepsize
  End Type Steepest_Descent

  Type Mod_Steepest_Descent
    Real(kind=RDbl)      :: max_stepsize, smear_factor
  End Type Mod_Steepest_Descent

  Type Mod_Velocity_Verlet
    Logical              :: init_velocities,checkstep
    Real(kind=RDbl)      :: start_temp !** starting temperature in Kelvin
    Real(kind=RDbl)      :: timestep   !** in picoseconds
    Real(kind=RDbl)      :: maxstep    !** in Angstroms
  End Type Mod_Velocity_Verlet

  Character(len=strLen), Parameter   :: stpdesc_idstring = 'STP_DESCENT'
  Character(len=strLen), Parameter   :: modstpdesc_idstring = &
      'MOD_STP_DESCENT'
  Character(len=strLen), Parameter   :: modvv_idstring = &
      'MOD_VELOCITY_VERLET'

Contains
  !--------------------------------------------------------------------
  ! Initializes the optimization step type from control file section
  ! Requires:  stepset -- optimization step pointer set
  !            ctrl_filename -- control filename to find init info
  !--------------------------------------------------------------------
  Subroutine optstep_init(stepset,ctrl_filename)
    Type(Optimization_Step), Intent(Out)     :: stepset
    Character(*), Intent(In)                 :: ctrl_filename    

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Stop

  End Subroutine optstep_init

  !--------------------------------------------------------------------
  ! Initializes the optimization step type from given strings and 
  ! numbers.  This is the simple version.
  ! Requires:  stepset -- optimization step pointer set
  !            idstring -- string identifying which step type
  !            nums -- array of reals containing parameters for step
  !--------------------------------------------------------------------
  Subroutine optstep_simpleinit(stepset,idstring,nums)
    Type(Optimization_Step), Intent(Out)       :: stepset
    Character(*), Intent(In)                   :: idstring
    Real(kind=RDbl), Dimension(:), Intent(In)  :: nums

    Integer          :: error

    !** nullify the pointer set
    Call optstep_nullset(stepset)

    Select Case (ToUpper(idstring))
    Case (stpdesc_idstring)
      Allocate(stepset%stpdesc, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call optstep_stpdescinit(stepset%stpdesc,nums)

    Case (modstpdesc_idstring)
      Allocate(stepset%modstpdesc, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call optstep_modstpdescinit(stepset%modstpdesc,nums)

    Case (modvv_idstring)
      Allocate(stepset%modvv, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call optstep_modvvinit(stepset%modvv,nums)

    Case Default
      Write(0,'(1x,2a,i4,a)') __FILE__,": ",__LINE__, & 
          " Could not interpret optimization step type: ",Trim(idstring)
      Stop
    End Select

  End Subroutine optstep_simpleinit

  !--------------------------------------------------------------------
  ! Nullifies the pointer set.
  ! Requires:  stepset -- optimization step pointer set
  !--------------------------------------------------------------------
  Subroutine optstep_nullset(stepset)
    Type(Optimization_Step), Intent(InOut)   :: stepset

    Nullify(stepset%stpdesc)
    Nullify(stepset%modstpdesc)
    Nullify(stepset%modvv)

  End Subroutine optstep_nullset

  !------------------------------------------------------------------------
  ! Make an optimization move based on the initialized pointer in set.
  ! NOTE: will need to split this routine when other move types are
  !       added later.  also switch to taking imodel as input?
  ! Requires:  stepset -- optimization step pointer set
  !            ffresults -- results from forcefield calculation on system
  !            species -- configuration
  !            simcell -- simulation cell information
  !            careful -- can request that a careful optimization be taken
  !------------------------------------------------------------------------
  Subroutine optstep_move(stepset,imodel,species,simcell,careful)
    Type(Optimization_Step), Intent(InOut)         :: stepset
    Type(Interaction_Model), Intent(In)            :: imodel
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Logical, Intent(In)                            :: careful

    If ((Associated(stepset%stpdesc)).Or.(Associated(stepset%modstpdesc))) Then
      Call optstep_stpdescsetup(stepset,imodel%results(1),species,simcell)

    Else If (Associated(stepset%modvv)) Then
      Call optstep_modvvmove(stepset%modvv,imodel,species,simcell,careful)

    End If

  End Subroutine optstep_move

  !------------------------------------------------------------------------
  ! Perform any follow-up calculations after an accepted move.  Thought I
  ! was going to need this, but no need yet.
  ! Requires:  stepset -- optimization step pointer set
  !            ffresults -- results from forcefield calculation on system
  !            species -- configuration
  !            simcell -- simulation cell information
  !------------------------------------------------------------------------
  Subroutine optstep_postmove(stepset,ffresults,species,simcell)
    Type(Optimization_Step), Intent(InOut)         :: stepset
    Type(Forcefield_Results), Intent(In)           :: ffresults
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell

    If (Associated(stepset%stpdesc)) Then
      !** nothing
    End If

    If (Associated(stepset%modstpdesc)) Then
      !** nothing
    End If

    If (Associated(stepset%modvv)) Then
      !** Call optstep_modvvpostmove(stepset%modvv,ffresults,species,simcell)
    End If

  End Subroutine optstep_postmove

  !------------------------------------------------------------------------
  ! Setup for a steepest descent optimization step, calls move for each
  ! species.  Do not pass stepset with non-steepest descent pointers 
  ! initialized.
  ! Requires:  stepset -- optimization step pointer set
  !            ffresults -- results from forcefield calculation on system
  !            species -- configuration
  !            simcell -- simulation cell information
  !------------------------------------------------------------------------
  Subroutine optstep_stpdescsetup(stepset,ffresults,species,simcell)
    Type(Optimization_Step), Intent(InOut)         :: stepset
    Type(Forcefield_Results), Intent(In)           :: ffresults
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell

    Integer                 :: i,j,error
    Integer                 :: spc,mol,atm,nspc,nmols,natms
    Type(VecType), Dimension(:,:), Allocatable   :: force

    nspc = Size(species,1)
    Do spc = 1,nspc
    
      !** Skip the FIXED species
      If (config_isfixed(species(spc))) Cycle
    
      !** Size the forces and magnitudes arrays
      natms = config_getnatoms(species,spc)
      nmols = config_getnmoles(species,spc)
      Allocate(force(natms,nmols), STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'force')    
    
      !** Extract the forces from the forcefield results
      Call storetop_fastforces(ffresults,spc,nmols,natms,force)

      !** Make the steepest descent move for each molecule of this species
      If (Associated(stepset%stpdesc)) Then
        Do mol = 1,nmols
          Call optstep_stpdescmove(stepset%stpdesc, &
              species(spc)%coords(:,mol)%rp,force(:,mol),natms)
        End Do
      Else If (Associated(stepset%modstpdesc)) Then
        Do mol = 1,nmols
          Call optstep_modstpdescmove(stepset%modstpdesc,&
              species(spc)%coords(:,mol)%rp,force(:,mol),natms)
        End Do
      Else
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Should not be here, no support for initialized pointer'
        Stop
      End If
    
      !** Apply PBCs to the new coordinates
      Call config_applypbcs(species,simcell,(/spc,0,0/))
    
      !** Deallocate the forces and magnitudes arrays
      Deallocate(force, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'force')    
    End Do

  End Subroutine optstep_stpdescsetup

  !------------------------------------------------------------------------
  ! Make a steepest descent optimization step for one molecule
  ! Requires:  stepinfo -- optimization step pointer set
  !            coords -- cartesian coordinates (atm,mol) for each atom
  !            force -- forces (atm,mol) on each atom
  !            natms -- number of atoms in each molecule
  !------------------------------------------------------------------------
  Subroutine optstep_stpdescmove(stepinfo,coords,force,natms)
    Type(Steepest_Descent), Intent(InOut)      :: stepinfo
    Type(VecType), Dimension(:), Intent(InOut) :: coords
    Type(VecType), Dimension(:), Intent(In)    :: force
    Integer, Intent(In)                        :: natms

    Integer                 :: mol,atm
    Real(kind=RDbl)         :: magforce,maxforce,factor
    Type(VecType)           :: disp
    
    maxforce = 0.0_RDbl
    !** Get the magnitudes of the force vectors and the maximum
    Do atm = 1,natms
      magforce = vector_getnorm(force(atm))
      If (magforce > maxforce) maxforce = magforce
    End Do
    
    !** Calculate the scaling factor so that max step size is obeyed
    !** disp = factor*F  =>  maxstep = factor*|F_max|  
    factor = stepinfo%max_stepsize/maxforce
    
    !** Create the displacement vectors and apply them
    Do atm = 1,natms
      disp = factor*force(atm)
      coords(atm) = coords(atm) + disp
    End Do

  End Subroutine optstep_stpdescmove

  !------------------------------------------------------------------------
  ! Make a MODIFIED steepest descent optimization step.  This is the same
  ! as the normal steepest descent move, but some pre-specified fraction of
  ! atomic forces on each molecule are smeared across the whole molecule.
  ! The idea is to help the molecule move more as a rigid body in response
  ! to external forces.  Additionally, the maximum step restriction is 
  ! applied to individual molecules, not whole species as previously.
  ! Requires:  stepinfo -- optimization step pointer set
  !            ffresults -- results from forcefield calculation on system
  !            force -- forces (atm) on each atom
  !            natms -- number of atoms in each molecule
  !------------------------------------------------------------------------
  Subroutine optstep_modstpdescmove(stepinfo,coords,force,natms)
    Type(Mod_Steepest_Descent), Intent(InOut)     :: stepinfo
    Type(VecType), Dimension(:), Intent(InOut)    :: coords
    Type(VecType), Dimension(:), Intent(In)       :: force
    Integer, Intent(In)                           :: natms

    Integer                 :: atm
    Real(kind=RDbl)         :: magforce,maxforce,factor
    Type(VecType)           :: disp,excess,sumforce,add
    Type(VecType), Dimension(natms)  :: forcecopy

    !** Copy the forces so we don't have to change them
    forcecopy = force

    !** Extract the set fraction of force from each atom's force vector
    Do atm = 1,natms
      excess = forcecopy(atm)
      excess = stepinfo%smear_factor*excess
      forcecopy(atm) = (1.0_RDbl - stepinfo%smear_factor)* &
          forcecopy(atm)
      sumforce = sumforce + excess
    End Do
    
    !** Calculate the smeared force to be added to each atom and do it
    add = sumforce/(1.0_RDbl*natms)
    Do atm = 1,natms
      forcecopy(atm) = forcecopy(atm) + add
    End Do
    
    maxforce = 0.0_RDbl
    !** Get the magnitudes of the force vectors and the maximum
    Do atm = 1,natms
      magforce = vector_getnorm(forcecopy(atm))
      If (magforce > maxforce) maxforce = magforce
    End Do
    
    !** Calculate the scaling factor so that max step size is obeyed
    !** disp = factor*F  =>  maxstep = factor*|F_max|  
    factor = 0.0_RDbl
    If (maxforce /= 0.0_RDbl) Then
      factor = stepinfo%max_stepsize/maxforce
    End If
    
    !** Create the displacement vectors and apply them
    Do atm = 1,natms
      disp = factor*forcecopy(atm)
      coords(atm) = coords(atm) + disp
    End Do

  End Subroutine optstep_modstpdescmove

  !-------------------------------------------------------------------------
  ! Make a modified Velocity Verlet move.  This is the move type recommended
  ! by Jonsson and coworkers for the Nudged Elastic Band method.
  ! Requires:  stepinfo -- modified Velocity Verlet parameters
  !            ffresults -- results from forcefield calculation on system
  !            species -- configuration
  !            simcell -- simulation cell information
  !            careful -- if True, check step size
  !-------------------------------------------------------------------------
  Subroutine optstep_modvvmove(stepinfo,imodel,species,simcell,careful)
    Type(Mod_Velocity_Verlet), Intent(InOut)       :: stepinfo
    Type(Interaction_Model), Intent(In)            :: imodel
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Logical, Intent(In)                            :: careful

    Integer                 :: i,j,error
    Integer                 :: spc,a,m,nspc,nmols,natms
    Logical                 :: fast,damped
    Real(kind=RDbl)         :: factor,modstep,stepsize,timestep
    Type(VecType), Dimension(:,:), Allocatable   :: force

    !** Initialize the velocities if required
    If (stepinfo%init_velocities) Then
      Call optstep_modvvinitvel(stepinfo,species,imodel%results(1))
    End If

    !** Get the number of species
    nspc = Size(species,1)

    !** Update the velocities while we still have the old accelerations
    Do spc = 1,nspc
    
      !** Skip the FIXED species
      If (config_isfixed(species(spc))) Cycle

      !** Get number of atoms and moles in this species
      natms = config_getnatoms(species,spc)
      nmols = config_getnmoles(species,spc)

      !** Update each molecule separately
      Do m = 1,nmols
        Call velverlet_halfIntegrateVel(species(spc)%coords(1:natms,m)%v, &
            species(spc)%afast(1:natms,m),stepinfo%timestep)
      End Do
    End Do

    !** Convert the forces into accelerations and store in species structure
    fast = .True.  !**HACK
    Call forcefield_putaccels2(imodel%ff(1),imodel%results(1),species,fast)

    !** Project the velocities in the direction of the forces (accelerations)
    Do spc = 1,nspc
      !** Skip the FIXED species
      If (config_isfixed(species(spc))) Cycle

      Call optstep_modvvproject(species,spc,fast)
    End Do

    !** Move each species if it's not fixed and update velocity
    Do spc = 1,nspc
      !** Skip the FIXED species
      If (config_isfixed(species(spc))) Cycle

      !** Get number of atoms and moles in this species
      natms = config_getnatoms(species,spc)
      nmols = config_getnmoles(species,spc)

      !** Move each molecule separately
      Do m = 1,nmols

        !** Perform maximum step checking if desired
        timestep = stepinfo%timestep
        If ((stepinfo%checkstep).Or.(careful)) Then
          damped = .False.
          Do a = 1,natms
            !** Estimate the step size
            stepsize = stepinfo%timestep*mag(species(spc)%coords(a,m)%v) + &
                (stepinfo%timestep**2)*mag(species(spc)%afast(a,m))

            !** Check step size, and get smallest modstep size for molecule
            If (stepsize > stepinfo%maxstep) Then
              damped = .True.
              modstep = stepinfo%timestep*(stepinfo%maxstep/stepsize)
!              Write(*,*) 'Damping integration ',a,m,stepinfo%maxstep/stepsize
              If (modstep < timestep) timestep = modstep
            End If
          End Do
          If (damped) Write(*,*) 'Damping molecule integration ',m,timestep
        End If

        !** Calculate the new positions for the principle coords
        Call velverlet_integratePos(species(spc)%coords(1:natms,m)%rp, &
            species(spc)%coords(1:natms,m)%v, &
            species(spc)%afast(1:natms,m),timestep)

        !** Generate simcell coordinates
        Call config_applypbcs(species,simcell,(/spc,m,0/))

        !** Finish the velocity calculation 
        Call velverlet_halfIntegrateVel(species(spc)%coords(1:natms,m)%v, &
            species(spc)%afast(1:natms,m),stepinfo%timestep)
      End Do
    End Do

  End Subroutine optstep_modvvmove

  !-------------------------------------------------------------------------
  ! Project out the component of the velocities that are not along the 
  ! acceleration vectors
  ! Requires:  species -- configuration
  !            spc -- species number
  !            fast -- acceleration label
  !-------------------------------------------------------------------------
  Subroutine optstep_modvvproject(species,spc,fast)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc
    Logical, Intent(In)                            :: fast

    Integer                 :: a,m,nmols,natms
    Real(kind=RDbl)         :: dotprod
    Type(VecType)           :: unitvec
    Type(VecType), Dimension(:,:), Pointer   :: accels

    !** Get the pointer to the correct accelerations
    Call config_getaccel(species,spc,fast,accels)

    !** Get number of atoms and moles in this species
    natms = config_getnatoms(species,spc)
    nmols = config_getnmoles(species,spc)

    !** Do the projection
    Do m = 1,nmols
      Do a = 1,natms
        unitvec = vector_getunitvec(accels(a,m))
        dotprod = unitvec*species(spc)%coords(a,m)%v
        If (dotprod > 0.0_RDbl) Then
          species(spc)%coords(a,m)%v = dotprod*unitvec
        Else
          species(spc)%coords(a,m)%v = 0.0_RDbl
        End If
      End Do
    End Do

    !** Nullify the pointer
    Nullify(accels)

  End Subroutine optstep_modvvproject

  !------------------------------------------------------------------------
  ! Get the maximum step size
  ! Requires:  stepset -- optimization step pointer set
  !------------------------------------------------------------------------
  Real(kind=RDbl) Function optstep_stepsize(stepset)
    Type(Optimization_Step), Intent(InOut)         :: stepset

    If (Associated(stepset%stpdesc)) Then
      optstep_stepsize = stepset%stpdesc%max_stepsize
    End If

    If (Associated(stepset%modstpdesc)) Then
      optstep_stepsize = stepset%modstpdesc%max_stepsize
    End If

  End Function optstep_stepsize

  !------------------------------------------------------------------------
  ! Set the maximum step size
  ! Requires:  stepset -- optimization step pointer set
  !            stepsize -- new step size
  !------------------------------------------------------------------------
  Subroutine optstep_setstep(stepset,stepsize)
    Type(Optimization_Step), Intent(InOut)         :: stepset
    Real(kind=RDbl), Intent(In)                    :: stepsize

    If (Associated(stepset%stpdesc)) Then
      stepset%stpdesc%max_stepsize = stepsize
    End If

    If (Associated(stepset%modstpdesc)) Then
      stepset%modstpdesc%max_stepsize = stepsize
    End If

    If (Associated(stepset%modvv)) Then
      stepset%modvv%timestep = stepsize
      stepset%modvv%init_velocities = .False.
    End If

  End Subroutine optstep_setstep

  !--------------------------------------------------------------------
  ! Initializes the steepest descent optimization step from numbers
  ! Requires:  stepinfo -- Steepest descent data structure
  !            nums -- array of reals containing parameters for step
  !--------------------------------------------------------------------
  Subroutine optstep_stpdescinit(stepinfo,nums)
    Type(Steepest_Descent), Intent(Out)        :: stepinfo
    Real(kind=RDbl), Dimension(:), Intent(In)  :: nums

    stepinfo%max_stepsize = nums(1)

  End Subroutine optstep_stpdescinit

  !--------------------------------------------------------------------
  ! Initializes the modified steepest descent optimization step 
  ! from numbers
  ! Requires:  stepinfo -- Modified steepest descent data structure
  !            nums -- array of reals containing parameters for step
  !--------------------------------------------------------------------
  Subroutine optstep_modstpdescinit(stepinfo,nums)
    Type(Mod_Steepest_Descent), Intent(Out)        :: stepinfo
    Real(kind=RDbl), Dimension(:), Intent(In)      :: nums

    stepinfo%max_stepsize = nums(1)
    stepinfo%smear_factor = nums(2)
    If ((stepinfo%smear_factor < 0.0_RDbl).Or. &
        (stepinfo%smear_factor > 1.0_RDbl)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Modified steepest descent smear factor must be between 0.0 and 1.0'
      Stop
    End If

  End Subroutine optstep_modstpdescinit

  !-------------------------------------------------------------------------
  ! Initializes the modified velocity verlet optimization step from numbers
  ! Requires:  stepinfo -- modified Velocity Verlet parameters
  !            nums -- array of reals containing parameters for step
  !-------------------------------------------------------------------------
  Subroutine optstep_modvvinit(stepinfo,nums)
    Type(Mod_Velocity_Verlet), Intent(Out)     :: stepinfo
    Real(kind=RDbl), Dimension(:), Intent(In)  :: nums

    stepinfo%init_velocities = .True.
    stepinfo%timestep = nums(1)

    !** If a second number is present, assume it's a maximum step
    If (Size(nums) > 1) Then
      stepinfo%checkstep = .True.
      stepinfo%maxstep = nums(2)
    Else
      stepinfo%maxstep = 0.01  !** Angstroms, used only when careful = True
    End If

    stepinfo%start_temp = 300.0_RDbl

  End Subroutine optstep_modvvinit

  !-------------------------------------------------------------------------
  ! Initialize the velocities for the modified velocity verlet type move.
  ! Basically, copy the forces into the velocity storage, then rescale
  ! the velocities to give some highish kinetic energy temperature.
  ! Requires:  stepinfo -- modified Velocity Verlet parameters
  !            ffresults -- results from forcefield calculation on system
  !            species -- configuration
  !            simcell -- simulation cell information
  !-------------------------------------------------------------------------
  Subroutine optstep_modvvinitvel(stepinfo,species,ffresults)
    Type(Mod_Velocity_Verlet), Intent(InOut)       :: stepinfo
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(Forcefield_Results), Intent(In)           :: ffresults

    Integer                 :: i,j,error
    Integer                 :: spc,mol,atm,nspc,nmols,natms
    Real(kind=RDbl)         :: factor,ke,desired_ke
    Type(VecType), Dimension(:,:), Allocatable   :: force

    nspc = Size(species)
    Do spc = 1,nspc

      !** Skip the FIXED species
      If (config_isfixed(species(spc))) Cycle

      !** Get number of atoms and moles in this species
      natms = config_getnatoms(species,spc)
      nmols = config_getnmoles(species,spc)

      !** Allocate space for the forces
      Allocate(force(natms,nmols), STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'force')    

      !** Extract the forces from the forcefield results
      Call storetop_fastforces(ffresults,spc,nmols,natms,force)

      !** Copy forces into velocity storage
      Do mol = 1,nmols
        species(spc)%coords(1:natms,mol)%v = force(1:natms,mol)
      End Do

      !** Get the scaling factor 
      ke = config_kineticEnergy(species,spc)
      If (Abs(ke) < 1.0e-8) Then
        Write(0,'(2a,i4,a,e14.4)') __FILE__,": ",__LINE__, &
            ' Kinetic Energy is unexpectedly small ',ke
        Stop
      End If
      desired_ke = 3.0_RDbl/2.0_RDbl*natms*nmols*stepinfo%start_temp*kjmole_kb
      factor = Sqrt(desired_ke/ke)

      !** Scale the velocities
      Do mol = 1,nmols
        Do atm = 1,natms
          species(spc)%coords(atm,mol)%v = &
              factor * species(spc)%coords(atm,mol)%v 
        End Do
      End Do

      !** Check it 
      ke = config_kineticEnergy(species,spc)*scaleke
      If ((ke - desired_ke) > 1.0e-2) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Ooops, velocity scaling did not work '
        Stop
      End If

      !** Deallocate space for the forces
      Deallocate(force, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'force')    
    End Do

    Return

    !** Dump velocity magnitudes to screen for checking
    Write(*,'(a)') 'spc,mol,atm, initial velocity (Ang/ps)'
    Do spc = 1,nspc
      !** Skip the FIXED species
      If (config_isfixed(species(spc))) Cycle

      !** Get number of atoms and moles in this species
      natms = config_getnatoms(species,spc)
      nmols = config_getnmoles(species,spc)

      Do mol = 1,nmols
        Do atm = 1,natms
          Write(*,'(3i3,f8.3)') spc,mol,atm,mag(species(spc)%coords(atm,mol)%v)
        End Do
      End Do

    End Do

  End Subroutine optstep_modvvinitvel

  !----------------------------------------------------------------------------
  ! Display the set of optimization steps structure
  ! Requires:  stepset -- optimization step pointer set
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine optstep_display(stepset,indent,unitno)
    Type(Optimization_Step), Intent(In)  :: stepset
    Integer, Intent(In)                  :: indent
    Integer, Optional, Intent(In)        :: unitno

    Integer                           :: unit
    Character(len=indent)             :: blank
    Character(len=strLen)             :: string

    blank = Repeat(' ',indent)
    
    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    If (Associated(stepset%stpdesc)) Then
      Write(unit,'(2a)') blank,'Steepest Descent step initialized'
      string = real2str(stepset%stpdesc%max_stepsize,5)
      Write(unit,'(a,2x,3a)') blank,'maximum step size: ', &
          Trim(string),' Angstroms'
    End If

    If (Associated(stepset%modstpdesc)) Then
      Write(unit,'(2a)') blank,'Modified Steepest Descent step initialized'
      string = real2str(stepset%modstpdesc%max_stepsize,5)
      Write(unit,'(a,2x,3a)') blank,'maximum step size: ', &
          Trim(string),' Angstroms'
      Write(unit,'(a,2x,a,f5.2,a)') blank,'molecular smear factor: ', &
          stepset%modstpdesc%smear_factor,' of atomic forces'
    End If

    If (Associated(stepset%modvv)) Then
      Write(unit,'(2a)') blank,'Modified Velocity Verlet step initialized'
      string = real2str(stepset%modvv%timestep,5)
      Write(unit,'(a,2x,3a)') blank,'timestep: ',Trim(string),' picoseconds'
    End If

  End Subroutine optstep_display

  !----------------------------------------------------------------------------
  ! Clean the structure
  ! Requires:  stepset -- optimization step pointer set
  !----------------------------------------------------------------------------
  Subroutine optstep_clean(stepset)
    Type(Optimization_Step), Intent(InOut)         :: stepset

    Integer         :: error

    If (Associated(stepset%stpdesc)) Then
      Deallocate(stepset%stpdesc, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'stpdesc')
    End If

    If (Associated(stepset%modstpdesc)) Then
      Deallocate(stepset%modstpdesc, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'modstpdesc')
    End If

    If (Associated(stepset%modvv)) Then
      Deallocate(stepset%modvv, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'modvv')
    End If

  End Subroutine optstep_clean

End Module optstep

