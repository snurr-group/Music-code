!------------------------------------------------------------------------------
! Gear 6th order integrator
!
! Needed Improvements:
! 1) a decent header
! 2) requires statements and more comments
! 3) reaches directly into imodel -> HACK (did this to get it running)
! 4) related to this, there are direct calls to intramolecular, use interact
!------------------------------------------------------------------------------

Module gear6

  Use defaults, Only: RDbl, strLen, kcalmole_kb, scalef, scalepe, kjmole_kb, &
      zero, scaleke, one, caltoj, MAX_SORBS
  Use simcell, Only: SimCell_Params
  Use config, Only: config_getconfig, AtMolCoords, config_getnatoms, &
      config_getnmoles, config_isfixed
  Use molecules, Only: molecules_getdof, molecules_dofIsOK, &
      molecules_getnatoms, molecules_AtomMass, &
      molecules_getatype, molecules_getnsorbs, molecules_getmass
  Use nemd, Only: NEMD_IS_ON, NEMDParams, nemd_updateNrg, nemd_getpe, &
      nemd_updateflux
  Use nosehoover, Only: nosehoover_getcoord, NHSorbates
  Use storestats, Only: storestats_updatenrg
  Use interact, Only: Interaction_Model, interact_hasint
  Use subinteract, Only: Subset_Interactions, subinteract_int 
  Use thermostats, Only: thermostats_init, thermostats_rescale, &
      thermostats_display, thermostats_getke, thermostats_getpe, &
      thermostats_T, ThermostatInfo, thermostats_isNH, &
      thermostats_updateenergy, thermostats_sampleCF
  Use intramolecular, Only: intramolecular_getpenalty, &   !HACK
      intramolecular_mconint, intramolecular_postadjust
  Use utils, Only: split, toint, toupper, tolower, allocErrDisplay, &
      checkandstop, cleanstring, deallocErrDisplay
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_display, vector_zerovec, mag
  Use fixedsorbate, Only: fixedsorbate_check

  Implicit None
  Save

  Private
  Public :: GearParams, GearHOT, GearSclrHOT, gear6_init, gear6_display, &
      gear6_integrate, gear6_initIntegrate, gear6_extraenergy, gear6_sampleCF

  Interface gear6_correct
    Module procedure gear6_vecCorrect
    Module procedure gear6_sclrCorrect
  End Interface

  Interface gear6_predict
    Module procedure gear6_vecPredict
    Module procedure gear6_sclrPredict
  End Interface

  Type GearHOT 
    Type(VecType), Dimension(:,:,:), Pointer :: hot
  End Type GearHOT

  Type GearSclrHOT
    Real(Kind=RDbl), Dimension(:), Pointer :: hot
  End Type GearSclrHOT

  ! initc, initcorr1, and initcorr2 are used by initintegrate routine
  Type GearParams
    Integer :: penaltyStep
    Type(ThermostatInfo), Pointer :: thermostat
    Real(kind=RDbl), Dimension(0:5) :: c
    Real(kind=RDbl), Dimension(0:5) :: corr1
    Real(kind=RDbl), Dimension(0:5) :: corr2
    Real(kind=RDbl), Dimension(0:5) :: initc, initcorr1, initcorr2
    Type(GearHOT), Dimension(:), Pointer :: drvs
    Type(GearSclrHOT), Dimension(:), Pointer :: nhDrvs
    Type(GearSclrHOT), Dimension(:,:), Pointer :: nhSpcDrvs
    Integer :: g   ! degrees of freedom in the system
    Integer, Dimension(:), Pointer :: nmoles_list
    Real(kind=RDbl), Dimension(:), Pointer :: spc_g 
    ! spc_g = degrees of freedom for each spc
  End Type GearParams

  Integer, Parameter :: MAX_CORRECTION = 2
  Real(kind=RDbl), Parameter :: NH_TOLERANCE = 1.0e-10_RDBL 

  Integer, Parameter :: GEAR6_INIT_SCALING_FACTOR = 1000

Contains

  !----------------------------------------------------------------------------
  ! Initializes the gear algorithm
  !----------------------------------------------------------------------------
  Subroutine gear6_init(g6,species,dt,ensemble,unitno,T)
    Type(GearParams), Intent(InOut) :: g6
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In) :: unitno
    Real(kind=RDbl), Intent(In) :: dt
    Character(*), Intent(In) :: ensemble
    Real(kind=RDbl), Intent(In), Optional :: T
    Integer :: error,nspc,spc,nfields, natoms,nmoles, dof
    Character(len=strLen) :: text
    Character(len=strLen), Dimension(strLen) :: params

    !** First, read any required information from the file
    Read(unitno,'(a)') text
    nfields = split(cleanstring(text),params)
    g6%penaltyStep = toint(params(1))

    nspc = molecules_getnsorbs()

    !** If it is NVT, read the thermostat info

    !** Make sure default is "no thermostat"
    Nullify(g6%thermostat)

    If (toupper(ensemble) == "NVT") Then
      Allocate(g6%thermostat, stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
            "Could not allocate memory for 'g6%thermostat'"
        Stop
      End If
      !** Call thermostats to read in the information
      Call thermostats_init(g6%thermostat,T,unitno)

      !** If the thermostat is Nose-Hoover, we need to store the 
      !** higher order terms for the 2 extended coordinates
      If (thermostats_isNH(g6%thermostat)) Then
        If (g6%thermostat%nh%species_wise) Then
          If (.Not.NEMD_IS_ON) THEN
            Write (*,*) "species wise thermo stat is currently only for NEMD"
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
          else
          endif

          Allocate(g6%nhSpcDrvs(nspc,2),Stat=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
          Do spc=1,nspc
            !** HOTs for s coordinate and thermal friction
            Allocate(g6%nhSpcDrvs(spc,1)%hot(5),g6%nhSpcDrvs(spc,2)%hot(5), &
                Stat=error)
            If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
            g6%nhSpcDrvs(spc,1)%hot(1:5) = zero
            g6%nhSpcDrvs(spc,2)%hot(1:5) = zero
          End Do
        Else
          Allocate(g6%nhDrvs(2),Stat=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
          !** HOTs for s coordinate and thermal friction
          Allocate(g6%nhDrvs(1)%hot(5),g6%nhDrvs(2)%hot(5), Stat=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
          g6%nhDrvs(1)%hot(1:5) = zero
          g6%nhDrvs(2)%hot(1:5) = zero
        Endif

      End If
    Endif


    !** Calculate the predictor and corrector constants
    Call gear6_calcTerms(g6,dt)

    !** Calculate the predictor and corrector constants for 
    !** initialization integrator
    Call gear6_calcInitTerms(g6,dt)



    !** Array which holds the dimensions of HOT etc. The values here 
    !** should match With , the nmoles in config. 
    Allocate(g6%nmoles_list(nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    g6%nmoles_list=0

    Allocate(g6%drvs(nspc),Stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Allocate memory for the higher order terms and zero them
    Do spc = 1, nspc
      !** Skip any fixed molecules
      If (config_isfixed(species(spc))) Cycle
      natoms = config_getnatoms(species,spc)
      nmoles = config_getnmoles(species,spc)
      g6%nmoles_list(spc)=nmoles
      Allocate(g6%drvs(spc)%hot(5,natoms,nmoles),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      g6%drvs(spc)%hot = VecType(0.0_RDbl)
    End Do

    Allocate(g6%spc_g(nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Calculate g = (total dof in system)
    g6%g = 0
    g6%spc_g = zero
    Do spc = 1, nspc
      !** Skip any fixed molecules
      If (config_isfixed(species(spc))) Cycle
      dof= molecules_getdof(spc)
      ! make sure dof is not junk
      If (molecules_dofIsOK(spc)) Then
        nmoles=config_getnmoles(species,spc)
        g6%g = g6%g + nmoles*dof
        g6%spc_g(spc) = Real(nmoles*dof,kind=RDbl)
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    End Do

  End Subroutine gear6_init

  !----------------------------------------------------------------------------
  ! Calculates the value of the constants and magic correctors given
  ! the time step.
  !----------------------------------------------------------------------------
  Subroutine gear6_calcTerms(g6,dt)
    Type(GearParams), Intent(InOut) :: g6
    Real(Kind=RDbl), Intent(In) :: dt
    Real(kind=RDbl), Dimension(0:5) :: magic1, magic2
    Integer :: i

    !** Assign values to the constants
    g6%c(0) = 1.0_RDbl
    Do i = 1,5
      g6%c(i) = g6%c(i-1)*dt/Real(i)
    End Do

    !** Calculate the correction factor using the magic terms
    magic1(0) = 95.0_RDbl/288.0_RDbl
    magic1(1) = 1.0_RDbl
    magic1(2) = 25.0_RDbl/24.0_RDbl
    magic1(3) = 35.0_RDbl/72.0_RDbl
    magic1(4) = 5.0_RDbl/48.0_RDbl
    magic1(5) = 1.0_RDbl/120.0_RDbl
    magic2(0) = 3.0_RDbl/20.0_RDbl
    magic2(1) = 251.0_RDbl/360.0_RDbl
    magic2(2) = 1.0_RDbl
    magic2(3) = 11.0_RDbl/18.0_RDbl
    magic2(4) = 1.0_RDbl/6.0_RDbl
    magic2(5) = 1.0_RDbl/60.0_RDbl
    Do i = 0,5
      g6%corr1(i) = magic1(i)*g6%c(1)/g6%c(i)
      g6%corr2(i) = magic2(i)*g6%c(2)/g6%c(i)
    End Do

  End Subroutine gear6_calcTerms

  !----------------------------------------------------------------------------
  ! Calculates the value of the constants and magic correctors given
  ! the time step. these will be used only during _initintegrate
  !----------------------------------------------------------------------------
  Subroutine gear6_calcInitTerms(g6,dt_orig)
    Type(GearParams), Intent(InOut) :: g6
    Real(Kind=RDbl), Intent(In) :: dt_orig

    Real(kind=RDbl), Dimension(0:5) :: magic1, magic2
    Real(kind=RDbl) :: dt 
    Integer :: i

    dt = dt_orig / GEAR6_INIT_SCALING_FACTOR

    !** Assign values to the constants
    g6%initc(0) = 1.0_RDbl
    Do i = 1,5
      g6%initc(i) = g6%initc(i-1)*dt/Real(i)
    End Do

    !** Calculate the correction factor using the magic terms
    magic1(0) = 95.0_RDbl/288.0_RDbl
    magic1(1) = 1.0_RDbl
    magic1(2) = 25.0_RDbl/24.0_RDbl
    magic1(3) = 35.0_RDbl/72.0_RDbl
    magic1(4) = 5.0_RDbl/48.0_RDbl
    magic1(5) = 1.0_RDbl/120.0_RDbl
    magic2(0) = 3.0_RDbl/20.0_RDbl
    magic2(1) = 251.0_RDbl/360.0_RDbl
    magic2(2) = 1.0_RDbl
    magic2(3) = 11.0_RDbl/18.0_RDbl
    magic2(4) = 1.0_RDbl/6.0_RDbl
    magic2(5) = 1.0_RDbl/60.0_RDbl

    Do i = 0,5
      g6%initcorr1(i) = magic1(i)*g6%initc(1)/g6%initc(i)
      g6%initcorr2(i) = magic2(i)*g6%initc(2)/g6%initc(i)
    End Do

  End Subroutine gear6_calcInitTerms


  !----------------------------------------------------------------------------
  ! Writes a sample of the require control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine gear6_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a,t30,a)') 'Integer', &
        '# Steps between calls to penalty function'
    Write(unitno,'(2a)') '# If NVT was specified, remove this line and ',&
        'include the following:'
    Call thermostats_sampleCF(unitno)

  End Subroutine gear6_sampleCF


  !----------------------------------------------------------------------------
  ! Checks whether the dimesnion of arrays in gr matches with that in species
  !----------------------------------------------------------------------------
  Subroutine gear6_checkandincr(g6,species)
    Type(GearParams), Intent(InOut)                :: g6
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer :: nspc, spc, oldsize, nmoles, natoms, error, i
    Type(VecType), Dimension(:,:,:), Allocatable :: tempptr

    !** g = (total dof in system)
    g6%g = 0

    nspc=Size(species,1)
    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle

      nmoles=config_getnmoles(species,spc)
      oldsize=g6%nmoles_list(spc)
      g6%g = g6%g + nmoles * molecules_getdof(spc)

      If (nmoles<=oldsize) Cycle

      !** Now we need to increase the size
      natoms = config_getnatoms(species,spc)      
      g6%nmoles_list(spc)=nmoles
      Allocate(tempptr(5,natoms,oldsize), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      ! copy to temperory storage array
      tempptr(1:5,1:natoms,1:oldsize)=g6%drvs(spc)%hot(1:5,1:natoms,1:oldsize)

      Deallocate(g6%drvs(spc)%hot,stat=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)

      Allocate(g6%drvs(spc)%hot(5,natoms,nmoles),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      ! copy values back from temperory storage array
      g6%drvs(spc)%hot(1:5,1:natoms,1:oldsize) = tempptr(1:5,1:natoms,1:oldsize)

      ! ** make new values zero
      g6%drvs(spc)%hot(1:5,1:natoms,(oldsize+1):nmoles) = vector_zerovec()

      Deallocate(tempptr, stat=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)

    End Do

  End Subroutine gear6_checkandincr

  !----------------------------------------------------------------------------
  ! This routine decides which (NVE, VR-NVT, NH-NVT) gear implimentation 
  ! to use
  !            success -- indicates whether integration was successful
  !----------------------------------------------------------------------------
  Subroutine gear6_integrate(subint, g6,simcell,species, &
      ensemble,dt,fast,step,success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(GearParams), Intent(InOut)                :: g6
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(Kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(In)                            :: fast
    Integer, Intent(In)                            :: step
    Logical, Intent(out)                           :: success

    Character(*), Intent(In)                       :: ensemble

    !** Check whether the dimensions of hot match with nmoles in config
    ! also reestimates the d.o.f, maybe this should not be called 
    ! during regular MD, waste of time
    Call gear6_checkandincr(g6, species)

    !** Select ensemble
    Select Case (toupper(ensemble))
    Case ("NVT")

      !** 
      If (NEMD_IS_ON) Then
        If (Associated(g6%thermostat%nh)) Then
          If (g6%thermostat%nh%species_wise) Then
            Call gear6_MultnhNEMD(subint, g6, simcell, species, dt, fast, &
                step,success)
          Else
            Call gear6_nhNEMD(subint, g6, simcell, species, dt, fast, &
                step,success)
          Endif
        Else
          Write(*,'(1x,2a,i4)') "NEMD as of now needs NH thermostat"
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        endif

        Return

      Endif

      !** NVT ensemble. Retrieve the thermostat and call the correct
      !** integrator
      If (Associated(g6%thermostat%vr)) Then
        Call gear6_vrIntegrate(subint,g6,simcell,species,dt,fast,step,success)
      Else If (Associated(g6%thermostat%nh)) Then
        Call gear6_nhIntegrate(subint,g6,simcell,species,dt,fast,step,success)
      End If




    Case ("NVE")
      !** NVE ensemble. Call the NVE integrator.
      Call gear6_nveIntegrate(subint,g6,simcell,species,dt,fast,step,success)

    End Select

  End Subroutine gear6_integrate

  !----------------------------------------------------------------------------
  ! This routine performs a single integration step of the 6th order
  ! Gear algorithm
  ! Requires:  g6 -- Gear integration parameters
  !            simcell -- definition of the simulation cell
  !            species -- the vital coordinates for the entire system
  !            dt -- time step ?
  !            fast -- 
  !            step -- ?
  !            success -- indicates whether integration was successful
  !----------------------------------------------------------------------------
  Subroutine gear6_nveIntegrate(subint, g6, simcell, species, &
      dt, fast, step,success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(GearParams), Intent(InOut)                :: g6
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(Kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(In)                            :: fast
    Integer, Intent(In)                            :: step
    Logical, Intent(out)                           :: success

    Integer                         :: nspc, natoms, nmoles
    Integer                         :: i, m, a, c
    Integer, Save                   :: steps          
    Logical                         :: mapflag, nrg_success
    Real(Kind=RDbl)                 :: u, sumobjf, uconst, TKelvin
    Type(AtMolCoords), Pointer      :: single

    !** Increment the steps counter
    steps = steps + 1
    success=.True.

    !** Get the number of species in the system
    nspc = molecules_getnsorbs()

    !** Loop through each species and do the prediction step
    Do i = 1,nspc

      natoms = molecules_getnatoms(i)
      nmoles = config_getnmoles(species,i)
      Call config_getconfig(species,i,single)

      !** Skip this species if it is fixed in place
      If (config_isfixed(single)) Cycle
      If ( fixedsorbate_check(i)) Cycle

      !** Predict the positions and velocities
      Do m = 1, nmoles
        Do a = 1, natoms
          Call gear6_predict(single%coords(a,m)%rp,single%coords(a,m)%v, &
              g6%drvs(i)%hot(1:5,a,m),g6%c,.True.)
          Call gear6_predict(single%coords(a,m)%r,single%coords(a,m)%v, &
              g6%drvs(i)%hot(1:5,a,m),g6%c)
          single%afast(a,m) = 0.0_RDbl
          single%aslow(a,m) = 0.0_RDbl
        End Do
      End Do

      !** Do any necessary post-integration adjustments
      !HACK! reaches into imodel
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i), &
          simcell,single)

    End Do

    !** Calculate the intermediate forces
    TKelvin = 200.0_RDbl ! dummy for HMD
    nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
    If (.Not.nrg_success) Then
      success=nrg_success
      Return
    Endif

    !** We check the penalty function after the correction and loop
    !** Correction-Constraints Loop
    Do c = 1, MAX_CORRECTION
      Do i = 1, nspc

        natoms = config_getnatoms(species,i)
        nmoles = config_getnmoles(species,i)
        Call config_getconfig(species,i,single)

        !** Skip this species if it is fixed in place
        If (config_isfixed(single)) Cycle
        If ( fixedsorbate_check(i)) Cycle

        !** Correct the positions and velocities using the new forces
        Do m = 1, nmoles
          Do a = 1, natoms
            If (fast) Then
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.True.)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2)
            Else
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                  g6%corr2,.True.)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                  g6%corr2)
            End If
          End Do
        End Do

        !** Call the constraints function
        If (.Not. interact_hasint(subint%imodel,i,'constraint')) Cycle

        uconst = 0.0_RDbl
        Do m = 1, nmoles
          If (Mod(steps,g6%penaltyStep) == 0) Then
            !** Call the penalty function
            Call intramolecular_getpenalty(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp,single%coords(1:natoms,m)%v,sumobjf)
          End If

          If (fast) Then
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp,single%coords(1:natoms,m)%v,u, &
                single%afast(1:natoms,m),.False.)
            uconst = u + uconst
          Else
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp,single%coords(1:natoms,m)%v,u, &
                single%aslow(1:natoms,m),.False.)
            uconst = u + uconst
          End If

        End Do

        Call storestats_updatenrg(subint%imodel%spcstats,i,'const',uconst)
        Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i), &
            simcell,single)

      End Do

#ifdef EXTRA_EVAL
      !** Increases the computational cost (says Shaji) and 
      !** doesn't improve energy conservation for Butane/sili/NVE MD
      nrg_success = subinteract_int(subint,species,simcell, &
          fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
      If (.Not.nrg_success) Then
        success=nrg_success
        Return
      Endif
#endif   

    End Do

  End Subroutine gear6_nveIntegrate


  !----------------------------------------------------------------------------
  ! This routine performs 8 (has to be more than 6) single integration 
  ! steps of the 6th order NVE-Gear algorithm with reduced timesteps so 
  ! that the hot's are initialized before the main integration. 
  ! It is always good to Call this init routine. It is tested with
  ! MD-GCMC and regular MD of butane
  ! Requires:  g6 -- Gear integration parameters
  !            simcell -- definition of the simulation cell
  !            species -- the vital coordinates for the entire system
  !            dt -- time step ?
  !            fast -- 
  !            step -- ?
  !            success -- indicates whether integration was successful
  !----------------------------------------------------------------------------
  Subroutine gear6_InitIntegrate(subint, g6, simcell, species, &
      dt, fast,success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(GearParams), Intent(InOut)                :: g6
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(Kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(In)                            :: fast
    Logical, Intent(out)                           :: success

    Integer                         :: nspc, natoms, nmoles, inititers
    Integer                         :: i, m, a, c, i_step, steps          
    Logical                         :: mapflag, nrg_success
    Real(Kind=RDbl)                 :: u, sumobjf, uconst, TKelvin
    Type(AtMolCoords), Pointer      :: single

    success=.True.
    inititers=8
    steps = 0

    !** Check whether the dimensions of hot match with nmoles in config
    ! also reestimates the d.o.f
    Call gear6_checkandincr(g6, species)

    !** Get the number of species in the system
    nspc = molecules_getnsorbs()

    Do i_step=1,inititers

      !** Loop through each species and do the prediction step
      Do i = 1,nspc

        natoms = molecules_getnatoms(i)
        nmoles = config_getnmoles(species,i)
        Call config_getconfig(species,i,single)

        !** Skip this species if it is fixed in place
        If (config_isfixed(single)) Cycle
        If ( fixedsorbate_check(i)) Cycle

        !** Predict the positions and velocities
        Do m = 1, nmoles
          Do a = 1, natoms
            Call gear6_predict(single%coords(a,m)%rp,single%coords(a,m)%v, &
                g6%drvs(i)%hot(1:5,a,m),g6%initc,.True.)
            Call gear6_predict(single%coords(a,m)%r,single%coords(a,m)%v, &
                g6%drvs(i)%hot(1:5,a,m),g6%initc)
            single%afast(a,m) = 0.0_RDbl
            single%aslow(a,m) = 0.0_RDbl
          End Do
        End Do

        !** Do any necessary post-integration adjustments
        !HACK! reaches into imodel
        Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i), &
            simcell,single)

      End Do


      !** Calculate the intermediate forces
      TKelvin = 200.0_RDbl ! dummy for HMD
      nrg_success = subinteract_int(subint,species,simcell, &
          fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
      If (.Not.nrg_success) Then
        success=nrg_success
        Return
      Endif

      !** We check the penalty function after the correction and loop
      !** Correction-Constraints Loop
      Do c = 1, MAX_CORRECTION
        Do i = 1, nspc

          natoms = config_getnatoms(species,i)
          nmoles = config_getnmoles(species,i)
          Call config_getconfig(species,i,single)

          !** Skip this species if it is fixed in place
          If (config_isfixed(single)) Cycle
          If ( fixedsorbate_check(i)) Cycle

          !** Correct the positions and velocities using the new forces
          Do m = 1, nmoles
            Do a = 1, natoms
              If (fast) Then
                Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                    g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                    g6%initcorr2,.True.)
                Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                    g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                    g6%initcorr2)
              Else
                Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                    g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                    g6%initcorr2,.True.)
                Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                    g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                    g6%initcorr2)
              End If
            End Do
          End Do

          !** Call the constraints function
          If (.Not. interact_hasint(subint%imodel,i,'constraint')) Cycle

          uconst = 0.0_RDbl
          Do m = 1, nmoles
            If (Mod(steps,g6%penaltyStep) == 0) Then
              !** Call the penalty function
              Call intramolecular_getpenalty(subint%imodel%ff(1)%iparams(i), &
                  single%coords(1:natoms,m)%rp,single%coords(1:natoms,m)%v,sumobjf)
            End If

            If (fast) Then
              Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                  single%coords(1:natoms,m)%rp,single%coords(1:natoms,m)%v,u, &
                  single%afast(1:natoms,m),.False.)
              uconst = u + uconst
            Else
              Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                  single%coords(1:natoms,m)%rp,single%coords(1:natoms,m)%v,u, &
                  single%aslow(1:natoms,m),.False.)
              uconst = u + uconst
            End If

          End Do

          Call storestats_updatenrg(subint%imodel%spcstats,i,'const',uconst)
          Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i), &
              simcell,single)

        End Do

#ifdef EXTRA_EVAL
        !** Increases the computational cost (says Shaji) and 
        !** doesn't improve energy conservation for Butane/sili/NVE MD
        nrg_success = subinteract_int(subint,species,simcell, &
            fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
        If (.Not.nrg_success) Then
          success=nrg_success
          Return
        Endif
#endif   

      End Do


    End Do

  End Subroutine gear6_InitIntegrate


  !----------------------------------------------------------------------------
  ! This routine performs a single integration step of the 6th order
  ! Gear algorithm using velocity rescaling to control temperature.
  !            success -- indicates whether integration was successful
  !----------------------------------------------------------------------------
  Subroutine gear6_vrIntegrate (subint, g6, simcell, species, &
      dt, fast, step, success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(GearParams), Intent(InOut)                :: g6
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(Kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(In)                            :: fast
    Integer, Intent(In)                            :: step
    Logical, Intent(out)                           :: success

    Integer                    :: nspc, i, m, a, natoms, nmoles, c
    Integer, Save              :: steps
    Logical                    :: mapflag,nrg_success
    Real(Kind=RDbl)            :: u, sumobjf, uconst, TKelvin
    Type(AtMolCoords), Pointer :: single

    steps = steps + 1
    success=.True.
    nspc = molecules_getnsorbs()

    Do i = 1,nspc
      natoms = molecules_getnatoms(i)
      nmoles = config_getnmoles(species,i)
      Call config_getconfig(species,i,single)

      !** Skip this species if it is fixed in place
      If (config_isfixed(single)) Cycle
      If ( fixedsorbate_check(i)) Cycle

      !** Perform the velocity rescaling first
      Call thermostats_rescale(subint%imodel,g6%thermostat,i, &
          single%coords(1:natoms,1:nmoles)%v,dt)

      !** Predict the positions and velocities
      Do m = 1, nmoles
        Do a = 1, natoms
          Call gear6_predict(single%coords(a,m)%rp,single%coords(a,m)%v, &
              g6%drvs(i)%hot(1:5,a,m),g6%c,.True.)
          Call gear6_predict(single%coords(a,m)%r,single%coords(a,m)%v, &
              g6%drvs(i)%hot(1:5,a,m),g6%c)
          single%afast(a,m) = 0.0_RDbl
          single%aslow(a,m) = 0.0_RDbl
        End Do
      End Do

      !** Do any necessary post-integration adjustments
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i), &
          simcell,single)

    End Do

    !** Calculate the intermediate forces
    TKelvin = 300.0_RDbl  !** HACK
    nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
    If (.Not.nrg_success) Then
      success=nrg_success
      Return
    Endif

    Do i = 1,nspc

      natoms = config_getnatoms(species,i)
      nmoles = config_getnmoles(species,i)
      Call config_getconfig(species,i,single)

      !** Skip this species if it is fixed in place
      If (config_isfixed(single)) Cycle
      If ( fixedsorbate_check(i)) Cycle

      !** Correction-Constraints Loop
      !** We check the penalty function after the correction and loop

      Do c = 1, MAX_CORRECTION

        !** Correct the positions and velocities using the new forces
        Do m = 1, nmoles
          Do a = 1, natoms
            If (fast) Then
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.True.)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2)
            Else
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                  g6%corr2,.True.)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                  g6%corr2)
            End If
          End Do
        End Do

        !** Call the constraints function
        If (.Not.interact_hasint(subint%imodel,i,'constraint')) Cycle
        uconst = 0.0_RDbl
        Do m = 1, nmoles
          If (Mod(steps,g6%penaltyStep) == 0) Then
            !** Call the penalty function
            Call intramolecular_getpenalty(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,sumobjf)
          End If
          If (fast) Then
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,u,single%afast(1:natoms,m),.False.)
            uconst = u + uconst
          Else
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,u,single%aslow(1:natoms,m),.False.)
            uconst = u + uconst
          End If
        End Do
        Call storestats_updatenrg(subint%imodel%spcstats,i,'const',uconst)

      End Do

    End Do

    !** Do any necessary post-integration adjustments
    Do i = 1,nspc
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i), &
          simcell,species(i))
    End Do

  End Subroutine gear6_vrIntegrate

  !----------------------------------------------------------------------------
  ! This routine performs a single integration step of the 6th order
  ! Gear algorithm using the Nose-Hoover thermostat. Here both the extended
  ! variables "s" and "Eta" are integrated using first-order predictor
  ! corrector methods
  !            success -- indicates whether integration was successful
  !----------------------------------------------------------------------------
  Subroutine gear6_nhIntegrate(subint, g6, simcell, species, dt, &
      fast, step, success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(GearParams), Intent(InOut)                :: g6
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(Kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(In)                            :: fast
    Integer, Intent(In)                            :: step
    Logical, Intent(out)                           :: success

    Integer          :: nspc, i, m, a, natoms, nmoles, c, nhstep 
    Integer, Save    :: steps
    Logical          :: mapflag,nrg_success
    Real(Kind=RDbl)  :: u, sumobjf, uconst, mass
    Real(Kind=RDbl)  :: mvsq   ! mass*vel**2 - for the Nose-Hoover xi term
    Real(Kind=RDbl)  :: cortf, cors ! correction term for Nose-Hoover coords
    Real(kind=RDbl)  :: ke, pe   ! Nose-Hoover hamiltonian
    Real(Kind=RDbl)  ::  TKelvin,beta
    Type(VecType)    :: correct
    Type(AtMolCoords), Pointer :: single
    Type(NHSorbates), Pointer :: nhsorb

    nhstep = 0
    mvsq   = 0.0_RDbl

    !** Beta has units of ps^2/Ang^2
    beta = 1.0_RDbl / (thermostats_T(g6%thermostat)*kcalmole_kb*scalef)

    steps = steps + 1
    success=.True.
    nspc = molecules_getnsorbs()

    Do i = 1,nspc
      natoms = molecules_getnatoms(i)
      nmoles = config_getnmoles(species,i)
      Call config_getconfig(species,i,single)

      !** Skip this species if it is fixed in place
      If (config_isfixed(single)) Cycle
      If ( fixedsorbate_check(i)) Cycle

      !** Predict the positions and velocities
      Do m = 1, nmoles
        Do a = 1, natoms
          Call gear6_predict(single%coords(a,m)%rp,single%coords(a,m)%v, &
              g6%drvs(i)%hot(1:5,a,m),g6%c,.True.)
          Call gear6_predict(single%coords(a,m)%r,single%coords(a,m)%v, &
              g6%drvs(i)%hot(1:5,a,m),g6%c)
          single%afast(a,m) = 0.0_RDbl
          single%aslow(a,m) = 0.0_RDbl

          !** Sum the momentum over mass
          !** This probably should be a routine in Nose-Hoover
          !** but since we have a loop here, we'll use it
          mass = molecules_AtomMass(i,a)
          mvsq = mvsq + single%coords(a,m)%v*single%coords(a,m)%v*mass

        End Do
      End Do

      !** Do any necessary post-integration adjustments
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i),&
          simcell,single)

    End Do

    !** Predict the values for the extended Nose-Hoover coordinates
    Call nosehoover_getcoord(g6%thermostat%nh,nhsorb)
    Call gear6_predict(nhsorb%s%r,nhsorb%s%v, &
        g6%nhDrvs(1)%hot(1:5),g6%c,.False.)
    Call gear6_predict(nhsorb%tf%r,nhsorb%tf%v, &
        g6%nhDrvs(2)%hot(1:5),g6%c,.False.)

    !** Calculate the intermediate forces
    !!    betaKcal = 1.0_RDbl / (thermostats_T(g6%thermostat)*kcalmole_kb)
    TKelvin = thermostats_T(g6%thermostat)
    nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
    If (.Not.nrg_success) Then
      success=nrg_success
      Return
    Endif

    !** Correction-Constraints Loop, MAX 2 corrections
    !** We check the penalty function after the correction and loop
    Do c = 1, MAX_CORRECTION
      Do i = 1, nspc

        natoms = config_getnatoms(species,i)
        nmoles = config_getnmoles(species,i)
        Call config_getconfig(species,i,single)

        !** Skip this species if it is fixed in place
        If (config_isfixed(single)) Cycle
        If ( fixedsorbate_check(i)) Cycle

        !** Correct the positions and velocities using the new forces
        Do m = 1, nmoles
          Do a = 1, natoms
            If (fast) Then
              correct = (single%afast(a,m) - single%coords(a,m)%v*nhsorb%tf%r- &
                  g6%drvs(i)%hot(2,a,m))
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.True.,correct)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.False.,correct)
            Else
              correct = (single%aslow(a,m)-single%coords(a,m)%v*nhsorb%tf%r- &
                  g6%drvs(i)%hot(2,a,m))
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                  g6%corr2,.True.,correct)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                  g6%corr2,.False.,correct)
            End If
          End Do
        End Do
      End Do

      !** correction factor for the thermal friction
      !** d(tfriction)/dt = (sum(momentum**2/mass) - dof*T*kb)/Qnose - d(tf)/dt
      cortf = ((mvsq - Real(g6%g,kind=RDbl)/beta)/(g6%thermostat%nh%Q*scalef) &
          - nhsorb%tf%v)

      !** correction factor for the s coordinate
      !** correct = xi*s - ds/dt
      cors = nhsorb%tf%r*nhsorb%s%r - nhsorb%s%v

      Call gear6_correct(nhsorb%s%r, nhsorb%s%v,g6%nhDrvs(1)%hot(1:5), &
          g6%nhDrvs(1)%hot(2),g6%corr1,.False.,cors)

      Call gear6_correct(nhsorb%tf%r,nhsorb%tf%v,g6%nhDrvs(2)%hot(1:5), &
          g6%nhDrvs(2)%hot(2),g6%corr1,.False.,cortf)

      !!      nhsorb%tf%v = (mvsq - Real(g6%g,kind=RDbl)/beta)/ &
      !!          (g6%thermostat%nh%Q*scalef)

      !** Check to make sure that we are within tolerance limits for
      !** the s coordinate and thermal friction
      If (Abs(nhsorb%tf%r*nhsorb%s%r - nhsorb%s%v) > NH_TOLERANCE) Then
        nhstep = nhstep + 1
        Cycle
      End If

      Do i = 1, nspc

        !** Skip this species if it is fixed in place
        If (config_isfixed(species(i))) Cycle
        If ( fixedsorbate_check(i)) Cycle

        !** Call the constraints function
        If (.Not.interact_hasint(subint%imodel,i,'constraint')) Cycle
        uconst = 0.0_RDbl
        nmoles = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)
        Call config_getconfig(species,i,single)
        Do m = 1, nmoles
          If (Mod(steps,g6%penaltyStep) == 0) Then
            !** Call the penalty function
            Call intramolecular_getpenalty(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,sumobjf)
          End If
          If (fast) Then
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,u,single%afast(1:natoms,m),.False.)
            uconst = u + uconst
          Else
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,u,single%aslow(1:natoms,m),.False.)
            uconst = u + uconst
          End If
        End Do
        Call storestats_updatenrg(subint%imodel%spcstats,i,'const',uconst)

      End Do
    End Do

    !** Do any necessary post-integration adjustments
    Do i = 1,nspc
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i), &
          simcell,species(i))
    End Do

    !** Update energies, In this routine we are propagating "s" not log(s)
    ke = 0.5_RDbl*g6%thermostat%nh%Q*nhsorb%tf%r**2 
    Call thermostats_updateEnergy(g6%thermostat,'ke',ke*scalepe)
    pe = Real(g6%g,kind=RDbl)*log(nhsorb%s%r)*thermostats_T(g6%thermostat) &
        *kjmole_kb
    Call thermostats_updateEnergy(g6%thermostat,'pe',pe)

  End Subroutine gear6_nhIntegrate



  !----------------------------------------------------------------------------
  ! This routine performs a single integration step of the 6th order
  ! Gear algorithm using the Nose-Hoover thermostat. Here both the extended
  ! variables "s" and "Eta" are integrated using first-order predictor
  ! corrector methods
  ! - This particular routine is for NEMD where an external field is 
  !   also applied
  ! - no more slow is allowed, only fast interactions are allowed
  !            success -- indicates whether integration was successful
  !----------------------------------------------------------------------------
  Subroutine gear6_nhNEMD(subint, g6, simcell, species, dt, &
      fast, step, success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(GearParams), Intent(InOut)                :: g6
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(Kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(In)                            :: fast
    Integer, Intent(In)                            :: step
    Logical, Intent(out)                           :: success

    Integer          :: nspc, spc, m, a, natoms, nmoles, c, nhstep 
    Integer, Save    :: steps
    Logical          :: mapflag,nrg_success, add_nemdforce
    Real(Kind=RDbl)  :: u, sumobjf, uconst, mass, ratio
    Real(Kind=RDbl)  :: mvsq   ! mass*vel**2 - for the Nose-Hoover xi term
    Real(Kind=RDbl)  :: cortf, cors ! correction term for Nose-Hoover coords
    Real(kind=RDbl)  :: ke, pe   ! Nose-Hoover hamiltonian
    Real(Kind=RDbl)  ::  TKelvin,beta, mult, max_v, first_atom_vel
    Type(VecType)    :: correct
    Type(VecType)    :: nemd_force
    Type(AtMolCoords), Pointer :: single
    Type(NHSorbates), Pointer :: nhsorb
    ! variable for checking magnitude ratios
    Logical, Save :: check_mags=.True.
    Real(Kind=RDbl),  Save  ::  mag_ratio=0.00
    Integer, Save :: mag_count=0

    nhstep = 0
    mvsq   = 0.0_RDbl


    !** Beta has units of ps^2/Ang^2
    beta = 1.0_RDbl / (thermostats_T(g6%thermostat)*kcalmole_kb*scalef)

    steps = steps + 1
    success=.True.
    nspc = molecules_getnsorbs()
    max_v=zero
    Do spc = 1,nspc
      natoms = molecules_getnatoms(spc)
      nmoles = config_getnmoles(species,spc)
      Call config_getconfig(species,spc,single)

      !** Skip this species if it is fixed in place
      If (config_isfixed(single)) Cycle

      !** Predict the positions and velocities
      Do m = 1, nmoles


        Do a = 1, natoms
          Call gear6_predict(single%coords(a,m)%rp,single%coords(a,m)%v, &
              g6%drvs(spc)%hot(1:5,a,m),g6%c,.True.)
          Call gear6_predict(single%coords(a,m)%r,single%coords(a,m)%v, &
              g6%drvs(spc)%hot(1:5,a,m),g6%c)
          single%afast(a,m) = 0.0_RDbl
          mass = molecules_AtomMass(spc,a)
          mvsq = mvsq + single%coords(a,m)%v*single%coords(a,m)%v*mass

        End Do

        ! update current maximum velocity based on first atom
        first_atom_vel=mag(single%coords(1,m)%v)
        If (first_atom_vel>max_v) max_v=first_atom_vel


      End Do

      !** Do any necessary post-integration adjustments
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(spc),&
          simcell,single)

    End Do

    !** Predict the values for the extended Nose-Hoover coordinates
    Call nosehoover_getcoord(g6%thermostat%nh,nhsorb)
    Call gear6_predict(nhsorb%s%r,nhsorb%s%v, &
        g6%nhDrvs(1)%hot(1:5),g6%c,.False.)
    Call gear6_predict(nhsorb%tf%r,nhsorb%tf%v, &
        g6%nhDrvs(2)%hot(1:5),g6%c,.False.)


    !** Calculate the intermediate forces
    !!    betaKcal = 1.0_RDbl / (thermostats_T(g6%thermostat)*kcalmole_kb)
    TKelvin = thermostats_T(g6%thermostat)
    nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
    If (.Not.nrg_success) Then
      success=nrg_success
      Return
    Endif

    ! potential nrg of color field, not really needed, but helpful in checking 
    ! hamiltonian conservation
    Call nemd_updateNrg(species)

    ! update flux values for reporting on the screen
    Call nemd_updateFlux(species,simcell,max_v*dt)

    !** Correction-Constraints Loop, MAX 2 corrections
    !** We check the penalty function after the correction and loop
    Do c = 1, MAX_CORRECTION
      Do spc = 1, nspc
        add_nemdforce=.False.
        If (NEMDParams%sorb==spc) add_nemdforce=.True.
        natoms = config_getnatoms(species,spc)
        nmoles = config_getnmoles(species,spc)
        Call config_getconfig(species,spc,single)

        !** Skip this species if it is fixed in place
        If (config_isfixed(single)) Cycle

        !** Correct the positions and velocities using the new forces
        Do m = 1, nmoles
          Do a = 1, natoms
            If (fast) Then

              If (c==1) Then
                !** Update accels with NEMD force also, so we add nemd accels 
                !** to regular accels, so this is needed only once when c==1
                ! units : mass = amu
                !         scalef = (amu. ang^2/ ps^2) / (Kcal/mol) 
                !         mult = (ang^2/ps^2)/(Kcal/mol)
                ! nemd_force*mult = ang/ps^2
                mass=molecules_AtomMass(spc,a)  
                mult=(scalef*one/mass)          ! 

                ! ** keep stats about NEMD force / regular force
                ! ** only done for first 50,000 iters
                If (check_mags.and.add_nemdforce) Then
                  !                  ratio= 100* mag((NEMDParams%force(a)*mult)) &
                  !                      / mag( single%afast(a,m))
                  !                  mag_ratio=(mag_ratio*mag_count+ratio)/(mag_count+1)
                  mag_count=mag_count+1
                  If (Mod(mag_count,10000)==0) Then
                    Write(*,*) "-------------- Ratio of NEMD/Reular accels --"
                    Write(*,'(a,f15.8)') " In Percentage", mag_ratio
                    If (mag_count==500000) check_mags=.False.
                  Endif
                Endif


                ! NEMDParams%force(a) = force from color field in Kcal/mol/ang.
                If (add_nemdforce) &
                    single%afast(a,m) = single%afast(a,m) + &
                    (NEMDParams%force(a)*mult)
              Endif

              correct = (single%afast(a,m) - &
                  single%coords(a,m)%v*nhsorb%tf%r- g6%drvs(spc)%hot(2,a,m))
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(spc)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.True.,correct)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(spc)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.False.,correct)
            Else
              Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
              Stop
            End If
          End Do
        End Do
      End Do


      !** correction factor for the thermal friction
      !** d(tfriction)/dt = (sum(momentum**2/mass) - dof*T*kb)/Qnose - d(tf)/dt
      cortf = ((mvsq - Real(g6%g,kind=RDbl)/beta)/(g6%thermostat%nh%Q*scalef) &
          - nhsorb%tf%v)

      !** correction factor for the s coordinate
      !** correct = xi*s - ds/dt
      cors = nhsorb%tf%r*nhsorb%s%r - nhsorb%s%v

      Call gear6_correct(nhsorb%s%r, nhsorb%s%v,g6%nhDrvs(1)%hot(1:5), &
          g6%nhDrvs(1)%hot(2),g6%corr1,.False.,cors)

      Call gear6_correct(nhsorb%tf%r,nhsorb%tf%v,g6%nhDrvs(2)%hot(1:5), &
          g6%nhDrvs(2)%hot(2),g6%corr1,.False.,cortf)

      !!      nhsorb%tf%v = (mvsq - Real(g6%g,kind=RDbl)/beta)/ &
      !!          (g6%thermostat%nh%Q*scalef)

      !** Check to make sure that we are within tolerance limits for
      !** the s coordinate and thermal friction
      If (Abs(nhsorb%tf%r*nhsorb%s%r - nhsorb%s%v) > NH_TOLERANCE) Then
        nhstep = nhstep + 1
        Cycle
      End If

      Do spc = 1, nspc

        !** Skip this species if it is fixed in place
        If (config_isfixed(species(spc))) Cycle

        !** Call the constraints function
        If (.Not.interact_hasint(subint%imodel,spc,'constraint')) Cycle
        uconst = 0.0_RDbl
        nmoles = config_getnmoles(species,spc)
        natoms = molecules_getnatoms(spc)
        Call config_getconfig(species,spc,single)
        Do m = 1, nmoles
          If (Mod(steps,g6%penaltyStep) == 0) Then
            !** Call the penalty function
            Call intramolecular_getpenalty(subint%imodel%ff(1)%iparams(spc), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,sumobjf)
          End If
          If (fast) Then
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(spc), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,u,single%afast(1:natoms,m),.False.)
            uconst = u + uconst
          Else
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
          End If
        End Do
        Call storestats_updatenrg(subint%imodel%spcstats,spc,'const',uconst)

      End Do
    End Do

    !** Do any necessary post-integration adjustments
    Do spc = 1,nspc
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(spc), &
          simcell,species(spc))
    End Do

    !** Update energies, In this routine we are propagating "s" not log(s)
    ke = 0.5_RDbl*g6%thermostat%nh%Q*nhsorb%tf%r**2 
    Call thermostats_updateEnergy(g6%thermostat,'ke',ke*scalepe)
    pe = Real(g6%g,kind=RDbl)*log(nhsorb%s%r)*thermostats_T(g6%thermostat) &
        *kjmole_kb
    Call thermostats_updateEnergy(g6%thermostat,'pe',pe)

  End Subroutine gear6_nhNEMD



  !----------------------------------------------------------------------------
  ! 
  ! This routine performs a single integration step of the 6th order
  ! Gear algorithm using 
  ! - one Nose-Hoover thermostat for each species
  ! - "s" is integrated using the second-order predictor corrector methods
  ! - This particular routine is for NEMD where an external field is 
  ! - no more "slow" is allowed, only fast interactions are allowed
  !            success -- indicates whether integration was successful
  ! - works with constant N only
  !----------------------------------------------------------------------------
  Subroutine gear6_MultNhNEMD(subint, g6, simcell, species, dt, &
      fast, step, success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(GearParams), Intent(InOut)                :: g6
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(Kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(In)                            :: fast
    Integer, Intent(In)                            :: step
    Logical, Intent(out)                           :: success

    Integer          :: nspc, spc, m, a, natoms, nmoles, c, nhstep 
    Integer, Save    :: steps
    Logical          :: mapflag,nrg_success, add_nemdforce
    Real(Kind=RDbl)  :: u, sumobjf, uconst, mass, ratio
    Real(Kind=RDbl),Dimension(MAX_SORBS)  :: mvsq   
    ! mass*vel**2 - for the Nose-Hoover xi term
    Real(Kind=RDbl)  :: cortf, acc_nh! correction term for Nose-Hoover coords
    Real(kind=RDbl)  :: ke, pe   ! Nose-Hoover hamiltonian
    Real(Kind=RDbl)  ::  TKelvin,beta, mult, max_v, first_atom_vel
    Type(VecType)    :: correct
    Type(VecType)    :: nemd_force
    Type(AtMolCoords), Pointer :: single
    Type(NHSorbates), Dimension(:),Pointer :: nhsorb
    ! variable for checking magnitude ratios
    Logical, Save :: check_mags=.True.
    Real(Kind=RDbl),  Save  ::  mag_ratio=0.00
    Integer, Save :: mag_count=0

    nhstep = 0
    mvsq   = zero

    !** Beta has units of ps^2/Ang^2
    beta = 1.0_RDbl / (thermostats_T(g6%thermostat)*kcalmole_kb*scalef)

    steps = steps + 1
    success=.True.
    nspc = molecules_getnsorbs()
    max_v=zero
    Do spc = 1,nspc
      natoms = molecules_getnatoms(spc)
      nmoles = config_getnmoles(species,spc)
      Call config_getconfig(species,spc,single)

      !** Skip this species if it is fixed in place
      If (config_isfixed(single)) Cycle

      !** Predict the positions and velocities
      Do m = 1, nmoles


        Do a = 1, natoms
          Call gear6_predict(single%coords(a,m)%rp,single%coords(a,m)%v, &
              g6%drvs(spc)%hot(1:5,a,m),g6%c,.True.)
          Call gear6_predict(single%coords(a,m)%r,single%coords(a,m)%v, &
              g6%drvs(spc)%hot(1:5,a,m),g6%c)
          single%afast(a,m) = 0.0_RDbl
          mass = molecules_AtomMass(spc,a)
          mvsq(spc) = mvsq(spc) + &
              single%coords(a,m)%v*single%coords(a,m)%v*mass

        End Do

        ! update current maximum velocity based on first atom
        first_atom_vel=mag(single%coords(1,m)%v)
        If (first_atom_vel>max_v) max_v=first_atom_vel


      End Do

      !** Do any necessary post-integration adjustments
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(spc),&
          simcell,single)

    End Do

    !** Predict the values for the extended Nose-Hoover coordinates
    !    Call nosehoover_getcoord(g6%thermostat%nh,nhsorb)
    nhsorb=>g6%thermostat%nh%spc
    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle
      Call gear6_predict(nhsorb(spc)%s%r, nhsorb(spc)%s%v, &
          g6%nhSpcDrvs(spc,1)%hot(1:5), g6%c,.False.)
    End Do

    !** Calculate the intermediate forces
    TKelvin = thermostats_T(g6%thermostat)
    nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
    If (.Not.nrg_success) Then
      success=nrg_success
      Return
    Endif

    ! potential nrg of color field, not really needed, but helpful in checking 
    ! hamiltonian conservation
    Call nemd_updateNrg(species)

    ! update flux values for reporting on the screen
    Call nemd_updateFlux(species,simcell,max_v*dt)

    !** Correction-Constraints Loop, MAX 2 corrections
    !** We check the penalty function after the correction and loop
    Do c = 1, MAX_CORRECTION
      Do spc = 1, nspc
        add_nemdforce=.False.
        If (NEMDParams%sorb==spc) add_nemdforce=.True.
        natoms = config_getnatoms(species,spc)
        nmoles = config_getnmoles(species,spc)
        Call config_getconfig(species,spc,single)

        !** Skip this species if it is fixed in place
        If (config_isfixed(single)) Cycle

        !** Correct the positions and velocities using the new forces
        Do m = 1, nmoles
          Do a = 1, natoms
            If (fast) Then

              If (c==1) Then
                !** Update accels with NEMD force also, so we add nemd accels 
                !** to regular accels, so this is needed only once when c==1
                ! units : mass = amu
                !         scalef = (amu. ang^2/ ps^2) / (Kcal/mol) 
                !         mult = (ang^2/ps^2)/(Kcal/mol)
                ! nemd_force*mult = ang/ps^2
                mass=molecules_AtomMass(spc,a)  
                mult=(scalef*one/mass)          ! 

                ! ** keep stats about NEMD force / regular force
                ! ** only done for first 50,000 iters
                If (check_mags.and.add_nemdforce) Then
                  ratio= 100* mag((NEMDParams%force(a)*mult)) &
                      / mag( single%afast(a,m))
                  mag_ratio=(mag_ratio*mag_count+ratio)/(mag_count+1)
                  mag_count=mag_count+1
                  If (Mod(mag_count,10000)==0) Then
                    Write(*,*) "-------------- Ratio of NEMD/Reular accels --"
                    Write(*,'(a,f15.8)') " In Percentage", mag_ratio
                    If (mag_count==500000) check_mags=.False.
                  Endif
                Endif


                ! NEMDParams%force(a) = force from color field in Kcal/mol/ang.
                If (add_nemdforce) &
                    single%afast(a,m) = single%afast(a,m) + &
                    (NEMDParams%force(a)*mult)
              Endif

              correct = (single%afast(a,m) - &
                  single%coords(a,m)%v*nhsorb(spc)%s%v- &
                  g6%drvs(spc)%hot(2,a,m))
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(spc)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.True.,correct)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(spc)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.False.,correct)
            Else
              Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
              Stop
            End If
          End Do
        End Do
      End Do


      !** Correct using the accelerations of "ln(s)" the extended system 
      !** coordinate
      Do spc=1,nspc
        If (config_isfixed(species(spc))) Cycle
        acc_nh = (mvsq(spc) - g6%spc_g(spc)/beta) / &
            (g6%thermostat%nh%Q*scalef)
        cortf=acc_nh-g6%nhSpcDrvs(spc,1)%hot(2)
        Call gear6_correct(nhsorb(spc)%s%r, nhsorb(spc)%s%v, &
            g6%nhSpcDrvs(spc,1)%hot(1:5), acc_nh,g6%corr2, .False., cortf)
      End Do

      Do spc = 1, nspc

        !** Skip this species if it is fixed in place
        If (config_isfixed(species(spc))) Cycle

        !** Call the constraints function
        If (.Not.interact_hasint(subint%imodel,spc,'constraint')) Cycle
        uconst = 0.0_RDbl
        nmoles = config_getnmoles(species,spc)
        natoms = molecules_getnatoms(spc)
        Call config_getconfig(species,spc,single)
        Do m = 1, nmoles
          If (Mod(steps,g6%penaltyStep) == 0) Then
            !** Call the penalty function
            Call intramolecular_getpenalty(subint%imodel%ff(1)%iparams(spc), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,sumobjf)
          End If
          If (fast) Then
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(spc), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,u,single%afast(1:natoms,m),.False.)
            uconst = u + uconst
          Else
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
          End If
        End Do
        Call storestats_updatenrg(subint%imodel%spcstats,spc,'const',uconst)

      End Do
    End Do

    !** Do any necessary post-integration adjustments
    Do spc = 1,nspc
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(spc), &
          simcell,species(spc))
    End Do

    !** Update energies, Note that here "s" refers to actually log(s) 
    !** in the original equations 
    ke=zero
    pe=zero
    Do spc = 1,nspc
      If (config_isfixed(species(spc))) Cycle
      ke = ke + (0.5_RDbl*g6%thermostat%nh%Q)*( (nhsorb(spc)%s%v)**2) 
      pe = pe + g6%spc_g(spc)*(nhsorb(spc)%s%r)* TKelvin* kjmole_kb
    End Do
    Call thermostats_updateEnergy(g6%thermostat,'ke',ke*scalepe)
    Call thermostats_updateEnergy(g6%thermostat,'pe',pe)

  End Subroutine gear6_MultNhNEMD



  !----------------------------------------------------------------------------
  ! This routine performs a single integration step of the 6th order
  ! Gear algorithm using the Nose-Hoover thermostat
  ! Here the NH extra variables are integrated using the second order methods
  ! predictor-corrector equations
  ! IT is not being currently used, If needed it can be substituted for 
  ! gear6_nhintegrate
  !            success -- indicates whether integration was successful
  !----------------------------------------------------------------------------
  Subroutine gear6_nhIntegrate2nd(subint, g6, simcell, species, dt,&
      fast, step, success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(GearParams), Intent(InOut)                :: g6
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(Kind=RDbl), Intent(InOut)                 :: dt
    Logical, Intent(In)                            :: fast
    Integer, Intent(In)                            :: step
    Logical, Intent(out)                           :: success

    Type(AtMolCoords), Pointer :: single
    Type(NHSorbates), Pointer :: nhsorb
    Integer :: nspc, i, m, a, natoms, nmoles, c
    Real(Kind=RDbl) :: u, sumobjf, uconst, mass
    Real(Kind=RDbl) :: mvsq   ! mass*vel**2 - for the Nose-Hoover xi term
    Real(Kind=RDbl) :: cortf, cors ! correction term for Nose-Hoover coords
    Real(kind=RDbl) :: ke, pe   ! Nose-Hoover hamiltonian
    Logical :: mapflag,nrg_success
    Integer, Save :: steps
    Real(Kind=RDbl) :: beta,acc_nh, TKelvin
    Type(VecType) :: correct

    mvsq   = 0.0_RDbl

    !** Beta has units of ps^2/Ang^2
    beta = 1.0_RDbl / (thermostats_T(g6%thermostat)*kcalmole_kb*scalef)

    steps = steps + 1
    success=.True.
    nspc = molecules_getnsorbs()

    Do i = 1,nspc

      natoms = molecules_getnatoms(i)
      nmoles = config_getnmoles(species,i)

      ! get the species pointer
      Call config_getconfig(species,i,single)

      ! Skip this species if it is fixed in place
      If (config_isfixed(single)) Cycle

      !** Predict the positions and velocities
      Do m = 1, nmoles
        Do a = 1, natoms
          !** Predicts the positions of the real atoms using taylor expnasion
          Call gear6_predict(single%coords(a,m)%rp,single%coords(a,m)%v, &
              g6%drvs(i)%hot(1:5,a,m),g6%c,.True.)
          Call gear6_predict(single%coords(a,m)%r,single%coords(a,m)%v, &
              g6%drvs(i)%hot(1:5,a,m),g6%c)
        End Do
      End Do

      !** Do any necessary post-integration adjustments, takes care of 
      !** ciccotti constraints
      Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i), &
          simcell,single)
    End Do

    !** Calculate the mvsq separately,
    mvsq=zero
    Do i = 1,nspc
      natoms = molecules_getnatoms(i)
      nmoles = config_getnmoles(species,i)
      Call config_getconfig(species,i,single)

      !** Skip this species if it is fixed in place
      If (config_isfixed(single)) Cycle

      !** Sum the momentum over mass
      !** This probably should be a routine in Nose-Hoover
      Do m = 1, nmoles
        Do a = 1, natoms
          mass = molecules_AtomMass(i,a)
          mvsq = mvsq + single%coords(a,m)%v*single%coords(a,m)%v*mass
        End Do
      End Do
    End Do

    !** Predict the values for the extended Nose-Hoover coordinates
    Call nosehoover_getcoord(g6%thermostat%nh,nhsorb)
    Call gear6_predict(nhsorb%s%r,nhsorb%s%v, g6%nhDrvs(1)%hot(1:5), &
        g6%c,.False.)

    !** Calculate the intermediate forces
    !!    betaKcal = 1.0_RDbl / (thermostats_T(g6%thermostat)*kcalmole_kb)
    TKelvin = thermostats_T(g6%thermostat)
    nrg_success = subinteract_int(subint,species,simcell, &
        fast,.True.,.False.,(/1000.0_RDbl,dt,TKelvin/),(/0,0,0/))
    If (.Not.nrg_success) Then
      success=nrg_success
      Return
    Endif

    !** Correction-Constraints Loop, Maximum 2 corrections
    !** We check the penalty function after the correction and loop
    Do c = 1, MAX_CORRECTION
!!$      temp_mvsq=zero
      Do i = 1, nspc

        natoms = config_getnatoms(species,i)
        nmoles = config_getnmoles(species,i)
        Call config_getconfig(species,i,single)

        !** Skip this species if it is fixed in place
        If (config_isfixed(single)) Cycle

        !** Correct the positions and velocities using the new forces
        Do m = 1, nmoles
          Do a = 1, natoms
            If (fast) Then
              correct = (single%afast(a,m) - single%coords(a,m)%v*nhsorb%s%v- &
                  g6%drvs(i)%hot(2,a,m))
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.True.,correct)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%afast(a,m), &
                  g6%corr2,.False.,correct)

!!$              !** Note temp_mvsq is not used currently.  Might need it for 
!!$              !** greater accuracy at a later time.
!!$              mass = atom_getmass(molecules_getatype(i,a))
!!$              temp_mvsq = temp_mvsq + single%coords(a,m)%v* &
!!$              single%coords(a,m)%v*mass
            Else
              correct = (single%aslow(a,m)-single%coords(a,m)%v*nhsorb%s%v- &
                  g6%drvs(i)%hot(2,a,m))
              Call gear6_correct(single%coords(a,m)%rp,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                  g6%corr2,.True.,correct)
              Call gear6_correct(single%coords(a,m)%r,single%coords(a,m)%v, &
                  g6%drvs(i)%hot(1:5,a,m),single%aslow(a,m), &
                  g6%corr2,.False.,correct)
            End If
          End Do
        End Do
      End Do

      !** Correct using the accelerations of "ln(s)" the extended system 
      !** coordinate
      acc_nh = (mvsq - Real(g6%g,kind=RDbl)/beta)/(g6%thermostat%nh%Q*scalef)
      cortf=acc_nh-g6%nhDrvs(1)%hot(2)
      Call gear6_correct(nhsorb%s%r, nhsorb%s%v,g6%nhDrvs(1)%hot(1:5), &
          acc_nh,g6%corr2,.False.,cortf)

      Do i = 1, nspc
        !** Skip this species if it is fixed in place
        If (config_isfixed(species(i))) Cycle

        !** Call the constraints function
        If (.Not.interact_hasint(subint%imodel,i,'constraint')) Cycle
        uconst = 0.0_RDbl
        nmoles = config_getnmoles(species,i)
        natoms = molecules_getnatoms(i)
        Call config_getconfig(species,i,single)
        Do m = 1, nmoles
          !** Call the penalty function for evansmorriss
          If (Mod(steps,g6%penaltyStep) == 0) Then
            !** Call the penalty function
            Call intramolecular_getpenalty(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,sumobjf)
          End If

          If (fast) Then
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,u,single%afast(1:natoms,m),.False.)
            uconst = u + uconst
          Else
            Call intramolecular_mconint(subint%imodel%ff(1)%iparams(i), &
                single%coords(1:natoms,m)%rp, &
                single%coords(1:natoms,m)%v,u,single%aslow(1:natoms,m),.False.)
            uconst = u + uconst
          End If

        End Do

        Call storestats_updatenrg(subint%imodel%spcstats, i,'const',uconst)
        Call intramolecular_postadjust(subint%imodel%ff(1)%iparams(i), &
            simcell,single)
      End Do
    End Do

    !** Update energies, Note that here "s" refers to actually log(s) 
    !** in the original equations 
    ke = (0.5_RDbl*g6%thermostat%nh%Q)*( (nhsorb%s%v)**2) 
    Call thermostats_updateEnergy(g6%thermostat,'ke',ke*scalepe)
    pe = Real(g6%g,kind=RDbl)*(nhsorb%s%r)*thermostats_T(g6%thermostat)* &
        kjmole_kb

    Call thermostats_updateEnergy(g6%thermostat,'pe',pe)

  End Subroutine gear6_nhIntegrate2nd


  !----------------------------------------------------------------------------
  ! Predicts the positions and optionally the velocities and HOTs
  !----------------------------------------------------------------------------
  Subroutine gear6_vecPredict(pos,vel,hot,c,onlyPos)
    Type(VecType), Intent(InOut) :: pos,vel  ! position, velocity
    Real(kind=RDbl), Dimension(:), Intent(In) :: c
    Type(VecType), Dimension(:), Intent(InOut) :: hot !Higher Order Terms
    Logical, Intent(In), Optional :: onlyPos
    Logical :: posOnly
    Integer :: i

    If (Present(onlyPos)) Then
      posOnly = onlyPos
    Else
      posOnly = .False.
    End If

    ! c indicies are offset by 1 when it is passed. Correct it.
    i = 1

    pos = pos + vel*c(1+i) + hot(2)*c(2+i) + hot(3)*c(3+i) + hot(4)*c(4+i) & 
        + hot(5)*c(5+i)

    If (.Not.posOnly) Then
      vel = vel + hot(2)*c(1+i) + hot(3)*c(2+i) + hot(4)*c(3+i) + hot(5)*c(4+i)
      hot(2) = hot(2) + hot(3)*c(1+i) + hot(4)*c(2+i) + hot(5)*c(3+i)
      hot(3) = hot(3) + hot(4)*c(1+i) + hot(5)*c(2+i)
      hot(4) = hot(4) + hot(5)*c(1+i)
    End If
  End Subroutine gear6_vecPredict

  !----------------------------------------------------------------------------
  ! Predicts the positions and optionally the velocities and HOTs
  !----------------------------------------------------------------------------
  Subroutine gear6_sclrPredict(pos,vel,hot,c,onlyPos)
    Real(kind=RDbl), Intent(InOut) :: pos,vel  ! position, velocity
    Real(kind=RDbl), Dimension(:), Intent(In) :: c
    Real(kind=RDbl), Dimension(:), Intent(InOut) :: hot !Higher Order Terms
    Logical, Intent(In), Optional :: onlyPos
    Logical :: posOnly
    Integer :: i

    If (Present(onlyPos)) Then
      posOnly = onlyPos
    Else
      posOnly = .False.
    End If

    ! c indicies are offset by 1 when it is passed. Correct it.
    i = 1

    pos = pos + vel*c(1+i) + hot(2)*c(2+i) + hot(3)*c(3+i) + hot(4)*c(4+i) & 
        + hot(5)*c(5+i)

    If (.Not.posOnly) Then
      vel = vel + hot(2)*c(1+i) + hot(3)*c(2+i) + hot(4)*c(3+i) + hot(5)*c(4+i)
      hot(2) = hot(2) + hot(3)*c(1+i) + hot(4)*c(2+i) + hot(5)*c(3+i)
      hot(3) = hot(3) + hot(4)*c(1+i) + hot(5)*c(2+i)
      hot(4) = hot(4) + hot(5)*c(1+i)
    End If
  End Subroutine gear6_sclrPredict

  !----------------------------------------------------------------------------
  ! Corrects the positions and optionally the velocities and HOTs
  !----------------------------------------------------------------------------
  Subroutine gear6_vecCorrect(pos,vel,hot,acc,mcorr,onlyPos,corr)
    Type(VecType), Intent(InOut) :: pos,vel,acc
    Type(VecType), Intent(InOut), Dimension(:) :: hot
    Real(Kind=Rdbl), Intent(In), Dimension(:) :: mcorr !Magic corrector terms
    Logical, Optional, Intent(In):: onlyPos
    Type(VecType), Intent(In), Optional :: corr
    Logical :: posOnly
    Type(VecType) :: correct
    Integer :: i,j

    If (Present(corr)) Then
      correct = corr
    Else
      correct = acc - hot(2)
    End If

    If (Present(onlyPos)) Then
      posOnly = onlyPos
    Else
      posOnly = .False.
    End If

    j = 1  ! correction constants index runs from 1:6, just to confuse

    pos = pos + correct*mcorr(0+j)

    If (.Not.posOnly) Then
      vel = vel + correct*mcorr(1+j)
      Do i = 2,5
        hot(i) = hot(i) + correct*mcorr(i+j)
      End Do
    End If

  End Subroutine gear6_vecCorrect

  !----------------------------------------------------------------------------
  ! Corrects the positions and optionally the velocities and HOTs
  !----------------------------------------------------------------------------
  Subroutine gear6_sclrCorrect(pos,vel,hot,acc,mcorr,onlyPos,corr)
    Real(kind=RDbl), Intent(InOut) :: pos,vel,acc
    Real(kind=RDbl), Intent(InOut), Dimension(:) :: hot
    Real(Kind=Rdbl), Intent(In), Dimension(:) :: mcorr !Magic corrector terms
    Logical, Optional, Intent(In):: onlyPos
    Real(Kind=RDbl), Intent(In), Optional :: corr
    Logical :: posOnly
    Real(Kind=RDbl) :: correct
    Integer :: i,j

    If (Present(corr)) Then
      correct = corr
    Else
      correct = acc - hot(2)
    End If

    If (Present(onlyPos)) Then
      posOnly = onlyPos
    Else
      posOnly = .False.
    End If

    j = 1  ! correction constants index runs from 1:6, just to confuse

    pos = pos + correct*mcorr(0+j)

    If (.Not.posOnly) Then
      vel = vel + correct*mcorr(1+j)
      Do i = 2,5
        hot(i) = hot(i) + correct*mcorr(i+j)
      End Do
    End If

  End Subroutine gear6_sclrCorrect

  !----------------------------------------------------------------------------
  ! Checks the constraints for the molecule
  !----------------------------------------------------------------------------
  Subroutine gear6_checkConst()

  End Subroutine gear6_checkConst

  !----------------------------------------------------------------------------
  ! Add any additional energy information to pe, ke, and their averages
  ! All energies are in returned kJ/mol
  !----------------------------------------------------------------------------
  Subroutine gear6_extraEnergy(g6,nrg,avgNrg,eType)
    Type(GearParams), Intent(In) :: g6
    Real(Kind=Rdbl), Intent(Out) :: nrg, avgNrg
    Character(*), Intent(In) :: eType

    !** Check to see if there is extra energy info in the thermostat
    If (Associated(g6%thermostat)) Then
      Select Case(tolower(eType))
      Case ('pe')
        nrg = thermostats_getpe(g6%thermostat,'inst')
        avgNrg = thermostats_getpe(g6%thermostat,'cavg')
      Case ('ke')
        nrg = thermostats_getke(g6%thermostat,'inst')
        avgNrg = thermostats_getke(g6%thermostat,'cavg')
      case default
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End Select
    Else
      nrg = 0.0_Rdbl
      avgNrg = 0.0_RDbl
    End If

    !** Add energies from color-field MD if any
    !** In color field energies are stored in kcal/mol; 
    If (NEMD_IS_ON) Then
      Select Case(tolower(eType))
      Case ('pe')
        nrg = nrg+ nemd_getpe('inst')*caltoj
        avgNrg = avgNrg + nemd_getpe('cavg')*caltoj
      Case ('ke')
      Case Default
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End Select
    Endif

  End Subroutine gear6_extraEnergy

  !----------------------------------------------------------------------------
  ! Displays the information about the Gear algorithm
  !----------------------------------------------------------------------------
  Subroutine gear6_display(g6,unitno)
    Type(GearParams), Intent(In) :: g6
    Integer, Intent(In) :: unitno

    Write(unitno,'(6x,a,i4)') "Steps between penalty function : ", &
        g6%penaltyStep
    If (Associated(g6%thermostat)) Then
      Call thermostats_display(g6%thermostat,unitno,8)
    End If
  End Subroutine gear6_display

End Module gear6

