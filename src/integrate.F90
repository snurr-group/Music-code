!------------------------------------------------------------------------------
! This module initializes and dishes out work to the correct integrator.
! It is called on a multiple timestep level and handles its own output.
! It contains the data structure "Integrate_Params" which contains
! both general integration parameters and those specific to the integrator
! being used.  The specific parameters take the form of a pointer to
! a data structure unique to the integrator.  Defining a pointer specifies
! which integrator will be used.  This module also contains routines for
! displaying the results of an on-going integration and for calculating 
! calculating quantities such as the kinetic temperature.
! 
! Needed Improvements:
! 1) better requires statements and comments
!------------------------------------------------------------------------------

Module integrate

  Use defaults, Only: RDbl,strLen,scaleke,kjmole_kb,zero,scalepe,TOTAL_INDEX
  Use file, Only: file_getunit, file_getname
  Use utils, Only: stripcmnt,toupper,split, int2str, allocErrDisplay
  Use datafile,Only: CONFILE, datafile_writeconfig
  Use vector, Only: VecType, Assignment(=), Operator(+), &
      Operator(-), Operator(*), Operator(/), mag
  Use velverlet, Only: VelVerletParams, velverlet_init, velverlet_integrate, &
      velverlet_extraenergy, velverlet_display, velverlet_sampleCF
  Use gear6, Only: GearParams, gear6_init, gear6_integrate, & 
      gear6_extraenergy, gear6_display, gear6_sampleCF, gear6_initIntegrate
  Use leapfrog, Only: LeapFrogParams, leapfrog_init, leapfrog_integrate, &
      leapfrog_extraenergy, leapfrog_display, leapfrog_sampleCF
  Use atom, Only: atom_getmass
  Use config, Only: AtMolCoords,config_getnmoles,config_getconfig, &
      config_isfixed, config_makeconsistent
  Use simcell, Only: SimCell_Params,simcell_pbc
  Use molecules, Only: molecules_getnsorbs, molecules_getnatoms, &
      molecules_getnthatom, molecules_getdof, molecules_getmass, &
      molecules_name, molecules_getnatomtypes
  Use storestats, Only: storestats_gettemp, &
      storestats_getintranrg,storestats_getke, storestats_updateEnergySS, &
      storestats_getnoncoul, storestats_getcoul, &
      storestats_updatetemp,storestats_updatenrg, storestats_updateIntraEnergy
  Use interact, Only: Interaction_Model, interact_hasint, &
      interact_dispDynInfo
  Use subinteract, Only: Subset_Interactions, subinteract_int, &
      subinteract_simpleupdate
  Use storetop, Only: storetop_display
  Use movie, Only: movie_init, movieInfo, movie_makemovie, movie_display 

  Implicit None
  Save

  Private
  Public :: Integrate_Params, integrate_init, &
      integrate_simdisplay,integrate_displayparams,integrate_sampleCF, &
      integrate_simdisplayExtra, integrate_getdt, integrate_integrate, &
      integrate_integrateMC

  Integer :: counterc
  !***** MOVIE STUFF
  Real(kind=RDbl) :: time !SHAJI_DEBUG...............
  !** Movie Parameters, for those who love to watch molecules
  Type(MovieInfo) :: movieParams
  Character(len=strLen), Parameter :: movie_tag = "Movie Information"
  Integer:: iter1

  Type Integrate_Params
    Real(kind=RDbl) :: simT      !sim temp
    Real(kind=RDbl) :: timeStep  !delta time, picoseconds
    Integer :: iterations        !no of itns
    Integer :: ioStep       
    Integer :: configStep
    Integer :: slowStep
    Character(len=strLen) :: ensemble
    Character(len=strLen) :: integrator
    Type(GearParams), Pointer      :: gear6
    Type(VelVerletParams), Pointer :: vverlet
    Type(LeapFrogParams), Pointer  :: lfrog
  End Type Integrate_Params

Contains

  !----------------------------------------------------------------------------
  ! Initializes the correct pointer for the given integrator
  !----------------------------------------------------------------------------
  Subroutine integrate_init(integParams, sorbates, filename)
    Type(Integrate_Params) :: integParams
    Type(AtMolCoords), Intent(In), Dimension(:) :: sorbates
    Character(*) ,Intent(in):: filename

    Integer :: error,j,unitno
    Character(len=255) :: text
    Character(len=strLen), Dimension(strLen) :: params

    unitno=file_getunit(filename)
    Read(unitno,*) integParams%iterations
    Read(unitno,*) integParams%timeStep
    Read(unitno,*) integParams%simT
    Read(unitno,*) integParams%ioStep
    Read(unitno,*) integParams%configStep
    counterc=0
    !** Get the ensemble for the move
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    j = split(text,params)
    integParams%ensemble = Trim( toupper( params(1) ) )

!!$    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__    
!!$    Write(*,*) "ensemble : ",integParams%ensemble ,Trim(toupper(params(1)))

    !** Get the type of integrator to use
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    j = split(text,params)
    integParams%integrator = Trim( toupper( params(1) )  )

    !** Initialize the integrator
    !LC    Write(0,*) "Initializing integrator: ",params(1)
    Select Case (integParams%integrator)

    Case ("VELOCITYVERLET")
      Allocate(integParams%vverlet,stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Nullify(integParams%gear6)
      Nullify(integParams%lfrog)
      If (Trim(integParams%ensemble)=="NVT") Then 
        Call velverlet_init(integParams%vverlet,integParams%ensemble &
            ,unitno,integParams%simT)
      Else If (Trim(integParams%ensemble)=="NVE") Then 
        Call velverlet_init(integParams%vverlet,integParams%ensemble,unitno)
      Else
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(2a)') "unable to identify ensemble: ",integParams%ensemble
        Stop
      End If

    Case ("LEAPFROG")
      Allocate(integParams%lfrog,stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Nullify(integParams%gear6)
      Nullify(integParams%vverlet)

      If (Trim(integParams%ensemble)=="NVT") Then 
        Call leapfrog_init(integParams%lfrog, sorbates, &
            integParams%ensemble, unitno, integParams%simT)
      Else If (Trim(integParams%ensemble)=="NVE") Then 
        Call leapfrog_init(integParams%lfrog, sorbates, &
            integParams%ensemble, unitno)
      Else
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(2a)') "unable to identify ensemble: ",integParams%ensemble
        Stop
      End If


    Case ("GEAR6")
      Allocate(integParams%gear6,stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Nullify(integParams%vverlet)
      Nullify(integParams%lfrog)
      If (Trim(integParams%ensemble)=="NVT") Then 
        Call gear6_init(integParams%gear6,sorbates,integParams%timeStep,&
            integParams%ensemble,unitno,integParams%simT)
      Else If (Trim(integParams%ensemble)=="NVE") Then 
        Call gear6_init(integParams%gear6,sorbates,integParams%timeStep,&
            integParams%ensemble,unitno)
      Else
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(2a)') "unable to identify ensemble: ",integParams%ensemble
        Stop
      End If

    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          "Don't know how to handle integrator ",Trim(integParams%integrator)
      Write(0,'(a)') "Allowed Types : VelocityVerlet, Gear6"
      Stop

    End Select

#ifdef DEBUG
    !** Initialize the movie
    Call movie_init(movieParams,file_getname(unitno),movie_tag)
    Call movie_display(movieParams,6,6)
    time=zero   ! for movie_display
    iter1=0
#endif

  End Subroutine integrate_init

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine integrate_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a,t30,a)') 'Integer', &
        '# Steps for this move type per iteration'
    Write(unitno,'(a,t30,a)') 'Real','# Time step, fs'
    Write(unitno,'(a,t30,a)') 'Real','# Temperature, K'
    Write(unitno,'(a,t30,a)') 'Integer','# Steps between writes to I/O'
    Write(unitno,'(a,t30,a)') 'Integer', &
        '# Steps between writes to configuration file'
    Write(unitno,'(a,t30,a)') '[NVE,NVT]','# Ensemble to simulate in'
    Write(unitno,'(a,t30,a)') '[VELOCITYVERLET,GEAR6]', &
        '# Integrator name'
    Write(unitno,'(a)')  &
        '# If GEAR6 was specified, remove this line and include:'
    Call gear6_sampleCF(unitno)
    Write(unitno,'(a)') '# End GEAR6 section (remove this line)'
    Write(unitno,'(a)') &
        '# If VELOCITYVERLET, remove this line and include:'
    Call velverlet_sampleCF(unitno)
    Write(unitno,'(a)') '# End VELOCITYVERLET section (remove this line)'

  End Subroutine integrate_sampleCF

  !----------------------------------------------------------------------------
  ! make an md-move, dump to io, repeat
  ! subint -- subint%temp will get energies of system after integrate.
  ! Note : If the main program decides to point to 
  !          subint%temp to imodel%ff%results then old energies will be 
  !          destroyed. In case of MD that's OK.
  ! Requires:  params -- integration move parameters
  !            subint -- subset interactions for species
  !            species -- the species coordinates
  !            simcell -- simulation cell information
  !            succ_flag -- a flag indicating success or failure
  !----------------------------------------------------------------------------
  Subroutine integrate_integrateMC(params, subint, species, simcell, &
      succ_flag)
    Type(Integrate_Params), Intent(InOut)          :: params
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Logical, Intent(out)                           :: succ_flag 

    Integer         :: i
    Logical         :: fast, updateFlag, success
    Real(kind=RDbl) :: ltime, dt

    fast = .True.
    ltime = 0.0_RDbl  !MC time
    updateFlag = .False.
    succ_flag = .True.

    !** initialize the hot fields of gear6
    If (Associated(params%gear6)) Then
      Call gear6_InitIntegrate(subint,params%gear6,simcell,species, &
          params%timeStep,fast,success)
    End If

    Do i = 1, params%iterations

      !** Set the default time step
      dt = params%timeStep

      Call integrate_OneStep(subint, params,simcell,species, &
          fast, i, dt, updateFlag,succ_flag)
      If (.Not. succ_flag) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            'Integration step failure -- ', &
            'probably unable to evaluate interactions'
        Return
      End If

      !** Update the simulation time
      ltime = ltime + dt

      !** Write to IO if required
      If ((Mod(i,params%ioStep) == 0).And.(params%iostep /= 0)) Then
        Call integrate_simdisplay(subint%imodel,params,species,i,ltime)
      End If

      !** Dump to the config file
      If ((params%configStep /= 0)) Then
        If (Mod(i,params%configStep) == 0) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "This needs more coding"
          Stop
        Endif
      End If

    End Do

    !** make sure generalized coordinates are consistent with %rp
    Call config_makeconsistent(species)

    !** make sure all energies are up to date. This is required
    ! because last step of gear_integrate is done after nrg
    ! calculation. (i.e. The particles are moved after calculating the
    ! energy. Also the constraints might alter positions a bit
    ! This check is only needed for -VERBOSE  option to work
    success = subinteract_int(subint,species,simcell, &
        .True.,.True.,.False.,(/1000.0_RDbl,dt,0.00_RDbl/),(/0,0,0/))
    If (.Not. success) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

  End Subroutine integrate_integrateMC

  !----------------------------------------------------------------------------
  ! Make a single Integrate move
  ! ** NOTE ** Fix calls to velverlet
  !----------------------------------------------------------------------------
  Subroutine integrate_OneStep(subint, params,simcell,sorbates,fast,&
      step,dt,updateFlag,success)
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Type(Integrate_Params), Intent(InOut)          :: params
    Type(AtMolCoords), Dimension(:)                :: sorbates
    Type(SimCell_Params), Intent(In)               :: simcell
    Logical, Intent(In)                            :: fast
    Integer, Intent(In)                            :: step
    Real(Kind=RDbl), Intent(InOut)                 :: dt       
    Logical, Intent(in)                            :: updateFlag
    Logical, Intent(out)                           :: success

    Integer                            :: i,j,natoms,nmoles,nsorbs
    Real(kind=RDbl)                    :: tmole, tatom, ke

    !** Initialize variables
    nsorbs = molecules_getnsorbs()

    !** Decide which integrator to use
    If (Associated(params%vverlet)) Then
      Call velverlet_integrate(subint, params%vverlet,simcell,sorbates,&
          params%ensemble,dt,fast,success)
    Else If (Associated(params%gear6)) Then
      Call gear6_integrate(subint,params%gear6,simcell,sorbates, &
          params%ensemble,dt,fast,step,success)
    Else If (Associated(params%lfrog)) Then
      Call leapfrog_integrate(subint,params%lfrog,simcell,sorbates, &
          params%ensemble,dt,fast,success)
    End If

    If (.Not.success) Return

    !** Loop over sorbates and perform necessary operations
    !** Apply pbc's, calculate and update temperatures
    Do i = 1,nsorbs
      !** Skip if this molecule has a FIXED configuration
      If (config_isfixed(sorbates(i))) Cycle

      !** Get the number of molecules and atoms
      natoms = molecules_getnatoms(i)
      nmoles = config_getnmoles(sorbates,i)

      !** Loop through the molecules and apply periodic boundary conditions
      Do j = 1,nmoles
        Call simcell_pbc(simcell,sorbates(i)%coords(1:natoms,j)%rp, &
            sorbates(i)%coords(1:natoms,j)%r)
      End Do

      !** this update can be done here or in the driving program.
      If (updateFlag) Then
!       Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Call subinteract_simpleupdate(subint,sorbates)
      End If

    End Do

#ifdef DEBUG
    !** Dump to movie file
    iter1 = iter1 + 1
    time = time + dt
    Call movie_makemovie(movieParams,sorbates,simcell,iter1,time)
#endif

  End Subroutine integrate_OneStep

  !----------------------------------------------------------------------------
  ! Do the MD simulation, i.e., make a move, dump to io, repeat
  ! imodel -- interaction model for forcefield, [also stores previous 
  !           system energies]
  ! subint -- subint%temp will get energies of system after integrate.
  ! Note : If the main program decides to point 
  !          subint%temp to imodel%ff%results then old energies will be 
  !          destroyed. In case of MD that's OK.
  !----------------------------------------------------------------------------
  Subroutine integrate_integrate(subint, params, simcell, species, &
      conf_file, updateFlag, time )
    Type(Integrate_Params), Intent(InOut)          :: params
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(CONFILE), Intent(in)                      :: conf_file
    Logical, Intent(In)                            :: updateFlag
    Type(Subset_Interactions), Intent(InOut)       :: subint    
    Real(kind=RDbl), Intent(InOut), Optional       :: time

    Integer         :: i
    Logical         :: fast,success
    Logical, Save   :: first_time = .True.
    Real(kind=RDbl) :: ltime, dt

    fast = .True.
    If (Present(time)) Then
      ltime = time
    Else
      ltime = 0.0_RDbl
    End If

    If (first_time) Then
      first_time = .False.
      !** initialize the hot fields of gear6
      If (Associated(params%gear6)) Then
        Call gear6_InitIntegrate(subint,params%gear6,simcell,species, &
            params%timeStep,fast,success)
      End If
    Endif
 
    Do i = 1, params%iterations

      !** Set the default time step
      dt = params%timeStep

      Call integrate_OneStep(subint, params,simcell,species, &
          fast, i, dt, updateFlag,success)
      If (.Not.success) Then
        Write(*,*) "problems encountered during integration"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif


      !** Update the simulation time
      ltime = ltime + dt

      !** Write to IO if required
      If ((Mod(i,params%ioStep) == 0).And.(params%iostep /= 0)) Then
        Call integrate_simdisplay(subint%imodel,params,species,i,ltime)
      End If

      !** Dump to the config file
      If ((params%configStep /= 0).And.(params%configStep /= 0)) Then
        If (Mod(i,params%configStep) == 0) Then
          counterc = counterc+1
          Call datafile_writeconfig(species,subint%imodel%spcstats,conf_file,ltime)
        End If
      End If
    End Do

    !----------- commented out for speed -------------------------------------
    !    !** make sure generalized coordinates are consistent with %rp
    !    Call config_makeconsistent(sorbates)
    !---- should be used in case this routine is used later by some md move --

    If (Present(time)) time = ltime

  End Subroutine integrate_integrate

  !----------------------------------------------------------------------------
  ! Calculate the kinetic temperature based either on the molecular or atomic 
  ! velocities
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function integrate_kineticTemp(molecule,vel,tType)
    Integer, Intent(In) :: molecule  ! This is the molecule type
    Type(VecType), Dimension(:,:), Intent(In) :: vel
    Character(*), Intent(In) :: tType
    Integer :: i,j,natoms,nmoles,atomt
    Real(kind=RDbl) :: ke, mmass, mmassi2,mass
    Type(VecType) :: vels

    nmoles = Size(vel,2)
    natoms = Size(vel,1)
    ke = 0.0_RDbl

    Select Case(tType)

    Case ('atomic')
      Do i = 1, nmoles
        Do j = 1, natoms
          atomt = molecules_getnthatom(molecule,j)
          mass = atom_getmass(atomt)
          ke = vel(j,i)*vel(j,i)*mass + ke
        End Do
      End Do
      integrate_kineticTemp = ke*scaleke / (kjmole_kb* &
          Real(molecules_getdof(molecule)*nmoles))

    Case ('molecular')
      mmass = molecules_getmass(molecule)
      mmassi2 = 1.0_RDbl / mmass / mmass
      Do i = 1,nmoles
        vels = 0.0_RDbl
        Do j = 1,natoms
          atomt = molecules_getnthatom(molecule,j)
          mass = atom_getmass(atomt)
          vels = vel(j,i)*mass + vels
        End Do
        ke = vels*vels*mmassi2 + ke
      End Do

      integrate_kineticTemp = ke*mmass*scaleke / (kjmole_kb*Real(3*nmoles))

    End Select

  End Function integrate_kineticTemp

  !----------------------------------------------------------------------------
  ! Calculate the total kinetic energy of the system
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function integrate_kineticEnergy(sorbates,sorb)
    Type(AtMolCoords), Intent(In), Dimension(:) :: sorbates
    Integer, Intent(In)                         :: sorb

    Integer         :: nmoles, natoms, atype, k, j
    Real(kind=RDbl) :: mass

    integrate_kineticEnergy = 0.0_RDbl
    nmoles = config_getnmoles(sorbates,sorb)
    natoms = molecules_getnatoms(sorb)

    If (.Not.Associated(sorbates(sorb)%coords)) Return
    Do k = 1, nmoles
      Do j = 1, natoms
        atype = molecules_getnthatom(sorb,j)
        mass = atom_getmass(atype)
        integrate_kineticEnergy = integrate_kineticEnergy + &
            sorbates(sorb)%coords(j,k)%v*sorbates(sorb)%coords(j,k)%v*mass
      End Do
    End Do

    integrate_kineticEnergy = 0.5_RDbl*integrate_kineticEnergy*scaleke

  End Function integrate_kineticEnergy

  !----------------------------------------------------------------------------
  ! Displays extra energy information for extended systems
  !----------------------------------------------------------------------------
  Subroutine integrate_extraEnergy(params,nrg,avgNrg,eType)
    Type(Integrate_Params), Intent(In) :: params
    Real(Kind=RDbl), Intent(Out) :: nrg, avgNrg
    Character(*), Intent(In) :: eType

    If (Associated(params%vverlet)) Then
      Call velverlet_extraEnergy(params%vverlet,nrg,avgNrg,eType)
    Else If (Associated(params%gear6)) Then
      Call gear6_extraEnergy(params%gear6,nrg,avgNrg,eType)
    End If

  End Subroutine integrate_extraEnergy

  !----------------------------------------------------------------------------
  ! Returns the value of the time step for the given integration move.
  !----------------------------------------------------------------------------
  Real(Kind=Rdbl) Function integrate_getdt(params)
    Type(Integrate_Params), Intent(In) :: params

    integrate_getdt = params%timeStep
  End Function integrate_getdt

  !----------------------------------------------------------------------------
  ! Calculate the center of mass velocity of the molecule
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function integrate_comVelocity(molecule,vel)
    Integer, Intent(In) :: molecule  
    Type(VecType), Dimension(:), Intent(In) :: vel

    Integer :: natoms, natypes, i 
    Integer, Dimension(Size(vel,1)) :: atypes
    Real(kind=RDbl) :: mass, comvel
    Type(VecType) :: comv

    !** Zero the velocity and mass
    comvel = 0.0_RDbl
    mass = 0.0_RDbl

    !** Get the number of atoms
    natoms = molecules_getnatoms(molecule)

    !** Get the atom types of the molecule
    natypes = molecules_getnatomtypes(molecule,atypes,.True.)

    !** Loop through the atoms and get the center of mass velocity
    Do i = 1, natoms
      comv = comv + vel(i)*atom_getmass(atypes(i))
      mass = mass + atom_getmass(atypes(i))
    End Do

    comv = comv/mass
    comvel = mag(comv)

    integrate_comVelocity = comvel

  End Function integrate_comVelocity

  !----------------------------------------------------------------------------
  ! Displays information about the integrator
  !----------------------------------------------------------------------------
  Subroutine integrate_displayparams(params,unitno)
    Type(Integrate_Params), Intent(In) :: params
    Integer, Intent(In) :: unitno

    Write(unitno,'(2x,2a)')    " Ensemble        : ", Trim(params%ensemble)
    If (params%ensemble=="NVT") Then
      Write(unitno,'(2x,a,f9.2)') " Temperature     : ",params%simT
    End If
    Write(unitno,'(2x,2a)') " Integrator Name : ", Trim(params%integrator)
    Write(unitno,'(2x,a,i9,10x,a,f12.5,a)') " No. Iterations  : ", &
        params%iterations, &
        " Time Step       : ",params%timeStep," ps"
    Write(unitno,'(2x,a,i9)') " slowStep        : ", &
        params%slowStep
    Write(unitno,'(2x,a,i9,10x,a,i9)')  " Writes to I/O   : ", &
        params%ioStep," Writes to Config: ",params%configStep
    If (Associated(params%vverlet)) Then
      Write(unitno,'(6x,a)') "Velocity Verlet Integrator Information"
      Call velverlet_display(params%vverlet,unitno)
    Else If (Associated(params%gear6)) Then
      Write(unitno,'(6x,a)') "6th Order Gear Integrator Information"
      Call gear6_display(params%gear6,unitno)
    Else If (Associated(params%lfrog)) Then
      Write(unitno,'(6x,a)') "Leap Frog Integrator Information"
      Call leapfrog_display(params%lfrog,unitno)
    End If

  End Subroutine integrate_displayparams

  !----------------------------------------------------------------------------
  ! Displays information about the MD simulation in progress
  !----------------------------------------------------------------------------
  Subroutine integrate_simdisplay(imodel,params,sorbates,step,time,unit)
    Type(Interaction_Model), Intent(InOut)      :: imodel
    Type(Integrate_Params), Intent(In)          :: params
    Type(AtMolCoords), Dimension(:), Intent(In) :: sorbates
    Integer, Intent(In)                         :: step
    Real(kind=RDbl), Intent(In)                 :: time
    Integer, Intent(In), Optional               :: unit

    Integer          :: i,j,nsorbs,moles,nmoles
    Real(kind=RDbl)  :: pe, ke, peavg, keavg, scaletime
    Real(kind=RDbl)  :: pei, kei, peavgi, keavgi
    Character(len=4) :: units
    Integer          :: dUnit  ! display unit number
    Character(len=strLen) :: outString, outString2

    If (.Not.Present(unit)) Then
      dUnit = 6
    Else
      dUnit = unit
    End If

    If (time < 0.1) Then
      units = "fs"
      scaletime = time*1000.0_RDbl
    Else If (time > 1000.0) Then
      units = "ns"
      scaletime = time/1000.0_Rdbl
    Else
      units = "ps"
      scaletime = time
    End If

    pe = 0.0_Rdbl
    ke = 0.0_Rdbl
    peavg = 0.0_Rdbl
    keavg = 0.0_Rdbl
    moles = 0

    nsorbs = molecules_getnsorbs()

    Write(dUnit,'(a)') " "
    outString = int2str(step)
    Write(dUnit,'(1x,2a,4x,a,f10.3,1x,a)') "MD Step: ",Trim(outString), &
        "Time ",scaletime,units
    Write(dUnit,'(1x,2a)') "----------------------------------------------", &
        "----------"
    Do i = 1,nsorbs
      If (config_isfixed(sorbates(i))) Cycle
      nmoles = config_getnmoles(sorbates,i)
      If (nmoles == 0) Cycle
      outString = molecules_name(i)
      Write(dUnit,'(1x,a,a)') Trim(outString), &
          " Information"
      Write(dUnit,'(3x,a16,4a15)') "Variable","Current","CumulAvg", &
          "Block","Std"
      Write(dUnit,'(3x,a16,4f15.3)') "Tmole          :", &
          storestats_gettemp(imodel%spcstats, i,'tmole','inst'), &
          storestats_gettemp(imodel%spcstats, i,'tmole','cavg'), &
          storestats_gettemp(imodel%spcstats, i,'tmole','block'), &
          storestats_gettemp(imodel%spcstats, i,'tmole','std')
      Write(dUnit,'(3x,a16,4f15.3)') "Tatom          :", &
          storestats_gettemp(imodel%spcstats, i,'tatom','inst'), &
          storestats_gettemp(imodel%spcstats, i,'tatom','cavg'), &
          storestats_gettemp(imodel%spcstats, i,'tatom','block'), &
          storestats_gettemp(imodel%spcstats, i,'tatom','std')
      Write(dUnit,'(3x,a16,4f15.6)') "Intramolecular :", &
          storestats_getintranrg(imodel%spcstats, i,'inst')/nmoles, &
          storestats_getintranrg(imodel%spcstats, i,'cavg')/nmoles, &
          storestats_getintranrg(imodel%spcstats, i,'block')/nmoles, &
          storestats_getintranrg(imodel%spcstats, i,'std')/nmoles
      Write(dUnit,'(3x,a16,4f15.6)') "Kinetic Energy :", &
          storestats_getke(imodel%spcstats, i,'inst')/nmoles, &
          storestats_getke(imodel%spcstats, i,'cavg')/nmoles, &
          storestats_getke(imodel%spcstats, i,'block')/nmoles, &
          storestats_getke(imodel%spcstats, i,'std')/nmoles 
      If (interact_hasint(imodel,i,'stretch')) Then
        Write(dUnit,'(3x,a16,4f15.6)') "Stretch PE     :", &
            storestats_getintranrg(imodel%spcstats, i,'inst','stretch')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'cavg','stretch')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'block','stretch')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'std','stretch')/nmoles
      End If
      If (interact_hasint(imodel,i,'constraint')) Then
        Write(dUnit,'(3x,a16,4f15.6)') "Constraint PE  :", &
            storestats_getintranrg(imodel%spcstats, i,'inst','constraint')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'cavg','constraint')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'block','constraint')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'std','constraint')/nmoles
      End If
      If (interact_hasint(imodel,i,'bending')) Then
        Write(dUnit,'(3x,a16,4f15.6)') "Bending PE     :", &
            storestats_getintranrg(imodel%spcstats, i,'inst','bending')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'cavg','bending')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'block','bending')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'std','bending')/nmoles
      End If
      If (interact_hasint(imodel,i,'torsion')) Then
        Write(dUnit,'(3x,a16,4f15.6)') "Torsional PE    :", &
            storestats_getintranrg(imodel%spcstats, i,'inst','torsion')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'cavg','torsion')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'block','torsion')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'std','torsion')/nmoles
      End If
      If (interact_hasint(imodel,i,'intrapair')) Then
        Write(dUnit,'(3x,a16,4f15.6)') "IntramolecPair :", &
            storestats_getintranrg(imodel%spcstats, i,'inst','intrapair')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'cavg','intrapair')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'block','intrapair')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'std','intrapair')/nmoles
      End If
      If (interact_hasint(imodel,i,'intracoul')) Then
        Write(dUnit,'(3x,a16,4f15.6)') "IntramolecCoul :", &
            storestats_getintranrg(imodel%spcstats, i,'inst','intracoul')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'cavg','intracoul')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'block','intracoul')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'std','intracoul')/nmoles
      End If
      pe = storestats_getintranrg(imodel%spcstats, i,'inst')/nmoles + pe
      peavg = storestats_getintranrg(imodel%spcstats, i,'cavg')/nmoles + peavg
      ke = storestats_getke(imodel%spcstats, i,'inst')/nmoles + ke
      keavg = storestats_getke(imodel%spcstats, i,'cavg')/nmoles + keavg
      moles = moles + nmoles
    End Do

    Do i = 1,nsorbs
      nmoles = config_getnmoles(sorbates,i)
      Do j = i, nsorbs
        If (i /= j) nmoles = nmoles + config_getnmoles(sorbates,j)
        If (config_isfixed(sorbates(i)).And.config_isfixed(sorbates(j))) Cycle
        outString = molecules_name(i)
        outString2= molecules_name(j)
        Write(dUnit,'(1x,4a)') "Energy for ",Trim(outString),"-", &
            Trim(outString2)
        Write(dUnit,'(3x,a16,4f15.6)') "Noncoulombic   :", &
            storestats_getnoncoul(imodel%spcstats, i,j,'inst')/nmoles, &
            storestats_getnoncoul(imodel%spcstats, i,j,'cavg')/nmoles, &
            storestats_getnoncoul(imodel%spcstats, i,j,'block')/nmoles, &
            storestats_getnoncoul(imodel%spcstats, i,j,'std')/nmoles
        Write(dUnit,'(3x,a16,4f15.6)') "Coulombic      :",&
            storestats_getcoul(imodel%spcstats, i,j,'inst')/nmoles, &
            storestats_getcoul(imodel%spcstats, i,j,'cavg')/nmoles, &
            storestats_getcoul(imodel%spcstats, i,j,'block')/nmoles, &
            storestats_getcoul(imodel%spcstats, i,j,'std')/nmoles
        pe = storestats_getnoncoul(imodel%spcstats, i,j,'inst')/nmoles + &
            storestats_getcoul(imodel%spcstats, i,j,'inst')/nmoles + pe
        peavg = storestats_getnoncoul(imodel%spcstats, i,j,'cavg')/nmoles + &
            storestats_getcoul(imodel%spcstats, i,j,'inst')/nmoles + peavg
      End Do
    End Do

    !** Get any extra energies from the integration routines
    !** e.g., Nose-Hoover energy, or others
    pei = 0.0_Rdbl
    kei = 0.0_RDbl
    peavgi = 0.0_RDbl
    keavgi = 0.0_Rdbl

    If (Associated(params%vverlet)) Then
      Call velverlet_extraEnergy(params%vverlet,pei,peavgi,'pe')
      Call velverlet_extraEnergy(params%vverlet,kei,keavgi,'ke')
    Else If (Associated(params%gear6)) Then
      Call gear6_extraEnergy(params%gear6,pei,peavgi,'pe')
      Call gear6_extraEnergy(params%gear6,kei,keavgi,'ke')
    Else If (Associated(params%lfrog)) Then
      Call leapfrog_extraEnergy(params%lfrog,pei,peavgi,'pe')
      Call leapfrog_extraEnergy(params%lfrog,kei,keavgi,'ke')
    End If

    If ((pei /= 0.0).And.(kei /= 0.0)) Then
      Write(dUnit,'(1x,a)') "Extra Integrator Energies (NH, etc)"
      Write(dUnit,'(3x,a16,2f15.6)') "Potential      :", &
          pei/moles,peavgi/moles
      Write(dUnit,'(3x,a16,2f15.6)') "Kinetic        :", &
          kei/moles,keavgi/moles
    End If

    pe = pe + pei/moles
    ke = ke + kei/moles
    peavg = peavg + peavgi/moles
    keavg = keavg + keavgi/moles

    !** Display the dynamic information from the imodel if necessary
    Call interact_dispDynInfo(imodel,1,dUnit)

    Write(dUnit,'(1x,a)') "System Energy Totals"
    Write(dUnit,'(3x,a16,2f15.6)') "Total Potential:",pe,peavg
    Write(dUnit,'(3x,a16,3f15.6)') "Total Kinetic  :",ke,keavg
    Write(dUnit,'(3x,a16,2f15.6)') "Total Energy   :",pe+ke,peavg+keavg
    Write(dUnit,'(3x,a16,f15.6)')  "COM Velocity   :", &
        integrate_comVelocity(1,sorbates(1)%coords(:,1)%v)

    Write(dUnit,'(1x,a)') "----------------------------------------------------"

  End Subroutine integrate_simdisplay

  !----------------------------------------------------------------------------
  ! Calls interact to display any extra energy statistics
  !----------------------------------------------------------------------------
  Subroutine integrate_simdisplayExtra(params,imodel,dUnit,indent)
    Type(Integrate_Params), Intent(In)         :: params
    Type(Interaction_Model), Intent(InOut)     :: imodel
    Integer, Intent(In)                        :: dUnit,indent

    Call interact_dispDynInfo(imodel,indent,dUnit)

  End Subroutine integrate_simdisplayExtra

  !----------------------------------------------------------------------------
  ! clean
  !----------------------------------------------------------------------------
  Subroutine integrate_cleanup

  End Subroutine integrate_cleanup

End Module integrate

