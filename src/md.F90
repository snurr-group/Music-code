!------------------------------------------------------------------------------
! Molecular Dynamics module
!------------------------------------------------------------------------------

Module md

  Use atom, Only: atom_getmass, atom_getname, atom_getsymbol
  Use auxmoveparams, Only: AuxMoveObjects
  Use config, Only: AtMolCoords, config_getnmoles, config_getnatoms, &
      config_getconfig, config_isfixed, config_kineticEnergy, &
      config_kineticTemp, config_getSpcCOMVel, config_getCOMVel
  Use datafile, Only: CONFILE
  Use defaults, Only: RDbl, strLen, lstrLen, kjmole_kb, scaleke, scalepe, &
      zero, dashedline, TUnit, velUnit, nrgUnit, MAX_SORBS
  Use emd, Only: emd_init, emd_display, emd_output, EMD_IS_ON
  Use file, Only: file_getunit, file_open
  Use gear6, Only: gear6_extraenergy
  Use interact, Only: Interaction_Model,interact_hasint,interact_initnrgs, &
      interact_totnrg
  Use subinteract, Only: Subset_Interactions, subinteract_simpleupdate, &
      subinteract_init, subinteract_int
  Use intramolecular, Only: intramolecular_getpenalty
  Use molecule, Only: MolecularParams
  Use molecules, Only: molecules_getnsorbs, molecules_name, &
      molecules_getnatoms, molecules_getnthatom, molecules_getnatomtypes, &
      molecules_getmass, molecules_getdof, molecules_dofIsOK, &
      molecules_isCoreShell
  Use moves,Only: Move_Params, moves_displayparams, moves_initparams, &
      moves_integrate, moves_sampleCF, moves_simdisplayExtra, moves_getparam
  Use nemd, Only: nemd_init, nemd_display, NEMD_IS_ON,nemd_simdisplay
  Use nbrlist, Only: Nbrlist_Params  , nbrlist_update

  Use random, Only: random_gaussian
  Use simcell, Only: SimCell_Params, simcell_pbc, simcell_getfillsorb
  Use stats, Only: stats_update
  Use storestats, Only: storestats_getintranrg, storestats_getke, &
      storestats_getnoncoul, storestats_gettemp, &
      storestats_updatetemp, storestats_updatenrg, &
      storestats_updateenergySS, storestats_getcoul, &
      storestats_updateIntraEnergy
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, &
      tolower, toupper, allocerrdisplay, int2str, real2str, checkandstop
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use velverlet, Only: velverlet_extraenergy
  Use fixedsorbate, Only: fixedsorbate_init
      
  Implicit None
  Save

  Private
  Public :: MDInfo, md_initParams, md_initEquil, md_dosim, md_display, &
      md_simdisplay, md_simsummary, md_sampleCF, md_getdt, md_displayinit

  Type MDInfo
    Real(kind=RDbl)        :: time
    Logical                :: sameSpeed,genVel
    Logical                :: movesUpdateFlag !** true-> do updates in move
    Character(len=lstrLen) :: velFile
    Integer                :: no_of_movetypes
    Type(Move_Params)                        :: eqParams !** for equilibration
    Type(Move_Params), Dimension(:), Pointer :: mvParams  !** for prodn run
    Type(Subset_Interactions) :: subint  !** points to energy 
                                       !** storage, interaction model etc 
    Type(AuxMoveObjects) :: auxmv
    Logical :: nbrlist_ON
    Type(Nbrlist_params),Pointer                      :: nbrlist
  End Type MDInfo

  Character(len=6), Parameter :: eunits = 'kJ/mol'
  Character(len=1), Parameter :: tunits = 'K'
  Character(len=strLen), Parameter :: md_tag = "MD Information"
  Character(len=lstrLen), Parameter :: mdequil_tag = &
      "MD Equilibration Information"

Contains

  !----------------------------------------------------------------------------
  ! Initializes the MD parameters from the control file. It reads in the
  ! information needed to run an MD simulation. This is what you want if
  ! you only want to read in the information from the md file.
  !----------------------------------------------------------------------------
  Subroutine md_initParams(mdParams,species,scell,ctrl_filename,tag)
    Type(MDInfo), Intent(InOut) :: mdParams
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Character(*), Intent(In) :: ctrl_filename
    Type(SimCell_Params), Intent(In)               :: scell
    Character(*), Intent(In) :: tag
    Character(len=strLen) :: thermostat, line, dof_origin

    Integer :: unitno, useMD, j, error, i, nsorbs, spc, dof
    Real(kind=RDbl) :: Tk
    Character(len=255) :: text
    Character(len=strLen), Dimension(strLen) :: params

    !** Check and initialize NEMD if necessary, to be done before moves
    Call nemd_init(ctrl_filename,species)


    !** Open the control file if it is not already open
    unitno = file_open(ctrl_filename)

    !** Find the MD section
    useMD = filesrchstr(unitno,tag,text,.True.)
    If (useMD /= 0) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      !** Read in the number of move types listed in the MD section
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      j = split(text,params)
      mdParams%no_of_movetypes = toint(params(1))

      !** Read in the location of the initial positions and velocities.
      !** This is either another MD run or "GENERATE", in which case
      !** we will need to run an equilibration simulation
      Read(unitno,'(a)') text
      text = stripcmnt(text)
      j = split(text,params)

      !** Find out if we generate an equilibrium and velocities
      !** or read them in from a file
      If (tolower(params(1)) == "generate") Then
        mdParams%genVel = .True.
      Else
        mdParams%genVel = .False.
        mdParams%velFile = params(1)
      End If


      !** Allocate space for the different move types listed in the
      !** control file.
      Allocate(mdParams%mvParams(mdParams%no_of_movetypes),stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            "Could not allocate memory for 'mdParams%mvparams'"
        Stop
      End If

      !** Read in all the MD move types
      Do i = 1, mdParams%no_of_movetypes
        !Read the blank line separating the move types
        Read(unitno,'(a)') line
        Call moves_initparams(mdParams%mvparams(i), mdParams%auxmv, &
            species,ctrl_filename)
      End Do

      !** this part coded only for one movetype
      If (mdParams%no_of_movetypes > 1) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, " need more code here "
        Stop
      End If

      thermostat=trim(toupper(md_getthermostat(mdParams%mvparams(1))))
      Write(*,*) "thermostat used is : ",thermostat 

      If (thermostat=="VELOCITYRESCALE")  Then
        mdParams%movesUpdateFlag = .True.
      Else
        mdParams%movesUpdateFlag = .False.
      Endif

      !** Get the number of molecule types
      nsorbs = molecules_getnsorbs()

    End If

    !** All MD calculations should be able to calculate their degrees of freedom --
    !** Other wise this check can be removed
    Do spc=1,nsorbs
      If (config_isfixed(species(spc))) Cycle
      dof = molecules_getdof(spc,dof_origin)
      If (molecules_dofIsOK(spc)) Then
        If (Trim(Toupper(dof_origin)) == 'MOLECULE_FILE') Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "MD is equipped with calculating dof, It might be erraneous to "
          Write(*,*) "specify it again in molecule file, for species : ",trim(molecules_name(spc))
          Write(*,*) " either recompile after commenting the following STOP"
          Write(*,*) " or remove DOF line from corresponding molecule file"
          Stop              !***** STOP ********
        End If
      Else
        Write(*,*) "wrong dof for molecule"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    End Do

    !** Check and initialize EMD for L_ij if necessary
    ! to be done after md ini
    Tk=mdParams%mvParams(1)%integrate%simT
    Call emd_init(ctrl_filename,species,scell,Tk)

    !** we dont want to run both emd and nemd simulatneously by mistake
    If ((NEMD_IS_ON).And.(EMD_IS_ON)) Then
      Write(*,*) "Decide whether you want nemd or emd and then modify ctrlfile"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Call fixedsorbate_init(ctrl_filename)

  End Subroutine md_initParams

  !----------------------------------------------------------------------------
  ! Initialize the MD parameters from the control file, and also runs
  ! the initial equilibration if necessary. This is probably what you want
  ! if you are going to be doing MD. It calls md_initParams to read stuff
  ! in from the control file.
  !----------------------------------------------------------------------------
  Subroutine md_initEquil(mdParams, imodel, scell, species, &
      ctrl_filename,tag, conf_file,unit)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(MDInfo), Intent(InOut)                    :: mdParams
    Type(SimCell_Params), Intent(In)               :: scell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Character(*), Intent(In)                       :: ctrl_filename
    Type(CONFILE), intent(In)                      :: conf_file
    Character(*), Intent(In)                       :: tag
    Integer, Intent(In), Optional                  :: unit

    Integer                :: unitno, useMD, j, i, nmoles, natoms, dUnit,nsorbs
    Logical                :: fast, recalc_flag, calc_accels, success, update
    Real(kind=RDbl)        :: objf,sum, dt, TKelvin
    Character(len=255)     :: text

    !** Set the default display unit
    dUnit = 6
    If (Present(unit)) dUnit = unit


    !** Set the genVel flag to false.
    mdParams%genVel = .False.

    !** Call md_initParams to read in the parameters from the control file.
    Call md_initParams(mdParams,species,scell,ctrl_filename,tag)



    !** initialize the subset interaction structure
    Call subinteract_init(mdParams%subint,imodel,'Molec_System', &
        'MD',(/0,0,0/),(/0,0,0/))

    !** check for neighborlist
    mdParams%nbrlist_ON=.False.
    If (Associated(imodel%ff(1)%nlist)) Then
      mdParams%nbrlist_ON=.True.
      mdParams%nbrlist=> imodel%ff(1)%nlist
      Call nbrlist_update(mdParams%nbrlist,species,scell)
    Endif



    !** Update all the stored interactions like accelerations before 
    !** start of integartion
    !** This is very important for comparing MD vs some other simulation
    !** Otherwise the first step of velocity verlet/Gear6 will give 
    !** wrong results
    !** don't update spc stats to maintain consistency with previous versions
    fast = .True.
    update = .False.  
    Call interact_initnrgs(imodel,fast,species,scell,update)

    !** Get the unitno for the control file
    unitno = file_open(ctrl_filename)

    !** Read the equilibration information if necessary
    If (mdParams%genVel) Then
      !** Search for the equilibration info
      useMD = filesrchstr(unitno,mdequil_tag,text,.True.)
      If (useMD == 0) Then
        Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
            " Could not find tag '",Trim(mdequil_tag),"' in file ", &
            ctrl_filename
        Stop
      End If

      Call moves_initparams(mdParams%eqParams, mdParams%auxmv, &
          species,ctrl_filename)
      Call moves_displayparams(mdParams%eqParams,9,dUnit)
    End If

    !** Generate the velocities if necessary
    If (mdParams%genVel) Then
      Write(6,'(1x)')
      Write(6,'(1x,a)') "Generating Velocities"
      Call md_genVel(imodel,scell,mdParams%eqParams,species)
    End If

    !** Check the constraints, if required
    !** This is HACK, should not reach directly into intramolecular
    nsorbs = molecules_getnsorbs()
    If (imodel%nff > 1) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          " More than one forcefield, cannot check constraints during MD init"
      Stop
    End If
    Do i = 1,nsorbs
      If (.Not. interact_hasint(imodel,i,'constraint')) Cycle
      natoms = config_getnatoms(species,i)
      nmoles = config_getnmoles(species,i)
      sum = 0.0_RDbl

      Do j = 1, nmoles
        Call intramolecular_getpenalty(imodel%ff(1)%iparams(i), &
            species(i)%coords(1:natoms,j)%rp, &
            species(i)%coords(1:natoms,j)%v,objf)
        sum = objf + sum
      End Do
      !MDEBUG
      If (i == 1) Then
        Write(*,*) "objf for ",Trim(molecules_name(i)),sum
      End If
      If (Abs(sum) > 1.0_RDbl) Then
        Write(0,'(4x,4a,f15.2)') "Sum of objf in penalty exceeds limit", &
            "for ",Trim(molecules_name(i))," : ",sum
        Stop
      End If
    End Do

    !** If we generated new velocities, we need to run 
    !** some MD steps to equilibrate
    If (mdParams%genVel) Then
      !Generate the equilibration
      Write(dUnit,'(1x)')
      Write(dUnit,'(1x,a)') "Running Equilibration"
      Call md_equilibrate(imodel,mdParams,scell,species,conf_file)
    End If

    !** Zero the time value
    mdParams%time = 0.0_Rdbl

  End Subroutine md_initEquil

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine md_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a)') md_tag
    Write(unitno,'(a,t30,a)') 'Integer','# Number of MD move types listed'
    Write(unitno,'(a,t30,a)') '[Generate,filename]', &
        '# Initial velocities, generate or read from filename'
    Write(unitno,'(t30,a)')  &
        '# Blank line (required), then for each MD move type'
    Call moves_sampleCF(unitno,'INTEGRATE')
    Write(unitno,'(a)') &
        '# If GENERATE velocities, remove this line and include:'
    Write(unitno,'(a)') dashedline
    Write(unitno,'(a)') trim(mdequil_tag)
    Call moves_sampleCF(unitno,'INTEGRATE')
    Write(unitno,'(a)') dashedline
    Write(unitno,'(a)') &
        '# End GENERATE section (remove this line)'

  End Subroutine md_sampleCF

  !----------------------------------------------------------------------------
  ! Sets the MD temperature according to the desired temperature from 
  ! a gaussian distribution. It correctly handles core/shell systems 
  ! as well.
  !----------------------------------------------------------------------------
  Subroutine md_genVel(imodel,simcell,mvParams,species,unit)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(Move_Params), Intent(In)                  :: mvParams
    Type(AtMolCoords), Intent(InOut), Dimension(:) :: species
    Integer, Intent(In), Optional                  :: unit

    Integer :: natoms, nmoles, i, j, k, dUnit, ntypes, error
    Logical :: isCoreShell  ! we set to TRUE if atom is part of core/shell pair
    Logical :: isCore       ! we set to TRUE if atom is a core
    Integer :: shellAtom    ! atom number of the shell corresponding to a core
    Integer, Allocatable, Dimension(:) :: atypes         ! array of atom TYPES
    Real(kind=RDbl)      :: ave, sigma, tg, beta, tmole, desired_T
    Real(Kind=RDbl)      :: comVelocity, tatom, imoles, vtemp, scale
    Type(VecType)        :: sum

    !** Assign the default display unit
    If (.Not.Present(unit)) Then
      dUnit = 6
    Else
      dUnit = unit
    End If

    !** Loop through each of the species and generate their velocities
    Do k = 1, Size(species)

      !** Update the user about our progress
      Write(dUnit,'(a,i5,2a)') __FILE__,__LINE__, &
          " Generating velocities for ",Trim(molecules_name(k))

      !** Is the molecule fixed? If so, skip it.
      If (config_isfixed(species(k))) Cycle

      !** get the number of molecules of this sorb type,
      !** and the number of atoms
      natoms = config_getnatoms(species,k)
      nmoles = config_getnmoles(species,k)
      If (nmoles==0) Cycle

      !** Set the desired temperature for this sorbate
      desired_T= mvParams%integrate%simT
      tg = desired_T*molecules_getdof(k)/3.0_RDbl

      !** Calculate the sigma parameter (distribution width) for
      !** the Gaussian velocity distribution 
      beta = 1.0_RDbl / (kjmole_kb/scaleke*tg)
      ave = 0.0_RDbl
      sigma = 1.0_RDbl / Sqrt(beta*molecules_getmass(k))
      sum = 0.0_RDbl

      !** Size the atom types array
      error = 0
      If (.Not.Allocated(atypes)) Then
        Allocate(atypes(molecules_getnatoms(k)),stat=error)
      Else
        If (Size(atypes,1) < molecules_getnatoms(k)) Then
          Deallocate(atypes,stat=error)
          If (error == 0) Then
            Allocate(atypes(molecules_getnatoms(k)),stat=error)
          Else
            Write(0,'(a,i5,a)') __FILE__,__LINE__, &
                ": Error deallocating atypes array"
            Stop
          End If
        End If
      End If
      !** Report allocation errors, if any
      If (error /= 0) Then
        Write(0,'(a,i5,a,i12)') __FILE__,__LINE__, &
            ": Error allocating atypes array of size ",molecules_getnatoms(k)
        Stop
      End If

      !** Get the atom types for this hack
      ntypes = molecules_getnatomtypes(k,atypes,.True.)

      If ((nmoles == 1).And.(simcell_getfillsorb(simcell) == k)) Then
        Write(0,'(a,i5,2a)') __FILE__,__LINE__,": ", &
            "Running genVel for the simulation cell fill sorb"
        !** Hack for siliflex
        !** Loop through and find the Silicon and Oxygen cores, assigning
        !** random velocites to them. For the shells, assign the same
        !** velocity as that assigned to the core.
        !** nmoles = 1
        j = 1

        !** Loop through all the atoms of the molecule
        Do i = 1, natoms

          !** If this is a core-shell pair, check if it is a core.
          !** If so, assign the core-shell a COM velocity. If it is
          !** a shell, cycle.
          isCoreShell = molecules_isCoreShell(k,i,isCore,shellAtom) 
          If (isCoreShell.And..Not.isCore) Cycle

          !** Assign a velocity to the atom
          species(k)%coords(i,j)%v = (/ random_gaussian(ave,sigma), &
              random_gaussian(ave,sigma), random_gaussian(ave,sigma) /)

          !** Add the velocity to the center-of-mass velocity sum
          sum = sum + species(k)%coords(i,j)%v*atom_getmass(atypes(i))

          !** If this is a core, assign the shell the same velocity
          If (isCore) Then 
            species(k)%coords(shellAtom,j)%v = species(k)%coords(i,j)%v
            sum = sum+species(k)%coords(i,j)%v*atom_getmass(atypes(shellAtom))
          End If
        End Do

        !** We are using the COM for the velocity distribution.
        imoles = 1.0_RDbl/molecules_getmass(k)
        vtemp = 0.0_RDbl

        !** Get the average velocity per atom
        sum = sum*imoles

        !** Hack for siliflex - subtract out the COM velocity
        Do i = 1, nmoles
          Do j = 1, natoms
            species(k)%coords(j,i)%v = species(k)%coords(j,i)%v - sum
            vtemp = vtemp + species(k)%coords(j,i)%v*species(k)%coords(j,i)%v
          End Do
        End Do

      Else
        !** Assign a single velocity to the entire molecule from a random
        !** Gaussian distribution
        If(nmoles>1) then
          !** EXPLAIN WHAT'S GOING ON
          Write(*,*) "all atoms of one molecule will be given same velocity"
          Write(*,*) "velocity of the center of mass of this species will be brought to zero"

          Do i = 1, nmoles
            species(k)%coords(1,i)%v = (/ random_gaussian(ave,sigma), &
                random_gaussian(ave,sigma), random_gaussian(ave,sigma) /)
            Do j = 2, natoms
              species(k)%coords(j,i)%v = species(k)%coords(1,i)%v 
            End Do

            sum = sum + species(k)%coords(1,i)%v
          End Do
        ElseIf(nmoles==1) then
          Write(*,*) "Each atoms of this one molecule will be given different velocity"
          Write(*,*) "velocity of the center of mass of this molecule will be brought to zero"
          Do j = 1, natoms
            species(k)%coords(j,1)%v = (/ random_gaussian(ave,sigma), &
                random_gaussian(ave,sigma), random_gaussian(ave,sigma) /)

          End Do
          sum = sum + config_getcomVel(species,k,1)
        End If

        imoles = 1.0_RDbl/Real(nmoles,kind=RDbl)
        sum = sum*imoles
        vtemp = 0.0_RDbl

        !        If (nmoles /= 1) Then
        Do i = 1, nmoles
          Do j = 1, natoms
            species(k)%coords(j,i)%v = species(k)%coords(j,i)%v - sum
          End Do
          vtemp = vtemp + species(k)%coords(1,i)%v*species(k)%coords(1,i)%v
        End Do
        !        End If
      End If

      !** Report RMS velocity for the molecule type
      vtemp = Sqrt(vtemp*imoles)
      Write(dUnit,'(a,f10.3)') 'Generated RMS velocity : ',vtemp

      !** Calculate the atomic and molecular temperatures
      tmole = config_kineticTemp(species,k,'molecular')
      tatom = config_kineticTemp(species,k,'atomic')

      Call storestats_updateTemp(imodel%spcstats, k, tmole, 'molecular')
      Call storestats_updateTemp(imodel%spcstats, k, tatom, 'atomic')

      comVelocity = md_comVelocity(species,k)

      !** Report the temperatures
      Write(dUnit,'(1x,2a)') 'Generated velocities for ',molecules_name(k)
      Write(dUnit,'(1x,a,f12.3,1x,a)') 'Molecular temperature : ',tmole, &
          trim(TUnit)
      Write(dUnit,'(1x,a,f12.3,1x,a)') 'Atomic temperature    : ',tatom, &
          trim(TUnit)
      Write(dUnit,'(2x,a,f12.3,1x,a)') '        COM velocity : ',comVelocity, &
          trim(velUnit)

      !** Scale parameter
      scale = (1.0_RDbl + (desired_T/tmole - 1.0_RDbl)*1.0_RDbl/4.0_RDbl) &
          **(0.5_RDbl)
      !      scale = (1.0_RDbl + (desired_T/tatom - 1.0_RDbl)*1.0_RDbl/4.0_RDbl) &
      !          **(0.5_RDbl)
      !MDEBUG
      Write(0,*) __FILE__,__LINE__," scale : ",scale

      Do i = 1, nmoles
        tg = desired_T*3.0_RDbl/3.0_RDbl
        !** Rescale the velocities based on the desired temperature
        Call md_rescaleVel2(tg,species(k)%coords(1:natoms,i)%v,tatom)
      End Do

      !** Calculate the new temperatures after rescaling
      tmole = config_kineticTemp(species,k,'molecular')
      tatom = config_kineticTemp(species,k,'atomic')

      Call storestats_updatetemp(imodel%spcstats, k, tmole, 'molecular')
      Call storestats_updatetemp(imodel%spcstats, k, tatom, 'atomic')

      !** Calculate the com velocity
      comVelocity = md_comVelocity(species,k)

      !** Report the new rescaled temperatures
      Write(dUnit,'(1x,a,f12.3,1x,a)') 'Rescaled molecular T  : ',tmole, &
          trim(TUnit)
      Write(dUnit,'(1x,a,f12.3,1x,a)') 'Rescaled atomic T     : ',tatom, &
          trim(TUnit)
      Write(dUnit,'(1x,a,f12.3,1x,a)') 'Rescaled COM velocity : ', &
          comVelocity, trim(velUnit)
      Write(dUnit,'(1x,a,f12.3,1x,a)') 'Kinetic Energy        : ', &
          config_kineticEnergy(species,k)/nmoles, Trim(nrgUnit) 

    End Do

  End Subroutine md_genVel

  !----------------------------------------------------------------------------
  ! Do the MD simulation, i.e., make a move, dump to io, repeat
  !----------------------------------------------------------------------------
  Subroutine md_doSim(mdParams,moveno,simcell,species,conf_file,time)
    Type(MDInfo), Intent(InOut)                    :: mdParams
    Integer,Intent(In)                             :: moveno
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(CONFILE), Intent(In)                      :: conf_file
    Real(kind=RDbl), Intent(InOut), Optional       :: time

    Integer                     :: i,j,nsorbs,natoms,nmoles
    Real(kind=RDbl)             :: tmole,tatom,ke,pot
    Type(AtMolCoords), Pointer  :: sorb

    !** Update nighborlist if necessary
    If (mdParams%nbrlist_ON) Then
      Call nbrlist_update(mdParams%nbrlist,species,simcell,time)
    End If

    !** Make the move(s)
    Call moves_integrate(mdParams%subint, mdParams%mvparams(moveno), simcell, &
        species, conf_file, mdParams%movesUpdateFlag, time)

    !** Do temperature and energy updates here
    If (.Not. mdParams%movesUpdateFlag) Then
      Call subinteract_simpleupdate(mdParams%subint,species)
    End If

    !** If you aim to calculate Onsager coefiicients then output velocities
    If (EMD_IS_ON) Then
      Call emd_output(species,time)
    Endif

  End Subroutine md_doSim

  !----------------------------------------------------------------------------
  ! Do the MD equilibration
  !----------------------------------------------------------------------------
  Subroutine md_equilibrate(imodel,mdParams,simcell,species,conf_file)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(MDInfo), Intent(InOut)                    :: mdParams
    Type(SimCell_Params), Intent(In)               :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type (CONFILE),intent(in)                      :: conf_file    

    Logical         :: updateFlag
    Real(kind=RDbl) :: time,pot

    updateFlag = .True.
    time = zero

    Call moves_integrate(mdParams%subint, mdParams%eqparams,simcell,species, &
        conf_file,.True.,time)

  End Subroutine md_equilibrate
!!$ Check below with rescale, generate options then delete this
!!$
!!$  !----------------------------------------------------------------------------
!!$  ! Calculate the center of mass velocity of the molecule
!!$  !----------------------------------------------------------------------------
!!$  Real(Kind=RDbl) Function md_comVelocity(molecule,vel)
!!$    Integer, Intent(In)                     :: molecule
!!$    Type(VecType), Dimension(:), Intent(In) :: vel
!!$
!!$    Integer :: natoms, natypes, i 
!!$    Integer, Dimension(Size(vel,1)) :: atypes
!!$    Real(kind=RDbl) :: mass, comvel
!!$    Type(VecType) :: comv
!!$
!!$    !** Zero the velocity and mass
!!$    comvel = 0.0_RDbl
!!$    mass = 0.0_RDbl
!!$    comv = 0.0_RDbl
!!$
!!$    !** Get the number of atoms
!!$    natoms = molecules_getnatoms(molecule)
!!$
!!$    !** Get the atom types of the molecule
!!$    natypes = molecules_getnatomtypes(molecule,atypes,.True.)
!!$
!!$    !** Loop through the atoms and get the center of mass velocity
!!$    Do i = 1, natoms
!!$      comv = comv + vel(i)*atom_getmass(atypes(i))
!!$      mass = mass + atom_getmass(atypes(i))
!!$    End Do
!!$
!!$    comv = comv/mass
!!$    comvel = mag(comv)
!!$
!!$    md_comVelocity = comvel
!!$
!!$  End Function md_comVelocity
!!$
!!$  !----------------------------------------------------------------------------
!!$  ! Calculate the center of mass velocity of the molecule Returns Vector !
!!$  !----------------------------------------------------------------------------
!!$  Type(VecType) Function md_comvVelocity(molecule,vel)
!!$    Integer, Intent(In)                     :: molecule
!!$    Type(VecType), Dimension(:), Intent(In) :: vel
!!$
!!$    Integer :: natoms, natypes, i 
!!$    Integer, Dimension(Size(vel,1)) :: atypes
!!$    Real(kind=RDbl) :: mass, comvel
!!$    Type(VecType) :: comv
!!$
!!$    !** Zero the velocity and mass
!!$    comvel = 0.0_RDbl
!!$    mass = 0.0_RDbl
!!$    comv = 0.0_RDbl
!!$
!!$    !** Get the number of atoms
!!$    natoms = molecules_getnatoms(molecule)
!!$
!!$    !** Get the atom types of the molecule
!!$    natypes = molecules_getnatomtypes(molecule,atypes,.True.)
!!$
!!$    !** Loop through the atoms and get the center of mass velocity
!!$    Do i = 1, natoms
!!$      comv = comv + vel(i)*atom_getmass(atypes(i))
!!$      mass = mass + atom_getmass(atypes(i))
!!$    End Do
!!$
!!$    comv = comv/mass
!!$    comvel = mag(comv)
!!$
!!$    md_comvVelocity = comv
!!$
!!$  End Function md_comvVelocity
!!$

  !----------------------------------------------------------------------------
  ! Calculate the magnitude of center of mass velocity of the molecule
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function md_comVelocity(species,spc)
    Integer, Intent(In)                     :: spc
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    md_comVelocity = mag(config_getSpcCOMVel(species,spc))
  End Function md_comVelocity


  !----------------------------------------------------------------------------
  ! Calculate the center of mass velocity of the molecule
  !----------------------------------------------------------------------------
  Type(VecType) Function md_comvVelocity(species,spc)
    Integer, Intent(In)                     :: spc
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    md_comvVelocity = config_getSpcCOMVel(species,spc)
  End Function md_comvVelocity

  !----------------------------------------------------------------------------
  ! Rescale the velocities according to the temperature
  !----------------------------------------------------------------------------
  Subroutine md_rescaleVel(simT,vel,temp)
    Type(VecType), Dimension(:,:), Intent(InOut) :: vel
    Real(kind=RDbl), Intent(In) :: temp, simT
    Real(kind=RDbl) :: scale
    Integer :: natoms, nmoles, i,j

    scale  = Sqrt(simT/temp)
    natoms = Size(vel,1)
    nmoles = Size(vel,2)

    Do i = 1, nmoles
      Do j = 1, natoms
        vel(j,i) = vel(j,i)*scale
      End Do
    End Do
  End Subroutine md_rescaleVel

  !----------------------------------------------------------------------------
  ! Rescale the velocities according to the temperature
  !----------------------------------------------------------------------------
  Subroutine md_rescaleVel2(simT,vel,temp)
    Type(VecType), Dimension(:), Intent(InOut) :: vel
    Real(kind=RDbl), Intent(In) :: temp, simT
    Real(kind=RDbl) :: scale
    Integer :: natoms, j

    scale  = Sqrt(simT/temp)
    natoms = Size(vel,1)

    Do j = 1, natoms
      vel(j) = vel(j)*scale
    End Do
  End Subroutine md_rescaleVel2

  !----------------------------------------------------------------------------
  ! Display information about the MD simulation
  !----------------------------------------------------------------------------
  Subroutine md_display(mdParams,unitno)
    Type(MDInfo), Intent(In) :: mdParams
    Integer, Intent(In), Optional :: unitno

    Integer :: dUnit  ! display unit number
    Integer :: i

    dUnit = 6
    If (Present(unitno))       dUnit = unitno

    Write(dUnit,'(a)') '-----------------------------------------------'
    Write(dUnit,'(a)') ' MD Simulation Information'
    Write(dUnit,'(a)') '-----------------------------------------------'
    Write(dUnit,'(4x,a,i6)') 'number of movetypes          : ',&
        mdParams%no_of_movetypes
    Write(dUnit,'(4x,a)') ' Different movetypes during prodn run: '
    Write(dUnit,'(4x,a)') '-----------------------------------------------'

    Do i = 1,mdParams%no_of_movetypes
      !** Indent by 9 spaces
      Call moves_displayparams(mdParams%mvparams(i),9,dUnit)
    End Do

    Write(dUnit,'(4x,a)') ' Equilibration details : '
    Write(dUnit,'(4x,a)') '-----------------------------------------------'
    If (mdParams%genVel) Then
      !** Indent by 9 spaces
      Call moves_displayparams(mdParams%eqParams,9,dUnit)
    Else
      Write(dUnit,'(4x,2a)') ' Velocities Taken from file : ',&
          Trim(mdParams%velFile)
    Endif
    If (NEMD_IS_ON) Then
      Write(dUnit,'(4x,a)') '-----------------------------------------------'
      Call nemd_display(dUnit, 6)
    Endif
    If (EMD_IS_ON) Then
      Write(dUnit,'(4x,a)') '-----------------------------------------------'
      Call emd_display(dUnit, 6)
    Endif
    Write(dUnit,*)
  End Subroutine md_display

  !----------------------------------------------------------------------------
  ! Displays some initial information about the MD system, specifically
  ! the temperatures, kinetic energy, and COM velocity. A useful check
  ! when starting from an equilibrated file, I think.
  !----------------------------------------------------------------------------
  Subroutine md_displayinit(species,dUnit,nspc)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In) :: dUnit
    Integer, Intent(In) :: nspc
    Character(nspc) :: spaces
    Integer :: i
    Real(Kind=RDbl) :: tmole, tatom
    Type(VecType) :: comVelocity
    spaces = ""
    Do i = 1, nspc
      spaces = spaces//" "
    End Do

    !** Report just a bit of information, like the temperature and the
    !** kinetic energy, both of which we can calculate here
    Do i = 1, molecules_getnsorbs()
      !** Temperatures
!      tmole = md_kineticTemp(i,species(i)%coords(:,:)%v,'molecular')
!      tatom = md_kineticTemp(i,species(i)%coords(:,:)%v,'atomic')

!      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      tmole = config_kineticTemp(species,i,'molecular')
      tatom = config_kineticTemp(species,i,'atomic')

      !** Calculate the com velocity (for first molecule)
      comVelocity = config_getCOMVel(species,i,1)

      !** Report the information
      Write(dUnit,'(2a)') spaces,dashedline
      Write(dUnit,'(2a)') spaces,"MD Start-up Kinetic Information"
      Write(dUnit,'(a)') spaces
      Write(dUnit,'(2a,i3,2a)') spaces,"Species ",i,": ",molecules_name(i)
      Write(dUnit,'(2x,2a,f10.3,1x,a)') spaces,"T(atom)  = ",tatom,Trim(TUnit)
      Write(dUnit,'(2x,2a,f10.3,1x,a)') spaces,"T(molec) = ",tmole,Trim(TUnit)
      Write(dUnit,'(2x,2a,f10.3,1x,a)') spaces,"Kinetic E= ", &
          config_kineticEnergy(species,i)/config_getnmoles(species,i), &
          Trim(nrgUnit)
      Write(dUnit,'(2x,2a,f10.3,1x,2a)') spaces,"COM Vel. = ", &
          comVelocity%comp,Trim(velUnit)," (first molecule only)"
    End Do
  End Subroutine md_displayinit
  
  !----------------------------------------------------------------------------
  ! Displays information about the MD simulation in progress
  ! Requires:  imodel -- interaction model
  !            mdParams -- MD parameters
  !            species -- coordinates for system
  !            step -- time step number
  !            time -- total elapsed simulation time in picoseconds
  !            unit -- optional unit to write to
  ! Needed Improvements:
  ! 1) The dumped energy quantities are confusing because they are averaged
  !    in strange ways.  Rework
  !----------------------------------------------------------------------------
  Subroutine md_simdisplay(imodel,mdParams,species,scell,step,time,unit)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(MDInfo), Intent(In)                       :: mdParams
    Type(AtMolCoords), Dimension(:), Intent(In)    :: species
    Type(SimCell_Params), Intent(In)            :: scell
    Integer, Intent(In)                            :: step
    Real(kind=RDbl), Intent(In)                    :: time
    Integer, Intent(In), Optional                  :: unit

    Integer                       :: i,j,nsorbs,moles,nmoles,dUnit
    Real(kind=RDbl)               :: pe, ke, peavg, keavg, scaletime
    Real(kind=RDbl)               :: pei, kei, peavgi, keavgi
    Real(Kind=RDbl),Dimension(MAX_SORBS)               :: comvels
    Character(len=4)              :: units
    Character(len=strLen)         :: string
    Character(len=80)         :: comvelstring

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

    nsorbs = molecules_getnsorbs()

    Write(dUnit,'(a)') " "
    string = int2str(step)
    Write(dUnit,'(1x,2a,4x,a,f10.3,1x,a)') "MD Step: ",Trim(string), &
        "Time ",scaletime,units
    Write(dUnit,'(1x,a)') "--------------------------------------------------"

    !** Loop through the species and dump PE and KE information for each
    Do i = 1,nsorbs
      nmoles = config_getnmoles(species,i)
      If (config_isfixed(species(i))) Cycle
      If (nmoles==0) Cycle

      Write(dUnit,'(1x,a,a)') Trim(molecules_name(i)), &
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
            storestats_getintranrg(imodel%spcstats, i,'block','stretch')/nmoles,&
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
        Write(dUnit,'(3x,a16,4f15.6)') "Tosional PE    :", &
            storestats_getintranrg(imodel%spcstats, i,'inst','torsion')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'cavg','torsion')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'block','torsion')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'std','torsion')/nmoles
      End If
      If (interact_hasint(imodel,i,'intrapair')) Then
        Write(dUnit,'(3x,a16,4f15.6)') "Intramolec Pair:", &
            storestats_getintranrg(imodel%spcstats, i,'inst','intrapair')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'cavg','intrapair')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'block','intrapair')/nmoles, &
            storestats_getintranrg(imodel%spcstats, i,'std','intrapair')/nmoles
      End If
      If (interact_hasint(imodel,i,'intracoul')) Then
        Write(dUnit,'(3x,a16,4f15.6)') "Intramolec Coul:", &
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

    !** Write the potential energies for the species pairs
    Do i = 1,nsorbs
      If (config_isfixed(species(i))) Cycle
      nmoles = config_getnmoles(species,i)
      If (nmoles==0) Cycle

      Do j = i, nsorbs
        If (config_isfixed(species(i)).And.config_isfixed(species(j))) Cycle
        If ((i /= j).And..Not.(config_isfixed(species(j)))) &
            nmoles = nmoles + config_getnmoles(species,j)
        If (nmoles==0) Cycle
        Write(dUnit,'(1x,4a)') "Energy for ",Trim(molecules_name(i)),"-", &
            Trim(molecules_name(j))
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
      End Do
    End Do

    !** zero some PE and KE counters
    pe = 0.0_Rdbl
    ke = 0.0_Rdbl
    peavg = 0.0_Rdbl
    keavg = 0.0_Rdbl
    moles = 0

    !** Collect the total potential and kinetic energies
    !** VERY IMPORTANT : these should be totals, not per mole basis
    !** If you calculate per mole, then multicomponent total energy 
    !** calculation will will get messed up
    Do i = 1,nsorbs
      If (config_isfixed(species(i))) Cycle

      nmoles = config_getnmoles(species,i)
      If (nmoles==0) Cycle
      moles = moles + nmoles
      pe = storestats_getintranrg(imodel%spcstats, i,'inst') + pe
      peavg = storestats_getintranrg(imodel%spcstats, i,'cavg') + peavg
      ke = storestats_getke(imodel%spcstats, i,'inst') + ke
      keavg = storestats_getke(imodel%spcstats, i,'cavg') + keavg

      Do j = i, nsorbs
        If (config_isfixed(species(i)).And.config_isfixed(species(j))) Cycle
        If ((i /= j).And..Not.(config_isfixed(species(j)))) &
            nmoles = nmoles + config_getnmoles(species,j)
        pe = storestats_getnoncoul(imodel%spcstats, i,j,'inst')+ &
            storestats_getcoul(imodel%spcstats, i,j,'inst') + pe
        peavg = storestats_getnoncoul(imodel%spcstats, i,j,'cavg') + &
            storestats_getcoul(imodel%spcstats, i,j,'cavg') + peavg
      End Do
    End Do

    !** Get any extra energies from the integration routines
    !** e.g., Nose-Hoover energy, or others
    pei = 0.0_Rdbl
    kei = 0.0_RDbl
    peavgi = 0.0_RDbl
    keavgi = 0.0_Rdbl

    If (Associated(mdParams%mvparams(1)%integrate%vverlet)) Then
      Call velverlet_extraEnergy( &
          mdParams%mvparams(1)%integrate%vverlet,pei,peavgi,'pe')
      Call velverlet_extraEnergy( &
          mdParams%mvparams(1)%integrate%vverlet,kei,keavgi,'ke')
    Else If (Associated(mdParams%mvparams(1)%integrate%gear6)) Then
      Call gear6_extraEnergy( &
          mdParams%mvparams(1)%integrate%gear6,pei,peavgi,'pe')
      Call gear6_extraEnergy( &
          mdParams%mvparams(1)%integrate%gear6,kei,keavgi,'ke')
    End If

    If ((pei /= 0.0).And.(kei /= 0.0)) Then
      Write(dUnit,'(1x,a)') "Extra Integrator Energies (NH, etc) for &
          & first move_type only"
      Write(dUnit,'(3x,a16,2f15.6)') "Potential      :", &
          pei/moles,peavgi/moles
      Write(dUnit,'(3x,a16,2f15.6)') "Kinetic        :", &
          kei/moles,keavgi/moles
    End If

    !** Call moves to dislay any extra information we're not privy to
    Call moves_simdisplayExtra(mdParams%mvParams(1),imodel,dUnit,1)


    ! get com vels
    comvels(1:nsorbs)=zero
    comvelstring=Repeat(" ",80)
    Do i=1,nsorbs
        If (config_isfixed(species(i))) Cycle
        If (config_getnmoles(species(i))==0) Cycle
        comvels(i)=mag(config_getSpcCOMVel(species,i))
        comvelstring=Trim(comvelstring)//"SPC-"//&
            Trim(int2str(i))//": "//Trim(real2str(comvels(i)))//", "
    End Do

    pe = (pe + pei)/moles
    ke = (ke + kei)/moles
    peavg = (peavg + peavgi)/moles
    keavg = (keavg + keavgi)/moles

    Write(dUnit,'(1x,a)') "System Energy Totals(per total moles)"
    Write(dUnit,'(3x,a16,2f15.6)') "Total Potential:",pe,peavg
    Write(dUnit,'(3x,a16,3f15.6)') "Total Kinetic  :",ke,keavg
    Write(dUnit,'(3x,a16,2f15.6)') "Total Energy   :",pe+ke,peavg+keavg
    Write(dUnit,'(3x,a)')  "COM Velocity (for each species) : "&
        //Trim(comvelstring)
    Write(dUnit,*) "Total # Moles ",moles
    Write(dUnit,'(1x,a)') "--------------------------------------------------"

    If (NEMD_IS_ON) Then
      Call nemd_simdisplay(species,scell,time,dUnit, 6)
      Write(dUnit,'(4x,a)') '-----------------------------------------------'
    Endif
 
  End Subroutine md_simdisplay

  !----------------------------------------------------------------------------
  ! Displays information about the MD simulation in progress
  !----------------------------------------------------------------------------
  Subroutine md_simsummary(imodel,mdParams,species,step,time,unit)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(MDInfo), Intent(In)                       :: mdParams
    Type(AtMolCoords), Dimension(:), Intent(In)    :: species
    Integer, Intent(In)                            :: step
    Real(kind=RDbl), Intent(In)                    :: time
    Integer, Intent(In), Optional                  :: unit

    Integer                    :: i,j,nsorbs,moles,nmoles,dUnit
    Real(kind=RDbl)            :: pe, ke, peavg, keavg, pei, totnrg
    Real(kind=RDbl)            :: kei, peavgi, keavgi, scaletime, intrae
    Character(len=4)           :: units
    Character(len=strLen)      :: string

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

    If (Present(unit)) Then
      dUnit = unit
    Else
      dUnit = 6
    End If

    pe = 0.0_Rdbl
    ke = 0.0_Rdbl
    peavg = 0.0_Rdbl
    keavg = 0.0_Rdbl
    moles = 0
    nsorbs = molecules_getnsorbs()

    Do i = 1, nsorbs
      If (config_isfixed(species(i))) Cycle
      nmoles = config_getnmoles(species,i)
      If (nmoles==0) Cycle
      peavg = storestats_getintranrg(imodel%spcstats, i,'cavg')/nmoles + peavg
      keavg = storestats_getke(imodel%spcstats, i,'cavg')/nmoles + keavg
      moles = moles + nmoles
      Do j = i, nsorbs
        If ((j /= i).And..Not.config_isfixed(species(j))) &
            nmoles = nmoles + config_getnmoles(species,j)
        peavg = storestats_getnoncoul(imodel%spcstats, i,j,'cavg')/nmoles + peavg
      End Do
    End Do

    !** Get any extra energies from the integration routines
    !** e.g., Nose-Hoover energy, or others
    pei = 0.0_Rdbl
    kei = 0.0_RDbl
    peavgi = 0.0_RDbl
    keavgi = 0.0_Rdbl

    If (Associated(mdParams%mvparams(1)%integrate%vverlet)) Then
      Call velverlet_extraEnergy( &
          mdParams%mvparams(1)%integrate%vverlet,pei,peavgi,'pe')
      Call velverlet_extraEnergy( &
          mdParams%mvparams(1)%integrate%vverlet,kei,keavgi,'ke')
    Else If (Associated(mdParams%mvparams(1)%integrate%gear6)) Then
      Call gear6_extraEnergy( &
          mdParams%mvparams(1)%integrate%gear6,pei,peavgi,'pe')
      Call gear6_extraEnergy( &
          mdParams%mvparams(1)%integrate%gear6,kei,keavgi,'ke')
    End If
!!$    Call integrate_extraEnergy(mdParams%mvParams%integParams,pei,peavgi,'pe')
!!$    Call integrate_extraEnergy(mdParams%mvParams%integParams,kei,keavgi,'ke')

    peavgi = peavgi/moles
    keavgi = keavgi/moles

    Write(dUnit,*) ""
    Write(dUnit,*) ""
    Write(dUnit,'(1x,a)') "-----------------------------------------------------"
    Write(dUnit,'(20x,a)')        "MD SUMMARY"
    Write(dUnit,'(12x,a,i10)')       "Total Steps: ",step
    Write(dUnit,'(12x,a,f10.3,1x,a)') "Total Time : ",scaletime,units

    Write(dUnit,'(1x,a)') "====================================================="
    Write(dUnit,'(15x,a)') "System Energy Totals"
    Write(dUnit,'(3x,a,t35,f15.3,1x,a)') "Average Potential Energy: ", &
        peavg,eunits
    Write(dUnit,'(3x,a,t35,f15.3,1x,a)') "Average Kinetic Energy: ", &
        keavg,eunits

    If ((pei /= 0.0).And.(kei /= 0.0)) Then
      Write(dUnit,'(3x,a,t35,f15.3,1x,a)') "Average Extended Pot Nrg: ", &
          peavgi,eunits
      Write(dUnit,'(3x,a,t35,f15.3,1x,a)') "Average Extended Kin Nrg: ", &
          keavgi,eunits
    End If

    Write(dUnit,'(3x,a,t35,f15.3,1x,a)') "Average Total Energy: ", &
        peavg+keavg+keavgi+peavgi,eunits

    totnrg = interact_totnrg(imodel,step,.False.) 
    Write(dUnit,'(3x,a,t35,f18.6,1x,a)') "Total System Energy (wo intra): ", &
        totnrg,eunits

    Write(dUnit,'(1x,a)') "====================================================="

    Do i = 1,nsorbs
      If (config_isfixed(species(i))) Cycle
      intrae = 0.0_RDbl
      nmoles = config_getnmoles(species,i)
      If (nmoles==0) Cycle
      Write(dUnit,'(1x,a,a)') Trim(molecules_name(i)), &
          " Summary Information"
      Write(dUnit,'(3x,a20,f15.3,1x,a)') "T(mole)      : ", &
          storestats_gettemp(imodel%spcstats, i,'tmole','cavg'), tunits
      Write(dUnit,'(3x,a20,f15.3,1x,a)') "T(atom)      : ", &
          storestats_gettemp(imodel%spcstats, i,'tatom','cavg'), tunits
      If (interact_hasint(imodel,i,'stretch')) Then
        Write(dUnit,'(3x,a20,f15.3,1x,a)') "Stretch Energy : ", &
            storestats_getintranrg(imodel%spcstats, i,'cavg','stretch')/nmoles, &
            eunits
        intrae = intrae + &
            storestats_getintranrg(imodel%spcstats, i,'cavg','stretch')
      End If
      If (interact_hasint(imodel,i,'constraint')) Then
        Write(dUnit,'(3x,a20,f15.3,1x,a)') "Cnstr Energy : ", &
            storestats_getintranrg(imodel%spcstats, i,'cavg', &
            'constraint')/nmoles,eunits
      End If
      If (interact_hasint(imodel,i,'bending')) Then
        Write(dUnit,'(3x,a20,f15.3,1x,a)') "Bending Energy : ", &
            storestats_getintranrg(imodel%spcstats, i,'cavg','bending')/nmoles, &
            eunits
        intrae = intrae + storestats_getintranrg(imodel%spcstats, i,'cavg', &
            'bending')
      End If
      If (interact_hasint(imodel,i,'torsion')) Then
        Write(dUnit,'(3x,a20,f15.3,1x,a)') "Torsion Energy : ", &
            storestats_getintranrg(imodel%spcstats, i,'cavg','torsion')/nmoles, &
            eunits
        intrae = intrae + storestats_getintranrg(imodel%spcstats, i,'cavg', &
            'torsion')
      End If
      If (interact_hasint(imodel,i,'intrapair')) Then
        Write(dUnit,'(3x,a20,f15.3,1x,a)') "Intra Pair Nrg : ", &
            storestats_getintranrg(imodel%spcstats, i,'cavg', &
            'intrapair')/nmoles, eunits
        intrae = intrae + storestats_getintranrg(imodel%spcstats, i,'cavg', &
            'intrapair')
      End If
     If (interact_hasint(imodel,i,'intracoul')) Then
       Write(dUnit,'(3x,a20,f15.3,1x,a)') "Intra Coul Nrg : ", &
           storestats_getintranrg(imodel%spcstats, i,'cavg', &
           'intracoul')/nmoles, eunits
       intrae = intrae + storestats_getintranrg(imodel%spcstats, i, &
           'cavg','intracoul')
     End If

     Write(dUnit,'(3x,a20,f15.3,1x,a)') "Total Intramolecular Energy : ", &
         intrae/nmoles, eunits

     Write(dUnit,'(3x,a20,f15.3,1x,a)') "Total Kinetic Energy : ", &
         storestats_getke(imodel%spcstats, i,'cavg')/nmoles, eunits

    End Do

    Write(dUnit,'(1x,a)') "====================================================="
    Write(dUnit,'(15x,a)') "Pairwise Interaction Energies"

    Do i = 1,nsorbs
      Do j = i, nsorbs
        nmoles = config_getnmoles(species,i)
        If (nmoles==0) Cycle
        If (config_isfixed(species(i)).And.config_isfixed(species(j))) Cycle
        If ((i /= j).And..Not.config_isfixed(species(j))) Then
          nmoles = config_getnmoles(species,j) + nmoles
        End If
        Write(dUnit,'(1x,4a)') "Energy for ",Trim(molecules_name(i)),"-", &
            Trim(molecules_name(j))
        string = int2str(nmoles)
        Write(dUnit,'(3x,a20,f15.3,2(1x,a),2a)') "Noncoul Inter Energy : ", &
            storestats_getnoncoul(imodel%spcstats, i,j,'cavg')/nmoles, &
            eunits, "(total of ",Trim(string)," mols)"
        Write(dUnit,'(3x,a20,f15.3,2(1x,a),2a)') "Coul Inter Energy : ", &
            storestats_getcoul(imodel%spcstats, i,j,'cavg')/nmoles, &
            eunits, "(total of ",Trim(string)," mols)"
            
      End Do
    End Do

    Write(dUnit,'(1x,a)') "-----------------------------------------------------"

  End Subroutine md_simsummary

  !----------------------------------------------------------------------------
  ! What the hell is this?
  !----------------------------------------------------------------------------
  Subroutine md_updatecoords
  End Subroutine md_updatecoords

  !----------------------------------------------------------------------------
  ! Returns the value of the time step used for the given MD simulation
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function md_getdt(mdparams,n)
    Type(MDInfo), Intent(In)      :: mdparams
    Integer, Intent(In), Optional :: n

    Integer               :: thisOne
    Real(Kind=RDbl)       :: val
    Integer, Dimension(2) :: intarray
    
    thisOne = 1
    If (Present(n)) thisOne = n

    Call moves_getparam(mdparams%mvParams(thisOne),"TIMESTEP",val,intarray)
    md_getdt = val
  End Function md_getdt
      
  !------------------------------------------------
  !returns the name of the thermostat used
  !------------------------------------------------
  Function md_getthermostat(mvparams)
    Character(len=strLen) :: md_getthermostat
    Type(Move_Params) ,Intent(in) :: mvparams
    Logical :: found
    md_getthermostat= "NULL"
    If (mvparams%integrate%ensemble=="NVT") Then
      found=.False.
      If (Associated(mvparams%integrate%gear6)) Then
        If (Associated(mvparams%integrate%gear6%thermostat%vr)) Then 
          md_getthermostat= "VELOCITYRESCALE"
          found=.True.
        Endif
        If (Associated(mvparams%integrate%gear6%thermostat%nh)) Then
          md_getthermostat= "NOSEHOOVER"
          found=.True.
        Endif
      Else If (Associated(mvparams%integrate%vverlet)) Then
        If (Associated(mvparams%integrate%vverlet%thermostat%vr)) Then
          md_getthermostat= "VELOCITYRESCALE"
          found=.True.
        Endif
      Else If (Associated(mvparams%integrate%lfrog)) Then
        If (Associated(mvparams%integrate%lfrog%thermostat%gauss)) Then
          md_getthermostat= "GAUSS"
          found=.True.
        Endif
      Endif
      If (.Not.found) Then
        Write(*,*) "Then may be your thermostat is not designed& 
            & to work; check it in ctrlfile"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      endif
    Endif
  End Function md_getthermostat
  
End Module md






