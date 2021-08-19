!------------------------------------------------------------------------------
! This 'interact' module serves to wrap the forcefield evaluations and the
! storage for the interaction results into a single module.  It contains
! the 'Interaction_Model' data type that defines the forcefield(s), 
! interactions storage and some statistics for a simulation.  
! 
! An equally important reason for wrapping all the interaction routines and
! update routines into one module is to provide universal support for 
! techniques that modify the system interactions or use multiple forcefield
! to determine system properties.  Examples of this class of simulation 
! techniques are Umbrella sampling and Voter's Hyperdynamics technique.
! Umbrella sampling uses results from one artifical forcefield to guide the 
! sampling of an ensemble, but computes the ensemble average using both
! the artificial and the "real" forcefield. 
!
! ALL FORCEFIELD EVALUATIONS AND STATISTICS UPDATES GO THROUGH THIS MODULE
!
! Important Routines:
!   interact_init -- initializes an interaction model
!   interact_updatestats -- updates the stored statistics
!   interact_initnrgs -- evaluates all interactions and stores statistics
!   interact_int -- evaluates interactions between subsets of the system
!   interact_checknrgs -- compares statistics storage to reevaluated energies
!   interact_changenmoles -- changes number of molecules in interaction storage
!
! Needed Improvements:
! 1) still lacks a completed control file section
! 2) modify forcefield and lower routines to check for maximum energy
!------------------------------------------------------------------------------

Module interact

  Use defaults, Only: RDbl,strLen,lstrLen,dashedline,zero,one,dbgflag, &
      NO_OF_INTRA_POTS,TOTAL_INDEX,caltoj,MAX_SORBS,kcalmole_kb,xlstrLen, &
      scalef, scalepe
  Use utils, Only: isfileopen,filesrchstr,stripcmnt,toupper,allocErrDisplay, &
      int2str,real2str,deallocErrDisplay,checkandstop,getdepth
  Use file, Only: file_open
  Use vector, Only: Assignment(=), mag
  Use molecules, Only: molecules_getnsorbs,molecules_getnatoms,molecules_name,&
      molecules_checkinit,molecules_gettype, molecules_getnthatom
  Use config, Only: AtMolCoords, config_isfixed, config_kineticenergy, &
      config_getnmoles, config_getnatoms, config_kineticTemp, &
      config_getnmoleslist, config_checkinit, config_changedof
  Use simcell, Only: SimCell_Params, simcell_getmolectype
  Use storestats, Only: Species_Stats, storestats_init, &
      storestats_incrAllIntraNrg,storestats_incrcoulnrg, &
      storestats_incrnoncoulnrg, storestats_incrkinnrg, &
      storestats_updateEnergySS, storestats_updatenrg, &
      storestats_updateintraenergy, storestats_updatetemp, &
      storestats_getallintra, storestats_getnoncoul, &
      storestats_getcoul, storestats_displaynrg, storestats_clean
  Use forcefield, Only: Forcefield_Info, forcefield_display, forcefield_init, &
      forcefield_initstore, forcefield_allint, forcefield_putaccels, &
      forcefield_hasint, forcefield_putaccels2, forcefield_msysint, &
      forcefield_conint, forcefield_listparams, forcefield_asint, &
      forcefield_clean, forcefield_getpotparameters, forcefield_boxinfo, &
      forcefield_setparams
  Use storetop, Only: Forcefield_Results,storetop_display,storetop_sum, &
      storetop_zero, storetop_fastzero, storetop_fastsum, storetop_init, &
      storetop_initcopy, storetop_copy, storetop_clean, storetop_totnrg, &
      storetop_update, storetop_extract, storetop_allnrgs, storetop_chklink, &
      storetop_changenmoles, storetop_chksizes, storetop_scalenrgs, &
      storetop_chkequal, storetop_fillsub, storetop_nmolesDecr
  Use storebase, Only: EnergyPlus, storebase_sumintra, storebase_copy, &
      storebase_clean, storebase_display, storebase_init, storebase_disp, &
      storebase_chkintra, storebase_chkequal, storebase_zero, &
      storebase_nrg, storebase_totnrg, storebase_initcopy
  Use boost, Only: BoostParams, boost_init, boost_calcboost, boost_sampleCF, &
      boost_display,boost_displayStats, boost_getBoostMolname, &
      boost_calcBoostHack
  Use umbrella, Only: Umbrella_Params,umbrella_init,umbrella_value, &
      umbrella_display,umbrella_clean

  Implicit None
  Save

  Private 
  Public :: Interaction_Model, interact_init, default_interact_tag, &
      interact_changenmoles, interact_int, interact_extract, &
      interact_updatestats, interact_hasint, interact_checknrgs, &
      interact_initnrgs, interact_dispDynInfo, &
      interact_usampling, interact_chkrecalc, interact_simpleint, &
      interact_totnrg, interact_display, interact_clean, &
      interact_getpotparameters, interact_listparams, &
      interact_chkmolnums, interact_resize, &
      interact_changeAllNmoles, interact_boxinfo, interact_setparams

  !** Stores all information about and from interaction model
  Type Interaction_Model
    Integer                        :: nderivs,nff
    Type(Species_Stats)            :: spcstats    
    Character(len=lstrLen)         :: keys
    Type(BoostParams), Pointer     :: boost
    Type(Umbrella_Params), Pointer :: usampling
    Type(Forcefield_Info), Dimension(:), Pointer     :: ff
    Type(Forcefield_Results), Dimension(:), Pointer  :: results

    !** Extra arrays for processing:  (temporary HACK?)
    Type(EnergyPlus), Dimension(:), Pointer    :: temp
    Type(EnergyPlus), Dimension(:), Pointer    :: intra
    Real(kind=RDbl), Dimension(:,:), Pointer   :: noncoul,coul,pscalings
  End Type Interaction_Model

  Character(len=strLen), Parameter :: default_interact_tag = &
      "Interaction Model"

Contains
  !---------------------------------------------------------------------
  ! Initializes the molecule-system data structure
  ! Requires:  imodel -- interaction model data type to initialize
  !            ctrl_filename -- the control filename to take info from
  !            simcell -- the simulation cell
  !            species -- species data structure
  !            simtype -- string defining simulation type
  !            opt_interact_tag -- optional control file ID string 
  ! Currently supported simulation types are:
  ! 'MC' -- energy-only forcefield evaluations
  ! 'MD','HMC' -- energy and force forcefield evaluations
  !---------------------------------------------------------------------
  Subroutine interact_init(imodel,ctrl_filename,simcell,species,simtype, &
      opt_interact_tag)
    Type(Interaction_Model), Intent(Out)            :: imodel
    Character(*), Intent(In)                        :: ctrl_filename
    Type(SimCell_Params), Intent(In)                :: simcell  
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Character(*), Intent(In)                        :: simtype
    Character(*), Intent(In),Optional               :: opt_interact_tag

    Integer                    :: i,unitno,lineno,error
    Integer                    :: nspc
    Logical                    :: flag,flag2
    Character(len=strLen)      :: keyword,interact_tag
    Character(len=xlstrLen)    :: full_line
    Character(len=strLen), Dimension(10)  :: ffid

    !** Set the tag
    interact_tag = default_interact_tag
    If (Present(opt_interact_tag)) interact_tag = opt_interact_tag

    !** Make sure the molecules are initialized
    If (.Not. molecules_checkinit()) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          " molecules must be initialized before the interaction model"
      Stop
    End If

    !** Make sure the configuration is initialized
    If (.Not. config_checkinit(species)) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          " configuration must be initialized before the interaction model"
      Stop
    End If

    !** Defaults
    imodel%nff = 1  
    imodel%keys = ''

    !** Get the number of species
    nspc = molecules_getnsorbs()

    !** Nullify the pointers in the structure
    Call interact_null(imodel)

    !** Open the ctrl_file if it is not already open and rewind
    unitno = file_open(ctrl_filename,110)
    Rewind(unit = unitno)    

    !** Look for the special interaction section in the control file
    lineno = filesrchstr(unitno,interact_tag,full_line)
    If (lineno == 0) Then   !** just use the basic forcefield, without frills
      ffid(1) = 'BASIC'
      imodel%keys = 'DEFAULT'
    Else

      !** Read and process keywords
      keyword = ''
      Do 
        Read(unitno,*) keyword
        If (Index(keyword,'-----') /= 0) Exit

        Write(imodel%keys,'(a,1x,a)') Trim(imodel%keys),Trim(keyword)

        Select Case (keyword)
        Case ('DEFAULT','NORMAL','SINGLEFF')
          Read(unitno,*) ffid(1)

        Case ('BOOST')
          Read(unitno,*) ffid(1)
          Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
              " WARNING: Boosting still being tested, some comment about NVT?"

          !** Initialize the boost parameters
          Allocate(imodel%boost, stat=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'imodel%boost')
          Call boost_init(imodel%boost,species,unitno)

          !** Initialize the dummy storage arrays, HACK?
          Allocate(imodel%intra(nspc), stat=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'imodel%intra')
          Do i = 1,nspc
            Call storebase_init(imodel%intra(i),1,.True.)
          End Do
          Allocate(imodel%noncoul(nspc,nspc), stat=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'imodel%noncoul')
          Allocate(imodel%coul(nspc,nspc), stat=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'imodel%coul')
          Allocate(imodel%pscalings(nspc,nspc), stat=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__, &
              'imodel%pscalings')
          imodel%noncoul = 0.0_RDbl
          imodel%coul = 0.0_RDbl
          imodel%pscalings = 0.0_RDbl

        Case ('UMBRELLA')
          Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
              " WARNING: Umbrella sampling still being developed"
          Stop

          !** Initialize the umbrella sampling parameters
          Allocate(imodel%usampling, stat=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
          imodel%nff = 2
          Call umbrella_init(imodel%usampling,ctrl_filename,species,ffid)

        Case Default
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              ' could not understand keyword: ',Trim(keyword)
          Stop
        End Select
      End Do
    End If

    !** Allocate space for the forcefields and the results
    Allocate(imodel%ff(imodel%nff), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(imodel%results(imodel%nff), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Initialize the forcefields
    Do i = 1,imodel%nff
      Call forcefield_init(imodel%ff(i),ffid(i),ctrl_filename,species,simcell)
    End Do

    !** If DOF is changed during forcefield_init, should be done in config too
    Call config_changedof(species)

    !** Identify simulation type and set storage parameters appropriately
    Select Case(ToUpper(Trim(simtype)))
    Case ('MD')
      imodel%nderivs = 1
    Case ('MC')
      imodel%nderivs = 0
    Case ('HMC')
      imodel%nderivs = 1
!!$      imodel%nderivs = 0  !** just for checking gcmc part of hmc
    Case Default
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Could not interpret simulation type ',Trim(simtype)
      Stop
    End Select

    !** Initialize the storage for the forcefield(s)
    Do i = 1,imodel%nff
      Call forcefield_initstore(imodel%results(i),imodel%ff(i),species, &
          imodel%nderivs,imodel%ff(i)%storelevel)
    End Do

    !** Initialize the temporary interaction structures for extractions
    Allocate(imodel%temp(imodel%nff), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1,imodel%nff
      Call storebase_initcopy(imodel%temp(i),imodel%results(i)%total)
    End Do

    !    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    !    Call storetop_display(imodel%results(1),.False.,2,6)
    !    stop

    !** Initialize species statistics
    If (Associated(imodel%usampling)) Then
      Call storestats_init(imodel%spcstats,.True.)
    Else
      Call storestats_init(imodel%spcstats,.False.)
    End If

  End Subroutine interact_init

  !---------------------------------------------------------------------
  ! Nullify the interaction model data type
  ! Requires:  imodel -- interaction model data type to initialize
  !---------------------------------------------------------------------
  Subroutine interact_null(imodel)
    Type(Interaction_Model), Intent(InOut)         :: imodel

    Nullify(imodel%boost)
    Nullify(imodel%usampling)
    Nullify(imodel%temp)
    Nullify(imodel%intra)
    Nullify(imodel%noncoul)
    Nullify(imodel%coul)
    Nullify(imodel%pscalings)

  End Subroutine interact_null

  !-------------------------------------------------------------------------
  ! Calculates interactions between two subsets of the system as perscribed
  ! by the interaction model and puts the results into a forcefield results
  ! data structure.  
  ! 
  ! Currently only does evaluations using the primary forcefield
  ! Requires:  imodel -- interaction model information
  !            results -- results of the evaluation
  !            species -- coordinate data structure 
  !            simcell -- simulation cell information
  !            fast -- logical flag indicating fast or slow evaluations
  !            recalc -- forces energy (re)calculation if True
  !            skip_intra -- skips intramolecular interaction if True
  !            calc_accels -- calculate accels and constraints if True
  !            auxparams -- additional parameters (maxnrg,timestep,T_Kelvin)
  !            subset1 -- 1st subset
  !            subset2 -- 2nd subset, optional, default is full system
  ! NOTE: currently setup NOT to be able to copy from storage.  It's not
  ! clear if all interactions or just the sums should be copied, so I've
  ! left it unimplemented.
  !-------------------------------------------------------------------------
  Logical Function interact_int(imodel,results,species,simcell,fast, &
      recalc,skip_intra,calc_accels,auxparams,subset1,subset2)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(Forcefield_Results), Intent(InOut)        :: results
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell  
    Logical, Intent(In)                            :: fast,recalc
    Logical, Intent(In)                            :: skip_intra,calc_accels
    Real(kind=RDbl), Dimension(:), Intent(In)      :: auxparams
    Integer, Dimension(3), Intent(In)              :: subset1
    Integer, Dimension(3), Intent(In), Optional    :: subset2

    Integer         :: depth1,depth2
    Logical         :: evaluated,success
    Real(kind=RDbl) :: nrg

    !** Determine if it's desired and possible to use stored interactions
    If (.Not. recalc) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' interact_int is not presently setup to copy from storage'
      Stop
    End If

    !** Zero the passed storage structure (including intramolecular)
    If ((results%nderivs == 1).And.(.Not. results%subsetonly)) Then
      Call storetop_fastzero(results,.False.)
    Else
      Call storetop_zero(results,.False.)
    End If

    !** Determine the depth of the subsets
    depth1 = getdepth(subset1)
    depth2 = 0
    If (Present(subset2)) depth2 = getdepth(subset2)

    !** Call forcefield to get specified subset1-subset2 interactions
    evaluated = .False.
    Select Case (depth1) 
    Case(0) !** full system -- subset2
      Select Case (depth2) 
      Case (0)  !** full system interactions
        interact_int = forcefield_allint(imodel%ff(1),results,species, &
            simcell,fast,.False.)
        evaluated = .True.
      End Select

    Case(2) !** molecule interacting with subset2 
      Select Case (depth2) 
      Case (0)  !** molecule -- system interactions
        interact_int = forcefield_msysint(imodel%ff(1),results,species, &
            subset1(1),subset1(2),simcell,fast,skip_intra,.False.)
        evaluated = .True.
      End Select

    End Select

    !** Check that interactions were evaluated
    If (.Not. evaluated) Then
      Write(0,'(1x,2a,i4,a,2i3)') __FILE__,' : ',__LINE__, &
          ' Requested subset-subset interactions not available for depths: ', &
          depth1,depth2
      Stop      
    End If

    !** Check that interaction evaluation was successful
    If (.Not. interact_int) Return

    !** Make sure the internal sums in the results are correct
    If ((results%nderivs == 1).And.(.Not. results%subsetonly)) Then
      Call storetop_fastsum(results)
    Else
      Call storetop_sum(results)
    End If

    !** Check the maximum energy criterium here, it would be best to
    !** do this incrementally in the lower routines
    nrg = storetop_totnrg(results,.True.)
    If (nrg > auxparams(1)) Then
      interact_int = .False.
      !      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      !      Write(*,*) 'Maximum energy criterium violated (max,current): ', &
      !          auxparams(1),nrg
      !      Return
    End If

    !** Convert the forces into accelerations and copy into species 
    If (calc_accels) Then
      Call forcefield_putaccels2(imodel%ff(1),results,species,fast)
    End If

    !** Do post-evaluation modifications as necessary
    !** Ideally, these would operate on forces, FIX and move above accel calcs
    If (Associated(imodel%boost)) Then
      Call interact_boost(imodel,results,species,auxparams)
    End If

    !** Get the constraint forces in acceleration form
    !** NOTE: this doesn't place the constraint energies properly, FIX
    success = forcefield_conint(imodel%ff(1),species,fast)
    Call checkandstop(success,__FILE__,__LINE__, &
        ' Could not calculate constraints')

  End Function interact_int

  !----------------------------------------------------------------------
  ! Reinitializes the storage portion of the interaction model structure.
  ! This is for use when the system size or number of species changes 
  ! drastically and the _changenmoles can't be used.
  ! Requires:  imodel -- interaction model data type to initialize
  !            ctrl_filename -- the control filename to take info from
  !            simcell -- the simulation cell
  !            species -- species data structure
  !----------------------------------------------------------------------
  Subroutine interact_resize(imodel,simcell,species)
    Type(Interaction_Model), Intent(InOut)       :: imodel
    Type(SimCell_Params), Intent(In)             :: simcell  
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species

    Integer                    :: i,unitno,lineno,error
    Integer                    :: nspc
    Logical                    :: flag,flag2
    Character(len=strLen)      :: keyword,interact_tag
    Character(len=xlstrLen)    :: full_line
    Character(len=strLen), Dimension(10)  :: ffid

    !** Get the number of species
    nspc = molecules_getnsorbs()

    !** Deallocate the results storage
    Call storestats_clean(imodel%spcstats)

    Do i = 1,imodel%nff
      Call storetop_clean(imodel%results(i))
    End Do
    Deallocate(imodel%results, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'results')

    !** Reallocate the results storage
    Allocate(imodel%results(imodel%nff), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Initialize the storage for the forcefield(s)
    Do i = 1,imodel%nff
      Call forcefield_initstore(imodel%results(i),imodel%ff(i),species, &
          imodel%nderivs,imodel%ff(i)%storelevel)
    End Do

    !    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    !    Call storetop_display(imodel%results(1),.False.,2,6)
    !    stop

    !** Initialize the temporary interaction structures for extractions
    Allocate(imodel%temp(imodel%nff), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1,imodel%nff
      Call storebase_initcopy(imodel%temp(i),imodel%results(i)%total)
    End Do

    !** Initialize species statistics
    If (Associated(imodel%usampling)) Then
      Call storestats_init(imodel%spcstats,.True.)
    Else
      Call storestats_init(imodel%spcstats,.False.)
    End If

  End Subroutine interact_resize

  !-------------------------------------------------------------------------
  ! Handle post-evaluation modification of the interaction results to 
  ! obtaining the boosting of the HyperMD method.  The accelerations
  ! must have already been calculated from results before this routine is 
  ! called.  THIS MAY BE BROKEN, please check
  ! Requires:  imodel -- interaction model information
  !            results -- results of the evaluation
  !            species -- coordinate data structure 
  !            auxparams -- additional parameters (maxnrg,timestep,T_Kelvin)
  !-------------------------------------------------------------------------
  Subroutine interact_boost(imodel,results,species,auxparams)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(Forcefield_Results), Intent(InOut)        :: results
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Real(kind=RDbl), Dimension(:), Intent(In)      :: auxparams

    Integer                :: a,error
    Integer                :: natoms,nmoles,spc
    Real(kind=RDbl)        :: betaKcal,dt,TKelvin
    Real(Kind=RDbl), Dimension(:), Allocatable, Save :: totalPE

    !** Extract from auxilliary parameters and calculate beta
    TKelvin = auxparams(3)
    betaKcal = 1.0_RDbl/(TKelvin*kcalmole_kb)
    dt = auxparams(2)

    !** Get the molecule type to boost and information about those molecules
    spc = molecules_gettype(boost_getBoostMolname(imodel%boost))
    nmoles = config_getnmoles(species,spc)
    natoms = molecules_getnatoms(spc)

    !** Size the totalPE array if this is the first time
    If (.Not. Allocated(totalPE)) Then
      Allocate(totalPE(nmoles),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'totalPE temp array')
    End If

    !** If there is only one molecule, we can just pass the spc-spc
    !** energy to boost_calcboost.  Otherwise, more complex (working?)
    If (nmoles == 1) Then
!!$      !MDEBUG
!!$      !** I don't really understand this storage stuff, so this is
!!$      !** an implementation to the best of my knowledge to get the
!!$      !** molecule's atom-species energies out.
!!$      Do a = 1,natoms
!!$        success = storetop_extract(imodel%results(1),.False.,.True.,.True., &
!!$            tmpNrg(a),(/1,1,a/))
!!$        Write(0,*) __FILE__,__LINE__,a,success,tmpNrg(a)%nrg
!!$      End Do

      !** Get the spc-spc energies for input to boost routines
      Call storetop_allnrgs(imodel%results(1),imodel%noncoul, &
          imodel%coul,imodel%intra)

      !** Call the boosting routine. If we boost, this will return 
      !** the boosted accelerations, timestep and potential scaling factors
      Call boost_calcBoost(imodel%boost,species,imodel%noncoul,imodel%coul, &
          betaKcal,dt,imodel%pscalings)

    Else
      !** In this case, we need the individual molecule-species energies.
      !** Extract the energies from the storage module. In this case,
      !** all the energies are summed together. I do this to keep the 
      !** force evaluations simple. Otherwise, separate calls with 
      !** separate forces might be required.  NOT READY
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ': Sorry, need more code before running Accelerated MD with more ', &
          'than 1 molecule'
      Stop

      !** Ok, with the TOTAL potential energy calculated for each 
      !** molecule, we can pass the list to boost_calcBoost to sort
      !** it all out.
      Call boost_calcBoost(imodel%boost,species,totalPE,betaKcal,dt, &
          imodel%pscalings)
    End If

    !** Scale the potentials in the storage according to the scaling factors
    Call storetop_scalenrgs(imodel%results(1),imodel%pscalings, &
        (/.True.,.True.,.True./))

  End Subroutine interact_boost

  !---------------------------------------------------------------------
  ! Calculates interactions for a given subset1-subset2 pair.  It does
  ! not use the storage structure to keep this fast and simple
  ! Requires:  imodel -- interaction model information
  !            ffout -- output from forcefield
  !            subset1 -- first system subset
  !            subset2 -- second system subset
  !            species -- coordinate data structure 
  !            simcell -- simulation cell information
  !            fast -- logical flag indicating fast or slow evaluations
  !---------------------------------------------------------------------
  Logical Function interact_simpleint(imodel,ffout,subset1,subset2,species, &
      simcell,fast,hot)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(EnergyPlus), Intent(InOut)                :: ffout
    Integer, Dimension(3), Intent(In)              :: subset1,subset2
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell  
    Logical, Intent(In)                            :: fast
    Real(Kind=RDbl), Dimension(:), Optional :: hot

    !** Zero the storage structure
    Call storebase_zero(ffout)

    !** Call forcefield to get all the interactions
    If ((getdepth(subset1) == 3).And.(getdepth(subset2) == 1)) Then
      If (.Not.Present(hot)) Then
        interact_simpleint = forcefield_asint(imodel%ff(1),ffout,species, &
            subset1(1),subset1(2),subset1(3),subset2(1),simcell,fast,.False.)
      Else
        interact_simpleint = forcefield_asint(imodel%ff(1),ffout,species, &
            subset1(1),subset1(2),subset1(3),subset2(1),simcell,fast,.False.,hot)
      End If
    Else
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' currently only setup to do atm-spc interactions, sorry'
      Stop
    End If

  End Function interact_simpleint

  !----------------------------------------------------------------------------
  ! Extracts a specified energy from the forcefield results storage using the
  ! interaction model.  Returns False if the quantity is not available.  The
  ! UNITS on the extracted quantities are kcal/mol and Angstroms.
  ! Requires:  imodel -- interaction model information
  !            iter -- iteration number
  !            intra -- flag indicating if intramolecular interactions sought
  !            ncoul -- flag indicating if Non-Coulombic interactions sought
  !            coul -- flag indicating if coulombic interactions sought
  !            info -- returned energy plus possible derivatives structure
  !            normfactor -- normalization factor, usually unity 
  !            subset1 -- 1st spc,mol,atm specifier
  !            subset2 -- 2nd spc,mol,atm specifier
  ! a "subset" is a 1D array (size > 0) with indices: species,molecule,atom
  ! Usage: the spc,molec or atm designators can also be set to zero.  This
  !        will tell the routine not descend to that level.  For example,
  !        spc1=1, molec1=0, spc2=2, molec2=0 => get overall interaction of 
  !        species 1 with species 2.
  !----------------------------------------------------------------------------
  Logical Function interact_extract(imodel,iter,intra,ncoul,coul,info, &
      normfactor,subset1,subset2)
    Type(Interaction_Model), Intent(InOut)       :: imodel
    Integer, Intent(In)                          :: iter
    Logical, Intent(In)                          :: intra,ncoul,coul
    Type(EnergyPlus), Intent(InOut)              :: info
    Real(kind=RDbl), Intent(Out)                 :: normfactor
    Integer, Dimension(3), Intent(In)            :: subset1 
    Integer, Dimension(3), Intent(In), Optional  :: subset2

    Integer            :: i,error
    Logical            :: success

    !** Defaults
    normfactor = one

    !** Extract the desired interactions for each forcefield
    Do i = 1,imodel%nff
      If (Present(subset2)) Then
        interact_extract = storetop_extract(imodel%results(i),intra,ncoul, &
            coul,imodel%temp(i),subset1,subset2)
      Else 
        interact_extract = storetop_extract(imodel%results(i),intra,ncoul, &
            coul,imodel%temp(i),subset1)
      End If
      If (.Not. interact_extract) Return
    End Do

    !** Process the interactions using the interaction model
    If (Associated(imodel%boost)) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' not able to do extractions with boost model yet'
      Stop
    Else If (Associated(imodel%usampling)) Then
      Call umbrella_value(imodel%usampling,imodel%temp,iter,info,normfactor)
    Else
      If (imodel%nff /= 1) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            " More than one forcefield, but no way to combine them"
        Stop
      End If

      Call storebase_copy(info,imodel%temp(1))
    End If

  End Function interact_extract

  !----------------------------------------------------------------------------
  ! Return the current total energy of the system using the _extract routine
  ! Requires:  imodel -- interaction model information
  !            iter -- iteration number
  !            intra -- if True, include intramolecular interactions in total
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function interact_totnrg(imodel,iter,intra)
    Type(Interaction_Model), Intent(InOut)       :: imodel
    Integer, Intent(In)                          :: iter
    Logical, Intent(In)                          :: intra

    Logical, SAVE           :: firsttime = .True.
    Logical                 :: success
    Real(kind=RDbl)         :: normfactor
    Type(EnergyPlus), SAVE  :: components

    !** Initialize the temporary energy structure for the first time
    If (firsttime) Then
      firsttime = .False.
      Call storebase_initcopy(components,imodel%results(1)%total)
    End If

    !** Use the extract routine
    success = interact_extract(imodel,iter,.True.,.True.,.True.,components, &
        normfactor,(/0,0,0/))
    Call checkandstop(success,__FILE__,__LINE__, &
        ' Could not extract total energy')

    !** Return the total energy, with or with intramolecular contributions
    If (intra) Then
      interact_totnrg = storebase_totnrg(components,.True.)
    Else
      interact_totnrg = storebase_nrg(components)
    End If

  End Function interact_totnrg

  !--------------------------------------------------------------------------
  ! Changes the number of molecules in the full system interaction 
  ! storage structure.  
  ! Requires:  imodel -- interaction model data type 
  !            spc -- species number for which to change nmoles 
  !            nmoles -- new number of molecules
  !            delmol -- optional molecule to delete
  !--------------------------------------------------------------------------
  Subroutine interact_changenmoles(imodel,spc,nmoles,delmol)
    Type(Interaction_Model), Intent(InOut)      :: imodel
    Integer, Intent(In)                         :: spc,nmoles
    Integer, Intent(In), Optional               :: delmol

    Integer         :: i

    If (Present(delmol)) Then
      Do i = 1,imodel%nff
        Call storetop_changenmoles(imodel%results(i),spc,nmoles,delmol)
      End Do
    Else
      Do i = 1,imodel%nff
        Call storetop_changenmoles(imodel%results(i),spc,nmoles)
      End Do
    End If

  End Subroutine interact_changenmoles



  !--------------------------------------------------------------------------
  ! Changes the number of molecules in the full system interaction 
  ! storage structure.  This makes changes for all spc, 
  ! Requires:  imodel -- interaction model data type 
  !            nmoles_list -- new number of molecules
  ! Note : but not tested for all cases, use with care -Shaji
  !--------------------------------------------------------------------------
  Subroutine interact_changeAllNmoles(imodel,nmoles_list)
    Type(Interaction_Model), Intent(InOut)      :: imodel
    Integer, Dimension(:), Intent(In)                       :: nmoles_list

    Integer         :: i,spc,nspc

    nspc=molecules_getnsorbs()

    Do spc=1,nspc
      Do i = 1,imodel%nff
        Call storetop_nmolesDecr(imodel%results(i),spc,nmoles_list(spc))
      End Do
    End Do

  End Subroutine interact_changeAllNmoles

  !----------------------------------------------------------------------
  ! Updates the species-species energy storage by pulling energies out
  ! of the storage structure and putting them in the statistics storage.
  ! NOTE: converted from kcal/mol to KJ/mol before storing here
  ! Requires:  imodel -- interaction model data type to initialize
  !            species -- species data structure
  !            mdflag -- true => also update temps etc for md
  ! need improvements:
  ! 1) reaches too far into imodel, hack, use storetop intermediate
  ! note: there is apparently something magical about the arguments of
  !       this subroutine in the optimized code, it runs much faster
  !       just the way it is.  this is why _updatestats is separate.
  !----------------------------------------------------------------------
  Subroutine interact_updatestats(imodel,species,mdflag)
    Type(Interaction_Model), Intent(InOut)      :: imodel
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Logical, Intent(In)                         :: mdflag

    Integer             :: i,j,nspc
    Real(kind=RDbl)     :: tmole,tatom,ke

    nspc = molecules_getnsorbs()

    If (Associated(imodel%usampling)) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Stop
    End If

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'before updating in interact_updateallstats'
    Call storetop_display(imodel%results(1),.False.,2,6)
#endif

    !** make sure the sums in the storage are internally consistent
    Call storetop_fastsum(imodel%results(1))  !assumes accurate mol-mol sums

    Do i = 1,nspc
      Do j = i,nspc
        !** Update upper triangle interpair coul- and non-coul energies
        Call storestats_updateEnergySS(imodel%spcstats, i,j,'noncoul',&
            imodel%results(1)%ncoul%ab(i,j)%total%nrg*calToJ)
        Call storestats_updateEnergySS(imodel%spcstats, i,j,'coul',&
            imodel%results(1)%coul%ab(i,j)%total%nrg*calToJ)
      End Do

      !** Update intramolecular energies
      Call storebase_sumintra(imodel%results(1)%intra(i)%total,.False.)
      Call storestats_updateIntraEnergy(imodel%spcstats,i, &
          imodel%results(1)%intra(i)%total%intranrg*calToJ)
    End Do

    !** Update temperatures, velocity and kinetic nrg statistics if desired
    If (mdflag) Then
      Do i = 1,nspc
        If (.Not. config_isfixed(species(i))) Then
          !** Calculate the quantities
          tmole = config_kineticTemp(species,i,'MOLECULAR')
          tatom = config_kineticTemp(species,i,'ATOMIC')
          ke = config_kineticEnergy(species,i)

          Call storestats_updateTemp(imodel%spcstats,i,tmole,'molecular')
          Call storestats_updateTemp(imodel%spcstats,i,tatom,'atomic')
          Call storestats_updatenrg(imodel%spcstats,i,'ke',ke)
        End If
      End Do
    End If

  End Subroutine interact_updatestats

  !----------------------------------------------------------------------------
  ! checks to see if a particular intramolecular interaction is turned on
  ! requires:  imodel -- interaction model information
  !            spc -- species number
  !            whichone -- string identifier for model type
  !----------------------------------------------------------------------------
  Logical Function interact_hasint(imodel,spc,whichone)
    Type(Interaction_Model), Intent(InOut)    :: imodel
    Integer, Intent(In)                       :: spc
    Character(*), Intent(In)                  :: whichone

    interact_hasint = forcefield_hasint(imodel%ff(1),spc,whichone)

  End Function interact_hasint

  !----------------------------------------------------------------------------
  ! Returns True if Umbrella sampling is turned on.
  ! requires:  imodel -- interaction model information
  !----------------------------------------------------------------------------
  Logical Function interact_usampling(imodel)
    Type(Interaction_Model), Intent(In)    :: imodel

    interact_usampling = .False.
    If (Associated(imodel%usampling)) Then
      interact_usampling = .True.
    End If

  End Function interact_usampling

  !---------------------------------------------------------------------------
  ! Initialize the species-species and intramolecular energy storage with
  ! an initial forcefield evaluation.  Used for Monte Carlo routines.
  ! requires:  imodel -- interaction model information
  !            fast -- logical flag indicating fast or slow evaluations
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            update -- if True, will update species statistics
  !---------------------------------------------------------------------------
  Subroutine interact_initnrgs(imodel,fast,species,simcell,update)
    Type(interaction_model), Intent(InOut)           :: imodel
    Logical, Intent(In)                              :: fast
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Type(simcell_params), Intent(In)                 :: simcell
    Logical, Intent(In)                              :: update

    Logical                          :: success,calc_accels
    Real(Kind=RDbl)                  :: dummydt,tkelvin

    !** Recalculate the entire system interactions
    dummydt = 1.0_Rdbl
    TKelvin = 300.0_Rdbl  !** HACK
    calc_accels = .False.
    If (imodel%nderivs > 0) calc_accels = .True.
    success = interact_int(imodel,imodel%results(1),species,&
        simcell,fast,.True.,.False.,calc_accels, &
        (/1000.0_RDbl,dummydt,TKelvin/),(/0,0,0/))

    Call checkandstop(success,__FILE__,__LINE__, &
        ' unexpected problem during forcefield evaluation')

#ifdef
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'INITIALIZATION results'
    Call storetop_display(imodel%results(1),.False.,2,6)
    !    Stop
#endif

    !** Update the spc-spc statistics
    If (update) Then
      Call interact_updatestats(imodel,species,.False.)
    End If

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call storestats_displaynrg(imodel%spcstats,2,6)
    Stop
#endif 

  End Subroutine interact_initnrgs

  !---------------------------------------------------------------------------
  ! This routine compares the stored species-species energies with those from
  ! a fresh calculation over the whole system.  It is essential for error 
  ! checking at the end of simulation types that only use partial system
  ! energy evaluations, such as monte carlo simulations.  This function 
  ! returns the RMS of the differences in individual species-species and
  ! intramolecular energies.
  ! Requires:  imodel -- interaction model information
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            indent -- no. of spaces from the left margin
  !            unit -- optional unit number, will dump info only if present
  !---------------------------------------------------------------------------
  Real(kind=RDbl) Function interact_checknrgs(imodel,fast,species,simcell, &
      indent,unitno)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Logical, Intent(In)                            :: fast
    Type(AtMolCoords), Dimension(:), Intent(In)    :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Integer, Intent(In)                            :: indent
    Integer, Intent(In), Optional                  :: unitno

    Integer                         :: i,j,nspc,unit
    Logical                         :: success,prob
    Real(kind = RDbl)               :: total_pot,diffsum,diff,tolerance
    Character(len=strlen)           :: names,string
    Character(len=20)               :: cutnames
    Character(len=indent)           :: blank
    Type(Forcefield_Results)        :: temp
    Type(EnergyPlus), Dimension(MAX_SORBS)             :: intra1,intra2
    Real(kind = RDbl), Dimension(MAX_SORBS,MAX_SORBS)  :: coul1,noncoul1
    Real(kind = RDbl), Dimension(MAX_SORBS,MAX_SORBS)  :: coul2,noncoul2

    blank = Repeat(' ',indent)
    nspc = molecules_getnsorbs()

    unit = 6
    If (Present(unitno)) unit = unitno

    tolerance = 1.0e-3_RDbl
    !    tolerance = 1.0e-6_RDbl
    interact_checknrgs = 0.0_RDbl

    !** Copy the full system storage structure
    Call storetop_initcopy(temp,imodel%results(1),(/0,0,0/))

    !** Zero the storage structure
    Call storetop_zero(temp,.False.)

    !** Call forcefield to get all the interactions
    success = forcefield_allint(imodel%ff(1),temp,species,simcell,fast,.False.)
    Call checkandstop(success,__FILE__,__LINE__, &
        ' Unexpected problem with full system forcefield calculation')

    !** Make sure the spc-spc level sums in the storage are internally consistent
    Call storetop_fastsum(temp)  

    !** Allocate the intramolecular energies
    Do i = 1,MAX_SORBS
      Call storebase_init(intra1(i),imodel%nderivs,.True.)
      Call storebase_init(intra2(i),imodel%nderivs,.True.)
    End Do

    !** Get the energies from the storage structure
    Call storetop_allnrgs(temp,noncoul2,coul2,intra2)

    !** Scale the energies from the storage
    Do i = 1,nspc
      Do j = 1,nspc
        noncoul2(i,j) = noncoul2(i,j)*calToJ
        coul2(i,j) = coul2(i,j)*calToJ
      End Do

      intra2(i)%intranrg = intra2(i)%intranrg*calToJ
      Call storebase_sumintra(intra2(i),.False.)
    End Do

    !** Compare the storage energies to the spc-spc statistics energies
    diffsum = 0.0_RDbl
    prob = .False.
    Do i = 1,nspc
      Call storestats_getallintra(imodel%spcstats,i,'inst',intra1(i))
      Call storebase_sumintra(intra1(i),.False.)

      Do j = 1,NO_OF_INTRA_POTS
        diff = intra1(i)%intranrg(j) - intra2(i)%intranrg(j)
        If (Abs(diff) > tolerance) Write(unit,*) blank,'PROBLEM: Intra ',i,j,diff
        If (Abs(diff) > tolerance) prob = .True.
        diffsum = diffsum + diff**2
      End Do

      Do j = i,nspc
        noncoul1(i,j) = storestats_getnoncoul(imodel%spcstats,i,j,'inst')
        diff = noncoul1(i,j) - noncoul2(i,j)
        If (Abs(diff) > tolerance) Write(unit,'(2a,2i3,3e14.4)') blank, &
            'PROBLEM: Non-Coul ',i,j,diff,noncoul1(i,j),noncoul2(i,j)
        If (Abs(diff) > tolerance) prob = .True.
        diffsum = diffsum + diff**2

        coul1(i,j) = storestats_getcoul(imodel%spcstats,i,j,'inst')
        diff = coul1(i,j) - coul2(i,j)
        If (Abs(diff) > tolerance) Write(unit,'(2a,2i3,3e14.4)') blank, &
            'PROBLEM: Coul ',i,j,diff,coul1(i,j),coul2(i,j)
        If (Abs(diff) > tolerance) prob = .True.
        diffsum = diffsum + diff**2
      End Do
    End Do

    interact_checknrgs = Sqrt(diffsum)

    !** Write to display if requested
    If (Present(unitno)) Then
      Write(unitno,'(2a)') blank,dashedline
      Write(unitno,'(2a)') blank,'Comparing Stored and Freshly&
          & Calculated Energies'
      Write(unitno,'(2a,5x,a22,3x,a15,2x,a)') blank,'Species','Energy', &
          'FRESH (kJ/mol)','STORED (kJ/mol)'

      !** Display Non-Coulombic and Coulombic information
      Do i = 1,nspc
        Do j = i,nspc
          Write(names,'(3a)') Trim(molecules_name(i)),'-',Trim(molecules_name(j))
          cutnames = names(1:20)
          Write(unitno,'(2x,2a,1x,a14,2x,e12.4,3x,e12.4)') blank,cutnames, &
              'Coulombic',coul2(i,j),coul1(i,j)
          Write(unitno,'(2x,2a,1x,a14,2x,e12.4,3x,e12.4)') blank,cutnames, &
              'Non-coulombic',noncoul2(i,j),noncoul1(i,j)
        End Do
      End Do

      !** Display intramolecular information
      Do i = 1,nspc
        Write(cutnames,'(a20)') molecules_name(i)
        Write(unitno,'(2x,2a,1x,a14,2x,e12.4,3x,e12.4)') blank,cutnames, &
            'Intramolecular', &
            intra2(i)%intranrg(TOTAL_INDEX),intra1(i)%intranrg(TOTAL_INDEX)
      End Do

      If (prob) Then
        string = real2str(interact_checknrgs,8)
        Write(*,'(3a)') blank,"Deviation between stored and newly &
            &calculated energies: ", Trim(string)        
        Call interact_chkrecalc(imodel,species,simcell,tolerance,indent,unit)
      End If

    End If

    !** Deallocate temporary storage and the intramolecular energies
    Call storetop_clean(temp)
    Call storebase_clean(intra1)
    Call storebase_clean(intra2)

    !** Stop for any errors that interact_chkrecalc may have missed
    If (prob) Then
      Write(*,'(3a)') blank,"Stopped due to interaction evaluation ", &
          "reproducibility problems"
      Stop
    End If

  End Function interact_checknrgs

  !---------------------------------------------------------------------------
  ! Copy the full system interaction structure and recalculate the full 
  ! system interactions.  Then compare these interactions in depth to the
  ! current full system interaction storage.  This routine can be used to
  ! fully check the partial system updates.
  ! Requires:  imodel -- interaction model information
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            tolerance -- tolerance for comparing old and new
  !            indent -- no. of spaces from the left margin
  !            unit -- optional unit number
  !---------------------------------------------------------------------------
  Subroutine interact_chkrecalc(imodel,species,simcell,tolerance,indent,unitno)
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(In)    :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Real(kind=RDbl), Intent(In)                    :: tolerance
    Integer, Intent(In)                            :: indent
    Integer, Intent(In), Optional                  :: unitno

    Integer                      :: unit
    Logical                      :: fast,success,equal
    Character(len=indent)        :: blank
    Character(len=strLen)        :: string
    Type(Forcefield_Results)     :: temp

    blank = Repeat(' ',indent)
    unit = 6
    If (Present(unitno)) unit = unitno

    !** Copy the full system storage structure
    Call storetop_initcopy(temp,imodel%results(1),(/0,0,0/))

    !** Zero the storage structure
    Call storetop_zero(temp,.False.)

    !** Call forcefield to get all the interactions
    fast = .True.
    success = forcefield_allint(imodel%ff(1),temp,species,simcell,fast,.False.)
    Call checkandstop(success,__FILE__,__LINE__, &
        ' problem with full system forcefield calculation')

    !** Compare the old and temporary forcefield storage structures
    !** reported differences are (current - recalculated)
    If (Present(unitno)) Then
      equal = storetop_chkequal(imodel%results(1),temp,tolerance,indent,unit)
    Else
      equal = storetop_chkequal(imodel%results(1),temp,tolerance)
    End If
    If (.Not. equal) Then
      Write(*,'(1x,2a,i4,a)') __FILE__,' : ',__LINE__, &
          ' ERROR Recalculation reproducibility test failed'
      Write(*,'(/a)') 'CORRECT interaction storage structure:'
      Call storetop_display(temp,.True.,2,6)
      Write(*,'(/a)') 'ACTUAL interaction storage structure:'
      Call storetop_display(imodel%results(1),.True.,2,6)
      Write(*,'(a)') 'STOPPED due to errors in recalculation reproducibility'
      Stop    
    End If

    !** Clean up the temporary storage structure
    Call storetop_clean(temp)

  End Subroutine interact_chkrecalc

  !----------------------------------------------------------------------------
  ! Returns information about the forcefield parameters governing the 
  ! interactions between two subsets of the system.  This information can
  ! then be used elsewhere.
  ! Requires:  imodel -- interaction model information
  !            subset1 -- 1st subset (specify only this for intra-only)
  !            subset2 -- 2nd subset
  !            intratype -- intramolecular interaction (uses usual keywords)
  !            nsets -- number of potential sets returned
  !            list -- list of atom numbers for interaction (set,1:Natoms)
  !            params -- array of strings containing type and parameters
  ! Currently works only with intramolecular potentials, needs to be 
  ! generalized to return any forcefield parameters!
  !----------------------------------------------------------------------------
  Subroutine interact_listparams(imodel,subset1,subset2,intratype,&
      nsets,list,params)
    Type(Interaction_Model), Intent(In)               :: imodel
    Integer, Dimension(:), Intent(In)                 :: subset1,subset2
    Character(*), Intent(In)                          :: intratype
    Integer, Intent(Out)                              :: nsets
    Integer, Dimension(:,:), Intent(Out)              :: list
    Character(len=lstrLen), Dimension(:), Intent(Out) :: params

    Call forcefield_listparams(imodel%ff(1),subset1,subset2,intratype, &
        nsets,list,params)

  End Subroutine interact_listparams


  !----------------------------------------------------------------------------
  ! sets forcefield parameters. this is required to change some
  ! parameter during a simulation. For example if we want to continuously
  ! change the strength of bond
  ! Requires:  imodel -- interaction model information
  !            subset -- speciefies spc numbers
  !            description -- bunchof strings
  !            alist -- lsit of atoms (needed only for intra)
  !            val -- the value
  !----------------------------------------------------------------------------
  Subroutine interact_setparams(imodel,subset,description, alist,val)
    Type(Interaction_Model), Intent(InOut)               :: imodel
    Integer, Dimension(:), Intent(In)                 :: subset
    Character(len=strlen),Dimension(:), Intent(In)    :: description
    Integer, Dimension(:), Intent(in)                 :: alist
    Real(kind=RDbl), Intent(in) :: val

    Call forcefield_setparams(imodel%ff(1),subset,description, alist,val)
  End Subroutine interact_setparams


  !------------------------------------------------------------------------------
  ! Tina added this, should make an interface with _listparams.
  ! Subroutine to get potential parameters, the array pot_params contains
  ! A, B, C, D, hicut, locut and returns zero if the parameters are not
  ! defined for the potential 
  !------------------------------------------------------------------------------
  Subroutine interact_getpotparameters(params, spc1,spc2,a1, a2, pot_params)
    Type(interaction_model), Intent(In)           :: params
    Integer, Intent(In)                           :: spc1, spc2
    Integer, Intent(In)                           :: a1, a2
    Real(kind = RDbl), Dimension(6), Intent(Out)  :: pot_params

    Call forcefield_getpotparameters(params%ff(1),spc1,spc2, a1, a2,pot_params)

  End Subroutine interact_getpotparameters

  !----------------------------------------------------------------------------
  ! Returns the box increments for a potential map if they are available.
  ! Requires:  imodel -- interaction model information
  !            species -- species data structure
  !            boxincrements -- x,y,z increments for potential map
  !            boxsteps -- number of boxes in x,y,z directions
  !----------------------------------------------------------------------------
  Subroutine interact_boxinfo(imodel,species,boxincrements,boxsteps)
    Type(Interaction_Model), Intent(In)           :: imodel
    Type(AtMolCoords), Dimension(:), Intent(In)   :: species
    Real(Kind=RDbl), Dimension(3), Intent(Out)    :: boxincrements
    Integer, Dimension(3), Intent(Out)            :: boxsteps

    Call forcefield_boxinfo(imodel%ff(1),species,boxincrements,boxsteps)

  End Subroutine interact_boxinfo

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to unit unitno
  ! Requires:  unitno -- unit number to write to
  !----------------------------------------------------------------------------
  Subroutine interact_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) '_sampleCF not finished'
    Stop

    Write(unitno,'(a)') &
        '# If INTEGRATE BOOST was specified, remove this line and include:'
    Call boost_sampleCF(unitno)
    Write(unitno,'(a)') '# End INTEGRATE BOOST section'

  End Subroutine interact_sampleCF

  !---------------------------------------------------------------------
  ! Handles the display of any dynamic quantities that might be useful
  ! to see during a simulation.
  ! Requires:  imodel -- interaction model structure
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !---------------------------------------------------------------------
  Subroutine interact_dispDynInfo(imodel,indent,unitno)
    Type(Interaction_Model), Intent(In)  :: imodel
    Integer, Intent(In)                  :: indent
    Integer, Optional, Intent(In)        :: unitno

    Integer                           :: unit
    Character(len=indent)             :: blank

    blank = Repeat(' ',indent)

    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    !** Display boost stuff
    If (Associated(imodel%boost)) Then
      Call boost_displayStats(imodel%boost,unit,1)
    End If

  End Subroutine interact_dispDynInfo

  !---------------------------------------------------------------------------
  ! Check the correspondence between the number of molecules of each species
  ! in the system with the numbers that the storage structure thinks it
  ! contains.  This is only used for debugging.
  ! Requires:  imodel -- interaction model structure
  !            species -- species data structure
  !            indent -- no. of spaces from the left margin
  !            unit -- unit number
  !---------------------------------------------------------------------------
  Subroutine interact_chkmolnums(imodel,species,indent,unitno)
    Type(Interaction_Model), Intent(In)         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: indent
    Integer, Intent(In)                         :: unitno

    Integer                  :: m,nspcs,spc,instorage
    Logical                  :: failed
    Character(len=indent)    :: blank
    Character(len=strLen)    :: string1,string2,string3
    Integer, Dimension(100)  :: list

    blank = Repeat(' ',indent)

    nspcs = molecules_getnsorbs()
    list = 0
    Call config_getnmoleslist(species,list)

    failed = .False.
    Do spc = 1,nspcs
      instorage = imodel%results(1)%nmoles_system(spc)
      If (instorage /= list(spc)) Then
        string1 = int2str(spc)
        string2 = int2str(instorage)
        string3 = int2str(list(spc))
        Write(unitno,'(7a)') blank,'nmoles agreement problem for species ', &
            Trim(string1),'   instorage, inconfig: ',Trim(string2), &
            ' ',Trim(string3)
        failed = .True.
      End If
    End Do

    If (failed) Then
      string1 = int2str(imodel%results(1)%nmoles_system(1:nspcs))
      string2 = int2str(list(1:nspcs))
      Write(unitno,'(3a)') blank,'in storage: ',Trim(string1)
      Write(unitno,'(3a)') blank,'in config:  ',Trim(string2)
    End If

  End Subroutine interact_chkmolnums

  !---------------------------------------------------------------------
  ! Displays the interaction structure
  ! Requires:  imodel -- interaction model structure
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !---------------------------------------------------------------------
  Subroutine interact_display(imodel,indent,unitno)
    Type(Interaction_Model), Intent(In)  :: imodel
    Integer, Intent(In)                  :: indent
    Integer, Optional, Intent(In)        :: unitno

    Integer                   :: i,unit
    Logical                   :: flag
    Character(len=indent)     :: blank,string

    blank = Repeat(' ',indent)

    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    Write(unitno,'(2a)') blank,dashedline
    Write(unitno,'(2a)') blank,'INTERACTION MODEL Information:'
    string = int2str(imodel%nderivs)
    Write(unitno,'(2x,3a)') blank,'number of derivatives: ',Trim(string)
    string = int2str(imodel%nff)
    Write(unitno,'(2x,3a)') blank,'number of forcefields: ',Trim(string)
    flag = Associated(imodel%boost)
    Write(unitno,'(2x,2a,l2)') blank,'boosting: ',flag
    flag = Associated(imodel%usampling)
    Write(unitno,'(2x,2a,l2)') blank,'umbrella sampling: ',flag

    !** Display the full details of the forcefield(s)
    Do i = 1,imodel%nff
      Call forcefield_display(imodel%ff(i),indent,unit)
    End Do

    !** Display the boost stuff if necessary
    If (Associated(imodel%boost)) Then
      Call boost_display(imodel%boost,unit,indent+2)
    End If

    !** Display the umbrella sampling stuff if necessary
    If (Associated(imodel%usampling)) Then
      Call umbrella_display(imodel%usampling,unit,indent+2)
    End If

  End Subroutine interact_display

  !---------------------------------------------------------------------
  ! Cleans the interation model structure
  ! Requires:  imodel -- interaction model structure
  !---------------------------------------------------------------------
  Subroutine interact_clean(imodel)
    Type(Interaction_Model), Intent(InOut)       :: imodel

    Integer            :: i,error

    Call storestats_clean(imodel%spcstats)

    If (Associated(imodel%boost)) Then
      Deallocate(imodel%boost, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'boost')    

      Do i = 1,Size(imodel%intra)
        Call storebase_clean(imodel%intra(i))
      End Do

      Deallocate(imodel%intra, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'intra')    
      Deallocate(imodel%noncoul, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'noncoul')    
      Deallocate(imodel%coul, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'coul')    
      Deallocate(imodel%pscalings, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'pscalings')    
    End If

    If (Associated(imodel%usampling)) Then
      Call umbrella_clean(imodel%usampling)
      Deallocate(imodel%usampling, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'usampling')
    End If

    Do i = 1,imodel%nff
      Call forcefield_clean(imodel%ff(i))
      Call storetop_clean(imodel%results(i))
    End Do

    Deallocate(imodel%ff, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'ff')
    Deallocate(imodel%results, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'results')

  End Subroutine interact_clean

End module interact





