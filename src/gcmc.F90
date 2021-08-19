!-----------------------------------------------------------------------------
! This module contains the main data structure and routines for doing 
! Grand Canonical Monte Carlo (GCMC) simulations.  The routines provided
! here are meant to be used by a driver that loops through the different
! simulation numbers, breaks them into chunks and calls the display routines.
!
! The selection of a species is handled at this level.  All species-specific
! move handling and information storage is handled by "gcmcmoves" which, in
! turn, uses the "moves" module to handle specific moves.
!
! Important routines are:
!    gcmc_init -- initializes the information from a control file
!    gcmc_dosim -- does a specified number of interations 
!    gcmc_displaystats -- displays intermediate statistics for simulation
!
! Needed Improvements:
! 1) remove dangerous hardcoding of 'fast' parameter
! 2) testing for more than one species
!-----------------------------------------------------------------------------

Module gcmc
  Use auxmoveparams, Only: AuxMoveObjects, auxmoveparams_init
  Use cavitylist, Only : cavitylist_displayCubeInfo, &
      cavitylist_display, cavitylist_displaySorbInfo
  Use defaults, Only: RDbl, strLen, dashedline,dashedline2, twopi, d_res_file, &
      d_con_file, hplanck, MAX_SORBS, NO_OF_INTRA_POTS, STRETCH_INDEX, &
      BENDING_INDEX,TORSION_INDEX,CONSTRAINT_INDEX,TOTAL_INDEX,Rgas, Nav, &
      scalepe, kcalmole_kb, lstrLen, dbgflag
  Use utils, Only: filesrchstr, stripcmnt, toint, toupper, int2str, real2str, &
      allocErrDisplay, deallocErrDisplay, readblank, checkandstop, split
  Use file, Only: file_getunit,file_gettype,file_open
  Use general, Only: genparams
  Use random, Only: rranf
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use eos, Only: eos_Models,eos_init,eos_settemperature,eos_setpressure,&
      eos_getfugacity,  eos_getConfInteg
  Use config, Only: AtMolCoords, config_writerestartfile, config_getnmoles, &
      config_display, config_isfixed
  Use simcell, Only: SimCell_Params, simcell_getvolume
  Use gcmcmoves, Only : GCMC_Move_Params,gcmcmoves_init,gcmcmoves_move, &
      gcmcmoves_displaystats, gcmcmoves_display, gcmcmoves_displaysimparams, &
      gcmcmoves_dispspcstats, gcmcmoves_updatenmolecs, gcmcmoves_beginsim, &
      gcmcmoves_copyfuglist, Fugacity_Params
  Use molecule, Only: MolecularParams
  Use molecules, Only: molecules_getmass, molecules_getnsorbs
  Use storestats, Only: Species_Stats, storestats_getcoul, &
      storestats_getnoncoul, storestats_displaynrg, storestats_displayavgs
  Use stats, Only: stats_update,stats_getvalue,stats_reset
  Use interact, Only: Interaction_Model,interact_checknrgs, &
      interact_initnrgs, interact_chkrecalc
  Use subinteract, Only: Subset_Interactions, subinteract_int, &
      subinteract_init, subinteract_oldnrg, &
      subinteract_newnrg, subinteract_chksubint
 
  Implicit None
  Save

  Private
  Public :: gcmc_init,gcmc_initdisplay,GCMC_Params,gcmc_dosim,&
      gcmc_getnosims,gcmc_displaystats,gcmc_displaysimparams, &
      gcmc_beginsim,gcmc_endsim, gcmc_temperature, gcmc_pressure, &
      gcmc_sampleCF,gcmc_chknrgs

  !** The tag marking the beginning of the GCMC section
  Character(len=strLen), Parameter    :: default_gcmc_tag = &
      "GCMC Information"

  !** Contains all the GCMC information for the whole simulation
  Type GCMC_Params
    Integer                :: nspc
    Integer                :: nsims        !* no. of simulations
    Integer                :: niterations
    Real(kind=RDbl)        :: tk           !* temperature in Kelvin
    Character(len=strLen)  :: eostag
    Type(EOS_Models)       :: eosparams

    !** gcmcspc are NOT indexed by species type.  They go from 1->nspc
    Type(GCMC_Move_Params), Dimension(:), Pointer   :: gcmcspc 

    !** Subset interaction storage, species-dependent, indexed by species number
    Type(Subset_Interactions), Dimension(:), Pointer :: subints

    Type(AuxMoveObjects),Pointer :: auxmv
    
    ! contains pressure list for all spc for all sims
    Type(Fugacity_Params), Dimension(:,:), Pointer   :: fuglist

  End Type GCMC_Params

Contains
  !-------------------------------------------------------------------
  ! Initializes the various GCMC parameters from the control file
  ! Requires: gcmcparams -- parameters for the GCMC simulation
  !           imodel -- interaction model information
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           ctrl_filename -- the name of the control file
  !           opt_gcmctag -- optional Tag to look for in control file
  !-------------------------------------------------------------------
  Subroutine gcmc_init(gcmcparams, imodel, species, simcell, ctrl_filename, &
      opt_gcmctag)
    Type(GCMC_Params), Intent(InOut)               :: gcmcparams
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Character(*), Intent(In)                       :: ctrl_filename
    Character(*), Optional, Intent(In)             :: opt_gcmctag

    Integer                :: unitno, gcmclineno, nspc, error, i, j,blocksize
    Integer                :: indent, spc
    Logical                :: fast
    Real(kind=RDbl)        :: volume
    Character(len=strLen)  :: eostag, tag, line, itype
    Character(len=lstrLen) :: errormsg

    If (Present(opt_gcmctag)) Then
      tag = opt_gcmctag
    Else
      tag = default_gcmc_tag
    End If

    !** Initialize the auxmv field. This contains other gcmc stuff 
    !** like cavity params, may not be required always. this will change 
    !** your position in ctrlfile
    Call auxmoveparams_init(gcmcparams%auxmv, species, simcell, ctrl_filename)
!    Call cavitylist_display(gcmcparams%auxmv%cavity, 5, 6)

    !** Open the ctrl_file if it is not opened
    unitno = file_open(ctrl_filename,110)
    Rewind(unitno)
    
    !** Find the GCMC section
    gcmclineno = filesrchstr(unitno, tag, line)
    If (gcmclineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", Trim(tag), " in the control file"
      Stop
    Endif
    
    Read(unitno, *) gcmcparams%niterations
    Read(unitno, *) gcmcparams%tk
    Read(unitno, '(a)') eostag
    gcmcparams%eostag = stripcmnt(eostag)
    Read(unitno, *) gcmcparams%nsims
    Read(unitno, *) blocksize
    Read(unitno, *) nspc
    gcmcparams%nspc = nspc !** number of species

    !** Warn about multiple species if necessary
    If (gcmcparams%nspc > 1) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' WARNING: GCMC has not been tested for multiple species'
    End If
    
    !** read the blankline at the end
    errormsg = 'must be a blank line before species parameters in GCMC section'
    Call readblank(unitno,__FILE__,__LINE__,errormsg)
    
    !** Allocate the gcmcparams%gcmcspc 
    Allocate(gcmcparams%gcmcspc(nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
    !** Set the GCMC move parameter information that should be the same
    !** for all moveable species such as temperature and no. of pts
    Do i = 1,nspc
      gcmcparams%gcmcspc(i)%tk   = gcmcparams%tk
      gcmcparams%gcmcspc(i)%npts = gcmcparams%nsims
      gcmcparams%gcmcspc(i)%blocksize = blocksize
    End Do

    !** Initialize the subset interactions for each moving species
    nspc = molecules_getnsorbs()   !** note changed nspc definition
    Allocate(gcmcparams%subints(nspc), STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'subints')
    Do spc = 1,nspc
      itype = 'MC'
      If (config_isfixed(species(spc))) itype = 'NULL'
      Call subinteract_init(gcmcparams%subints(spc),imodel,'Molec_System', &
          itype,(/spc,1,0/),(/0,0,0/))
    End Do
    
    !** The rest of the stuff will be read and initialized by the 
    !** gcmcmove init
    Call gcmcmoves_init(gcmcparams%gcmcspc, species, simcell, ctrl_filename, &
       gcmcparams%auxmv)
!    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!    Write(*,*)   gcmcparams%auxmv%cavity%ncubes, Size(gcmcparams%auxmv%cavity%sorbs)

    !** Initialize the parameters for the chosen equation of state 
    Rewind(unitno)
    Call eos_init(gcmcparams%eosparams, ctrl_filename, gcmcparams%eostag, &
        gcmclineno)
    
    !** Generate the fugacitylist 
    volume = simcell_getvolume(simcell)
    Call gcmc_genfuglist(gcmcparams, volume)
    
    Close(unit=unitno)
    
    !** Update all the stored interactions
    fast = .True.
    Call interact_initnrgs(imodel,fast,species,simcell,.True.)

    Write(*,*) "Energies at the beginning of GCMC simulation: "
    Call storestats_displaynrg(imodel%spcstats,3,6)

  End Subroutine gcmc_init
 
  !---------------------------------------------------------------------------
  ! Performs "ninterations" of the specified GCMC simulation.
  ! We need to ensure that insertions and deletions are done with equal 
  ! probability for each species type. The algorithm used here is 
  ! similar to the one used by Karavias and Myers (Mol. Sim. (1991) v8, p.51)
  ! 1) Pick a species with equal probability.
  ! 2) Decide whether we want to attempt an
  !    insertion, deletion, translation or rotation.
  ! Enforcement of these constraints are handled in gcmcmoves (init)
  ! Requires: gcmcparams -- parameters for the GCMC simulation
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           simno -- the simulation number
  !---------------------------------------------------------------------------
  Subroutine gcmc_dosim(gcmcparams,species,simcell,simno)
    Type(GCMC_Params), Intent(InOut)                 :: gcmcparams
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Type(SimCell_Params), Intent(In)                 :: simcell
    Integer, Intent(In)                              :: simno

    Integer               :: iter,nmoles,spc
    Logical               :: success

    !** Make sure that the structure has been initialized
    If (.Not. Associated(gcmcparams%gcmcspc)) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " GCMC Species must be initialized before starting simulation"
      Stop
    End If

    Do iter = 1,gcmcparams%niterations
      !** Pick a species randomly
      spc = Int(rranf()*gcmcparams%nspc) + 1

      !** Pick a move type randomly and execute it
      success = gcmcmoves_move(gcmcparams%gcmcspc(spc),gcmcparams%subints, &
          species,simcell,simno)

      !** Do updating of quantities to be averaged
      Call gcmcmoves_updatenmolecs(gcmcparams%gcmcspc,species)
    End Do
#ifdef DEBUG
    If (dbgflag) Then
      Call cavitylist_displaycubeInfo(gcmcparams%auxmv%cavity, "cube.xyz", "cube", "Ne ")
      Call cavitylist_displaySorbInfo(gcmcparams%auxmv%cavity, species, "sorb.xyz", "sorb", "Ar ")
      
    Endif
#endif
  End Subroutine gcmc_dosim

  !------------------------------------------------------------------------
  ! Fill in the rest of the fields in the fugacity list
  ! Requires: gcmcparams -- GCMC simulation parameters
  !           volume -- volume of the simulation cell
  ! NOTE: this routine should not directly access gcmcspc params, rewrite
  !------------------------------------------------------------------------
  Subroutine gcmc_genfuglist(gcmcparams, volume)
    Type(GCMC_Params), Intent(InOut)  :: gcmcparams
    Real(kind=RDbl), Intent(In)       :: volume

    Integer        :: i, j, spc, nsims, nspc, error
    Real(kind=RDbl):: tk, pp, fug, B, sivolume, mass, murti, Lambda, ratio, Zig

    !** Set the temperature and the pressure of the different
    !** species
    tk = gcmcparams%tk
    Call eos_settemperature(gcmcparams%eosparams, tk)

    nsims= gcmcparams%nsims
    nspc=molecules_getnsorbs() ! total number, NOT size of gcmcspc

    !** array for fugacities
    Allocate(gcmcparams%fuglist(nsims, nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'fuglist')

    !** Get the fugacity
    Do i=1, nsims

      !** Set the pressure of each component
      Do j=1, gcmcparams%nspc
        spc = gcmcparams%gcmcspc(j)%spc

        !** pp is in kPa
        pp = gcmcparams%gcmcspc(j)%fuglist(i)%pressure
!        Write(*,*) spc, pp

        Call eos_setpressure(gcmcparams%eosparams, spc, pp)
      End Do

      !** Get the fugacity of each component
      Do j=1, gcmcparams%nspc
        spc = gcmcparams%gcmcspc(j)%spc
        fug = eos_getfugacity(gcmcparams%eosparams, spc) ! fug [=] kPa
        gcmcparams%gcmcspc(j)%fuglist(i)%fugacity = fug

        !** Get the excess chemical potential (B in Adams Notation)
        sivolume = volume*1.0e-30*Nav ! convert to m^3/mole
        B = Log(fug*1.0e3*sivolume/(Rgas*tk))

       !** This is Z/Omega, got from ideal parameters
       Zig = eos_getConfInteg(gcmcparams%eosparams, j)
!       Write(*,*) B, Zig
       ratio = log(Zig) 

       !** we need (PV)/(RTZ)
       B = B - ratio
!        Write(0,'(2a,i4,a,f16.2,6f16.10)') __FILE__,": ",__LINE__, &
!        " B ", ratio,B,Zig

        gcmcparams%gcmcspc(j)%fuglist(i)%B = B 

        !** Get mu the chemical potential
        !** Calculate the DeBroglie wavelength
        mass = molecules_getmass(spc)
        mass = mass*1.0e-3              ! convert to kg
        Lambda = hplanck/Sqrt(twopi*mass/Nav*Rgas/Nav*tk)
        murti = (B - Log(sivolume/Nav/Lambda**3))
        gcmcparams%gcmcspc(j)%fuglist(i)%murti = murti
!        Write(*,*) fug, B, ratio, murti
      End Do
    End Do

    !** copy the values to local array also
    Do i=1,nsims
      Do j=1, gcmcparams%nspc
        spc = gcmcparams%gcmcspc(j)%spc
        gcmcparams%fuglist(i,spc)%pressure = &
            gcmcparams%gcmcspc(j)%fuglist(i)%pressure
        gcmcparams%fuglist(i,spc)%fugacity = &
            gcmcparams%gcmcspc(j)%fuglist(i)%fugacity
        gcmcparams%fuglist(i,spc)%murti = &
            gcmcparams%gcmcspc(j)%fuglist(i)%murti
        gcmcparams%fuglist(i,spc)%B = &
            gcmcparams%gcmcspc(j)%fuglist(i)%B
      End Do
    End Do

  End Subroutine gcmc_genfuglist

  !--------------------------------------------
  ! Gets the total no. of simulations
  !--------------------------------------------
  Integer Function gcmc_getnosims(gcmcparams)
    Type(GCMC_Params), Intent(In)  :: gcmcparams

    gcmc_getnosims = gcmcparams%nsims

  End Function gcmc_getnosims

  !--------------------------------------------------------------------------
  ! Check the energy calculation and updating system that the GCMC uses.  If
  ! unitno is specified, it will produce output.  Otherwise, output will 
  ! only be produced if there is an error.
  ! Requires:  gcmcparams -- GCMC simulation parameters
  !            imodel -- interaction model information
  !            species -- species data structure  
  !            simcell -- the simulation cell information
  !            simno -- simulation number
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional unit number for display
  !--------------------------------------------------------------------------
  Subroutine gcmc_chknrgs(gcmcparams,imodel,species,simcell,indent,unitno)
    Type(GCMC_Params), Intent(InOut)               :: gcmcparams
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Integer, Intent(In)                            :: indent
    Integer, Intent(In), Optional                  :: unitno

    Integer                      :: spc
    Logical                      :: dump
    Real(kind=RDbl)              :: tolerance
    Character(len=indent)        :: blank

    tolerance = 1.0e-6_RDbl
    blank = Repeat(' ',indent)
    dump = .False.
    If (Present(unitno)) dump = .True.

    !** Check reproducibility of full system forcefield calculation
    !** This recalculates the full interaction and compares them in 
    !** detail to the stored interactions, which come from updating
    If (dump) Then
      Write(unitno,'(a)') 'Checking recalculation reproducibility'
    End If
    If (dump) Then
      Call interact_chkrecalc(imodel,species,simcell,tolerance,indent,unitno)
    Else
      Call interact_chkrecalc(imodel,species,simcell,tolerance,indent)
    End If

    !** Check interact's consistency in molec-molec interactions
    If (dump) Then
      Write(unitno,'(a)') 'Checking extraction/recalculation reproducibility'
    End If
    Do spc = 1,molecules_getnsorbs()
      If (config_isfixed(species(spc))) Cycle
      If (dump) Then
        Call subinteract_chksubint(gcmcparams%subints(spc),spc,species, &
            simcell,indent,unitno)
      Else
        Call subinteract_chksubint(gcmcparams%subints(spc),spc,species, &
            simcell,indent)
      End If
    End Do

  End Subroutine gcmc_chknrgs

  !--------------------------------------------
  ! Gets the simulation temperature
  !--------------------------------------------
  Real(kind=RDbl) Function gcmc_temperature(gcmcparams, simno)
    Type(GCMC_Params), Intent(In)  :: gcmcparams
    Integer, Intent(in) :: simno ! as of now irrelevant
    gcmc_temperature = gcmcparams%tk
  End Function gcmc_temperature

  !--------------------------------------------
  ! Gets the pressure of first gcmc-species
  !--------------------------------------------
  Real(kind=RDbl) Function gcmc_pressure(gcmcparams, simno)
    Type(GCMC_Params), Intent(In)  :: gcmcparams
    Integer, Intent(in) :: simno ! as of now irrelevant
    gcmc_pressure = gcmcparams%gcmcspc(1)%fuglist(simno)%pressure
  End Function gcmc_pressure
  
  !--------------------------------------------------------------------------
  ! Display the GCMC statistics for the current simulation
  ! Requires:  gcmcparams -- GCMC simulation parameters
  !            imodel -- interaction model information
  !            species -- species data structure  
  !            simcell -- the simulation cell information
  !            simno -- simulation number
  !            indent -- no. of spaces from the left margin
  !            optunit -- optional unit number for display
  !--------------------------------------------------------------------------
  Subroutine gcmc_displaystats(gcmcparams,imodel,species,simcell,simno,indent, &
      optunit)
    Type(GCMC_Params), Intent(InOut)               :: gcmcparams
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Integer, Intent(In)                            :: simno,indent
    Integer, Intent(In), Optional                  :: optunit

    Integer                      :: unitno,spc
    Logical                      :: fast
    Real(kind=RDbl)              :: nrg_devn
    Character(len=indent)        :: blank
    Character(len=strLen)        :: string

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    !** Write the simulation no.
    Write(unitno, '(2a)') blank, dashedline
    Write(unitno, '(2a)') blank, "The GCMC Stats:"
    Call gcmcmoves_displaystats(gcmcparams%gcmcspc,simno,indent+2,unitno)

    !** Check storage versus fresh, full-system evaluation of interactions
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      !** the stats values
      Write(unitno,'(2x,2a)') blank, &
          "----------- Nrg stats (CAVG) Values-----------"
      string = "CAVG"
      Call storestats_displayavgs(imodel%spcstats,string,indent+2,unitno)

      fast = .True.
      nrg_devn = interact_checknrgs(imodel,fast,species,simcell,indent+2,unitno)
      string = real2str(nrg_devn,6)
      Write(unitno,'(4x,3a)') blank,"Deviation between stored and newly &
          &calculated energies: ", Trim(string)
    End If

  End Subroutine gcmc_displaystats

  !------------------------------------------------------------------
  ! Display the gcmc simulation parameters for a specific simulation.
  ! Useful for summarizing parameters at the beginning of a new sim
  ! Requires: gcmcparams -- GCMC simulation parameters
  !           simno -- the simulation number
  !           indent -- no. of spaces from the left margin
  !           optunit -- optional unit number for display
  !-----------------------------------------------------------------
  Subroutine gcmc_displaysimparams(gcmcparams,simno,indent,optunit)
    Type(GCMC_Params), Intent(In) :: gcmcparams
    Integer, Intent(In)           :: simno,indent
    Integer, Optional, Intent(In) :: optunit

    Integer                 :: i, unitno, funit
    Character(len=indent)   :: blank
    Character(len=strLen)   :: restartfile, configfile

    blank = Repeat(' ',indent)    

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    !** Get the configuration and restart file names
    Call file_gettype(d_con_file,configfile,funit)
    Call file_gettype(d_res_file,restartfile,funit)

    !** Write the simulation no.
    Write(unitno, '(2a)') blank,"The GCMC Simulation Parameters:"
    Write(unitno, '(a,2x,a,i4)') blank,"Simulation Number : ", simno
    Write(unitno, '(a,2x,2a)') blank,"Configuration file: ", &
        Trim(configfile)
    Write(unitno, '(a,2x,2a)') blank, "Restart file      : ", &
        Trim(restartfile)
    Call gcmcmoves_displaysimparams(gcmcparams%gcmcspc,simno,indent+2,unitno)

  End Subroutine gcmc_displaysimparams

  !-------------------------------------------------------------------------
  ! Handles the beginning of GCMC simulations
  ! 1) Write feedback to screen
  ! Requires: gcmcparams -- GCMC simulation parameters
  !           simno -- simulation number
  !           species -- species data structure  
  !           simcell -- the simulation cell information
  !           indent -- no. of spaces from the left margin
  !           unit -- display unit number
  !-------------------------------------------------------------------------
  Subroutine gcmc_beginsim(gcmcparams,simno,species,simcell,indent,unit)
    Type(GCMC_Params), Intent(InOut)         :: gcmcparams
    Integer, Intent(In)                      :: simno
    Type(AtMolCoords), Dimension(:), Pointer :: species
    Type(SimCell_Params), Intent(In)         :: simcell
    Integer, Intent(In)                      :: indent,unit

    Integer                           :: spc

    !** Reset the counters and stats for each species
    Do spc = 1,gcmcparams%nspc
      Call gcmcmoves_beginsim(gcmcparams%gcmcspc(spc),simno)
    End Do
    
    !** copy the list of fugacities from gcmc to gcmcmoves_objects
    Do spc = 1,gcmcparams%nspc
      Call gcmcmoves_copyfuglist(gcmcparams%gcmcspc(spc),&
          gcmcparams%fuglist,simno)
    End Do

  End Subroutine gcmc_beginsim

  !-------------------------------------------------------------------------
  ! Handles the ending of GCMC simulations
  ! 1) Check the accumulated and fresh species-species and intra energies
  ! 2) Write feedback to screen
  ! 3) Optionally write to a restartfile (if stopfile /= '')
  ! 4) halt the program only if desired
  ! Step (4) can be useful for getting a starting file for subsequent MD runs.
  ! Requires:  gcmcparams -- GCMC simulation parameters
  !            simno -- simulation number
  !            imodel -- interaction model information
  !            species -- species data structure  
  !            simcell -- the simulation cell information
  !            stopfile -- name of restart file
  !            indent -- no. of spaces from the left margin
  !            unit -- display unit number
  !            stopflag -- optional flag to signal program halt
  !-------------------------------------------------------------------------
  Subroutine gcmc_endsim(gcmcparams,simno,imodel,species,simcell,stopfile, &
      indent,unit,stopflag)
    Type(GCMC_Params), Intent(In)            :: gcmcparams
    Integer, Intent(In)                      :: simno
    Type(Interaction_Model), Intent(InOut)   :: imodel
    Type(AtMolCoords), Dimension(:), Pointer :: species
    Type(SimCell_Params), Intent(In)         :: simcell
    Character(*), Intent(In)                 :: stopfile
    Integer, Intent(In)                      :: indent,unit
    Logical, Intent(In), Optional            :: stopflag 

    Integer                           :: spc
    Logical                           :: killflag,fast
    Real(kind = RDbl)                 :: rmsdev
    Character(len=indent)             :: blank
    Character(len=strLen)             :: string

    blank = Repeat(' ',indent)    

    If (Present(stopflag)) Then
      killflag = stopflag
    Else
      killflag = .False.
    End If

    !** check the updated and fresh energies for agreement
    fast = .True.
    rmsdev = interact_checknrgs(imodel,fast,species,simcell,0)
    If (rmsdev > 1.0E-2_RDbl) Then
      Write(0,'(1x,2a,i4,a)') __FILE__,": ",__LINE__, &
          ' ERROR: rms deviation between stored and'
      Write(0,'(5x,a,f8.3)') &
          'fresh energies is greater than 0.01 kJ/mol at ',rmsdev
      If (.Not.killflag) Stop
    Else If (rmsdev > 1.0E-6_RDbl) Then
      Write(0,'(1x,2a,i4,a,e14.4)') __FILE__,": ",__LINE__, &
          ' WARNING: rms deviation between stored and'
      Write(0,'(5x,a,e14.4)') 'fresh energies is significant: ',rmsdev
    End If

    !** Feedback 
    Write(unit,'(2a)') blank,dashedline
    string = int2str(simno)
    Write(unit,'(3a,4x,a,e9.3)') blank,'Ending GCMC simulation number ', &
        Trim(string),'Nrg dev check: ',rmsdev

    Do spc = 1,gcmcparams%nspc
      Call gcmcmoves_dispspcstats(gcmcparams%gcmcspc(spc),simno,indent+2,unit)
    End Do
    
    If (Trim(stopfile) /= '') Then
      Write(unit,*) "Writing a final restart file into  ", Trim(stopfile)
      Call config_writerestartfile(species, Trim(stopfile))
    End If

    Write(unit,'(a)') dashedline    
    
    !** Halts the program if desired
    If (killflag) Then
      Write(0,*) "Exiting after writing the restartfile"
      Stop
    Endif
    
  End Subroutine gcmc_endsim

  !-----------------------------------------------------------------------
  ! Display the gcmc initialization parameters, ie full information
  ! Requires: gcmcparams -- GCMC simulation parameters
  !           indent -- no. of spaces from the left margin
  !           optunit -- optional unit number for display
  !-----------------------------------------------------------------------
  Subroutine gcmc_initdisplay(gcmcparams,indent,optunit)
    Type(GCMC_Params), Intent(in)         :: gcmcparams
    Integer, Intent(in)                   :: indent
    Integer, Optional, Intent(in)         :: optunit

    Integer               :: i, unitno
    Character(len=indent) :: blank

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    Write(unitno,'(2a)') blank, dashedline
    Write(unitno,'(2a)') blank, "The GCMC Parameters Section:"
    Write(unitno,'(a,2x,a,i6)') blank,"No. of iterations  : ", &
        gcmcparams%niterations
    Write(unitno,'(a,2x,a,i6)') blank, "No. of simulations : ", &
        gcmcparams%nsims
    Write(unitno,'(a,2x,a,f8.3)') blank, "Temperature(K)     : ", &
        gcmcparams%tk
    Write(unitno,'(a,2x,a,i6)') blank, "No. of species     : ",gcmcparams%nspc
    Call gcmcmoves_display(gcmcparams%gcmcspc,indent+2,unitno)

  End Subroutine gcmc_initdisplay


  !----------------------------------------------------------------------------
  ! Writes a sample section of the control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine gcmc_sampleCF(unitno)
    Integer, Intent(In) :: unitno
    
    Write(unitno,'(a)') "---- "//Trim(default_gcmc_tag)//" ----"
    Write(unitno,'(a,t30,a)') 'Integer','# Number of iterations per gcmc move'
    Write(unitno,'(a,t30,a)') 'Real',   '# Temperature K '
    Write(unitno,'(a,t30,a)') 'Character','# Eqn of State Tag'
    Write(unitno,'(a,t30,a)') 'Integer','# Number of Simulations for isotherm'
    Write(unitno,'(a,t30,a)') 'Integer','# Statistics block size '
    Write(unitno,'(a,t30,a)') 'Integer','# Number of species '
    Write(unitno,'(t30,a)') '# Blank line (required)'
    Write(unitno,*) "MORE TO FOLLOW HERE"
  End Subroutine gcmc_sampleCF


End Module gcmc







