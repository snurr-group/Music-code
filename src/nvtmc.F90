!-----------------------------------------------------------------------------
! This module contains the main data structure and routines for doing 
! Canonical Monte Carlo (NVTMC) simulations.  The routines provided
! here are meant to be used by a driver that loops through the different
! simulation numbers, breaks them into chunks and calls the display routines.
!
! The selection of a species is handled at this level.  All species-specific
! move handling and information storage is handled by "nvtmcmoves" which, in
! turn, uses the "moves" module to handle specific moves.
!
! Important routines are:
!    nvtmc_init -- initilizes the information from a control file
!    nvtmc_dosim -- does a specified number of interations 
!    nvtmc_displaystats -- displays intermediate statistics for simulation
!
! Needed Improvements:
! 1) remove dangerous hardcoding of 'fast' parameter
! 2) testing for more than one species
!-----------------------------------------------------------------------------

Module nvtmc
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl, strLen, dashedline,dashedline2, twopi, d_res_file, &
      d_con_file, hplanck, MAX_SORBS, Rgas, Nav, scalepe, kcalmole_kb, lstrLen, &
      dbgflag
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, readblank, int2str, toreal, deallocerrdisplay, isdigit, &
      real2str, checkandstop 
  Use file, Only: file_getunit,file_gettype,file_open
  Use general, Only: genparams
  Use random, Only: rranf
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use eos, Only: eos_Models,eos_init,eos_settemperature,eos_setpressure,&
      eos_getfugacity
  Use config, Only: AtMolCoords, config_writerestartfile, config_getnmoles, &
      config_display, config_isfixed, config_getnmoleslist
  Use simcell, Only: SimCell_Params, simcell_getvolume
  Use nvtmcmoves, Only: NVTMC_Move_Params,nvtmcmoves_init,nvtmcmoves_move, &
      nvtmcmoves_displaystats, nvtmcmoves_display, nvtmcmoves_displaysimparams, &
      nvtmcmoves_dispspcstats,nvtmcmoves_clean
  Use molecule, Only: MolecularParams
  Use molecules, Only: molecules_getmass, molecules_getnsorbs
  Use storestats, Only: Species_Stats, storestats_getcoul, &
      storestats_getnoncoul, storestats_displaynrg
  Use stats, Only: stats_update,stats_getvalue,stats_reset
  Use interact, Only: Interaction_Model,interact_checknrgs,interact_initnrgs
  Use subinteract, Only: Subset_Interactions, subinteract_init, &
      subinteract_int, subinteract_oldnrg, subinteract_newnrg
  Use storetop, Only: storetop_display
  Use storesym, Only: storesym_extract,storesym_sum
  Use storebase, Only: EnergyPlus, storebase_disp

  Implicit None
  Save

  Private
  Public :: nvtmc_init,nvtmc_initdisplay,NVTMC_Params,nvtmc_dosim,&
      nvtmc_getnosims,nvtmc_displaystats,nvtmc_displaysimparams, &
      nvtmc_beginsim,nvtmc_endsim,nvtmc_sampleCF,nvtmc_clean

  !** The tag marking the beginning of the NVTMC section
  Character(len=strLen), Parameter    :: default_nvtmc_tag = &
      "NVTMC Information"

  !** Contains all the NVTMC information for the whole simulation
  Type NVTMC_Params
    Integer                :: nspc
    Integer                :: nsims        !* no. of simulations
    Integer                :: niterations

    Character(len=strLen)  :: eostag,tempfile
    Type(EOS_Models)       :: eosparams

    !** lists of temperature (in Kelvin) and rti = 1.0/(Rgas*T_Kelvin)
    Real(kind=RDbl), Dimension(:), Pointer :: templist,rtilist

    !** nvtmcspc are NOT indexed by species type.  They go from 1->nspc
    Integer, Dimension(:), Pointer                   :: spclist
    Type(NVTMC_Move_Params), Dimension(:), Pointer   :: nvtmcspc 

    !** Subset interaction storage, species-dependent, indexed by species number
    Type(Subset_Interactions), Dimension(:), Pointer :: subints

    !** common objects like cavity lists can be stored here
    Type(AuxMoveObjects),Pointer :: auxmv

  End Type NVTMC_Params

Contains
  !-------------------------------------------------------------------
  ! Initializes the various NVTMC parameters from the control file
  ! Requires: nvtmcparams -- parameters for the NVTMC simulation
  !           imodel -- interaction model information
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           ctrl_filename -- the name of the control file
  !           opt_nvtmctag -- optional Tag to look for in control file
  !-------------------------------------------------------------------
  Subroutine nvtmc_init(nvtmcparams, imodel, species, simcell, ctrl_filename, &
      opt_nvtmctag)
    Type(NVTMC_Params), Intent(InOut)              :: nvtmcparams
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Character(*), Intent(In)                       :: ctrl_filename
    Character(*), Optional, Intent(In)             :: opt_nvtmctag

    Integer                :: error,i,j,spc
    Integer                :: unitno,nvtmclineno,nspc,blocksize
    Logical                :: fast
    Real(kind=RDbl)        :: volume
    Character(len=strLen)  :: eostag, tag, line, temperatureline, itype
    Character(len=lstrLen) :: errormsg

    !** Pad the 'tag' to distinguish it from other incidents of the name
    If (Present(opt_nvtmctag)) Then
      Write(tag,'(2a)') '- ',Trim(opt_nvtmctag)
    Else
      Write(tag,'(2a)') '- ',Trim(default_nvtmc_tag)
    End If
    
    !** Open the ctrl_file if it is not opened
    unitno = file_open(ctrl_filename)
    Rewind(unitno)
    
    !** Find the NVTMC section
    nvtmclineno = filesrchstr(unitno, tag, line)
    If (nvtmclineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag '",Trim(tag),"' in the control file"
      Stop
    End If
    
    Read(unitno,*) line
    line = stripcmnt(line)
    nvtmcparams%niterations = toint(line,'number of NVTMC steps')
    Read(unitno,*) nvtmcparams%nsims
    Read(unitno,*) temperatureline
    Call nvtmc_gettemperatures(nvtmcparams,temperatureline)

    Read(unitno,'(a)') eostag
    nvtmcparams%eostag = stripcmnt(eostag)

    Read(unitno,*) blocksize
    Read(unitno,*) nspc    
    nvtmcparams%nspc = nspc

    !** Read the blankline at the end
    errormsg = 'must be a blank line before species parameters in NVTMC section'
    Call readblank(unitno,__FILE__,__LINE__,errormsg)
    
    !** Allocate the nvtmcparams%nvtmcspc 
    Allocate(nvtmcparams%nvtmcspc(nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
    !** Set the NVTMC move parameter information that should be the same
    !** for all moveable species such as temperature and no. of pts
    Do i = 1,nspc
      nvtmcparams%nvtmcspc(i)%blocksize = blocksize
    End Do

    !** Initialize the subset interactions for each moving species
    nspc = molecules_getnsorbs()   !** note changed nspc definition
    Allocate(nvtmcparams%subints(nspc), STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'subints')
    Do spc = 1,nspc
      itype = 'MC'
      If (config_isfixed(species(spc))) itype = 'NULL'
      Call subinteract_init(nvtmcparams%subints(spc),imodel,'Molec_System', &
          itype,(/spc,1,0/),(/0,0,0/))
    End Do
    
    !** The rest of the stuff will be read and initialized by the 
    !** nvtmcmove init
    Call nvtmcmoves_init(nvtmcparams%nvtmcspc,species,simcell,ctrl_filename, &
        nvtmcparams%auxmv )

    !** Allocate species number list
    Allocate(nvtmcparams%spclist(nvtmcparams%nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    
    Do spc = 1,nvtmcparams%nspc
      nvtmcparams%spclist = nvtmcparams%nvtmcspc%spc  !HaCK
    End Do
    
    !** Initialize the parameters for the chosen equation of state 
    Rewind(unitno)
    Call eos_init(nvtmcparams%eosparams, ctrl_filename, nvtmcparams%eostag, &
        nvtmclineno)
        
    Close(unit=unitno)
    
    !** Update all the energies stored in molecules
    fast = .True.
    Call interact_initnrgs(imodel,fast,species,simcell,.True.)

#ifdef DEBUG
    !** Just some feedback
    Do i = 1,nspc
      Do j = i,nspc
        Write(*,'(a,2I4,2f18.6)') "Species:",i,j, &
            molecules_getcoul(i,j,'inst'),molecules_getnoncoul(i,j,'inst')
      End Do
    End Do
#endif
    
  End Subroutine nvtmc_init
 
  !---------------------------------------------------------------------------
  ! Performs "ninterations" of the specified NVTMC simulation.
  ! We need to ensure that insertions and deletions are done with equal 
  ! probability for each species type. The algorithm used here is 
  ! similar to the one used by Karavias and Myers (Mol. Sim. (1991) v8, p.51)
  ! 1) Pick a species with equal probability.
  ! 2) Decide whether we want to attempt an
  !    insertion, deletion, translation or rotation.
  ! Enforcement of these constraints are handled in nvtmcmoves (init)
  ! Requires: nvtmcparams -- parameters for the NVTMC simulation
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           simno -- the simulation number
  !---------------------------------------------------------------------------
  Real(kind=RDbl) Function nvtmc_dosim(nvtmcparams,species,simcell,simno)
    Type(NVTMC_Params), Intent(InOut)                :: nvtmcparams
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Type(SimCell_Params), Intent(In)                 :: simcell
    Integer, Intent(In)                              :: simno

    Integer               :: iter,nmoles,spc,indx
    Logical               :: success

    Integer, Dimension(20)   :: molcount

    !** Zero the acceptance counter
    nvtmc_dosim = 0.0_RDbl

    !** Make sure that the structure has been initialized
    If (.Not. Associated(nvtmcparams%nvtmcspc)) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " NVTMC Species must be initialized before starting simulation"
      Stop
    End If

    Do iter = 1,nvtmcparams%niterations
      !** Pick a species randomly, making sure it has molecules
      Do 
        indx = Int(rranf()*nvtmcparams%nspc) + 1
        spc = nvtmcparams%spclist(indx)
        nmoles = config_getnmoles(species,spc)
        If (nmoles /= 0) Exit
      End Do

      !** Pick a move type randomly and execute it
      success = nvtmcmoves_move(nvtmcparams%nvtmcspc(indx), &
          nvtmcparams%rtilist(simno),nvtmcparams%subints,species,simcell,simno)

      If (success) Then
        nvtmc_dosim = nvtmc_dosim + 1.0_RDbl
      End If
    End Do

    !** Calculate the acceptance ratio and return it
    nvtmc_dosim = nvtmc_dosim/nvtmcparams%niterations

  End Function nvtmc_dosim

  !----------------------------------------------------------------
  ! Parses the temperature line to set the temperatures of the
  ! sorbate type.  It can either do the temperatures at equal
  ! intervals or read them from a file.
  !---------------------------------------------------------------
  Subroutine nvtmc_gettemperatures(nvtmcparams, temperatureline)
    Type(NVTMC_Params), Intent(inout)  :: nvtmcparams
    Character(*), Intent(in)           :: temperatureline

    Integer                :: nfields, npts, i, error
    Real(kind=RDbl)        :: startT, endT, tincr
    Character(len=strLen)  :: stripline
    Character(len=strLen), Dimension(10)   :: strfields

    !** allocate memory for the temperature lists
    Allocate(nvtmcparams%templist(nvtmcparams%nsims), STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'templist')

    stripline = stripcmnt(temperatureline)

    !** Check the no. of fields.  If the no. of fields is more than one then
    !** we need to either generate the temperatues at equal intervals, read
    !** them from a file or just read a single number.  
    nfields = split(stripline, strfields, ",")
    If (nfields == 1) Then
      If (isdigit(strfields(1)(1:1))) Then  !** just use this number
        startT = toreal(strfields(1))
        Call nvtmc_gentemperatures(nvtmcparams,startT,startT)
      Else        !** Read from a file
        nvtmcparams%tempfile = strfields(1)
        Call nvtmc_readtemperatures(nvtmcparams, nvtmcparams%tempfile)
      End If

    Else   !** Generate the temperatures
      startT = toreal(strfields(1))
      endT   = toreal(strfields(2))
      Call nvtmc_gentemperatures(nvtmcparams, startT, endT)
    End If

    Call nvtmc_getinvtemps(nvtmcparams)
    
  End Subroutine nvtmc_gettemperatures

  !----------------------------------------------------------------------
  ! Generate the rti = 1.0/(Rgas*T_Kelvin) values from the temperatures
  !----------------------------------------------------------------------
  Subroutine nvtmc_getinvtemps(nvtmcparams)
    Type(NVTMC_Params), Intent(InOut)  :: nvtmcparams

    Integer     :: i, error

    If (.Not. Associated(nvtmcparams%templist)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          "  Temperature list must first be associated "
      Stop
    End If

    Allocate(nvtmcparams%rtilist(nvtmcparams%nsims), STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'rtilist')

    Do i = 1,nvtmcparams%nsims
      nvtmcparams%rtilist(i) = 1.0_RDbl/(kcalmole_kb*nvtmcparams%templist(i))
    End Do

  End Subroutine nvtmc_getinvtemps

  !-----------------------------------------------------
  ! This routine reads in the temperature values from the
  ! file "filename"
  !-----------------------------------------------------
  Subroutine nvtmc_readtemperatures(nvtmcparams, filename)
    Type(NVTMC_Params), Intent(inout)  :: nvtmcparams    
    Character(*), Intent(in)           :: filename

    Integer      :: unitno, error, npts, i
    
    !** Open the file
    unitno = isfileopen(filename)
    If (unitno < 0) Then
      unitno = file_getunit(filename)
      Open(unit = unitno, file=filename, status='old', IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            "  Could not open file ", trim(filename)
        Stop
      End If
    End If

    !** Read the no. of points
    Read(unitno, *) npts

    !** Check that the no. of points is the same as specified in the
    !** control file
    If (npts /= nvtmcparams%nsims) Then
      Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
          " The no. of points in the file ", Trim(filename), &
          " does not match that in the control file"
      Stop
    End If
    Read(unitno, *) (nvtmcparams%templist(i), i=1, npts)
    Close(unitno)
    Return
  End Subroutine nvtmc_readtemperatures

  !----------------------------------------------------
  ! This routine chooses between a linear or a log
  ! scale to increment/decrement the temperature based
  ! on the temperature endpoint values
  !----------------------------------------------------
  Subroutine nvtmc_gentemperatures(nvtmcparams, startT, endT)
    Type(NVTMC_Params), Intent(inout)  :: nvtmcparams
    Real(kind=RDbl), Intent(in)        :: startT, endT

    Integer          :: npts, i
    Real(kind=RDbl)  :: tincr

    npts   = nvtmcparams%nsims
    If (Abs(Log10(endT/startT)) > Log10(20.0)) Then
      ! use a log scale
      If (npts /= 1) Then
        tincr  = Exp(Log(endT/startT)/(npts - 1.0_RDbl))
      End If
      Do i=1, npts
        nvtmcparams%templist(i) = startT*(tincr**(i-1))
      End Do
    Else
      If (npts /= 1) Then
        tincr  = (endT - startT)/(npts - 1.0_RDbl)
      End If
      Do i=1, npts
        nvtmcparams%templist(i) = startT + (i - 1.0_RDbl)*tincr
      End Do
    End If
    Return
  End Subroutine nvtmc_gentemperatures

  !--------------------------------------------
  ! Gets the total no. of simulations
  !--------------------------------------------
  Integer Function nvtmc_getnosims(nvtmcparams)
    Type(NVTMC_Params), Intent(In)  :: nvtmcparams

    nvtmc_getnosims = nvtmcparams%nsims

  End Function nvtmc_getnosims
  
  !---------------------------------------------------------------------------
  ! Display the NVTMC statistics for the current simulation
  ! Requires:  nvtmcparams -- NVTMC simulation parameters
  !            imodel -- interaction model information
  !            species -- species data structure  
  !            simcell -- the simulation cell information
  !            simno -- simulation number
  !            indent -- no. of spaces from the left margin
  !            optunit -- optional unit number for display
  !---------------------------------------------------------------------------
  Subroutine nvtmc_displaystats(nvtmcparams,imodel,species,simcell,simno, &
      indent,optunit)
    Type(NVTMC_Params), Intent(InOut)              :: nvtmcparams
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Integer, Intent(In)                            :: simno,indent
    Integer, Intent(In), Optional                  :: optunit

    Integer               :: i, unitno
    Integer               :: spc,m
    Logical               :: fast,skip_intra_flag,success
    Real(kind=RDbl)       :: nrg_devn,nrg_stored,nrg_recalc,max_devn
    Character(len=indent) :: blank
    Character(len=strLen) :: restartfile,configfile,string
    Type(EnergyPlus)      :: nrgs

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    !** Write the simulation no.
    Write(unitno, '(2a)') blank, dashedline
    Write(unitno, '(2a)') blank, "The NVTMC Stats:"
    Call nvtmcmoves_displaystats(nvtmcparams%nvtmcspc,simno,indent+2,unitno)

    If (Trim(genparams%displaymode) == "VERBOSE") Then
      fast = .True.
      nrg_devn = interact_checknrgs(imodel,fast,species,simcell,indent+4,unitno)
      string = real2str(nrg_devn,6)
      Write(unitno,'(4x,3a)') blank,"Deviation between stored and newly &
          &calculated energies: ", Trim(string)
    End If

#ifdef DEBUG
    !** Check interact's consistency in molec-molec interactions
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(unitno,'(2a)') blank,'Checking storage vs recalc reproducibility'
      Write(unitno,'(a,2a5,3a12)') blank,'spc','mol','stored','recalc','diff'
      fast = .True.
      skip_intra_flag = .False.

      Do spc = 1,molecules_getnsorbs()
        Do m = 1,config_getnmoles(species,spc)
          If (config_isfixed(species(spc))) Cycle

          !** Get the stored value
          success = interact_int(nvtmcparams%nvtmcspc(spc)%subint,species, &
              simcell,fast,.False.,skip_intra_flag,1000.0_RDbl,(/spc,m/))
          Call checkandstop(success,__FILE__,__LINE__,' nrg calc')
          nrg_stored = interact_oldnrg(nvtmcparams%nvtmcspc(spc)%subint)

          !** Get the recalculated value
          success = interact_int(nvtmcparams%nvtmcspc(spc)%subint,species, &
              simcell,fast,.True.,skip_intra_flag,1000.0_RDbl,(/spc,m/))
          Call checkandstop(success,__FILE__,__LINE__,' nrg calc')
          nrg_recalc = interact_newnrg(nvtmcparams%nvtmcspc(spc)%subint)

          nrg_devn = Abs(nrg_stored - nrg_recalc)
          If (nrg_devn > max_devn) max_devn = nrg_devn
          If (nrg_devn > 1.0e-8_RDbl) Then
            Write(unitno,'(a,2i5,3f12.5)') blank,spc,m,nrg_stored, &
                nrg_recalc,nrg_devn
          End If
        End Do
      End Do

      If (max_devn > 1.0e-7_RDbl) Then
        Write(0,'(1x,2a,i4,a)') __FILE__,": ",__LINE__, &
            ' ERROR: storage versus recalc shows unacceptable reproducibility'
        Stop

        !** debugging stuff
        Call storetop_display(imodel%ff%results,.True.,2,6)
        Call storesym_sum(imodel%ff%results%ncoul)
        success = storesym_extract(imodel%ff%results%ncoul,nrgs, &
            (/1,7,0/),(/0,0,0/))
        Write(*,*) success,Trim(storebase_disp(nrgs))
        Stop
      End If
    End If
#endif

  End Subroutine nvtmc_displaystats

  !------------------------------------------------------------------
  ! Display the nvtmc simulation parameters for a specific simulation.
  ! Useful for summarizing parameters at the beginning of a new sim
  ! Requires: nvtmcparams -- NVTMC simulation parameters
  !           simno -- the simulation number
  !           indent -- no. of spaces from the left margin
  !           optunit -- optional unit number for display
  !-----------------------------------------------------------------
  Subroutine nvtmc_displaysimparams(nvtmcparams,simno,indent,optunit)
    Type(NVTMC_Params), Intent(In) :: nvtmcparams
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
    Write(unitno, '(2a)') blank,"The NVTMC Simulation Parameters:"
    Write(unitno, '(a,2x,a,i4)') blank,"Simulation Number : ", simno
    Write(unitno, '(a,2x,2a)') blank,"Configuration file: ", &
        Trim(configfile)
    Write(unitno, '(a,2x,2a)') blank, "Restart file      : ", &
        Trim(restartfile)
    Call nvtmcmoves_displaysimparams(nvtmcparams%nvtmcspc,simno,indent+2,unitno)

  End Subroutine nvtmc_displaysimparams

  !-------------------------------------------------------------------------
  ! Handles the beginning of NVTMC simulations
  ! 1) Write feedback to screen
  ! Requires: nvtmcparams -- NVTMC simulation parameters
  !           simno -- simulation number
  !           species -- species data structure  
  !           simcell -- the simulation cell information
  !           indent -- no. of spaces from the left margin
  !           unit -- display unit number
  !-------------------------------------------------------------------------
  Subroutine nvtmc_beginsim(nvtmcparams,simno,species,simcell,indent,unit)
    Type(NVTMC_Params), Intent(InOut)            :: nvtmcparams
    Integer, Intent(In)                          :: simno
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell
    Integer, Intent(In)                          :: indent,unit

    Integer                           :: spc
    Character(len=indent)             :: blank

    blank = Repeat(' ',indent)    

    !** Reset the counters for each species
    Do spc = 1,nvtmcparams%nspc
!LC      Call stats_reset(nvtmcparams%nvtmcspc(spc)%nmoles)  !shouldn't be here
    End Do

  End Subroutine nvtmc_beginsim

  !-------------------------------------------------------------------------
  ! Handles the ending of NVTMC simulations
  ! 1) Check the accumulated and fresh species-species and intra energies
  ! 2) Write feedback to screen
  ! 3) Optionally write to a restartfile (if stopfile /= '')
  ! 4) halt the program only if desired
  ! Step (4) can be useful for getting a starting file for subsequent MD runs.
  ! Requires: nvtmcparams -- NVTMC simulation parameters
  !           simno -- simulation number
  !           species -- species data structure  
  !           simcell -- the simulation cell information
  !           stopfile -- name of restart file
  !           indent -- no. of spaces from the left margin
  !           unit -- display unit number
  !           stopflag -- optional flag to signal program halt
  !-------------------------------------------------------------------------
  Subroutine nvtmc_endsim(nvtmcparams,imodel,simno,species,simcell,stopfile, &
      indent,unit,stopflag)
    Type(NVTMC_Params), Intent(In)               :: nvtmcparams
    Integer, Intent(In)                          :: simno
    Type(Interaction_Model), Intent(InOut)       :: imodel
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell
    Character(*), Intent(In)                     :: stopfile
    Integer, Intent(In)                          :: indent,unit
    Logical, Intent(In), Optional                :: stopflag 

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
    Else If (rmsdev > 1.0E-4_RDbl) Then
      Write(0,'(1x,2a,i4,a,e14.4)') __FILE__,": ",__LINE__, &
          ' WARNING: rms deviation between stored and'
      Write(0,'(5x,a,e14.4)') 'fresh energies is significant: ',rmsdev
    End If

    !** Feedback 
    Write(unit,'(2a)') blank,dashedline
    string = int2str(simno)
    Write(unit,'(3a,4x,a,e9.3)') blank,'Ending NVTMC simulation number ', &
        Trim(string),'Nrg dev check: ',rmsdev

    !** Commented out because there's nothing to report
!    Do spc = 1,nvtmcparams%nspc
!      Call nvtmcmoves_dispspcstats(nvtmcparams%nvtmcspc(spc), &
!          simno,indent+2,unit)
!    End Do
    
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
    
  End Subroutine nvtmc_endsim

  !-----------------------------------------------------------------------
  ! Display the nvtmc initialization parameters, ie full information 
  ! Requires: nvtmcparams -- NVTMC simulation parameters
  !           indent -- no. of spaces from the left margin
  !           optunit -- optional unit number for display
  !-----------------------------------------------------------------------
  Subroutine nvtmc_initdisplay(nvtmcparams,indent,optunit)
    Type(NVTMC_Params), Intent(In)        :: nvtmcparams
    Integer, Intent(In)                   :: indent
    Integer, Optional, Intent(In)         :: optunit

    Integer               :: i, unitno
    Character(len=indent) :: blank
    Character(len=strLen) :: string

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    Write(unitno,'(2a)') blank, dashedline
    Write(unitno,'(2a)') blank, "The NVTMC Parameters Section:"
    string = int2str(nvtmcparams%niterations)
    Write(unitno,'(a,2x,2a)') blank,"Number of iterations: ",Trim(string)
    string = int2str(nvtmcparams%nsims)
    Write(unitno,'(a,2x,2a)') blank, "Number of simulations: ",Trim(string)

    Write(unitno,'(a,2x,a)',ADVANCE='NO') blank, "Temperature(s) (K): "
    Do i = 1,nvtmcparams%nsims
      string = real2str(nvtmcparams%templist(i),4)
      Write(unitno,'(a,2x)',ADVANCE='NO') Trim(string)
    End Do
    Write(unitno,*)

    string = int2str(nvtmcparams%nspc)
    Write(unitno,'(a,2x,2a)') blank, "Number of species: ", Trim(string)
        
    Call nvtmcmoves_display(nvtmcparams%nvtmcspc,indent+2,unitno)

  End Subroutine nvtmc_initdisplay

  !----------------------------------------------------------------------------
  ! Writes a sample section of the control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine nvtmc_sampleCF(unitno)
    Integer, Intent(In) :: unitno
    
    Write(unitno,'(a)') "---- "//Trim(default_nvtmc_tag)//" ----"
    Write(unitno,'(a,t30,a)') 'Integer','# Number of iterations per NVTMC move'
    Write(unitno,'(a,t30,a)') 'Integer','# Number of simulation points'
    Write(unitno,'(a,t30,a)') 'Real',&
        '# Temperature(s) (K) (can specific range here: lo, hi)'
    Write(unitno,'(a,t30,a)') 'Character', &
        '# Tag for the equation of state (NULL = Ideal Gas)'
    Write(unitno,'(a,t30,a)') 'Integer','# Statistics block size '
    Write(unitno,'(a,t30,a)') 'Integer','# Number of species '
    Write(unitno,'(t30,a)') '# Blank line (required)'
    Write(unitno,*) "MORE TO FOLLOW HERE"
  End Subroutine nvtmc_sampleCF

  !----------------------------------------------------------------------
  ! Cleans up the various pointers once we are done with this structure.
  !----------------------------------------------------------------------
  Subroutine nvtmc_clean(nvtmcparams)
    Type(NVTMC_Params), Intent(InOut)        :: nvtmcparams

    Integer    :: i,nspc,error

    Deallocate(nvtmcparams%templist, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'templist')

    nspc = Size(nvtmcparams%nvtmcspc, 1)
    Do i = 1,nspc
      Call nvtmcmoves_clean(nvtmcparams%nvtmcspc(i))
    End Do

    Deallocate(nvtmcparams%nvtmcspc, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'nvtmcspc')

  End Subroutine nvtmc_clean

End Module nvtmc
 
