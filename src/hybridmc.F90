!-----------------------------------------------------------------------------
! This module contains the main data structure and routines for doing 
! Hybrid Monte Carlo (GCMC) simulations.  The routines provided
! here are meant to be used by a driver that loops through the different
! simulation numbers, breaks them into chunks and calls the display routines.
! The simulations can be either GCMC [ constant T,V, P(1..nsims)], or 
! NVT [constant N, V, T(1..nsims)]
! The selection of a species is handled at this level.  All species-specific
! move handling and information storage is handled by "gcmcmoves" which, in
! turn, uses the "moves" module to handle specific moves.
! movelike integration are called thru hybridmoves
!
! Important routines are:
!    hybridmc_init -- initializes the information from a control file
!    hybridmc_dosim -- does a specified number of iterations 
!    hybridmc_displaystats -- displays intermediate statistics for simulation
!
! Needed Improvements:
! 1) remove dangerous hardcoding of 'fast' parameter
!-----------------------------------------------------------------------------
Module hybridmc
  Use auxmoveparams, Only: AuxMoveObjects, auxmoveparams_init
  Use defaults, Only: RDbl, strLen, dashedline,dashedline2, twopi, d_res_file, &
      d_con_file, hplanck, Rgas, Nav, one, zero, dbgflag
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toreal, &
      toupper, allocErrDisplay, cleanstring, real2str, deallocErrDisplay
  Use general, Only: genparams
  Use molecules, Only: molecules_getmass,molecules_getnsorbs
  Use random, Only: rranf
  Use file, Only: file_open, file_getunit, file_gettype
  Use config, Only: AtMolCoords, config_writerestartfile, config_isfixed
  Use datafile, Only: CONFILE
  Use interact, Only: Interaction_Model, &
      interact_checknrgs,interact_initnrgs,interact_chkrecalc
  Use subinteract, Only: Subset_Interactions,subinteract_init, &
      subinteract_chksubint
  Use simcell, Only: SimCell_Params, simcell_getvolume
  Use eos, Only: eos_Models,eos_init,eos_settemperature,eos_setpressure,&
      eos_getfugacity
  Use gcmcmoves, Only: GCMC_Move_Params,gcmcmoves_init,gcmcmoves_move, &
      gcmcmoves_displaystats, gcmcmoves_updatenmolecs, &
      gcmcmoves_display, gcmcmoves_displaysimparams, gcmcmoves_settemplist, &
      gcmcmoves_beginsim
  Use hybridmoves, Only: Hybrid_Move_Params,hybridmoves_init,&
      hybridmoves_integrate, hybridmoves_displaystats
!!$      hybridmoves_displayNrgs,hybridmoves_displayparams
  Use storestats, Only: storestats_displaynrg

 !SDEBUG
  Use cavitylist, Only: cavitylist_displayCubeinfo, cavitylist_displaySorbinfo 
 !SDEBUG

  Implicit None
  Save

  Private
  Public :: hybridmc_init,hybridmc_initdisplay,HYBRIDMC_Params, &
      hybridmc_dosim, hybridmc_chknrgs, &
      hybridmc_getnosims,hybridmc_displaystats,hybridmc_displaysimparams, &
      hybridmc_nooflibsorbs, hybridmc_libUpdtFreq, hybridmc_beginsim, &
      hybridmc_pressure, hybridmc_temperature

  !** The tag marking the beginning of the GCMC section
  Character(len=strLen), Parameter    :: default_hybridmc_tag = &
      "HYBRID GCMC Information"

  Type HYBRIDMC_Params
    !** gcmcsorbs are NOT indexed by sorbate type.  They go from 1-nsorbs
    Type(GCMC_Move_Params), Dimension(:), Pointer   :: gcmcspc 

    Type(Hybrid_Move_Params), Pointer               :: hybrid
    Character(len=strLen)  :: eostag, hybridmc_section_tag, gcmc_section_tag
    Character(len=strLen)  :: lib_update_tag

    Type(EOS_Models)       :: eosparams
    Integer                :: nsorbs
    Integer                :: nsims            ! no. of simulations
    Integer                :: niterations
    Integer                :: blocksize
    Integer                :: gcmc_hybrid_ratio

    !** list of temperatures for the simulations,dimension(nsims)
    !** since temperature is common for all sorbates it is stored here
    !** pressure is different for each components, so that's in gcmcmoves
    Real(kind=Rdbl), Dimension(:), Pointer :: templist

    !** Subset interaction storage, spc-dependent, indexed by species number
    Type(Subset_Interactions), Dimension(:), Pointer :: subints

    Type(AuxMoveObjects),Pointer :: auxmv

  End Type HYBRIDMC_Params

Contains

  !-------------------------------------------------------------------
  ! Initializes the various GCMC parameters from the control file
  ! Note, hybrid-gcmc sims are where librray updates also are done.
  ! Requires:  gcmcparams -- parameters for the GCMC simulation
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            ctrl_filename -- the name of the control file
  !            opt_gcmctag -- optional Tag to look for in control file
  !            opt_update_tag -- optional Tag to be returned for 
  !-------------------------------------------------------------------
  Subroutine hybridmc_init(params, imodel, species, simcell, &
      ctrl_filename, opt_hybridtag, opt_update_tag)
    Type(HYBRIDMC_Params)                          :: params
    Type(Interaction_Model)                        :: imodel
    Type(AtMolCoords), Dimension(:), Intent(inout) :: species
    Type(SimCell_Params), Intent(in)               :: simcell
    Character(*), Intent(in)                       :: ctrl_filename
    Character(*), Optional, Intent(in)             :: opt_hybridtag
    Character(*), Optional, Intent(out)            :: opt_update_tag

    Integer   :: unitno, hybridmclineno, n_gc_sorbs, error, i
    Integer   :: blocksize,nfields
    Integer   :: lineno, ios, spc, nspc
    Logical   :: gcmcmoves_yes, hybridmcmoves_yes, fast, rewind_flag
    Real(kind=RDbl)           :: volume, gcmc_hybrid_ratio
    Character(len=strLen)     :: eostag, tag, line, text, itype
    Character(len=strLen)     :: gcmc_section_tag, hybridmc_section_tag
    Character(len=strLen), Dimension(strLen) :: fields
    Character(len=strLen), Dimension(2) :: two_strings 

    !** Set the tag
    tag = default_hybridmc_tag
    If (Present(opt_hybridtag)) tag = opt_hybridtag

    !** Initialize the auxmv field. This contains other gcmc stuff 
    !** like cavity params, may not be required always. this will change 
    !** your position in ctrlfile
    Call auxmoveparams_init(params%auxmv, species, simcell, ctrl_filename)

    !** Open the ctrl_file if it is not opened
    unitno=file_open(ctrl_filename,110)
    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
        " intializing Hybrid GCMC section"
    !** Find the hybrid-GCMC section
    rewind_flag=.True.
    hybridmclineno = filesrchstr(unitno, tag, line, rewind_flag)
    If (hybridmclineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ",Trim(tag)," in the control file"
      Stop
    End If

    !** Read the hybrid-GCMC section
    Read(unitno, *) params%niterations
    Read(unitno, *) params%nsims
    Read(unitno, '(a)') eostag
    params%eostag = cleanstring(eostag)
    Read(unitno, *) params%blocksize ! used for initializing stats 
    ! blocksize for gcmcmoves 

    !** See whether only one temperature is specified or is there a range
    Read(unitno,'(a)') text
    Call hybridmc_getTempList(params, text)

    !** Read specification for GCMC-type moves
    Read(unitno, '(a)') text
    text=cleanstring(text)
    nfields=split(text, fields, ",")
    gcmcmoves_yes=.True.
    If (cleanstring(fields(1)) == "NO_GCMC_MOVES" )  Then
      gcmcmoves_yes=.False.
      params%nsorbs=0
    Else
      n_gc_sorbs=toint(fields(2),ios)
      If (ios/=0) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      If (n_gc_sorbs==0) Then
        gcmcmoves_yes=.False.
      Else
        gcmcmoves_yes=.True.
        gcmc_section_tag=cleanstring(fields(1))
        params%gcmc_section_tag=gcmc_section_tag
        params%nsorbs=n_gc_sorbs
      Endif
    Endif

    !** Read specification for HMC-type moves
    Read(unitno, '(a)') text
    text=cleanstring(text)
    nfields=split(text, fields, ",")
    hybridmcmoves_yes=.True.
    If (cleanstring(fields(1))=="NO_HYBRIDMC_INTEGRATION")  Then
      hybridmcmoves_yes=.False.
      If (gcmcmoves_yes) Then
        gcmc_hybrid_ratio=1
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End If
    Else
      gcmc_hybrid_ratio=toint(fields(2), ios )
      If ( (gcmc_hybrid_ratio<1).Or.(ios/=0) ) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Unable to parce hybridmc tag section of hybridgcmc secion"
        stop
      Else
        hybridmcmoves_yes=.True.
        hybridmc_section_tag=cleanstring(fields(1))
        params%hybridmc_section_tag=hybridmc_section_tag
      End If
    End If

    params%gcmc_hybrid_ratio=gcmc_hybrid_ratio

    !** Define Tag for library-update-section
    If (Present(opt_update_tag)) Then
      Read(unitno, *) text
      opt_update_tag=cleanstring(text)
      params%lib_update_tag=opt_update_tag
    Else
      params%lib_update_tag="NO_UPDATE"
    End If

    params%nsorbs = n_gc_sorbs

    !** Allocate space for GCMC-type moves if desired
    If (gcmcmoves_yes) Then
      Allocate(params%gcmcspc(n_gc_sorbs), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      !** Set the GCMC move parameter information that should be the same
      !** for all moveable sorbates such as temperature and no. of pts
      Do i=1, n_gc_sorbs
        params%gcmcspc(i)%npts = params%nsims
        params%gcmcspc(i)%blocksize = params%blocksize
      End Do
    Else
      Nullify(params%gcmcspc)
    End If

    !** Allocate space for HMC-type moves if desired
    If (hybridmcmoves_yes) Then
      Allocate(params%hybrid,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Else
      Nullify(params%hybrid)
    End If

    !** Initialize the subset interactions for each moving species
    nspc = molecules_getnsorbs()   !** note changed nspc definition
    Allocate(params%subints(nspc), STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'subints')
    Do spc = 1,nspc
      itype = 'MC'
      If (config_isfixed(species(spc))) itype = 'NULL'
      Call subinteract_init(params%subints(spc),imodel,'Molec_System', &
          itype,(/spc,1,0/),(/0,0,0/))
    End Do

    !** Initialize the GCMC parameters, if desired
    If (Associated(params%gcmcspc)) Then
      !** Find the gcmc section and initialize gcmcmoves 
      rewind_flag=.True.
      two_strings(1)=gcmc_section_tag
      two_strings(2)="SECTION_TAG"
      lineno=filesrchstr(unitno, two_strings , line, rewind_flag)
      If (lineno==0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not find the tags: --TAG :--, --"//&
            Trim(gcmc_section_tag)//"-- in the control file"
        Stop
      Else
        !** allows for different temperatures, not yet functional
        !** need to be fixed, linsert libraries are temperature dependent. 
        Call gcmcmoves_settemplist(params%gcmcspc, params%templist)
        Call gcmcmoves_init(params%gcmcspc, species, simcell, &
            ctrl_filename, params%auxmv)
      End If
    End If

    !** Initialize the HMC parameters, if desired
    If (Associated(params%hybrid)) Then

      !** Find the hybridmc section and initialize
      rewind_flag=.True.
      two_strings(1)=hybridmc_section_tag
      two_strings(2)="SECTION_TAG"

      lineno=filesrchstr(unitno, two_strings, line, rewind_flag)
      If (lineno==0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not find the tags: --TAG :--, --"//&
            Trim(hybridmc_section_tag)//"-- in the control file"
        Stop
      Else
        params%hybrid%npts = params%nsims
        params%hybrid%blocksize = params%blocksize
        Call hybridmoves_init(params%hybrid, imodel, species, &
            simcell, params%templist, ctrl_filename, params%auxmv )
      End If
    End If

    !** Initialize the parameters for the chosen equation of state 
    Rewind(unitno)
    Call eos_init(params%eosparams, ctrl_filename, &
        params%eostag, hybridmclineno)

    !** Generate the fugacitylist 
    volume = simcell_getvolume(simcell)
    Call hybridmc_genfuglist(params, volume)
    Close(unit=unitno)

    !** Update all the stored interactions
    fast = .True.
    Call interact_initnrgs(imodel,fast,species,simcell,.True.)

    Write(*,*) "Energies at the beginning of HMC simulation: "
    Call storestats_displaynrg(imodel%spcstats,3,6)

  End Subroutine hybridmc_init

  !------------------------------------------------
  ! Fill in the fields in templist
  !------------------------------------------------
  Subroutine hybridmc_getTempList(params,text)
    Type(HYBRIDMC_Params), Intent(inout)  :: params
    Character(*),Intent(inout) :: text
    Character(len=strLen),Dimension(strLen) :: fields 
    Integer :: nfields, nTemps, error, i
    Real(kind=RDbl) :: lowtemp, hightemp
    text=trim(stripcmnt(text))
    nfields=split(text, fields, ",")
    nTemps=nfields
    
    Allocate(params%templist(params%nsims), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)


    If (nfields==1) Then
      lowtemp=toreal(fields(1))
      hightemp=lowtemp
      params%templist=lowtemp

    Elseif (nfields==2) Then
      !** check whether both high and low are same
      lowtemp=toreal(fields(1))
      hightemp=toreal(fields(2))
      If (Abs((Abs(hightemp/lowtemp)-one))<1.0e-7) Then
        nfields=1
        nTemps=1
        params%templist=lowtemp
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Need more coding here, to generate temperature list"
      Endif

    Else
      ! we are given a list of temperatures, the number of these should 
      ! correspond to the number of sims
      If (nTemps==params%nsims) Then
        Do i=1,nfields
          params%templist(i)=toreal(fields(i))
        End Do
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) " the number of simulations and number of temperatures ",&
            " given in hybridmc section do not match"
        Stop
      Endif
    Endif
  End Subroutine hybridmc_getTempList

  !------------------------------------------------
  ! Fill in the rest of the fields in the fuglist
  !------------------------------------------------
  Subroutine hybridmc_genfuglist(params, volume)
    Type(HYBRIDMC_Params), Intent(inout)  :: params
    Real(kind=RDbl), Intent(in)       :: volume

    Integer    :: i, j, sorbtype
    Real(kind=RDbl)   :: tk, pp, fug, B, sivolume, mass, murti, Lambda
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__    

    !** Get the fugacity
    Do i=1, params%nsims

      !** Set the temperature and the pressure of the different
      !** species
      tk = params%templist(i)
      Call eos_settemperature(params%eosparams, tk)

      If (associated(params%gcmcspc)) Then

        ! we are doing gcmcspc, so need fugacity
        ! Set the pressure of each component
        Do j=1, params%nsorbs
          sorbtype = params%gcmcspc(j)%spc
          pp = params%gcmcspc(j)%fuglist(i)%pressure
          ! pp is in kPa
          Call eos_setpressure(params%eosparams, sorbtype, pp)
        End Do



        ! Get the fugacity of each component
        Do j=1, params%nsorbs
          sorbtype = params%gcmcspc(j)%spc
          fug = eos_getfugacity(params%eosparams, sorbtype) ! fug [=] kPa
          params%gcmcspc(j)%fuglist(i)%fugacity = fug

          ! Get the excess chemical potential (B in Adams Notation)
          sivolume = volume*1.0e-30*Nav ! convert to m^3/mole
          B = Log(fug*1.0e3*sivolume/(Rgas*tk))
          params%gcmcspc(j)%fuglist(i)%B = B 

          ! Get mu the chemical potential
          ! Calculate the DeBroglie wavelength
          mass = molecules_getmass(sorbtype)
          mass = mass*1.0e-3              ! convert to kg
          Lambda = hplanck/Sqrt(twopi*mass/Nav*Rgas/Nav*tk)
          murti = (B - Log(sivolume/Nav/Lambda**3))
          params%gcmcspc(j)%fuglist(i)%murti = murti
        End Do




      Else

        ! this might be an NVT simulation no ned to call eos

      Endif
    End Do
  End Subroutine hybridmc_genfuglist

  !--------------------------------------------------------------------------
  ! Check the energy calculation and updating system that the HMC uses.  If
  ! unitno is specified, it will produce output.  Otherwise, output will 
  ! only be produced if there is an error.
  ! Requires:  params -- HMC simulation parameters
  !            imodel -- interaction model information
  !            species -- species data structure  
  !            simcell -- the simulation cell information
  !            simno -- simulation number
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional unit number for display
  !--------------------------------------------------------------------------
  Subroutine hybridmc_chknrgs(params,imodel,species,simcell,indent,unitno)
    Type(HYBRIDMC_Params), Intent(InOut)           :: params
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
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) spc
      If (dump) Then
        Call subinteract_chksubint(params%subints(spc),spc,species, &
            simcell,indent,unitno)
      Else
        Call subinteract_chksubint(params%subints(spc),spc,species, &
            simcell,indent)
      End If
    End Do

  End Subroutine hybridmc_chknrgs

  !--------------------------------------------
  ! Gets the total no. of simulations
  !--------------------------------------------
  Integer Function hybridmc_getnosims(params)
    Type(HYBRIDMC_Params), Intent(in)  :: params
    hybridmc_getnosims =  params%nsims
  End Function hybridmc_getnosims

  !--------------------------------------------
  ! Gets the library update frequency
  !--------------------------------------------
  Integer Function hybridmc_libUpdtFreq(params)
    Type(HYBRIDMC_Params), Intent(in)  :: params

    If (Associated(params%gcmcspc)) Then
      hybridmc_libUpdtFreq=0
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

  End Function hybridmc_libUpdtFreq

  !--------------------------------------------
  ! Gets the total no. of simulations
  !--------------------------------------------
  Integer Function hybridmc_nooflibsorbs(params)
    Type(HYBRIDMC_Params), Intent(in)  :: params
    hybridmc_nooflibsorbs = params%nsorbs 
  End Function hybridmc_nooflibsorbs

  !-------------------------------------------
  ! Does one iteration of hybrid-GCMC simulation
  !-------------------------------------------
  Subroutine hybridmc_dosim(params, species, simcell, conf_file,simno)
    Type(HYBRIDMC_Params), Intent(inout)  :: params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: species
    Type(SimCell_Params), Intent(inout)  :: simcell
    Type(CONFILE) ,Intent(in) :: conf_file
    Integer, Intent(in) :: simno
    Integer             :: iter, nmoles, moveno, spc
    Integer             :: i, j, nsorbs,sorbno,iseed,integ_itn
    Real(kind=RDbl)     :: fugacity
    Character(len=strLen)  :: moveType
    Logical :: success, debugflag
    debugflag=.False.

    !** Do the iterations at the set conditions
    Do iter = 1, params%niterations

      !** first do GCMC-type moves 
      If (Associated(params%gcmcspc)) Then
        If (debugflag) Then
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Write(*,*) 'BEGINNING GCMC-type moves'
        Endif
        Do integ_itn=1,params%gcmc_hybrid_ratio
          spc = Int(rranf()*params%nsorbs) + 1

          !** Pick a move type randomly and execute it
          success =  gcmcmoves_move(params%gcmcspc(spc),params%subints,species,&
              simcell,simno)

          !** Do updating of quantities to be averaged
          Call gcmcmoves_updatenmolecs(params%gcmcspc,species)

        End Do
        If (debugflag) Then
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Write(*,*) 'END GCMC-type moves'
        Endif
      End If

      If (Associated(params%hybrid)) Then
        If (debugflag) Then
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Write(*,*) 'BEGINNING Integration-type move'
        Endif

        !** the insert/delet moves were done above
        !** now it is time to do configurational move
        Call hybridmoves_integrate(params%hybrid, species, simcell, simno)

        !** Number averages are stored in gcmcspc. This call will not 
        !** make much difference if integrations are not done very often.
        !** It affects only that stats variable used for display during the run
      If (associated(params%gcmcspc)) &
          Call gcmcmoves_updatenmolecs(params%gcmcspc,species)
 
     
        If (debugflag) Then
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Write(*,*) 'END Integration-type move'
        Endif
      End If

    End Do

  End Subroutine hybridmc_dosim


  !-------------------------------------------------------------------------
  ! Handles the beginning of HGCMC simulations
  ! Requires:  params -- HybridMC parameters
  !           simno -- simulation number
  !            species -- full species data structure
  !            scell -- simulation cell
  !            indent -- number of spaces from left margin
  !            optunit -- optional output unit
  !-------------------------------------------------------------------------
  Subroutine hybridmc_beginsim(params,simno,species,scell,indent,unitno)
    Type(HYBRIDMC_Params), Intent(InOut)     :: params
    Integer, Intent(In)                      :: simno
    Type(AtMolCoords), Dimension(:), Pointer :: species
    Type(SimCell_Params), Intent(In)         :: scell
    Integer, Intent(In)                      :: indent,unitno

    Integer                           :: spc

    If (Associated(params%gcmcspc)) Then
      !** Reset the counters and stats for each species
      Do spc = 1,params%nsorbs
        Call gcmcmoves_beginsim(params%gcmcspc(spc),simno)
      End Do
    Endif
  End Subroutine hybridmc_beginsim


  !--------------------------------------------
  ! Gets the simulation temperature
  !--------------------------------------------
  Real(kind=RDbl) Function hybridmc_temperature(params, simno)
    Type(HYBRIDMC_Params), Intent(In)  :: params
    Integer, Intent(in) :: simno 
    hybridmc_temperature = params%templist(simno)
  End Function hybridmc_temperature

  !--------------------------------------------
  ! Gets the pressure of first gcmc-species
  !--------------------------------------------
  Real(kind=RDbl) Function hybridmc_pressure(params, simno)
    Type(HYBRIDMC_Params), Intent(In)  :: params
    Integer, Intent(in) :: simno 
    hybridmc_pressure=zero ! default for nvt-hmc
    If (associated(params%gcmcspc)) Then
      hybridmc_pressure = params%gcmcspc(1)%fuglist(simno)%pressure
    Endif
  End Function hybridmc_pressure

  
  !-----------------------------------------------------------
  ! Display the HMC statistics for the current simulation
  ! Requires:  params -- HybridMC parameters
  !            species -- full species data structure
  !            imodel -- interaction model
  !            scell -- simulation cell
  !            simno -- simulation number
  !            indent -- number of spaces from left margin
  !            optunit -- optional output unit
  !-----------------------------------------------------------
  Subroutine hybridmc_displaystats(params, species, imodel, scell, &
      simno, indent, optunit)
    Type(HYBRIDMC_Params), Intent(in)              :: params
    Type(AtMolCoords), Dimension(:), Intent(Inout) :: species
    Type(Interaction_Model), Intent(inout)         :: imodel
    Type(SimCell_Params), Intent(in)               :: scell
    Integer, Intent(in)                            :: simno,indent
    Integer, Optional, Intent(in)                  :: optunit

    Integer               :: i, unitno
    Logical               :: fast
    Real(kind=RDbl)       :: nrg_devn
    Character(len=indent) :: blank
    Character(len=strLen) :: restartfile, configfile, string

    blank = Repeat(" ",indent)
    unitno = 6
    If (Present(optunit)) unitno = optunit

    !** Write the simulation no.
    Write(unitno, '(2a)') blank, dashedline
    Write(unitno, '(2a)') blank, "The HYBRID GCMC Stats:"
    If (Associated(params%gcmcspc)) Then
      Call gcmcmoves_displaystats(params%gcmcspc, simno, indent+2, unitno)
    Else 
      Write(unitno,'(2x,2a)') blank, "--- No GCMC Types of moves ---"
    Endif

    If (Associated(params%hybrid)) Then
      Call hybridmoves_displaystats(params%hybrid,species,simno,indent+2,unitno)
    Else 
      Write(unitno,'(2x,2a)') blank, "--- No HMC Types of moves ---"
    End If
    
    !** Check storage versus fresh, full-system evaluation of interactions
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      fast = .True.
      nrg_devn = interact_checknrgs(imodel,fast,species,scell,indent+2,unitno)
      string = real2str(nrg_devn,6)
      Write(unitno,'(4x,4a)') blank,"Deviation between stored and newly ", &
          "calculated energies: ", Trim(string)
    End If

  End Subroutine hybridmc_displaystats

  !------------------------------------------------------------------
  ! Display the hybrid-gcmc simulation parameters for the given 
  ! simulation number
  ! Requires:  params -- HybridMC parameters
  !            simno -- simulation number
  !            indent -- number of spaces from left margin
  !            optunit -- optional output unit
  !-----------------------------------------------------------------
  Subroutine hybridmc_displaysimparams(params, simno, indent, optunit)
    Type(HYBRIDMC_Params), Intent(In) :: params
    Integer, Intent(In)               :: simno,indent
    Integer, Optional, Intent(In)     :: optunit

    Integer               :: i, unitno, funit
    Character(len=indent) :: blank
    Character(len=strLen) :: restartfile, configfile

    blank = Repeat(' ',indent)
    unitno = 6
    If (Present(optunit)) unitno = optunit

    !** Get the configuration and restart file names
    Call file_gettype(d_con_file, configfile, funit)
    Call file_gettype(d_res_file, restartfile,funit)

    !** Write the simulation no.
    Write(unitno, '(2a)') blank, dashedline
    Write(unitno, '(2a)') blank, "The HYBRID GCMC Simulation Parameters:"
    Write(unitno, '(a,2x,a,i4)') blank, "Simulation Number   : ", simno
    Write(unitno, '(a,2x,2a)') blank, "Configuration file  : ", &
        Trim(configfile)
    Write(unitno, '(a,2x,2a)') blank, "Restart file        : ", &
        Trim(restartfile)
    Write(*,*) blank,"ADD MORE CODE HERE, for what?"

    If (Associated(params%gcmcspc)) Then
      Call gcmcmoves_displaysimparams(params%gcmcspc,simno,indent+2,unitno)
    End If

!!$    Call hybridmoves_displaysimparams(params%hybrid,simno,indent+2,unitno)

  End Subroutine hybridmc_displaysimparams

  !-----------------------------------------------------------------
  ! Display the gcmc initialization parameters
  ! Requires:  params -- HybridMC parameters
  !            indent -- number of spaces from left margin
  !            optunit -- optional output unit
  !-----------------------------------------------------------------
  Subroutine hybridmc_initdisplay(params, indent, optunit)
    Type(HYBRIDMC_Params), Intent(in) :: params
    Integer, Intent(in)               :: indent
    Integer, Optional, Intent(in)     :: optunit

    Integer               :: i, unitno
    Character(len=indent) :: blank
    Character(len=strLen) :: restartfile, configfile

    blank = Repeat(' ',indent)
    unitno = 6
    If (Present(optunit)) unitno = optunit
    
    Write(unitno, '(2a)') blank, dashedline
    Write(unitno, '(2a)') blank, "The hybrid-GCMC Parameters Section:"
    Write(unitno, '(a,2x,a,i6)')    blank,"No. of iterations   : ", &
        params%niterations
    Write(unitno, '(a,2x,a,i6)')    blank, "No. of simulations : ", &
        params%nsims
    Write(unitno, '(a,2x,2a)')      blank, "EOS Tag            : ", &
        Trim(params%eostag)
    Write(unitno, '(a,2x,a,i6)')    blank, "stats blocksize    : ", &
        params%blocksize
    Write(unitno, '(a,2x,a,2f12.3)')blank, "Temp.  Range (K)   : ", & 
        params%templist(1), params%templist(params%nsims)

    !** Display GCMC moves
    If (Associated(params%gcmcspc)) Then
      Write(unitno, '(a,2x,a,i6)')  blank, "No. of gcmc-sorbs  : ", &
          params%nsorbs
      Call gcmcmoves_display(params%gcmcspc, indent+2, unitno)
    Else
      Write(unitno, '(a,4x,a )') blank, "-- No gcmc moves will be used --"
    End If

    !** Display Hybrid moves
    If (Associated(params%hybrid)) Then
      Write(unitno, '(a,2x,2a)')    blank, "hybrid section tag : ", &
        Trim(params%hybridmc_section_tag)      
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "need hybridmoves display code here"
!!$    Call hybridmoves_displayparams(params%hybrid, indent+4, unitno)
    Else
      Write(unitno,'(a,4x,a)') blank, "No hybridmc moves will be used"
    End If
    
    Write(unitno,'(a,4x,2a)') blank, "lib update      tag : ", &
        Trim(params%lib_update_tag)      

 End Subroutine hybridmc_initdisplay

End Module hybridmc







