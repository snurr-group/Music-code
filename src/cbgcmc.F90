!-----------------------------------------------------------------------------
! This module contains the main data structure and routines for doing 
! Configurational biased Grand Canonical Monte Carlo (GCMC) simulations.  
! The routines provided here are meant to be used by a driver that loops 
! through the different  simulation numbers, breaks them into chunks and 
! calls the display routines. This is very similar to gcmc.F90 in its structure
! This can work only with the branchedcoords and brmoves modules, 
! other gen coords can not be used
! The selection of a species is handled at this level.  All species-specific
! move handling and information storage is handled by "cbgcmoves" which, in
! turn, uses the "moves" module to handle specific moves.
!
! Important routines are:
!    cbgcmc_init -- initilizes the information from a control file
!    cbgcmc_dosim -- does a specified number of interations 
!    cbgcmc_displaystats -- displays intermediate statistics for simulation
!
! Needed Improvements:
! 1) remove dangerous hardcoding of 'fast' parameter
!-----------------------------------------------------------------------------


Module cbgcmc
  Use brmoves, Only: brmoves_initTrialArrays
  Use cbgcmoves, Only : CBGCMC_Move_Params,cbgcmoves_init,cbgcmoves_move, &
      cbgcmoves_beginsim, cbgcmoves_nvtmove, &
      cbgcmoves_displaystats, cbgcmoves_display, &
      cbgcmoves_displaysimparams, cbgcmoves_updatenmolecs, &
      cbgcmoves_dispspcstats
  Use cbgcaccel, Only : CBGC_Accelerator, cbgcaccel_init, cbgcaccel_dosim 
  Use config, Only: AtMolCoords, config_writerestartfile, config_getnmoles,&
      config_isfixed
  Use defaults, Only: RDbl, strLen, dashedline,dashedline2, twopi, d_res_file &
      ,d_con_file, hplanck, MAX_SORBS, NO_OF_INTRA_POTS, STRETCH_INDEX,&
      BENDING_INDEX,TORSION_INDEX,CONSTRAINT_INDEX,TOTAL_INDEX,&
      INTRAPAIR_INDEX,Rgas, Nav, scalepe, kcalmole_kb, one,zero
  Use eos, Only : eos_Models,eos_init,eos_settemperature,eos_setpressure,&
      eos_getfugacity,eos_getConfInteg
  Use file, Only: file_getunit,file_gettype,file_open
  Use forcefield, Only : forcefield_getssint, forcefield_checknrgs, forcefield_initnrgs
  Use frame, Only : Frame_Library_Params
  Use general, Only: genparams
  Use intramolecular, Only : intramolecular_hasint
  Use molecule, Only : MolecularParams
  Use molecules, Only: molecules_getmass, molecules_getcoul &
      ,molecules_getnoncoul, molecules_getintranrg, molecules_getnsorbs, &
      molecules_updateEnergySS, molecules_getpointer, molecules_updateenergy, &
      molecules_name, molecules_updateIntraEnergy, molecules_displaynrg
  Use random,Only :rranf
  Use smap,Only:Smap_Params,smap_init
  Use stats, Only : Statistics, stats_init,stats_update,stats_getvalue, &
      stats_display 
  Use simcell, Only: SimCell_Params, simcell_getvolume
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, int2str, readblank
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag

  Implicit None
  Save

  Private
  Public :: cbgcmc_init,cbgcmc_initdisplay,CBGCMC_Params,cbgcmc_dosim,&
      cbgcmc_getnosims,cbgcmc_displaystats,cbgcmc_displaysimparams , &
      cbgcmc_beginsim

  ! The tag marking the beginning of the CBGCMC section
  Character(len=strLen), Parameter    :: default_cbgcmc_tag = &
      "CBGCMC Information"
  ! The tag marking acceleration tag
  Character(len=strLen), Parameter    :: cbgc_accel_tag = &
      "CBGCMC Acceleration Info"

  Type CBGCMC_Params
    ! cbgcmcsorbs are NOT indexed by sorbate type.  They go from 1-nsorbs
    Type(CBGCMC_Move_Params), Dimension(:), Pointer   :: chain_sorbs 
    Type(CBGCMC_Move_Params), Dimension(:), Pointer   :: rigid_sorbs !we can 
                                       ! also have methane Type of molecules 
    Character(len=strLen)  :: eostag
    Type(EOS_Models)       :: eosparams
    Real(kind=RDbl)        :: tk        ! temperature in Kelvin
    Real(kind=RDbl)        :: max_nrg   ! maximum allowed energy 
    Integer    :: nsorbs
    Integer    :: nsims            ! no. of simulations
    Integer    :: niterations
    Type(Smap_Params), Pointer :: smap
    Character(len=strLen)      :: smapname
    Real(kind=RDbl),Dimension(:,:),Pointer :: noncoul,coul,intra
    Type(SimCell_Params)       :: scell
    Type(Statistics),Dimension(:,:),Pointer       :: intrastats
    Type(Statistics)          :: totNrg
    Type(CBGC_Accelerator),Pointer    :: accel
  End Type CBGCMC_Params

Contains
  !-------------------------------------------------------------------
  ! Initializes the various CBGCMC parameters from the control file
  ! Requires: cbgcparams -- parameters for the GCMC simulation
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           ctrl_filename -- the name of the control file
  !           opt_cbgcmctag -- optional Tag to look for in control file
  !-------------------------------------------------------------------
  Subroutine cbgcmc_init(cbgcmcparams, sorbates, simcell, ctrl_filename, &
      opt_cbgcmctag)
    Type(CBGCMC_Params)         :: cbgcmcparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(inout)  :: simcell
    Character(*), Intent(in)  :: ctrl_filename
    Character(*), Optional, Intent(in) :: opt_cbgcmctag

    Character(len=strLen)     :: eostag, tag, line
    Integer   :: unitno, cbgcmclineno, nsorbs, nmtypes, error, i, j,blocksize
    Integer   :: disk_samp_num, sphere_samp_num, NpsiSS, NthetaSS, NphiDS
    Real(kind=RDbl)     :: volume,pot
    Character(len=strLen) :: tempstr
    Character(len=2*strLen) :: errormsg
    Logical :: fast, pcalcflag, mapflag

    If (Present(opt_cbgcmctag)) Then
      tag = opt_cbgcmctag
    Else
      tag = default_cbgcmc_tag
    End If

    !** Open the ctrl_file if it is not opened
    unitno = file_open(ctrl_filename)
    Rewind(unitno)
    
    !** Find the CBGCMC section
    Write(0,*) "Looking for CBGMC section in ctrlfile"
    cbgcmclineno = filesrchstr(unitno, tag, line)
    If (cbgcmclineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif
    Write(0,*) "Going to initialize CBGMC section "
    Read(unitno, *) cbgcmcparams%niterations   
    Read(unitno, *) cbgcmcparams%tk         ! Temperature
    Read(unitno, '(a)') eostag             
    cbgcmcparams%eostag = stripcmnt(eostag)
    Read(unitno, *) cbgcmcparams%nsims
    Read(unitno, *) blocksize               ! block size for stats
    Read(unitno, *) nsorbs                  ! number of chain molecules
    cbgcmcparams%nsorbs = nsorbs

    ! number of trials in disk sampling, and sphere sampling 
    Read(unitno,*) disk_samp_num, sphere_samp_num
    NpsiSS=Int(Sqrt(sphere_samp_num*one))+1
    NthetaSS=NpsiSS
    Write(*,*) "Actual number of sphere sampling points used : ", &
        NpsiSS**2
    NPhiDS=disk_samp_num

    !** Read the maximum energy to be used during cbgcmc insertion/deletions
    !** If the enrgy goes above this , attempt will be stopped
    Read(unitno,*) cbgcmcparams%max_nrg

    ! some sort og hacking here, but no other way out..
    Call brmoves_initTrialArrays( NpsiSS, NthetaSS, NPhiDS)

    !** read the blankline at the end
    errormsg = 'must be a blank line before species parameters in CBGCMC &
        & section'
    Call readblank(unitno,__FILE__,__LINE__,errormsg)

    !** Allocate the cbgcmcparams%chain_sorbs
    Allocate(cbgcmcparams%chain_sorbs(nsorbs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Set the CBGCMC move parameter information that should be the same
    !** for all moveable sorbates such as temperature and no. of pts
    Do i=1, nsorbs
      cbgcmcparams%chain_sorbs(i)%tk        = cbgcmcparams%tk
      cbgcmcparams%chain_sorbs(i)%npts      = cbgcmcparams%nsims
      cbgcmcparams%chain_sorbs(i)%blocksize = blocksize
      cbgcmcparams%chain_sorbs(i)%max_nrg   = cbgcmcparams%max_nrg
    End Do

    !** The rest of the stuff will be read and initialized by the 
    !** cbgcmcmove init
    Write(*,*) " Going to initialize individual cbgcmc sections"
    Call cbgcmoves_init(cbgcmcparams%chain_sorbs, sorbates, simcell, &
        ctrl_filename)
    Write(*,*) " Finished initializing individual cbgcmc sections"

    !** Initialize the parameters for the chosen equation of state 
    Rewind(unitno)
    Call eos_init(cbgcmcparams%eosparams, ctrl_filename, cbgcmcparams%eostag, &
        cbgcmclineno)
   
    !** Generate the fugacitylist 
    cbgcmcparams%scell=simcell
    volume = simcell_getvolume(simcell)
    Call cbgcmc_BranchGenfuglist(sorbates, cbgcmcparams, volume)
    
    nmtypes=molecules_getnsorbs()
    Allocate(cbgcmcparams%noncoul(nmtypes,nmtypes),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"noncoul")
    Allocate(cbgcmcparams%coul(nmtypes,nmtypes),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"coul")
    Allocate(cbgcmcparams%intra(nmtypes,NO_OF_INTRA_POTS),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"intra")
    Allocate(cbgcmcparams%intrastats(nmtypes,NO_OF_INTRA_POTS),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"intrastats")
    Do i=1,nmtypes
      Do j=1,NO_OF_INTRA_POTS
        tempstr="# intra-"
        If (j==STRETCH_INDEX) tempstr=Trim(tempstr)//"STRH"
        If (j==BENDING_INDEX) tempstr=Trim(tempstr)//"BEND"
        If (j==TORSION_INDEX) tempstr=Trim(tempstr)//"TORS"
        If (j==TOTAL_INDEX) tempstr=Trim(tempstr)//"TOTL"
        If (j==INTRAPAIR_INDEX) tempstr=Trim(tempstr)//"INTR"
        
        Call stats_init(cbgcmcparams%intrastats(i,j),tempstr,1000, &
            .False.,"f14.4")
      End Do
    End Do
    Call stats_init(cbgcmcparams%totnrg,"TotalBrNrg",1000,.False.)
    
    !** Update all nrgs based on initial configs
    Call forcefield_initnrgs(sorbates, simcell)
    Call molecules_displaynrg(6)

    Call cbgcmc_initAccel(cbgcmcparams, sorbates, ctrl_filename)

  End Subroutine cbgcmc_init


  !------------------------------------------------
  ! If needed , initializes the accelerator
  !------------------------------------------------
  Subroutine cbgcmc_initAccel(params, sorbates, ctrlfile)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(CBGCMC_Params), Intent(inout)  :: params
    Character(len=strLen), Intent(in)   :: ctrlfile
    
    Integer :: unitno, accel_line, error
    Character(len=strLen) :: line
    unitno=file_open(ctrlfile)
    accel_line= filesrchstr(unitno, cbgc_accel_tag, line)
    
    If (accel_line==0) Then
      !** no acceleration asked for
      Nullify(params%accel)
    Else
      Allocate(params%accel, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call cbgcaccel_init(params%accel, params%chain_sorbs, sorbates, unitno)
    Endif
  End Subroutine cbgcmc_initAccel


  
  !------------------------------------------------
  ! Generate the fugacity of Ideal Chains 
  ! Need to evaluate the Z- configuration integral of the chain here
  !------------------------------------------------
  Subroutine cbgcmc_BranchGenfuglist(sorbates, cbgcmcparams, volume)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(CBGCMC_Params), Intent(inout)  :: cbgcmcparams
    Real(kind=RDbl), Intent(in)       :: volume

    Integer    :: i, j, k, sorbtype, found
    Real(kind=RDbl)   :: tk, pp, fug, B, sivolume, mass, murti, Lambda
    Real(kind=RDbl)   :: energy,ratio,Zig

    !** Set the temperature and the pressure of the different
    !** species
    tk = cbgcmcparams%tk
    Call eos_settemperature(cbgcmcparams%eosparams, tk)

    !** Get the fugacity
    Do i=1, cbgcmcparams%nsims
      ! Set the pressure of each component
      Do j=1, cbgcmcparams%nsorbs
        sorbtype = cbgcmcparams%chain_sorbs(j)%sorbtype
        pp = cbgcmcparams%chain_sorbs(j)%fuglist(i)%pressure
        ! pp is in kPa
        Call eos_setpressure(cbgcmcparams%eosparams, sorbtype, pp)
      End Do

      ! Get the fugacity of each component
      Do j=1, cbgcmcparams%nsorbs

        !** This is not fugacity, its ideal fugacity
       sorbtype = cbgcmcparams%chain_sorbs(j)%sorbtype
       fug = eos_getfugacity(cbgcmcparams%eosparams, sorbtype) ! fug [=] kPa
       cbgcmcparams%chain_sorbs(j)%fuglist(i)%fugacity = fug

       ! Get the excess chemical potential (B in Adams Notation)
       sivolume = volume*1.0e-30*Nav ! convert to m^3/mole
       B = Log(fug*1.0e3*sivolume/(Rgas*tk))

       !** This is Z/Omega, got from ideal parameters
       Zig = eos_getConfInteg(cbgcmcparams%eosparams, sorbtype)

       ratio=log(Zig) 
       !** we need (PV)/(RTZ)
       B = B - ratio
        Write(0,'(2a,i4,a,f16.2,6f16.10)') __FILE__,": ",__LINE__, &
        " B ", ratio,B
        cbgcmcparams%chain_sorbs(j)%fuglist(i)%B = B

        ! Get mu the chemical potential
        ! Calculate the DeBroglie wavelength
        mass = molecules_getmass(sorbtype)
        mass = mass*1.0e-3              ! convert to kg
        Lambda = hplanck/Sqrt(twopi*mass/Nav*Rgas/Nav*tk)
        murti = (B - Log(sivolume/Nav/Lambda**3))
        cbgcmcparams%chain_sorbs(j)%fuglist(i)%murti = murti
      End Do
    End Do
  End Subroutine cbgcmc_BranchGenfuglist


  !--------------------------------------------
  ! Gets the total no. of simulations
  !--------------------------------------------
  Integer Function cbgcmc_getnosims(cbgcmcparams)
    Type(CBGCMC_Params), Intent(in)  :: cbgcmcparams
    cbgcmc_getnosims = cbgcmcparams%nsims
  End Function cbgcmc_getnosims
  
  !------------------------------------------------------------------
  ! Calls one of the dosim routines
  !------------------------------------------------------------------
  Subroutine cbgcmc_dosim(cbgcmcparams, sorbates, simcell, simno)
    Type(CBGCMC_Params), Intent(inout)  :: cbgcmcparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(in)  :: simcell
    Integer, Intent(in) :: simno
    If (Associated(cbgcmcparams%accel)) Then
      Call cbgcmc_dosimAccel(cbgcmcparams, sorbates, simcell, simno)
    Else
      Call cbgcmc_dosimNormal(cbgcmcparams, sorbates, simcell, simno)
    Endif
  End Subroutine cbgcmc_dosim

  !---------------------------------------------------------------------------
  ! Performs "ninterations" of the specified CBGCMC simulation.
  ! We need to ensure that insertions and deletions are done with equal 
  ! probability for each species type. 
  ! This is shaji's home made algorithm, and should be checked 
  ! thorouhly before being used anywhere else 
  ! ( it allows for number of insertions and deleteions to be different )
  ! 1) Picks a species with equal probability.
  ! 2) Decide whether we want to attempt an
  !    insertion, deletion, translation or rotation.
  ! Enforcement of these constraints are handled in cbgcmoves
  ! Requires: cbgcparams -- parameters for the GCMC simulation
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           simno -- the simulation number
  !---------------------------------------------------------------------------
  Subroutine cbgcmc_dosimNormal(cbgcmcparams, sorbates, simcell, simno)
    Type(CBGCMC_Params), Intent(inout)  :: cbgcmcparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(in)  :: simcell
    Integer, Intent(in) :: simno
    Integer             :: iter, nmoles, moveno, sorbtype
    Integer             :: i, j, nsorbs,sorbno,iseed
    Real(kind=RDbl)     :: fugacity
    Character(len=strLen)  :: moveType

    Logical :: mapflag, fast , success
    Integer :: k 
    Real    :: noncoulpot   

    !** Do the iterations at the set conditions
    Do iter = 1, cbgcmcparams%niterations
      ! Make sure that cbgcmcsorbs is assigned, nsorbs should be > 0
      If (cbgcmcparams%nsorbs>0) Then
        sorbtype = Int(rranf()*cbgcmcparams%nsorbs) + 1
        !** Pick a move type randomly and execute it
        success =  cbgcmoves_move(cbgcmcparams%chain_sorbs(sorbtype), &
            sorbates, simcell,simno)

        !** Do updating of quantities to be averaged
        Call cbgcmoves_updatenmolecs(cbgcmcparams%chain_sorbs,sorbates)
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    End Do

  End Subroutine cbgcmc_dosimNormal


  !---------------------------------------------------------------------------
  ! WITH ACCELERATOR
  ! Performs "ninterations" of the specified CBGCMC simulation.
  ! We need to ensure that insertions and deletions are done with equal 
  ! probability for each species type. 
  ! This is shaji's home made algorithm, and should be checked 
  ! thorouhly before being used anywhere else 
  ! ( it allows for number of insertions and deleteions to be different )
  ! 1) Picks a species with equal probability.
  ! 2) Decide whether we want to attempt an
  !    insertion, deletion, translation or rotation.
  ! Enforcement of these constraints are handled in cbgcmoves
  ! Requires: cbgcparams -- parameters for the GCMC simulation
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           simno -- the simulation number
  !---------------------------------------------------------------------------
  Subroutine cbgcmc_dosimAccel(cbgcmcparams, sorbates, simcell, simno)
    Type(CBGCMC_Params), Intent(inout)  :: cbgcmcparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(in)  :: simcell
    Integer, Intent(in) :: simno
    Integer             :: iter, nmoles, moveno, sorbtype
    Integer             :: nsorbs,sorbno,iseed
    Real(kind=RDbl)     :: fugacity
    Character(len=strLen)  :: moveType

    Logical :: mapflag, fast , success
    Integer :: k ,i
    Real    :: randnum

    !** Do the iterations at the set conditions
    Do iter = 1, cbgcmcparams%niterations
      ! Make sure that cbgcmcsorbs is assigned, nsorbs should be > 0
      sorbtype = Int(rranf()*cbgcmcparams%nsorbs) + 1
!!$      If (rranf()<cbgcmcparams%accel%gc_factor)  Then

        If (cbgcmcparams%nsorbs>0) Then
          !** Prepare the accelerator
          Call cbgcaccel_dosim(cbgcmcparams%accel, sorbates,simcell,&
              cbgcmcparams%chain_sorbs, sorbtype, simno)

          !** Do updating of quantities to be averaged
          Call cbgcmoves_updatenmolecs(cbgcmcparams%chain_sorbs,sorbates)

          !** do some nvt here too, but use all the sorbtypes
          !** This might not be very microscopically reversible
          Do i=1,cbgcmcparams%accel%sorbs(sorbtype)%nvt_simnum
            sorbtype = Int(rranf()*cbgcmcparams%nsorbs) + 1
            success =  cbgcmoves_nvtmove(cbgcmcparams%chain_sorbs(sorbtype), &
                sorbates, simcell,simno)
          End Do
          
        Else

          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
          
        Endif
!!$      Else
!!$        !** Pick a move type randomly and execute it
!!$        success =  cbgcmoves_nvtmove(cbgcmcparams%chain_sorbs(sorbtype), &
!!$            sorbates, simcell,simno)
!!$      Endif
    End Do
    
  End Subroutine cbgcmc_dosimAccel
  
  
  !-------------------------------------------------------------------------
  ! Handles the beginning of CBGCMC simulations
  ! 1) Write feedback to screen
  ! Requires: cbgcmcparams -- CBGCMC simulation parameters
  !           simno -- simulation number
  !           sorbates -- species data structure  
  !           simcell -- the simulation cell information
  !           indent -- no. of spaces from the left margin
  !-------------------------------------------------------------------------
  Subroutine cbgcmc_beginsim(cbgcmcparams,simno,sorbates,simcell,indent)
    Type(CBGCMC_Params), Intent(InOut)         :: cbgcmcparams
    Integer, Intent(In)                      :: simno
    Type(AtMolCoords), Dimension(:), Pointer :: sorbates
    Type(SimCell_Params), Intent(In)         :: simcell
    Integer, Intent(In)                      :: indent

    Integer                           :: sorbno
    Character(len=indent)             :: blank

    blank = Repeat(' ',indent)    

    !** Reset the counters for each species
    Do sorbno = 1,cbgcmcparams%nsorbs
      Call cbgcmoves_beginsim(cbgcmcparams%chain_sorbs(sorbno), &
          sorbates, simno)
    End Do
    
  End Subroutine cbgcmc_beginsim

  !-------------------------------------------------------------------------
  ! Handles the ending of CBGCMC simulations
  ! 1) Check the accumulated and fresh species-species and intra energies
  ! 2) Write feedback to screen
  ! 3) Optionally write to a restartfile (if stopfile /= '')
  ! 4) halt the program only if desired
  ! Step (4) can be useful for getting a starting file for subsequent MD runs.
  ! Requires: cbgcmcparams -- CBGCMC simulation parameters
  !           simno -- simulation number
  !           sorbates -- species data structure  
  !           simcell -- the simulation cell information
  !           stopfile -- name of restart file
  !           indent -- no. of spaces from the left margin
  !           unit -- display unit number
  !           stopflag -- optional flag to signal program halt
  !-------------------------------------------------------------------------
  Subroutine cbgcmc_endsim(cbgcmcparams,simno,sorbates,simcell,stopfile, &
      indent, unit,stopflag)
    Type(CBGCMC_Params), Intent(In)          :: cbgcmcparams
    Integer, Intent(In)                      :: simno
    Type(AtMolCoords), Dimension(:), Pointer :: sorbates
    Type(SimCell_Params), Intent(In)         :: simcell
    Character(*), Intent(In)                 :: stopfile
    Integer, Intent(In)                      :: indent,unit
    Logical, Intent(In), Optional            :: stopflag 

    Integer                           :: spc
    Logical                           :: killflag
    Real(kind = RDbl)                 :: rmsdev
    Character(len=indent)             :: blank

    blank = Repeat(' ',indent)    

    If (Present(stopflag)) Then
      killflag = stopflag
    Else
      killflag = .False.
    End If

    !** check the updated and fresh energies for agreement
    rmsdev = forcefield_checknrgs(sorbates,simcell,6)
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

    Call molecules_displaynrg(6)

    !** Feedback 
    Write(unit,'(2a)') blank,dashedline
    Write(unit,'(3a,4x,a,e9.3)') blank,'Ending GCMC simulation number ', &
        Trim(int2str(simno)),'Nrg dev check: ',rmsdev

    Do spc = 1,cbgcmcparams%nsorbs
      Call cbgcmoves_dispspcstats(cbgcmcparams%chain_sorbs(spc),simno, &
          indent+2,unit)
    End Do
    
    If (Trim(stopfile) /= '') Then
      Write(unit,*) "Writing a final restart file into  ", Trim(stopfile)
      Call config_writerestartfile(sorbates, Trim(stopfile))
    End If

    Write(unit,'(a)') dashedline    
    
    !** Halts the program if desired
    If (killflag) Then
      Write(0,*) "Exiting after writing the restartfile"
      Stop
    Endif
    
  End Subroutine cbgcmc_endsim

  !-----------------------------------------------------------
  ! Display the CBCGCMC statistics for the current simulation
  ! Requires: cbgcmcparams -- CBGCMC simulation parameters
  !           sorbates -- species data structure  
  !           nspc -- no. of spaces from the left margin
  !           optunit -- optional unit number for display
  !-----------------------------------------------------------
  Subroutine cbgcmc_displaystats(cbgcmcparams, sorbates, simno, nspc, optunit)
    Type(CBGCMC_Params), Intent(inout) :: cbgcmcparams
    Integer, Intent(in)   :: simno, nspc
    Integer, Optional, Intent(in) :: optunit
    Type(AtMolCoords), Dimension(:), Intent(Inout) :: sorbates
    
    Character(len=nspc) :: spc
    Character(len=strLen) :: restartfile, configfile
    Integer    :: i, unitno, funit
    
    spc = Repeat(' ', nspc)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "The CBGCMC Stats:"
    Call cbgcmoves_displaystats(cbgcmcparams%chain_sorbs, simno, nspc+4, &
        unitno)
    If (Trim(genparams%displaymode)=="VERBOSE") Then
      Call cbgcmc_displayNrgs(cbgcmcparams, sorbates,nspc+4, unitno)
    Endif
  End Subroutine cbgcmc_displaystats

  !----------------------------------------------------------------------
  ! Dispalys Nrg's after estimating them freshly
  ! and also after obtaining  from molecules
  !----------------------------------------------------------------------
  Subroutine cbgcmc_displayNrgs(params,sorbates,nspc,unitno)
    Type(CBGCMC_Params), Intent(inout) :: params
    Integer, Intent(in)   :: nspc,unitno
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: sorbates

    Character(len=nspc) :: spc
    Logical           :: mapflag,pcalcFlag,fastFlag    
    Real(kind=RDbl)   :: pot,totalnrg
 !SDEBUG
 !   Real(kind=RDbl) ::totalcoulnrg,totalnoncoulnrg,totalintranrg
 !SDEBUG
    Integer :: i,mtype1,mtype2,nmtypes,nmoles,nmoles1,nmoles2,totalmoles
    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do    

    fastFlag=.true.
    pcalcFlag=.true.
    Call forcefield_getssint(sorbates,params%scell,fastFlag,pot,mapflag,&
        params%noncoul,params%coul,params%intra,pcalcFlag)
    
    nmtypes=molecules_getnsorbs()
    !** display coul and non coul
    
    Write(unitno,*) spc, dashedline
    Write(unitno,*) "Energies calculated Fresh vs Stored"
    totalmoles=0
    totalnrg=zero

    Write(unitno,*) spc//" Pairs        Fresh-NCOUL ,   Fresh-COUL"
    Do mtype1=1,nmtypes
      
      nmoles1=config_getnmoles(sorbates,mtype1)
      If (config_isfixed(sorbates(mtype1))) nmoles1 = 0
      totalmoles=totalmoles+nmoles1
      
      Do mtype2=mtype1,nmtypes
        nmoles2=config_getnmoles(sorbates,mtype2)
        If (config_isfixed(sorbates(mtype2))) nmoles2 = 0
        
        nmoles=nmoles1+nmoles2
        If (mtype2==mtype1) nmoles=nmoles1
        If (nmoles==0) Cycle
        Write(*,'(a,i3,a,i3,3x,2f15.4)') spc,mtype1,"-",mtype2, &
            params%noncoul(mtype1,mtype2)*scalepe/nmoles, &
            params%coul(mtype1,mtype2)*scalepe/nmoles

        totalnrg=totalnrg+params%noncoul(mtype1,mtype2)*scalepe+&
            params%coul(mtype1,mtype2)*scalepe
 !       totalcoulnrg=totalcoulnrg+params%coul(mtype1,mtype2)*scalepe
 !       totalnoncoulnrg=totalnoncoulnrg+params%noncoul(mtype1,mtype2)*scalepe

        Write(unitno,'(a20,4f10.4)') spc//"Stored-NCOUL : ", &
            molecules_getnoncoul(mtype1,mtype2,'inst')/nmoles, &
            molecules_getnoncoul(mtype1,mtype2,'cavg')/nmoles, &
            molecules_getnoncoul(mtype1,mtype2,'block')/nmoles, &
            molecules_getnoncoul(mtype1,mtype2,'std')/nmoles
        Write(unitno,'(a20,4f10.4)') spc//"Stored-COUL : ", &
            molecules_getcoul(mtype1,mtype2,'inst')/nmoles, &
            molecules_getcoul(mtype1,mtype2,'cavg')/nmoles, &
            molecules_getcoul(mtype1,mtype2,'block')/nmoles, &
            molecules_getcoul(mtype1,mtype2,'std')/nmoles
        
        Write(unitno,*) 
      End Do
    End Do

    Do mtype1=1,nmtypes
      If (config_isfixed(sorbates(mtype1))) Cycle
      nmoles=config_getnmoles(sorbates,mtype1)
      Write(*,'(a,i4,a,i4)') "Intra nrg for molecule: ", &
          mtype1, "Nmoles ", nmoles
      If (nmoles==0) Cycle

      Call stats_update(params%intrastats(mtype1,STRETCH_INDEX),&
          params%intra(mtype1,STRETCH_INDEX)*scalepe/nmoles)
      Call stats_display(params%intrastats(mtype1,STRETCH_INDEX),&
          2,unitno)

      Call stats_update(params%intrastats(mtype1,BENDING_INDEX),&
          params%intra(mtype1,BENDING_INDEX)*scalepe/nmoles)
      Call stats_display(params%intrastats(mtype1,BENDING_INDEX),&
          2,unitno)

      Call stats_update(params%intrastats(mtype1,TORSION_INDEX),&
          params%intra(mtype1,TORSION_INDEX)*scalepe/nmoles)
      Call stats_display(params%intrastats(mtype1,TORSION_INDEX),&
          2,unitno)

      Call stats_update(params%intrastats(mtype1,INTRAPAIR_INDEX),&
          params%intra(mtype1,INTRAPAIR_INDEX)*scalepe/nmoles)
      Call stats_display(params%intrastats(mtype1,INTRAPAIR_INDEX),&
          2,unitno)

      Call stats_update(params%intrastats(mtype1,TOTAL_INDEX),&
          params%intra(mtype1,TOTAL_INDEX)*scalepe/nmoles)
      Call stats_display(params%intrastats(mtype1,TOTAL_INDEX),&
          2,unitno)
      
      Write(unitno,'(a14,4f14.4)') "Stored-Intra : ", &
          molecules_getintranrg(mtype1,'inst')/nmoles, &
          molecules_getintranrg(mtype1,'cavg')/nmoles, &
          molecules_getintranrg(mtype1,'block')/nmoles, &
          molecules_getintranrg(mtype1,'std')/nmoles
      Write(unitno,'(a)') dashedline
      Write(unitno,*) 
      

      totalnrg=totalnrg+params%intra(mtype1,TOTAL_INDEX)*scalepe
  !    totalintranrg=totalintranrg+params%intra(mtype1,TOTAL_INDEX)*scalepe
    End do


!      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!    Write(*,*) totalcoulnrg
!    Write(*,*) totalnoncoulnrg
!    Write(*,*) totalintranrg
!    Write(*,*) totalnrg
!    Write(*,*) pot*scalepe
!      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    If (totalmoles>0) Then
      Write(unitno,'(a,f18.5)') "total avg nrg calculated with current &
          &config: ",totalnrg/totalmoles
      Write(unitno,'(a,i4)') " Number of moles : ", totalmoles
!      Write(unitno,'(a,e18.5)') "pot nrg calculated with current &
!          &config: ", pot/totalmoles
    Else
      Write(unitno,'(a,i4)') " Number of moles : ", totalmoles
    Endif

  End Subroutine cbgcmc_displayNrgs
  
  !------------------------------------------------------------------
  ! Display the cbgcmc simulation parameters for simulation no. "simno"
  ! The display is sent to the optional unit no. "optunit".  "nspc"
  ! is the no. of spaces to leave from the left margin
  !-----------------------------------------------------------------
  Subroutine cbgcmc_displaysimparams(cbgcmcparams, simno, nspc, optunit)
    Type(CBGCMC_Params), Intent(in) :: cbgcmcparams

    Integer, Intent(in) :: nspc
    Integer, Intent(in) :: simno
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Character(len=strLen) :: restartfile, configfile
    Integer    :: i, unitno, funit
    
    spc = Repeat(' ',nspc)
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    !** Get the configuration and restart file names
    Call file_gettype(d_con_file, configfile, funit)
    Call file_gettype(d_res_file, restartfile,funit)

    !** Write the simulation no.
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "The CBGCMC Simulation Parameters:"
    Write(unitno, '(a,2x,a,i4)') spc, "Simulation Number   : ", simno
    Write(unitno, '(a,2x,2a) ') spc, "Configuration file  : ", &
        Trim(configfile)
    Write(unitno, '(a,2x,2a) ') spc, "Restart file        : ", &
        Trim(restartfile)

    Call cbgcmoves_displaysimparams(cbgcmcparams%chain_sorbs,simno,nspc+2,unitno)
  End Subroutine cbgcmc_displaysimparams


  !---------------------------------------------
  ! Display the cbgcmc initialization parameters
  !---------------------------------------------
  Subroutine cbgcmc_initdisplay(cbgcmcparams, nspc, optunit)
    Type(CBGCMC_Params), Intent(in) :: cbgcmcparams
    Integer, Intent(in) :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer    :: i, unitno
    
    spc = Repeat(' ',nspc)
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "The CBGCMC Parameters Section:"
    Write(unitno, '(a,2x,a,i6)') spc,"No. of iterations  : ", &
        cbgcmcparams%niterations
    Write(unitno, '(a,2x,a,i6)') spc, "No. of simulations : ", cbgcmcparams%nsims
    Write(unitno, '(a,2x,a,f8.3)') spc, "Temperature(K)     : ", cbgcmcparams%tk
    Write(unitno, '(a,2x,a,i6)') spc, "No. of sorbates    : ",cbgcmcparams%nsorbs
    Call cbgcmoves_display(cbgcmcparams%chain_sorbs, nspc+2, unitno)
  End Subroutine cbgcmc_initdisplay

End Module cbgcmc







