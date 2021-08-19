!-------------------------------------------------------------------
! This module contains the different data structures
! and routines for doing the hybrid monte carlo integration moves
! It has an mcmove object inside, which is used for conducting the monte 
! carlo move
!-------------------------------------------------------------------

Module hybridmoves
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl, strLen, MAX_EXCLSITES,MAX_SORBS,kcalmole_kb ,&
      kjmole_kb, zero, NO_OF_INTRA_POTS,TOTAL_INDEX,STRETCH_INDEX, &
      BENDING_INDEX, TORSION_INDEX,INTRAPAIR_INDEX, lstrLen, &
      CONSTRAINT_INDEX,scalepe, scaleke, dashedline, INCR_SIZE, MAX_ATOMS
  Use file, Only: file_getunit, file_gettype, file_open
  Use atom, Only: atom_getmass
  Use config, Only: Atmolcoords,config_getnmoles,config_getnatoms,&
      config_allocfields,config_checkandincr,config_setnmoles,&
      config_copymolec,config_kineticEnergy, config_isfixed, config_reverseVels
  Use datafile, Only: CONFILE
  Use gear6, Only: gear6_extraenergy
  Use gcmodels, Only: GeneralizedCoords,gcmodels_getcom
  Use interact, Only: Interaction_Model, interact_checknrgs
  Use subinteract, Only: Subset_Interactions, subinteract_int, &
      subinteract_init, subinteract_update
  Use mcmoves, Only: MC_Move_Params,mcmoves_init, &
      mcmoves_perturb, mcmoves_getbasictag, mcmoves_display, &
      mcmoves_stats, mcmoves_resetstats, mcmoves_checkandincr
  Use molecule, Only: MolecularParams
  Use molecules, Only: molecules_getnatoms, &
      molecules_getnsorbs, molecules_getnthatom, &
      molecules_gettype, molecules_getgcmodeltype, molecules_getatype, &
      molecules_getpointer, molecules_name
  Use random, Only: rranf
  Use randomvec, Only: randomvec_gaussian
  Use stats, Only: Statistics,stats_setvalue,stats_getvalue, &
      stats_update, stats_display
  Use simcell, Only: SimCell_Params, simcell_pbc
  Use storetop, Only: storetop_sum, storetop_display
  Use storebase, Only:  storebase_copy  
  Use utils, Only: isfileopen, filesrchstr, toupper, split, &
      stripcmnt, toint , allocErrDisplay, cleanstring
  Use velverlet, Only: velverlet_extraenergy
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)

  Implicit None
  Save

  Private
  Public :: Hybrid_Move_Params,hybridmoves_init,hybridmoves_integrate, &
      hybridmoves_displaystats

  Type Hybrid_Move_Params
    Type(MC_Move_Params)                           :: mcmove
    Integer                                        :: npts
    Integer                                        :: blocksize 
    !** Thermophysical Properties
    Real(kind=RDbl),Dimension(:), Pointer          :: templist
    Real(kind=RDbl),Dimension(:), Pointer          :: rti !* 1.0/(Rgas*tk)
    
    !** Subset interaction storage, for the whole system 
    !** I had to make this an array so I can pass to mcmoves_perturb
    Type(Subset_Interactions), Dimension(1) :: subint
    
    !** Type of velocity generation 
    Logical :: boltzmann_vels, keep_vels

  End Type Hybrid_Move_Params

Contains

  !----------------------------------------------------------------
  ! Initializes the various Hybridmoves parameters from the control
  ! file "ctrl_filename".  
  ! allocates all the dummy arrays in Hybrid_Move_Params
  !----------------------------------------------------------------
  Subroutine hybridmoves_init(params, imodel, species, scell, templist, &
      ctrl_filename,auxmv)
    Type(Hybrid_Move_Params),Intent(inout)           :: params
    Type(Interaction_Model), Intent(in)              :: imodel 
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: species
    Type(Simcell_Params), Intent(in)                 :: scell
    Character(*), Intent(in)                         :: ctrl_filename
    Real(kind=RDbl) ,Dimension(:), Intent(in)        :: templist
    Type(AuxMoveObjects), Pointer           :: auxmv 

    Integer                   :: unitno, nsorbs, i, error,nfields,lineno
    Integer                   :: nmoles, j,nspecies, npts
    Character(len=strLen)             :: sorbname,blankline,ens_string
    Character(len=2*strLen)           :: text,line
    Character(len=strLen),Dimension(MAX_SORBS) :: fields

    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) "Initializing hybridmoves section"

    !** unit.no of control file, 
    !** Its already Open and pointing to the hybrid moves section
    unitno = file_getunit(ctrl_filename)

    !** hybridmoves are for the whole species structure
    !** make sure that we are reading the correct parameters
    Read(unitno, *) sorbname
    sorbname=Trim(stripcmnt(sorbname))
    If(sorbname=="ALLSORBS") Then

      ! ** figure out how to get the velocities prior to the integration
      Read(unitno,'(a)') text
      Select Case(cleanstring(text))
      Case ("BOLTZMANN") 
        params%boltzmann_vels=.True.
        params%keep_vels=.False.
        ens_string="NVT"
      Case ("KEEP_VELOCITIES")
        params%boltzmann_vels=.False.
        params%keep_vels=.True.
        ens_string="NVE_VELS"
      Case default
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End Select

      npts=params%npts

      Allocate( params%templist(npts) , params%rti(npts), STAT=error )
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      params%templist(1:npts) = templist(1:npts)

      !** all force/nrg calculations are going to be done in Kcal/mole 
      !** kcalmole_kb is the "R" value in the units : kcal/mol/K
      params%rti(1:npts) = 1.0_RDbl/( kcalmole_kb * params%templist(1:npts) )

      !** initialize the subset interaction structure
      Call subinteract_init(params%subint(1),imodel,'System_System', &
          'HMC',(/0,0,0/),(/0,0,0/))

      !** get the move-params used for the integartion , should be an MD move
      Call mcmoves_init(params%mcmove, ctrl_filename, scell, species, &
          ens_string, auxmv)

    Else
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,&
          "could not find sorbname",sorbname
      Stop
    End If

  End Subroutine hybridmoves_init

  !-------------------------------------------------------------------
  ! Decides which hybridmc integration move to make and performs it.
  ! Requires:  hybridparams -- Hybrid move structure
  !            species -- configuration data structure
  !            scell -- simulation cell
  !            simno -- simulation number
  !-------------------------------------------------------------------
  Subroutine hybridmoves_integrate(hybridparams, species, scell, &
      simno)
    Type(Hybrid_Move_Params), Intent(InOut)        :: hybridparams
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(Simcell_Params), Intent(In)               :: scell
    Integer, Intent(In)                            :: simno

    !** Identify the move type
    If (hybridparams%boltzmann_vels) Then
      Call hybridmoves_integrate_boltz(hybridparams, species, &
          scell, simno)
    Else If (hybridparams%keep_vels) Then
      Call hybridmoves_integrate_vels(hybridparams, species, scell, &
          simno)

    Else
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

  End Subroutine hybridmoves_integrate


  !-------------------------------------------------------------------
  ! Does an Hybird Monte Carlo Integration move. 
  ! All moves from here are accepted/ejected by mcmoves. These are NVE moves 
  ! Requires:  params -- Hybrid move structure
  !            species -- configuration data structure
  !            scell -- simulation cell
  !            simno -- simulation number
  !-------------------------------------------------------------------
  Subroutine hybridmoves_integrate_boltz(params, species, &
      scell, simno)
    Type(Hybrid_Move_Params), Intent(inout)        :: params
    Type(AtMolCoords), Dimension(:), Intent(inout) :: species
    Type(Simcell_Params), Intent(in)               :: scell
    Integer, Intent(In)                            :: simno

    Logical                :: nrgCalcFlag
    Integer                :: i, j, nsorbs, nmoles
    Logical                :: fast, success, reCalcFlag, knrgflag
    Real(kind=RDbl)        :: max_nrg, dt, Tk     

    !** HACKed that all interactions are FAST
    fast=.True.
    dt=zero  ! just lazy to write the mcmoves_getfileds routine
    Tk=params%templist(simno)

    nsorbs= molecules_getnsorbs()

    !** make sure all arrays and storage structures are correctly sized
    !** resizes mcmove%subint%temp (which is same as params%subint)
    Call mcmoves_checkandincr(params%mcmove,params%subint(1),species,(/0,0,0/))

    !** assign Boltzman Velocites to species, get the total kin nrg back.
    Call hybridmoves_boltzVel(params,species,Tk)

    !** Initialize the forces , we need this to make sure that the 
    !** on all atoms are set before starting the integration
    recalcFlag=.True.
    max_nrg=1000000.00_RDbl
    success = subinteract_int(params%subint(1),species,scell,fast, &
      recalcFlag,.False.,(/max_nrg, dt, Tk/),(/0,0,0/))

    !** Update the forces and intra nrgs, things which entered the system 
    !** through regular gcmc(insert) moves might not have their 
    !** intra-nrgs, forces etc up to date
    Call subinteract_update(params%subint,(/1/))

    !** copy the results into current_sum
    Call storetop_sum(params%subint(1)%temp_try)
    Call storebase_copy(params%subint(1)%now_sum, &
        params%subint(1)%temp_try%total)

    If (.Not. success) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Interaction evaluation unsuccessful for current config. '
      Write(*,*) "check max_nrg values also, it may be set too low. "
      Stop    
    End If

    !** HMC needs KE in acceptance criterion
    knrgflag=.True.  
    success = mcmoves_perturb(params%mcmove, (/0,0,0/), params%rti(simno), &
        params%subint, species, scell, knrgflag) 

  End Subroutine hybridmoves_integrate_boltz

  !----------------------------------------------------------------------------
  ! Assigns boltzman weighted velocities to each atom of the species
  ! separately. and returns the total Kin Energy of the system.
  !----------------------------------------------------------------------------
  Subroutine hybridmoves_boltzVel(params,species , tk)
    Type(Hybrid_Move_Params), Intent(In) :: params
    Type(AtMolCoords), Intent(InOut), Dimension(:) :: species
    Real(kind=RDbl), Intent(in) :: tk

    Real(kind=RDbl) :: kinNrg

    Integer :: natoms, nmoles, mtype,atype,a,m 
    Real(kind=RDbl) :: avgVelocity
    Real(kind=RDbl),Dimension(MAX_ATOMS) :: sigma,mass

    kinNrg=zero
    
    Do mtype = 1, Size(species)
      If (config_isfixed(species(mtype))) Cycle

      natoms = config_getnatoms(species,mtype)

      !** check whether "sigma" array and "mass" array are big enough
      If (natoms>MAX_ATOMS) Then
        Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,&
            "MAX_ATOMS should be increased"
        Stop
      End If

      !** store the STD Devn in velocity for each atom, in the molecule
      !** and their massses
      Do a = 1, natoms
        atype = molecules_getnthatom( mtype, a )
        mass(a) = atom_getmass(atype)

        !** we need momentum such that it has sigma^2= m K T
        !** "scaleke" rqd for converting Kin Nrg to KJ/mole 
        sigma(a) = Sqrt(  ( kjmole_kb * tk ) / ( scaleke * mass(a) )  )

      End Do
      
      !** actual number of molecules in the sorb'th sorbate
      nmoles = config_getnmoles(species,mtype)

      avgVelocity = zero
      Do m = 1, nmoles
        Do a = 1, natoms

          !** assign random gaussian velocities
          species(mtype)%coords(a,m)%v = &
              randomvec_gaussian( avgVelocity , sigma(a))

          !** update total kinetic energy
          kinNrg = kinNrg + mass(a) *  & 
              ( species(mtype)%coords(a,m)%v * species(mtype)%coords(a,m)%v )
        End Do
      End Do
      
    End Do
    
    !** returns in KJ/mol this is calculated only for debugging purposes
    ! REMOVE LATER
    kinNrg = kinNrg * scaleke / 2

  End Subroutine hybridmoves_boltzVel



  !-------------------------------------------------------------------
  ! Does an Integration move. All moves from here are expected to be 
  ! accepted by mcmoves. These are NVE moves 
  ! Requires:  params -- Hybrid move structure
  !            species -- configuration data structure
  !            scell -- simulation cell
  !            conf_file -- configuration file structure
  !            simno -- simulation number
  !-------------------------------------------------------------------
  Subroutine hybridmoves_integrate_vels(params, species, scell, &
      simno)
    Type(Hybrid_Move_Params), Intent(InOut)        :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(Simcell_Params), Intent(in)               :: scell
    Integer, Intent(In)                            :: simno

    Integer                :: i, j, nsorbs, nmoles,outunit
    Logical                :: recalcFlag, calcAccels, fast, success, knrgflag
    Real(KIND=RDBL)        :: max_nrg, dt, Tk
    Real(kind=RDbl)        :: nrg_devn
    Character(len=lstrLen) :: filename
    
    !** HACKed that all interactions are FAST
    fast=.True.
    dt=zero  ! just lazy to write the mcmoves_getfileds routine
    Tk=zero ! lazy again

    !** make sure all arrays and storage structures are correctly sized
    Call mcmoves_checkandincr(params%mcmove,params%subint(1),species,(/0,0,0/))
    
    !** go through all molecules and reverse their velocities with a 
    !** probability of 0.5
    Call hybridmoves_reverseVels(params, species)

!!$
!!$    filename = 'beforeint.txt'
!!$    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!!$    Write(*,*) 'Writing storage to: ',Trim(filename)
!!$    outunit = file_open(filename,0)
!!$    Call storetop_display(params%subint%temp,.False.,2,outunit)
!!$    Close(unit=outunit)

    !** Initialize the forces , we need this to make sure that the 
    !** on all atoms are set before starting the integration
    recalcFlag=.True.
    max_nrg = 1000000.00_RDbl
    success = subinteract_int(params%subint(1),species,scell,fast, &
      recalcFlag,.False.,(/max_nrg,dt,Tk/),(/0,0,0/))

    !** Update the forces and intra nrgs, things which entered the system 
    !** through regular gcmc(insert) moves might not have their 
    !** intra-nrgs, forces etc up to date
    Call subinteract_update(params%subint,(/1/))

#ifdef DEBUG
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) "BEFORE INTEGRATE"
    nrg_devn = interact_checknrgs(params%subint%imodel,fast,species,scell,&
        8 , 6)
#endif
    If (.Not. success) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Interaction evaluation unsuccessful for current config. '
      Write(*,*) "check max_nrg values also, it may be set too low. "
      Stop    
    End If

    !** copy the results into now_sum
    Call storetop_sum(params%subint(1)%temp_try)
    Call storebase_copy(params%subint(1)%now_sum, &
        params%subint(1)%temp_try%total)

    ! knrg will chnage during these moves , but no need to use it in 
    ! acceptance criterion
    knrgflag=.False.  
    success = mcmoves_perturb(params%mcmove, (/0,0,0/), params%rti(simno), &
        params%subint, species, scell, knrgflag) 

#ifdef DEBUG
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) "AFTER INTEGRATE"
    nrg_devn = interact_checknrgs(params%subint%imodel,fast,species,scell,&
        8 , 6)
#endif

    !** This is nve_vels ensemble. All moves will be accepted in mcmoves.
    !** This check is only to make sure that nothing goes wrong in mcmoves 
    If (.Not.success) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif

  End Subroutine hybridmoves_integrate_vels

  !-------------------------------------------------------------
  ! Display the statistics of the hybridmoves
  !-------------------------------------------------------------
  Subroutine hybridmoves_displaystats(params, species, simno, nspc, optunit)
    Type(Hybrid_Move_Params), Intent(in) :: params
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(in)   :: simno, nspc
    Integer, Optional, Intent(in) :: optunit
    Character(len=nspc) :: spaces
    Integer     :: i, j, unitno
    Real(kind=RDbl)         :: integratio

    spaces = Repeat(" ",nspc)
    unitno = 6  
    If (Present(optunit))       unitno = optunit

    Call mcmoves_stats(params%mcmove, nspc+2, unitno)

  End Subroutine hybridmoves_displaystats

  !--------------------------------------------------------------------
  !Checks and make sures that all coords are big enough 
  !--------------------------------------------------------------------
  Subroutine hybridmoves_reverseVels(params,species)
    Type(Hybrid_Move_Params),Intent(inout)         :: params
    Type(AtMolCoords), Dimension(:), Intent(inout) :: species

    Integer :: nmoles,nsorbs, spc, m

    nsorbs=molecules_getnsorbs()
    Do spc = 1,nsorbs
      If (config_isfixed(species(spc))) Cycle
      nmoles   = config_getnmoles(species, spc)
      Do m=1,nmoles
        If (rranf()>0.5000_RDbl) Then
          Call config_reverseVels(species, spc, m)
        Endif
      End Do
    End Do

  End Subroutine hybridmoves_reverseVels

End Module hybridmoves


!********* OLD ROUTINES SOME OF THEM CAN BE DELETED, SOME HAS TO BE
!********* CHECKED WITH THE NEW CODE AND REINSERTED ABOVE
!!$  !----------------------------------------------------
!!$  ! All energies are stored in molecules, gets it from there, sums it
!!$  ! and returns the total potential energy of the system
!!$  !---------------------------------------------------- 
!!$  Subroutine hybridmoves_getPotNrg(params,PotNrg)
!!$    Type(Hybrid_Move_Params), Intent(in) :: params
!!$    Real(kind=RDbl),Intent(out) :: PotNrg
!!$
!!$    Real(kind=RDbl) :: pe,intra,coul,noncoul
!!$    Integer:: nsorbs,i,j
!!$
!!$    pe=zero
!!$    intra=zero
!!$    nsorbs=molecules_getnsorbs()
!!$
!!$    Do i=1,nsorbs
!!$      If (hybridmoves_checkDummy(params,i)) Cycle
!!$      !** intramolecular potential nrgs strecthing,bending etc..
!!$      intra=intra+molecules_getintranrg(i,'inst')
!!$      coul=zero
!!$      noncoul=zero
!!$      Do j=i,nsorbs
!!$      If (hybridmoves_checkDummy(params,j)) Cycle
!!$        !** coul and noncoul interactions.
!!$        pe=pe+molecules_getnoncoul(i,j,'inst')
!!$        noncoul=noncoul+molecules_getnoncoul(i,j,'inst')
!!$        pe=pe+molecules_getcoul(i,j,'inst')
!!$        coul=coul+molecules_getcoul(i,j,'inst')
!!$
!!$      End Do
!!$
!!$    End Do
!!$
!!$    PotNrg=pe+intra
!!$
!!$  End Subroutine hybridmoves_getPotNrg
!!$
!!$
!!$  !----------------------------------------------------------------------------
!!$  ! Displays the enrgies during the hybridmoves, displays all intramolecular
!!$  ! components separately, adopted from Marty's MD-display routine
!!$  !----------------------------------------------------------------------------
!!$  Subroutine hybridmoves_displayNrgs(params,species, unit)
!!$    Type(Hybrid_Move_Params), Intent(In) :: params
!!$    Type(AtMolCoords), Dimension(:), Intent(In) :: species
!!$
!!$    Integer :: i,j,nsorbs,moles,nmoles,nmoles_i,nmoles_j
!!$    Real(kind=RDbl) :: pe, ke, peavg, keavg, scaletime
!!$    Real(kind=RDbl) :: pei, kei, peavgi, keavgi
!!$    Character(len=4) :: units
!!$    Integer, Intent(In), Optional :: unit
!!$    Integer :: dUnit  ! display unit number
!!$    Logical :: xtraInfo
!!$
!!$    If (.Not.Present(unit)) Then
!!$      dUnit = 6
!!$    Else
!!$      dUnit = unit
!!$    End If
!!$
!!$    pe = 0.0_Rdbl
!!$    ke = 0.0_Rdbl
!!$    peavg = 0.0_Rdbl
!!$    keavg = 0.0_Rdbl
!!$    moles = 0
!!$
!!$    nsorbs = molecules_getnsorbs()
!!$
!!$    Write(dUnit,'(a)') " "
!!$
!!$    Write(dUnit,'(1x,2a)') "----------------------------------------------", &
!!$        "----------"
!!$    Write(dUnit,*) "Hybrid Display"
!!$    Do i = 1,nsorbs
!!$      nmoles = config_getnmoles(species,i)
!!$      If (nmoles == 0) Cycle
!!$      If (config_isfixed(species(i))) Cycle
!!$      Write(dUnit,'(1x,a,a)') Trim(molecules_name(i)), &
!!$          " Information"
!!$      Write(dUnit,'(3x,a12,4a15)') "Variable","Current","CumulAvg", &
!!$          "Block","Std"
!!$      Write(dUnit,'(3x,a15,4f15.6)') "Intra Energy : ", &
!!$          molecules_getintranrg(i,'inst')/nmoles, &
!!$          molecules_getintranrg(i,'cavg')/nmoles, &
!!$          molecules_getintranrg(i,'block')/nmoles, &
!!$          molecules_getintranrg(i,'std')/nmoles
!!$      Write(dUnit,'(3x,a15,4f15.6)') "Total Kin Nrg: ", &
!!$          molecules_getke(i,'inst')/nmoles, &
!!$          molecules_getke(i,'cavg')/nmoles, &
!!$          molecules_getke(i,'block')/nmoles, &
!!$          molecules_getke(i,'std')/nmoles 
!!$      If (intramolecular_hasint(i,'stretch')) Then
!!$        Write(dUnit,'(3x,a15,4f15.6)') "Stretch Energy : ", &
!!$            molecules_getintranrg(i,'inst','stretch')/nmoles, &
!!$            molecules_getintranrg(i,'cavg','stretch')/nmoles, &
!!$            molecules_getintranrg(i,'block','stretch')/nmoles, &
!!$            molecules_getintranrg(i,'std','stretch')/nmoles
!!$      End If
!!$      If (intramolecular_hasint(i,'constraint')) Then
!!$        Write(dUnit,'(3x,a15,4f15.6)') "Cnstr Energy  : ", &
!!$            molecules_getintranrg(i,'inst','constraint')/nmoles, &
!!$            molecules_getintranrg(i,'cavg','constraint')/nmoles, &
!!$            molecules_getintranrg(i,'block','constraint')/nmoles, &
!!$            molecules_getintranrg(i,'std','constraint')/nmoles
!!$      End If
!!$      If (intramolecular_hasint(i,'bending')) Then
!!$        Write(dUnit,'(3x,a15,4f15.6)') "Bending Energy: ", &
!!$            molecules_getintranrg(i,'inst','bending')/nmoles, &
!!$            molecules_getintranrg(i,'cavg','bending')/nmoles, &
!!$            molecules_getintranrg(i,'block','bending')/nmoles, &
!!$            molecules_getintranrg(i,'std','bending')/nmoles
!!$      End If
!!$      If (intramolecular_hasint(i,'torsion')) Then
!!$        Write(dUnit,'(3x,a15,4f15.6)') "Tosion Energy : ", &
!!$            molecules_getintranrg(i,'inst','torsion')/nmoles, &
!!$            molecules_getintranrg(i,'cavg','torsion')/nmoles, &
!!$            molecules_getintranrg(i,'block','torsion')/nmoles, &
!!$            molecules_getintranrg(i,'std','torsion')/nmoles
!!$      End If
!!$      If (intramolecular_hasint(i,'intrapair')) Then
!!$        Write(dUnit,'(3x,a15,4f15.6)') "1-4 LJ Energy : ", &
!!$            molecules_getintranrg(i,'inst','intrapair')/nmoles, &
!!$            molecules_getintranrg(i,'cavg','intrapair')/nmoles, &
!!$            molecules_getintranrg(i,'block','intrapair')/nmoles, &
!!$            molecules_getintranrg(i,'std','intrapair')/nmoles
!!$      End If
!!$      pe = molecules_getintranrg(i,'inst')/nmoles + pe
!!$      peavg = molecules_getintranrg(i,'cavg')/nmoles + peavg
!!$      ke = molecules_getke(i,'inst')/nmoles + ke
!!$      keavg = molecules_getke(i,'cavg')/nmoles + keavg
!!$      moles = moles + nmoles
!!$    End Do
!!$
!!$    !** Interaction energies
!!$    Do i = 1,nsorbs
!!$      Do j = i, nsorbs
!!$       nmoles_i = config_getnmoles(species,i)
!!$       if (config_isfixed(species(i))) nmoles_i = 0
!!$
!!$       nmoles_j = config_getnmoles(species,j)
!!$       if (config_isfixed(species(j))) nmoles_j = 0
!!$
!!$       !** Averaging will be done by the total number of moles in the i-j pair
!!$       nmoles=nmoles_i
!!$       If (i /= j) nmoles = nmoles + nmoles_j 
!!$
!!$        If (nmoles == 0) Cycle
!!$        Write(dUnit,'(1x,5a,i5)') "NonCoulEnergy for ",&
!!$            Trim(molecules_name(i)),"-", Trim(molecules_name(j)) &
!!$            ,"Total moles : ",nmoles
!!$        Write(dUnit,'(3x,a15,4f15.6)') "Inter Energy : ", &
!!$            molecules_getnoncoul(i,j,'inst')/nmoles, &
!!$            molecules_getnoncoul(i,j,'cavg')/nmoles, &
!!$            molecules_getnoncoul(i,j,'block')/nmoles, &
!!$            molecules_getnoncoul(i,j,'std')/nmoles
!!$        pe = molecules_getnoncoul(i,j,'inst')/nmoles + pe
!!$        peavg = molecules_getnoncoul(i,j,'cavg')/nmoles + peavg
!!$      End Do
!!$    End Do
!!$
!!$    !** Get any extra energies from the integration routines
!!$    !** e.g., Nose-Hoover energy, or others
!!$    pei = 0.0_Rdbl
!!$    kei = 0.0_RDbl
!!$    peavgi = 0.0_RDbl
!!$    keavgi = 0.0_Rdbl
!!$
!!$    If (Associated(params%mvparams%IntegParams%vverlet)) Then
!!$      Call velverlet_extraEnergy( &
!!$          params%mvparams%IntegParams%vverlet,pei,peavgi,'pe')
!!$      Call velverlet_extraEnergy( &
!!$          params%mvparams%IntegParams%vverlet,kei,keavgi,'ke')
!!$    Else If (Associated(params%mvparams%IntegParams%gear6)) Then
!!$      Call gear6_extraEnergy( &
!!$          params%mvparams%IntegParams%gear6,pei,peavgi,'pe')
!!$      Call gear6_extraEnergy( &
!!$          params%mvparams%IntegParams%gear6,kei,keavgi,'ke')
!!$    End If
!!$
!!$    If ((pei /= 0.0).And.(kei /= 0.0)) Then
!!$      Write(dUnit,'(1x,a)') "Extra Integrator Energies (NH, etc) for &
!!$          & first move_type only"
!!$      Write(dUnit,'(3x,a15,2f15.6)') "Potential Nrg: ",pei/moles,peavgi/moles
!!$      Write(dUnit,'(3x,a15,2f15.6)') "Kinetic Energ: ",kei/moles,keavgi/moles
!!$    End If
!!$
!!$    if (moles>0) then
!!$    pe = pe + pei/moles
!!$    ke = ke + kei/moles
!!$    peavg = peavg + peavgi/moles
!!$    keavg = keavg + keavgi/moles
!!$    endif
!!$    Write(dUnit,'(1x,a)') "System Energy Totals : current, Cumulative "
!!$    If(moles==0) Then
!!$      Write(dUnit,*) "No Molecules at this moment"
!!$    Else
!!$      Write(dUnit,'(3x,a15,2f15.6)') "Total Pot Nrg: ",pe,peavg
!!$      Write(dUnit,'(3x,a15,3f15.6)') "Total Kin Nrg: ",ke,keavg
!!$      Write(dUnit,'(3x,a15,2f15.6)') "Total Energy : ",pe+ke,peavg+keavg
!!$    Endif
!!$    Write(dUnit,'(1x,a)') &
!!$        "--------------------------------------------------------"
!!$
!!$  End Subroutine hybridmoves_displayNrgs
