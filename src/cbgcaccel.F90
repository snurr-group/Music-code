!-----------------------------------------------------------------------------
! This module contains the data structure and routines for accelerating 
! Configurational biased Grand Canonical Monte Carlo (GCMC) simulations.
!-----------------------------------------------------------------------------

Module cbgcaccel
  Use branchedcoords, Only : branchedcoords_toxyz, branchedcoords_assignNode, &
      branchedcoords_dumpNode
  Use config, Only: AtMolCoords, config_getnmoles, config_getnatoms, &
      config_setnmoles, config_checkandincr, config_delmolec
  Use cbgcmoves, Only : CBGCMC_Move_Params, cbgcmoves_move, &
      cbgcmoves_updatenmolecs, cbgcmoves_nvtmove
  Use defaults, Only: RDbl, strLen, dashedline,dashedline2, twopi, d_res_file &
      ,d_con_file, hplanck, MAX_SORBS, NO_OF_INTRA_POTS, STRETCH_INDEX,&
      BENDING_INDEX,TORSION_INDEX,CONSTRAINT_INDEX,TOTAL_INDEX,&
      INTRAPAIR_INDEX,Rgas, Nav, scalepe, kcalmole_kb, one,zero, &
      default_MAX_EXP, MAX_ATOMS
  Use file, Only: file_getunit,file_gettype,file_open
  Use general, Only: genparams
  Use interact, Only: Molecule_System_Energies, interact_msyscopy, &
      interact_msystotnrg, interact_initmsys, interact_msysupdate
  Use molecule, Only : MolecularParams
  Use molecules, Only: molecules_getmass, molecules_getcoul &
      ,molecules_getnoncoul, molecules_getintranrg, molecules_getnsorbs, &
      molecules_updateEnergySS, molecules_getpointer, molecules_updateenergy, &
      molecules_name, molecules_updateIntraEnergy, molecules_displaynrg, molecules_getnatoms
  Use moves, Only : moves_insert, moves_delete
  Use mcmoveset, Only: Move_Set, mcmoveset_getmoveno
  Use random,Only :rranf
  Use stats, Only : Statistics, stats_init,stats_update,  stats_getvalue,stats_display 
  Use simcell, Only: SimCell_Params , simcell_pbc
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, int2str, readblank
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag

  Implicit None
  Save

  Private
  Public :: cbgcaccel_init, CBGC_Accelerator, sorb_accelerator, &
      Delete_Storage, Insert_Storage, cbgcaccel_dosim
 
  Type CBGC_Accelerator
    Integer                             :: nsorbs
    Type(sorb_accelerator),Dimension(:),Pointer :: sorbs
    Integer                :: max_attempts    
    Real(kind=RDbl)        :: gc_factor
    Type(Delete_Storage),Dimension(:), Pointer   :: del_store
    Type(Insert_storage), Dimension(:), Pointer   :: ins_store
    Type(Molecule_System_Energies) :: newnrgs, oldnrgs
  End Type CBGC_Accelerator
  
  Type sorb_accelerator
    Integer                :: N_insert, max_Trials, nvt_simnum
    Real(kind=RDbl)        :: nvtratio, insdelratio,inswt,delwt
    Integer                :: currenttrials,current_success
  End Type sorb_accelerator

  Type Delete_Storage
    Real(kind=RDbl),Dimension(:), Pointer :: bias,expfac
    Type(Molecule_System_Energies), Dimension(:), Pointer :: nrgs
  End Type Delete_Storage

  Type Insert_Storage
    Integer, Dimension(:) , Pointer       :: index
    Real(kind=RDbl),Dimension(:), Pointer :: bias,expfac
    Real(kind=RDbl),Dimension(:,:,:), Pointer :: posns
    Type(Molecule_System_Energies), Dimension(:), Pointer :: nrgs
  End Type Insert_Storage

 !SDEBUG ! maximum allowed molecules per sorb
  Integer :: Max_Molecs =100
  Real(kind=RDbl) , Dimension(MAX_ATOMS) :: angle1, angle2
  Integer,Dimension(MAX_ATOMS) :: num
  Real(kind=RDbl) BOLTZ_DEL_MAX, BOLTZ_INS_MAX
  Integer :: mul_fac=1000
 !SDEBUG
  Contains

  !------------------------------------------------
  ! Initializes the accelerator
  !------------------------------------------------
  Subroutine cbgcaccel_init(params, chainsorbs, sorbates, unitno)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(CBGCMC_Move_Params), Dimension(:), Intent(inout) :: chainsorbs
    Type(CBGC_Accelerator), Intent(inout)  :: params
    Integer, Intent(in)                    :: unitno

    Integer :: natoms, nsorbs,error,i, j, N_insert, max_Trials, sorbtype

    Write(*,*) " Initializing the cbgcmc acceleration session"
    !** number of cbgcmc sorbs
    Read(unitno,*) nsorbs
    params%nsorbs=nsorbs

    Allocate(params%sorbs(nsorbs),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-Accelerator")

    Allocate(params%del_store(nsorbs),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-Accelerator")

    Allocate(params%ins_store(nsorbs),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-Accelerator")

    Read(unitno,*) params%max_attempts
    
    ! temperory arrays for passing nrgs
    Call interact_initmsys(params%newnrgs)
    Call interact_initmsys(params%oldnrgs)

    Do i=1, nsorbs
      Read(unitno,*) params%sorbs(i)%nvt_simnum
      Read(unitno,*) N_insert, max_Trials
      params%sorbs(i)%N_insert=N_insert
      params%sorbs(i)%max_Trials=max_Trials
      Read(unitno,*) params%sorbs(i)%nvtratio
      Read(unitno,*) params%sorbs(i)%insdelratio
      params%sorbs(i)%inswt= params%sorbs(i)%insdelratio/ &
          (1+params%sorbs(i)%insdelratio)     
      params%sorbs(i)%delwt= 1-params%sorbs(i)%inswt

      !** The insertion storage
      sorbtype=chainsorbs(i)%sorbtype
      natoms=config_getnatoms(sorbates,sorbtype)
      Allocate(params%ins_store(i)%index(max_Trials),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-accelerator")
      Allocate(params%ins_store(i)%posns(N_insert,natoms,3),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-accelerator")
      Allocate(params%ins_store(i)%bias(N_insert),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-accelerator")
      Allocate(params%ins_store(i)%expfac(N_insert),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-accelerator")
      Allocate(params%ins_store(i)%nrgs(N_insert),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-accelerator")
      Do j=1,N_insert
        Call interact_initmsys(params%ins_store(i)%nrgs(j))
      End Do

      !** The deletion storage
      Allocate(params%del_store(i)%bias(Max_Molecs),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-accelerator")
      Allocate(params%del_store(i)%expfac(Max_Molecs),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-accelerator")
      Allocate(params%del_store(i)%nrgs(Max_Molecs),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbgc-accelerator")

      Do j=1,Max_Molecs
        Call interact_initmsys(params%del_store(i)%nrgs(j))
      End Do

    End Do

  End Subroutine cbgcaccel_init

  !------------------------------------------------
  ! Attempts - deletion of all molecule
  !          - insertion of many configurations
  !  and then stores these values, Accelerator will later use these 
  !  results repeatedly for insertions and deletions
  !------------------------------------------------
  Subroutine cbgcaccel_prepare(params, sorbates, chainsorbs, cbsorb, simno)
    Type(CBGC_Accelerator), Intent(inout)  :: params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(CBGCMC_Move_Params), Dimension(:), Intent(inout) :: chainsorbs
    Integer, Intent(in)                    :: cbsorb, simno

    Integer :: j,sorbtype, moveno, molec, nmoles, totalconfs
    Real(kind=RDbl) :: biasfactor,ufinal,min_nrg
    Logical :: nrg_calc_flag

    BOLTZ_DEL_MAX=zero
    BOLTZ_INS_MAX=zero

    !** Attempt deletion of all molecules
    sorbtype=chainsorbs(cbsorb)%sorbtype
    nmoles=config_getnmoles(sorbates,sorbtype)
    moveno=mcmoveset_getmoveno(chainsorbs(cbsorb)%moveset,"DELETE")
      Do molec=1,nmoles
        Call moves_delete(sorbates, sorbtype, molec, &
            chainsorbs(cbsorb)%mcmoves(moveno)%mvparams,biasfactor, &
            params%oldnrgs)
        Call cbgcaccel_copydelete(params, sorbates, sorbtype, molec, &
            biasfactor, params%oldnrgs, cbsorb, chainsorbs(cbsorb), simno)
      End Do

    !** Attemp insertion of molecules
      sorbtype=chainsorbs(cbsorb)%sorbtype
      nmoles=config_getnmoles(sorbates, sorbtype)
      molec= nmoles+1
      !** Check if we need to increment the size of the various config arrays
      Call config_checkandincr(sorbates, sorbtype, nmoles)

      Call config_setnmoles(sorbates, sorbtype, molec)

      moveno=mcmoveset_getmoveno(chainsorbs(cbsorb)%moveset,"INSERT")
      min_nrg=zero
      totalconfs=0
      !** default trials
      params%sorbs(cbsorb)%currenttrials=params%sorbs(cbsorb)%max_Trials
      params%sorbs(cbsorb)%current_success=0
      Do j=1,params%sorbs(cbsorb)%max_Trials
        Call moves_insert(sorbates, sorbtype, molec, &
            chainsorbs(cbsorb)%mcmoves(moveno)%mvparams, biasfactor, &
            nrg_calc_flag, params%newnrgs)
        If (nrg_calc_flag) Then
          ufinal=interact_msystotnrg(params%newnrgs)
          If (ufinal>zero) nrg_calc_flag=.False.
          If (ufinal<min_nrg) min_nrg=ufinal
        Endif

        If (nrg_calc_flag) Then
          totalconfs=totalconfs+1
          Call cbgcaccel_copyinsert(params, sorbates, sorbtype, molec,&
              biasfactor, params%newnrgs, cbsorb, chainsorbs(cbsorb), &
              simno ,totalconfs)
          params%ins_store(cbsorb)%index(j)=totalconfs
          If (totalconfs==params%sorbs(cbsorb)%N_insert) Then
            params%sorbs(cbsorb)%currenttrials=j
            Exit
          Endif
        Else
          params%ins_store(cbsorb)%index(j)=0
        Endif
      End Do
      params%sorbs(cbsorb)%current_success=totalconfs
    !** setting actual moles back to correct
    Call config_setnmoles(sorbates, sorbtype, molec-1)

!!$    Write(*,*) "max, del, ins",BOLTZ_DEL_MAX, BOLTZ_INS_MAX

  End Subroutine cbgcaccel_prepare



  !---------------------------------------------------------------------
  ! This uses the stored values for attempting insertions and deletions
  !--------------------------------------------------------------------
  Subroutine cbgcaccel_dosim(params, sorbates, scell, chainsorbs, cbgcsorb, &
      simno)
    Type(CBGC_Accelerator), Intent(inout)  :: params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(In) :: scell
    Type(CBGCMC_Move_Params), Dimension(:), Intent(inout) :: chainsorbs
    Integer, Intent(in)                    :: cbgcsorb, simno

    Integer :: total_trials, max_attempts, max_allowed, trials, k, sorbno
    Logical :: do_nvt, success
    total_trials=0

    !** This is in millions
    max_attempts=params%max_attempts
    Do
!!$ exit     Write(*,*) "Preparing for simulation , sorb :", cbgcsorb
      !** prepare storage lists
      Call cbgcaccel_prepare(params, sorbates, chainsorbs, cbgcsorb, simno)    
      
      max_allowed=max_attempts - total_trials

!!$       Write(*,*) "Preparing for attempts, max=", max_allowed
!!$      !** making attempts using the above storage
      Call cbgcaccel_makeattempts( params, sorbates, scell, chainsorbs, &
          cbgcsorb, simno, trials, max_allowed, do_nvt)
    

      If (do_nvt) Then
!!$        Write(*,*) "doing NVT"
        !** doing for all the sorbs
        Do k=1,params%sorbs(cbgcsorb)%nvt_simnum
          sorbno=Int(rranf()*Size(chainsorbs,1))+1
        !** Pick a move type randomly and execute it
        success =  cbgcmoves_nvtmove(chainsorbs(sorbno), &
            sorbates, scell,simno)
        End do
      Endif
!!$      Write(*,*) "Trials conducted ", trials
!!$      total_trials=total_trials + trials
!!$
!!$      If (total_trials >= max_attempts ) Exit
!!$
!!$      !** keep looping till we finisg required number of attempts
      Exit

    End Do

    
  End Subroutine cbgcaccel_dosim



  !---------------------------------------------------------------------
  ! This uses the stored values for attempting insertions and deletions
  !--------------------------------------------------------------------
  Subroutine cbgcaccel_makeattempts(params, sorbates, scell, chainsorbs, &
      cbgcsorb, simno, trials, maxtrials, do_nvt)
    Type(CBGC_Accelerator), Intent(inout)  :: params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(In) :: scell
    Type(CBGCMC_Move_Params), Dimension(:), Intent(inout) :: chainsorbs
    Integer, Intent(in)                    :: cbgcsorb, simno
    Integer , Intent(out) :: trials
    Integer , Intent(in) ::  maxtrials
    Logical, Intent(out) :: do_nvt
    Integer :: i, mil, sorbtype,nconfs,confnum,index, nmoles, molec, del_attempt,ins_attempt, del_num, moveno, tottrials
    Logical :: exitflag,no_ins
    Real(kind=RDbl) :: inswt
    do_nvt=.False.
    trials=0
    del_attempt=0
    ins_attempt=0

    exitflag=.False.

    sorbtype=chainsorbs(cbgcsorb)%sorbtype
    nconfs=params%sorbs(cbgcsorb)%current_success
    tottrials=params%sorbs(cbgcsorb)%currenttrials

    inswt=params%sorbs(cbgcsorb)%inswt
    no_ins=.False.
    If (nconfs<1) no_ins=.True.
    nmoles=config_getnmoles(sorbates,sorbtype)


    If (.Not.no_ins) Then
      !** we have to try both insertion aand deletions

      !** mul_fac=1000, is used to avoid excessive rranf usage
      Do i=1,maxtrials*mul_fac
        !** actual nvt ratio= nvtratio/mulfac
        If (rranf()<params%sorbs(cbgcsorb)%nvtratio) Then
          exitflag=.True.
          do_nvt=.True.
        Endif

        !** Note that maxtrials was read in millions
        !** Our aim is to do maxtrials*million trials
        Do mil=1,(1000000)/mul_fac
          
          !** Ins or Del
          If (rranf()< inswt) Then

            confnum=Int(rranf()*tottrials)+1
            index=params%ins_store(cbgcsorb)%index(confnum)

            !** Only use good configs
            If (index>0) Then 
              If (rranf()<params%ins_store(cbgcsorb)%expfac(index)) Then
                Call cbgcaccel_acceptInsert(sorbates, scell, sorbtype, &
                    cbgcsorb, params, index )
                moveno=mcmoveset_getmoveno(chainsorbs(cbgcsorb)%moveset, &
                    "INSERT")
                chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%succ=&
                    chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%succ+1
                chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%att=&
                    chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%att+1

                
                exitflag=.True.
                Exit
              Endif
            Endif

          Else
            !** select molec
            molec=Int(rranf()*nmoles)+1
            If (rranf()<params%del_store(cbgcsorb)%expfac(molec)) Then
              Call cbgcaccel_acceptDelete(sorbates, sorbtype, cbgcsorb, params, &
                  molec)
                moveno=mcmoveset_getmoveno(chainsorbs(cbgcsorb)%moveset, &
                    "DELETE")
                chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%succ=&
                    chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%succ+1
                chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%att=&
                    chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%att+1
              exitflag=.True.
              Exit
            Endif
          Endif
        End Do
        If (exitflag) exit
      End Do
      Write(*,*) "Molec", sorbtype, "Percentage Finished ",i/&
          (one * maxtrials*mul_fac )
    Else

      !** only deletions have to be tried, so we assume that we are trying 
      ! all those insertions but since the energies are high , they are
      ! anyway getting deleted
      ! del_num is in In millions
      del_num=Int(1000000/params%sorbs(cbgcsorb)%insdelratio+0.5)
      If (del_num/mul_fac<2) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Else
!!$        Write(*,*) "del_num",del_num
      Endif
      Do i=1,maxtrials*mul_fac
        If (rranf()<params%sorbs(cbgcsorb)%nvtratio) Then
          exitflag=.True.
          do_nvt=.True.
        Endif
        Do mil=1,del_num/mul_fac
          !** select molec
          molec=Int(rranf()*nmoles)+1
          If (rranf()<params%del_store(cbgcsorb)%expfac(molec)) Then
            Call cbgcaccel_acceptDelete(sorbates, sorbtype, cbgcsorb, params, &
                molec)
                moveno=mcmoveset_getmoveno(chainsorbs(cbgcsorb)%moveset, &
                    "DELETE")
                chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%succ=&
                    chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%succ+1
                chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%att=&
                    chainsorbs(cbgcsorb)%mcmoves(moveno)%stats%att+1
            exitflag=.True.
            Exit
          Endif
        End Do
        If (exitflag) Exit
      End Do
      Write(*,*) "Molec", sorbtype, "Percentage Finished ",i/&
          (one *maxtrials*mul_fac) 

!!$      Write(*,*) "Finished ",i
    Endif

  End Subroutine cbgcaccel_makeattempts



  !---------------------------------------------------------------------
  ! Does mundane things to be done with an acceptance
  ! Requires : -sorbates : original positions structure
  !            -params   : accel params (includes storage
  !            -index    : index in the storage aray of the new conf
  !---------------------------------------------------------------------
  Subroutine cbgcaccel_acceptInsert(sorbates, scell, sorbtype, cbgcsorb, &
      params, index)
    Type(CBGC_Accelerator), Intent(in)  :: params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(In) :: scell
    Integer, Intent(in)                    :: index, sorbtype, cbgcsorb


    Integer :: nmoles, molec, natoms, moveno
    ! params%oldnrgs will not be used anyway
    ! This update only takes care of inst energy, all averages will be skewed
    Call interact_msysupdate(sorbtype, &
        params%oldnrgs, params%ins_store(cbgcsorb)%nrgs(index), 0 )

    !** Now copy the actual coordinates
    nmoles=config_getnmoles(sorbates,sorbtype)
    molec=nmoles+1
    Write(*,*) "Accepting insertion at molec", molec


    Call config_setnmoles(sorbates,sorbtype,molec)
    natoms=molecules_getnatoms(sorbtype)
    !** first atom x,y,z
    sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos%comp(1)= &
        params%ins_store(cbgcsorb)%posns(index,1,1)
    sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos%comp(2)= &
        params%ins_store(cbgcsorb)%posns(index,1,2)
    sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos%comp(3)=&
        params%ins_store(cbgcsorb)%posns(index,1,3)
    angle1(1:natoms-1)=params%ins_store(cbgcsorb)%posns(index,2:natoms,1)
    angle2(1:natoms-1)=params%ins_store(cbgcsorb)%posns(index,2:natoms,2)
    Call branchedcoords_AssignNode(sorbates(sorbtype)%gcoords(molec)%branchedcoords%root, angle1, angle2, num)

    Call branchedcoords_toxyz( &
        sorbates(sorbtype)%gcoords(molec)%branchedcoords, & 
        sorbates(sorbtype)%coords(1:natoms,molec)%rp) 
    Call simcell_pbc(scell, sorbates(sorbtype)%coords(1:natoms,molec)%rp, &
        sorbates(sorbtype)%coords(1:natoms, molec)%r, &
        sorbates(sorbtype)%coords(1:natoms, molec)%cr)
    
  End Subroutine cbgcaccel_acceptInsert



  !---------------------------------------------------------------------
  ! Does mundane things to be done with an acceptance
  ! Requires : -sorbates : original positions structure
  !            -params   : accel params (includes storage
  !            -molec    : molec to be deleted 
  !---------------------------------------------------------------------
  Subroutine cbgcaccel_acceptDelete(sorbates, sorbtype, cbgcsorb, &
      params, molec)
    Type(CBGC_Accelerator), Intent(in)  :: params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Integer, Intent(in)                    :: molec, sorbtype, cbgcsorb

    Write(*,*) "Accepting Deletion at molec", molec

    ! params%newnrgs will not be used anyway
    ! This update only takes care of inst energy, all averages will be skewed
    Call interact_msysupdate(sorbtype, &
        params%del_store(cbgcsorb)%nrgs(molec), params%newnrgs, 2)

    !** remove a molecule and adjust the array; also sets nmoles=nmoles-1
    Call config_delmolec(sorbates, sorbtype, molec)
    
  End Subroutine cbgcaccel_acceptDelete




  !-------------------------------------------------------------------------
  ! Copies the bias, energies,  expfac of trial deletions into
  ! the storage arrays
  !-------------------------------------------------------------------------
  Subroutine cbgcaccel_copydelete(params, sorbates, sorbtype, molec, &
            biasfactor, oldnrgs, cbgcsorbno, chainsorb , simno)
    Type(CBGC_Accelerator), Intent(inout)            :: params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(CBGCMC_Move_Params), Intent(inout)          :: chainsorb
    Integer, Intent(in)                        :: sorbtype, molec, simno
    Integer, Intent(in)                        :: cbgcsorbno
    Type(Molecule_System_Energies),Intent(in)  :: oldnrgs
    Real(kind=RDbl), Intent(in)                :: biasfactor

    Real(kind=RDbl) :: uinit, rti, B, utemp, boltz, factor
    Integer         :: nmoles,i

    uinit = interact_msystotnrg(oldnrgs)
    rti=chainsorb%rti
    B= chainsorb%fuglist(simno)%B
    utemp=rti*uinit -B

    ! store biasfactor
    params%del_store(cbgcsorbno)%bias(molec)   = biasfactor

    factor=biasfactor*params%sorbs(cbgcsorbno)%insdelratio
    utemp = utemp + log(factor)

    !** Avoid overflow/underflow errors during exponentiation
    If (utemp > default_MAX_EXP) utemp = default_MAX_EXP
    If (utemp < -default_MAX_EXP) utemp = -default_MAX_EXP
    nmoles=config_getnmoles(sorbates, sorbtype)

    boltz = Exp(utemp)*Real(nmoles)
    params%del_store(cbgcsorbno)%expfac(molec) = boltz 
    Call interact_msyscopy(oldnrgs,params%del_store(cbgcsorbno)%nrgs(molec))
    If (boltz>BOLTZ_DEL_MAX) BOLTZ_DEL_MAX=boltz
    
  End Subroutine cbgcaccel_copydelete


  !-------------------------------------------------------------------------
  ! Copies the coords, bias, energies,  expfac of trial insertions into
  ! the storage arrays
  !-------------------------------------------------------------------------
  Subroutine cbgcaccel_copyinsert(params, sorbates, sorbtype, molec, &
            biasfactor, newnrgs, cbgcsorbno, chainsorb , simno, confnum)
    Type(CBGC_Accelerator), Intent(inout)            :: params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(CBGCMC_Move_Params), Intent(inout)          :: chainsorb
    Integer, Intent(in)                        :: sorbtype, molec, simno    
    Integer, Intent(in)                        :: cbgcsorbno
    Type(Molecule_System_Energies),Intent(in)  :: newnrgs
    Real(kind=RDbl), Intent(in)                :: biasfactor
    Integer, Intent(in)                        :: confnum

    Real(kind=RDbl) :: ufinal, rti, B, utemp, boltz, factor
    Integer         :: natoms, nmoles, i, atomnum

    ufinal = interact_msystotnrg(newnrgs)
    rti=chainsorb%rti
    B= chainsorb%fuglist(simno)%B

    utemp = -rti * ufinal + B

    ! store biasfactor
    params%ins_store(cbgcsorbno)%bias(confnum)   = biasfactor

    factor = biasfactor/params%sorbs(cbgcsorbno)%insdelratio

    utemp = utemp + log(biasfactor)

    !** Avoid overflow/underflow errors during exponentiation
    If (utemp > default_MAX_EXP) utemp = default_MAX_EXP
    If (utemp < -default_MAX_EXP) utemp = -default_MAX_EXP
    nmoles=config_getnmoles(sorbates, sorbtype)

    boltz = Exp(utemp)/Real(molec)
    If (boltz>BOLTZ_INS_MAX) BOLTZ_INS_MAX=boltz
    params%ins_store(cbgcsorbno)%expfac(confnum) = boltz
 
    Call interact_msyscopy(newnrgs,params%ins_store(cbgcsorbno)%nrgs(confnum))

    natoms=config_getnatoms(sorbates, sorbtype)

    !** first atom x,y,z
    params%ins_store(cbgcsorbno)%posns(confnum,1,1)=sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos%comp(1)
    params%ins_store(cbgcsorbno)%posns(confnum,1,2)=sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos%comp(2)
    params%ins_store(cbgcsorbno)%posns(confnum,1,3)=sorbates(sorbtype)%gcoords(molec)%branchedcoords%startpos%comp(3)

    Call branchedcoords_DumpNode(sorbates(sorbtype)%gcoords(molec)%branchedcoords%root, angle1, angle2, num)

    params%ins_store(cbgcsorbno)%posns(confnum,2,1)=sorbates(sorbtype)%gcoords(molec)%branchedcoords%root%baPc
    params%ins_store(cbgcsorbno)%posns(confnum,2,2)=sorbates(sorbtype)%gcoords(molec)%branchedcoords%root%taC


    Do atomnum=3,natoms
    params%ins_store(cbgcsorbno)%posns(confnum,atomnum,1)=angle1(atomnum-1)
    params%ins_store(cbgcsorbno)%posns(confnum,atomnum,2)=angle2(atomnum-1)
    End Do

  End Subroutine cbgcaccel_copyinsert


End Module cbgcaccel
