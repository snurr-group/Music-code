!-------------------------------------------------------------------------------
! This module also contains a data type (Subset_Interactions) designed to 
! efficiently encapsulate sets of interactions and their storage.  
!
! For many Monte Carlo methods, only part of the system is perturbed.
! In such a case, it is wasteful to reevaluate all the system
! interactions.  This module can store temporary interactions between
! two subsets of the system while a move is being considered.  If the
! move is accepted, this partial set of interactions can be then copied
! into the full interaction storage.  In the case where the storage
! isn't detailed enough to be updated simply by copying replacement
! values in, two temporary storage structures are necessary.  One
! contains the current configuration interactions for that subset and
! one contains the new interactions.  If the move is accepted, the current
! temporary storage is subtracted from the full and the new temporary
! storage is added to the full.  One such example is when single molecules
! are moved, but only SPC-detail storage is used.
!
! ALL MOVE ROUTINES ROUTE THEIR INTERACTION CALCULATIONS THROUGH THIS MODULE
!
! Important Routines:
!   subinteract_init -- initializes an subset interaction data type
!   subinteract_update -- updates the full interactions from the subsets
!   subinteract_updateall -- updates the full interactions and stored statistics
!   subinteract_int -- evaluates interactions between subsets of the system
!-------------------------------------------------------------------------------

Module subinteract

  Use defaults, Only: RDbl,strLen,lstrLen,xlstrLen,zero,one,dbgflag, &
      kcalmole_kb
  Use utils, Only: toupper,getdepth,checkandstop,allocErrDisplay,findint
  Use molecules, Only: molecules_getnatoms,molecules_gettype
  Use config, Only: AtMolCoords, config_getnmoles, config_isfixed
  Use simcell, Only: SimCell_Params, simcell_getmolectype
  Use storebase, Only: EnergyPlus, storebase_chkintra, storebase_copy, &
      storebase_init, storebase_clean, storebase_subtract, storebase_display
  Use storetop, Only: Forcefield_Results, storetop_update, storetop_copy, &
      storetop_setmap, storetop_chksizes, storetop_initcopy, storetop_fillsub, &
      storetop_extract, storetop_zero, storetop_display, storetop_clean, &
      storetop_subtract, storetop_add, storetop_delpossible, storetop_sum
  Use interact, Only: Interaction_Model, interact_init, interact_chkmolnums, &
      interact_int, interact_changenmoles, interact_updatestats

  Implicit None
  Save

  Private 
  Public :: Subset_Interactions, subinteract_init, &
      subinteract_int, subinteract_copy, subinteract_zerotemp, &
      subinteract_oldnrg, subinteract_newnrg, subinteract_changenmoles, &
      subinteract_update, subinteract_updateall, subinteract_chksubint, &
      subinteract_simpleupdate, subinteract_chksizes, subinteract_clean, &
      subinteract_display, subinteract_fullint, subinteract_chkmolnums

  !** Stores path to and interactions between a subset of the system
  !** and every other species of the system.  Temporary partial storage 
  !** structures are: temp_now (current configuration) and temp_try (trial)
  !** Update modes are:  0 -> no updates, straight use of storage
  !**                    1 -> replacement mode updates
  !**                    2 -> subtract/add mode updates
  Type Subset_Interactions
    Character(len=strLen)              :: id
    Type(Interaction_Model), Pointer   :: imodel
    Integer                            :: depth1,depth2
    Integer, Dimension(3)              :: subset1,subset2
    Logical                            :: calc_accels,skip_intra
    Integer                            :: update_mode
    Type(EnergyPlus)                   :: now_sum,try_sum
    Type(Forcefield_Results), Pointer  :: temp_now,temp_try
  End Type Subset_Interactions

Contains
  !---------------------------------------------------------------------------
  ! Initializes the subset interaction data structure.  If necessary, it will
  ! create temporary storage structures for subset1--subset2 interactions.
  ! Note that the interaction model must first be initialized.
  ! Requires:  subinfo -- subset interaction structure
  !            imodel -- interaction model information
  !            label -- string label for new structure
  !            sim_type -- identifies simulation type for updating decisions
  !            subset1 -- 1st subset
  !            subset2 -- 2nd subset
  !---------------------------------------------------------------------------
  Subroutine subinteract_init(subinfo,imodel,label,sim_type,subset1,subset2)
    Type(Subset_Interactions), Intent(Out)       :: subinfo
    Type(Interaction_Model), Intent(In), Target  :: imodel
    Character(*), Intent(In)                     :: label,sim_type
    Integer, Dimension(3), Intent(In)            :: subset1,subset2

    Integer           :: i,nspc,error,nderivs,detail

    !** Nullify the temperary storage structure pointer
    Nullify(subinfo%temp_try)
    Nullify(subinfo%temp_now)

    !** Store basic information about interaction model
    subinfo%id = label
    subinfo%imodel => imodel
    subinfo%update_mode = 1

    !** Get the depth for each subset and store the subsets
    subinfo%depth1 = getdepth(subset1)    
    subinfo%depth2 = getdepth(subset2)    
    subinfo%subset1 = subset1
    subinfo%subset2 = subset2

    !** Set defaults
    subinfo%calc_accels = .False.
    subinfo%skip_intra = .False.

    !** Modify defaults according to specific simulation type needs
    Select Case(ToUpper(Trim(sim_type)))
    Case ('MD')
      subinfo%update_mode = 0
      subinfo%calc_accels = .True.
      !** Point the temporary storage to the full storage
      subinfo%temp_try => subinfo%imodel%results(1)

    Case ('NEB')
      subinfo%update_mode = 0
      subinfo%calc_accels = .False.
      !** Point the temporary storage to the full storage
      subinfo%temp_try => subinfo%imodel%results(1)

    Case ('MC')
      Allocate(subinfo%temp_try,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Case ('HMC')
      subinfo%calc_accels = .True.
      Allocate(subinfo%temp_try,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Case ('NULL')
      subinfo%update_mode = 0
      subinfo%calc_accels = .False.
      !** Point the temporary storage to the full storage
      subinfo%temp_try => subinfo%imodel%results(1)

    Case Default
      Write(0,'(1x,2a,i4,2a)') __FILE__,' : ',__LINE__, &
          ' Unable to understand simulation type specification: ',Trim(sim_type)
      Stop
    End Select
    
    !** Allocate the temporary interaction storage if it will be used
    If (subinfo%update_mode > 0) Then

      !** Determine the type of updating needed if necessary.  The
      !** subtract/add mode will be used if the perturbed subset 
      !** is smaller than the full storage detail
      detail = imodel%results(1)%storelevel - 1
      If (detail < (Max(subinfo%depth1,subinfo%depth2))) subinfo%update_mode = 2

      If (subinfo%depth2 /= 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__,' : ',__LINE__, &
            ' Cannot initialize temporary storage with subset2 /= full system'
        Stop
      End If
    
      !** Initialize the two sums, even if not required
      nderivs = subinfo%imodel%results(1)%nderivs
      Call storebase_init(subinfo%now_sum,nderivs,.True.)
      Call storebase_init(subinfo%try_sum,nderivs,.True.)

      !** Initialize the temporary interactions storage
      Select Case(subinfo%depth1)
      Case (0)
        If (subinfo%depth2 == 0) Then
          subinfo%update_mode = 1
          Call storetop_initcopy(subinfo%temp_try,imodel%results(1),(/0,0,0/))
        Else
          Write(0,'(1x,2a,i4,a,i2)') __FILE__,' : ',__LINE__, &
              ' Cannot initialize temporary storage with depth1 = 0 /= depth2'
          Stop
        End If

      Case (2)
        If (subinfo%update_mode > 1) Then
          !** In this case we need a full copy of the storage in order
          !** to store mol-sys interactions in a way that they can be used
          !** for updating the 'full' storage.
          Call storetop_initcopy(subinfo%temp_try,imodel%results(1),(/0,0,0/))
        Else
          Call storetop_initcopy(subinfo%temp_try,imodel%results(1), &
              subinfo%subset1)
          Call storetop_fillsub(subinfo%temp_try,subinfo%subset1)
          Call storetop_setmap(subinfo%temp_try,.True.,subinfo%subset1)
        End If

      Case Default
        Write(0,'(1x,2a,i4,a,i2)') __FILE__,' : ',__LINE__, &
            ' Cannot initialize temporary storage with subset1 depth ', &
            subinfo%depth1
        Stop
      End Select

      !** Create the temp_now structure if needed for subtract/add update mode
      If (subinfo%update_mode > 1) Then
        Allocate(subinfo%temp_now,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
        Call storetop_initcopy(subinfo%temp_now,subinfo%temp_try,(/0,0,0/))
      End If
      
    End If

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'original structure:'
    Call storetop_display(imodel%results(1),.False.,2,6)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'temporary mol-spc structure:'
    Write(*,*) 'for subset1, subset2: ',subinfo%subset1,subinfo%subset2
    Call storetop_display(subinfo%temp_try,.False.,2,6)
    Stop
#endif

  End Subroutine subinteract_init

  !----------------------------------------------------------------------
  ! Copies the subset interaction structures
  ! Requires:  dest -- subset interactions structure (already init)
  !            orig -- original subset interactions structure
  !----------------------------------------------------------------------
  Subroutine subinteract_copy(dest,orig)
    Type(Subset_Interactions), Intent(Out)  :: dest
    Type(Subset_Interactions), Intent(In)   :: orig

    Logical          :: success

    dest%id = orig%id
    dest%imodel => orig%imodel
    dest%depth1 = orig%depth1
    dest%depth2 = orig%depth2
    dest%calc_accels = orig%calc_accels
    dest%update_mode = orig%update_mode
    dest%skip_intra = orig%skip_intra

    !** Copy or link interaction storage
    If (orig%update_mode > 0) Then
      Call storebase_copy(dest%now_sum,orig%now_sum)
      Call storebase_copy(dest%try_sum,orig%try_sum)
      If (.Not. storetop_copy(dest%temp_try,orig%temp_try)) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
            ' could not copy "try" subset interactions'
        Stop
      End If
      If (.Not. storetop_copy(dest%temp_now,orig%temp_now)) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
            ' could not copy "now" subset interactions'
        Stop
      End If
    Else
      dest%temp_try => orig%temp_try
      Nullify(dest%temp_now)
    End If
    
  End Subroutine subinteract_copy

  !---------------------------------------------------------------------------
  ! Calculates or retrieves interactions between subsets of the system and 
  ! stores the results for later use in the designated temporary storage 
  ! structure.
  ! Requires:  subinfo -- subset interaction structure
  !            species -- coordinate data structure 
  !            simcell -- simulation cell information
  !            fast -- logical flag indicating fast or slow evaluations
  !            recalc -- True => forces energy recalculation
  !            skip_intra -- flag allows intra calc to be skipped
  !            auxparams -- additional parameters (maxnrg,timestep,T_Kelvin)
  !            subset1 -- 1st subset
  !            subset2 -- 2nd subset
  !---------------------------------------------------------------------------
  Logical Function subinteract_int(subinfo,species,simcell,fast, &
      recalc,skip_intra,auxparams,subset1,subset2)
    Type(Subset_Interactions), Intent(InOut)       :: subinfo
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell  
    Logical, Intent(In)                            :: fast,recalc,skip_intra
    Real(kind=RDbl), Dimension(:), Intent(In)      :: auxparams
    Integer, Dimension(3), Intent(In)              :: subset1
    Integer, Dimension(3), Intent(In), Optional    :: subset2

    Logical      :: intra

    !** Default
    subinteract_int = .False.

    !** Store the subset information
    subinfo%subset1 = subset1
    If (Present(subset2)) subinfo%subset2 = subset2

    !** Set the parameters
    If (subinfo%update_mode == 1) Then
      Call storetop_setmap(subinfo%temp_try,.True.,subinfo%subset1)
    End If
    subinfo%skip_intra = skip_intra

    !** If needed, make sure the subset storage is sized properly
    Call subinteract_chksizes(subinfo,subinfo%subset1,recalc)

    !** Either extract results from storage or pass to reevaluation
    If (.Not. recalc) Then
      intra = (.Not. skip_intra)
      subinteract_int = storetop_extract(subinfo%imodel%results(1), &
          intra,.True.,.True.,subinfo%now_sum,subinfo%subset1,&
          subinfo%subset2)

      !** Return now if extraction was successful
      If (subinteract_int) Return

      !** Otherwise, do evaluation with the 'now' temporary storage
      subinteract_int = interact_int(subinfo%imodel,subinfo%temp_now,species,&
          simcell,fast,.True.,subinfo%skip_intra,subinfo%calc_accels, &
          auxparams,subinfo%subset1,subinfo%subset2)

      !** Copy the total into now_sum if this was a current configuration
      Call storebase_copy(subinfo%now_sum,subinfo%temp_now%total)

    Else 
      !** Do evaluation with 'try' temporary storage
      subinteract_int = interact_int(subinfo%imodel,subinfo%temp_try,species,&
          simcell,fast,.True.,subinfo%skip_intra,subinfo%calc_accels, &
          auxparams,subinfo%subset1,subinfo%subset2)

      !** Copy the total summed interactions information from temp storage
      If (subinfo%update_mode > 0) Then
        Call storebase_copy(subinfo%try_sum,subinfo%temp_try%total)
      End If

    End If

  End Function subinteract_int

  !---------------------------------------------------------------------------
  ! In rare cases it is necessary to reevaluate all the system interactions.
  ! This routine serves this purpose.
  ! Requires:  subinfo -- subset interaction structure
  !            species -- coordinate data structure 
  !            simcell -- simulation cell information
  !            fast -- logical flag indicating fast or slow evaluations
  !            auxparams -- additional parameters (maxnrg,timestep,T_Kelvin)
  !---------------------------------------------------------------------------
  Logical Function subinteract_fullint(subinfo,species,simcell,fast,auxparams)
    Type(Subset_Interactions), Intent(InOut)       :: subinfo
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell  
    Logical, Intent(In)                            :: fast
    Real(kind=RDbl), Dimension(:), Intent(In)      :: auxparams

    !** Default
    subinteract_fullint = .False.

    !** Do evaluation with the full storage structure
    !** HACK, has to reach into imodel, create routine in interact for this
    subinteract_fullint = interact_int(subinfo%imodel, &
        subinfo%imodel%results(1),species,simcell,fast,.True.,.False., &
        subinfo%calc_accels,auxparams,(/0,0,0/))

  End Function subinteract_fullint

  !---------------------------------------------------------------------
  ! Returns the total energy for the existing system subset 
  ! Requires:  subinfo -- subset interactions structure
  !            no_intra -- tells to skip adding intra nrgs
  !---------------------------------------------------------------------
  Real(kind=RDbl) Function subinteract_oldnrg(subinfo, no_intra)
    Type(Subset_Interactions), Intent(In)    :: subinfo
    Logical, Intent(In)                      :: no_intra

    subinteract_oldnrg = subinfo%now_sum%nrg
    If (no_intra) Return

    !** Add the intramolecular portion to the returned quantity
    If (storebase_chkintra(subinfo%now_sum)) Then
      subinteract_oldnrg = subinteract_oldnrg + Sum(subinfo%now_sum%intranrg)
    End If
    
#ifdef DEBUG
    string = storebase_disp(subinfo%now_sum)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**debug_insert
    Write(*,*) 'existing: ',Trim(string)
#endif

  End Function subinteract_oldnrg

  !---------------------------------------------------------------------
  ! Returns the total energy for the temporary system subset 
  ! Requires:  subinfo -- subset interactions structure
  !            no_intra -- tells to skip adding intra nrgs
  !---------------------------------------------------------------------
  Real(kind=RDbl) Function subinteract_newnrg(subinfo,no_intra)
    Type(Subset_Interactions), Intent(In)     :: subinfo
    Logical, Intent(In)                       :: no_intra

    subinteract_newnrg = subinfo%try_sum%nrg 
    If (no_intra) Return

    !** Add the intramolecular portion to the returned quantity
    If (storebase_chkintra(subinfo%try_sum)) Then
      subinteract_newnrg = subinteract_newnrg + sum(subinfo%try_sum%intranrg)
    End If

#ifdef DEBUG
    string = storebase_disp(subinfo%try_sum)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(0,*) 'temporary: ',Trim(string)
#endif

  End Function subinteract_newnrg

  !---------------------------------------------------------------------------
  ! Zero the temporary storage structure(s)
  ! Requires:  subinfo -- Subset Interaction data type for specific species
  !            subset -- the subset for which the temp structure is intended
  !---------------------------------------------------------------------------
  Subroutine subinteract_zerotemp(subinfo)
    Type(Subset_Interactions), Intent(InOut)       :: subinfo

    Call storetop_zero(subinfo%temp_try,.False.)
    If (Associated(subinfo%temp_now)) Then
      Call storetop_zero(subinfo%temp_now,.False.)
    End If

  End Subroutine subinteract_zerotemp

  !--------------------------------------------------------------------------
  ! Changes the number of molecules in the full system interaction 
  ! storage structure.  Essentially just a cover for interact_changenmoles.
  ! This function returns False if it was required to make a deletion in
  ! the storage and was unable due to the form of the storage.  In this case,
  ! a full system interaction evaluation will be required.
  ! Requires:  subinfo -- subset interactions structure
  !            spc -- species number for which to change nmoles 
  !            nmoles -- new number of molecules
  !            delmol -- optional molecule to delete
  !--------------------------------------------------------------------------
  Logical Function subinteract_changenmoles(subinfo,spc,nmoles,delmol)
    Type(Subset_Interactions), Intent(InOut)    :: subinfo
    Integer, Intent(In)                         :: spc,nmoles
    Integer, Intent(In), Optional               :: delmol

    Logical            :: reevaluate

    !** set default
    subinteract_changenmoles = .True.
    reevaluate =  .False.

    If (Present(delmol)) Then
      !** Subtract the temp_now structure from the full storage if it exists
      If (subinfo%update_mode == 2) Then
        !** Determine if a full system reevaluation must be done
        reevaluate = (.Not. storetop_delpossible(subinfo%imodel%results(1)))

        If (.Not. reevaluate) Then
          Call storetop_subtract(subinfo%imodel%results(1),subinfo%temp_now)
        End If
      End If

      !** Change molecule position or sizing in storage structure
      If (reevaluate) Then
        Call interact_changenmoles(subinfo%imodel,spc,nmoles)
        subinteract_changenmoles = .False.
      Else
        Call interact_changenmoles(subinfo%imodel,spc,nmoles,delmol)
      End If
    Else
      Call interact_changenmoles(subinfo%imodel,spc,nmoles)
    End If

  End Function subinteract_changenmoles

  !----------------------------------------------------------------------------
  ! Copies the results from the %temp_try storage of subset to full system
  ! storage.  Operates in two alternative update modes.  It either
  ! subtracts the current ('now') interactions and adds trial interactions
  ! or it simply replaces the old interactions with the new.  Replacement
  ! is simpliest, but can only be done with storage that is more detailed
  ! than the subset being changed.  The types of perturbations done to the
  ! system are specified in the list of species numbers.  Examples:
  ! (/-1,1/) -- perturbation move on species 1 (added and subtracted)
  ! (/-2/) -- deletion move on species 2
  ! (/1/) -- insertion move on species 1
  ! (/-1,-3,2/) -- deletions for species 1,3 insertion of species 2
  !                (example is the situation during idchange moves)
  ! The 'species' here doesn't necessarily refer to a full species move.
  ! It is assumed that all of the subset interactions have the same full 
  ! system storage structure(s).
  ! Requires:  subints -- subset interactions structures
  !            spcs -- list of species, can be negative
  !----------------------------------------------------------------------------
  Subroutine subinteract_update(subints,spcs)
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Integer, Dimension(:), Intent(In)                      :: spcs

    Logical               :: success
    Integer               :: i,spc,nlisted

    !** Return now if the full and temporary storage are identical
    If (subints(1)%update_mode == 0) Return

    !** Get number of listed species
    nlisted = Size(spcs)

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'passed species number indicators: ',spcs
    Write(*,*) 'BEFORE sub-update: full forcefield storage:'
    Call storetop_display(subints(1)%imodel%results(1),.True.,2,6)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Do i = 1,Size(spcs)
      If (spcs(i) == 0) Cycle
      Write(*,*) 'going to update full storage with subset'
      Write(*,*) 'BEFORE sub-update: temp forcefield storage:'
      spc = Abs(spcs(i))
      Write(*,*) 'skip_intra flag: ',subints(spc)%skip_intra
      Write(*,*) 'subset1: ',subints(spc)%subset1
      If ((spcs(i) < 0).And.(subints(spc)%update_mode == 2)) Then
        Write(*,*) 'SUBTRACTING temp_now'
        Call storetop_display(subints(spc)%temp_now,.True.,2,6)
      Else
        Write(*,*) 'ADDING or REPLACING temp_try'
        Call storetop_display(subints(spc)%temp_try,.True.,2,6)
      End If
    End Do
#endif

    !** Update either in replacement or subtract/add mode
    If (subints(1)%update_mode == 1) Then
      !** Replace parts of full storage with the temporary storage
      Do i = 1,nlisted
        If (spcs(i) <= 0) Cycle
        spc = spcs(i)
        success = storetop_update(subints(1)%imodel%results(1), &
            subints(spc)%temp_try,subints(spc)%subset1,subints(spc)%skip_intra)
        Call checkandstop(success,__FILE__,__LINE__, &
            ' could not successfully update the storage ',(/spc/))
      End Do

    Else
      !** Subtract the 'now' storage or add the 'try' storage for each listed
      Do i = 1,nlisted
        If (spcs(i) < 0) Then
          spc = -spcs(i)          
          Call storetop_subtract(subints(1)%imodel%results(1), &
              subints(spc)%temp_now)
        Else If (spcs(i) > 0) Then
          spc = spcs(i)
          Call storetop_add(subints(1)%imodel%results(1), &
              subints(spc)%temp_try)
        Else   !** (== 0)
          Cycle
        End If

        !** Also copy the storage sums from temp to current
        Call storebase_copy(subints(spc)%now_sum,subints(spc)%try_sum)
      End Do

    End If

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'AFTER sub-update: full forcefield storage:'
    Call storetop_display(subints(1)%imodel%results(1),.True.,2,6)
#endif

  End Subroutine subinteract_update

  !-------------------------------------------------------------------------
  ! Update the forcefield results storage using a subset interaction
  ! structure.  This is done when a move is accepted The routine also
  ! updates the statistics
  !   The 2nd subinfo structure is passed when the subtract/add update mode
  ! is in use and when the 'try' and 'now' temporary structures are stored
  ! in different subset interaction structures.  An example is the situation 
  ! during idchange moves.
  ! Requires:  subints -- subset interactions structures
  !            species -- species data structure
  !            spcs -- list of species, can be negative
  !            copysubint -- if True, will copy subset interactions to full
  !-------------------------------------------------------------------------
  Subroutine subinteract_updateall(subints,species,spcs,copysubint)
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(In)            :: species
    Integer, Dimension(:), Intent(In)                      :: spcs
    Logical, Intent(In)                                    :: copysubint

    Integer             :: i,spc

    !** Copy from temp storage to full system storage 
    If (copysubint) Call subinteract_update(subints,spcs)      
    
    !** Update the statistics for each (+) species
    Do i = 1,Size(spcs)
      !** Skip a 2nd call for perturbed species (X,-X), but not the deleted (-X)
      spc = Abs(spcs(i))
      If ((spcs(i) <= 0).And.(findint(spcs,spc) /= 0)) Cycle
      Call interact_updatestats(subints(spc)%imodel,species, &
          subints(spc)%calc_accels)
    End Do

  End Subroutine subinteract_updateall

  !-------------------------------------------------------------------------
  ! A special, simple version of the update routines for full system MD-type
  ! moves.  It simply updates the statistics because it can be assumed
  ! that there is no temporary storage.  I had to make this separate because
  ! the MD-type moves don't have arrays of subset interaction types.
  ! Requires:  subint -- subset interactions structure
  !            species -- species data structure
  !-------------------------------------------------------------------------
  Subroutine subinteract_simpleupdate(subint,species)
    Type(Subset_Interactions), Intent(InOut)        :: subint
    Type(AtMolCoords), Dimension(:), Intent(In)     :: species

    !** Update the statistics
    Call interact_updatestats(subint%imodel,species,subint%calc_accels)

  End Subroutine subinteract_simpleupdate

  !---------------------------------------------------------------------------
  ! Make sure the temporary interaction storage is sized properly
  ! Requires:  subinfo -- Subset Interaction data type for specific species
  !            subset -- the subset for which the temp structure is intended
  !            recalc -- same as in subinteract_int call, True => 'now' calc
  ! HACK, write a routine in interact that does this and call here
  !---------------------------------------------------------------------------
  Subroutine subinteract_chksizes(subinfo,subset,recalc)
    Type(Subset_Interactions), Intent(InOut)       :: subinfo
    Integer, Dimension(3), Intent(In)              :: subset
    Logical, Intent(In)                            :: recalc

    Select Case(subinfo%update_mode)
    Case(0)   !** no temporary storage, do nothing
      Return

    Case(1)   !** replacement mode
      Call storetop_chksizes(subinfo%temp_try,subinfo%imodel%results(1),subset)

    Case(2)   !** subtract/add mode
      If (recalc) Then
        Call storetop_chksizes(subinfo%temp_try, &
            subinfo%imodel%results(1),(/0,0,0/))
      Else
        Call storetop_chksizes(subinfo%temp_now, &
            subinfo%imodel%results(1),(/0,0,0/))
      End If
    End Select

  End Subroutine subinteract_chksizes

  !---------------------------------------------------------------------------
  ! Check the ability of the interaction storage code to return the same
  ! mol-system energies from extraction from the whole storage and from
  ! evaluation of only a subset interaction.
  ! Requires:  subinfo -- Subset Interaction data type for specific species
  !            spc -- species number
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            indent -- no. of spaces from the left margin
  !            unit -- optional unit number, will dump info only if present
  !---------------------------------------------------------------------------
  Subroutine subinteract_chksubint(subinfo,spc,species,simcell,indent,unitno)
    Type(Subset_Interactions), Intent(InOut)       :: subinfo
    Integer, Intent(In)                            :: spc
    Type(AtMolCoords), Dimension(:), Intent(Inout) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Integer, Intent(In)                            :: indent
    Integer, Intent(In), Optional                  :: unitno

    Integer               :: m
    Logical               :: fast,success,skip_intra_flag
    Character(len=indent) :: blank
    Character(len=strLen) :: string
    Real(kind=RDbl)       :: nrg_devn,nrg_recalc,nrg_stored,max_devn,maxnrg
    Type(EnergyPlus)      :: nrgs

    blank = Repeat(' ',indent)
    max_devn = 0.0_RDbl
    maxnrg = 1000.0_RDbl

    !** Set flags governing calculation
    fast = .True.
    skip_intra_flag = .False.

    If (Present(unitno)) Then
      Write(unitno,'(2a)') blank,'Checking storage vs recalc reproducibility'
      Write(unitno,'(a,2a5,3a12)') blank,'spc','mol','stored','recalc','diff'
    End If

    !** Loop through species to check all molecule-system interactions
    Do m = 1,config_getnmoles(species,spc)
      If (config_isfixed(species(spc))) Cycle
    
      !** Get the stored value using extraction
      success = subinteract_int(subinfo,species,simcell,fast, &
          .False.,skip_intra_flag,(/maxnrg/),(/spc,m,0/))
      Call checkandstop(success,__FILE__,__LINE__,' nrg calc')
      nrg_stored = subinteract_oldnrg(subinfo,skip_intra_flag)
    
      !** Get the recalculated value
      success = subinteract_int(subinfo,species,simcell,fast, &
          .True.,skip_intra_flag,(/maxnrg/),(/spc,m,0/))
      Call checkandstop(success,__FILE__,__LINE__,' nrg calc')
      nrg_recalc = subinteract_newnrg(subinfo,skip_intra_flag)
    
      nrg_devn = Abs(nrg_stored - nrg_recalc)
      If (nrg_devn > max_devn) max_devn = nrg_devn
!      If (nrg_devn > 1.0e-8_RDbl) Then
        If (Present(unitno)) Then
          Write(unitno,'(a,2i5,3f12.5,a)') blank,spc,m,nrg_stored, &
              nrg_recalc,nrg_devn,' COMPARISON'
        End If
!      End If
    End Do

    !** Stop if the deviation is large enough
    If (max_devn > 1.0e-7_RDbl) Then
      Write(0,'(1x,2a,i4,a)') __FILE__,": ",__LINE__, &
          ' ERROR: storage versus recalc shows unacceptable reproducibility'
      Stop
      Call storetop_display(subinfo%imodel%results(1),.True.,2,6)
    End If

  End Subroutine subinteract_chksubint

  !---------------------------------------------------------------------------
  ! Check the correspondence between the number of molecules of each species
  ! in the system with the numbers that the storage structure thinks it
  ! contains.  This is only used for debugging.
  ! Requires:  subinfo -- Subset Interaction data type for specific species
  !            species -- species data structure
  !            indent -- no. of spaces from the left margin
  !            unit -- unit number
  !---------------------------------------------------------------------------
  Subroutine subinteract_chkmolnums(subinfo,species,indent,unitno)
    Type(Subset_Interactions), Intent(In)       :: subinfo
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: indent
    Integer, Intent(In)                         :: unitno

    Call interact_chkmolnums(subinfo%imodel,species,indent,unitno)

  End Subroutine subinteract_chkmolnums

  !---------------------------------------------------------------------
  ! Display the subset interaction structure
  ! Requires:  subinfo -- subset interaction structure
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !---------------------------------------------------------------------
  Subroutine subinteract_display(subinfo,indent,unitno)
    Type(Subset_Interactions), Intent(In)  :: subinfo
    Integer, Intent(In)                    :: indent
    Integer, Optional, Intent(In)          :: unitno

    Integer                           :: unit,i
    Character(len=indent)             :: blank

    blank = Repeat(' ',indent)
    
    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    Write(unitno,'(4a)') blank,'The "',Trim(subinfo%id), &
        '" subset interaction information:'
    Write(unitno,'(3a)') blank,'Using interaction model: ', &
        Trim(subinfo%imodel%keys)
    Write(unitno,'(2a,i1,a,i1)') blank,'Level 1: ',subinfo%depth1, &
        '    Level 2: ',subinfo%depth2
    Write(unitno,'(2a,l2)') blank,'Calculate Accelerations? ',subinfo%calc_accels
    Write(unitno,'(2a,l2)') blank,'Update mode is ',subinfo%update_mode

    Write(unitno,'(2a,3i5)') blank,'Subset 1 (spc,mol,atm): ', &
        (subinfo%subset1(i),i=1,3)
    Write(unitno,'(2a,3i5)') blank,'Subset 2 (spc,mol,atm): ', &
        (subinfo%subset2(i),i=1,3)
    Write(unitno,'(2a,f8.3)') blank,'Current Subset Energy:   ', &
        subinfo%now_sum%nrg
    Write(unitno,'(2a,f8.3)') blank,'Temporary Subset Energy: ', &
        subinfo%try_sum%nrg
    
    If (subinfo%update_mode > 0) Then
      Write(unitno,'(2a)') blank,'Temporary Storage:'
      Call storetop_display(subinfo%temp_try,.False.,indent+2,unitno)

      If (Associated(subinfo%temp_now)) Then
        Write(unitno,'(2a)') blank,'Temporary "now" Storage:'
        Call storetop_display(subinfo%temp_now,.False.,indent+2,unitno)
      End If
    End If

  End Subroutine subinteract_display

  !---------------------------------------------------------------------
  ! Cleans the subset interaction structure
  ! Requires:  subinfo -- subset interactions structure
  !---------------------------------------------------------------------
  Subroutine subinteract_clean(subinfo)
    Type(Subset_Interactions), Intent(InOut)       :: subinfo

    Nullify(subinfo%imodel)
    Call storebase_clean(subinfo%now_sum)
    Call storebase_clean(subinfo%try_sum)
    Call storetop_clean(subinfo%temp_try)
    If (Associated(subinfo%temp_now)) Call storetop_clean(subinfo%temp_now)

  End Subroutine subinteract_clean

End module subinteract


