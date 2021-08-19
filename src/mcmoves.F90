!------------------------------------------------------------------------------
! This module handles the acceptance/rejection and implementation of BASIC
! Monte Carlo move types.  Current "basic" move types are:
!  1) Perturbation of current molecule (rotation, translation, geometry change)
!  2) Insertion of a new molecule
!  3) Deletion of an existing molecule
!  4) Identity change
! Within these basic types, there are different varieties such as biased, etc. 
!
! Any basic Monte Carlo move type that could possibly be used by more 
! than one simulation type should be contained in, or accessed through 
! this module.  This module uses the "moves" module to actually change 
! the state of the system.  The "moves" module returns all the information
! necessary to accept or reject a move, this includes energies.
!
! The classic Metropolis solution for the transition matrix plus the addition
! of biasing is used to accept or reject the configurations:
! P_mn = min[1, (alpha_nm/alpha_mn)*(rho_n/rho_m)],   Where:
!    P_mn = acceptance probility for moving from state m -> state n
!    rho_n = probability of being at state n, according to ensemble statistics
!    alpha_mn = probability of attempting move from state m -> state n
!    (alpha_nm/alpha_mn) = bias factor, returned from moves
! 
! The module also contains a data type for auxilliary rejection parameters.
! This data type and the routines associated with it could be relocated to
! their own module if necessary.
!
! Important routines and move types provided here are:
!    mcmoves_init -- handles initilization 
!    mcmoves_perturb -- handles Metropolis-accepted perturbations on
!                       existing molecules, i.e. rotations/translations
!    mcmoves_insert -- handles insertions
!    mcmoves_delete -- handles deletions 
!    mcmoves_volchange -- handles volume changes (doesn't exist yet)
!    mcmoves_idchange -- handles identity changes
!    mcmoves_stats -- display one line of intermediate statistics
!    mcmoves_display -- displays information about the move type
!
! For all move types, thermophysical properties and the "mcmove" data
! type are passed into the function and the logical return value
! indicates the success or failure of the move.  Attempt and success
! counters are kept within the "mcmove" data type.
!------------------------------------------------------------------------------

Module mcmoves
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl,strLen,MAX_EXCLSITES,MAX_SORBS,caltoj,one, &
      dashedline,zero,default_MAX_EXP,dbgflag,lstrLen
  Use utils, Only: checkandstop
  Use file, Only: file_getunit,file_open
  Use config, Only: Atmolcoords,config_getnmoles,config_getnatoms,&
      config_allocfields,config_checkandincr,config_setnmoles,&
      config_delmol,config_config2xyz,config_copysubset, config_dumpmol, &
      config_isfixed, config_subsysKE
  Use subinteract, Only: Subset_Interactions, subinteract_changenmoles, &
      subinteract_oldnrg, subinteract_newnrg, subinteract_updateall, &
      subinteract_zerotemp, subinteract_chksizes, subinteract_fullint
  Use moves, Only: Move_Params,Move_stats,moves_initparams,moves_gettag, &
      moves_makemove,moves_displayparams,moves_adjustparams,moves_movestats, &
      moves_restore,moves_movestatsdisplay,moves_dyninfo, moves_setMaxnrg, &
      moves_checkoppmoves, moves_postadjust, &
      moves_getparam, moves_checkIntra
  Use gcmodels, Only: gcmodels_getcom, gcmodels_dumprestartinfo
  Use molecules, Only: molecules_getnatoms, molecules_getnsorbs, molecules_name
  Use storestats, Only: Species_Stats, storestats_incrnoncoulnrg, &
      storestats_incrAllIntraNrg,storestats_incrcoulnrg
  Use storetop, Only: storetop_display, storetop_chksizes, storetop_zero
  Use random, Only: rranf
  Use simcell, Only: SimCell_Params,simcell_pbc
  Use smap, Only: Smap_Params,smap_init,smap_getPosSiteType
  Use storebase, Only: storebase_display
  Use utils, Only: toupper,split,stripcmnt,toint,allocErrDisplay,int2str, &
      cleanstring,real2str
  Use vector, Only: VecType, Assignment(=)
  
  Implicit None
  Save

  Private
  Public :: MC_move_Params,AuxReject_Params,mcmoves_init,mcmoves_getbasictag, &
      mcmoves_perturb,mcmoves_insert,mcmoves_delete,mcmoves_idchange, &
      moves_restore,mcmoves_initaux, mcmoves_stats,mcmoves_auxrejectdisplay,&
      mcmoves_display,mcmoves_setMaxNrg,mcmoves_resetStats, &
      mcmoves_checkOppMoves,mcmoves_postadjust,mcmoves_checkandincr

  !** These are auxilliary parameters currently used as additional
  !** criteria for move rejection.  They could be in their own module
  Type AuxReject_Params
    Integer                           :: nexclsites    
    Integer, Dimension(:), Pointer    :: exclsitetype
    Character(len=strLen)             :: smapname
    Type(Smap_Params), Pointer        :: smap
  End Type AuxReject_Params

  Type MC_Move_Params
    Integer               :: spc           !* molecule type
    Character(len=strLen) :: mvTag         !* full name of move type
    Character(len=strLen) :: basic_tag     !* basic name of move type
    Type(Move_Stats)      :: stats      
    Type(Move_Params)     :: mvparams 
    Integer               :: ensemble
    Type(AtMolCoords),Dimension(:), Pointer     :: oldcoords
    Integer, Dimension(3)                       :: oldcoords_subset
    Type(AuxReject_Params), Pointer    :: auxparams
    Logical                            :: use_intra
  End Type MC_Move_Params

  !** Integers corresponding to different ensembles
  !** I like to use these integers rather than using string like "NVT"
  Integer, Parameter :: RANDOM_ENS = -1
  Integer, Parameter :: NVE_VELS_ENS = 0
  Integer, Parameter :: NVT_ENS = 1
  Integer, Parameter :: MuVT_ENS = 2
  
Contains
  !---------------------------------------------------------------------
  ! Initializes the MC move parameters from the control file
  ! Requires:  mcmove -- general MC move parameters 
  !            filename -- the name of the control file
  !            spc -- species number
  !            simcell -- the simulation cell information
  !            species -- species data structure
  !            ens_string -- string indicating the ensemble
  !            auxparams -- optional auxillaiary params for rejection
  ! Note : Pass "spc" only if species -specific move
  !---------------------------------------------------------------------
  Subroutine mcmoves_init(mcmove, filename, simcell, species, ens_string, &
      auxmv, spc, auxparams )
    Type(MC_Move_Params), Intent(InOut)           :: mcmove
    Character(*), Intent(In)                      :: filename
    Type(SimCell_Params), Intent(In)              :: simcell
    Type(AtMolCoords), Dimension(:), Intent(In)   :: species
    Character(*), Intent(In)                      :: ens_string
    Type(AuxMoveObjects), Pointer           :: auxmv 
    Integer, Intent(In), Optional                 :: spc
    Type(AuxReject_Params), Pointer, Optional     :: auxparams

    Integer            :: nsorbs,error, i

    !** Initialize the move parameters
    If (Present(spc)) Then
      mcmove%spc = spc
      Call moves_initparams(mcmove%mvparams, auxmv, simcell,spc,filename, &
          mcmove%basic_tag,species)
    Else
      mcmove%spc = 1
      Call moves_initparams(mcmove%mvparams, auxmv, species,filename, &
          mcmove%basic_tag)
    Endif

    !** Set the ensemble
    Select Case(cleanstring(ens_string))
    Case("RANDOM")
      mcmove%ensemble = RANDOM_ENS
    Case("NVT")
      mcmove%ensemble = NVT_ENS
    Case("NVE_VELS")
      mcmove%ensemble = NVE_VELS_ENS
    Case("MuVT")
      mcmove%ensemble = MuVT_ENS
    Case Default
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' Could not interpret passed ensemble type string "', &
          Trim(cleanstring(ens_string)),'"'
      Stop      
    End Select

    nsorbs = Size(species,1)
    mcmove%mvTag = moves_gettag(mcmove%mvparams)

    !** Allocate oldcoords based on whether system move or species-
    !** specific move. This might not be used in all movetypes
    If (Trim(mcmove%mvTag) == "INTEGRATE") Then
      mcmove%oldcoords_subset=(/ 0,0,0/)
      Allocate(mcmove%oldcoords(nsorbs), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Do i = 1,nsorbs
        !** default amount of memory to store 1 molecule
        Call config_allocfields(mcmove%oldcoords(i),i,1)
        Call config_setnmoles(mcmove%oldcoords(i),1)
      End Do
    Else
      mcmove%oldcoords_subset=(/1,1,0/)
      Allocate(mcmove%oldcoords(1), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      Call config_allocfields(mcmove%oldcoords(1),spc,1)
      Call config_setnmoles(mcmove%oldcoords(1),1)
    End If

!LC    Call moves_displayparams(mcmove%mvparams,5,6)

    !** set all stat counter to zero
    Call mcmoves_resetStats(mcmove)

    Nullify(mcmove%auxparams)
    If (Present(auxparams)) Then
      mcmove%auxparams => auxparams
    End If

    ! checks whether we need intra nrgs in acceptance criteria or not
    mcmove%use_intra=.False.
    Call moves_checkIntra(mcmove%mvparams, mcmove%use_intra)

  End Subroutine mcmoves_init

  !-----------------------------------------------------------------------
  ! Retrieves the basic tag name from the data structure
  ! Requires: mcmove -- general move type parameters
  !-----------------------------------------------------------------------
  Function mcmoves_getbasictag(mcmove)
    Character(len=strLen)              :: mcmoves_getbasictag
    Type(MC_Move_Params), Intent(In)   :: mcmove

    mcmoves_getbasictag = mcmove%basic_tag

  End Function mcmoves_getbasictag

  !----------------------------------------------------------------------------
  ! Performs the updates necessary after a successful move.  Note that
  ! input energies are kcal/mol and stored energies are kJ/mol.  Use the
  ! 'copysubint' flag to update the full system interaction from the subset 
  ! interaction ONLY if an insertion or perturbation move was accepted.
  ! For deletion moves, the subset interactions aren't filled in because
  ! the subset interactions contains information for a molecule that has been
  ! removed.  The types of perturbations done to the
  ! system are specified in the list of species numbers.  Examples:
  ! (/-1,1/) -- perturbation move on species 1 (added and subtracted)
  ! (/-2/) -- deletion move on species 2
  ! (/1/) -- insertion move on species 1
  ! (/-1,-3,2/) -- deletions for species 1,3 insertion of species 2
  !                (example is the situation during idchange moves)
  ! The 'species' here doesn't necessarily refer to a full species move.
  ! Requires:  mcmove -- general move type parameters
  !            subints -- subset interactions for all species
  !            species -- species data structure
  !            spcs -- species numbers to update, can be negative, see above
  !            copysubint -- if True, will copy subset interactions to full
  !----------------------------------------------------------------------------
  Subroutine mcmoves_accept(mcmove,subints,species,spcs,copysubint)
    Type(MC_Move_Params), Intent(InOut)                    :: mcmove
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(In)            :: species
    Integer, Dimension(:), Intent(In)                      :: spcs
    Logical, Intent(In)                                    :: copysubint

    !** Update the full forcefield structure with subset and update statistics
    Call subinteract_updateall(subints,species,spcs,copysubint)

    !** Update the success counter
    mcmove%stats%succ = mcmove%stats%succ + 1

  End Subroutine mcmoves_accept

  !-----------------------------------------------------------------------
  ! Performs the updates necessary after an UNsuccessful move
  ! Requires:  mcmove -- general move type parameters
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !-----------------------------------------------------------------------
  Subroutine mcmoves_reject(mcmove,subints,species)
    Type(MC_Move_Params), Intent(InOut)                    :: mcmove
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(In)            :: species

    !** Update the statistics
    Call subinteract_updateall(subints,species,(/mcmove%spc/),.False.)

  End Subroutine mcmoves_reject

  !------------------------------------------------------------------------
  ! Attempts the move.  In other words, it calculates the acceptance
  ! probability, generates a random number and returns True if the
  ! move should be accepted.
  ! P_mn = min[1, (alpha_nm/alpha_mn)*(rho_n/rho_m)]
  !        min[1, (alpha_nm/alpha_mn)*prefactor*exp(addxtra - deltaE/RT)]
  !        min[1, biasfactor*prefactor*exp(addxtra - deltaE*rti)]
  ! Requires:  mcmove -- general move type parameters
  !            deltaE -- energy difference for Boltzmann factor
  !            rti -- 1.0/(Rgas*tk)
  !            biasfactor -- the bias factor
  !            addxtra -- extra term to add to exponential
  !            prefactor -- extra factor to multiply by
  ! Note: definition of biasfactor - factor used to remove the bias towards 
  ! the newly generated config. For example, if new posn A was selected 
  ! with probability Exp(-U_A)/ \sigma Exp(-U_i) from a set of N positions, 
  ! Then the biasfactor is N* \sigma Exp(-U_i) /Exp(-U_A)
  ! If its the deletion move then biasfactor is inverse of the above.
  !------------------------------------------------------------------------
  Logical Function mcmoves_attempt(mcmove,deltaE,rti,biasfactor, &
      addxtra,prefactor)
    Type(MC_Move_Params), Intent(In)       :: mcmove
    Real(kind=RDbl), Intent(In)            :: deltaE,rti,biasfactor
    Real(kind=RDbl), Intent(In), Optional  :: addxtra,prefactor

    Real(kind=RDbl)        :: boltz,utemp

    !** Accept immediately if this is a random ensemble
    If (mcmove%ensemble == RANDOM_ENS) Then
      mcmoves_attempt = .True.
      Return
    End If
      
    !** Acceptance/Rejection Calculation
    utemp = -rti * deltaE
    utemp = utemp + Log(biasfactor)

    !** Add the extra term if it's present
    If (Present(addxtra)) utemp = utemp + addxtra
    
    !** Avoid overflow/underflow errors during exponentiation
    If (utemp > default_MAX_EXP) utemp = default_MAX_EXP   
    If (utemp < -default_MAX_EXP) utemp = -default_MAX_EXP

    !** Calculate Boltzmann factor
    boltz  = Exp(utemp)

    !** Modify Boltzmann factor if necessary
    If (Present(prefactor)) boltz = boltz*prefactor

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,'(a,e14.4)') 'delta energy: ',deltaE
    Write(*,'(a,e14.4)') 'utemp: ',utemp
    Write(*,'(a,e14.4)') 'Acceptance probability: ',boltz
#endif

    !** Generate random number and accept or reject move
    If (boltz > rranf()) Then
      mcmoves_attempt = .True.
    Else
      mcmoves_attempt = .False.
    End If

  End Function mcmoves_attempt

  !-------------------------------------------------------------------
  ! Handles any move type that involves a perturbation of an existing
  ! molecules coordinates.  It uses the simple Metropolis acceptance
  ! criteria and therefore works for rotation and translation moves.
  ! Recently it has been upgraded to deal with perturbation of part of system
  ! rather that just the molecule
  ! Requires: mcmove -- general move type parameters
  !           molec -- molecule to attempt deletion on
  !           rti -- 1.0/(Rgas*tk)
  !           subints -- subset interactions for each species
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           knrgflag -- if .True. include kinetic_energy in the acceptance
  !                       criteria 
  !-------------------------------------------------------------------
  Logical Function mcmoves_perturb(mcmove, subset, rti, subints, species, &
      simcell, knrgflag)
    Type(MC_Move_Params), Intent(InOut)                    :: mcmove
    Integer, Dimension(:), Intent(in)                      :: subset
    Real(kind=RDbl), Intent(In)                            :: rti  
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Logical, Intent(in)                                    :: knrgflag

    Integer                :: nmoles,natoms,spc, molec
    Logical                :: succ_flag, accept, no_intra
    Real(kind=RDbl)        :: biasfactor, old_ke, new_ke, delta_ke
    Real(kind=RDbl)        :: uinit,ufinal,deltaU,boltz,utemp,ratio
    Character(len=strLen)  :: comment,string1,string2,string3

    !** Set default
    mcmoves_perturb = .False.

    If (knrgflag) old_ke = config_subsysKE(species, subset)

    !** Store the old configuration
    Call config_copysubset(species, mcmove%oldcoords, subset, &
        mcmove%oldcoords_subset)

    !** Update the attempt statistics
    mcmove%stats%att = mcmove%stats%att + 1

    !** Make the move perturbation using the moves module
    succ_flag = moves_makemove(mcmove%mvparams,subset,subints,species,simcell, &
        biasfactor)

    If (knrgflag) new_ke = config_subsysKE(species, subset)

    !** Do what we need to do based on ensemble
    Select Case(mcmove%ensemble)
    Case(NVE_VELS_ENS)
      !** NVE ensemble, always accept
      accept = .True.
      If (.Not.(succ_flag)) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
            ' Something wrong with NVE ensemble move'
        Stop
      End If

    Case(NVT_ENS, MuVT_ENS, RANDOM_ENS)
      !** Check the criteria that could lead to automatic rejection
      accept = .True.
      If (.Not.(succ_flag)) Then
        !** Check if the energy calc. was successful.  If not, reject move
        accept = .False.
      Else If (Associated(mcmove%auxparams)) Then
        !** Check the auxillary rejection criteria and reject if appropriate
        !** note that this is equivalent to assigning a very high energy
        If (mcmoves_auxreject(mcmove%auxparams,subset,species)) &
            accept = .False.
      End If

      !** If preliminary check say the config is OK, apply Metropolis criteria
      If (accept) Then
        ! we could use mcmove%use_intra, 
        ! but need to test for translate,rotate,integrate cases
        no_intra = .not.mcmove%use_intra
        ufinal = subinteract_newnrg(subints(mcmove%spc),no_intra)
        uinit  = subinteract_oldnrg(subints(mcmove%spc),no_intra)

        !** Acceptance/Rejection Calculation
        deltaU = ufinal - uinit
        If (knrgflag) Then
          !** change kinetic energy to kcals
          delta_ke = (new_ke - old_ke)/calToj
          deltaU = deltaU + delta_ke
        End If

        accept = mcmoves_attempt(mcmove,deltaU,rti,biasfactor)

#ifdef DEBUG
        string1 = real2str(ufinal,12)
        string2 = real2str(uinit,12) 
        string3 = real2str(deltaU,12)
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,'(6a)') 'perturb energies(final-init): ' &
            ,Trim(string1),' - ',Trim(string2),' = ',Trim(string3)
#endif
      End If

    Case Default
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' Could not interpret ensemble number during mcmove: ',mcmove%ensemble
      Stop

    End Select

    mcmoves_perturb = accept
    If (accept) Then
      Call mcmoves_accept(mcmove,subints,species, &
          (/-mcmove%spc,mcmove%spc/),.True.)
    Else
      Call mcmoves_reject(mcmove,subints,species)
!!$
!!$      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!!$      Write(*,*) subset,mcmove%oldcoords_subset

      !** Restore the original coordinates
      Call config_copysubset(mcmove%oldcoords, species, &
          mcmove%oldcoords_subset, subset)
    End If

    !** Adjust the jump length based on the acceptance ratio
    Call moves_adjustParams(mcmove%mvparams,mcmove%stats)

  End Function mcmoves_perturb

  !-------------------------------------------------------------------
  ! Handles selection and acceptance or rejection of INSERTIONS
  ! Requires: mcmove -- general move type parameters
  !           molec -- RETURNED new molecule number
  !           rti -- 1.0/(Rgas*tk)
  !           B -- Adam's B parameter (contains fugacity)
  !           subints -- subset interactions for each species
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           insdelratio -- ratio of insertions /deletion
  !-------------------------------------------------------------------
  Logical Function mcmoves_insert(mcmove,molec,rti,B,subints,species,simcell, &
      insdelratio)
    Type(MC_Move_Params), Intent(InOut)                    :: mcmove
    Integer, Intent(Out)                                   :: molec
    Real(kind=RDbl), Intent(In)                            :: rti,B
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Optional, Intent(In)                  :: insdelratio

    Integer                :: nmoles,natoms,spc
    Logical                :: succ_flag, accept, no_intra
    Real(kind=RDbl)        :: ufinal,boltz,utemp,biasfactor

    mcmoves_insert = .False.
    spc = mcmove%spc
    nmoles = config_getnmoles(species,spc)
    natoms = config_getnatoms(species,spc)

    !** Update the attempt statistics
    mcmove%stats%att = mcmove%stats%att + 1
    accept = .True.

    !** Get the new molecule number for insertion
    molec = nmoles + 1
    Call config_setnmoles(species, spc, molec)

    !** Check if we need to increment the size of the configuration storage
    !** increases if we are less by 1 
    Call config_checkandincr(species,spc,nmoles)

    !** Increment the interaction storage by one molecule
    succ_flag = subinteract_changenmoles(subints(mcmove%spc),spc,nmoles+1)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

    !** Make the move perturbation using the moves module
    succ_flag = moves_makemove(mcmove%mvparams,(/spc,molec,0/),subints,&
        species,simcell,biasfactor)

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    If (.Not. succ_flag) Write(*,'(a)') 'insertion move unsucessful'
#endif

    !** First pass, accept stays True if there were no big problems
    If (.Not. succ_flag) Then 
      accept = .False.

    Else If (Associated(mcmove%auxparams)) Then
      !** Check the auxillary rejection criteria and reject if appropriate
      If (mcmoves_auxreject(mcmove%auxparams,(/spc,molec,0/),species)) Then
        accept = .False.
      End If
    End If

    !** Second pass, calculate the acceptance probability and make attempt
    If (accept) Then
      no_intra =  .Not.(mcmove%use_intra) !** 
      ufinal = subinteract_newnrg(subints(mcmove%spc), no_intra )

#ifdef DEBUG
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,'(a,f8.3)') 'insert energy: ',ufinal
#endif

      !** Adjust biasfactor for different ins/del ratio
      If (present(insdelratio)) biasfactor = biasfactor/insdelratio

      !** Acceptance/Rejection Calculation
      ! p = [1, 1/(N+1) exp(B - del V / RT) ]
      ! See Snurr et al. JPC (1993) vol 97, pages 13742-13752
      ! Note that as far as V is the external potential and the Ideal 
      ! EOS is used hybrid-GCMC also has same acceptance criteria 
      accept = mcmoves_attempt(mcmove,ufinal,rti,biasfactor, &
          B,1.0_RDbl/Real(molec))
    End If

    !** Finalize the move
    If (accept) Then
      Call mcmoves_accept(mcmove,subints,species,(/mcmove%spc/),.True.)
      mcmoves_insert = .True.
    Else
      Call mcmoves_reject(mcmove,subints,species)

      !** Decrease number of molecules in configuration and interaction storage
      Call config_setnmoles(species,spc,molec-1)
      succ_flag = subinteract_changenmoles(subints(mcmove%spc),spc,nmoles)
      Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')
    End If

  End Function mcmoves_insert

  !-------------------------------------------------------------------
  ! Handles selection and acceptance or rejection of DELETIONS
  ! Requires: mcmove -- general move type parameters
  !           mol -- molecule to attempt deletion on
  !           rti -- 1.0/(Rgas*tk)
  !           B -- Adam's B parameter (contains fugacity)
  !           subints -- subset interactions for each species
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !-------------------------------------------------------------------
  Logical Function mcmoves_delete(mcmove,mol,rti,B,subints,species,simcell, &
      insdelratio)
    Type(MC_Move_Params), Intent(InOut)                    :: mcmove
    Integer, Intent(In)                                    :: mol
    Real(kind=RDbl), Intent(In)                            :: rti,B
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Optional ,Intent(in)                  :: insdelratio

    Integer                :: nmoles,natoms,spc
    Logical                :: succ_flag,accept,no_intra
    Real(kind=RDbl)        :: uinit,boltz,utemp,biasfactor

    mcmoves_delete = .False.
    spc = mcmove%spc
    nmoles = config_getnmoles(species,spc)
    natoms = config_getnatoms(species,spc)

    !** Update the attempt statistics
    mcmove%stats%att = mcmove%stats%att + 1

    !** Make the move perturbation using the moves module
    succ_flag = moves_makemove(mcmove%mvparams,(/spc,mol,0/),subints,&
        species,simcell,biasfactor)

    If (.Not. succ_flag) Then 
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' Unexpected move failure during deletion attempt'
      Stop      
    End If

    !** Get the energy of the deleted molecule
    no_intra = (.Not. mcmove%use_intra)
    uinit = subinteract_oldnrg(subints(mcmove%spc), no_intra) 
    utemp = -1.0_RDbl*uinit

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,'(a,f8.3)') 'delete energy: ',uinit
#endif

    !** Adjust for different ins/del ratio
    If (Present(insdelratio)) biasfactor = biasfactor*insdelratio

    !** Acceptance/Rejection Calculation
    !** p = min [1, N*exp(-B - del V / RT) ]
    !** See Snurr et al. JPC (1993) vol 97, pages 13742-13752
    accept = mcmoves_attempt(mcmove,utemp,rti,biasfactor, &
        -1.0_RDbl*B, 1.0_RDbl*nmoles)
          
    !** Finalize the move
    If (accept) Then
      !** Decrease number of molecules in configuration and interaction storage
      Call config_delmol(species,spc,mol)
      succ_flag = subinteract_changenmoles(subints(mcmove%spc),spc,nmoles-1,mol)
      If (.Not. succ_flag) Then
        succ_flag = subinteract_fullint(subints(mcmove%spc),species,simcell, &
            .True.,(/10000.0_RDbl/))
        Call checkandstop(succ_flag,__FILE__,__LINE__,'full system eval failed')
      End If

      !** Update statistics
      Call mcmoves_accept(mcmove,subints,species,(/-mcmove%spc/),.False.)
      mcmoves_delete = .True.
    Else
      Call mcmoves_reject(mcmove,subints,species)
    End If

  End Function mcmoves_delete

  !-------------------------------------------------------------------
  ! Handles selection and acceptance or rejection of identity changes
  ! Requires:  mcmove -- general move type parameters
  !            subset -- definition of subset of which to change ID
  !            rti -- 1.0/(Rgas*tk)
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            B_term -- ratio of fugacities for gcmc 
  ! NOTE: it takes a subset by only operates on a molecule
  !-------------------------------------------------------------------
  Logical Function mcmoves_idchange(mcmove, subset, rti, subints, &
      species, simcell, B_term_array)
    Type(MC_Move_Params), Intent(InOut)                    :: mcmove
    Integer, Dimension(3), Intent(In)                      :: subset
    Real(kind=RDbl), Intent(In)                            :: rti  
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Dimension(:), Intent(In),Optional     :: B_term_array

    Integer                :: nmoles1,nmoles2,natoms,spc,i,newspc
    Logical                :: succ_flag,accept,no_intra
    Real(kind=RDbl)        :: ufinal,uinit,deltaE,boltz,biasfactor,realnum
    Real(kind=RDbl)        :: molratio, B_term
    Integer, Dimension(3)  :: spcs
    Character(len=strLen)  :: string1,string2,string3

    mcmoves_idchange = .False.
    spc = subset(1)
    nmoles1 = config_getnmoles(species,spc)
    natoms = config_getnatoms(species,spc)

    !** Update the attempt statistics
    mcmove%stats%att = mcmove%stats%att + 1
    accept = .True.

    !** Make the move perturbation using the moves module
    succ_flag = moves_makemove(mcmove%mvparams,subset,subints,&
        species,simcell,biasfactor)

    !** First pass, accept stays True if there were no big problems
    If (.Not. succ_flag) Then 
      accept = .False.

    Else If (Associated(mcmove%auxparams)) Then
      !** Check the auxillary rejection criteria and reject if appropriate
      If (mcmoves_auxreject(mcmove%auxparams,subset,species)) Then
        accept = .False.
      End If
    End If

    !** Second pass, calculate the acceptance probability and make attempt
    If (accept) Then
      !** Get the changed species number, following subinteract_update convention
      Call moves_getparam(mcmove%mvparams,'SPCS',realnum,spcs)

      !** Get energies and use intramolecular energy if it's there
      no_intra = .False.  
      uinit = 0.0_RDbl
      ufinal = 0.0_RDbl
      Do i = 1,Size(spcs)
        If (spcs(i) < 0) Then
          uinit = uinit + subinteract_oldnrg(subints(-spcs(i)),no_intra)
          !** remove neg numbers so _update doesn't sub, done in _changenmoles
          spcs(i) = 0   
        Else If (spcs(i) > 0) Then
          ufinal = ufinal + subinteract_newnrg(subints(spcs(i)),no_intra)
          newspc = spcs(i)
        End If
      End Do
      deltaE = ufinal - uinit

#ifdef DEBUG
      string1 = real2str(ufinal,12)
      string2 = real2str(uinit,12)
      string3 = real2str(deltaE,12)
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,'(6a)') 'idchange energies(final-init): ' &
            ,Trim(string1),' - ',Trim(string2),' = ',Trim(string3)
#endif

      !** Acceptance/Rejection Calculation
      !** B_term = B of new - B of old
        If (Present(B_term_array)) Then
        B_term = B_term_array(newspc) - B_term_array(mcmove%spc)
        nmoles2 = config_getnmoles(species,newspc)
        molratio = Real(nmoles1)/Real((nmoles2))
        !** gcmc type of move
        accept = mcmoves_attempt(mcmove,deltaE,rti,biasfactor,B_term,molratio)
      Else
        !** some other ensemble
        accept = mcmoves_attempt(mcmove,deltaE,rti,biasfactor)
      End If

    End If

    If (accept) Then
      !** Erase specified molecule's interaction storage
      succ_flag = subinteract_changenmoles(subints(subset(1)),subset(1), &
          nmoles1-1,subset(2))
      Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

      Call mcmoves_accept(mcmove,subints,species,spcs,.True.)
      mcmoves_idchange = .True.
    Else
      Call mcmoves_reject(mcmove,subints,species)
      Call moves_restore(mcmove%mvparams,subset,subints,species,simcell)
    End If

  End Function mcmoves_idchange

  !---------------------------------------------------------------------
  ! Increses the size of oldcoords and storage(temp) so that the subset
  ! values can be copied into these.  Note that this routine has an 
  ! expensive call to storetop_zero, so should not be called very often.
  ! Requires:  mcmove -- general move type parameters
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            subset -- the system subset for this move 
  !---------------------------------------------------------------------
  Subroutine mcmoves_checkandincr(mcmove,subint,species,subset)
    Type(MC_Move_Params), Intent(InOut)              :: mcmove
    Type(Subset_Interactions), Intent(InOut)         :: subint
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Integer, Dimension(3), Intent(In)                :: subset

    Integer               :: nspc, nmoles, spc

    If (subset(1) == 0) Then   !** Subset is whole system
      nspc = Size(species)
      If (Size(mcmove%oldcoords) /= nspc) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(a)') 'dimensions of species and oldcoords do not match'
        Stop
      End If

      !** Make sure the configuration storage is sized properly
      Do spc = 1,nspc
        If (.Not.(config_isfixed(species(spc)))) Then
          nmoles = config_getnmoles(species,spc)

          !** Check if we need to increment the size of the config storage
          Call config_checkandincr(mcmove%oldcoords , spc, nmoles, .False. )
        End If
      End Do

      !** Make sure the interaction storage is sized properly
      Call subinteract_chksizes(subint,subset,.False.)

      !** Zero the temporary storage structure(s)
      Call subinteract_zerotemp(subint)

    Else If (subset(1) > 0) Then   !** Subset is less than whole system

      !** check molecule number, nothing much to be done here
      If (subset(2) < 1) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(a)') 'copying one whole species not yet implemented'
        Stop
      End If

    Else
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

  End Subroutine mcmoves_checkandincr

  !-------------------------------------------------------------------
  ! Initializes the auxilliary criteria from file
  ! Requires: params -- auxilliary rejection params 
  !           filename -- filename to read from
  !           simcell -- the simulation cell information
  !-------------------------------------------------------------------
  Subroutine mcmoves_initaux(params,filename,simcell)
    Type(AuxReject_Params), Pointer                :: params
    Character(*), Intent(In)                       :: filename
    Type(SimCell_Params), Intent(In)               :: simcell

    Integer                           :: i,unitno,error
    Character(len=strLen)             :: dummy_string,line
    Character(len=strLen), Dimension(MAX_EXCLSITES) :: fields

    unitno = file_getunit(filename)

    If (.Not. Associated(params)) Then
      Allocate(params, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'params')
    End If
    Nullify(params%exclsitetype)
    Nullify(params%smap)

    !** Get sitemap name
    Read(unitno,*) dummy_string
    params%smapname = stripcmnt(dummy_string)

    !** Initialize the sitemap if it is not null
    If (toupper(params%smapname) /= "NULL") Then
      Call smap_init(params%smap,simcell,params%smapname)
      !** Read the line containing the list of sites to be excluded
      Read(unitno,'(a)') dummy_string
      line = stripcmnt(dummy_string)
    Else
      Nullify(params%smap)
    End If

    !** Define the sites to exclude
    If (Associated(params%smap)) Then
      If (Trim(toupper(line)) == 'NULL') Then
        params%nexclsites = 0
      Else
        params%nexclsites = split(line, fields, ",")
        Allocate(params%exclsitetype(params%nexclsites), STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'exclsitetype')
        Do i = 1,params%nexclsites
          params%exclsitetype(i) = toint(fields(i))
        End Do
      End If
    Else
      params%nexclsites = 0
    End If

  End Subroutine mcmoves_initaux

  !-------------------------------------------------------------------
  ! Some Monte Carlo moves can be stopped while attepmting itself, if we 
  ! know that energy has exceeded a particular limit, especially 
  ! cbgcmc Type of moves. That value (the threshold) is set from here 
  ! Requires:  mcmove --  mc move params
  !            max_nrg --  value of maximum energy (kcal/mol)
  !-------------------------------------------------------------------
  Subroutine mcmoves_setMaxNrg(mcmove, max_nrg)
    Type(MC_Move_Params), Intent(InOut)            :: mcmove
    Real(kind=RDbl) , Intent(in)                   :: max_nrg

    Call moves_setMaxNrg(mcmove%mvparams,max_nrg)

  End Subroutine mcmoves_setMaxNrg

  !-------------------------------------------------------------------
  ! Uses the auxilliary criteria to determine if the move should be
  ! reject outright.  Currently used only for rejecting moves into
  ! excluded sites.
  ! Requires: params -- auxilliary rejection params
  !           spc -- species number
  !           molec -- molecule number
  !           species -- species data structure
  !-------------------------------------------------------------------
  Logical Function mcmoves_auxreject(params,subset,species)
    Type(AuxReject_Params), Pointer                :: params
    Integer, Dimension(3), Intent(In)                            :: subset 
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species

    Integer                :: s,a,sitetype, spc, molec,natoms
    Type(VecType)          :: com,avec

    mcmoves_auxreject = .False.

    If (.Not. Associated(params)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " ERROR: Parameters not associated"
      Stop
    End If
    If (subset(1)==0) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Else
      spc=subset(1)
    Endif
    If (subset(2)==0) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Else
      molec=subset(2)
    Endif

    If (Associated(params%smap)) Then
      If (Associated(species(spc)%gcoords(molec)%branchedcoords)) Then
        ! check each atoms
        natoms=molecules_getnatoms(spc)
        Do a=1,natoms
          avec = species(spc)%coords(a,molec)%r
          sitetype = smap_getPosSiteType(params%smap,avec)

          !** Reject move if the molecule is in an excluded site
          Do s = 1,params%nexclsites
            If (sitetype == params%exclsitetype(s)) Then
              mcmoves_auxreject = .True.          
              Return
            End If
          End Do
        End Do
      Else
        com = gcmodels_getcom(species(spc)%gcoords(molec))
        sitetype = smap_getPosSiteType(params%smap,com)

        !** Reject move if the molecule is in an excluded site
        Do s = 1,params%nexclsites
          If (sitetype == params%exclsitetype(s)) Then
            mcmoves_auxreject = .True.          
            Return
          End If
        End Do
      Endif


    End If

  End Function mcmoves_auxreject

  !-----------------------------------------------------------------------
  ! Displays the contents of the auxilliary rejection data type, displays
  ! nothing if the critical fields are not allocated
  ! Requires: params -- auxilliary rejection data type
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !-----------------------------------------------------------------------
  Subroutine mcmoves_auxrejectdisplay(params,indent,unit)
    Type(AuxReject_Params), Intent(In)  :: params
    Integer, Intent(In)                 :: indent,unit

    Character(len=indent)       :: blank
    Character(len=strLen)       :: display

    blank = Repeat(' ',indent)

    If (Associated(params%exclsitetype)) Then
      Write(unit,'(2a)') blank,"Auxilliary move rejection information:"
      Write(unit,'(2a,i2)') blank," number of excluded sites: ", &
          params%nexclsites
      Write(display,'(a,i1,a)') '(2a,',Size(params%exclsitetype,1),'i3)'
      Write(unit,display) blank," excluded sites: ", &
          params%exclsitetype(1:Size(params%exclsitetype,1))
      Write(unit,'(3a)') blank," Sitemap filename: ", &
          Trim(params%smapname)
    End If

  End Subroutine mcmoves_auxrejectdisplay

  !-----------------------------------------------------------------------
  ! Displays the information in the "mcmove" data structure
  ! Requires: mcmove -- general move type parameters
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !-----------------------------------------------------------------------
  Subroutine mcmoves_display(mcmove,indent,unit)
    Type(MC_Move_Params), Intent(In)    :: mcmove
    Integer, Intent(In)                 :: indent,unit

    Character(len=indent)       :: blank
    Character(len=strLen)       :: display,molec_name

    blank = Repeat(' ',indent)

    Write(unit,'(3a)') blank,"Basic Move type:    ",Trim(mcmove%basic_tag)
    Write(unit,'(3a)') blank,"Specific Move type: ",Trim(mcmove%mvTag)

    Call moves_displayparams(mcmove%mvparams,indent+2,unit)
    Call mcmoves_auxrejectdisplay(mcmove%auxparams,indent+2,unit)

    Write(unit,'(2a)') blank,Trim(dashedline)

  End Subroutine mcmoves_display

  !-----------------------------------------------------------------------
  ! Single line display information giving the move statistics 
  ! Requires: mcmove -- general move type parameters
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !-----------------------------------------------------------------------
  Subroutine mcmoves_stats(mcmove,indent,unit)
    Type(MC_Move_Params), Intent(In)    :: mcmove
    Integer, Intent(In)                 :: indent,unit

    Character(len=indent)       :: blank
    Character(len=lstrLen)      :: display,string1,string2

    blank = Repeat(' ',indent)

    string1 = int2str(indent+27)
    Write(display,'(3a)') "(5a,t",Trim(string1),",a,2x,a)"

    string1 = moves_movestats(mcmove%stats)
    string2 = moves_dyninfo(mcmove%mvparams) 
    Write(unit,display) blank,Trim(mcmove%basic_tag)," (", &
        Trim(mcmove%mvTag),") ",Trim(string1),Trim(string2)

  End Subroutine mcmoves_stats

  !---------------------------------------------------------------------------
  ! Checks whether opposite moves like insert-delete are properly initialized
  ! Also copies some relevant fields from insert to rotate/translate/integrate
  ! The need for this : While performing a ROTATE move, you also might need
  ! to know what cavity object is being used by INSERT. This can be achieved 
  ! by repeating cavitybias information in control file in the rotate
  ! section also.  But we want to do it with minimal controlfile
  ! changes, and without increasing the total number of movetypes. 
  ! So this subroutine seems to be the best solution.
  ! Requires: mcmove -- general move type parameters
  !           species -- coordinates
  !---------------------------------------------------------------------------
  Subroutine mcmoves_checkOppMoves(params, species, auxmv, ctrlfilename)
    Type(MC_Move_Params), Dimension(:), Intent(Inout)    :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)       :: species
    Type(AuxMoveObjects),Pointer           :: auxmv 
    Character(len=strLen), Intent(In)                    :: ctrlfilename

    Integer           :: i, nmoves, ins_n, del_n, trans_n, rot_n
    Integer           :: integ_n, xform_n
    Integer           :: insnum, delnum, transnum, rotnum, integnum, xformnum
    Integer           :: icnum, ic_n

    !** Get nmoves and zero the counters for each move type
    nmoves=Size(params,1)
    insnum=0
    delnum=0
    transnum=0
    rotnum=0
    integnum=0
    xformnum=0

    !** Increment the counters based on the basic tag
    Do i = 1,nmoves
      Select Case(params(i)%basic_tag)
      Case("INSERT")
        insnum=insnum+1
        ins_n=i
      Case("DELETE")
        delnum=delnum+1
        del_n=i
      Case("TRANSLATE")      
        transnum=transnum+1
        trans_n=i
      Case("ROTATE")      
        rotnum=rotnum+1
        rot_n=i
      Case("TRANSFORM")      
        xformnum=xformnum+1
        xform_n=i
      Case("INTEGRATE")
        integnum=integnum+1
        integ_n=i
      Case("IDCHANGE")
        icnum=icnum+1
        ic_n=i
      Case Default
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            " Could not interpert basic tag: ",Trim(params(i)%basic_tag)
        Stop
      End Select
    End Do

    !** Stop with error message if there is more than one insert or delete
    If ((insnum > 1).Or.(delnum > 1)) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(1x,a)') 'More than one INSERT or DELETE move ', &
          'detected in control file'
      Stop
    Endif

    !** Note again ; may change position in currently open ctrlfile
    !** make sure cavity bias parameters are equal for all movetypes
    If ((insnum+delnum)==2) Then
      Call moves_checkOppMoves(params(ins_n)%mvparams, &
          params(del_n)%mvparams, species, auxmv, ctrlfilename)
    End If

  End Subroutine mcmoves_checkOppMoves

  !-----------------------------------------------------------------------
  ! Updates things like cavity bias objects 
  ! Requires:  params -- general move type parameters
  !            species -- the coordinates
  !            spc, molec  -- on which move was performed
  !            success -- whether accepted or not
  !-----------------------------------------------------------------------
  Subroutine mcmoves_postadjust(params,species,spc,molec,success)
    Type(MC_Move_Params), Intent(Inout)             :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Integer, Intent(In)                             :: spc, molec
    Logical, Intent(In)                             :: success

    If (success) Then
      Call moves_postadjust(params%mvparams, species, spc, molec )
    Else
      Return
    End If

  End Subroutine mcmoves_postadjust

  !-----------------------------------------------------------------------
  ! Reset all counter to zero 
  ! Requires: mcmove -- general move type parameters
  !-----------------------------------------------------------------------
  Subroutine mcmoves_resetStats(mcmove)
    Type(MC_Move_Params), Intent(InOut)    :: mcmove

    mcmove%stats%att = 0.0_RDbl
    mcmove%stats%succ = 0.0_RDbl
    mcmove%stats%nrg_succ = 0.0_RDbl

  End Subroutine mcmoves_resetStats

End Module mcmoves


