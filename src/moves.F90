!----------------------------------------------------------------------------
! This is a "proxy module" acting as a common interface between all possible
! movetypes.  The movetypes interfaced by this module all perform a change
! in the system state and return any quantities, including new energies that
! are necessary to calculate the acceptance of this new system configruation.
! The acceptance or rejection of these new states is not handled by this 
! module, but instead at a higher level.  Most of the moves operate only
! on a single molecules.  Those that handle more, are called "common moves"
!
! The interfaces group the specific move types by basic move type.  For
! example, "moves_translate" feeds to the fixed or random translate routines
!
! The "basic" move types are:
!   translate -- rigid translation of a molecule
!   rotate -- rigid rotation of a molecule
!   insert -- insertion of a molecule into system
!   delete -- deletion of a molecule from system
!   transform -- coordinate perturbation that changes internal shape
!   idchange -- change the identity of a molecule
!
! Currently available fundamental move types are:
!   fixed translation - move the molecule a specified distance
!   random translation - move the molecule a bounded random distance
!   fixed rotation - rotate the molecule by a specified increment
!   random rotation - rotate the molecule by bounded random increment
!   reinsertion - reinserts an existing molecule with random orientation
!   random insertion - insertion of randomly generated molecule
!   biased insertion - energy-biased insertion of random molecule
!   configurational biased insertion - configurational biased insertion
!   library insertion - insertion from library
!   random deletion - deletion of molecule
!   biased deletion - energy-biased deletion of molecule
!   configurational biased deletion - configurational biased deletion
!   cut and regrow - deletion of molecule part with biased regrowth
!   eigenmode perturb - perturbation along molecular non-rot/trans eigenmodes
!   integration - molecular dynamics integration of molecule in time
!
! Needed Improvements:
! 1) finish the "requires" comments
! 2) consider grouping all translate, rotate, insert, delete or transform
!    pointers into distinct types which would be stored in these routines.
!    This would make things a little less crowded and messy here.  I just
!    did this for the transform moves, use as example.
!----------------------------------------------------------------------------

Module moves
  Use auxmoveparams, Only: AuxMoveObjects
  Use brmoves, Only: brmoves_checkOppMoves
  Use config, Only: AtMolCoords
  Use datafile, Only: CONFILE
  Use defaults, Only: RDbl,strLen,lstrLen,dashedline,zero,one,TOTAL_INDEX
  Use delete, Only: Delete_Params, BDelete_Params, CBDelete_Params, &
      delete_init, delete_move, delete_display, &
      delete_setmaxnrg, delete_idstrings, delete_postadjust, delete_checkIntra
  Use file, Only: file_getunit
  Use molecules, Only: molecules_getgcmodeltype
  Use simcell, Only: SimCell_Params
  Use utils, Only: isfileopen,filesrchstr,stripcmnt,toupper,allocErrDisplay, &
      int2str,real2str,findstr
  Use vector, Only: VecType,Assignment(=),Operator(+),Operator(-), &
      Operator(*),Operator(/)
  Use idchange, Only: IDchange_Moves,idchange_init,idchange_move, &
      idchange_idstrings,idchange_restore,idchange_display,idchange_clean, &
      idchange_spcs, idchange_dyninfo, idchange_postadjust
  Use insert, Only: Insert_Params, BInsert_Params, CBInsert_Params, &
      insert_init, insert_display, &
      insert_setmaxnrg, insert_idstrings, insert_postadjust, insert_move, &
       insert_checkIntra
  Use integrate, Only :Integrate_Params,integrate_init,integrate_integrate,&
      integrate_displayparams, integrate_sampleCF, integrate_simdisplayExtra, &
      integrate_getdt, integrate_integrateMC
  Use interact, Only: Interaction_Model
  Use subinteract, Only: Subset_Interactions
  Use rotate, Only: Rotate_Params, rotate_init, rotate_adjustdelta, &
      rotate_display, rotate_dyninfo, rotate_idstrings, &
      rotate_postadjust, rotate_move, rotate_setmaxnrg
  Use translate, Only: Translate_Params, translate_init, translate_display, &
      translate_dyninfo, translate_adjustdelta, translate_idstrings, &
      translate_postadjust, translate_move, translate_setmaxnrg
  Use transform, Only: Transform_Moves,transform_init,transform_move, &
      transform_setmaxnrg,transform_idstrings,transform_adjustdelta, &
      transform_dyninfo,transform_display,transform_clean, &
      transform_postadjust
  Use rigidmoves, Only: rigidmoves_postadjust, rigidmoves_checkoppmoves
  
  Implicit None
  Save

  Private 
  Public :: Move_Params,Move_stats,moves_initparams,moves_makemove, &
      moves_integrate,moves_adjustparams,moves_gettag,moves_sampleCF, &
      moves_setmaxnrg,moves_movestats,moves_movestatsdisplay, &
      moves_dyninfo,moves_simdisplayExtra,moves_getparam,moves_restore, &
      moves_checkOppMoves, moves_postadjust, &
      moves_displayparams,moves_clean, moves_checkIntra

  !** att=number of attempts,succ=number of successes 
  !** nrg_succ=number of nrg calculation successe
  Type Move_Stats
    Integer       :: att
    Integer       :: succ,nrg_succ
  End Type Move_Stats
  
  Type Move_Params
    Character(len=strLen)              :: tag        ! specific move-name
    Character(len=strLen)              :: basic_tag  ! general name
    Character(len=strLen)              :: modeltype

    Type(Transform_Moves),  Pointer    :: transform    !general transformation
    Type(Integrate_Params),  Pointer   :: integrate    !Integration

    Type(Translate_Params), Pointer    :: translate    !general translation
    Type(Rotate_Params), Pointer       :: rotate       !general rotattion 

    Type(Insert_Params), Pointer       :: insert       !general insertion
    Type(Delete_Params),  Pointer      :: delete       !general deletion

    Type(IDchange_Moves), Pointer      :: idchange     !general ID change
  End Type Move_params

  Interface moves_initparams
    Module Procedure moves_initparams
    Module Procedure moves_initCommonMoves
  End Interface

Contains
  !---------------------------------------------------------------------
  ! Initializes "mvparams" for the given "sorbtype".  This routine is
  ! called when a move type section of the control file is encountered.
  ! It forms a link to the individual move type module initializations.  
  ! It also labels each move type with a basic tag, which is then used
  ! by the calling routines to decide how the moves should be called.
  ! Requires: mvparams -- parameter set to be initialized
  !           simcell -- the simulation cell information
  !           sorbtype -- the species type
  !           filename -- the filename to get the information from
  !           basic_tag -- the basic move type
  !           sorbates -- the simulation species
  !---------------------------------------------------------------------
  Subroutine moves_initparams(mvparams, auxmv, simcell, sorbtype, & 
      filename, basic_tag, sorbates) 
    Type(Move_Params), Intent(inout)          :: mvparams
    Type(AuxMoveObjects),Pointer           :: auxmv 
    Type(SimCell_Params), Intent(in)          :: simcell
    Integer, Intent(in)                       :: sorbtype
    Character(*), Intent(in)                  :: filename
    Character(len=strLen), Intent(Out)        :: basic_tag
    Type(AtMolCoords), Dimension(:), Optional, Intent(In) :: sorbates 

    Integer                   :: error,unitno,index
    Character(len=strLen)     :: modeltype,move_type,line
    Logical :: init_success

    init_success=.False.

    !** Nullify all the pointer instances in mvparams
    Call moves_nullifyAllFields(mvparams)

    !** get the gcmodel type
    modeltype = molecules_getgcmodeltype(sorbtype)
    mvparams%modeltype=Trim(toupper(modeltype))

    unitno = file_getunit(filename)

    Read(unitno,'(a)') line
    move_type=toupper(Trim(stripcmnt(line)))

    !** Check to see if it's a transform move and initialize
    index = findstr(transform_idstrings,move_type)
    If (index /= 0) Then
      Allocate(mvparams%transform, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call transform_init(mvparams%transform,move_type,modeltype, &
          filename,sorbtype,sorbates,simcell,auxmv)
      basic_tag = 'TRANSFORM' 
      mvparams%tag = transform_idstrings(index)
      mvparams%basic_tag = Trim(basic_tag)
      init_success=.True.
      Return
    End If

    !** Check to see if it's an insert move and initialize
    index = findstr(insert_idstrings,move_type)
    If (index /= 0) Then
      Allocate(mvparams%insert, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call insert_init(mvparams%insert,move_type,modeltype, &
          filename,sorbtype,sorbates,simcell,auxmv)
      basic_tag = 'INSERT' 
      mvparams%tag = insert_idstrings(index)
      mvparams%basic_tag = Trim(basic_tag)
      init_success=.True.
      Return
    End If

    !** Check to see if it's an delete move and initialize
    index = findstr(delete_idstrings,move_type)
    If (index /= 0) Then
      Allocate(mvparams%delete, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call delete_init(mvparams%delete, move_type,modeltype, &
          filename,sorbtype,sorbates,simcell,auxmv)
      basic_tag = 'DELETE' 
      mvparams%tag = delete_idstrings(index)
      mvparams%basic_tag = Trim(basic_tag)
      init_success=.True.
      Return
    End If

    !** Check to see if it's an translate move and initialize
    index = findstr(translate_idstrings,move_type)
    If (index /= 0) Then
      Allocate(mvparams%translate, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call translate_init(mvparams%translate, move_type,modeltype, &
          filename,sorbtype,sorbates,simcell,auxmv)
      basic_tag = 'TRANSLATE' 
      mvparams%tag = translate_idstrings(index)
      mvparams%basic_tag = Trim(basic_tag)
      init_success=.True.
      Return
    End If

    !** Check to see if it's an rotate move and initialize
    index = findstr(rotate_idstrings,move_type)
    If (index /= 0) Then
      Allocate(mvparams%rotate, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call rotate_init(mvparams%rotate, move_type,modeltype, &
          filename,sorbtype,sorbates,simcell,auxmv)
      basic_tag = 'ROTATE' 
      mvparams%tag = rotate_idstrings(index)
      mvparams%basic_tag = Trim(basic_tag)
      init_success=.True.
      Return
    End If

    !** Check to see if it's an identity change move and initialize
    index = findstr(idchange_idstrings,move_type)
    If (index /= 0) Then
      Allocate(mvparams%idchange, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call idchange_init(mvparams%idchange,move_type,&
          filename,sorbtype,sorbates,simcell,auxmv)
      basic_tag = 'IDCHANGE' 
      mvparams%tag = idchange_idstrings(index)
      mvparams%basic_tag = Trim(basic_tag)
      init_success = .True.
      Return
    End If

    If (.Not.init_success) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' Could not get parameters for move_type: "', Trim(move_type),'"'
      Stop
    Endif

  End Subroutine moves_initparams

  !--------------------------------------------------------------
  ! Initializes "mvparams" for the whole system moves 
  !--------------------------------------------------------------
  Subroutine moves_initCommonMoves(mvparams, auxmv, sorbates,filename, opt_tag)
    Type(Move_Params), Intent(inout)            :: mvparams
    Type(AuxMoveObjects),Intent(inout)           :: auxmv 
    Type(AtMolCoords), Intent(In), Dimension(:) :: sorbates
    Character(*), Intent(in)                    :: filename
    Character(*), Intent(out), Optional         :: opt_tag 

    Integer                   :: error,unitno
    Character(len=strLen)     :: move_type,line,basic_tag
    
    !** Nullify all the pointer instances in mvparams
    Call moves_nullifyAllFields(mvparams)
    unitno=file_getunit(filename)
    Read(unitno,'(a)') line

    move_type=toupper(Trim(stripcmnt(line)))
    
    !** Call the appropriate move-init routine based on the move_type
    Select Case (Trim(toupper(move_type)))
    Case('INTEGRATE')
      Allocate(mvparams%integrate, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Call integrate_init(mvparams%integrate,sorbates,filename)
      basic_tag = 'INTEGRATE'
      mvparams%tag = 'INTEGRATE'
      If (Present(opt_tag)) opt_tag=basic_tag
       
    Case Default
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not get parameters for move_type ", Trim(move_type)
      Stop
         
    End Select

  End Subroutine moves_initCommonMoves

  !----------------------------------------------------------------------
  ! Actually perform the move, pick the correct module and pass the 
  ! move parameters and system state to it.
  ! Requires:  params -- generalized move parameters
  !            subset -- system subset to which to apply move
  !            subints -- subset interactions for each species
  !            species -- the species coordinates
  !            simcell -- simulation cell information
  !            biasfactor -- the returned bias factor
  !----------------------------------------------------------------------
  Logical Function moves_makemove(params,subset,subints,species,simcell, &
      biasfactor)
    Type(Move_Params), Intent(InOut)                       :: params
    Integer, Dimension(3), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Intent(Out)                           :: biasfactor

    Logical         :: succ_flag
    
    !** Identify the basic move type and call appropriate module
    If (Associated(params%transform)) Then
      succ_flag = transform_move(params%transform,subset,species, &
          simcell,biasfactor,subints(subset(1)))

    Else If (Associated(params%translate)) Then
      succ_flag = translate_move(params%translate,  subset, species, &
          simcell,biasfactor,subints(subset(1)))

    Else If (Associated(params%rotate)) Then
      succ_flag = rotate_move(params%rotate, subset, species, simcell, &
          biasfactor,subints(subset(1)))
      biasfactor = one

    Else If (Associated(params%insert)) Then
      succ_flag = insert_move(params%insert, subset, species, simcell, &
          biasfactor,subints(subset(1)))

    Else If (Associated(params%delete)) Then
      succ_flag = delete_move(params%delete, subset, species, simcell, &
          biasfactor,subints(subset(1)))

    Else If ((Associated(params%integrate))) Then
      Call integrate_integrateMC(params%integrate, subints(1), &
          species,simcell,succ_flag)
      biasfactor = one

    Else If (Associated(params%idchange)) Then
      succ_flag = idchange_move(params%idchange,subset,subints,species, &
          simcell,biasfactor)

    Else
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Could not identify associated move pointer"
      Stop
    End If

    !** set the success status of the move
    moves_makemove = succ_flag

  End Function moves_makemove

  !----------------------------------------------------------------------
  ! Reverse a move.  Can be used to restore the state of a system after
  ! an unsuccesful move is made.
  ! Requires:  params -- generalized move parameters
  !            subset -- system subset to which to apply move
  !            subints -- subset interactions for each species
  !            species -- the species coordinates
  !            simcell -- simulation cell information
  !----------------------------------------------------------------------
  Subroutine moves_restore(params,subset,subints,species,simcell)
    Type(Move_Params), Intent(InOut)                       :: params
    Integer, Dimension(3), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell

    Logical           :: succ_flag

    !** Identify the basic move type and call appropriate module
    If (Associated(params%idchange)) Then
      Call idchange_restore(params%idchange,subset,subints,species,simcell)
          
    Else
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Move type not supported in moves_restore "
      Stop
    End If

  End Subroutine moves_restore

  !----------------------------------------------------------------------
  ! Calls the integrate routine in the integrate module
  !----------------------------------------------------------------------
  Subroutine moves_integrate(subint, mvparams, simcell, sorbates, conf_file, &
      updateFlag, time)
    Type(Subset_Interactions), Intent(InOut)          :: subint
    Type(Move_Params), Intent(InOut)                :: mvparams
    Type(SimCell_Params), Intent(In)                :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: sorbates
    Type (CONFILE),intent(in)                       :: conf_file          
    Real(kind=RDbl),Intent(inout)                   :: time
    Logical,Intent(in)                              :: updateFlag

    Call integrate_integrate(subint, mvparams%integrate, simcell, &
        sorbates, conf_file, updateFlag, time)

  End Subroutine moves_integrate

  !--------------------------------------------------------------
  ! Gets the move tag, useful for identifying move type
  !--------------------------------------------------------------
  Function moves_gettag(mvparams)
    Character(len=strLen)               :: moves_gettag
    Type(Move_Params), Intent(In)       :: mvparams

    moves_gettag = mvparams%tag
  End Function moves_gettag

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file to the unit unitno
  !----------------------------------------------------------------------------
  Subroutine moves_sampleCF(unitno,type)
    Integer, Intent(In) :: unitno
    Character(*), Intent(in) :: type

    Integer        :: index

    Select Case(toupper(type))
    Case('TRANSLATE')
      Write(unitno,'(a,t30,a)') 'TRANSLATE','# Translate move type'
      Write(unitno,'(a)') 'Lazy Programmer Fault: Routine does not exist'
      
    Case('ROTATE')
      Write(unitno,'(a,t30,a)') 'ROTATE','# Rotate move type'
      Write(unitno,'(a)') 'Lazy Programmer Fault: Routine does not exist'

    Case('INSERT')
      Write(unitno,'(a,t30,a)') 'INSERT','# Insert move type'
      Write(unitno,'(a)') 'Lazy Programmer Fault: Routine does not exist'

    Case('DELETE')
      Write(unitno,'(a,t30,a)') 'DELETE','# Delete move type'
      Write(unitno,'(a)') 'Lazy Programmer Fault: Routine does not exist'

    Case ('INTEGRATE')
      Write(unitno,'(a,t30,a)') 'INTEGRATE','# Integrate move type'
      Call integrate_sampleCF(unitno)

    Case('IDCHANGE')
      Write(unitno,'(a,t30,a)') 'IDCHANGE','# Identity change move type'
      Write(unitno,'(a)') 'Lazy Programmer Fault: Routine does not exist'
      
    Case Default
      !** check to see if it's a transform move
      index = findstr(transform_idstrings,type)
      If (index /= 0) Then
        Write(unitno,'(a,t30,a)') 'PERTURB','# Perturb move type'
        Write(unitno,'(a)') 'Lazy Programmer Fault: Routine does not exist'
        Stop
      Else
        Write(unitno,'(2a)') 'Unable to identify move tag: ',Trim(Type)
      End If

    End Select

  End Subroutine moves_sampleCF

  !----------------------------------------------------------------------
  ! For setting the maximum allowed energy during a move.
  ! Useful for cbgcmc
  ! Requires:  mvparams -- generalized move parameters
  !            max_nrg -- new maximum energy
  !----------------------------------------------------------------------
  Subroutine moves_setMaxNrg(mvparams,max_nrg)
    Type(Move_Params), Intent(InOut)            :: mvparams
    Real(kind=RDbl), Intent(In)                 :: max_nrg

    If (Associated(mvparams%transform)) Then
      Call transform_setmaxnrg(mvparams%transform,max_nrg)

    Else If (Associated(mvparams%delete)) Then
      Call delete_setmaxnrg(mvparams%delete, max_nrg)

    Else If (Associated(mvparams%insert)) Then
      Call insert_setmaxnrg(mvparams%insert,max_nrg)

    Else If (Associated(mvparams%translate)) Then
      Call translate_setmaxnrg(mvparams%translate,max_nrg)

    Else If (Associated(mvparams%rotate)) Then
      Call rotate_setmaxnrg(mvparams%rotate,max_nrg)

    Else
      Write(0,*) "This movetype does not require max_nrg"
    Endif

  End Subroutine moves_setMaxNrg

  !----------------------------------------------------------------------
  ! For setting the maximum allowed energy during a move.
  ! Useful for cbgcmc
  ! Requires:  mvparams -- generalized move parameters
  !            max_nrg -- new maximum energy
  !----------------------------------------------------------------------
  Subroutine moves_checkIntra(mvparams,use_intra)
    Type(Move_Params), Intent(In)            :: mvparams
    Logical , Intent(InOut)                 :: use_intra

    If (Associated(mvparams%transform)) Then
      use_intra=.True.
    Else If (Associated(mvparams%delete)) Then
      Call delete_checkIntra(mvparams%delete, use_intra)
    Else If (Associated(mvparams%insert)) Then
      Call insert_checkIntra(mvparams%insert, use_intra)
    Else If (Associated(mvparams%translate)) Then
       use_intra=.False.
    Else If (Associated(mvparams%rotate)) Then
      use_intra=.False.
    Else If (Associated(mvparams%integrate)) Then
      use_intra=.True.
    Else If (Associated(mvparams%idchange)) Then
      ! may not be necessary for some idchange types 
      use_intra=.True.
    Else
      Write(0,*) "This movetype does not require check_intra?"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__      
      Stop
    Endif

  End Subroutine moves_checkIntra


  !----------------------------------------------------------------------
  ! Nullifies all fields
  !----------------------------------------------------------------------
  Subroutine moves_nullifyAllFields(mvparams)
    Type(Move_Params), Intent(inout)                    :: mvparams

    Nullify(mvparams%transform)
    Nullify(mvparams%translate)
    Nullify(mvparams%rotate)
    Nullify(mvparams%insert)
    Nullify(mvparams%delete)
    Nullify(mvparams%integrate)
    Nullify(mvparams%idchange)

  End Subroutine moves_nullifyAllFields

  !----------------------------------------------------------------------
  ! Calls the various adjust routines
  !----------------------------------------------------------------------
  Subroutine moves_adjustParams(mvparams,stats)
    Type(Move_Params), Intent(InOut)        :: mvparams
    Type(Move_Stats), Intent(In)            :: stats

    Real(kind=RDbl)                 :: ratio

    ratio = stats%succ/(stats%att*one)

    If (Associated(mvparams%translate)) Then
      Call translate_adjustdelta(mvparams%translate, ratio)

    Else If (Associated(mvparams%rotate)) Then
      Call rotate_adjustdelta(mvparams%rotate,ratio)

    Else If (Associated(mvparams%transform)) Then
      Call transform_adjustdelta(mvparams%transform,ratio)

    Else
      !** this routine might be called redundantly, that's OK.
      !   Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__," Nothing Associated"
      !   Stop
    End If

  End Subroutine moves_adjustParams

  !----------------------------------------------------------------------------
  ! This is a general subroutine that may be used to return information about
  ! the various move types, such as the time step for MD or the temperature
  ! for NVTMC. Supplied with one of a predetermined list of keywords, it will
  ! return the proper value.  When returning species numbers, the optional
  ! integer array must be present.
  ! Requires:  mvparams -- generalized move parameters
  !            key -- string keyword
  !            val -- returned real value
  !            intarray -- optional returned integer array
  !----------------------------------------------------------------------------
  Subroutine moves_getparam(mvparams,key,val,intarray)
    Type(Move_params), Intent(In)       :: mvparams
    Character(*), Intent(In)            :: key
    Real(Kind=RDbl), Intent(Out)        :: val
    Integer, Dimension(:), Intent(Out)  :: intarray

    Logical          :: nomatch

    !** True if the keyword was recognized, but no appropriate pointer associated
    nomatch = .False.

    !** Check to see which key to use
    Select Case (Trim(toupper(key)))
    Case ("TIMESTEP")   !** Returns the timestep
      If (Associated(mvparams%integrate)) Then
        val = integrate_getdt(mvparams%integrate)
        Return
      Else
        nomatch = .True.
      End If

    Case ("SPCS")     !** Returns the species number(s) for ID change moves
      If (Associated(mvparams%idchange)) Then
        Call idchange_spcs(mvparams%idchange,intarray)
        Return
      Else
        nomatch = .True.
      End If

    Case Default
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' Unable to recognized move parameter keyword: ',Trim(key)
      Stop
    End Select

    If (nomatch) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' Matched keyword, but no appropriate pointer initialized'
      Stop
    End If

  End Subroutine moves_getparam

  !--------------------------------------------------------------
  ! Display the information in the passed move parameters
  ! Requires: mvparams -- general move type parameters
  !           unit -- unit to write into
  !           indent -- number of spaces to indent
  !--------------------------------------------------------------
  Subroutine moves_displayparams(mvparams,indent,unitno)
    Type(Move_Params), Intent(In)  :: mvparams
    Integer, Intent(In)            :: indent,unitno

    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    If (Associated(mvparams%transform)) &
      Call transform_display(mvparams%transform,indent,unitno)
    If (Associated(mvparams%translate)) &
      Call translate_display(mvparams%translate,indent,unitno)
    If (Associated(mvparams%rotate)) &
      Call rotate_display(mvparams%rotate,indent,unitno)
    If (Associated(mvparams%insert)) &
      Call insert_display(mvparams%insert,indent,unitno)
    If (Associated(mvparams%delete)) &
      Call delete_display(mvparams%delete,indent,unitno)
    If (Associated(mvparams%idchange))   &
        Call idchange_display(mvparams%idchange,indent,unitno)
    If (Associated(mvparams%integrate))   &
        Call integrate_displayparams(mvparams%integrate,unitno)
    
  End Subroutine moves_displayparams

  !----------------------------------------------------------------------------
  ! Displays extra information that the calling routines may not know about
  !----------------------------------------------------------------------------
  Subroutine moves_simdisplayExtra(params,imodel,unitno,indent)
    Type(Move_Params), Intent(In)              :: params
    Type(Interaction_Model), Intent(InOut)     :: imodel
    Integer, Intent(In)                        :: unitno, indent

    If (Associated(params%integrate)) Then
      Call integrate_simdisplayExtra(params%integrate,imodel,unitno,indent)
    Else
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          ' No recognized pointers associated'
      Stop
    End If

  End Subroutine moves_simdisplayExtra

  !--------------------------------------------------------------
  ! Display the information in the passed move statistics type
  ! Requires: mvstats -- move statistics type
  !           unit -- unit to write into
  !           indent -- number of spaces to indent
  !--------------------------------------------------------------
  Subroutine moves_movestatsdisplay(mvstats,unitno,indent)
    Type(Move_Stats), Intent(In)                 :: mvstats
    Integer, Intent(In)                          :: unitno,indent

    Real(kind=RDbl)             :: ratio
    Character(len=indent)       :: blank
    Character(len=strLen)       :: display

    blank = Repeat(' ',indent)

    If (mvstats%att > 0) Then
      ratio = (mvstats%succ*1.0_RDbl)/mvstats%att
    Else
      ratio = 0.0_RDbl
    End If

    display = 'f7.5'
    If (ratio < 1.0E-4_RDbl) display = 'e12.6'

    Write(unitno,'(a,2(a,i10,4x))') blank,"Attempts: ",mvstats%att, &
        "Successes: ",mvstats%succ
    Write(unitno,'(2a,i10,5x,a)') blank,"Energy calc. successes: ", &
        mvstats%nrg_succ,' broken?'
    display = "(2a,"//Trim(Adjustl(display))//")"
    Write(unitno,display) blank,"Acceptance Ratio: ",ratio

  End Subroutine moves_movestatsdisplay

  !----------------------------------------------------------------------
  ! Single line display information for the passed move statistics type
  ! Requires: mvstats -- move statistics type
  !----------------------------------------------------------------------
  Function moves_movestats(mvstats)
    Character(len=lstrLen)          :: moves_movestats
    Type(Move_Stats), Intent(In)    :: mvstats

    Real(kind=RDbl)             :: ratio
    Character(len=strLen)       :: string1,string2,string3

    If (mvstats%att/=0) Then
      ratio = (mvstats%succ*1.0_RDbl)/mvstats%att
    Else
      ratio = zero
    End If

    string1 = real2str(ratio,5)
    string2 = int2str(mvstats%succ)
    string3 = int2str(mvstats%att)
    Write(moves_movestats,'(2a,3x,5a)') "Acc. Ratio: ", &
        Trim(string1),"(",Trim(string2)," of ",Trim(string3),")"

  End Function moves_movestats

  !----------------------------------------------------------------------
  ! Single line display information for dynamic move parameters
  ! Requires: mvparams -- move parameters type
  !----------------------------------------------------------------------
  Function moves_dyninfo(mvparams)
    Character(len=lstrLen)          :: moves_dyninfo
    Type(Move_Params), Intent(In)   :: mvparams

    moves_dyninfo = ''

    If (Associated(mvparams%translate)) Then
      moves_dyninfo = translate_dyninfo(mvparams%translate)

    Else If (Associated(mvparams%rotate)) Then
      moves_dyninfo = rotate_dyninfo(mvparams%rotate)

    Else If (Associated(mvparams%insert)) Then
!      moves_dyninfo = insert_dyninfo(mvparams%insert)

    Else If (Associated(mvparams%delete)) Then
!      moves_dyninfo = delete_dyninfo(mvparams%delete)

    Else If (Associated(mvparams%transform)) Then
      moves_dyninfo = transform_dyninfo(mvparams%transform)

    Else If (Associated(mvparams%integrate)) Then
!      moves_dyninfo = integrate_dyninfo(mvparams%integrate)

    Else If (Associated(mvparams%idchange)) Then
      moves_dyninfo = idchange_dyninfo(mvparams%idchange)

    Else
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          ' No recognized pointers associated for dynamic info display'
      Stop
    End If
        
  End Function moves_dyninfo

  !----------------------------------------------------------------------
  ! Check whether two move types are opposite (insert and delete), and 
  ! make sure their biases also are properly set.  They must use the 
  ! same bmap and cavitylist objects.
  ! Requires:  mv1 -- move type #1 parameters
  !            mv2 -- move type #2 parameters
  !            species -- configuration data structure
  !            ctrlfilename -- control filename
  ! Note : here we will have to bypass insert/delete and directly call 
  ! the rigidmoves or bmoves routines
  ! Hmmm.. : one way of overcoming this problem might be let 'insert' 
  ! "Use" 'delete' since its the opposite move of 'insert'
  !--------------------------------------------------------------
  Subroutine moves_checkOppMoves(mv1,mv2,species,auxmv, ctrlfilename)
    Type(Move_Params), Intent(In), Target          :: mv1, mv2
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(AuxMoveObjects), Pointer           :: auxmv 
    Character(len=strLen), Intent(In)              :: ctrlfilename

    Logical                     :: Fatal_error, Warning_flag
    Type(Move_Params), Pointer  :: ins,del

    !** Set defaults
    Fatal_error  = .True.
    Warning_flag = .False.    !** not currently used

    !** Identify the insertion and deletion moves
    If (Trim(mv1%basic_tag) == "INSERT") Then
      ins => mv1
      del => mv2
    Else If (Trim(mv2%basic_tag) == "INSERT") Then
      ins => mv2
      del => mv1
    Else
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(2x,a)') 'one of the moves must be INSERT'
      Stop
    End If

    !** Make sure the other one is a delete move
    If (del%basic_tag /= "DELETE") Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(2x,2a)') 'An INSERT and a non-DELETE type move were passed: ', &
          Trim(del%basic_tag)
      Stop      
    End If

    !** Do the checking based on the insertion move type
    Select Case(Trim(ins%tag))
    Case("RINSERT")
      If ((del%tag) == "RDELETE") Then
        Fatal_error = .False.
      End If

    Case("BINSERT")
      If ((del%tag) == "BDELETE") Then
        If (Associated(ins%insert%biased%rigid)) Then
          Call rigidmoves_checkOppMoves(ins%insert%biased%rigid, &
              del%delete%biased%rigid,species, auxmv, ctrlfilename)
          Fatal_error = .False.
        End If
      End If
    
    Case("LINSERT")
      If ((del%Tag) == "BDELETE") Then
        If (Associated(ins%insert%lib%rigid)) Then
          Write(*,*) Associated(ins%insert%lib%rigid), &
              Associated(del%delete%biased%rigid)
          Call rigidmoves_checkOppMoves(ins%insert%lib%rigid, &
              del%delete%biased%rigid,species, auxmv, ctrlfilename)
          Fatal_error = .False.
        End If
      End If
    
    Case("CBINSERT")
      If ((del%Tag) == "CBDELETE") Then
        If (Associated(ins%insert%cbiased%branched)) Then
          Call brmoves_checkOppMoves(ins%insert%cbiased%branched, &
              del%delete%cbiased%branched, ctrlfilename)
          Fatal_error = .False.
        End If
      End If
    End Select

    !** Give the 'fatal' error if necessary
    If (Fatal_error) Then 
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(2x,3a)') 'INSERT and DELETE moves are of different types', &
          Trim(ins%tag), Trim(del%tag)
      Stop      
    End If

    !** Just give the 'warning' if necessary (not used for now)
    If (Warning_Flag) Then 
      Write(*,*) " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
      Write(*,*) " XXXXXXXXXXX  WARNING        XXXXXXXXXXXXX"
      Write(*,*) " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
      Write(*,*) " Looks like wrong movetypes were passed here "
      Write(*,*) "types : "//Trim(mv1%tag)//" "//Trim(mv2%tag)
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    End If

  End Subroutine moves_checkOppMoves

  !----------------------------------------------------------------------
  ! Adjusts the cavity bias object after a move has been accepted
  !----------------------------------------------------------------------
  Subroutine moves_postadjust(params, sorbates, sorb, molec)
    Type(Move_Params), Intent(InOut)               :: params
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Integer , Intent(in)                           :: sorb, molec

    Select Case(Trim(params%basic_tag))
    Case ("ROTATE")
      Call rotate_postadjust(params%rotate, sorbates, sorb, molec)
    Case ("TRANSLATE")
      Call translate_postadjust(params%translate, sorbates, sorb, molec)
    Case ("TRANSFORM")
      Call transform_postadjust(params%transform, sorbates, sorb, molec)
    Case ("DELETE")
      Call delete_postadjust(params%delete, sorbates, sorb, molec)
    Case ("INSERT")
      Call insert_postadjust(params%insert, sorbates, sorb, molec) 
    Case ("IDCHANGE")
      Call idchange_postadjust(params%idchange, sorbates)
    Case default
    End Select

  End Subroutine moves_postadjust

  !--------------------------------------------------------------
  ! Clean the move parameters
  ! Requires: mvparams -- general move type parameters
  !--------------------------------------------------------------
  Subroutine moves_clean(mvparams)
    Type(Move_Params), Intent(InOut)  :: mvparams

    If (Associated(mvparams%transform)) Then
      Call transform_clean(mvparams%transform)
    Else If (Associated(mvparams%idchange)) Then
      Call idchange_clean(mvparams%idchange)
    End If

    !** most of these routines don't exist and are therefore commented out
#if COMMENTEDOUT
    If (Associated(mvparams%FTransParams))  &
        Call translate_clean(mvparams%FTransParams)
    If (Associated(mvparams%RTransParams))  &
        Call translate_clean(mvparams%RTransParams)
    If (Associated(mvparams%FRotParams))  &
        Call rotate_clean(mvparams%FRotParams)
    If (Associated(mvparams%RRotParams))  &
        Call rotate_clean(mvparams%RRotParams)
    If (Associated(mvparams%BInsParams))  &
        Call insert_clean(mvparams%BInsParams)
    If (Associated(mvparams%CBInsParams))  &
        Call insert_clean(mvparams%CBInsParams)
    If (Associated(mvparams%LInsParams))  &
        Call insert_clean(mvparams%LInsParams)
    If (Associated(mvparams%RInsParams))  &
        Call insert_clean(mvparams%RInsParams)
    If (Associated(mvparams%BDelParams))   &
        Call delete_clean(mvparams%BDelParams)
    If (Associated(mvparams%CBDelParams))   &
        Call delete_clean(mvparams%CBDelParams)
    If (Associated(mvparams%RDelParams))   &
        Call delete_clean(mvparams%RDelParams)
    If (Associated(mvparams%IntegParams))   &     
        Call integrate_clean(mvparams%IntegParams)
#endif

  End Subroutine moves_clean

End Module moves
