!----------------------------------------------------------------------------
! This module contains the Move Set data structure and its routines.
! Contains the information necessary to pick a random move from a set 
! of moves.  Included are also routines for initializing the structure,
! reading the weights from a file line, initializing just the tags, getting
! tag information and picking a move randomly from the set.
!
! NOTE-1 : We allow for insert and delete to have different weights, 
!          This should be taken care of in the monte carlo acceptance 
!          criteria during insertions/deletions 
!----------------------------------------------------------------------------

Module mcmoveset

  Use defaults, Only: RDbl,strLen,zero,kcalmole_kb,one
  Use utils, Only: filesrchstr,toupper,split,stripcmnt,toreal,toint, &
      allocErrDisplay
  Use file, Only: file_open
  Use random, Only: rranf
       
  Implicit None
  Save

  Private
  Public :: Move_Set,mcmoveset_initms,mcmoveset_initmstags, &
      mcmoveset_setmstag,mcmoveset_readmswts,mcmoveset_getmoveno, &
      mcmoveset_getmstag,mcmoveset_normms,mcmoveset_pickmove, &
      mcmoveset_displayms,mcmoveset_getmswt, mcmoveset_setInsDelRatio, &
      mcmoveset_insdelratio

  Type Move_Set
    Integer                                      :: nmoves
    Real(kind=RDbl), Dimension(:), Pointer       :: wt,cumwt
    Character(len=strLen), Dimension(:), Pointer :: tag
    !** : number of insertions/ number of deletions, see NOTE-1 of header
    Real(kind=RDbl)                              :: ins_del_ratio
  End Type Move_Set

  Interface mcmoveset_getmswt
    Module Procedure mcmoveset_getmswt_number
    Module Procedure mcmoveset_getmswt_tag
  End Interface

Contains

  !---------------------------------------------------------------
  ! Initializes the Move Set structure
  ! Requires: moveset -- move set data structure
  !           nmoves -- number of moves
  !           weights -- array of weights (can be unnormalized)
  !           tags -- array of tags
  !---------------------------------------------------------------
  Subroutine mcmoveset_initms(moveset,nmoves,weights,tags)
    Type(Move_Set), Intent(InOut)               :: moveset
    Integer, Intent(In)                         :: nmoves
    Real(kind=RDbl), Dimension(:), Intent(In)   :: weights
    Character(len=strLen), Dimension(:), Intent(In), Optional :: tags

    Integer                 :: error
    Logical                 :: flag

    moveset%nmoves = nmoves

    flag = .False.
    If (Present(tags)) Then
      If (Size(tags,1) < nmoves) flag = .True.
    End If
    If (Size(weights,1) < nmoves) flag = .True.

    If (flag) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Passed array(s) are not large enough"
      Stop
    End If

    !** allocate space
    Allocate(moveset%wt(nmoves), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'moveset%wt')

    !** nullify so that moveset%tag can be allocated at a point when we come 
    !** across tags
    Nullify(moveset%tag)
    
    If (Present(tags)) Then
      Allocate(moveset%tag(nmoves), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'moveset%tag')
    End If

    !** make the assignments
    moveset%wt(1:nmoves) = weights(1:nmoves)
    If (Present(tags)) moveset%tag(1:nmoves) = tags(1:nmoves)

    !** normalize the weights and get cumulative weights
    Call mcmoveset_normms(moveset)

  End Subroutine mcmoveset_initms

  !---------------------------------------------------------------
  ! Initializes just the tags in the Move Set structure
  ! Requires: moveset -- move set data structure
  !           tags -- array of tags
  !---------------------------------------------------------------
  Subroutine mcmoveset_initmstags(moveset,tags)
    Type(Move_Set), Intent(InOut)                   :: moveset
    Character(len=strLen), Dimension(:), Intent(In) :: tags

    Integer                 :: error

    If (.Not. Associated(moveset%wt)) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Passed Move Set weights must already be allocated"
      Stop      
    End If

    If (Size(tags,1) < moveset%nmoves) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Passed array(s) are not large enough"
      Stop
    End If

    !** allocate space
    Allocate(moveset%tag(moveset%nmoves), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'moveset%wt')

    !** make the assignments
    moveset%tag(1:moveset%nmoves) = tags(1:moveset%nmoves)

  End Subroutine mcmoveset_initmstags

  !---------------------------------------------------------------
  ! Initializes just ONE tag in the Move Set structure
  ! Requires: moveset -- move set data structure
  !           moveno -- move number
  !           tag -- tag
  !---------------------------------------------------------------
  Subroutine mcmoveset_setmstag(moveset,moveno,tag)
    Type(Move_Set), Intent(InOut)        :: moveset
    Integer, Intent(In)                  :: moveno
    Character(len=strLen), Intent(In)    :: tag

    Integer                 :: error

    If (.Not. Associated(moveset%tag)) Then
      Allocate(moveset%tag(moveset%nmoves), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'moveset%wt')
    End If

    !** make the assignment
    moveset%tag(moveno) = tag

  End Subroutine mcmoveset_setmstag

  !---------------------------------------------------------------
  ! Gets the weights for each move_type from a file
  ! Requires: moveset  -- move set data structure
  !           filename -- file to read from
  !           nmoves   -- no of moves (already read from ctrlfile)
  !---------------------------------------------------------------
  Subroutine mcmoveset_readmswts(moveset,filename,nmoves)
    Type(Move_Set), Intent(InOut)            :: moveset
    Character(*), Intent(In)                 :: filename
    Integer, Intent(in)                      :: nmoves

    Integer                   :: i,nfields,unitno

    Character(len=strlen)     :: weightline
    Real(kind=RDbl), Dimension(strLen)        :: weights
    Character(len=strLen), Dimension(strLen)  :: strfields

    unitno = file_open(filename,110)

    !** Split the fields by commas
    Read(unitno,'(a)') weightline

    weightline = Trim(stripcmnt(weightline))
    nfields = split(weightline, strfields, ",") 

    If (nfields /= nmoves) Then
      Write(*,*) "Possibly something wrong in ctrlfile. no of moves and &
          & move weights dont match"
      Write(*,'(a,i4,a,i4)') "moves = ",nmoves,"  weights = ",nfields
      Write(*,*) " Check the moves section"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Do i = 1,nfields
      weights(i) = toreal(strfields(i))
    End Do

    Call mcmoveset_initms(moveset,nfields,weights)

  End Subroutine mcmoveset_readmswts

  !---------------------------------------------------------------
  ! Gets the tag given the move number
  ! Requires: moveset -- move set data structure
  !           moveno -- move number
  !---------------------------------------------------------------
  Function mcmoveset_getmstag(moveset,moveno)
    Character(len=strLen)                   :: mcmoveset_getmstag
    Type(Move_Set), Intent(In)              :: moveset
    Integer, Intent(In)                     :: moveno

    If (moveno > moveset%nmoves) Then
      mcmoveset_getmstag = moveset%tag(moveno)
    Else
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Passed move number too large"
      Stop      
    End If

  End Function mcmoveset_getmstag

  !---------------------------------------------------------------
  ! Gets the weight given the move number
  ! Requires: moveset -- move set data structure
  !           moveno -- move number
  !---------------------------------------------------------------
  Real(kind=RDbl) Function mcmoveset_getmswt_number(moveset,moveno)
    Type(Move_Set), Intent(In)              :: moveset
    Integer, Intent(In)                     :: moveno

    If (moveno <= moveset%nmoves) Then
      mcmoveset_getmswt_number = moveset%wt(moveno)
    Else
      Write(0,*) __FILE__," : ",__LINE__, &
          " Passed move number: ", moveno," is too large"
      Write(0,*) "actual number of moves", moveset%nmoves
      Stop      
    End If

  End Function mcmoveset_getmswt_number

  !---------------------------------------------------------------
  ! Gets the weight given the tag, returns -1.0 upon failure
  ! Requires: moveset -- move set data structure
  !           tag -- tag to get weight for
  !---------------------------------------------------------------
  Real(kind=RDbl) Function mcmoveset_getmswt_tag(moveset,tag)
    Type(Move_Set), Intent(In)     :: moveset
    Character(*), Intent(In)       :: tag

    Integer           :: moveno

    moveno = mcmoveset_getmoveno(moveset,tag)
    If (moveno == 0) Then
      mcmoveset_getmswt_tag = -1.0_RDbl
    Else
      mcmoveset_getmswt_tag = moveset%wt(moveno)
    End If

  End Function mcmoveset_getmswt_tag

  !---------------------------------------------------------------
  ! Gets the move number given the tag, returns zero at failure
  ! Requires: moveset -- move set data structure
  !           tag -- tag to search for
  !---------------------------------------------------------------
  Integer Function mcmoveset_getmoveno(moveset,tag)
    Type(Move_Set), Intent(In)           :: moveset
    Character(*), Intent(In)             :: tag

    Integer                 :: i

    mcmoveset_getmoveno = 0

    Do i = 1,moveset%nmoves
      If (Toupper(Trim(moveset%tag(i))) == Toupper(Trim(tag))) Then
        mcmoveset_getmoveno = i
      End If
    End Do

  End Function mcmoveset_getmoveno

  !---------------------------------------------------------------
  ! Normalize the weights and create the cumulative weights
  ! Requires: moveset -- move set data structure
  !---------------------------------------------------------------
  Subroutine mcmoveset_normms(moveset)
    Type(Move_Set), Intent(InOut)     :: moveset

    Integer                   :: i,error
    Real(kind=RDbl)           :: cum_sum

    If (.Not. Associated(moveset%wt)) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Passed Move Set weights must already be allocated"
      Stop      
    End If

    Nullify(moveset%cumwt)
    Allocate(moveset%cumwt(moveset%nmoves), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'moveset%cumwt')

    !** normalize the weights
    moveset%wt = moveset%wt/Sum(moveset%wt(1:moveset%nmoves))

    !** calculate the cumulative weights from the weights
    cum_sum = zero

    Do i = 1,moveset%nmoves
      cum_sum = cum_sum + moveset%wt(i)
      moveset%cumwt(i) = cum_sum
    End Do

  End Subroutine mcmoveset_normms

  !---------------------------------------------------------------
  ! Picks a move number randomly, returns the tag optionally
  ! Requires: moveset -- move set data structure
  !           tags -- optional returned tag
  !---------------------------------------------------------------
  Integer Function mcmoveset_pickmove(moveset,tag)
    Type(Move_Set), Intent(In)                   :: moveset
    Character(len=strLen), Intent(Out), Optional :: tag

    Integer                 :: i,moveno
    Real(kind=RDbl)         :: prob

    prob = rranf()
    mcmoveset_pickmove = 1
    Do i = 1,moveset%nmoves
      If (prob > moveset%cumwt(i)) Then
        mcmoveset_pickmove = mcmoveset_pickmove + 1
      Else 
        Exit
      End If
    End Do

    If (Present(tag)) Then
      tag = moveset%tag(mcmoveset_pickmove)
    End If

  End Function mcmoveset_pickmove

  !---------------------------------------------------------------
  ! Display the moveset structure
  ! Requires: moveset -- move set data structure
  !---------------------------------------------------------------
  Subroutine mcmoveset_displayms(moveset,indent,unit)
    Type(Move_Set), Intent(In)     :: moveset
    Integer, Intent(In)            :: indent,unit

    Integer                        :: i
    Character(len=indent)          :: blank

    blank = Repeat(' ',indent)

    Write(unit,'(2a,i3)') blank,"Move Set number of moves: ",moveset%nmoves

    If (Associated(moveset%wt)) Then
      Write(unit,'(2x,2a,15(f6.4,2x))') blank,'weights: ',&
          moveset%wt(1:moveset%nmoves)
    End If
    If (Associated(moveset%cumwt)) Then
      Write(unit,'(2x,2a,15(f6.4,2x))') blank,'cum weights: ',&
          moveset%cumwt(1:moveset%nmoves)
    End If
    If (Associated(moveset%tag)) Then
      Write(unit,'(2x,2a,15(a,2x))') blank,'tags: ',&
          (Trim(moveset%tag(i)),i=1,moveset%nmoves)
    End If

  End Subroutine mcmoveset_displayms


  !------------------------------------------------------------------
  ! Sets the ratio of insertions vs deletions. This will be used in 
  ! the monte carlo acceptanece criteria. 
  ! Requires: params -- move set data structure
  !------------------------------------------------------------------
  Subroutine mcmoveset_setInsDelRatio(params)
    Type(Move_Set), Intent(InOut)     :: params

    Integer             :: ins_move, del_move,moveno
    Logical             :: found_ins, found_del, error_flag
    Real(kind=RDbl)     :: ins_del_ratio,ins_wt,del_wt

    found_ins=.false.
    found_del=.false.
    error_flag=.false.
    Do moveno=1,params%nmoves
      If (params%tag(moveno)=="INSERT") Then
        !** make sure one move was not already found previously
        !** If this is second move then exit and give error mesage
        If (found_ins) Then
          found_ins=.False.
          error_flag=.True.
          Exit
        Else
          found_ins=.True.
          ins_move=moveno
        Endif
      Elseif(params%tag(moveno)=="DELETE") Then
        ! same here, avoiding more than one delete move
        If (found_del) Then
          found_del=.False.
          error_flag=.true.
          Exit
        Else
          found_del=.True.
          del_move=moveno
        Endif
      Endif
    End Do

    If (.Not.(found_ins.And.found_del)) error_flag=.True.

    If (error_flag) Then
      Write(*,*) "WARNING :--------------------------------------------"
      Write(*,*) "This Routine should not have been called in the -"
      Write(*,*) "absence of insertion and deletion moves. Check whether -"
      Write(*,*) "move_tags are set. Also check that there are only one -"
      Write(*,*) " deletion and insertion movetype each"
      Write(*,*) "found_ins = ",found_ins,"; found_del = ", found_del,";"
      Write(*,*) "Setting ins_del_ratio=1"
      Write(*,*) "WARNING :--------------------------------------------"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      ins_del_ratio=one
      params%ins_del_ratio=ins_del_ratio
      Return
    End If
    ins_wt=mcmoveset_getmswt(params,ins_move)
    del_wt=mcmoveset_getmswt(params,del_move)
    ins_del_ratio=ins_wt/del_wt

!    Write(*,*) "Insertion/Deletion ratio calculated :",ins_del_ratio

    params%ins_del_ratio = ins_del_ratio

  End Subroutine mcmoveset_setInsDelRatio

  !--------------------------------------------------------------
  ! Returns the ratio of [insert attempts]/[delete attempts]
  ! Requires: params -- move set data structure
  !--------------------------------------------------------------
  Real(kind=RDbl) Function mcmoveset_insdelratio(params)
    Type(Move_Set), Intent(InOut)     :: params

    mcmoveset_insdelratio=params%ins_del_ratio
  End Function mcmoveset_insdelratio

End Module mcmoveset
