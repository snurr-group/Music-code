!----------------------------------------------------------------------------
! This module contains the different data structures and routines for 
! handling the CBGCMC move types.  The data type "CBGCMC_Move_Params" contains
! the simulation information for a single species type.  Move-specific
! information and handling are done using the "moves" module through the
! "mcmoves" module.
!
! Important routines are:
!    cbgcmoves_init -- handles initilization for each species
!    cbgcmoves_move -- picks a movetype and executes
!    cbgcmoves_displaystats -- displays intermediate statistics for species
!----------------------------------------------------------------------------
Module cbgcmoves

  Use branchedcoords, Only: BranchedMolec, Nodetype   
  Use brmoves, Only: BranchedCBIns_Params
  Use config,Only:Atmolcoords,config_getnmoles,config_getnatoms,&
      config_allocfields,config_checkandincr,config_setnmoles,&
      config_delmolec,config_copymolec
  Use defaults, Only: RDbl, strLen, MAX_EXCLSITES, MAX_SORBS,kcalmole_kb , &
      zero,caltoj, one, dashedline, dashedline2,NO_OF_INTRA_POTS,TOTAL_INDEX, &
      default_MAX_EXP, INITIAL_SIZE
  Use file, Only: file_getunit, file_gettype,file_open
  Use gcmodels, Only: GeneralizedCoords,gcmodels_getcom
  Use insert, Only: CBInsert_Params  
  Use intramolecular, Only: intramolecular_getint
  Use mcmoveset, Only: Move_Set,mcmoveset_initms,mcmoveset_readmswts, &
      mcmoveset_pickmove,mcmoveset_displayms,mcmoveset_setmstag, &
      mcmoveset_getmswt, mcmoveset_setInsDelRatio, mcmoveset_insdelratio, &
      mcmoveset_getmoveno
  Use mcmoves, Only: MC_Move_Params,AuxReject_Params,mcmoves_init, &
      mcmoves_initaux,mcmoves_perturb, &
      mcmoves_insert,mcmoves_delete, &
      mcmoves_getbasictag,mcmoves_display,mcmoves_auxrejectdisplay, &
      mcmoves_stats, mcmoves_resetStats, mcmoves_setMaxNrg
  Use molecules, Only: molecules_name, molecules_getnsorbs, & 
      molecules_gettype, molecules_getnatoms, molecules_getgcmodeltype, &
      molecules_getatype,molecules_incrcoulnrg,molecules_incrnoncoulnrg, &
      molecules_incrintranrg

  Use moves,Only: Move_Params,Move_stats, moves_initparams,&
      moves_insert,moves_delete,moves_rotate, moves_translate, &
      moves_transform, moves_displayparams,moves_adjustparams, &
      moves_setMaxNrg 
  Use random,Only :rranf
  Use stats,Only: Statistics, stats_setvalue, stats_getvalue, stats_getcavg, &
      stats_update, stats_display, stats_init
  Use simcell, Only: SimCell_Params, simcell_pbc
  Use ssnoncoul,Only: ssnoncoul_msdriver
  Use sscoul,Only: sscoul_msdriver
  Use transform, Only: Regrow_Params  
  Use utils, Only: isfileopen, filesrchstr, toupper, split, stripcmnt, &
      toreal, toint, allocErrDisplay
  Use vector, Only: VecType, Assignment(=)

  Implicit None
  Save

  Private
  Public :: CBGCMC_Move_Params, Fugacity_Params,cbgcmoves_init, & 
!!      cbgcmoves_insert, cbgcmoves_delete,cbgcmoves_rot_trans, &
!!      cbgcmoves_transform, &
      cbgcmoves_move, cbgcmoves_updatenmolecs, &
      cbgcmoves_displaystats,cbgcmoves_display,cbgcmoves_displaysimparams, &
      cbgcmoves_dispspcstats,cbgcmoves_beginsim, cbgcmoves_nvtmove
  !!    cbgcmoves_getrandommovetype, 
  !!    

  Type Fugacity_Params
    Real(kind=RDbl)   :: pressure
    Real(kind=RDbl)   :: murti       ! Chemical Potential/RT
    Real(kind=RDbl)   :: fugacity
    Real(kind=RDbl)   :: B        ! Adams Notation
  End Type Fugacity_Params

  Type CBGCMC_Move_Params
    Integer                   :: sorbtype
    Integer                   :: npts
    Integer                   :: nunitcells
    Type(Statistics)          :: nmoles
    Integer                   :: blocksize ! Block size for collecting stats

    !** Thermophysical Properties
    Real(kind=RDbl)           :: tk
    Real(kind=RDbl)           :: rti   ! 1.0/(Rgas*tk)
    Character(len=strLen)     :: pressurefile
    Type(Fugacity_Params), Dimension(:), Pointer   :: fuglist

    !** Move information: tags can be "INSERT","DELETE","ROTATE" or "TRANSLATE"
    Integer                   :: no_of_movetypes
    Type(Move_Set)                            :: moveset
    Type(MC_Move_Params),Dimension(:),Pointer :: mcmoves
    Type(AuxReject_Params), Pointer           :: auxparams

    Real(kind=RDbl)                           :: max_nrg 

  End Type CBGCMC_Move_Params

  !** maximum to be used as the exponent 
  !** in exponentials, during insert, delete, translate and rotate moves
  Real(kind=RDbl) , Parameter :: EXP_MAX=default_MAX_EXP

Contains
  !----------------------------------------------------------------
  ! Initializes the various CBGCMC parameters from the control file
  ! Requires: cbgcmcsorbs -- parameters for individual species
  !           sorbates -- species data structure
  !           simcell -- the simulation cell information
  !           ctrl_filename -- the name of the control file
  !----------------------------------------------------------------
  Subroutine cbgcmoves_init(cbgcmcsorbs, sorbates, simcell, ctrl_filename)
    Type(CBGCMC_Move_Params), Dimension(:), Intent(inout) :: cbgcmcsorbs
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(in)  :: simcell
    Character(*), Intent(in)          :: ctrl_filename

    Real(kind=RDbl)      :: tk
    Real(kind=RDbl), Dimension(3) :: ell
    Real(kind=RDbl)      :: inswt,delwt
    Integer         :: unitno, nsorbs, i, sorbtype, npts, error, blocksize,l
    Integer         :: nmoles, nfields, j,no_of_movetypes,nsorbates
    Character(len=strLen) :: sorbname, pressureline, gcmodeltype
    Character(len=strLen) :: line, exclusionsites
    Character(len=strLen), Dimension(MAX_EXCLSITES) :: fields
    Type(NodeType), Pointer    :: node  

    !** from the control file read no of sorbates
    unitno = file_getunit(ctrl_filename)
    nsorbs = Size(cbgcmcsorbs, 1)

    !** initialise the parameters for each adsorbate
    Do i=1, nsorbs

      Read(unitno, *) sorbname
      Write(*,*) "Initializing cbgcmoves section for sorbate "//Trim(sorbname)
      sorbname=Trim(stripcmnt(sorbname))
      sorbtype = molecules_gettype(sorbname)
      If (sorbtype == 0) Then
        Write(0,'(1x,2a,i4,1x,2a)') __FILE__," : ",__LINE__, &
            Trim(sorbname), " is not one of the sorbate types in molecule list"
        Stop
      End If

      cbgcmcsorbs(i)%sorbtype = sorbtype
      
      ! make sure that sorbates has enough memoery, sorbates fields will be 
      ! accessed during brmoves_init
      nmoles = config_getnmoles(sorbates, sorbtype)
      Call config_checkandincr(sorbates, sorbtype, nmoles)
      
      ! auxillary rejection criteria like sitemaps etc..
      Nullify(cbgcmcsorbs(i)%auxparams)

      !** Initialize the pressure, fugacity etc.
      npts = cbgcmcsorbs(i)%npts

      Allocate(cbgcmcsorbs(i)%fuglist(npts), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"fuglist")

      !** different pressure values in the isotherm
      Read(unitno,'(a)') pressureline
      Call cbgcmoves_getpressures(cbgcmcsorbs(i), pressureline)      

      !** Read the auxilliary parameters for rejections
      Call mcmoves_initaux(cbgcmcsorbs(i)%auxparams,ctrl_filename,simcell)

      !** Check whether we are using the correct gcmodel type
      gcmodeltype = molecules_getgcmodeltype(sorbtype)
      If (Trim(toupper(gcmodeltype))/="BRANCHED") Then
        Write(*,*) "The gcmodel :"//Trim(gcmodeltype)//" cant be used for &
            & cbgcmc type of moves, need -BRANCHED gcmodel"
        Stop
      Endif

      tk  = cbgcmcsorbs(i)%tk

      !** Read the line containing number of movetypes
      Read(unitno,*) no_of_movetypes
      cbgcmcsorbs(i)%no_of_movetypes=no_of_movetypes
      Allocate(cbgcmcsorbs(i)%mcmoves(no_of_movetypes), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'mvparams')
      Write(*,*) "Reading the move-weights line"
      Call mcmoveset_readmswts(cbgcmcsorbs(i)%moveset,ctrl_filename, &
          no_of_movetypes)
      Write(*,*) "Finished Reading the move-weights line"

      Do j=1,no_of_movetypes
        Write(*,'(a,i3)') "Initializing Movetype : ",j
        Call mcmoves_init(cbgcmcsorbs(i)%mcmoves(j),ctrl_filename,sorbtype, &
            simcell, sorbates,cbgcmcsorbs(i)%auxparams)
        Call mcmoves_setMaxNrg(cbgcmcsorbs(i)%mcmoves(j), &
            cbgcmcsorbs(i)%max_nrg)
        Call mcmoveset_setmstag(cbgcmcsorbs(i)%moveset,j, &
            mcmoves_getbasictag(cbgcmcsorbs(i)%mcmoves(j)))
      End Do

      !** Calculate Beta (1/(Rgas*tk))
      cbgcmcsorbs(i)%rti = 1.0_RDbl/(kcalmole_kb*tk)

      !** Check the input move types (insert = delete, existence)
      inswt = mcmoveset_getmswt(cbgcmcsorbs(i)%moveset,'INSERT') 
      If (inswt <= zero) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Move set for ',Trim(sorbname),' must contain an INSERT move'
!!        Stop
      End If

      delwt = mcmoveset_getmswt(cbgcmcsorbs(i)%moveset,'DELETE') 
      If (delwt <= zero) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Move set for ',Trim(sorbname),' must contain a DELETE move'
!!        Stop
      End If

      If (Abs((inswt/delwt)-one) > 1.0e-6) Then
        Write(*,'(1x,3a,i4)') "WARNING : ",__FILE__," : ",__LINE__
        Write(0,*) "WARNING : different insertion / deletion rates "
        Write(0,*) "WARNING : Make sure that this is correct for your Case "
        !** allow for having different insertion and deletion rates
      End If

      !** Ratio of insertio attempts/ delete attempts
      Call mcmoveset_setInsDelRatio(cbgcmcsorbs(i)%moveset)
      
      If (mcmoveset_getmswt(cbgcmcsorbs(i)%moveset,'TRANSLATE') <= zero) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Move set for ',Trim(sorbname),' must contain a TRANSLATE move'
!!        Stop
      End If
 !SDEBUG
      Write(*,*) __FILE__," : ",__LINE__,"Check config_alloc"

      !** Read the blank separating line
      If (i /= nsorbs) Then 
        Read(unitno,*)
      End If

    End Do !** loop over all the sorbates

  End Subroutine cbgcmoves_init

  !------------------------------------------------------
  ! Initialize the stats variable used for display during the simulation 
  !--------------------------s----------------------------
  Subroutine cbgcmoves_beginsim(cbgcmcsorb, sorbates,simno)
    Type(CBGCMC_Move_Params), Intent(inout)  :: cbgcmcsorb
    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Integer,Intent(in) :: simno

    !** set stats variables, attempt/success values
    Call cbgcmoves_initRunStats(cbgcmcsorb,sorbates(cbgcmcsorb%sorbtype) )

  End Subroutine cbgcmoves_beginsim

  !------------------------------------------------------
  ! Initializ the stats variable used for display during the simulation 
  !------------------------------------------------------
  Subroutine cbgcmoves_initRunStats(cbgcmcsorb, sorbate)
    Type(CBGCMC_Move_Params), Intent(inout)  :: cbgcmcsorb
    Type(AtMolCoords), Intent(in) :: sorbate

    Integer :: blocksize,nmoles,i

    blocksize = cbgcmcsorb%blocksize

    nmoles = config_getnmoles(sorbate)

    Call stats_init(cbgcmcsorb%nmoles, &
        "# Moles (ins, blk, cum, stdd)", blocksize, "f7.2")

    Call stats_setvalue(cbgcmcsorb%nmoles, nmoles*one)
    Do i=1,cbgcmcsorb%no_of_movetypes
      Call mcmoves_resetStats(cbgcmcsorb%mcmoves(i))
    End Do


  End Subroutine cbgcmoves_initRunStats



  !------------------------------------------------------
  ! Set the simulation temperature for each sorbate type
  !------------------------------------------------------
  Subroutine cbgcmoves_settemp(cbgcmcsorb, tk)
    Type(CBGCMC_Move_Params), Intent(inout)  :: cbgcmcsorb
    Real(kind=RDbl), Intent(in)   :: tk
    cbgcmcsorb%tk = tk
    Return
  End Subroutine cbgcmoves_settemp

  !------------------------------------------------------
  ! Set the number of simulation points.
  !------------------------------------------------------
  Subroutine cbgcmoves_setnpts(cbgcmcsorb, npts)
    Type(CBGCMC_Move_Params), Intent(inout)  :: cbgcmcsorb
    Integer, Intent(in)   :: npts

    cbgcmcsorb%npts = npts
    Return
  End Subroutine cbgcmoves_setnpts

  !---------------------------------------------------------
  ! Parses the fugacity line to set the fugacities of the
  ! sorbate type.  It can either do the fugacities at equal
  ! intervals or read them from a file.
  !---------------------------------------------------------
  Subroutine cbgcmoves_getpressures(cbgcmcsorb, pressureline)
    Type(CBGCMC_Move_Params), Intent(inout) :: cbgcmcsorb
    Character(*), Intent(in)  :: pressureline
    Real(kind=RDbl) :: startp, endp
    Integer         :: i, nfields, npts, convfields, error
    Character(len=strLen)     :: strfields(10), filename, str

    ! Split the fields by commas
    nfields = split(pressureline, strfields, ",") 
    ! If nfields is 1, it's possible that it is not a file but
    ! the fugacity we want IF npts is 1.
    If ((nfields == 1).And.(cbgcmcsorb%npts == 1)) Then
      ! Check if the line is an integer 
      startp = toint(strfields(1),error)
      If (error /= 0) Then
        ! It's a file. Read it.
        nfields = split(pressureline, strfields, " ")
        cbgcmcsorb%pressurefile = strfields(1)
        Call cbgcmoves_readpressures(cbgcmcsorb)
      Else
        ! It's a number. Store it.
        cbgcmcsorb%fuglist(1)%pressure = startp
      End If
    Else If (nfields == 1) Then
      ! We need to read the points from a file
      ! Get the filename
      nfields = split(pressureline, strfields, " ")
      cbgcmcsorb%pressurefile = strfields(1)
      Call cbgcmoves_readpressures(cbgcmcsorb)
    Else
      ! We need to generate the range
      startp = toreal(strfields(1))

      ! Get the last field after separating it from the rest of the line
      str = strfields(2)
      nfields = split(str, strfields, " ")
      endp = toreal(strfields(1))

      ! Fill in the cbgcmcsorbs pressure data
      Call cbgcmoves_genpressures(cbgcmcsorb, startp, endp)
    End If
  End Subroutine cbgcmoves_getpressures


!!$  !---------------------------------------------------------
!!$  ! gets the weights each move_type from the control file
!!$  !---------------------------------------------------------
!!$  Subroutine cbgcmoves_getmvWeights(cbgcmcsorb, unitno)
!!$    Type(CBGCMC_Move_Params), Intent(inout) :: cbgcmcsorb
!!$    Integer, Intent(in)                   :: unitno
!!$    Character(len=strlen)              :: weightline
!!$    Integer         :: i, nfields
!!$    Character(len=strLen)     :: strfields(cbgcmcsorb%no_of_movetypes &
!!$        )
!!$    Real(kind=RDbl) :: sum,cum_sum
!!$
!!$    ! Split the fields by commas
!!$    Read(unitno,'(a)') weightline
!!$    weightline=Trim( stripcmnt(weightline) )
!!$    nfields = split(weightline, strfields, ",") 
!!$    If (nfields >=cbgcmcsorb%no_of_movetypes) Then
!!$      sum=zero
!!$      Do i=1,cbgcmcsorb%no_of_movetypes
!!$        cbgcmcsorb%mvweights(i)=toreal(strfields(i))
!!$        sum=sum+cbgcmcsorb%mvweights(i)
!!$      End Do
!!$      !** convert weights to cumulative weights and normalise weights
!!$      cum_sum=zero
!!$      Do i=1,cbgcmcsorb%no_of_movetypes
!!$        cum_sum=cum_sum+cbgcmcsorb%mvweights(i)
!!$        cbgcmcsorb%mvweights(i)=cum_sum/sum
!!$      End Do
!!$    Else
!!$      Write(0,*) " No of weights given is less than ", &
!!$          cbgcmcsorb%no_of_movetypes 
!!$      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__," not able&
!!$          & to understand weights line: ",trim(weightline)
!!$      Stop
!!$    End If
!!$    Write(*,*) "gccmsorb%mvweights : ", cbgcmcsorb%mvweights 
!!$  End Subroutine cbgcmoves_getmvWeights

!!$  !---------------------------------------------------------
!!$  ! sets a Tag corresponding to cbgcmcsorb%mvparams(moveno)
!!$  !---------------------------------------------------------
!!$  Subroutine cbgcmoves_setmvTag(cbgcmcsorb, moveno)
!!$    Type(CBGCMC_Move_Params), Intent(inout) :: cbgcmcsorb
!!$    Integer, Intent(in)                   :: moveno
!!$    Integer         :: i, nfields
!!$    Character(len=strLen)     :: strfields(cbgcmcsorb%no_of_movetypes &
!!$        )
!!$
!!$    ! Split the fields by commas
!!$    If (Associated(cbgcmcsorb%mvparams(moveno)%BInsParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="INSERT"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%CBInsParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="INSERT"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%LInsParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="INSERT"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%RInsParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="INSERT"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%FTransParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="TRANSLATE"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%RTransParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="TRANSLATE"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%FRotParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="ROTATE"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%RRotParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="ROTATE"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%RegrowParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="CUT_REGROW"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%BDelParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="DELETE"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%CBDelParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="DELETE"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%RDelParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="DELETE"
!!$    Else If (Associated(cbgcmcsorb%mvparams(moveno)%IntegParams)) Then
!!$      cbgcmcsorb%mvTags(moveno)="INTEGRATE"
!!$    Else
!!$      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,"not able&
!!$          &to decide a tag for the given movetype"
!!$      Stop
!!$    End If
!!$  End Subroutine cbgcmoves_setmvTag


  !---------------------------------------------------------
  ! Read the fugacities from the file "filename"
  !---------------------------------------------------------
  Subroutine cbgcmoves_readpressures(cbgcmcsorb)
    Type(CBGCMC_Move_Params), Intent(inout) :: cbgcmcsorb

    Integer     :: unitno, sorbtype, lineno, npts, i, error
    Character(len=strLen)    :: line, sorbname, filename

    !** Open the file if not already open
    filename = cbgcmcsorb%pressurefile
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

    !** Find the sorbate name in the file
    sorbtype = cbgcmcsorb%sorbtype
    sorbname = molecules_name(sorbtype)
    lineno = filesrchstr(unitno, Trim(sorbname), line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4, 4a)') __FILE__," : ",__LINE__, &
          " Could not find the sorbate ", Trim(sorbname), " in the file ",&
          Trim(filename)
      Stop
    End If

    !** Read the no. of points
    Read(unitno, *) npts
    ! Check that the no. of points is the same as that specfied
    ! in the control file
    If (npts /= cbgcmcsorb%npts) Then
      Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
          " The no. of points in the file ", Trim(filename), &
          " does not match that in the control file"
      Stop
    End If
    Read(unitno, *) (cbgcmcsorb%fuglist(i)%pressure, i=1, npts)
    Close(unit=unitno)

  End Subroutine cbgcmoves_readpressures

  !-----------------------------------------------------------
  ! Generates the actual fugacity list from the given intervals
  ! contains the start fugacity, end fugacity and no. of points
  !--------------------------------------------------------
  Subroutine cbgcmoves_genpressures(cbgcmcsorb, startp, endp)
    Type(CBGCMC_Move_Params), Intent(inout) :: cbgcmcsorb
    Real(kind=RDbl), Intent(in)           :: startp, endp
    Real(kind=RDbl)     :: pressincr 
    Integer             :: i, npts

    ! Check to see what kind of a scale we are going to use
    ! A log scale is used if the ending fugacity is greater
    ! than the starting fugacity by 2 orders of magnitude
    npts = cbgcmcsorb%npts
    If (Abs(Log10(endp/startp)) > 2) Then
      !use a log scale
      If (npts /= 1) Then
        pressincr = Exp(Log(endp/startp)/(npts-1.0_RDbl))
      Endif
      cbgcmcsorb%fuglist(1)%pressure = startp
      Do i=1, npts-1
        cbgcmcsorb%fuglist(i+1)%pressure = startp*(pressincr**i)
      End Do
    Else
      If (npts /= 1) Then
        pressincr = (endp - startp)/(npts - 1.0_RDbl)
      End If
      cbgcmcsorb%fuglist(1)%pressure = startp
      Do i=1, npts-1
        cbgcmcsorb%fuglist(i+1)%pressure = startp + pressincr*i
      End Do
    Endif
    Return
  End Subroutine cbgcmoves_genpressures


  !--------------------------------------------------------------------
  ! Pick a move and make it using the mcmoves module.  Returns true if
  ! move was successful.
  ! Requires: cbgcmcspc -- parameters for individual species
  !           moveno -- the move number to implement
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           simno -- the simulation number
  !--------------------------------------------------------------------
  Logical Function cbgcmoves_move(cbgcmcspc,species,simcell,simno)
    Type(CBGCMC_Move_Params), Intent(InOut)           :: cbgcmcspc
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Type(SimCell_Params), Intent(In)                :: simcell
    Integer, Intent(In)                             :: simno

    Integer                   :: nmoles,spc,moveno,molec
    Character(len=strLen)     :: movename
    Real(kind=RDbl)           :: insdelratio

    spc = cbgcmcspc%sorbtype
    moveno = mcmoveset_pickmove(cbgcmcspc%moveset,movename)
    nmoles = config_getnmoles(species,cbgcmcspc%sorbtype)

    Select Case(Trim(movename))
    Case("INSERT")
      insdelratio=mcmoveset_insdelratio(cbgcmcspc%moveset)
      cbgcmoves_move = mcmoves_insert(cbgcmcspc%mcmoves(moveno),molec, &
          cbgcmcspc%rti,cbgcmcspc%fuglist(simno)%B,species,simcell, &
          insdelratio)

    Case("DELETE")
      molec = Int(rranf()*nmoles) + 1

      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
      insdelratio=mcmoveset_insdelratio(cbgcmcspc%moveset)
        cbgcmoves_move = mcmoves_delete(cbgcmcspc%mcmoves(moveno),molec, &
            cbgcmcspc%rti,cbgcmcspc%fuglist(simno)%B,species,simcell,&
            insdelratio)
      Else
        cbgcmoves_move = .False.
      End If

    Case("TRANSLATE")
      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
        molec = Int(rranf()*nmoles) + 1
        cbgcmoves_move = mcmoves_perturb(cbgcmcspc%mcmoves(moveno),molec, &
            cbgcmcspc%rti,species,simcell)
      Else
        cbgcmoves_move = .False.
      End If

    Case("ROTATE")
      molec = Int(rranf()*nmoles) + 1

      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
      cbgcmoves_move = mcmoves_perturb(cbgcmcspc%mcmoves(moveno),molec, &
          cbgcmcspc%rti,species,simcell)
      Else
        cbgcmoves_move = .False.
      End If

    Case ("TRANSFORM")
        molec = Int(rranf()*nmoles) + 1

      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
        cbgcmoves_move = mcmoves_perturb(cbgcmcspc%mcmoves(moveno), &
            molec, cbgcmcspc%rti,species,simcell)
      Else
        cbgcmoves_move = .False.
      End If
      
    Case Default
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Could not Identify movename : ", Trim(movename)
      Stop
    End Select

#ifdef DEBUG  !** very useful debugging feedback
    Write(*,'(2a,t34,a,i4,t55,a,l4)') 'GCMC feedback: move = ', &
        Trim(movename),'   molecule = ',molec,'    Accepted? ',gcmcmoves_move
    Write(*,*)
#endif

  End Function cbgcmoves_move



  !--------------------------------------------------------------------
  ! Does only nvt moves. This routine does not use mcmoveset to pick 
  ! the movetype randomly
  ! This is special to simulations with cbgcmc-accelerator on.
  ! Requires: cbgcmcspc -- parameters for individual species
  !           moveno -- the move number to implement
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           simno -- the simulation number
  !--------------------------------------------------------------------
  Logical Function cbgcmoves_nvtmove(cbgcmcspc,species,simcell,simno)
    Type(CBGCMC_Move_Params), Intent(InOut)           :: cbgcmcspc
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Type(SimCell_Params), Intent(In)                :: simcell
    Integer, Intent(In)                             :: simno

    Integer                   :: nmoles,spc,moveno,molec
    Character(len=strLen)     :: movename
    Real(kind=RDbl)           :: insdelratio

    spc = cbgcmcspc%sorbtype
    nmoles = config_getnmoles(species,cbgcmcspc%sorbtype)

    If (rranf()>0.5) then
      !** Make sure that the no. of molecules is not zero
      moveno = mcmoveset_getmoveno(cbgcmcspc%moveset,"TRANSLATE")
      If (nmoles /= 0) Then
        molec = Int(rranf()*nmoles) + 1
        cbgcmoves_nvtmove = mcmoves_perturb(cbgcmcspc%mcmoves(moveno),molec, &
            cbgcmcspc%rti,species,simcell)
      Else
        cbgcmoves_nvtmove = .False.
      End If
    Else
      moveno = mcmoveset_getmoveno(cbgcmcspc%moveset,"TRANSFORM")
      molec = Int(rranf()*nmoles) + 1

      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
        cbgcmoves_nvtmove = mcmoves_perturb(cbgcmcspc%mcmoves(moveno), &
            molec, cbgcmcspc%rti,species,simcell)
      Else
        cbgcmoves_nvtmove = .False.
      End If
    Endif

  End Function cbgcmoves_nvtmove






  !-----------------------------------------------------------------
  ! Update the number of moles for all species
  ! Requires: gcmcspc -- array of GCMC species parameters
  !           species -- species data structure
  !-----------------------------------------------------------------
  Subroutine cbgcmoves_updatenmolecs(params,species)
    Type(CBGCMC_Move_Params), Dimension(:), Intent(InOut) :: params
    Type(AtMolCoords), Dimension(:), Intent(In)         :: species

    Integer            :: sorbno,nmoles

    Do sorbno = 1,Size(params,1)
      nmoles = config_getnmoles(species,params(sorbno)%sorbtype)
      Call stats_update(params(sorbno)%nmoles,nmoles*one)      
    End Do

  End Subroutine cbgcmoves_updatenmolecs

  !----------------------------------------------------
  ! Cleans up the various pointers once we are done with
  ! this structure.
  !---------------------------------------------------- 
  Subroutine cbgcmoves_cleanup(cbgcmcsorbs)
    Type(CBGCMC_Move_Params), Dimension(:), Intent(inout) :: cbgcmcsorbs

    Integer    :: i, nsorbs, error

    nsorbs = Size(cbgcmcsorbs, 1)
    Do i=1, nsorbs
      Deallocate(cbgcmcsorbs(i)%fuglist, STAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,a, i4)') __FILE__," : ",__LINE__, &
            "Could not deallocate 'fugacitylist' for sorbate", i
      End If
      Stop
    End Do
  End Subroutine cbgcmoves_cleanup

  !-------------------------------------------------------------
  ! Display the simulation parameters
  !-------------------------------------------------------------
  Subroutine cbgcmoves_displaystats(cbgcmcsorbs, simno, indent, optunit)
    Type(CBGCMC_Move_Params), Dimension(:), Intent(in) :: cbgcmcsorbs
    Integer, Intent(in)   :: simno, indent
    Integer, Optional, Intent(in) :: optunit
    Character(len=indent) :: blank
    Integer     :: i, j, unitno, nsorbs, sorbtype 
    Character(len=strLen)   :: molecname
    Real(kind=RDbl)         :: insertratio, deleteratio,transratio,rotratio
    Real(kind=RDbl)         :: transfratio

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    nsorbs = Size(cbgcmcsorbs, 1)
    Do i=1, nsorbs
      sorbtype = cbgcmcsorbs(i)%sorbtype
      molecname = molecules_name(sorbtype)
      Write(unitno,'(4a,f7.2,a)') blank,"Species Name: ",Trim(molecname), &
          " at ",cbgcmcsorbs(i)%fuglist(simno)%pressure," kPa"

      Do j=1,cbgcmcsorbs(i)%no_of_movetypes
        Call mcmoves_stats(cbgcmcsorbs(i)%mcmoves(j), indent+2, unitno)
      End Do
      Call stats_display(cbgcmcsorbs(i)%nmoles, indent+2)

      
!!$      If (cbgcmcsorbs(i)%insert%att /= 0) Then
!!$       insertratio = cbgcmcsorbs(i)%insert%succ/(cbgcmcsorbs(i)%insert%att*1.0)
!!$      Else
!!$        insertratio = 0.0_RDbl
!!$      End If
!!$      If (cbgcmcsorbs(i)%delete%att /= 0) Then
!!$        deleteratio = cbgcmcsorbs(i)%delete%succ/(cbgcmcsorbs(i)%delete%att*1.0)
!!$      Else
!!$        deleteratio = 0.0_RDbl
!!$      End If
!!$      If (cbgcmcsorbs(i)%trans%att /= 0) Then
!!$        transratio = cbgcmcsorbs(i)%trans%succ/(cbgcmcsorbs(i)%trans%att*1.0)
!!$      Else
!!$        transratio = 0.0_RDbl
!!$      End If
!!$      If (cbgcmcsorbs(i)%rot%att /= 0) Then
!!$        rotratio = cbgcmcsorbs(i)%rot%succ/(cbgcmcsorbs(i)%rot%att*1.0)
!!$      Else
!!$        rotratio = 0.0_RDbl
!!$      End If
!!$      If (cbgcmcsorbs(i)%transform%att /= 0) Then
!!$        transfratio = cbgcmcsorbs(i)%transform%succ/ &
!!$            (cbgcmcsorbs(i)%transform%att*1.0)
!!$      Else
!!$        transfratio = 0.0_RDbl
!!$      End If
!!$
!!$
!!$      Write(unitno, '(a,2a)') spc,"Sorbate Name        : ",  Trim(molecname)
!!$      Write(unitno, '(a,2x,a,e11.5,1x,a,i11,a,i11,a)') &
!!$          spc,"Ratio (tranform)  : ", transfratio,   &
!!$          "(",cbgcmcsorbs(i)%transform%succ,"/",&
!!$          cbgcmcsorbs(i)%transform%att,")"
!!$      Write(unitno, '(a,2x,a,e11.5,1x,a,i11,a,i11,a)') &
!!$          spc,"Ratio (insert)    : ", insertratio,   &
!!$          "(",cbgcmcsorbs(i)%insert%succ,"/",cbgcmcsorbs(i)%insert%att,")"
!!$      Write(unitno, '(a,2x,a,e11.5,1x,a,i11,a,i11,a)') &
!!$          spc,"Ratio (delete)    : ", deleteratio,   &
!!$          "(",cbgcmcsorbs(i)%delete%succ,"/",cbgcmcsorbs(i)%delete%att,")"
!!$      Write(unitno, '(a,2x,a,e11.5,1x,a,i11,a,i11,a)') &
!!$          spc,"Ratio (translate) : ", transratio,   &
!!$          "(",cbgcmcsorbs(i)%trans%succ,"/",cbgcmcsorbs(i)%trans%att,")"
!!$      Write(unitno, '(a,2x,a,e11.5,1x,a,i11,a,i11,a)') &
!!$          spc,"Ratio (rotate)    : ", rotratio,   &
!!$          "(",cbgcmcsorbs(i)%rot%succ,"/",cbgcmcsorbs(i)%rot%att,")"
!!$      Write(unitno,'(a,2x,a,a,i11,a,i11,a)') &
!!$          spc,"Nrg Calc success, ins : ",  "(",&
!!$          cbgcmcsorbs(i)%insert%nrg_succ,"/",cbgcmcsorbs(i)%insert%att,")"
!!$      Call stats_display(cbgcmcsorbs(i)%nmoles, nspc+2)
!!$

    End Do

  End Subroutine cbgcmoves_displaystats

  !-------------------------------------------------------------
  ! Display the simulation parameters
  !-------------------------------------------------------------
  Subroutine cbgcmoves_displaysimparams(cbgcmcsorbs, simno, nspc, optunit)
    Type(CBGCMC_Move_Params), Dimension(:), Intent(in) :: cbgcmcsorbs
    Integer, Intent(in)   :: simno
    Integer, Intent(in)   :: nspc
    Integer, Optional, Intent(in) :: optunit
    Character(len=nspc) :: spc
    Integer     :: i, j, unitno, nsorbs, sorbtype 
    Character(len=strLen)   :: molecname

    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    nsorbs = Size(cbgcmcsorbs, 1)
    Write(unitno,'(2a)') spc, dashedline
    Write(unitno,'(2a)') spc, "CBGCMC Sorbate Parameters:"

    !** if not associated nsorbs=0, then does not go into the loop
    Do i=1, nsorbs
      sorbtype = cbgcmcsorbs(i)%sorbtype
      molecname = molecules_name(sorbtype)
      Write(unitno, '(a,2x,a)')  spc, dashedline2
      Write(unitno, '(a,2x,2a)') spc,"Sorbate Name        : ",  Trim(molecname)
      Write(unitno, '(a,4x, 4a11)') spc,"Press.(kPa)", "Fugacity", "murti", &
          "B(Adams)"
      Write(unitno, '(a,4x,4f11.3)') spc, &
          cbgcmcsorbs(i)%fuglist(simno)%pressure, &
          cbgcmcsorbs(i)%fuglist(simno)%fugacity, &
          cbgcmcsorbs(i)%fuglist(simno)%murti, &
          cbgcmcsorbs(i)%fuglist(simno)%B
    End Do
    Write(unitno, '(a,2x,a)')  spc, dashedline2
  End Subroutine cbgcmoves_displaysimparams

  !--------------------------------------------------------
  ! Display the cbgcmc sorbate information
  !--------------------------------------------------------
  Subroutine cbgcmoves_display(cbgcmcsorbs, nspc, optunit)
    Type(CBGCMC_Move_Params), Dimension(:), Intent(in) :: cbgcmcsorbs
    Integer, Intent(in)   :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=strLen)   :: molecname
    Character(len=nspc) :: blank
    Integer     :: i, j, unitno, nsorbs, sorbtype 

    blank = Repeat(' ',nspc)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    nsorbs = Size(cbgcmcsorbs, 1)
    Write(unitno,'(2a)') blank, dashedline
    Write(unitno,'(2a)') blank, "CBGCMC Moves Parameters:"

    Do i=1, nsorbs
      sorbtype = cbgcmcsorbs(i)%sorbtype
      molecname = molecules_name(sorbtype)
      Write(unitno, '(a,2x,a)')  blank, dashedline
      Write(unitno, '(a,2x,2a)') blank,"Sorbate Name        : ",  &
          Trim(molecname)

      Call mcmoves_auxrejectdisplay(cbgcmcsorbs(i)%auxparams,nspc+2,unitno)


      Write(unitno,'(2x,2a)') blank,"Thermophysical simulation points:"
      Write(unitno,'(a,4a10)') blank,"Pressure", "Fugacity", "murti", &
          "B(Adams)"
      Do j=1, cbgcmcsorbs(i)%npts
        Write(unitno, '(a,4x,4f10.3)') blank,&
            cbgcmcsorbs(i)%fuglist(j)%pressure, &
            cbgcmcsorbs(i)%fuglist(j)%fugacity, &
            cbgcmcsorbs(i)%fuglist(j)%murti, &
            cbgcmcsorbs(i)%fuglist(j)%B
      End Do
      Write(unitno, '(a,2x,a)')  blank, dashedline
    End Do

    !** Write the move-params
    Write(unitno,'(2a)') blank, "Allowed Move Types and Parameters:"
    Write(unitno,'(2a)') blank,dashedline
    Do i=1, nsorbs
      Write(unitno, '(a,2x,2a)') blank,"Sorbate Name        : ",  &
          Trim(molecname)
      Do j=1,cbgcmcsorbs(i)%no_of_movetypes
        Call mcmoves_display(cbgcmcsorbs(i)%mcmoves(j),nspc+2,unitno)
      end do
      Write(unitno,'(2a)') blank,dashedline
    End Do

  End Subroutine cbgcmoves_display

  !----------------------------------------------------------------------
  ! Calculates the energy of the molecule "molec" with all other molecules
  ! It returns the total energy in "pot" and the sorbate-sorbate 
  ! coulombic and non-coulombic energies in "coulnrg"  and "noncoulnrg"
  ! The return value of the function is true or false depending on whether
  ! the calculation was successful
  !----------------------------------------------------------------------
  Logical Function cbgcmoves_getnrg(sorbates, simcell, nsorbs, sorbtype, &
      molec, natoms, fast, pot, coulnrg, noncoulnrg)
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(in)  :: simcell
    Integer, Intent(in)    :: nsorbs, sorbtype, molec, natoms
    Logical, Intent(in)    :: fast
    Real(kind=RDbl), Intent(out)   :: pot
    Real(kind=RDbl), Dimension(:), Intent(out):: coulnrg
    Real(kind=RDbl), Dimension(:), Intent(out):: noncoulnrg

    Integer    :: i
    Logical    :: mapflag
    Real(kind=RDbl)    :: coulpot, noncoulpot

    !** Set the default value of function to false
    cbgcmoves_getnrg = .False.

    !** Do the energy calculation
    coulnrg = 0.0_RDbl
    noncoulnrg = 0.0_RDbl
    pot = 0.0_RDbl
    Do i=1, nsorbs
      ! Get the non-coulombic energy and store it
      noncoulpot = 0.0_RDbl
      mapflag = ssnoncoul_msdriver(sorbates, simcell, sorbtype, molec, &
          i, fast, noncoulpot, sorbates(sorbtype)%afast(1:natoms, molec), &
          sorbates(i)%afast)
      If (.NOT. mapflag) Return 
      noncoulnrg(i) = noncoulpot

      ! Get the coulombic energy and store it
      coulpot = 0.0_RDbl
      mapflag = sscoul_msdriver(sorbates, simcell, sorbtype, molec,  &
          i, fast, coulpot, sorbates(sorbtype)%afast(1:natoms, molec), &
          sorbates(i)%afast)
      If (.NOT. mapflag) Return
      coulnrg(i) = coulpot

      ! Update the total energy
      pot = pot + noncoulpot + coulpot
    End Do

    !** If it made it this far then the calcn. was successful
    cbgcmoves_getnrg = .True.

  End Function cbgcmoves_getnrg

  !-----------------------------------------------------------------------
  ! Display one line statistics for each species (fugacity, loading)
  ! Requires: gcmcspc -- GCMC species parameters
  !           simno -- simulation number
  !           indent -- no. of spaces from the left margin
  !           unitno -- optional display unit number
  !-----------------------------------------------------------------------
  Subroutine cbgcmoves_dispspcstats(params,simno,indent,optunit)
    Type(CBGCMC_Move_Params), Intent(In)    :: params
    Integer, Intent(In)                   :: simno,indent
    Integer, Optional, Intent(In)         :: optunit

    Integer                      :: unit
    Real(kind=RDbl)              :: loading
    Character(len=indent)        :: blank
    Character(len=strLen)        :: molecname

    unit = 6
    If (Present(optunit)) unit = optunit
    blank = Repeat(' ',indent)    

    molecname = molecules_name(params%sorbtype)
    loading = stats_getcavg(params%nmoles)/(params%nunitcells*one)

    Write(unit,'(3a,2(2x,a,f6.1))') blank,Trim(molecname),':', &
        'pres. (kPa): ',params%fuglist(simno)%pressure, &
        'loading (molec/uc): ',loading

  End Subroutine cbgcmoves_dispspcstats


!!$  !----------------------------------------------------------------------
!!$  ! Get a Random Move Type based on the weights of differenet
!!$  ! move types in cbgcmcsorb, and returns the number of the move and its Tag
!!$  !----------------------------------------------------------------------  
!!$  Subroutine cbgcmoves_getRandomMoveType(cbgcmcsorb,moveno,moveType)
!!$    Type(CBGCMC_Move_Params),Intent(in) :: cbgcmcsorb
!!$    Integer, Intent(out)   :: moveno
!!$    Character(*), Intent(out)   :: moveType
!!$    Real(kind=RDbl) :: prob
!!$    Integer :: i
!!$    prob=rranf()
!!$    moveno=1
!!$    Do i=1,cbgcmcsorb%no_of_movetypes
!!$      If (prob>cbgcmcsorb%mvWeights(i)) Then
!!$        moveno=moveno+1
!!$      Else 
!!$        Exit
!!$      Endif
!!$    End Do
!!$    moveType = Trim(cbgcmcsorb%mvTags(moveno))
!!$  End Subroutine cbgcmoves_getRandomMoveType
End Module cbgcmoves

!!$
!!$
!!$
!!$  !-------------------------------------------
!!$  ! Does the CBGCMC insertion move
!!$  !-------------------------------------------
!!$  Subroutine cbgcmoves_insert(cbgcmcsorb, sorbates, simcell, moveno,simno)
!!$    Type(CBGCMC_Move_Params), Intent(inout):: cbgcmcsorb
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Type(SimCell_Params), Intent(in)  :: simcell
!!$    Integer, Intent(in)    :: moveno,simno
!!$
!!$    Integer                :: index, sorbtype, nmoles, molec, natoms, i
!!$    Integer                :: nsorbs, j, sitetype
!!$    Real(kind=RDbl)        :: uinit, coulpot, noncoulpot, ufinal, utemp
!!$    Real(kind=RDbl)        :: biasfactor, B, boltz, totnrg,change
!!$    Type(GeneralizedCoords), Pointer     :: gcoords
!!$    Type(VecType), Dimension(:), Pointer :: xyzcoords
!!$    Type(VecType)          :: com
!!$    Logical                :: fast, mapflag, nrg_calc_flag
!!$    Real(kind=RDbl), Dimension(MAX_SORBS):: coulnrg
!!$    Real(kind=RDbl), Dimension(MAX_SORBS):: noncoulnrg
!!$    Real(kind=RDbl), Dimension(NO_OF_INTRA_POTS):: potlist_init,potlist_final
!!$    Real(kind=RDbl)        :: intra
!!$    Real(kind=RDbl), Dimension(2,size(sorbates))        :: nrg_arr 
!!$    !** Assumed that all interactions are fast
!!$
!!$    fast = .True.
!!$
!!$    !** Update the move stats
!!$    cbgcmcsorb%insert%att = cbgcmcsorb%insert%att + 1
!!$
!!$    !** Get the current no. of molecules of the sorbate type and other
!!$    !** sundry details
!!$    nsorbs = molecules_getnsorbs()
!!$    sorbtype = cbgcmcsorb%sorbtype
!!$    nmoles = stats_getvalue(cbgcmcsorb%nmoles)
!!$    natoms = config_getnatoms(sorbates, sorbtype)
!!$
!!$    !** Check if we need to increment the size of the various config arrays
!!$    Call config_checkandincr(sorbates, sorbtype, nmoles)
!!$
!!$    !** Check everything is alright with the storage arrays
!!$    Call cbgcmoves_checkIncrStorage(cbgcmcsorb,nmoles)
!!$
!!$    !** Get the molecule no. for insertion and its initial energy
!!$    molec = nmoles + 1
!!$    uinit = 0.0_RDbl
!!$    nrg_calc_flag = .True.
!!$
!!$    !** set the sorbates number 
!!$    Call config_setnmoles(sorbates, sorbtype, molec)
!!$
!!$    !** attempt insertion, calculate the bias factor, energies
!!$    !** ufinal is the energy of the inserted molecule, includes 
!!$    !** noncoul, coul and intra
!!$    !** If nrg_calc_flag=.false. then it means either there was error 
!!$    !** While nrg calculation or that this nrg is very high and will 
!!$    !** have anyway zero acceptance probability 
!!$!!!    Call moves_copyBiasStorage(cbgcmcsorb%mvparams(moveno), &
!!$!!!        cbgcmcsorb%tempb)
!!$    Call moves_insert(sorbates,sorbtype,molec, &
!!$        cbgcmcsorb%mvparams(moveno), biasfactor,nrg_calc_flag, intra, ufinal, &
!!$        nrg_arr)
!!$
!!$    If( .NOT. nrg_calc_flag ) Then 
!!$      Call config_setnmoles(sorbates, sorbtype, molec-1)
!!$      Call stats_update(cbgcmcsorb%nmoles, (molec-1)*1.0_RDbl) 
!!$      Return 
!!$    End If
!!$    cbgcmcsorb%insert%nrg_succ = cbgcmcsorb%insert%nrg_succ + 1    
!!$
!!$    If(Associated(sorbates(sorbtype)%gcoords(molec)%branchedcoords))Then 
!!$      utemp = -cbgcmcsorb%rti * ufinal
!!$      B   = cbgcmcsorb%fuglist(simno)%B
!!$      utemp = utemp + B + Log(biasfactor/Real(molec))
!!$      If (utemp>EXP_MAX) utemp=EXP_MAX
!!$      If (utemp<(-EXP_MAX)) utemp=-EXP_MAX
!!$      boltz = Exp(utemp)/cbgcmcsorb%ins_del_ratio
!!$
!!$      If(boltz > rranf()) Then 
!!$        totnrg=zero
!!$        Do i=1, nsorbs
!!$          ! Get the non-coulombic energy and store it
!!$          !!          mapflag = ssnoncoul_msdriver(sorbates, simcell, sorbtype, molec, &
!!$          !!              i, fast, noncoulpot, sorbates(sorbtype)%afast(1:natoms, molec), &
!!$          !!              sorbates(i)%afast)
!!$          !!          noncoulnrg(i) = noncoulpot
!!$          noncoulnrg(i)=nrg_arr(2,i)
!!$
!!$          coulpot = 0.0_RDbl
!!$          !!          mapflag = sscoul_msdriver(sorbates, simcell, sorbtype, molec, &
!!$          !!              i, fast, coulpot, sorbates(sorbtype)%afast(1:natoms, molec), &
!!$          !!              sorbates(i)%afast)
!!$          !!          coulnrg(i) = coulpot
!!$          coulnrg(i) = nrg_arr(1,i)
!!$          !!          Write(*,*) "Comparing Values from moves and Insert "
!!$          !!          Write(*,*) noncoulpot,nrg_arr(2,i),i,mapflag
!!$          Call molecules_incrcoulnrg(sorbtype, i, coulnrg(i)*caltoj)
!!$          Call molecules_incrnoncoulnrg(sorbtype, i, noncoulnrg(i)*caltoj)
!!$
!!$        End Do
!!$
!!$        Call stats_update(cbgcmcsorb%nmoles, molec*1.0_RDbl)
!!$        cbgcmcsorb%insert%succ = cbgcmcsorb%insert%succ + 1
!!$
!!$        !** Update the total intramolecular energy of sorbtype
!!$        Call molecules_incrIntraNrg(sorbtype,intra*caltoj)
!!$!!!        Call brmoves_copyMolecStorage(cbgcmcsorb%tempb, &
!!$!!!            cbgcmcsorb%insb(molec))
!!$        Return 
!!$      Else    
!!$        Call config_setnmoles(sorbates, sorbtype, molec-1)
!!$        Call stats_update(cbgcmcsorb%nmoles, (molec-1)*1.0_RDbl) 
!!$        Return 
!!$      End If
!!$
!!$    End If
!!$
!!$    Call simcell_pbc(simcell, sorbates(sorbtype)%coords(1:natoms,molec)%rp, &
!!$        sorbates(sorbtype)%coords(1:natoms, molec)%r, &
!!$        sorbates(sorbtype)%coords(1:natoms, molec)%cr)
!!$
!!$
!!$  End Subroutine cbgcmoves_insert
!!$
!!$  !-----------------------------------------
!!$  ! Does the CBGCMC deletion move
!!$  !-----------------------------------------
!!$  Subroutine cbgcmoves_delete(cbgcmcsorb, sorbates, simcell,moveno,simno)
!!$    Type(CBGCMC_Move_Params), Intent(inout) :: cbgcmcsorb
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Type(SimCell_Params), Intent(in)  :: simcell
!!$    Integer, Intent(in)    :: moveno,simno
!!$
!!$    Integer                :: nsorbs, natoms, sorbtype, nmoles, molec, i
!!$    Real(kind=RDbl)        :: biasfactor
!!$
!!$    Logical                :: fast, mapflag
!!$    Real(kind=RDbl)        :: uinit, coulpot, noncoulpot, ufinal, utemp
!!$    Real(kind=RDbl)        :: B, boltz,intra, totnrg
!!$    Real(kind=RDbl), Dimension(MAX_SORBS):: coulnrg
!!$    Real(kind=RDbl), Dimension(MAX_SORBS):: noncoulnrg
!!$    Real(kind=RDbl) ,Dimension(2,Size(sorbates,1)) :: nrg_arr
!!$
!!$    !** Get the sorbate type and other sundry values
!!$    nsorbs = molecules_getnsorbs()
!!$    sorbtype = cbgcmcsorb%sorbtype
!!$    natoms   = config_getnatoms(sorbates, sorbtype)
!!$    nmoles   = stats_getvalue(cbgcmcsorb%nmoles)
!!$
!!$    !** Make sure there are molecules to delete
!!$    If (nmoles == 0) Return
!!$
!!$    !** Increment the no. of deletion attempts
!!$    cbgcmcsorb%delete%att = cbgcmcsorb%delete%att + 1
!!$
!!$    !** Pick a molecule 
!!$    molec = Int(rranf()*nmoles) + 1
!!$
!!$    !** Get the bias weight of the molecule,calls the delete move
!!$    Call moves_delete(sorbates,sorbtype,molec, &
!!$        cbgcmcsorb%mvparams(moveno), biasfactor,intra,uinit,nrg_arr)
!!$
!!$    !** these nrg calcs should always be success
!!$    cbgcmcsorb%delete%nrg_succ = cbgcmcsorb%delete%nrg_succ + 1
!!$
!!$    If(Associated(sorbates(sorbtype)%gcoords(molec)%branchedcoords))Then 
!!$      utemp = cbgcmcsorb%rti * uinit
!!$      B     = cbgcmcsorb%fuglist(simno)%B
!!$      utemp = utemp - B +log(biasfactor*Real(nmoles))
!!$      If (utemp>EXP_MAX) utemp=EXP_MAX
!!$      If (utemp<(-EXP_MAX)) utemp=-EXP_MAX
!!$      boltz = Exp(utemp)*cbgcmcsorb%ins_del_ratio
!!$      If (boltz > rranf()) Then
!!$        nmoles = nmoles - 1
!!$        fast = .True.
!!$        noncoulpot = 0.0_RDbl
!!$        coulpot = 0.0_RDbl
!!$        totnrg=zero
!!$        mapflag = .True.
!!$        !!      Write(*,*) "Attemted deletioin of molecule ", molec
!!$
!!$
!!$        Do i=1, nsorbs
!!$          !!          mapflag = ssnoncoul_msdriver(sorbates, simcell, sorbtype, molec, &
!!$          !!              i, fast, noncoulpot, sorbates(sorbtype)%afast(1:natoms, molec), &
!!$          !!              sorbates(i)%afast)
!!$          !!      Write(*,*) "compare delete nrgs ", nrg_arr(2,i), noncoulpot 
!!$
!!$          !!          noncoulnrg(i) = noncoulpot
!!$          noncoulnrg(i) = nrg_arr(2,i)
!!$!!!          mapflag = sscoul_msdriver(sorbates, simcell, sorbtype, molec, &
!!$!!!              i, fast, coulpot, sorbates(sorbtype)%afast(1:natoms, molec), &
!!$!!!              sorbates(i)%afast)
!!$!!!          coulnrg(i) = coulpot
!!$          coulnrg(i) = nrg_arr(1,i)
!!$          Call molecules_incrcoulnrg(sorbtype, i, -coulnrg(i)*caltoj)
!!$          Call molecules_incrnoncoulnrg(sorbtype, i, -noncoulnrg(i)*caltoj)
!!$        End Do
!!$
!!$        !** Substract the intra nrg of deleted molecule from total
!!$        Call molecules_incrIntraNrg(sorbtype,(-intra*caltoj))
!!$        Call stats_update(cbgcmcsorb%nmoles, nmoles*1.0_RDbl)
!!$        cbgcmcsorb%delete%succ = cbgcmcsorb%delete%succ + 1
!!$        Call config_delmolec(sorbates, sorbtype, molec)
!!$        Return 
!!$      Else   
!!$        Call stats_update(cbgcmcsorb%nmoles, nmoles*1.0_RDbl)
!!$        Return 
!!$      End If
!!$    End If
!!$  End Subroutine cbgcmoves_delete
!!$
!!$  !-------------------------------------------------------------------
!!$  ! Does the translation/rotation/Integration move
!!$  !-------------------------------------------------------------------
!!$  Subroutine cbgcmoves_rot_trans(cbgcmcsorb, sorbates, simcell,moveno)
!!$    Type(CBGCMC_Move_Params), Intent(inout) :: cbgcmcsorb
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Type(SimCell_Params), Intent(in)  :: simcell
!!$    Integer ,Intent(in) :: moveno
!!$
!!$    Integer                :: i, j, nsorbs, sorbtype, nmoles, molec, natoms
!!$    Integer                :: sitetype
!!$    Logical                :: fast, nrgflag
!!$    Real(kind=RDbl)        :: coulpot, noncoulpot,ratio
!!$    Real(kind=RDbl)        :: uinit, ufinal, deltaU, boltz, utemp
!!$    Real(kind=RDbl), Dimension(MAX_SORBS):: oldcoulnrg, newcoulnrg
!!$    Real(kind=RDbl), Dimension(MAX_SORBS):: oldnoncoulnrg, newnoncoulnrg
!!$    Type(VecType)          :: com
!!$
!!$    !** Get the sorbate type and no. of molecules
!!$    sorbtype = cbgcmcsorb%sorbtype
!!$    nsorbs   = molecules_getnsorbs()
!!$    nmoles   = config_getnmoles(sorbates, sorbtype)
!!$    natoms   = config_getnatoms(sorbates, sorbtype)
!!$
!!$    !** Make sure that the no. of molecules is not zero
!!$    If (nmoles == 0) Return
!!$
!!$    !** Update the move stats
!!$    If (cbgcmcsorb%mvTags(moveno)=="TRANSLATE") Then
!!$      cbgcmcsorb%trans%att = cbgcmcsorb%trans%att + 1
!!$    Else If (cbgcmcsorb%mvTags(moveno)=="ROTATE") Then
!!$      cbgcmcsorb%rot%att = cbgcmcsorb%rot%att + 1
!!$    Endif
!!$
!!$    !** Pick a molecule for perturbing
!!$    molec = Int(rranf()*nmoles) + 1
!!$
!!$    !** Get the initial energy
!!$    fast = .True.
!!$    nrgflag = cbgcmoves_getnrg(sorbates, simcell, nsorbs, sorbtype, molec, &
!!$        natoms, fast, uinit, oldcoulnrg, oldnoncoulnrg)
!!$    !** Save the old coordinates
!!$    Call config_copymolec(sorbates(sorbtype), molec, cbgcmcsorb%oldcoords, 1)
!!$
!!$    !** Do the perturbation
!!$    If (cbgcmcsorb%mvTags(moveno)=="TRANSLATE") Then
!!$      Call moves_translate(sorbates,sorbtype,molec,&
!!$          cbgcmcsorb%mvparams(moveno))
!!$    Else If (cbgcmcsorb%mvTags(moveno)=="ROTATE") Then
!!$      Call moves_rotate(sorbates,sorbtype,molec,&
!!$          cbgcmcsorb%mvparams(moveno))
!!$    Else
!!$      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
!!$          " Could not get movetype"
!!$      Stop
!!$    Endif
!!$
!!$    !** Update the other arrays in sorbates
!!$    Call simcell_pbc(simcell, sorbates(sorbtype)%coords(1:natoms,molec)%rp, &
!!$        sorbates(sorbtype)%coords(1:natoms, molec)%r, &
!!$        sorbates(sorbtype)%coords(1:natoms, molec)%cr)
!!$
!!$    !** Check the siting if the smap is not "null"
!!$    If (Associated(cbgcmcsorb%smap)) Then
!!$      com = gcmodels_getcom(sorbates(sorbtype)%gcoords(molec))
!!$      sitetype = smap_getPosSiteType(cbgcmcsorb%smap, com)
!!$      Do i=1, cbgcmcsorb%nexclsites
!!$        If (sitetype == cbgcmcsorb%exclsitetype(i)) Then
!!$          ! Restore the original coordinates
!!$          Call config_copymolec(cbgcmcsorb%oldcoords, 1, &
!!$              sorbates(sorbtype), molec)
!!$
!!$          ! Update the pairwise energies irrespective of whether we
!!$          ! accepted the move or not
!!$          Do j=1, nsorbs
!!$            Call molecules_incrcoulnrg(sorbtype, j, 0.0_Rdbl)
!!$            Call molecules_incrnoncoulnrg(sorbtype, j, 0.0_Rdbl)
!!$          End Do
!!$          Call molecules_incrIntraNrg(sorbtype,zero)
!!$          Return
!!$        End If
!!$      End Do
!!$    End If
!!$
!!$    !** Get the final energy
!!$    nrgflag = cbgcmoves_getnrg(sorbates, simcell, nsorbs, sorbtype, molec, &
!!$        natoms, fast, ufinal, newcoulnrg, newnoncoulnrg)
!!$    !** Check if the energy calc. was successful.  If not
!!$    !** get out of here after resetting the molecule
!!$    If (.Not. (nrgflag)) Then
!!$      Call config_copymolec(cbgcmcsorb%oldcoords, 1, sorbates(sorbtype), molec)
!!$      Return
!!$    End If
!!$
!!$    !** Accept/Reject
!!$    deltaU = ufinal - uinit
!!$    utemp  = -cbgcmcsorb%rti * deltaU
!!$    If (utemp > 30) utemp = 30
!!$    If (utemp < -30) utemp = -30
!!$    boltz  = Exp(utemp)
!!$
!!$    If (boltz > rranf()) Then
!!$      ! Update the pairwise energies irrespective of whether we
!!$      ! accepted the move or not
!!$      Do i=1, nsorbs
!!$        Call molecules_incrcoulnrg(sorbtype, i,    &
!!$            (newcoulnrg(i)-oldcoulnrg(i))*caltoj)
!!$        Call molecules_incrnoncoulnrg(sorbtype, i, &
!!$            (newnoncoulnrg(i)-oldnoncoulnrg(i))*caltoj)
!!$      End Do
!!$
!!$      ! Update the move stats
!!$      If (cbgcmcsorb%mvTags(moveno)=="TRANSLATE") Then
!!$        cbgcmcsorb%trans%succ = cbgcmcsorb%trans%succ + 1
!!$      Else If (cbgcmcsorb%mvTags(moveno)=="ROTATE") Then
!!$        cbgcmcsorb%rot%succ = cbgcmcsorb%rot%succ + 1
!!$      Endif
!!$    Else
!!$      ! Restore the original coordinates
!!$      Call config_copymolec(cbgcmcsorb%oldcoords, 1, sorbates(sorbtype), molec)
!!$
!!$      ! Update the pairwise energies irrespective of whether we
!!$      ! accepted the move or not
!!$      Do i=1, nsorbs
!!$        Call molecules_incrcoulnrg(sorbtype, i, 0.0_Rdbl)
!!$        Call molecules_incrnoncoulnrg(sorbtype, i, 0.0_Rdbl)
!!$      End Do
!!$    End If
!!$
!!$    !** Adjust the jump length based on the acceptance ratio
!!$    If (cbgcmcsorb%mvTags(moveno)=="TRANSLATE") Then
!!$      ratio= cbgcmcsorb%trans%succ / ( cbgcmcsorb%trans%att*one )
!!$      Call moves_adjustParams(cbgcmcsorb%mvparams(moveno), ratio)
!!$    Else If (cbgcmcsorb%mvTags(moveno)=="ROTATE") Then
!!$      ratio= cbgcmcsorb%rot%succ / ( cbgcmcsorb%rot%att*one )
!!$      Call moves_adjustParams(cbgcmcsorb%mvparams(moveno), ratio)
!!$    Endif
!!$
!!$  End Subroutine cbgcmoves_rot_trans
!!$
!!$  !-------------------------------------------------------------------
!!$  ! Transforms a molecule into a different shape ( cut and regrow in case of
!!$  ! linear alkanes )
!!$  !-------------------------------------------------------------------
!!$  Subroutine cbgcmoves_transform(cbgcmcsorb, sorbates, simcell,moveno)
!!$    Type(CBGCMC_Move_Params), Intent(inout) :: cbgcmcsorb
!!$    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
!!$    Type(SimCell_Params), Intent(in)  :: simcell
!!$    Integer ,Intent(in) :: moveno
!!$
!!$    Integer                :: i, nsorbs, sorbtype, nmoles, molec, natoms
!!$    Integer                :: sorb,atomno,sitetype, siteno
!!$    Logical                :: fast, nrgflag, nrg_calc_flag
!!$    Real(kind=RDbl)        :: coulpot, noncoulpot,ratio
!!$    Real(kind=RDbl)        :: uinit, ufinal, deltaU, boltz, utemp, biasfactor
!!$    Real(kind=RDbl),Dimension(2,NO_OF_INTRA_POTS)         :: intra
!!$    Real(kind=RDbl), Dimension(2,2,Size(sorbates,1))      :: nrg_arr
!!$    Real(kind=RDbl), Dimension(MAX_SORBS):: oldcoul, newcoul
!!$    Real(kind=RDbl), Dimension(MAX_SORBS):: oldnoncoul, newnoncoul
!!$    Type(VecType)          :: com
!!$
!!$    !** Get the sorbate type and no. of molecules
!!$    sorbtype = cbgcmcsorb%sorbtype
!!$    nsorbs   = molecules_getnsorbs()
!!$    nmoles   = config_getnmoles(sorbates, sorbtype)
!!$    natoms   = config_getnatoms(sorbates, sorbtype)
!!$    fast=.True.
!!$
!!$    ! Make sure that 1) the no. of molecules is not zero, 
!!$    !                2) there are atleast two atoms in the molecule
!!$    If (nmoles == 0) Return
!!$    If (natoms<2) Return
!!$
!!$    ! Update the move stats
!!$    If (cbgcmcsorb%mvTags(moveno)=="CUT_REGROW") Then
!!$      cbgcmcsorb%transform%att = cbgcmcsorb%transform%att + 1
!!$    Else 
!!$      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$      stop
!!$    Endif
!!$
!!$    ! Pick a molecule and node for perturbing
!!$    molec = Int(rranf()*nmoles) + 1
!!$
!!$    ! the atomno represents the atom that will be left unperturbed
!!$    ! the chain will be regrown by selecting new positions for 
!!$    ! (atomno+1) to (natoms)
!!$    atomno = Int(rranf()*(natoms-1)) + 1
!!$
!!$    !** Get the initial energy, note that this uinit does not have 
!!$    ! the intramolecular energy of the molecule that will have to be 
!!$    ! obtained during the transform move
!!$
!!$    fast = .True.
!!$!!!    nrgflag = cbgcmoves_getnrg(sorbates, simcell, nsorbs, sorbtype, molec, &
!!$!!!        natoms, fast, uinit, oldcoul, oldnoncoul)
!!$!!!    
!!$!!!    If (.Not.nrgflag) Then
!!$!!!      Write(*,*) "Problem with energy of existing molecule"
!!$!!!      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$!!!      Stop
!!$!!!    Endif
!!$
!!$    !** Save the old coordinates
!!$    Call config_copymolec(sorbates(sorbtype), molec, cbgcmcsorb%oldcoords, 1)
!!$
!!$    !** Do the perturbation. Note that here intra is an array, 
!!$    !** intra(1)=initial nrg, intra(2)=nrg of New configuration
!!$    If (cbgcmcsorb%mvTags(moveno)=="CUT_REGROW") Then
!!$!!!     Write(*,*) "Attempting tranform at molec, atomno: ", molec, atomno
!!$      Call moves_transform( sorbates, sorbtype, molec, atomno, &
!!$          cbgcmcsorb%mvparams(moveno), biasfactor, nrg_calc_flag, ufinal, &
!!$          intra, nrg_arr)
!!$    Else
!!$      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
!!$          " Could not get movetype"
!!$      Stop
!!$    Endif
!!$
!!$
!!$    !** Check the siting if the smap is not "null"
!!$    If (Associated(cbgcmcsorb%smap)) Then
!!$      Write(*,*) "This might not be very appropriate for long chain alkanes"
!!$      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$      Stop
!!$      com = gcmodels_getcom(sorbates(sorbtype)%gcoords(molec))
!!$      sitetype = smap_getPosSiteType(cbgcmcsorb%smap, com)
!!$      Do siteno=1, cbgcmcsorb%nexclsites
!!$        If (sitetype == cbgcmcsorb%exclsitetype(siteno)) Then
!!$          ! Restore the original coordinates
!!$          Call config_copymolec(cbgcmcsorb%oldcoords, 1, &
!!$              sorbates(sorbtype), molec)
!!$
!!$          ! Update the pairwise energies irrespective of whether we
!!$          ! accepted the move or not
!!$          Do sorb=1, nsorbs
!!$            Call molecules_incrcoulnrg(sorbtype, sorb, 0.0_Rdbl)
!!$            Call molecules_incrnoncoulnrg(sorbtype, sorb, 0.0_Rdbl)
!!$          End Do
!!$          Return
!!$        End If
!!$      End Do
!!$    End If
!!$
!!$    !** Check if the energy calc. was successful.  If not
!!$    !** get out of here after resetting the molecule
!!$    If (nrg_Calc_flag) Then
!!$      cbgcmcsorb%transform%nrg_succ = cbgcmcsorb%transform%nrg_succ + 1
!!$    Else
!!$      Call config_copymolec(cbgcmcsorb%oldcoords, 1, sorbates(sorbtype), molec)
!!$      Return
!!$    End If
!!$
!!$    uinit=zero
!!$    ufinal=zero
!!$    Do i=1,nsorbs
!!$      uinit=uinit+nrg_arr(1,1,i)+nrg_arr(1,2,i)
!!$      ufinal=ufinal+nrg_arr(2,1,i)+nrg_arr(2,2,i)
!!$    End Do
!!$    uinit=uinit+intra(1,TOTAL_INDEX)
!!$    ufinal=ufinal+intra(2,TOTAL_INDEX)
!!$    !!    Write(*,'(a,2f16.5)') " Checking again :Initial, Final",uinit,ufinal
!!$    !!    Write(*,'(a,e16.6)') "Bias is ", biasfactor
!!$
!!$    !** Accept/Reject
!!$    deltaU = ufinal - uinit
!!$    utemp  = -cbgcmcsorb%rti * deltaU
!!$
!!$    If (utemp > EXP_MAX) utemp = EXP_MAX
!!$    If (utemp < -EXP_MAX) utemp = -EXP_MAX
!!$    boltz  = Exp(utemp)*biasfactor
!!$
!!$    If (boltz > rranf()) Then
!!$      ! add new energies
!!$      Do sorb=1, nsorbs
!!$        Call molecules_incrcoulnrg(sorbtype, sorb, &
!!$            (nrg_arr(2,1,sorb)-nrg_arr(1,1,sorb))*caltoj)
!!$        Call molecules_incrnoncoulnrg(sorbtype, sorb, &
!!$            (nrg_arr(2,2,sorb)-nrg_arr(1,2,sorb))*caltoj)
!!$      End Do
!!$
!!$      !** Register the change in total intra energy
!!$      Call molecules_incrIntraNrg(sorbtype, &
!!$          ( intra(2,TOTAL_INDEX)-intra(1,TOTAL_INDEX) )*caltoj )
!!$      ! Update the move stats
!!$      cbgcmcsorb%transform%succ = cbgcmcsorb%transform%succ + 1
!!$    Else
!!$      ! Restore the original coordinates
!!$      Call config_copymolec(cbgcmcsorb%oldcoords, 1, sorbates(sorbtype), molec)
!!$
!!$      ! update with old energies, for ensemble average
!!$      Do sorb=1, nsorbs
!!$        Call molecules_incrcoulnrg(sorbtype, sorb, 0.0_Rdbl)
!!$        Call molecules_incrnoncoulnrg(sorbtype, sorb, 0.0_Rdbl)
!!$      End Do
!!$      Call molecules_incrIntraNrg(sorbtype,zero)
!!$
!!$    End If
!!$
!!$  End Subroutine cbgcmoves_transform
!!$
