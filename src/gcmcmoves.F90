!----------------------------------------------------------------------------
! This module contains the different data structures and routines for 
! handling the GCMC move types.  The data type "GCMC_Move_Params" contains
! the simulation information for a single species type.  Move-specific
! information and handling are done using the "moves" module through the
! "mcmoves" module.
!
! Important routines are:
!    gcmcmoves_init -- handles initilization for each species
!    gcmcmoves_move -- picks a movetype and executes
!    gcmcmoves_displaystats -- displays intermediate statistics for species
!----------------------------------------------------------------------------

Module gcmcmoves
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl,strLen,MAX_EXCLSITES,zero,dashedline,dashedline2, &
      kcalmole_kb,one,lstrlen,dbgflag
  Use utils, Only: filesrchstr,toupper,split,stripcmnt,toreal,toint, &
      allocErrDisplay,readblank,real2str,deallocErrDisplay, cleanstring, int2str
  Use file, Only: file_open
  Use general, Only: genparams
  Use config, Only: AtMolcoords,config_getnmoles,config_config2xyz, &
      config_changeflex,config_dumpmol
  Use simcell, Only: SimCell_Params,simcell_getnuc
  Use molecules, Only: molecules_name,molecules_gettype, & 
      molecules_getgcmodeltype,molecules_getnatoms,molecules_changedof, &
      molecules_getnsorbs
  Use mcmoveset, Only: Move_Set,mcmoveset_initms,mcmoveset_readmswts, &
      mcmoveset_pickmove,mcmoveset_displayms,mcmoveset_setmstag, &
      mcmoveset_getmswt, mcmoveset_setInsDelratio, mcmoveset_insdelratio
  Use mcmoves, Only: MC_Move_Params,AuxReject_Params,mcmoves_init, &
      mcmoves_initaux,mcmoves_perturb, mcmoves_insert,mcmoves_delete, &
      mcmoves_getbasictag,mcmoves_display,mcmoves_auxrejectdisplay, &
      mcmoves_stats, mcmoves_resetstats,mcmoves_checkOppMoves, &
      mcmoves_postadjust, mcmoves_idchange
  Use random, Only: rranf,random_getnewseed
  Use stats, Only: Statistics,stats_setvalue,stats_getcavg,stats_update,&
      stats_display,stats_init, stats_reset
  Use subinteract, Only: Subset_Interactions, subinteract_chkmolnums
  Use storetop, Only: storetop_display
      
  Implicit None
  Save

  Private
  Public :: GCMC_Move_Params,Fugacity_Params,gcmcmoves_init,gcmcmoves_move, &
      gcmcmoves_display,gcmcmoves_displaystats,gcmcmoves_displaysimparams, &
      gcmcmoves_dispspcstats,gcmcmoves_updatenmolecs, gcmcmoves_beginsim, &
      gcmcmoves_settemplist,  gcmcmoves_copyfuglist

  Type Fugacity_Params
    Real(kind=RDbl)   :: pressure
    Real(kind=RDbl)   :: murti       ! Chemical Potential/RT
    Real(kind=RDbl)   :: fugacity
    Real(kind=RDbl)   :: B           ! Adams Notation
  End Type Fugacity_Params

  !** All move parameters necessary for each species
  Type GCMC_Move_Params
    Integer                   :: spc  !* normal index of species 
    Integer                   :: npts
    Integer                   :: nunitcells
    Type(Statistics)          :: nmoles
    Integer                   :: blocksize !* Block size for collecting stats

    !** Thermophysical Properties
    Real(kind=RDbl)           :: tk
    Real(kind=RDbl)           :: rti       !* 1.0/(Rgas*tk)
    Character(len=strLen)     :: pressurefile
    Type(Fugacity_Params), Dimension(:), Pointer   :: fuglist
    Type(Fugacity_Params), Dimension(:), Pointer   :: single_sim_fug_arr
    Real(kind=RDbl), Dimension(:), Pointer         :: single_sim_B_arr

    !** Move information: tags can be "INSERT","DELETE","ROTATE" or "TRANSLATE"
    Integer                                   :: no_of_movetypes
    Type(Move_Set)                            :: moveset
    Type(MC_Move_Params),Dimension(:),Pointer :: mcmoves
    Type(AuxReject_Params), Pointer           :: auxparams

  End Type GCMC_Move_Params

  ! ** maximum allowed pressure points for gcmc simulations
  Integer, Parameter :: MAX_PRESSURE_POINTS=50
  
Contains
  !----------------------------------------------------------------
  ! Initializes the various GCMC parameters from the control file
  ! Requires: gcmcspc -- parameters for individual species
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           ctrl_filename -- the name of the control file
  !----------------------------------------------------------------
  Subroutine gcmcmoves_init(gcmcspc,species,simcell,ctrl_filename, auxmv)
    Type(GCMC_Move_Params), Dimension(:), Intent(InOut) :: gcmcspc
    Type(AtMolCoords), Dimension(:), Intent(InOut)      :: species
    Type(SimCell_Params), Intent(In)                    :: simcell
    Character(*), Intent(In)                            :: ctrl_filename
    Type( AuxMoveObjects),Pointer :: auxmv

    Integer                       :: unitno,i,j,m,error,blocksize,nfields
    Integer                       :: spc,nspc,npts,nmoles,no_of_movetypes
    Integer                       :: size_spc
    Real(kind=RDbl)               :: tk,inswt,delwt
    Character(len=strLen)         :: spcname,gcmodeltype
    Character(len=MAX_PRESSURE_POINTS*strLen)         :: pressureline
    Character(len=strLen), Dimension(MAX_EXCLSITES) :: fields

    !** from the control file read no of species
    unitno = file_open(ctrl_filename,110)
    nspc = Size(gcmcspc, 1)
    size_spc=molecules_getnsorbs()
    !** Initialize the parameters for each adspecies
    Do i = 1,nspc
      gcmcspc(i)%nunitcells = simcell_getnuc(simcell)

      Nullify(gcmcspc(i)%auxparams)

      Read(unitno, *) spcname
      spcname = Trim(stripcmnt(spcname))
      Write(*,'(2a)') 'Initializing GCMC moves for: ',spcname
      spc = molecules_gettype(spcname)
      If (spc == 0) Then
        Write(0,'(1x,2a,i4,1x,2a)') __FILE__," : ",__LINE__, &
            Trim(spcname), " is not one of the species types in molecule list"
        Stop
      End If

      gcmcspc(i)%spc = spc

      !** Initialize the pressure, fugacity etc.
      npts = gcmcspc(i)%npts
      Allocate(gcmcspc(i)%fuglist(npts), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'fuglist')
      Allocate(gcmcspc(i)%single_sim_fug_arr(size_spc), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'fuglist')
      Allocate(gcmcspc(i)%single_sim_B_arr(size_spc), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'fuglist')


!      Write(*,*) "Reading the pressure line from ctrlfile"
      !** Different pressure values in the isotherm
      Read(unitno,'(a)') pressureline
      Call gcmcmoves_getpressures(gcmcspc(i),pressureline,unitno)      

      !** Read the auxilliary parameters for rejections
      Call mcmoves_initaux(gcmcspc(i)%auxparams,ctrl_filename,simcell)

      gcmodeltype = molecules_getgcmodeltype(spc)
      tk  = gcmcspc(i)%tk

      !** Read the line containing number of movetypes
      Read(unitno,*) no_of_movetypes
      gcmcspc(i)%no_of_movetypes = no_of_movetypes

!** IS THIS NECESSARY? It breaks hgcmc
!!$ I think dof should be explicitly specified in molecules, or let the 
!!$ program calculate it automatically. If dof is not specified in molecule 
!!$ file the we assume interally_flexible. GCMC should not assume dof=6. 
!!$ That will break linsert move
!!$      !** with only 4 move types, assume we're without TRANSFORM, change DOF
!!$      If (no_of_movetypes == 4) Then
!!$        Call molecules_changedof(spc,6,'MODIFIED_BY_GCMC')
!!$        Call config_changeflex(species,spc,.False.)
!!$      End If

      !** allocate memory 
      Allocate(gcmcspc(i)%mcmoves(no_of_movetypes), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'mcmoves')

      Call mcmoveset_readmswts(gcmcspc(i)%moveset,ctrl_filename, &
          no_of_movetypes)

      !** Read and initialize the different movetypes, set the tags
      Do j = 1,no_of_movetypes
        Call mcmoves_init(gcmcspc(i)%mcmoves(j),ctrl_filename, &
            simcell, species, "MuVT", auxmv, spc, gcmcspc(i)%auxparams)
        Call mcmoveset_setmstag(gcmcspc(i)%moveset,j, &
            mcmoves_getbasictag(gcmcspc(i)%mcmoves(j)))
      End Do

!LC      Call mcmoveset_displayms(gcmcspc(i)%moveset,0,6)

      !** Check the input move types 
      inswt = mcmoveset_getmswt(gcmcspc(i)%moveset,'INSERT') 
      If (inswt <= zero) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Move set for ',Trim(spcname),' must contain an INSERT move'
!        Stop
      End If

      delwt = mcmoveset_getmswt(gcmcspc(i)%moveset,'DELETE') 
      If (delwt <= zero) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Move set for ',Trim(spcname),' must contain a DELETE move'
!        Stop
      End If

      !** make sure that we have a transform move if movetypes > 4
      If (no_of_movetypes > 4) Then
        If (mcmoveset_getmswt(gcmcspc(i)%moveset,'TRANSFORM') <= zero) Then
          Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
              ' Move set for ',Trim(spcname),' should usually contain a TRANSFORM move'
!          Stop
        End If
      End If

      If (Abs((inswt/delwt)-one) > 1.0e-6) Then
        Write(*,'(1x,3a,i4)') "WARNING : ",__FILE__," : ",__LINE__
        Write(0,*) "WARNING : different insertion / deletion rates "
        Write(0,*) "WARNING : Make sure that this is correct for your Case "
      End If

      !** Ratio of insertio attempts/ delete attempts
      Call mcmoveset_setInsDelRatio(gcmcspc(i)%moveset)

#ifdef FULLCHECK
      If (mcmoveset_getmswt(gcmcspc(i)%moveset,'TRANSLATE') <= zero) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Move set for ',Trim(spcname),' must contain a TRANSLATE move'
        Stop
      End If

      If (molecules_getnatoms(spc) > 1) Then
        If (mcmoveset_getmswt(gcmcspc(i)%moveset,'ROTATE') <= zero) Then
          Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
              ' Move set for multiatomic ',Trim(spcname), &
              ' must contain a ROTATE move'
          Stop
        End If
      End If
#endif

      !** Calculate Beta (1/(Rgas*tk))
      gcmcspc(i)%rti = 1.0_RDbl/(kcalmole_kb*tk)

      !** Initialize the statistics variables
      blocksize = gcmcspc(i)%blocksize
      nmoles = config_getnmoles(species, spc)
      Call stats_init(gcmcspc(i)%nmoles, &
          "Number of Molecules (ins, blk, cum, stdd)",blocksize,.False.,"f7.2")
      Call stats_setvalue(gcmcspc(i)%nmoles, nmoles*1.0_RDbl)
      
      !** Read the blank separating line
      If (i /= nspc) Then 
        Call readblank(unitno,__FILE__,__LINE__)
      End If
      
    End Do
    
    !** Make sure that opposite moves are properly initialized
    !** Some of these might lead to changing the current-position 
    !** in the control file. So this part should not be taken inside 
    !** the above Do-loop
    Do i = 1,Size(gcmcspc)
      Call mcmoves_checkOppMoves(gcmcspc(i)%mcmoves,species, auxmv, &
          ctrl_filename)
    End Do

  End Subroutine gcmcmoves_init

  !--------------------------------------------------------------------
  ! Pick a move and make it using the mcmoves module.  Returns true if
  ! move was successful.
  ! Requires:  gcmcspc -- parameters for individual species
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            simno -- the simulation number
  !--------------------------------------------------------------------
  Logical Function gcmcmoves_move(gcmcspc,subints,species,simcell,simno)
    Type(GCMC_Move_Params), Intent(InOut)                  :: gcmcspc
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Integer, Intent(In)                                    :: simno

    Integer                   :: i,nmoles,spc,moveno,molec,dumpunit
    Real(kind=RDbl)           :: insdelratio, B_term
    Character(len=strLen)     :: movename
    Character(len=lstrLen)    :: string,filename
    Integer, Dimension(3)     :: subset
    Logical                   :: success, knrgflag

    !** Pick the move type 
    spc = gcmcspc%spc
    nmoles = config_getnmoles(species,gcmcspc%spc)
    moveno = mcmoveset_pickmove(gcmcspc%moveset,movename)
    success = .False.

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call subinteract_chkmolnums(subints(1),species,2,6)
    Write(*,'(a,i5)') 'GCMC iteration: ',genparams%currentiteration
    Write(*,'(2a)') 'GCMC feedback: selected move = ',Trim(movename)
    string = random_getnewseed()
    Write(*,'(2a)') 'GCMC feedback: rranf iseed =   ',Trim(string)
    Write(*,'(a,i3)') 'GCMC feedback: spc = ',spc
    molec = 0
#endif

    Select Case(Trim(movename))
    Case("INSERT")
      insdelratio = mcmoveset_insdelratio(gcmcspc%moveset)
      success = mcmoves_insert(gcmcspc%mcmoves(moveno),molec, &
          gcmcspc%rti,gcmcspc%fuglist(simno)%B,subints,species,simcell,&
          insdelratio)

    Case("DELETE")
      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
        molec = Int(rranf()*nmoles) + 1
        insdelratio = mcmoveset_insdelratio(gcmcspc%moveset)
        success = mcmoves_delete(gcmcspc%mcmoves(moveno),molec, &
            gcmcspc%rti,gcmcspc%fuglist(simno)%B,subints,species, &
            simcell,insdelratio)
      End If

    Case("TRANSLATE", "ROTATE", "TRANSFORM")
      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
        molec = Int(rranf()*nmoles) + 1
        knrgflag = .False.  !** these moves dont change knrg
        success = mcmoves_perturb(gcmcspc%mcmoves(moveno),&
            (/spc,molec,0/), gcmcspc%rti,subints,species,simcell,knrgflag)
      End If

    Case("IDCHANGE")
      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
        molec = Int(rranf()*nmoles) + 1
        success = mcmoves_idchange(gcmcspc%mcmoves(moveno),&
            (/spc,molec,0/), gcmcspc%rti, subints, species, simcell, &
            gcmcspc%single_sim_B_arr)
      End If

    Case Default
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Could not get movename"
      Stop
    End Select
    
    gcmcmoves_move = success
    
    !** Note : check whether "molec" is defined carefully above
    Call mcmoves_postadjust(gcmcspc%mcmoves(moveno),species,spc,molec,success) 

#ifdef DEBUG  !** very useful debugging feedback
    Write(*,'(2a,t34,a,2i3,t55,a,l4)') 'GCMC feedback: move = ', &
        Trim(movename),'   spc,mol =',spc,molec,'    Accepted? ',gcmcmoves_move
    Write(*,*) 'number of molecules: ',config_getnmoles(species,spc)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Call subinteract_chkmolnums(subints(1),species,2,6)
    Write(*,*)
#endif
#ifdef EXTRA
    string = int2str(genparams%currentiteration)
    Write(filename,'(2a)') 'config',Trim(string)
    dumpunit =  file_open(filename)
    Call config_dumpmol(species,2,1,2,dumpunit)
    Close(unit=dumpunit)
#endif


  End Function gcmcmoves_move

  !-----------------------------------------------------------------
  ! Update the number of moles for all species
  ! Requires:  gcmcspc -- array of GCMC species parameters
  !            species -- species data structure
  !-----------------------------------------------------------------
  Subroutine gcmcmoves_updatenmolecs(gcmcspc,species)
    Type(GCMC_Move_Params), Dimension(:), Intent(InOut) :: gcmcspc
    Type(AtMolCoords), Dimension(:), Intent(In)         :: species

    Integer            :: spc,nmoles

    Do spc = 1,Size(gcmcspc,1)
      nmoles = config_getnmoles(species,gcmcspc(spc)%spc)
      Call stats_update(gcmcspc(spc)%nmoles,nmoles*one)      
    End Do

  End Subroutine gcmcmoves_updatenmolecs
  
  !------------------------------------------------------
  ! Set the simulation temperature for each species type
  !------------------------------------------------------
  Subroutine gcmcmoves_settemp(gcmcspc, tk)
    Type(GCMC_Move_Params), Intent(inout)  :: gcmcspc
    Real(kind=RDbl), Intent(in)   :: tk
    gcmcspc%tk = tk
    Return
  End Subroutine gcmcmoves_settemp
  
  !------------------------------------------------------
  ! Set the number of simulation points
  !------------------------------------------------------
  Subroutine gcmcmoves_setnpts(gcmcspc, npts)
    Type(GCMC_Move_Params), Intent(inout)  :: gcmcspc
    Integer, Intent(in)   :: npts

    gcmcspc%npts = npts
    Return
  End Subroutine gcmcmoves_setnpts



  !------------------------------------------------------------------------
  ! Fill in the rest of the fields in the temperature list
  ! Requires: gcmcparams -- GCMC simulation parameters
  !           templist -- volume of the simulation cell
  ! NOTE: this routine does not work the way it should 
  !  just copies first field of templist for now
  ! this should replace _settemp
  !------------------------------------------------------------------------
  Subroutine gcmcmoves_settemplist(gcmcparams, templist)
    Type(GCMC_Move_Params), Dimension(:),Intent(InOut)  :: gcmcparams
    Real(kind=RDbl), Dimension(:)     :: templist 
    Integer           :: spc
    Do spc=1,Size(gcmcparams)
      gcmcparams(spc)%tk=templist(1)
    End Do
  End Subroutine gcmcmoves_settemplist

  !---------------------------------------------------------
  ! Parses the fugacity line to set the fugacities of the
  ! species type.  It can either do the fugacities at equal
  ! intervals or read them from a file.
  ! Works only with pressureline<strLen
  !---------------------------------------------------------
  Subroutine gcmcmoves_getpressures(gcmcspc, pressureline, unitno)
    Type(GCMC_Move_Params), Intent(inout) :: gcmcspc
    Character(*), Intent(in)              :: pressureline
    Integer                   :: unitno

    Integer                   :: i,npts,error,convfields,maxnfields
    Integer                   :: nfields1,nfields2,nfields3, count,spc
    Real(kind=RDbl)           :: startp, endp
    Character(len=strLen)     :: filename
    Character(len=2*strLen)   :: str1,str2,str3
    Character(len=strLen), Dimension(:), Pointer :: strfields1 ,strfields2, &
        strfields3 

    spc = gcmcspc%spc

    !** Make sure there is enough memory before calling split
    maxnfields = Max(strLen, gcmcspc%npts)

    Allocate(strfields1(maxnfields),strfields2(maxnfields),&
        strfields3(maxnfields),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__) 

    !** Split the fields by commas
    str1=cleanstring(pressureline)
    nfields1 = split(str1, strfields1, ",") 
    nfields2=0
    nfields3=0   

    ! only first strLen characters of each line are read, so read more lines
    ! if necessary, maximum 3 such lines allowed, 
    If (Adjustl(Trim(strfields1(nfields1)))=="+") Then
      nfields1=nfields1-1 ! exclude the "+" sign
      If (nfields1>=(gcmcspc%npts)) Then
        Write(*,*) "remove + sign from pressure line; wrong pressureline"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      Read(unitno,'(a)') str2
      str2=cleanstring(str2)
      nfields2 = split(str2, strfields2, ",") 
      If (Adjustl(Trim(strfields2(nfields2)))=="+") Then
        ! note : nfields = npts +1 ( to include "+" sign at end)
        nfields2=nfields2-1

        If ((nfields2+nfields1)>=(gcmcspc%npts)) Then
          Write(*,*) "wrong 2nd pressureline"
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
        Read(unitno,'(a)') str3
        str3=cleanstring(str3)
        ! note : last line, so no "+" here
        nfields3 = split(str3, strfields3, ",") 
      Endif
    Endif

!    Write(*,*) nfields1,  nfields2,  nfields3

    If ((gcmcspc%npts==1) .And. (nfields1>1)) Then
      Write(*,*) "Number of points is only one. Still more than one pressure"
      Write(*,*) "specifeid in ctrlfile. This could be erraneous"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    !** If nfields is 1, it's possible that it is not a file but
    !** the fugacity we want IF npts is 1.
    If ((nfields1 == 1).And.(gcmcspc%npts == 1)) Then

      !** Check if the line is an integer
      startp = toreal(strfields1(1),error)
      If (error /= 0) Then

        !** It's a file. Read it.
        nfields1 = split(pressureline, strfields1, " ")
        gcmcspc%pressurefile = strfields1(1)
        Call gcmcmoves_readpressures(gcmcspc)

      Else
        !** It's a number. Store it.
        gcmcspc%fuglist(1)%pressure = startp

      End If

    Else If (nfields1 == 1) Then
      !** We need to read the points from a file, Get the filename
      nfields1 = split(pressureline, strfields1, " ")
      gcmcspc%pressurefile = strfields1(1)
      Call gcmcmoves_readpressures(gcmcspc)

    Else If (nfields1==2) Then
      !** We need to generate the range
      startp = toreal(strfields1(1))

      !** Get the last field after separating it from the rest of the line
      str1 = strfields1(2)
      nfields1 = split(str1, strfields1, " ")
      endp = toreal(strfields1(1))

      !** Fill in the gcmcspc pressure data
      Call gcmcmoves_genpressures(gcmcspc, startp, endp)

    Elseif (nfields1==gcmcspc%npts) Then
      !** we have all pressure points listed there
      Do i=1,gcmcspc%npts
        gcmcspc%fuglist(i)%pressure=toreal(Trim(strfields1(i)))
      Enddo
    Elseif ((nfields1+nfields2+nfields3)==gcmcspc%npts) Then
      count=0
      If (nfields1>0) Then
        !** we have all pressure points listed there
        Do i=1,nfields1
          count=count+1
          gcmcspc%fuglist(count)%pressure=toreal(Trim(strfields1(i)))
        Enddo
      Endif
      If (nfields2>0) Then
        Do i=1,nfields2
          count=count+1
          gcmcspc%fuglist(count)%pressure=toreal(Trim(strfields2(i)))
        Enddo
      Endif
      If (nfields3>0) Then
        Do i=1,nfields3
          count=count+1
          gcmcspc%fuglist(count)%pressure=toreal(Trim(strfields3(i)))
        Enddo
      Endif
    Else 
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "could not parse the pressureline, try using a pressurefile?"
      Write(*,*) "pressureline passed : ",  str1
      Write(*,*) "nfields, gcmcspc%npts",nfields1, gcmcspc%npts
      Stop
    End If

    !** Free some paltry memory
    DeAllocate(strfields1,strfields2,strfields3,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__) 

  End Subroutine gcmcmoves_getpressures

  !---------------------------------------------------------
  ! Read the fugacities from the file "filename"
  !---------------------------------------------------------
  Subroutine gcmcmoves_readpressures(gcmcspc)
    Type(GCMC_Move_Params), Intent(InOut) :: gcmcspc

    Integer                :: unitno, spc, lineno, npts, i
    Character(len=strLen)  :: line, spcname, filename

    !** Open the file if not already open
    filename = gcmcspc%pressurefile
    unitno = file_open(filename,110)

    !** Find the species name in the file
    spc = gcmcspc%spc
    spcname = molecules_name(spc)
    lineno = filesrchstr(unitno, Trim(spcname), line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4, 4a)') __FILE__," : ",__LINE__, &
          " Could not find the species ", Trim(spcname), " in the file ",&
          Trim(filename)
      Stop
    End If

    !** Read the no. of points
    Read(unitno, *) npts

    !** Check that the no. of points is the same as that specfied
    !** in the control file
    If (npts /= gcmcspc%npts) Then
      Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
          " The no. of points in the file ", Trim(filename), &
          " does not match that in the control file"
      Stop
    End If
    Read(unitno, *) (gcmcspc%fuglist(i)%pressure, i=1, npts)

    Close(unit=unitno)

  End Subroutine gcmcmoves_readpressures

  !-----------------------------------------------------------
  ! Generates the actual fugacity list from the given intervals
  ! contains the start fugacity, end fugacity and no. of points
  !--------------------------------------------------------
  Subroutine gcmcmoves_genpressures(gcmcspc, startp, endp)
    Type(GCMC_Move_Params), Intent(inout) :: gcmcspc
    Real(kind=RDbl), Intent(in)           :: startp, endp
    Real(kind=RDbl)     :: pressincr 
    Integer             :: i, npts

    ! Check to see what kind of a scale we are going to use
    ! A log scale is used if the ending fugacity is greater
    ! than the starting fugacity by 2 orders of magnitude
    npts = gcmcspc%npts
    If (Abs(Log10(endp/startp)) > 2) Then
      !use a log scale
      If (npts /= 1) Then
        pressincr = Exp(Log(endp/startp)/(npts-1.0_RDbl))
      Endif
      gcmcspc%fuglist(1)%pressure = startp
      Do i=1, npts-1
        gcmcspc%fuglist(i+1)%pressure = startp*(pressincr**i)
      End Do
    Else
      If (npts /= 1) Then
        pressincr = (endp - startp)/(npts - 1.0_RDbl)
      End If
      gcmcspc%fuglist(1)%pressure = startp
      Do i=1, npts-1
        gcmcspc%fuglist(i+1)%pressure = startp + pressincr*i
      End Do      
    Endif
    Return
  End Subroutine gcmcmoves_genpressures

  !----------------------------------------------------------------------
  ! Cleans up the various pointers once we are done with this structure.
  !----------------------------------------------------------------------
  Subroutine gcmcmoves_cleanup(gcmcspc)
    Type(GCMC_Move_Params), Dimension(:), Intent(inout) :: gcmcspc

    Integer    :: i, nspc, error
    
    nspc = Size(gcmcspc, 1)
    Do i = 1,nspc
      Deallocate(gcmcspc(i)%fuglist, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'fuglist')   
    End Do

  End Subroutine gcmcmoves_cleanup

  !-------------------------------------------------------------
  ! Display the intermediate statistics from the run
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !-------------------------------------------------------------
  Subroutine gcmcmoves_displaystats(gcmcspc, simno, indent, optunit)
    Type(GCMC_Move_Params), Dimension(:), Intent(In) :: gcmcspc
    Integer, Intent(In)                              :: simno,indent
    Integer, Optional, Intent(In)                    :: optunit

    Integer                 :: i, j, unitno, nspc, spc 
    Character(len=indent)   :: blank
    Character(len=strLen)   :: molecname,string
    Real(kind=RDbl)         :: insertratio,deleteratio,transratio,rotratio

    blank = Repeat(' ',indent)
    unitno = 6
    If (Present(optunit)) unitno = optunit

    nspc = Size(gcmcspc, 1)
    Do i = 1,nspc
      spc = gcmcspc(i)%spc
      molecname = molecules_name(spc)
      string = real2str(gcmcspc(i)%fuglist(simno)%pressure,6)
      Write(unitno,'(3a,4x,3a)') blank,"Species Name: ",Trim(molecname), &
          "Pressure: ",Trim(string)," kPa"

      Do j = 1,gcmcspc(i)%no_of_movetypes
        Call mcmoves_stats(gcmcspc(i)%mcmoves(j),indent+2,unitno)
      End Do

      Call stats_display(gcmcspc(i)%nmoles,indent+2)

    End Do

  End Subroutine gcmcmoves_displaystats

  !-------------------------------------------------------------------------
  ! resets the stats, and counter to zero
  ! Requires: gcmcparams -- GCMC simulation parameters
  !           simno -- simulation number
  !-------------------------------------------------------------------------
  Subroutine gcmcmoves_beginsim(gcmcspc,simno)
    Type(GCMC_Move_Params), Intent(InOut)    :: gcmcspc
    Integer, Intent(In)                      :: simno

    Integer :: i

    Call stats_reset(gcmcspc%nmoles)

    Do i = 1,gcmcspc%no_of_movetypes
      Call mcmoves_resetStats(gcmcspc%mcmoves(i))
    End Do

  End Subroutine gcmcmoves_beginsim


  !-------------------------------------------------------------------------
  ! copies the fugacities of other compounds to this gcmcparams
  ! fuglist dimensions are (1:nsims, 1:nsorbs)
  ! Requires: gcmcparams -- GCMC simulation parameters
  !           simno -- simulation number
  !-------------------------------------------------------------------------
  Subroutine gcmcmoves_copyfuglist(gcmcspc,fuglist,simno)
    Type(GCMC_Move_Params), Intent(InOut)    :: gcmcspc
    Type(Fugacity_Params), Dimension(:,:), Pointer   :: fuglist
    Integer, Intent(In)                      :: simno

    Integer :: nspc,spc
    nspc=molecules_getnsorbs()
    Do spc = 1,nspc
      gcmcspc%single_sim_fug_arr(spc)=fuglist(simno,spc)
      gcmcspc%single_sim_B_arr(spc)=fuglist(simno,spc)%B
    End Do

  End Subroutine gcmcmoves_copyfuglist


  !-------------------------------------------------------------
  ! Display the simulation parameters
  ! Requires: gcmcspc -- GCMC species parameters
  !           simno -- simulation number
  !           indent -- no. of spaces from the left margin
  !           unitno -- optional display unit number
  !-------------------------------------------------------------
  Subroutine gcmcmoves_displaysimparams(gcmcspc,simno,indent,optunit)
    Type(GCMC_Move_Params), Dimension(:), Intent(In) :: gcmcspc
    Integer, Intent(In)                              :: simno,indent
    Integer, Optional, Intent(In)                    :: optunit

    Integer                 :: i, j, unit, nspc, spc 
    Character(len=indent)   :: blank
    Character(len=strLen)   :: molecname

    blank = Repeat(' ',indent)
    
    If (Present(optunit)) Then
      unit = optunit
    Else
      unit = 6
    End If
    
    nspc = Size(gcmcspc, 1)
    Write(unit,'(2a)') blank, dashedline
    Write(unit,'(2a)') blank, "GCMC Species Parameters:"

    !** if not associated nspc=0, then does not go into the loop
    Do i = 1,nspc
      spc = gcmcspc(i)%spc
      molecname = molecules_name(spc)
      Write(unit, '(2a)')  blank, dashedline2
      Write(unit, '(2x,3a)') blank,"Species Name: ", &
          Trim(molecname)
      Write(unit, '(a,4x, 4a11)') blank,"Press.(kPa)", "Fugacity", "murti", &
          "B(Adams)"
      Write(unit, '(a,4x,4f11.3)') blank, &
          gcmcspc(i)%fuglist(simno)%pressure, &
          gcmcspc(i)%fuglist(simno)%fugacity, &
          gcmcspc(i)%fuglist(simno)%murti, &
          gcmcspc(i)%fuglist(simno)%B
    End Do
    Write(unit, '(2a)')  blank, dashedline2

  End Subroutine gcmcmoves_displaysimparams

  !-----------------------------------------------------------------------
  ! Display one line statistics for each species (fugacity, loading)
  ! Requires: gcmcspc -- GCMC species parameters
  !           simno -- simulation number
  !           indent -- no. of spaces from the left margin
  !           unitno -- optional display unit number
  !-----------------------------------------------------------------------
  Subroutine gcmcmoves_dispspcstats(gcmcspc,simno,indent,optunit)
    Type(GCMC_Move_Params), Intent(In)    :: gcmcspc
    Integer, Intent(In)                   :: simno,indent
    Integer, Optional, Intent(In)         :: optunit

    Integer                      :: unit
    Real(kind=RDbl)              :: loading
    Character(len=indent)        :: blank
    Character(len=strLen)        :: molecname,string1,string2

    unit = 6
    If (Present(optunit)) unit = optunit
    blank = Repeat(' ',indent)    

    molecname = molecules_name(gcmcspc%spc)
    loading = stats_getcavg(gcmcspc%nmoles)/(gcmcspc%nunitcells*one)
    string1 = real2str(gcmcspc%fuglist(simno)%pressure,5)
    string2 = real2str(loading,5)
    Write(unit,'(3a,2(2x,2a))') blank,Trim(molecname),':','P (kPa): ', &
        Trim(string1),'  loading (molec/uc): ',Trim(string2)

  End Subroutine gcmcmoves_dispspcstats

  !--------------------------------------------------------
  ! Display the gcmc species information
  ! Requires: gcmcspc -- GCMC species parameters
  !           indent -- no. of spaces from the left margin
  !           unitno -- optional display unit number
  !--------------------------------------------------------
  Subroutine gcmcmoves_display(gcmcspc,indent,optunit)
    Type(GCMC_Move_Params), Dimension(:), Intent(In) :: gcmcspc
    Integer, Intent(In)                              :: indent
    Integer, Optional, Intent(In)                    :: optunit

    Integer                 :: i, j, unit, nspc, spc 
    Character(len=indent)   :: blank
    Character(len=strLen)   :: molecname

    blank = Repeat(' ',indent)
    
    If (Present(optunit)) Then
      unit = optunit
    Else
      unit = 6
    End If

    nspc = Size(gcmcspc, 1)
    Write(unit,'(2a)') blank, dashedline
    Write(unit,'(2a)') blank, "GCMC Moves Parameters:"

    Do i = 1,nspc
      spc = gcmcspc(i)%spc
      molecname = molecules_name(spc)
      Write(unit, '(2a)')  blank, dashedline
      Write(unit, '(3a)') blank,"Species Name: ", &
          Trim(molecname)
      Call mcmoves_auxrejectdisplay(gcmcspc(i)%auxparams,indent+2,unit)

      Write(unit,'(2x,2a)') blank,"Thermophysical simulation points:"
      Write(unit,'(a,4a10)') blank,"Pressure", "Fugacity", "murti", &
          "B(Adams)"
      Do j = 1,gcmcspc(i)%npts
        Write(unit, '(a,2e10.3,2f10.3)') blank,&
            gcmcspc(i)%fuglist(j)%pressure, &
            gcmcspc(i)%fuglist(j)%fugacity, gcmcspc(i)%fuglist(j)%murti, &
            gcmcspc(i)%fuglist(j)%B
      End Do

      !** Write the move-params
!LC      Write(unit,'(2x,2a)') blank,dashedline
      Write(unit,'(2x,2a)') blank,"Allowed Move Types and Parameters:"
!LC      Write(unit,'(2x,2a)') blank,dashedline

      Do j = 1,gcmcspc(i)%no_of_movetypes
        Call mcmoves_display(gcmcspc(i)%mcmoves(j),indent+4,unit)
      End Do
      If (i /= nspc) Write(unit,'(2a)') blank,dashedline
    End Do
   
  End Subroutine gcmcmoves_display
  
End Module gcmcmoves

!!$
!!$
!!$  !---------------------------------------------------------
!!$  ! Parses the fugacity line to set the fugacities of the
!!$  ! species type.  It can either do the fugacities at equal
!!$  ! intervals or read them from a file.
!!$  !---------------------------------------------------------
!!$  Subroutine gcmcmoves_getpressures(gcmcspc, pressureline)
!!$    Type(GCMC_Move_Params), Intent(inout) :: gcmcspc
!!$    Character(*), Intent(in)              :: pressureline
!!$    
!!$    Real(kind=RDbl) :: startp, endp
!!$    Integer         :: i, nfields, npts, convfields, error
!!$    Character(len=strLen)     :: fields(strLen), filename, str
!!$    
!!$    !** Split the fields by commas
!!$    nfields = split(pressureline, fields, ",") 
!!$    
!!$    !** Previously ther was an option here for reading the pressure from a file
!!$    !** It does not simplify anything and makes post-code calculations messy
!!$    !** So NO MORE use of pressurefile
!!$    
!!$    !** If nfields is 1, it's possible that npts is 1.
!!$    If ((nfields == 1).And.(gcmcspc%npts == 1)) Then
!!$      
!!$      !** It's a number. Store it.
!!$      gcmcspc%fuglist(1)%pressure = toreal(fields(1))
!!$      
!!$    Else If (nfields==2) Then
!!$      !** We need to generate the range
!!$      startp = toreal(fields(1))
!!$      
!!$      !** Get the last field after separating it from the rest of the line
!!$      str = fields(2)
!!$      nfields = split(str, fields, " ")
!!$      endp = toreal(fields(1))
!!$       
!!$      !** Fill in the gcmcspc pressure data
!!$      Call gcmcmoves_genpressures(gcmcspc, startp, endp)
!!$      
!!$    Elseif (nfields==gcmcspc%npts) Then
!!$      !** we have all pressure points listed there
!!$      Do i=1,gcmcspc%npts
!!$        gcmcspc%fuglist(i)%pressure=toreal(Trim(fields(i)))
!!$      Enddo
!!$    Else 
!!$      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$      Write(*,*) "could not parse the pressureline"
!!$      
!!$    End If
!!$    
!!$  End Subroutine gcmcmoves_getpressures
!!$
!!$  !---------------------------------------------------------
!!$  ! Read the fugacities from the file "filename"
!!$  !---------------------------------------------------------
!!$  Subroutine gcmcmoves_readpressures(gcmcspc)
!!$    Type(GCMC_Move_Params), Intent(InOut) :: gcmcspc
!!$
!!$    Integer                :: unitno, spc, lineno, npts, i
!!$    Character(len=strLen)  :: line, spcname, filename
!!$
!!$    !** Open the file if not already open
!!$    filename = gcmcspc%pressurefile
!!$    unitno = file_open(filename,110)
!!$
!!$    !** Find the species name in the file
!!$    spc = gcmcspc%spc
!!$    spcname = molecules_name(spc)
!!$    lineno = filesrchstr(unitno, Trim(spcname), line)
!!$    If (lineno == 0) Then
!!$      Write(0,'(1x,2a,i4, 4a)') __FILE__," : ",__LINE__, &
!!$          " Could not find the species ", Trim(spcname), " in the file ",&
!!$          Trim(filename)
!!$      Stop
!!$    End If
!!$
!!$    !** Read the no. of points
!!$    Read(unitno, *) npts
!!$
!!$    !** Check that the no. of points is the same as that specfied
!!$    !** in the control file
!!$    If (npts /= gcmcspc%npts) Then
!!$      Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
!!$          " The no. of points in the file ", Trim(filename), &
!!$          " does not match that in the control file"
!!$      Stop
!!$    End If
!!$    Read(unitno, *) (gcmcspc%fuglist(i)%pressure, i=1, npts)
!!$
!!$    Close(unit=unitno)
!!$
!!$  End Subroutine gcmcmoves_readpressures
!!$
!!$  !-----------------------------------------------------------
!!$  ! Generates the actual fugacity list from the given intervals
!!$  ! contains the start fugacity, end fugacity and no. of points
!!$  !--------------------------------------------------------
!!$  Subroutine gcmcmoves_genpressures(gcmcspc, startp, endp)
!!$    Type(GCMC_Move_Params), Intent(inout) :: gcmcspc
!!$    Real(kind=RDbl), Intent(in)           :: startp, endp
!!$    Real(kind=RDbl)     :: pressincr 
!!$    Integer             :: i, npts
!!$
!!$    ! Check to see what kind of a scale we are going to use
!!$    ! A log scale is used if the ending fugacity is greater
!!$    ! than the starting fugacity by 2 orders of magnitude
!!$    npts = gcmcspc%npts
!!$    If (Abs(Log10(endp/startp)) > 2) Then
!!$      !use a log scale
!!$      If (npts /= 1) Then
!!$        pressincr = Exp(Log(endp/startp)/(npts-1.0_RDbl))
!!$      Endif
!!$      gcmcspc%fuglist(1)%pressure = startp
!!$      Do i=1, npts-1
!!$        gcmcspc%fuglist(i+1)%pressure = startp*(pressincr**i)
!!$      End Do
!!$    Else
!!$      If (npts /= 1) Then
!!$        pressincr = (endp - startp)/(npts - 1.0_RDbl)
!!$      End If
!!$      gcmcspc%fuglist(1)%pressure = startp
!!$      Do i=1, npts-1
!!$        gcmcspc%fuglist(i+1)%pressure = startp + pressincr*i
!!$      End Do      
!!$    Endif
!!$    Return
!!$  End Subroutine gcmcmoves_genpressures
