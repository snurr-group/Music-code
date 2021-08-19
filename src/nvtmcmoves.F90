!----------------------------------------------------------------------------
! This module contains the different data structures and routines for 
! handling the NVTMC move types.  The data type "NVTMC_Move_Params" contains
! the simulation information for a single species type.  Move-specific
! information and handling are done using the "moves" module through the
! "mcmoves" module.
!
! Important routines are:
!    nvtmcmoves_init -- handles initilization for each species
!    nvtmcmoves_move -- picks a movetype and executes
!    nvtmcmoves_displaystats -- displays intermediate statistics for species
!----------------------------------------------------------------------------

Module nvtmcmoves
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl,strLen,MAX_EXCLSITES,zero,dashedline,dashedline2, &
      kcalmole_kb,one,lstrlen
  Use utils, Only: filesrchstr,toupper,split,stripcmnt,toreal,toint, &
      allocErrDisplay,readblank,real2str,deallocerrdisplay
  Use file, Only: file_open
  Use config, Only: AtMolcoords,config_getnmoles,config_config2xyz
  Use simcell, Only: SimCell_Params,simcell_getnuc
  Use molecules, Only: molecules_name,molecules_gettype,molecules_getnatoms
  Use mcmoveset, Only: Move_Set,mcmoveset_initms,mcmoveset_readmswts, &
      mcmoveset_pickmove,mcmoveset_displayms,mcmoveset_setmstag, &
      mcmoveset_getmswt
  Use mcmoves, Only: MC_Move_Params,AuxReject_Params,mcmoves_init, &
      mcmoves_initaux,mcmoves_perturb, mcmoves_insert,mcmoves_delete, &
      mcmoves_getbasictag,mcmoves_display,mcmoves_auxrejectdisplay, &
      mcmoves_stats, mcmoves_idchange
  Use random, Only: rranf,random_getnewseed
  Use stats, Only: Statistics,stats_setvalue,stats_getcavg,stats_update,&
      stats_display,stats_init, stats_reset
  Use subinteract, Only: Subset_Interactions,subinteract_init
      
  Implicit None
  Save

  Private
  Public :: NVTMC_Move_Params,nvtmcmoves_init,nvtmcmoves_clean, &
      nvtmcmoves_move, nvtmcmoves_display,nvtmcmoves_displaystats, &
      nvtmcmoves_displaysimparams, nvtmcmoves_dispspcstats

  !** All move parameters necessary for each species
  Type NVTMC_Move_Params
    Integer                   :: spc  !* normal index of species 
    Integer                   :: blocksize !* Block size for collecting stats

    !** Thermophysical Properties
    Integer                   :: nunitcells
    Real(kind=RDbl)           :: volume,nmoles

    !** Move information: tags can be "INSERT","DELETE","ROTATE" or "TRANSLATE"
    Integer                                   :: no_of_movetypes
    Type(Move_Set)                            :: moveset
    Type(MC_Move_Params),Dimension(:),Pointer :: mcmoves
    Type(AuxReject_Params), Pointer           :: auxparams
  End Type NVTMC_Move_Params
  
Contains
  !----------------------------------------------------------------
  ! Initializes the various NVTMC parameters from the control file
  ! Requires:  nvtmcspc -- parameters for individual species
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            ctrl_filename -- the name of the control file
  !----------------------------------------------------------------
  Subroutine nvtmcmoves_init(nvtmcspc,species,simcell,ctrl_filename, auxmv)
    Type(NVTMC_Move_Params), Dimension(:), Intent(InOut) :: nvtmcspc
    Type(AtMolCoords), Dimension(:), Intent(InOut)       :: species
    Type(SimCell_Params), Intent(In)                     :: simcell
    Character(*), Intent(In)                             :: ctrl_filename
    Type( AuxMoveObjects),Pointer :: auxmv

    Integer                       :: unitno,i,j,error,blocksize,nfields
    Integer                       :: spc,nspc,npts,nmoles,no_of_movetypes
    Real(kind=RDbl)               :: tk,inswt,delwt
    Character(len=strLen)         :: spcname,pressureline
    Character(len=lstrLen)        :: msg
    Character(len=strLen), Dimension(MAX_EXCLSITES) :: fields
    
    !** from the control file read no of species
    unitno = file_open(ctrl_filename)
    nspc = Size(nvtmcspc, 1)
     
!SDEBUG    Write(*,'(a,i4,a,i4)') __FILE__,__LINE__,": number of species ",nspc

    !** Initialize the parameters for each species
    Do i = 1,nspc
      nvtmcspc(i)%nunitcells = simcell_getnuc(simcell)

      Nullify(nvtmcspc(i)%auxparams)

      Read(unitno, *) spcname
      spcname = Trim(stripcmnt(spcname))
      Write(*,'(2a)') 'Initializing NVTMC moves for: ',spcname
      spc = molecules_gettype(spcname)
      If (spc == 0) Then
        Write(0,'(1x,2a,i4,1x,2a)') __FILE__," : ",__LINE__, &
            Trim(spcname), " is not one of the species types in control file"
        Stop
      End If

      nvtmcspc(i)%spc = spc

      !** Read the auxilliary parameters for rejections
      Call mcmoves_initaux(nvtmcspc(i)%auxparams,ctrl_filename,simcell)

      !** Read the line containing number of movetypes
      Read(unitno,*) no_of_movetypes
      nvtmcspc(i)%no_of_movetypes = no_of_movetypes

      !** allocate memory 
      Allocate(nvtmcspc(i)%mcmoves(no_of_movetypes), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'mcmoves')

      Call mcmoveset_readmswts(nvtmcspc(i)%moveset,ctrl_filename,no_of_movetypes)

      !** Read and initialize the different movetypes, set the tags
      Do j = 1,no_of_movetypes
        Call mcmoves_init(nvtmcspc(i)%mcmoves(j),ctrl_filename, &
            simcell,species,"NVT",auxmv,spc,nvtmcspc(i)%auxparams)
        Call mcmoveset_setmstag(nvtmcspc(i)%moveset,j, &
            mcmoves_getbasictag(nvtmcspc(i)%mcmoves(j)))
      End Do

!LC      Call mcmoveset_displayms(nvtmcspc(i)%moveset,0,6)

      !** Check the input move types 
      If (mcmoveset_getmswt(nvtmcspc(i)%moveset,'TRANSLATE') <= zero) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Move set for ',Trim(spcname),' must contain a TRANSLATE move'
!        Stop
      End If
      If (molecules_getnatoms(spc) > 1) Then
        If (mcmoveset_getmswt(nvtmcspc(i)%moveset,'ROTATE') <= zero) Then
          Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
              ' Move set for multiatomic ',Trim(spcname), &
              ' must contain a ROTATE move'
!          Stop
        End If
      End If

#if ALL
      !** Initialize the statistics variables
      blocksize = nvtmcspc(i)%blocksize
      nmoles = config_getnmoles(species, spc)
      Call stats_init(nvtmcspc(i)%nmoles, &
          "Number of Molecules (ins, blk, cum, stdd)", blocksize,.False., "f7.2")
      Call stats_setvalue(nvtmcspc(i)%nmoles, nmoles*1.0_RDbl)
#endif

      !** Read the blank separating line
      If (i /= nspc) Then 
        Write(msg,'(2a)') 'just finished reading NVT MC moves for species: ', &
            Trim(spcname)
        Call readblank(unitno,__FILE__,__LINE__,msg)
      End If
    End Do

  End Subroutine nvtmcmoves_init

  !--------------------------------------------------------------------
  ! Pick a move and make it using the mcmoves module.  Returns true if
  ! move was successful.
  ! Requires:  nvtmcspc -- parameters for individual species
  !            rti -- 1.0/(Rgas*T_Kelvin)
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            simno -- the simulation number
  !--------------------------------------------------------------------
  Logical Function nvtmcmoves_move(nvtmcspc,rti,subints,species,simcell,simno)
    Type(NVTMC_Move_Params), Intent(InOut)                 :: nvtmcspc
    Real(kind=RDbl), Intent(In)                            :: rti 
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Integer, Intent(In)                                    :: simno

    Integer                   :: nmoles,spc,moveno,molec
    Character(len=strLen)     :: movename
    Character(len=lstrLen)    :: string

    spc = nvtmcspc%spc
    moveno = mcmoveset_pickmove(nvtmcspc%moveset,movename)
    nmoles = config_getnmoles(species,nvtmcspc%spc)

#ifdef DEBUG
    Write(*,'(2a)') 'NVTMC feedback: selected move = ',Trim(movename)
    string = random_getnewseed()
    Write(*,'(2a)') 'NVTMC feedback: rranf iseed =   ',Trim(string)
#endif

    Select Case(Trim(movename))
    Case("TRANSLATE", "ROTATE", "TRANSFORM")
      molec = Int(rranf()*nmoles) + 1

      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
        nvtmcmoves_move = mcmoves_perturb(nvtmcspc%mcmoves(moveno),&
            (/spc,molec,0/),rti,subints,species,simcell,.False.)
      Else
        nvtmcmoves_move = .False.
      End If

    Case("IDCHANGE")
      molec = Int(rranf()*nmoles) + 1

      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
        nvtmcmoves_move = mcmoves_idchange(nvtmcspc%mcmoves(moveno),&
            (/spc,molec,0/),rti,subints,species,simcell)
      Else
        nvtmcmoves_move = .False.
      End If

    Case Default
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            " Could not interpret movename: ",Trim(movename)
      Stop
    End Select

#ifdef DEBUG  !** very useful debugging feedback
    Write(*,'(2a,t34,2(a,i3),t60,a,l4)') 'NVTMC feedback: move = ', &
        Trim(movename),' mol = ',molec,' spc = ',spc, &
        '    Accepted? ',nvtmcmoves_move
    Write(*,*)
#endif

  End Function nvtmcmoves_move

  !----------------------------------------------------------------------
  ! Cleans up the various pointers once we are done with this structure.
  !----------------------------------------------------------------------
  Subroutine nvtmcmoves_clean(nvtmcspc)
    Type(NVTMC_Move_Params), Intent(InOut) :: nvtmcspc

    Integer    :: i,error

    Deallocate(nvtmcspc%auxparams, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'auxparams')

    Do i = 1,Size(nvtmcspc%mcmoves)
!LC      Call mcmoves_clean()   !routine doesn't exist at the moment
    End Do

    Deallocate(nvtmcspc%mcmoves, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'mcmoves')

  End Subroutine nvtmcmoves_clean

  !-------------------------------------------------------------
  ! Display the intermediate statistics from the run
  !-------------------------------------------------------------
  Subroutine nvtmcmoves_displaystats(nvtmcspc, simno, indent, optunit)
    Type(NVTMC_Move_Params), Dimension(:), Intent(In) :: nvtmcspc
    Integer, Intent(In)                              :: simno,indent
    Integer, Optional, Intent(In)                    :: optunit

    Integer                 :: i, j, unit, nspc, spc 
    Character(len=indent)   :: blank
    Character(len=strLen)   :: molecname
    Real(kind=RDbl)         :: insertratio,deleteratio,transratio,rotratio

    blank = Repeat(' ',indent)
    
    If (Present(optunit)) Then
      unit = optunit
    Else
      unit = 6
    End If

    nspc = Size(nvtmcspc, 1)
    Do i = 1,nspc
      spc = nvtmcspc(i)%spc
      molecname = molecules_name(spc)

      Write(unit,'(3a)') blank,"Species Name: ",Trim(molecname)
      Do j = 1,nvtmcspc(i)%no_of_movetypes
        Call mcmoves_stats(nvtmcspc(i)%mcmoves(j),indent+2,unit)
      End Do
    End Do

  End Subroutine nvtmcmoves_displaystats

  !-------------------------------------------------------------
  ! Display the simulation parameters
  ! Requires: nvtmcspc -- NVTMC species parameters
  !           simno -- simulation number
  !           indent -- no. of spaces from the left margin
  !           unitno -- optional display unit number
  !-------------------------------------------------------------
  Subroutine nvtmcmoves_displaysimparams(nvtmcspc,simno,indent,optunit)
    Type(NVTMC_Move_Params), Dimension(:), Intent(In) :: nvtmcspc
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
    
    nspc = Size(nvtmcspc, 1)
    Write(unit,'(2a)') blank, dashedline
    Write(unit,'(2a)') blank, "NVTMC Species Parameters:"

    !** if not associated nspc=0, then does not go into the loop
    Do i = 1,nspc
      spc = nvtmcspc(i)%spc
      molecname = molecules_name(spc)
      Write(unit, '(2a)')  blank, dashedline2
      Write(unit, '(2x,3a)') blank,"Species Name: ", &
          Trim(molecname)
    End Do
    Write(unit, '(2a)')  blank, dashedline2

  End Subroutine nvtmcmoves_displaysimparams

  !-----------------------------------------------------------------------
  ! Display one line statistics for each species (fugacity, loading)
  ! Requires: nvtmcspc -- NVTMC species parameters
  !           simno -- simulation number
  !           indent -- no. of spaces from the left margin
  !           unitno -- optional display unit number
  !-----------------------------------------------------------------------
  Subroutine nvtmcmoves_dispspcstats(nvtmcspc,simno,indent,optunit)
    Type(NVTMC_Move_Params), Intent(In)    :: nvtmcspc
    Integer, Intent(In)                   :: simno,indent
    Integer, Optional, Intent(In)         :: optunit

    Integer                      :: unit
    Real(kind=RDbl)              :: realnumber
    Character(len=indent)        :: blank
    Character(len=strLen)        :: molecname

    unit = 6
    If (Present(optunit)) unit = optunit
    blank = Repeat(' ',indent)    

    molecname = molecules_name(nvtmcspc%spc)
    Write(unit,'(2a)') blank,Trim(molecname)

  End Subroutine nvtmcmoves_dispspcstats

  !--------------------------------------------------------
  ! Display the nvtmc species information
  !--------------------------------------------------------
  Subroutine nvtmcmoves_display(nvtmcspc,indent,optunit)
    Type(NVTMC_Move_Params), Dimension(:), Intent(In) :: nvtmcspc
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

    nspc = Size(nvtmcspc, 1)
    Write(unit,'(2a)') blank, dashedline
    Write(unit,'(2a)') blank, "NVTMC Moves Parameters:"

    Do i = 1,nspc
      spc = nvtmcspc(i)%spc
      molecname = molecules_name(spc)
      Write(unit, '(2a)')  blank, dashedline
      Write(unit, '(3a)') blank,"Species Name: ", &
          Trim(molecname)
      Call mcmoves_auxrejectdisplay(nvtmcspc(i)%auxparams,indent+2,unit)

      !** Write the move-params
!LC      Write(unit,'(2x,2a)') blank,dashedline
      Write(unit,'(2x,2a)') blank,"Allowed Move Types and Parameters:"
!LC      Write(unit,'(2x,2a)') blank,dashedline

      Do j = 1,nvtmcspc(i)%no_of_movetypes
        Call mcmoves_display(nvtmcspc(i)%mcmoves(j),indent+4,unit)
      End Do
      If (i /= nspc) Write(unit,'(2a)') blank,dashedline
    End Do
   
  End Subroutine nvtmcmoves_display
  
End Module nvtmcmoves


 
