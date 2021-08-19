!*************** NOTE ******************
!This is just a copy of hybrid-gc which was doing only
!NVT part, has to be cleaned up
!----------------------------------------------------
! This module contains the different data structures
! and routines for doing the GCMC move types.
!----------------------------------------------------
Module hybridnvt

  Use defaults, Only: RDbl, strLen, dashedline,dashedline2, twopi, d_res_file &
      ,d_con_file, hplanck, Rgas, Nav
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay
  Use file, Only: file_getunit,file_gettype,file_open
  Use random,Only :rranf
  Use eos,Only : eos_Models,eos_init,eos_settemperature,eos_setpressure,&
      eos_getfugacity
  Use config, Only: AtMolCoords, config_writerestartfile
  Use simcell, Only: SimCell_Params, simcell_getvolume
!  Use gcmcmoves, Only : GCMC_Move_Params,gcmcmoves_init,gcmcmoves_insert, &
      gcmcmoves_delete, gcmcmoves_nvt, gcmcmoves_displaystats, &
      gcmcmoves_display ,gcmcmoves_displaysimparams
  Use hybridmoves, Only : Hybrid_Move_Params,hybridmoves_init,&
      hybridmoves_integrate, hybridmoves_displaystats,  &
      hybridmoves_displayNrgs,hybridmoves_displayparams
  Use molecules, Only: molecules_getmass 
  Use smap,Only:Smap_Params,smap_init
  Use datafile,Only:CONFILE

  Implicit None
  Save

  Private
  Public :: hybrid_init,hybridnvt_initdisplay,HYBRIDNVTMC_Params, &
      hybridnvt_dosim,&
      hybridnvt_getnosims,hybridnvt_displaystats,hybridnvt_displaysimparams


  ! The tag marking the beginning of the GCMC section
  Character(len=strLen), Parameter    :: default_hybridnvt_tag = &
      "HYBRID GCMC Information"

  Type HYBRIDNVTMC_Params
    ! gcmcsorbs are NOT indexed by sorbate type.  They go from 1-nsorbs
    Type(GCMC_Move_Params), Dimension(:), Pointer   :: gcmcsorbs 
    Type(Hybrid_Move_Params),Pointer                 :: hybrid
    Character(len=strLen)  :: eostag
    Type(EOS_Models)       :: eosparams
    Real(kind=RDbl)        :: tk   ! temperature in Kelvin
    Integer    :: nsorbs
    Integer    :: nsims            ! no. of simulations
    Integer    :: niterations
    Type(Smap_Params), Pointer :: smap
    Character(len=strLen)      :: smapname
  End Type HYBRIDNVTMC_Params

Contains
  !----------------------------------------------------------
  ! Initializes the various GCMC parameters from the control
  ! file "ctrl_file"
  !----------------------------------------------------------
  Subroutine hybridnvt_init(hybridnvtparams, sorbates, simcell, ctrl_filename, &
      opt_hybridtag)
    Type(HYBRIDNVTMC_Params)         :: hybridnvtparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(in)  :: simcell
    Character(*), Intent(in)  :: ctrl_filename
    Character(*), Optional, Intent(in) :: opt_hybridtag
    Character(len=strLen)     :: eostag, tag, line
    Integer   :: unitno, hybridnvtlineno, nsorbs, error, i, blocksize
    Real(kind=RDbl)     :: volume

    If (Present(opt_hybridtag)) Then
      tag = opt_hybridtag
    Else
      tag = default_hybridnvt_tag
    End If

    !** Open the ctrl_file if it is not opened
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
!      Open(file=ctrl_filename, unit=unitno)
      Call file_open(ctrl_filename,unitno)
    Endif
    Rewind(unitno)
    
    !** Find the GCMC section
    hybridnvtlineno = filesrchstr(unitno, tag, line)
    If (hybridnvtlineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif

    Read(unitno, *) hybridnvtparams%niterations
    Read(unitno, *) hybridnvtparams%tk
    Read(unitno, '(a)') eostag
    hybridnvtparams%eostag = stripcmnt(eostag)
    Read(unitno, *) hybridnvtparams%nsims
    Read(unitno, *) blocksize
    Read(unitno, *) hybridnvtparams%smapname
    Read(unitno, *) nsorbs
    hybridnvtparams%nsorbs = nsorbs

    Nullify(hybridnvtparams%hybrid)
    Read(unitno, *) line
    If (Trim(stripcmnt(line))=="YES") Then 
      Allocate(hybridnvtparams%hybrid, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Endif
    
    Read(unitno, *)
    
    !** Allocate the hybridnvtparams%gcmcsorbs 
    Allocate(hybridnvtparams%gcmcsorbs(nsorbs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Initialize the sitemap if it is not null
    If (toupper(hybridnvtparams%smapname) /= "NULL") Then
      Call smap_init(hybridnvtparams%smap, simcell, hybridnvtparams%smapname)
    Else
      Nullify(hybridnvtparams%smap)
    End If

    !** Set the GCMC move parameter information that should be the same
    !** for all moveable sorbates such as temperature and no. of pts
    Do i=1, nsorbs
      hybridnvtparams%gcmcsorbs(i)%tk   = hybridnvtparams%tk
      hybridnvtparams%gcmcsorbs(i)%npts = hybridnvtparams%nsims
      hybridnvtparams%gcmcsorbs(i)%blocksize = blocksize
      hybridnvtparams%gcmcsorbs(i)%smap      => hybridnvtparams%smap
    End Do
    
    !** The rest of the stuff will be read and initialized by the 
    !** gcmcmove init
    Call gcmcmoves_init(hybridnvtparams%gcmcsorbs, sorbates, simcell, &
        ctrl_filename)
    
    If (Associated(hybridnvtparams%hybrid)) &
        !reads the blnakline before the params
        Read(unitno,*)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__    
  Call hybridmoves_init( hybridnvtparams%hybrid, sorbates, &
        hybridnvtparams%tk, ctrl_filename )
          Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    !** Initialize the parameters for the chosen equation of state 
    Rewind(unitno)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Call eos_init(hybridnvtparams%eosparams, ctrl_filename, &
        hybridnvtparams%eostag, hybridnvtlineno)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__   
    !** Generate the fugacitylist 
    volume = simcell_getvolume(simcell)
    Call hybridnvt_genfuglist(hybridnvtparams, volume)
    Close(unit=unitno)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
  End Subroutine hybridnvt_init

  !------------------------------------------------
  ! Fill in the rest of the fields in the fuglist
  !------------------------------------------------
  Subroutine hybridnvt_genfuglist(hybridnvtparams, volume)
    Type(HYBRIDNVTMC_Params), Intent(inout)  :: hybridnvtparams
    Real(kind=RDbl), Intent(in)       :: volume

    Integer    :: i, j, sorbtype
    Real(kind=RDbl)   :: tk, pp, fug, B, sivolume, mass, murti, Lambda
    
    !** Set the temperature and the pressure of the different
    !** species
    tk = hybridnvtparams%tk
    Call eos_settemperature(hybridnvtparams%eosparams, tk)

    !** Get the fugacity
    Do i=1, hybridnvtparams%nsims
      ! Set the pressure of each component
      Do j=1, hybridnvtparams%nsorbs
        sorbtype = hybridnvtparams%gcmcsorbs(j)%sorbtype
        pp = hybridnvtparams%gcmcsorbs(j)%fuglist(i)%pressure
        ! pp is in kPa
        Call eos_setpressure(hybridnvtparams%eosparams, sorbtype, pp)
      End Do

      ! Get the fugacity of each component
      Do j=1, hybridnvtparams%nsorbs
        sorbtype = hybridnvtparams%gcmcsorbs(j)%sorbtype
        fug = eos_getfugacity(hybridnvtparams%eosparams, sorbtype) ! fug [=] kPa
        hybridnvtparams%gcmcsorbs(j)%fuglist(i)%fugacity = fug

        ! Get the excess chemical potential (B in Adams Notation)
        sivolume = volume*1.0e-30*Nav ! convert to m^3/mole
        B = Log(fug*1.0e3*sivolume/(Rgas*tk))
        hybridnvtparams%gcmcsorbs(j)%fuglist(i)%B = B 

        ! Get mu the chemical potential
        ! Calculate the DeBroglie wavelength
        mass = molecules_getmass(sorbtype)
        mass = mass*1.0e-3              ! convert to kg
        Lambda = hplanck/Sqrt(twopi*mass/Nav*Rgas/Nav*tk)
        murti = (B - Log(sivolume/Nav/Lambda**3))
        hybridnvtparams%gcmcsorbs(j)%fuglist(i)%murti = murti
      End Do
    End Do
  End Subroutine hybridnvt_genfuglist


  !--------------------------------------------
  ! Gets the total no. of simulations
  !--------------------------------------------
  Integer Function hybridnvt_getnosims(hybridnvtparams)
    Type(HYBRIDNVTMC_Params), Intent(in)  :: hybridnvtparams
    hybridnvt_getnosims = hybridnvtparams%nsims
  End Function hybridnvt_getnosims

  !-------------------------------------------
  ! Does the actual GCMC simulation
  !-------------------------------------------
  Subroutine hybridnvt_dosim(hybridnvtparams, sorbates, simcell, conf_file,simno)
    Type(HYBRIDNVTMC_Params), Intent(inout)  :: hybridnvtparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(inout)  :: simcell
    Type(CONFILE) ,Intent(in) :: conf_file
    Integer, Intent(in) :: simno
    Integer             :: iter, nmoles, moveno, sorbtype
    Integer             :: i, j, nsorbs,sorbno,iseed
    Real(kind=RDbl)     :: fugacity
    Character(len=strLen)  :: moveType

    !** Do the iterations at the set conditions
    Do iter = 1, hybridnvtparams%niterations
      !** We need to ensure that insertions and deletions
      !** are done with equal probability for each sorbate
      !** type. The algorithm used here is similar to the
      !** one used by Myers.
      !** 1) Pick a sorbate with equal probability.
      !** 2) Decide whether we want to attempt an
      !**    insertion or deletion.




!      sorbtype=1
      !** Pick sorbate type
      !** Pick a movetype and corresponding moveno
      ! Make sure that gcmcsorbs is assigned, nsorbs should be > 0
      If (hybridnvtparams%nsorbs>0) Then
        sorbtype = Int(rranf()*hybridnvtparams%nsorbs) + 1

        !** call up gcmcmoves and ask which movetype should be selected
        Call gcmcmoves_getRandomMoveType(hybridnvtparams%gcmcsorbs(sorbtype)&
            ,moveno,movetype)

      Endif

      movetype=Trim(movetype)

!SHAJI----- HARDCODED FOR DEBUGGING
      movetype="INTEGRATE"
      Select Case(movetype)

      Case("INSERT")
        Call gcmcmoves_insert(hybridnvtparams%gcmcsorbs(sorbtype), sorbates, &
            simcell, moveno,simno)
      Case("DELETE")
        ! Insert
        Call gcmcmoves_delete(hybridnvtparams%gcmcsorbs(sorbtype), sorbates, &
            simcell, moveno,simno)
      Case("TRANSLATE")
        Call gcmcmoves_nvt(hybridnvtparams%gcmcsorbs(sorbtype), sorbates, &
            simcell,moveno)          
      Case("ROTATE")
        Call gcmcmoves_nvt(hybridnvtparams%gcmcsorbs(sorbtype), sorbates, &
            simcell,moveno)
      Case("INTEGRATE")

       Call hybridmoves_integrate(hybridnvtparams%hybrid, sorbates, &
            simcell,conf_file)     

      Case Default
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Could not get movetype"
        Stop
      End Select

    End Do
  End Subroutine hybridnvt_dosim
  
  !-----------------------------------------------------------
  ! Display the GCMC statistics for the current simulation
  !-----------------------------------------------------------
  Subroutine hybridnvt_displaystats(hybridnvtparams, sorbates,nspc, optunit)
    Type(HYBRIDNVTMC_Params), Intent(in) :: hybridnvtparams
    Integer, Intent(in)   :: nspc
    Integer, Optional, Intent(in) :: optunit
    Type(AtMolCoords), Dimension(:), Intent(In) :: sorbates

    Character(len=nspc) :: spc
    Character(len=strLen) :: restartfile, configfile
    Integer    :: i, unitno, funit
    
    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    !** Write the simulation no.
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "The HYBRID GCMC Stats:"
    Call gcmcmoves_displaystats(hybridnvtparams%gcmcsorbs, nspc+4, unitno)

    If (Associated(hybridnvtparams%hybrid)) Then
      Write(*,'(1x,2a,i4)') __FILE__," : Calling hybrid display at",__LINE__
      Call hybridmoves_displaystats(hybridnvtparams%hybrid, sorbates,nspc+4, &
          unitno)
    Endif
    
  End Subroutine hybridnvt_displaystats

  !------------------------------------------------------------------
  ! Display the gcmc simulation parameters for simulation no. "simno"
  ! The display is sent to the optional unit no. "optunit".  "nspc"
  ! is the no. of spaces to leave from the left margin
  !-----------------------------------------------------------------
  Subroutine hybridnvt_displaysimparams(hybridnvtparams, simno, nspc, optunit)
    Type(HYBRIDNVTMC_Params), Intent(in) :: hybridnvtparams

    Integer, Intent(in) :: nspc
    Integer, Intent(in) :: simno
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Character(len=strLen) :: restartfile, configfile
    Integer    :: i, unitno, funit
    
    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

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
    Write(unitno, '(2a)') spc, "The HYBRID GCMC Simulation Parameters:"
    Write(unitno, '(a,2x,a,i4)') spc, "Simulation Number   : ", simno
    Write(unitno, '(a,2x,2a) ') spc, "Configuration file  : ", &
        Trim(configfile)
    Write(unitno, '(a,2x,2a) ') spc, "Restart file        : ", &
        Trim(restartfile)

    Call gcmcmoves_displaysimparams(hybridnvtparams%gcmcsorbs,simno,nspc+2,unitno)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__    
    End Subroutine hybridnvt_displaysimparams

  !---------------------------------------------
  ! Display the gcmc initialization parameters
  !---------------------------------------------
  Subroutine hybridnvt_initdisplay(hybridnvtparams, nspc, optunit)
    Type(HYBRIDNVTMC_Params), Intent(in) :: hybridnvtparams
    Integer, Intent(in) :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer    :: i, unitno
    
    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "The GCMC Parameters Section:"
    Write(unitno, '(a,2x,a,i6)') spc,"No. of iterations  : ", &
        hybridnvtparams%niterations
    Write(unitno, '(a,2x,a,i6)') spc, "No. of simulations : ", hybridnvtparams%nsims
    Write(unitno, '(a,2x,a,f8.3)') spc, "Temperature(K)     : ", hybridnvtparams%tk
    Write(unitno, '(a,2x,a,i6)') spc, "No. of sorbates    : ",hybridnvtparams%nsorbs
!    Call gcmcmoves_display(hybridnvtparams%gcmcsorbs, nspc+2, unitno)
    If (Associated(hybridnvtparams%hybrid)) then
            Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Call hybridmoves_displayparams(hybridnvtparams%hybrid, nspc+4, unitno)
   endif
    
 End Subroutine hybridnvt_initdisplay

End Module hybridnvt







