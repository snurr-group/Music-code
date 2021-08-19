!----------------------------------------------------
! This module contains the different data structures
! and routines for doing the GCMC move types.
!----------------------------------------------------
Module hybridgc

  Use defaults, Only: RDbl, strLen, dashedline,dashedline2, twopi, d_res_file &
      ,d_con_file, hplanck, Rgas, Nav
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay
  Use file              ,Only : file_open, file_getunit, file_gettype
  Use random,Only :rranf
  Use eos,Only : eos_Models,eos_init,eos_settemperature,eos_setpressure,&
      eos_getfugacity
  Use config, Only: AtMolCoords, config_writerestartfile
  Use simcell, Only: SimCell_Params, simcell_getvolume
  Use gcmcmoves, Only : GCMC_Move_Params,gcmcmoves_init,gcmcmoves_insert, &
      gcmcmoves_delete, gcmcmoves_nvt, gcmcmoves_displaystats, &
      gcmcmoves_display ,gcmcmoves_displaysimparams, &
      gcmcmoves_getrandommovetype
  Use hybridmoves, Only : Hybrid_Move_Params,hybridmoves_init,&
      hybridmoves_integrate, hybridmoves_multFFinteg, &
      hybridmoves_displaystats,  &
      hybridmoves_displayNrgs,hybridmoves_displayparams
  Use molecules, Only: molecules_getmass 
  Use smap,Only:Smap_Params,smap_init
  Use datafile,Only:CONFILE

  Implicit None
  Save

  Private
  Public :: hybridgc_init,hybridgc_initdisplay,HYBRIDGCMC_Params, &
      hybridgc_dosim,&
      hybridgc_getnosims,hybridgc_displaystats,hybridgc_displaysimparams

  ! The tag marking the beginning of the GCMC section
  Character(len=strLen), Parameter    :: default_hybridgc_tag = &
      "HYBRID GCMC Information"

  Type HYBRIDGCMC_Params
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
    Integer                    :: integ_to_insert_ratio 
    Logical                    :: multFF
  End Type HYBRIDGCMC_Params

  !** we need to look for this tag to see of more than one forcefield exists
  Character(len=2*strLen) :: multFF_check_tag="More Than One INTRA FF" 

Contains
  !----------------------------------------------------------
  ! Initializes the various GCMC parameters from the control
  ! file "ctrl_file"
  !----------------------------------------------------------
  Subroutine hybridgc_init(hybridgcparams, sorbates, simcell, ctrl_filename, &
      opt_hybridtag)
    Type(HYBRIDGCMC_Params)         :: hybridgcparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(in)  :: simcell
    Character(*), Intent(in)  :: ctrl_filename
    Character(*), Optional, Intent(in) :: opt_hybridtag
    Character(len=strLen)     :: eostag, tag, line
    Character(len=strLen), Dimension(strLen) :: fields
    Integer   :: unitno, hybridgclineno, nsorbs, error, i, blocksize,nfields
    Integer   :: lineno
    Real(kind=RDbl)     :: volume

    If (Present(opt_hybridtag)) Then
      tag = opt_hybridtag
    Else
      tag = default_hybridgc_tag
    End If

    !** Open the ctrl_file if it is not opened
    unitno=file_open(ctrl_filename)

    !** Check whether simple hybrid gcmc OR multiple forcefield hybridgcmc
    lineno = filesrchstr(unitno, multFF_check_Tag, line,.True.)    
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4, 4a)') __FILE__," : ",__LINE__, &
          " Could not find the sting ", multFF_check_Tag, " in the file ",&
          Trim(ctrl_filename)
      Write(*,*) "So we will be just using one forcefield for &
          &integration and acceptance"
      hybridgcparams%multFF=.False.
    Else
      Write(*,*) "We will be using one forcefield for integration and &
          & another for acceptance"
      hybridgcparams%multFF=.True.
    Endif
    Rewind(unitno)
    
    !** Find the GCMC section
    hybridgclineno = filesrchstr(unitno, tag, line)
    If (hybridgclineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif

    Read(unitno, *) hybridgcparams%niterations
    Read(unitno, *) hybridgcparams%tk
    Read(unitno, '(a)') eostag
    hybridgcparams%eostag = stripcmnt(eostag)
    Read(unitno, *) hybridgcparams%nsims
    Read(unitno, *) blocksize
    Read(unitno, *) hybridgcparams%smapname
    Read(unitno, *) nsorbs
    hybridgcparams%nsorbs = nsorbs

    Nullify(hybridgcparams%hybrid)

    Read(unitno, '(a)') line
    Write(*,*) "line ",line
    nfields=split(Trim(stripcmnt(line)),fields,",")
    Write(*,*) nfields
    Write(*,*) "fields are fields",(1),fields(2)
    If ( Trim(fields(1) ) == "YES" ) Then 
      Allocate(hybridgcparams%hybrid, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Endif
    hybridgcparams%integ_to_insert_ratio = toint(fields(2))
    Write(*,'(a,i4)') " Ratio of Integration to Insert/Delete : ", &
        hybridgcparams%integ_to_insert_ratio
    Read(unitno, *)
    
    !** Allocate the hybridgcparams%gcmcsorbs 
    Allocate(hybridgcparams%gcmcsorbs(nsorbs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Initialize the sitemap if it is not null
    If (toupper(hybridgcparams%smapname) /= "NULL") Then
      Call smap_init(hybridgcparams%smap, simcell, hybridgcparams%smapname)
    Else
      Nullify(hybridgcparams%smap)
    End If

    !** Set the GCMC move parameter information that should be the same
    !** for all moveable sorbates such as temperature and no. of pts
    Do i=1, nsorbs
      hybridgcparams%gcmcsorbs(i)%tk   = hybridgcparams%tk
      hybridgcparams%gcmcsorbs(i)%npts = hybridgcparams%nsims
      hybridgcparams%gcmcsorbs(i)%blocksize = blocksize
      hybridgcparams%gcmcsorbs(i)%smap      => hybridgcparams%smap
    End Do
    
    !** The rest of the stuff will be read and initialized by the 
    !** gcmcmove init
    Call gcmcmoves_init(hybridgcparams%gcmcsorbs, sorbates, simcell, &
        ctrl_filename)
    
    If (Associated(hybridgcparams%hybrid)) &
        !reads the blnakline before the params
        Read(unitno,*)
      Call hybridmoves_init( hybridgcparams%hybrid, sorbates, &
        hybridgcparams%tk, ctrl_filename )
    
    !** Initialize the parameters for the chosen equation of state 
    Rewind(unitno)
    Call eos_init(hybridgcparams%eosparams, ctrl_filename, &
        hybridgcparams%eostag, hybridgclineno)
   
    !** Generate the fugacitylist 
    volume = simcell_getvolume(simcell)
    Call hybridgc_genfuglist(hybridgcparams, volume)
    Close(unit=unitno)
!    should have an energy int here
  End Subroutine hybridgc_init

  !------------------------------------------------
  ! Fill in the rest of the fields in the fuglist
  !------------------------------------------------
  Subroutine hybridgc_genfuglist(hybridgcparams, volume)
    Type(HYBRIDGCMC_Params), Intent(inout)  :: hybridgcparams
    Real(kind=RDbl), Intent(in)       :: volume

    Integer    :: i, j, sorbtype
    Real(kind=RDbl)   :: tk, pp, fug, B, sivolume, mass, murti, Lambda
    
    !** Set the temperature and the pressure of the different
    !** species
    tk = hybridgcparams%tk
    Call eos_settemperature(hybridgcparams%eosparams, tk)

    !** Get the fugacity
    Do i=1, hybridgcparams%nsims
      ! Set the pressure of each component
      Do j=1, hybridgcparams%nsorbs
        sorbtype = hybridgcparams%gcmcsorbs(j)%sorbtype
        pp = hybridgcparams%gcmcsorbs(j)%fuglist(i)%pressure
        ! pp is in kPa
        Call eos_setpressure(hybridgcparams%eosparams, sorbtype, pp)
      End Do

      ! Get the fugacity of each component
      Do j=1, hybridgcparams%nsorbs
        sorbtype = hybridgcparams%gcmcsorbs(j)%sorbtype
        fug = eos_getfugacity(hybridgcparams%eosparams, sorbtype) ! fug [=] kPa
        hybridgcparams%gcmcsorbs(j)%fuglist(i)%fugacity = fug

        ! Get the excess chemical potential (B in Adams Notation)
        sivolume = volume*1.0e-30*Nav ! convert to m^3/mole
        B = Log(fug*1.0e3*sivolume/(Rgas*tk))
        hybridgcparams%gcmcsorbs(j)%fuglist(i)%B = B 

        ! Get mu the chemical potential
        ! Calculate the DeBroglie wavelength
        mass = molecules_getmass(sorbtype)
        mass = mass*1.0e-3              ! convert to kg
        Lambda = hplanck/Sqrt(twopi*mass/Nav*Rgas/Nav*tk)
        murti = (B - Log(sivolume/Nav/Lambda**3))
        hybridgcparams%gcmcsorbs(j)%fuglist(i)%murti = murti
      End Do
    End Do
  End Subroutine hybridgc_genfuglist


  !--------------------------------------------
  ! Gets the total no. of simulations
  !--------------------------------------------
  Integer Function hybridgc_getnosims(hybridgcparams)
    Type(HYBRIDGCMC_Params), Intent(in)  :: hybridgcparams
    hybridgc_getnosims = hybridgcparams%nsims
  End Function hybridgc_getnosims

  !-------------------------------------------
  ! Does the actual GCMC simulation
  !-------------------------------------------
  Subroutine hybridgc_dosim(hybridgcparams, sorbates, simcell, conf_file,simno)
    Type(HYBRIDGCMC_Params), Intent(inout)  :: hybridgcparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(inout)  :: simcell
    Type(CONFILE) ,Intent(in) :: conf_file
    Integer, Intent(in) :: simno
    Integer             :: iter, nmoles, moveno, sorbtype
    Integer             :: i, j, nsorbs,sorbno,iseed,integ_itn
    Real(kind=RDbl)     :: fugacity
    Character(len=strLen)  :: moveType

    !** Do the iterations at the set conditions
    Do iter = 1, hybridgcparams%niterations
      !** We need to ensure that insertions and deletions
      !** are done with equal probability for each sorbate
      !** The algorithm is as follows ( assumes nosrbs=1 )
      !** 1) Pick a sorbate with equal probability.
      Do integ_itn=1,hybridgcparams%integ_to_insert_ratio
        ! Make sure that gcmcsorbs is assigned, nsorbs should be > 0
        If (hybridgcparams%nsorbs>0) Then
          sorbtype = Int(rranf()*hybridgcparams%nsorbs) + 1
          
          !** call up gcmcmoves and ask which movetype should be selected
          Call gcmcmoves_getRandomMoveType(hybridgcparams%gcmcsorbs(sorbtype)&
              ,moveno,movetype)
        Endif
        
        movetype=Trim(movetype)
        Select Case(movetype)
        Case("INSERT")
          Call gcmcmoves_insert(hybridgcparams%gcmcsorbs(sorbtype), sorbates, &
              simcell, moveno,simno)
        Case("DELETE")
          ! Insert
          Call gcmcmoves_delete(hybridgcparams%gcmcsorbs(sorbtype), sorbates, &
              simcell, moveno,simno)
        Case Default
          Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
              " Could not get movetype"
          Stop
        End Select
      End Do
      
      !SDEBUG      
      !SDEBUG
      !** the insert/delet moves were done above
      !** now it is time yo do configurational move
      If (hybridgcparams%multFF) Then
        Call hybridmoves_multFFinteg(hybridgcparams%hybrid, sorbates, &
            simcell,conf_file)     
      Else
        Call hybridmoves_integrate(hybridgcparams%hybrid, sorbates, &
            simcell,conf_file)     
      Endif
      !SDEBUG
      !SDEBUG
      
    End Do
    
  End Subroutine hybridgc_dosim
  
  !-----------------------------------------------------------
  ! Display the GCMC statistics for the current simulation
  !-----------------------------------------------------------
  Subroutine hybridgc_displaystats(hybridgcparams, sorbates,nspc, optunit)
    Type(HYBRIDGCMC_Params), Intent(in) :: hybridgcparams
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
    Call gcmcmoves_displaystats(hybridgcparams%gcmcsorbs, nspc+4, unitno)
    If (Associated(hybridgcparams%hybrid)) Then
      Call hybridmoves_displaystats(hybridgcparams%hybrid, sorbates,nspc+4, &
          unitno)
    Endif
    
  End Subroutine hybridgc_displaystats

  !------------------------------------------------------------------
  ! Display the gcmc simulation parameters for simulation no. "simno"
  ! The display is sent to the optional unit no. "optunit".  "nspc"
  ! is the no. of spaces to leave from the left margin
  !-----------------------------------------------------------------
  Subroutine hybridgc_displaysimparams(hybridgcparams, simno, nspc, optunit)
    Type(HYBRIDGCMC_Params), Intent(in) :: hybridgcparams

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

    Call gcmcmoves_displaysimparams(hybridgcparams%gcmcsorbs,simno,nspc+2,unitno)
    End Subroutine hybridgc_displaysimparams

  !---------------------------------------------
  ! Display the gcmc initialization parameters
  !---------------------------------------------
  Subroutine hybridgc_initdisplay(hybridgcparams, nspc, optunit)
    Type(HYBRIDGCMC_Params), Intent(in) :: hybridgcparams
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
        hybridgcparams%niterations
    Write(unitno, '(a,2x,a,i6)') spc, "No. of simulations : ", hybridgcparams%nsims
    Write(unitno, '(a,2x,a,f8.3)') spc, "Temperature(K)     : ", hybridgcparams%tk
    Write(unitno, '(a,2x,a,i6)') spc, "No. of sorbates    : ",hybridgcparams%nsorbs
    Call gcmcmoves_display(hybridgcparams%gcmcsorbs, nspc+2, unitno)
    If (Associated(hybridgcparams%hybrid)) then
            Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Call hybridmoves_displayparams(hybridgcparams%hybrid, nspc+4, unitno)
   endif
    
 End Subroutine hybridgc_initdisplay

End Module hybridgc







