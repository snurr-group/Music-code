!-----------------------------------------------------------------------------
! This module handles generating random configurations under a given
! energy cutoff using Monte-Carlo type moves.
!-----------------------------------------------------------------------------

Module randconf

  Use defaults, Only: RDbl, strLen, dashedline,dashedline2, twopi, d_res_file, &
      d_con_file, hplanck, MAX_SORBS, Rgas, Nav, calToj, kcalmole_kb, lstrLen, &
      dbgflag
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, readblank, int2str, toreal, deallocerrdisplay, isdigit, &
      real2str, checkandstop 
  Use file, Only: file_getunit,file_gettype,file_open
  Use general, Only: genparams
  Use random, Only: rranf
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use config, Only: AtMolCoords, config_writerestartfile, config_getnmoles, &
      config_display, config_isfixed, config_getnmoleslist
  Use simcell, Only: SimCell_Params, simcell_getvolume
  Use molecule, Only: MolecularParams
  Use molecules, Only: molecules_getmass, molecules_getnsorbs, molecules_name, &
      molecules_gettype
  Use mcmoveset, Only: Move_Set,mcmoveset_readmswts, &
      mcmoveset_pickmove,mcmoveset_setmstag
  Use mcmoves, Only: MC_Move_Params,AuxReject_Params,mcmoves_init, &
      mcmoves_initaux,mcmoves_perturb,mcmoves_getbasictag,mcmoves_display, &
      mcmoves_setMaxNrg,mcmoves_stats,mcmoves_auxrejectdisplay
  Use random, Only: rranf,random_getnewseed
  Use interact, Only: Interaction_Model,interact_checknrgs,interact_initnrgs
  Use subinteract, Only: Subset_Interactions, subinteract_init, &
      subinteract_int, subinteract_oldnrg, subinteract_newnrg

  Implicit None
  Save

  Private
  Public :: randconf_init,RandConf_Params,randconf_dosim,randconf_endsim, &
      randconf_samplecf,randconf_displaystats,randconf_initdisplay, &
      randconf_clean,randconf_display

  !** The tag marking the beginning of the Random configuration control section
  Character(len=strLen), Parameter    :: default_randconf_tag = &
      "Random Configuration Generator"

  !** Contains all the necessary information
  Type RandConf_Params
    Integer                :: spc,niterations
    Real(kind=RDbl)        :: max_nrg

    !** Move information: tags can be "INSERT","DELETE","ROTATE" or "TRANSLATE"
    Integer                                   :: no_of_movetypes
    Type(Move_Set)                            :: moveset
    Type(MC_Move_Params),Dimension(:),Pointer :: mcmoves
    Type(AuxReject_Params), Pointer           :: auxparams

    !** Subset interaction storage, species-dependent, indexed by species number
    Type(Subset_Interactions), Dimension(:), Pointer :: subints
  End Type RandConf_Params

Contains
  !-------------------------------------------------------------------
  ! Initializes the various parameters from the control file
  ! Requires:  params -- random configuration parameters
  !            imodel -- interaction model information
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            ctrl_filename -- the name of the control file
  !            opt_randconftag -- optional Tag to look for in control file
  !-------------------------------------------------------------------
  Subroutine randconf_init(params, imodel, species, simcell, ctrl_filename, &
      opt_randconftag)
    Type(RandConf_Params), Intent(InOut)           :: params
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Character(*), Intent(In)                       :: ctrl_filename
    Character(*), Optional, Intent(In)             :: opt_randconftag

    Integer                :: error,i,j,spc,no_of_movetypes
    Integer                :: unitno,randconflineno,nspc,blocksize
    Logical                :: fast
    Real(kind=RDbl)        :: volume
    Character(len=strLen)  :: spcname, tag, line
    Character(len=lstrLen) :: errormsg

    !** Pad the 'tag' to distinguish it from other incidents of the name
    If (Present(opt_randconftag)) Then
      Write(tag,'(2a)') '- ',Trim(opt_randconftag)
    Else
      Write(tag,'(2a)') '- ',Trim(default_randconf_tag)
    End If
    
    !** Open the ctrl_file if it is not opened
    unitno = file_open(ctrl_filename)
    
    !** Find the section
    randconflineno = filesrchstr(unitno, tag, line, .True.)
    If (randconflineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag '",Trim(tag),"' in the control file"
      Stop
    End If
    
    !** Read the number of moves/step
    Read(unitno,*) line
    line = stripcmnt(line)
    params%niterations = toint(line,'number of random config steps')

    !** Read the maximum energy and convert to kcal/mol
    Read(unitno,*) line
    line = stripcmnt(line)
    params%max_nrg = toreal(line,'maximum energy')
    params%max_nrg = params%max_nrg/calToj

    !** Read the blankline at the end
    errormsg = 'must be a blank line before species parameters'
    Call readblank(unitno,__FILE__,__LINE__,errormsg)

    !** Read species-specific move type information
    Read(unitno, *) spcname
    spcname = Trim(stripcmnt(spcname))
    Write(*,'(2a)') 'Initializing random configuration moves for: ',spcname
    spc = molecules_gettype(spcname)
    If (spc == 0) Then
      Write(0,'(1x,2a,i4,1x,2a)') __FILE__," : ",__LINE__, &
          Trim(spcname), " is not one of the species types in control file"
      Stop
    End If
    
    params%spc = spc
    
    !** Read the auxilliary parameters for rejections
    Nullify(params%auxparams)
    Call mcmoves_initaux(params%auxparams,ctrl_filename,simcell)
    
    !** Read the line containing number of movetypes
    Read(unitno,*) no_of_movetypes
    params%no_of_movetypes = no_of_movetypes
    
    !** Allocate memory 
    Allocate(params%mcmoves(no_of_movetypes), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'mcmoves')

    !** Read the weightings for each move type
    Call mcmoveset_readmswts(params%moveset,ctrl_filename,no_of_movetypes)
    
    !** Read and initialize the different movetypes, set the tags
    Do j = 1,no_of_movetypes
      Call mcmoves_init(params%mcmoves(j),ctrl_filename, &
          simcell,species,"RANDOM",spc,params%auxparams)
      Call mcmoveset_setmstag(params%moveset,j, &
          mcmoves_getbasictag(params%mcmoves(j)))
    End Do

    !** Set maximum energy for each move type
    Do j = 1,no_of_movetypes
      Call mcmoves_setMaxNrg(params%mcmoves(j), params%max_nrg)
    End Do

    !** Initialize the subset interactions array
    nspc = molecules_getnsorbs()   
    Allocate(params%subints(nspc), STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'subints')
    
    !** Initialize the subset interactions for species
    Call subinteract_init(params%subints(spc),imodel,'Molec_System', &
        'MC',(/spc,1,0/),(/0,0,0/))
        
    Close(unit=unitno)
    
    !** Update all the energies stored in molecules
    fast = .True.
    Call interact_initnrgs(imodel,fast,species,simcell,.True.)

#ifdef DEBUG
    !** Just some feedback
    Do i = 1,nspc
      Do j = i,nspc
        Write(*,'(a,2I4,2f18.6)') "Species:",i,j, &
            molecules_getcoul(i,j,'inst'),molecules_getnoncoul(i,j,'inst')
      End Do
    End Do
#endif
    
  End Subroutine randconf_init
 
  !---------------------------------------------------------------------------
  ! Performs "ninterations" of the specified RANDCONF simulation.
  ! Requires:  params -- random configuration parameters
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !---------------------------------------------------------------------------
  Real(kind=RDbl) Function randconf_dosim(params,species,simcell)
    Type(RandConf_Params), Intent(InOut)             :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Type(SimCell_Params), Intent(In)                 :: simcell

    Integer               :: iter,nmoles,spc,indx
    Logical               :: success

    Integer, Dimension(20)   :: molcount

    !** Zero the acceptance counter
    randconf_dosim = 0.0_RDbl

    Do iter = 1,params%niterations

      !** Pick a move type randomly and execute it
      success = randconf_move(params,params%subints,species,simcell)

      If (success) Then
        randconf_dosim = randconf_dosim + 1.0_RDbl
      End If
    End Do

    !** Calculate the acceptance ratio and return it
    randconf_dosim = randconf_dosim/params%niterations

  End Function randconf_dosim

  !--------------------------------------------------------------------
  ! Pick a move and make it using the mcmoves module.  Returns true if
  ! move was successful.
  ! Requires:  params -- random configuration parameters
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !--------------------------------------------------------------------
  Logical Function randconf_move(params,subints,species,simcell)
    Type(RandConf_Params), Intent(InOut)                   :: params
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell

    Integer                   :: nmoles,spc,moveno,molec
    Character(len=strLen)     :: movename
    Character(len=lstrLen)    :: string

    spc = params%spc
    moveno = mcmoveset_pickmove(params%moveset,movename)
    nmoles = config_getnmoles(species,spc)

#ifdef DEBUG
    Write(*,'(2a)') 'RandConf feedback: selected move = ',Trim(movename)
    string = random_getnewseed()
    Write(*,'(2a)') 'RandConf feedback: rranf iseed =   ',Trim(string)
#endif

    Select Case(Trim(movename))
    Case("TRANSLATE", "ROTATE", "TRANSFORM")
      molec = Int(rranf()*nmoles) + 1

      !** Make sure that the no. of molecules is not zero
      If (nmoles /= 0) Then
        randconf_move = mcmoves_perturb(params%mcmoves(moveno),&
            (/spc,molec,0/),0.0_RDbl,subints,species,simcell,.False.)
      Else
        randconf_move = .False.
      End If

    Case Default
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            " Could not interpret movename: ",Trim(movename)
      Stop
    End Select

#ifdef DEBUG  !** very useful debugging feedback
    Write(*,'(2a,t35,2(a,i3),t60,a,l4)') 'RandConf feedback: move = ', &
        Trim(movename),' mol = ',molec,' spc = ',spc, &
        '    Accepted? ',randconf_move
    Write(*,*)
#endif

  End Function randconf_move

  !-----------------------------------------------------------------------
  ! Display the randconf initialization parameters, ie full information 
  ! Requires:  params -- random configuration parameters
  !            indent -- no. of spaces from the left margin
  !            optunit -- optional unit number for display
  !-----------------------------------------------------------------------
  Subroutine randconf_initdisplay(params,indent,optunit)
    Type(RandConf_Params), Intent(In)     :: params
    Integer, Intent(In)                   :: indent
    Integer, Optional, Intent(In)         :: optunit

    Integer               :: i,j,unitno
    Character(len=indent) :: blank
    Character(len=strLen) :: molecname,string

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    Write(unitno,'(2a)') blank, dashedline
    Write(unitno,'(2a)') blank, "The Random Configuration Parameters Section:"
    string = int2str(params%niterations)
    Write(unitno,'(a,2x,2a)') blank,"Number of iterations: ",Trim(string)

    Write(unitno,'(2x,2a)') blank, dashedline
    Write(unitno,'(2x,2a)') blank, "Moves Parameters:"

    molecname = molecules_name(params%spc)
    Write(unitno,'(2x,2a)') blank, dashedline
    Write(unitno,'(2x,3a)') blank,"Species Name: ", &
        Trim(molecname)
    Call mcmoves_auxrejectdisplay(params%auxparams,indent+4,unitno)
    
    !** Write the move-params
    Write(unitno,'(4x,2a)') blank,"Allowed Move Types and Parameters:"
    Do j = 1,params%no_of_movetypes
      Call mcmoves_display(params%mcmoves(j),indent+6,unitno)
    End Do

    Write(unitno,'(2a)') blank,dashedline

  End Subroutine randconf_initdisplay

  !---------------------------------------------------------------------------
  ! Display the Random Configuration statistics for the current simulation
  ! Requires:  params -- random configuration parameters
  !            imodel -- interaction model information
  !            species -- species data structure  
  !            simcell -- the simulation cell information
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional unit number for display
  !---------------------------------------------------------------------------
  Subroutine randconf_displaystats(params,imodel,species,simcell, &
      indent,unitno)
    Type(RandConf_Params), Intent(In)              :: params
    Type(Interaction_Model), Intent(InOut)         :: imodel
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Integer, Intent(In)                            :: indent
    Integer, Intent(In), Optional                  :: unitno

    Integer                 :: i, j, unit, nspc, spc 
    Logical                 :: fast
    Real(kind=RDbl)         :: nrg_devn
    Character(len=indent)   :: blank
    Character(len=strLen)   :: molecname,string

    blank = Repeat(' ',indent)
    
    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    !** Write the simulation no.
    Write(unit, '(2a)') blank, dashedline
    Write(unit, '(2a)') blank, "The NVTMC Stats:"

    spc = params%spc
    molecname = molecules_name(spc)
    
    Write(unit,'(3a)') blank,"Species Name: ",Trim(molecname)
    Do j = 1,params%no_of_movetypes
      Call mcmoves_stats(params%mcmoves(j),indent+2,unit)
    End Do

    If (Trim(genparams%displaymode) == "VERBOSE") Then
      fast = .True.
      nrg_devn = interact_checknrgs(imodel,fast,species,simcell,indent+4,unit)
      string = real2str(nrg_devn,6)
      Write(unitno,'(4x,3a)') blank,"Deviation between stored and newly &
          &calculated energies: ", Trim(string)
    End If

  End Subroutine randconf_displaystats

  !------------------------------------------------------------------
  ! Display the randconf simulation parameters.  Useful for 
  ! summarizing parameters at the beginning of a new sim
  ! Requires:  params -- random configuration parameters
  !            indent -- no. of spaces from the left margin
  !            optunit -- optional unit number for display
  !-----------------------------------------------------------------
  Subroutine randconf_display(params,indent,optunit)
    Type(RandConf_Params), Intent(In) :: params
    Integer, Intent(In)               :: indent
    Integer, Optional, Intent(In)     :: optunit

    Integer                 :: i, unitno, funit
    Character(len=indent)   :: blank
    Character(len=strLen)   :: restartfile, configfile, molecname

    blank = Repeat(' ',indent)    

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    !** Get the configuration and restart file names
    Call file_gettype(d_con_file,configfile,funit)
    Call file_gettype(d_res_file,restartfile,funit)

    !** Write the simulation no.
    Write(unitno, '(2a)') blank,"The Random Configuration Simulation Parameters:"
    Write(unitno, '(a,2x,2a)') blank,"Configuration file: ", &
        Trim(configfile)
    Write(unitno, '(a,2x,2a)') blank, "Restart file      : ", &
        Trim(restartfile)
    molecname = molecules_name(params%spc)
    Write(unitno, '(2x,3a)') blank,"Species Name: ",Trim(molecname)

  End Subroutine randconf_display

  !-------------------------------------------------------------------------
  ! Handles the ending of RANDCONF simulations
  ! 1) Check the accumulated and fresh species-species and intra energies
  ! 2) Write feedback to screen
  ! 3) Optionally write to a restartfile (if stopfile /= '')
  ! 4) halt the program only if desired
  ! Step (4) can be useful for getting a starting file for subsequent MD runs.
  ! Requires: params -- RANDCONF simulation parameters
  !           species -- species data structure  
  !           simcell -- the simulation cell information
  !           stopfile -- name of restart file
  !           indent -- no. of spaces from the left margin
  !           unit -- display unit number
  !           stopflag -- optional flag to signal program halt
  !-------------------------------------------------------------------------
  Subroutine randconf_endsim(params,imodel,species,simcell,stopfile, &
      indent,unit,stopflag)
    Type(RandConf_Params), Intent(In)            :: params
    Type(Interaction_Model), Intent(InOut)       :: imodel
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell
    Character(*), Intent(In)                     :: stopfile
    Integer, Intent(In)                          :: indent,unit
    Logical, Intent(In), Optional                :: stopflag 

    Integer                           :: spc
    Logical                           :: killflag,fast
    Real(kind = RDbl)                 :: rmsdev
    Character(len=indent)             :: blank
    Character(len=strLen)             :: string

    blank = Repeat(' ',indent)    

    If (Present(stopflag)) Then
      killflag = stopflag
    Else
      killflag = .False.
    End If

    !** check the updated and fresh energies for agreement
    fast = .True.
    rmsdev = interact_checknrgs(imodel,fast,species,simcell,0)
    If (rmsdev > 1.0E-2_RDbl) Then
      Write(0,'(1x,2a,i4,a)') __FILE__,": ",__LINE__, &
          ' ERROR: rms deviation between stored and'
      Write(0,'(5x,a,f8.3)') &
          'fresh energies is greater than 0.01 kJ/mol at ',rmsdev
      If (.Not.killflag) Stop
    Else If (rmsdev > 1.0E-4_RDbl) Then
      Write(0,'(1x,2a,i4,a,e14.4)') __FILE__,": ",__LINE__, &
          ' WARNING: rms deviation between stored and'
      Write(0,'(5x,a,e14.4)') 'fresh energies is significant: ',rmsdev
    End If

    !** Feedback 
    Write(unit,'(2a)') blank,dashedline
    Write(unit,'(3a,e9.3)') blank,'Ending random config simulation ', &
        'Nrg dev check: ',rmsdev
    
    If (Trim(stopfile) /= '') Then
      Write(unit,*) "Writing a final restart file into  ", Trim(stopfile)
      Call config_writerestartfile(species, Trim(stopfile))
    End If

    Write(unit,'(a)') dashedline    
    
    !** Halts the program if desired
    If (killflag) Then
      Write(0,*) "Exiting after writing the restartfile"
      Stop
    Endif
    
  End Subroutine randconf_endsim

  !----------------------------------------------------------------------------
  ! Writes a sample section of the control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine randconf_sampleCF(unitno)
    Integer, Intent(In) :: unitno
    
    Write(unitno,'(a)') "---- "//Trim(default_randconf_tag)//" ----"
    Write(unitno,'(a,t30,a)') 'Integer','# Number of iterations per move'
    Write(unitno,'(a,t30,a)') 'Real',&
        '# Maximum energy (kJ/mol)'
    Write(unitno,*) "species-specific move type information"

  End Subroutine randconf_sampleCF

  !----------------------------------------------------------------------
  ! Cleans up the various pointers once we are done with this structure.
  !----------------------------------------------------------------------
  Subroutine randconf_clean(params)
    Type(RandConf_Params), Intent(InOut)        :: params

    Integer    :: i,nspc,error

    Deallocate(params%auxparams, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'auxparams')

    Do i = 1,Size(params%mcmoves)
!LC      Call mcmoves_clean()   !routine doesn't exist at the moment
    End Do

    Deallocate(params%mcmoves, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'mcmoves')

  End Subroutine randconf_clean

End Module randconf
