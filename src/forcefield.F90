!------------------------------------------------------------------------------
! This module serves as a front-end for all forcefield calculations.  Its 
! data structure contains all the parameters for a full system forcefield
! evaluation and for storing the information from the evaluation.  Calls
! from this module are separated into species-species calculations and
! the intramolecular calculations that are done on a single species.
!
! Example qualitative call hierarchies:            
!   forcefield -> ssdriver -> ssbasic -> pairmodel -> lj/buck/coul/etc
!   forcefield -> ssdriver -> ssmap ->-> interpolate 
!   forcefield -> intramolecular -> bsmodel/bbmodel/tormodel/ipmodel/conmodel
! 
! Note that many of the calls require that the fast/slow options be set.
! This allows only those interactions labeled "fast" or "slow" to be evaluated
! The feature will be useful for multiple timestep MD or other splitting
! algorithms.  Must of the calls also require the simcell to be passed
! since many require periodic boundary condition information.
!
! A note concerning units: all returned energies are in kcal/mol
!
! Needed Improvements:
! 1) strongly consider removing forcefield storage from data type
!------------------------------------------------------------------------------

Module forcefield

  Use defaults, Only: RDbl,strLen,scalepe,dashedline,zero,d_ctrl_file, &
      d_aa_file,d_ss_file,lstrLen,dbgflag,xlstrlen
  Use utils, Only: split, filesrchstr, toupper, real2str, findstr, &
      allocErrDisplay,deallocErrDisplay, checkandstop
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use file, Only: file_gettype, file_getunit, file_settag, file_open
  Use atom, Only: atom_checkinit, atoms
  Use molecules, Only: molecules_checkinit, molecules_name, &
      molecules_getnsorbs, molecules_getfilename, molecules_getnatoms
  Use nbrlist, Only: nbrlist_init, Nbrlist_Params  
  Use simcell, Only: SimCell_Params, simcell_getmolectype, simcell_checkinit
  Use pairmodel, Only: pairmodel_clean
  Use config, Only: AtMolCoords, config_isfixed, config_getnmoles, &
      config_getnatoms, config_accelptr, config_conv2accel, config_getaccel, &
      config_dump, config_getnmoleslist, config_getfixed
  Use intramolecular, Only: IntramolecularInfo, intramolecular_init, &
      intramolecular_initstore, intramolecular_molint, intramolecular_spcint, &
      intramolecular_sampleCF, intramolecular_clean, intramolecular_display, &
      intramolecular_hasint, intramolecular_mconint, &
      intramolecular_intrainfo, intramolecular_setparams
  Use sssum, Only: SSSumParams, sssum_initCalc
  Use ssdriver, Only: SpcSpc_Params, ssdriver_init, ssdriver_initstore, &
      ssdriver_ssint, ssdriver_display, ssdriver_clean, ssdriver_link, &
      ssdriver_ison, ssdriver_msint, ssdriver_asint, ssdriver_listparams, &
      ssdriver_getpotparameters, ssdriver_boxinfo
  Use store, Only: Store_Level_Pair,store_display,store_disp
  Use storebase, Only: EnergyPlus,storebase_disp, storebase_display
  Use storetop, Only: Forcefield_Results,storetop_init,storetop_clean, &
      storetop_getforces, storetop_fastforces, storetop_ptrs, storetop_display, &
      storetop_level
  Use storesym, Only: store_levels,storesym_display

  Implicit None

  Private
  Public :: Forcefield_Info, forcefield_init, forcefield_initstore, &
      forcefield_allint, forcefield_clean,forcefield_display, &
      forcefield_putaccels, forcefield_hasint, forcefield_putaccels2, &
      forcefield_msysint, forcefield_asint, forcefield_conint, &
      forcefield_listparams, forcefield_storelevel, &
      forcefield_getpotparameters, forcefield_boxinfo, forcefield_setparams

  Character(len=strLen), Parameter :: default_forcefield_tag = & 
      "Forcefield Information"

  Type Forcefield_Info
    Integer                                           :: nspc
    Character(len=strLen)                             :: id  
    Character(len=3)                                  :: storelevel
    Type(SpcSpc_Params), Dimension(:,:), Pointer      :: ncparams,cparams
    Type(IntramolecularInfo), Dimension(:), Pointer   :: iparams
    Type(Nbrlist_params),Pointer                      :: nlist
  End Type Forcefield_Info

Contains  
  !----------------------------------------------------------------------------
  ! This routine initializes a specified forcefield
  ! Requires:  ff -- forcefield data structure to initialize
  !            id -- forcefield identifier (present in control file)
  !            ctrl_filename -- the control filename to take info from
  !            simcell -- the simulation cell
  !            opt_tag -- optionatl, alternative flag to identify info section
  !----------------------------------------------------------------------------
  Subroutine forcefield_init(ff,id,ctrl_filename,species,simcell,opt_tag)
    Type(Forcefield_Info), Intent(Out)     :: ff
    Character(*), Intent(In)               :: id
    Character(*), Intent(In)               :: ctrl_filename
    Type(AtMolCoords), Dimension(:), Intent(in)  :: species
    Type(SimCell_Params), Intent(In)       :: simcell  
    Character(*), Optional, Intent(In)     :: opt_tag

    Integer                :: ios,error,unitno,lineno
    Integer                :: spc1,spc2,nspc,nfields
    Logical                :: foundit
    Character(len=strLen)  :: tag,srchstr,storelevel
    Character(len=strLen)  :: aa_filename,ss_filename,intra_filename
    Character(len=lstrLen) :: line
    Character(len=strLen), Dimension(strLen) :: fields

    If (Present(opt_tag)) Then
      tag = opt_tag
    Else
      tag = default_forcefield_tag
    End If

    !** Check initialization of necessary modules
    error = 0
    If (.NOT.atom_checkinit(atoms)) error = 1
    If (.NOT.molecules_checkinit()) error = 1
    If (.NOT.simcell_checkinit(simcell)) error = 1
    If (error /= 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Atom, Molecule, Simcell structures ', &
          'must be initialized before forcefield is initialized'
      Stop
    End If

    !** Open the ctrl_file if it is not already open and rewind
    unitno = file_open(ctrl_filename,110)
    Rewind(unit = unitno)

    !** Find the forcefield sections, look for the one matching 'id'
    foundit = .False.
    lineno = 1
    Do While((lineno /= 0).And.(.Not. foundit))
      lineno = filesrchstr(unitno, tag, srchstr)
      Read(unitno,*,IOSTAT=ios) line
      If (ios /= 0) Then
        foundit = .False.
        Exit
      End If
      If (Index(line,id) /= 0) foundit = .True.
    End Do
    If (.Not. foundit) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(1x,3a)') 'WARNING: Could not find the tag "', Trim(tag),'"'
      Write(0,'(2x,4a)') 'in the control file: ', Trim(ctrl_filename), &
          '  with ID: ',Trim(id)
      Write(0,'(2x,3a)') 'Assuming should use first Forcefield section in file'

      !** Find the first occurance of the forcefield section tag
      lineno = filesrchstr(unitno,tag,srchstr,.True.)
    End If

    ff%id = Trim(id)

    !** Read the storage level if it's there and atm-atm filename. 
    Read(unitno,*) line
    nfields = split(line,fields)
    If (findstr(store_levels,fields(1)) == 0) Then
      !      storelevel = 'SPC'     !** Otherwise assume 'SPC'
      storelevel = 'MOL'      !** Otherwise assume 'MOL' (SPC doesn't work yet)
      aa_filename = line
    Else
      storelevel = line
      Read(unitno,*) aa_filename
    End If

    !** make sure the storage level makes sense
    Select Case (ToUpper(Trim(storelevel)))
    Case ('ATM')
      ff%storelevel = 'ATM'
    Case ('MOL')
      ff%storelevel = 'MOL'
    Case ('SPC')
      ff%storelevel = 'SPC'
    Case ('SYS')
      ff%storelevel = 'SPC'
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Sorry, system-level storage option not yet available, use MOL or ATM'
      Stop
    Case Default
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Specified storage level must be: "ATM", "MOL", "SPC" or "SYS"'
      Stop
    End Select

    !** Read the species-species interaction filenames
    Call file_settag(aa_filename,d_aa_file)
    Read(unitno,*) ss_filename
    Call file_settag(ss_filename,d_ss_file)

    !** Read the intramolecular interaction filename OR keyword
    Read(unitno,*) intra_filename

    !** store the number of species
    nspc = molecules_getnsorbs()
    ff%nspc = nspc

    !** Allocate space for intramolecular interactions
    nspc = molecules_getnsorbs()
    Allocate(ff%iparams(nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ff%iparams')

    !** Initialize the intramolecular interaction information
    Select Case (toupper(Trim(intra_filename)))
    Case ("INMOLECULE")
      Do spc1 = 1,nspc
        intra_filename = molecules_getfilename(spc1)
        Call intramolecular_init(ff%iparams(spc1),spc1,intra_filename,simcell)
      End Do

    Case ("HERE")
      Do spc1 = 1,nspc
        intra_filename = ctrl_filename
        Call intramolecular_init(ff%iparams(spc1),spc1,intra_filename,simcell)
      End Do

    Case Default  
      Do spc1 = 1,nspc
        Call intramolecular_init(ff%iparams(spc1),spc1,intra_filename,simcell)
      End Do
    End Select

    Write(*,*) "Finished Initializing Intramolecular interactions"

    !** Allocate the arrays for the coulombic and non-coulombic interactions
    Allocate(ff%cparams(nspc,nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ff%cparams')
    Allocate(ff%ncparams(nspc,nspc),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ff%ncparams')


    !** initialize neighborlist info from ctrlfile
    Nullify(ff%nlist)
    Call nbrlist_init(ff%nlist,species, simcell, ctrl_filename)

    Write(*,*) "Initializing Species-Species interactions...."

    !** Initialize the species-species interaction information
    Do spc1 = 1,nspc
      Do spc2 = spc1,nspc
        Call ssdriver_init(ff%ncparams(spc1,spc2), &
            .False.,spc1,spc2,ss_filename,simcell)
        Call ssdriver_init(ff%cparams(spc1,spc2), &
            .True.,spc1,spc2,ss_filename,simcell)
        If (Associated(ff%nlist))  Then
          If (Associated(ff%ncparams(spc1,spc2)%basic)) Then
            ff%ncparams(spc1,spc2)%basic%nlist=>ff%nlist
          Endif
        Endif


      End Do
    End Do
    Write(*,*) "Finished Initializing Species-Species interactions."
    !** mirror the spc-spc interactions in the lower triangle
    Do spc1 = 1,nspc
      Do spc2 = (spc1+1),nspc
        Call ssdriver_link(ff%ncparams(spc2,spc1),ff%ncparams(spc1,spc2))
        Call ssdriver_link(ff%cparams(spc2,spc1),ff%cparams(spc1,spc2))
      End Do
    End Do



    Write(*,*) "Finished Initializing Forcefield interactions."
  End Subroutine forcefield_init

  !----------------------------------------------------------------------------
  ! This routine initializes the storage for a specified forcefield.  It's
  ! in the forcefield module because the forcefield understands how to 
  ! structure the storage.
  ! Requires:  storage -- forcefield results storage structure
  !            ff -- forcefield data structure 
  !            species -- species data structure
  !            nderivs -- number of derivatives in storage
  !            storelevel -- string indicating desired storage level
  !----------------------------------------------------------------------------
  Subroutine forcefield_initstore(storage,ff,species,nderivs,storelevel)
    Type(Forcefield_Results), Intent(InOut)      :: storage
    Type(Forcefield_Info), Intent(In)            :: ff
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: nderivs
    Character(len=3), Intent(In)                 :: storelevel

    Integer                            :: nmoles1,nmoles2,spc1,spc2,intralevel
    Integer, Dimension(Size(species))  :: nmoleslist
    Character(len=xlstrLen)            :: string

    !** Get the number of molecules of each species
    Call config_getnmoleslist(species,nmoleslist)

    !** Allocate the top of the storage structure
    Call storetop_init(storage,ff%nspc,nmoleslist,storelevel,nderivs)

    !** Set the number of molecules of any species to one if it's zero
    Do spc1 = 1,Size(species)
      If (nmoleslist(spc1) == 0) nmoleslist(spc1) = 1
    End Do

    !** Allocate the species-based storage
    Do spc1 = 1,ff%nspc
      nmoles1 = nmoleslist(spc1)

      !** Intramolecular storage, use MOL-level or ATM-level
      intralevel = Max(3,storetop_level(storelevel))
      storage%intraon(spc1) = intramolecular_hasint(ff%iparams(spc1))
      Call intramolecular_initstore(ff%iparams(spc1),storage%intra(spc1), &
          spc1,nderivs,nmoles1,intralevel)
#ifdef DEBUG
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      string = storebase_disp(storage%intra(spc1)%total)
      Write(*,*) spc1,Trim(string)
      Call store_disp(storage%intra(spc1),.False.,2,6)
#endif

      Do spc2 = spc1,ff%nspc
        nmoles2 = nmoleslist(spc2)

        !** NON-Coulombic storage
        storage%ncoul%on(spc1,spc2) = ssdriver_ison(ff%ncparams(spc1,spc2))
        !        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        !        Write(*,*) 'non-coul ',spc1,spc2,storage%ncoul%on(spc1,spc2)
        Call ssdriver_initstore(ff%ncparams(spc1,spc2),storage%ncoul, &
            spc1,nmoles1,spc2,nmoles2,nderivs,storelevel)

        !** Coulombic storage
        storage%coul%on(spc1,spc2) = ssdriver_ison(ff%cparams(spc1,spc2))
        !        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        !        Write(*,*) 'coul ',spc1,spc2,storage%coul%on(spc1,spc2)
        Call ssdriver_initstore(ff%cparams(spc1,spc2),storage%coul, &
            spc1,nmoles1,spc2,nmoles2,nderivs,storelevel)
      End Do

    End Do

    Return

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call storetop_display(storage,.False.,2,6)

    Stop

  End Subroutine forcefield_initstore

  !----------------------------------------------------------------------------
  ! This routine initializes any parameters, values, or variables that need
  ! to be updated or calculated each time BEFORE calling forcefield drivers.
  ! Requires: ff -- forcefield parameters structure
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           fast -- the logical flag that toggles between fast/slow evals
  !----------------------------------------------------------------------------
  Subroutine forcefield_initCalc(ff,species,simcell,fast)
    Type(Forcefield_Info), Intent(In)            :: ff
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell  
    Logical, Intent(In)                          :: fast

    Logical                          :: sum
    Integer                          :: i, nspcs, spc1, spc2
    Integer, Dimension(Size(species,1)**2,2)         :: sumSpcs
    Type(SSSumParams), Dimension(Size(species,1)**2) :: sumparams

    !** initialize counter i 
    i = 0
    sumSpcs = 0

    !** number of spc types
    nspcs = molecules_getnsorbs()

    !** Ugh, silly loop. Is there a better way to do this?
    Do spc1 = 1, nspcs
      Do spc2 = spc1, nspcs
        If (Associated(ff%cparams(spc1,spc2)%sum)) Then
          sum = .True.
          i = i+1
          sumSpcs(i,1) = spc1
          sumSpcs(i,2) = spc2
          sumparams(i) = ff%cparams(spc1,spc2)%sum
        End If
      End Do
    End Do

    !** Call the initCalc routine in summation
    If (sum) Call sssum_initCalc(sumparams(1:i),sumSpcs(1:i,1:2), &
        species,simcell,fast)

  End Subroutine forcefield_initCalc

  !----------------------------------------------------------------------------
  ! This routine evaluates ALL the inter- and intra-molecular interactions 
  ! in the forcefield.  It is setup so that it only evaluates the fast or 
  ! the slow components.  It does NOT zero the storage structure before it 
  ! starts.  Returns false if it was unable to make a full evaluation.
  ! Requires: ff -- forcefield parameters structure
  !           ffout -- output from forcefield
  !           species -- species data structure
  !           simcell -- the simulation cell information
  !           fast -- the logical flag that toggles between fast/slow evals
  !           pcalc -- optional logical flag for allowing brute force map calcs
  !----------------------------------------------------------------------------
  Logical Function forcefield_allint(ff,ffout,species,simcell,fast,pcalc)
    Type(Forcefield_Info), Intent(InOut)         :: ff
    Type(Forcefield_Results), Intent(InOut)      :: ffout
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell  
    Logical, Intent(In)                          :: fast
    Logical, Intent(In), Optional                :: pcalc

    Integer                        :: spc,spc1,spc2,nspc,simcellspc
    Logical                        :: pcalcflag,allintflag,success

    !** set defaults
    If (.Not.Present(pcalc)) Then 
      pcalcflag = .False.
    Else
      pcalcflag = pcalc
    End If
    forcefield_allint = .False.

    !** Avoid double counting between molecules of same species
    allintflag = .True.

    nspc = molecules_getnsorbs()

    !** Calculate any parameters that may be needed to handle the
    !** species-species interactions
    Call forcefield_initCalc(ff,species,simcell,fast)

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Dumping incoming ncoul structure'
    spc = file_open('ffin.txt')
    Call storesym_display(ffout%ncoul,.False.,2,spc)
    Close(unit=spc)
    !    stop
#endif    

    !** Get the new quantities from the species-species interactions
    !** These cover interactions that can be calculated as pairwise additive
    Do spc1 = 1,nspc
      Do spc2 = spc1,nspc

        !** Evaluate the NON-Coulombic interactions
        success = ssdriver_ssint(ff%ncparams(spc1,spc2),ffout%ncoul, &
            spc1,spc2,species,simcell,fast,pcalcflag,allintflag)
        If (.Not.success) Return

        !** Evaluate the Coulombic interactions
        success = ssdriver_ssint(ff%cparams(spc1,spc2),ffout%coul, &
            spc1,spc2,species,simcell,fast,pcalcflag,allintflag)
        If (.Not.success) Return

      End Do
    End Do

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Dumping ncoul structure'
    spc = file_open('ff.txt')
    Call storesym_display(ffout%ncoul,.False.,2,spc)
    Close(unit=spc)
    !    stop
#endif    

    !** Get the "special" molecule, i.e., the one that fills the simcell
    simcellspc = simcell_getmolectype(simcell)

    !** Do the INTRAmolecular forcefield calculations
    Do spc = 1,nspc
      !** Cycle if this species type is fixed
      If (config_isfixed(species(spc))) Cycle

      !** Do the evaluations
      If (spc == simcellspc) Then
        Call intramolecular_spcint(ff%iparams(spc),ffout%intra(spc), &
            spc,species(spc),fast,simcell)
      Else
        Call intramolecular_spcint(ff%iparams(spc),ffout%intra(spc), &
            spc,species(spc),fast)
      End If
    End Do

    !** if we've made it this far, all is ok
    forcefield_allint = .True.

  End Function forcefield_allint

  !----------------------------------------------------------------------------
  ! This routine evaluates a single MOLECULE-SYSTEM interaction.  It is 
  ! setup so that it only evaluates the fast or the slow components.  It does 
  ! NOT zero the storage structure before it starts.  Returns false if it was 
  ! unable to make a full evaluation.
  ! Requires: ff -- forcefield parameters structure
  !           ffout -- output from forcefield
  !           species -- species data structure
  !           spc -- species number 
  !           mol -- molecule number 
  !           simcell -- the simulation cell information
  !           fast -- the logical flag that toggles between fast/slow evals
  !           skip_intra -- True => skip intramolecular interactions
  !           pcalc -- optional logical flag for allowing brute force map calcs
  !----------------------------------------------------------------------------
  Logical Function forcefield_msysint(ff,ffout,species,spc,mol,simcell, &
      fast,skip_intra,pcalc)
    Type(Forcefield_Info), Intent(InOut)         :: ff
    Type(Forcefield_Results), Intent(InOut)      :: ffout
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: spc,mol
    Type(SimCell_Params), Intent(In)             :: simcell  
    Logical, Intent(In)                          :: fast,skip_intra
    Logical, Intent(In),Optional                 :: pcalc

    Integer                          :: spc2,nspc,natoms,simcellspc,junk
    Logical                          :: pcalcflag,success
    Integer, Dimension(2)            :: indices
    Type(Store_Level_Pair), Dimension(:), Pointer  :: storage_ptr

    !** set defaults
    If (.Not.Present(pcalc)) Then 
      pcalcflag = .False.
    Else
      pcalcflag = pcalc
    End If
    forcefield_msysint = .False.

    nspc = molecules_getnsorbs()

#ifdef DEBUG
    !** dump forcefield storage structure to file for debugging
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    junk = file_open('ffstructure.txt')
    Call storetop_display(ffout,.False.,2,junk)
    Close(unit=junk)
#endif

    !** Calculate any parameters that may be needed to handle the
    !** species-species interactions
    Call forcefield_initCalc(ff,species,simcell,fast)

    !** Get the new interactions from the molecules-species interactions
    !** These cover interactions that can be calculated as pairwise additive
    Do spc2 = 1,nspc

      !** Evaluate the NON-Coulombic interactions
      !      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      !      Write(*,*) 'NON-coul: ',spc,spc2
      success = ssdriver_msint(ff%ncparams(spc,spc2),ffout%ncoul, &
          spc,mol,spc2,species,simcell,fast,pcalcflag)
      If (.Not. success) Return

      !** Evaluate the Coulombic interactions
      !      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      !      Write(*,*) 'coul: ',spc,spc2
      success = ssdriver_msint(ff%cparams(spc,spc2),ffout%coul, &
          spc,mol,spc2,species,simcell,fast,pcalcflag)
      If (.Not. success) Return

    End Do

    !** Get the "special" molecule, i.e., the one that fills the simcell
    simcellspc = simcell_getmolectype(simcell)

    !** Do the INTRAmolecular forcefield calculations
    If ((.Not. skip_intra).And.(intramolecular_hasint(ff%iparams(spc)))) Then 
      If (.Not. config_isfixed(species(spc))) Then

        natoms = molecules_getnatoms(spc)

        !** Get the appropriate pointer
        Call storetop_ptrs(ffout,.True.,.False.,(/spc,mol,0/),(/0,0,0/), &
            indices,storage_ptr)

        !** Do the evaluations
        If (spc == simcellspc) Then
          success = intramolecular_molint(ff%iparams(spc), &
              storage_ptr(indices(1)),species(spc)%coords(1:natoms,mol)%rp, &
              fast,simcell)
        Else
          success = intramolecular_molint(ff%iparams(spc), &
              storage_ptr(indices(1)),species(spc)%coords(1:natoms,mol)%rp,fast)
        End If

        If (.Not. success) Then
          Write(0,'(2a,i5,a,2i3)') __FILE__,":",__LINE__, &
              " Intramolecular evaluation unsuccessful (spc,mol): ",spc,mol 
          Stop
        End If

      End If
    End If

    !** if we've made it this far, all is ok
    forcefield_msysint = .True.

  End Function forcefield_msysint

  !----------------------------------------------------------------------------
  ! This routine evaluates a single ATOM-SPECIES interaction.  It is 
  ! setup so that it only evaluates the fast or the slow components.  It does 
  ! NOT zero the storage structure before it starts.  Returns false if it was 
  ! unable to make a full evaluation.
  ! Requires:  ff -- forcefield parameters structure
  !            ffout -- output from forcefield
  !            species -- species data structure
  !            spc1 -- first species number
  !            mol1 -- first molecule number
  !            atm1 -- first atom number
  !            spc2 -- second species number
  !            simcell -- the simulation cell information
  !            fast -- the logical flag that toggles between fast/slow evals
  !            pcalc -- allow brute force map calculations if True
  !----------------------------------------------------------------------------
  Logical Function forcefield_asint(ff,ffout,species,spc1,mol1,atm1,spc2, &
      simcell,fast,pcalc,hot)
    Type(Forcefield_Info), Intent(InOut)         :: ff
    Type(EnergyPlus), Intent(InOut)              :: ffout
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: spc1,mol1,atm1,spc2
    Type(SimCell_Params), Intent(In)             :: simcell  
    Logical, Intent(In)                          :: fast
    Logical, Intent(In), Optional                :: pcalc
    Real(Kind=RDbl), Dimension(:), Optional :: hot

    Integer           :: nspc,natoms
    Logical           :: pcalcflag,success

    !** set defaults
    If (.Not.Present(pcalc)) Then 
      pcalcflag = .False.
    Else
      pcalcflag = pcalc
    End If
    forcefield_asint = .False.

    !** Calculate any parameters that may be needed to handle the
    !** species-species interactions
    Call forcefield_initCalc(ff,species,simcell,fast)

    !** Evaluate the NON-Coulombic interactions
    If (.Not.Present(hot)) Then
      success = ssdriver_asint(ff%ncparams(spc1,spc2),ffout, &
          spc1,mol1,atm1,spc2,species,simcell,fast,pcalcflag)
    Else
      success = ssdriver_asint(ff%ncparams(spc1,spc2),ffout, &
          spc1,mol1,atm1,spc2,species,simcell,fast,pcalcflag,hot)
    End If
    If (.Not. success) Return

    !** Evaluate the Coulombic interactions
    If (.Not.Present(hot)) Then
      success = ssdriver_asint(ff%cparams(spc1,spc2),ffout, &
          spc1,mol1,atm1,spc2,species,simcell,fast,pcalcflag)
    Else
      success = ssdriver_asint(ff%cparams(spc1,spc2),ffout, &
          spc1,mol1,atm1,spc2,species,simcell,fast,pcalcflag,hot)
    End If
    If (.Not. success) Return

    !** if we've made it this far, all is ok
    forcefield_asint = .True.

  End Function forcefield_asint

  !----------------------------------------------------------------------------
  ! This routine handles the evaluation of constraint interactions.
  ! Requires: ff -- forcefield parameters structure
  !           species -- species data structure
  !           fast -- the logical flag that toggles between fast/slow evals
  !----------------------------------------------------------------------------
  Logical Function forcefield_conint(ff,species,fast)
    Type(Forcefield_Info), Intent(In)              :: ff
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Logical, Intent(In)                            :: fast

    Integer              :: m,spc,nmoles,natoms
    Real(kind=RDbl)      :: uconst,u

    uconst = 0.0_RDbl
    Do spc = 1,molecules_getnsorbs()
      If (config_isfixed(species(spc))) Cycle

      If (.Not. intramolecular_hasint(ff%iparams(spc),'CONSTRAINT')) Cycle

      nmoles = config_getnmoles(species,spc)
      natoms = molecules_getnatoms(spc)

      Do m = 1,nmoles
        If (fast) Then
          Call intramolecular_mconint(ff%iparams(spc), &
              species(spc)%coords(1:natoms,m)%rp, &
              species(spc)%coords(1:natoms,m)%v, &
              u,species(spc)%afast(1:natoms,m),.True.)
          uconst = u + uconst
        Else
          Call intramolecular_mconint(ff%iparams(spc), &
              species(spc)%coords(1:natoms,m)%rp, &
              species(spc)%coords(1:natoms,m)%v, &
              u,species(spc)%aslow(1:natoms,m),.True.)
          uconst = u + uconst
        End If
      End Do

    End Do

    forcefield_conint = .True.

  End Function forcefield_conint

  !----------------------------------------------------------------------------
  ! Take gradient information from the storage structure and put it into the
  ! acceleration arrays in the species structure for use in MD simulations.
  ! The call to storetop_getforces extracts the gradients and the call to 
  ! config_putaccels converts the gradients to accelerations and puts them
  ! into the species arrays.
  ! Requires:  ff -- forcefield data structure
  !            results -- forcefield evaluation results
  !            species -- species data structure
  !            fast -- True => treat Fast interactions, False => treat slow
  !----------------------------------------------------------------------------
  Subroutine forcefield_putaccels(ff,results,species,fast)
    Type(Forcefield_Info), Intent(In)               :: ff
    Type(Forcefield_Results), Intent(In)            :: results
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Logical, Intent(In)                             :: fast

    Integer                                 :: spc,natoms,nmoles
    Type(VecType), Dimension(:,:), Pointer  :: accels

    Do spc = 1,ff%nspc
      If (config_isfixed(species(spc))) Cycle
      nmoles = config_getnmoles(species,spc)
      natoms = config_getnatoms(species,spc)

      !** put the gradients into the species structure
      Call config_getaccel(species,spc,fast,accels)
      accels = VecType(0.0_Rdbl)
      Call storetop_getforces(results,spc,nmoles,natoms,accels)

#ifdef DEBUG
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Call config_dump(species,spc,2,6)
      stop
#endif
      !** convert to accelerations
      Call config_conv2accel(species(spc),spc,fast)
    End Do

  End Subroutine forcefield_putaccels

  !----------------------------------------------------------------------------
  ! An attempt to speed up _putaccels by using routines that directly access
  ! the storage structure.  Helps alot.
  ! Requires:  ff -- forcefield data structure
  !            results -- the results are copied from here to species
  !            species -- species data structure
  !            fast -- True => treat Fast interactions, False => treat slow
  !----------------------------------------------------------------------------
  Subroutine forcefield_putaccels2(ff,results,species,fast)
    Type(Forcefield_Info), Intent(In)               :: ff
    Type(Forcefield_Results), Intent(In)            :: results
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Logical, Intent(In)                             :: fast

    Integer                                 :: spc,natoms,nmoles,a,m
    Type(VecType), Dimension(:,:), Pointer  :: accels

    Do spc = 1,ff%nspc

      If (config_isfixed(species(spc))) Cycle
      nmoles = config_getnmoles(species,spc)
      If (nmoles==0) Cycle
      natoms = config_getnatoms(species,spc)

      !** Get the species acceleration storage and zero it
      Nullify(accels)
      Call config_getaccel(species,spc,fast,accels)
      accels = VecType(0.0_Rdbl)

      !** Put the forces into the species structure
      Call storetop_fastforces(results,spc,nmoles,natoms,accels)

#ifdef DEBUG
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Do m = 1,nmoles
        Write(*,*) 'back from storetop_fastforces (spc,mol): ',spc,m
        Do a = 1,natoms
          Write(*,*) accels(a,m)
        End Do
      End Do
#endif
#ifdef DEBUG
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Call config_dump(species,spc,2,6)
      !      Stop
#endif

      !** Convert to accelerations
      Call config_conv2accel(species(spc),spc,fast)
      Nullify(accels)
    End Do

  End Subroutine forcefield_putaccels2


  !----------------------------------------------------------------------------
  ! Checks to see if a particular intramolecular interaction is turned on
  ! Requires:  ff -- forcefield data structure
  !            spc -- species number
  !            whichOne -- string identifier for model type
  !----------------------------------------------------------------------------
  Logical Function forcefield_hasint(ff,spc,whichOne)
    Type(Forcefield_Info), Intent(In) :: ff
    Integer, Intent(In)               :: spc
    Character(*), Intent(In)          :: whichOne

    forcefield_hasint = intramolecular_hasint(ff%iparams(spc),whichOne)

  End Function forcefield_hasint

  !----------------------------------------------------------------------------
  ! Checks to see if a particular intramolecular interaction is turned on
  ! Requires:  ff -- forcefield data structure
  !            spc -- species number
  !            whichOne -- string identifier for model type
  !----------------------------------------------------------------------------
  Function forcefield_storelevel(ff,depth)
    Character(len=3)                  :: forcefield_storelevel
    Type(Forcefield_Info), Intent(In) :: ff
    Integer, Intent(Out), Optional    :: depth

    forcefield_storelevel = ff%storelevel

    If (Present(depth)) Then
      Select Case (ff%storelevel)
      Case ('ATM')
        depth = 3
      Case ('MOL')
        depth = 2
      Case ('SPC')
        depth = 1
      Case Default
        Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
            ' Unexpected storelevel string: ',ff%storelevel
        Stop      
      End Select
    End If

  End Function forcefield_storelevel

  !----------------------------------------------------------------------------
  ! Returns information about the forcefield parameters governing the 
  ! interactions between two subsets of the system.  This information can
  ! then be used elsewhere.
  ! Requires:  ff -- forcefield data structure
  !            subset1 -- 1st subset (specify only this for intra-only)
  !            subset2 -- 2nd subset
  !            intratype -- intramolecular interaction (uses usual keywords)
  !            nsets -- number of potential sets returned
  !            list -- list of atom numbers for interaction (set,1:Natoms)
  !            params -- array of strings containing type and parameters
  ! Currently works only with intramolecular potentials, needs to be 
  ! generalized to return any forcefield parameters!
  !----------------------------------------------------------------------------
  Subroutine forcefield_listparams(ff,subset1,subset2,intratype,nsets, &
      list,params)
    Type(Forcefield_Info), Intent(In)                 :: ff
    Integer, Dimension(:), Intent(In)                 :: subset1,subset2
    Character(*), Intent(In)                          :: intratype
    Integer, Intent(Out)                              :: nsets
    Integer, Dimension(:,:), Intent(Out)              :: list
    Character(len=lstrLen), Dimension(:), Intent(Out) :: params

    !** ONLY intramolecular!
    Call intramolecular_intrainfo(ff%iparams(subset1(1)), &
        intratype,nsets,list,params)    

    !    Call ssdriver_listparams(ff%ncparams(),subset1,subset2,nsets,list,params)
    !    Call ssdriver_listparams(ff%cparams(),subset1,subset2,nsets,list,params)

  End Subroutine forcefield_listparams

  !----------------------------------------------------------------------------
  ! sets forcefield parameters. this is required to change some
  ! parameter during a simulation. For example if we want to continuously
  ! change the strength of bond
  ! Requires:  imodel -- interaction model information
  !            subset -- speciefies spc numbers
  !            description -- bunchof strings
  !            alist -- lsit of atoms (needed only for intra)
  !            val -- the value
  !----------------------------------------------------------------------------
  Subroutine forcefield_setparams(ff,subset,description, alist,val)
    Type(Forcefield_Info), Intent(InOut)                 :: ff
    Integer, Dimension(:), Intent(In)                 :: subset
    Character(len=strlen),Dimension(:), Intent(In)    :: description
    Integer, Dimension(:), Intent(in)                 :: alist
    Real(kind=RDbl), Intent(in) :: val

    Integer :: ndescr
    ndescr=Size(description)
    If (ndescr<1) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    Select Case(Trim(description(1)))
    Case("INTRA")
      If (ndescr<3) Then
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
      Endif
      Call intramolecular_setparams(ff%iparams(subset(1)),&
          description(2:3), alist,val)
    Case default
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
    End Select

!    Call forcefield_setparams(imodel%ff(1),description, list,val)
  End Subroutine forcefield_setparams


  !----------------------------------------------------------------------------
  ! Returns the box increments for a potential map if they are available.
  ! Requires:  ff -- forcefield data structure
  !            species -- species data structure
  !            boxincrements -- x,y,z increments for potential map
  !            boxsteps -- number of boxes in x,y,z directions
  !----------------------------------------------------------------------------
  Subroutine forcefield_boxinfo(ff,species,boxincrements,boxsteps)
    Type(Forcefield_Info), Intent(In)               :: ff
    Type(AtMolCoords), Dimension(:), Intent(In)     :: species
    Real(Kind=RDbl), Dimension(3), Intent(Out)      :: boxincrements
    Integer, Dimension(3), Intent(Out)              :: boxsteps

    Integer                 :: nfixed,spc2,nspc
    Integer, Dimension(10)  :: list

    nfixed = config_getfixed(species,list)
    nspc = molecules_getnsorbs()

    !** Loop over species pair possibilities until one works
    Do spc2 = 1,nspc
      If (spc2 == list(1)) Cycle
      If (ssdriver_boxinfo(ff%ncparams(list(1),spc2), &
          boxincrements,boxsteps)) Return
    End Do

  End Subroutine forcefield_boxinfo

  !----------------------------------------------------------------------------
  ! This routine displays the forcefield
  ! Requires: ff -- forcefield data structure
  !           indent -- no. of spaces from the left margin
  !           unit -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine forcefield_display(ff,indent,unitno)
    Type(Forcefield_Info), Intent(In)   :: ff
    Integer, Intent(In)                 :: indent
    Integer, Optional, Intent(In)       :: unitno

    Integer               :: unit,spc1,spc2,nspc
    Character(len=indent) :: blank
    Character(len=strLen) :: filename

    blank = Repeat(' ',indent)

    unit = 6
    If (Present(unitno)) unit = unitno

    nspc = Size(ff%iparams)
    Call file_gettype(d_ss_file,filename)

    Write(unitno,'(2a)') blank,dashedline
    Write(unitno,'(4a)') blank,'The "',Trim(ff%id),'" FORCEFIELD information:'

    Write(unit,'(2a)') blank,dashedline

    !** display Non-coulombic spc-spc parameters
    Write(unit,'(2a)') blank,&
        'NON-COULOMBIC SPECIES-SPECIES pair parameters:'
    Write(unit,'(3a)') blank,'Information taken from: ',Trim(filename)
    Do spc1 = 1,nspc
      Do spc2 = spc1,nspc
        Write(unit,'(2x,2a)') blank,dashedline
        Call ssdriver_display(ff%ncparams(spc1,spc2),spc1,spc2,indent+2,unit)
      End Do
    End Do

    Write(unit,'(2a)') blank,dashedline

    !** display coulombic spc-spc parameters
    Write(unit,'(2a)') blank,'COULOMBIC SPECIES-SPECIES pair parameters:'
    Write(unit,'(3a)') blank,'Information taken from: ',Trim(filename)

    Do spc1 = 1,nspc
      Do spc2 = spc1,nspc
        Write(unit,'(2x,2a)') blank,dashedline
        Call ssdriver_display(ff%cparams(spc1,spc2),spc1,spc2,indent+2,unit)
      End Do
    End Do

    Write(unit,'(2a)') blank,dashedline
    Write(unit,'(2a)') blank,'INTRAMOLECULAR information for each species '

    !** display intramolecular parameters
    Do spc1 = 1,nspc
      Write(unit,'(2x,2a)') blank,dashedline
      Write(unit,'(2x,4a)') blank,'Intramolecular information for ', &
          Trim(molecules_name(spc1)),':'
      Call intramolecular_display(ff%iparams(spc1),spc1,indent+4,unit)
    End Do

    Write(unit,'(a)') dashedline

  End Subroutine forcefield_display

  !----------------------------------------------------------------------------
  ! This routine cleans the forcefield
  ! Requires: ff -- forcefield data structure
  !----------------------------------------------------------------------------
  Subroutine forcefield_clean(ff)
    Type(Forcefield_Info), Intent(InOut)   :: ff

    Integer                :: error,spc1,spc2,nspc

    nspc = molecules_getnsorbs()

    !** Deallocate the internals
    Do spc1 = 1,nspc
      Do spc2 = 1,nspc
        Call ssdriver_clean(ff%ncparams(spc1,spc2))
        Call ssdriver_clean(ff%cparams(spc1,spc2))
      End Do
      Call intramolecular_clean(ff%iparams(spc1))
    End Do

    !** Deallocate the pointers themselves
    Deallocate(ff%ncparams, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'ncparams')
    Deallocate(ff%cparams, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'cparams')
    Deallocate(ff%iparams, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'iparams')

  End Subroutine forcefield_clean

  !---------------------------------------------------------------------------
  ! Tina added
  ! Subroutine to get potential parameters, the array pot_params contains
  ! A, B, C, D, hicut, locut and returns zero if the parameters are not
  ! defined for the potential 
  !------------------------------------------------------------------------------
  Subroutine forcefield_getpotparameters(params,spc1,spc2,a1,a2,pot_params)
    Type(forcefield_info), Intent(in) :: params
    Integer, Intent(IN)              :: a1, a2
    Integer, Intent(IN)              :: spc1, spc2
    Real(kind = Rdbl), Dimension(6)  :: pot_params

    Call ssdriver_getpotparameters(params%ncparams(spc1,spc2),a1,a2,pot_params)

  End Subroutine forcefield_getpotparameters

End Module forcefield





