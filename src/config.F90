!-----------------------------------------------------------------------------
! This module contains the data structures for the positions and the 
! velocities of the molecules and the routines that operate on them.  It also 
! handles the initialization of the species data structure.
!
! Important routines:
!   config_defaultinit -- just allocates the species arrays to sa default size
!   config_init -- reads initialization information from control on a species
!                  by species basis and gets the initial configurations
!
! Need Improvements:
! 1) makes some direct usage of objects in atom and molecules = hack
! 2) lacks requires statements in most cases
! 3) should be expanded to include variable charges on individual atoms
! 4) consider moving some of the specific configuration init routines out
!-----------------------------------------------------------------------------

Module config

  Use defaults, Only: RDbl, strLen, d_ctrl_file, INITIAL_SIZE, INCR_SIZE, &
      scalef, d_res_file, dashedline, dashedline2,scaleke,MAX_SORBS, &
      kjmole_kb, xlstrLen, lstrLen, one, zero, degToRad, MAX_DOF
  Use utils, Only: isfileopen, filesrchstr, toupper, split, stripcmnt, &
      toint, allocErrDisplay, DeallocErrDisplay, getdepth, combine, int2str
  Use file, Only: file_getunit, file_gettype,file_open
  Use general, Only: general_getfiledesc, general_getcurrentiter, &
      general_getcurrenttime
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), IntVecType, vector_display, vector_getnorm, &
      mag, vector_angle
  Use matrix, Only: MatrixType, matrix_b2feuler, Operator(*), Assignment(=), &
      matrix_rotAroundVec
  Use atom, Only: Atomic_Params, atom_getmass, atom_getsymbol
  Use molecules, Only: molecules_name, molecules_getnsorbs, molecules_getcom, & 
      molecules_gettype, molecules_getnatoms, molecules_getfilename, &
      molecules_getatype, molecules_getnthatom, molecules_getcomdefn, &
      molecules_getdof, molecules_getmass, molecules_getinvmasses, &
      molecules_isCoreShell, molecules_getcharges, molecules_getcharge, &
      molecules_atomsymbol, molecules_isflexible, molecules_AtomMass
  Use gcmodels, Only: GeneralizedCoords, gcmodels_initcoords, &
      gcmodels_readrestartfile, Assignment(=), gcmodels_dumprestartinfo, &
      gcmodels_setgencoords, gcmodels_toxyz, gcmodels_fromxyz, &
      gcmodels_changeflex, gcmodels_changedof, gcmodels_setfrom_rp, &
      gcmodels_clean, gcmodels_getgencoords
  Use simcell, Only: SimCell_Params, simcell_getnuc, simcell_pbc, &
      simcell_geteff, simcell_maptosimcell,simcell_getvolume
  Use readoldgcmc, Only: OldGCMCSystem,readoldgcmc_init,readoldgcmc_chkspc, &
      readoldgcmc_nmoles,readoldgcmc_molecstate
  Use readstruc, Only: Structure,readstruc_xform,readstruc_clean,readstruc_set
  Use match, Only: match_compare
  Use visxyz, Only: visxyz_dump, XYZ_Entry, visxyz_make

  Implicit None
  Save

  Private
  Public :: AtMolCoords, RectCoords, default_config_tag, config_display, &
      config_getnmoles, config_defaultinit, config_initdefaultsize, &
      config_getrp, config_getr, config_getMolecCOM, config_getatype, &
      config_checkinit, config_getnatoms, config_getconfig, config_getaccel, &
      config_init, config_writerestartfile, config_allocfields, &
      config_checkandincr,config_setnmoles, config_delmol, config_applypbcs, &
      config_copymolec,config_kineticEnergy, config_isfixed, config_makegrid, &
      config_totalvecs,config_sampleCF,config_getsourcefile, &
      config_getsourcetype, config_zeroforces, config_conv2force, & 
      config_genfromparent, config_initmolecule, config_initFixed, &
      config_config2xyz, config_setfromgcoords, config_changeflex, &
      config_getMolecCTR, config_dump, config_dumpmol, &
      config_kineticTemp, config_accelptr, config_conv2accel, &
      config_getnmoleslist, config_displaycoords, config_getTotMoles, &
      config_nsubsetatoms, config_getmolr, config_getmolq, config_subset2xyz, &
      config_setxyz, config_copysubset, config_rotatexyz, config_getq, &
      config_makeconsistent, config_changedof, config_rotateAroundCom, &
      config_getdistance, config_getangle, config_gettorangle,config_changerp,&
      config_reverseVels, config_subsysKE, config_makeimages, &
      config_initcopy, config_getspcdensity, &
      config_config2struc, config_changefixed, config_simpleinit, &
      config_getCOMVel, config_getSpcCOMVel, config_idrebuild, &
      config_setfromrp, config_getfixed, config_clean, config_xyz2config, &
      config_xyz2molec
  Character(len=strLen), Parameter :: default_config_tag = & 
      "Configuration Initialization"

  Type RectCoords
    Type(VecType)    :: r    ! sim. cell coords.
    Type(VecType)    :: rp   ! principal coords.
    Type(IntVecType) :: cr   ! continuation vectors
    Type(VecType)    :: v    ! velocities
  End Type RectCoords

  !** The coordinates themselves, 2D arrays are defined as (atom:molecule)
  Type AtMolCoords
    Logical                                   :: isinit  !* True if initialized
    Logical                                   :: fixed   !* True if FIXED
    Integer                                   :: natoms
    Integer                                   :: nmoles
    Character(len=strLen)                     :: sourcetype
    Character(len=strLen)                     :: filename
    Type(VecType), Dimension(:,:), Pointer    :: afast  !"fast" accl.
    Type(VecType), Dimension(:,:), Pointer    :: aslow  !"slow" accl.
    Type(RectCoords), Dimension(:,:), Pointer :: coords
    Type(GeneralizedCoords), Dimension(:), Pointer :: gcoords
  End Type AtMolCoords

  Interface config_getnmoles
    Module Procedure config_getnmoles_array
    Module Procedure config_getnmoles_single
  End Interface

  Interface config_getnatoms
    Module Procedure config_getnatoms_array
    Module Procedure config_getnatoms_single
  End Interface

  Interface config_setnmoles
    Module Procedure config_setnmoles_array
    Module Procedure config_setnmoles_single
  End Interface

  Interface config_init
    Module Procedure config_init
    Module Procedure config_defaultinit
  End Interface

Contains
  !--------------------------------------------------------------------
  ! Initializes the "sorbates" data structure for "nsorbates" with default
  ! sizes , does not use the start-up files (just allocates enough memory)
  ! No actual positions stored
  !--------------------------------------------------------------------
  Subroutine config_defaultinit(sorbates)
    Type(AtMolCoords), Dimension(:), Pointer :: sorbates

    Integer    :: nsorbates,sorb,error

    !** Allocate the sorbates structure
    nsorbates = molecules_getnsorbs()  ! number of molecule types
    Allocate(sorbates(nsorbates), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"sorbates")
    
    !** Initialize the sourcetype field
    Do sorb = 1, nsorbates
      sorbates(sorb)%sourcetype = "DFLT_SIZE"
      sorbates(sorb)%filename = "NULL"
      sorbates(sorb)%isinit = .False.
    End Do
    
    !** Allocate memory for co-ordinates
    Do sorb = 1, nsorbates
      
      !** Set the initialization flag
      sorbates(sorb)%isinit = .True.
      
      If (config_isfixed(sorbates(sorb))) Then
        !There is no configuration for this molecule. Skip it.
        !This will probably only be used for the zeolite.
        sorbates(sorb)%natoms = 0
        sorbates(sorb)%nmoles = 0
        Nullify(sorbates(sorb)%aslow)
        Nullify(sorbates(sorb)%afast)
        Nullify(sorbates(sorb)%coords)
        Nullify(sorbates(sorb)%gcoords)
      Else
        ! allocate defaultsize
        Call config_initdefaultsize(sorbates, sorb)
      End If
    End Do

  End Subroutine config_defaultinit
  
  !--------------------------------------------------------------------
  ! Initializes the "sorbates" data structure for "nsorbates" given
  ! the "ctrl_file" and the optional control file tag "optconfig_tag"
  ! Requires: sorbates -- species data structure
  !           simcell -- simulation cell information
  !           ctrl_filename -- control filename to find init info
  !           optconfig_tag -- optional tag to look for in control file
  !--------------------------------------------------------------------
  Subroutine config_init(species, simcell, ctrl_filename, optconfig_tag)
    Type(AtMolCoords), Dimension(:)          :: species
    Type(SimCell_Params), Intent(In)         :: simcell
    Character(*), Intent(In)                 :: ctrl_filename
    Character(*), Optional, Intent(In)       :: optconfig_tag

    Integer                :: unitno, lineno, ios, i, j
    Integer                :: nspc,spcsread,spc
    Character(len=strLen)  :: filename, tag, spcname, srchstr, sourcetype

    !** Pad the 'tag' to distinguish it from other incidents of the name
    If (Present(optconfig_tag)) Then
      Write(tag,'(2a)') '- ',Trim(optconfig_tag)
    Else
      Write(tag,'(2a)') '- ',Trim(default_config_tag)
    End If

    !** Open the ctrl_file if it is not already open
    unitno = file_open(ctrl_filename,110)

    !** Find the configuration initialization section
    lineno = filesrchstr(unitno, tag, srchstr, .True.)
    If (lineno == 0) Then
      Write(0,'(1x, 4a)') 'Could not find the tag "', Trim(tag), &
          '" in the control file ', d_ctrl_file
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

    !** Set the initialization flags for each species
    nspc = molecules_getnsorbs()  ! number of molecule types
    Do i = 1,nspc
      species(i)%isinit = .False.
    End Do

    !** Read the species type and the source information
    Do i = 1,nspc
      spcsread = 0
      Read(unitno, *, IOSTAT=ios) spcname

      !** If ios /= 0 then we haven't found at least one of the speciess.  
      !** Loop through to discover which one hasn't been found
      If (ios /= 0) Then
        Do j = 1,nspc
          spcname = molecules_name(j)
          If (.Not. species(j)%isinit) Then
            Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
                " Could not find configuration information for ", &
                Trim(spcname)
            Stop
          End If
        End Do
      End If
      
      spc = molecules_gettype(spcname)
      If (spc == 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__,' : ',__LINE__, &
            ' Undefined species name specified in config section: ', &
            Trim(spcname)
        Stop      
      End If

      !** Get the information from one of many sources      
      Read(unitno, *, IOSTAT=ios) sourcetype, filename
      species(spc)%sourcetype = sourcetype
      species(spc)%filename = filename
      Call config_initspc(species,spc,simcell)
    End Do

    Close(unitno)

  End Subroutine config_init

  !--------------------------------------------------------------------
  ! Initializes a single species once the source identifier string and 
  ! a filename to read from are set.
  ! Requires:  species -- species data structure
  !            spc -- species number
  !            simcell -- simulation cell information
  !--------------------------------------------------------------------
  Subroutine config_initspc(species,spc,simcell)
    Type(AtMolCoords), Dimension(:)          :: species
    Integer, Intent(In)                      :: spc
    Type(SimCell_Params), Intent(In)         :: simcell

    Character(len=strLen)             :: spcname,filename

    !** Set the initialization flag
    species(spc)%isinit = .True.

    !** Set the fixed coordinate flag to False
    species(spc)%fixed = .False.
    
    !** Now go and get the initial configuration
    filename = species(spc)%filename
    spcname = molecules_name(spc)
    Select Case (ToUpper(Trim(species(spc)%sourcetype)))
    Case ("EQUILFILE")
      Call config_read_equilfile(species,simcell,spcname,filename)
    Case ("RESTARTFILE")
      Call config_readrestartfile(species,simcell,spc,filename)
    Case ("CRASHFILE")
      Call config_readcrashfile(species,simcell,spcname,filename)
    Case ("MEPRESTARTFILE")
      Call config_initmeprestart(species,spc,filename)
    Case ("CUSTOM")
      Call config_initcustomrestart(species,spc,filename)
    Case ("GCMC")
      Call config_initdefaultsize(species,spc)
    Case ("OLDGCMC")
      Call config_readoldgcmc(species,simcell,spcname,filename)
    Case ("DFLT_SIZE","DEFAULT")
      Call config_initdefaultsize(species,spc)
    Case ("MOLECULE")
      !** Just read the coordinates from molecule
      Call config_initmolecule(species,spc)
    Case ("FIXED")
      !** The coordinates are fixed and the molecule is never moved
      Call config_initFixed(species,spc)
    Case Default
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(1x,2a)')  &
          "Don't know how to initialize configurations for ", &
          species(spc)%sourcetype
      Write(0,'(1x,2a)') &
         "Options are: EQUILFILE RESTARTFILE CRASHFILE MEPRESTARTFILE ", &
         "CUSTOM GCMC OLDGCMC DFLT_SIZE MOLECULE FIXED" 
      Stop
    End Select

  End Subroutine config_initspc

  !----------------------------------------------------------------------
  ! Initialize the species configurations using string inputs.  Assumes
  ! that information for all of the species will come from a single file
  ! Requires:  species -- species data structure
  !            simcell -- simulation cell information
  !            filetype -- string defining file type to read from
  !            filename -- filename with information for whole species
  !----------------------------------------------------------------------
  Subroutine config_simpleinit(species,simcell,filetype,filename)
    Type(AtMolCoords), Dimension(:)          :: species
    Type(SimCell_Params), Intent(In)         :: simcell
    Character(*), Intent(In)                 :: filetype,filename

    Integer                   :: spc,nspc
    Character(len=strLen)     :: spcname

    !** Set the initialization flags for each species
    nspc = molecules_getnsorbs()  ! number of molecule types
    Do spc = 1,nspc
      species(spc)%isinit = .False.
    End Do

    !** Read the species type and the source information
    Do spc = 1,nspc
      spcname = molecules_name(spc)

      species(spc)%sourcetype = filetype
      species(spc)%filename = filename
      Call config_initspc(species,spc,simcell)
    End Do

  End Subroutine config_simpleinit

  !----------------------------------------------------------------------------
  ! Writes information about the control file contents to unit unitno
  !----------------------------------------------------------------------------
  Subroutine config_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a,t30,a)') default_config_tag,'# then for each molecule:'
    Write(unitno,'(a,t30,a)') 'Character','# Molecule name'
    Write(unitno,'(a,t30,2a)') 'Character Character', &
        '# Source type (Equilfile, Restartfile, Crashfile, GCMC, Molecule,', &
        'Meprestartfile, Fixed, Noconfig), Filename (or NULL)'
    Write(unitno,'(t30,a)') '# Blank line (required)'

  End Subroutine config_sampleCF

  !-------------------------------------------------------------------------
  ! Allocate the various fields for the coordinates of "sorbate"
  ! It takes the no. of atoms in "natoms" and the no. of molecules "nmoles"
  ! Also initialize the fields to zero.  
  ! NOTE: the number of molecules is NOT set in this routine, this must be
  !       done with config_setnmoles
  ! Requires: species -- single species data structure
  !           spc -- species number
  !           nmoles -- number of molecules to size arrays to
  !           init_gc -- optional flag to control init. of gen. coordinates
  !-------------------------------------------------------------------------
  Subroutine config_allocfields(species, spc, nmoles, init_gc)
    Type(AtMolCoords), Intent(InOut) :: species
    Integer, Intent(In)              :: spc
    Integer, Intent(In)              :: nmoles
    Logical, Optional, Intent(In)    :: init_gc

    Integer               :: error, natoms, i, m, a
    Logical               :: gc, use_ref
    Type(VecType), Dimension(:), Allocatable  :: ref_struc

    !** Check to see if we need to initialize the generalized coords
    gc = .True.
    If (Present(init_gc)) gc = init_gc

    natoms = molecules_getnatoms(spc)

    species%natoms = natoms
    
    Allocate(species%coords(natoms, nmoles), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,a,i5,a,2(a,i5))') __FILE__,__LINE__," : ", &
          " Could not allocate 'species%coords' of size ",natoms," x ",nmoles
      Stop
    End If

    !** Allocate the generalized coordinate array and allocate
    !** memory for the appropriate generalized coordinate model
    If (gc) Then 
      Allocate(species%gcoords(nmoles), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'%gcoords')

      use_ref = .False.

! i think i have done it alright, let's see
!      If (molecules_getdof(spc) > 6) Then
      If (molecules_isflexible(spc)) Then
        natoms = molecules_getnatoms(spc)
        Allocate(ref_struc(natoms), STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ref_struc')
        Call molecules_getcomdefn(spc,ref_struc)
        use_ref = .True.
      End If

      Do i = 1,nmoles
        If (use_ref) Then
          Call gcmodels_initcoords(species%gcoords(i),spc,ref_struc)
        Else
          Call gcmodels_initcoords(species%gcoords(i),spc)
        End If
      End Do

      If (use_ref) Then
        Deallocate(ref_struc, STAT=error)
        If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'ref_struc')
      End If
    End If

    !** Allocate the fast acceleration structure
    Allocate(species%afast(natoms, nmoles), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Could not allocate 'species%afast'"
      Stop
    End If

    !** Allocate the slow acceleration structure
    Allocate(species%aslow(natoms, nmoles), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Could not allocate 'speciess%aslow'"
      Stop
    End If

    !** Initialize all the fields 
    Do m = 1,nmoles
      Do a = 1,natoms
        species%coords(a,m)%rp = 0.0_RDbl
        species%coords(a,m)%r  = 0.0_RDbl
        species%coords(a,m)%cr = 0
        species%coords(a,m)%v  = 0.0_RDbl
        species%afast(a,m) = 0.0_RDbl
        species%aslow(a,m) = 0.0_RDbl
      End Do
    End Do

  End Subroutine config_allocfields

  !-------------------------------------------------------------------------
  ! Initialize and optionally copy one configuration into another.  
  ! Note that the image array size must already be set and be equal to
  ! that of the 'old' configuration.
  ! Requires:  image -- species data structure to be initialize
  !            old -- species number
  !            copy -- if True, will also copy the configuration
  !-------------------------------------------------------------------------
  Subroutine config_initcopy(image,old,copy)
    Type(AtMolCoords), Dimension(:), Intent(Out)   :: image
    Type(AtMolCoords), Dimension(:), Intent(In)    :: old
    Logical, Intent(In)                            :: copy

    Integer     :: spc,mol,atm,i,j
    Integer     :: natms,nmols,nspc

    !** Check the overall shape of the structures
    nspc = Size(old,1)
    If (nspc /= Size(image,1)) Then
      Write(0,'(1x,a,i5,a)') __FILE__,__LINE__," : ", &
          ' image and old array sizes do not match'
      Stop
    End If    

    !** Loop through the individual species and allocate the natms,nmols
    Do spc = 1,nspc
      natms = old(spc)%natoms
      nmols = old(spc)%nmoles
      If (natms /= molecules_getnatoms(spc)) Then
        Write(0,'(1x,a,i5,a)') __FILE__,__LINE__," : ", &
            ' Number of atoms in old configuration is inconsistent with species'
        Stop
      End If

      If (Associated(old(spc)%gcoords)) Then
        Call config_allocfields(image(spc),spc,nmols,.True.)
      Else
        Call config_allocfields(image(spc),spc,nmols,.False.)
      End If

      !** set the number of moles (not done in _allocfields)
      image(spc)%nmoles = old(spc)%nmoles
    End Do

    !** Copy the characteristics fields
    Do spc = 1,nspc
      image(spc)%isinit = old(spc)%isinit
      image(spc)%fixed = old(spc)%fixed
      image(spc)%sourcetype = old(spc)%sourcetype
      image(spc)%filename = old(spc)%filename
    End Do

    If (.Not. copy) Return

    !** Optionally, copy the coordinates themselves
    Do spc = 1,nspc
      natms = old(spc)%natoms
      nmols = old(spc)%nmoles
      Do mol = 1,nmols
        Do atm = 1,natms
          image(spc)%afast(atm,mol) = old(spc)%afast(atm,mol)
          image(spc)%aslow(atm,mol) = old(spc)%aslow(atm,mol)
          image(spc)%coords(atm,mol) = old(spc)%coords(atm,mol)
        End Do
        image(spc)%gcoords(mol) = old(spc)%gcoords(mol)
      End Do
    End Do

  End Subroutine config_initcopy

  !---------------------------------------------------------------------
  ! Allocate the memory for the MEPRESTART File.  Although the file
  ! contains a coordinates for each path in a plane the sorbates structure
  ! only contains one set of coordinates and it cycles through the 
  ! different sets.
  !---------------------------------------------------------------------
  Subroutine config_initmeprestart(sorbates, sorbtype, filename)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Integer, Intent(in)      :: sorbtype
    Character(*), Intent(in) :: filename

    !** Allocate the memory for storing the triplets
    Call config_allocfields(sorbates(sorbtype), sorbtype, 1)
    sorbates(sorbtype)%natoms     = molecules_getnatoms(sorbtype)
    sorbates(sorbtype)%nmoles     = 1
    sorbates(sorbtype)%sourcetype = "MEPRESTARTFILE"
    sorbates(sorbtype)%filename   = filename

  End Subroutine config_initmeprestart

  !-----------------------------------------------------------------------
  ! Allocate the memory for the CUSTOM File.  This applies to any format
  ! that is routine specific.  That is the file contains a format that
  ! is custom designed for solving a specific problem.
  !-----------------------------------------------------------------------
  Subroutine config_initcustomrestart(sorbates, sorbtype, filename)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Integer, Intent(in)      :: sorbtype
    Character(*), Intent(in) :: filename

    !** Allocate the memory for storing the triplets
    sorbates(sorbtype)%natoms = 0
    sorbates(sorbtype)%nmoles = 0
    sorbates(sorbtype)%sourcetype = "CUSTOM"
    sorbates(sorbtype)%filename   = filename
  End Subroutine config_initcustomrestart
    
  !----------------------------------------------------------
  ! This just allocates an initial default amount of space
  ! for the different arrays in config. 
  !----------------------------------------------------------
  Subroutine config_initdefaultsize(sorbates, sorbtype)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Integer, Intent(in)    :: sorbtype

    Integer    :: natoms

    !** Get the memory for the different fields with initial space for
    !** molecules set to "INITIAL_SIZE"
    Call config_allocfields(sorbates(sorbtype), sorbtype, INITIAL_SIZE)

    !** Allocate a default amount of memory for the different
    !** arrays in sorbates.  Later on if the arrays are too
    !** small their size is incremented
    natoms = molecules_getnatoms(sorbtype)
    sorbates(sorbtype)%natoms = natoms
    sorbates(sorbtype)%nmoles = 0

  End Subroutine config_initdefaultsize

  !----------------------------------------------------------------------------
  ! The configuration is initialized using the atom coordinates that are 
  ! stored in the molecule file. This will be useful when using a flexible
  ! zeolites or other sorbant. Maybe useful for debugging. Who knows, some
  ! day maybe it can clean your dog.
  !----------------------------------------------------------------------------
  Subroutine config_initmolecule(sorbates,sorb)
    Implicit None
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Integer, Intent(IN)                            :: sorb

    Integer :: natoms, error, i
    Type(VecType), Dimension(:), Allocatable :: vecs

    !** Get the number of atoms in the molecule
    natoms = molecules_getnatoms(sorb)
    sorbates(sorb)%natoms = natoms

    !** Set the number of molecules to 1
    sorbates(sorb)%nmoles = 1

    !** Allocate the fields of the sorbate array based on the number
    !** of molecules and atoms.
    Call config_allocfields(sorbates(sorb),sorb,sorbates(sorb)%nmoles)

    !** Allocate our temporary position vector array. This will
    !** hold the atom coordinates obtained from molecules.
    Allocate(vecs(natoms),stat=error)
    If (error /=0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not allocate vecs"
      Stop
    End If

    !** Get the coordinates from the molecule
    Call molecules_getcomdefn(sorb,vecs)

    !** Copy the coodinates from the molecule atom array to 
    !** the new sorbates array
    Do i = 1,sorbates(sorb)%nmoles
      sorbates(sorb)%coords(1:natoms,i)%r = vecs
      sorbates(sorb)%coords(1:natoms,i)%rp = vecs
    End Do

    !** Deallocate the temporary vecs array
    Deallocate(vecs,stat=error)
    If (error/=0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not deallocate vecs"
      Stop
    End If

    !** Warn about this routine
    !MDEBUG
    Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
        " MOLECULE config type is not yet complete. Please check."

  End Subroutine config_initmolecule

  !----------------------------------------------------------------------------
  ! Initialize a FIXED type of sorbate configuration. This type
  ! obtains the coordinates from the molecule structure and is
  ! never moved.
  !----------------------------------------------------------------------------
  Subroutine config_initFixed(sorbates,sorb)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: sorbates
    Integer, Intent(In) :: sorb

    Integer :: natoms, error, i
    Type(VecType), Dimension(:), Allocatable :: vecs
    
    !** Get the number of atoms from molecule.
    natoms = molecules_getnatoms(sorb)
    sorbates(sorb)%natoms = natoms

    !** Allocate our temporary position vector array. This will
    !** hold the atom coordinates obtained from molecules.
    Allocate(vecs(natoms),stat=error)
    If (error /=0) Then
      Write(0,'(2a,i4,2a,i4)') __FILE__,":",__LINE__, &
          " Could not allocate vecs ", "natoms is-",natoms
      Stop
    End If

    !** Set the number of molecules, which is 1 for the time
    !** being. It should be whatever is required.
    sorbates(sorb)%nmoles = 1
    
    !** Allocate the sorbates array. This initializes all the
    !** fields of sorbates(sorb) like coords
    Call config_allocfields(sorbates(sorb),sorb,sorbates(sorb)%nmoles)

    !** Set the fixed flag to true. This means the sorbate is 
    !** never moved.
    sorbates(sorb)%fixed = .True.

    !** Get the coordinates from the molecules module
    !** and copy them into the sorbates array.
    Call molecules_getcomdefn(sorb,vecs)

    Do i = 1,sorbates(sorb)%nmoles
      sorbates(sorb)%coords(1:natoms,i)%r = vecs
      sorbates(sorb)%coords(1:natoms,i)%rp = vecs
    End Do

    !** Deallocate the temporary vecs array
    Deallocate(vecs,stat=error)
    If (error/=0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not deallocate vecs"
      Stop
    End If

  End Subroutine config_initFixed

  !--------------------------------------------------------
  ! Check if the initial configuration for all "sorbates"
  ! has been initialized
  !--------------------------------------------------------
  Logical Function config_checkinit(sorbates)
    Type(AtMolCoords), Dimension(:), Intent(In), Target :: sorbates

    Type(AtMolCoords), Dimension(:), Pointer :: sorbptr
    
    Integer     :: i

    sorbptr => sorbates
    If (.Not. Associated(sorbptr)) Then
      config_checkinit = .False.
      Return
    End If

    Do i = 1,Size(sorbates, 1)
      If (.Not. (sorbates(i)%isinit)) Then
        config_checkinit = .False.
        Return
      End If
    End Do

    config_checkinit = .True.

  End Function config_checkinit

  !----------------------------------------------------------------------------
  ! Creates a 3D grid of points inside the simulation cell simcell given
  ! the desired spacing gridspace for the sorbate sorb. It returns the 
  ! coordinates grid and the number of points in each direction npnts.
  !----------------------------------------------------------------------------
  Subroutine config_makegrid(simcell,sorb,gridspace,grid,npnts,fund)
    Type(SimCell_Params), Intent(In) :: simcell
    Real(Kind=RDbl), Intent(In)      :: gridspace
    Type(AtMolCoords), Intent(InOut) :: grid
    Integer, Dimension(3)            :: npnts
    Integer, Intent(In)              :: sorb
    Logical, Optional, Intent(In)    :: fund

    Real(Kind=RDbl), Dimension(3)    :: eff, step, coord
    Integer :: pnts, error
    Integer :: m, x, y, z
    Logical :: fundOnly

    !** Check to see if the optional arguement fund is passed.
    !** If fund is true, the grid is only created for the 
    !** fundamental cell of the simcell.
    fundOnly = .False.
    If (Present(fund)) fundOnly = fund
    
    !** Warn if this is used for a molecule with more than one
    !** atom
    If (molecules_getnatoms(sorb) > 1) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " WARNING: only creates coords for a 1 atom molecule"
      Stop
    End If
      
    !** Get the effective edge lengths of the simulation cell
    eff = simcell_geteff(simcell,fundOnly)

    !** Calculate the number of points
    npnts = Int(eff/gridspace) + 1
    pnts = npnts(1)*npnts(2)*npnts(3)

    !** The number of points is the number of molecules
    grid%nmoles = pnts

    !** The number of atoms in the molecule is the number of atoms here
    grid%natoms = molecules_getnatoms(sorb)

    !** Allocate the sorbates array by calling the config init routines
!    Call config_allocfields(grid,sorb,pnts,.True.)
    Allocate(grid%coords(1,pnts),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not allocate memory for grid%coords"
      Stop
    End If

    !** Figure out the step size in each direction
    step = eff/Real(npnts-1,Kind=RDbl)

    !** Loop through the points and store the coordinates in grid
    m = 0
    Do z = 1, npnts(3)
      coord(3) = Real(z-1,kind=RDbl)*step(3)
      Do y = 1, npnts(2)
        coord(2) = Real(y-1,kind=RDbl)*step(2)
        Do x = 1, npnts(1)
          coord(1) = Real(x-1,kind=RDbl)*step(1)
          m = m + 1
          grid%coords(1,m)%r = coord
          grid%coords(1,m)%rp = coord
        End Do
      End Do
    End Do

    !MDEBUG
    Write(0,*) "All coordinates generated!"

  End Subroutine config_makegrid

  !----------------------------------------------------------------------------
  ! Returns true if the sorbate coordinates are fixed and are not to be changed
  !----------------------------------------------------------------------------
  Logical Function config_isfixed(sorbate)
    Type(AtMolCoords), Intent(In) :: sorbate
    config_isfixed = sorbate%fixed
  End Function config_isfixed

  !-----------------------------------------------------------
  ! Set the no. of molecules of sorbate type "sorbtype" to
  ! "nmoles"
  !-----------------------------------------------------------
  Subroutine config_setnmoles_array(sorbates, sorbtype, nmoles)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Integer, Intent(in)    :: sorbtype, nmoles
    
    sorbates(sorbtype)%nmoles = nmoles
  End Subroutine config_setnmoles_array

  !-----------------------------------------------------------
  ! Set the no. of molecules of sorbate "nmoles"
  !-----------------------------------------------------------
  Subroutine config_setnmoles_single(sorbate,nmoles)
    Type(AtMolCoords), Intent(inout) :: sorbate
    Integer, Intent(in)    :: nmoles
    sorbate%nmoles = nmoles
  End Subroutine config_setnmoles_single

  !-------------------------------------------------------------------------
  ! Sets the generalized coordinates of a molecule in the configuration.
  ! It also takes the new generalized coordinates and uses them to generate
  ! the rp coordinates.  Note that the interpretation of the flat array of 
  ! new coordinates is dependent on the generalized coordinate model.
  ! Requires: species -- the species data structure
  !           molec -- the molecule number 
  !           spc -- the species number
  !           genarray -- an array of reals containing new coordinates
  !-------------------------------------------------------------------------
  Subroutine config_setfromgcoords(species,molec,spc,genarray)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,molec
    Real(kind=RDbl), Dimension(:), Intent(In)      :: genarray

    Call gcmodels_setgencoords(species(spc)%gcoords(molec),genarray)
    Call gcmodels_toxyz(species(spc)%gcoords(molec), &
        species(spc)%coords(:,molec)%rp)

  End Subroutine config_setfromgcoords

  !-------------------------------------------------------------------------
  ! Takes the generalized coordinates from one molecule of one species and
  ! uses that information to build a molecule of another species.
  ! Requires:  species1 -- 1st species data structure
  !            mol1 -- molecule number in 1st species
  !            species2 -- 1st species data structure
  !            mol2 -- molecule number in 2nd species to create
  !-------------------------------------------------------------------------
  Subroutine config_idrebuild(species1,mol1,species2,mol2)
    Type(AtMolCoords), Intent(In)    :: species1
    Integer, Intent(In)              :: mol1
    Type(AtMolCoords), Intent(InOut) :: species2
    Integer, Intent(In)              :: mol2

    Real(kind=RDbl), Dimension(MAX_DOF)      :: genarray

    !** Get the generalized coordinates of molecule 1
    genarray = gcmodels_getgencoords(species1%gcoords(mol1))

    ! as of now first 6 degrees of freedom are being used, for anything more 
    ! modify this part.
    !** Replace molecule 2's generalized coordinates and rebuild
    Call gcmodels_setgencoords(species2%gcoords(mol2),genarray(1:6))
    Call gcmodels_toxyz(species2%gcoords(mol2),species2%coords(:,mol2)%rp)

  End Subroutine config_idrebuild

  !-------------------------------------------------------------------------
  ! Changes the 'fixed' status of a species
  ! Requires:  species -- the species data structure
  !            spc -- the species number
  !            fixed -- new value for the fixed flag
  !-------------------------------------------------------------------------
  Subroutine config_changefixed(species,spc,fixed)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc
    Logical, Intent(In)                            :: fixed

    species(spc)%fixed = fixed

  End Subroutine config_changefixed

  !-------------------------------------------------------------------------
  ! Gets the number of fixed species and their species numbers
  ! Requires:  species -- the species data structure
  !            list -- the list of fixed species numbers
  !-------------------------------------------------------------------------
  Integer Function config_getfixed(species,list)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Dimension(:), Intent(Out)          :: list

    Integer       :: spc

    config_getfixed = 0
    Do spc = 1,Size(species)
      If (species(spc)%fixed) Then
        config_getfixed = config_getfixed + 1
        If (config_getfixed > Size(list)) Then
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              ' Passed integer array too small'
          Stop
        End If
        list(config_getfixed) = spc
      End If
    End Do

  End Function config_getfixed

  !-------------------------------------------------------------------------
  ! Changes the flexibility status of a molecule.  If no molecule number is
  ! passed than the flexibility of all the allocated generalized coordinates
  ! of a particularly species is changed.
  ! Requires: species -- the species data structure
  !           spc -- the species number
  !           flexflag -- new value for internal flexibility
  !           molec -- the molecule number, optional
  !-------------------------------------------------------------------------
  Subroutine config_changeflex(species,spc,flexflag,molec)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc
    Logical, Intent(In)                            :: flexflag
    Integer, Intent(In), Optional                  :: molec

    Integer           :: m

    If (Present(molec)) Then
      Call gcmodels_changeflex(species(spc)%gcoords(molec),flexflag)
    Else
      Do m = 1,Size(species(spc)%gcoords)
        Call gcmodels_changeflex(species(spc)%gcoords(m),flexflag)
      End Do
    End If

  End Subroutine config_changeflex

  !-------------------------------------------------------------------------
  ! This routine checks the size of various arrays against "nmoles" and 
  ! reallocates the arrays with a size which are INCR_SIZE bigger than 
  ! the previous size.  INCR_SIZE is defined in the defaults module.
  ! It will resize the arrays if the available storage space is less than
  ! or equal to the number of molecules specified.
  ! Requires: species -- the species data structure
  !           spc -- the species number
  !           nmoles -- current number of molecules
  !           optPreserve -- logical flag which specifies, whether old values 
  !                       should be preserved or not
  !-------------------------------------------------------------------------
  Subroutine config_checkandincr(species, spc, nmoles, optPreserve)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc, nmoles
    Logical, Intent(in), Optional                  :: optPreserve

    Logical :: preserve
    
    Integer                :: nsize, m, newsize
    Type(AtMolCoords)      :: sorbatetempptr
    
    preserve = .True.
    If (Present(optPreserve)) preserve = optPreserve

    !** Check if the species array is big enough
    !** Get the current number of molecules 
    nsize = Size(species(spc)%coords, 2)
    
    !** Check if we need to reallocate stuff
    If (nmoles < nsize) Then
      Return
    Else If (nmoles >= nsize) Then
      newsize = nsize + INCR_SIZE

      If (nmoles > newsize) Then
        Write(*,*) "nmoles too large, try increasing INCR_SIZE in defaults "
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End If

      !** HERE here HERE
      sorbatetempptr = species(spc)
      Nullify(species(spc)%coords)
      Nullify(species(spc)%gcoords)
      Nullify(species(spc)%afast)
      Nullify(species(spc)%aslow)

      !** Reallocate the fields (does not set %nmoles)
      Call config_allocfields(species(spc), spc, newsize)

      !** Copy old values if desired
      If (preserve) Then
        Do m = 1,nmoles
          If (m <= nsize) Then
            Call config_copymolec(sorbatetempptr, m, species(spc), m)
          Else
            Write(*,*) "---- Warning -------"
            Write(*,*) "Can not copy molecule- ", m,&
                "There was no space allocated for it before "
            Write(*,*) "Probably this should have been called with"
            Write(*,*) "preserve=.False."
          End If
        End Do
      End If

    End If

  End Subroutine config_checkandincr

  !-------------------------------------------------------------------
  ! Copies a specified subset of the configuration from one structure
  ! to another.
  ! Requires:  orig -- origin configuration structure
  !            dest -- destination configuration structure
  !            osubset -- subset to copy from origin
  !            dsubset -- subset to place information in destination
  !-------------------------------------------------------------------
  Subroutine config_copysubset(orig,dest,osubset,dsubset)
    Type(AtMolCoords), Dimension(:), Intent(In)    :: orig
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: dest
    Integer, Dimension(3), Intent(In)              :: osubset,dsubset

    Integer    :: spc,mol,atm,dspc,dmol,datm
    Integer    :: nspcs,nmols,natms,depth

    !** Check that the subset depths match
    depth = getdepth(osubset)
    If (depth /= getdepth(dsubset)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Subset depths do not match, cannot copy'
      Stop
    End If

    !** Loop over the origin species number and copy by species
    nspcs = Size(orig)
    Do spc = 1,nspcs
      !** Skip this species if not pertinent
      If ((osubset(1) /= 0).And.(osubset(1) /= spc)) Cycle
      dspc = dsubset(1)
      If (dspc == 0) dspc = spc

      !** Get number of molecules and atoms in the origin species
      nmols = orig(spc)%nmoles
      natms = orig(spc)%natoms

      If (depth > 1) Then   !** copy molecules, atoms individually
        Do mol = 1,nmols
          !** Skip this molecule if not pertinent
          If ((osubset(2) /= 0).And.(osubset(2) /= mol)) Cycle
          dmol = dsubset(2)
          If (dmol == 0) dmol = mol

          !** Check number of atoms if we're copying whole atoms
          If ((osubset(3) == 0).And.(dest(dspc)%natoms /= natms)) Then
            Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
                ' number of atoms does not match, cannot copy'
            Stop
          End If

          Do atm = 1,natms
            !** Skip this atom if not pertinent and get destination atm number
            If ((osubset(3) /= 0).And.(osubset(3) /= atm)) Cycle
            datm = dsubset(3)
            If (datm == 0) datm = atm

            dest(dspc)%coords(datm,dmol) = orig(spc)%coords(atm,mol)
            dest(dspc)%afast(datm,dmol) = orig(spc)%afast(atm,mol)
            dest(dspc)%aslow(datm,dmol) = orig(spc)%aslow(atm,mol)
          End Do  !** loop over atoms

          !** Copy generalized coordinates for molecule if possible
          If (depth < 3) Then
            dest(dspc)%gcoords(dmol) = orig(spc)%gcoords(mol)
          Else
            Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
                ' WARNING: cannot copy portions of molecular generalized coords'
            Stop
          End If

        End Do  !** loop over molecules

      Else   !** copy the whole species in bulk to speed process
        !** Set nmoles at destination to be the same as the origin
        dest(dspc)%nmoles = nmols

        !** Make sure there's enough space for the molecules
        If (Size(dest(dspc)%coords,2) < nmols) Then
          Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              ' number of molecules does not match, cannot copy'
          Write(0,'(2x,a,i4)') 'destination nmolecules = ',dest(dspc)%nmoles
          Write(0,'(2x,a,i4)') 'origin nmolecules = ',nmols
          Write(0,'(2x,a,3i4)') 'destination subset = ',dsubset
          Write(0,'(2x,a,3i4)') 'origin subset = ',osubset
          Stop
        End If

        !** check equivalence of the number of atoms
        If (dest(dspc)%natoms /= natms) Then
          Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              ' number of atoms does not match, cannot copy'
          Write(0,'(2x,a,i4)') 'destination nmolecules = ',dest(dspc)%natoms
          Write(0,'(2x,a,i4)') 'origin nmolecules = ',natms
          Stop
        End If

        !** Copy all fields, using operator overloading for gcoords
        dest(dspc)%coords(1:natms,1:nmols) = orig(spc)%coords(1:natms,1:nmols)
        dest(dspc)%afast(1:natms,1:nmols) = orig(spc)%afast(1:natms,1:nmols) 
        dest(dspc)%aslow(1:natms,1:nmols) = orig(spc)%aslow(1:natms,1:nmols) 
        Do mol = 1,nmols
          dest(dspc)%gcoords(mol) = orig(spc)%gcoords(mol)
        End Do
      End If

    End Do  !** loop over species

  End Subroutine config_copysubset

#ifdef OLDER_COPYSUBSET   !** remove soon
  !-------------------------------------------------------------------
  ! Copies the coordinates of from_conf to to_conf, only those parts 
  ! specified by the from_subset are copied
  ! -- Very ugly looking routine -- 
  ! NOTE : please notice that coords(index1, index2) are indexed so 
  !        that index1 is atom number and index2 is molecule number
  !-------------------------------------------------------------------
  Subroutine config_copysubset(from_conf , to_conf, from_subset, to_subset)
    Type(AtMolCoords), Dimension(:), Intent(In)    :: from_conf 
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: to_conf
    Integer, Dimension(3), Intent(in) :: from_subset
    Integer, Dimension(3), Intent(in) :: to_subset

    Integer    :: nmoles, natoms, nspc, to_m, to_spc, m, i

    nspc= Size(from_conf)

    If ( (from_subset(3)/=0) .Or. (to_subset(3)/=0 ) ) Then
      Write(*,*) " No atoms specific copying yet"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    !** more checks
    If ( (from_subset(1)<0) .Or. (from_subset(1)>nspc) ) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    !** loop thru all species and check where copying is necessary
    Do i=1,nspc

      If ( (from_subset(1)==0) .Or. (from_subset(1)==i) ) Then

        !** copy this species
        nmoles=from_conf(i)%nmoles
        natoms=from_conf(i)%natoms

        !** find the corresponding to_spc
        If (to_subset(1)/=0) Then
          If (from_subset(1)==0) Then
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
          End If
          to_spc=to_subset(1)
        Else
          to_spc=i
        End If

        If (from_subset(2)==0)  Then
          !** Copy all the molecules
          to_conf(to_spc)%nmoles = nmoles
          to_conf(to_spc)%natoms = natoms

          !** check array sizes again, see NOTE in header
          If ( (Size(to_conf(to_spc)%coords,2)<nmoles) .Or. &
              (Size(to_conf(to_spc)%coords,1)/=natoms ) ) Then
            Write(*,*) "sizes dont match in copysubset"
            Write(*,*) "nmoles-to : ", Size(to_conf(to_spc)%coords,2), &
                "nmoles from : " ,nmoles
            Write(*,*) "natoms-to : ", Size(to_conf(to_spc)%coords,1), &
                "natoms from : " ,natoms
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
          Else
            to_conf(to_spc)%coords(1:natoms,1:nmoles) = &
                from_conf(i)%coords(1:natoms, 1:nmoles)
            to_conf(to_spc)%afast(1:natoms,1:nmoles) = &
                from_conf(i)%afast(1:natoms, 1:nmoles)
            to_conf(to_spc)%aslow(1:natoms,1:nmoles) = &
                from_conf(i)%aslow(1:natoms, 1:nmoles)
            Do m=1,nmoles
              !** Operator overloading here, all fields of gcoords are copied
              to_conf(to_spc)%gcoords(m) = from_conf(i)%gcoords(m) 
            End Do
          End If

        Else
          !** need to copy only one molecule
          to_m=to_subset(2)

          !** see note in header
          If ( (Size(to_conf(to_spc)%coords,2)<to_m) .Or. &
              (Size(to_conf(to_spc)%coords,1)/=natoms ) ) Then
            Write(*,*) from_subset, to_subset, natoms, to_spc, to_m
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop

          Else
            ! ** we have to copy only one molecule
            to_conf(to_spc)%coords(1:natoms,to_m) = &
                from_conf(i)%coords(1:natoms, from_subset(2) )
            to_conf(to_spc)%afast(1:natoms,to_m) = &
                from_conf(i)%afast(1:natoms, from_subset(2) )
            to_conf(to_spc)%aslow(1:natoms,to_m) = &
                from_conf(i)%aslow(1:natoms, from_subset(2) )

            !** Operator overloading here, all fields of gcoords are copied
            to_conf(to_spc)%gcoords(to_m) = &
                from_conf(i)%gcoords(from_subset(2)) 
          End If

        End If !** end of If subset(2)==0 
      End If !** end of 'If' for checking whether this species is to be copied
    End Do !** end of loop over spc

  End Subroutine config_copysubset
#endif

  !-------------------------------------------------------------------
  ! Rotates the coordinates of molec by given angles around ORIGIN
  ! If (rotateVels) Then rotates the velocities also
  !-------------------------------------------------------------------
  Subroutine config_rotatexyz(coords, molec, theta, phi, psi, rotateVels)
    Type(AtMolCoords), Intent(InOut) :: coords
    Integer, Intent(In)              :: molec
    Real(kind=RDbl), Intent(In)      :: theta, phi, psi 
    Logical, Intent(in)              :: rotateVels

    Integer            :: natoms, a
    Type(MatrixType)   :: b2f_eulermat

    !** Get the transformation matrix
    b2f_eulermat = matrix_b2feuler(phi, theta, psi)
    
    !** Check the no. of atoms 
    natoms = config_getnatoms(coords)

    If (rotateVels) Then
      Do a = 1,natoms
        coords%coords(a,molec)%rp = b2f_eulermat * coords%coords(a,molec)%rp 
        coords%coords(a,molec)%v = b2f_eulermat * coords%coords(a,molec)%v
      End Do
    Else
      Do a = 1,natoms
        coords%coords(a,molec)%rp = b2f_eulermat * coords%coords(a,molec)%rp 
      End Do
    End If

  End Subroutine config_rotatexyz

  !-------------------------------------------------------------------
  ! Rotates the coordinates of molec around the unit vector given by phi, theta
  ! The amount of rotation is psi 
  ! If (rotateVels) Then rotates the velocities also
  !-------------------------------------------------------------------
  Subroutine config_rotateAroundCOM( sorbates, sorb , molec, theta, phi, &
      psi, rotateVels)
    Type(AtMolCoords), Dimension(:), Intent(Inout)    :: sorbates
    Integer, Intent(In)              :: sorb, molec
    Real(kind=RDbl), Intent(In)      :: theta, phi, psi 
    Logical, Intent(in)              :: rotateVels

    Type(VecType) :: com, tempvec
    Type(MatrixType)   :: transMatrix
    Integer    :: natoms, a
    
    ! find COM
    com=config_getMolecCOM(sorbates, sorb, molec)

    !** Get the transformation matrix
    transMatrix = matrix_rotAroundVec(phi, theta, psi)
    
    !** Check the no. of atoms 
    natoms = config_getnatoms(sorbates,sorb)

    If (rotateVels) Then
      Do a =1,natoms
        tempvec = sorbates(sorb)%coords(a,molec)%rp - com
        sorbates(sorb)%coords(a,molec)%rp = com + transMatrix * tempvec
        sorbates(sorb)%coords(a,molec)%v = transMatrix * &
            sorbates(sorb)%coords(a,molec)%v
      End Do
    Else
      Do a =1,natoms
        tempvec = sorbates(sorb)%coords(a,molec)%rp - com
        sorbates(sorb)%coords(a,molec)%rp = com + transMatrix * tempvec
       End Do
    Endif

  End Subroutine config_rotateAroundCOM



  !-------------------------------------------------------------------
  ! Reverses the velocities of a given molecule 
  !-------------------------------------------------------------------
  Subroutine config_reverseVels( sorbates, sorb , molec)
    Type(AtMolCoords), Dimension(:), Intent(Inout)    :: sorbates
    Integer, Intent(In)              :: sorb, molec

    Integer    :: natoms, a

    !** Check the no. of atoms 
    natoms = config_getnatoms(sorbates,sorb)
    Do a =1,natoms
      sorbates(sorb)%coords(a,molec)%v = (-one) * &
          sorbates(sorb)%coords(a,molec)%v
    End Do

  End Subroutine config_reverseVels


  !-------------------------------------------------------------------
  ! Copies the coordinates of "molec1" in "coords1" into the "molec2"
  ! location of "coords2"
  ! Requires:  coords1 -- full species storage with molecule to copy
  !            molec1 -- 1st molecule number
  !            coords2 -- full species storage to copy molecule into
  !            molec2 -- 2nd molecule number
  !            opt_copy_a -- do not copy accelerations also if False
  !-------------------------------------------------------------------
  Subroutine config_copymolec(coords1,molec1,coords2,molec2,opt_copy_a)
    Type(AtMolCoords), Intent(In)    :: coords1
    Type(AtMolCoords), Intent(InOut) :: coords2
    Integer, Intent(In)              :: molec1, molec2
    Logical, Intent(in), Optional    :: opt_copy_a
    
    Integer    :: natoms
    Logical    :: copy_accels

    natoms = coords1%natoms
    coords2%natoms = natoms
    copy_accels = .True.
    If (Present(opt_copy_a)) copy_accels = opt_copy_a 

    !** Copy the coordinates
    coords2%coords(1:natoms,molec2) = coords1%coords(1:natoms,molec1)

    !** Use operator overloading to copy generalized coordinates
    coords2%gcoords(molec2) = coords1%gcoords(molec1)

    !** Copy accelerations if desired
    If (copy_accels) Then
      coords2%afast(1:natoms, molec2) = coords1%afast(1:natoms, molec1)
      coords2%aslow(1:natoms, molec2) = coords1%aslow(1:natoms, molec1)
    End If

  End Subroutine config_copymolec

#ifdef OBSOLETE
  !-----------------------------------------------------------------------
  ! This routine deletes molecule "molec" of type "sorbno" and compresses
  ! the list afterward.
  !-----------------------------------------------------------------------
  Subroutine config_delmolec(sorbates, sorbno, molecno)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Integer, Intent(in)  :: sorbno, molecno
    Integer     :: nmoles, natoms, m

    nmoles = sorbates(sorbno)%nmoles
    natoms = sorbates(sorbno)%natoms

    !** Copy the xyz coords
    Do m = molecno, nmoles-1
      Call config_copymolec(sorbates(sorbno), m+1, sorbates(sorbno), m)
    End Do

    nmoles = nmoles - 1
    sorbates(sorbno)%nmoles = nmoles
  End Subroutine config_delmolec
#endif

  !-----------------------------------------------------------------------
  ! This routine deletes one molecule of the specified species type.  It
  ! Then moves the last molecule in the list to the empty position.
  ! Requires:  species -- full system configuration
  !            spc -- species number
  !            mol -- molecule number
  !-----------------------------------------------------------------------
  Subroutine config_delmol(species,spc,mol)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol

    Integer     :: nmoles

    !** Copy the molecule at the end of the list into the vacant position
    nmoles = species(spc)%nmoles
    Call config_copymolec(species(spc),nmoles,species(spc),mol)

    !** Decrement the count to effectively erase the last molecule
    species(spc)%nmoles = nmoles - 1

  End Subroutine config_delmol

  !----------------------------------------------------------------------------
  ! This function returns the configuration info of type sorb
  !----------------------------------------------------------------------------
  Subroutine config_getconfig(sorbates,sorb,sorbpntr)
    Type(AtMolCoords), Intent(In), Dimension(:), Target :: sorbates
    Integer, Intent(In)                                 :: sorb
    Type(AtMolCoords), Pointer                          :: sorbpntr

    sorbpntr => sorbates(sorb)

  End Subroutine config_getconfig

  !---------------------------------------------------------------------
  ! This function returns the source where the configuration information
  ! for the sorbate type "sorbtype" was obtained from
  !---------------------------------------------------------------------
!!$  Character(len=strLen) Function config_getsourcetype(sorbates, sorbtype)
  Function config_getsourcetype(sorbates, sorbtype)
    Character(len=strLen) :: config_getsourcetype
    Type(AtMolCoords), Intent(In), Dimension(:) :: sorbates
    Integer, Intent(In) :: sorbtype
    
    config_getsourcetype = sorbates(sorbtype)%sourcetype
  End Function config_getsourcetype

  !---------------------------------------------------------------------
  ! Returns the filename for reading the initial configuration of
  ! "sorbtype" as stored in "sorbates"
  !---------------------------------------------------------------------
!!$  Character(len=strLen) Function config_getsourcefile(sorbates, sorbtype)
  Function config_getsourcefile(sorbates, sorbtype)
    Character(len=strLen) :: config_getsourcefile
    Type(AtMolCoords), Intent(In), Dimension(:) :: sorbates
    Integer, Intent(In) :: sorbtype
    
    config_getsourcefile = sorbates(sorbtype)%filename
  End Function config_getsourcefile

  !----------------------------------------------------------------------------
  ! This Subroutine returns a pointer to the coordinates of sorbate sorb
  !----------------------------------------------------------------------------
  Subroutine config_getcoord(sorbates,sorb,pntr)
    Type(AtMolCoords), Dimension(:), Intent(In), Target :: sorbates
    Type(RectCoords), Dimension(:,:), Pointer :: pntr
    Integer, Intent(In) :: sorb

    pntr => sorbates(sorb)%coords
  End Subroutine config_getcoord

  !----------------------------------------------------------------------------
  ! Apply periodic boundary conditions to the parent coordinates to regenerate
  ! the simulation cell coordinates and the continuation vectors.
  ! Requires:  species -- species data structure
  !            simcell -- simulation cell information
  !            subset -- system subset to which to apply PBCs
  ! NOTE: I'm worried that this may slow things down unnecessarily
  !----------------------------------------------------------------------------
  Subroutine config_applypbcs(species,simcell,subset)
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Type(SimCell_Params), Intent(In)                 :: simcell
    Integer, Dimension(3), Intent(In)                :: subset

    Integer             :: depth,natoms,nmoles,mol

    depth = getdepth(subset)

    If (depth == 3) Then
      Call simcell_pbc(simcell, &
          species(subset(1))%coords(subset(3),subset(2))%rp, &
          species(subset(1))%coords(subset(3),subset(2))%r, &
          species(subset(1))%coords(subset(3),subset(2))%cr)    

    Else If (depth == 2) Then
      natoms = molecules_getnatoms(subset(1))
      Call simcell_pbc(simcell, &
          species(subset(1))%coords(1:natoms,subset(2))%rp, &
          species(subset(1))%coords(1:natoms,subset(2))%r, &
          species(subset(1))%coords(1:natoms,subset(2))%cr)    

    Else If (depth == 1) Then
      natoms = molecules_getnatoms(subset(1))
      nmoles = config_getnmoles(species,subset(1))
      Do mol = 1,nmoles
        Call simcell_pbc(simcell, &
            species(subset(1))%coords(1:natoms,mol)%rp, &
            species(subset(1))%coords(1:natoms,mol)%r, &
            species(subset(1))%coords(1:natoms,mol)%cr)    
      End Do

    Else
      Write(0,'(1x,2a,i4,a,i2)') __FILE__," : ",__LINE__, &
          ' not prepared for this depth ',depth
      Stop
    End If

  End Subroutine config_applypbcs

  !----------------------------------------------------------------------------
  ! This function returns the principal coordinates of an atom
  ! Requires:  species -- full species data structure
  !            subset -- atom specification in subset form
  !----------------------------------------------------------------------------
  Type(VecType) Function config_getrp(species,subset)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Dimension(3), Intent(In)           :: subset
     
    config_getrp = species(subset(1))%coords(subset(3),subset(2))%rp

  End Function config_getrp

  !----------------------------------------------------------------------------
  ! This function returns the simulation cell coordinates of an atom
  ! Requires:  species -- full species data structure
  !            subset -- atom specification in subset form
  !----------------------------------------------------------------------------
  Type(VecType) Function config_getr(species,subset)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Dimension(3), Intent(In)           :: subset
     
    config_getr = species(subset(1))%coords(subset(3),subset(2))%r

  End Function config_getr

  !----------------------------------------------------------------------------
  ! Returns the sim cell coordinates of a molecule
  ! Requires:  species -- full configuration information
  !            spc -- species number
  !            mol -- molecule number
  !            coords -- returned coordinates
  !----------------------------------------------------------------------------
  Subroutine config_getmolr(species,spc,mol,coords)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: spc,mol
    Type(VecType), Dimension(:), Intent(Out)    :: coords
     
    coords = species(spc)%coords(:,mol)%r

  End Subroutine config_getmolr

  !----------------------------------------------------------------------------
  ! Returns the charges of one atom in a molecule
  ! Requires:  species -- full configuration information
  !            spc -- species number
  !            mol -- molecule number
  !            atm -- atom number
  !            charge -- returned charge
  ! NOTE: will need to be changed when we start storing charges in species
  !       does not use species (that's why it doesn't use 'mol')
  !----------------------------------------------------------------------------
  Subroutine config_getq(species,spc,mol,atm,charge)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: spc,mol,atm
    Real(kind=RDbl), Intent(Out)                :: charge
     
    charge = molecules_getcharge(spc,atm)

  End Subroutine config_getq

  !----------------------------------------------------------------------------
  ! Returns the charges of each atom in a molecule
  ! Requires:  species -- full configuration information
  !            spc -- species number
  !            mol -- molecule number
  !            charges -- returned charges on each atom
  ! NOTE: will need to be changed when we start storing charges in species
  !       does not use species (that's why it doesn't use 'mol')
  !----------------------------------------------------------------------------
  Subroutine config_getmolq(species,spc,mol,charges)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: spc,mol
    Real(kind=RDbl), Dimension(:), Intent(Out)  :: charges
     
    Call molecules_getcharges(spc,charges)

  End Subroutine config_getmolq

  !----------------------------------------------------------------------------
  ! Returns a pointer to the velocities for sorbate sorb
  !----------------------------------------------------------------------------
  Subroutine config_getv(sorbates,sorb,vpntr)
    Type(AtMolCoords), Dimension(:), Intent(In), Target :: sorbates
    Integer, Intent(In) :: sorb
    Type(VecType), Dimension(:,:), Pointer :: vpntr
    Integer :: natoms, nmoles

    If (Associated(sorbates(sorb)%coords)) Then
      natoms = config_getnatoms(sorbates,sorb)
      nmoles = config_getnmoles(sorbates,sorb)
      vpntr => sorbates(sorb)%coords(1:natoms,1:nmoles)%v
    Else
      nullify(vpntr)
    End If
  End Subroutine config_getv

  !----------------------------------------------------------------------------
  ! This function returns the type of an atom
  ! Requires: species -- species data structure
  !           spc -- species number
  !           molec -- molecule number
  !           atom -- atom number
  ! NOTE: doesn't current use 'molec', but may need to later
  !       doesn't current use species
  !----------------------------------------------------------------------------
  Integer Function config_getatype(species, spc, molec, atom)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: spc, molec, atom

    config_getatype = molecules_getatype(spc, atom)
  End Function config_getatype

  !----------------------------------------------------
  ! This function returns the no. of atoms of sorbate
  ! "sorbno"
  !----------------------------------------------------
  Integer Function config_getnatoms_array(sorbates, sorbno)
    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Integer, Intent(in)  :: sorbno

    config_getnatoms_array = sorbates(sorbno)%natoms
  End Function config_getnatoms_array

  !----------------------------------------------------
  ! This function returns the no. of atoms of sorbate
  ! "sorbno"
  !----------------------------------------------------
  Integer Function config_getnatoms_single(sorbate)
    Type(AtMolCoords), Intent(In) :: sorbate

    config_getnatoms_single = sorbate%natoms

  End Function config_getnatoms_single

  !----------------------------------------------------------------------------
  ! Counts the number of atoms in a given system subset
  ! Requires:  species -- species data structure
  !            subset -- subset of system coordinates (spc,mol,atm)
  !----------------------------------------------------------------------------
  Integer Function config_nsubsetatoms(species,subset)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Dimension(3), Intent(In)           :: subset

    Integer         :: spc,natoms

    config_nsubsetatoms = 0

    !** Check to see if it's just a single atom
    If (subset(3) /= 0) Then
      config_nsubsetatoms = 1
      Return
    End If

    Do spc = 1,Size(species)
      If ((subset(1) /= 0).And.(subset(1) /= spc)) Cycle

      natoms = molecules_getnatoms(spc)
      If (subset(2) == 0) Then  !** all molecules
        config_nsubsetatoms = config_nsubsetatoms + natoms*species(spc)%nmoles
      Else                      !** just one molecule
        config_nsubsetatoms = config_nsubsetatoms + natoms
      End If
    End Do

  End Function config_nsubsetatoms

  !---------------------------------------------------------------
  ! This routine returns a list of the number of molecules in the
  ! full species structure
  ! Requires:  species -- full species structure
  !---------------------------------------------------------------
  Subroutine config_getnmoleslist(species,list)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Dimension(:), Intent(Out)          :: list

    Integer          :: spc

    Do spc = 1,Size(species)
      list(spc) = species(spc)%nmoles
    End Do

  End Subroutine config_getnmoleslist

  !----------------------------------------------------
  ! This function returns the no. of molecules of sorbate
  ! "sorbno"
  !----------------------------------------------------
  Integer Function config_getnmoles_array(sorbates, sorbno)
    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Integer, Intent(in)  :: sorbno

    config_getnmoles_array = sorbates(sorbno)%nmoles
  End Function config_getnmoles_array

  !----------------------------------------------------
  ! This function returns the total no. of molecules in the system, 
  ! excluding the sorbates that are fixed 
  !----------------------------------------------------
  Integer Function config_getTotMoles(sorbates )
    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Integer :: sorb
    config_getTotMoles=0
    Do sorb=1,Size(sorbates,1)
      If (.Not.(config_isfixed(sorbates(sorb)))) &
          config_getTotMoles=    &
          config_getTotMoles+config_getnmoles(sorbates,sorb)
    End Do
  End Function config_getTotMoles

  !----------------------------------------------------
  ! This function returns the no. of molecules of the
  ! sorbate object sorbate passed
  !----------------------------------------------------
  Integer Function config_getnmoles_single(sorbate)
    Type(AtMolCoords), Intent(in) :: sorbate

    config_getnmoles_single = sorbate%nmoles
  End Function config_getnmoles_single

  !-----------------------------------------------------------------------
  ! This routine gets the centers of mass of all molecules of one species
  ! using the PRINCIPAL coordinates.
  ! Requires:  species -- species data structure
  !            spc -- species number
  !            com -- output array of COM vectors
  !-----------------------------------------------------------------------
  Subroutine config_getSorbCOM(species,spc,com)
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: spc
    Type(VecType), Intent(Out), Dimension(:)     :: com

    Integer    :: nmoles, m

    nmoles = species(spc)%nmoles
    Do m = 1,nmoles
      com(m) = config_getMolecCOM(species,spc,m)
    End Do

  End Subroutine config_getSorbCOM

  !-----------------------------------------------------------------------
  ! Calculate the center of mass (COM) of a single molecule using the 
  ! ** PRINCIPAL ** coordinates.
  ! Requires:  species -- species data structure
  !            spc -- species number
  !            molec -- molecule number
  !-----------------------------------------------------------------------
  Type(VecType) Function config_getMolecCOM(species,spc,molec)
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: spc,molec

    Integer          :: a,natoms
    Real(kind=RDbl)  :: mass,mass_sum

    natoms = species(spc)%natoms
    config_getMolecCOM = 0.0_RDbl
    mass_sum = molecules_getmass(spc)

    Do a = 1,natoms
      mass = atom_getmass(molecules_getnthatom(spc,a))
      config_getMolecCOM = config_getMolecCOM + &
          (mass*species(spc)%coords(a,molec)%rp)
    End Do
    config_getMolecCOM = config_getMolecCOM/mass_sum

  End Function config_getMolecCOM

  !-----------------------------------------------------------------------
  ! Calculate the geometric center of a single molecule using the 
  ! PRINCIPAL coordinates.
  ! Requires:  species -- species data structure
  !            spc -- species number
  !            molec -- molecule number
  !-----------------------------------------------------------------------
  Type(VecType) Function config_getMolecCTR(species,spc,molec)
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: spc,molec

    Integer          :: a,natoms
    Type(VecType)    :: ctr
    
    natoms = species(spc)%natoms
    ctr = 0.0_RDbl
    Do a = 1, natoms
      ctr = ctr + species(spc)%coords(a,molec)%rp
    End Do
    config_getMolecCTR = ctr/natoms

  End Function config_getMolecCTR


  !-----------------------------------------------------------------------
  ! Calculate the velocity of center of mass , of one molecule
  ! Requires:  species -- species data structure
  !            spc -- species number
  !            molec -- molecule number
  !-----------------------------------------------------------------------
  Type(VecType) Function config_getCOMVel(species,spc,molec)
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: spc,molec

    Integer          :: a,natoms
    Real(kind=RDbl)  :: mass,mass_sum
    Type(VecType)    :: comvel


    natoms = species(spc)%natoms
    config_getCOMVel = zero
    mass_sum = molecules_getmass(spc)

    Do a = 1,natoms
      mass = atom_getmass(molecules_getnthatom(spc,a))
      config_getCOMVel = config_getCOMVel + &
          (mass*species(spc)%coords(a,molec)%v)
    End Do
    config_getCOMVel = config_getCOMVel/mass_sum


  End Function config_getCOMVel

  !-----------------------------------------------------------------------
  ! Calculate the velocity of center of mass , of one species.
  ! This is the instantaneous average velocity of that species
  ! Requires:  species -- species data structure
  !            spc -- species number
  !-----------------------------------------------------------------------
  Type(VecType) Function config_getSpcCOMVel(species,spc)
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: spc

    Integer          :: a,natoms,m,nmoles
    Real(kind=RDbl)  :: mass,mass_sum
    Type(VecType)    :: velsum


    natoms = species(spc)%natoms
    nmoles = config_getnmoles(species,spc)

    velsum = zero


    ! calculate \sigma m_i v_i where i runs over each atom of species
    Do a = 1,natoms
      mass = molecules_AtomMass(spc,a)
      Do m=1,nmoles
        velsum = velsum + (mass*species(spc)%coords(a,m)%v)
      End Do
    End Do

    mass_sum = nmoles * (molecules_getmass(spc))
    config_getSpcCOMVel = velsum/mass_sum

  End Function config_getSpcCOMVel


  !-----------------------------------------------------------------------
  ! Calculate the Density of one species.
  ! in molecules/ang.^3
  ! Requires:  species -- species data structure
  !            scell -- simcell for volume calcn
  !            spc -- species number
  !-----------------------------------------------------------------------
  Real(Kind=RDbl)Function config_getSpcDensity(species,scell,spc)
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)   :: scell
    Integer, Intent(In)                          :: spc
    config_getSpcDensity=config_getnmoles(species,spc) &
        /simcell_getvolume(scell)
  End Function config_getSpcDensity

  !----------------------------------------------------------
  ! Zero the "fast" or "slow" forces of the sorbate "sorbate"
  ! depending on which is true
  !----------------------------------------------------------
  Subroutine config_zeroforces(sorbate, fast)
    Type(AtMolCoords), Intent(inout) :: sorbate
    Logical, Intent(in)              :: fast

    Integer   :: natoms, nmoles, a, m
    natoms = sorbate%natoms
    nmoles = sorbate%nmoles

    If (fast) Then
      Do a=1, natoms
        Do m=1, nmoles
          sorbate%afast(a,m) = 0.0_RDbl
        EndDo
      End Do
    Else
      Do a=1, natoms
        Do m=1, nmoles
          sorbate%aslow(a,m) = 0.0_RDbl
        EndDo
      End Do
    End If
  End Subroutine config_zeroforces

  !----------------------------------------------------------
  ! This routine takes the array of accelerations which are
  ! stored in units of ang/ps^2 and converts them to forces
  ! in the units of kcal/ang.
  !----------------------------------------------------------
  Subroutine config_conv2force(sorbate, sorbtype, fast)
    Type(AtMolCoords), Intent(inout) :: sorbate
    Integer, Intent(in)              :: sorbtype
    Logical, Intent(in)              :: fast

    Integer   :: natoms, nmoles, a, m
    Real(kind=RDbl) :: atom_mass

    natoms = sorbate%natoms
    nmoles = sorbate%nmoles

    If (fast) Then
      Do m=1, nmoles
        Do a=1, natoms
          atom_mass = atom_getmass(molecules_getatype(sorbtype, a))
          sorbate%afast(a,m) = sorbate%afast(a,m)*atom_mass/scalef
        EndDo
      End Do
    Else
      Do m=1, nmoles
        Do a=1, natoms
          atom_mass = atom_getmass(molecules_getatype(sorbtype, a))
          sorbate%aslow(a,m) = sorbate%aslow(a,m)*atom_mass/scalef
        EndDo
      End Do
    End If
  End Subroutine config_conv2force

  !----------------------------------------------------------------------------
  ! This routine converts forces in kcal/mol/Ang in the species structure to 
  ! accelerations stored in units of ang/ps^2.
  ! Requires:  species -- single species data structure 
  !            spc -- species number
  !            fast -- True => modify Fast interactions, False => only slow
  !----------------------------------------------------------------------------
  Subroutine config_conv2accel(species,spc,fast)
    Type(AtMolCoords), Intent(InOut) :: species
    Integer, Intent(In)              :: spc
    Logical, Intent(In)              :: fast

    Integer                   :: natoms,nmoles,a,m
    Real(kind=RDbl), Dimension(species%natoms) :: invmass

    !** get number of molecules and inverse masses for each atom in species
    nmoles = species%nmoles
    natoms = species%natoms
    Call molecules_getinvmasses(spc,invmass)

    !** Do the conversion
    If (fast) Then
      Do m = 1,nmoles
        Do a = 1,natoms
          species%afast(a,m) = species%afast(a,m)*invmass(a)*scalef
        EndDo
      End Do
    Else
      Do m = 1,nmoles
        Do a = 1,natoms
          species%aslow(a,m) = species%aslow(a,m)*invmass(a)*scalef
        EndDo
      End Do
    End If

  End Subroutine config_conv2accel

  !----------------------------------------------------------------------------
  ! Get a pointer to the correct 2D array of accelerations.  Can be used to
  ! update the accelerations.
  ! Requires:  species -- species data structure 
  !            spc -- species number
  !            fast -- True => treat Fast interactions, False => treat slow
  !----------------------------------------------------------------------------
  Function config_accelptr(species,spc,fast)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: spc
    Logical, Intent(In)                         :: fast
    Type(VecType), Dimension(:,:), Pointer      :: config_accelptr

    If ((.Not.Associated(species(spc)%afast)).Or. &
        (.Not.Associated(species(spc)%aslow))) Then
      Nullify(config_accelptr)
      Return
    End If

    If (fast) Then
      config_accelptr => species(spc)%afast
    Else
      config_accelptr => species(spc)%aslow
    End If

  End Function config_accelptr

  !----------------------------------------------------------------------------
  ! Get a pointer to the correct 2D array of accelerations.  Can be used to
  ! update the accelerations.  Same as config_accelptr, but a subroutine
  ! Requires:  species -- species data structure 
  !            spc -- species number
  !            fast -- True => treat Fast interactions, False => treat slow
  !----------------------------------------------------------------------------
  Subroutine config_getaccel(sorbates,sorb,fast,accel)
    Type(AtMolCoords), Intent(In), Dimension(:) :: sorbates
    Integer, Intent(In)                         :: sorb
    Logical, Intent(In)                         :: fast
    Type(VecType), Dimension(:,:), Pointer      :: accel
    
    If ((.Not.Associated(sorbates(sorb)%afast)).Or. &
        (.Not.Associated(sorbates(sorb)%aslow))) Then
      Nullify(accel)
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " neither acceleration array is initialized"
      Stop
      Return
    End If

    If (fast) Then
      accel => sorbates(sorb)%afast
    Else
      accel => sorbates(sorb)%aslow
    End If

  End Subroutine config_getaccel


  !----------------------------------------------------------------------------
  ! Zero the accelerations in the species array and insert new from the 
  ! storage.  Note that we are assuming here that the storage structure
  ! contains forces (= -gradients).
  ! Requires:  species -- species data structure 
  !            grads -- 2D array of gradient vectors (atom,molec)
  !            fast -- True => treat Fast interactions, False => treat slow
  !----------------------------------------------------------------------------
  Subroutine config_putaccels(species,spc,grads,fast)
    Type(AtMolCoords), Intent(InOut)          :: species
    Integer, Intent(In)                       :: spc
    Type(VecType), Dimension(:,:), Intent(In) :: grads
    Logical, Intent(In)                       :: fast

    Integer            :: a,m

    !** put gradients(forces?) into species for each atom
    If (fast) Then
      Do m = 1,Size(grads,2)
        Do a = 1,Size(grads,1)
          species%afast(a,m) = grads(a,m)
        End Do
      End Do
    Else
      Do m = 1,Size(grads,2)
        Do a = 1,Size(grads,1)
          species%aslow(a,m) = grads(a,m)
        End Do
      End Do
    End If

    !** convert to accelerations
    Call config_conv2accel(species,spc,fast)

  End Subroutine config_putaccels

  !------------------------------------------------------------------------
  ! This routine takes the coordinate structures and generates
  ! the simulation cell coordinates (%r) and the continuation vectors 
  ! (%cr) from the principal coordinates.  Operates on a single species
  !------------------------------------------------------------------------
  Subroutine config_genfromparent(simcell,species)
    Type(SimCell_Params), Intent(In)   :: simcell
    Type(AtMolCoords), Intent(InOut)   :: species

    Integer                          :: a,m

    Do m = 1,Size(species%coords,2)
      Do a = 1,Size(species%coords,1)
        species%coords(a,m)%r = simcell_maptosimcell(simcell, &
            species%coords(a,m)%rp,species%coords(a,m)%cr)
      End Do
    End Do

  End Subroutine config_genfromparent

  !----------------------------------------------------------------------------
  ! This function returns the TOTAL number of vectors of all sorbates, i.e.,
  ! sorbates(1)%nmoles*sorbates(1)%natoms*sorbates(2)%nmoles*...
  !----------------------------------------------------------------------------
  Integer Function config_totalvecs(sorbates)
    Type(AtMolCoords), Dimension(:), Intent(In) :: sorbates
    Integer :: i,j

    j = 0

    Do i = 1, Size(sorbates,1)
      j = j + (sorbates(i)%nmoles)*(sorbates(i)%natoms)
    End Do
    config_totalvecs = j
  End Function config_totalvecs

  !----------------------------------------------------------------------------
  ! Read the fields of sorbate "sorbno" from the restart file
  !----------------------------------------------------------------------------
  Subroutine config_readrestartfile(sorbates, simcell, sorbno, filename)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(SimCell_Params), Intent(in)               :: simcell
    Integer, Intent(in)                            :: sorbno
    Character(len=strLen), Intent(in)              :: filename

    Integer                  :: i, unitno, natoms, a, m, temp_index, lo, hi
    Integer                  :: lineno, nfields, restart_nmoles, restart_natoms
    Logical                  :: found, need_update
    Character(len=strLen)    :: molecname, srchstr, comment_line
    Character(len=lstrLen)   :: line, portion
    Character(len=strLen), Dimension(20) :: fields
    
    !** Open the file if not already open
    unitno = file_open(filename,110)
    Rewind(unitno)

    !** Get some sundry information
    molecname = molecules_name(sorbno)
    natoms    = molecules_getnatoms(sorbno)

    !** Find the string "_Molecule_Name_:"
    srchstr = "_MOLECULE_NAME_:"
    found = .False.

    !** temp_index is used to avoid an infinte loop here 
    !** Assumes that only MAX_SORBS molecules are there in the file 
    Do temp_index=1,MAX_SORBS
      lineno = filesrchstr(unitno, srchstr, line, .false.)
      nfields = split(line, fields, " ")
      If (fields(2) == molecname) Then
        found = .True.
        Exit
      End If
    End Do

    If (.Not. (found)) Then
      Write(0,'(1x,2a,i4, 4a)') __FILE__," : ",__LINE__, &
          " Could not find the molecule ", Trim(molecname), &
          " in restart file ", filename
      Stop
    End If

    !** Read the natoms and nmoles
    Read(unitno,*) restart_nmoles
    Read(unitno,*) restart_natoms
    If (restart_natoms /= natoms) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " The number of atoms do not match "
      Stop
    End If
    sorbates(sorbno)%natoms = natoms
    sorbates(sorbno)%nmoles = restart_nmoles
    Call config_allocfields(sorbates(sorbno), sorbno, restart_nmoles)
    
    !** Read the principal coordinates and fill remainder of position fields
    Read(unitno,*) comment_line      !** Comment line
    nfields =  split(comment_line, fields)
    If (toupper(fields(1)) /= "_PRINCIPAL") Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Not reading from a RESTART file"
      Stop
    End If
    Do m = 1,restart_nmoles
      a = 0
      Do
        !** Read the line, strip comments, then break into vectors
        Read(unitno,'(a)') line
        line = stripcmnt(line)
        nfields = split(line,fields)
        If (Mod(nfields,3) /= 0) Then
          Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
              " Number of entries on restart file line not divisible by three"
          Write(0,'(2x,2a)') "Line is: ",Trim(line)
          Stop          
        End If
        Do i = 1,Int(nfields/3)
          lo = 1 + 3*(i-1)
          hi = lo + 2
          a = a + 1
          portion = combine(fields(lo:hi))
          Read(portion,*) sorbates(sorbno)%coords(a,m)%rp
        End Do
        
        If (a == natoms) Exit
      End Do

    End Do

    !** Read the velocities
    Read(unitno,*)       !** Comment line
    Do m = 1,restart_nmoles
      a = 0
      Do
        !** Read the line, strip comments, then break into vectors
        Read(unitno,'(a)') line
        line = stripcmnt(line)
        nfields = split(line,fields)
        If (Mod(nfields,3) /= 0) Then
          Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
              " Number of entries on restart file line not divisible by three"
          Write(0,'(2x,2a)') "Line is: ",Trim(line)
          Stop          
        End If
        Do i = 1,Int(nfields/3)
          lo = 1 + 3*(i-1)
          hi = lo + 2
          a = a + 1
          portion = combine(fields(lo:hi))
          Read(portion,*) sorbates(sorbno)%coords(a,m)%v
        End Do
        
        If (a == natoms) Exit
      End Do
    End Do

    !** Read the generalized coordinates
    Read(unitno,*)       !** Comment line
    Do m = 1,restart_nmoles
      Call gcmodels_readrestartfile(sorbates(sorbno)%gcoords(m), sorbno, &
          unitno,need_update)
      If (need_update) Then
        Call gcmodels_toxyz(sorbates(sorbno)%gcoords(m),&
            sorbates(sorbno)%coords(1:natoms,m)%rp)
      Endif
    End Do

    Do m = 1,restart_nmoles
      Call simcell_pbc(simcell, sorbates(sorbno)%coords(1:natoms,m)%rp, &
          sorbates(sorbno)%coords(1:natoms, m)%r, &
          sorbates(sorbno)%coords(1:natoms, m)%cr)
    End Do

  End Subroutine config_readrestartfile

  !-------------------------------------------------------------------------
  ! Updates coordinates from an xyz file
  !-------------------------------------------------------------------------
  Subroutine config_xyz2config(species,scell,filename,success)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(Simcell_Params) :: scell
    Character(len=strLen), Intent(in)              :: filename
    Logical, Intent(out) :: success    
    Integer       :: nsorbs, spc
    nsorbs = Size(species)
    Do spc = 1,nsorbs
      If (config_isfixed(species(spc))) Cycle
        Call config_xyz2spc(species, scell, spc, filename,success)
        If (.Not.success) Return
    End Do
  End Subroutine config_xyz2config




  !-------------------------------------------------------------------------
  ! Updates coordinates from an xyz file
  !-------------------------------------------------------------------------
  Subroutine config_xyz2spc(species,scell,spc,filename,success)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(Simcell_Params) :: scell
    Integer, Intent(in)                            :: spc
    Character(len=strLen), Intent(in)              :: filename
    Logical, Intent(out) :: success    

    Integer       :: nmoles, mol
    nmoles = config_getnmoles(species, spc)
    Do mol = 1,nmoles
      Call config_xyz2molec(species, scell, spc, mol, filename,success)
      If (.Not.success) Return
    End Do
  End Subroutine config_xyz2spc


  !----------------------------------------------------------------------------
  ! Read the fields of species "sorbno" from the xyzfile
  !----------------------------------------------------------------------------
  Subroutine config_xyz2molec(species, simcell, spc, mol, filename, success)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: species
    Type(SimCell_Params), Intent(in)               :: simcell
    Integer, Intent(in)                            :: spc, mol
    Character(len=strLen), Intent(in)              :: filename
    Logical, Intent(out) :: success    
    Integer                  :: i, unitno, natoms, a, err
    Integer                  :: lineno, nfields, xyz_natoms
    Character(len=strLen)    :: molecname, srchstr, comment_line
    Character(len=lstrLen)   :: line, portion
    Character(len=strLen), Dimension(20) :: fields
    
    !** Open the file if not already open
    unitno = file_open(filename,110)
    success=.True.
    
    !** Get some sundry information
    molecname = molecules_name(spc)
    natoms    = molecules_getnatoms(spc)

    !** temp_index is used to avoid an infinte loop here 
    !** Assumes that only MAX_SORBS molecules are there in the file 

    Read(unitno,*,IOSTAT=err) xyz_natoms
    If (err/=0) Then
      success=.false.
      Return
    Endif
    If (xyz_natoms /= natoms) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " The number of atoms do not match in xyz file "
      Stop
    End If
    species(spc)%natoms = natoms

    !** Read the principal coordinates and fill remainder of position fields
    Read(unitno,*) comment_line      !** Comment line

    Do a=1,natoms 
      !** Read the line, strip comments, then break into vectors
        Read(unitno,'(a)') line
        line = stripcmnt(line)
        nfields = split(line,fields)
        If (nfields<4) Then
          Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
              " Probelm reading xyz line "
          Write(0,'(2x,2a)') "Line is: ",Trim(line)
          Stop          
        End If
        portion = combine(fields(2:4))
        Read(portion,*) species(spc)%coords(a,mol)%rp
      End Do
      ! we are reading xyz so set vels=0
      species(spc)%coords(a,mol)%v=zero

      Call simcell_pbc(simcell, species(spc)%coords(1:natoms,mol)%rp, &
          species(spc)%coords(1:natoms, mol)%r, &
          species(spc)%coords(1:natoms, mol)%cr)

    End Subroutine config_xyz2molec


  !---------------------------------------------------------------
  ! Write the  coordinates and energies to the restart/crash file
  ! "filename" for all sorbates.
  !---------------------------------------------------------------
  Subroutine config_writerestartfile(sorbates, optfilename)
    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Character(*), Intent(in), Optional        :: optfilename 

    Integer                 :: i, unitno, rsunit, error, natoms, nmoles
    Integer                 :: a, m, nsorbs
    Character(len=xlstrLen) :: string
    Character(len=strLen)   :: molecname, filename

    !** Get the filename if it is not passed a parameter
    If (Present(optfilename)) Then
      filename = optfilename
      rsunit = file_getunit(filename)
    Else
      ! Get the filename from the file tag "d_res_file"
      Call file_gettype(d_res_file, filename, rsunit)
      If (rsunit == 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " restart file tag not set"
        Stop
      End If
    End If
    
    !** Open the file if not already open
    unitno = isfileopen(filename)
    If (unitno < 0) Then
      unitno = rsunit
      Open(unit = unitno, file=filename, status='unknown', IOSTAT=error)
      If (error /= 0) Then
        Close(unitno)
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            "  Could not open file ", trim(filename)
        Stop
      End If
    End If
    
    Write(unitno,*) Trim(general_getfiledesc())
    Write(unitno,*) general_getcurrentiter()
    Write(unitno,*) general_getcurrenttime()

    nsorbs = molecules_getnsorbs()
    Do i=1, nsorbs
      molecname = molecules_name(i)
      nmoles = sorbates(i)%nmoles
      natoms = sorbates(i)%natoms
      Write(unitno,'(2a)') "_MOLECULE_NAME_: ", molecname
      Write(unitno,'(i20, a)') nmoles, " # Nmoles"
      Write(unitno,'(i20, a)') natoms, " # Natoms"
      Write(unitno,*) "_Principal Coordinates Only_"
      Do m = 1, nmoles
        Do a = 1, natoms
          string = vector_display(sorbates(i)%coords(a,m)%rp, "f13.5")
          Write(unitno,'(a)') Trim(string)
        End do
      End Do

      Write(unitno,*) "_Velocities_"
      Do m = 1, nmoles
        Do a = 1,natoms
          string = vector_display(sorbates(i)%coords(a,m)%v, "f13.5")
          Write(unitno,'(a)') Trim(string)
          End do
      End Do

      Write(unitno,*) "_Generalized Coordinates_"
      Do m = 1,nmoles
        Call gcmodels_dumprestartinfo(sorbates(i)%gcoords(m), i, unitno)
      End Do
    End Do

    !** All done here
    Close(unitno)
    Return
  End Subroutine config_writerestartfile



  !----------------------------------------------------------------------------
  ! Returns the total kinetic energy of the subsystem 
  ! Returns energy in KJ
  ! Requires:  species -- species data structure
  !            spc -- species number
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function config_subsysKE(species,subsys)
    Type(AtMolCoords), Intent(In), Dimension(:) :: species 
    Integer, Dimension(3)                         :: subsys

    Integer         :: spc, nsorbs
    Real(kind=RDbl) :: sumKE

    !** Zero the kinetic energy.
    sumKE = zero 

    nsorbs=Size(species,1)
    If (subsys(1)/=0) Then
      Write(*,*) "wrong subset", subsys
      ! ** right now we can only return full system KE
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Else
      Do spc=1,nsorbs
        If (config_isfixed(species(spc))) Cycle
        sumKe=sumKE+config_kineticEnergy(species,spc)
      Enddo
    Endif

    config_subsysKE=sumKE

  End Function config_subsysKE


  !----------------------------------------------------------------------------
  ! Calculate the total kinetic energy of molecules "sorb"
  ! Returns in KJ
  ! Requires:  species -- species data structure
  !            spc -- species number
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function config_kineticEnergy(sorbates,mType)
    Type(AtMolCoords), Intent(In), Dimension(:) :: sorbates
    Integer, Intent(In)                         :: mType

    Integer         :: nmoles, natoms, atype, a, m, partner
    Logical         :: isCore
    Real(kind=RDbl) :: mass, pmass
    Type(vecType)   :: v

    !** Zero the kinetic energy.
    config_kineticEnergy = 0.0_RDbl

    !** Check if the molecule is fixed. If so, there is no kinetic energy.
    If (config_isfixed(sorbates(mType))) Return

    !** Better check to ensure these coordinates exist!
    If (.Not.Associated(sorbates(mType)%coords)) Return

    !** Get the number of molecules and atoms
    nmoles = config_getnmoles(sorbates,mType)
    natoms = molecules_getnatoms(mType)

    Do m = 1, nmoles
      Do a = 1, natoms
        atype = molecules_getnthatom(mType,a)
        mass = atom_getmass(atype)
        !** Check to see if this is a core-shell pair. If so, we
        !** must add the COM velocity, not the core and shell 
        !** velocities individually.
        If (molecules_isCoreShell(mType,a,isCore,partner)) Then
          If (.Not.isCore) Cycle
          !** We need the COM velocity
          pMass = atom_getmass(molecules_getnthatom(mType,partner))
          v = (sorbates(mType)%coords(a,m)%v*mass +  &
              sorbates(mType)%coords(partner,m)%v*pMass) / (mass + pMass)
          mass = mass + pMass
        Else
          v = sorbates(mType)%coords(a,m)%v
        End If
        config_kineticEnergy = config_kineticEnergy + v*v*mass
      End Do
    End Do

    !** What is the unit here-> Returns energy in KJ
    config_kineticEnergy = 0.5_RDbl*config_kineticEnergy*scaleke

  End Function config_kineticEnergy

  !----------------------------------------------------------------------------
  ! Calculate the total kinetic energy of molecules "sorb"
  ! Requires:  species -- species data structure
  !            spc -- species number
  !            tType -- string identifying temperature type
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function config_kineticTemp(species,spc,tType)
    Type(AtMolCoords), Intent(In), Dimension(:) :: species
    Integer, Intent(In)                         :: spc
    Character(*), Intent(In)                    :: tType

    Integer               :: a,m,natoms,nmoles,atomt,partner,dof
    Logical               :: isCore
    Real(kind=RDbl)       :: ke,mmass,mmassi2,mass,pMass
    Type(VecType)         :: vels,v,vatom,vpartner

    natoms = molecules_getnatoms(spc)
    nmoles = config_getnmoles(species,spc)
    If (nmoles==0) Then
      config_kineticTemp = zero
      Return
    Endif

    ke = 0.0_RDbl

    !** Get the temperature based on the passed string
    Select Case(ToUpper(tType))
    Case ('ATOMIC')
      Do m = 1,nmoles
        Do a = 1,natoms
          !** Get the atom type and mass
          atomt = molecules_getnthatom(spc,a)
          mass = atom_getmass(atomt)

          !** Check to see if this is a core-shell pair. If so, we
          !** must add the COM velocity, not the core and shell 
          !** velocities individually.
          If (molecules_isCoreShell(spc,a,isCore,partner)) Then
            If (.Not.isCore) Cycle

            !** We need the COM velocity
            pMass = atom_getmass(molecules_getnthatom(spc,partner))
            vatom = species(spc)%coords(a,m)%v
            vpartner = species(spc)%coords(partner,m)%v
            v = (vatom*mass + vpartner*pMass) / (mass + pMass)
            mass = mass + pMass
          Else
            v = species(spc)%coords(a,m)%v
          End If

          !** Sum the kinetic energy
          ke = ke + v*v*mass 
        End Do
      End Do

      dof = molecules_getdof(spc)
      config_kineticTemp = ke*scaleke / (kjmole_kb*Real(dof*nmoles))

    Case ('MOLECULAR')
      mmass = molecules_getmass(spc)
      mmassi2 = 1.0_RDbl / mmass / mmass

      Do m = 1, nmoles
        vels = 0.0_RDbl
        Do a = 1, natoms

          !** Get the atom type and mass
          atomt = molecules_getnthatom(spc,a)
          mass = atom_getmass(atomt)

          !** Check for a core-shell pair
          If (molecules_isCoreShell(spc,a,isCore,partner)) Then
            If (.Not.isCore) Cycle

            !** We need the COM velocity
            pMass = atom_getmass(molecules_getnthatom(spc,partner))
            vatom = species(spc)%coords(a,m)%v
            vpartner = species(spc)%coords(partner,m)%v
            v = (vatom*mass + vpartner*pMass) / (mass + pMass)
            mass = mass + pMass
          Else
            v = species(spc)%coords(a,m)%v
          End If

          vels = v*mass + vels

        End Do
        ke = vels*vels*mmassi2 + ke
      End Do
      config_kineticTemp = ke*mmass*scaleke / (kjmole_kb*Real(3*nmoles))

    End Select

  End Function config_kineticTemp

  !-------------------------------------------------------
  ! Read equilibrium file generated by the old MD code
  !-------------------------------------------------------
  Subroutine config_readcrashfile(sorbates,simcell,sorbname,filename)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: sorbates
    Type(SimCell_Params), Intent(In) :: simcell
    Character(*), Intent(In) :: sorbname
    Character(*), Intent(In) :: filename

    Integer :: nsorbs, i, j, unitno, maxatoms, totalmoles, lastindex, &
        error, natoms
    Character(len=25) :: molecname
    Integer, Dimension(:), Allocatable :: nmoles, index1, index2
    Real(kind=RDbl), Dimension(:,:), Allocatable :: rx, ry, rz, rxp, ryp, &
        rzp, vx, vy, vz
    Integer, Dimension(:,:), Allocatable ::  cxp, cyp, czp

    !** Check if the file is open otherwise open it
    unitno = isfileopen(filename)
    If (unitno < 0) Then
      unitno = file_getunit(filename)
      Open(file=filename, status = 'old', unit=unitno, form='unformatted')
    End If

    !** First, we need to know how many sorbate types we should have
    nsorbs = molecules_getnsorbs()

    !** I'm not sure how this will fare in the new zeolite-less code
    !** Without the zeolite, we'll have to read in each of the 
    !** sorbates and fill in the info as we go along. Ugh.
!!$    !** Check to see if any are zeolites, which the old code does not track
!!$    j = nsorbs
!!$    Do i = 1,j
!!$      If (zeolite_iszeolite(i)) nsorbs = nsorbs - 1
!!$    End Do

    Allocate(nmoles(nsorbs),index2(nsorbs),index1(nsorbs),STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, & 
          " Could not allocate 'nmoles' and 'molecname'"
      Stop
    End If

    !** Now read in the sorbate names and the number of molecules for each
    ! read junk first
    Read(unitno) i,j
    maxatoms = 0
    totalmoles = 0
    lastindex = 0
    Do i = 1,nsorbs
      Read(unitno) molecname
      j = molecules_gettype(molecname)
      If (j == 0) Then
        Write(0,'(1x,2a,i4,5a)') __FILE__," : ",__LINE__, & 
            " Could not find ",Trim(molecname),&
            " listed specified in crashfile ", Trim(filename), &
            " in molecules structure."
        Stop
      End If
      Read(unitno) nmoles(j)
      index1(j) = lastindex + 1              
      index2(j) = index1(j) + nmoles(j) - 1
      lastindex = index2(j)
      totalmoles = totalmoles + nmoles(j)
      If (maxatoms < molecules_getnatoms(i)) maxatoms = molecules_getnatoms(i)

      !** Allocate the memory for all the fields of sorbates(sorbtype)
      Call config_allocfields(sorbates(j), j, nmoles(j))
      sorbates(j)%nmoles = nmoles(j)
    End Do

    !** Now we allocate the old data structure type. 
    Allocate(rx(maxatoms,totalmoles),ry(maxatoms,totalmoles), &
        rz(maxatoms,totalmoles),rxp(maxatoms,totalmoles), &
        ryp(maxatoms,totalmoles), rzp(maxatoms,totalmoles), &
        vx(maxatoms,totalmoles),vy(maxatoms,totalmoles), &
        vz(maxatoms,totalmoles),cxp(maxatoms,totalmoles), &
        cyp(maxatoms,totalmoles),czp(maxatoms,totalmoles), stat=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, & 
          " Could not allocate memory to read crash file."
      Stop
    End If

    !** Read the file
    Read(unitno) rx, ry, rz, rxp, ryp, rzp, cxp, cyp, czp, vx, vy, vz

    !** Copy the data into the sorbates structure
    j = molecules_gettype(sorbname)
    natoms = molecules_getnatoms(j)
    sorbates(j)%coords(1:natoms,1:nmoles(j))%r%comp(1) =  &
        Real(rx(1:natoms,index1(j):index2(j)),kind=RDbl)
    sorbates(j)%coords(1:natoms,1:nmoles(j))%r%comp(2) =  &
        Real(ry(1:natoms,index1(j):index2(j)),kind=RDbl)
    sorbates(j)%coords(1:natoms,1:nmoles(j))%r%comp(3) =  &
        Real(rz(1:natoms,index1(j):index2(j)),kind=RDbl)
    sorbates(j)%coords(1:natoms,1:nmoles(j))%rp%comp(1) = &
        rxp(1:natoms,index1(j):index2(j))
    sorbates(j)%coords(1:natoms,1:nmoles(j))%rp%comp(2) = &
        ryp(1:natoms,index1(j):index2(j))
    sorbates(j)%coords(1:natoms,1:nmoles(j))%rp%comp(3) = &
        rzp(1:natoms,index1(j):index2(j))
    sorbates(j)%coords(1:natoms,1:nmoles(j))%cr%comp(1) = &
        cxp(1:natoms,index1(j):index2(j))
    sorbates(j)%coords(1:natoms,1:nmoles(j))%cr%comp(2) = &
        cyp(1:natoms,index1(j):index2(j))
    sorbates(j)%coords(1:natoms,1:nmoles(j))%cr%comp(3) = &
        czp(1:natoms,index1(j):index2(j))
    sorbates(j)%coords(1:natoms,1:nmoles(j))%v%comp(1) = &
        vx(1:natoms,index1(j):index2(j))
    sorbates(j)%coords(1:natoms,1:nmoles(j))%v%comp(2) = &
        vy(1:natoms,index1(j):index2(j))
    sorbates(j)%coords(1:natoms,1:nmoles(j))%v%comp(3) = &
        vz(1:natoms,index1(j):index2(j))

    !** Apply periodic boundry conditions
    Do i = 1, nmoles(j)
      Call simcell_pbc(simcell, &
          sorbates(j)%coords(1:natoms,i)%rp, &
          sorbates(j)%coords(1:natoms,i)%r, &
          sorbates(j)%coords(1:natoms,i)%cr)
    End Do

    Deallocate(rx,ry,rz,rxp,ryp,rzp,cxp,cyp,czp,vx,vy,vz,index1, &
        index2,nmoles,stat=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, & 
          " Could not deallocate memory from reading crash file."
      Stop
    End If

  End Subroutine config_readcrashfile

  !-------------------------------------------------------
  ! Read equilibrium file generated by one of the
  ! Monte Carlo routines.
  !-------------------------------------------------------
  Subroutine config_read_equilfile(sorbates, simcell, sorbname, filename, &
      sznrg)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(SimCell_Params), Intent(in) :: simcell
    Character(*), Intent(in) :: sorbname
    Character(*), Intent(in)  :: filename
    Real(kind=RDbl), Intent(Out), Optional :: sznrg

    ! **** Local Variables
    Character(len=2*strLen):: line
    Character(len=strLen)  :: molec_name
    Character(len=strLen), Dimension(10) :: fields

    Integer, Parameter  :: MAX_ATOMS = 10000
    Integer             :: j, k, nmoles, natoms
    Integer             :: unitno, sorbtype, lineno, nfields
    Integer             :: found
    Real(kind=Rdbl), Dimension(3, MAX_ATOMS)  :: atmcoords
    Real(kind=RDbl)        :: sorb_zeonrg
    Integer :: x,y,z
    
    ! Check if the file is open otherwise open it
    unitno = isfileopen(filename)
    If (unitno < 0) Then
      unitno = file_getunit(filename)
      Open(file=filename, status = 'old', unit=unitno)
    End If

    Rewind(unitno)

    ! Read the first few lines of stuff
    Read(unitno,*) line    ! This is the description line
    Read(unitno,*) x,y,z   ! This is the number of unit cells.
    ! Check the unit cell lengths
    If ((x /= simcell_getnuc(simcell,'x')).Or. &
        (y /= simcell_getnuc(simcell,'y')).Or. &
        (z /= simcell_getnuc(simcell,'z'))) Then
      Write(0,'(1x,2a)') "WARNING: Number of unit cells defined in the ", &
          "equilibrium"
      Write(0,'(1x,a)') "file does not match the number in the control file."
      Write(0,'(3x,a,3i3)') "Equil file : ",x,y,z
      Write(0,'(3x,a,3i3)') "Cntrl file : ",simcell_getnuc(simcell,'x'), &
          simcell_getnuc(simcell,'y'), simcell_getnuc(simcell,'z')
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

    ! Find the sorbate name in the file
    found = 0
    Do While (found == 0)
      lineno = filesrchstr(unitno, sorbname, line)
      If (lineno == 0) Then
        Write(0,'(1x,2a,i4, 4a)') __FILE__," : ",__LINE__, &
            " Could not find the sorbate ", sorbname, " in the file ", & 
            filename
        Stop
      End If
      line = stripcmnt(line)
      nfields = split(line, fields)
      If (fields(1) /= '') Then
        found = 1
      End If
    End Do

    ! Split the input line to get the sorbate name and the no. of moles
    nfields  = split(line, fields)
    molec_name = fields(1)
    nmoles   = toint(fields(2))
    If (trim(molec_name) /= trim(sorbname)) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Sorbate ", trim(sorbname), " not in the equil file"
      Stop
    End If

    ! Get the sorbate type for this molecule
    sorbtype = molecules_gettype(molec_name)
    natoms = molecules_getnatoms(sorbtype)
    sorbates(sorbtype)%natoms = natoms
    sorbates(sorbtype)%nmoles = nmoles

    ! Read in the sorbate zeolite energy
    Read(unitno,*) sorb_zeonrg
    !MDEBUG
    Write(*,*) "Sorbate-zeolite energy read in : ",sorb_zeonrg
    ! Assign sznrg if it is present
    If (Present(sznrg)) sznrg = sorb_zeonrg

    !** Allocate the memory for all the fields of sorbates(sorbtype)
    Call config_allocfields(sorbates(sorbtype), sorbtype, nmoles)

    !** Read in the coordinates from the equilibrium file
    Do j=1, nmoles
      Read(unitno,*) (atmcoords(1,k), atmcoords(2,k), &
          atmcoords(3,k), k=1, natoms)
      Do k=1, natoms
        sorbates(sorbtype)%coords(k,j)%rp = atmcoords(1:3, k)
      End Do
      Call simcell_pbc(simcell, sorbates(sorbtype)%coords(1:natoms,j)%rp, &
          sorbates(sorbtype)%coords(1:natoms, j)%r, &
          sorbates(sorbtype)%coords(1:natoms, j)%cr)
    Enddo

    Close(unit=unitno)
  End Subroutine config_read_equilfile
  
  !----------------------------------------------------------------------------
  ! Just reads in and returns the total energy, pairwise energy, and 
  ! temperature stored in the equilibration file for checking
  !----------------------------------------------------------------------------
  Subroutine config_read_equilfile_nrg(filename,ssnrg,totnrg,T)
    Character(*), Intent(In)     :: filename
    Real(kind=RDbl), Intent(Out) :: ssnrg, totnrg, T

    ! **** Local Variables
    Character(len=2*strLen):: line
    Integer                :: unitno
    Integer                :: x,y,z
    
    ! Check if the file is open otherwise open it
    unitno = isfileopen(filename)
    If (unitno < 0) Then
      unitno = file_getunit(filename)
      Open(file=filename, status = 'old', unit=unitno)
    End If

    Rewind(unitno)

    ! Read the first few lines of stuff
    Read(unitno,*) line    ! This is the description line
    Read(unitno,*) x,y,z   ! This is the number of unit cells.
    ! Read in the temperature for equilibration
    Read(unitno,*) T
    ! Read in the total energy
    Read(unitno,*) totnrg
    ! Read in the sorb-sorb energy
    Read(unitno,*) ssnrg

    ! That's all folks
    Close(unitno)
  End Subroutine config_read_equilfile_nrg

  !----------------------------------------------------------------------------
  ! Reads in a species configuration from a restart file made using Louie's 
  ! old GCMC code.
  ! Requires: species -- species data structure
  !           simcell -- simulation cell information
  !           spcname -- species name
  !           filename -- restart filename
  !----------------------------------------------------------------------------
  Subroutine config_readoldgcmc(species,simcell,spcname,filename)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Character(*), Intent(In)                       :: spcname
    Character(*), Intent(In)                       :: filename

    Integer                          :: spc,m,natoms,error
    Real(kind=RDbl)                  :: theta,phi,psi
    Type(OldGCMCSystem)              :: oldsystem
    Type(VecType)                    :: com
    Type(VecType), Dimension(:), Pointer   :: xyzcoords

    !** initialize the old molecule reference configurations
    Call readoldgcmc_init(oldsystem,'molecule_defns',filename)

    !** make sure the reference structures match
    spc = molecules_gettype(spcname)
    Call readoldgcmc_chkspc(spc)

    !** Allocate the memory for all the fields of species(spc)
    species(spc)%nmoles = readoldgcmc_nmoles(spc,oldsystem)
    Call config_allocfields(species(spc),spc,species(spc)%nmoles)

    natoms = molecules_getnatoms(spc)

    !** allocate memory for the temporary xyz coordinates
    Allocate(xyzcoords(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'xyzcoords')    

    !** set the new molecular states
    Do m = 1,species(spc)%nmoles
      Call gcmodels_initcoords(species(spc)%gcoords(m),spc)
      Call readoldgcmc_molecstate(oldsystem,m,spc,theta,phi,psi,com,xyzcoords)

#if DEBUG
      Do a = 1,natoms
        Write(*,*) a,Trim(vector_display(xyzcoords(a),'f8.3'))
      End Do
#endif

      !** unfortunately, can't use the above given eulerian angles because
      !** they are from another definition, regenerate angles from xyz here
      Call gcmodels_fromxyz(species(spc)%gcoords(m),xyzcoords,spc)

#if DEBUG
      !** check
      Call gcmodels_toxyz(species(spc)%gcoords(m),xyzcoords)
      Do a = 1,natoms
        Write(*,*) a,Trim(vector_display(xyzcoords(a),'f8.3'))
      End Do
      stop
#endif

      species(spc)%coords(1:natoms,m)%rp = xyzcoords

      Call simcell_pbc(simcell, species(spc)%coords(1:natoms,m)%rp, &
          species(spc)%coords(1:natoms,m)%r, &
          species(spc)%coords(1:natoms,m)%cr)
    End Do

    !** deallocate memory for the temporary xyz coordinates
    Deallocate(xyzcoords, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'xyzcoords')    

  End Subroutine config_readoldgcmc

  !---------------------------------------------------------
  ! Displays the coordinates of sorbate type "sorbno"
  !---------------------------------------------------------
  Subroutine config_displaycoords(sorbates, sorbno, optunit)
    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Integer, Intent(in)    :: sorbno
    Integer, Intent(in), Optional :: optunit

    Integer     :: nmoles, natoms, a, m, unitno 
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    Endif

    nmoles = config_getnmoles(sorbates, sorbno)
    natoms = config_getnatoms(sorbates, sorbno)
    Do m=1, nmoles
      Write(unitno, '(2x, a, i6)') "Molecule: ", m
      Do a=1, natoms
        Write(unitno,'(3x,i4, 2a24, a12)') a,  &
            Trim(vector_display(sorbates(sorbno)%coords(a,m)%rp, "f8.3")), &
            Trim(vector_display(sorbates(sorbno)%coords(a,m)%r, "f8.3")), &
            Trim(vector_display(sorbates(sorbno)%coords(a,m)%cr, "i4"))
      End Do
    End Do
  End Subroutine config_displaycoords

  !---------------------------------------------------------
  ! Displays the velocities of sorbate type "sorbno"
  !---------------------------------------------------------
  Subroutine config_displayvelocity(sorbates, sorbno, optunit)
    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Integer, Intent(in)    :: sorbno
    Integer, Intent(in), Optional :: optunit

    Integer     :: nmoles, natoms, a, m, unitno 
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    Endif

    nmoles = config_getnmoles(sorbates, sorbno)
    natoms = config_getnatoms(sorbates, sorbno)
    Do m=1, nmoles
      Write(unitno, '(2x, a, i6)') "Molecule: ", m
      Do a=1, natoms
        Write(unitno,'(3x,i4,a)') a,  &
            Trim(vector_display(sorbates(sorbno)%coords(a,m)%v, "e8.2"))
      End Do
    End Do
  End Subroutine config_displayvelocity

  !----------------------------------------------------------------------
  ! HACK for now
  ! Requires:  sorbates -- species data structures
  !            spc -- desired species number
  !            indent -- no. of spaces from the left margin
  !            optunit -- optional output unit
  !----------------------------------------------------------------------
  Subroutine config_dump(species,spc,indent,optunit)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: spc,indent
    Integer, Intent(In), Optional               :: optunit

    Integer                    :: m,natoms,nmoles,unitno
    Character(len=indent)      :: blank

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    nmoles = config_getnmoles(species,spc)
    natoms = config_getnatoms(species,spc)

    Do m = 1,nmoles
      Write(unitno, '(2a,i6)') blank,"Molecule: ", m
      Call config_dumpmol(species,spc,m,indent,unitno)
    End Do

  End Subroutine config_dump

  !----------------------------------------------------------------------
  ! Dumps information for one molecule to a specified unit or the screen
  ! Requires:  species -- species data structures
  !            spc -- desired species number
  !            mol -- desired molecule number
  !            indent -- no. of spaces from the left margin
  !            optunit -- optional output unit
  !----------------------------------------------------------------------
  Subroutine config_dumpmol(species,spc,mol,indent,optunit)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In)                         :: spc,mol,indent
    Integer, Intent(In), Optional               :: optunit

    Integer                    :: a,natoms,unitno
    Character(len=indent)      :: blank
    Character(len=lstrLen)     :: string

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    natoms = config_getnatoms(species,spc)

    Write(unitno,'(a,2i3)') 'Position vectors for spc,mol: ',spc,mol
    Do a = 1,natoms
      string = vector_display(species(spc)%coords(a,mol)%r,'f12.6')
!      string = vector_display(species(spc)%coords(a,mol)%rp,'f12.6')
      Write(unitno,'(a,i4,3x,a)') blank,a,Trim(string)
    End Do

    Write(unitno,'(a,2i3)') 'Acceleration vectors for spc,mol: ',spc,mol
    Do a = 1,natoms
!      string = vector_display(species(spc)%afast(a,mol),'e18.8')
      string = vector_display(species(spc)%afast(a,mol),'f14.6')
      Write(unitno,'(a,i4,3x,a)') blank,a,Trim(string)
    End Do

  End Subroutine config_dumpmol
  
  !----------------------------------------------------------------------
  ! The routine writes a brief overview of the configuration structure
  ! to the specified unit.
  ! Requires:  species -- species data structures
  !            optunit -- optional output unit
  !----------------------------------------------------------------------
  Subroutine config_display(species, optunit)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer, Intent(In), Optional               :: optunit

    Integer                 :: spcno, nspcs
    Integer                 :: nmoles, natoms, unitno 
    Character(len=strLen)   :: fmt
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    Endif

    nspcs = molecules_getnsorbs()
    fmt = "f8.3"

    Write(unitno, '(a)') dashedline
    Write(unitno, '(a)') "Configuration Section "
    Write(unitno, '(a, i6)') "No. of species : ", nspcs
    Write(unitno, '(a)') dashedline2
    Do spcno = 1, nspcs
      nmoles = species(spcno)%nmoles
      natoms = species(spcno)%natoms
      Write(unitno, '(2x, 2a)')  "Species name : ", molecules_name(spcno)
      Write(unitno, '(4x, a, i6)')  "No. of molecules : ", nmoles
      Write(unitno, '(4x, a, i4)')  "No. of atoms     : ", natoms
      Write(unitno, '(4x, 2a)')     "Initial Config. Source : ", &
          Trim(species(spcno)%sourcetype)
      Write(unitno, '(4x, 2a)')      "Source Filename: ", &
          Trim(species(spcno)%filename)
      Write(unitno, '(2x, a)') dashedline2
    End Do
    
  End Subroutine config_display

  !----------------------------------------------------------------------------
  ! Dumps the configuration in species to an .xyz file
  ! Requires:  species -- species data structure
  !            filename -- output file name
  !            pbcs -- flag indicating use of %rp (w/o PBCs) or %r (w/PBCs)
  !            comment -- optional comment for .xyz file
  !----------------------------------------------------------------------------
  Subroutine config_config2xyz(species,filename,pbcs,comment)
    Type(AtMolCoords), Dimension(:), Intent(In)    :: species
    Character(*), Intent(In)                       :: filename
    Logical, Intent(In)                            :: pbcs
    Character(*), Optional                         :: comment

    Integer                                 :: n,error,size,spc,high
    Integer                                 :: natoms,nmoles
    Type(XYZ_Entry), Dimension(:), Pointer  :: entries

    !** set the comment
    If (.NOT.(Present(comment))) Then
      comment = 'Created by MUSIC code'
    End If

    !** get the number of entries necessary
    size = 0
    Do spc = 1,molecules_getnsorbs()
      nmoles = config_getnmoles(species,spc)
      natoms = config_getnatoms(species,spc)
      size = size + nmoles*natoms
    End Do

    !** Size the entries array
    Allocate(entries(size),STAT=error)
    If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Initialize the entries array
    entries = XYZ_Entry('  ',VecType(0.0_RDbl))

    !** Get the species entries
    n = config_setxyz(species,(/0,0,0/),entries,pbcs)

    !** Dump the file
    Call visxyz_dump(entries,filename,'f8.3',comment)

    !** Deallocate the entries array
    Deallocate(entries,STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__)
    
  End Subroutine config_config2xyz

  !----------------------------------------------------------------------------
  ! Converts the single-species configuration into a Structure data type
  ! Requires:  struc -- resulting structure data type
  !            species -- species data structure
  !            spc -- species number
  !----------------------------------------------------------------------------
  Subroutine config_config2struc(struc,species,spc)
    Type(Structure), Intent(Out)     :: struc
    Type(AtMolCoords), Intent(In)    :: species
    Integer, Intent(In)              :: spc

    Integer                    :: a,atm,mol
    Type(VecType), Dimension(species%natoms*species%nmoles)    :: atomcoords
    Character(len=2), Dimension(species%natoms*species%nmoles) :: elements

    !** Copy the coordinates and element symbols
    a = 0
    Do atm = 1,species%natoms
      Do mol = 1,species%nmoles
        a = a + 1
        atomcoords(a) = species%coords(atm,mol)%rp
        elements(a) = molecules_atomsymbol(spc,atm)
      End Do
    End Do

    !** Set the structure
    Call readstruc_set(struc,atomcoords,elements)
    
  End Subroutine config_config2struc

  !----------------------------------------------------------------------------
  ! Dumps one or two subsets of the configuration in species to an .xyz file
  ! Requires:  species -- species data structure
  !            filename -- output file name
  !            comment -- comment for .xyz file
  !            pbcs -- flag indicating use of %rp (w/o PBCs) or %r (w/PBCs)
  !            subset1 -- first subset
  !            subset2 -- optional second subset
  !----------------------------------------------------------------------------
  Subroutine config_subset2xyz(species,filename,comment,pbcs,subset1,subset2)
    Type(AtMolCoords), Dimension(:), Intent(In)    :: species
    Character(*), Intent(In)                       :: filename,comment
    Logical, Intent(In)                            :: pbcs
    Integer, Dimension(3), Intent(In)              :: subset1
    Integer, Dimension(3), Intent(In), Optional    :: subset2

    Integer                                   :: n,error
    Integer                                   :: natoms,unitno
    Type(XYZ_Entry), Dimension(:), Pointer    :: entries

    !** Get the number of entries necessary
    natoms = config_nsubsetatoms(species,subset1)
    If (Present(subset2)) Then
      natoms = natoms + config_nsubsetatoms(species,subset2)
    End If

    !** Size the entries array
    Allocate(entries(natoms),STAT=error)
    If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__)

    !** initialize the entries array
    entries = XYZ_Entry('  ',VecType(0.0_RDbl))

    !** Get the individual entries for each subset
    n = config_setxyz(species,subset1,entries,pbcs)
    If (Present(subset2)) Then
      n = config_setxyz(species,subset2,entries(n+1:),pbcs)
    End If

    !** Dump the file
    Call visxyz_dump(entries,filename,'f8.3',comment)
    unitno = file_getunit(filename)
    Close(unit = unitno)

    !** Deallocate the entries array
    Deallocate(entries,STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__)
    
  End Subroutine config_subset2xyz

  !----------------------------------------------------------------------------
  ! This subroutine puts a subset of the species coordinates into the entries
  ! of an XYZ array.  Returns the number of entries created.  It dumps the
  ! simulation cell coordinates.
  ! Requires:  species -- species data structure
  !            subset -- system subset
  !            entries -- output xyz entries
  !            pbcs -- flag indicating use of %rp (w/o PBCs) or %r (w/PBCs)
  !----------------------------------------------------------------------------
  Integer Function config_setxyz(species,subset,entries,pbcs)
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Dimension(3), Intent(In)            :: subset
    Type(XYZ_Entry), Dimension(:), Intent(Out)   :: entries
    Logical, Intent(In)                          :: pbcs

    Integer          :: n,atype,arraysize,depth
    Integer          :: lospc,hispc,spc
    Integer          :: lomol,himol,mol
    Integer          :: loatm,hiatm,atm

    arraysize = Size(entries)
    depth = getdepth(subset)
    config_setxyz = 0
    n = 0

    !** Set the species number limits
    lospc = 1
    hispc = Size(species)
    If (depth > 0) Then
      lospc = subset(1)
      hispc = subset(1)
    End If

    !** Create the entries
    Do spc = lospc,hispc
      !** Set the molecule number limits
      lomol = 1
      himol = config_getnmoles(species,spc)
      If (depth > 1) Then
        lomol = subset(2)
        himol = subset(2)
      End If

      Do mol = lomol,himol
        !** Set the atom number limits
        loatm = 1
        hiatm = config_getnatoms(species,spc)
        If (depth > 2) Then
          loatm = subset(3)
          hiatm = subset(3)
        End If

        Do atm = loatm,hiatm
          n = n + 1
          If (n > arraysize) Then
            Write(0,'(2a,i4,a,2i5)') __FILE__,": ",__LINE__, &
                 '  Passed entries array is too small ',n,arraysize
            Stop                      
          End If 

          atype = config_getatype(species,spc,mol,atm)
          If (pbcs) Then
            entries(n) = visxyz_make(config_getr(species,(/spc,mol,atm/)), &
                atom_getsymbol(atype))
          Else
            entries(n) = visxyz_make(config_getrp(species,(/spc,mol,atm/)), &
                atom_getsymbol(atype))
          End If
        End Do

      End Do
    End Do

    config_setxyz = n

  End Function config_setxyz

  !----------------------------------------------------------------------------
  ! Dumps the positions, velocities, and forces for the given molecule
  ! of sorbate sorb to unitno
  !----------------------------------------------------------------------------
  Subroutine config_dumpmole(sorbate,molecule,caption,unitno)
    Type(AtMolCoords), Intent(In) :: sorbate
    Integer, Intent(In) :: molecule
    Integer, Intent(In), Optional :: unitno
    Integer :: unit, i, natoms
    Character(*), Intent(In), Optional :: caption

    If (.Not.Present(unitno)) Then
      unit = 0
    Else
      unit = unitno
    End If

    If (Present(caption)) Then
      Write(unit,'(a)') caption
    End If
    
    natoms = sorbate%natoms

    Do i = 1,natoms
      Write(unit,'(a,i2,a,i2,a,a)') "Princ. Pos for mol# ",molecule, &
          " atom# ",i," : ", &
          vector_display(sorbate%coords(i,molecule)%rp,'f15.4')
    End Do

    Do i = 1,natoms
      Write(unit,'(a,i2,a,i2,a,a)') "PBC Position for mol# ",molecule, &
          " atom# ",i," : ", &
          vector_display(sorbate%coords(i,molecule)%r,'f15.4')
    End Do

    Do i = 1,natoms
      Write(unit,'(a,i2,a,i2,a,a)') "Velocities for mol# ",molecule,&
          " atom# ",i," : ", &
          vector_display(sorbate%coords(i,molecule)%v,'f15.4')
    End Do
    
    If (Associated(sorbate%afast)) Then
      Do i = 1,natoms
        Write(unit,'(a,i2,a,i2,a,a)') " Fast Force on mol# ",molecule, &
            " atom# ",i," : ", &
            vector_display(sorbate%afast(i,molecule),'f15.4')
      End Do
    End If
    If (Associated(sorbate%aslow)) Then
      Do i = 1,natoms
        Write(unit,'(a,i2,a,i2,a,a)') " Slow Force on mol# ",molecule, &
            " atom# ",i," : ", &
            vector_display(sorbate%aslow(i,molecule),'f15.4')
      End Do
    End If
  End Subroutine config_dumpmole

  !-------------------------------------------------------------------------
  ! Makes sure that gccords and coords are consistent with each other. To 
  ! avoid any problems that may have happened due to wrong restartfile 
  ! Requires:  species -- full species data structure
  !-------------------------------------------------------------------------
  Subroutine config_makeconsistent(species)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species

    Integer       :: nsorbs, nmoles, spc, m, dof

    nsorbs = Size(species)
    Do spc = 1,nsorbs
      If (config_isfixed(species(spc))) Cycle
      dof = molecules_getdof(spc)
      If (dof < 5) Then
        Cycle
      Elseif (dof==5) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "This is the case of diatomic linear molecule,contact shaji"
        Write(*,*) "shaji@northwestern.edu"
        stop
      Endif

      nmoles = config_getnmoles(species, spc)

      Do m = 1,nmoles
        Call config_setfromrp(species, spc, m)
      End Do
    End Do

  End Subroutine config_makeconsistent

  !-------------------------------------------------------------------------
  ! Set the generalized coordinates for a molecule from the current 
  ! parent coordinates.
  ! Requires:  species -- full species data structure
  !            spc -- species number
  !            mol -- molecule number
  !-------------------------------------------------------------------------
  Subroutine config_setfromrp(species,spc,mol)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer, Intent(In)                            :: spc,mol

    Integer          :: natoms
    Type(VecType)    :: com

    natoms = config_getnatoms(species, spc)

    !** Calculate the COM from the principal coords
    com = config_getMolecCOM(species,spc,mol)

    Call gcmodels_setfrom_rp(species(spc)%gcoords(mol), spc, &
        species(spc)%coords(1:natoms,mol)%rp,com)

  End Subroutine config_setfromrp

  !-------------------------------------------------------------------------
  ! brings the %coords%rp fields close to origin. This helps in
  ! maintaining the significant digits while doing ideal gas
  ! simulations. Otherwise molecules tend to fly away
  ! Requires: species -- full species data structure
  !-------------------------------------------------------------------------
  Subroutine config_changerp(species,scell)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(Simcell_Params) :: scell

    Integer       :: nsorbs, nmoles, natoms, spc, m, a
    Type(VecType) :: com,newcom

    nsorbs = Size(species)
    Do spc = 1,nsorbs
      If (config_isfixed(species(spc))) Cycle
      nmoles = config_getnmoles(species, spc)
      natoms = config_getnatoms(species, spc)
      Do m = 1,nmoles

        !** from principal coords
        com = config_getMolecCOM(species, spc, m)

        !** do pbc on com
        Call simcell_pbc(scell, com, newcom )


        !** find positions relative to com and shift
        Do a=1,natoms
          species(spc)%coords(a,m)%rp=species(spc)%coords(a,m)%rp - com
          species(spc)%coords(a,m)%rp=species(spc)%coords(a,m)%rp + newcom
        End Do

        If (associated(species(spc)%gcoords(m)%rigidcoords)) Then
          If (species(spc)%gcoords(m)%rigidcoords%internally_flexible) Then
            Call gcmodels_setfrom_rp(species(spc)%gcoords(m), spc, &
                species(spc)%coords(1:natoms,m)%rp, newcom)
          Else
            species(spc)%gcoords(m)%rigidcoords%com=newcom
          Endif
        Endif

        Call simcell_pbc(scell, species(spc)%coords(1:natoms,m)%rp, &
            species(spc)%coords(1:natoms, m)%r, &
            species(spc)%coords(1:natoms, m)%cr)

      End Do
    End Do

  End Subroutine config_changerp

  !-------------------------------------------------------------------------
  ! Change dof to what is present in molecules
  ! Important for mc where gcoords are used
  !-------------------------------------------------------------------------
  Subroutine config_changedof(species )
    Type(AtMolCoords), Dimension(:),Intent(InOut) :: species
    Integer :: nsorbs, nmoles, spc, m, dof

    nsorbs = Size(species)
    Do spc = 1,nsorbs
      If (config_isfixed(species(spc))) Cycle
      dof = molecules_getdof(spc)

      !** for al fields
      nmoles = Size(species(spc)%gcoords)
      Do m = 1,nmoles
        Call gcmodels_changedof(species(spc)%gcoords(m), dof)
        If (molecules_isflexible(spc)) Then
          Call gcmodels_changeflex(species(spc)%gcoords(m), .True.)
        End If
      End Do
    End Do
    
  End Subroutine config_changedof
  
  !-------------------------------------------------------------------------
  ! get distance between two atoms , uses %rp
  !-------------------------------------------------------------------------
  Real(kind=RDbl) Function config_getdistance(species,spc1,m1,a1,spc2,m2,a2)
    Type(AtMolCoords), Dimension(:),Pointer :: species

    Integer            :: spc1, m1, a1, spc2, m2, a2

    config_getdistance = mag(species(spc1)%coords(a1,m1)%rp - &
        species(spc2)%coords(a2,m2)%rp )

  End Function config_getdistance
  
  !-------------------------------------------------------------------------
  ! get angle  between three atoms of one molecule in radian
  ! note : maybe we should generalise it to deal with any 3 atoms
  !-------------------------------------------------------------------------
  Real(kind=RDbl) Function config_getangle(species , spc, m, a1, a2, a3)
    Type(AtMolCoords), Dimension(:), Pointer :: species

    Integer          :: spc, m, a1, a2, a3

    config_getangle = vector_angle(species(spc)%coords(a1,m)%rp, &
        species(spc)%coords(a2,m)%rp,species(spc)%coords(a3,m)%rp)

    !** vector_angle gives back in degrees, convert back to radians 
    config_getangle=config_getangle*degToRad

  End Function config_getangle
  
  !-------------------------------------------------------------------------
  ! get angle  between three atoms of one molecule in radian
  ! look at cosexpansion.F90 for explanation
  !-------------------------------------------------------------------------
  Real(kind=RDbl) Function config_gettorangle(species , spc, m, a1, a2, a3, a4)
    Type(AtMolCoords), Dimension(:),Pointer :: species
    Integer :: spc, m, a1, a2, a3, a4

    Real(kind=RDbl)        :: a_a,a_b,a_c,b_b,b_c,c_c
    Real(kind=RDbl)        :: T1,T2,T3, T2T3_invroot
    Real(kind=RDbl)        :: cos, cos_mag
    Type(VecType)   :: a,b,c
    
    !**calculate a b c, these represent the vectors along the bonds
    a = species(spc)%coords(a2,m)%rp - species(spc)%coords(a1,m)%rp
    b = species(spc)%coords(a3,m)%rp - species(spc)%coords(a2,m)%rp 
    c = species(spc)%coords(a4,m)%rp - species(spc)%coords(a3,m)%rp  

    !**calculate the dot products 
    a_a=a*a
    a_b=a*b
    a_c=a*c
    b_b=b*b
    b_c=b*c
    c_c=c*c

    !**calculate T1,T2,T3, term that will have to used later
    T1= a_b*b_c - a_c*b_b
    T2= a_a*b_b - a_b*a_b
    T3= b_b*c_c - b_c*b_c

    T2T3_invroot = one/Sqrt(T2*T3)

    !** When theta is zero atom.1 and atom.4 are farthest from each other.
    !**cos(theta) is defined so that the above condn on theta is satisfied
    cos=-(T1*T2T3_invroot)
    
    cos_mag=Abs(cos)
    If ( (cos_mag>0.99999999_RDbl) .And. (cos_mag<1.00000001_RDBl) )Then
      cos = 0.999999999_RDbl * cos
    Endif
    
    config_gettorangle=Acos(cos)

  End Function config_gettorangle

  !--------------------------------------------------------------------
  ! Interpolates between two configuration end points to create new
  ! intermediate configurations.  If a species is marked as fixed, then
  ! this routine will also check to make sure that the two end points
  ! have the same coordinates for that species.
  ! Note that the image structures must already be sized
  ! Requires:  startpt -- start point species data structure
  !            endpt -- end point species data structure
  !            image -- intermediate species data structures to fill
  !            simcell -- simulation cell information
  !--------------------------------------------------------------------
  Subroutine config_makeimages(startpt,endpt,image,simcell)
    Type(AtMolCoords), Dimension(:), Intent(In)       :: startpt
    Type(AtMolCoords), Dimension(:), Intent(In)       :: endpt
    Type(AtMolCoords), Dimension(:,:), Intent(InOut)  :: image
    Type(SimCell_Params), Intent(In)                  :: simcell

    Integer               :: i,j,atm,mol,spc
    Integer               :: nspc,natms,nimages,nomatch
    Type(VecType)         :: disp
    Type(Structure)       :: struc1,struc2
    Character(len=strLen) :: string

    nspc = Size(startpt)
    nimages = Size(image,1)

    !** Check the end points for equilvalent FIXED species coordinates
    Do spc = 1,nspc
      If (startpt(spc)%fixed) Then
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) 'checking equivalence of FIXED species number ',spc
        Call config_config2struc(struc1,startpt(spc),spc)
        Call config_config2struc(struc2,endpt(spc),spc)
        nomatch = match_compare(struc1,struc2,1.0e-6_RDbl,.False.)
        If (nomatch /= 0) Then
          string = int2str(nomatch)
          Write(0,'(1x,2a,i4,3a,i2)') __FILE__,": ",__LINE__, &
              '  WARNING: Could not match at least atom ',Trim(string), &
              ' of species ',spc
          Write(0,'(2x,a)') 'between the two end-points'
!          nomatch = match_compare(struc1,struc2,1.0e-3_RDbl,.True.)
        End If
        Call readstruc_clean(struc1)
        Call readstruc_clean(struc2)
      End If
    End Do

    !** Size the images
    Do i = 1,nimages
      Call config_initcopy(image(i,:),startpt,.False.)
    End Do

    !** Perform the interpolation for each image
    Do i = 1,nimages
      Do spc = 1,nspc
        !** Simply copy the coordinates from the start point if it's FIXED
        If (startpt(spc)%fixed) Then
          image(i,spc)%coords = startpt(spc)%coords
          Cycle
        End If

        !** Otherwise, create each image along atomic vectors between endpts
        natms = startpt(spc)%natoms
        Do mol = 1,startpt(spc)%nmoles
          Do atm = 1,natms
            disp = endpt(spc)%coords(atm,mol)%rp - &
                startpt(spc)%coords(atm,mol)%rp
            disp = disp*(i * 1.0_Rdbl/(nimages + 1))
            image(i,spc)%coords(atm,mol)%rp = &
                startpt(spc)%coords(atm,mol)%rp + disp
          End Do

          !** Create the simulation cell coordinates using PBCs
          Call simcell_pbc(simcell, &
              image(i,spc)%coords(1:natms,mol)%rp, &
              image(i,spc)%coords(1:natms,mol)%r, &
              image(i,spc)%coords(1:natms,mol)%cr)
        End Do

      End Do
    End Do

  End Subroutine config_makeimages

  !----------------------------------------------------------------------
  ! Clean the configuration structure
  ! Requires:  species -- single species data structure
  !----------------------------------------------------------------------
  Subroutine config_clean(species)
    Type(AtMolCoords), Intent(InOut) :: species

    Integer       :: error,i

    Deallocate(species%afast, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    Deallocate(species%aslow, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    Deallocate(species%coords, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    Do i = 1,Size(species%gcoords)
      Call gcmodels_clean(species%gcoords(i))
    End Do

  End Subroutine config_clean
  
End Module config

