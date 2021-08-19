!-------------------------------------------------------------------------
! This module handles the interface of MUSIC with GULP.  All Input/Output
! (IO) calls go through this module.  It handles the creation of GULP
! input files, the calling of GULP and the translation of the output into
! something that MUSIC understands.  For now, it is setup only to do
! single-point calculations.  However, it is written to be easily 
! generalized to optimizations.  It is important to note that GULP handles
! forcefield parameters very differently and that it is not always 
! possible to have the same forcefield in GULP.
!
! The creation of a GULP control file is managed in two different ways,
! depending on the initialization information.  
! METHOD 1: The filename of a file containing the forcefield to be used
!           is passed into the init routine.  This file must contain a 
!           'species' section that specifies the atom symbol used in the
!           rest of the forcefield input, the core/shell type and the
!           charge.  Example:
!              species
!              si core  4.00000
!              o2 shel -2.81753
!           When the gulp atom types are initialized, this section will
!           be read and prepared for usage when dumping coordinates.  The
!           contents of the entire file will be dumped to the end of new
!           GULP control files.
!
! METHOD 2: MuSiC atom numbers and string-based parameters can alternatively
!           be passed into the initialization.  They will then be used to
!           build the GULP atom types and forcefield parameter input.  This
!           is a riskier procedure.  Be sure to check the resulting input
!           files.  May not be working.
!
! Needed Improvements: 
! 1) forcefield dumping should be improved, add 4-bodies and interpairs
!-------------------------------------------------------------------------

Module gulpio

  Use defaults, Only: strLen,lstrLen,RDbl,zeroTolerance,xlstrLen
  Use file, Only: file_open
  Use general, Only: genparams
  Use vector, Only: VecType,vector_display,Assignment(=),Operator(+), &
      Operator(-),Operator(*),Operator(/)
  Use utils, Only: toint,filesrchstr,split,toupper,toreal,allocerrdisplay, &
      int2str,firstchars,digits,deallocerrdisplay
  Use readgulp, Only: GULPsystem,GULPconfig,readgulp_getsystem, &
      readgulp_config2xyz,readgulp_display,readgulp_nrg,readgulp_forces
  Use simcell, Only: SimCell_Params,simcell_getell,simcell_angles,simcell_display
  Use atom, Only: atom_getsymbol,atom_getcharge,atom_getntypes,atom_getname, &
      atom_iontype
  Use molecules, Only: molecules_getnatoms,molecules_getcharge,molecules_name, &
      molecules_getnsorbs,molecules_gettype,molecules_getnatomtypes
  Use bsmodel, Only: STRETCH_KEY
  Use bbmodel, Only: BENDING_KEY
  Use tormodel, Only: TORSION_KEY
  Use ipmodel, Only: INTRAPAIR_KEY
#ifdef NAGCOMPILER
  Use f90_unix_proc, Only: system
#endif

  Implicit None

  Private
  Public :: GULP_Atomic_Params,GULP_Atomic_Info,GULP_Forcefield_Input, &
      gulpio_init,gulpio_snglpt,gulpio_display,gulpio_clean

  Character(len=strLen)   :: gulpio_outfile = 'gulp_output'

  !** GULP atoms must have a unique charge, they differ from atoms structure
  Type GULP_Atomic_Params
    Character(len=strLen)   :: atom_name    ! to match atoms structure
    Integer                 :: atom_index   ! to match atoms structure
    Character(len=2)        :: orig_element ! atomic symbol in atoms structure
    Character(len=2)        :: element      ! atomic symbol
    Character(len=4)        :: symbol       ! GULP atomic symbol (element+number)
    Integer                 :: index        ! the number from above
    Character(len=5)        :: id           ! 'core' or 'shell'
    Real(kind=RDbl)         :: charge
  End Type GULP_Atomic_Params

  Type GULP_Atomic_Info
    Integer                 :: natypes = 0     ! number of atom types
    Type(GULP_Atomic_Params), Dimension(:), Pointer  :: atype
  End Type GULP_Atomic_Info

  !** Holds either a filename which contains GULP control file parameter lines
  !** or parameters with which this module can construct the needed lines
  Type GULP_Forcefield_Input
    Character(len=lstrLen)          :: gulp_params_file
    Integer                         :: nsets
    Logical                         :: qok
    Type (GULP_Atomic_Info)         :: gulp_atoms
    Integer, Dimension(2)           :: intraoff_spc
    Type(SimCell_Params), Pointer   :: simcell  !** hack used with intramol int
    Character(len=lstrLen), Dimension(:), Pointer  :: params
    Integer, Dimension(:,:), Pointer               :: list
  End Type GULP_Forcefield_Input

  Type GULP_Forcefield_Info
    Character(len=StrLen)   :: spec_line
    Character(len=lStrLen)  :: params_line
  End Type GULP_Forcefield_Info

Contains

  !-----------------------------------------------------------------------
  ! Initialize the GULP Forcefield input data type.  
  ! Requires:  specs -- data type to initialize
  !            filename -- filename containing GULP FF spec lines or blank
  !            intraoff_spc -- list of species number to be without intra
  !            simcell -- simulation cell information
  !            list -- list of atom numbers for interaction (set,1:Natoms)
  !            params -- array of strings containing type and parameters
  !-----------------------------------------------------------------------
  Subroutine gulpio_init(specs,filename,intraoff_spc,simcell,list,params)
    Type(GULP_Forcefield_Input), Intent(Out)                    :: specs
    Character(len=lstrLen), Intent(In)                          :: filename
    Integer, Dimension(2), Intent(In)                           :: intraoff_spc
    Type(SimCell_Params), Intent(In), Target                    :: simcell
    Integer, Dimension(:,:), Intent(In), Optional               :: list
    Character(len=lstrLen), Dimension(:), Intent(Out), Optional :: params

    Integer         :: error,i,nsets

    If ((Present(params)).And.(filename /= '')) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " if params passed, filename must be empty"
      Stop      
    End If

    specs%nsets = 0
    specs%gulp_params_file = filename
    specs%intraoff_spc = intraoff_spc
    specs%simcell => simcell
    Nullify(specs%params)
    Nullify(specs%list)

    If (Present(params)) Then
      !** Initialize the GULP atom types
      Call gulpio_initatominfo(specs%gulp_atoms)

      !** Allocate storage space for the parameters
      nsets = Size(params)
      specs%nsets = nsets
      Allocate(specs%params(nsets),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'input_specs%params')
      Allocate(specs%list(nsets,3),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'input_specs%list')
      
      !** Copy the atom list and parameters into their new storage
      Do i = 1,nsets
        specs%list(i,:) = list(i,:)
        specs%params(i) = params(i)
      End Do
    Else
      !** Initialize the GULP atom types from filename
      Call gulpio_initatominfo(specs%gulp_atoms,specs%gulp_params_file)
    End If

  End Subroutine gulpio_init

  !-----------------------------------------------------------------------
  ! Performs a single point energy calculation using GULP.  Will
  ! attempt to generate the potentials specs for GULP from the existing
  ! MUSIC information if the filename is not present in gulp_specs
  ! Requires:  gulp_specs -- specifications for forcefield input
  !            qok -- flag indicating if GULP 'qok' flag should be used
  !            simcell -- simulation cell data structure
  !            natoms -- number of atoms in arrays
  !            coords -- coordinates of each atom
  !            atypes -- atom type of each atom
  !            charges -- charge of each atom
  !            subsets -- subset identification of each atom
  !            nrg -- resulting energy from GULP
  !            forces -- optional resulting forces on each atom
  !-----------------------------------------------------------------------
  Logical Function gulpio_snglpt(gulp_specs,qok,simcell,natoms, &
      coords,atypes,charges,subsets,nrg,forces)
    Type(GULP_Forcefield_Input), Intent(InOut)         :: gulp_specs
    Logical, Intent(In)                                :: qok
    Type(SimCell_Params), Intent(In)                   :: simcell  
    Integer, Intent(In)                                :: natoms
    Type(VecType), Dimension(:), Intent(In)            :: coords
    Integer, Dimension(:), Intent(In)                  :: atypes
    Real(kind=RDbl), Dimension(:), Intent(In)          :: charges
    Integer, Dimension(:,:), Intent(In)                :: subsets
    Real(kind=RDbl), Intent(Out)                       :: nrg  
    Type(VecType), Dimension(:), Intent(Out), Optional :: forces

    Integer                        :: unit,error,i,j
    Logical                        :: rflag
    Real(kind=RDbl)                :: number
    Character(len=lstrLen)         :: command,inputfile
#ifndef NAGCOMPILER
    Integer                        :: system,error_flag
#endif
    Type(GULPsystem)               :: gulp_system

    !** set defaults
    gulpio_snglpt = .False.
    nrg = 0.0_RDbl

    !** give feedback if -VERBOSE flag is used
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(*,'(1x,a)') 'Making external call to GULP:'
    End If

    !** create the input file
    inputfile = 'gulp_control_file'
    Call gulpio_snglptfile(inputfile,gulp_specs,simcell,coords(1:natoms), &
        atypes(1:natoms),charges(1:natoms),subsets(1:natoms,:),qok)

    !** set the command line
    Write(command,'(4a)') 'gulp < ',Trim(inputfile),' > ',Trim(gulpio_outfile)

    !** run GULP
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(*,'(1x,2a,i4,6a)') __FILE__," : ",__LINE__,' Running GULP ',&
          '(input: ',Trim(inputfile),') (output: ',Trim(gulpio_outfile),')...'
    End If
#ifndef NAGCOMPILER
    error_flag = system(command)
    If (error_flag < 0) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Error executing GULP code, please check output file: ", &
          Trim(gulpio_outfile)
      Stop
    End If
#endif
#ifdef NAGCOMPILER
    Call system(command)
#endif

    !** Analyze the output
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(*,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' Analyzing GULP output (',Trim(gulpio_outfile),')...'
    End If
    Call readgulp_getsystem(gulpio_outfile,gulp_system)
    nrg = readgulp_nrg(gulp_system,1)
    gulpio_snglpt = .True.

    !** Also extract the gradients, if desired
    If (Present(forces)) Then
      Call readgulp_forces(gulp_system,1,forces)
    End If

  End Function gulpio_snglpt

  !---------------------------------------------------------------------------
  ! Creates a GULP control file for a single point calculation.  Will
  ! attempt to generate the potentials specs for GULP from the existing
  ! MUSIC information if 'gulp_pots' filename is not present.
  ! Requires:  filename -- file to create
  !            gulp_specs -- specifications for forcefield input
  !            simcell -- simulation cell data structure
  !            coords -- coordinates of each atom
  !            atypes -- atom type of each atom
  !            charges -- charge of each atom
  !            subsets -- subset identification of each atom
  !            qok -- optional flag indicating if GULP qok flag should be used
  ! NOTE: the gulp_atom information must be initialized first
  !---------------------------------------------------------------------------
  Subroutine gulpio_snglptfile(filename,gulp_specs,simcell,coords, &
      atypes,charges,subsets,qok)
    Character(*), Intent(In)                     :: filename
    Type(GULP_Forcefield_Input), Intent(InOut)   :: gulp_specs
    Type(SimCell_Params), Intent(In)             :: simcell  
    Type(VecType), Dimension(:), Intent(In)      :: coords
    Integer, Dimension(:), Intent(In)            :: atypes
    Real(kind=RDbl), Dimension(:), Intent(In)    :: charges
    Integer, Dimension(:,:), Intent(In)          :: subsets
    Logical, Intent(In), Optional                :: qok

    Integer                        :: unit
    Character(len=xlstrLen)        :: keywords

    !** create the input file
    unit = file_open(filename)

    !** set default keywords: single point, constant volume, no symmetry
    keywords = 'single conv nosymm'

    !** add keywords if necessary
    If (Present(qok)) Then
      If (qok) Write(keywords,'(2a)') Trim(keywords),' qok'
    End If
    If (gulp_specs%intraoff_spc(1) /= 0) Then
      Write(keywords,'(2a)') Trim(keywords),' molecule'
    End If

    !** write keywords to file
    Write(unit,'(a)') Trim(keywords)
    Write(unit,*) 

    !** write coordinates to file
    Call gulpio_wrtcoords(filename,gulp_specs,simcell,coords,atypes,charges)

    !** Write forcefield to file
    Call gulpio_wrtffield(filename,gulp_specs,atypes,charges,subsets)

    Close(unit=unit)

  End Subroutine gulpio_snglptfile

  !-----------------------------------------------------------------------
  ! Dumps the MUSIC coordinate structure and simulation cell information
  ! into a GULP input file using the "fractional" format
  ! Requires:  filename -- GULP output filename
  !            specs -- GULP forcefield input data
  !            simcell -- simulation cell data structure
  !            coords -- coordinates of each atom
  !            atypes -- atom type of each atom
  !            charges -- charge of each atom
  !-----------------------------------------------------------------------
  Subroutine gulpio_wrtcoords(filename,specs,simcell,coords,atypes,charges)
    Character(*), Intent(In)                     :: filename
    Type(GULP_Forcefield_Input), Intent(InOut)   :: specs
    Type(SimCell_Params), Intent(In)             :: simcell  
    Type(VecType), Dimension(:), Intent(In)      :: coords
    Integer, Dimension(:), Intent(In)            :: atypes
    Real(kind=RDbl), Dimension(:), Intent(In)    :: charges

    Integer                        :: unit,error,i,a,m,spc
    Integer                        :: nspc,natoms
    Logical                        :: fractional
    Character(len=5)               :: id
    Character(len=4)               :: symbol
    Type(VecType)                  :: pt
    character(len=lstrLen)         :: string
    Real(kind=RDbl), Dimension(3)  :: ell,ellinv

    unit = file_open(filename,110)

    ell = simcell_getell(simcell)
    ellinv = 1.0_RDbl/ell

    !** write the simcell info (angs => angle format)
    Write(unit,'(a)') 'cell angs'  
    Write(unit,'(3f12.8,4x,3f12.8)') ell,simcell_angles(simcell,.True.)
    Write(unit,*) 
!    Write(unit,'(a)') 'spacegroup'
!    Write(unit,'(a)') '1'
    Write(unit,'(a)') 'origin 0'        !origin type 0
    Write(unit,*) 

    !** Actually dump the coordinates to file
    natoms = Size(coords)
    fractional = .False.
    If (fractional) Then   !** dump as fractional coordinates (not useful?)
      Write(unit,'(a)') 'fractional'
      Write(unit,'(i8)') natoms
      
      Do a = 1,natoms
        Do i = 1,3
          pt%comp(i) = coords(a)%comp(i)*ellinv(i)
        End Do
  
        !** Get the GULP atomic symbol and id
        Call gulpio_getatominfo(atypes(a),charges(a),specs%gulp_atoms,symbol,id)
  
        !** Write to file
        string = vector_display(pt,'f12.8')
        Write(unit,'(a,3x,a,3x,a)') Trim(symbol),Trim(id),Trim(string)
      End Do

    Else   !** dump as cartesian coordinates
      Write(unit,'(a)') 'cartesian angs'
      Do a = 1,natoms
        !** Get the GULP atomic symbol and id
        Call gulpio_getatominfo(atypes(a),charges(a),specs%gulp_atoms,symbol,id)
  
        !** Write to file
        string = vector_display(coords(a),'f12.8')
        Write(unit,'(a,3x,a,3x,a)') Trim(symbol),Trim(id),Trim(string)
      End Do
    End If

    Write(unit,*) 

  End Subroutine gulpio_wrtcoords

  !-----------------------------------------------------------------------
  ! Dumps the MUSIC forcefield information into a GULP input file.  Will
  ! attempt to generate the potentials specs for GULP from the existing
  ! MUSIC information if 'gulp_pots' filename is not present.
  ! Requires:  filename -- GULP output filename
  !            gulp_specs -- specifications for forcefield input
  !            atypes -- atom type of each atom
  !            charges -- charge of each atom
  !            subsets -- subset identification of each atom
  ! NOTE: the gulp_atom information must be initialized first
  !-----------------------------------------------------------------------
  Subroutine gulpio_wrtffield(filename,gulp_specs,atypes,charges,subsets)
    Character(*), Intent(In)                     :: filename
    Type(GULP_Forcefield_Input), Intent(InOut)   :: gulp_specs
    Integer, Dimension(:), Intent(In)            :: atypes
    Real(kind=RDbl), Dimension(:), Intent(In)    :: charges  
    Integer, Dimension(:,:), Intent(In)          :: subsets
    
    Integer                  :: atype,unit,readunit,ios
    Real(kind=RDbl)          :: charge
    Character(len=2)         :: symbol
    Character(len=strLen)    :: id
    Character(len=xlstrLen)  :: line

    unit = file_open(filename,110)

    !** write the 'species' information
    Write(unit,'(a)') 'species'  
    Do atype = 1,gulp_specs%gulp_atoms%natypes
      symbol = gulp_specs%gulp_atoms%atype(atype)%symbol
      id = gulp_specs%gulp_atoms%atype(atype)%id
      charge = gulp_specs%gulp_atoms%atype(atype)%charge
      Write(unit,'(a,3x,a,3x,f9.5)') Trim(symbol),Trim(id),charge
    End Do
    Write(unit,*)

    !** Add the section to adjust covalent radii if necessary
    Call gulpio_adjustradii(gulp_specs,unit)

    !** write the potential parameter information
    If (Associated(gulp_specs%params)) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " WARNING, please check correctness of GULP forcefield in ", &
          Trim(filename)
      Call gulpio_wrt2body(filename,gulp_specs,atypes,charges,subsets)
      Call gulpio_wrtintrapair(filename,gulp_specs,atypes,charges,subsets)
      Call gulpio_wrt3body(filename,gulp_specs,atypes,charges,subsets)
!      Call gulpio_wrt4body(filename,gulp_specs,atypes,charges,subsets)
    Else
      readunit = file_open(gulp_specs%gulp_params_file,110)
      Do 
        Read(readunit,'(a)',IOSTAT=ios) line
        If (ios /= 0) Exit
        Write(unit,'(a)') Trim(line)
      End Do
      Close(unit=readunit)
    End If

  End Subroutine gulpio_wrtffield

  !-----------------------------------------------------------------------
  ! Dumps the GULP forcefield section that adjusts element covalent radii.
  ! Adjust the covalent radii if we're using the 'molecule' keyword.
  ! We may want to make sure that some species are not defined as 
  ! molecules and therefore do not have their intramolecular charge
  ! interactions substracted out.  This is important for shell-models
  ! Requires:  gulp_specs -- specifications for forcefield input
  !            unit -- unit to write to, if necessary
  ! NOTE: the gulp_atom information must be initialized first
  !-----------------------------------------------------------------------
  Subroutine gulpio_adjustradii(gulp_specs,unit)
    Type(GULP_Forcefield_Input), Intent(In)      :: gulp_specs
    Integer, Intent(In)                          :: unit
    
    Integer                             :: spc,nelements,natypes,a,i
    Logical                             :: match
    Character(len=2)                    :: symbol
    Integer, Dimension(100)             :: alist
    Character(len=2), Dimension(100)    :: element_list

    !** leave if this is not appropriate
    If (gulp_specs%intraoff_spc(1) == 0) Return

    !** Add elements that are in the species that should maintain
    !** their intramolecular charge interactions
    Do spc = 1,molecules_getnsorbs()
      If (gulp_specs%intraoff_spc(1) == spc) Cycle
      If (gulp_specs%intraoff_spc(2) == spc) Cycle

      !** get atom types in this species type
      natypes = molecules_getnatomtypes(spc,alist,.False.)

      !** put these in a list if they aren't repeats
      nelements = 0
      Do i = 1,natypes
        symbol = atom_getsymbol(alist(i))
        match = .False.
        Do a = 1,nelements
          If (element_list(a) == symbol) match = .True.
        End Do
        If (match) Cycle
        nelements = nelements + 1
        element_list(nelements) = symbol
      End Do
    End Do

#ifdef DEBUG
    Do a = 1,nelements
      Write(*,*) a,element_list(a)
    End Do
#endif

    !** Check to make sure we didn't add any elements that are in
    !** the intramolecular-charge-free species
    Do spc = 1,molecules_getnsorbs()
      If ((gulp_specs%intraoff_spc(1) == spc).Or. &
          (gulp_specs%intraoff_spc(2) == spc)) Then
        Do a = 1,nelements
          Do i = 1,gulp_specs%gulp_atoms%natypes
            If (element_list(a) == gulp_specs%gulp_atoms%atype(i)%element) Then 
              Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
                  ' Could not zero covalent radius of element ',element_list(a)
              Stop
            End If
          End Do
        End Do
      End If
    End Do

    !** Actually do the writing to the file
    Write(unit,'(a)') 'element'
    Do a = 1,nelements
      Write(unit,'(a,2x,a,f8.3)') 'covalent',element_list(a),0.0_RDbl
    End Do
    Write(unit,'(a)') 'end'
    Write(unit,*)      

  End Subroutine gulpio_adjustradii

  !-------------------------------------------------------------------------
  ! Dumps the MUSIC forcefield intra-pair information into a GULP input file
  ! Requires:  filename -- GULP output filename
  !            specs -- specifications for forcefield input
  !            atypes -- atom type of each atom
  !            charges -- charge of each atom
  !            subsets -- subset identification of each atom
  ! NOTE: the gulp_atom information must be initialized first
  !-----------------------------------------------------------------------
  Subroutine gulpio_wrtintrapair(filename,specs,atypes,charges,subsets)
    Character(*), Intent(In)                     :: filename
    Type(GULP_Forcefield_Input), Intent(InOut)   :: specs
    Integer, Dimension(:), Intent(In)            :: atypes
    Real(kind=RDbl), Dimension(:), Intent(In)    :: charges  
    Integer, Dimension(:,:), Intent(In)          :: subsets

    Integer                               :: i,a1,a2,param_spc,natoms
    Integer                               :: nfields,nnew_sets
    Logical                               :: test_flag
    Type(GULP_Forcefield_Info)            :: trial
    Character(len=4)                      :: asym1,asym2
    Character(len=5)                      :: id1,id2
    Character(len=strLen), Dimension(100)      :: fields
    Type(GULP_Forcefield_Info), Dimension(100) :: fflist

    natoms = Size(atypes)
    nnew_sets = 0
    
    !** handle intrapair potentials
    Do i = 1,specs%nsets
      !** split the parameters line and skip if not pertinent
      nfields = split(specs%params(i),fields)
      If (fields(1) /= INTRAPAIR_KEY) Cycle
      param_spc = molecules_gettype(fields(2))
      If (param_spc == 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Unable to associate species name with number ',Trim(fields(2))
        Stop
      End If

      !** convert molecular atom numbers to indices in 'atypes' and 'charges'
      a1 = gulpio_idsubset(subsets,param_spc,specs%list(i,1))
      If (a1 == 0) Cycle
      a2 = gulpio_idsubset(subsets,param_spc,specs%list(i,2))
      If (a2 == 0) Cycle

      !** get the GULP atom IDs
      Call gulpio_getatominfo(atypes(a1),charges(a1), &
          specs%gulp_atoms,asym1,id1)
      Call gulpio_getatominfo(atypes(a2),charges(a2), &
          specs%gulp_atoms,asym2,id2)

      !** interpret the potential set parameters
      Select Case (Trim(fields(3)))  !* see ipmodel_displayParams
      Case('Buck')
        !** decide if it's a core-shell or normal harmonic
        If ((id1 == 'shell').Or.(id2 == 'shell')) Then
          trial%spec_line = 'buckingham kcal'
          Write(trial%params_line,'(7(a,2x),f8.3)') asym1,id1,asym2,id2, &
              Trim(fields(5)),Trim(fields(8)),Trim(fields(11)),10.0_RDbl
        Else 
          trial%spec_line = 'buckingham bond kcal'
          Write(trial%params_line,'(7(a,2x))') asym1,id1,asym2,id2, &
              Trim(fields(5)),Trim(fields(8)),Trim(fields(11))
        End If
      
      Case Default
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not interpret ID string: ",Trim(fields(3))
        Stop
      End Select

      !** compare with old forcefield lines and add if it's new
      test_flag = gulpio_matchff(trial,nnew_sets,fflist)
    End Do

    !** dump parameters to file
    Call gulpio_dumpff(filename,nnew_sets,fflist)

  End Subroutine gulpio_wrtintrapair

  !-----------------------------------------------------------------------
  ! Dumps the MUSIC forcefield 2-BODY information into a GULP input file
  ! Requires:  filename -- GULP output filename
  !            specs -- specifications for forcefield input
  !            atypes -- atom type of each atom
  !            charges -- charge of each atom
  !            subsets -- subset identification of each atom
  ! NOTE: the gulp_atom information must be initialized first
  !-----------------------------------------------------------------------
  Subroutine gulpio_wrt2body(filename,specs,atypes,charges,subsets)
    Character(*), Intent(In)                     :: filename
    Type(GULP_Forcefield_Input), Intent(InOut)   :: specs
    Integer, Dimension(:), Intent(In)            :: atypes
    Real(kind=RDbl), Dimension(:), Intent(In)    :: charges  
    Integer, Dimension(:,:), Intent(In)          :: subsets
    
    Integer                               :: i,a1,a2,param_spc,natoms
    Integer                               :: nfields,nnew_sets
    Logical                               :: test_flag
    Type(GULP_Forcefield_Info)            :: trial
    Character(len=4)                      :: asym1,asym2
    Character(len=5)                      :: id1,id2
    Character(len=strLen), Dimension(100)      :: fields
    Type(GULP_Forcefield_Info), Dimension(100) :: fflist

    natoms = Size(atypes)
    nnew_sets = 0
    
    !** handle bond stretch potentials
    Do i = 1,specs%nsets
      !** split the parameters line and skip if not pertinent
      nfields = split(specs%params(i),fields)
      If (fields(1) /= STRETCH_KEY) Cycle
      param_spc = molecules_gettype(fields(2))
      If (param_spc == 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Unable to associate species name with number ',Trim(fields(2))
        Stop
      End If

      !** convert molecular atom numbers to indices in 'atypes' and 'charges'
      a1 = gulpio_idsubset(subsets,param_spc,specs%list(i,1))
      If (a1 == 0) Cycle
      a2 = gulpio_idsubset(subsets,param_spc,specs%list(i,2))
      If (a2 == 0) Cycle

      !** get the GULP atom IDs
      Call gulpio_getatominfo(atypes(a1),charges(a1), &
          specs%gulp_atoms,asym1,id1)
      Call gulpio_getatominfo(atypes(a2),charges(a2), &
          specs%gulp_atoms,asym2,id2)
    
      !** interpret the potential set parameters
      Select Case (Trim(fields(3)))
      Case('Harmonic')
        !** decide if it's a core-shell or normal harmonic
        If ((asym1 == asym2).And.(id1 /= id2)) Then
          trial%spec_line = 'spring kcal'
          Write(trial%params_line,'(a,3x,a)') asym1,Trim(fields(4))
        Else
          trial%spec_line = 'harmonic bond kcal'
          Write(trial%params_line,'(5(a,2x))') asym1,id1,asym2,id2, &
              Trim(fields(4))
        End If
    
      Case Default
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not interpret ID string: ",Trim(fields(3))
        Stop
      End Select

      !** compare with old forcefield lines and add if it's new    
      test_flag = gulpio_matchff(trial,nnew_sets,fflist)
    End Do

    !** dump parameters to file
    Call gulpio_dumpff(filename,nnew_sets,fflist)

  End Subroutine gulpio_wrt2body

  !-------------------------------------------------------------------------
  ! Dumps the MUSIC forcefield 3-BODY information into a GULP input file
  ! Requires:  filename -- GULP output filename
  !            specs -- specifications for forcefield input
  !            atypes -- atom type of each atom
  !            charges -- charge of each atom
  !            subsets -- subset identification of each atom
  ! NOTE: the gulp_atom information must be initialized first
  !-------------------------------------------------------------------------
  Subroutine gulpio_wrt3body(filename,specs,atypes,charges,subsets)
    Character(*), Intent(In)                     :: filename
    Type(GULP_Forcefield_Input), Intent(InOut)   :: specs
    Integer, Dimension(:), Intent(In)            :: atypes
    Real(kind=RDbl), Dimension(:), Intent(In)    :: charges  
    Integer, Dimension(:,:), Intent(In)          :: subsets
    
    Integer                               :: i,a1,a2,a3,param_spc,natoms
    Integer                               :: nfields,nnew_sets
    Logical                               :: test_flag
    Type(GULP_Forcefield_Info)            :: trial
    Character(len=4)                      :: asym1,asym2,asym3
    Character(len=5)                      :: id1,id2,id3
    Character(len=strLen), Dimension(100)      :: fields
    Type(GULP_Forcefield_Info), Dimension(100) :: fflist

    natoms = Size(atypes)
    nnew_sets = 0

    !** handle 3-body potentials
    Do i = 1,specs%nsets
      !** split the parameters line and skip if not pertinent
      nfields = split(specs%params(i),fields)
      If (fields(1) /= BENDING_KEY) Cycle    
      param_spc = molecules_gettype(fields(2))
      If (param_spc == 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Unable to associate species name with number ',Trim(fields(2))
        Stop
      End If

      !** convert molecular atom numbers to indices in 'atypes' and 'charges'
      a1 = gulpio_idsubset(subsets,param_spc,specs%list(i,1))
      If (a1 == 0) Cycle
      a2 = gulpio_idsubset(subsets,param_spc,specs%list(i,2))
      If (a2 == 0) Cycle
      a3 = gulpio_idsubset(subsets,param_spc,specs%list(i,3))
      If (a3 == 0) Cycle

      !** get the GULP atom IDs
      Call gulpio_getatominfo(atypes(a1),charges(a1), &
          specs%gulp_atoms,asym1,id1)
      Call gulpio_getatominfo(atypes(a2),charges(a2), &
          specs%gulp_atoms,asym2,id2)
      Call gulpio_getatominfo(atypes(a3),charges(a3), &
          specs%gulp_atoms,asym3,id3)
    
      !** interpret the potential set parameters
      Select Case (Trim(fields(3)))  !* see bbmodel_displayParams
      Case('Harmonica')
        !** decide if it's a core-shell or normal harmonic
        If ((id1 == 'shell').Or.(id2 == 'shell').Or.(id3 == 'shell')) Then
          trial%spec_line = 'three-body kcal'
          Write(trial%params_line,'(8(a,2x),3f8.3)') asym1,id1,asym2,id2, &
              asym3,id3,Trim(fields(4)),Trim(fields(5)),2.5,2.5,4.0
        Else 
          trial%spec_line = 'three-body bond kcal'
          Write(trial%params_line,'(8(a,2x),3f8.3)') asym1,id1,asym2,id2, &
              asym3,id3,Trim(fields(4)),Trim(fields(5)),2.5,2.5,4.0
        End If
    
      Case Default
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not interpret ID string: ",Trim(fields(3))
        Stop
      End Select

      !** compare with old forcefield lines and add if it's new    
      test_flag = gulpio_matchff(trial,nnew_sets,fflist)
    End Do

    !** Dump parameters to file
    Call gulpio_dumpff(filename,nnew_sets,fflist)

  End Subroutine gulpio_wrt3body

  !-------------------------------------------------------------------------
  ! Returns the index in the array of subsets that matches a certain species
  ! and atom number.  Will return zero if it can't find a match.
  ! Requires:  subset -- array of subsets
  !            spc -- species number
  !            atm -- atom number
  !-----------------------------------------------------------------------
  Integer Function gulpio_idsubset(subsets,spc,atm)
    Integer, Dimension(:,:), Intent(In)  :: subsets
    Integer, Intent(In)                  :: spc,atm

    Integer         :: i

    gulpio_idsubset = 0

    Do i = 1,Size(subsets)
      If (subsets(i,1) /= spc) Cycle
      If (subsets(i,3) /= atm) Cycle

      !** if we're still in loop, we've found a match
      gulpio_idsubset = i
      Return
    End Do

  End Function gulpio_idsubset

  !-----------------------------------------------------------------------
  ! Dumps the information in the GULP forcefield structure to file
  ! Requires: filename -- GULP output filename
  !           nsets -- number of forcefield terms in list
  !           fflist -- list of forcefield terms
  !-----------------------------------------------------------------------
  Subroutine gulpio_dumpff(filename,nsets,fflist)
    Character(*), Intent(In)                             :: filename
    Integer, Intent(In)                                  :: nsets
    Type(GULP_Forcefield_Info), Dimension(:), Intent(In) :: fflist

    Integer                :: unit,i

    unit = file_open(filename,110)
    Write(unit,*) 

    If (nsets == 0) Return

    !** dump the parameters to the file
    Write(unit,'(a)') Trim(fflist(1)%spec_line)
    Write(unit,'(a)') Trim(fflist(1)%params_line)
    Do i = 2,nsets
      If (fflist(i)%spec_line /= fflist(i-1)%spec_line) Then
        Write(unit,'(a)') Trim(fflist(i)%spec_line)
      End If
      Write(unit,'(a)') Trim(fflist(i)%params_line)
    End Do

  End Subroutine gulpio_dumpff

  !-----------------------------------------------------------------------
  ! Handles the compilation of forcefield parameters.  Checks for a match
  ! in the existing list.  If no match is found, inserts the new info set
  ! into the list.
  ! Requires:  trial -- new info set to match
  !            nsets -- current number of sets
  !            list -- list of existing sets
  !-----------------------------------------------------------------------
  Logical Function gulpio_matchff(trial,nsets,list)
    Type(GULP_Forcefield_Info), Intent(In)                  :: trial
    Integer, Intent(InOut)                                  :: nsets
    Type(GULP_Forcefield_Info), Dimension(:), Intent(InOut) :: list

    Integer        :: i

    gulpio_matchff = .False.
    If (nsets > 0) Then
      Do i = 1,nsets
        If ((Trim(trial%spec_line) == Trim(list(i)%spec_line)).And. &
           (Trim(trial%params_line) == Trim(list(i)%params_line))) Then
          gulpio_matchff = .True.
          Exit
        End If
      End Do
    End If 

    If (.Not. gulpio_matchff) Then
      nsets = nsets + 1
      list(nsets)%spec_line = trial%spec_line
      list(nsets)%params_line = trial%params_line
    End If

  End Function gulpio_matchff

  !--------------------------------------------------------------------------
  ! Initialize the GULP atom types.  If the optional filename is passed, then
  ! it reads the GULP 'species' section and builds an atom type list based
  ! on this specifications.  This helps ensure that the coordinates that
  ! are dumped into the GULP control file use matching atom symbols.
  ! Requires:  gulp_atoms -- GULP atomic information structure
  !            filename -- optional filename to build type info from
  !--------------------------------------------------------------------------
  Subroutine gulpio_initatominfo(gulp_atoms,filename)
    Type(GULP_Atomic_Info), Intent(Out)   :: gulp_atoms
    Character(*), Intent(In), Optional    :: filename

    Integer                               :: error,lineno,unit,a,i
    Integer                               :: nfields
    Logical                               :: iontype_match
    Real(kind=RDbl)                       :: atype_charge
    Character(len=2)                      :: atype_element
    Character(len=5)                      :: atype_id
    Character(len=strLen)                 :: element,index
    character(len=lstrLen)                :: line       
    Character(len=strLen), Dimension(20)  :: fields

    !** Allocate space 
    Allocate(gulp_atoms%atype(1000),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'gulp_atom%atype')
    gulp_atoms%natypes = 0

    If (Present(filename)) Then
      unit = file_open(filename,110)
      lineno = fileSrchStr(unit,'species',line,.True.)
      If (lineno == 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not find 'species' section in: ",Trim(filename)
        Stop        
      End If

      !** read and interpret the atomic species lines
      a = 0
      Do 
        Read(unit,'(a)') line
        nfields = split(line,fields)
        If (nfields /= 3) Exit
        a = a + 1

        !** store the values on the line
        gulp_atoms%atype(a)%symbol = fields(1)
        index = digits(gulp_atoms%atype(a)%symbol)
        If (index == '') Then
          gulp_atoms%atype(a)%index = 0
        Else
          gulp_atoms%atype(a)%index = toint(index)
        End If
        element = firstchars(fields(1))
        gulp_atoms%atype(a)%element = element 
        Select Case(Toupper(fields(2)(1:4)))
        Case('CORE')
          gulp_atoms%atype(a)%id = 'core'
        Case('SHEL')
          gulp_atoms%atype(a)%id = 'shell'
        Case Default
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              " Could not interpret core/shell specification ",Trim(fields(2))
          Stop        
        End Select
        gulp_atoms%atype(a)%charge = toreal(fields(3))

        !** Match to a MuSiC atom type
        gulp_atoms%atype(a)%atom_name = ''
        Do i = 1,atom_getntypes()
          atype_element = atom_getsymbol(i)
          atype_id = atom_iontype(i)
          atype_charge = atom_getcharge(i)
          iontype_match = .False.
          If ((ToUpper(gulp_atoms%atype(a)%id) == ToUpper(atype_id)).Or. &
              ((ToUpper(gulp_atoms%atype(a)%id) == 'CORE').And. &
              (ToUpper(atype_id) == 'NONE'))) iontype_match = .True.

          !** Check element, iontype and charge for match
          If ((ToUpper(element) == ToUpper(atype_element)).And. &
              (iontype_match).And. &
              (gulp_atoms%atype(a)%charge == atype_charge)) Then
            gulp_atoms%atype(a)%atom_name = atom_getname(i)
            gulp_atoms%atype(a)%atom_index = i
            gulp_atoms%atype(a)%orig_element = atom_getsymbol(i)
            gulp_atoms%natypes = gulp_atoms%natypes + 1
            Exit
          End If
        End Do

        If (gulp_atoms%atype(a)%atom_name == '') Then
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              " Could not find MuSiC atom type match for ",Trim(line)
          Write(0,'(1x,a)') "Please make sure one of the atom types matches"
          Write(*,*) ToUpper(element),ToUpper(gulp_atoms%atype(a)%id), &
              gulp_atoms%atype(a)%charge
          Stop                
        End If

      End Do
    End If

  End Subroutine gulpio_initatominfo

  !--------------------------------------------------------------------------
  ! Creates new atomic symbols for use in the GULP forcefield statements.
  ! For example, GULP distinguishes O, O1, and O2.  This routine ensures that
  ! each MUSIC atom that has a different name, charge (from molecules) and
  ! element symbol.  It maintains the gulp_atoms database. 
  ! Requires: atype -- atom type
  !           charge -- atomic charge
  !           gulp_atoms -- GULP atomic information structure
  !           element -- returned atomic symbol
  !           id -- 'shell' or 'core' designation
  ! The gulp_atoms%natypes variable MUST have a value before this routine
  ! is called.  Set gulp_atoms%natypes = 0 before this routine is called
  ! for the first time.
  !--------------------------------------------------------------------------
  Subroutine gulpio_getatominfo(atype,charge,gulp_atoms,symbol,id)
    Integer, Intent(In)                          :: atype
    Real(kind=RDbl)                              :: charge
    Type(GULP_Atomic_Info), Intent(InOut)        :: gulp_atoms
    Character(len=4), Intent(Out)                :: symbol
    Character(len=5), Intent(Out)                :: id

    Integer                 :: error,match,a,num
    Character(len=2)        :: new_element
    Character(len=4)        :: newsymbol
    Character(len=strLen)   :: aelement,aname,string

    !** Allocate space if this is the first time
    If (gulp_atoms%natypes == 0) Then
      Allocate(gulp_atoms%atype(1000),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'gulp_atom%atype')
    End If

    !** get atom type, charge, name and atomic symbol of the atom
    aname = atom_getname(atype)
    aelement = atom_getsymbol(atype)

    !** look for a match in the recorded types
    match = 0
    Do a = 1,gulp_atoms%natypes
      !** check name, atomic symbol and charge for match
!LC      Write(*,*) 'name: ','(',Trim(aname),')(', &
!LC          Trim(gulp_atoms%atype(a)%atom_name),')'
      If (Trim(aname) /= Trim(gulp_atoms%atype(a)%atom_name)) Cycle
!LC      Write(*,*) 'element: ',aelement,gulp_atoms%atype(a)%element
      If (Trim(aelement) /= Trim(gulp_atoms%atype(a)%orig_element)) Cycle
!LC      Write(*,*) 'charge: ',charge,gulp_atoms%atype(a)%charge
      If (Abs(charge - gulp_atoms%atype(a)%charge) > zeroTolerance) Cycle

      !** match found
      match = a
!LC      Write(*,*) 'MATCH'
      Exit
    End Do

    !** No match, create a new symbol
    If (match == 0) Then
      If (gulp_atoms%natypes >= Size(gulp_atoms%atype)) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Not enough storage space for atom types"
        Stop        
      End If

      !** find the largest index of this atomic number
      num = 0
      Do a = 1,gulp_atoms%natypes
        If (Trim(gulp_atoms%atype(a)%element) == Trim(aelement)) Then
          If (gulp_atoms%atype(a)%index > num) num = gulp_atoms%atype(a)%index
        End If
      End Do

      !** create the new GULP atom type
      num = num + 1
      gulp_atoms%natypes = gulp_atoms%natypes + 1
      gulp_atoms%atype(gulp_atoms%natypes)%index = num
      gulp_atoms%atype(gulp_atoms%natypes)%atom_name = aname
      gulp_atoms%atype(gulp_atoms%natypes)%charge = charge

      !** make the element symbol and id
      Select Case(ToUpper(aelement))
      Case('OS')
        new_element = 'O'
        gulp_atoms%atype(gulp_atoms%natypes)%id = 'shell'
      Case('OC')
        new_element = 'O'
        gulp_atoms%atype(gulp_atoms%natypes)%id = 'core'
      Case Default
        new_element = aelement
        gulp_atoms%atype(gulp_atoms%natypes)%id = 'core'
      End Select 

!LC      Write(*,*) aelement,symbol

      gulp_atoms%atype(gulp_atoms%natypes)%orig_element = aelement
      gulp_atoms%atype(gulp_atoms%natypes)%element = new_element
      string = int2str(num)
      Write(newsymbol,'(2a)') Trim(new_element),Trim(string)
      gulp_atoms%atype(gulp_atoms%natypes)%symbol = newsymbol
      match = gulp_atoms%natypes

    End If

    symbol = gulp_atoms%atype(match)%symbol
    id = gulp_atoms%atype(match)%id

  End Subroutine gulpio_getatominfo

  !----------------------------------------------------------------------------
  ! Display the parameters
  ! Requires:  specs -- specifications for forcefield input
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine gulpio_display(specs,indent,unitno)
    Type(GULP_Forcefield_Input), Intent(In)  :: specs
    Integer, Intent(In)                      :: indent
    Integer, Intent(In), Optional            :: unitno

    Integer                     :: unit
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)
    unit = 6
    If (Present(unitno)) unit = unitno

    If (Associated(specs%list)) Then
      Write(unit,'(2a)') blank,'GULP used with MuSiC-based parameters'
    Else
      Write(unit,'(3a)') blank,'GULP used with parameters from file: ', &
          Trim(specs%gulp_params_file)
    End If

  End Subroutine gulpio_display

  !----------------------------------------------------------------------------
  ! Clean the parameters
  ! Requires:  specs -- specifications for input
  !----------------------------------------------------------------------------
  Subroutine gulpio_clean(specs)
    Type(GULP_Forcefield_Input), Intent(InOut)        :: specs

    Integer         :: error

    If (Associated(specs%simcell)) Then
      Deallocate(specs%simcell, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'simcell')
    End If
    If (Associated(specs%params)) Then
      Deallocate(specs%params, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'params')
    End If
    If (Associated(specs%list)) Then
      Deallocate(specs%list, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'list')
    End If
    If (Associated(specs%gulp_atoms%atype)) Then
      Deallocate(specs%gulp_atoms%atype, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'atypes')
    End If

  End Subroutine gulpio_clean


End Module gulpio


