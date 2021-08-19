!------------------------------------------------------------------------
! This module handles aspects of the eigenmode perturbation move.  The
! move makes perturbations along a select set of eigenmodes from the 
! reference configuration stored in molecules.  It requires a Hessian as
! input.  The eigenvalues associated with the rotational and translational 
! eigenmodes of this Hessian are zeroed.  Consequently, this move type 
! is capable of perturbing only the internal degrees of freedom, assuming 
! that the correct modes are choosen.
!------------------------------------------------------------------------

Module emperturb

  Use defaults, Only: RDbl, strLen, lstrLen, MAX_ATOMS
  Use utils, Only:  filesrchstr, split, toupper, allocErrDisplay, real2str, &
      int2str, deallocErrDisplay, str2seq, stripcmnt, combine
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), &
      Operator(-),vector_display
  Use config, Only: AtMolCoords
  Use simcell, Only: SimCell_Params,simcell_pbc
  Use molecules, Only: molecules_getcomdefn,molecules_getlinewtag, &
      molecules_getfilename,molecules_getnatoms,molecules_atomsymbol, &
      molecules_name
  Use hessian, Only: HessianInfo,hessian_removert,hessian_eigeninfo, &
      hessian_getmode,hessian_dispmode,hessian_dispevalues
  Use readhess, Only: Hessian_Specs,readhess_init,readhess_display, &
      readhess_clean
  Use subinteract, Only: Subset_Interactions,subinteract_int
  Use random, Only: rranf
  Use rigidcoords, Only: rigidcoords_changeref,rigidcoords_toxyz
  Use visxyz, Only: XYZ_Entry,visxyz_make

  Implicit None
  Save

  Private
  Public :: EMPerturb_Params,emperturb_init,emperturb_move,emperturb_display, &
      emperturb_clean

  Type EMPerturb_Params
    Integer                       :: spc,natoms,nmodes2perturb,emode_dim
    Real(kind=RDbl)               :: max_disp,max_nrg
    Type(HessianInfo)             :: hessian        !* Hessian to be used
    Type(Hessian_Specs)           :: hessian_specs  !* Hessian input method
    Type(VecType), Dimension(:), Pointer     :: ref_struc !* reference structure
    Integer, Dimension(:), Pointer           :: emodes    !* indices in Hessian
    Real(kind=RDbl), Dimension(:), Pointer   :: evalues
    Real(kind=RDbl), Dimension(:,:), Pointer :: evectors  !* in columns
  End Type EMPerturb_Params

Contains
  !---------------------------------------------------------------------------
  ! Initializes Eigenmode Perturbation move for the given species type
  ! Requires:  params -- EM perturb parameters
  !            spc -- species number
  !            species -- species data structure
  !            filename -- filename where initialization data can be found
  !---------------------------------------------------------------------------
  Subroutine emperturb_init(params,spc,species,filename)
    Type(EMPerturb_Params), Intent(InOut)       :: params
    Integer, Intent(In)                         :: spc
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Character(*), Intent(In)                    :: filename            

    Integer                                  :: i,mode
    Integer                                  :: unit,nfields,nchunks,error
    Character(len=lstrLen)                   :: string,specs_string
    Character(len=lstrLen)                   :: line,tag,xyzfilename
    Character(len=strLen), Dimension(strLen) :: fields,chunks
    Integer, Dimension(MAX_ATOMS*3)          :: nums
    Type(XYZ_Entry), Dimension(:), Pointer   :: xyzcoords

    !** Set the easy stuff
    params%spc = spc
    params%natoms = molecules_getnatoms(spc)
    params%max_nrg = 10000.0_RDbl  !** maximum energy HACK for now

    !** Get the reference coordinates from MuSiC species storage
    Allocate(params%ref_struc(params%natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Call molecules_getcomdefn(spc,params%ref_struc)

    !** Open file for reading
    unit = file_open(filename,110)

    !** Get the modes to perturb
    Read(unit,'(a)') line
    params%nmodes2perturb = str2seq(stripcmnt(line),nums)
    Allocate(params%emodes(params%nmodes2perturb), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    params%emodes(1:params%nmodes2perturb) = nums(1:params%nmodes2perturb)
    
    !** Get the maximum displacement along any eigenvector from reference
    Read(unit,*) params%max_disp

    !** Get the location of the hessian specifications
    Read(unit,*) line

    !** Split the line
    nfields = split(line,fields)

    !** Get the hessian location specifications
    !** Specs line will look like:   Hessian_Info: SPECIFICATIONS
    Select Case(Toupper(Trim(fields(1))))
    Case('INMOLECULE')  !** the specs line is in molecule file

      tag = 'Hessian_Info:'
      If (.Not. molecules_getlinewtag(spc,tag,line)) Then
        Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
            ' Could not find "',Trim(tag),'" tag in molecule file: ', &
            Trim(molecules_getfilename(spc))
        Stop
      End If

      nfields = split(line,fields)
      specs_string = fields(2)

    Case('HERE')  !** the specs line is current line here in control file
      specs_string = fields(2)
      
    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Could not identify keyword as hessian specification: ', &
          Trim(fields(1))
      Stop
    End Select

    !** Get the hessian and compare species against its generating structure
    nchunks = split(specs_string,chunks,'@')
    If (nchunks == 2) Then
      string = combine(chunks(1:2),'@')
      Call readhess_init(params%hessian_specs,params%hessian,string,spc)
    Else If (nchunks == 3) Then
      string = combine(chunks(1:2),'@')
      Call readhess_init(params%hessian_specs,params%hessian,string, &
          spc,chunks(3))
    Else 
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Wrong format for Hessian location specification: ', &
          Trim(fields(2))
      Write(0,'(2x,2a)') ' Format has 2 or 3 @-separated fields and is ', &
          'HESSIAN_FILETYPE@FILENAME[:atom_numbers][@COORDS_FILENAME]'
      Stop
    End If

    !** Remove rotational and translation components and diagonalize the Hessian
    Call hessian_removert(params%hessian,params%ref_struc)
    Call hessian_eigeninfo(params%hessian)

    !** Allocate memory for the eigenmodes to be perturbed
    Allocate(params%evalues(params%nmodes2perturb), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    params%emode_dim = params%natoms*3
    Allocate(params%evectors(params%emode_dim,params%nmodes2perturb), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
#ifdef DEBUG
    Allocate(xyzcoords(params%natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
#endif

    !** Get and store the desired eigenmodes
    Do i = 1,params%nmodes2perturb
      mode = params%emodes(i)
      Call hessian_getmode(params%hessian,mode,params%evalues(i), &
          params%evectors(1:params%emode_dim,i))

#ifdef DEBUG
      !** make xyz movies of the eigenmodes
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      xyzfilename = 'eigenmode'//Trim(int2str(mode))//'_movie.xyz'
      Do a = 1,params%natoms
        xyzcoords(a) = visxyz_make(params%ref_struc(a), &
            molecules_atomsymbol(spc,a))
      End Do
      Call hessian_dispmode(params%hessian,xyzcoords,mode,xyzfilename)
#endif
    End Do

  End Subroutine emperturb_init

  !---------------------------------------------------------------------------
  ! Performs a the internal transformation by displacing the molecular 
  ! geometry along one or more non-rotational, non-translational eigenmodes.
  ! Requires:  params -- eigenmode perturbation parameters
  !            spc -- species number
  !            molec -- molecule number
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias factor associated with move
  !            succ_flag -- returns true if move and energy evaluation went ok
  !            subint -- subset interactions (molecule-system)
  !---------------------------------------------------------------------------
  Subroutine emperturb_move(params,spc,molec,species,simcell, &
      biasfactor,succ_flag,subint)
    Type(EMPerturb_Params), Intent(In)             :: params
    Integer, Intent(In)                            :: spc,molec
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Real(kind=RDbl), Intent(Out)                   :: biasfactor
    Logical, Intent(Out)                           :: succ_flag
    Type(Subset_Interactions), Intent(InOut)       :: subint

    Integer                                  :: mode,atom,idx
    Logical                                  :: fast,skip_intra_flag
    Type(VecType), Dimension(params%natoms)  :: newref_struc
    Real(kind=RDbl), Dimension(params%emode_dim,params%nmodes2perturb) :: vecs

    succ_flag = .False.
    biasfactor = 1.0_RDbl  !** currently, no biasing

    !** get the old energies
    fast = .True.
    skip_intra_flag = .False.

    succ_flag = subinteract_int(subint,species,simcell,fast, &
      .False.,skip_intra_flag,(/params%max_nrg/),(/spc,molec,0/))
    If (.Not. succ_flag) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "Problem with nrg calculation of an existing structure"
      Stop
    End If

    !** Pick eigenmodes
    vecs = params%evectors

    !** Pick displacements along eigenmode(s) for perturbation
    Do mode = 1,params%nmodes2perturb
      vecs(:,mode) = vecs(:,mode) * (rranf()*params%max_disp)
    End Do

    !** Copy the reference structure
    newref_struc = params%ref_struc

    !** Perform the displacements on the new reference structure
    Do atom = 1,params%natoms
      idx = 1+(atom-1)*3 
      Do mode = 1,params%nmodes2perturb
        newref_struc(atom) = newref_struc(atom) + vecs(idx:idx+2,mode)
      End Do
    End Do

    !** Insert this new structure into the generalized coordinates
    Call rigidcoords_changeref(species(spc)%gcoords(molec)%rigidcoords,&
        newref_struc)

    !** Build the cartesian coordinates from the generalized coordinates
    Call rigidcoords_toxyz(species(spc)%gcoords(molec)%rigidcoords,&
        species(spc)%coords(1:params%natoms,molec)%rp)

    !** Apply periodic boundary conditions before forcefield evaluation
    Call simcell_pbc(simcell, &
        species(spc)%coords(1:params%natoms,molec)%rp, &
        species(spc)%coords(1:params%natoms,molec)%r, &
        species(spc)%coords(1:params%natoms,molec)%cr)

    !** Get the new energy, hardcoded to fast, don't skip intra
    succ_flag = subinteract_int(subint,species,simcell,fast, &
      .True.,skip_intra_flag,(/params%max_nrg/),(/spc,molec,0/))

  End Subroutine emperturb_move

  !---------------------------------------------------------------------------
  ! Display the move information for the Eigenmode perturbation move
  ! Requires:  params -- eigenmode perturbation parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine emperturb_display(params,indent,unit)
    Type(EMPerturb_Params), Intent(In)       :: params
    Integer, Intent(In)                      :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2

    blank = Repeat(' ',indent)

    Write(unit,'(3a)') blank,'Eigenmode perturbation move information for ', &
        Trim(molecules_name(params%spc))
    Call readhess_display(params%hessian_specs,indent+2,6)
    string1 = int2str(params%nmodes2perturb)
    string2 = real2str(params%max_disp,5)
    Write(unit,'(6a)') blank,'perturbing ', &
        Trim(string1),' eigenmode(s) with a max displacement of ', &
        Trim(string2),' Angstroms'
    Do i = 1,params%nmodes2perturb
      string1 = int2str(params%emodes(i))
      string2 = real2str(params%evalues(i),6)
      Write(unit,'(6a)') blank,'  mode ',Trim(string1), &
          ': eigenvalue = ', Trim(string2),' kJ/molAng^2'
    End Do

  End Subroutine emperturb_display

  !---------------------------------------------------------------------------
  ! Cleans the move parameter structure
  ! Requires:  params -- eigenmode perturbation parameters
  !---------------------------------------------------------------------------
  Subroutine emperturb_clean(params)
    Type(EMPerturb_Params), Intent(InOut)       :: params

    Integer            :: error

    Deallocate(params%ref_struc, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'ref_struc')
    Call readhess_clean(params%hessian_specs)
    Deallocate(params%emodes, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'emodes')
    Deallocate(params%evalues, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'evalues')
    Deallocate(params%evectors, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'evectors')

  End Subroutine emperturb_clean

End Module emperturb
