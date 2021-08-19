!-------------------------------------------------------------------------
! This module perturbs pre-defined internal rotation angles of a molecule.
! It can be used, for example, to perturb the rotation of a methyl group 
! or any other group in a molecule.  It functions by making rotations on 
! the rigid molecule reference structure.
!
! The module's operation is organized to support rotation perturbations
! with the axis defined from:
! 1) a bond, using: 
!      BOND@[bond_atom1,bond_atom2]  [rotset_atoms] 
! 2) perpendicular to a plane formed by 3-atom angle, using: 
!      ANGLE@[atom1,atom2,atom3]  [rotset_atoms] 
! 3) a torsion angle (somehow, not yet determined)
!      TORSION@[atom1,atom2,atom3,atom4]  [rotset_atoms] (not yet supported)
!------------------------------------------------------------------------

Module irotperturb

  Use defaults, Only: RDbl, strLen, lstrLen, MAX_ATOMS, twoPi
  Use utils, Only:  filesrchstr, split, toupper, allocErrDisplay, real2str, &
      int2str, deallocErrDisplay, str2seq, stripcmnt, toint
  Use file, Only: file_open
  Use random, Only: random_gaussian
  Use vector, Only: VecType, Assignment(=), Operator(+), &
      Operator(-),vector_display,vector_getnorm,vector_crossprod
  Use matrix, Only: MatrixType, Assignment(=), Operator(*),matrix_transpose, &
      matrix_display
  Use config, Only: AtMolCoords
  Use simcell, Only: SimCell_Params,simcell_pbc
  Use molecules, Only: molecules_getcomdefn,molecules_getlinewtag, &
      molecules_getfilename,molecules_getnatoms,molecules_getcoords, &
      molecules_origxform,molecules_atomsymbol,molecules_name
  Use subinteract, Only: Subset_Interactions,subinteract_int
  Use random, Only: rranf
  Use rigidcoords, Only: rigidcoords_changeref,rigidcoords_toxyz, &
      rigidcoords_getref
  Use visxyz, Only: XYZ_Entry,visxyz_make

  Implicit None
  Save

  Private
  Public :: IRotPerturb_Params,IRotSet,irotperturb_init,irotperturb_move, &
      irotperturb_adjust,irotperturb_dyninfo, &
      irotperturb_display,irotperturb_clean

  Real(kind=RDbl), Parameter        :: MIN_SIGMA = 1.0e-2_RDbl
  Real(kind=RDbl), Parameter        :: MAX_SIGMA = 10.0_RDbl

  !** Contains information for a group of atoms to be rotated
  Type IRotSet
    Integer                         :: natoms,ndefnatms
    Character(len=strLen)           :: axisdefn
    Character(len=lstrLen)          :: atomseq
    Integer, Dimension(:), Pointer  :: atoms
    Integer, Dimension(4)           :: defnatm
  End Type IRotSet

  !** Contains internal rotation information for a whole species-type
  Type IRotPerturb_Params
    Logical                              :: scale_sigma
    Integer                              :: spc,nrotsets,natoms
    Real(kind=RDbl)                      :: sigma,max_nrg
    Type(IRotSet), Dimension(:), Pointer :: set
  End Type IRotPerturb_Params

Contains
  !---------------------------------------------------------------------------
  ! Initializes Internal Rotation Perturbation move for the given species type
  ! Requires: params -- internal rotation perturb parameters
  !           spc -- species number
  !           filename -- control filename with initialization data ready
  !---------------------------------------------------------------------------
  Subroutine irotperturb_init(params,spc,filename)
    Type(IRotPerturb_Params), Intent(InOut)     :: params
    Integer, Intent(In)                         :: spc
    Character(*), Intent(In)                    :: filename            

    Integer                                  :: unit,ios,error,natoms,i,j
    Integer                                  :: nfields,nchunks,nbits
    Integer                                  :: molunit, scale_sigma
    Character(len=lstrLen)                   :: line,tag,xyzfilename
    Character(len=strLen), Dimension(strLen) :: fields,chunks,bits
    Integer, Dimension(MAX_ATOMS*3)          :: nums

    !** Set easy parameters
    params%spc = spc
    natoms = molecules_getnatoms(spc)
    params%natoms = natoms
    params%max_nrg = 10000.0_RDbl  !** maximum energy HACK for now

    !** Get the starting sigma 
    unit = file_open(filename,110)
    Read(unit,*,IOSTAT=ios) params%sigma, scale_sigma
    If (ios /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Expected to read: "gaussian_sigma, [1 or 0]" from control file'
      Stop
    End If

    !** Set the flag to decide whether we scale the jump lengths or not
    params%scale_sigma = .False.
    If (scale_sigma == 1) Then
      params%scale_sigma = .True.
    End If

    !** Get the parameters from the molecule file
    tag = 'RotationSet_Info:'
    If (.Not. molecules_getlinewtag(spc,tag,line,molunit)) Then
      Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
          ' Could not find "',Trim(tag),'" tag in molecule file: ', &
          Trim(molecules_getfilename(spc))
      Stop
    End If

    !** Get the number of rotation sets
    Read(molunit,'(a)') line
    nfields = split(line,fields)    
    params%nrotsets = toint(fields(1))
    Allocate(params%set(params%nrotsets), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)      

    !** Get the rotation sets, format should be one of the following:
    !** BOND@[bond_atom1,bond_atom2]  [rotset_atoms] 
    !** ANGLE@[atom1,atom2,atom3]  [rotset_atoms] 
    !** TORSION@[atom1,atom2,atom3,atom4]  [rotset_atoms] (not yet supported)
    !** the rotation set of atoms is specified as a sequence e.g. 1,2-6,8,9
    Do i = 1,params%nrotsets
      Read(molunit,'(a)') line
      nfields = split(line,fields)    

      !** Get the rotation axis definition
      nchunks = split(fields(1),chunks,'@')    
      params%set(i)%axisdefn = ToUpper(chunks(1))
      
      !** Interpret the keyword defining the rotation axis
      Select Case(params%set(i)%axisdefn)
      Case('BOND')
        params%set(i)%ndefnatms = 2
      Case('ANGLE')
        params%set(i)%ndefnatms = 3
      Case('TORSION')
        params%set(i)%ndefnatms = 4
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' TORSION option for internal rotations not yet supported'
        Stop
      Case Default
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' Could not understand internal rotation axis definition: ', &
            Trim(params%set(i)%axisdefn)
        Stop
      End Select

      !** Get the rotation axis definition
      nbits = split(chunks(2),bits,',')    
      If (nbits /= params%set(i)%ndefnatms) Then
        Write(0,'(2a,i4,a,2i3)') __FILE__,": ",__LINE__, &
            ' Number of atoms in irot definition does not match that expected', &
            nbits,params%set(i)%ndefnatms
        Stop
      End If
      Do j = 1,params%set(i)%ndefnatms
        params%set(i)%defnatm(j) = toint(bits(j), &
            'error reading internal rotation axis definition')
      End Do
      
      !** Get the atom number sequence to specify group to rotate
      params%set(i)%atomseq = fields(2)
      params%set(i)%natoms = str2seq(fields(2),nums)

      !** store the atom number sequence
      Allocate(params%set(i)%atoms(params%set(i)%natoms), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)      
      params%set(i)%atoms = nums(1:params%set(i)%natoms)
    End Do

  End Subroutine irotperturb_init

  !---------------------------------------------------------------------------
  ! Performs all the internal rotation transformations on a single molecule. 
  ! Also calculates initial and final energies.
  ! Requires:  params -- internal rotation perturbation parameters
  !            spc -- species number
  !            molec -- molecule number
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias factor associated with move
  !            succ_flag -- returns true if move and energy evaluation went ok
  !            subint -- subset interactions (molecule-system)
  !---------------------------------------------------------------------------
  Subroutine irotperturb_move(params,spc,molec,species,simcell, &
      biasfactor,succ_flag,subint)
    Type(IRotPerturb_Params), Intent(In)           :: params
    Integer, Intent(In)                            :: spc,molec
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In)               :: simcell
    Real(kind=RDbl), Intent(Out)                   :: biasfactor
    Logical, Intent(Out)                           :: succ_flag
    Type(Subset_Interactions), Intent(InOut)       :: subint

    Integer                                  :: i,a
    Logical                                  :: fast,skip_intra_flag
    Type(VecType), Dimension(params%natoms)  :: newref_struc

    succ_flag = .False.
    biasfactor = 1.0_RDbl  !** currently, no biasing

    !** Get the old energies
    fast = .True.
    skip_intra_flag = .False.
    succ_flag = subinteract_int(subint,species,simcell,fast, &
      .False.,skip_intra_flag,(/params%max_nrg/),(/spc,molec,0/))
    If (.Not. succ_flag) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "Problem with nrg calculation of an existing structure"
      Stop
    End If

    !** Copy the reference structure to local variable
    Call rigidcoords_getref(species(spc)%gcoords(molec)%rigidcoords,newref_struc)

    !** Perform the rotation(s) to get the new reference structure
    Do i = 1,params%nrotsets
      Call irotperturb_irot(params%set(i),params%sigma,newref_struc)
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

    !** get the new energy, hardcoded to fast, don't skip intra
    succ_flag = subinteract_int(subint,species,simcell,fast, &
      .True.,skip_intra_flag,(/params%max_nrg/),(/spc,molec,0/))

  End Subroutine irotperturb_move

  !---------------------------------------------------------------------------
  ! Performs the internal rotation transformation on one of the rotation sets
  ! in a single molecule.
  ! Requires:  setparams -- parameters for a single internal rotation set
  !            sigma -- sigma for gaussian distribution (move size)
  !            coords -- coordinates of entire molecule, indexed by atom num
  !---------------------------------------------------------------------------
  Subroutine irotperturb_irot(setparams,sigma,coords)
    Type(IRotSet), Intent(In)                  :: setparams
    Real(kind=RDbl), Intent(In)                :: sigma
    Type(VecType), Dimension(:), Intent(InOut) :: coords

    Integer                    :: i,a,atom
    Real(kind=RDbl)            :: rotangle
    Type(VecType)              :: axisvec,basecoord
    Type(MatrixType)           :: rotmtx

    !** Get the rotation axis
    axisvec = irotperturb_getaxis(setparams,coords,basecoord)

    !** Pick a random number from a Gaussian distribution centered at zero
    rotangle = random_gaussian(0.0_RDbl,sigma)
    rotangle = Mod(rotangle,twoPi)

    !** Get the rotation matrix
    Call irotperturb_buildmtx(axisvec,rotangle,rotmtx)

    !** Translate the coordinates back to some base and rotate
    Do a = 1,setparams%natoms
      atom = setparams%atoms(a)
      coords(atom) = coords(atom) - basecoord
      coords(atom) = rotmtx*coords(atom)
      coords(atom) = coords(atom) + basecoord
    End Do

  End Subroutine irotperturb_irot

  !---------------------------------------------------------------------------
  ! Calculates the rotation axis given the definition and the coordinates
  ! Requires:  setparams -- parameters for a single internal rotation set
  !            coords -- coordinates of entire molecule, indexed by atom num
  !            base -- position vector forming 'base' of rotation axis
  !---------------------------------------------------------------------------
  Type(VecType) Function irotperturb_getaxis(setparams,coords,base)
    Type(IRotSet), Intent(In)                  :: setparams
    Type(VecType), Dimension(:), Intent(InOut) :: coords
    Type(VecType), Intent(Out)                 :: base

    Type(VecType)       :: bvec1,bvec2

    !** Construct the vector defining the rotation axis
    Select Case(setparams%axisdefn)
    Case('BOND')
      irotperturb_getaxis = coords(setparams%defnatm(2)) - &
          coords(setparams%defnatm(1))
      base = coords(setparams%defnatm(1))

    Case('ANGLE')
      bvec1 = coords(setparams%defnatm(1)) - coords(setparams%defnatm(2))
      bvec2 = coords(setparams%defnatm(3)) - coords(setparams%defnatm(2))
      irotperturb_getaxis = vector_crossprod(bvec1,bvec2)
      base = coords(setparams%defnatm(2))

    Case('TORSION')
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' TORSION option not yet supported'
      Stop
    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Could not understand internal rotation axis definition: ', &
          Trim(setparams%axisdefn)
      Stop
    End Select

  End Function irotperturb_getaxis

  !---------------------------------------------------------------------------
  ! Builds the matrix needed to perform the rotation about an axis.  Too bad
  ! I don't remember how it works, my notes are 4000 miles away.  Probably
  ! it's a composition of three matrices, using similarity transform to rotate
  ! a standard rotation matrix to some arbitary basis.  One rotation transform
  ! rotates the coordinate system to coincide with the bond axis and another
  ! rotates around the bond axis.
  ! Requires:  rotaxis -- the vector defining the rotation axis
  !            rotangle -- rotation angle in radians
  !            rotmtx -- returned rotation matrix
  !---------------------------------------------------------------------------
  Subroutine irotperturb_buildmtx(rotaxis,rotangle,rotmtx)
    Type(VecType), Intent(In)      :: rotaxis
    Real(kind=RDbl), Intent(In)    :: rotangle
    Type(MatrixType), Intent(Out)  :: rotmtx

    Real(kind=RDbl)      :: length,xylength,singo22
    Real(kind=RDbl)      :: cosa,cosa2,cos2a,sina,sina2,sin2a
    Real(kind=RDbl)      :: cosb,cosb2,cos2b,sinb,sinb2,sin2b
    Real(kind=RDbl)      :: cosg,cosg2,cos2g,sing,sing2,sin2g

    length = vector_getnorm(rotaxis)
    xylength = Sqrt(rotaxis%comp(1)*rotaxis%comp(1) + &
        rotaxis%comp(2)*rotaxis%comp(2))
    
    !** ALPHA portion
    cosa = rotaxis%comp(1)/xylength
    sina = rotaxis%comp(2)/xylength
    cosa2 = cosa*cosa
    sina2 = sina*sina
    cos2a = 1.0e0 - 2*sina2
    sin2a = 2*sina*cosa
    
    !** BETA portion
    cosb = xylength/length
    sinb = -rotaxis%comp(3)/length
    cosb2 = cosb*cosb
    sinb2 = sinb*sinb
    cos2b = 1.0e0 - 2*sinb2
    sin2b = 2*sinb*cosb
    
    !** GAMMA portion
    cosg = Cos(rotangle)
    sing = Sin(rotangle)
    cosg2 = cosg*cosg
    sing2 = sing*sing
    cos2g = 1.0e0 - 2*sing2
    sin2g = 2*sing*cosg
    
    singo22 = Sin(rotangle/2.0e0)   !**sin^2(gamma/2)
    singo22 = singo22*singo22

    rotmtx%comp(1,1) = cosg*sina2 + cosa2*(cosb2 + cosg*sinb2)
    rotmtx%comp(1,2) = cosb2*sin2a*singo22 - sinb*sing
    rotmtx%comp(1,3) = cosb*(cosa*(-1.0e0 + cosg)*sinb - sina*sing)
    
    rotmtx%comp(2,1) = cosb2*sin2a*singo22 + sinb*sing
    rotmtx%comp(2,2) = cosa2*cosg + sina2*(cosb2 + cosg*sinb2)
    rotmtx%comp(2,3) = cosb*((-1.0e0 + cosg)*sina*sinb + cosa*sing)
    
    rotmtx%comp(3,1) = cosb*(cosa*(-1.0e0 + cosg)*sinb + sina*sing)
    rotmtx%comp(3,2) = cosb*((-1.0e0 + cosg)*sina*sinb - cosa*sing)
    rotmtx%comp(3,3) = cosb2*cosg + sinb2

  End Subroutine irotperturb_buildmtx

  !-------------------------------------------------------------------------
  ! Adjust the sigma of the Gaussian jump size distribution to try to 
  ! keep the acceptance ratio close to 50%
  ! Requires:  params -- internal rotation perturbation parameters
  !            ratio -- current acceptance ratio
  !-------------------------------------------------------------------------
  Subroutine irotperturb_adjust(params,ratio)
    Type(IRotPerturb_Params), Intent(InOut)   :: params
    Real(kind=RDbl), Intent(In)               :: ratio

    If (.Not. params%scale_sigma) Return

    !** These values should not be hard coded here, should go into control file
    If (ratio < 0.49) Then
      params%sigma = Max(params%sigma*0.95_RDbl, MIN_SIGMA)
    Else If (ratio > 0.51) Then
      params%sigma = Min(params%sigma*1.05_RDbl, MAX_SIGMA)
    End If

  End Subroutine irotperturb_adjust 

  !----------------------------------------------------------------------
  ! Single line display information for dynamic Translation parameters
  ! Requires:  params -- internal rotation perturbation parameters
  !----------------------------------------------------------------------
  Function irotperturb_dyninfo(params)
    Character(len=strLen)                  :: irotperturb_dyninfo
    Type(IRotPerturb_Params), Intent(In)   :: params

    Write(irotperturb_dyninfo,'(a,f6.3)') 'Disp. Sigma: ',&
        params%sigma

  End Function irotperturb_dyninfo

  !---------------------------------------------------------------------------
  ! Display the move information from the internal rotation parameters
  ! Requires:  params -- irot parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine irotperturb_display(params,indent,unit)
    Type(IRotPerturb_Params), Intent(In)       :: params
    Integer, Intent(In)                        :: indent,unit

    Integer                     :: i,j
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2

    blank = Repeat(' ',indent)

    Write(unit,'(4a)') blank,'Internal rotation perturbation move ', &
        'information for ',Trim(molecules_name(params%spc))
    Write(unit,'(2a,f8.3)') blank,'sigma for Gaussian distribution: ', &
        params%sigma
    Do i = 1,params%nrotsets
      string1 = int2str(i)
      string2 = int2str(params%set(i)%natoms)
      Write(unit,'(6a)') blank,'Rotation Set: ',Trim(string1),' with ', &
          Trim(string2),' atoms'
      Write(unit,'(a,2x,2a)') blank,'Axis definition type: ', &
          Trim(params%set(i)%axisdefn)
      Write(unit,'(a,2x,2a)') blank,'atoms in rotation group: ', &
          Trim(params%set(i)%atomseq)
    End Do

  End Subroutine irotperturb_display

  !---------------------------------------------------------------------------
  ! Cleans the move parameter structure
  ! Requires:  params -- internal rotation perturbation parameters
  !---------------------------------------------------------------------------
  Subroutine irotperturb_clean(params)
    Type(IRotPerturb_Params), Intent(InOut)       :: params

    Integer            :: i,error

    Do i = 1,params%nrotsets
      Deallocate(params%set(i)%atoms, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'set atoms')
    End Do

    Deallocate(params%set, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'emodes')

  End Subroutine irotperturb_clean


End Module irotperturb
