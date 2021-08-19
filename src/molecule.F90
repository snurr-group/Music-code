!-----------------------------------------------------------------------
! This module contains the data type that holds individual molecule 
! information.  Molecules are entities that are made up of a collection 
! of one or more bonded atoms.  Molecule specific information includes:
!    reference coordinates (from molecule file)
!    mass, charge, stereochemistry
!    local coordinate system 
!    principle moments of inertia
!    nearest neighbor bonded connections for each ATOM
!    nearest neighbor bonded connections for each atom TYPE
!    energy statistics based on sorbate type (should be in molecules?)
!-----------------------------------------------------------------------

Module molecule

  Use defaults, Only: RDbl, strLen, zeroTolerance, lstrLen, MAX_ATOMS, one, &
      zero, Echarge
  Use utils, Only: split, stripcmnt, findint, tolower, getpath, filesrchstr, &
      toint, isfileopen, toupper, toreal, filesrchwocomment, real2str, &
      allocErrDisplay, int2str, deAllocErrDisplay
  Use file, Only: file_getunit
  Use vector, Only: VecType, mag, vector_getnormsq, vector_display, &
      vector_getcomp, vector_iscollinear, vector_crossprod,&
      Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use matrix, Only: MatrixType, matrix_getinv, Assignment(=), Operator(*), &
      matrix_transpose, matrix_display, matrix_identity
  Use atom, Only: atom_gettypename, atom_getmass, atom_nbonds, &
      atom_getntypes, atom_getsymbol, atom_getname, atom_invmass, atom_isIon
  Use localcoords, Only: LocalCoordInfo,localcoords_init, &
      localcoords_display
  Use fundcell, Only: Fundamental_Cell
  Use stats, Only: Statistics, stats_init, stats_update
  Use connects, Only: MAX_CONNECTIONS, &
      AtomConnections, connects_nconnected, connects_display, &
      connects_initarray, connects_initmaxcon, connects_connectatoms, &
      connects_generateAll, ConnectMatrix, &
      connects_initconnection, connects_displayall, &
      connects_regenerateAll, connects_atomsxaway
  Use visxyz, Only: XYZ_Entry, visxyz_dump

  Implicit None
  Save

  Private
  Public :: AtomInfo, MolecularParams, molecule_init, &
      molecule_getname, molecule_ischarged, molecule_getfilename, &
      molecule_getbodyaxes, molecule_getnatoms, molecule_getcharge, &
      molecule_getnthatom, molecule_getdof, molecule_getatomcoords, &
      molecule_getnatomtypes, molecule_getmass, molecule_AtomMass, &
      molecule_changedof, molecule_getmoi, molecule_getcomdefn, &
      molecule_getgcmodeltype, molecule_calcdipole, &
      molecule_display, molecule_connectinit, molecule_netcharge, &
      molecule_getnthstereo, molecule_expand, molecule_molec2xyz, &
      molecule_reconnect, molecule_getneighbors, molecule_origxform, &
      molecule_sampleCF, molecule_isCoreShell, molecule_getrefcom, &
      molecule_getcom, molecule_checkconnect, molecule_getconnects

  External tred2, tql2
 
  !** NAG strLen workaround
  Integer :: dummyMolecule = strLen
  
  Integer, Parameter    :: MAX_SETATOMS = 50
 
  Type AtomInfo
    Character(len=strLen)        :: name
    Integer                      :: atype
    Type(VecType)                :: r
    Real(kind=RDbl)              :: q0       ! Charge
    Integer                      :: set,type
    Character                    :: stereo   ! stereochemistry
                                             ! 'S','R', or 'N' for none
  End Type AtomInfo
  
  Type MolecularParams
    Character(len=strLen)                :: name
    Character(len=strLen)                :: file  
    Character(len=strLen)                :: geometry
    Logical                              :: zeolite,charged
    Real(kind=RDbl)                      :: mass
    Real(kind=RDbl),Dimension(:),Pointer :: atom_mass_list
    Integer                              :: natoms, natypes

    !** list of types : atypes(natypes), 
    !** list of types of all atoms : allatypes(natoms)
    Integer, Dimension(:), Pointer       :: atypes, allatypes

    Integer                              :: dof     !Degrees of freedom 
    Character(len=strLen)                :: dof_origin
    Type(AtomInfo),Dimension(:), Pointer :: atoms    

    !** The name of the coordinate system in which the atom coords are stored
    Character(len=strLen) :: coordsys

    !** The TYPE of the generalized coordinates for this molecule
    Character(len=strLen)                :: gcmodeltype

    !** The axes of the body or local coordinates system
    Type(LocalCoordInfo)          :: bodyaxes
    Real(kind=RDbl), Dimension(3) :: pMoI  ! Principle Moments of Inertia

    !** information on the transformation that took the molecule file
    !** coordinates to the stored reference configuration
    Type(VecType)     :: trans_xform    ! molfile + trans_xform -> ref (1st)
    Type(MatrixType)  :: rot_xform      ! rot_xform*molfile -> ref (2nd)
    Logical           :: rotateflag, translateflag

    !** Connection information
    !** connect stores the connection information for each ATOM in the molecule
    Type(AtomConnections), Dimension(:), Pointer :: connect

    !** connections stores the connection information for atom-pair TYPES
    Type(ConnectMatrix) :: connections
    
    !** Stores the statistics for each sorbate type
    Type(Statistics)   :: tmole      ! Temp. of each species
    Type(Statistics)   :: tatom                
    Type(Statistics)   :: u_sz       ! Sorbate-zeolite nrg
  End Type MolecularParams
  
  Interface getnthatom
    Module procedure molecule_getnthatom
  End Interface 

Contains


  !----------------------------------------------------------------------------
  ! FUNCTION getnatomtypes
  ! REQUIRES molecule (type MolecularParams), fully initialized
  ! OPTIONAL atypes (Integer array or arbitrary size), allTypes (logical)
  ! DESCRIPTION
  ! Gets the number of atom TYPES in a sorbate. Optionally returns an array
  ! 'atypes' with the atom types in it. If allTypes is true, it returns the
  ! atom type for each atom. Otherwise, only the overall types are returned.
  !----------------------------------------------------------------------------
  Integer Function molecule_getnatomtypes(molecule,atypes,allTypes)

    Type(MolecularParams) :: molecule
    Integer, Dimension(:), Optional :: atypes
    Logical, Optional :: allTypes

    Logical :: all

    molecule_getnatomtypes=molecule%natypes
    If (Present(atypes)) Then
      all=.False.
      If (Present(allTypes)) all=allTypes
      If (all) Then
        atypes(1:molecule%natoms)=molecule%allatypes
      Else
        atypes(1:molecule%natypes)=molecule%atypes
      Endif
    Endif

  End Function molecule_getnatomtypes
 


  !----------------------------------------------------------------------------
  ! FUNCTION getnatypesInit
  ! REQUIRES molecule (type MolecularParams)
  ! OPTIONAL atypes (Integer array or arbitrary size), allTypes (logical)
  ! DESCRIPTION
  ! Gets the number of atom TYPES in a sorbate. Optionally returns an array
  ! 'atypes' with the atom types in it. If allTypes is true, it returns the
  ! atom type for each atom. Otherwise, only the overall types are returned.
  !**********************************************************************
  ! This is reasonably expensive, while dealing with a large molecule
  ! like a zeolite, use molecules_getnatomtypes forfaster access
  ! this routine should be called only during _init
  !***********************************************************************
  !----------------------------------------------------------------------------
  Integer Function molecule_getnatypesInit(molecule,atypes,allTypes)
 
    Type(MolecularParams) :: molecule
    Integer :: i,atomt1,k,j
    Integer, Dimension(:), Optional :: atypes
    Logical, Optional :: allTypes
 
    Integer, Dimension(molecule%natoms) :: alist
    Logical :: all

    If (Present(allTypes)) Then
      all = allTypes
    Else
      all = .False.
    End If

    k = 0
    alist = 0
 
    Do i = 1,molecule%natoms
      atomt1 = molecule%atoms(i)%atype
      If (Present(atypes).And.all) atypes(i) = atomt1
      j = findint(alist,atomt1)
      If (j == 0) Then
        k = k+1
        alist(k) = atomt1
        If (Present(atypes).And..Not.all) Then
          atypes(k) = atomt1
        End If
      End If
    End Do
 
    molecule_getnatypesInit = k
    
  End Function molecule_getnatypesInit
 
  !--------------------------------------------
  ! Gets the number of atoms in one sorbate
  !--------------------------------------------
  Integer Function molecule_getnatoms(molecule)
    Type(MolecularParams) :: molecule
 
    molecule_getnatoms = molecule%natoms
  End Function molecule_getnatoms

  !----------------------------------------------------------------------------
  ! Gets the charge atom atom "atomno" of molecule "molecule"
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function molecule_getcharge(molecule, atomno)
    Type(MolecularParams) :: molecule
    Integer, Intent(in)   :: atomno
    
    molecule_getcharge = molecule%atoms(atomno)%q0
  End Function molecule_getcharge

  !----------------------------------------------------------------------------
  ! Gets the charge atom atom "atomno" of molecule "molecule"
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function molecule_calcdipole(molecule)
    Type(MolecularParams) :: molecule
    Integer   :: a
    Real(kind=RDbl) :: dipole
    Type(VecType) :: dvec
    dvec=zero
    Do a=1, molecule%natoms
     dvec=dvec+(molecule%atoms(a)%r*molecule%atoms(a)%q0)
    End Do 

    dipole=mag(dvec)
    ! convert to coulomb-meter
    dipole=dipole*Echarge*1.0e-10

    ! convert to debye
    dipole=dipole/3.33e-30
    molecule_calcdipole=dipole
  End Function molecule_calcdipole
 
  !-----------------------------------------------
  ! Get the molecule name
  !-----------------------------------------------
!!$  Character(len=strLen) Function molecule_getname(molecule)
  Function molecule_getname(molecule)
    Character(len=strLen) :: molecule_getname
    Type(MolecularParams) :: molecule
 
    molecule_getname = molecule%name
  End Function molecule_getname

  !----------------------------------------------------------------------------
  ! Get the molecule's mass
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function molecule_getmass(molecule)
    Type(MolecularParams), Intent(In) :: molecule
    molecule_getmass = molecule%mass
  End Function molecule_getmass

  !----------------------------------------------------------
  ! This routine gets the mass of species-spc atom-a 
  !----------------------------------------------------------
  Real(kind=RDbl) Function molecule_AtomMass(molec, a)
    Type(MolecularParams),Intent(in)           :: molec
    Integer, Intent(In)                        :: a
    molecule_AtomMass=molec%atom_mass_list(a)
  End Function molecule_AtomMass

  !-----------------------------------------------
  ! Get the molecule filename
  !-----------------------------------------------
  Function molecule_getfilename(molecule)
    Character(len=strLen) :: molecule_getfilename
    Type(MolecularParams) :: molecule
 
    molecule_getfilename = molecule%file
  End Function molecule_getfilename
 
  !---------------------------------------------------------------
  ! Get the type of the generalized coordinates used
  ! to describe this molecule
  !---------------------------------------------------------------
!!$  Character(len=strLen) Function molecule_getgcmodeltype(molecule)
  Function molecule_getgcmodeltype(molecule)
    Character(len=strLen) :: molecule_getgcmodeltype
    Type(MolecularParams) :: molecule
    
    molecule_getgcmodeltype = molecule%gcmodeltype
    Return
  End Function molecule_getgcmodeltype

  !----------------------------------------------------------------------------
  ! Get the Center of Mass coordinates of the molecule
  ! The original routine was a function, but was causing compiler
  ! segmentation faults.
  !----------------------------------------------------------------------------
  Subroutine molecule_getcomdefn(molecule,com)
    Type(MolecularParams), Intent(In)          :: molecule
    Type(VecType), Dimension(:), Intent(InOut) :: com

    Integer   :: i, natoms
    
    natoms = molecule%natoms
    Do i=1, natoms
      com(i) = molecule%atoms(i)%r
    End Do
    com(1:molecule%natoms) = molecule%atoms(1:natoms)%r
  End Subroutine molecule_getcomdefn

  !-----------------------------------------------------
  ! Get the name of the molecule
  !-----------------------------------------------------
!!$  Character(len=strLen) Function molecule_name(molecule)
  Function molecule_name(molecule)
    Character(len=strLen) :: molecule_name
    Type(MolecularParams), Intent(in) :: molecule
    
    molecule_name = molecule%name
  End Function molecule_name

  !----------------------------------------------------------------------------
  ! Gets the number of degrees of freedom for the molecule
  !----------------------------------------------------------------------------
  Integer Function molecule_getdof(molecule,dof_origin)
    Type(MolecularParams), Intent(In) :: molecule
    Character(*), Optional, Intent(Out) :: dof_origin
    molecule_getdof = molecule%dof
    If (Present(dof_origin)) Then
      dof_origin = molecule%dof_origin
    End If
  End Function molecule_getdof

  !----------------------------------------------------------------------------
  ! Change the number of degrees of freedom for the molecule
  !----------------------------------------------------------------------------
  Subroutine molecule_changedof(molecule,newdof,new_origin)
    Type(MolecularParams), Intent(InOut) :: molecule
    Integer, Intent(In)                  :: newdof
    Character(*), Intent(In), Optional   :: new_origin

    molecule%dof = newdof
    If (Present(new_origin)) molecule%dof_origin = new_origin
  End Subroutine molecule_changedof

  !----------------------------------------------------------------------------
  ! Returns the type of the nth atom in molecule
  !----------------------------------------------------------------------------
  Integer Function molecule_getnthatom(molecule,atom)
    Type(MolecularParams), Intent(in) :: molecule
    Integer, Intent(In) :: atom
    molecule_getnthatom = molecule%atoms(atom)%atype
  End Function molecule_getnthatom

  !---------------------------------------------------------------------------
  ! Returns the stereochemistry of the nth atom in molecule
  !---------------------------------------------------------------------------
  Character Function molecule_getnthstereo(molecule,atom)
    Type(MolecularParams), Intent(in) :: molecule
    Integer, Intent(in)               :: atom
    molecule_getnthstereo = molecule%atoms(atom)%stereo
  End Function molecule_getnthstereo


  !----------------------------------------------------------------------------
  ! Returns the coordinates of atom number atom in molecule
  !----------------------------------------------------------------------------
  Type(VecType) Function molecule_getatomcoords(molecule,atom)
    Type(MolecularParams), Intent(in) :: molecule
    Integer, Intent(In) :: atom
    molecule_getatomcoords = molecule%atoms(atom)%r
  End Function molecule_getatomcoords

  !--------------------------------------------------------------------------
  ! Returns the information about the local coordinate system of the molecule
  !--------------------------------------------------------------------------
  Type(LocalCoordInfo) Function molecule_getbodyaxes(molecule)
    Type(MolecularParams), Intent(In) :: molecule

    molecule_getbodyaxes = molecule%bodyaxes

  End Function molecule_getbodyaxes

  !------------------------------------------------------------------
  ! Check to see if the atoms of the molecule "molecule" are charged
  !------------------------------------------------------------------
  Logical Function molecule_ischarged(molecule)
    Type(MolecularParams) :: molecule

    Integer     :: i
        
    Do i=1, molecule%natoms
      If (Abs(molecule%atoms(i)%q0) > 1.0e-5_RDbl) Then
        molecule_ischarged = .True.
        Return
      End If
    End Do
    molecule_ischarged = .False.
  End Function molecule_ischarged

  !------------------------------------------------------------------
  ! Calculates the net charge on the molecule
  !------------------------------------------------------------------
  Real(kind=RDbl) Function molecule_netcharge(molecule)
    Type(MolecularParams) :: molecule

    Integer     :: i

    molecule_netcharge = 0.0_RDbl
        
    Do i = 1,molecule%natoms
!LC      Write(*,*) i,molecule%atoms(i)%q0
      molecule_netcharge = molecule_netcharge + molecule%atoms(i)%q0
    End Do

  End Function molecule_netcharge

  !----------------------------------------------------------------------------
  ! Returns true if the atom is part of a core-shell. Optionally
  ! returns the corresponding atom of the core/shell pair.
  !----------------------------------------------------------------------------
  Logical Function molecule_isCoreShell(molecule,atom,isCore,partner)
    Type(MolecularParams), Intent(In) :: molecule
    Integer, Intent(In) :: atom
    Logical, Intent(Out):: isCore
    Integer, Intent(Out), Optional :: partner

    Integer, Dimension(molecule%natoms) :: neighbors
    Integer :: nneighbors, i 
    Logical :: isPartnerCore, isPartnerIon

    !** Default value
    molecule_isCoreShell = .False.
    isCore = .False.

    !** Call atom_isIon to check the atom's status
    molecule_isCoreShell = atom_isIon(molecule%atoms(atom)%atype,isCore)

    !** If partner isn't here, just return
    If (.Not.Present(partner)) Return

    !** If the molecule isn't a core or shell, return
    If (.Not.molecule_isCoreShell) Return

    !** Get the atom connections
    nneighbors = molecule_getneighbors(molecule,atom,neighbors)

    !** Loop through the neighbors and find the core or shell
    Do i = 1, nneighbors
      isPartnerIon = atom_isIon(molecule%atoms(neighbors(i))%atype, &
          isPartnerCore)
      !** Check if the partner is an ion shell. If so, accept and return.
      If (isCore.And.isPartnerIon.And..Not.isPartnerCore) Then
        partner = neighbors(i)
        Exit
        !** Check if the partner is an ion core. If so, accept and return.
      Else If (.Not.isCore.And.isPartnerIon.And.isPartnerCore) Then
        partner = neighbors(i)
        Exit
      End If
    End Do

  End Function molecule_isCoreShell

  !---------------------------------------------------------------
  ! Initialize the molecule information
  ! Requires: molecule -- the molecule's data structure
  !           molecName -- molecule name
  !           filename -- file to take definition from
  !---------------------------------------------------------------
  Subroutine molecule_init(molecule,molecName,filename)
    Type(MolecularParams), Intent(InOut)     :: molecule
    Character(*), Intent(In)                 :: molecName,filename

    Integer                                  :: i,j,ios,natoms, natypes
    Integer                                  :: unit,lineno,nfields,error
    Real(kind=RDbl) :: avg_charge
    Logical                                  :: isCoreShell, isCore
    Character(len=strLen)                    :: string
    Character(len=255)                       :: line,srchstr,filewpath,key
    Character(len=255)                       :: stripLine
    Character(len=strLen), Dimension(strLen) :: fields

    !** get the filename and unit
    molecule%file = filename
    molecule%zeolite = .False.
    molecule%charged = .False.
    unit = file_getunit(filename)
    Write(filewpath,'(2a)') Trim(getpath('MOLSDIR')),Trim(filename)

    !** NOTE: it's important that the output shows where the molecule came from
    Write(*,'(4x,4a)') 'Reading molecule information for ',Trim(MolecName), &
        ' from: ',Trim(filewpath)

    !** Open the molecule file for reading
    Open(unit=unit, file=trim(filewpath), status='old', IOSTAT=ios)
    If (ios /= 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          '  Could not open file: ',Trim(filename)
      Stop
    End If

    !** find the molecule in the file
    srchstr = 'Molecule_Name:'
    lineno = filesrchstr(unit, srchstr, line)
    If (lineno == 0) Then
      Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
          ' Could not find search string: ',Trim(srchstr), &
          ' in file ', Trim(filewpath)
      Stop
    End If

    nfields = split(line,fields)
    molecule%name = fields(2)
    If (molecName /= molecule%name) Then
      Write(0,'(2a,i5,3a)') __FILE__,":",__LINE__, &
          " Requested molecule to initialize (",Trim(molecName),")"
      Write(0,'(2x,3a)') "does not match name in the molecule file (", &
          Trim(molecule%name),")"
      Stop
    End If

    !** Interpret other fields if they're there
    If (nfields > 2) Then
      If (fields(3) == 'CHARGED') Then
        molecule%charged = .True.
      End If
    End If

    !** Initialize the atom coordinates and mass
    Call molecule_coordinit(molecule,unit)

    !** Set stereochemistries
    Call molecule_stereoinit(molecule,unit)

    !** Build connections
    Nullify(molecule%connect)
    Call molecule_connectinit(molecule,unit)

    !** Default degrees of freedom
    !** If the molecule contains core-shell pairs, remove the shells
    !** since they are tied to the cores.
    natoms = 0
    Do i = 1, molecule%natoms
      isCoreShell = molecule_isCoreShell(molecule,i,isCore)
      If (isCoreShell.And..Not.isCore) Cycle
      natoms = natoms + 1
    End Do
    molecule%dof = 3*natoms
    molecule%dof_origin = '3N_Calculation'

    !** Check to see if dof is specified in the file
    key = "Molecule_DOF"
    lineno = filesrchwocomment(unit,key,line,.True.)

    If (lineno /= 0) Then

      !**Found specified dof
      stripLine = stripcmnt(line)
      j = split(stripLine,fields)
      molecule%dof = toint(fields(2),'molecule DOF read')
      molecule%dof_origin = 'Molecule_File'

    End If

    !** Check for charge neutrality
    avg_charge=(Abs(molecule_netcharge(molecule)))/natoms
    If (avg_charge > zeroTolerance) Then
      string = real2str(molecule_netcharge(molecule),4)
      Write(0,'(1x,2a,i5,5a)') __FILE__," :",__LINE__, &
          ' WARNING: ',Trim(molecule%name),' is not charge neutral (',&
          Trim(string),')'
      If (.Not. molecule%charged) Then
        Write(0,'(2x,2a)') '"CHARGED" must be used on the Molecule_Name line ',&
            'if charged molecules are to be employed'
        Stop
      End If
    End If

    Close(unit)
    !    Call molecule_molec2xyz(molecule,'test.xyz','just a test')

    ! initialize some arrays that are stored for faster access
    ! this avoids calls to expensive molecule_getnatomtypes() frequently
    natypes =  molecule_getnatypesInit(molecule)

    molecule%natypes=natypes
    Allocate(molecule%atypes(natypes), molecule%allatypes(natoms), &
        STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)


    natypes =  molecule_getnatypesInit(molecule,molecule%atypes, .False. )
    natypes =  molecule_getnatypesInit(molecule,molecule%allatypes, .True. )

    !** fill in atom-mass-list array
    Allocate(molecule%atom_mass_list(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1, natoms
      molecule%atom_mass_list(i)=atom_getmass(molecule%allatypes(i))
    End Do    

  End Subroutine molecule_init

  !---------------------------------------------------------------
  ! Nullify the pointers in the molecule data type
  ! Requires: molecule -- the molecule's data structure
  !---------------------------------------------------------------
  Subroutine molecule_null(molecule)
    Type(MolecularParams), Intent(InOut)     :: molecule

    Nullify(molecule%atoms)
    Nullify(molecule%connect)

  End Subroutine molecule_null

  !---------------------------------------------------------------------------
  ! Expands the molecule by the offset. This produces n EXACT copies of the
  ! molecule and creates one larger molecule. It can only be used currently
  ! to expand a molecule whose COM coordinates are in cartesian coords.
  ! Requires:  spc -- species number
  !            offset -- direction to offset new copies
  !            n -- number of copies
  !            reconnect -- flag if the connects should be remade
  !            fcell -- optional fundamental cell if PBCs are to be used
  !---------------------------------------------------------------------------
  Subroutine molecule_expand(molecule,offset,n,reconnect,fcell)
    Type(MolecularParams), Intent(InOut)           :: molecule
    Type(VecType), Intent(In)                      :: offset
    Integer, Intent(In)                            :: n
    Logical, Intent(In)                            :: reconnect   
    Type(Fundamental_Cell), Intent(In), Optional   :: fcell    

    Integer               :: i, error, a, newa, newatoms, oldatoms
    Type(VecType)         :: moveby
    Type(MolecularParams) :: old

    !** If n = 0, then just return
    If (n == 0) Return

    !** First, we need to make a copy of the original molecule
    oldatoms=molecule%natoms
    old%natoms = molecule%natoms
    Call molecule_null(old)
    Call molecule_copyAtoms(molecule,old)

    Allocate(old%allatypes(oldatoms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"atomtypes")
    old%allaTypes(1:oldatoms)=molecule%allaTypes(1:oldatoms)


    !** Deallocate the old molecule atoms, atomtypes arrays
    Deallocate(molecule%atoms,molecule%allatypes,stat=error)
    If (error/=0) Call deAllocErrDisplay(__FILE__,__LINE__,"atoms,atomtypes")

    !** Set the new number of atoms
    newatoms = oldatoms*(n+1)

    !** Allocate space for the new atoms in the original molecule
    Allocate(molecule%atoms(newatoms),molecule%allatypes(newatoms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"atoms,atomtypes")

    !** Copy the original atoms back into the newly expanded atom list
    Call molecule_copyAtoms(old,molecule)
    molecule%allaTypes(1:oldatoms)=old%allaTypes(1:oldatoms)

    !** Set the new number of atoms
    molecule%natoms = newatoms

    !** Set the new molecule mass
    molecule%mass = molecule%mass*(n+1)

    !** Loop through the number of copies to make
    Do i = 2, n+1
      !** This is the amount we are offsetting the new atoms by
      moveby = offset*Real(i-1,kind=RDbl)

      !** Loop through the atoms to change the coordinates
      Do a = 1, oldatoms
        newa = oldatoms*(i-1)
        molecule%atoms(a+newa)%name  = old%atoms(a)%name
        molecule%atoms(a+newa)%atype = old%atoms(a)%atype
        molecule%atoms(a+newa)%q0    = old%atoms(a)%q0
        molecule%atoms(a+newa)%set   = old%atoms(a)%set
        molecule%atoms(a+newa)%type  = old%atoms(a)%type
        molecule%atoms(a+newa)%r     = old%atoms(a)%r+moveby
        molecule%atoms(a+newa)%stereo= old%atoms(a)%stereo

        molecule%allaTypes(a+newa) = old%allaTypes(a)

      End Do
    End Do

    !** Correct the degrees of freedom
    Call molecule_changedof(molecule,(n+1)*molecule%dof)

    !** Generate the local coordinate axes in the body coordinate system
    !** and determine the geometry of the molecule if the molecule has
    !** generalized coordinate of type "rigid".
    !** This routine also ensures that the coordinates read in have their
    !** center-of-mass at the origin and are expressed in terms of the
    !** principle axes of rotation
    molecule%trans_xform = 0.0_RDbl
    molecule%rot_xform = matrix_identity()
    If (Trim(toupper(molecule%gcmodeltype)) == "RIGID") Then
       Call molecule_genprincipleaxes(molecule)
       Call molecule_genaxescoeff(molecule)
    End If

    !** Get rid of the molecule copy
    Deallocate(old%atoms,old%allatypes,stat=error)
    If (error/=0) Call deAllocErrDisplay(__FILE__,__LINE__,"atoms,atypes")

    !** Rebuild the connections, if requested. Otherwise return.
    If (.Not.reconnect) Return
    If (Present(fcell)) Then
      Call molecule_reconnect(molecule,fcell)
    Else
      Call molecule_reconnect(molecule)
    End If

  End Subroutine molecule_expand

  !----------------------------------------------------------------------------
  ! Makes an exact copy of the original molecule mold's atoms in
  ! the new molecule mnew atom list
  !----------------------------------------------------------------------------
  Subroutine molecule_copyAtoms(mold,mnew,start)
    Type(MolecularParams), Intent(In)    :: mold
    Type(MolecularParams), Intent(InOut) :: mnew
    Integer, Intent(In), Optional        :: start

    Integer :: error, lo, i

    !** See if we are to copy them into a different position 
    !** in the new molecule. The offset is given in start
    lo = 0
    If (Present(start)) lo = start

    !** Allocate space for the new atoms, if neccessary
    If (.Not. Associated(mnew%atoms)) Then
      Allocate(mnew%atoms(mnew%natoms),stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
            " Could not allocate space for mnew"
        Stop
      End If
    End If

    !** Loop through the atoms and copy the information
    Do i = 1, mnew%natoms
      mnew%atoms(lo+i)%name  = mold%atoms(i)%name
      mnew%atoms(lo+i)%atype = mold%atoms(i)%atype
      mnew%atoms(lo+i)%r     = mold%atoms(i)%r
      mnew%atoms(lo+i)%q0    = mold%atoms(i)%q0
      mnew%atoms(lo+i)%set   = mold%atoms(i)%set
      mnew%atoms(lo+i)%type  = mold%atoms(i)%type
    End Do
  End Subroutine molecule_copyAtoms


  !------------------------------------------------------------------------
  ! This routine gets the center of mass of a molecule using the reference
  ! coordinates from the molecule data structure
  ! Requires:  molecule -- single molecule data structure
  !------------------------------------------------------------------------
  Type(VecType) Function molecule_getrefcom(molecule)
    Type(MolecularParams), Intent(In) :: molecule

    Integer                   :: a
    Real(kind=RDbl)           :: mass_sum,mass

    mass_sum = 0.0_RDbl
    molecule_getrefcom = 0.0_RDbl

    Do a = 1,molecule_getnatoms(molecule)
      mass = atom_getmass(molecule_getnthatom(molecule,a))
      mass_sum = mass_sum + mass
      molecule_getrefcom = molecule_getrefcom + molecule%atoms(a)%r*mass
    End Do
    mass_sum = 1.0_RDbl/mass_sum
    molecule_getrefcom = molecule_getrefcom*mass_sum

  End Function molecule_getrefcom

  !-------------------------------------------------------------------------
  ! This routine gets the center of mass of a molecule using passed 
  ! coordinates and masses from the molecule data structure.  Note that
  ! the ORDER OF THE PASSED COORDINATES MUST MATCH that in the molecule 
  ! structure.
  ! Requires: molecule -- single molecule data structure
  !           xyzcoords -- xyz coordinates for each atom
  !------------------------------------------------------------------------
  Type(VecType) Function molecule_getcom(molecule,xyzcoords)
    Type(MolecularParams), Intent(In)        :: molecule
    Type(VecType), Dimension(:), Intent(In)  :: xyzcoords

    Integer                   :: a
    Real(kind=RDbl)           :: mass_sum,mass

    mass_sum = 0.0_RDbl
    molecule_getcom = 0.0_RDbl
    Do a = 1,molecule_getnatoms(molecule)
      mass = atom_getmass(molecule_getnthatom(molecule,a))
      mass_sum = mass_sum + mass
      molecule_getcom = molecule_getcom + xyzcoords(a)*mass
    End Do
    mass_sum = 1.0_RDbl/mass_sum
    molecule_getcom = molecule_getcom*mass_sum

  End Function molecule_getcom

  !------------------------------------------------------------------------
  ! This routine sets up a coordinate system at the center-of-mass of the
  ! molecule along the principle axes of the body (only for rigid bodies).
  ! It also expresses the center-of-mass coordinates in the principle 
  ! coordinate system of the body.
  ! Requires: molecule -- the molecule's data structure
  !------------------------------------------------------------------------
  Subroutine molecule_genprincipleaxes(molecule)
    Type(MolecularParams), Intent(InOut) :: molecule

    Integer                         :: i, j, a, mass, natoms, atype, ierr
    Type(VecType)                   :: rcom
    Real(kind=RDbl)                 :: rsq, cmp1, cmp2, moi
    Real(kind=RDbl), Dimension(3,3) :: mInertia, eigvec
    Real(kind=RDbl), Dimension(3)   :: eigval, subdiagel
    Type(VecType), Dimension(molecule%natoms)  :: atomcoords
    Character(len=lstrLen)          :: string

    natoms = molecule%natoms
    If (natoms < 3) Return
    
    !** Translate the molecule such that center-of-mass is at the origin
    rcom = molecule_getrefcom(molecule)
    molecule%trans_xform = (-1.0_RDbl)*rcom

    If (molecule%translateflag) Then   
      Do i = 1,natoms
        molecule%atoms(i)%r = molecule%atoms(i)%r + molecule%trans_xform
      End Do

      !** Give feeback
      If (mag(rcom) > zeroTolerance) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            " WARNING: Molecule '",Trim(molecule%name), &
            "' was not centered at origin, now shifted"
        Write(0,'(a)') &
            "  Note that the NOTRANSFORM flag will prevent this shifting"
        string = vector_display(rcom, "f12.6")
        Write(0,'(2x,a,a40)') "Former center-of-mass: ",Trim(string)
      End If
    Else
      molecule%trans_xform = (/0.0_RDbl,0.0_RDbl,0.0_RDbl/)
    End If

    Do i = 1,natoms
      atomcoords(i) = molecule%atoms(i)%r
    End Do

    !** Find the moments of inertia
    Do i = 1,3
      Do j = 1,3
        moi = 0.0_RDbl
        Do a = 1,natoms
          atype = molecule%atoms(a)%atype
          mass  = atom_getmass(atype)
          If (i == j) Then
            !** Do the diagonal terms
            rsq   = vector_getnormsq(atomcoords(a))
            cmp1  = vector_getcomp(atomcoords(a), i)
            moi   = moi + mass*(rsq - cmp1*cmp1)
          Else
            !** Do the off-diagonal terms
            cmp1  = vector_getcomp(atomcoords(a), i)
            cmp2  = vector_getcomp(atomcoords(a), j)
            moi   = moi - mass*cmp1*cmp2
          End If
        End Do
        mInertia(i,j) = moi
      End Do
    End Do

    !** Diagonalize the moment of inertia matrix and get the eigenvalues
    !** and eigenvectors.  The eigenvectors are the principle axes
    Call tred2(3, 3, mInertia, eigval, subdiagel, eigvec)
    Call tql2(3, 3, eigval, subdiagel, eigvec, ierr)
    If (ierr /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not diagonlize the moment-of-inertia matrix"
      Stop
    End If
    
    !** Set the Principle Moments of Inertia
    molecule%pMoI = eigval

    !** Set the transformation matrix
    molecule%rot_xform = eigvec
    molecule%rot_xform = matrix_transpose(molecule%rot_xform)

    !** Rotate the coordinate system to coincide with the principal axes
    If (molecule%rotateflag) Then
      Do a = 1,natoms
        molecule%atoms(a)%r = molecule%rot_xform*molecule%atoms(a)%r
      End Do
      Write(0,'(1x,2a,i4,4a)') __FILE__," : ",__LINE__, &
          " WARNING: Molecule '",Trim(molecule%name), &
          "' has been realigned to place"
      Write(0,'(2x,a)') "the cartesian axes along the moments of inertia"
      Write(0,'(2x,a)') &
          "Note that the NOROTATE or NOTRANSFORM flags will prevent this"
    Else
      molecule%rot_xform = matrix_identity()
    End If

  End Subroutine molecule_genprincipleaxes

  !---------------------------------------------------------------------------
  ! This routine returns the principal moments-of-inertia of "molecule"
  ! This is probably only useful for a rigid molecule
  !---------------------------------------------------------------------------
  Function molecule_getmoi(molecule)
    Type(MolecularParams), Intent(In) :: molecule
    Real(kind=RDbl), Dimension(3)     :: molecule_getmoi
    
    molecule_getmoi = molecule%pMoI
  End Function molecule_getmoi

  !---------------------------------------------------------------------------
  ! This routine returns the transformations that originally took the 
  ! coordinates in the molecule file to those in the molecule data structure
  ! Requires: molecule -- the molecule data structure
  !           trans_xform -- translation vector
  !           rot_xform -- 3x3 rotation matrix
  ! The orginal transform was:
  !   ref_vector = rot_xform*(molfile_vector + trans_xform)
  !---------------------------------------------------------------------------
  Subroutine molecule_origxform(molecule,trans_xform,rot_xform)
    Type(MolecularParams), Intent(In)     :: molecule
    Type(VecType), Intent(Out)            :: trans_xform
    Type(MatrixType), Intent(Out)         :: rot_xform
    
    trans_xform = molecule%trans_xform
    rot_xform = molecule%rot_xform

  End Subroutine molecule_origxform

  !---------------------------------------------------------------------------
  ! This routine sets up the local coordinates system by getting coefficients
  ! which when multiplied by certain atom vectors give the x, y and z axes.
  ! This routine stores the vectors and their respective coefficients for
  ! the three different coordinate axes.  These coefficients and vectors
  ! are later used for the purpose of calculating generalized eulerian angles
  ! from cartesian coordinates.  This routine also defines the basic geometry
  ! of the molecule, 'Spherical', 'Linear' or default
  ! Requires: molecule -- single molecule data structure
  !---------------------------------------------------------------------------
  Subroutine molecule_genaxescoeff(molecule)
    Type(MolecularParams), Intent(InOut) :: molecule

    !** Initialize the local coordinates
    Call localcoords_init(molecule%bodyaxes,molecule%natoms, &
        molecule%atoms%r,molecule_getrefcom(molecule),molecule%geometry)

  End Subroutine molecule_genaxescoeff

  !----------------------------------------------------------------------------
  ! Initialize the molecule coordinate information
  ! Requires: molecule -- the molecule's data structure
  !           unit -- unit to take coordinates from
  !----------------------------------------------------------------------------
  Recursive Subroutine molecule_coordinit(molecule,unit)
    Type(MolecularParams), Intent(InOut) :: molecule
    Integer, Intent(In)                  :: unit

    Integer                       :: uno, natoms, nparams
    Integer                       :: i,usecoord,error,n,index,atype
    Real(kind=RDbl)               :: mass
    Character(len=255)            :: text, stripText
    Character(len=strLen)         :: key
    Character(len=strLen), Dimension(strLen) :: params
 
    key = "Coord_Info"

    usecoord = filesrchstr(unit,key,text,.True.)
    
    If (usecoord == 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ': molecule file is missing ',trim(key)
      Stop
    End If

    stripText = stripcmnt(text)
    n = split(stripText,params)

    Select Case (Trim(tolower(params(2))))

    Case ('listed')
      !** Store the coordinate system name
      molecule%coordsys = Trim(params(3))
      molecule%gcmodeltype = Trim(params(4))
      text = adjustl(Trim(params(5)))

      !** this means molecule will be translated so that com coincides 
      !   with coord system origin
      molecule%translateflag = .True.

      If (Trim(ToUpper(text)) == "NOROTATE") Then
        molecule%rotateflag = .False.
      Else If (Trim(ToUpper(text)) == "NOTRANSFORM") Then
        molecule%rotateflag = .False.
        molecule%translateflag = .False.
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) " Molecule will be rotated so that the principle axes "
        Write(*,*) "  will coincide with coordinate system axes"
        molecule%rotateflag=.True.
      End If

      !** Read in atom information 
      Read(unit,'(a)') text
      stripText = stripcmnt(text)
      n = split(stripText,params)
      natoms = toint(params(1),'molecule natoms read')
      If (natoms > MAX_ATOMS) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' WARNING: number of atoms in molecule exceeds MAX_ATOMS'
        Write(0,'(1x,a)') 'change default if you encounter problems'
        Write(0,'(1x,a,i4)') 'natoms = ',natoms
        Write(0,'(1x,a,i4)') 'MAX_ATOMS = ',MAX_ATOMS
      End If
      molecule%natoms = natoms

      !** Allocate the atoms array
      Allocate(molecule%atoms(molecule%natoms),STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' unable to allocate molecule%atoms'
        Stop
      End If

      mass = 0.0_RDbl
      Do i = 1,natoms
        Read(unit,'(a)') text
        stripText = stripcmnt(text)
        nparams = split(stripText,params)

        index = toint(params(1),'molecule coordinate line read')
        molecule%atoms(i)%r%comp(1) = toreal(params(2))
        molecule%atoms(i)%r%comp(2) = toreal(params(3))
        molecule%atoms(i)%r%comp(3) = toreal(params(4))
        molecule%atoms(i)%name = trim(params(5))
        molecule%atoms(i)%q0 = toreal(params(6))
        molecule%atoms(i)%set = toint(params(7),'set number conversion')
        molecule%atoms(i)%type = toint(params(8),'type number conversion')

        atype = atom_gettypename(molecule%atoms(i)%name)
        If (atype == 0) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              ': atom name does not exist in atoms structure: ',&
              Trim(molecule%atoms(i)%name)
          Stop          
        End If
        molecule%atoms(i)%atype = atype
        
        !** Calculate a molecular weight
        mass = mass + atom_getmass(atype)
      End Do
      molecule%mass = mass

    Case Default
      !** Assume that it is a filename. Try to open it.
      uno = isfileopen(Trim(tolower(params(2))))
      If (uno < 0) Then
        uno = file_getunit(trim(params(2)))
        Open(uno,file=Trim(params(2)),iostat=error)
        If (error /= 0) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              " : Cannot open connection file ",trim(params(2))
          Stop
        End If
      End If
      Call molecule_coordinit(molecule,uno)
      Close(uno)

    End Select
    

    !** Generate the local coordinate axes in the body coordinate system
    !** and determine the geometry of the molecule if the molecule has
    !** generalized coordinate of type "rigid".
    !** This routine also ensures that the coordinates read in have their
    !** center-of-mass at the origin and are expressed in terms of the
    !** principle axes of rotation
    molecule%trans_xform = 0.0_RDbl
    molecule%rot_xform = matrix_identity()
    If (Trim(toupper(molecule%gcmodeltype)) == "RIGID") Then
       Call molecule_genprincipleaxes(molecule)
       Call molecule_genaxescoeff(molecule)
    End If

  End Subroutine molecule_coordinit
  
  !---------------------------------
  !  Initialize steroechemistries
  !---------------------------------
  Subroutine molecule_stereoinit(molecule,unitno)
    
    Type(MolecularParams) :: molecule
    Integer, Intent(IN)   :: unitno
    
    Integer :: i,j,usecon,num_stereo
    Character(len=255) :: text, stripText
    Character(len=strLen), Dimension(strLen) :: params
    Character(len=strLen) :: key
    
    ! Initialize all stereochemistries to 'N'
    Do i = 1,molecule%natoms
      molecule%atoms(i)%stereo = 'N'
    End Do
    
    key = "Stereo_Info"
    
    usecon = filesrchstr(unitno,key,text,.True.)
    
    If (usecon == 0) Return
    
    j = split(text,params)
    
    Select Case (Trim(tolower(params(2))))    
       
    Case ('listed')
       
      Read(unitno,'(a)') text
      stripText = stripcmnt(text)
      j = split(stripText,params)
       
      num_stereo = toint(params(1),'number of stereo centers in molecule read')
       
      Do i=1,num_stereo
        Read(unitno,'(a)') text
        j = split(text,params)
        
        ! Hopefully this assignment will truncate string to first character
        ! if necessary
        molecule%atoms(toint(params(1)))%stereo = toupper(Trim(params(2)))
      End Do
       
    Case Default
       ! Stereochemistries remain initialized at 'N'
       
    End Select
    
  End Subroutine molecule_stereoinit
  
  !----------------------------------------------------------------------------
  ! Initialize the connections
  !----------------------------------------------------------------------------
   Recursive Subroutine molecule_connectinit(molecule,unitno)
    Type(MolecularParams) :: molecule
    Integer, Intent(In)   :: unitno

    Integer                    :: i,i2,j,usecon,error
    Integer                    :: atom1, atom2, atype1, atype2, natypes
    Integer                    :: nconnect, uno   
    Real(kind=RDbl)            :: length
    Character(len=strLen)      :: key
    Character(len=255)         :: text, stripText
    Character(len=strLen), Dimension(strLen)  :: params
    Type(VecType), Dimension(molecule%natoms) :: coords

    Write(*,*) __FILE__,__LINE__,' Generating connections'

    key = "Connect_Info"
!!$    usecon = filesrchwocomment(unitno,key,text,.True.)
    ! find the line with the key
    usecon = filesrchstr(unitno,key,text,.True.)
    ! make sure the line was not commented out
    If (trim(stripcmnt(text))=="") usecon=0
    If (usecon == 0) Return

    !** Initialize the connection matrix for this molecule
    natypes = atom_getntypes()
    Call connects_initarray(molecule%connections,natypes)

    !** Allocate the connect array
    Allocate(molecule%connect(molecule%natoms), STAT=error) 
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          " : Error allocating array 'connect'"
    End If

    !** initialize the connections
    molecule%connect(1:molecule%natoms)%nconnect = 0
    Do i = 1, molecule%natoms
      Call connects_initmaxcon(molecule%connect(i), &
          atom_nbonds(molecule%atoms(i)%atype))
    End Do

    ! ** dont consider if it is comented out
    text=stripcmnt(text)
    j = split(text,params)

    Select Case (Trim(tolower(params(2))))
    Case ('generate')

      !MDEBUG
      Write(*,*) "Generating Connections for ",molecule%name
      Read(unitno,'(a)') text
      j = split(text,params)
      nconnect = toint(params(1))
      Do i = 1, molecule%natoms
        coords(i) = molecule%atoms(i)%r
      End Do
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__

      !** Read in and store the connection types. Loop through each
      !** connection type
      Do i = 1, nconnect
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__

        !** Read the line, strip out comments, and split the line
        Read(unitno,'(a)') text
        stripText = stripcmnt(text)
        j = split(stripText,params)

        !** Get the atom types of each atom and the required length
        atype1 = atom_gettypename(Trim(params(1)))
        atype2 = atom_gettypename(Trim(params(2)))

        If ((atype1==0).Or.(atype2==0)) Then
          Write(*,*) __FILE__,__LINE__
          Write(*,'(3a)') "WRONG PAIRS: ", params(1), params(2)
          Write(*,*) "wrong atom types in the connect section of molecule"
          stop
        Endif
        length = toreal(params(3))

        If (j > 4) Then
          !** An extra parameter is specified, which includes
          !** things like the shellcore designation
          Call connects_initconnection(molecule%connections, &
              atype1, atype2, length, toreal(params(4)), Trim(params(5)))
        Else If (j > 3) Then
          ! user specified a tolerance in cf
          Call connects_initconnection(molecule%connections, &
              atype1, atype2, length, toreal(params(4)))
        Else
          Call connects_initconnection(molecule%connections, &
              atype1, atype2, length)
        End If
      End Do

      !** Actually generate the connections list
      !MDEBUG
      !LC       Write(*,*) __FILE__,__LINE__," Calling generateconnections"
      Call connects_generateAll(molecule%connect, &
          molecule%atoms(1:molecule%natoms)%atype, coords, molecule%connections)

!!$       !MDEBUG
!!$       Write(*,*) __FILE__,__LINE__," Displaying connections: "
!!$          Call connects_displayall(molecule%connect,6)

    Case ('listed')

      Read(unitno,'(a)') text
      j = split(text,params)
      text = stripcmnt(text)
      nconnect = toint(params(1))

      ! Reading in bond lengths for various types of connections
      Do i = 1, nconnect
        Read(unitno,'(a)') text
        text = stripcmnt(text)
        j = split(text,params)
        atype1 = atom_gettypename(trim(params(1)))
        atype2 = atom_gettypename(trim(params(2)))
        length = toreal(params(3))
        Call connects_initconnection(molecule%connections, &
            atype1, atype2, length)
      End Do

      ! Now read in the connectivities immediately following
      ! There should be same number of entries as there are atoms
      Do i = 1, molecule%natoms
        Read(unitno,'(a)') text
        j = split(text,params)
        If (j>1) Then
          atom1 = toint(params(1))
          atype1 = molecule_getnthatom(molecule,atom1)
          Do i2 = 2,j
            atom2 = toint(params(i2))
            atype2 = molecule_getnthatom(molecule,atom2)
            Call connects_connectatoms(molecule%connect,atom1,atom2,&
                molecule%connections%bonds(atype1,atype2)%length)
          End Do
        End If
      End Do

    Case Default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) " Calling molecule_connectinit recursively"
      Write(*,*) "from file  named : ", Trim(params(2))


      !** we assume this must be a file. Open it and read it.
      uno = isfileopen(Trim(tolower(params(2))))
      If (unitno < 0) Then
        uno = file_getunit(trim(params(2)))
        Open(uno,file=Trim(params(2)),iostat=error)
        If (error /= 0) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              " : Cannot open connection file ",trim(params(2))
          Stop
        End If
      End If

      !** Recursion
      Call molecule_connectinit(molecule,uno)
      Close(uno)

    End Select

  End Subroutine molecule_connectinit

  !----------------------------------------------------------------------------
  ! Regenerates the molecules's connections based on the 
  ! information already stored in the molecules's connections 
  ! Requires:  spc -- species number
  !            fcell -- optional fundamental cell if PBCs are to be used
  !----------------------------------------------------------------------------
  Subroutine molecule_reconnect(molecule,fcell)
    Type(MolecularParams), Intent(InOut)           :: molecule
    Type(Fundamental_Cell), Intent(In), Optional   :: fcell    

    Integer :: error,i

    !** make sure it's really necessary first
    If (.Not. Associated(molecule%connect)) Return

    !** Initialize the connections by first dumping the
    !** original connection list
    Deallocate(molecule%connect,stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          " Could not deallocate molecule%connect"
      Stop
    End If

    !** Now reallocate it based on the number of atoms
    Allocate(molecule%connect(molecule%natoms),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i6,a,i9)') __FILE__,":",__LINE__, &
          " Could not allocate molecule%connect of size ",molecule%natoms
      Stop
    End If
      
    !** Initialize each entry in molecule%connect
    molecule%connect(1:molecule%natoms)%nconnect = 0
    Do i = 1, molecule%natoms
      Call connects_initmaxcon(molecule%connect(i), &
          atom_nbonds(molecule%atoms(i)%atype))
    End Do
      
    !** Regenerate the connections first
    If (Present(fcell)) Then
      Call connects_regenerateAll(molecule%connect, &
          molecule%atoms(1:molecule%natoms)%atype, &
          molecule%atoms(1:molecule%natoms)%r, molecule%connections, fcell)
    Else
      Call connects_regenerateAll(molecule%connect, &
          molecule%atoms(1:molecule%natoms)%atype, &
          molecule%atoms(1:molecule%natoms)%r, molecule%connections)
    End If

  End Subroutine molecule_reconnect

  !----------------------------------------------------------------------------
  ! Check the connections for the molecule to make sure they are complete
  !----------------------------------------------------------------------------
  Logical Function molecule_checkconnect(molecule)
    Type(MolecularParams), Intent(In) :: molecule
    Integer :: error, i, nbonds, nconnect

    !** Set the default value to be TRUE.
    molecule_checkconnect = .True.
    error = 0

    !** Loop through each atom and compare the maximum number of atom
    !** connections and the actual number of connections.
    Do i = 1, molecule%natoms

      !** Number of bonds required by the atom type
      nbonds = atom_nbonds(molecule%atoms(i)%atype)

      !** Cycle if there are no bonds for this atom. Stop if this is
      !** a molecule with more than 1 atom but the atom has no connections
      !** required.
      If ((nbonds == 0).And.(molecule%natoms > 1)) Then
        Write(0,'(a,i5,a)') __FILE__,__LINE__, &
            " : Error, molecule with more than 1 atom has an atom type with &
            &zero required bonds."
        Write(0,'(5x,3a,i5)') "molecule : ",Trim(molecule%name), &
            "    atom number : ",i
        Stop
      Else If (nbonds == 0) Then
        Cycle
      End If

      !** Actual number of bonds created by connects
      nconnect = connects_nconnected(molecule%connect(i))

      !** Compare the numbers
      If (nbonds /= nconnect) Then
        molecule_checkconnect = .False.
        Write(0,'(a,i5,3a,i2,a,i2,a)') 'Atom ',i,' of ', &
            Trim(molecule_getname(molecule)),' only has ',nconnect,' of ', &
            nbonds,' required bonds.'
        Write(0,'(a)') &
            'Please check the specified bond length and/or atom types.'
        Write(0,'(a)') 'Error info:'
        Write(0,'(2a)') 'Atom Type :',atom_getname(molecule%atoms(i)%atype)
        Write(0,'(2a)') 'Connected to :',connects_display(molecule%connect(i))
      End If
    End Do
    If (.Not.molecule_checkconnect) Then
      Write(0,'(a,i4,a)') __FILE__,__LINE__,' : Error in number of bonds'
      Stop
    End If
  End Function molecule_checkconnect

  !----------------------------------------------------------------------------
  ! Return a list of the neighbors given the molecule and atom number a
  !----------------------------------------------------------------------------
  Integer Function molecule_getneighbors(molecule,a,nlist)
    Type(MolecularParams), Intent(In) :: molecule
    Integer, Intent(In)               :: a
    Integer, Dimension(:), Intent(Out):: nlist

    !** Call connects to get the neighbors
    molecule_getneighbors = connects_atomsxaway(molecule%connect,a,1,nlist)

  End Function molecule_getneighbors

  !----------------------------------------------------------------------------
  ! Return the pointer to the molecule atom connection list
  ! Requires:  molecule -- the molecule's data structure
  !            connections -- returned pointer to the connections
  !----------------------------------------------------------------------------
  Subroutine molecule_getconnects(molecule,connections)
    Type(MolecularParams), Intent(In)            :: molecule
    Type(AtomConnections), Dimension(:), Pointer :: connections

    connections => molecule%connect

  End Subroutine molecule_getconnects

  !----------------------------------------------------------------------------
  ! Dump the fields of the "molecules" structure
  ! Requires: molecule -- molecule data structure
  !           unit -- unit number to dump into
  !           indent -- indentation from left margin
  !----------------------------------------------------------------------------
  Subroutine molecule_display(molecule,unit,indent)
    Type(MolecularParams), Intent(In)   :: molecule
    Integer, Intent(In)                 :: unit,indent

    Integer                             :: atom
    Character(len=indent)               :: blank
    Character(len=strLen)               :: string

    blank = Repeat(' ',indent)
    
    Write(unit,'(3a)') blank,'Molecule Name       : ', &
        Trim(molecule%name)
    string = real2str(molecule%mass,5)
    Write(unit,'(3a)') blank,'Molecule Mass       : ', &
        Trim(string)
    Write(unit,'(3a)') blank,'Molecule File       : ', & 
        Trim(molecule%file)
    string = int2str(molecule%dof)
    Write(unit,'(3a)') blank,'Degrees of freedom  : ', &
        Trim(string)
    Write(unit,'(3a)') blank,'DOF origin          : ', &
        Trim(molecule%dof_origin)
    Write(unit,'(2a,3i4)') blank,'body axes atom nos  : ', &
        molecule%bodyaxes%atomnos(1),molecule%bodyaxes%atomnos(2), &
        molecule%bodyaxes%atomnos(3)
    string = int2str(molecule%natoms)
    Write(unit,'(3a)') blank,'Number of atoms     : ', & 
        Trim(string)
    string = real2str(molecule_netcharge(molecule),4)
    Write(unit,'(3a)') blank,'Net Charge          : ', & 
        Trim(string)
    string = real2str(molecule_calcdipole(molecule))
    Write(unit,'(3a)') blank,'Dipole Moment(debye): ', & 
        Trim(string)
    string = vector_display(molecule%trans_xform,'f10.5')
    Write(unit,'(3a)') blank,'Translation xform   : ', & 
        Trim(string)
    Write(unit,'(2a)') blank,'Rotational xform    : '
    Call matrix_display(molecule%rot_xform,'f10.5',unit)

    If (molecule%natoms <= 10) Then
      Do atom = 1,molecule%natoms
        Call molecule_display_atom(molecule,unit,atom)
      End Do
    Else
      Write(unit,'(2a)') blank,'Atom list omitted (natoms > 10)'
    End If

    Write(unit,'(2a)')      blank,'Connection Information : '
    Call molecule_connectdisplay(molecule,unit,indent+2)

  End Subroutine molecule_display

  !----------------------------------------------------------------------------
  ! Dump the information for a single atom
  !----------------------------------------------------------------------------
  Subroutine molecule_display_atom(molecule,unit,atom)
    Type(MolecularParams)    :: molecule
    Integer, Intent(IN)      :: unit,atom
    Integer                  :: j

    Write(unit,'(6x,2a)')      'Atom Name           : ', &
        Trim(molecule%atoms(atom)%name)
    Write(unit,'(6x,a,i4)')    'Atom Type           : ', & 
        molecule%atoms(atom)%atype
    Write(unit,'(6x,a,3f8.4)')'Atom Position       : ', & 
        (molecule%atoms(atom)%r%comp(j),j=1,3)
    Write(unit,'(6x,a,f8.4)')  'Equilibrium charge  : ', & 
        molecule%atoms(atom)%q0
    Write(unit,'(6x,a,2i4)')   'MD Set, type        : ', & 
        molecule%atoms(atom)%set,molecule%atoms(atom)%type
    Return
  End Subroutine molecule_display_atom

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine molecule_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    !** Warn this isn't done
    Write(unitno,'(a,i6,a)') __FILE__,__LINE__, &
        ": Sorry, sample molecule file section not complete"
    
  End Subroutine molecule_sampleCF

  !----------------------------------------------------------------------------
  ! Dump an xyz file for a molecule
  !----------------------------------------------------------------------------
  Subroutine molecule_molec2xyz(molecule,filename,opt_cmnt)
    Type(MolecularParams), Intent(In)    :: molecule
    Character(*), Intent(In)             :: filename
    Character(*), Optional               :: opt_cmnt

    Integer                              :: i,atype
    Character(len=strLen)                :: display_format, comment
    Type(XYZ_Entry), Dimension(molecule%natoms)   :: entries

    !set the comment
    If (Present(opt_cmnt)) Then
      comment=opt_cmnt
    Else
      comment = 'Created by MUSIC code'
    End If

    Do i = 1,molecule%natoms
      entries(i)%r = molecule%atoms(i)%r
      atype = molecule%atoms(i)%atype
      entries(i)%symbol = atom_getsymbol(atype)
    End Do

    display_format = 'f12.5'
    Call visxyz_dump(entries,filename,display_format,comment)

    Return
  End Subroutine molecule_molec2xyz


  !-----------------------------------------------
  ! Display the information about connections
  !-----------------------------------------------
  Subroutine molecule_connectdisplay(molecule,unit,indent)
    Type(MolecularParams) :: molecule
    Integer, Intent(In)   :: unit,indent

    Integer                :: i
    Character(len=indent)  :: blank
    Character(len=lstrLen) :: string

    blank = Repeat(' ',indent)

    If (Associated(molecule%connect)) Then
      Write(unit,'(3a)') blank,"Connection Information for ",molecule%name
      If (Size(molecule%connect,1) < 50) Then
        Do i = 1,Size(molecule%connect)
          string = connects_display(molecule%connect(i))
          Write(unit,'(2a,i5,2a)') blank,"Atom ",i,": ", &
              Trim(string)
        End Do
      Else
        Write(unit,'(3a)') blank,"More than 50 atoms, connection information" &
            ," will not be displayed"
      End If
    Else
      Write(unit,'(3a)') blank,"No connectivity information for ", &
          Trim(molecule%name)
    End If
  End Subroutine molecule_connectdisplay


End Module molecule








