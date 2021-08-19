!-----------------------------------------------------------------------
! This module contains the data type that stores the information about 
! the molecule types on a multiple species level.  The type that it 
! defines and operates on encompasses all of the molecule type 
! definitions.
! In many cases this module is simply a front-end for the routines in
! the molecule module.  It serves the purpose of actually storing the
! molecular definitions.
!-----------------------------------------------------------------------

Module molecules

  Use defaults, Only: RDbl, strLen, lstrLen, d_ctrl_file, d_ss_file, &
      dashedline, scalepe
  Use utils, Only: split, filesrchstr, stripcmnt, toupper, findint, &
      getlineswpair, allocErrDisplay, getpath, int2str, real2str, cleanstring
  Use file, Only: file_gettype, file_settag, file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use matrix, Only: MatrixType
  Use fundcell, Only: Fundamental_Cell
  Use connects, Only: AtomConnections
  Use localcoords, Only: LocalCoordInfo
  Use molecule, Only: MolecularParams, molecule_init, molecule_getname, &
      molecule_ischarged, molecule_getfilename, &
      molecule_getbodyaxes, molecule_getnatoms, molecule_getcharge, &
      molecule_getnthatom, molecule_getdof, molecule_getatomcoords, &
      molecule_getnatomtypes, molecule_getmass, molecule_changedof, &
      molecule_getgcmodeltype, molecule_getmoi, molecule_getcomdefn, &
      molecule_display, molecule_expand, molecule_netcharge, &
      molecule_molec2xyz, molecule_reconnect, molecule_origxform, &
      molecule_sampleCF, molecule_isCoreShell, molecule_getcom, &
      molecule_getrefcom, molecule_getneighbors, molecule_checkconnect, &
      molecule_getconnects, molecule_atomMass
  Use atom, Only: atom_checkinit, atom_invmass, atoms, atom_getsymbol, &
      atom_getname, atom_getmass

  Implicit None
  Save
  
  Private
      !  ------------ very basic functions and init/display -------
  Public :: Molecule_List, molecs, molecules_gettype, molecules_getmolecule, &
      molecules_getpointer, molecules_getnsorbs, molecules_getfilename, &
      molecules_checkinit, molecules_init, molecules_exists, &
      molecules_display, molecules_dump2xyz, molecules_sampleCF, &
      !  ------------ name/ mass / com/ dof  -------
      molecules_name, molecules_getmass, molecules_getcomdefn, &
      molecules_getcom, molecules_getdof, molecules_changedof, &
      molecules_dofIsOK,  &
      ! ------ coords/gcmodel -----------
      molecules_getgcmodeltype, molecules_getcoords, molecules_isflexible, &
      !  ------------ atoms-related functions  -------
      molecules_getnatoms, molecules_getnthatom, molecules_getatype, &
      molecules_getnatomtypes, molecules_getinvmasses, molecules_atomsymbol, &
      molecules_getatomcoords, molecules_getmaxnatoms, molecules_getaname, &
      molecules_getneighbors,  molecules_getmasses, molecules_AtomMass, & 
      ! ------- related to Moment of Inertia ---------
      molecules_getmoi, molecules_getbodyaxes, molecules_origxform, & 
      ! ------- charged molecules ----------
      molecules_getcharge, molecules_getcharges, molecules_netcharge, &
      molecules_ischarged, &
      ! ---------- connections ----------------
      molecules_expand, molecules_reconnect, molecules_getconnects, & 
      ! ------ I dont know where these fit ----------
      molecules_isCoreShell, molecules_getlinewtag

  !** NAG compiler workaround for strLen
  Integer :: dummyMolecules = strLen

  Type Molecule_List
    Integer                                        :: nmolec_types 
    Type(MolecularParams), Dimension(:), Pointer   :: mo
  End Type Molecule_List                           

  Type(Molecule_List), Target                      :: molecs

  Interface getmolecule
    Module procedure molecules_getmolecule
  End Interface

  Character(strLen), Parameter :: molecules_idstring = 'Molecule Types'

Contains
  !----------------------------------------------------------------------------
  ! Initialize the molecule information and returns the number of molecule
  ! types initialized.
  !----------------------------------------------------------------------------
  Integer Function molecules_init()
    Integer                               :: i,j,ios,error
    Integer                               :: unit,lineno,spc
    Character(len=lstrLen)                :: line,srchstr,filename
    Character(len=strLen)                 :: molecName, molecFile

    !** Make sure that the "atoms" structure has been initialized
    If (.NOT. (atom_checkinit(atoms))) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          "Initialize the atoms structure before initializing the molecules"
      Stop
    End If

    !** Read the species-species section of the control file
    Call file_gettype(d_ctrl_file,filename,unit)
    Open(unit=unit, file=filename, status='unknown', IOSTAT=ios)
    If (ios /= 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
           'Could not open file ',trim(filename)
      Stop
    End If 
    
    srchstr = molecules_idstring
    lineno = filesrchstr(unit, srchstr, line)
    If (lineno == 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
           ' Could not find search string '//Trim(srchstr), &
           ' in file '//Trim(filename)
      Stop
    End If

    Read(unit,*) molecs%nmolec_types
    molecules_init = molecs%nmolec_types
    Allocate(molecs%mo(molecs%nmolec_types),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'molecs%mo')

    Read(unit,'(a)') line 

    !** Read the individual molecule names and files
    Do spc = 1,molecs%nmolec_types
      Read(unit,*) molecName
      Read(unit,*) molecFile
      !** Read the blank line if we are not reading the last
      !** molecule
      If (spc /= molecs%nmolec_types) Then
        Read(unit,'(a)') line
      End If

      Call molecule_init(molecs%mo(spc),molecName,molecFile)
    End Do

    !** Check here to make sure that molecule names are not duplicated
    !** Probably need it in atoms too
    If ((molecs%nmolec_types)>1) Then
      Do i = 1,molecs%nmolec_types
        Do j = i+1,molecs%nmolec_types
          If (molecule_getname(molecs%mo(i))==&
              molecule_getname(molecs%mo(j))) Then
            Write(*,*) "The molecule - ",Trim(molecs%mo(i)%name),&
                "is duplicated in the list in ctrlfile"
            Stop
          End If
        End Do
      End Do
    End If
    
    Close(unit = unit)
  End Function molecules_init

  !----------------------------------------------------------------------------
  ! Writes a sample of the control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine molecules_sampleCF(unitno)
    Integer, Intent(In) :: unitno
    
    Write(unitno,'(a)') molecules_idstring
    Write(unitno,'(a,t30,a)') 'Integer','# Number of molecule types listed'
    Write(unitno,'(a,t30,a)') 'Character', &
        '# Molecule-molecule interaction file'
    Write(unitno,'(t30,a)') '# Blank line (required), then for each molecule:'
    Write(unitno,'(a,t30,a)') 'Character','# Name of molecule'
    Write(unitno,'(a,t30,a)') 'Character','# Molecule file filename'
    Call molecule_sampleCF(unitno)
    Write(unitno,'(t30,a)') '# Blank line (required)'
  End Subroutine molecules_sampleCF

  !-------------------------------------------------------
  ! Checks if the atoms of the molecule "spc" have charge
  !-------------------------------------------------------
  Logical Function molecules_ischarged(spc)
    Integer, Intent(in)     :: spc

    molecules_ischarged = molecule_ischarged(molecs%mo(spc))
    Return
  End Function molecules_ischarged

  !----------------------------------------------------------------------------
  ! Returns true if the the atom is part of a core-shell. Optionally
  ! returns the atom number of the corresponding core or shell.
  !----------------------------------------------------------------------------
  Logical Function molecules_isCoreShell(spc,atom,isCore,partner)
    Integer, Intent(In) :: spc, atom
    Logical, Intent(Out):: isCore
    Integer, Intent(Out), Optional :: partner

    If (Present(partner)) Then
      molecules_isCoreShell =molecule_isCoreShell(molecs%mo(spc),atom,isCore, &
          partner)
    Else
      molecules_isCoreShell =molecule_isCoreShell(molecs%mo(spc),atom,isCore)
    End If
  End Function molecules_isCoreShell



  !----------------------------------------------------------------------------
  ! Returns true if the the molecule dof is alright, dont call frequently
  !----------------------------------------------------------------------------
  Logical Function molecules_dofIsOK(spc)
    Integer, Intent(In) :: spc

    Integer :: natoms,dof

    natoms=molecules_getnatoms(spc)
    dof=molecules_getdof(spc)

    molecules_dofIsOK=.True.
    If (natoms>2) Then
      If ((dof<6).Or.(dof>3*natoms))     molecules_dofIsOK=.False.  
    Elseif(natoms==2) Then
      If ((dof>6).Or.(dof<5)) molecules_dofIsOK=.False.  
    Elseif(natoms==1) Then
      If (dof/=3) molecules_dofIsOK=.False.  
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
  End Function molecules_dofIsOK

  !----------------------------------------------------------------------------
  ! Returns true if the the molecule is internally flexible 
  !----------------------------------------------------------------------------
  Logical Function molecules_isflexible(spc,opt_dof)
    Integer, Intent(In) :: spc
    Integer, Intent(Out), Optional :: opt_dof
    Integer :: dof, natoms
    Logical ::  isflexible 

    ! many checks here are unnecessary, just extra safety
    isflexible=.False.
    dof=molecules_getdof(spc)
    natoms=molecules_getnatoms(spc)

    If (.Not.(molecules_dofIsOK(spc))) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif


    If (natoms>2) Then
      If (dof>6) isflexible=.True.
    Elseif(natoms==2) Then
      If (dof==6) isflexible=.True.
    Elseif(natoms==1) Then
      isflexible=.False.
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    molecules_isflexible=isflexible
    If (Present(opt_dof)) opt_dof=dof
  End Function molecules_isflexible

  !----------------------------------------
  ! Returns the molecule DATA for molec n
  !----------------------------------------
  Type(MolecularParams) Function molecules_getmolecule(n)

    Integer, Intent(In) :: n
    molecules_getmolecule = molecs%mo(n)

  End Function molecules_getmolecule

  !----------------------------------------------------------------------------
  ! Returns the neighboring atoms given atom number a for molecule spc
  ! Requires:  spc -- species number
  !            atm -- atom number
  !            nlist -- returned list of neighboring atom number
  !----------------------------------------------------------------------------
  Integer Function molecules_getneighbors(spc,atm,nlist)
    Integer, Intent(In) :: spc
    Integer, Intent(In) :: atm
    Integer, Dimension(:), Intent(Out) :: nlist
    
    molecules_getneighbors = molecule_getneighbors(molecs%mo(spc),atm,nlist)

  End Function molecules_getneighbors

  !----------------------------------------------------------------------------
  ! Return the pointer to the molecule atom connection list
  ! Requires:  spc -- species number
  !            connections -- returned pointer to the connections
  !----------------------------------------------------------------------------
  Subroutine molecules_getconnects(spc,connections)
    Integer, Intent(In)                          :: spc
    Type(AtomConnections), Dimension(:), Pointer :: connections

    Call molecule_getconnects(molecs%mo(spc),connections)

  End Subroutine molecules_getconnects

  !-----------------------------------------------------
  ! Returns the molecule filename for molecule type n
  !-----------------------------------------------------
  Function molecules_getfilename(n)
    Character(len=strLen) :: molecules_getfilename
    Integer, Intent(In) :: n
    molecules_getfilename = molecule_getfilename(molecs%mo(n))
  End Function molecules_getfilename

  !-------------------------------------------------------------------------
  ! Returns the details of the body coordinate system
  !-------------------------------------------------------------------------
  Type(LocalCoordInfo) Function molecules_getbodyaxes(spc)
    Integer, Intent(in) :: spc
    
    molecules_getbodyaxes = molecule_getbodyaxes(molecs%mo(spc))
  End Function molecules_getbodyaxes

  !----------------------------------------------------------------------------
  ! Gets the number of species
  !----------------------------------------------------------------------------
  Integer Function molecules_getnsorbs()
    molecules_getnsorbs = molecs%nmolec_types
  End Function molecules_getnsorbs

  !----------------------------------------------------------------------------
  ! Gets the number of atoms in one species
  ! I added a check for spc being proper, this might increase cpu time ?
  !----------------------------------------------------------------------------
  Integer Function molecules_getnatoms(spc)
    Integer, Intent(In)              :: spc

    If ((spc<1) .Or. (spc>Size(molecs%mo,1))) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "unreasonable species number passed here", spc
      Stop
    Else
      molecules_getnatoms = molecule_getnatoms(molecs%mo(spc))
    End If

  End Function molecules_getnatoms

  !----------------------------------------------------------------------------
  ! Get the maximum number of atoms in the collection of molecule types
  ! Requires:  skipspc -- species numbers to skip (could be fixed species)
  !----------------------------------------------------------------------------
  Integer Function molecules_getmaxnatoms(skipspc)
    Integer, Dimension(:), Intent(In)    :: skipspc

    Integer        :: spc,natoms

    molecules_getmaxnatoms = 0

    Do spc = 1,molecules_getnsorbs()
      If (findint(skipspc,spc) /= 0) Cycle
      natoms = molecules_getnatoms(spc)
      If (natoms > molecules_getmaxnatoms) molecules_getmaxnatoms = natoms
    End Do

  End Function molecules_getmaxnatoms

  !----------------------------------------------------------------------------
  ! Gets the charge of atom "atomno" of species type "spc"
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function molecules_getcharge(spc, atomno)
    Integer, Intent(IN)              :: spc
    Integer, Intent(in)              :: atomno

    molecules_getcharge = molecule_getcharge(molecs%mo(spc), atomno)
  End Function molecules_getcharge

  !----------------------------------------------------------------------------
  ! Gets the nth atom type of species spc
  !----------------------------------------------------------------------------
  Integer Function molecules_getnthatom(spc,nthAtom)
    Integer, Intent(In) :: spc,nthAtom
    molecules_getnthatom = molecule_getnthatom(molecs%mo(spc),nthAtom)
  End Function molecules_getnthatom

  !----------------------------------------------------------------------------
  ! Gets the number of degrees of freedom for molecule type spc
  !----------------------------------------------------------------------------
  Integer Function molecules_getdof(spc,dof_origin)
    Integer, Intent(In)            :: spc
    Character(*), Optional         :: dof_origin
    molecules_getdof = molecule_getdof(molecs%mo(spc),dof_origin)
  End Function molecules_getdof

  !----------------------------------------------------------------------------
  ! Change the number of degrees of freedom for the molecule of type spc
  !----------------------------------------------------------------------------
  Subroutine molecules_changedof(spc,newdof,new_origin)
    Integer, Intent(In)                :: spc,newdof
    Character(*), Intent(In), Optional :: new_origin

    If (Present(new_origin)) Then
      Call molecule_changedof(molecs%mo(spc),newdof,new_origin)
    Else 
      Call molecule_changedof(molecs%mo(spc),newdof)
    End If
  End Subroutine molecules_changedof

  !----------------------------------------------------------------------------
  ! Gets the coordinates of the nth atom type of species spc
  !----------------------------------------------------------------------------
  Type(VecType) Function molecules_getatomcoords(spc,nthAtom)
    Integer, Intent(In) :: spc,nthAtom
    molecules_getatomcoords = molecule_getatomcoords(molecs%mo(spc),nthAtom)
  End Function molecules_getatomcoords

  !----------------------------------------------------------------------------
  ! Gets the coordinates and elements of an entire molecule
  ! Requires: spc -- species number
  !           atomcoords -- array of coordinate vectors
  !           elements -- array of elemental symbols
  !----------------------------------------------------------------------------
  Subroutine molecules_getcoords(spc,atomcoords,elements)
    Integer, Intent(In)                         :: spc
    Type(VecType), Dimension(:), Intent(Out)    :: atomcoords
    Character(len=2), Dimension(:), Intent(Out) :: elements

    Integer                 :: atom,atype,natoms

    natoms = molecule_getnatoms(molecs%mo(spc))

    If (Size(atomcoords,1) < natoms) Then
      Write(0,'(2a,i4,a,2i4)') __FILE__,": ",__LINE__, &
          ' Supplied atomcoords array too small to contain all atoms', &
          Size(atomcoords,1),natoms
      Stop      
    End If
    If (Size(elements,1) < natoms) Then
      Write(0,'(2a,i4,a,2i4)') __FILE__,": ",__LINE__, &
          ' Supplied elements array too small to contain all atoms', &
          Size(elements,1),natoms
      Stop      
    End If

    !** get the coordinates
    Do atom = 1,natoms
      atomcoords(atom) = molecule_getatomcoords(molecs%mo(spc),atom)
      atype = molecule_getnthatom(molecs%mo(spc),atom)
      elements(atom) = atom_getsymbol(atype)
    End Do

  End Subroutine molecules_getcoords

  !----------------------------------------------------------------------------
  ! Assign pointer molecpntr to point at molecule number spc
  !----------------------------------------------------------------------------
  Subroutine molecules_getpointer(spc,molecpntr)
    Type(MolecularParams), Pointer :: molecpntr
    Integer, Intent(In) :: spc
    molecpntr => molecs%mo(spc)
  End Subroutine molecules_getpointer

  !----------------------------------------------------------------------------
  ! Checks whether a molecule exists with the given name 
  ! Works similar to _gettype
  !----------------------------------------------------------------------------
  Logical Function molecules_exists(name)
    Character(*), Intent(IN)         :: name
    Integer                          :: spc

    molecules_exists = .False. 
    Do spc = 1,molecules_getnsorbs()
      If (cleanstring(name) == cleanstring(molecule_getname(molecs%mo(spc)))) &
          Then
        molecules_exists = .True. 
        Return
      Endif
    End Do
  End Function molecules_exists

  !----------------------------------------------------------------------------
  ! Use the molecule name to get the molecule number
  ! Returns zero if it could not find out the molecule type
  !----------------------------------------------------------------------------
  Integer Function molecules_gettype(name)
    Character(*), Intent(IN)         :: name
    Integer                          :: spc

    molecules_gettype = 0
    Do spc = 1,molecules_getnsorbs()
      If (trim(name) == trim(molecule_getname(molecs%mo(spc)))) Then
        molecules_gettype = spc
        Return
      Endif
    End Do
  End Function molecules_gettype

  !----------------------------------------------------------------------------
  ! Get the atom symbol given the species number and atom number
  ! Requires: spc -- species number
  !           atom -- atom number
  !----------------------------------------------------------------------------
  Function molecules_atomsymbol(spc,atom)
    Character(len=2)          :: molecules_atomsymbol
    Integer, Intent(In)       :: spc,atom

    molecules_atomsymbol = atom_getsymbol(molecules_getatype(spc,atom))

  End Function molecules_atomsymbol

  !----------------------------------------------------------------------------
  ! Get the atom name given the species number and atom number
  ! Requires: spc -- species number
  !           atom -- atom number
  !----------------------------------------------------------------------------
  Function molecules_getaname(spc,atom)
    Character(len=strLen)          :: molecules_getaname
    Integer, Intent(In)            :: spc,atom

    molecules_getaname = atom_getname(molecules_getatype(spc,atom))

  End Function molecules_getaname

  !----------------------------------------------------------------------------
  ! Returns the number of atom types in the molecule number spc. Optionally
  ! returns a list of atom types in alist.
  !----------------------------------------------------------------------------
  Integer Function molecules_getnatomtypes(spc,alist,allTypes)
    Integer, Intent(In) :: spc
    Integer, Dimension(:), Intent(Out), Optional :: alist
    Logical, Optional, Intent(In) :: allTypes
    Logical :: all

    If (Present(allTypes)) Then
      all = allTypes
    Else
      all = .False.
    End If

    If (Present(alist)) Then
      molecules_getnatomtypes = molecule_getnatomtypes(molecs%mo(spc),alist,&
          all)
    Else
      molecules_getnatomtypes = molecule_getnatomtypes(molecs%mo(spc))
    End If
  End Function molecules_getnatomtypes

  !----------------------------------------------------------------------------
  ! Returns the atom charges for molecule number spc
  !----------------------------------------------------------------------------
  Subroutine molecules_getcharges(spc,clist)
    Integer, Intent(In) :: spc
    Real(Kind=RDbl), Dimension(:), Intent(Out) :: clist
    Integer :: natoms

    natoms = molecule_getnatoms(molecs%mo(spc))

    !** WARNING HACK HACK HACK
    !** Should call molecule to get the charges
    clist = Reshape(molecs%mo(spc)%atoms(1:natoms)%q0,(/natoms/))
  End Subroutine molecules_getcharges

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
  Subroutine molecules_expand(spc,offset,n,reconnect,fcell)
    Integer, Intent(In)                            :: spc
    Type(VecType), Intent(In)                      :: offset
    Integer, Intent(In)                            :: n
    Logical, Intent(In)                            :: reconnect   
    Type(Fundamental_Cell), Intent(In), Optional   :: fcell    

    If (Present(fcell)) Then
      Call molecule_expand(molecs%mo(spc),offset,n,reconnect,fcell)
    Else
      Call molecule_expand(molecs%mo(spc),offset,n,reconnect)
    End If

  End Subroutine molecules_expand

  !----------------------------------------------------------------------------
  ! Regenerates the connection list for molecule spc and then regenerates
  ! the intramolecular interaction lists
  ! Requires:  spc -- species number
  !            fcell -- optional fundamental cell if PBCs are to be used
  !----------------------------------------------------------------------------
  Subroutine molecules_reconnect(spc,fcell)
    Integer, Intent(In)                            :: spc
    Type(Fundamental_Cell), Intent(In), Optional   :: fcell    

    If (Present(fcell)) Then
      Call molecule_reconnect(molecs%mo(spc),fcell)
    Else
      Call molecule_reconnect(molecs%mo(spc))
    End If

  End Subroutine molecules_reconnect

  !----------------------------------------------
  ! Get the mass of the molecule
  !----------------------------------------------
  Real(kind=RDbl) Function molecules_getmass(spc)
    Integer, Intent(in)  :: spc
    molecules_getmass = molecule_getmass(molecs%mo(spc))
    Return
  End Function molecules_getmass

  !----------------------------------------------
  ! Get the net charge of a specific molecule
  !----------------------------------------------
  Real(kind=RDbl) Function molecules_netcharge(spc)
    Integer, Intent(in)  :: spc

    molecules_netcharge = molecule_netcharge(molecs%mo(spc))

  End Function molecules_netcharge

  !----------------------------------------------------------------------------
  ! Get the atom type given the atom number in the molecule
  !----------------------------------------------------------------------------
  Integer Function molecules_getatype(spc,atom)
    Integer, Intent(In)              :: spc,atom

    molecules_getatype = molecs%mo(spc)%atoms(atom)%atype
  End Function molecules_getatype

  !----------------------------------------------------------------------------
  ! Get the type of the generalized coordinates used
  ! to describe this molecule
  !----------------------------------------------------------------------------
!!$  Character(len=strLen) Function molecules_getgcmodeltype(spcno)
  Function molecules_getgcmodeltype(spcno)
    Character(len=strLen) :: molecules_getgcmodeltype
    Integer, Intent(in)  :: spcno
    
    molecules_getgcmodeltype = molecule_getgcmodeltype(molecs%mo(spcno))
    Return
  End Function molecules_getgcmodeltype

  !----------------------------------------------------------------------
  ! This routine returns the principal moments-of-inertia of species type
  ! "spcno". This is probably only useful for a rigid molecule
  !----------------------------------------------------------------------
  Function molecules_getmoi(spcno)
    Integer, Intent(in)  :: spcno
    Real(kind=RDbl), Dimension(3) :: molecules_getmoi
    
    molecules_getmoi = molecule_getmoi(molecs%mo(spcno))
  End Function molecules_getmoi

  !----------------------------------------------------------------------
  ! This routine returns the original transformation that took the 
  ! molecule file coordinates to the current coordinates in 
  !----------------------------------------------------------------------
  Subroutine molecules_origxform(spc,trans_xform,rot_xform)
    Integer, Intent(In)              :: spc
    Type(VecType), Intent(Out)       :: trans_xform
    Type(MatrixType), Intent(Out)    :: rot_xform
    
    Call molecule_origxform(molecs%mo(spc),trans_xform,rot_xform)
        
  End Subroutine molecules_origxform

  !-----------------------------------------------------------------------
  ! Returns TRUE if all the connections are complete, FALSE otherwise.
  !-----------------------------------------------------------------------
  Subroutine molecules_checkconnect(spc,allOk)
    Integer, Intent(In) :: spc
    Logical, Intent(Out) :: allOk

    !MDEBUG
    Write(0,*) __FILE__,__LINE__

    allOk = .True.

    !** Call the molecule_checkconnect function to check all the connections
    allOk = molecule_checkconnect(molecs%mo(spc))
  End Subroutine molecules_checkconnect

  !-----------------------------------------------------------------------
  ! Get the Reference coordinate definition of each molecule type.  The 
  ! original function was causing compiler faults. Changed to a subroutine
  ! Requires: spc -- species number
  !           com -- returned reference coordinates
  !-----------------------------------------------------------------------
  Subroutine molecules_getcomdefn(spc,com)
    Integer, Intent(In)                        :: spc
    Type(VecType), Dimension(:), Intent(InOut) :: com
    
    Call molecule_getcomdefn(molecs%mo(spc),com)

  End Subroutine molecules_getcomdefn

  !-----------------------------------------------------------------------
  ! Get the center of mass of a molecule using either the stored reference
  ! coordinates or, optionally, a new set of xyz coordinates.  Note that
  ! the ORDER OF THE PASSED COORDINATES MUST MATCH that in the molecule 
  ! structure.
  ! Requires: spc -- species number
  !           xyzcoords -- optional xyz coordinates for each atom
  !-----------------------------------------------------------------------
  Type(VecType) Function molecules_getcom(spc,xyzcoords)
    Integer, Intent(In)                                :: spc
    Type(VecType), Dimension(:), Intent(In), Optional  :: xyzcoords
    
    If (Present(xyzcoords)) Then
      molecules_getcom = molecule_getcom(molecs%mo(spc),xyzcoords)
    Else
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      molecules_getcom = molecule_getrefcom(molecs%mo(spc))
    End If

  End Function molecules_getcom

  !----------------------------------------------------------------------------
  ! Get the name of a specific species type
  !----------------------------------------------------------------------------
!!$  Character(len=strLen) Function molecules_name(spc)
  Function molecules_name(spc)
    Character(len=strLen) :: molecules_name
    Integer, Intent(IN)              :: spc
    
    molecules_name = molecule_getname(molecs%mo(spc))
  End Function molecules_name

  !----------------------------------------------------------
  ! This routine gets an array containing the inverse masses
  ! of each atom in a specified molecule type
  !----------------------------------------------------------
  Subroutine molecules_getinvmasses(type,massinv)
    Integer, Intent(In)                        :: type
    Real(kind=RDbl), Dimension(:), Intent(Out) :: massinv

    Integer                   :: a

    Do a = 1,molecules_getnatoms(type)
      massinv(a) = atom_invmass(molecules_getatype(type,a))
    End Do

  End Subroutine molecules_getinvmasses

  !----------------------------------------------------------
  ! This routine gets an array containing the masses
  ! of each atom in a specified molecule type
  !----------------------------------------------------------
  Subroutine molecules_getmasses(type,mass)
    Integer, Intent(In)                        :: type
    Real(kind=RDbl), Dimension(:), Intent(Out) :: mass

    Integer                   :: a

    Do a = 1,molecules_getnatoms(type)
      mass(a) = atom_getmass(molecules_getatype(type,a))
    End Do

  End Subroutine molecules_getmasses

  !----------------------------------------------------------
  ! This routine gets the mass of species-spc atom-a 
  !----------------------------------------------------------
  Real(kind=RDbl) Function molecules_AtomMass(spc, a)
    Integer, Intent(In)                        :: spc, a
    molecules_AtomMass=molecule_AtomMass(molecs%mo(spc),a)
  End Function molecules_AtomMass

  !--------------------------------------------------------------------
  ! This routine gets a line in one of the molecule files whose first
  ! string matches a given tag.  Function returns true if successful.
  ! Requires:  spc -- species number
  !            tag -- tag to search for
  !            line -- full line containing the tag
  !            unit -- optional unit number of open file
  ! NOTE: if the optional unit number is not passed, the function will
  ! close the file when it is done.
  !--------------------------------------------------------------------
  Logical Function molecules_getlinewtag(spc,tag,line,unit)
    Integer, Intent(In)                        :: spc
    Character(*), Intent(In)                   :: tag
    Character(len=lstrLen), Intent(Out)        :: line
    Integer, Intent(Out), Optional             :: unit

    Integer                                  :: lineno,unitno,nfields
    Character(len=lstrLen)                   :: filewpath
    Character(len=strLen), Dimension(strLen) :: fields

    molecules_getlinewtag = .False.

    !** Open the molecule file
    Write(filewpath,'(2a)') Trim(getpath('MOLSDIR')),&
        Trim(molecules_getfilename(spc))
    unitno = file_open(filewpath,110)
    
    !** Find the tag in the file
    Rewind(unit=unitno)
    Do While (.Not. molecules_getlinewtag)
      lineno = filesrchstr(unitno,tag,line,.False.)
      If (lineno > 0) Then
        nfields = split(line,fields)
  
        !** Make sure the tag is the first string present on line
        If (Toupper(Trim(fields(1)(1:2))) == Toupper(Trim(tag(1:2)))) Then
          molecules_getlinewtag = .True.
        End If
      End If
    End Do

    If (Present(unit)) Then
      unit = unitno
    Else
      Close(unit=unitno)
    End If

  End Function molecules_getlinewtag

  !----------------------------------------------------------------------------
  ! Checks to see if the molecules have been initialized
  !----------------------------------------------------------------------------
  Logical Function molecules_checkinit()
    If (Associated(molecs%mo)) Then
      molecules_checkinit = .TRUE.
    Else
      molecules_checkinit = .FALSE.
    End If
    Return
  End Function molecules_checkinit

  !----------------------------------------------------------------------------
  ! Given a molecule number, it dumps the molecule coordinates to an xyz 
  ! type file name filename
  !----------------------------------------------------------------------------
  Subroutine molecules_dump2xyz(molec,filename,comment)
    Integer, Intent(In)         :: molec
    Character(*), Intent(In)    :: filename
    Character(*), Optional      :: comment

    If (Present(comment)) Then
      Call molecule_molec2xyz(molecs%mo(molec),filename,comment)
    Else
      Call molecule_molec2xyz(molecs%mo(molec),filename)
    End If

  End Subroutine molecules_dump2xyz

  !-------------------------------------------------------
  ! Dump the fields of the "molecules" structure
  !-------------------------------------------------------
  Subroutine molecules_display(unit)
    Integer              :: i,unit
    
    Write(unit,'(a)') dashedline
    Write(unit,'(a)') 'The MOLECULES structure:'
    Write(unit,'(2x,a,i5)') 'No. of molecule types: ', molecs%nmolec_types
    Do i = 1,molecs%nmolec_types
      Write(unit,'(a)') dashedline
      Call molecule_display(molecs%mo(i),unit,2)
    End Do
    Write(unit,'(a)') dashedline
  End Subroutine molecules_display

  !----------------------------------------------------------------------------
  ! Cleanup the molecules allocations
  !----------------------------------------------------------------------------
  Subroutine molecules_cleanup()
    Integer :: error

    Deallocate(molecs%mo,STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' unable to deallocate molecs%mo'
      Stop
    End If
  End Subroutine molecules_cleanup

End Module Molecules
