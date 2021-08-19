!------------------------------------------------------------------------------
! This module handles molecule connectivity.  It provides structures and their
! associated routines to determine the type of connection and the connectivity
! of a structure.  
!
! Two important types of basic data structures are provided here.  One stores
! information on a TYPE of connection (ConnectionParams) between two atoms.
! The other type stores the atom numbers and lengths of the connected atoms 
! (AtomConnections).  
!
! The first type (ConnectionParams) are wrapped into a 2D array (ConnectMatrix)
! that is indexed by atom-type.  This array then allows one to access basic
! information on the bonds between two atom types.
! Important routines:
!   connects_initarray -- initialize the ConnectMatrix shape
!   connects_initconnection -- initialize a single atype-atype entry
!
! The second type (AtomConnections) are usually bundled together in an array 
! with one index for each atom in a structure.  Collectively, they then hold
! all the information about the connectivity of a structure.  In a set, these
! data structure can be recursively searched to determine paths between two
! atoms or any other connectivity information.
! Important routines:
!   connects_initmaxcon -- initialize single AtomConnections data structure
!   connects_connectatoms -- given the array, make an atom-atom connection 
!   connects_allneighbors -- get a list of all atom number within given ncon
!   connects_makeTypeList -- returns lists of connected atoms with given atypes
!   connects_makeXList -- lists of all connectivity sets of a given length
!
! Needed Improvements
! 1) May fail for ring structures, fix the commented out routines
! 2) remove connects_regenerateAll after making sure it's unneeded
! 3) consider splitting the connection type routines into another module
!------------------------------------------------------------------------------

Module connects

  Use defaults, Only: RDbl, lstrLen, strLen
  Use utils, Only: findint, findgt, toupper, allocerrdisplay, deallocerrdisplay
  Use file, Only: file_getunit
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use fundcell, Only: Fundamental_Cell, fundcell_snglMinImage

  Implicit None
  Save

  Private
  Public :: MAX_CONNECTIONS, AtomConnections, AtomList, ConnectionParams, &
      ConnectMatrix, ConnectPaths, connects_initarray, &
      connects_initconnection, connects_initmaxcon, connects_connectatoms, &
      connects_initcopy, connects_nconnected, connects_nthcon, &
      connects_allneighbors, connects_atomsxaway, &
      connects_makeXList, connects_makeTypeList, &
      connects_generateAll, connects_regenerateAll, &
      connects_display, connects_displayall, connects_cleanup

  Integer, Parameter :: MAX_CONNECTIONS = 6
  Integer, Parameter :: MAX_NEIGHBORS = 20
  Real(kind=RDbl), Parameter :: CONNECT_TOL = 0.01_RDbl
  Character(strLen), Parameter :: shellcore_idstring = "SHELLCORE"

  !** This type stores the information about an individual connection
  !** type, including bond length, connection tolerance, the min and
  !** max bond length based on the tolerance. Init indicates if the
  !** bond type is initialized, and shellCore indicates a shell-core
  !** bond type
  Type ConnectionParams
    Real(kind=RDbl) :: length
    Real(kind=RDbl) :: tol
    Real(kind=RDbl) :: maxlength
    Real(kind=RDbl) :: minlength
    Logical         :: init
    Logical         :: shellCore
  End Type ConnectionParams

  !** This data type is used to store connect types indexed by ATOM TYPE
  Type ConnectMatrix
    Type(ConnectionParams), Dimension(:,:), Pointer :: bonds
  End Type ConnectMatrix

  !** This is a single connection to an individual atom. atom is 
  !** the atom number, and length is the equilibrium distance 
  !** between the two. length is NOT updated as coordinates in the
  !** simulation change.
  Type AtomList
    Integer         :: atom
    Real(kind=RDbl) :: length
  End Type AtomList

  !** This is the list of connections, usually indexed by
  !** atom number. nconnect is the number of atoms connected to
  !** this atom, and the clist contains the connected atoms and
  !** their distances.
  Type AtomConnections
    Integer                               :: nconnect
    Type(AtomList), Dimension(:), Pointer :: clist
  End Type AtomConnections

  !** This is a list of paths, 'path' gives atom number and 
  !** is indexed: pathno,atomno
  Type ConnectPaths
    Integer                          :: npaths
    Integer, Dimension(:,:), Pointer :: path
  End Type ConnectPaths

Contains
  !----------------------------------------------------------------------------
  ! Initialize the size of the atype-indexed connections array
  ! Requires:  conmatrix -- Connection Matrix structure to initialize
  !            natypes -- number of atom types
  !----------------------------------------------------------------------------
  Subroutine connects_initarray(conmatrix,natypes)
    Type(ConnectMatrix), Intent(InOut) :: conmatrix
    Integer, Intent(In)                :: natypes

    Integer            :: error, i, j

    !** Size the matrix based on the number of atom types. This way,
    !** we can look up bond pairs based solely upon their atom types.
    Allocate(conmatrix%bonds(natypes,natypes),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Nullify the pointers in the bonds matrix. Later on, 
    !** we will check to see if the pointers are associated
    !** for building connections
    Do i = 1,natypes
      Do j = 1,natypes
        conmatrix%bonds(i,j)%init = .False.
      End Do
    End Do

  End Subroutine connects_initarray

  !----------------------------------------------------------------------------
  ! Stores the connection in the connection matrix given the atypes, 
  ! length, and (optionally) the tolerance
  ! Requires:  conmatrix -- Connection Matrix structure 
  !            atype1 -- 1st atom type
  !            atype2 -- 2nd atom type
  !            length -- connection length
  !            tolerance -- length tolerance
  !            extraInfo -- string containing extra information
  !----------------------------------------------------------------------------
  Subroutine connects_initconnection(conmatrix,atype1,atype2,length, &
      tolerance,extraInfo)
    Type(ConnectMatrix), Intent(InOut)    :: conmatrix
    Integer, Intent(In)                   :: atype1, atype2
    Real(Kind=RDbl), Intent(In)           :: length
    Real(Kind=RDbl), Intent(In), Optional :: tolerance
    Character(*), Intent(In), Optional    :: extraInfo

    Real(Kind=RDbl)             :: tol

    !** If tolerance is passed, use it. Otherwise, default tolerance is 0.
    If (Present(tolerance)) Then
      tol = tolerance
    Else
      tol = 0.0_RDbl
    End If

    !** If extraInfo is passed, check to see what it contains
    conmatrix%bonds(atype1,atype2)%shellCore = .False.
    If (Present(extraInfo)) Then
      Select Case (toupper(extraInfo))

      Case (shellcore_idstring)
        conmatrix%bonds(atype1,atype2)%shellCore = .True.
      Case Default
        Write(0,'(2a,i6,2a)') __FILE__,":",__LINE__, &
            " Ignoring unsupported flag ",extraInfo
      End Select
    End If

    !** Add the pair to the matrix
    conmatrix%bonds(atype1,atype2)%init = .True.
    conmatrix%bonds(atype2,atype1)%init = .True.

    !** Store the length, tolerance, and max and min lengths based on tol
    conmatrix%bonds(atype1,atype2)%length = length
    conmatrix%bonds(atype1,atype2)%tol = tol
    conmatrix%bonds(atype1,atype2)%maxlength = length + tol
    conmatrix%bonds(atype1,atype2)%minlength = length - tol

    !** Make sure the complementary pair points to the same data
    conmatrix%bonds(atype2,atype1)%length = length
    conmatrix%bonds(atype2,atype1)%tol = tol
    conmatrix%bonds(atype2,atype1)%maxlength = length + tol
    conmatrix%bonds(atype2,atype1)%minlength = length - tol
    conmatrix%bonds(atype2,atype1)%shellCore = &
        conmatrix%bonds(atype1,atype2)%shellCore

#ifdef DEBUG
    Write(*,*) atype1,atype2,' length: ',conmatrix%bonds(atype1,atype2)%length
    Write(*,*) atype1,atype2,' tolerance: ',conmatrix%bonds(atype1,atype2)%tol
    Write(*,'(2x,2i3,a,2f8.5)') atype1,atype2,' Max, Min connection lengths :', &
        conmatrix%bonds(atype1,atype2)%maxlength, &
        conmatrix%bonds(atype1,atype2)%minlength
#endif

  End Subroutine connects_initconnection

  !----------------------------------------------------------------------------
  ! Initialize the number of connections for each entry
  ! Requires:  connect -- single Atom Connections data structure
  !            maxcon -- maximum number of connections
  !----------------------------------------------------------------------------
  Subroutine connects_initmaxcon(connect,maxcon)
    Type(AtomConnections), Intent(InOut) :: connect
    Integer, Intent(In)                  :: maxcon

    Integer             :: error

    Allocate(connect%clist(maxcon),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"connect%clist")

    connect%clist(1:maxcon)%atom = 0
    connect%clist(1:maxcon)%length = 0.0_RDbl
    connect%nconnect = 0

  End Subroutine connects_initmaxcon

  !-------------------------------------------------------------------
  ! Given the full array of atom connectivity data structure, place a
  ! connection between two given atom numbers.  Note that it makes two
  ! additions, one under each atom's connectivity storage.
  ! Requires:  connect -- Array of Atom Connections for the structure
  !            atom1 -- 1st atom number
  !            atom2 -- 2nd atom number
  !            length -- distance between the two atoms
  !-------------------------------------------------------------------
  Subroutine connects_connectatoms(connect,atom1,atom2,length)
    Type(AtomConnections), Dimension(:), Intent(InOut) :: connect
    Integer, Intent(In)                                :: atom1,atom2
    Real(kind=RDbl), Intent(In)                        :: length

    Character(len=lstrLen)               :: string

    If (.Not. connects_iscon(connect(atom1),atom2)) Then
      If (connect(atom1)%nconnect == Size(connect(atom1)%clist)) Then
        Write(0,'(2a,i4,a,i6)') __FILE__,": ",__LINE__, &
            " : Max connections reached ",Size(connect(atom1)%clist)
        Write(0,'(2x,a)') 'Error information'
        Write(0,'(2x,a,i5,a,i5)') &
            'Attempted to connect atoms ',atom1,' and ',atom2
        Write(0,'(2x,a,f10.3)') 'Connection length :',length
        string = connects_display(connect(atom1))
        Write(0,'(2x,a,i5,2a)') 'Atom ',atom1,' connections: ',Trim(string)
        string = connects_display(connect(atom2))
        Write(0,'(2x,a,i5,2a)') 'Atom ',atom2,' connections: ',Trim(string)
        Stop
      End If
      
      connect(atom1)%nconnect = connect(atom1)%nconnect+1
      connect(atom1)%clist(connect(atom1)%nconnect)%atom  = atom2
      connect(atom1)%clist(connect(atom1)%nconnect)%length  = length
    End If
    
    If (.Not. connects_iscon(connect(atom2),atom1)) Then
      If (connect(atom2)%nconnect == Size(connect(atom2)%clist)) Then
        Write(0,'(2a,i4,a,i6)') __FILE__,": ",__LINE__, &
            " : Max connections reached ",Size(connect(atom2)%clist)
        Write(0,'(2x,a)') 'Error information'
        Write(0,'(2x,a,i5,a,i5)') &
            'Attempted to connect atoms ',atom1,' and ',atom2
        Write(0,'(2x,a,f10.3)') 'Connection length :',length
        string = connects_display(connect(atom1))
        Write(0,'(2x,a,i5,2a)') 'Atom ',atom1,' connections: ',Trim(string)
        string = connects_display(connect(atom2))
        Write(0,'(2x,a,i5,2a)') 'Atom ',atom2,' connections: ',Trim(string)
        Write(0,'(2x,a,i5)') 'Number of existing connections: ', &
            connect(atom1)%nconnect
        Stop
      End If

      connect(atom2)%nconnect = connect(atom2)%nconnect+1
      connect(atom2)%clist(connect(atom2)%nconnect)%atom  = atom1
      connect(atom2)%clist(connect(atom2)%nconnect)%length  = length
    End If
    
  End Subroutine connects_connectatoms

  !---------------------------------------------------------------------------
  ! Initialize and make a copy of a single atom connection structure
  ! Requires:  image -- new single Atom Connections data structure
  !            old -- single Atom Connections data structure to copy
  !---------------------------------------------------------------------------
  Subroutine connects_initcopy(image,old)
    Type(AtomConnections), Intent(Out) :: image
    Type(AtomConnections), Intent(In)  :: old

    Integer      :: i

    Call connects_initmaxcon(image,Size(old%clist))

    image%nconnect = old%nconnect
    Do i = 1,old%nconnect
      image%clist(i)%atom = old%clist(i)%atom
      image%clist(i)%length = old%clist(i)%length
    End Do

  End Subroutine connects_initcopy

  !--------------------------------------------------------------------------
  ! Checks to see if an atom is connected to another.  
  ! Requires:  connect -- atom connections structure for one atom
  !            atom -- atom number to look for in the connection list
  !--------------------------------------------------------------------------
  Logical Function connects_iscon(connect,atom)
    Type(AtomConnections), Intent(In) :: connect
    Integer, Intent(In)               :: atom

    Integer :: i

    Do i = 1,connect%nconnect
      If (connect%clist(i)%atom == atom) Then
        connects_iscon = .True.
        Return 
      End If
    End Do

    connects_iscon = .False.

  End Function connects_iscon

  !--------------------------------------------------------------
  ! Return the number of atoms connected to atom
  ! Requires:  connect -- single Atom Connections data structure
  !--------------------------------------------------------------
  Integer Function connects_nconnected(connect)
    Type(AtomConnections), Intent(In) :: connect

    connects_nconnected = connect%nconnect

  End Function connects_nconnected

  !--------------------------------------------------------------------
  ! Returns the AtomList structure for the nth connection
  ! Requires:  connect -- single Atom Connections data structure
  !            n -- connection number for which to return info
  !            output -- returned AtomList structure
  !--------------------------------------------------------------------
  Subroutine connects_nthcon(connect,n,output)
    Type(AtomConnections), Intent(In) :: connect
    Integer, Intent(In)               :: n
    Type(AtomList), Intent(Out)       :: output

    output = connect%clist(n)
    
  End Subroutine connects_nthcon

  !----------------------------------------------------------------------------
  ! Returns a list of all the neighbors of atom within a given number of 
  ! connections xaway away and their distance away.  This routine is very 
  ! similar to atomsxaway execpt that it will not keep track of the exact 
  ! path, but only all the atoms within the given natom cutoff.
  ! Requires:  connect -- Array of Atom Connections for the structure
  !            atom -- atom number at one end of chain
  !            xaway -- maximum number of atoms away to include in output
  !            conlist -- a list of atom numbers connected within xaway atoms
  !----------------------------------------------------------------------------
  Integer Function connects_allneighbors(connect,atom,conlist,xaway)
    Type(AtomConnections), Intent(In), Dimension(:) :: connect
    Integer, Intent(In)                             :: atom, xaway
    Integer, Dimension(:), Intent(Out)              :: conlist

    Integer                             :: i, npals
    Integer, Dimension(Size(conlist,1)) :: pals

    !** Init the npals counter
    npals = 0
    
    !** Call the recursive function
    i = connects_genneighbors(connect,atom,atom,xaway,pals,npals)

    !** Copy pals to conlist
    conlist(1:npals) = pals(1:npals)

    !** Return the number of neighbors.
    connects_allneighbors = npals

  End Function connects_allneighbors

  !----------------------------------------------------------------------------
  ! Recursive part of _allneighbors.  Given an atom number, it steps 
  ! recursively along all of its listed connections and accumulates unique
  ! atom numbers within a given connectivity number.
  ! Requires:  connect -- Array of Atom Connections for the structure
  !            atom -- current atom to branch from
  !            last -- last atom number visited
  !            xaway -- maximum number of atoms away to include in output
  !            pals -- a list of atom numbers connected within xaway atoms
  !            npals -- number of neighbors found so far
  !----------------------------------------------------------------------------
  Recursive Integer Function connects_genneighbors(connect,atom,last, &
      xaway,pals,npals) Result(nneighbors)
    Type(AtomConnections), Intent(In), Dimension(:) :: connect
    Integer, Intent(In)                             :: atom, last, xaway
    Integer, Dimension(:), Intent(InOut)            :: pals
    Integer, Intent(InOut)                          :: npals

    Integer            :: next, n, i

    !** Check for the ending condition, i.e., xaway = 0
    If (xaway == 0) Then
      !** We've reached the furthest we want to go
      nneighbors = npals
      Return
    End If

    !** Go through the list
    Do i = 1, connect(atom)%nconnect
      !** next atom to check out
      next = connect(atom)%clist(i)%atom
      If (next /= last) Then
        !** It's a new atom, add it to the list
        npals = npals+1
        pals(npals) = next
        n = connects_genneighbors(connect,next,atom,xaway-1,pals,npals)
      End If
    End Do
      
    nneighbors = npals

  End Function connects_genneighbors

  !------------------------------------------------------------------------
  ! Find all atoms x connections away, and return the paths to these atoms
  ! Requires:  connect -- Array of Atom Connections for the structure
  !            atom -- current atom to branch from
  !            xaway -- maximum number of atoms away to include in output
  !            pals -- a list of atom numbers connected AT xaway atoms
  !            paths -- optional returned set of paths
  !------------------------------------------------------------------------
  Integer Function connects_atomsxaway(connect,atom,xaway,pals,paths)
    Type(AtomConnections), Intent(In), Dimension(:) :: connect
    Integer, Intent(In)                             :: xaway,atom
    Integer, Dimension(:)                           :: pals
    Type(ConnectPaths), Intent(InOut), Optional     :: paths

    Integer               :: npals
    Integer, Dimension(1) :: lastpath

    !** Call the recursive function to find the closest atoms to atom
    npals = 0

    If (Present(paths)) Then
      !** Initialize the last path by adding the first atom to it
      lastpath(1) = atom

      !** Initialize the paths structure
      paths%path = 0
      paths%npaths = 0
      npals = connects_getatomsxawayp(connect,atom,atom,xaway,pals,npals, &
          lastpath,paths)
    Else
      npals = connects_getatomsxaway(connect,atom,atom,xaway,pals,npals)
    End If

    connects_atomsxaway = npals

  End Function connects_atomsxaway

  !------------------------------------------------------------------------
  ! Recurse through the tree to find connections at the atom number that
  ! are at a given number of connections away.
  ! Requires:  connect -- Array of Atom Connections for the structure
  !            atom -- current atom to branch from
  !            lastatom -- last atom visited
  !            xaway -- maximum number of atoms away to include in output
  !            pals -- a list of atom numbers connected AT xaway atoms
  !            npals -- number of pals
  !------------------------------------------------------------------------
  Recursive Function connects_getatomsxaway(connect,atom,lastatom,xaway, &
      pals,npals) Result(nfound)
    Type(AtomConnections), Intent(In), Dimension(:) :: connect
    Integer, Intent(In)                             :: atom, lastatom, xaway
    Integer, Dimension(:), Intent(Out)              :: pals
    Integer, Intent(InOut)                          :: npals

    Integer :: i, nfound, nextatom

    If (xaway == 0) Then    !** We're done
      npals = npals + 1
      pals(npals) = atom
      nfound = npals
      Return
    End If

    !** search the tree
    Do i = 1, connect(atom)%nconnect
      nextatom = connect(atom)%clist(i)%atom
      If (nextatom /= lastatom) Then
        nfound = connects_getatomsxaway(connect,nextatom,atom, &
            xaway-1,pals,npals)
      End If
    End Do

    nfound = npals

  End Function connects_getatomsxaway

  !------------------------------------------------------------------------
  ! Recurse through the tree to find connections at the atom number that
  ! are at a given number of connections away.  Same as _getatomsxaway, 
  ! but this routine also returns the paths.
  ! Requires:  connect -- Array of Atom Connections for the structure
  !            atom -- current atom to branch from
  !            lastatom -- last atom visited
  !            xaway -- maximum number of atoms away to include in output
  !            pals -- a list of atom numbers connected AT xaway atoms
  !            npals -- number of pals
  !            paths -- returned set of paths
  !------------------------------------------------------------------------
  Recursive Function connects_getatomsxawayp(connect,atom,lastatom,xaway, &
      pals,npals,lastpath,paths) Result(nfound)
    Type(AtomConnections), Intent(in), Dimension(:) :: connect
    Integer, Intent(In)                             :: atom, lastatom, xaway
    Integer, Dimension(:), Intent(Out)              :: pals
    Integer, Intent(InOut)                          :: npals
    Integer, Dimension(:), Intent(In)               :: lastpath
    Type(ConnectPaths), Intent(InOut)               :: paths

    Integer :: i,nfound,nextatom

    If (xaway == 0) Then   !** We're done
      npals = npals + 1
      pals(npals) = atom
      nfound = npals

      !** Add the path from the original atom to the final to paths array
      paths%npaths = paths%npaths+1
      paths%path(paths%npaths,1:Size(lastpath,1)) = lastpath
      Return
    End If

    !** search the tree
    Do i = 1, connect(atom)%nconnect
      nextatom = connect(atom)%clist(i)%atom
      If (nextatom /= lastatom) Then
        !** When the variables are passed to the recursive function
        !** the new atom, nextatom, is added to the lastpath array.
        !** This way, when it returns here we still have the original
        !** path
        nfound = connects_getatomsxawayp(connect,nextatom,atom, &
            xaway-1,pals,npals,(/lastpath,nextatom/),paths)
      End If
    End Do

    nfound = npals

  End Function connects_getatomsxawayp

  !----------------------------------------------------------------------------
  ! Generates the connections for an entire structure given the connection
  ! matrix.  Uses PBCs if the fundamental cell is passed.
  ! Requires:  connect -- Array of Atom Connections for the structure
  !            atypes -- a list of atom types for each atom in structure
  !            coords -- position vector for each atom
  !            conmatrix -- Connection Matrix structure 
  !            fcell -- optional fundamental cell if PBCs are to be used
  !----------------------------------------------------------------------------
  Subroutine connects_generateAll(connect,atypes,coords,conmatrix,fcell)
    Type(AtomConnections), Dimension(:), Intent(InOut)  :: connect
    Integer, Dimension(:), Intent(In)                   :: atypes
    Type(VecType), Dimension(:), Intent(In)             :: coords
    Type(ConnectMatrix), Intent(In)                     :: conmatrix
    Type(Fundamental_Cell), Intent(In), Optional        :: fcell    

    Integer           :: i, j
    Real(kind=RDbl)   :: sep
    Type(VecType)     :: sepvec

    Do i = 1,Size(coords,1)
      Do j = i,Size(coords,1)

        !** Check to see if the two atom types should be
        !** connected in the connection matrix
        If (conmatrix%bonds(atypes(i),atypes(j))%init) Then

          !** Get the separation vector with or without PBCs
          If (Present(fcell)) Then          
            Call fundcell_snglMinImage(fcell,coords(i),coords(j),sepvec)
          Else
            sepvec = coords(i) - coords(j)
          End If

          sep = mag(sepvec)
          
          !** Check to see if the atoms are within the given length
          If ((sep <= conmatrix%bonds(atypes(i),atypes(j))%maxlength).And. &
              (sep >= conmatrix%bonds(atypes(j),atypes(i))%minlength)) Then
            !** Connect the two atoms
            Call connects_connectatoms(connect,i,j,mag(coords(i)-coords(j)))
          End If
        End If

      End Do
    End Do
    
  End Subroutine connects_generateAll

  !----------------------------------------------------------------------------
  ! Regenerates the connections for a given molecule based on the parameters
  ! already stored in the connect structure.  Uses PBCs if the fundamental 
  ! cell is passed.
  !
  ! We've broken the OOP paradigm here by including minimum image calculations
  ! but I can't figure out a better way to do it. --Marty
  !
  ! Requires:  connect -- Array of Atom Connections for the structure
  !            atypes -- a list of atom types for each atom in structure
  !            coords -- position vector for each atom
  !            conmatrix -- Connection Matrix structure 
  !            fcell -- optional fundamental cell if PBCs are to be used
  !----------------------------------------------------------------------------
  Subroutine connects_regenerateAll(connect,atypes,coords,conmatrix,fcell)
    Type(AtomConnections), Dimension(:), Intent(InOut)  :: connect
    Integer, Dimension(:), Intent(In)                   :: atypes
    Type(VecType), Dimension(:)                         :: coords
    Type(ConnectMatrix), Intent(In)                     :: conmatrix
    Type(Fundamental_Cell), Intent(In), Optional        :: fcell    

    Integer           :: i, j, k
    Real(kind=RDbl)   :: sep
    Type(VecType)     :: sepvec

    !** Loop through the cooridinates to find the bonds
    Do i = 1, Size(coords,1)
      Do j = i, Size(coords,1)

        !** Check to see if the two atom types should be
        !** connected in the connection matrix
        If (conmatrix%bonds(atypes(i),atypes(j))%init) Then
          
          !** Get the separation vector with or without PBCs
          If (Present(fcell)) Then          
            Call fundcell_snglMinImage(fcell,coords(i),coords(j),sepvec)
          Else
            sepvec = coords(i) - coords(j)
          End If

          !** Calculate the separation magnitude
          sep = mag(sepvec)

          !** Check to see if the atoms are within the given length
          If ((sep <= conmatrix%bonds(atypes(i),atypes(j))%maxlength).And. &
              (sep >= conmatrix%bonds(atypes(j),atypes(i))%minlength)) Then

            !** Connect the two atoms
            Call connects_connectatoms(connect,i,j,sep)
          End If
        End If
      End Do
    End Do

  End Subroutine connects_regenerateAll

  !----------------------------------------------------------------------------
  ! Given the connections array, this routine returns a list of all atoms
  ! that are xaway and their paths in between.
  ! Requires:  connect -- Array of Atom Connections for the structure
  !            xaway -- the number of atoms distant to be returned
  !            conlist -- a list of the paths, for example:
  !   With xaway = 2, the routine returns conlist(path,1:xaway).
  !   Here path is the index of the path and the numbers in the entries 
  !   1:xaway are the atom numbers that make up the entire path.
  !----------------------------------------------------------------------------
  Integer Function connects_makeXList(connect,xaway,conlist)
    Type(AtomConnections), Dimension(:), Intent(In) :: connect
    Integer, Intent(In)                             :: xaway
    Integer, Dimension(:,:), Intent(Out)            :: conlist

    Integer                             :: i, n, error
    Integer                             :: nmatch, nxaway, oldmatch
    Integer, Dimension(Size(connect,1)) :: atomlist, all
    Type(ConnectPaths)                  :: paths

    !** Initialize the counter
    nmatch = 0

    !** Allocate temporary storage for the paths
    Allocate(paths%path(Size(connect,1),xaway+1),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'paths')    

    !** Loop through the list to find the desired atom connections
    Do i = 1, Size(connect,1)

      !** First find atoms xaway
      nxaway = connects_atomsxaway(connect,i,xaway,atomlist,paths)

      !** Remove any paths that are invalid, i.e., paths whose last
      !** atom is > i.
      n = findgt(paths%path(1:nxaway,xaway+1),i,all)

      !** If there are no valid paths, then cycle
      If (n == 0) Cycle

      !** Copy the returned paths into conlist
      oldmatch = nmatch
      nmatch = nmatch + n

      !** If we're looking for atoms xaway, the TOTAL number of atoms
      !** is xaway+1. For example, for the connection list of atoms 2
      !** places away, the list is atom1 atom2 atom3
      conlist(oldmatch+1:nmatch,1:xaway+1)=paths%path(all(1:n),1:xaway+1)

    End Do

    !** Free the temporary path memory
    Deallocate(paths%path,stat=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'paths')    

    connects_makeXList = nmatch

  End Function connects_makeXList

  !-----------------------------------------------------------------------------
  ! Given the array of atom connections, and a list of atypes to match, this 
  ! subroutine returns an array of lists containing connected atoms in the 
  ! structure that match the given atom types.
  ! 
  ! This routine is useful for generating, say, a list of atom triplets for 
  ! bond bending interactions. It assumes that the type list is SYMMETRIC, 
  ! unless symmFlag is passed as .False.
  ! Requires:  connect -- Array of Atom Connections data structures
  !            atypes -- a list of all atom types in structure
  !            typeSpec -- an ordered list of atom types to match
  !            conlist -- output array of atom nums in connection (listno,index)
  !            symmFlag -- optional specification atype sequence symmetry
  !-----------------------------------------------------------------------------
  Integer Function connects_makeTypeList(connect,atypes,typeSpec,conlist,&
      symmFlag)
    Type(AtomConnections), Dimension(:), Intent(In) :: connect
    Integer, Dimension(:), Intent(In)               :: atypes
    Integer, Dimension(:), Intent(In)               :: typeSpec
    Integer, Dimension(:,:), Intent(Out)            :: conlist
    Logical, Intent(In), Optional                   :: symmFlag

    Integer        :: xaway, i, j, ntypes, k
    Integer        :: nmatch, ngroups
    Logical        :: add, symmetric
    Integer, Dimension(Size(conlist,1),Size(conlist,2)) :: possibilities

    !** Set the symmetric flag
    symmetric = .True.
    If (Present(symmFlag)) symmetric = symmFlag

    !** Initialize the number of matches
    nmatch = 0 

    !** The size of the atom neighbors. This is the number of connections
    !** away to look
    xaway = Size(typeSpec,1)-1
    ntypes = Size(typeSpec,1)

    !** warn if the number of connections is greater than 3, we
    !** haven't checked its validity yet
    If (xaway > 3) Then
      Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
          " Warning, this has not been verified for more than 3 atom groups"
    End If

    !** Make a list of all the possibilities
    ngroups = connects_makeXList(connect,xaway,possibilities)

    !** Loop through all the groups and see which match
    Do i = 1, ngroups

      add = .False.

      !** Check the match in the foward direction (typeSpec(1) typeSpec(2) ...)
      If (atypes(possibilities(i,1)) == typeSpec(1)) Then
        !** Ok, first one matches
        add = .True.

        !** Loop through the rest to check if they match
        Do j = 2, ntypes
          If (atypes(possibilities(i,j)) /= typeSpec(j)) Then
            add = .False.
            Exit
          End If
        End Do

      !** Check the match in the reverse (typeSpec(end) typeSpec(end-1) ...)
      Else If (atypes(possibilities(i,1)) == typeSpec(ntypes)) Then
        If (.Not.symmetric) Cycle

        !** Ok, last matches
        add = .True.

        !** Start the second counter that goes from 1 to ntypes
        k = 1

        !** Loop through the rest to check if they match
        Do j = ntypes-1, 1, -1
          k = k+1
          If (atypes(possibilities(i,k)) /= typeSpec(j)) Then
            add = .False.
            Exit
          End If
        End Do

        !** Flip it around
        possibilities(i,1:ntypes) = possibilities(i,ntypes:1:-1)

      Else 
        !** Doesn't match either first or last.
        Cycle

      End If

      !** If add is still true, then add it to the conlist array
      If (add) Then
        nmatch = nmatch + 1
        conlist(nmatch,1:ntypes) = possibilities(i,1:ntypes)
      End If

    End Do

    !** Return the number of matches 
    connects_makeTypeList = nmatch

  End Function connects_makeTypeList

!!$  !--------------------------------------------------
!!$  ! Driver for the recursive getpath function
!!$  !--------------------------------------------------
!!$  Integer Function connects_getpathdriver(connect,atom1,atom2,path)
!!$    Integer, Intent(In) :: atom1, atom2
!!$    Type(AtomConnections), Dimension(:) :: connect
!!$    Integer, Dimension(:) :: path
!!$    Integer :: pathlen
!!$
!!$    pathlen = 0
!!$
!!$    pathlen=connects_getpath(connect,atom1,atom2,path,pathlen)
!!$
!!$    connects_getpathdriver = pathlen
!!$
!!$  End Function connects_getpathdriver
!!$
!!$  !---------------------------------------------------------------------------
!!$  ! Return the atoms between atom1 and atom2. If there are two paths, 
!!$  ! it returns the shortest UNLESS the flag long is passed. 
!!$  ! It does not currently check to see if it has already passed a given 
!!$  ! atom, so it's a BAD IDEA to use this for ringed molecules.
!!$  !---------------------------------------------------------------------------
!!$  Recursive Function connects_getpath(connect,atom1,atom2,path,pathlen) &
!!$      Result(pathl)
!!$    Type(AtomConnections), Dimension(:)  :: connect
!!$    Integer, Intent(In)                  :: atom1,atom2
!!$    Integer, Dimension(:), Intent(InOut) :: path
!!$    Integer, Intent(InOut)               :: pathlen
!!$
!!$    Integer :: i,nextatom,pathl
!!$
!!$    pathlen = pathlen + 1
!!$    path(pathlen) = atom1 
!!$    pathl = 0
!!$
!!$    If (atom1 == atom2) Then
!!$      Return
!!$    End If
!!$
!!$    !** search the tree
!!$    Do i = 1, connect(atom1)%nconnect
!!$      nextatom = connect(atom1)%clist(i)%atom
!!$      If (nextatom == atom2) Then
!!$        pathlen = pathlen+1
!!$        path(pathlen) = atom2
!!$        !** We're done
!!$        pathl = pathlen
!!$        Exit
!!$      Else If (nextatom /= path(pathlen-1)) Then
!!$        If (findint(path(1:pathlen),nextatom) /= 0) Cycle
!!$        pathl = connects_getpath(connect,nextatom,atom2,path,pathlen)
!!$        If (pathl /= 0) Then
!!$          exit
!!$        End If
!!$      End If
!!$    End Do
!!$
!!$    If (pathl == 0) Then
!!$      path(pathlen) = 0
!!$      pathlen = pathlen-1
!!$    End If
!!$
!!$  End Function connects_getpath
!!$
!!$  !---------------------------------------------------------------------------
!!$  ! Generates the connections for an atom type-atom type pair
!!$  !---------------------------------------------------------------------------
!!$  Subroutine connects_regenerateconnection(connect,atom1,atom2,atomt1, &
!!$      atomt2,length)
!!$    
!!$    Type(AtomConnections), Dimension(:), Intent(InOut) :: connect
!!$    Integer, Intent(In) :: atom1, atom2, atomt1, atomt2
!!$    Real(kind=RDbl), Intent(In) :: length
!!$    Real(kind=RDbl) :: blength, tolerance
!!$    Real(kind=RDbl) :: maxlen, minlen
!!$
!!$    !** Get the tolerance and bond length for the pair
!!$    tolerance = connections(atomt1,atomt2)%tol
!!$    blength = connections(atomt1,atomt2)%length
!!$
!!$    !No bond length recorded. Exit.
!!$    If (blength == 0.0_RDbl) Return
!!$
!!$    !Calculate the upper and lower bound on the bond length
!!$    maxlen = blength*(1.0_RDbl+tolerance)
!!$    minlen = blength*(1.0_RDbl-tolerance)
!!$
!!$    !** Compare the lengths and make the connection if necessary
!!$    If ((length <= maxlen).And.(length >= minlen).And. &
!!$        .Not.connects_iscon(connect(atom1),atom2)) Then
!!$      ! Make the connection
!!$      Call connects_connectatoms(connect,atom1,atom2,length)
!!$    End If
!!$
!!$  End Subroutine Connects_regenerateconnection
!!$
!!$  !----------------------------------------------------
!!$  ! Returns an atom-atom pair list of connections
!!$  !----------------------------------------------------
!!$  Integer Function connects_makepairlist(connect,pairlist)
!!$
!!$    Type(AtomConnections), Dimension(:), Intent(In) :: connect
!!$    Integer, Dimension(:,:) :: pairlist
!!$    Integer, Dimension(Size(pairlist,1),2) :: pairindex
!!$
!!$    Integer :: i,j,k,atom1,atom2
!!$
!!$    pairindex = 0
!!$    pairlist = 0
!!$    k = 0
!!$
!!$    ! This would work better with a linked list or tree for sorting
!!$    ! However, I want to keep this a bit simplier for the time being.
!!$    ! We traverse the connect array, which is ordered from 1 to natoms,
!!$    ! so pairlist will be ordered with pairs including 1 first, 2 second,
!!$    ! etc. We keep track of the location of each section of atoms, i.e.,
!!$    ! where do the pairs with 1 begin and end, the pairs with 2, etc.
!!$
!!$    ! Put the first pair in the list
!!$    k = k+1
!!$    pairlist(k,1) = Min(1,connect(1)%clist(1)%atom)
!!$    pairlist(k,2) = Max(1,connect(1)%clist(1)%atom)
!!$    pairindex(1,1) = k
!!$    pairindex(1,2) = k
!!$
!!$    Do i = 1,Size(connect)
!!$      Do j = 1,connect(i)%nconnect
!!$        ! Get the pair
!!$        atom1 = Min(i,connect(i)%clist(j)%atom)
!!$        atom2 = Max(i,connect(i)%clist(j)%atom)
!!$        ! Now check to see if it's in the array yet        
!!$        If (pairindex(atom1,1) == 0) Then
!!$          ! Just add it to the array since we don't have a
!!$          ! starting and ending location stored
!!$          k = k+1
!!$          pairlist(k,1) = atom1
!!$          pairlist(k,2) = atom2
!!$          pairindex(atom1,1) = k
!!$          pairindex(atom1,2) = k
!!$        Else If (findint(pairlist(pairindex(atom1,1):& 
!!$            pairindex(atom1,2),2),atom2) == 0) Then
!!$          ! The pair isn't in the array. Add it.
!!$          k = k+1
!!$          pairlist(k,1) = atom1
!!$          pairlist(k,2) = atom2
!!$          pairindex(atom1,2) = k
!!$        End If
!!$      End Do
!!$    End Do
!!$
!!$    connects_makepairlist = k
!!$
!!$  End Function connects_makepairlist
!!$
!!$  !---------------------------------------------
!!$  ! Return the number of branch points
!!$  !---------------------------------------------
!!$  Integer Function connects_nbranches(connect,branchindx)
!!$    Type(AtomConnections), Dimension(:), Intent(In) :: connect
!!$    Integer, Dimension(:), Optional, Intent(InOut) :: branchindx
!!$
!!$    Integer :: i, nbranch
!!$
!!$    nbranch = 0
!!$
!!$    Do i = Lbound(connect,1),Ubound(connect,1)
!!$      If (connect(i)%nconnect > 2) Then
!!$        nbranch = nbranch+1
!!$        If (Present(branchindx)) Then
!!$          branchindx(nbranch) = i
!!$        End If
!!$      End If
!!$    End Do
!!$    connects_nbranches = nbranch
!!$  End Function connects_nbranches
!!$
!!$  !---------------------------------------------------------------------------
!!$  ! Check to see if there are any rings
!!$  !---------------------------------------------------------------------------
!!$  Integer Function connects_isrings(connect,atom,rings)
!!$    Type(AtomConnections), Dimension(:), Intent(In) :: connect
!!$    Integer, Intent(In) :: atom
!!$    Integer, Dimension(:), Intent(Out) :: rings
!!$    Integer :: nbranches
!!$    Integer, Dimension(size(connect)) :: branches
!!$    Logical :: isend
!!$
!!$    isend = .False.
!!$    rings = 0
!!$    
!!$    !** Get the branches
!!$    nbranches = connects_nbranches(connect,branches)
!!$    
!!$!    If (nbranches /= 0) Then
!!$!      Do i = 1,nbranches
!!$!        Do j = 1,connect(branches(i))%nconnect
!!$!          ctemp = connect
!!$!          ctemp(branches(i))%nconnect = ctemp(branches(i))%nconnect-1
!!$!          ctemp(branches(i))%clist(j) = &
!!$!              connect(branches(i))%clist(connect(branches(i))%nconnect)
!!$!          pathlen = connects_getpathdriver(ctemp,branches(i), &
!!$!              ctemp(branches(i))%clist(j)%atom)
!!$!        End Do
!!$!      End Do
!!$!    End If
!!$
!!$    connects_isrings = 0
!!$
!!$  End Function connects_isrings
!!$
!!$  !---------------------------------------------------------------------------
!!$  ! Get the rings
!!$  !---------------------------------------------------------------------------
!!$  Integer Function connects_getringsdriver(connect,atom,rings)
!!$    Type(AtomConnections), Dimension(:), Intent(In) :: connect
!!$    Integer, Intent(In) :: atom
!!$    Integer, Dimension(:), Optional :: rings
!!$    Integer, Dimension(Size(connect)) :: ringstemp, branches, visited
!!$    Integer :: nrings, atomn, ringz, j,nbranches
!!$
!!$    atomn = atom
!!$    j = 0
!!$    nrings = 0
!!$    visited = 0
!!$    ringstemp = 0
!!$
!!$    !** Get the branches
!!$    nbranches = connects_nbranches(connect,branches)
!!$
!!$    !** Find the ends and rings
!!$    ringz = connects_getrings(connect,atomn,atomn,visited, &
!!$        nrings,ringstemp)
!!$
!!$    !** Return just the rings
!!$!    Do i = 1,nrings
!!$!      If (connect(ringstemp(i))%nconnect < 2) Then
!!$!        ringz = ringz-1
!!$!      Else If (Present(rings)) Then
!!$!        j = j+1
!!$!        rings(j) = ringstemp(i)
!!$!      End If
!!$!    End Do
!!$    rings = ringstemp
!!$    connects_getringsdriver = ringz
!!$  End Function connects_getringsdriver
!!$
!!$  !---------------------------------------------------------------------------
!!$  ! Find ends and rings
!!$  !---------------------------------------------------------------------------
!!$  Recursive Integer Function connects_getrings(connect,currentatom,lastatom, &
!!$      visited,nrings,rings)
!!$
!!$    Type(AtomConnections), Dimension(:), Intent(In) :: connect
!!$    Integer, Intent(InOut) :: currentatom,lastatom,nrings
!!$    Integer, Dimension(:), Optional :: rings
!!$    Integer, Dimension(:), Intent(InOut) :: visited
!!$    Integer :: nextatom,i,ringz
!!$
!!$    ringz = 0
!!$
!!$    !** Update visited list
!!$    visited(currentatom) = visited(currentatom) + 1
!!$    If (visited(currentatom) > 1) Then
!!$      !this is a ring atom
!!$      ringz = ringz+1
!!$      nrings = nrings+1
!!$      rings(nrings) = currentatom
!!$      connects_getrings = ringz
!!$      Return
!!$    End If
!!$
!!$    Do i = 1,connect(currentatom)%nconnect
!!$!      nextatom = connect(currentatom)%clist(i)%atom
!!$!      If ((nextatom == lastatom).And. &
!!$!          (connect(currentatom)%nconnect == 1)) Then
!!$!        !** Found an end!
!!$!        nrings = nrings + 1
!!$!        ringz = ringz+1
!!$!        rings(nrings) = currentatom
!!$      If (nextatom == lastatom) Then
!!$        cycle
!!$      Else
!!$        !** Continue to search
!!$        ringz = connects_getrings(connect,nextatom,currentatom,visited, &
!!$            nrings,rings)+ringz
!!$      End If      
!!$    End Do
!!$
!!$    connects_getrings = ringz
!!$
!!$  End Function connects_getrings
!!$
!!$  !---------------------------------------------------------------------------
!!$  ! Get the ends of the branches
!!$  !---------------------------------------------------------------------------
!!$  Integer Function connects_getendsdriver(connect,atom,ends)
!!$    Type(AtomConnections), Dimension(:), Intent(In) :: connect
!!$    Integer, Intent(In) :: atom
!!$    Integer, Dimension(:), Optional :: ends
!!$    Integer, Dimension(Size(connect)) :: endstemp,visited
!!$    Integer :: nends, atomn, endz, j,i
!!$
!!$    atomn = atom
!!$    j = 0
!!$    nends = 0
!!$    visited = 0
!!$
!!$    !** Get the branches
!!$    nbranches = connects_nbranches(connect,branches)
!!$    !MDEBUG
!!$    Write(*,*) "nbranches, branches :",nbranches,branches
!!$
!!$
!!$    !** Find the ends and rings
!!$    endz = connects_getends(connect,atomn,atomn,visited, &
!!$        nends,endstemp)
!!$
!!$    !** Return just the ends
!!$    Do i = 1,nends
!!$      If (connect(endstemp(i))%nconnect > 1) Then
!!$        endz = endz-1
!!$      Else If (Present(ends)) Then
!!$        j = j+1
!!$        ends(j) = endstemp(i)
!!$      End If
!!$    End Do
!!$
!!$    connects_getendsdriver = endz
!!$  End Function connects_getendsdriver
!!$
!!$  !---------------------------------------------------------------------------
!!$  ! Find ends and rings
!!$  !---------------------------------------------------------------------------
!!$  Recursive Integer Function connects_getends(connect,currentatom,lastatom, &
!!$      visited,nends,ends)
!!$
!!$    Type(AtomConnections), Dimension(:), Intent(In) :: connect
!!$    Integer, Intent(InOut) :: currentatom,lastatom,nends
!!$    Integer, Dimension(:), Optional :: ends
!!$    Integer, Dimension(:), Intent(InOut) :: visited
!!$    Integer :: nextatom,i,endz
!!$
!!$    endz = 0
!!$
!!$    !** Update visited list
!!$    visited(currentatom) = visited(currentatom) + 1
!!$    If (visited(currentatom) > 1) Then
!!$      !this is a ring atom
!!$      endz = endz+1
!!$      nends = nends+1
!!$      ends(nends) = currentatom
!!$      connects_getends = endz
!!$      Return
!!$    End If
!!$
!!$    Do i = 1,connect(currentatom)%nconnect
!!$      nextatom = connect(currentatom)%clist(i)%atom
!!$      If ((nextatom == lastatom).And. &
!!$          (connect(currentatom)%nconnect == 1)) Then
!!$        !** Found an end!
!!$!        nends = nends + 1
!!$!        endz = endz+1
!!$!        ends(nends) = currentatom
!!$      Else If (nextatom == lastatom) Then
!!$        cycle
!!$      Else
!!$        !** Continue to search
!!$        endz = connects_getends(connect,nextatom,currentatom,visited, &
!!$            nends,ends)+endz
!!$      End If      
!!$    End Do
!!$
!!$    connects_getends = endz
!!$
!!$  End Function connects_getends

  !--------------------------------------------------------------------
  ! Generates coords for a gnuplot file
  ! Requires:  connect -- Array of Atom Connections data structures
  !            coords -- position vector for each atom
  !            filename -- file name to dump coordinates into
  !--------------------------------------------------------------------
  Subroutine connects_connect2plot(connect,coords,filename)
    Type(AtomConnections), Dimension(:), Intent(In) :: connect
    Type(VecType), Dimension(:), Intent(In)         :: coords
    Character(*), Intent(In)                        :: filename

    Integer                             :: i,unitno,error,npairs
    Character(len=strLen)               :: frmt
    Integer, Dimension(Size(connect),2) :: pairlist

    !** Make a pair list
!!$    npairs = connects_makepairlist(connect,pairlist)
    npairs = connects_makeXList(connect,2,pairlist)

    unitno = file_getunit(filename)
    Open(unitno,file=filename,iostat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a,a)') __FILE__,": ",__LINE__, &
          " : Error opening file ",trim(filename)
      Stop
    End If
    frmt = "(1x,i3,2x,e10.4,2x,e10.4,2x,e10.4)"

    Do i = 1, npairs
      Write(unitno,frmt) pairlist(i,1),coords(pairlist(i,1))
      Write(unitno,frmt) pairlist(i,1),coords(pairlist(i,2))
    End Do

    Close(unitno)

#ifdef NOTWORKING
    !** Get the number of branches
!    j = connects_getnbranches(connect,branchpnts)

    If (nbranches == 0) Then
      !** easy straight chain
      j = connects_getpathdriver(connect,Lbound(connect,1), &
          Ubound(connect,1),path)
      Do i = 1,j
        Write(unitno,frmt) path(i),coords(path(i))
      End Do
    Else
      Do i = 1,nbranches
        !** Get the ends of the branches
!        nends = connects_getends(connect,branchpnts)
        !** Get a path
        j = connects_getpathdriver(connect,Lbound(connect,1),branchpnts(i),path)
        Do k = 1,j
          Write(unitno,frmt) path(k),coords(path(k))
        End Do
        Write(unitno,'(a)') ''
        Write(unitno,'(a)') ''
      End Do
    End If
#endif

  End Subroutine connects_connect2plot

  !-----------------------------------------------------------------------
  ! Display the content for a single atom
  ! Requires:  connect -- single Atom Connections data structure
  !-----------------------------------------------------------------------
  Function connects_display(connect)
    Character(len=lstrLen)            :: connects_display
    Type(AtomConnections), Intent(In) :: connect

    Integer               :: i
    Character(len=strLen) :: strformat

    If (connect%nconnect <= 0) Then
      Write(connects_display,'(a)') 'No connections'
      Return
    End If

    Write(strformat,'(a,i2,a,i2,a)') "(1x,",connect%nconnect,"i5,1x,", & 
        connect%nconnect,"f8.3)"
    Write(connects_display,strformat) &
        (connect%clist(i)%atom,i=1,connect%nconnect), &
        (connect%clist(i)%length,i=1,connect%nconnect)

  End Function connects_display

  !--------------------------------------------------------------------
  ! Display all the information to unitio
  ! Requires:  connect -- Array of Atom Connections data structures
  !            unit -- unit number to dump to
  !--------------------------------------------------------------------
  Subroutine connects_displayall(connect,unit)
    Type(AtomConnections), Dimension(:), Intent(In) :: connect
    Integer, Intent(In)                             :: unit
    
    Integer                  :: j
    Character(len=lstrLen)   :: string

    Write(unit,'(2a)') "Connect information (Atom, Atoms ", &
        "Connected To, Bond Lengths)"
    Do j = 1,Size(connect)
      string = connects_display(connect(j))
      Write(unit, '(a,i3,2a)') "Atom ",j,": ",Trim(string)
    End Do

  End Subroutine connects_displayall

  !---------------------------------------------------------------------------
  ! Displays the connection information matrix
  ! Requires:  conmatrix -- Connection Matrix structure to initialize
  !            unit -- unit number to dump to
  !---------------------------------------------------------------------------
  Subroutine connects_displayArray(conmatrix,unit)
    Type(ConnectMatrix), Intent(In) :: conmatrix
    Integer, Intent(In)             :: unit

    Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
        '  routine not finished'
    Stop

  End Subroutine connects_displayArray

  !----------------------------------------------------------------------
  ! Cleanup a single atom connection type
  ! Requires:  connect -- single Atom Connections data structure
  !----------------------------------------------------------------------
  Subroutine connects_cleanup(connects)
    Type(AtomConnections), Intent(InOut) :: connects

    Integer :: error

    If (Associated(connects%clist)) Then
      Deallocate(connects%clist,stat=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,"connect%clist")    
    End If

  End Subroutine connects_cleanup

End Module connects



