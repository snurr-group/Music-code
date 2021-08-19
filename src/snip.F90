!----------------------------------------------------------------------------
! This module handles the creation of a 'host' and a 'cluster' from a single
! structure.
! 
! Needed Improvements:
! 1) how to generalize to structures other than zeolites?
!----------------------------------------------------------------------------

Module snip

  Use defaults, Only: strLen, RDbl, dashedline, zero, radTodeg
  Use utils, Only: split, stripcmnt, filesrchstr, toupper, allocerrdisplay, &
      deallocerrdisplay, condenseint, subintlist, unionintlist
  Use file, Only: file_getunit,file_open
  Use vector, Only: VecType, vector_getnorm, Assignment(=), Operator(+), &
      Operator(-), Operator(*), Operator(/), vector_angle, vector_getdist, &
      vector_getunitvec, vector_getplanenorm
  Use fundcell, Only: Fundamental_Cell, fundcell_angles, fundcell_init, &
      fundcell_getell, fundcell_slant, fundcell_unslant, fundcell_latticevec, &
      fundcell_isortho, fundcell_snglMinImage, fundcell_minwidth, &
      fundcell_display
  Use readstruc, Only: Structure, readstruc_natoms, readstruc_subcopy, &
      readstruc_recenter, readstruc_xform, readstruc_initcopy, &
      readstruc_visxyz, readstruc_dumpcar
  Use connects, Only: AtomConnections, ConnectMatrix, connects_generateAll, &
      connects_initarray, connects_initconnection, connects_allneighbors, &
      connects_initmaxcon

  Implicit None
  Save

  Private
  Public  :: snip_zeo
  
  Real(kind=RDbl), Parameter        :: len_AlOH = 0.9628_RDbl
  Real(kind=RDbl), Parameter        :: len_SiOH = 0.9666_RDbl
  Real(kind=RDbl), Parameter        :: len_maxSiO = 1.70_RDbl
  Real(kind=RDbl), Parameter        :: len_maxAlO = 1.80_RDbl

Contains

  !----------------------------------------------------------------------------
  ! Takes in the full structure, the fundamental cell that defines the
  ! PBCS and the desired cluster size and returns the cluster.  Passing in
  ! a pre-determined connectivity array will speed the process.
  ! Requires:  full -- 'host' structure
  !            fcell -- fundamental cell structure
  !            tatom -- the T-atom number, around which the cluster is centered
  !            csize -- cluster size
  !            method -- specifies method for partial shell identification
  !            cluster -- returned 'cluster' structure (uninitialized)
  !            connects -- optional connectivity array 
  !----------------------------------------------------------------------------
  Subroutine snip_zeo(full,fcell,tatom,csize,method,cluster,connects)
    Type(Structure), Intent(In)            :: full
    Type(Fundamental_Cell), Intent(In)     :: fcell    
    Integer, Intent(In)                    :: tatom,csize
    Character(*), Intent(In)               :: method
    Type(Structure), Intent(Out)           :: cluster
    Type(AtomConnections), Dimension(:), Intent(In), Target, Optional :: connects

    Integer                     :: natoms,error,a,ncatoms,xaway
    Logical, Save               :: firsttime = .True.
    Type(ConnectMatrix), Save   :: conmatrix
    Type(VecType)               :: shift
    Type(Structure)             :: copy
    Integer, Dimension(100)                       :: catoms
    Integer, Dimension(:), Allocatable            :: atypes
    Type(AtomConnections), Dimension(:), Pointer  :: connectivity

    !** Get the number of atoms in the structure
    natoms = readstruc_natoms(full)

    !** Establish the connectivity matrix with 5 atypes (Si,O,Al,H,Ge)
    If (firsttime) Then
      firsttime = .False.
      Call connects_initarray(conmatrix,6)
      Call connects_initconnection(conmatrix,1,2,1.7_RDbl,0.2_RDbl)  !** Si-O
      Call connects_initconnection(conmatrix,2,3,1.7_RDbl,0.2_RDbl)  !** Al-O
      Call connects_initconnection(conmatrix,2,4,1.0_RDbl,0.2_RDbl)  !** O-H
      Call connects_initconnection(conmatrix,4,4,1.0_RDbl,0.1_RDbl)  !** H-H
      Call connects_initconnection(conmatrix,4,5,1.1_RDbl,0.2_RDbl)  !** C-H
      Call connects_initconnection(conmatrix,5,2,1.7_RDbl,0.2_RDbl)  !** Ge-O
    End If

    !** Construct the array of atom types
    Allocate(atypes(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do a = 1,natoms
      Select Case(ToUpper(full%elements(a)))
      Case ('SI')
        atypes(a) = 1
      Case ('O')
        atypes(a) = 2
      Case ('AL')
        atypes(a) = 3
      Case ('H')
        atypes(a) = 4
      Case ('GE')
        atypes(a) = 5
      Case ('C')
        atypes(a) = 6
      Case Default
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' Could not interpret element type: ',Trim(full%elements(a))
        Stop
      End Select
    End Do

    !** Copy and recenter the full structure around the Tatom
    Call readstruc_initcopy(copy,full)
    Call readstruc_recenter(copy,tatom,fcell,shift)
    shift = shift*(-1.0_RDbl)
    Call readstruc_xform(copy,shift)

    !** Establish the connectivity, either from scratch or from input
    If (Present(connects)) Then
      connectivity => connects
    Else
      Allocate(connectivity(natoms), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Do a = 1,natoms
        Call connects_initmaxcon(connectivity(a),4)
      End Do
      Call connects_generateAll(connectivity,atypes,copy%coords, &
          conmatrix,fcell)
    End If

    !** Get the atom numbers in the cluster
    Select Case(csize)
    Case (1)
      ncatoms = snip_cluster1(copy,fcell,tatom,method,catoms,connectivity)
    Case (2)
      ncatoms = snip_cluster2(copy,fcell,tatom,catoms,connectivity)
    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Could not interpret cluster size parameter: ',csize
      Stop
    End Select

    !** Remove repeats from the list of cluster atoms
    ncatoms = condenseint(catoms(1:ncatoms))

    !** Copy select atom numbers to give the cluster
    cluster = readstruc_subcopy(copy,catoms(1:ncatoms))

    !** Terminate the cluster with hydrogen atoms
    Call snip_termzeo(cluster)

    !** Clean the atypes and the connectivity
    Deallocate(atypes, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)          
    If (Present(connects)) Then
      Nullify(connectivity)
    Else
      Deallocate(connectivity, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)      
    End If

  End Subroutine snip_zeo

  !----------------------------------------------------------------------------
  ! Performs the operations necessary to pick atoms for the 1.5 T-atom 
  ! coordination sphere cluster.  This means all the T-atoms 2 atoms away
  ! from the central T-atom plus those 4 atoms away in the direction of
  ! smallest curvature.
  ! Requires:  full -- 'host' structure
  !            fcell -- fundamental cell structure
  !            tatom -- the T-atom number, around which the cluster is centered
  !            border_select -- string containing selector for method
  !            catoms -- array containing cluster atom numbers
  !            connects -- initialized connectivity array 
  !----------------------------------------------------------------------------
  Integer Function snip_cluster1(full,fcell,tatom,border_select,catoms,connects)
    Type(Structure), Intent(In)            :: full
    Type(Fundamental_Cell), Intent(In)     :: fcell    
    Integer, Intent(In)                    :: tatom
    Character(*), Intent(In)               :: border_select
    Integer, Dimension(:), Intent(InOut)   :: catoms
    Type(AtomConnections), Dimension(:), Intent(In) :: connects

    Integer                     :: i,j,n,a
    Integer                     :: natoms,ncatoms,xaway,oatom
    Real(kind=RDbl)             :: maxangle,angle,mindist,dist,mindot
    Type(VecType)               :: vec1,vec2,vec3
    Integer, Dimension(2)       :: pair,notpair
    Integer, Dimension(4)       :: tatoms,termatoms,nneighbors,dummy
    Integer, Dimension(20)      :: fromborder,results
    Integer, Dimension(4,10)    :: tneighbors
    Type(VecType), Dimension(4) :: normal

    !** Get the core cluster, up to 2 atoms away
    xaway = 2
    ncatoms = connects_allneighbors(connects,tatom,catoms,xaway)

!    snip_cluster1 = ncatoms
!    return

    !** Get the atom numbers for the four border T-atoms
    natoms = 0
    Do a = 1,ncatoms
      If ((ToUpper(full%elements(catoms(a))) == 'SI').Or. &
          (ToUpper(full%elements(catoms(a))) == 'GE')) Then
        natoms = natoms + 1
        tatoms(natoms) = catoms(a)
      End If
    End Do
    If (natoms /= 4) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Could not find the four border Si atoms'
      Write(0,*) tatoms(1:natoms)
      Stop
    End If

    !** Add the T-atom to the list
    ncatoms = ncatoms + 1
    catoms(ncatoms) = tatom

    !** Get the neighbors up to 2 away for each of the border T-atoms
    xaway = 2
    Do i = 1,4
      natoms = connects_allneighbors(connects,tatoms(i), &
          tneighbors(i,:),xaway)
      natoms = subintlist(tneighbors(i,1:natoms),catoms(1:ncatoms))
      nneighbors(i) = natoms

      !** Find the atom that makes the smallest angle with the previous Tatoms
      termatoms(i) = snip_pickdirect(full,tatom,tatoms(i), &
          tneighbors(i,1:natoms))
    End Do

    !** Select the border T-atoms based on given criteria
    Select Case(ToUpper(border_select))
    Case('MAXANGLE')
      !** Find the two "terminal" atoms that make the largest angle with Tatom
      maxangle = 0.0_RDbl
      Do i = 1,4
        vec1 = full%coords(termatoms(i))
        Do j = i+1,4
          vec3 = full%coords(termatoms(j))
          angle = vector_angle(vec1,full%coords(tatom),vec3)
  
          !** Get the atom pairs that have the two maximum angles
          If (angle > maxangle) Then
            maxangle = angle
            pair(1) = i
            pair(2) = j
          End If
        End Do
      End Do

    Case('MINDIST')
      !** Find the two "terminal" atoms that have the shortest separation distance
      mindist = 1000.0_RDbl
      Do i = 1,4
        vec1 = full%coords(termatoms(i))
        Do j = i+1,4
          vec2 = full%coords(termatoms(j))
          dist = vector_getdist(vec1,vec2)
  
          !** Get the atom pairs that have the min distance
          If (dist < mindist) Then
            mindist = dist
            notpair(1) = i
            notpair(2) = j
          End If
        End Do
      End Do
  
      !** Get the Tatom numbers of the other (non-min distance) pair
      dummy = (/1,2,3,4/)
      n = subintlist(dummy,notpair)
      If (n /= 2) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' Could not eliminate atoms properly'
        Stop
      End If
      pair = dummy(1:n)
    
    Case('INPLANE')
      !** Select the border atoms based on the min angle between the
      !** plane normals formed by the center T-atom, oxygen and border T-atom.
      Do i = 1,4
        !** Get all atoms one away from border atom
        fromborder = 0
        n = connects_allneighbors(connects,tatoms(i),fromborder,1)

        !** Find oxygen atom connecting the two T-atoms
        n = unionintlist(catoms(1:ncatoms),fromborder(1:n),results)
        If (n > 1) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              ' Not able to find oxygen atom uniquely'
          Write(0,*) results(1:n)
          Stop
        Else
          oatom = results(1)
        End If

        !** Get the normal to the plane
        normal(i) = vector_getplanenorm((/full%coords(tatom), &
            full%coords(oatom),full%coords(tatoms(i))/))
      End Do

      !** Find the pair of border T-atoms with the smallest angle between planes
      mindot = 1000.0_RDbl
      Do i = 1,4
        Do j = i+1,4        
          angle = normal(i)*normal(j)
          If (angle < -0.5_RDbl) angle = angle + 1.0_RDbl
          If (Abs(angle) < mindot) Then
            pair(1) = i
            pair(2) = j
            mindot = Abs(angle)
          End If
        End Do
      End Do

    Case('GELABELED')
      !** Find the two Ge-labeled atoms 
      pair = 0
      Do i = 1,4
        n = tatoms(i)
        If (ToUpper(full%elements(n)) == 'GE') Then
          If (pair(1) == 0) Then
            pair(1) = i
          Else
            pair(2) = i
          End If
        End If
      End Do

      If ((pair(1) == 0).Or.(pair(2) == 0)) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Could not find two Germanium atoms in terminal set'
        Stop
      End If

    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Could not interpret border atom selection string: ', &
          Trim(border_select)
      Stop
    End Select

    !** Add the unique neighbors of these two maximum angle terminal atoms
    Do i = 1,2
      n = nneighbors(pair(i))
      catoms(ncatoms+1:ncatoms+n) = tneighbors(pair(i),1:n)
      ncatoms = ncatoms + n
    End Do

    snip_cluster1 = ncatoms

  End Function snip_cluster1

  !----------------------------------------------------------------------------
  ! Performs the operations necessary to pick atoms for the 2.0 T-atom 
  ! coordination sphere cluster.  This means all the T-atoms 4 atoms away
  ! from the central T-atom.
  ! Requires:  full -- 'host' structure
  !            fcell -- fundamental cell structure
  !            tatom -- the T-atom number, around which the cluster is centered
  !            catoms -- array containing cluster atom numbers
  !            connects -- initialized connectivity array 
  !----------------------------------------------------------------------------
  Integer Function snip_cluster2(full,fcell,tatom,catoms,connects)
    Type(Structure), Intent(In)            :: full
    Type(Fundamental_Cell), Intent(In)     :: fcell    
    Integer, Intent(In)                    :: tatom
    Integer, Dimension(:), Intent(InOut)   :: catoms
    Type(AtomConnections), Dimension(:), Intent(In) :: connects

    Integer                     :: i,j,n,a
    Integer                     :: natoms,ncatoms,xaway
    Real(kind=RDbl)             :: maxangle,angle,mindist,dist
    Type(VecType)               :: vec1,vec2,vec3
    Integer, Dimension(2)       :: pair,notpair
    Integer, Dimension(4)       :: tatoms,termatoms,nneighbors,dummy
    Integer, Dimension(4,10)    :: tneighbors

    !** Get the cluster, up to 4 atoms away, easy
    xaway = 4
    ncatoms = connects_allneighbors(connects,tatom,catoms,xaway)

    !** Add the central T-atom to the list
    ncatoms = ncatoms + 1
    catoms(ncatoms) = tatom

    snip_cluster2 = ncatoms

  End Function snip_cluster2

  !----------------------------------------------------------------------------
  ! Picks an atom number from a list that gives the largest angle with two
  ! pre-specified atoms.  Only consider 'SI' atoms for now.
  ! Requires:  full -- 'host' structure
  !            endatom -- forms one end of the 3-atom angle
  !            midatom -- forms the middle of the 3-atom angle
  !            catoms -- array containing cluster atom numbers
  !----------------------------------------------------------------------------
  Integer Function snip_pickdirect(full,endatom,midatom,atoms)
    Type(Structure), Intent(In)         :: full
    Integer, Intent(In)                 :: endatom,midatom
    Integer, Dimension(:), Intent(In)   :: atoms

    Integer               :: i,n,natoms,a,ncatoms,xaway
    Real(kind=RDbl)       :: angle,maxangle
    Type(VecType)         :: vec1,vec2,vec3

    natoms = Size(atoms)
    maxangle = 0.0_RDbl
    snip_pickdirect = 0

    vec1 = full%coords(endatom)
    vec2 = full%coords(midatom)
    Do a = 1,natoms
      If (ToUpper(full%elements(atoms(a))) == 'SI') Then
        vec3 = full%coords(atoms(a))
!        Write(*,*) atoms(a)
!        Write(*,*) vec1
!        Write(*,*) vec2
!        Write(*,*) vec3
        angle = vector_angle(vec1,vec2,vec3)
        If (angle > maxangle) Then
          snip_pickdirect = atoms(a)
          maxangle = angle
        End If
      End If
    End Do

  End Function snip_pickdirect

  !----------------------------------------------------------------------------
  ! Takes in the cluster structure, finds the terminating Si atoms, and 
  ! replaces these atoms with hydrogen atoms at the correct O-H distance.
  ! Also returns a list of the terminating hydrogen atoms
  ! Requires:  full -- 'host' structure
  !            fcell -- fundamental cell structure
  !            tatom -- the T-atom number, around which the cluster is centered
  !            csize -- cluster size
  !            cluster -- returned 'cluster' structure (uninitialized)
  !            connects -- optional connectivity array 
  !----------------------------------------------------------------------------
  Subroutine snip_termzeo(cluster,termatoms)
    Type(Structure), Intent(InOut)                 :: cluster
    Integer, Dimension(:), Intent(Out), Optional   :: termatoms

    Integer          :: i,j,nneighbors,neighbor,tatom,nterm
    Real(kind=RDbl)  :: dist
    Type(VecType)    :: shift
    
    nterm = 0

    !** Loop through atoms in cluster and look for ones to terminate
    Do i = 1,cluster%natoms

      If (ToUpper(cluster%elements(i)) /= 'SI') Cycle
      nneighbors = 0

      !** Check this atom's neighbors
      Do j = 1,cluster%natoms      
        If (i == j) Cycle
        dist = vector_getdist(cluster%coords(i),cluster%coords(j))
        If (dist < len_maxSiO) Then
          nneighbors = nneighbors + 1
          neighbor = j
        End If
        If (nneighbors > 1) Exit
      End Do

      !** Switch this atom to a terminating hydrogen if it meets the criteria
      If (nneighbors == 1) Then
        If (ToUpper(cluster%elements(neighbor)) == 'O') Then
          If (Present(termatoms)) Then
            nterm = nterm + 1
            termatoms(nterm) = i
            If (nterm > Size(termatoms)) Then
              Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
                  ' passed termatoms array too small'
              Stop
            End If 
          End If

          !** Switch to H and begin to build shift vector
          cluster%elements(i) = 'H'
          If (cluster%fftypes_on) Then
            cluster%fftype(i) = 'h1'
          End If
          shift = cluster%coords(i) - cluster%coords(neighbor)
          shift = vector_getunitvec(shift)
          
          !** Get the nearest T-atom type
          tatom = 0
          Do j = 1,cluster%natoms      
            If ((j == neighbor).Or.(j == i)) Cycle
            dist = vector_getdist(cluster%coords(j),cluster%coords(neighbor))
            If (dist < len_maxAlO) Then
              tatom = j
              Exit
            End If
          End Do
          If (tatom == 0) Then
            Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
                ' Could not find Tatom in H-O-Tatom sequence'
            Stop
          End If

          !** Reposition the H terminating atom based on the T-atom type
          Select Case(ToUpper(cluster%elements(tatom)))
          Case ('SI')
            cluster%coords(i) = cluster%coords(neighbor) + shift*len_SiOH
          Case ('GE')
            cluster%coords(i) = cluster%coords(neighbor) + shift*len_SiOH
          Case ('AL')
            cluster%coords(i) = cluster%coords(neighbor) + shift*len_AlOH
          Case ('H')
            !** Do nothing, don't want to terminate an acidic proton
          Case Default
            Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
                ' Could not identify Tatom in H-O-Tatom sequence: ', &
                ToUpper(cluster%elements(tatom))
            Stop
          End Select

        End If
      End If
    End Do

  End Subroutine snip_termzeo


End Module snip
