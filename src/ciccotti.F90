!--------------------------------------------------------------------
! This module contains the data structure and the routines necessary
! for the implementation of Ciccotti, Ferrario and Ryckaert's rigid
! body constraints (Mol. Phys. vol. 47 pg. 1253-1264).  These are used 
! to handle rigid units such as linear, planar or tetrahedral sub-units.  
! If there are too many atoms in one of these units then bond length 
! constraints alone cannot keep the unit rigid.  We get around this 
! problem by reducing the number of atoms in the unit.  Extra atoms 
! are handled as sort of fictitious atoms, they don't get moved 
! directly by the integrator.  All the forces that they feel are mapped 
! appropriately onto the real atoms.  Amit's thesis gives a formulation
! of this method using the Edberg, Evans and Morriss bond constraints.
! You will also need the notes for this algorithm if you want to make
! any sense out of this module.  This algorithm does not affect potentials.
!
! Terminology:
!  The rigid units are called "units" and the group of them for one 
!  is called a "set".  Each molecule can have an arbitrary number 
!  of units.  Unit (0) is the set of atoms representing the atoms that
!  are not part of a rigid unit.  Each unit, other than unit (0), contains
!  both primary and secondary atoms.  The primary atoms are those that
!  actually get moved by the integrator.  Forces on the secondary atoms
!  are mapped onto the primary.  Indices for the primary and secondary
!  atoms are i,j,k and alpha,beta respectively.
!--------------------------------------------------------------------

Module ciccotti

  Use defaults, Only: MAX_ATOMS,RDbl,strLen, lstrLen
  Use general, Only: genparams
  Use vector, Only: VecType, mag, vector_getcomp, vector_getnormsq, &
      Assignment(=), Operator(+), Operator(-), Operator(*), Operator(/), &
      vector_display
  Use utils
  Use file, Only: file_getunit
  Use visxyz, Only: XYZ_Entry,visxyz_make,visxyz_dump
  Use matrixops, Only: matrixops_ludcmp, matrixops_lubksb

  Implicit None
  Save

  Private
  Public :: ciccotti_init, ciccotti_mapforces, CiccottiConnectInfo, &
      CiccottiRigidUnitInfo, CiccottiRigidSetInfo, ciccotti_display, &
      ciccotti_regenSet, ciccotti_getdofred

  Interface ciccotti_mapforces
    Module Procedure ciccotti_mapallforces
    Module Procedure ciccotti_mapunitforces
  End Interface

  Interface ciccotti_display
    Module Procedure ciccotti_displayset
    Module Procedure ciccotti_displayunit
  End Interface

  Type CiccottiConnectInfo
    Integer :: atom1,atom2   !atom numbers of the length constraint
    Integer :: idx1,idx2     !ijk indices used in the stored matrices
    Real(kind=RDbl) :: length
  End Type CiccottiConnectInfo

  Type CiccottiRigidUnitInfo
    Integer                                    :: natoms,np,ns
    Integer, Dimension(MAX_ATOMS)              :: patom,satom
    Real(kind=RDbl), Dimension(:), Pointer     :: patom_invmass
    Real(kind=RDbl), Dimension(:,:), Pointer   :: Calphai,Aalphabeta,Ainv
    Real(kind=RDbl), Dimension(:,:,:), Pointer :: Sijk,Balphajk
    Integer                                    :: ncon
    Type(CiccottiConnectInfo), Dimension(:), Pointer  :: con
  End Type CiccottiRigidUnitInfo

  Type CiccottiRigidSetInfo
    Integer                              :: nunits,totalncon
    Integer                              :: nmobile,natoms
    Integer                              :: nprimary   ! No. of primary atoms
    Integer                              :: nsecondary ! No. of second. atoms
    Type(CiccottiRigidUnitInfo), Dimension(:), Pointer :: unit
  End Type CiccottiRigidSetInfo

Contains

  !----------------------------------------------------------------------------
  ! Initialize the set information and construct matrices.  This
  ! routine should only be used after the atom numbers are fixed
  ! in "molecule".
  ! called with: rigidsetptr -- empty data type to contain info
  !              natoms -- number of atoms in molecule
  !              atomcoords -- list of atomic position vectors
  !              invmasslist -- list of inverse atomic masses
  !              setlist -- list of set numbers for each atom
  !              typelist -- list of type numbers for each atom
  !              atom1list -- list of the first atoms in the constraints
  !              atom2list -- list of the second atoms in the constraints
  !              lengthlist -- list of bond lengths
  ! (later we could have this module generate its own set and type
  ! numbers, but for how, we'll just import them from the molecule)
  !----------------------------------------------------------------------------
  Subroutine ciccotti_init(rigidsetptr,natoms,atomcoords,invmasslist, &
      setlist,typelist,atom1list,atom2list,lengthlist)
    Type(CiccottiRigidSetInfo), Pointer        :: rigidsetptr
    Integer, Intent(In)                        :: natoms
    Type(VecType), Dimension(:), Intent(In)    :: atomcoords
    Real(kind=RDbl), Dimension(:), Intent(In)  :: invmasslist
    Integer, Dimension(:), Intent(In)          :: setlist
    Integer, Dimension(:), Intent(In)          :: typelist
    Integer, Dimension(:), Intent(In)          :: atom1list,atom2list
    Real(Kind=RDbl), Dimension(:), Intent(In)  :: lengthlist

    Integer                                    :: error,a,u,i,j,type
    Integer                                    :: nunitatoms,np,ns
    Integer                                    :: match1,match2
    Logical                                    :: success,failure
    Integer, Dimension(MAX_ATOMS)              :: patom,satom
    Integer, Dimension(0:MAX_ATOMS)            :: nsubcon
    Type(VecType), Dimension(Size(atomcoords))             :: regencoords
    Type(CiccottiConnectInfo), Dimension(Size(atom1list))  :: con,subcon
    Type(CiccottiRigidUnitInfo), Pointer       :: dummyptr
    Character(len=strLen)                      :: tempstrg

    !** Copy the atomic coordinates, to be partially overwritten
    regencoords = atomcoords

    !** Do some counting
    satom = 0
    patom = 0
    rigidsetptr%nmobile = 0
    rigidsetptr%nprimary = 0
    rigidsetptr%nsecondary = 0
    rigidsetptr%nunits = 0
    rigidsetptr%natoms = natoms
    Do a = 1,natoms
      If (setlist(a) == 0) Then
        rigidsetptr%nmobile = rigidsetptr%nmobile + 1
      Else If (setlist(a) > 0) Then
        If (setlist(a) > rigidsetptr%nunits) Then 
          rigidsetptr%nunits = setlist(a)
        End If
        If (typelist(a) == 1) Then
          rigidsetptr%nprimary = rigidsetptr%nprimary + 1
          rigidsetptr%nmobile = rigidsetptr%nmobile + 1
        Else If (typelist(a) == 2) Then
          rigidsetptr%nsecondary = rigidsetptr%nsecondary + 1
        Else
          Write(0,'(2a,i4,a,i2)') __FILE__,":",__LINE__, &
              " type number not reasonable: ",typelist(a)
          Stop 
        End If
      Else 
        Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
            " set number less than zero"
        Stop
      End If
    End Do

    !** Create a temporary depository of all the constraint info
    rigidsetptr%totalncon = Size(atom1list)
    Do i = 1,rigidsetptr%totalncon
      con(i)%atom1 = atom1list(i)
      con(i)%atom2 = atom2list(i)
      con(i)%idx1 = 0
      con(i)%idx2 = 0
      con(i)%length = lengthlist(i)
    End Do

    !** Allocate space for the set information
    If (Associated(rigidsetptr%unit)) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " rigidset%unit has already been allocated"
      Stop
    Else
      Allocate(rigidsetptr%unit(0:rigidsetptr%nunits),stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
            " Could not allocate rigidset%unit"
        Stop
      End If
    End If

    !** Deal with initialization of individual rigid units
    failure = .False.
    Do u = 0,rigidsetptr%nunits
      nunitatoms = 0
      np = 0
      ns = 0
      Do a = 1,natoms
        If (setlist(a) == u) Then
          type = typelist(a)
          If (u == 0) type = 1  !put all unit 0 atoms in patom
          nunitatoms = nunitatoms + 1
          If (type == 1) Then
            np = np + 1
            patom(np) = a
          Else If (type == 2) Then
            ns = ns + 1
            satom(ns) = a
          Else
            Write(0,'(2a,i4,a,i2,a,i2)') __FILE__,":",__LINE__, &
                " ERROR: atom ",a," has unknown type ",typelist(a)
            Stop
          End If
        End If
      End Do

      !** check to be sure that none of the secondary atoms are bonded
!!$      Write(0,*) __FILE__,__LINE__,": Connection check skipped"
      Do i = 1,rigidsetptr%totalncon
        If (findint(satom,con(i)%atom1) /= 0) Then
          Write(0,'(2a,i4,a,i2,a)') __FILE__,":",__LINE__, &
              " ERROR: a secondary atom (",con(i)%atom1,") cannot be bonded"
          Write(0,'(a)') 'atom1, atom2, length'
          Do j = 1,rigidsetptr%totalncon
            Write(0,'(2i4,f8.3)') con(j)%atom1,con(j)%atom2,con(j)%length
          End Do
          Stop
        End If
      End Do
    
      !** make the sub-constraint data structure, want to find
      !** all constraints that are in this single rigid unit
      nsubcon(u) = 0  !number of constraints in this rigid unit
      Do i = 1,rigidsetptr%totalncon
        match1 = 0  !primary atom index that corresponds to first atom
        match2 = 0  !primary atom index that corresponds to second atom
        Do j = 1,np
          If (patom(j) == con(i)%atom1) Then
            match1 = j
          Else If (patom(j) == con(i)%atom2) Then
            match2 = j
          End If
        End Do
        If ((match1 > 0).And.(match2 > 0)) Then
          nsubcon(u) = nsubcon(u) + 1
          subcon(nsubcon(u)) = con(i)
          subcon(nsubcon(u))%idx1 = match1
          subcon(nsubcon(u))%idx2 = match2
          !** place constraints between units into unit 0, but give error
        Else If ((u == 0).And.((match1 > 0).Or.(match2 > 0))) Then
          nsubcon(u) = nsubcon(u) + 1
          subcon(nsubcon(u)) = con(i)
          subcon(nsubcon(u))%idx1 = match1
          subcon(nsubcon(u))%idx2 = match2
          Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
              ' ERROR: bond constraints cannot extend out of rigid units'
          Write(0,'(6a)') 'Problem is in unit ',Trim(int2str(u)), &
              ', there is a bond between atoms ',Trim(int2str(con(i)%atom1)), &
              ' and ',Trim(int2str(con(i)%atom2))
          Stop
        End If
      End Do

      tempstrg=int2str(u)
      If ((u > 0).And.(ns == 0)) Then
        Write(0,'(2a,i4,3a)') __FILE__,":",__LINE__, &
            ' rigid unit ',Trim(tempstrg),' has no secondary atoms'
        Write(0,'(a)') 'please put these atoms into unit "0"'
        Stop
      End If

      dummyptr => rigidsetptr%unit(u)
      success = ciccotti_initRigidUnit(dummyptr,np,ns,patom,satom, &
          atomcoords,invmasslist,subcon(1:nsubcon(u)),regencoords)
      If (.Not. success) failure = .True.
!LC      Call ciccotti_displayunit(dummyptr,6,2)

    End Do

    If (failure) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Fatal problem with regeneration of coordinates"
      Write(0,'(2x,a)') "Dumping regenerated set.  You may be able to use"
      Write(0,'(2x,a)') "this set of coordinates instead, please look."
      Call ciccotti_visunits(rigidsetptr,regencoords,.False.)
!      Call ciccotti_visunits(rigidsetptr,atomcoords,.False.)
      Stop
    End If

!LC    Call ciccotti_visunits(rigidsetptr,atomcoords,.True.)

    !** check the division of the constraints into the rigid units
    If (Sum(nsubcon(0:rigidsetptr%nunits)) /= rigidsetptr%totalncon) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " The sum of constraints in each rigid unit does not agree"
      Write(0,'(a)') 'unit, nconstraints_in_unit'
      Do j = 0,rigidsetptr%nunits
        Write(0,'(2i4)') j,nsubcon(j)
      End Do
      Write(0,'(a,i4)') 'ncon in molec: ',rigidsetptr%totalncon
      Write(0,'(a,i4)') 'sum:           ',Sum(nsubcon(0:rigidsetptr%nunits))
      Stop
    End If

  End Subroutine ciccotti_init

  !------------------------------------------------------------------------
  ! Get the reduction in degrees of freedom due to the rigid units
  ! This reduction is equivalent to 3 times the number of secondary atoms
  ! Requires: rigidsetptr -- set of rigid units (WHOLE molec)
  !------------------------------------------------------------------------
  Integer Function ciccotti_getdofred(rigidsetptr)
    Type(CiccottiRigidSetInfo), Pointer      :: rigidsetptr

    Integer       :: unit

    ciccotti_getdofred = 0  

    Do unit = 1,rigidsetptr%nunits
      ciccotti_getdofred = ciccotti_getdofred + 3*rigidsetptr%unit(unit)%ns
    End Do

  End Function ciccotti_getdofred


  !---------------------------------------------------------------
  ! Map the forces on the secondary atoms back to the primary
  ! atoms.  Then, zero the forces on the secondary atoms
  ! called with: rigidsetptr -- set of rigid units (WHOLE molec)
  !              coords -- one molecule's atomic position vectors
  !              v -- one molecule's velocity vectors
  !              f -- one molecule's force vectors  (correct?)
  !---------------------------------------------------------------
  Subroutine ciccotti_mapallforces(rigidsetptr,coords,v,f)
    Type(CiccottiRigidSetInfo), Pointer             :: rigidsetptr
    Type(VecType), Dimension(:), Intent(In)         :: coords,v
    Type(VecType), Dimension(:), Intent(InOut)      :: f
 
    Integer       :: a,unit

    Do unit = 1,rigidsetptr%nunits
      Call ciccotti_mapunitforces(rigidsetptr%unit(unit),coords,v,f)
    End Do

    Return

    Write(*,*)
    Write(*,*) __FILE__,': Finished mapping forces on molecule '
    Write(*,'(2a,i4,a)') __FILE__,": ",__LINE__, &
        "  Dumping new forces and velocities"
    Do a = 1,rigidsetptr%natoms
      Write(*,'(1x,i10,a,3e14.4)') a,'  Force:    ',f(a)
      Write(*,'(1x,i10,a,3e14.4)') a,'  Velocity: ',v(a)
    End Do

  End Subroutine ciccotti_mapallforces


  !---------------------------------------------------------------
  ! Map the forces on the secondary atoms back to the primary
  ! atoms.  Then, zero the forces on the secondary atoms
  ! called with: rs -- single rigid unit info pointer
  !              coords -- one molecule's atomic position vectors
  !              v -- one molecule's velocity vectors
  !              f -- one molecule's force vectors  (correct?)
  !---------------------------------------------------------------
  Subroutine ciccotti_mapunitforces(rs,coords,v,f)
    Type(CiccottiRigidUnitInfo), Intent(In)          :: rs
    Type(VecType), Dimension(:), Intent(In)          :: coords,v
    Type(VecType), Dimension(:), Intent(InOut)       :: f

    Integer                                     :: i,j,k
    Integer                                     :: alpha,beta
    Integer                                     :: satom,patom
    Integer                                     :: jk,jatom,katom
    Real(Kind=RDbl)                             :: imass,velocity2
    Type(VecType), Dimension(rs%np,rs%ncon)     :: scriptr
    Type(VecType), Dimension(rs%np)             :: scriptf
    Type(VecType), Dimension(rs%ns)             :: tvec
    Type(VecType), Dimension(rs%ncon)           :: bondvec
    Type(VecType)                               :: Fjk,vec
    Real(Kind=RDbl), Dimension(rs%ncon,rs%ncon) :: A
    Real(Kind=RDbl), Dimension(rs%ncon)         :: b
    Real(Kind=RDbl)                             :: nperm  
    Integer, Dimension(rs%ncon)                 :: perm   

#if DEBUG
    Write(*,*) 'ciccotti_mapunitforces, input coords, velocities, forces: '
    Do i = 1,Size(f)
      Write(*,'(i2,2x,a)') i,Trim(vector_display(coords(i),'f10.6'))
      Write(*,'(4x,a)') Trim(vector_display(v(i),'f10.6'))
      Write(*,'(4x,a)') Trim(vector_display(f(i),'f10.6'))
    End Do                                             
#endif

    !** Calculate the Talpha vector for the secondary atoms
    Do alpha = 1,rs%ns
      satom = rs%satom(alpha)
      tvec(alpha) = (-1.0d0)*f(satom)
      Do i = 1,rs%np
        patom = rs%patom(i)
        vec = rs%Calphai(alpha,i)*f(patom)
        tvec(alpha) = tvec(alpha) + vec
      End Do
    End Do

    !** Calculate the script force vectors
    scriptf = VecType(0.0_RDbl)
    Do i = 1,rs%np
       
      patom = rs%patom(i)
      imass = rs%patom_invmass(i)
      vec = 0.0_RDbl
      Do alpha = 1,rs%ns
        Do beta = 1,rs%ns
          vec = vec + rs%Calphai(alpha,i)*rs%Ainv(alpha,beta)*tvec(beta)
        End Do
      End Do

      scriptf(i) = f(patom) - vec*imass
   
#if DEBUG 
      Write(*,'(7x,a,i4)') 'primary atom ',i
      Write(*,'(4x,3f14.3,a)') f(patom),' (original force on atom)'
      Write(*,'(4x,3f14.3,a)') vec*((-1.0d0)*imass), &
          ' (contributed force from secondaries)'
      Write(*,'(4x,3f14.3,a)') scriptf(i),' (weighted resultant)'
      Write(*,*) 'inverse mass: ',imass
#endif
    End Do

    !** zero the forces on the secondary atoms, don't need them anymore
    Do alpha = 1,rs%ns
      f(rs%satom(alpha)) = 0.0_Rdbl
    End Do

    !** DEAL WITH CONSTRAINT FORCES HERE
    !** Loop through the molecules of this sorbate type
    !** NOTE: we need to construct the system Ax=b 
    !**       (x is the vector of constraint "forces" that we're looking for)
    Do jk = 1,rs%ncon    !loop over the constraints on the primary atoms
    
      !** Construct the constraint bond vectors, Here: (Rjk = Rj - Rk)
      !**   jk = the connection number in the data structure
      jatom = rs%con(jk)%atom1    !these are ATOM indices
      katom = rs%con(jk)%atom2

      j = rs%con(jk)%idx1    
      k = rs%con(jk)%idx2

      bondvec(jk) = coords(jatom) - coords(katom)

      !** Construct the script R's
      !**    script Rijk = 2*Sijk*Rjk
      !**    NOTE: scriptr is 2D rather than 3D because the possible jk
      !**          combinations are restricted by bond connectivity
      Do i = 1,rs%np
        scriptr(i,jk) = (2.0d0*rs%Sijk(i,j,k))*bondvec(jk)
      End Do

    End Do

    !** Construct the A matrix (using operator overloading here)
    !** Note: The Rjk vectors change between rows 
    !**       The i's in the scriptRijk vectors change between rows
    !**       The jk's in the scriptRijk vectors change between columns
    Do jk = 1,rs%ncon      !loop through the jk pairs  (i.e. the rows)
      Do i = 1,rs%ncon     !loop through the columns
        j = rs%con(jk)%idx1    
        k = rs%con(jk)%idx2
        A(jk,i) = bondvec(jk)*(scriptr(k,i) - scriptr(j,i))
      End Do
    End Do

    !** Construct the b vector
    !**   jk = the connection number in the data structure
    Do jk = 1,rs%ncon
      j = rs%con(jk)%idx1    
      k = rs%con(jk)%idx2
      jatom = rs%con(jk)%atom1    !these are ATOM indices
      katom = rs%con(jk)%atom2

      Fjk = scriptf(j) - scriptf(k)
      velocity2 = vector_getnormsq(v(katom) - v(jatom))
     
      !** Finally, evaluate the b-vector itself
      b(jk) = (-1.0d0)*(bondvec(jk)*Fjk + velocity2)
    End Do

    !** Solve the system of equations (solution x is written into b)
    Call matrixops_ludcmp(A,rs%ncon,rs%ncon,perm,nperm)
    Call matrixops_lubksb(A,rs%ncon,rs%ncon,perm,b)

    !** Finally, map the secondary forces onto the primary atoms
    Do i = 1,rs%np
      patom = rs%patom(i)
      vec = 0.0_RDbl
      Do jk = 1,rs%ncon
        vec = vec + scriptr(i,jk)*b(jk)
      End Do
     
!LC      Write(*,'(4x,i4,3e14.4)') i,vec
!LC      Fjk = (scriptf(i) - vec) - f(patom)    !additional mapped forces
!LC      Write(*,'(10x,a,i4,3e14.4)') 'Additional force on atom ',i,Fjk   

      f(patom) = scriptf(i) - vec
    End Do

#if DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'ciccotti_mapunitforces, output forces: '
    Do i = 1,Size(f)
      Write(*,*) i,vector_display(f(i),'f8.3')
    End Do                                                
    Stop
#endif

  End Subroutine ciccotti_mapunitforces


  !---------------------------------------------------------------
  ! Regenerate the positions and velocities of all the secondary 
  ! atoms in the molecule
  ! called with: rigidsetptr -- set of rigid units (WHOLE molec)
  !              coords -- one molecule's atomic position vectors
  !              v -- one molecule's velocity vectors
  !---------------------------------------------------------------
  Subroutine ciccotti_regenSet(rsptr,coords,v)
    Type(CiccottiRigidSetInfo), Pointer             :: rsptr
    Type(VecType), Dimension(:), Intent(InOut)      :: coords,v
 
    Integer       :: s,alpha,satom,patom,i

    Do s = 1,rsptr%nunits
      Do alpha = 1,rsptr%unit(s)%ns
        satom = rsptr%unit(s)%satom(alpha)
        coords(satom) = 0.0_RDbl
        v(satom) = 0.0_RDbl
      
        Do i = 1,rsptr%unit(s)%np
          patom = rsptr%unit(s)%patom(i)
          coords(satom) = coords(satom) + &
              rsptr%unit(s)%Calphai(alpha,i)*coords(patom)
          v(satom) = v(satom) + rsptr%unit(s)%Calphai(alpha,i)*v(patom)
        End Do

#if UCSTUFF      
         !*** Generate the positions of the atoms in the unit
         !*** cell from rxp, ryp and rzp coordinates.
         !*** We want the greatest integer less than rxp/ell(1)...
         sorbates%cxp(satom, molec) = floor(sorbates%rxp(satom,molec)/ell(1))
         sorbates%cyp(satom, molec) = floor(sorbates%ryp(satom,molec)/ell(2))
         sorbates%czp(satom, molec) = floor(sorbates%rzp(satom,molec)/ell(3))
         sorbates%rx(satom,molec) = sorbates%rxp(satom,molec) - 
     &                              sorbates%cxp(satom,molec)*ell(1)
         sorbates%ry(satom,molec) = sorbates%ryp(satom,molec) - 
     &                              sorbates%cyp(satom,molec)*ell(2)
         sorbates%rz(satom,molec) = sorbates%rzp(satom,molec) - 
     &                              sorbates%czp(satom,molec)*ell(3)
#endif
      
      End Do    !loop over secondary atoms
    End Do

!    Write(*,*) __FILE__,': Finished regenerating'

  End Subroutine ciccotti_regenSet


  !---------------------------------------------------------------
  ! Initialize the set information and construct matrices.  
  ! Returns False if the generation was unsuccesful.
  ! This allows the calling routine to act.
  ! called with:  rigidunitptr -- empty data type to contain info
  !               np -- number of primary atoms
  !               ns -- number of secondary atoms
  !               patom -- atom indices of the primary atoms
  !               satom -- atom indices of the secondary atoms
  !               atomcoords -- list of atomic position vectors
  !               invmasslist -- list of inverse atomic masses
  !               con -- a list of constraint connections 
  !               regencoords -- list of regenerated coordinates 
  !---------------------------------------------------------------
  Logical Function ciccotti_initRigidUnit(rigidunitptr,np,ns,patom, &
      satom,atomcoords,invmasslist,con,regencoords)
    Type(CiccottiRigidUnitInfo), Pointer       :: rigidunitptr
    Integer, Intent(In)                        :: np
    Integer, Intent(In)                        :: ns
    Integer, Dimension(:), Intent(In)          :: patom
    Integer, Dimension(:), Intent(In)          :: satom
    Type(VecType), Dimension(:), Intent(In)    :: atomcoords
    Real(kind=RDbl), Dimension(:), Intent(In)  :: invmasslist
    Type(CiccottiConnectInfo), Dimension(:), Intent(In)  :: con
    Type(VecType), Dimension(:), Intent(Out)   :: regencoords

    Integer              :: error,alpha,a
    Logical              :: success

    !** basic non-matrix assignments
    ciccotti_initRigidUnit = .True.
    rigidunitptr%patom = patom
    rigidunitptr%satom = satom
    rigidunitptr%np = np
    rigidunitptr%ns = ns
    rigidunitptr%natoms = np + ns
    rigidunitptr%ncon = Size(con)

    !** Allocate the primary atom inverse mass structure
    Allocate(rigidunitptr%patom_invmass(rigidunitptr%np), STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not allocate rigidunitptr%patom_invmass"
      Stop
    End If
    Do a = 1,rigidunitptr%np
      rigidunitptr%patom_invmass(a) = invmasslist(rigidunitptr%patom(a))
    End Do

    !** Allocate the constraint structure
    Allocate(rigidunitptr%con(rigidunitptr%ncon), STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not allocate rigidunitptr%con"
      Stop
    End If
    rigidunitptr%con = con

    !** Nullify the internal pointers
    Nullify(rigidunitptr%Calphai)
    Nullify(rigidunitptr%Aalphabeta)
    Nullify(rigidunitptr%Ainv)
    Nullify(rigidunitptr%Balphajk)
    Nullify(rigidunitptr%Sijk)

    If (rigidunitptr%ns == 0) Then  !**this is the zero unit
      Return
    End If

    !** allocate the Calphai matrix
    Allocate(rigidunitptr%Calphai(ns,np), STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not allocate rigidunitptr%Calphai"
      Stop
    End If

    !** make the Calphai matrix   
    rigidunitptr%Calphai = 0.0d0
    Do alpha = 1,ns
      success = ciccotti_makeCalphai(rigidunitptr,alpha, &
          atomcoords,regencoords)
      If (.Not. success) ciccotti_initRigidUnit = .False.
    End Do

    !** Allocate memory for the Aalphabeta matrix and its inverse
    Allocate(rigidunitptr%Aalphabeta(ns,ns), STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not allocate rigidunitptr%Aalphabeta"
      Stop
    End If
    Allocate(rigidunitptr%Ainv(ns,ns), STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not allocate rigidunitptr%Ainv"
      Stop
    End If

    !** make the Aalphabeta matrix and its inverse
    Call ciccotti_makeAalphabeta(rigidunitptr,invmasslist)

    !** allocate memory for the Balphajk matrix 
    Allocate(rigidunitptr%Balphajk(ns,np,np), STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not allocate rigidunitptr%Balphajk"
      Stop
    End If

    !** make the Balphajk matrices
    Call ciccotti_makeBalphajk(rigidunitptr,invmasslist)

    !** allocate memory for the Sijk matrix 
    Allocate(rigidunitptr%Sijk(np,np,np), STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not allocate rigidunitptr%Sijk"
      Stop
    End If

    !** make the Sijk matrices
    Call ciccotti_makeSijk(rigidunitptr,invmasslist)

    !** display the whole structure in its current state
!    Call ciccotti_displayunit(rigidunitptr,6,2)

  End Function ciccotti_initRigidUnit

  !--------------------------------------------------------------------------
  ! Build a single vector in the Calphai matrix.  This is the one that 
  ! allows generation of the secondary atoms based on the coordinates
  ! of the primary atoms.  Returns False if the generation was unsuccesful.
  ! This allows the calling routine to act.
  ! called with:  ptr -- pointer to rigid unit data structure
  !               alpha -- the secondary for which to build conversion vector
  !               atomcoords -- list of atomic position vectors
  !               newcoords -- regenerated coordinates (as a test)
  ! indices in Calphai run from 1->ns 
  !--------------------------------------------------------------------------
  Logical Function ciccotti_makeCalphai(ptr,alpha,atomcoords,newcoords)
    Type(CiccottiRigidUnitInfo), Pointer       :: ptr
    Integer, Intent(In)                        :: alpha
    Type(VecType), Dimension(:), Intent(In)    :: atomcoords
    Type(VecType), Dimension(:), Intent(Out)   :: newcoords

    Integer                                       :: ns,np,dim,i,j,aindex
    Type(VecType)                                 :: refcoord,alphacoord
    Type(VecType)                                 :: jcoord
    Real(kind=RDbl)                               :: nperm,norm
    Integer, Dimension(ptr%np-1)                  :: perm
    Real(kind=RDbl), Dimension(ptr%np-1)          :: b
    Real(kind=RDbl), Dimension(ptr%np-1,ptr%np-1) :: A
    Character(len=lstrLen)                        :: tempstrg

    ciccotti_makeCalphai = .True.
    ns = ptr%ns
    np = ptr%np
    dim = np - 1  !dimension of molecule

    refcoord = atomcoords(ptr%patom(1))
    alphacoord = atomcoords(ptr%satom(alpha))
    aindex = ptr%satom(alpha)

    !** set-up system of equations (Ax=b) and solve
    Do i = 1,dim
      b(i) = alphacoord%comp(i) - refcoord%comp(i)
      Do j = 2,(dim+1)
        jcoord = atomcoords(ptr%patom(j)) 
        A(i,(j-1)) = jcoord%comp(i) - refcoord%comp(i)
      End Do
    End Do

    !** Solve the system of equations
    Call matrixops_ludcmp(A,dim,dim,perm,nperm)
    Call matrixops_lubksb(A,dim,dim,perm,b)

    !** Translate the information from the solution vector into Calphai's
    ptr%Calphai(alpha,1) = 1.0_RDbl

    Do i = 2,(dim + 1)
      ptr%Calphai(alpha,1) = ptr%Calphai(alpha,1) - b(i-1)
      ptr%Calphai(alpha,i) = b(i-1)
    End Do

    !** Check the normalization constraint
    norm = Sum(ptr%Calphai(alpha,1:np))   !check this command on website
    If (Abs(norm - 1.0d0) > 1.0d-5) Then
      Write(0,'(2a,i4,a,i2)') __FILE__,":",__LINE__, &
          "problem with normalization of Calphai ",alpha
      Stop
    End If

    !** Regenerate the reference alpha atom coordinates as a check
    !** Problems with this regeneration can indicate that the constraints
    !** are poorly specified, for example, not divided well into the units
    newcoords(aindex) = 0.0_RDbl
    Do i = 1,np
      newcoords(aindex) = newcoords(aindex) + &
          ptr%Calphai(alpha,i)*atomcoords(ptr%patom(i))
    End Do
    norm = mag(newcoords(aindex) - alphacoord)

!LC    Write(0,'(2a,i4,a,i2,a,e10.4)') __FILE__,":",__LINE__, &
!LC        '  Secondary atom: ',alpha,'  coordinate regeneration error: ',norm

    If (Abs(norm) > 1.0d-3) Then  !NOTE: TOLERANCE HERE
      Write(0,'(2a,i4,a,2i3)') __FILE__,":",__LINE__, &
          " ERROR problem with regeneration of alpha coord ", &
          alpha,ptr%satom(alpha)
      tempstrg = vector_display(alphacoord,'f12.5')
      Write(0,'(2a)') 'Original:    ',Trim(tempstrg)
      tempstrg = vector_display(newcoords(aindex),'f12.5')
      Write(0,'(2a)') 'Regenerated: ',Trim(tempstrg)
      Write(0,'(a,f12.5)') 'Norm(diff):  ',norm
      Write(0,'(a)') 'Consider increasing the number of primary atoms, up to 4'
      ciccotti_makeCalphai = .False.      
      Return

    Else If (Abs(norm) > 1.0d-6) Then
      Write(0,'(2a,i4,a,i3)') __FILE__,":",__LINE__, &
          " WARNING potential problem with regeneration of alpha coord ", alpha
      Write(0,'(4x,a,i2,a,f8.5)') 'for atom ',ptr%satom(alpha), &
          ': abs(new_coord - old_coord) = ',norm
    End If

  End Function ciccotti_makeCalphai

  !--------------------------------------------------------------------------
  ! Make the Aalphabeta matrix and its inverse
  ! called with: ptr -- pointer to rigid unit data structure
  !              invmasslist -- list of inverse atomic masses
  ! indices in Aalphabeta run from 1->ns 
  !--------------------------------------------------------------------------
  Subroutine ciccotti_makeAalphabeta(ptr,invmasslist)

    Type(CiccottiRigidUnitInfo), Pointer       :: ptr
    Real(kind=RDbl), Dimension(:), Intent(In)  :: invmasslist

    Integer                                    :: ns,np,alpha,beta,i,j
    Real(kind=RDbl)                            :: nperm,massinv
    Integer, Dimension(ptr%ns)                 :: perm
    Real(kind=RDbl), Dimension(ptr%ns,ptr%ns)  :: A,Asave,b

    ns = ptr%ns
    np = ptr%np

    !** make the Aalphabeta matrix
    ptr%Aalphabeta = 0.0_RDbl
    Do alpha = 1,ns
      Do beta = 1,ns

!        Write(*,*) 'A alpha,beta ',alpha,beta

        Do i = 1,np
          massinv = invmasslist(ptr%patom(i))
          ptr%Aalphabeta(alpha,beta) = ptr%Aalphabeta(alpha,beta) &
              + ptr%Calphai(alpha,i)*ptr%Calphai(beta,i)*massinv

!        Write(*,*) i,ptr%Calphai(alpha,i),ptr%Calphai(beta,i),massinv
        End Do

        If (alpha == beta) Then
          massinv = invmasslist(ptr%satom(alpha))
          ptr%Aalphabeta(alpha,beta) = ptr%Aalphabeta(alpha,beta) + massinv
        End If
   
!       Write(*,*) 'A alpha,beta = ',ptr%Aalphabeta(alpha,beta)

      End Do
    End Do

    !** setup an identity matrix and copy the Aalphabeta matrix
    A = ptr%Aalphabeta
    Asave = ptr%Aalphabeta
    Do i = 1,ns
      Do j = 1,ns
        b(i,j) = 0
      End Do
      b(i,i) = 1
    End Do

    !** Perform LU decomposition on the Aalphabeta matrix
    Call matrixops_ludcmp(A,ns,ns,perm,nperm)

    !** Find the inverse by solving an equation for each column
    Do j = 1,ns
      Call matrixops_lubksb(A,ns,ns,perm,b(1:ns,j))
    End Do

    !** store the inverse 
    ptr%Ainv = b

  End Subroutine ciccotti_makeAalphabeta


  !--------------------------------------------------------------------------
  ! Make the Balphajk matrix
  ! called with: ptr -- pointer to rigid unit data structure
  !              atomcoords -- list of atomic position vectors
  !              invmasslist -- list of inverse atomic masses
  ! indices in Balphajk run from 1->ns (alpha) and 1->np (j,k)
  !--------------------------------------------------------------------------
  Subroutine ciccotti_makeBalphajk(ptr,invmasslist)

    Type(CiccottiRigidUnitInfo), Pointer       :: ptr
    Real(kind=RDbl), Dimension(:), Intent(In)  :: invmasslist

    Integer                                    :: alpha,j,k
    Real(kind=RDbl)                            :: invmassj,invmassk
    Character(strLen)                          :: form

    form = achar(iachar('1')+ptr%np-1)
    form = '(1x,2i4,'//trim(form)//'f8.3)'

    ptr%Balphajk = 0.0_RDbl
    
!    Write(*,'(1x,a)') 'The Balphajk matrix is:'
    
    Do alpha = 1,ptr%ns
       Do j = 1,ptr%np
          invmassj = invmasslist(ptr%patom(j))
    
          Do k = 1,ptr%np
             invmassk = invmasslist(ptr%patom(k))
             ptr%Balphajk(alpha,j,k) = ptr%Calphai(alpha,j)*invmassj &
                 - ptr%Calphai(alpha,k)*invmassk
          EndDo
    
!          Write(*,form) alpha,j,(ptr%Balphajk(alpha,j,k),k=1,np) 
       EndDo
    EndDo

  End Subroutine ciccotti_makeBalphajk


  !--------------------------------------------------------------------------
  ! Make the Sijk matrices
  ! called with: ptr -- pointer to rigid unit data structure
  !              invmasslist -- list of inverse atomic masses
  ! indices in Sijk run from 1->np
  !--------------------------------------------------------------------------
  Subroutine ciccotti_makeSijk(ptr,invmasslist)

    Type(CiccottiRigidUnitInfo), Pointer       :: ptr
    Real(kind=RDbl), Dimension(:), Intent(In)  :: invmasslist

    Integer                                    :: alpha,beta,i,j,k
    Real(kind=RDbl)                            :: invmass
!!$    Character(strLen)                          :: form

!!$    form = achar(iachar('1')+ptr%np-1)
!!$    form = '(1x,2i4,'//trim(form)//'f8.3)'

    ptr%Sijk = 0.0_RDbl
   
!!$   Write(*,'(1x,a)') 'The Sijk matrix is:'
   
    Do i = 1,ptr%np
      invmass = invmasslist(ptr%patom(i))
   
      Do j = 1,ptr%np
        Do k = 1,ptr%np

          If (i == j) Then
            ptr%Sijk(i,j,k) = ptr%Sijk(i,j,k) + 1.0d0
          End If
   
          If (i == k) Then
            ptr%Sijk(i,j,k) = ptr%Sijk(i,j,k) - 1.0d0
          End If
                  
          Do alpha = 1,ptr%ns
            Do beta = 1,ptr%ns
              ptr%Sijk(i,j,k) = ptr%Sijk(i,j,k) &
                  - ptr%Calphai(alpha,i)*ptr%Ainv(alpha,beta)* &
                  ptr%Balphajk(beta,j,k)
            End Do
          End Do
   
          ptr%Sijk(i,j,k) = ptr%Sijk(i,j,k)*invmass
   
        End Do
!        Write(*,form) i,j,(ptr%Sijk(i,j,k),k=1,ptr%np) 
      End Do
    End Do

  End Subroutine ciccotti_makeSijk

  !-------------------------------------------------------------------------
  ! Visualize the rigid unit information for a molecule
  ! called with: rigidset -- rigid set pointer
  !              atomcoords -- position vector array indexed by atom number
  !              units -- flag indicated if the units should be separated
  !-------------------------------------------------------------------------
  Subroutine ciccotti_visunits(rigidsetptr,atomcoords,units)
    Type(CiccottiRigidSetInfo), Pointer        :: rigidsetptr
    Type(VecType), Dimension(:), Intent(In)    :: atomcoords
    Logical, Intent(In)                        :: units

    Integer                                    :: u,a,n,disp,atom
    Character(len=strLen)                      :: filename
    Type(XYZ_Entry), Dimension(rigidsetptr%natoms)  :: entries

    !** Dump either the straight coordinates labeled as C or separated in units
    If (units) Then
      n = 0
      disp = 0
      Do u = 0,rigidsetptr%nunits
        disp = disp + 1
        Do a = 1,rigidsetptr%unit(u)%np
          n = n + 1
          atom = rigidsetptr%unit(u)%patom(a)
          entries(n) = visxyz_make(atomcoords(atom),'',disp)
        End Do
        disp = disp + 1
        Do a = 1,rigidsetptr%unit(u)%ns
          n = n + 1
          atom = rigidsetptr%unit(u)%satom(a)
          entries(n) = visxyz_make(atomcoords(atom),'',disp)
        End Do
      End Do
    Else
      Do a = 1,rigidsetptr%natoms
        entries(a) = visxyz_make(atomcoords(a),'C ')
      End Do
    End If

    filename = 'units.xyz'
    Write(*,'(3a)') __FILE__,': Dumping rigid units visualization to: ', &
        Trim(filename)
    Call visxyz_dump(entries,filename,'f14.8','Ciccotti Rigid Units')

  End Subroutine ciccotti_visunits


  !---------------------------------------------------------------
  ! Display the rigid set information
  ! called with: rigidset -- the appropriate data type for a SET
  !---------------------------------------------------------------
  Subroutine ciccotti_displayset(rigidset,unit,indent)
    Type(CiccottiRigidSetInfo), Intent(In)     :: rigidset
    Integer, Intent(In)                        :: unit, indent

    Integer                     :: i,s
    Character(len=indent)       :: blank

    blank = ""
    Do i = 1, indent
      blank = blank//" "
    End Do

    Write(unit,'(2a)') blank,"Ciccotti Rigid Set structure information:"
    Write(unit,'(2a,i2)') blank," number of rigid units: ",rigidset%nunits
    Write(unit,'(2a,i2)') blank," number of atoms: ",rigidset%natoms
    Write(unit,'(2a,i2)') blank," number of mobile atoms: ",rigidset%nmobile
    Write(unit,'(2a,i2)') blank," number of primary atoms: ",rigidset%nprimary
    Write(unit,'(2a,i2)') blank," number of secondary atoms: ", &
        rigidset%nsecondary

    If (Associated(rigidset%unit)) Then
      Do s = 1,rigidset%nunits
        Write(unit,'(2a,i2,a)') blank," Rigid Unit ",s," information:"
        Call ciccotti_displayunit(rigidset%unit(s),unit,indent+3)
      End Do
    Else
      Write(unit,'(2a)') blank,"Individual units not associated "
    End If

  End Subroutine ciccotti_displayset

  !---------------------------------------------------------------
  ! Display the rigid UNIT information
  ! called with: rigidunit -- the appropriate data type for a UNIT
  !---------------------------------------------------------------
  Subroutine ciccotti_displayunit(rigidunit,unit,indent)
    Type(CiccottiRigidUnitInfo), Intent(In)     :: rigidunit
    Integer, Intent(In)                        :: unit, indent

    Integer                     :: i,width,a,alpha
    Character(len=indent)       :: blank
    Character(len=strLen)       :: display,nmbr,string

    blank = ""
    Do i = 1, indent
      blank = blank//" "
    End Do

    Write(unit,'(2a)') blank,"Ciccotti Rigid Unit information:"

    !** Display the basics
    Write(unit,'(2a,i2)') blank," number of atoms: ",rigidunit%natoms
    Write(unit,'(2a,i2)') blank," number of primary atoms: ",rigidunit%np
    width = rigidunit%np
    Write(nmbr,'(i6)') width

    display = "(2a,"//Trim(Adjustl(nmbr))//"i3)"
    Write(unit,display) blank," primary atoms: ",rigidunit%patom(1:width)
    Write(unit,'(2a,i2)') blank," number of secondary atoms: ",rigidunit%ns
    width = rigidunit%ns
    Write(nmbr,'(i6)') width

    display = "(2a,"//Trim(Adjustl(nmbr))//"i3)"
    Write(unit,display) blank," secondary atoms: ",rigidunit%satom(1:width)

    If (Trim(genparams%displaymode) == "VERBOSE") Then
      !** Display the Calphai matrix
      If(Associated(rigidunit%Calphai)) Then
        Write(unit,'(2a)') blank," Calphai Matrix: "
        width = Size(rigidunit%Calphai,2)
        Write(nmbr,'(i6)') width
        display = "(a,"//Trim(nmbr)//"f8.3)"
        Do a = 1,rigidunit%ns
          Write(unit,'(2a,i2,a,i2)') blank," Calphai for secondary atom: ", &
              a," atomno: ",rigidunit%satom(a)
          Write(unit,display) blank,rigidunit%Calphai(a,1:width)
        End Do
      Else
        Write(unit,'(2a)') blank," Calphai matrix not associated "
      End If

      !** Display the Aalphabeta matrix
      If (Associated(rigidunit%Aalphabeta)) Then
        Write(unit,'(2a)') blank," Aalphabeta Matrix: "
        width = Size(rigidunit%Aalphabeta,2)
        Write(nmbr,'(i6)') width
        display = "(a,"//Trim(nmbr)//"f8.3)"
        Do a = 1,rigidunit%ns
          Write(unit,display) blank,rigidunit%Aalphabeta(a,1:width)
        End Do
      Else
        Write(unit,'(2a)') blank," Aalphabeta matrix not associated "
      End If
  
      !** Display the Ainv matrix
      If (Associated(rigidunit%Ainv)) Then
        Write(unit,'(2a)') blank," Ainv Matrix: "
        width = Size(rigidunit%Ainv,2)
        Write(nmbr,'(i6)') width
        display = "(a,"//Trim(nmbr)//"f8.3)"
        Do a = 1,rigidunit%ns
          Write(unit,display) blank,rigidunit%Ainv(a,1:width)
        End Do
      Else
        Write(unit,'(2a)') blank," Ainv matrix not associated "
      End If
  
      !** Display the Balphajk matrix
      If (Associated(rigidunit%Balphajk)) Then
        Do alpha = 1,rigidunit%ns
          Write(unit,'(2a,i2)') blank," Balphajk Matrix for alpha= ",alpha
          width = Size(rigidunit%Balphajk,3)
          Write(nmbr,'(i6)') width
          display = "(a,"//Trim(nmbr)//"f8.3)"
          Do a = 1,rigidunit%np
            Write(unit,display) blank,rigidunit%Balphajk(alpha,a,1:width)
          End Do
        End Do
      Else
        Write(unit,'(2a)') blank," Balphajk matrix not associated "
      End If
  
      !** Display the Sijk matrix
      If (Associated(rigidunit%Sijk)) Then
        Do i = 1,rigidunit%np
          string = int2str(i)
          Write(unit,'(3a)') blank," Sijk Matrix for i= ",Trim(string)
          display = "(a,"//Trim(int2str(Size(rigidunit%Sijk,3)))//"f8.3)"
          Do a = 1,rigidunit%np
            Write(unit,display) blank,rigidunit%Sijk(i,a,1:width)
          End Do
        End Do
      Else
        Write(unit,'(2a)') blank," Sijk matrix not associated "
      End If
  
      Write(unit,'(2a,i2)') blank," number of constraints: ",rigidunit%ncon
      If (Associated(rigidunit%con)) Then
        Do i = 1,rigidunit%ncon
          Write(unit,'(2a,2(i3,a,i3,a),f10.3)') blank,"Atom pair distance: ",& 
              rigidunit%con(i)%atom1,"(",rigidunit%con(i)%idx1,")", &
              rigidunit%con(i)%atom2,"(",rigidunit%con(i)%idx2,")", &
              rigidunit%con(i)%length
        End Do
      Else
        Write(unit,'(2a)') blank," Constraints not associated "
      End If
    End If

  End Subroutine ciccotti_displayunit


  !----------------------------------------------------------------------------
  ! Cleanup the data structure
  !----------------------------------------------------------------------------
  Subroutine ciccotti_cleanup(rigidset)
    Type(CiccottiRigidSetInfo)               :: rigidset
    Integer                                  :: error,i

    Do i = 1,rigidset%nunits
      Deallocate(rigidset%unit(i)%Calphai,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'rigidset%unit%Calphai'"
        Stop
      End If

      Deallocate(rigidset%unit(i)%Aalphabeta,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'rigidset%unit%Aalphabeta'"
        Stop
      End If

      Deallocate(rigidset%unit(i)%Sijk,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'rigidset%unit%Sijk'"
        Stop
      End If

      Deallocate(rigidset%unit(i)%Ainv,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'rigidset%unit%Ainv'"
        Stop
      End If

      Deallocate(rigidset%unit(i)%Balphajk,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'rigidset%unit%Balphajk'"
        Stop
      End If

      Deallocate(rigidset%unit,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'rigidset%unit'"
        Stop
      End If

    End Do

  End Subroutine ciccotti_cleanup


End Module ciccotti
