!------------------------------------------------------------------------------
! This module contains the data structure defining the generalized coordinates 
! of a rigid molecule.  Namely, the 3 center of mass coordinates and the 3
! eulerian angles.  It also includes the reference configuration definition.
! If the corresponding logical flag is positive, the individual molecule 
! will be treated as a rigid body with internal degrees of freedom.  If the
! molecule is internally flexible, then the necessary reference configuration
! will be taken from the molecules structure for the specific species type
! during initialization.  All outputs reflect the possibility of internal
! flexibility, for example, the reference configuration is dumped to the
! restartfile if and only if the molecule is internally flexible.
! 
! Important routines: 
!   rigidcoords_init -- initializes structure for one molecule
!   rigidcoords_toxyz -- generates toxyz coordinates based from the gen. coords
! 
! Needed Improvements:
! 1) _fromxyz and _getmetrictensor will not work properly for flexible case
! 2) _getgenforce may not work, should be checked.
!------------------------------------------------------------------------------
Module rigidcoords

  Use defaults, Only: RDbl, strLen, lstrLen,twopi,zerotolerance,MAX_DOF,zero, dbgflag
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toreal,&
      getinvangle, allocerrdisplay, ToUpper, deallocerrdisplay
  Use vector, Only: VecType, vector_crossprod, vector_display, Assignment(=), &
      Operator(+), Operator(-), Operator(*), Operator(/), vector_getnorm
  Use matrix, Only: MatrixType, matrix_b2feuler, matrix_Db2fDTheta, &
      matrix_Db2fDPhi, matrix_Db2fDPsi, Operator(*), Assignment(=)
  Use simcell, Only: SimCell_Params,simcell_pbc
  Use molecules, Only: molecules_getbodyaxes, molecules_getdof, &
      molecules_getmass, molecules_getmoi, molecules_getnthatom, &
      molecules_getcomdefn, molecules_getnatoms, molecules_getcom
  Use localcoords, Only: LocalCoordInfo, localcoords_getaxes

  Implicit None
  Save

  Private
  Public :: RigidMolec, rigidcoords_displaystr, rigidcoords_getcom, & 
      rigidcoords_getgencoords, rigidcoords_initcoords, rigidcoords_clean, &
      rigidcoords_gcoords_eq_vectype, rigidcoords_toxyz, rigidcoords_fromxyz, &
      rigidcoords_setgencoords, rigidcoords_getmetrictensor, &
      rigidcoords_pbc, rigidcoords_readrestartfile, rigidcoords_changeflex, &
      rigidcoords_dumprestartinfo, rigidcoords_convstr, rigidcoords_getref, &
      rigidcoords_getgenforce, rigidcoords_copy, rigidcoords_changeref, &
      rigidcoords_setfrom_rp, rigidcoords_changedof

  Type RigidMolec
    Logical           :: internally_flexible  
    Integer           :: dof                           !* degrees of freedom
    Type(VecType)     :: com                           !* center of mass
    Real(kind=RDbl)   :: theta                         !* Range 0 to pi
    Real(kind=RDbl)   :: phi, psi                      !* Range 0 to 2*pi
    Type(VecType), Dimension(:), Pointer  :: ref_struc !* reference structure
  End Type RigidMolec
    
Contains
  !-------------------------------------------------------------------------
  ! Initializes the model params, with or without optional reference 
  ! structure.  If the optional structure isn't present it will get the
  ! reference structure from molecules using the species number.
  ! Requires: gcoords -- generalized rigid coordinates for one molecule
  !           spc -- species type
  !           ref_struc -- optional reference structure
  !-------------------------------------------------------------------------
  Subroutine rigidcoords_initcoords(gcoords,spc,ref_struc)
    Type(RigidMolec), Pointer              :: gcoords
    Integer, Intent(In)                    :: spc
    Type(VecType), Dimension(:), Optional  :: ref_struc

    Integer        :: error, dof, natoms

    dof = molecules_getdof(spc)

    Allocate(gcoords, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'gcoords')

    gcoords%com = 0.0_RDbl
    gcoords%theta = 0.0_RDbl
    gcoords%phi = 0.0_RDbl
    gcoords%psi = 0.0_RDbl
    gcoords%dof = dof

    Nullify(gcoords%ref_struc)

    !** allocate memory for reference structure
    natoms = molecules_getnatoms(spc)
    Allocate(gcoords%ref_struc(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ref_struc')

    If (Present(ref_struc)) Then
      !** check number of atoms
      natoms = molecules_getnatoms(spc)
      If (Size(ref_struc) /= natoms) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
            " number of atoms between molecule defn and passed defn do not match"
        Stop
      End If
      gcoords%internally_flexible = .True.
      gcoords%ref_struc = ref_struc

      If (gcoords%dof <= 6) Then
        If (.Not.(natoms==2)) Then 
          Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
              " ERROR: DOF <= 6, but the molecule is internally flexible"
          Stop
        Endif
      End If

    Else
      gcoords%internally_flexible = .False.
      Call molecules_getcomdefn(spc,gcoords%ref_struc)

      If (gcoords%dof > 6) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
            " ERROR: DOF is greater than 6, but the molecule is fully rigid"
        Write(0,'(5x,a)') 'Consider using the "None" generalized coordinates'
        Write(0,'(5x,a)') 'if only molecular dynamics is being done'
        Stop
      End If
    End If

  End Subroutine rigidcoords_initcoords

  !-------------------------------------------------------------------------
  ! Changes the internal flexibility flag 
  ! Requires: gcoords -- generalized rigid coordinates for one molecule
  !           flexflag -- new value for internal flexibility
  !-------------------------------------------------------------------------
  Subroutine rigidcoords_changeflex(gcoords,flexflag)
    Type(RigidMolec), Pointer              :: gcoords
    Logical, Intent(In)                    :: flexflag

    gcoords%internally_flexible = flexflag

  End Subroutine rigidcoords_changeflex


  !-------------------------------------------------------------------------
  ! Changes the dof
  ! Requires: gcoords -- generalized rigid coordinates for one molecule
  !           dof -- new value for degrees of freedom 
  !-------------------------------------------------------------------------
  Subroutine rigidcoords_changedof(gcoords,dof)
    Type(RigidMolec), Pointer              :: gcoords
    Integer, Intent(In)                    :: dof 
    gcoords%dof = dof
  End Subroutine rigidcoords_changedof
  

  !-------------------------------------------------------------------------
  ! Set the center-of-mass of the generalized coordinates
  ! Requires: gcoords -- generalized rigid coordinates for one molecule
  !           vec -- center-of-mass
  !-------------------------------------------------------------------------
  Subroutine rigidcoords_gcoords_eq_vectype(gcoords, vec)
    Type(RigidMolec), Pointer       :: gcoords
    Type(VecType), Intent(in)       :: vec
    
    gcoords%com = vec
  End Subroutine rigidcoords_gcoords_eq_vectype

  !-------------------------------------------------------------------------
  ! Copies the fields of rigidcoords
  ! Requires: to_coord -- destination generalized rigid coordinates 
  !           from_coord -- origin generalized rigid coordinates 
  !-------------------------------------------------------------------------
  Subroutine rigidcoords_copy(to_coord, from_coord)
    Type(RigidMolec), Pointer       :: from_coord, to_coord

    If (Associated(to_coord)) Then
      to_coord%internally_flexible = from_coord%internally_flexible
      to_coord%com = from_coord%com
      to_coord%theta = from_coord%theta
      to_coord%psi = from_coord%psi
      to_coord%phi = from_coord%phi
      to_coord%dof = from_coord%dof
      to_coord%ref_struc = from_coord%ref_struc

#if REMOVE_LATER
      !** copy the reference structure too, if it's there
      If (Associated(from_coord%ref_struc)) Then
        If (Associated(to_coord%ref_struc)) Then
          !** shouldn't need this, could remove for speed
          If (Size(to_coord%refstruc) /= Size(to_coord%ref_struc)) Then 
            Deallocate(to_coord%refstruc, STAT=error)
            If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ref_struc')
            Allocate(to_coord%ref_struc(Size(from_coord%ref_struc)), STAT=error)
            If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ref_struc')
          End If
        Else 
          Allocate(to_coord%ref_struc(Size(from_coord%ref_struc)), STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ref_struc')
        End If
        to_coord%ref_struc = from_coord%ref_struc
      Else
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Reference structure is not associated, it should be"
        Stop
      End If
#endif

    Else
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " New structure must first be initialized"
      Stop
    End If

  End Subroutine rigidcoords_copy

  !----------------------------------------------------------------------------
  ! Change the reference structure.  This is what makes this "coordinate
  ! system" internally flexible.
  ! Requires: gcoords -- generalized rigid coordinates 
  !           ref_struc -- new reference structure
  ! NOTE: because there are no constraints on this new reference structure,
  ! there is no guarantee that the new reference structure's COM will really
  ! be at the origin.
  !----------------------------------------------------------------------------
  Subroutine rigidcoords_changeref(gcoords,ref_struc)
    Type(RigidMolec), Pointer                 :: gcoords
    Type(VecType), Dimension(:), Intent(In)   :: ref_struc

    Integer         :: natoms

    natoms = Size(gcoords%ref_struc)

    gcoords%ref_struc = ref_struc(1:natoms)

  End Subroutine rigidcoords_changeref

  !----------------------------------------------------------------------------
  ! Get the  reference structure.  
  ! Requires: gcoords -- generalized rigid coordinates 
  !           ref_struc -- returned reference structure
  !----------------------------------------------------------------------------
  Subroutine rigidcoords_getref(gcoords,ref_struc)
    Type(RigidMolec), Pointer                 :: gcoords
    Type(VecType), Dimension(:), Intent(Out)  :: ref_struc

    Integer         :: natoms

    natoms = Size(gcoords%ref_struc)
    ref_struc(1:natoms) = gcoords%ref_struc

  End Subroutine rigidcoords_getref

  !----------------------------------------------------------------------------
  ! Generate the xyz coordinates from the generalized coordinates.
  ! Requires: gcoords -- generalized rigid coordinates 
  !           xyzcoords -- output xyz coordinates
  !----------------------------------------------------------------------------
  Subroutine rigidcoords_toxyz(gcoords,xyzcoords)
    Type(RigidMolec), Pointer                 :: gcoords
    Type(VecType), Dimension(:), Intent(Out)  :: xyzcoords

    Type(MatrixType)   :: b2f_eulermat
    Integer            :: natoms, i
    
    !** Check the no. of atoms and return if there is only one
    natoms = Size(gcoords%ref_struc,1)
    If (natoms == 1) Then
      xyzcoords(1) = gcoords%com
      Return
    End If

    !** Get the transformation matrix
    b2f_eulermat = matrix_b2feuler(gcoords%phi, gcoords%theta, gcoords%psi)

    !** Do the transformation
    Do i = 1,natoms
      xyzcoords(i) = gcoords%com + b2f_eulermat * gcoords%ref_struc(i)
    End Do

  End Subroutine rigidcoords_toxyz

  !----------------------------------------------------------------------------
  ! This routine takes an array of xyz coordinates (xyzcoords) and generates 
  ! generalized coordinates from them using the previously defined coordinate
  ! system set up in the molecules structure.  
  ! Requires: gcoords -- output generalized rigid coordinates 
  !           xyzcoords -- array of xyz coordinates
  !           spc -- species number
  ! Needed Improvements:
  ! 1) only works when the reference structure matches the molecules structure
  !----------------------------------------------------------------------------
  Subroutine rigidcoords_fromxyz(gcoords, xyzcoords, spc)
    Type(RigidMolec), Pointer               :: gcoords
    Type(VecType), Dimension(:), Intent(In) :: xyzcoords
    Integer, Intent(In)                     :: spc

    Real(kind=RDbl)              :: phi, sinphi, cosphi
    Real(kind=RDbl)              :: theta, costheta, sintheta
    Real(kind=RDbl)              :: psi, sinpsi, cospsi
    Type(VecType)                :: rcom, xf, yf, zf, e1, e2, e3
    Type(LocalCoordInfo)         :: lcoords

    If (gcoords%internally_flexible) Then
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Routine does not work for internally flexible molecules"
      Stop    
    End If

    !** Get the center of mass from the xyz coordinates
    rcom = molecules_getcom(spc,xyzcoords)
    !LC    Write(*,*) 'Center of mass: ',vector_display(rcom,'f8.3')
    If (molecules_getnatoms(spc)<3) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "diatomic/monoatomic molecule,contact shaji"
      Write(*,*) "shaji@northwestern.edu"
      stop
    Endif

    !** Generate in the fixed coordinate system what the position of 
    !** the body axes would be.  The coordinates of the body axes
    !** as expressed in the fixed coordinates are given by "xf, yf, zf"
    lcoords = molecules_getbodyaxes(spc)
    Call localcoords_getaxes(lcoords,xyzcoords,rcom,xf,yf,zf)

    !** July2001, I don't understand why xf,yf,zf aren't normalized
    !** since the definition of the cos (= AdotB/(|A||B|)) that is used 
    !** below assumes it.  I'm going to normalize them here -LC, works now
    xf = xf/vector_getnorm(xf)
    yf = yf/vector_getnorm(yf)
    zf = zf/vector_getnorm(zf)

    !** Now generate cos(theta), sin(theta) by looking at the zth component
    !** of zf.  This follows from the fact that the element matrix_b2f(3,3)
    !** is costheta.  So, the body axis (0,0,1) transforms to (*,*,costheta)
    !** where the "*" could be any number.
    e1 = (/1.0_Rdbl, 0.0_RDbl, 0.0_Rdbl/)
    e2 = (/0.0_Rdbl, 1.0_RDbl, 0.0_Rdbl/)
    e3 = (/0.0_Rdbl, 0.0_RDbl, 1.0_Rdbl/)
    costheta = e3*zf
    theta    = Acos(costheta)      ! theta is between 0 and pi
    sintheta = Sin(theta)

    If (Abs(sintheta) > zeroTolerance) Then
      !** Similarly for phi.  The xth component of zf gives "sintheta*sinphi"
      !** and the yth component of zf gives "-sintheta*cosphi"
      sinphi = (e1*zf)/sintheta
      cosphi = (e2*zf)/(-sintheta)
      Call getinvangle(sinphi, cosphi, phi, .False.)

      !** Finally, the zth component of xf corresponds to "sinpsi*sintheta"
      !** and the zth component of yf corresponds to "cospsi*sintheta"
      sinpsi = (e3*xf)/sintheta
      cospsi = (e3*yf)/sintheta
      If (dbgflag) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) sinpsi, cospsi, sintheta
        Write(*,*) e3, xf, yf

      Endif

      Call getinvangle(sinpsi, cospsi, psi, .False.)

    Else
      !** Need an alternate way to find the eulerian angles since sintheta is
      !** zero.  In this case the rotations are simply two rotation about
      !** the z-axis.  Hence, it is similar to doing one rotation with an
      !** angle of phi+psi about z-axis.  So, we will set the angle "phi" to
      !** be that single rotation and "psi" to zero
      If (costheta > 0.0_RDbl) Then
        cosphi =  (e1*xf)
        sinphi = -(e1*yf)
        Call getinvangle(sinphi, cosphi, phi, .False.)
        psi = 0.0_RDbl
      Else
        cosphi = (e1*xf)
        sinphi = (e1*yf)
        Call getinvangle(sinphi, cosphi, phi, .False.)
        psi = 0.0_RDbl
      End If
    End If

    gcoords%com = rcom
    gcoords%theta = theta
    gcoords%phi   = phi
    gcoords%psi   = psi

  End Subroutine rigidcoords_fromxyz

  !---------------------------------------------------------------------------
  ! Gets the center-of-mass of the molecule described by "gcoords"
  ! Requires: gcoords -- output generalized rigid coordinates 
  !---------------------------------------------------------------------------
  Type(VecType) Function rigidcoords_getcom(gcoords)
    Type(RigidMolec), Pointer   :: gcoords

    rigidcoords_getcom = gcoords%com
  End Function rigidcoords_getcom

  !---------------------------------------------------------------------------
  ! Return the generalized coordinates as an array
  ! Requires:  gcoords -- output generalized rigid coordinates 
  !            refstruc -- optional returned reference structure
  !---------------------------------------------------------------------------
  Subroutine rigidcoords_getgencoords(gcoords,coords,refstruc)
    Type(RigidMolec), Pointer                           :: gcoords
    Real(kind=RDbl), Dimension(6), Intent(Out)          :: coords
    Type(VecType), Dimension(:), Intent(Out), Optional  :: refstruc

    Integer       :: natoms

    coords(1:3) = gcoords%com
    coords(4:6) = (/gcoords%theta, gcoords%phi, gcoords%psi/)

    !** Return the reference structure if necessary
    If (Present(refstruc)) Then
      natoms = Size(gcoords%ref_struc)
      If (natoms > Size(refstruc)) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
            ' Passed array is too small to contain reference structure'
        Stop
      End If
      refstruc(1:natoms) = gcoords%ref_struc
    End If

  End Subroutine rigidcoords_getgencoords

  !---------------------------------------------------------------------
  ! Set the generalized coordinates from the passed array.  Any number
  ! of array values beyond 3 is optional, see definition below.
  ! Requires: gcoords -- single rigid molecule generalized coords
  !           genarray -- flattened array containing new coords
  !             genarray(1:3) -- new center of mass (x,y,z)
  !             genarray(3:6) -- theta,phi,psi
  !             genarray(7) -- degrees of freedom
  !             genarray(8:8+natoms*3) -- used to set reference coords
  !---------------------------------------------------------------------
  Subroutine rigidcoords_setgencoords(gcoords,genarray)
    Type(RigidMolec), Pointer                  :: gcoords
    Real(kind=RDbl), Dimension(:), Intent(In)  :: genarray

    Integer            :: n,a,natoms,arraysize

    arraysize = Size(genarray)
    Select Case (arraysize)
    Case (0:2)
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " 'genarray' has too few values"
      Stop
    Case (3)
      gcoords%com = genarray(1:3)
    Case (4:6)
      gcoords%com   = genarray(1:3)
      gcoords%theta = genarray(4)
      gcoords%phi   = genarray(5)
      gcoords%psi   = genarray(6)
    Case (7)
      gcoords%com   = genarray(1:3)
      gcoords%theta = genarray(4)
      gcoords%phi   = genarray(5)
      gcoords%psi   = genarray(6)
      gcoords%dof   = Int(genarray(7))
    Case (8:)
      gcoords%com   = genarray(1:3)
      gcoords%theta = genarray(4)
      gcoords%phi   = genarray(5)
      gcoords%psi   = genarray(6)
      gcoords%dof   = Int(genarray(7))

      !** Just return the 1st seven for now
      Return

      Write(0,'(1x,2a,i4,a,i4)') __FILE__," : ",__LINE__, &
          " Size(genarray) >= 8 has not been tested ",arraysize
      Stop

      n = arraysize - 8
      If (Mod(arraysize,3) /= 0) Then
        Write(0,'(1x,2a,i4)') __FILE__,' : ',__LINE__, &
            ' size beyond first 8 should be divisible by 3'
        Stop
      End If

      natoms = n/3
      Do a = 1,natoms,2
        Write(*,*) a,a+2
        gcoords%ref_struc(a) = genarray(a:a+2)
      End Do
      
    Case Default
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " could not understand passed 'genarray'"
      Stop
    End Select

  End Subroutine rigidcoords_setgencoords

  !---------------------------------------------------------------------
  ! Returns the forces in generalized coordinates "genforce"
  ! given the forces in xyz coordinates on each atom in "xyzforce"
  ! Requires: gcoords -- single rigid molecule generalized coords
  !           xyzforce -- forces on each atom
  !           genforce -- forces generalized to the 6 rigid DOF
  !---------------------------------------------------------------------
  Subroutine rigidcoords_getGenForce(gcoords, xyzforce, genforce)
    Type(RigidMolec), Pointer                  :: gcoords
    Type(VecType), Dimension(:), Intent(In)    :: xyzforce
    Real(kind=RDbl), Dimension(:), Intent(Out) :: genforce

    Integer                        :: natoms, i, comcoord, dof
    Real(kind=RDbl)                :: phi, theta, psi
    Real(kind=RDbl), Dimension(3)  :: f
    Type(MatrixType)               :: Db2fDTheta, Db2fDPhi, Db2fDPsi

    If (gcoords%internally_flexible) Then
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Routine might not work for internally flexible molecules, check"
      Stop    
    End If
    
    !** Initialize some stuff
    natoms = Size(xyzforce,1)
    dof    = gcoords%dof
    genforce = 0.0_RDbl

    !** Get the center of mass derivatives
    Do i=1, natoms
      f = xyzforce(i)
      Do comcoord = 1, 3
        genforce(comcoord) = genforce(comcoord) + f(comcoord)
      End Do
    End Do

    !** If the molecule is spherical we are done
    If (dof == 3) Then
      Return
    End If

    !** Get the derivatives wrt to the Eulerian angles
    phi   = gcoords%phi
    theta = gcoords%theta
    psi   = gcoords%psi
    Db2fDTheta = matrix_Db2fDTheta(phi, theta, psi)
    Do i=1, natoms
      genforce(4) = genforce(4) + xyzforce(i)*(Db2fDTheta*gcoords%ref_struc(i))
    End Do

    !** Get the derivatives wrt to the Eulerian angles
    Db2fDPhi = matrix_Db2fDPhi(phi, theta, psi)
    Do i=1, natoms
      genforce(5) = genforce(5) + xyzforce(i)*(Db2fDPhi*gcoords%ref_struc(i))
    End Do

    !** Get the derivatives wrt to the Eulerian angles
    Db2fDPsi = matrix_Db2fDPsi(phi, theta, psi)
    Do i=1, natoms
      genforce(6) = genforce(6) + xyzforce(i)*(Db2fDPsi*gcoords%ref_struc(i))
    End Do
    
  End Subroutine rigidcoords_getGenForce

  !---------------------------------------------------------------------------
  ! Returns the metric tensor for the sorbate "sorbtype".  The metric tensor
  ! is defined by Randy in J. Phys. Chem., Vol. 98, 1994, pp. 11948--11961
  ! "ierr" is non-zero when metric tensor is singular which happens when 
  ! theta is zero.
  ! Requires: gcoords -- single rigid molecule generalized coords
  !           spc -- species number
  !           metricTensor -- the metric tensor
  !           ierr -- error integer flag
  !---------------------------------------------------------------------------
  Subroutine rigidcoords_getMetricTensor(gcoords, spc, metricTensor, ierr)
    Type(RigidMolec), Pointer                    :: gcoords
    Integer, Intent(In)                          :: spc
    Real(kind=RDbl), Dimension(:,:), Intent(Out) :: metricTensor
    Integer, Intent(Out)                         :: ierr

    Integer              :: i, dof
    Real(kind=RDbl)      :: mass
    Real(kind=Rdbl)      :: sintheta, costheta, sinphi, cosphi, sinpsi, cospsi
    Real(kind=RDbl), Dimension(3) :: pMoI

    If (gcoords%internally_flexible) Then
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Routine probably is not correct for internally flexible molecules"
      Stop    
    End If

    dof = molecules_getdof(spc)
    
    !** For a spherical molecule it is just the mass multiplied
    !** by a unit matrix
    metricTensor = 0.0_Rdbl
    If (dof == 3) Then
      mass = molecules_getmass(spc)
      Do i=1, dof
        metricTensor(i,i) = mass
      End Do
    Else
      mass = molecules_getmass(spc)
      pMoi = molecules_getmoi(spc)
      Do i=1, 3
        metricTensor(i,i) = mass
      End Do
      sintheta = Sin(gcoords%theta)
      costheta = Cos(gcoords%theta)
      sinphi   = Sin(gcoords%phi)
      cosphi   = Cos(gcoords%phi)
      sinpsi   = Sin(gcoords%psi)
      cospsi   = Cos(gcoords%psi)
      !** Check if theta is close to zero
      If (Abs(sintheta) < 1.0e-4_RDbl) Then
        ierr = 1
        Return
      End If

      !** theta
      metricTensor(4,4) = pMoI(1)*cospsi*cospsi + pMoI(2)*sinpsi*sinpsi
      metricTensor(4,5) = pMoI(1)*sintheta*sinpsi*cospsi &
                        - pMoI(2)*sintheta*sinpsi*cospsi
      metricTensor(4,6) = 0.0_RDbl
      !** phi              
      metricTensor(5,4) = metricTensor(4,5)
      metricTensor(5,5) = pMoI(1)*sintheta*sintheta*sinpsi*sinpsi &
                        + pMoI(2)*sintheta*sintheta*cospsi*cospsi &
                        + pMoI(3)*costheta*costheta
      metricTensor(5,6) = pMoI(3)*costheta
      !** psi              
      metricTensor(6,4) = metricTensor(4,6)
      metricTensor(6,5) = metricTensor(5,6)
      metricTensor(6,6) = pMoI(3)
    End If

  End Subroutine rigidcoords_getMetricTensor

  !-----------------------------------------------------------------------
  ! This routine does "periodic boundary conditions" on the generalized
  ! coordinates.  E.g. center-of-mass is translated to the simulation cell
  ! and all angles are between -pi and pi
  ! Requires: gcoords -- single rigid molecule generalized coords
  !           simcell -- simulation cell information
  !-----------------------------------------------------------------------
  Subroutine rigidcoords_pbc(gcoords, simcell)
    Type(RigidMolec), Pointer        :: gcoords
    Type(Simcell_Params), Intent(in) :: simcell
    
    Type(VecType)   :: pbccom
    
    Call simcell_pbc(simcell, gcoords%com, pbccom)
    gcoords%com = pbccom
    gcoords%phi = gcoords%phi - Anint(gcoords%phi/twopi)*twopi
    gcoords%psi = gcoords%psi - Anint(gcoords%psi/twopi)*twopi

    !** Change theta to a value between -pi and pi.  Ideally we want
    !** theta to be between 0 and pi but that probably requires us to 
    !** change phi and psi. For now this will have to Do
    gcoords%theta = gcoords%theta - Anint(gcoords%theta/twopi)*twopi

  End Subroutine rigidcoords_pbc

  !-----------------------------------------------------------------------
  ! Read the generalized coordinates from the restart file
  ! Requires: gcoords -- single rigid molecule generalized coords
  !           unitno -- unit to read from
  !-----------------------------------------------------------------------
  Subroutine rigidcoords_readrestartfile(gcoords, unitno)
    Type(RigidMolec), Pointer      :: gcoords
    Integer, Intent(in)            :: unitno

    Integer                        :: a, error, nfields
    Real(kind=Rdbl)                :: x, y, z, check
    Character(len=lstrLen)          :: flexi,junk
    Character(len=strLen),Dimension(strLen) :: fields

    !** Read in the flexibility flag
    Read(unitno,'(a)') flexi
    flexi = stripcmnt(flexi)
    nfields = split(flexi,fields)

    !** Check to see if we are reading an old restart file
    !** An older restart file will just have x,y,z,theta,phi,psi
    !** and no other information.
    check = toreal(fields(1),error)
    If ((error == 0).And.(nfields == 6)) Then
      !** We are reading an old file. Use the old format
      x = toreal(fields(1))
      y = toreal(fields(2))
      z = toreal(fields(3))
      gcoords%theta = toreal(fields(4))
      gcoords%phi = toreal(fields(5))
      gcoords%psi = toreal(fields(6))

      gcoords%com = (/ x, y, z /)
      
      Return
    End If

    !** If we made it this far, it must not be an old file. Read in 
    !** using the new format
    If (ToUpper(Trim(fields(1))) == 'INTERNALLY_FLEXIBLE') Then
      gcoords%internally_flexible = .True.
    Else
      gcoords%internally_flexible = .False.
    End If

    Read(unitno,*) x,y,z,gcoords%theta,gcoords%phi,gcoords%psi,gcoords%dof
    gcoords%com = (/x, y, z/)

    !** consider pulling reference coordinates from molecules here

    If (.Not. gcoords%internally_flexible) Return 

    Read(unitno,*) junk
    If (ToUpper(Trim(junk)) /= 'REFERENCE_COORDINATES') Then
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Expected reference coordinates ID string"
      Stop    
    End If

    If (.Not. Associated(gcoords%ref_struc)) Then
      Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Expected the reference coordinates to be associated"
      Stop    
    End If

    !** read the references coordinates
    Do a = 1,Size(gcoords%ref_struc)
      Read(unitno,*) x, y, z
      gcoords%ref_struc(a) = (/x, y, z/)
    End Do

  End Subroutine rigidcoords_readrestartfile

  !-----------------------------------------------------------------------
  ! Write out the generalized coordinates to the restart file.  This will
  ! Dump the reference coordinates only if the molecule is internally 
  ! flexible.
  ! Requires: gcoords -- single rigid molecule generalized coords
  !           unitno -- unit to dump into
  !-----------------------------------------------------------------------
  Subroutine rigidcoords_dumprestartinfo(gcoords, unitno)
    Type(RigidMolec), Pointer      :: gcoords
    Integer, Intent(In)            :: unitno

    Integer                        :: a
    Character(len=lstrLen)         :: string

    If (gcoords%internally_flexible) Then
      Write(unitno,'(a)') 'INTERNALLY_FLEXIBLE'
!      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
!          " WARNING: untested restartfile dumps with internal flexibility"
    Else
      Write(unitno,'(a)') 'FULLY_RIGID'
    End If

    string = vector_display(gcoords%com, "f13.5")
    Write(unitno,'(a, 3f13.5,i5)') Trim(string), &
        gcoords%theta, gcoords%phi, gcoords%psi, gcoords%dof

    If (.Not. gcoords%internally_flexible) Return

    Write(unitno,'(a)') 'REFERENCE_COORDINATES'

    !** dump the references coordinates
    Do a = 1,Size(gcoords%ref_struc)
      string = vector_display(gcoords%ref_struc(a),'f13.5')
      Write(unitno,'(a)') Trim(string)
    End Do

  End Subroutine rigidcoords_dumprestartinfo


  !-----------------------------------------------------------------------
  ! This routine returns the values of the generalized coordinates
  ! as a string
  ! Requires: gcoords -- single rigid molecule generalized coords
  !           strformat -- display parameters for COM, example: 'f12.5'
  !-----------------------------------------------------------------------
!!$  Character(len=lstrLen) Function rigidcoords_displaystr(gcoords, strformat)
  Function rigidcoords_displaystr(gcoords, strformat)
    Character(len=lstrLen)       :: rigidcoords_displaystr
    Type(RigidMolec), Pointer    :: gcoords
    Character(*), Intent(in)     :: strformat

    Character(len=strLen)        :: fmt

    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
        " WARNING: _displaystr only dumps COM and Eulerian angles"

    Write(fmt, '(3a)') "(a,3",Trim(strformat),")"
    Write(rigidcoords_displaystr, fmt) Trim(vector_display(gcoords%com, &
        strformat)),gcoords%theta, gcoords%phi, gcoords%psi

  End Function rigidcoords_displaystr

  !-------------------------------------------------------------------------
  ! Converts a string to generalized coordinates.  This will ensure that
  ! the order of the generalized coordinates is maintained when we dump
  ! strings and read strings using the routines rigidgc_displaystr, 
  ! rigidgc_convstr
  !-------------------------------------------------------------------------
  Subroutine rigidcoords_convstr(gcoords, str)
    Type(RigidMolec), Pointer   :: gcoords
    Character(*), Intent(in)    :: str

    Integer               :: nfields
    Character(len=strLen), Dimension(10) :: fields
    Real(kind=RDbl)       :: x, y, z

    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
        " WARNING: _convstr only gets COM and Eulerian angles"
    
    nfields = split(str, fields, " ")
    If (nfields > 6) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " This does not look like rigid generalized coordinates"
      Stop
    End If

    x = toreal(fields(1))
    y = toreal(fields(2))
    z = toreal(fields(3))
    gcoords%com = (/x, y, z/)
    gcoords%theta = toreal(fields(4))
    gcoords%phi   = toreal(fields(5))
    gcoords%psi   = toreal(fields(6))
  End Subroutine rigidcoords_convstr

  !-----------------------------------------------------------------------
  ! Clean a single molecule's generalized rigid coordinates
  ! Requires: gcoords -- single rigid molecule generalized coords
  !-----------------------------------------------------------------------
  Subroutine rigidcoords_clean(gcoords)
    Type(RigidMolec), Pointer      :: gcoords

    Integer              :: error

    Deallocate(gcoords%ref_struc, STAT=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__,'ref_struc')
    
  End Subroutine rigidcoords_clean

  !-----------------------------------------------------------------------
  ! Sets the generalized coordinates given a set of xyz coordinates
  ! Requires:  gcoords -- generalized rigid coordinates for one molecule
  !            spc -- species number
  !            rp -- set of xyz coordinate vectors for molecule
  !            com -- molecular center of mass
  !-----------------------------------------------------------------------
  Subroutine rigidcoords_setfrom_rp(gcoords, spc, rp, com)
    Type(RigidMolec), Pointer                    :: gcoords
    Integer, Intent(In)                          :: spc
    Type(VecType), Dimension(:), Intent(In)      :: rp 
    Type(VecType), Intent(In)                    :: com

    Integer              :: atm

    !** Set COM
    gcoords%com = com

    !** Set the angles and reference structure if appropriate
    If (.Not. gcoords%internally_flexible) Then
      Call rigidcoords_fromxyz(gcoords, rp, spc)

    Else
      Do atm = 1,Size(gcoords%ref_struc)
        gcoords%ref_struc(atm) = rp(atm) - com
      End Do

      !** Zero the Euler angles
      gcoords%theta = zero
      gcoords%phi = zero
      gcoords%psi = zero
    End If

  End Subroutine rigidcoords_setfrom_rp


End Module rigidcoords
