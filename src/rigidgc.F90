!---------------------------------------------------
! This module contains the data structure defining
! the generalized coordinates of a rigid molecule.
! Namely, the 3 center of mass coordinates and the
! 3 eulerian angles.  It also has coordinates that
! go between xyz and generalized coords.
!---------------------------------------------------
Module rigidgc

  Use defaults, Only: RDbl, MAX_DOF, lstrLen, strLen, zeroTolerance, Pi, twopi
  Use file, Only:
  Use random, Only: rranf
  Use utils, Only: getcom, split, toreal, getinvangle
  Use vector, Only: VecType, vector_crossprod, vector_display, Assignment(=), &
      Operator(+), Operator(-), Operator(*), Operator(/)
  Use matrix, Only: MatrixType, matrix_b2feuler, matrix_Db2fDTheta, &
      matrix_Db2fDPhi, matrix_Db2fDPsi, Operator(*), Assignment(=)
  Use bmap, Only: bmap_getbiasindex, bmap_getncubelets, bmap_getbiaswt, &
      bmap_getbiaspt, bmap_getcellwt
  Use simcell, Only: simcell_getcellorigin, simcell_maptouc, SimCell_Params, &
      simcell_pbc
  Use molecule, Only: LocalCoordInfo
  Use molecules, Only: molecules_getdof, molecules_getbodyaxes, &
      molecules_getmass, molecules_getmoi, molecules_getcomdefn
  Use rigidmcparams, Only: RigidMolecGCMCParams, RigidMolecNVTParams

  Implicit None
  Save

  Private
  Public :: RigidMolec, rigidgc_displaystr, rigidgc_getcom, &
      rigidgc_getgencoords, rigidgc_initcoords, &
      rigidgc_gcoords_eq_vectype, rigidgc_toxyz, rigidgc_fromxyz, &
      rigidgc_setgencoords, rigidgc_getmetrictensor, &
      insert, delete, translate, rotate, rigidgc_pbc, &
      rigidgc_readrestartfile, rigidgc_dumprestartinfo, rigidgc_convstr, &
      rigidgc_displaycoords, rigidgc_getgenforce

  Type RigidMolec
    Type(VecType)     :: com
    Real(kind=RDbl)   :: theta       ! Range 0 to pi
    Real(kind=RDbl)   :: phi, psi    ! Range 0 to 2*pi
    Integer           :: dof         ! # of degrees of freedom
  End Type RigidMolec
  
!!$  Interface init
!!$    Module Procedure rigidgc_initcoords
!!$  End Interface

  Interface insert
    Module Procedure rigidgc_insert
  End Interface

  Interface delete
    Module Procedure rigidgc_delete
  End Interface

  Interface translate
    Module Procedure rigidgc_randomtranslate
    Module Procedure rigidgc_translate
  End Interface

  Interface rotate
    Module Procedure rigidgc_randomrotate
  End Interface

!!$  Interface rigidgc_display
!!$    Module Procedure rigidgc_displaycoords
!!$    Module Procedure rigidgc_displaystr
!!$  End Interface
    
Contains
  !----------------------------------------------
  ! Initializes the model params
  !----------------------------------------------

  Subroutine rigidgc_initcoords(gcoords, sorbtype)

    Type(RigidMolec), Pointer  :: gcoords
    Integer, Intent(in)        :: sorbtype

    Integer  :: error, dof

    dof = molecules_getdof(sorbtype)

    Allocate(gcoords, STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          "Could not allocate 'gcoords'"
      Stop
    Endif

    gcoords%com = 0.0_RDbl
    gcoords%theta = 0.0_RDbl
    gcoords%phi = 0.0_RDbl
    gcoords%psi = 0.0_RDbl
    gcoords%dof = dof

  End Subroutine rigidgc_initcoords

  !-------------------------------------------------------------
  ! Set the center-of-mass of the generalized coordinates
  !-------------------------------------------------------------
  Subroutine rigidgc_gcoords_eq_vectype(gcoords, vec)
    Type(RigidMolec), Pointer       :: gcoords
    Type(VecType), Intent(in)       :: vec
    
    gcoords%com = vec
  End Subroutine rigidgc_gcoords_eq_vectype

  !--------------------------------------------------------------
  ! Generate the xyz coordinates "xyzcoords" from the generalized
  ! coordinates "gcoords"
  !--------------------------------------------------------------
  Subroutine rigidgc_toxyz(gcoords, xyzcoords, com_defn)
    Type(RigidMolec), Pointer               :: gcoords
    Type(VecType), Dimension(:), Intent(out):: xyzcoords
    Type(VecType), Dimension(:), Intent(in) :: com_defn

    Type(MatrixType)   :: b2f_eulermat
    Integer            :: natoms, i
    
    !** Check the no. of atoms and return if there is only one
    natoms = Size(com_defn,1)
    If (natoms == 1) Then
      xyzcoords(1) = gcoords%com
      Return
    End If

    !** Get the transformation matrix
    b2f_eulermat = matrix_b2feuler(gcoords%phi, gcoords%theta, gcoords%psi)
    
    !** Do the transformation for each atom
    Do i=1, natoms
      xyzcoords(i) = gcoords%com + b2f_eulermat * com_defn(i)
    End Do
  End Subroutine rigidgc_toxyz

  !------------------------------------------------------------
  ! This routine takes an array of xyz coordinates (xyzcoords)
  ! and generates generalized coordinates from them
  !------------------------------------------------------------
  Subroutine rigidgc_fromxyz(gcoords, xyzcoords, sorbtype)
    Type(RigidMolec), Pointer               :: gcoords
    Type(VecType), Dimension(:), Intent(in) :: xyzcoords
    Integer, Intent(in)     :: sorbtype
    
    Type(VecType) :: rcom, xf, yf, zf, e1, e2, e3
    Type(VecType), Dimension(3) :: atvec
    Type(LocalCoordInfo) :: lcoord
    Integer         :: i, atom1, atom2
    Real(kind=RDbl) :: theta, costheta, sintheta, psi, sinpsi, cospsi
    Real(kind=RDbl) :: phi, sinphi, cosphi
    
    !** Get the center of mass from the xyz coordinates
    rcom = getcom(xyzcoords)
    
    !** Generate in the fixed coordinate system what the position of 
    !** the body axes would be.  The coordinates of the body axes
    !** as expressed in the fixed coordinates is giveny by "xf, yf, zf"
    lcoord = molecules_getbodyaxes(sorbtype)
    atom1  = lcoord%atomnos(1)
    atom2  = lcoord%atomnos(2)
    atvec(1) = xyzcoords(atom1) - rcom
    atvec(2) = xyzcoords(atom2) - rcom
    atvec(3) = vector_crossprod(atvec(1),atvec(2))
    xf = 0.0_RDbl
    yf = 0.0_RDbl
    zf = 0.0_RDbl
    Do i=1, 3
      xf  = xf + atvec(i)*lcoord%xcoeff(i)
      yf  = yf + atvec(i)*lcoord%ycoeff(i)
      zf  = zf + atvec(i)*lcoord%zcoeff(i)
    End Do

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
  End Subroutine rigidgc_fromxyz

  !------------------------------------------------------------
  ! Gets the center-of-mass of the molecule described by "gcoords"
  !------------------------------------------------------------
  Type(VecType) Function rigidgc_getcom(gcoords)
    Type(RigidMolec), Pointer   :: gcoords

    rigidgc_getcom = gcoords%com
  End Function rigidgc_getcom

  !-----------------------------------------------------------
  ! Return the generalized coordinates as an array
  !-----------------------------------------------------------
  Function rigidgc_getgencoords(gcoords)
    Type(RigidMolec), Pointer           :: gcoords
    Real(kind=RDbl), Dimension(MAX_DOF) :: rigidgc_getgencoords
    
    rigidgc_getgencoords(1:3) = gcoords%com
    rigidgc_getgencoords(4:6) = (/gcoords%theta, gcoords%phi, gcoords%psi/)
    Return
  End Function rigidgc_getgencoords

  !------------------------------------------------------------------
  ! Set the generalized coordinates from the values in "genarray"
  ! "genarray" can have 3 elements in which case only the COM is set.
  ! If the number of values >= 6 then the eulerian angles are set.
  !------------------------------------------------------------------
  Subroutine rigidgc_setgencoords(gcoords, genarray)
    Type(RigidMolec), Pointer                  :: gcoords
    Real(kind=RDbl), Dimension(:), Intent(in)  :: genarray

    If (Size(genarray) >= 3) Then
      gcoords%com   = genarray(1:3)
      If (Size(genarray) >= 6 ) Then
        gcoords%theta = genarray(4)
        gcoords%phi   = genarray(5)
        gcoords%psi   = genarray(6)
      End If
    Else
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " 'genarray' has too few values"
      Stop
    End If

  End Subroutine rigidgc_setgencoords

  !---------------------------------------------------------------
  ! Returns the forces in generalized coordinates "genforce"
  ! given the forces in xyz coordinates on each atom in "xyzforce"
  !---------------------------------------------------------------
  Subroutine rigidgc_getGenForce(gcoords, xyzforce, comdefn, genforce)
    Type(RigidMolec), Pointer               :: gcoords
    Type(VecType), Dimension(:), Intent(in) :: xyzforce
    Type(VecType), Dimension(:), Intent(in) :: comdefn
    Real(kind=RDbl), Dimension(:), Intent(out) :: genforce

    Integer         :: natoms, i, comcoord, dof
    Real(kind=RDbl) :: phi, theta, psi
    Real(kind=RDbl), Dimension(3)  :: f
    Type(MatrixType)   :: Db2fDTheta, Db2fDPhi, Db2fDPsi
    
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
      genforce(4) = genforce(4) + xyzforce(i)*(Db2fDTheta*comdefn(i))
    End Do

    !** Get the derivatives wrt to the Eulerian angles
    Db2fDPhi = matrix_Db2fDPhi(phi, theta, psi)
    Do i=1, natoms
      genforce(5) = genforce(5) + xyzforce(i)*(Db2fDPhi*comdefn(i))
    End Do

    !** Get the derivatives wrt to the Eulerian angles
    Db2fDPsi = matrix_Db2fDPsi(phi, theta, psi)
    Do i=1, natoms
      genforce(6) = genforce(6) + xyzforce(i)*(Db2fDPsi*comdefn(i))
    End Do
    
  End Subroutine rigidgc_getGenForce

  !---------------------------------------------------------------------------
  ! Returns the metric tensor for the sorbate "sorbtype".  The metric tensor
  ! is defined by Randy in J. Phys. Chem., Vol. 98, 1994, pp. 11948--11961
  ! "ierr" is non-zero when metric tensor is singular which happens when 
  ! theta is zero.
  !---------------------------------------------------------------------------
  Subroutine rigidgc_getMetricTensor(gcoords, sorbtype, metricTensor, ierr)
    Type(RigidMolec), Pointer                    :: gcoords
    Integer, Intent(in) :: sorbtype
    Real(kind=RDbl), Dimension(:,:), Intent(out) :: metricTensor
    Integer, Intent(out)  :: ierr

    Integer              :: i, dof
    Real(kind=RDbl)      :: mass
    Real(kind=RDbl), Dimension(3) :: pMoI
    Real(kind=Rdbl)      :: sintheta, costheta, sinphi, cosphi, sinpsi, cospsi

    dof = molecules_getdof(sorbtype)
    
    !** For a spherical molecule it is just the mass multiplied
    !** by a unit matrix
    metricTensor = 0.0_Rdbl
    If (dof == 3) Then
      mass = molecules_getmass(sorbtype)
      Do i=1, dof
        metricTensor(i,i) = mass
      End Do
    Else
      mass = molecules_getmass(sorbtype)
      pMoi = molecules_getmoi(sorbtype)
      Do i=1, 3
        metricTensor(i,i) = mass
      End Do
      sintheta = Sin(gcoords%theta)
      costheta = Cos(gcoords%theta)
      sinphi   = Sin(gcoords%phi)
      cosphi   = Cos(gcoords%phi)
      sinpsi   = Sin(gcoords%psi)
      cospsi   = Cos(gcoords%psi)
      ! Check if theta is close to zero
      If (Abs(sintheta) < 1.0e-4_RDbl) Then
        ierr = 1
        Return
      End If

      ! theta
      metricTensor(4,4) = pMoI(1)*cospsi*cospsi + pMoI(2)*sinpsi*sinpsi
      metricTensor(4,5) = pMoI(1)*sintheta*sinpsi*cospsi &
                        - pMoI(2)*sintheta*sinpsi*cospsi
      metricTensor(4,6) = 0.0_RDbl
      ! phi              
      metricTensor(5,4) = metricTensor(4,5)
      metricTensor(5,5) = pMoI(1)*sintheta*sintheta*sinpsi*sinpsi &
                        + pMoI(2)*sintheta*sintheta*cospsi*cospsi &
                        + pMoI(3)*costheta*costheta
      metricTensor(5,6) = pMoI(3)*costheta
      ! psi              
      metricTensor(6,4) = metricTensor(4,6)
      metricTensor(6,5) = metricTensor(5,6)
      metricTensor(6,6) = pMoI(3)
    Endif
  End Subroutine rigidgc_getMetricTensor
  
  !-----------------------------------------------------------
  ! Generates a center-of-mass of the molecule according to a
  ! bias map and the using random eulerian angles generates the 
  ! xyz coordinates.  It also returns the bias of the insertion
  ! in "biasfactor"
  !-----------------------------------------------------------
  Subroutine rigidgc_insert(gcoords, xyzcoords, gcparams, biasfactor)
    Type(RigidMolec), Pointer                   :: gcoords
    Type(VecType), Dimension(:), Intent(inout)  :: xyzcoords
    Type(RigidMolecGCMCParams), Intent(in) :: gcparams
    Real(kind=RDbl), Intent(out)   :: biasfactor

    Integer              :: natoms, i, index, nx, ny, nz
    Integer              :: icellx,icelly,icellz, ncubelets
    Real(kind=RDbl)      :: costheta
    Type(VecType)        :: com

    !** Get the position of the center of mass using the bias map
    !** 1.) Pick a point in the unit cell according to the biasmap
    index = bmap_getbiasindex(gcparams%bmap)
    ncubelets  = bmap_getncubelets(gcparams%bmap)
    biasfactor = 1.0_RDBl/bmap_getbiaswt(gcparams%bmap,index)/ncubelets
    com   = bmap_getbiaspt(gcparams%bmap, index)

    !** 2.) Select one of the unit cells in the simulation cell
    icellx = Int(rranf()*gcparams%simcell%nx)
    icelly = Int(rranf()*gcparams%simcell%ny)
    icellz = Int(rranf()*gcparams%simcell%nz)

    !** 3.) Translate the pt. to the appropriate unit cell
    com = com + simcell_getcellorigin(gcparams%simcell, icellx, icelly, icellz)
    gcoords%com = com
    gcoords%theta = 0.0_RDbl
    gcoords%phi   = 0.0_RDbl
    gcoords%psi   = 0.0_RDbl

    !** If no. of atoms is 1 we are done
    If (Size(xyzcoords, 1) == 1) Then
      xyzcoords(1) = com
      Return
    End If

    !** Generate the eulerian angles at random
    ! Sample theta from a cosine distribution.
    costheta = 1 - rranf() * 2.0
    gcoords%theta = Acos(costheta)
    gcoords%phi = pi - rranf()*twopi
    gcoords%psi = pi - rranf()*twopi

    !** Generate the xyz coordinates
    Call rigidgc_toxyz(gcoords, xyzcoords, gcparams%com_defn)
  End Subroutine rigidgc_insert

  !-----------------------------------------------------------------
  ! This routine gets the bias-weight for deleting a molecule
  !-----------------------------------------------------------------
  Subroutine rigidgc_delete(gcoords, gcparams, biasfactor)
    Type(RigidMolec), Pointer              :: gcoords
    Type(RigidMolecGCMCParams), Intent(in) :: gcparams
    Real(kind=RDbl), Intent(out)   :: biasfactor    

    Integer         :: index, ncubelets
    Real(kind=RDbl) :: wt
    Real(kind=RDbl), Dimension(3)  :: uccom
    
    !** Map the position to the unit cell. Also convert the 
    !** vector type to an array
    uccom = simcell_maptouc(gcparams%simcell, gcoords%com)

    !** Get the weigth of cubelet where the COM lies
    wt = bmap_getcellwt(gcparams%bmap, uccom)

    !** Calculate the biasfactor
    ncubelets  = bmap_getncubelets(gcparams%bmap)
    biasfactor = ncubelets*wt

    Return
  End Subroutine rigidgc_delete

  !-----------------------------------------------------------------
  ! This routine shifts the center of mass of a molecule
  ! with generalized coordinates "gcoords" by amount "deltacom" and 
  ! generates the "xyzcoords"
  !----------------------------------------------------------------
  Subroutine rigidgc_translate(gcoords, xyzcoords, sorbtype, deltacom)
    Type(RigidMolec), Pointer                   :: gcoords
    Type(VecType), Dimension(:), Intent(inout)  :: xyzcoords
    Integer, Intent(in)       :: sorbtype
    Type(VecType), Intent(in) :: deltacom

    !** Displace the molecule
    gcoords%com = gcoords%com + deltacom

    !** Generate the new coordinates
    Call rigidgc_toxyz(gcoords, xyzcoords, molecules_getcomdefn(sorbtype))
  End Subroutine rigidgc_translate

  !-------------------------------------------------------------
  ! Translates the center-of-mass of the molecule and generates
  ! xyzcoords.
  !-------------------------------------------------------------
  Subroutine rigidgc_randomtranslate(gcoords, xyzcoords, gcparams)
    Type(RigidMolec), Pointer                   :: gcoords
    Type(VecType), Dimension(:), Intent(inout)  :: xyzcoords
    Type(RigidMolecNVTParams), Intent(in) :: gcparams

    Type(VecType)       :: deltadisp, com
    Real(kind=RDbl)     :: deltatrans, xdelta, ydelta, zdelta
    
    deltatrans = gcparams%deltatrans
    If (gcparams%constrained) Then
      !** Do a coordinate transformation to the constrained plane
      !** system where the normal is along the x-axis
      com = gcoords%com
      com = gcparams%constrtrans*com

      !** Get the displacement in the y-z plane
      ydelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      zdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      deltadisp = (/0.0_RDbl, ydelta, zdelta/)
      com = com + deltadisp

      !** Tranform back to the original coordinate system
      gcoords%com = gcparams%constrinvtrans*com
    Else
      xdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      ydelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      zdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      deltadisp = (/xdelta, ydelta, zdelta/)
      gcoords%com = gcoords%com + deltadisp
    End If
    
    !** Generate the new coordinates
    Call rigidgc_toxyz(gcoords, xyzcoords, gcparams%com_defn)
  End Subroutine rigidgc_randomtranslate

  !-------------------------------------------------------------
  ! Rotate the molecule by generating new eulerian angles and 
  ! generates xyzcoords.
  !-------------------------------------------------------------
  Subroutine rigidgc_randomrotate(gcoords, xyzcoords, gcparams)
    Type(RigidMolec), Pointer                   :: gcoords
    Type(VecType), Dimension(:), Intent(inout)  :: xyzcoords
    Type(RigidMolecNVTParams), Intent(in) :: gcparams

    Integer             :: natoms
    Type(VecType)       :: deltadisp, com
    Real(kind=RDbl)     :: deltarot, oldphi, phi, oldpsi, psi
    Real(kind=RDbl)     :: oldcostheta, costheta

    !** Check the no. of atoms and return if no. of atoms is 1
    natoms = Size(xyzcoords, 1)
    If (natoms == 1) Return

    deltarot = gcparams%deltarot

    !** Sample phi and psi uniformly between -pi,pi but theta
    !** needs to be sampled from a cosine distribution
    !** See Allen and Tildesley, Pg. 132-133, "Molecular Liquids"
    phi = gcoords%phi + ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
    gcoords%phi = phi - Anint(phi/twopi)*twopi

    psi = gcoords%psi + ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
    gcoords%psi = psi - Anint(psi/twopi)*twopi

    costheta = Cos(gcoords%theta) + ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
    costheta = costheta - Anint(costheta/2.0_RDbl)*2.0_RDbl
    gcoords%theta = Acos(costheta)
    
    !** Generate the new coordinates
    Call rigidgc_toxyz(gcoords, xyzcoords, gcparams%com_defn)
  End Subroutine rigidgc_randomrotate

  !-----------------------------------------------------------------------
  ! This routine does "periodic boundary conditions" on the generalized
  ! coordinates.  E.g. center-of-mass is translated to the simulation cell
  ! and all angles are between -pi and pi
  !-----------------------------------------------------------------------
  Subroutine rigidgc_pbc(gcoords, simcell)
    Type(RigidMolec), Pointer        :: gcoords
    Type(Simcell_Params), Intent(in) :: simcell
    
    Type(VecType)   :: pbccom
    
    Call simcell_pbc(simcell, gcoords%com, pbccom)
    gcoords%com = pbccom
    gcoords%phi = gcoords%phi - Anint(gcoords%phi/twopi)*twopi
    gcoords%psi = gcoords%psi - Anint(gcoords%psi/twopi)*twopi
    !Change theta to a value between -pi and pi.  Ideally we want
    !theta to be between 0 and pi but that probably requires us to 
    !change phi and psi. For now this will have to Do
    gcoords%theta = gcoords%theta - Anint(gcoords%theta/twopi)*twopi
  End Subroutine rigidgc_pbc

  !----------------------------------------------------------------
  ! Read the generalized coordinates from the restart file
  !----------------------------------------------------------------
  Subroutine rigidgc_readrestartfile(gcoords, unitno)
    Type(RigidMolec), Pointer      :: gcoords
    Integer, Intent(in)            :: unitno

    Real(kind=Rdbl)                :: x, y, z

    Read(unitno,*) x, y, z, gcoords%theta, gcoords%phi, gcoords%psi
    gcoords%com = (/x, y, z/)
  End Subroutine rigidgc_readrestartfile

  !-----------------------------------------------------------
  ! Write out the generalized coordinates to the restart file
  !-----------------------------------------------------------
  Subroutine rigidgc_dumprestartinfo(gcoords, unitno)
    Type(RigidMolec), Pointer      :: gcoords
    Integer, Intent(in)            :: unitno

    Write(unitno,'(a, 3f13.5)') Trim(vector_display(gcoords%com, "f13.5")), &
        gcoords%theta, gcoords%phi, gcoords%psi
  End Subroutine rigidgc_dumprestartinfo


  !----------------------------------------------------------------
  ! This routine returns the values of the generalized coordinates
  ! as a string
  !----------------------------------------------------------------
!!$  Character(len=lstrLen) Function rigidgc_displaystr(gcoords, strformat)
  Function rigidgc_displaystr(gcoords, strformat)
    Character(len=lstrLen) :: rigidgc_displaystr
    Type(RigidMolec), Pointer    :: gcoords
    Character(*), Intent(in)     :: strformat

    Character(len=strLen)        :: fmt

    Write(fmt, '(3a)') "(a,3",Trim(strformat),")"
   Write(rigidgc_displaystr, fmt) Trim(vector_display(gcoords%com, &
        strformat)),gcoords%theta, gcoords%phi, gcoords%psi
  End Function rigidgc_displaystr

  !-------------------------------------------------------------------------
  ! Converts a string to generalized coordinates.  This will ensure that
  ! the order of the generalized coordinates is maintained when we dump
  ! strings and read strings using the routines rigidgc_displaystr, 
  ! rigidgc_convstr
  !-------------------------------------------------------------------------
  Subroutine rigidgc_convstr(gcoords, str)
    Type(RigidMolec), Pointer   :: gcoords
    Character(*), Intent(in)    :: str

    Integer               :: nfields, i
    Character(len=strLen), Dimension(10) :: fields
    Real(kind=RDbl)       :: x, y, z
    
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
  End Subroutine rigidgc_convstr

  !------------------------------------------------------
  ! Writes out the various parameters and the coordinate
  ! information to the unit "unit"
  !------------------------------------------------------
  Subroutine rigidgc_displaycoords(gcoords, xyzcoords, unit)
    Type(RigidMolec), Pointer               :: gcoords
    Type(VecType), Dimension(:), Intent(in) :: xyzcoords
    Integer, Optional, Intent(in)   :: unit
    Integer   :: unitno, i, natoms
    
    If (Present(unit)) Then
      unitno = unit
    Else
      unitno = 6
    End If
    Write(unitno,'(1x,2a)') 'Center of Mass: ', &
         vector_display(gcoords%com, "f8.3")
    Write(unitno,'(1x,a,3f8.3)') 'Eulerian Angles: ', gcoords%phi, &
        gcoords%theta, gcoords%psi
    natoms = Size(xyzcoords,1)
    ! Dump out the xyzcoords in a RASMOL format with atomtype "C"
    Write(unitno,'(1x, i6)') natoms
    Write(unitno,*)
    Do i = 1,natoms
      Write(unitno,'(1x, 2a)') "C", vector_display(xyzcoords(i), "f8.3")
    End Do

  End Subroutine rigidgc_displaycoords

End Module rigidgc
