!------------------------------------------------------------------
! This module contains the parameters for each plane making up the 
! minimum energy path
!-----------------------------------------------------------------
Module mepplanes
  Use defaults
  Use utils
  Use file
  Use nvtmc
  Use stats
  Use molecules
  Use config
  Implicit None
  Save

  !** The basic idea is that a number of paths cross the planes.  At
  !** each plane there is a set of "initialcoords" which are
  !** the final coordinates of the previous plane translated by an
  !** amount giving by the inter-plane spacing.  These initial 
  !** coords are then annealed giving the "finalcoords" in the plane
  !** for each path.
  Type PATH_Params
    Type(GeneralizedCoords), Dimension(:), Pointer  :: initialcoords
    Type(GeneralizedCoords), Dimension(:), Pointer  :: finalcoords
    Real(kind=RDbl)    :: totnrg
    Integer            :: nsorbs
    Type(Statistics), Dimension(MAX_SORBS) :: avgx,avgy,avgz,theta,phi,psi
  End Type PATH_Params

  !** This structure contains the parameter for each plane.  
  !** Each plane typically will have a number of paths going
  !** through it.  These paths will either be distinct or
  !** coincident
  Type PLANE_Params
    Integer            :: csorbtype  ! The index of the constrained molecule
                                     ! as stored in the sorbates array
    Integer            :: csorbindex ! The index of the constrained molecule
                                     ! as stored in the mc params
    Type(VecType)      :: initpoint
    Type(VecType)      :: normal
    Type(PATH_Params), Dimension(:), Pointer :: path
    Integer            :: npaths
  End Type PLANE_Params

Contains
  !----------------------------------------------------------------
  ! Initialize the parameters for an MEP simulation
  !----------------------------------------------------------------
  Subroutine mepplanes_init(planes, sorbates, mcparams, unitno)
    Type(PLANE_Params), Dimension(:), Intent(out) :: planes
    Type(AtMolCoords), Dimension(:), Intent(inout)  :: sorbates
    Type(NVTMC_Params), Intent(in) :: mcparams
    Integer, Intent(in) :: unitno

    Character(len=strLen)  :: normalline, spacingline, initptline
    Character(len=strLen)  :: sorbname, csorbname
    Integer                :: error, i, j, sorbtype, nplanes, npaths
    Integer                :: csorbtype, csorbindex
    Type(VecType), Dimension(Size(planes,1)) :: normals
    Real(kind=RDbl)  :: ptx, pty, ptz

    nplanes = Size(planes, 1)

    !** Get the no. of iterations in each plane, constraint sorbate name
    !** and the normals to each plane
    Read(unitno,*) csorbname
    Read(unitno,*) npaths
    Read(unitno,'(a)') normalline
    Call mepplanes_parsenormals(normals, normalline)

    csorbtype = molecules_gettype(csorbname)
    Do i=1, nplanes
      planes(i)%csorbtype = csorbtype
      planes(i)%csorbindex = nvtmc_getsorbindex(mcparams, csorbtype)
      planes(i)%npaths = npaths
      planes(i)%normal = normals(i)
      Call mepplanes_initpaths(planes(i), mcparams)
    End Do

    !** Set the initial coordinates for all the paths in the first plane
    Call mepplanes_setfirstplanepts(planes(1), sorbates, mcparams, initptline)

  End Subroutine mepplanes_init

  !--------------------------------------------------------------------
  ! Allocate and initialize the generalized coordinates that store the 
  ! coordinates of the minimum for each path in the plane "plane" 
  !--------------------------------------------------------------------
  Subroutine mepplanes_initpaths(plane, mcparams)
    Type(PLANE_Params), Intent(inout) :: plane
    Type(NVTMC_Params), Intent(in)    :: mcparams

    Integer   :: i, j, npaths, nsorbs, error, sorbtype

    !** Initialize all the paths
    npaths = plane%npaths
    Allocate(plane%path(npaths), Stat=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Could not allocate memory for 'paths'"
      Stop
    End If

    nsorbs = nvtmc_getnsorbs(mcparams)
    Do i=1, plane%npaths
      plane%path(i)%nsorbs = nsorbs
      plane%path(i)%totnrg = 0.0_RDbl

      Allocate(plane%path(i)%initialcoords(nsorbs), Stat=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not allocate memory for 'initialcoords'"
        Stop
      End If

      Allocate(plane%path(i)%finalcoords(nsorbs), Stat=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not allocate memory for 'finalcoords'"
        Stop
      End If
      Do j=1, nsorbs
        sorbtype = nvtmc_getsorbtype(mcparams, j)
        Call gcmodels_initcoords(plane%path(i)%initialcoords(j), sorbtype)
        Call gcmodels_initcoords(plane%path(i)%finalcoords(j), sorbtype)
        Call stats_init(plane%path(i)%avgx(j), "x position :", 10, &
            .False., "f8.3")
        Call stats_init(plane%path(i)%avgy(j), "y position :", 10, &
            .False., "f8.3")
        Call stats_init(plane%path(i)%avgz(j), "z position :", 10, &
            .False., "f8.3")
        Call stats_init(plane%path(i)%theta(j),"theta (rad):", 10, &
            .False., "f8.3")
        Call stats_init(plane%path(i)%phi(j),  "phi   (rad):", 10, &
            .False., "f8.3")
        Call stats_init(plane%path(i)%psi(j),  "psi   (rad):", 10, &
            .False., "f8.3")
      End Do
    End Do
  End Subroutine mepplanes_initpaths

  !----------------------------------------------------------------------
  ! There are two options for getting the points for the first plane.
  ! i) Get the coordinates of one triplet from the "restart" file and
  !    duplicate the coordinates for all the paths.
  !ii) Get the coordinates of all paths from a "meprestart" file.
  !----------------------------------------------------------------------
  Subroutine mepplanes_setfirstplanepts(plane, sorbates, mcparams, initptline)
    Type(PLANE_Params), Intent(inout) :: plane
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(NVTMC_Params), Intent(in)    :: mcparams
    Character(*), Intent(in)          :: initptline

    Character(len=strLen)  :: sourcetype, sourcefilename
    Integer                :: i, j, nfields, sorbtype, nsorbs, npaths
    Real(kind=RDbl), Dimension(3)  :: initcoords
    Type(VecType)          :: initcom

    !** Get some sundry information
    nsorbs = nvtmc_getnsorbs(mcparams)
    npaths  = plane%npaths

    !** For each source type find the source of the initial configuration
    !** and accordingly set the paths.
    Do i=1, nsorbs
      sorbtype   = nvtmc_getsorbtype(mcparams, i)
      sourcetype = config_getsourcetype(sorbates, sorbtype)
      
      If (sourcetype == "MEPRESTARTFILE") Then
        sourcefilename = config_getsourcefile(sorbates, sorbtype)
        Call mepplanes_readmeprestartfile(plane, sorbates, mcparams, &
            sorbtype, sourcefilename)

      Else If (sourcetype == "RESTARTFILE") Then
        Do j=1, npaths
          plane%path(j)%initialcoords(i) = sorbates(sorbtype)%gcoords(1)
        End Do
        
      End If
    End Do
  End Subroutine mepplanes_setfirstplanepts

  !----------------------------------------------------------------
  ! This routine reads the initial generalized coordinates from the
  ! given by "filename".  Also initialize the sorbates structure.
  !----------------------------------------------------------------
  Subroutine mepplanes_readmeprestartfile(plane, sorbates, mcparams, &
      sorbno, filename)
    Type(PLANE_Params), Intent(inout) :: plane
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(NVTMC_Params), Intent(in) :: mcparams
    Integer, Intent(in)            :: sorbno
    Character(*), Intent(in)       :: filename

    Character(len=strLen), Dimension(10) :: fields
    Character(len=strLen)          :: srchstr, molecname
    Character(len=lstrLen)         :: line
    Integer     :: unitno, error, natoms, lineno, nfields
    Integer     :: i, mep_nmoles, mep_natoms, sorbindx
    Logical     :: found

    !** Open the file if not already open
    unitno = isfileopen(filename)
    If (unitno < 0) Then
      unitno = file_getunit(filename)
      Open(unit = unitno, file=filename, status='old', IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            "  Could not open file ", trim(filename)
        Stop
      End If
    End If

    !** Get some sundry information
    molecname = molecules_name(sorbno)
    natoms    = molecules_getnatoms(sorbno)

    !** Find the string "_Molecule_Name_:"
    srchstr = "_MOLECULE_NAME_:"
    found = .False.
    Do 
      lineno = filesrchstr(unitno, srchstr, line, .False.)
      nfields = split(line, fields, " ")
      If (fields(2) == molecname) Then
        found = .True.
        Exit
      End If
    End Do
    
    If (.Not. (found)) Then
      Write(0,'(1x,2a,i4, 4a)') __FILE__," : ",__LINE__, &
          " Could not find the molecule ", Trim(molecname), &
          " in restart file ", filename
      Stop
    End If

    !** Read the natoms and nmoles
    Read(unitno,*) mep_nmoles
    Read(unitno,*) mep_natoms
    If (mep_natoms /= natoms) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " The number of atoms do not match "
      Stop
    End If

    !** Read in the coordinates from the mep restart file
    Read(unitno,*)       ! Comment line
    Do i=1, plane%npaths
      !** Read the generalized coordinates
      sorbindx = nvtmc_getsorbindex(mcparams, sorbno)
      Call gcmodels_readrestartfile(plane%path(i)%initialcoords(sorbindx), &
          unitno)
    End Do
  End Subroutine mepplanes_readmeprestartfile

  !----------------------------------------------------------------
  ! This routine writes the generalized coordinates to the file
  ! given by "filename".
  !----------------------------------------------------------------
  Subroutine mepplanes_writemeprestartfile(plane,sorbates,pathno,optfilename)
    Type(PLANE_Params), Intent(inout) :: plane
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Integer, Intent(in)    :: pathno
    Character(*), Optional, Intent(in) :: optfilename

    Character(len=strLen)    :: filename
    Character(len=strLen)    :: molecname
    Integer     :: unit, unitno, error, natoms, nmoles, nsorbs
    Integer     :: i, m, sorbindx

    !** Get the filename if it is not passed a parameter
    If (Present(optfilename)) Then
      filename = optfilename
      unit = file_getunit(filename)
    Else
      Call file_gettype(d_res_file, filename, unit)
    End If
    
    !** Open the file if not already open
    unitno = isfileopen(filename)
    If (unitno < 0) Then
      unitno = unit
      Open(unit = unitno, file=filename, IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            "  Could not open file ", trim(filename)
        Stop
      End If
    End If

    Write(unitno,*) Trim(general_getfiledesc())
    Write(unitno,*) general_getcurrentiter()
    Write(unitno,*) general_getcurrenttime()

    nsorbs = plane%path(pathno)%nsorbs
    Do i=1, nsorbs
      molecname = molecules_name(i)
      nmoles = sorbates(i)%nmoles
      natoms = sorbates(i)%natoms
      Write(unitno,'(2a)') "_MOLECULE_NAME_: ", molecname
      Write(unitno,'(i20, a)') nmoles, " # Nmoles"
      Write(unitno,'(i20, a)') natoms, " # Natoms"

      Write(unitno,*) "_Generalized Coordinates_"
      Do m=1, pathno
        Call gcmodels_dumprestartinfo(plane%path(m)%finalcoords(i), unitno)
      End Do
    End Do
    
    !** All done here
    Close(unitno)
    Return
  End Subroutine mepplanes_writemeprestartfile

  !----------------------------------------------------------------
  ! Parses the normal line to set the normals to the planes
  ! If the value is specified in the control file then all the
  ! planes have the same normal otherwise it reads it from a file.
  !---------------------------------------------------------------
  Subroutine mepplanes_parsenormals(normals, normalline)
    Type(VecType), Dimension(:), Intent(inout) :: normals
    Character(*), Intent(in)           :: normalline

    Character(len=strLen)  :: stripline
    Character(len=strLen), Dimension(10)   :: strfields
    Integer                :: nfields, npts, i, nplanes
    Real(kind=RDbl)        :: n1, n2, n3

    nplanes = Size(normals, 1)
    stripline = stripcmnt(normalline)
    
    !** Check the no. of fields.  If the no. of fields is more than one then
    !** all the planes have the same normal, otherwise read
    !** them from a file.  The fields are separated by a comma
    nfields = split(stripline, strfields, ",")
    
    If (nfields == 1) Then
      ! Read from a file
      Call mepplanes_readnormals(normals, strfields(1))
    Else
      ! Generate the normals
      Do i=1, nplanes
        n1 = toreal(strfields(1))
        n2 = toreal(strfields(2))
        n3 = toreal(strfields(3))
        normals(i) = (/n1, n2, n3/)
      End Do
    End If
  End Subroutine mepplanes_parsenormals

  !-----------------------------------------------------------
  ! This routine reads in the values of the normals from the
  ! file "filename"
  !-----------------------------------------------------------
  Subroutine mepplanes_readnormals(normals, filename)
    Type(VecType), Dimension(:), Intent(inout) :: normals
    Character(*), Intent(in)           :: filename

    Integer               :: unitno, error, npts, i, nplanes
    Real(kind=RDbl)       :: n1, n2, n3
    
    nplanes = Size(normals, 1)

    !** Open the file
    unitno = isfileopen(filename)
    If (unitno < 0) Then
      unitno = file_getunit(filename)
      Open(unit = unitno, file=filename, status='old', IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            "  Could not open file ", trim(filename)
        Stop
      End If
    End If

    !** Read the no. of points
    Read(unitno, *) npts

    !** Check that the no. of points is the same as specified in the
    !** control file
    If (npts /= nplanes) Then
      Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
          " The no. of points in the file ", Trim(filename), &
          " does not match that in the control file"
      Stop
    End If
    Do i=1, npts
      Read(unitno, *) n1, n2, n3
      normals(i) = (/n1, n2, n3/)
    End Do
    Close(unitno)
    Return
  End Subroutine mepplanes_readnormals
  
  !-----------------------------------------------------------------
  ! Nudges the final coords back into the plane if they have moved
  ! out due to the machine round-off
  !-----------------------------------------------------------------
  Subroutine mepplanes_setfinalcoords(sorbates, mcparams, plane, pathno)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(NVTMC_Params), Intent(in)   :: mcparams
    Type(PLANE_Params), Intent(inout):: plane
    Integer, Intent(in)              :: pathno

    Integer          :: i, nsorbs, sorbtype, csorbtype
    Type(VecType)    :: currentpt, normal, displace, initpoint
    Type(VecType), Dimension(1) :: xyzcoords
    Real(kind=RDbl)  :: dist
    
    !** Get some sundry information
    csorbtype = plane%csorbtype
    currentpt = config_getrp(sorbates, (/csorbtype, 1, 1/))
    normal    = plane%normal
    nsorbs    = plane%path(pathno)%nsorbs
    initpoint = gcmodels_getcom(plane%path(pathno)%initialcoords(csorbtype))

    !** Adjust the final position of the constrained molecule so that
    !** it lies in the plane
    If (.Not. isinplane(currentpt, initpoint, normal)) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " The constrained molecule is not in the correct plane"
      Stop
    Else
      ! We are most probably a little bit off.  Lets nudge it.
      dist = getplanedist(currentpt, normal, initpoint, normal)
      displace = normal*dist
      Call gcmodels_translate(sorbates(csorbtype)%gcoords(1), &
          xyzcoords(1:1), csorbtype, displace)
    End If
    
    !** Copy the generalized coordinates into finalcoords
    Do i=1, nsorbs
      sorbtype = nvtmc_getsorbtype(mcparams, i)
      plane%path(pathno)%finalcoords(i) = sorbates(sorbtype)%gcoords(1)
    End Do
  End Subroutine mepplanes_setfinalcoords


  !-------------------------------------------------------------
  ! Copies the final generalized coordinates stored in "plane1"
  ! to the initial coordinates of "plane2"
  !-------------------------------------------------------------
  Subroutine mepplanes_copyfinalgcoords(plane1, plane2)
    Type(PLANE_Params), Intent(in) :: plane1
    Type(PLANE_Params), Intent(out) :: plane2

    Integer    :: npaths, i, j

    npaths = Min(plane1%npaths, plane2%npaths)
    Do i=1, npaths
      Do j=1, plane2%path(i)%nsorbs
        plane2%path(i)%initialcoords(j) = plane1%path(i)%finalcoords(j)
      End Do
    End Do
  End Subroutine mepplanes_copyfinalgcoords

  !----------------------------------------------------------------------
  ! This routine the displaces the generalized coordinates of the
  ! constrained molecule by amount "spacing" along the vector "dispvect"
  !----------------------------------------------------------------------
  Subroutine mepplanes_displace(plane, spacing, dispvect)
    Type(PLANE_Params), Intent(inout) :: plane
    Real(kind=RDbl), Intent(in)       :: spacing
    Type(VecType), Intent(in)         :: dispvect

    Integer        :: npaths, i, csorbtype, csorbindex
    Type(VecType)  :: displacement
    Type(VecType), Dimension(1) :: xyzcoords

    npaths = plane%npaths
    csorbtype  = plane%csorbtype
    csorbindex = plane%csorbindex
    displacement = dispvect * spacing
    Do i=1, npaths
      Call gcmodels_translate(plane%path(i)%initialcoords(csorbindex), &
          xyzcoords(1:1), csorbtype, displacement)
    End Do
  End Subroutine mepplanes_displace
  
  !--------------------------------------------------------------
  ! Get the no. of paths
  !--------------------------------------------------------------
  Integer Function mepplanes_getnopaths(plane)
    Type(PLANE_Params), Intent(in)    :: plane
    mepplanes_getnopaths = plane%npaths
    Return
  End Function mepplanes_getnopaths

  !--------------------------------------------------------------
  ! Get the normal to a plane
  !--------------------------------------------------------------
  Type(VecType) Function mepplanes_getnormal(plane)
    Type(PLANE_Params), Intent(in)    :: plane
    mepplanes_getnormal = plane%normal
    Return
  End Function mepplanes_getnormal


!!$  !----------------------------------------------------------
!!$  ! Display the mep parameters for the simulation at plane
!!$  ! "planeno"
!!$  !----------------------------------------------------------
!!$  Subroutine mepplanes_displaysimparams(mepparams, planeno, nspc, optunit)
!!$    Type(MEPPLANES_Params), Intent(in) :: mepparams
!!$    Integer, Intent(in)          :: planeno
!!$    Integer, Intent(in)          :: nspc
!!$    Integer, Optional, Intent(in) :: optunit
!!$
!!$    Character(len=nspc) :: spc
!!$    Integer    :: i, unitno, sorbtype
!!$    
!!$    spc = ''
!!$    Do i=1, nspc
!!$      spc = spc//' '
!!$    End Do
!!$
!!$    If (Present(optunit)) Then
!!$      unitno = optunit
!!$    Else
!!$      unitno = 6
!!$    End If
!!$    
!!$    sorbtype = mepparams%constrsorbtype
!!$    Write(unitno, '(2a)') spc, dashedline
!!$    Write(unitno, '(2a)') spc, "The MEP Simulation Section:"
!!$    Write(unitno, '(a,2x,a,i6)') spc,"Plane No.     : ", planeno
!!$    Write(unitno, '(a,2x,2a)')   spc,"Initial Point : ", &
!!$        Trim(display(mepparams%constrcoords, "f8.3"))
!!$    Write(unitno, '(a,2x, 2a)') spc, "Normal        : ", &
!!$        Trim(display(mepparams%normals(planeno), "f8.3"))
!!$    Write(unitno, '(a,2x,a)') spc, dashedline
!!$    Write(unitno, '(a,2x,a)') spc, "The Initial Monte Carlo Parameters:"
!!$    Call nvtmc_displaysimparams(mepparams%initmcparams, nspc+2, unitno)
!!$    Write(unitno, '(a,2x,a)') spc, dashedline
!!$    Write(unitno, '(a,2x,a)') spc, "The Annealing Parameters:"
!!$    Call nvtmc_displaysimparams(mepparams%annealparams, nspc+2, unitno)
!!$    
!!$  End Subroutine mepplanes_displaysimparams


  !---------------------------------------------
  ! Display the mep parameters
  !---------------------------------------------
  Subroutine mepplanes_display(planes, nspc, optunit)
    Type(PLANE_Params), Dimension(:), Intent(in) :: planes
    Integer, Intent(in) :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer    :: i, unitno, nplanes
    
    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    nplanes = Size(planes, 1)

    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "The MEP PLANE Parameters Section:"
    Write(unitno, '(a, 2x, a8, a16, 8x, a12)') spc, &
        "Plane No.", "Normal", "Iters/plane"
    Do i=1, nplanes
      Write(unitno, '(a,2x,i8, a24, i8)') spc, i, &
          Trim(vector_display(planes(i)%normal, "f8.3")), planes(i)%npaths
    End Do
  End Subroutine mepplanes_display

End Module mepplanes

