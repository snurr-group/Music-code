!------------------------------------------------------------------
! This is the driver module for getting the Free Energies By
! Putting Planes that are tangent to a given Minimum Energy Path
!-----------------------------------------------------------------
Module fnrg
  Use defaults
  Use utils
  Use vector
  Use file
  Use stats
  Use moves
  Use random
  Use forcefield
  Use stats
  Implicit None
  Save

!!$  Private
!!$  Public :: FNRG_SORB_Params, fnrg_init 

  Integer, Parameter   :: MAX_PATH_PTS = 500
  Integer, Parameter   :: MAX_DIMS = 30
  Integer, Parameter   :: MEP_HEADER_LINES = 19
  Integer, Parameter   :: MAX_MOVES = 2

  !** The tag marking the beginning of the Minimum Energy Path section
  Character(len=strLen), Parameter  :: default_fnrg_tag = &
      "Free Energy Path Info"  

  !** Parameters for individual sorbate types
  Type FNRG_SORB_Params
    Character(len=strLen) :: sorbname
    Integer               :: sorbtype
    Integer               :: dof
    Integer               :: nmoves
    Type(Move_Params), Dimension(MAX_MOVES) :: mvparams
    Type(Move_Stats), Dimension(MAX_MOVES)  :: accratio
  End Type FNRG_SORB_Params

  !** The box parameters defines a box bounded by 6 planes
  !** Each plane is defined by a point and a plane
  Type Box_Params
    Real(kind=RDbl), Dimension(3) :: lb, ub !lower bound, upper bound
    Logical, Dimension(3)  :: pbc ! Whether to do pbc in x, y, z directions
  End Type Box_Params

  !** The overall free energy parameters
  Type FNRG_Params
    Integer                :: nplanes
    Integer                :: startplaneno
    Integer                :: max_planes
    Integer                :: niterations
    Integer                :: dumpfreq
    Logical                :: reversepath
    Integer                :: nsorbs
    Type(FNRG_SORB_Params), Dimension(:), Pointer :: sorbparams
    Real(kind=RDbl), Dimension(MAX_DIMS) :: begconfig
    Integer                :: nopts1, nopts2
    Real(kind=RDbl), Dimension(:,:), Pointer :: origpath
    Integer                :: ncpts1, ncpts2
    Real(kind=RDbl), Dimension(:,:), Pointer :: cmppath
    Real(kind=RDbl), Dimension(:,:), Pointer :: normals
    Real(kind=RDbl)        :: spacing
    Real(kind=RDbl)        :: lastenergy
    Real(kind=RDbl)        :: temp
    Real(kind=RDbl)        :: rti
    Character(len=strLen)  :: pathfile
    Character(len=strLen)  :: basefile
    Integer                :: currentplane
    Integer                :: dof
    Integer                :: totaldims
    Integer                :: cnstrsorbno    ! The sorb no. to be constrained
    Integer                :: cnstrsorbtype  ! The sorbtype to be constrained
    Type(Box_Params)       :: cnstrBox
  End Type FNRG_Params

Contains
  !----------------------------------------------------------------
  ! Initialize the parameters for an MEP simulation
  !----------------------------------------------------------------
  Subroutine fnrg_init(feparams, simcell, sorbates, ctrl_filename, opt_fnrgtag)
    Type(FNRG_Params), Intent(out)   :: feparams
    Type(Simcell_Params), Intent(in) :: simcell
    Type(AtMolCoords), Dimension(:), Intent(inout), Target :: sorbates
    Character(*), Intent(in) :: ctrl_filename
    Character(*), Optional, Intent(in) :: opt_fnrgtag

    Character(len=strLen)    :: tag, line, nvtmc_tag, normalline, spacingline
    Character(len=strLen)    :: sorbname
    Integer          :: unitno, lineno, error, i, j, sorbtype, dof, nmoves
    Integer          :: reverseflag, natoms, nfields, flag
    Real(kind=RDbl)  :: ptx, pty, ptz
    Character(len=strLen)  :: gcmodeltype, junk, pt_norm_line, constrname
    Type(VecType)    :: pt, norm

    If (Present(opt_fnrgtag)) Then
      tag = opt_fnrgtag
    Else
      tag = default_fnrg_tag
    End If

    !** Open the ctrl_file if it is not opened
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
      Open(file=ctrl_filename, unit=unitno)
    Endif
    
    !** Find the Free Energy section
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif

    !** Read the rest of the stuff
    Read(unitno,*) feparams%nplanes
    Read(unitno,*) feparams%startplaneno
    Read(unitno,*) feparams%niterations
    Read(unitno,*) feparams%dumpfreq

    !** Get the spacing between the planes
    Read(unitno,*) feparams%spacing

    !** Read the path file
    feparams%reversepath = .False.
    Read(unitno,*) feparams%pathfile
    Read(unitno,*) reverseflag
    If (reverseflag /= 0) Then
      feparams%reversepath = .True.
    End If

    !** Read the basefilename, MC section tag and no. of sorbates
    Read(unitno,*) feparams%basefile
    Read(unitno,*) feparams%temp
    Read(unitno,*) feparams%nsorbs
    Read(unitno,*)
    feparams%rti = 1.0_Rdbl/(feparams%temp*kjmole_kb) ! 1/(kj/mole K)
    
    !** Initialize the different sorbate parameters
    Allocate(feparams%sorbparams(feparams%nsorbs), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,  &
          " Could not allocate memory for 'feparams%sorbparams'"
      Stop
    End If
    dof = 0
    natoms = 0
    Do i=1, feparams%nsorbs
      Read(unitno,*) feparams%sorbparams(i)%sorbname
      sorbtype = molecules_gettype(feparams%sorbparams(i)%sorbname)
      feparams%sorbparams(i)%dof = molecules_getdof(sorbtype)
      dof = dof + molecules_getdof(sorbtype)
      natoms = natoms + molecules_getnatoms(sorbtype)
      feparams%sorbparams(i)%sorbtype = sorbtype
      
      ! Initialize the different move types
      Read(unitno,*) feparams%sorbparams(i)%nmoves
      Do j=1, feparams%sorbparams(i)%nmoves
        Call moves_initparams(feparams%sorbparams(i)%mvparams(j), &
            simcell, sorbtype, ctrl_filename)
      End Do

      ! Read the inter-sorbate spacing
      If (i /= feparams%nsorbs) Then
        Read(unitno, *)
      End If
    End Do
    feparams%dof = dof
    feparams%totaldims = 3*natoms

    !** Read the constrain parameters
    Read(unitno,*) 
    Read(unitno,*) constrname
    Do i=1, 3
      Read(unitno,*) flag, feparams%cnstrBox%lb(i), feparams%cnstrBox%ub(i)
      If (flag == 1) Then
        feparams%cnstrBox%pbc(i) = .True.
      Else
        feparams%cnstrBox%pbc(i) = .False.
      End If
    End Do
    Do i=1, feparams%nsorbs
      If (feparams%sorbparams(i)%sorbname == constrname) Then
        feparams%cnstrsorbno = i
        feparams%cnstrsorbtype = feparams%sorbparams(i)%sorbtype
        Exit
      End If
    End Do
    If (i > feparams%nsorbs) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find sorbate ", Trim(constrname), " in control file"
      Stop
    End If
    
    !** This routine sets the path.  It does two things, 1) sets up the
    !** initial configuration of the system and 2) generates the constraint
    !** normals along the way
    Call fnrg_setpath(feparams, sorbates, simcell)

    !** Set the default planeno to a negative number indicating it has
    !** not been initialized
    feparams%currentplane = -1
  End Subroutine fnrg_init

  !-------------------------------------------------------------------
  ! 1) Sets the initial configuration
  ! 2) Also generates the constraint normals from tangents to
  ! points along the path read from the control file
  !-------------------------------------------------------------------
  Subroutine fnrg_setpath(feparams, sorbates, simcell)
    Type(FNRG_Params), Intent(inout) :: feparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(in) :: simcell 

    Integer    :: i, j, pathunit, npts1, npts2, dof, nfields, lineno
    Integer    :: totalpts, natoms, ncpts1, ncpts2, max_planes, nsorbs, error
    Integer    :: sorbtype, max_pts
    Character(len=strLen) :: line, tag1, junk1, junk2
    Character(len=strLen), Dimension(10) :: fields
    Real(kind=RDbl), Dimension(MAX_PATH_PTS, feparams%dof) :: path
    Real(kind=RDbl), Dimension(MAX_PATH_PTS, feparams%dof) :: cmppath
    Real(kind=RDbl), Dimension(feparams%dof) :: temp
    Real(kind=RDbl)     :: norm

    !** Get some sundry information
    dof = feparams%dof
    nsorbs = feparams%nsorbs

    !** Open the file with the path
    pathunit = isfileopen(feparams%pathfile)
    If (pathunit < 0) Then
      pathunit = file_getunit(feparams%pathfile)
      Open(unit=pathunit, file=feparams%pathfile)
    End If
    Rewind(pathunit)

    !** Read in the path from the "pathfile"
    ! Find the no. of points from the minima to the saddle-point
    tag1 = "_NPTS1"
    lineno = filesrchstr(pathunit, tag1, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Did not find field ", Trim(tag1), " in MEP header"
      Stop
    End If
    nfields = split(line, fields)
    npts1 = toint(fields(3))
    Read(pathunit, *) junk1, junk2, npts2
    ! Now read the path coordinates from MIN1 -> TS
    Read(pathunit, *) junk1, (path(1, j), j=1, dof)
    Do i=2, npts1
      Read(pathunit, *) (path(i, j), j=1, dof)
    End Do

    ! Read the transition state coordinates
    Read(pathunit, *) junk1, (path(npts1+1, j), j=1, dof)

    ! Now read the path coordinates from TS -> MIN2
    Do i=1, npts2-1
      Read(pathunit, *) (path(npts1+1+i, j), j=1, dof)
    End Do
    Read(pathunit, *) junk1, (path(npts1+1+i, j), j=1, dof)
    totalpts = npts1 + 1 + npts2

    !** Check the flag from the control file to see if we need
    !** to reverse the direction of the path
    If (feparams%reversepath) Then
      Do i=1, totalpts/2
        temp = path(i, :)
        path(i,:) = path(totalpts-i+1, :)
        path(totalpts-i+1, :) = temp
      End Do
    End If

    !** Save the original path
    Allocate(feparams%origpath(totalpts, dof), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Could not allocate 'feparams%path'"
      Stop
    End If
    Do i=1, totalpts
      feparams%origpath(i,1:dof) = path(i, 1:dof)
    End Do
    feparams%nopts1 = npts1
    feparams%nopts2 = npts2

    !** Break up the path into the spacing specified in the control
    !** file.  For now only look at the displacement of the centers-of-mass
    !** of the different sorbates
    ! Compress the points from MIN1 -> TS
    ncpts1 = fnrg_compressPath(feparams, path(1:npts1, :), cmppath)
    ! Copy the transition state
    cmppath(ncpts1+1, :) = path(npts1+1,:) 
    ! Compress the points from TS -> MIN2
    ncpts2 = fnrg_compressPath(feparams, path(npts1+2:totalpts, :), &
        cmppath(ncpts1+2:MAX_PATH_PTS, :))
    max_pts = ncpts1 + ncpts2 + 1
    max_planes = max_pts - 1 ! One less than the total no. of points
    feparams%max_planes = max_planes

    !** Save the compressed path
    Allocate(feparams%cmppath(max_pts, dof), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Could not allocate 'feparams%cmppath'"
      Stop
    End If
    Do i=1, max_pts
      feparams%cmppath(i,1:dof) = cmppath(i, 1:dof)
    End Do
    feparams%ncpts1 = ncpts1
    feparams%ncpts2 = ncpts2

    !** Make sure that the no. of planes specified is less than the
    !** total no. of points possible
    If (feparams%nplanes > (max_planes)) Then
      Write(0,'(1x,2a,i4,a,i4)') __FILE__," : ",__LINE__, &
          " Too many planes in control file.  Max. no.: ", max_planes
      Stop
    End If

    !** If the no. of planes is 0 then set it to max_planes
    If (feparams%nplanes == 0) Then
      feparams%nplanes = max_planes
    End If

    !** Initialize the sorbates structure
    Do i=1, feparams%nsorbs
      sorbtype = feparams%sorbparams(i)%sorbtype
      natoms   = molecules_getnatoms(sorbtype)
      Call config_allocfields(sorbates(sorbtype), sorbtype, 1)
      sorbates(sorbtype)%natoms = natoms
      sorbates(sorbtype)%nmoles = 1
      sorbates(sorbtype)%filename = feparams%pathfile
      sorbates(sorbtype)%sourcetype = "Steepest Descent Path"
    End Do
    Call fnrg_xToSorb(feparams, sorbates, simcell, path(1,:))
    
    !** Generate the normals
    Allocate(feparams%normals(feparams%nplanes, feparams%dof), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate memrory for 'dispvect'"
      Stop
    End If
    Do i=1, feparams%nplanes
      feparams%normals(i,:) = cmppath(i+1,:) - cmppath(i,:) 
      norm = arrnorm(feparams%normals(i,:))
      feparams%normals(i,:) = feparams%normals(i,:)/norm
!!$      Write(*,'(a,15f7.2)') "i  ", (cmppath(i+1,j), j=1, 15)
!!$      Write(*,'(a,15f7.2)') "i-1", (cmppath(i,j), j=1, 15)
!!$      Write(*,'(a,15f7.2)') "n",   (feparams%normals(i,j), j=1, 15)
    End Do
  End Subroutine fnrg_setpath


  !----------------------------------------------------------------------
  ! This function takes a path and the spacing between planes and gives
  ! the compressed path in "cpath".  It returns the no. of points in the
  ! compressed path
  !----------------------------------------------------------------------
  Integer Function fnrg_compressPath(feparams, path, cpath)
    Type(FNRG_Params), Intent(in) :: feparams
    Real(kind=RDbl), Dimension(:,:), Intent(in)  :: path
    Real(kind=RDbl), Dimension(:,:), Intent(out) :: cpath

    Integer     :: i, j, npts, ncomppts, nsorbs
    Real(kind=RDbl) :: spacing, max_diff
    Type(VecType), Dimension(feparams%nsorbs) :: sorb_coms1, sorb_coms2
    
    !** Initialize some variables and store the first point of the path
    npts = Size(path, 1)
    ncomppts = 1
    spacing = feparams%spacing
    nsorbs  = feparams%nsorbs
    cpath(ncomppts, :) = path(1, :)

    Call fnrg_getcoms(feparams, path(1,:), sorb_coms1)
    Do i=2, npts
      Call fnrg_getcoms(feparams, path(i,:), sorb_coms2)

      ! Find the maximum displacement
      max_diff = -1.0_RDbl     ! Some negative number
      Do j=1, nsorbs
        If (mag(sorb_coms1(j)-sorb_coms2(j)) > max_diff) Then
          max_diff = mag(sorb_coms1(j)-sorb_coms2(j))
        End If
      End Do

      ! Copy the configuration which is the required distance away or
      ! if it is the last point of the path
      If (max_diff > spacing .Or. i==npts) Then
        ncomppts = ncomppts + 1
        cpath(ncomppts, :) = path(i, :)
        ! Update the current centers-of-mass
        sorb_coms1 = sorb_coms2
      End If
    End Do
    
    fnrg_compressPath = ncomppts
    Return
  End Function fnrg_compressPath

  !------------------------------------------------------------------------
  ! Get the centers-of-mass of the sorbates as stored in array "x"
  !------------------------------------------------------------------------
  Subroutine fnrg_getcoms(feparams, x, coms)
    Type(FNRG_Params), Intent(in) :: feparams
    Real(kind=RDbl), Dimension(:), Intent(in)   :: x
    Type(VecType), Dimension(:), Intent(out) :: coms

    Integer :: begindex, i, dof, sorbtype

    begindex = 1
    Do i=1, feparams%nsorbs
      sorbtype = feparams%sorbparams(i)%sorbtype
      dof = molecules_getdof(sorbtype)
      coms(i) = x(begindex:begindex+2)
      begindex = begindex + dof
    End Do
  End Subroutine fnrg_getcoms

  !------------------------------------------------------------------------
  ! Change from the "sorbates" representation to the "x"
  ! representation. It takes the generalized coordinates 
  ! from the "sorbates" structure and puts them in the "x" array.
  !------------------------------------------------------------------------
  Subroutine fnrg_sorbToX(feparams, sorbates, x)
    Type(FNRG_Params), Intent(in) :: feparams
    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Real(kind=RDbl), Dimension(:), Intent(out)  :: x

    Integer  :: nvar, sorbno, dof, begidx, endidx, nsorbs, sorbtype
    Real(kind=RDbl), Dimension(feparams%dof) :: gencoords

    nsorbs = feparams%nsorbs

    begidx = 1
    Do sorbno=1, nsorbs
      sorbtype = feparams%sorbparams(sorbno)%sorbtype
      dof = molecules_getdof(sorbtype)
      endidx = begidx + dof - 1
      gencoords = gcmodels_getgencoords(sorbates(sorbtype)%gcoords(1))
      x(begidx:endidx) = gencoords(1:dof)
      begidx = endidx + 1
    Enddo
  End Subroutine fnrg_sorbToX

  !-----------------------------------------------------------------
  ! Change from the array representation "x" to the "sorbates"
  ! representation
  !-----------------------------------------------------------------
  Subroutine fnrg_xToSorb(feparams, sorbates, simcell, x)
    Type(FNRG_Params), Intent(in) :: feparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(in)               :: simcell
    Real(kind=RDbl), Dimension(:), Intent(in)  :: x

    Integer  :: nvar, sorbno, dof, begidx, endidx, natoms, sorbtype, nsorbs
    Real(kind=RDbl), Dimension(MAX_DOF) :: gencoords

    nsorbs = feparams%nsorbs

    begidx = 1
    Do sorbno=1, nsorbs
      sorbtype = feparams%sorbparams(sorbno)%sorbtype
      natoms = molecules_getnatoms(sorbtype)
      dof = molecules_getdof(sorbtype)
      endidx = begidx + dof - 1

      Call gcmodels_setgencoords(sorbates(sorbtype)%gcoords(1), &
          x(begidx:endidx)) 

      begidx = endidx + 1

      !** Generate the xyz principal coordinates
      Call gcmodels_toxyz(sorbates(sorbtype)%gcoords(1), &
          sorbates(sorbtype)%coords(1:natoms,1)%rp, &
          molecules_getcomdefn(sorbtype))
      
      !** Generate the other coordinates
      Call simcell_pbc(simcell, sorbates(sorbtype)%coords(1:natoms,1)%rp, &
          sorbates(sorbtype)%coords(1:natoms, 1)%r, &
          sorbates(sorbtype)%coords(1:natoms, 1)%cr)    
    Enddo
  End Subroutine fnrg_xToSorb

  !--------------------------------------------------------------
  ! Get the no. of simulations
  !--------------------------------------------------------------
  Integer Function fnrg_getendingplane(feparams)
    Type(FNRG_Params), Intent(in)    :: feparams
    fnrg_getendingplane = feparams%nplanes
    Return
  End Function fnrg_getendingplane

  !--------------------------------------------------------------
  ! Get the no. of simulations
  !--------------------------------------------------------------
  Integer Function fnrg_getstartingplane(feparams)
    Type(FNRG_Params), Intent(in)    :: feparams
    fnrg_getstartingplane = feparams%startplaneno
    Return
  End Function fnrg_getstartingplane

  !-------------------------------------------------------------------------
  ! This routine moves to the next plane and initializes the transformation
  !  matrix based on the normal
  !-------------------------------------------------------------------------
  Subroutine fnrg_movetoplane(feparams, sorbates, simcell, planeno)
    Type(FNRG_Params), Intent(inout) :: feparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(inout) :: simcell
    Integer, Intent(in)           :: planeno

    Integer        :: dof, i, j, nsorbs
    Real(kind=RDbl), Dimension(feparams%dof) :: normal, prod
    Real(kind=RDbl), Dimension(feparams%dof, feparams%dof) :: vecs1, vecs2
    Real(kind=RDbl), Dimension(feparams%dof, feparams%dof) :: tempvec
    Real(kind=RDbl):: pot
    Logical        :: fast, mapflag

    !** Initialize some stuff
    dof    = feparams%dof
    nsorbs = feparams%nsorbs
    normal = feparams%normals(planeno, 1:dof)
    vecs1  = 0.0_RDbl
    vecs2  = 0.0_RDbl
    Do i=1, dof
      vecs1(i,i) = 1.0_RDbl
    End Do
    
    !** Set the current plane to be the next plane
    feparams%currentplane = planeno

    !** Initialize the no. of attempts and accepted move counters
    Do i=1, nsorbs
      feparams%sorbparams(i)%accratio(1:MAX_MOVES)%succ = 0
      feparams%sorbparams(i)%accratio(1:MAX_MOVES)%att = 0
    End Do
    
    !** Move the starting position to the minimum in the new plane
    feparams%begconfig(1:dof) = feparams%cmppath(planeno,1:dof)

    !** Calculate the initial energy
    Call fnrg_xToSorb(feparams, sorbates, simcell, feparams%begconfig)
    fast = .False.
    pot  = 0.0_RDbl
    Call forcefield_getssint(sorbates, simcell, fast, pot, mapflag)
    If (.Not. mapflag) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Initial energy calculation in untabulated region"
      Stop
    End If
    feparams%lastenergy = pot*caltoj
    
  End Subroutine fnrg_movetoplane

  !-----------------------------------------------------------------------
  !  Does the periodic boundary conditions in the direction bounded by
  !  "lobound" and "upbound".  The new point is returned in "pbcpt"
  !-----------------------------------------------------------------------
  Subroutine fnrg_pbc(simcell, sorbates, sorbtype, box)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Integer, Intent(in) :: sorbtype
    Type(Simcell_Params), Intent(in):: simcell
    Type(Box_Params), Intent(in) :: box

    Real(kind=RDbl) :: lb, ub, edgelen, temp
    Real(kind=RDbl), Dimension(3) :: arrpt
    Integer         :: i, indx
    Logical         :: pbc

    !** Use the principal cell coordinates for our purpose
    arrpt = sorbates(sorbtype)%coords(1,1)%rp
    pbc   = .False.

    Do i=1, 3
      !** Do we want to do pbc on this dimension
      If (.Not. box%pbc(i)) Cycle
      
      !** Get the lower bound, upper bound and the edge length
      lb = box%lb(i)
      ub = box%ub(i)
      If (lb > ub) Then
        ! Swap the lower bound with the upper bound
        temp = lb
        lb = ub
        ub = temp
      End If
      edgelen = ub - lb

      !** Translate this coordinate system to one where the box has its
      !** origin at (0, 0, 0).
!!$      Write(*,'(a,3f8.3)') 'arrpt(i) ', arrpt(i), lb, arrpt(i) - lb
      arrpt(i) = arrpt(i) - lb

      !** Now do the periodic boundary condition as with 'simcell'
      indx = Floor(arrpt(i)/edgelen)
      arrpt(i) = arrpt(i) - indx*edgelen
      If (indx /= 0) Then
!!$        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
!!$            " Doing pbc"
!!$        Write(*,*) 'indx ', indx

        !** Translate it back to the original coordinate system
        arrpt(i) = arrpt(i) + lb

!!$        Write(*,'(a,3f8.3)') 'pbc pt.   ', arrpt(1:3)
        pbc = .True.
      End If
    End Do

    !** Update the coordinates if we need be
    If (pbc) Then
      sorbates(sorbtype)%coords(1,1)%rp = arrpt
      
      !** Now generate the principal coords
      !** We assume that the continuation vectors do not change
      Call simcell_pbc(simcell, sorbates(sorbtype)%coords(1,1)%rp, &
          sorbates(sorbtype)%coords(1,1)%r, sorbates(sorbtype)%coords(1,1)%cr)
      Call gcmodels_setgencoords(sorbates(sorbtype)%gcoords(1), arrpt)
    End If
  End Subroutine fnrg_pbc

  !----------------------------------------------------------------------
  ! Do the constrained MC move
  !----------------------------------------------------------------------
  Subroutine fnrg_cnstrMCmove(feparams, sorbates, simcell, sorbtype, &
      movetype, planeno, df)
    Type(FNRG_Params), Intent(inout)  :: feparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_params), Intent(inout) :: simcell
    Integer, Intent(in) :: sorbtype, movetype
    Integer, Intent(in) :: planeno
    Real(kind=RDbl), Dimension(:), Intent(inout) :: df

    Integer :: i, nsorbs, sorbno, nmoves, pcoord, dof, natoms
    Integer :: sorbdof, accept, attempts
    Real(kind=RDbl), Dimension(feparams%dof) :: currentpos, dispvect
    Real(kind=RDbl), Dimension(feparams%dof) :: normal
    Real(kind=RDbl), Dimension(feparams%dof,feparams%dof) :: metricT
    Real(kind=RDbl), Dimension(feparams%dof) :: newdf, savedf, dispdf 
    Real(kind=RDbl) :: comp, deltaE, boltz, rti, newpot, oldpot, disp, ratio
    Real(kind=RDbl) :: ojac, newjac, sintheta, utemp
    Logical         :: fast, mapflag

    !** Get some sundry stuff
    dof    = feparams%dof
    nsorbs = feparams%nsorbs
    rti    = feparams%rti
    oldpot = feparams%lastenergy
    savedf = df
    normal = feparams%normals(planeno,:)
    accept   = feparams%sorbparams(sorbtype)%accratio(movetype)%succ
    attempts = feparams%sorbparams(sorbtype)%accratio(movetype)%att

    !** Do the random displacement
    natoms = molecules_getnatoms(sorbtype)
    !** changed this August 15, 2001 LC, from a call to moves_nvt, correct?
    Call moves_translate(sorbates, sorbtype, 1, &
        feparams%sorbparams(sorbtype)%mvparams(movetype))
    Call simcell_pbc(simcell, sorbates(sorbtype)%coords(1:natoms,1)%rp, &
        sorbates(sorbtype)%coords(1:natoms, 1)%r, &
        sorbates(sorbtype)%coords(1:natoms, 1)%cr)
    attempts = attempts + 1
    feparams%sorbparams(sorbtype)%accratio(movetype)%att  =  attempts

    !** Also make sure that the penetrant stays within the box
    !** Use periodic boundary conditions to map it back into the box
    If (sorbtype == feparams%cnstrsorbtype) Then
      Call fnrg_pbc(simcell, sorbates, feparams%cnstrsorbtype, &
          feparams%cnstrBox)
    End If
    
    !** Convert to the array representation of sorbates
    Call fnrg_sorbToX(feparams, sorbates, df)

    !** Project the displacement vector on to the plane
    dispdf = df - savedf
    df     = df - (Dot_product(dispdf, normal)*normal)

    !** Calculate the energy
    Call fnrg_xToSorb(feparams, sorbates, simcell, df)
    fast = .False.
    Call forcefield_getssint(sorbates, simcell, fast, newpot, mapflag)
    newpot = newpot * caltoj

    !** Check the mapflag
    If (.Not. mapflag) Then
      df = savedf
      Call fnrg_xToSorb(feparams, sorbates, simcell, df)
      Return
    End If

    !** Decide to accept or reject this move
    deltaE = newpot - oldpot   ! kj
    ! Prevent overflow and underflow
    utemp = -deltaE*rti
    If (utemp > 30.0_RDbl)  utemp =  30.0_RDbl
    If (utemp < -30.0_RDbl) utemp = -30.0_RDbl
!!$    Write(*,'(a,4e16.7)') '**** deltaE ', deltaE, Exp(utemp), boltz
    boltz  = Exp(utemp)

    If (rranf() > boltz) Then
      ! Don't accept the move
      df = savedf
      Call fnrg_xToSorb(feparams, sorbates, simcell, df)
    Else
      ! Accept the move
      accept = accept + 1
      feparams%lastenergy = newpot
    End If
    
    !** Rescale the displacement to get an acceptance ratio of about 0.5
    ratio = accept/(attempts*1.0_RDbl)
    feparams%sorbparams(sorbtype)%accratio(movetype)%succ =  accept
    Call moves_adjustparams(feparams%sorbparams(sorbtype)%mvparams(movetype),&
        ratio)
  End Subroutine fnrg_cnstrMCmove

  !--------------------------------------------------------------
  ! Does the actual moves.  For each plane it creates a file with
  ! the ratio of the configuration integral and the delta F.  It
  ! also creates a file for each sorbate with the center-of-mass
  ! position.
  !--------------------------------------------------------------
  Subroutine fnrg_dosim(feparams, sorbates, simcell, planeno)
    Type(FNRG_Params), Intent(inout)    :: feparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(inout)  :: simcell
    Integer, Intent(in)  :: planeno
    
    Integer    :: i, j, iter, sorbno, sorbtype, movetype, ierr
    Integer    :: sorbdof, currentplane, nsorbs, nfields, fnrgunit, nmoves
    Integer    :: att, acc, dumpfreq
    Integer, Dimension(feparams%nsorbs) :: posunits
    Character(len=strLen)         :: molecname
    Character(len=lstrLen)        :: avgposstr
    Character(len=strLen), Dimension(10) :: fields
    Real(kind=RDbl)  :: dist, mag, ratio, deltaE, rti, pot, newpot, spacing
    Real(kind=RDbl)  :: utemp
    Real(kind=RDbl), Dimension(feparams%dof) :: df, startingdf, dispdf, normal
    Logical          :: mapflag, fast
    Type(Statistics) :: deltaF, configratio
    Type(GeneralizedCoords), Dimension(feparams%nsorbs) :: pbcgcoords

    !** Check that we have actually moved to the right plane
    If (feparams%currentplane /= planeno) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Please move to the correct plane before doing the minimization"
      Stop
    End If
    
    !** Get some sundry information
    currentplane = feparams%currentplane
    nsorbs       = feparams%nsorbs
    startingdf   = feparams%begconfig
    normal       = feparams%normals(planeno, :)
    rti          = feparams%rti
    spacing      = feparams%spacing
    dumpfreq     = feparams%dumpfreq
    Call fnrg_sorbToX(feparams, sorbates, df)
    Call stats_init(deltaF, "Free Energy Change", 1000,.False.)
    Call stats_init(configratio, "Ratio of Config. Integ.", 1000,.False.)
    Do i=1, nsorbs
      sorbtype = feparams%sorbparams(i)%sorbtype
      Call gcmodels_initcoords(pbcgcoords(i), sorbtype)
    End Do

    !** Open the output files
    Call fnrg_openoutfiles(feparams, posunits, fnrgunit, planeno)

    !** Do the simulation for the current plane
    Do iter = 1, feparams%niterations
      !** Do the constrained monte-carlo move
      ! Pick a sorbtype
      sorbno   = Int(rranf()*nsorbs) + 1
      sorbtype = feparams%sorbparams(sorbno)%sorbtype
      sorbdof  = molecules_getdof(sorbtype)
      ! One MC cycle consists of "sorbdof" moves
      Do i=1, sorbdof
        nmoves = feparams%sorbparams(sorbno)%nmoves
        ! Pick a move type
        movetype = Int(rranf()*nmoves) + 1

        Call fnrg_cnstrMCmove(feparams, sorbates, simcell, sorbtype, &
            movetype, planeno, df)
      End Do
      ! Get the current energy in kJ/mol
      pot = feparams%lastenergy

      !** Check if the point is still in the plane
      dispdf = df - startingdf
      mag    = Dot_product(dispdf, normal)
      If (mag > 1.0e-7) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Displaced point is no longer in the plane"
        Stop
      End If
      
!!$      Write(*,'(a,15f8.3)') 'df ', df
!!$      Do i=1, 3
!!$        Write(*,'(a,i5,2x,a)') 'sorbate ', i, &
!!$            Trim(gcmodels_displaystr(sorbates(i)%gcoords(1), "f8.3"))
!!$      End Do
      
      !** Displace the current configuration along the normal
      dispdf = df + spacing*normal

      !** Calculate the change in energy
      Call fnrg_xToSorb(feparams, sorbates, simcell, dispdf)
      fast = .False.
      Call forcefield_getssint(sorbates, simcell, fast, newpot, mapflag)
      newpot = newpot * caltoj
      If (.Not. mapflag) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
         ' Displaced point in untabulated region '
        Stop
      End If

      !** Gather the statistics
      deltaE = newpot - pot
      utemp = -deltaE*rti
      Call stats_update(configratio, Exp(utemp))

      !** Restore the point back to the original plane
      Call fnrg_XToSorb(feparams, sorbates, simcell, df)

      !** Display the move stats
      If (Mod(iter, dumpfreq) == 0) Then
        Write(*,*) 'Iteration no. :', iter
        !** Write the move stats
        Do i=1, nsorbs
          sorbtype = feparams%sorbparams(i)%sorbtype
          nmoves = feparams%sorbparams(i)%nmoves
          sorbtype = feparams%sorbparams(i)%sorbtype
          nmoves   = feparams%sorbparams(i)%nmoves
          Write(*,*) 'Sorbate Name  : ', Trim(molecules_name(sorbtype))
          Do j=1, nmoves
            acc = feparams%sorbparams(i)%accratio(j)%succ
            att = feparams%sorbparams(i)%accratio(j)%att
            If (att /= 0) Then
              ratio = acc/(att*1.0_RDbl)
            Else
              ratio = 0.0_RDbl
            End If
            Write(*,'(2x,a)') 'Acceptance Ratio(s): '
            Write(*,'(2x,a,i2,a,2x,2a,f5.2)') "Move Type ", j, " :", &
                Trim(makestr(acc, "/", att)), "=", ratio
            Call moves_displayparams(feparams%sorbparams(i)%mvparams(j), &
                sorbtype)
          End Do
        End Do
      End If

      !** Write to the various output files
      If (Mod(iter, dumpfreq) == 0) Then
        ! Write the ratio and delta free energy
        Write(fnrgunit,'(i8,4e12.4,f9.4)') iter, stats_getvalue(configratio), &
            stats_getblock(configratio), stats_getcavg(configratio), &
            stats_getstd(configratio), &
            -1.0/rti*Log(stats_getcavg(configratio))
        ! Write the center-of-mass position
        Do i=1, nsorbs
          sorbtype = feparams%sorbparams(i)%sorbtype
          pbcgcoords(i) = sorbates(sorbtype)%gcoords(1)
          Call gcmodels_pbc(pbcgcoords(i), simcell)
          Write(posunits(i), '(i8, 5x, a)') iter, &
              Trim(gcmodels_displaystr(pbcgcoords(i),"f8.3"))
        End Do
      End If
    End Do     ! End of iters loop

    !** Close the output files
    Close(fnrgunit)
    Do i=1, nsorbs
      Close(posunits(i))
    End Do
  End Subroutine fnrg_dosim

  !----------------------------------------------------------------
  ! Open the files to save the last positions of the different
  ! sorbates and the energies
  !----------------------------------------------------------------
  Subroutine fnrg_openoutfiles(feparams, posunits, fnrgunit, planeno)
    Type(FNRG_Params), Intent(in)    :: feparams
    Integer, Dimension(:), Intent(inout)  :: posunits
    Integer, Intent(inout)   :: fnrgunit
    Integer, Intent(in)      :: planeno

    Integer    :: i, j, sorbtype
    Integer    :: currentplane, nsorbs, simstart
    Character(len=strLen)  :: posfile, molecname, fnrgfile
    Character(len=strLen)  :: basefile

    !** Get some sundry information
    nsorbs    = feparams%nsorbs
    simstart  = general_getsimstart()
    basefile  = feparams%basefile

    !** Open the free energy file
    fnrgfile = Trim(basefile)//"."//"fnrg"
    fnrgfile = genfilename(fnrgfile, planeno)
    fnrgunit = file_getunit(fnrgfile)
    Open(unit=fnrgunit, file=fnrgfile)
    Write(fnrgunit, '(a,a5,6x,a,7x,a,7x,a,7x,a,6x,a)') "#", &
        "iter", "inst.", "block", "cum", "stddev", "deltaF(kJ)"
    
    !** Open the position files
    Do i=1, nsorbs
      sorbtype   = feparams%sorbparams(i)%sorbtype
      molecname  = molecules_name(sorbtype)
      posfile = Trim(basefile)//"."//molecname
      posfile = genfilename(posfile, planeno)
      posunits(i) = file_getunit(posfile)
      Open(unit=posunits(i), file=posfile)
    End Do
  End Subroutine fnrg_openoutfiles

!!$  !---------------------------------------------------------
!!$  ! This routine writes the current energy to "unitno"
!!$  !---------------------------------------------------------
!!$  Subroutine fnrg_dumpenergy(feparams, planeno, pathno, unitno)
!!$    Type(FNRG_Params), Intent(inout) :: feparams
!!$    Integer, Intent(in) :: planeno, pathno
!!$    Integer, Intent(in) :: unitno
!!$
!!$    Integer, Save :: currentplane = 0
!!$    Integer       :: nsorbs, i, j, simstart
!!$    Real(kind=RDbl)    ::   coulnrg, noncoulnrg, totpairnrg, totnrg
!!$    
!!$    nsorbs    = molecules_getnsorbs()
!!$    simstart  = general_getsimstart()
!!$
!!$    !** Check if this the first call for this plane. If so dump
!!$    !** the header
!!$    If (planeno /= currentplane) Then
!!$      currentplane = planeno
!!$      Write(unitno, '(a,9x)', Advance='No') "#"
!!$      Do i=1, nsorbs-1
!!$        Do j=i+1, nsorbs
!!$          Write(unitno,'(5x,i1,a,i1,4x)', Advance='No') i,"-",j
!!$        End Do
!!$      End Do
!!$      Write(unitno,'(a12)') "Total(kJ)"
!!$    End If
!!$
!!$    !** Dump the energies
!!$    totnrg = 0.0_RDbl
!!$    Write(unitno,'(2i5)', Advance='No') planeno+simstart-1, pathno
!!$    Do i=1, nsorbs-1
!!$      Do j=i+1, nsorbs
!!$        noncoulnrg = molecules_getnoncoul(i,j,"inst")
!!$        coulnrg = molecules_getcoul(i,j,"inst")
!!$        totpairnrg = (coulnrg + noncoulnrg)
!!$        totnrg = totnrg + totpairnrg
!!$        feparams%planes(planeno)%path(pathno)%totnrg = totnrg
!!$        Write(unitno,'(f12.4)', Advance='No') totpairnrg
!!$      End Do
!!$    End Do
!!$    Write(unitno,'(f12.4)') totnrg
!!$  End Subroutine fnrg_dumpenergy

  !----------------------------------------------------------
  ! Display the mep parameters for the simulation at plane
  ! "planeno"
  !----------------------------------------------------------
  Subroutine fnrg_displaysimparams(feparams, planeno, nspc, optunit)
    Type(FNRG_Params), Intent(in) :: feparams
    Integer, Intent(in)          :: planeno
    Integer, Intent(in)          :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer    :: i, unitno, sorbtype, dof
    
    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    !**
    dof = feparams%dof

    !**
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a,i5)') spc, &
        "The Free Energy Simulation Section For Plane No. : ", planeno
    Write(unitno, '(a,2x,a)') spc, dashedline

    ! Write the starting point in this plane
    Write(unitno, '(a,2x,a)')   spc, "Starting Point in this Plane :"
    Do i=1, dof
      Write(unitno, '(a,2x,f7.3)', Advance='No') spc, &
          feparams%begconfig(i)
    End Do
    Write(unitno,*)

    ! Write the normal
    Write(unitno, '(a,2x,a)') spc, "Normal :"
    Do i=1, dof
      Write(unitno, '(a,2x,f5.2)', Advance='No') spc, &
          feparams%normals(planeno, i)
    End Do
    Write(unitno,*)
    
    ! Write the starting energy
    Write(unitno, '(a,2x,a,f12.4)') spc, "Starting Energy(kJ) : ", &
        feparams%lastenergy
  End Subroutine fnrg_displaysimparams


  !----------------------------------------------------------------------
  ! Display the overall free energy calculation parameters
  !----------------------------------------------------------------------
  Subroutine fnrg_display(feparams, nspc, optunit)
    Type(FNRG_Params), Intent(in) :: feparams
    Integer, Intent(in) :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Character, Dimension(3) :: dir
    Integer    :: i, j, unitno, totalpts, sorbtype
    
    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "The Free Energy Parameters Section:"
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a,i6)')   spc,"Max. No. of planes      : ", &
        feparams%max_planes
    Write(unitno, '(a,2x,a,i6)')   spc,"No. of planes           : ", &
        feparams%nplanes
    Write(unitno, '(a,2x,a,i6)')   spc,"No. of iterations/plane : ", &
        feparams%niterations
    Write(unitno, '(a,2x,a,f8.3)') spc,"Inter-plane spacing     : ", &
        feparams%spacing
    Write(unitno, '(a,2x,a,i4)') spc,  "Total Degrees of Freedom: ", &
        feparams%dof
    Write(unitno, '(a,2x,a,f8.2)') spc,  "Temperature             : ", &
        feparams%temp
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno,'(a,2x,a)')  spc, "Constrained Box Parameters"
    Write(unitno, '(a,2x,a)') spc, dashedline
    dir = (/'x', 'y', 'z'/)
    Do i=1, 3
      If (feparams%cnstrBox%pbc(i)) Then
        Write(unitno,'(a,4x,a, a, 2f8.3)')  spc, dir(i), " lbound, ubound :", &
            feparams%cnstrBox%lb(i), feparams%cnstrBox%ub(i)
      Else
        Write(unitno,'(a,4x,a,a)')  spc, dir(i), ": no bounds"
      End If
    End Do

    !** The parameters for each sorbate
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "The Sorbate Parameters Information"    
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a,i4)') spc, "No. of sorbates          : ", &
        feparams%nsorbs
    Do i=1, feparams%nsorbs
      sorbtype = feparams%sorbparams(i)%sorbtype
      Write(unitno,'(a,2x,2a)') spc, "Sorbate Name: ", &
          Trim(molecules_name(sorbtype))
      Write(unitno,'(a,4x,a,i3)') spc, "No. of Move Types: ", &
          feparams%sorbparams(i)%nmoves
      Do j=1, feparams%sorbparams(i)%nmoves
        Call moves_displayparams(feparams%sorbparams(i)%mvparams(j), sorbtype)
      End Do
    End Do
    Write(unitno, '(a)') dashedline    

    !** The original and compressed path information
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "The Path Information"    
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a,a)')  spc, "Pathfile                : ", &
        Trim(feparams%pathfile)
    Write(unitno, '(a,2x,a)') spc, "Original Path"
    Write(unitno, '(a,2x,a,i5)') spc, "No. of points (Min1-TS) : ", &
        feparams%nopts1
    Write(unitno, '(a,2x,a,i5)') spc, "No. of points (TS-Min2) : ", &
        feparams%nopts2
    totalpts = feparams%nopts1+1+feparams%nopts2
    Do i=1, totalpts
      If (i==1) Write(unitno, '(a)', Advance='No') 'MIN_1: '
      If (i==feparams%nopts1+1) Write(unitno, '(a)', Advance='No') 'TS: '
      If (i==totalpts) Write(unitno, '(a)', Advance='No') 'MIN_2: '
      Do j=1, feparams%dof
        Write(unitno, '(f6.2)', Advance='No') feparams%origpath(i, j)
      End Do
      Write(unitno,*)
    End Do

    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "The Compressed Path Information"    
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "Compressed Path"
    Write(unitno, '(a,2x,a,i5)') spc, "No. of points (Min1-TS) : ", &
        feparams%ncpts1
    Write(unitno, '(a,2x,a,i5)') spc, "No. of points (TS-Min2) : ", &
        feparams%ncpts2
    totalpts = feparams%ncpts1+1+feparams%ncpts2
    Do i=1, totalpts
      If (i==1) Write(unitno, '(a)', Advance='No') 'MIN_1: '
      If (i==feparams%ncpts1+1) Write(unitno, '(a)', Advance='No') 'TS: '
      If (i==totalpts) Write(unitno, '(a)', Advance='No') 'MIN_2: '

      Do j=1, feparams%dof
        Write(unitno, '(f6.2)', Advance='No') feparams%cmppath(i, j)
      End Do
      Write(unitno,*)
    End Do

    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "The Normals Information"    
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a, 2x, a9, 20x, a)') spc, "Plane No.", "Normal"
    Do i=1, feparams%nplanes
      Write(unitno, '(a, 2x, i5)', Advance='No') spc, i
      Call dispvec(feparams%normals(i,1:feparams%dof), "f7.3")
      Write(unitno, *)
    End Do
  End Subroutine fnrg_display
End Module fnrg



