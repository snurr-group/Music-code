!------------------------------------------------------------------
! This is the driver module for getting the Minimum Energy Path
! using simulated annealing.
!-----------------------------------------------------------------
Module mep
  Use defaults
  Use utils
  Use file
  Use nvtmc
  Use stats
  Use mepplanes
  Implicit None
  Save

  !** The tag marking the beginning of the Minimum Energy Path section
  Character(len=strLen), Parameter  :: default_mep_tag = &
      "Minimum Energy Path Info"  

  !** This structure stores the lastcoords of each sorbate type at
  !** each point(plane) along a path.  Also stores the total energy
  !** of the last configuration
  Type MEP_PATH
    Type(GeneralizedCoords), Dimension(:), Pointer :: gcoords
    Integer            :: nsorbs
    Real(kind=RDbl)    :: totnrg
  End Type MEP_PATH
  
  Type MEP_Params
    Character(len=strLen)    :: initmctag     
    Type(NVTMC_Params)       :: initmcparams
    Character(len=strLen)    :: annealtag
    Type(NVTMC_Params)       :: annealparams
    Integer                  :: nplanes
    Type(PLANE_Params), Dimension(:), Pointer :: planes
    Real(kind=RDbl), Dimension(:), Pointer    :: spacing
    Type(VecType), Dimension(:), Pointer      :: dispvect
    Character(len=strLen)    :: spacingfile
    Character(len=strLen)    :: basefile
    Integer                  :: currentplane
    Integer                  :: csorbtype
    Integer                  :: csorbindex
    Type(VecType), Pointer   :: constrcoords
    ! The averages are stored in the same index order as they are stored
    ! in the annealparams data structure
    Integer                  :: nsorbs
  End Type MEP_Params

Contains
  !----------------------------------------------------------------
  ! Initialize the parameters for an MEP simulation
  !----------------------------------------------------------------
  Subroutine mep_init(mepparams, simcell, sorbates, ctrl_filename, opt_meptag)
    Type(MEP_Params), Intent(out)  :: mepparams
    Type(Simcell_Params), Intent(in) :: simcell
    Type(AtMolCoords), Dimension(:), Intent(inout), Target :: sorbates
    Character(*), Intent(in) :: ctrl_filename
    Character(*), Optional, Intent(in) :: opt_meptag

    Character(len=strLen)    :: tag, line, nvtmc_tag, normalline, spacingline
    Character(len=strLen)    :: sorbname
    Integer          :: unitno, lineno, error, i, j, sorbtype
    Real(kind=RDbl)  :: ptx, pty, ptz
    Character(len=strLen)    :: gcmodeltype, junk

    If (Present(opt_meptag)) Then
      tag = opt_meptag
    Else
      tag = default_mep_tag
    End If

    !** Open the ctrl_file if it is not opened
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
      Open(file=ctrl_filename, unit=unitno)
    Endif
    
    !** Find the MEP section
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif

    !** Read the rest of the stuff
    Read(unitno,*) mepparams%nplanes

    !** Get the spacing between the planes and the displacement
    !** vector to go from one plane to the next
    Allocate(mepparams%spacing(mepparams%nplanes), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate memrory for 'spacing'"
      Stop
    End If
    Allocate(mepparams%dispvect(mepparams%nplanes), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate memrory for 'dispvect'"
      Stop
    End If
    Read(unitno,'(a)') spacingline
    Call mep_getspacing(mepparams, spacingline)

    Read(unitno,*) mepparams%initmctag
    Read(unitno,*) mepparams%annealtag
    Read(unitno,*) mepparams%basefile
    Read(unitno,*)

    !** Initialize the different monte carlo simulations
    nvtmc_tag = mepparams%initmctag
    Call nvtmc_init(mepparams%initmcparams, simcell, ctrl_filename, &
        nvtmc_tag, lineno)
    nvtmc_tag = mepparams%annealtag
    Call nvtmc_init(mepparams%annealparams, simcell, ctrl_filename, &
        nvtmc_tag, lineno)
    ! Set the file position to the appropriate line number
    Call skiplines(unitno, lineno+6)

    !** Allocate memory for the each plane and initialize it
    Allocate(mepparams%planes(mepparams%nplanes), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate memrory for 'planes'"
      Stop
    End If
    Call mepplanes_init(mepparams%planes, sorbates, &
        mepparams%annealparams, unitno)

    !** Set the constraint sorbate type and index
    mepparams%csorbtype  = mepparams%planes(1)%csorbtype
    mepparams%csorbindex = mepparams%planes(1)%csorbindex
    mepparams%nsorbs     = mepparams%planes(1)%path(1)%nsorbs

    !** Set the default planeno to a negative number indicating it has
    !** not been initialized
    mepparams%currentplane = -1
  End Subroutine mep_init

  !----------------------------------------------------------------
  ! Parses the spacing line to set the spacing between the planes
  ! If the value is specified in the control file then all the
  ! planes have the same spacing between them otherwise the spacing
  ! between planes 1,2 is in mepparams%spacing(1) and so on
  !---------------------------------------------------------------
  Subroutine mep_getspacing(mepparams, spacingline)
    Type(MEP_Params), Intent(inout)  :: mepparams
    Character(*), Intent(in)         :: spacingline

    Character(len=strLen)  :: stripline
    Character(len=strLen), Dimension(10)   :: strfields
    Integer                :: nfields, npts, i
    Real(kind=RDbl)        :: spacing, dx, dy, dz

    stripline = stripcmnt(spacingline)
    nfields = split(stripline, strfields, ",")

    !** Check to see if the number of fields is greater than one.
    !** If it is read from the filename specified
    If (nfields == 1) Then
      ! Read from a file
      mepparams%spacingfile = strfields(1)
      Call mep_readspacing(mepparams, mepparams%spacingfile)
    Else
      ! Generate the spacing and the displacement normal
      spacing = toreal(strfields(1))
      dx      = toreal(strfields(2))
      dy      = toreal(strfields(3))
      dz      = toreal(strfields(4))
      Do i=1, mepparams%nplanes
        mepparams%spacing(i) = spacing
        mepparams%dispvect(i) = (/dx, dy, dz/)
      End Do
    End If
  End Subroutine mep_getspacing

  !-----------------------------------------------------------
  ! This routine reads in the values of the normals from the
  ! file "filename"
  !-----------------------------------------------------------
  Subroutine mep_readspacing(mepparams, filename)
    Type(MEP_Params), Intent(inout)   :: mepparams    
    Character(*), Intent(in)           :: filename

    Integer               :: unitno, error, npts, i
    Real(kind=RDbl)       :: spacing, dx, dy, dz
    
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
    If (npts /= mepparams%nplanes) Then
      Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
          " The no. of points in the file ", Trim(filename), &
          " does not match that in the control file"
      Stop
    End If
    Do i=1, npts
      Read(unitno, *) spacing, dx, dy, dz
      mepparams%spacing(i) = spacing
      mepparams%dispvect(i) = (/dx, dy, dz/)
    End Do
    Close(unitno)
    Return
  End Subroutine mep_readspacing

  !--------------------------------------------------------------
  ! Get the no. of simulations
  !--------------------------------------------------------------
  Integer Function mep_getnoplanes(mepparams)
    Type(MEP_Params), Intent(in)    :: mepparams
    mep_getnoplanes = mepparams%nplanes
    Return
  End Function mep_getnoplanes

  !--------------------------------------------------------------
  ! This routine sets the constraint plane to the plane specified
  ! by "simno" and moves the molecules to that plane
  !--------------------------------------------------------------
  Subroutine mep_movetoplane(mepparams, sorbates, simcell, planeno)
    Type(MEP_Params), Intent(inout)  :: mepparams
    Type(AtMolCoords), Intent(inout), Dimension(:) :: sorbates
    Type(SimCell_Params), Intent(inout) :: simcell
    Integer, Intent(in)              :: planeno

    Integer              :: i, constrindex, csorbtype, sorbtype
    Real(kind=RDbl)      :: spacing, avgx, avgy, avgz, pot
    Type(VecType)        :: normal, planept, displace, currentpt, avgpos
    Logical              :: fast, mapflag
    Type(VecType)        :: dispvect

    !** Set the current planeno
    mepparams%currentplane = planeno

    !** Set the constrain plane based on the simulation number
    constrindex = mepparams%csorbindex
    normal      = mepplanes_getnormal(mepparams%planes(planeno))
    Call mcparams_setconstrnormal( &
        mepparams%initmcparams%nvtmcsorbs(constrindex)%gcnvtparams, normal)
    Call mcparams_setconstrnormal( &
        mepparams%annealparams%nvtmcsorbs(constrindex)%gcnvtparams, normal)
    
    !** Do not do anything else if it is the first plane
    If (planeno == 1) Then
      Return
    End If

    !** Copy the initial coords of the currentplane from the previous
    !** plane
    Call mepplanes_copyfinalgcoords(mepparams%planes(planeno-1), &
        mepparams%planes(planeno))
    
    !** Displace the constrained molecule to the current plane
    spacing  = mepparams%spacing(planeno-1)
    dispvect = mepparams%dispvect(planeno-1)
    Call mepplanes_displace(mepparams%planes(planeno), spacing, dispvect)
  End Subroutine mep_movetoplane

  !------------------------------------------------------------------------
  ! This function sets the initial coordinates of the current path
  ! and stores them in the "sorbates" structure from the "planes"
  ! field within "mepparams".  The "sorbates" structure is used for
  ! the annealing. "ierr" is set to 1 if we run into an untabulated region
  !------------------------------------------------------------------------
  Subroutine mep_movetopath(mepparams,sorbates,simcell,planeno,pathno,ierr)
    Type(MEP_Params), Intent(inout)  :: mepparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(SimCell_Params), Intent(inout) :: simcell
    Integer, Intent(in)              :: planeno, pathno
    Integer, Intent(out)             :: ierr
    
    Integer        :: i, natoms, sorbtype, nextpath, triedpaths, npaths
    Real(kind=RDbl)         :: pot
    Type(GeneralizedCoords) :: gcoords
    Logical                 :: fast, mapflag
    
    ierr = 0
    nextpath = pathno
    triedpaths = 0
    npaths = mepparams%planes(planeno)%npaths
    Do
      !** Generate the xyz coordinates from the generalized coordinates store
      !** in the initialcoords of the appropriate plane and path.
      Do i=1, mepparams%nsorbs
        sorbtype = nvtmc_getsorbtype(mepparams%annealparams, i)
        natoms   = molecules_getnatoms(sorbtype)
        
        ! Copy the generalized coordinates
        sorbates(sorbtype)%gcoords(1) = &
            mepparams%planes(planeno)%path(nextpath)%initialcoords(i)
        
        ! Generate the xyz principal coordinates
        Call gcmodels_toxyz(sorbates(sorbtype)%gcoords(1), &
          sorbates(sorbtype)%coords(1:natoms,1)%rp, &
          molecules_getcomdefn(sorbtype))

        ! Generate the other coordinates
        Call simcell_pbc(simcell, sorbates(sorbtype)%coords(1:natoms,1)%rp, &
            sorbates(sorbtype)%coords(1:natoms, 1)%r, &
            sorbates(sorbtype)%coords(1:natoms, 1)%cr)
      End Do
    
      !** Calculate the initial energy of the system and see if we
      !** are still in the tabulated region.  If not then try picking
      !** the coordinates of a different path.  We cycle through all
      !** all the paths until we either exhaust all paths or find one
      fast = .False.
      Call forcefield_getssint(sorbates, simcell, fast, pot, mapflag)
      If (mapflag) Then
        Return
      Else
        Write(0,'(1x,2a,i4, a, i4)') __FILE__," : ",__LINE__, &
            " We are in an untabulated region for path :", nextpath
        ierr = nextpath - pathno
        nextpath = Mod(nextpath, npaths)
        nextpath = nextpath + 1
        triedpaths = triedpaths + 1
      End If
      If (triedpaths == npaths) Exit
    End Do
    
    !** We did not find any paths
    Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
        " Could not proceed further with any paths", &
        " All paths lead into the walls"
    Stop
  End Subroutine mep_movetopath


  !--------------------------------------------------------------
  ! Does the actual moves
  !--------------------------------------------------------------
  Subroutine mep_dosim(mepparams, sorbates, simcell)
    Type(MEP_Params), Intent(inout)    :: mepparams
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(inout)  :: simcell

    Integer    :: i, j, pathno, nsims, simno, sorbtype, ierr
    Integer, Dimension(MAX_SORBS) :: lastposunits
    Integer    :: currentplane, nsorbs, nfields, nrgunit
    Type(VecType)                 :: currentpt, displace
    Character(len=strLen)         :: molecname
    Character(len=lstrLen)        :: avgposstr
    Character(len=strLen), Dimension(10) :: fields
    Real(kind=RDbl)     :: dist
    Type(GeneralizedCoords) :: gcoords

    !** Check that we have actually moved to the right plane
    If (mepparams%currentplane < 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Please move to the correct plane before doing the minimization"
      Stop
    End If

    !** Get some sundry information
    currentplane = mepparams%currentplane
    nsorbs       = mepparams%nsorbs

    !** Open the output files
    Call mep_openoutfiles(mepparams, lastposunits, nrgunit)
    
    !** Do the simulation repeatedly in each plane
    !** and get the average final position
    Do pathno=1, mepparams%planes(currentplane)%npaths !No. of iterations/plane
      !** Move to the next path in the current plane
      Call mep_movetopath(mepparams, sorbates, simcell, currentplane, &
          pathno, ierr)
      ! Check if we ran into a wall following this path.  If so discard
      ! it and continue with the next path
      If (ierr /= 0) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
            " Discarding path number ", pathno
        Cycle
      End If
      
      !** Randomize the initial placement of the molecules by doing 
      !** a short NVT Monte Carlo and write the crash file
      Call nvtmc_dosim(mepparams%initmcparams, sorbates, simcell, 1)

      !** Do the simulated annealing
      nsims = nvtmc_getnosims(mepparams%annealparams)
      Do simno = 1, nsims   ! Iterations at different temperatures
        Call nvtmc_initcounters(mepparams%annealparams, simno)
        Call nvtmc_dosim(mepparams%annealparams, sorbates, simcell, simno)
        Call nvtmc_displaystats(mepparams%annealparams, simno, 0)
      End Do

      !** Make sure the lastposition for each plane is in the constrained
      !** plane.  It is possible that the point is slightly off due to
      !** the finite machine precision, so we also want to "nudge" it
      !** back into the plane
      Call mepplanes_setfinalcoords(sorbates, mepparams%annealparams, &
          mepparams%planes(currentplane), pathno)

      !** Write the crash file after the annealing
      Call mepplanes_writemeprestartfile(mepparams%planes(currentplane), &
          sorbates, pathno)

      !** Get the final position of the molecules. Write the
      !** final position to the "lastpos" file and the energies
      Do i=1, nsorbs
        sorbtype = nvtmc_getsorbtype(mepparams%annealparams, i)
        avgposstr = gcmodels_display( &
            mepparams%planes(currentplane)%path(pathno)%finalcoords(i), "f8.3")
        Write(lastposunits(i), '(i6, a54)') pathno, Trim(avgposstr)
      End Do
      Call mep_dumpenergy(mepparams, currentplane, pathno, nrgunit)
    End Do   ! No. of paths/plane
    
    !** Close the files
    Do i=1, nsorbs
      Close(unit=lastposunits(i))
    End Do
  End Subroutine mep_dosim

  !----------------------------------------------------------------
  ! Open the files to save the last positions of the different
  ! sorbates and the energies
  !----------------------------------------------------------------
  Subroutine mep_openoutfiles(mepparams, lastposunits, nrgunit)
    Type(MEP_Params), Intent(in)    :: mepparams
    Integer, Dimension(:), Intent(inout)  :: lastposunits
    Integer, Intent(inout)                :: nrgunit

    Integer    :: i, j, sorbtype
    Integer    :: currentplane, nsorbs, simstart
    Type(VecType)                 :: normal
    Character(len=strLen)         :: lastposbasefile, lastposfile
    Character(len=strLen)         :: energyfile, molecname

    !** Get some sundry information
    currentplane = mepparams%currentplane
    nsorbs       = mepparams%nsorbs
    simstart     = general_getsimstart()

    !** Generate the correct file extensions
    lastposbasefile = Trim(mepparams%basefile)//".lastpos"
    energyfile      = Trim(mepparams%basefile)//".nrg"
    
    energyfile  = genfilename(energyfile, simstart+currentplane-1)
    nrgunit     = file_getunit(energyfile)
    Open(unit=nrgunit, file=energyfile)
    ! Open the one file for each sorbate type to store the final positions
    ! in a plane
    Do i=1, nsorbs
      sorbtype   = nvtmc_getsorbtype(mepparams%annealparams, i)
      molecname  = molecules_name(sorbtype)
      lastposfile = Trim(lastposbasefile)//"."//molecname
      lastposfile = genfilename(lastposfile, simstart+currentplane-1)
      lastposunits(i) = file_getunit(lastposfile)
      Open(unit=lastposunits(i), file=lastposfile)
    End Do
  End Subroutine mep_openoutfiles

  !---------------------------------------------------------
  ! This routine writes the current energy to "unitno"
  !---------------------------------------------------------
  Subroutine mep_dumpenergy(mepparams, planeno, pathno, unitno)
    Type(MEP_Params), Intent(inout) :: mepparams
    Integer, Intent(in) :: planeno, pathno
    Integer, Intent(in) :: unitno

    Integer, Save :: currentplane = 0
    Integer       :: nsorbs, i, j, simstart
    Real(kind=RDbl)    ::   coulnrg, noncoulnrg, totpairnrg, totnrg
    
    nsorbs    = molecules_getnsorbs()
    simstart  = general_getsimstart()

    !** Check if this the first call for this plane. If so dump
    !** the header
    If (planeno /= currentplane) Then
      currentplane = planeno
      Write(unitno, '(a,9x)', Advance='No') "#"
      Do i=1, nsorbs-1
        Do j=i+1, nsorbs
          Write(unitno,'(5x,i1,a,i1,4x)', Advance='No') i,"-",j
        End Do
      End Do
      Write(unitno,'(a12)') "Total(kJ)"
    End If

    !** Dump the energies
    totnrg = 0.0_RDbl
    Write(unitno,'(2i5)', Advance='No') planeno+simstart-1, pathno
    Do i=1, nsorbs-1
      Do j=i+1, nsorbs
        noncoulnrg = molecules_getnoncoul(i,j,"inst")
        coulnrg = molecules_getcoul(i,j,"inst")
        totpairnrg = (coulnrg + noncoulnrg)
        totnrg = totnrg + totpairnrg
        mepparams%planes(planeno)%path(pathno)%totnrg = totnrg
        Write(unitno,'(f12.4)', Advance='No') totpairnrg
      End Do
    End Do
    Write(unitno,'(f12.4)') totnrg
  End Subroutine mep_dumpenergy

  !----------------------------------------------------------
  ! Display the mep parameters for the simulation at plane
  ! "planeno"
  !----------------------------------------------------------
  Subroutine mep_displaysimparams(mepparams, planeno, nspc, optunit)
    Type(MEP_Params), Intent(in) :: mepparams
    Integer, Intent(in)          :: planeno
    Integer, Intent(in)          :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer    :: i, unitno, sorbtype
    
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
    Write(unitno, '(2a)') spc, "The MEP Simulation Section:"
    Write(unitno, '(a,2x,a,i6)') spc,"Plane No.     : ", planeno
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "The Initial Monte Carlo Parameters:"
    Call nvtmc_displaysimparams(mepparams%initmcparams, nspc+2, unitno)
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "The Annealing Parameters:"
    Call nvtmc_displaysimparams(mepparams%annealparams, nspc+2, unitno)
    
  End Subroutine mep_displaysimparams


  !---------------------------------------------
  ! Display the mep parameters
  !---------------------------------------------
  Subroutine mep_display(mepparams, nspc, optunit)
    Type(MEP_Params), Intent(in) :: mepparams
    Integer, Intent(in) :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer    :: i, unitno
    
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
    Write(unitno, '(2a)') spc, "The MEP Parameters Section:"
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a,i6)') spc,"No. of simulations        : ", &
        mepparams%nplanes
    Write(unitno, '(a,2x,a,i4)') spc, "No. of sorbates           : ", &
        mepparams%nsorbs
    Write(unitno, '(a,2x,a,a)') spc, "Constrained sorbate       : ", &
        Trim(molecules_name(mepparams%csorbtype))
    Write(unitno, '(a, 2x, a9, a24, 8x, a8)') spc, &
        "Plane No.", "Disp. Vector", "Spacing"
    Do i=1, mepparams%nplanes
      Write(unitno, '(a, 2x, i5, a28, 8x, f8.3)') spc, i, &
          Trim(vector_display(mepparams%dispvect(i), "f8.3")),  &
          mepparams%spacing(i)
    End Do
    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "The Planes Section:"
    Call mepplanes_display(mepparams%planes, nspc+2, unitno)

    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "The Initial Monte Carlo Section:"
    Call nvtmc_display(mepparams%initmcparams, nspc+2, unitno)

    Write(unitno, '(a,2x,a)') spc, dashedline
    Write(unitno, '(a,2x,a)') spc, "The Annealing Section:"
    Call nvtmc_display(mepparams%annealparams, nspc+2, unitno)
    
  End Subroutine mep_display

End Module mep

