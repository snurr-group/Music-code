Module stdescent
  Use defaults, Only: lstrLen, strLen, RDbl, MAX_DOF, &
      caltoj, MAX_SORBS, dashedline
  Use utils, Only: arrnorm, comb, isfileopen, filesrchstr, toupper, &
      split, toint, toreal, multarrvec, multvecvec, dispvec
  Use matrix, Only: invertarray
  Use file, Only: file_getunit
  Use config, Only: AtMolCoords, config_allocfields
  Use molecules, Only: molecules_gettype, molecules_getnsorbs, &
      molecules_getnatoms, molecules_name, molecules_getdof
  Use gcmodels, Only: gcmodels_display, gcmodels_pbc
  Use mcparams, Only:
  Use bakersubs, Only: bakersubs_getenergy, bakersubs_calchess, &
      bakersubs_calcmetrictensor, bakersubs_calcgrad, bakersubs_xtosorb
  Use sort, Only: SortArray_Type, sort_init, sort_heapsort
  Use smap, Only: Smap_Params, smap_getPosSiteType, smap_display, smap_init
  Use search, Only: SearchArray_Type, MAX_SRCHKEYS, search_binarysearch, &
      search_init
  Use histogram, Only: Histogram_Params, histogram_init, histogram_update, &
      histogram_display
  Use simcell, Only: SimCell_Params, simcell_pbc
  Use vector, Only: VecType, mag

  Implicit None
  Save

  Private
  Public :: STDESCENT_Params, stdescent_init, MININFO_Params, &
      stdescent_dosim, stdescent_display

  Integer, Parameter  :: MAX_PTS = 40000
  Integer, Parameter  :: MAX_NRGINTER = 15

  Character(len=strLen)    :: default_stdescent_tag = &
      "Baker Minimum Energy Path"

  Type PATHINFO_Params
    Real(kind=RDbl), Dimension(MAX_DOF)      :: df
    Real(kind=RDbl), Dimension(MAX_NRGINTER) :: nrg
  End Type PATHINFO_Params

  Type MININFO_Params
    Real(kind=RDbl), Dimension(MAX_DOF) :: df
    Real(kind=RDbl), Dimension(MAX_NRGINTER) :: nrg
    Integer                  :: recordno       
    Integer                  :: sitetype
  End Type MININFO_Params

  Type STDESCENT_Params
    Logical                :: histflag
    Logical                :: descent
    Character(len=strLen)  :: inbasefile
    Character(len=strLen)  :: outbasefile
    Real(kind=RDbl)        :: minActNrg
    Real(kind=RDbl)        :: stepsize
    Real(kind=RDbl)        :: gradtol
    Real(kind=RDbl)        :: energytol
    Real(kind=RDbl)        :: mintol
    Real(kind=RDbl)        :: minspacing
    Real(kind=RDbl)        :: minminnrg, maxminnrg
    Real(kind=RDbl)        :: minsdptnrg, maxsdptnrg
    Integer                :: maxsteps
    Real(kind=RDbl)        :: binsize
    Integer                :: nsorbs
    Integer, Dimension(:), Pointer :: mtypes
    Integer                :: sitesorbtype
    Integer                :: dof
    Integer                :: sortindx   ! indx in df(1:dof) on which to sort
    Integer                :: nrgnfields ! No. of energy fields
    Integer                :: nminpts
    Integer                :: nsdlpts
    Type(Smap_Params), Pointer :: smap
    Character(len=strLen)      :: smapfile
    Character(len=strLen)      :: sitemolec
    Type(MININFO_Params), Dimension(:), Pointer :: minima
    Type(MININFO_Params), Dimension(:), Pointer :: sdlpts
    Type(SearchArray_Type), Dimension(:), Pointer :: minsrch
  End Type STDESCENT_Params

  Interface display
    Module Procedure stdescent_display
  End Interface


Contains
  !------------------------------------------------------------------
  ! Initialize the parameters for getting the Minimum Energy Paths
  ! using the Baker's algorithm.  Get the parameters from the control
  ! file "ctrl_filename" and store them in the "stdescentsparams" object
  !------------------------------------------------------------------
  Subroutine stdescent_init(stdparams, sorbates, simcell, ctrl_filename, &
      opt_stdescenttag)
    Type(STDESCENT_Params), Intent(out)  :: stdparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(SimCell_Params), Intent(inout) :: simcell
    Character(*), Intent(in)            :: ctrl_filename
    Character(*), Intent(in),Optional:: opt_stdescenttag

    Integer                 :: ierr, recordno, ptno, sysnsorbs, sitesorbtype
    Real(kind=RDbl)         :: pot, tk, totnrg
    Character(len=strLen)   :: tag, sorbname, sourcetype, filename, gcmodeltype
    Character(len=lstrLen)  :: line
    Integer                 :: lineno, unitno, i, j, pertlineno, natoms
    Integer                 :: bakerlineno, sorbtype, dof, error, npts, index
    Integer                 :: nsorbs, ins_sorbs, sortindx, nrgnfields
    Logical                 :: mapflag
    Integer                 :: descflag, makehist
    Real(kind=RDbl), Dimension(MAX_DOF) :: df
    
    If (Present(opt_stdescenttag)) Then
      tag = opt_stdescenttag
    Else
      tag = default_stdescent_tag
    End If
    
    !** Open the ctrl_file if it is not opened
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
      Open(file=ctrl_filename, unit=unitno)
    Endif
    
    !** Find the Steepest Descent section
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif
    
    !** Read the rest of the stuff
    Read(unitno,*) descflag
    If (descflag == 1) Then
      stdparams%descent = .True.
    Else
      stdparams%descent = .False.
    End If
    Read(unitno,*) makehist
    If (makehist == 1) Then
      stdparams%histflag = .True.
    Else
      stdparams%histflag = .False.
    End If
    Read(unitno,*) stdparams%minActNrg
    stdparams%minActNrg = stdparams%minActNrg
    Read(unitno,*) stdparams%inbasefile
    Read(unitno,*) stdparams%outbasefile
    Read(unitno,*) stdparams%stepsize
    Read(unitno,*) stdparams%gradtol
    Read(unitno,*) stdparams%energytol
    Read(unitno,*) stdparams%mintol
    Read(unitno,*) stdparams%minspacing
    Read(unitno,*) stdparams%maxsteps
    Read(unitno,*) stdparams%binsize
    Read(unitno,*) stdparams%smapfile
    Read(unitno,*) stdparams%sitemolec
    Read(unitno,*) nsorbs
    Read(unitno,*)
    sitesorbtype = molecules_gettype(stdparams%sitemolec)
    stdparams%sitesorbtype = sitesorbtype
    stdparams%nsorbs = nsorbs
    ! Calculate the no. of energy fields. This should be the no. of ways
    ! we can take sets of two out of the total no. of sorbates plus one
    ! for the total energy
    sysnsorbs   = molecules_getnsorbs()
    stdparams%nrgnfields =  comb(sysnsorbs, 2) + 1
    dof = 0
    
    !** Allocate space for the molecule types
    Allocate(stdparams%mtypes(nsorbs),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i6,a,i4)') __FILE__,":",__LINE__, &
          " Could not allocate stdparams%mtypes of size ",nsorbs
      Stop
    End If

    Do i=1, nsorbs
      Read(unitno,*) sorbname
      sorbtype = molecules_gettype(sorbname)
      stdparams%mtypes(i) = sorbtype
      If (sorbtype == sitesorbtype) Then
        stdparams%sortindx = dof+1
      End If
      dof = dof + molecules_getdof(sorbtype)
    End Do
    stdparams%dof = dof
    stdparams%nsorbs = nsorbs

    !** Make sure that atleast one of the running flags(descent, histogram)
    !** is set to be true
    If ((.Not. stdparams%histflag) .And. (.Not. stdparams%descent)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Please set either the descent flag and/or the histogram flag to 1"
      Stop
    End If
    
    !** Read the sitemap if it is not null
    If (toupper(stdparams%smapfile) /= 'NULL') Then
      Call smap_init(stdparams%smap, simcell, stdparams%smapfile)
    End If

    !** Initialize the config structure.  We are using a CUSTOM
    !** restart file
    Do i=1, nsorbs
      sorbtype = stdparams%mtypes(i)
      Call config_allocfields(sorbates(sorbtype), sorbtype, 1)
      sorbates(sorbtype)%natoms     = molecules_getnatoms(sorbtype)
      sorbates(sorbtype)%nmoles     = 1
      sorbates(sorbtype)%sourcetype = "CUSTOM"
      sorbates(sorbtype)%filename   = Trim(stdparams%inbasefile)//".sdpt"
    End Do

    !** Sort the two input files, containing the saddle points and minima, 
    !** based on their energies and generate ".srtd" files
    Call stdescent_initsites(stdparams, ".sdpt", stdparams%sdlpts, &
        stdparams%minsdptnrg, stdparams%maxsdptnrg, stdparams%nsdlpts)
    Call stdescent_initsites(stdparams, ".min", stdparams%minima, &
        stdparams%minminnrg, stdparams%maxminnrg, stdparams%nminpts)

    !** Initialize the minsearch array
    npts = stdparams%nminpts
    sortindx = stdparams%sortindx
    nrgnfields = stdparams%nrgnfields
    Allocate(stdparams%minsrch(npts), Stat=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate 'stdparams%minsrch'"
      Stop
    End If
    ! Initialize the search array by giving it a list of values
    Call search_init(stdparams%minima(1:npts)%nrg(nrgnfields), &
        stdparams%minima(1:npts)%df(sortindx), &
        stdparams%minima(1:npts)%df(sortindx+1), &
        stdparams%minima(1:npts)%df(sortindx+2), stdparams%minsrch, npts)
  End Subroutine stdescent_init

  !-------------------------------------------------------------------------
  ! This routine initializes the minima and the transition state sites.
  ! It reads the sites from the file inbasefile//ext.  It creates
  ! the sorted array "pts" and returns the no. of pts in the array("npts").
  ! It also creates a file inbasefile//".srtd"//.ext
  !-------------------------------------------------------------------------
  Subroutine stdescent_initsites(stdparams, ext, srtdpts, min, max, npts)
    Type(STDESCENT_Params), Intent(inout) :: stdparams
    Character(*), Intent(in)  :: ext
    Type(MININFO_Params), Dimension(:), Pointer   :: srtdpts
    Real(kind=RDbl), Intent(out)       :: min, max
    Integer, Intent(out)               :: npts

    Character(len=strLen)  :: infilename, outfilename, molecname
    Character(len=strLen)  :: commentline, sorbname
    Character(len=strLen), Dimension(10) :: fields1, fields2
    Character(len=2*lstrLen):: sorblist
    Integer, Dimension(MAX_SORBS) :: inunitno
    Integer                :: outunitno, iterno, molecdof, sitetype, index
    Integer                :: nsorbs, nfields1, nfields2, sorbtype, nrgfields
    Integer                :: nrgunit, oursorbtype, begdof, enddof, sortindx
    Type(MININFO_Params), Dimension(MAX_PTS) :: inppts
    Type(SortArray_Type), Dimension(MAX_PTS) :: sortarray
    Type(SearchArray_Type), Dimension(MAX_PTS) :: srcharray
    Integer     :: i, j, dof, error, recordno, sitesorbtype, mtype
    Real(kind=RDbl), Dimension(MAX_DOF)      :: df
    Real(kind=RDbl), Dimension(MAX_NRGINTER) :: nrg
    Real(kind=RDbl), Dimension(MAX_PTS, MAX_NRGINTER) :: nrglist
    Type(VecType)   :: com
    
    !** Initialize some stuff
    dof         = stdparams%dof
    nsorbs      = stdparams%nsorbs
    nrgfields   = stdparams%nrgnfields
    sitesorbtype = stdparams%sitesorbtype
    sortindx     = stdparams%sortindx

    !** Open the various input, output and energy files
    !** of this file we can sort the coordinates files
    ! Energy file
    infilename = Trim(stdparams%inbasefile)//Trim(ext)//".nrg"
    nrgunit = file_getunit(infilename)
    Open(unit = nrgunit, file=infilename, status='old', IOSTAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          "  Could not open file ", Trim(infilename)
      Stop
    End If
    ! Open the input files
    Do i=1, nsorbs
      mtype = stdparams%mtypes(i)
      molecname = molecules_name(mtype)
      infilename  = &
          Trim(stdparams%inbasefile)//Trim(ext)//"."//Trim(molecname)
          
      inunitno(i) = file_getunit(infilename)
      Open(unit = inunitno(i), file=infilename, status='old', IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            "  Could not open file ", Trim(infilename)
        Stop
      End If
    End Do
    ! Open the output file
    outfilename = &
        Trim(stdparams%inbasefile)//".srtd"//Trim(ext)
    outunitno = file_getunit(outfilename)
    Open(unit = outunitno, file=outfilename, IOSTAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          "  Could not open file ", Trim(outfilename)
      Stop
    End If

    !** Read the energy file first
    npts = 0
    ! Read the file line which maps the sorbate names and sorbate types
    ! We want to make sure our sorbate name and sorbate types have the
    ! same correspondence to the file we are reading
    Read(nrgunit,'(a)') sorblist
    nfields1 = split(sorblist, fields1, ":")
    Do i=1, nfields2
      nfields2 = split(fields1(i), fields2, "=")
      sorbtype = toint(fields2(1))
      sorbname = fields2(2)
      oursorbtype = molecules_gettype(sorbname)
      If (oursorbtype /= sorbtype) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
            " The order of sorbates do not match those in the input file"
        Stop
      End If
    End Do

    Read(nrgunit,*) ! Read the column headers
    Do
      Read(nrgunit, *, IOSTAT=error) iterno, (nrg(j), j=1, nrgfields)
      If (error /= 0) Exit
      If (npts == MAX_PTS) Then
        Write(0,'(1x,2a,i4,a,i7,a)') __FILE__," : ",__LINE__, &
            " Reading the first ", MAX_PTS, " only"
        Exit
      End If
      npts = npts + 1
      nrglist(npts, 1:nrgfields) = nrg(1:nrgfields)

      ! Find the minimum and the maximum energies
      If (npts == 1) Then
        min = nrg(nrgfields)
        max = nrg(nrgfields)
      Else
        If (nrg(nrgfields) < min) Then
          min = nrg(nrgfields)
        End If
        If (nrg(nrgfields) > max) Then
          max = nrg(nrgfields)
        End If
      End If
    End Do

    !** Allocate the coordinates array based on the energy read in
    Allocate(srtdpts(npts), Stat=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate 'srtdpts'"
      Stop
    End If

    !** Read the various input files
    ! Combine the positions from the various output files
    begdof = 1
    Do i=1, nsorbs
      mtype = stdparams%mtypes(i)
      molecdof = molecules_getdof(mtype)
      Read(inunitno(i), '(a)') commentline
      enddof = begdof + molecdof - 1
      Do j=1, npts
        Read(inunitno(i),*, Iostat=error) recordno, df(1:molecdof)
        If (error /= 0) Then
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              " Error reading inputfile for sorbate ", &
              Trim(molecules_name(mtype))
          Write(0, '(a,i6)') "No. of points read: ", j
          Stop
        End If
        inppts(j)%df(begdof:enddof) = df(1:molecdof)
        inppts(j)%recordno = j
        If (sitesorbtype == i) Then
          com = df(1:3)
          inppts(j)%sitetype = smap_getPosSiteType(stdparams%smap, com)
        End If
      End Do
      begdof = enddof + 1
    End Do
    ! Set the energy
    Do j=1, npts
      inppts(j)%nrg(1:MAX_NRGINTER) = nrglist(j,1:MAX_NRGINTER)
    End Do

    !** Sort the input based on the "pot" field and then on the basis 
    !** of the position of molecule used for siting
    Call sort_init(inppts(1:npts)%nrg(nrgfields), &
        inppts(1:npts)%df(sortindx), inppts(1:npts)%df(sortindx+1), &
        inppts(1:npts)%df(sortindx+2), sortarray, npts)
    Call sort_heapsort(sortarray(1:npts), npts)
    Do i=1, npts
      srtdpts(i) = inppts(sortarray(i)%index)
    End Do

    !** Write the sorted output
    Write(outunitno,'(a)') commentline
    Do i=1, npts
      Write(outunitno,'(i6,i7)', Advance='No') i, srtdpts(i)%recordno
      Write(outunitno, '(i4)', Advance='No') srtdpts(i)%sitetype
      Do j=1, dof
        Write(outunitno, '(f9.4)', Advance='No') srtdpts(i)%df(j)
      End Do
      Write(outunitno, '(f13.7)') srtdpts(i)%nrg(nrgfields)
    End Do

    !** Finish up
    Close(nrgunit)
    Do i=1, nsorbs
      Close(inunitno(i))
    End Do
    Close(outunitno)
  End Subroutine stdescent_initsites

  !------------------------------------------------------------------------
  ! This routine reads the input file containing the saddlepts.  It
  ! gets the generalized coordinates of the saddlepoint in the array
  ! "df" and returns the record number and the total energy in "recordno"
  ! and "totnrg" respectively. "ierr" is set to a non-zero number when
  ! the end-of-file or any other error in reading the file occurs
  !-----------------------------------------------------------------------
  Subroutine stdescent_readinputfile(stdparams, df, recordno, totnrg, ierr)
    Type(STDESCENT_Params), Intent(in)        :: stdparams
    Integer, Intent(out)                      :: recordno
    Real(kind=RDbl), Dimension(:), Intent(out):: df
    Real(kind=RDbl), Intent(out)              :: totnrg
    Integer, Intent(out)        :: ierr
    
    Integer    :: unitno, error, dof
    Character(len=strLen) :: filename
    Logical, Save         :: firsttime = .True.

    !** If this the first time here open the file and read the file
    !** description
    If (firsttime) Then
      firsttime = .False.
      filename = Trim(stdparams%inbasefile)//".srtd.sdpt"
      unitno = isfileopen(filename)
      If (unitno < 0) Then
        unitno = file_getunit(filename)
        Open(unit = unitno, file=filename, status='old', IOSTAT=error)
        If (error /= 0) Then
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              "  Could not open file ", Trim(filename)
          Stop
        End If
      End If
      Rewind(unitno)
      ! Read the first line which is a comment line
      Read(unitno,*)
    End If

    !** Read the line and if we have reached the end close the file
    ierr = 0
    dof = stdparams%dof
    Read(unitno, *, IOSTAT=error) recordno, df(1:dof), totnrg
    If (error /= 0) Then
      ! We have reached the end of file
      ierr = error
      Close(unitno)
    End If
  End Subroutine stdescent_readinputfile

  !------------------------------------------------------------------------
  ! Reduce the no. of points
  !------------------------------------------------------------------------
  Integer Function stdescent_reducepts(pathinfo,  npts, dof, minspacing)
    Type(PATHINFO_Params), Dimension(:), Intent(inout) :: pathinfo
    Integer, Intent(in) :: npts, dof
    Real(kind=RDbl), Intent(in) :: minspacing

    Integer      :: i, j, lastconfig, rednpts
    Real(kind=RDbl) :: maxdiff
    Type(PATHINFO_Params), Dimension(Size(pathinfo, 1)) :: redpts
    
    lastconfig = 1
    rednpts = 0
    Do i=2, npts-1
      !** Get the largest change in the degrees-of-freedom between two
      !** consecutive configurations
      maxdiff = Abs(pathinfo(lastconfig)%df(1) - pathinfo(i)%df(1))
      Do j=2, dof
        If (maxdiff < Abs(pathinfo(lastconfig)%df(j) - pathinfo(i)%df(j))) Then
          maxdiff = Abs(pathinfo(lastconfig)%df(j) - pathinfo(i)%df(j))
        End If
      End Do
      
      !** Check if maxdiff is less than minspacing if not add the configuration
      If (maxdiff > minspacing) Then
        lastconfig = i
        rednpts = rednpts + 1
        redpts(rednpts) = pathinfo(i)
      End If
    End Do
    !** Add the last point since we want to keep the minimum in our list
    rednpts = rednpts + 1
    redpts(rednpts) = pathinfo(npts)
    
    !** Finally copy the temp array back to pathinfo
    pathinfo(1:rednpts) = redpts(1:rednpts)
    stdescent_reducepts = rednpts
  End Function stdescent_reducepts

  !-------------------------------------------------------------------
  ! This routine returns the vector along which we take the first step
  ! from the saddle point
  !-------------------------------------------------------------------
  Subroutine stdescent_firststepvec(stdparams, sorbates, simcell,df,vec,ierr)
    Type(STDESCENT_Params), Intent(in)             :: stdparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(inout)            :: simcell
    Real(kind=RDbl), Dimension(:), Intent(in) :: df
    Real(kind=RDbl), Dimension(:), Intent(out):: vec
    Integer, Intent(out)   :: ierr

    Integer      :: i, nsorbs, dof, negvals
    Real(kind=RDbl), Dimension(stdparams%dof) :: eigval, subdiag
    Real(kind=RDbl), Dimension(stdparams%dof,stdparams%dof):: eigvec, hess
    Real(kind=RDbl), Dimension(stdparams%dof,stdparams%dof):: G, Ginv, Prod
    Real(kind=RDbl), Dimension(stdparams%dof,stdparams%dof):: GinvH, GinvG

    !** Get some sundry information
    nsorbs = stdparams%nsorbs
    dof    = stdparams%dof
    ierr   = 0
    
    !** Get the hessian
    Call bakersubs_calcHess(sorbates,simcell,stdparams%mtypes,df,hess,ierr)
    If (ierr /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Problem in calculating the hessian"
      Stop
    Endif
    
    !** Get the eigenvalues of the Hessian matrix
    Call tred2(dof, dof, hess, eigval, subdiag, eigvec)
    Call tql2(dof, dof, eigval, subdiag, eigvec, ierr)
    If (ierr /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Could not get the eigenvalues for the first step"
      Stop
    Endif
    
    !** The vector is along the eigenvector corresponding to the negative
    !** eigenvalue
    ! Make sure only one eigenvector is negative i.e. we are actually on
    ! a saddle point
    negvals = 0
    Do i=1, dof
      If (eigval(i) < 0.0_RDbl) Then
        negvals = negvals + 1
      Endif
    End Do
    If (negvals /= 1) Then
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " We are not on a first order saddle point. ", &
          " Please check the following eigenvalues and eigenvectors"
      ! It is possible that due to machine precision it finds more than one 
      ! negative eigenvalue.  Since we are only storing 6 decimal values in the
      ! input file this could lead to some error in the calculation.
      Do i=1, dof
        Write(0,'(a, i5, g17.8)')    'eigval  :', i, eigval(i)
        Write(0,'(a)', Advance='No') 'eigvect :'
        Call dispvec(eigvec(1:dof, i), "g15.6", 0)
        Write(0,*)
      End Do
    Endif
    
    ! eigval has the values is ascending order so just take the eigenvector
    ! stored in the first column
    vec = eigvec(1:dof, 1)
  End Subroutine stdescent_firststepvec

  !------------------------------------------------------------------------
  ! Check if we are near a previously tabulated minimum.  In its current
  ! form it does NOT necessarily return the correct minimum i.e. it does
  ! not return the index with the same position but only the SAME ENERGY
  !------------------------------------------------------------------------
  Logical Function stdescent_foundminima(stdparams, simcell, energy, df, index)
    Use search
    Type(STDESCENT_Params), Intent(in)       :: stdparams
    Type(Simcell_Params), Intent(in)         :: simcell
    Real(kind=RDbl), Intent(in)              :: energy
    Real(kind=RDbl), Dimension(:), Intent(in):: df
    Integer, Intent(out)                     :: index
    
    Integer        :: begindx, endindx, npts, nrgnfields, niters
    Real(kind=RDbl):: totnrg1, totnrg2, mintol, nrgdiff
    Real(kind=RDbl), Dimension(MAX_SRCHKEYS) :: srchkey
    Real(kind=Rdbl):: dxpos, xpos1, xpos2
    Logical        :: lowerindxclosest
    Type(VecType)  :: com, pbccom

    !** Get some sundry information
    begindx = stdparams%sortindx
    npts    = stdparams%nminpts
    mintol  = stdparams%mintol
    nrgnfields = stdparams%nrgnfields
    endindx = begindx + 2

    ! Construct the search key from the energy and the position
    srchkey(1)   =  energy
    com          = df(begindx:endindx) 
    Call simcell_pbc(simcell, com, pbccom)
    srchkey(2:4) = pbccom 
    
    !** Do the search for the minima with coordinates given by df
    Call search_binarysearch(stdparams%minsrch, srchkey, npts, index)

    !** The search gives us an "index" where the "srchkey" would fit.
    !** How do we know whether the group search key belongs to is at
    !** "index" or "index-1" or none?  We can compare the srch energy
    !** with those at index-1 and index and see which one is closest
    !** and whether that is within a certain tolerance.
    lowerindxclosest = .False.
    If (index == 1) Then
      nrgdiff = Abs(stdparams%minima(index)%nrg(nrgnfields)-energy)
      lowerindxclosest = .False.
    Else If (index > npts) Then
      nrgdiff = Abs(stdparams%minima(index-1)%nrg(nrgnfields)-energy)
      lowerindxclosest = .True.
    Else
      totnrg1 = stdparams%minima(index-1)%nrg(nrgnfields)
      totnrg2 = stdparams%minima(index)%nrg(nrgnfields)
      If (Abs(totnrg1-energy) > Abs(totnrg2-energy)) Then
        nrgdiff = Abs(totnrg2-energy)
        lowerindxclosest = .False.
      Else
        nrgdiff = Abs(totnrg1-energy)
        lowerindxclosest = .True.
      End If
    End If

    ! Check if the nrdiff is within the tolerance if not we have not 
    ! found a minimum
    If (nrgdiff > mintol) Then
      stdescent_foundminima = .False.
      Return
    End If

    !** We have found a minimum.  Need to make sure it is the right one
    !** If the lower index is closer in energy than we have a round-off
    !** problem.
    If (lowerindxclosest) Then
      ! We probably have a round-off problem.
      srchkey(1) = stdparams%minima(index-1)%nrg(nrgnfields)
      Call search_binarysearch(stdparams%minsrch, srchkey, npts, index)
      stdescent_foundminima = .True.
    End If
    stdescent_foundminima = .True.
    Return
  End Function stdescent_foundminima

  !------------------------------------------------------------------------
  ! This routine does the steepest descent. 
  ! ierr = 1 if we descend into an untabulated region
  ! ierr = 2 if the energy does not go down with every step or we have a new
  !          minimum
  ! ierr = 3 if a minimum was not found in the max. no. of steps
  ! ierr = 4 if the eqs. of motion diverge because theta is close to zero
  ! ierr = 5 if we have a NEW minimum
  !-----------------------------------------------------------------------
  Subroutine stdescent_descend(stdparams, sorbates, simcell, df, &
      firstvec, pot, pathinfo, nsteps, index, ierr)
    Type(STDESCENT_Params), Intent(in)             :: stdparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(inout)            :: simcell
    Real(kind=RDbl), Dimension(:), Intent(inout)   :: df
    Real(kind=RDbl), Dimension(:), Intent(in) :: firstvec
    Real(kind=RDbl), Intent(out)              :: pot
    Type(PATHINFO_Params), Dimension(:)       :: pathinfo
    Integer, Intent(out)                      :: nsteps
    Integer, Intent(out)                      :: index
    Integer, Intent(out)                      :: ierr
    
    Integer            :: i, nsorbs, dof, maxsteps, nrgnfields
    Real(kind=RDbl)    :: stepsize, dtau, dssq, gradtol, lastpot
    Real(kind=RDbl), Dimension(stdparams%dof, stdparams%dof) ::  G, Ginv
    Real(kind=RDbl), Dimension(stdparams%dof)   ::  grad, GinvGrad, GGinvGrad
    Real(kind=RDbl), Dimension(stdparams%dof)   ::  dq

    !** Get some sundry information and initialize some stuff
    stepsize = stdparams%stepsize
    nsorbs   = stdparams%nsorbs
    dof      = stdparams%dof
    maxsteps = stdparams%maxsteps
    gradtol  = stdparams%gradtol
    nrgnfields = stdparams%nrgnfields
    nsteps   = 0
    ierr     = 0

    !** Get the initial energy
    Call bakersubs_calcgrad(sorbates,simcell,stdparams%mtypes, &
        df,pot,grad,ierr)
    pot = pot*caltoj
    If (ierr /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " We are in unknown territory"
      ierr = 1
      Return
    End If
    lastpot = pot

    !** Move a bit along the vector
    Do i=1, dof
      df(i) = df(i) + firstvec(i)*stepsize
    End Do

    !** Get the gradient and save the first point
    Call bakersubs_calcgrad(sorbates, simcell, stdparams%mtypes, &
        df, pot, grad, ierr)
    pot = pot*caltoj
    If (ierr /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " We are in unknown territory"
      ierr = 1
      Return
    End If
    nsteps = nsteps + 1
    pathinfo(nsteps)%df(1:dof) = df(1:dof)
    pathinfo(nsteps)%nrg(nrgnfields) = pot
    
    !** Start the descent
    Do
      ! Check the ending condition.  There are many ways to exit this 
      ! darn loop. 
      If (arrnorm(grad) < gradtol) Then
        ! We have found a minimum.  Check if it pretabulated
        If (stdescent_foundminima(stdparams, simcell, pot, df, index)) Then
          Return
        Else
          ierr = 5
          Return
        End If
      Else If (nsteps >= maxsteps) Then
        ! We did not find a minimum in the maximum no. of iterations
        ierr = 3
        Return
      Else If (pot > lastpot) Then
        ! We may be circling the minimum.  Check if we are near one
        If (stdescent_foundminima(stdparams, simcell, pot, df, index)) Then
          Return
        Else
          ! We have either found a new minimum or we actually have a strange
          ! case of the energy increasing even though we are doing a descent
          Write(0,*) ' ****** Maybe a new minima ********'
          Write(0,*) ' energy (kj): ', pot
          Write(0,*) ' df :', df(1:dof)
          Write(0,*) ' gradient norm :', arrnorm(grad)
          Write(0,*) ' nearest index :', index
          ierr = 2
          Return
        End If
      End If

      ! Set the lastpot to make sure we are going down in energy
      lastpot = pot

      ! Increment the step counter
      nsteps = nsteps + 1
      
      ! Get the Metric Tensor and its inverse
      Call bakersubs_calcMetricTensor(sorbates, simcell, &
          stdparams%mtypes, df, G, ierr)
      If (ierr /= 0) Then
        ! theta is close to zero and the equations of motion cannot be
        ! solved in this coordinate system.
        ierr = 4
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) 'theta is close to zero'
        Return
      End If
      Call invertarray(G, Ginv)
!!$      Call disparray(G, "MetricTensor: ", "f8.3")

      
      ! Calculate the stepsize. This is done according to the method
      ! described in Randy's paper.  The idea is that we actually set
      ! the approx length moved (ds) and from that back out the 
      ! "time"(dtau). See eqns. (13, 16) of Randy's paper
      Call multarrvec(Ginv, grad, GinvGrad)
      Call multarrvec(G, GinvGrad, GGinvGrad)
      Call multvecvec(GGinvGrad, GinvGrad, dssq)
      dtau = stepsize/Sqrt(dssq)

      ! Get the displacement vector dq
      dq = GinvGrad*dtau

      ! Move the system by "dq".  The "dq" we have solved for with a 
      ! positive value of dtau corresponds to (q_n - q_n+1)
      ! Therefore, q_n+1 = q_n - dq
      Do i=1, dof
        df(i) = df(i) - dq(i)
      End Do

      ! Recalculate the energy and the gradient
      Call bakersubs_calcgrad(sorbates, simcell, stdparams%mtypes, &
          df, pot, grad, ierr)
      pot = pot*caltoj
      If (ierr /= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " We are in unknown territory"
        ierr = 1
        Return
      End If

      ! Save the path coordinates and the total energy
      pathinfo(nsteps)%df(1:dof) = df(1:dof)
      pathinfo(nsteps)%nrg(nrgnfields) = pot
    End Do
  End Subroutine stdescent_descend

  
  !-------------------------------------------------------------------
  ! Do the actual steepest descent using Fukui's IRC approach
  !-------------------------------------------------------------------
  Subroutine stdescent_DoSim(stdparams, sorbates, simcell)
    Type(STDESCENT_Params), Intent(inout)           :: stdparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(SimCell_Params), Intent(inout)            :: simcell

    Integer        :: nsorbs, i, j, index, dof, statsunit, nsucc, tsrecordno
    Integer        :: recordno, ierr, nsdlpts, nrgnfields, sitesorbtype
    Integer        :: nsteps1, nsteps2, sitetypeTS, sitetype1, sitetype2
    Integer        :: min1index, min2index, begindx, endindx
    Integer, Dimension(MAX_SORBS) :: lastposunits
    Integer, Dimension(10) :: errorlog
    Character(len=strLen)  :: outfilename, outbasefile
    Character(len=2*lstrLen) :: sorborder
    Real(kind=RDbl)        :: totnrg, sdptpot, min1pot, min2pot, minspacing
    Real(kind=RDbl)        :: minminnrg, minActNrg
    Real(kind=RDbl)        :: lastenergy, energytol, min1TSnrg, min2TSnrg
    Real(kind=RDbl), Dimension(stdparams%dof)   :: df, pdf, ndf, grad, eigval
    Real(kind=RDbl), Dimension(stdparams%dof)   :: descentvec
    Real(kind=RDbl), Dimension(stdparams%dof,stdparams%dof):: eigvec, hess
    Type(PATHINFO_Params),Dimension(stdparams%maxsteps):: pathinfo1, pathinfo2
    Logical         :: mapflag, done
    Type(VecType)   :: com1, com2, comTS
    
    !** Get some sundry information
    nsorbs = stdparams%nsorbs
    dof    = stdparams%dof
    outbasefile = stdparams%outbasefile
    lastenergy = 0.0_RDbl
    energytol = stdparams%energytol
    nsdlpts   = stdparams%nsdlpts
    sitesorbtype = stdparams%sitesorbtype
    begindx   = stdparams%sortindx
    endindx   = begindx + 2
    nrgnfields = stdparams%nrgnfields  ! Last field has the total energy
    minspacing = stdparams%minspacing
    minminnrg  = stdparams%minminnrg
    minActNrg  = stdparams%minActNrg
    errorlog = 0
    nsucc    = 0

    !** Generate a string holding the order of sorbates
    sorborder = ""
    Do i=1, stdparams%nsorbs
      sorborder = Trim(sorborder)//" "// &
          Trim(molecules_name(stdparams%mtypes(i)))
    End Do

    !** Open the output files only if are doing the descent
    !** It is possible that we are only doing the histogram calculation
    !** on an existing "cxns" file
    If (stdparams%descent) Then
      outfilename = Trim(outbasefile)//".cxns"
      statsunit = file_getunit(outfilename)
      Open(unit=statsunit, file=outfilename, status='new')
    End If

    !** Do the descent on the unique stationary points
    recordno   = 0
    tsrecordno = 0
    done = .False.
    Do
      !** It is possible we just want to do the histogram calculation
      If (.Not. stdparams%descent) Exit

      !** Read the input file.  The points from which to descend are chosen
      !** based on two criteria: 1) The energy difference with the last 
      !** point and 2) whether the total energy - lowest minimum energy
      !** is greater than a given cutoff specified in the control file.
      ierr = 0
      Do
        tsrecordno = tsrecordno + 1
        If (tsrecordno > nsdlpts) Then
          done = .True.
          Exit
        End If
        df = stdparams%sdlpts(tsrecordno)%df
        totnrg = stdparams%sdlpts(tsrecordno)%nrg(nrgnfields)
        If (Abs(lastenergy - totnrg) > energytol &
            .And. totnrg > minminnrg+minActNrg) Then
          Exit
        End If
      End Do
      
      If (done) Then
        ! We are done
        Exit
      End If
      lastenergy = totnrg
      Write(*,'(a)') '-------------------------------------------------------'
      Write(*,'(a, i5)') 'Processing Record No. : ', tsrecordno

      !** Make sure that the energies match
      Call bakersubs_xToSorb(sorbates, simcell, df, stdparams%mtypes)
      mapflag = bakersubs_getenergy(sorbates, simcell, &
          stdparams%mtypes, sdptpot)
      sdptpot = sdptpot*calToj
      Write(*,'(a,e10.4)') "Transition State Energy Calculated(kJ): ", sdptpot
      If ( (.Not. mapflag) .Or.  &
          Abs((sdptpot - totnrg)/sdptpot) > 1.0e-2_RDbl) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Error in energy calculation"
        Write(0,*) "Energy Read in: ", totnrg
        Write(0,*) "Energy Calculated: ", sdptpot
        Stop
      Endif

      !** Get some preliminary information
      sitetypeTS = stdparams%sdlpts(tsrecordno)%sitetype
      comTS      = df(begindx:endindx)

      !** Get the first step vector
      Call stdescent_firststepvec(stdparams, sorbates, simcell, df, &
          descentvec, ierr)
      If (ierr /= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not take the first step"
        Cycle
      End If

      !** Do the steepest descent in the positive direction
      pdf = df
      Call stdescent_descend(stdparams, sorbates, simcell, pdf, &
          descentvec, min1pot, pathinfo1, nsteps1, min1index, ierr)
      Write(*,*) 'ierr ', ierr
      ! Check if the search was successful
      If (ierr /= 0) Then
        errorlog(ierr) = errorlog(ierr) + 1
        ! We still want to keep the path even though it did not find a 
        ! matching minimum.  We are assuming that for ierr=2 we are close
        ! to a minimum and that is why the energy is increasing.
        If (ierr /= 2 .And. ierr /= 5) Then
          Cycle
        End If
      End If

      min1TSnrg = sdptpot - min1pot
      ! Get the site type and reduce the no. of points by deleting points
      ! that are closer than "minspacing" 
      com1 = pathinfo1(nsteps1)%df(begindx:endindx)
      sitetype1 = smap_getPosSiteType(stdparams%smap, com1)
      nsteps1 = &
          stdescent_reducepts(pathinfo1(1:nsteps1), nsteps1, dof, minspacing)
      
      !** Do the steepest descent in the negative direction
      descentvec = -1.0_RDbl*descentvec
      ndf = df
      Call stdescent_descend(stdparams, sorbates, simcell, ndf, &
          descentvec, min2pot, pathinfo2, nsteps2, min2index, ierr)
      Write(*,*) 'ierr ', ierr
      ! Check if the search was successful
      If (ierr /= 0) Then
        errorlog(ierr) = errorlog(ierr) + 1
        ! We still want to keep the path even though it did not find a 
        ! matching minimum.  We are assuming that for ierr=2 we are close
        ! to a minimum and that is why the energy is increasing.
        If (ierr /= 2 .And. ierr /= 5) Then
          Cycle
        End If
      End If
      min2TSnrg = sdptpot - min2pot
      ! Get the site type and reduce the no. of points by deleting points
      ! that are too close together
      com2 = pathinfo2(nsteps2)%df(begindx:endindx)
      sitetype2 = smap_getPosSiteType(stdparams%smap, com2)
      nsteps2 = &
          stdescent_reducepts(pathinfo2(1:nsteps2), nsteps2, dof, minspacing)
      
      !** Write out stuff to the output file
      recordno = recordno + 1
      Write(statsunit,'(a)') dashedline    
      Write(statsunit,'(a,i7)') "_RECORD_NO.        :", recordno
      Write(statsunit,'(2a)')   "_SORB_ORDER        :", Trim(sorborder)
      Write(statsunit,'(a,i7)') "_MOVING_SORBINDX   :", begindx 
      Write(statsunit,'(a,i7)') "_TS_RECORD_NO.     :", tsrecordno
      Write(statsunit,'(a)', Advance='No')    "_TS_DF             :"
      Call dispvec(df, "f8.3", statsunit)
      Write(statsunit, *)
      Write(statsunit,'(a,i7)') "_MIN1_RECORD_NO.   :", min1index
      Write(statsunit,'(a)', Advance='No')    "_MIN1_DF           :"
      Call dispvec(pathinfo1(nsteps1)%df(1:dof), "f8.3", statsunit)
      Write(statsunit,*)
      Write(statsunit,'(a,i7)') "_MIN2_RECORD_NO.   :", min2index
      Write(statsunit,'(a)', Advance='No')    "_MIN2_DF           :"
      Call dispvec(pathinfo2(nsteps2)%df(1:dof), "f8.3", statsunit)
      Write(statsunit,*)
      Write(statsunit,'(a,f10.4)') "_DIST_MOVED1_TS :", mag(comTS-com1)
      Write(statsunit,'(a,f10.4)') "_DIST_MOVED2_TS :", mag(comTS-com2)
      Write(statsunit,'(a,i3)') "_SITE_TYPES(1)     :", sitetype1
      Write(statsunit,'(a,i3)') "_SITE_TYPES(TS)    :", sitetypeTS
      Write(statsunit,'(a,i3)') "_SITE_TYPES(2)     :", sitetype2
      Write(statsunit,'(a,f10.4)') "_1_TS_ENERGY       :", min1TSnrg
      Write(statsunit,'(a,f10.4)') "_2_TS_ENERGY       :", min2TSnrg
      Write(statsunit,'(a,i7)') "_NPTS1             :", nsteps1
      Write(statsunit,'(a,i7)') "_NPTS2             :", nsteps2
      ! Write the coordinates of min1
      Write(statsunit,'(a)', Advance='No') "_MIN1:"
      Do i=nsteps1, 1, -1
        Do j=1, dof
          Write(statsunit,'(f8.3)', Advance='No') pathinfo1(i)%df(j)
        End Do
        Write(statsunit,*)
      End Do
      ! Write the coordinates of the Transition State
      Write(statsunit,'(a)', Advance='No') "_TS_:"
      Do j=1, dof
        Write(statsunit, '(f8.3)', Advance='No') df(j)
      End Do
      Write(statsunit,*)
      ! Write the coordinates of min2
      Do i=1, nsteps2
        If (i == nsteps2) Then
          Write(statsunit,'(a)', Advance='No') "_MIN2:"
        End If
        Do j=1, dof
          Write(statsunit,'(f8.3)', Advance='No') pathinfo2(i)%df(j)
        End Do
        Write(statsunit,*)
      End Do
    End Do

    !** Print out the error log
    Do i=1, 5
      Write(*,'(a,i3,a,i5)') "No. of errors of type : ", i, " = ", errorlog(i)
    End Do
    Close(statsunit)

    !** If the histogram flag is set then make it
    If (stdparams%histflag) Then
      Call stdescent_makehistogram(stdparams)
    End If
  End Subroutine stdescent_DoSim

  !----------------------------------------------------------
  ! This routine makes a histogram of activation energies 
  !----------------------------------------------------------
  Subroutine stdescent_makehistogram(stdparams)
    Type(STDESCENT_Params), Intent(inout)  :: stdparams

    Type HEADER_INFO
      Integer         :: npts1, npts2
      Real(kind=RDbl) :: actnrg1, actnrg2
    End Type HEADER_INFO

    Character(len=strLen) :: outbasefile, cxnsfile, junk1, junk2, histfile
    Character(len=strLen) :: line
    Integer     :: cxnsunit, i, error, histunit, lineno, nrecords, temp
    Integer     :: npts1, npts2, nfields, lastnpts1, lastnpts2
    Integer     :: nhistelems
    Character(len=strLen), Dimension(10) :: fields
    Real(kind=RDbl) :: lowval, highval, actnrg1, actnrg2 
    Real(kind=RDbl) :: lastactnrg1, lastactnrg2
    Type(Histogram_Params) :: histobj
    Type(Header_Info), Dimension(MAX_PTS)   :: header, tempheader
    Type(SortArray_Type), Dimension(MAX_PTS):: sortarray
    
    outbasefile = stdparams%outbasefile
    cxnsfile    = Trim(outbasefile)//".cxns"
    cxnsunit = file_getunit(cxnsfile)
    Open(unit=cxnsunit, file=cxnsfile, IOSTAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Could not open file ", Trim(cxnsfile)
      Stop
    End If
    Rewind(cxnsunit)
    
    !** Initialize the histogram
    ! The lowest values is 0
    lowval = 0.0_RDbl
    ! The highest possible
    highval = (stdparams%maxsdptnrg - stdparams%minminnrg)
    Call histogram_init(histobj, lowval, highval, stdparams%binsize)
    
    !** Read the records in the ".cxns" file and store them in an array
    nrecords  = 0
    Do
      ! Find the separator line "------------------"
      lineno = filesrchstr(cxnsunit, '1_TS_ENERGY', line)
      If (lineno == 0) Then
        ! We have reached the end
        Exit
      End If
      nfields = split(line, fields, ":")
      actnrg1 = toreal(fields(2))
      Read(cxnsunit, *) junk1, junk2, actnrg2
      Read(cxnsunit, *) junk1, junk2, npts1
      Read(cxnsunit, *) junk1, junk2, npts2
      
      ! Keep track of the no. of records read
      nrecords = nrecords + 1
      If (Mod(nrecords, 100) == 0) Then
        Write(*,*) "No. of records read ", nrecords
      End If
      
      ! Store the values in the array
      If (npts1 < npts2) Then
        header(nrecords)%actnrg1 = actnrg1
        header(nrecords)%npts1   = npts1
        header(nrecords)%actnrg2 = actnrg2
        header(nrecords)%npts2   = npts2
      Else
        header(nrecords)%actnrg1 = actnrg2
        header(nrecords)%npts1   = npts2
        header(nrecords)%actnrg2 = actnrg1
        header(nrecords)%npts2   = npts1
      End If
    End Do

    !** Sort the header array on the basis of the first activation energy
    !** and then the second activation energy
    Call sort_init(header(1:nrecords)%actnrg1, header(1:nrecords)%actnrg2, &
        sortarray, nrecords)
    Call sort_heapsort(sortarray, nrecords)
    tempheader = header
    Do i=1, nrecords
      header(i) = tempheader(sortarray(i)%index)
    End Do

    !** Now create the histogram after sifting out the duplicate entries
    nhistelems = 0
    lastnpts1 = 0
    lastnpts2 = 0
    lastactnrg1 = 0.0_RDbl
    lastactnrg2 = 0.0_RDbl
    Do i=1, nrecords
      actnrg1 = header(i)%actnrg1
      actnrg2 = header(i)%actnrg2
      npts1   = header(i)%npts1
      npts2   = header(i)%npts2
!!$      Write(*,'(a,2f9.4, 2i4)', Advance='No') &
!!$          'e1, e2, n1, n2 ', actnrg1,actnrg2,npts1,npts2
      
      ! Make sure that we actually have a unique path.  This is done
      ! by comparing npts1 and npts2 with lastnpts1 and lastnpts2.  If the
      ! the nos. match within +/- 2 pts
      If (Abs(lastnpts1-npts1)<2 .And. &
          Abs(lastnpts2-npts2)<2 .And. &
          Abs(lastactnrg1-actnrg1)<1.0e-2_RDbl .And. &
          Abs(lastactnrg2-actnrg2)<1.0e-2_RDbl) Then
!!$        Write(*,*) 'skip '
        Cycle
      Else
        lastactnrg1 = actnrg1
        lastnpts1 = npts1
        lastactnrg2 = actnrg2
        lastnpts2 = npts2
!!$        Write(*,*) 
      End If

      ! Set the lastpts      
      ! We could have round-off error for very small paths (e.g. 1 pt paths)
      ! making the activation energy negative.  This would crash the histogram
      If (actnrg1 > 0.0_RDbl) Then
        Call histogram_update(histobj, actnrg1)
      End If
      If (actnrg2 > 0.0_RDbl) Then
        nhistelems = nhistelems + 1
        Call histogram_update(histobj, actnrg2)
      End If
    End Do
    Close(cxnsunit)
    
    !** Write the no. of elements
    Write(*,*) 'Total No. of Records Read :', nrecords
    Write(*,*) 'No. of histogram elements :', nhistelems
    
    !** Write the histogram to a file
    histfile = Trim(outbasefile)//".hist"
    histunit = file_getunit(histfile)
    Open(unit=histunit, file=histfile)
    Call histogram_display(histobj, histunit)
    Close(histunit)
  End Subroutine stdescent_makehistogram

  !---------------------------------------------------------
  ! This routine writes the current energy to "unitno"
  !---------------------------------------------------------
  Subroutine stdescent_dumpstats(sorbates, simcell, mtypes, unitno, pot)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(in)               :: simcell
    Integer, Intent(in), Dimension(:) :: mtypes
    Integer, Intent(in)  :: unitno
    Real(kind=RDbl), Intent(in) :: pot

    Integer               :: i, j, nsorb, mtype, nsorbs
    Integer, Save         :: iterno
    Character(len=strLen) :: avgposstr, molecname
    Logical, Save         :: firstpass = .True.

    !** Number of stdescent sorb types
    nsorbs = Size(mtypes,1)
    
    !** Check if this the first call for this plane. If so dump
    !** the header
    If (firstpass) Then
      firstpass = .False.
      Write(unitno, '(a,4x)', Advance='No') "#"
      Do j=1, nsorbs
        mtype = mtypes(j)
        molecname = molecules_name(mtype)
        Write(unitno,'(a27,23x)', Advance='No') Trim(molecname)
      End Do
      Write(unitno,'(a13)') "Total(kJ)"
      iterno = 0
    End If
    
    !** Dump the generalized coordinates of the stationary point
    iterno = iterno + 1
    Write(unitno, '(i5)', Advance='No') iterno
    Do i=1, nsorbs
      mtype = mtypes(i)
      Call gcmodels_pbc(sorbates(mtype)%gcoords(1), simcell)
      avgposstr = gcmodels_display(sorbates(mtype)%gcoords(1), "f8.3")
      Write(unitno, '(a50)', Advance='No') Trim(avgposstr)
    End Do
    Write(unitno, '(f13.7)') pot
  End Subroutine stdescent_dumpstats

  !----------------------------------------------------------------
  ! Display the minimum energy path calculation parameters 
  !----------------------------------------------------------------
  Subroutine stdescent_display(stdparams, nspc, optunit)
    Type(STDESCENT_Params), Intent(in) :: stdparams
    Integer, Intent(in)     :: nspc
    Integer, Intent(in), Optional :: optunit

    Character(len=nspc)  :: spc
    Integer              :: unitno, i
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "Steepest Descent Parameters"
    Write(unitno, '(2a)') spc, dashedline
    If (stdparams%descent) Then
      Write(unitno, '(a,2x,a)') spc, "Doing Steepest Descent     : "
    Else
      Write(unitno, '(a,2x,a)') spc, "Reading from Existing File : "
    End If
    If (stdparams%histflag) Then
      Write(unitno, '(a,2x,a)') spc, "Making Histogram           : "
      Write(unitno, '(a,4x,a,f8.3)') spc, "Histogram Bin Size         : ", &
          stdparams%binsize
    Else
      Write(unitno, '(a,2x,a)') spc, "Histogram Not Generated    : "
    End If

    Write(unitno, '(a,2x,2a)') spc, "Input File        : ", &
        stdparams%inbasefile
    Write(unitno, '(a,2x,2a)') spc, "Output File       : ", &
        stdparams%outbasefile
    Write(unitno, '(a,2x,a, i7)') spc, "No. of Minima     : ", &
        stdparams%nminpts
    Write(unitno, '(a,2x,a, 2f13.7)') spc, "Range of Minima   : ", &
        stdparams%minminnrg, stdparams%maxminnrg
    Write(unitno, '(a,2x,a, i7)') spc, "No. of Saddle Pts  : ", &
        stdparams%nsdlpts
    Write(unitno, '(a,2x,a, 2f13.7)') spc, "Range of Saddle-Pts: ", &
        stdparams%minsdptnrg, stdparams%maxsdptnrg
    Call smap_display(stdparams%smap, nspc+2, unitno)
  End Subroutine stdescent_display

End Module stdescent
