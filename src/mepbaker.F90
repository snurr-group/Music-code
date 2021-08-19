Module mepbaker

  Use baker, Only: baker_isSaddleSearch, baker_display, baker_baker, &
      baker_init, BAKER_Params
  Use bakersubs, Only: bakersubs_getenergy, bakersubs_sorbtox, &
      bakersubs_calcgrad, bakersubs_xtosorb
  Use config, Only: AtMolCoords, config_getnatoms, config_getsourcefile, &
      config_getsourcetype, config_setnmoles, config_isfixed, config_getnmoles
  Use defaults, Only: RDbl, twopi, pi, lstrLen, strLen, dashedline, &
      MAX_ATOMS, MAX_DOF, MAX_SORBS
  Use file, Only: file_getunit
  Use gcmc, Only: GCMC_Params, gcmc_init, gcmc_dosim, gcmc_displaysimparams
  Use gcmodels, Only: gcmodels_getcom, gcmodels_display, gcmodels_toxyz, &
      gcmodels_readrestartfile, gcmodels_pbc, gcmodels_setgencoords
  Use molecules, Only: molecules_getnatoms, molecules_getcomdefn, &
      molecules_name, molecules_getdof, molecules_getgcmodeltype, &
      molecules_gettype, molecules_getcoul, molecules_getnoncoul, &
      molecules_getnsorbs, molecules_displaynrg, molecules_updateenergySS
  Use nvtmc, Only: NVTMC_Params, nvtmc_display, nvtmc_dosim, nvtmc_init
  Use random, Only: rranf, random_getiseed
  Use simcell, Only: Simcell_Params, simcell_pbc
  Use utils, Only: stripcmnt, split, filesrchstr, isfileopen
  Use vector, Only: VecType, vector_display
  
  Implicit None
  Save

  Private
  Public :: default_mepbaker_tag, MEPBAKER_Params, Initconfig_Params, &
      Sphereins_Params, mepbaker_init, mepbaker_dosim, mepbaker_display

  Character(len=strLen)    :: default_mepbaker_tag = &
      "Baker Minimum Energy Path"

  Character(len=strLen), Dimension(0:2), Parameter :: mepbaker_inittype = &
      (/"NVTMC              ", &
        "Random Insertion   ", &
        "Spherical Insertion"/)

  Integer, Parameter :: MAX_INSCENTERS = 10
  Private MAX_INSCENTERS

  Type SPHEREINS_Params
    Type(VecType)    :: center
    Real(kind=RDbl)  :: maxradius
    Real(kind=RDbl)  :: minradius
  End Type SPHEREINS_Params

  !** This data structure keeps all the parameters/structures we need
  !** to generate the initial configurations for the minimization
  Type INITCONFIG_Params
    Character(len=strlen)   :: ctrltag
    Integer                 :: taglineno
    Integer                 :: pertflag
    ! Generalized Coord Params for insertion
    Type(GCMC_Params)       :: insparams
    Real(kind=RDbl)         :: tk
    Integer                 :: simno
    ! Params for doing the NVT Monte Carlo
    Type(NVTMC_Params)      :: nvtmcparams
    ! Params for doing the Spherical Insertion
    Integer                 :: nsphinsparams
    Type(SPHEREINS_Params), Dimension(MAX_INSCENTERS)  :: sphinsparams
  End Type INITCONFIG_Params
  
  !** Note that the arrays in this data type are indexed by the Baker
  !** molec number, NOT molecule type. If you want the molecule TYPE, look 
  !** it up in mtypes
  Type MEPBAKER_Params
    Integer                :: niterations
    Integer                :: lastiterunit
    Integer                :: nsucc   !No. of searches that succeeded
    Integer                :: nsorbs
    Integer                :: startiter
    Character(len=strLen)  :: bakertag
    Type(Baker_Params)     :: bakerparams
    Character(len=strLen)  :: basefile, debugfile
    Integer                :: dof
    Real(kind=RDbl)        :: initucut
    Type(INITCONFIG_Params), Dimension(:), Pointer :: initconf
    Integer, Dimension(:), Pointer :: mtypes
  End Type MEPBAKER_Params

  Type MININFO_Params
    Real(kind=RDbl), Dimension(MAX_DOF) :: df
    Real(kind=RDbl)                     :: pot
    Integer                             :: sitetype
  End Type MININFO_Params

  Interface display
    Module Procedure mepbaker_display
  End Interface

Contains
  !------------------------------------------------------------------
  ! Initialize the parameters for getting the Minimum Energy Paths
  ! using the Baker's algorithm.  Get the parameters from the control
  ! file "ctrl_filename" and store them in the "mepbakersparams" object
  !------------------------------------------------------------------
  Subroutine mepbaker_init(mepparams, sorbates, simcell, ctrl_filename, &
      opt_mepbakertag)
    Type(MEPBAKER_Params), Intent(out)  :: mepparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(SimCell_Params), Intent(inout) :: simcell
    Character(*), Intent(in)            :: ctrl_filename
    Character(*), Intent(in),Optional:: opt_mepbakertag

    Real(kind=RDbl)         :: pot, tk
    Character(len=strLen)   :: tag, sorbname, sourcetype, filename, gcmodeltype
    Character(len=strLen)   :: debugfile
    Character(len=lstrLen)  :: line
    Integer                 :: lineno, unitno, i, j, pertlineno, natoms, mtype
    Integer                 :: bakerlineno, dof, error
    Integer                 :: nsorbs, ins_sorbs, lastiterunit
    Logical                 :: mapflag
    Type(VecType), Dimension(MAX_ATOMS) :: comdefn

    If (Present(opt_mepbakertag)) Then
      tag = opt_mepbakertag
    Else
      tag = default_mepbaker_tag
    End If
    
    !** Open the ctrl_file if it is not opened
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
      Open(file=ctrl_filename, unit=unitno)
    Endif
    Rewind(unitno)
    
    !** Find the MEPBAKER section
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif
    
    !** Read the rest of the stuff
    Read(unitno,*) mepparams%niterations
    Read(unitno,*) mepparams%startiter
    Read(unitno,'(a)') mepparams%bakertag
    mepparams%bakertag = stripcmnt(mepparams%bakertag)
    Read(unitno,*) mepparams%basefile
    Read(unitno,*) mepparams%initucut
    Read(unitno,*) nsorbs
    mepparams%nsorbs = nsorbs
    mepparams%nsucc = 0

    !** Allocate space for the mtypes and initconf array
    Allocate(mepparams%mtypes(mepparams%nsorbs),stat=error)
    If (error /=0) Then
      Write(0,'(2a,i6,a,i4)') __FILE__,":",__LINE__, &
          " Could not allocate mepparams%mtypes of size ",mepparams%nsorbs
      Stop
    End If

    Allocate(mepparams%initconf(mepparams%nsorbs),stat=error)
    If (error /=0) Then
      Write(0,'(2a,i6,a,i4)') __FILE__,":",__LINE__, &
          " Could not allocate mepparams%initconf of size ",mepparams%nsorbs
      Stop
    End If

    !** Initialize the initial configuration generation parameters
    Call mepbaker_initgenconfig(mepparams, sorbates, simcell, lineno, &
        ctrl_filename, unitno)

    !** Initialize the Baker's algorithm parameters
    bakerlineno = lineno + 2
    Rewind(unitno)
    Call baker_init(mepparams%bakerparams,ctrl_filename,mepparams%bakertag, &
        bakerlineno)
    mepparams%bakerparams%dof = mepparams%dof

    !** Set up a debug file
    debugfile = Trim(mepparams%basefile)//".dbx"  ! In case the run crashes
    lastiterunit = file_getunit(debugfile)        ! we can look at this file
    mepparams%lastiterunit = lastiterunit
    mepparams%debugfile    = debugfile
    mepparams%bakerparams%lastiterunit = lastiterunit

    !** Initialize the starting configuration.
    !** It may be read in from a RESTARTFILE in which case we are all set.
    !** If the initial configuration has to be taken from an MEPRESTARTFILE
    !** we have to do that here
    Do i=1, nsorbs
      mtype = mepparams%mtypes(i)
      sourcetype = config_getsourcetype(sorbates, mtype)
      filename   = config_getsourcefile(sorbates, mtype)
      natoms     = config_getnatoms(sorbates, mtype)
      If (sourcetype == 'MEPRESTARTFILE') Then
        Call mepbaker_readmeprestartfile(sorbates, mtype, filename)

        Call molecules_getcomdefn(mtype,comdefn)

        ! Generate the xyz principal coordinates
        Call gcmodels_toxyz(sorbates(mtype)%gcoords(1), &
            sorbates(mtype)%coords(1:natoms,1)%rp, comdefn)
        
        ! Generate the other coordinates
        Call simcell_pbc(simcell, sorbates(mtype)%coords(1:natoms,1)%rp, &
            sorbates(mtype)%coords(1:natoms, 1)%r, &
            sorbates(mtype)%coords(1:natoms, 1)%cr)
      End If
    End Do

    ! Calculate the energy of the initial configuration
    mapflag =  bakersubs_getenergy(sorbates, simcell, mepparams%mtypes, pot)
    If ( .Not. mapflag) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Initial configuration in untabulated region"
      Stop
    End If

    Call molecules_displaynrg(2,6)
  End Subroutine mepbaker_init

  !------------------------------------------------------------------------
  ! Initialize the parameters that are used to generate the initial
  ! configuration.  "lineno" is the line number in the control file where
  ! the entire init section begins. "unitno" is the LUN of the control file
  !------------------------------------------------------------------------
  Subroutine mepbaker_initgenconfig(mepparams,sorbates,simcell,sectionlineno, &
      ctrl_filename, unitno)
    Type(MEPBAKER_Params), Intent(inout)  :: mepparams
    Type(Simcell_Params), Intent(in)   :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: sorbates
    Integer, Intent(in)    :: sectionlineno
    Character(len=strLen), Intent(in) :: ctrl_filename
    Integer, Intent(in)    :: unitno

    Integer                 :: i, j, dof, currlineno, taglineno
    Integer                 :: nsorbs, sorbtype
    Real(kind=RDbl)         :: tk, x1, y1, z1, maxradius, minradius
    Character(len=strLen)   :: tag, sorbname, gcmodeltype
    Character(len=lstrLen)  :: line

    !** Get the no. of sorbates
    nsorbs = mepparams%nsorbs

    !** Get the initialization parameters for the different sorbates
    dof = 0
    currlineno = sectionlineno + 10
    Do i=1, nsorbs
      Read(unitno,*)
      Read(unitno,*) sorbname
      mepparams%mtypes(i) = molecules_gettype(trim(sorbname))
      If (mepparams%mtypes(i) == 0) Then
        Write(0,'(2a,i5,2a)') __FILE__,":",__LINE__, &
            " Unable to recognize species type ",trim(sorbname)
        Stop
      End If
      If (config_isfixed(sorbates(mepparams%mtypes(i)))) Then
        Write(0,'(2a,i5,3a)') __FILE__,":",__LINE__, &
            " Unable to use molecule ",Trim(sorbname), &
            " since it is a FIXED type"
        Stop
      End If
      Read(unitno,*) mepparams%initconf(i)%pertflag
      Read(unitno,'(a)') tag
      mepparams%initconf(i)%ctrltag = Trim(stripcmnt(tag))
      mepparams%initconf(i)%taglineno = currlineno
      sorbtype = molecules_gettype(sorbname)
      If (sorbtype /= i) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
            " The sorbates to be minimized should be listed first in", &
            " the control file and they should have the same order in", &
            " the initialization section"
        Stop
      End If
      currlineno = currlineno + 4
      dof = dof + molecules_getdof(mepparams%mtypes(i))
    End Do
    mepparams%dof = dof

    !** Initialize the technique that will be used to generate the initial
    !** configuration.  We have the option of using NVTMC to perturb an
    !** initial configuration or use GCMC style insertion to generate a
    !** random configuration
    Do i=1, nsorbs
      Rewind(unitno)
      tag = mepparams%initconf(i)%ctrltag
      taglineno = mepparams%initconf(i)%taglineno
      If (mepparams%initconf(i)%pertflag == 0) Then
        ! Use Monte Carlo
        Call nvtmc_init(mepparams%initconf(i)%nvtmcparams, sorbates, simcell, &
            ctrl_filename, tag, taglineno)

      Else If (mepparams%initconf(i)%pertflag == 1) Then
        ! Use Random Insertion
        ! Find the tag that marks the section for the bias tags
        currlineno = filesrchstr(unitno, tag, line)
        If (currlineno == 0 .Or. currlineno == taglineno) Then
          Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
              " Could not find the section '", Trim(tag), "'"
          Write(0,'(1x,2a)') &
              "Either it is not defined in the control file or", &
              " it is defined AFTER the section referring to it"
          Stop
        End If

        !** Find the string
        Read(unitno,*) tk
        gcmodeltype = molecules_getgcmodeltype(mepparams%mtypes(i))
        Call gcmc_init(mepparams%initconf(i)%insparams,sorbates, &
            simcell, ctrl_filename, tag)
        mepparams%initconf(i)%simno = 1

      Else If (mepparams%initconf(i)%pertflag == 2) Then
        ! Use a spherical random insertion
        currlineno = filesrchstr(unitno, tag, line)
        If (currlineno == 0 .Or. currlineno == taglineno) Then
          Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
              " Could not find the section '", Trim(tag), "'"
          Write(0,'(1x,2a)') &
              "Either it is not defined in the control file or", &
              " it is defined AFTER the section referring to it"
          Stop
        End If
        Read(unitno,*) mepparams%initconf(i)%nsphinsparams
        Write(*,*) mepparams%initconf(i)%nsphinsparams
        Do j=1,  mepparams%initconf(i)%nsphinsparams
          Read(unitno,*) x1, y1, z1
          Read(unitno,*) minradius, maxradius
          mepparams%initconf(i)%sphinsparams(j)%center = (/x1, y1, z1/)
          mepparams%initconf(i)%sphinsparams(j)%maxradius = maxradius
          mepparams%initconf(i)%sphinsparams(j)%minradius = minradius
        End Do
      End If
    End Do
  End Subroutine mepbaker_initgenconfig

  !-----------------------------------------------------------------------
  ! This routine reads the initial generalized coordinates from the
  ! given by "filename".  This format is refered to as MEPRESTARTFILE
  ! format. The sorbate structure is already allocated in "config.f90"
  !-----------------------------------------------------------------------
  Subroutine mepbaker_readmeprestartfile(sorbates, sorbno, filename)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
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
    Rewind(unitno)

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
    Call gcmodels_readrestartfile(sorbates(sorbno)%gcoords(1),sorbno,unitno)
  End Subroutine mepbaker_readmeprestartfile

  !------------------------------------------------------------------------
  ! This subroutine generates the intial configuration using either
  ! a NVTMC simulation to perturb an initial configuration or a GCMC type
  ! insertion to generate an initial configuration in the case of a single
  ! sorbate
  !------------------------------------------------------------------------
  Subroutine mepbaker_genconfig(mepparams, sorbates, simcell)
    Type(MEPBAKER_Params), Intent(inout)  :: mepparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(inout)   :: simcell

    Integer     :: i, nsorbs, natoms, centerno, ncenters, moveno, j, mtype
    Real(kind=RDbl) :: biasfactor, maxradius, minradius, radius, pot
    Real(kind=RDbl) :: theta, phi, x, y, z
    Real(kind=RDbl) :: costheta, sintheta, cosphi, sinphi
    Real(kind=RDbl), Dimension(6) :: gencoordarray
    Type(VecType)   :: com, com1, com2, center, origin_com
    Logical         :: mapflag
    Type(VecType), Dimension(MAX_ATOMS) :: comdefn

    !** Number of sorbates involved in the mepbaker routine
    nsorbs = mepparams%nsorbs

    !** Do the perturbation kind of moves first since this is most time
    !** consuming
    Do i=1, nsorbs
      mtype = mepparams%mtypes(i)
      If (mepparams%initconf(i)%pertflag == 0) Then
        com1 = gcmodels_getcom(sorbates(mtype)%gcoords(1))
        Call nvtmc_dosim(mepparams%initconf(i)%nvtmcparams, &
            sorbates, simcell, 1)
        com2 = gcmodels_getcom(sorbates(mtype)%gcoords(1))
!!$        Write(*,*) 'moved sorbate ', i, ' by : ', mag(com1-com2)
      End If
    End Do

    !** Do the insertion kind of moves first and check for the energy
    !** of the system as we insert.
    i = 1
    mtype = mepparams%mtypes(i)
    Do 
      !** Don't mess with fixed position sorbates
      If (config_isfixed(sorbates(mtype))) Then
        i = i+1
      End If
      If (i > nsorbs) Exit
      If (mepparams%initconf(i)%pertflag == 0) Then
        ! This is already done
        i = i+1
        Cycle
      Else If (mepparams%initconf(i)%pertflag == 1) Then
        !** We need to insert one molecule in here. Reset the number of
        !** sorbates to zero for this sorb type.
        Call config_setnmoles(sorbates,mtype,0)
        !** Do an insertion
        Do
          Call gcmc_dosim(mepparams%initconf(i)%insparams, sorbates, simcell, &
              mepparams%initconf(i)%simno)
          If (config_getnmoles(sorbates(mtype)) == 1) Exit
        End Do

      Else If (mepparams%initconf(i)%pertflag == 2) Then        
        !** Do an insertion inside a sphere
        ! Pick a center
        ncenters = mepparams%initconf(i)%nsphinsparams
        centerno = Int(rranf()*ncenters) + 1

        ! Pick a radius 
        maxradius = mepparams%initconf(i)%sphinsparams(centerno)%maxradius
        minradius = mepparams%initconf(i)%sphinsparams(centerno)%minradius
        center     = mepparams%initconf(i)%sphinsparams(centerno)%center
        radius = rranf()*(maxradius-minradius) + minradius
        ! Pick one angle from a cosine distribution
        costheta = (rranf()*2.0_Rdbl) - 1.0_RDbl
        If (costheta > 1.0_RDbl) Then
          theta = 0.0_RDbl;
        Else If (costheta < -1.0_RDbl) Then
          theta = pi;
        Else
          theta = Acos(costheta);
        End If
        costheta = Cos(theta);
        sintheta = Sin(theta);
        ! Pick the second angle from 0-2*pi
        phi = rranf()*twopi
        sinphi = Sin(phi)
        cosphi = Cos(phi)
        ! Generate the x,y,z coordinates 
        x = radius*sintheta*cosphi
        y = radius*sintheta*sinphi
        z = radius*costheta
        origin_com = (/x, y, z/)
        ! Displace the system to the center
        com = origin_com + center

        ! Set the coordinates in the sorbates structure
        gencoordarray(1:3) = com
        gencoordarray(4:6) = 0.0_RDbl
        Call gcmodels_setgencoords(sorbates(mtype)%gcoords(1), gencoordarray)
        natoms = config_getnatoms(sorbates, mtype)
        Call molecules_getcomdefn(mtype,comdefn)
        Call gcmodels_toxyz(sorbates(mtype)%gcoords(1), &
            sorbates(mtype)%coords(1:natoms, 1)%rp, comdefn)
        
        ! Set the total no. of moles and update the other arrays in sorbates
        Call config_setnmoles(sorbates, mtype, 1)
        Call simcell_pbc(simcell, sorbates(mtype)%coords(1:natoms,1)%rp, &
            sorbates(mtype)%coords(1:natoms, 1)%r, &
            sorbates(mtype)%coords(1:natoms, 1)%cr)
      End If
      mapflag = bakersubs_getenergy(sorbates, simcell, mepparams%mtypes, pot)
      If (mapflag) Then
        If (pot <=  mepparams%initucut) i=i+1
      End If
    End Do
  End Subroutine mepbaker_genconfig

  
  !----------------------------------------------------
  ! Do the actual search using Baker's algorithm.
  !----------------------------------------------------
  Subroutine mepbaker_DoSim(mepparams, sorbates, simcell)
    Type(MEPBAKER_Params), Intent(inout)           :: mepparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(SimCell_Params), Intent(inout)            :: simcell

    Integer        :: nsorbs, sorb1, sorb2, i, j, ierr, iterno, nrgunit
    Integer        :: nmoles, inssorbtype, dof, minno, junk, statsunit, nsucc
    Integer        :: lastiterunit, error, startiter, nstep, mtype
    Integer, Dimension(MAX_SORBS) :: lastposunits
    Integer, Dimension(10)        :: errorlog
    Character(len=strLen)         :: outfilename, basefile, debugfile
    Real(kind=RDbl), Dimension(mepparams%dof)   :: savex, x, grad, eigval
    Real(kind=RDbl), Dimension(mepparams%dof,mepparams%dof):: eigvec, hess
    Real(kind=RDbl) :: pot, totnrg
    Logical         :: mapflag
    Character(len=lstrlen)   :: ext, avgposstr

    !** Get some sundry information
    nsorbs = mepparams%nsorbs
    dof    = mepparams%dof
    errorlog = 0
    minno    = 0
    nsucc    = 0
    basefile = mepparams%basefile
    debugfile = mepparams%debugfile
    lastiterunit = mepparams%lastiterunit
    startiter    = mepparams%startiter

    !** Open the output files
    If (baker_isSaddleSearch(mepparams%bakerparams)) Then
      ext = ".sdpt"
    Else
      ext = ".min"
    End If
    Call mepbaker_openoutfiles(mepparams, ext, lastposunits, nrgunit)

    !** Save the old configuration.  We want to perturb this 
    !** configuration every time we iterate. At this point we don't
    !** necessarily have the xyz coordinates in sorbates but that does
    !** not matter since "savex" only has the generalized coordinates
    Call bakersubs_sorbtoX(sorbates, savex, mepparams%mtypes)
    Write(0,'(a)', Advance='No') 'initial x: '
    Do i=1, dof
      Write(0,'(f6.2)', Advance='No') savex(i)
    End Do
    Write(0,*)

    !** Do a bunch of iterations
    Do iterno = startiter, startiter+mepparams%niterations-1
      Open(unit=lastiterunit, file=debugfile, status='unknown', IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
            " Could not open debug file ", Trim(debugfile)
      End If
      Rewind(lastiterunit)
      Write(lastiterunit,*) 'iteration no. : ', iterno
      Write(lastiterunit,*) 'iseed : ', random_getiseed()
      
      !** Restore the configuration to the initial configuration.  This
      !** will prevent the "blocking molecules" from wandering off
      Call bakersubs_xToSorb(sorbates, simcell, savex, mepparams%mtypes)

      !** Generate/Perturb the initial configuration depending on
      !** the value of mepparams%pertflag
      Do
        Call mepbaker_genconfig(mepparams, sorbates, simcell)
        mapflag = bakersubs_getenergy(sorbates, simcell, mepparams%mtypes, pot)
        If (mapflag) Then
          If (pot <= mepparams%initucut) Exit
        End If
      End Do
      
      !** Change the coordinate representation to "x". Do the baker's search
      !** and then go back to the sorbates representation from the new "x"
      Call bakersubs_sorbtoX(sorbates, x, mepparams%mtypes)
      Write(lastiterunit,*) 'x             : ', x(1:dof)
      Call baker_baker(mepparams%bakerparams, sorbates, simcell, &
          mepparams%mtypes, x, pot, grad, eigval, eigvec, hess, nstep, ierr)
      Call bakersubs_xToSorb(sorbates, simcell, x, mepparams%mtypes)


      !** Close the debug file
      Close(unit=lastiterunit, IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Error closing debugfile"
      End If

      !** Print out the iteration number for my peace of mind
      If (Mod(iterno, 10) == 0) Then
        Write(*,'(a, i5)') 'Iteration No. : ', iterno
      End If

      !** Check if the search was successful
      If (ierr /= 0) Then
        errorlog(ierr) = errorlog(ierr) + 1
        Cycle
      Else
        Write(*,'(a, i5)') ' No. of steps:', nstep
      End If
      
      !** Update the success counter
      nsucc = nsucc + 1

      !** Get the final position of the molecules. Write the
      !** final position to the "lastpos" file and the energies
      Call mepbaker_dumpenergy(iterno, nrgunit, totnrg, mepparams%mtypes)
      Do i=1, nsorbs
        mtype = mepparams%mtypes(i)
        Call gcmodels_pbc(sorbates(mtype)%gcoords(1), simcell)
        avgposstr = gcmodels_display(sorbates(mtype)%gcoords(1), "f9.4")
        Write(lastposunits(i), '(i6, a54, f12.4)') iterno, &
            Trim(avgposstr), totnrg
      End Do
    End Do
    
    !** Print out the error log
    Do i=1, 4
      Write(*,'(a,i3,a,i5)') "No. of error of type : ", i, " = ", errorlog(i)
    End Do

    !** Finish up
    Do i=1, nsorbs
      Close(unit=lastposunits(i))
    End Do
    Close(unit=nrgunit)
    mepparams%nsucc = nsucc
    Return
  End Subroutine mepbaker_DoSim

  !------------------------------------------------------------------------
  ! Opens the various position and energy output file by appending the
  ! appropriate extensions
  !------------------------------------------------------------------------
  Subroutine mepbaker_openoutfiles(mepparams, ext, lastposunits, nrgunit)
    Type(MEPBAKER_Params), Intent(in)    :: mepparams
    Character(*), Intent(in)             :: ext
    Integer, Dimension(:), Intent(inout)  :: lastposunits
    Integer, Intent(inout)                :: nrgunit
    
    Integer    :: i, j, sorbtype, nsorbs, mtype
    Type(VecType)                 :: normal
    Character(len=strLen)         :: lastposbasefile, lastposfile
    Character(len=strLen)         :: energyfile, molecname
    
    !** Get some sundry information
    nsorbs       = mepparams%nsorbs

    !** Generate the correct file extensions
    lastposbasefile = Trim(mepparams%basefile)//Trim(ext)
    energyfile      = Trim(mepparams%basefile)//Trim(ext)//".nrg"
    
    nrgunit     = file_getunit(energyfile)
    Open(unit=nrgunit, file=energyfile)
    ! Open the one file for each sorbate type to store the final positions
    ! in a plane
    Do i=1, nsorbs
      mtype = mepparams%mtypes(i)
      molecname  = molecules_name(mtype)
      lastposfile = Trim(lastposbasefile)//"."//molecname
      lastposunits(i) = file_getunit(lastposfile)
      Open(unit=lastposunits(i), file=lastposfile)
      Write(lastposunits(i),'(2a)') "Generalized Coordinates for ", molecname
    End Do
  End Subroutine mepbaker_openoutfiles


  !-------------------------------------------------------------------
  ! Write the various energy components. It returns the total energy
  ! in "totnrg"
  !-------------------------------------------------------------------
  Subroutine mepbaker_dumpenergy(iterno, unitno, totnrg, mtypes)
    Integer, Intent(in) :: iterno
    Integer, Intent(in) :: unitno
    Real(kind=RDbl), Intent(out) :: totnrg
    Integer, Intent(In), Dimension(:) :: mtypes

    Integer, Save :: currentplane = 0
    Integer       :: nsorbs, i, j, mtype1, mtype2
    Real(kind=RDbl)    :: coulnrg, noncoulnrg, totpairnrg
    Logical, Save      :: firsttime = .True.
    Character(len=strLen) :: molecname
    
    !** Number of baker molecules
    nsorbs    = Size(mtypes,1)

    !** Check if this the first call for this plane. If so dump
    !** the header
    If (firsttime) Then
      firsttime = .False.
      ! Write which number corresponds to which sorbate type.  This
      ! will be useful if the files from here are used by other programs.
      Write(unitno, '(a)', Advance='No') "#"
      Do i=1, molecules_getnsorbs()
        mtype1 = i
        molecname = molecules_name(mtype1)
        Write(unitno,'(i2,3a)', Advance='No') i,'=',Trim(molecname)
        If (i /= nsorbs) Then
          Write(unitno,'(a)', Advance='No') ':'
        End If
      End Do
      Write(unitno,*)

      ! Now write the headings of the different columns as numbers
      Write(unitno, '(a,8x)', Advance='No') "#"
!!$      Do i=1, nsorbs-1
!!$        Do j=i+1, nsorbs
      Do i = 1, molecules_getnsorbs()
        Do j = i, molecules_getnsorbs()
          Write(unitno,'(5x,i1,a,i1,5x)', Advance='No') i,"-",j
        End Do
      End Do
      Write(unitno,'(a12)') "Total(kJ)"
    End If

    !** Dump the energies
    totnrg = 0.0_RDbl
    Write(unitno,'(2i5)', Advance='No') iterno
!!$    Do i=1, nsorbs-1
!!$      Do j=i+1, nsorbs
    Do i = 1, molecules_getnsorbs()
      Do j = i, molecules_getnsorbs()
        mtype1 = i
        mtype2 = j
        noncoulnrg = molecules_getnoncoul(mtype1,mtype2,"inst")
        coulnrg = molecules_getcoul(mtype1,mtype2,"inst")
        totpairnrg = (coulnrg + noncoulnrg)
        totnrg = totnrg + totpairnrg
        Write(unitno,'(f14.6)', Advance='No') totpairnrg
      End Do
    End Do
    Write(unitno,'(f14.6)') totnrg
  End Subroutine mepbaker_dumpenergy


!!$  !---------------------------------------------------------------------
!!$  ! Reads the file of stationary points dumped by the routine
!!$  ! "mepbaker_doSim" and generates a file of unique stationary points
!!$  ! using the total energy as the criterion.
!!$  !---------------------------------------------------------------------
!!$  Subroutine mepbaker_genuniq(mepparams)
!!$    Type(MEPBAKER_Params), Intent(in)   :: mepparams
!!$
!!$    Integer                 :: outunit, inunit, npts, error, dof, i, j
!!$    Integer                 :: junk, sitetype, diffunit
!!$    Character(len=strLen)   :: basefile, infilename, outfilename, difffilename
!!$    Type(MinInfo_Params), Dimension(mepparams%niterations)::stpts, sortedstpts
!!$    Type(SortArray_Type), Dimension(mepparams%niterations)::sortarray
!!$    Real(kind=RDbl)         :: energydiff, parttol, pot
!!$    Real(kind=RDbl), Dimension(6) :: df
!!$    Real(kind=Rdbl), Dimension(mepparams%niterations) :: diffarray
!!$
!!$    !** Some sundries
!!$    basefile = mepparams%basefile
!!$    dof      = mepparams%dof
!!$    parttol  = mepparams%parttol
!!$
!!$    !** Open the input file
!!$    If (baker_isSaddleSearch(mepparams%bakerparams)) Then
!!$      infilename = Trim(basefile)//".sdpt"
!!$    Else
!!$      infilename = Trim(basefile)//".min"
!!$    End If
!!$    inunit = file_getunit(infilename)
!!$    Open(unit=inunit, file=infilename, status='old')
!!$    
!!$    !** Open the output file
!!$    If (baker_isSaddleSearch(mepparams%bakerparams)) Then
!!$      outfilename = Trim(basefile)//".uniq.sdpt"
!!$    Else
!!$      outfilename = Trim(basefile)//".uniq.min"
!!$    End If
!!$    outunit = file_getunit(outfilename)
!!$    Open(unit=outunit, file=outfilename)
!!$
!!$    !** Open the file that contains the sorted differences of
!!$    !** the consecutive energies.
!!$    If (baker_isSaddleSearch(mepparams%bakerparams)) Then
!!$      difffilename = Trim(basefile)//".diff.sdpt"
!!$    Else
!!$      difffilename = Trim(basefile)//".diff.min"
!!$    End If
!!$    diffunit = file_getunit(difffilename)
!!$    Open(unit=diffunit, file=difffilename)
!!$
!!$    !** Read the input file
!!$    npts = 0
!!$    Read(inunit,*)  ! Read the comment line
!!$    Do 
!!$      Read(inunit,*, Iostat=error) junk, df, pot
!!$      If (error /= 0) Exit
!!$      npts = npts + 1
!!$      If (npts > mepparams%niterations) Then
!!$        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
!!$            " Maximum number of iterations exceeded specified number"
!!$        Stop
!!$      End If
!!$      stpts(npts)%df(1:dof) = df(1:dof)
!!$      stpts(npts)%pot = pot
!!$    End Do
!!$
!!$    !** Sort the stationary points on the basis of energy
!!$    Call sort_init(stpts(1:npts)%pot, sortarray, npts)
!!$    Call sort_heapsort(sortarray, npts)
!!$    Do i=1, npts
!!$      sortedstpts(i) = stpts(sortarray(i)%index)
!!$    End Do
!!$
!!$    !** Go down the sorted list and based on the difference in energy
!!$    !** of consequtive rows partition the array into dissimilar energies
!!$    sitetype = 1
!!$    sortedstpts(1)%sitetype = sitetype
!!$    Do i=2, npts
!!$      energydiff = Abs(sortedstpts(i)%pot - sortedstpts(i-1)%pot)
!!$      diffarray(i-1) = -energydiff
!!$      If (energydiff > parttol) Then
!!$        ! Start a new site type
!!$        sitetype = sitetype + 1
!!$      End If
!!$      sortedstpts(i)%sitetype = sitetype
!!$    End Do
!!$    Write(*,*) 'Total No. of sites : ', sitetype
!!$
!!$    !** Write out what we have so far
!!$    Do i=1, npts
!!$      Write(outunit,'(i5)', Advance='No') i
!!$      Do j=1, dof
!!$        Write(outunit,'(f8.3)', Advance='No') sortedstpts(i)%df(j)
!!$      End Do
!!$      Write(outunit,'(f13.7, i5)') sortedstpts(i)%pot, sortedstpts(i)%sitetype
!!$    End Do
!!$
!!$    !** Also write out the list of sorted differences
!!$    ! This gives us an easy way to estimate how the tolerance in the 
!!$    ! partitioning will change the number of unique states
!!$    Call sort_init(diffarray, sortarray, npts-1)
!!$    Call sort_heapsort(sortarray, npts-1)
!!$    Do i=1, npts-1
!!$      Write(diffunit,*) i, sortarray(i)%key
!!$    End Do
!!$
!!$    !** Finish up here
!!$    Close(inunit)
!!$    Close(outunit)
!!$    Close(diffunit)
!!$  End Subroutine mepbaker_genuniq


  !---------------------------------------------------------------
  ! Do the finite difference test to see whether the analytical
  ! derivatives are correct.  The finite differencing is done with
  ! a step-size of h
  !-----------------------------------------------------------------
  Subroutine mepbaker_FiniteDiff(mepparams, sorbates, simcell, h)
    Type(MEPBAKER_Params), Intent(inout)           :: mepparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(SimCell_Params), Intent(inout)            :: simcell
    Real(kind=RDbl), Intent(in)                    :: h

    Integer        :: nsorbs, sorb1, sorb2, i, j, err
    Real(kind=RDbl), Dimension(mepparams%dof)   :: x, grad
    Real(kind=RDbl) :: pot, potorig
    Logical         :: pcalcflag, mapflag

    !** Get some sundry information
    nsorbs = mepparams%nsorbs
    
    !** Change the representation of the coordinates to a single array
    Call bakersubs_sorbtoX(sorbates, x, mepparams%mtypes)
    
    !** Calculate the generalized gradients
    Call bakersubs_calcgrad(sorbates,simcell,mepparams%mtypes,x,pot,grad,err)

    !** Change X to do the finite differencing
    mapflag = bakersubs_getenergy(sorbates, simcell, mepparams%mtypes, potorig)
    Write(*,'(/,a,e10.3)') 'Step Size: ', h
    Write(*,'(a7, a8, 2a12)') "Var.", "Value", "Grad", "FinDiff"
    Do i=1, mepparams%dof
      x(i) = x(i) + h
      Call bakersubs_xToSorb(sorbates, simcell, x, mepparams%mtypes)
      mapflag = bakersubs_getenergy(sorbates, simcell, mepparams%mtypes, pot)
      x(i) = x(i) - h
      Call bakersubs_xToSorb(sorbates, simcell, x, mepparams%mtypes)
      Write(*,'(i5, 3f12.6)') i, x(i), grad(i), (pot-potorig)/h
    End Do
  End Subroutine mepbaker_FiniteDiff

  !----------------------------------------------------------------
  ! Display the minimum energy path calculation parameters 
  !----------------------------------------------------------------
  Subroutine mepbaker_display(mepparams, nspc, optunit)
    Type(MEPBAKER_Params), Intent(in) :: mepparams
    Integer, Intent(in)     :: nspc
    Integer, Intent(in), Optional :: optunit

    Character(len=nspc)  :: spc
    Integer              :: unitno, i, j, nsorbs, mtype
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    nsorbs = mepparams%nsorbs
    
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "Minimum Energy Paths using Baker's Algorithm"
    Write(unitno, '(2a)') spc, dashedline
    Call baker_display(mepparams%bakerparams, nspc, unitno)
    Write(unitno, '(2a,i5)') spc, "No. of sorbates: ", nsorbs
    Do i=1, nsorbs
      mtype = mepparams%mtypes(i)
      Write(unitno, '(3a)') spc, "Molecule Name : ",Trim(molecules_name(mtype))
      Write(unitno, '(3a)') spc, "Init Method   : ", &
          mepbaker_inittype(mepparams%initconf(i)%pertflag)
      If (mepparams%initconf(i)%pertflag == 0) Then
        Call nvtmc_display(mepparams%initconf(i)%nvtmcparams, nspc+2, unitno)
      Else If (mepparams%initconf(i)%pertflag == 1) Then
        Call gcmc_displaysimparams(mepparams%initconf(i)%insparams, &
            mepparams%initconf(i)%simno,nspc+2,unitno)
      Else
        Write(unitno, '(a,a,i4)') spc, "No. of centers: ", &
            mepparams%initconf(i)%nsphinsparams
        Do j=1, mepparams%initconf(i)%nsphinsparams
          Write(unitno, '(a,2x,a,a)') spc, "Center: ", Trim(vector_display( &
              mepparams%initconf(i)%sphinsparams(j)%center, "f8.3"))
          Write(unitno, '(a,2x,a,f8.3)') spc, "Min. Radius: ", &
              mepparams%initconf(i)%sphinsparams(j)%minradius
          Write(unitno, '(a,2x,a,f8.3)') spc, "Max. Radius: ", &
              mepparams%initconf(i)%sphinsparams(j)%maxradius
        End Do
      End If
    End Do
  End Subroutine mepbaker_display
End Module mepbaker
