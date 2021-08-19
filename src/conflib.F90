!---------------------------------------------------------------------------
! This module contains the different data structures and routines for
! working with a library (reservoir) of configurations.  The aim is to
! store energy, positions, velocities etc of atoms of a molecule. The
! library will typically represent an ensemble (NVT, NVE, etc.). For
! example it could have different intra-molecular configurations of a
! big molecule like octane.  It also contains basic routines for using
! that library
!
! NOTE: This module uses datafile, species etc which will not be of 
! the same number/type as the main simulation.  Beware..
!---------------------------------------------------------------------------

Module conflib

  Use utils, Only: stripcmnt,isfileopen,toupper,tolower, &
      filesrchstr,allocErrDisplay,deallocErrDisplay,filetoStringArr, &
      cleanstring
  Use defaults, Only: MAX_ATOMS,MAX_SORBS,RDbl,strLen,lstrLen, &
      d_ctrl_file,zero,TOTAL_INDEX,NO_OF_INTRA_POTS
  Use file, Only: file_open, file_getunit, file_settag
  Use random, Only: rranf
  Use atom, Only: atom_getmass, atom_init, atom_display
  Use config, Only: AtMolcoords, config_init, config_display, &
      config_allocfields, config_clean
  Use datafile, Only: CONFILE, datafile_initin, datafile_gencontrolfile, &
      datafile_readconfig, datafile_getnconfs, datafile_readOneMolec, &
      datafile_readOneSpc, datafile_close, datafile_clean, datafile_hasvelocity

  Use general, Only: general_getnoofconfigs, general_init, general_display
  Use molecules, Only: molecules_init, molecules_display, &
      molecules_getnthatom, molecules_getnsorbs,molecs, &
      molecules_name, molecules_gettype, molecules_getnatoms, &
      molecules_getmasses
  Use simcell, Only: SimCell_Params, simcell_init, simcell_display
  Use vector, Only: VecType,Assignment(=),Operator(*),Operator(-), &
      Operator(+), Operator(/), vector_zerovec

  Implicit None
  Save

  Private
  Public :: ConfLibrary, SingleConfig, conflib_basicInit, conflib_dataInit, &
      conflib_getrandomconfig, conflib_clean
  
  !** A single entry in the configuration library
  Type SingleConfig
    Integer                       :: spc      !* species number
    Logical                       :: energies !* True if energies available
    Logical                       :: centered !* COM and net momentum zero if T
    Logical                       :: vels     !* velocities also stored
    Type(VecType)                 :: r0, v0   !* position velocity of com
    Type(VecType), Dimension(:), Pointer :: r !* atomic positions
    Type(VecType), Dimension(:), Pointer :: v !* atomic velocities

    Real(kind=RDbl)                        :: nonintra_nrg
    Real(kind=RDbl), Dimension(:), Pointer :: intra_nrg
    !** nonintra_nrg is the summed noncoul+coul energy 
    !** intra_nrg is the summed intramolecular energy 
    
    Real(kind=RDbl)    :: weight
  End Type SingleConfig 

  !** This defines the Library-Object, a collection of configurations
  Type ConfLibrary
    !** Flags telling what is to be stored
    Logical :: stores_posns, stores_vels, stores_intranrg, chiral_image
    Logical :: centered !* COM and net momentum zero if T

    !** Some info about the source of this library (index by spc number)
    !** spc_origname is the name of the species in the datafile
    Character(len=strlen), Dimension(:), Pointer :: sourcefiles
    Character(len=strlen), Dimension(:), Pointer :: spc_orignames

    Integer                                   :: nconfs  !* number of configs
    Type(SingleConfig), Dimension(:), Pointer :: confs
  End Type ConfLibrary

Contains

  !------------------------------------------------------------------------
  ! Initialises a library for storing the configurations.  It only 
  ! intializes the molecule name other details are handled elsewhere.
  ! Actual memroy allocations and reading of the datafile are done by
  ! conflib_dataInit.  The species name(s) in the data file may be different
  ! that the species name(s) of the library.
  ! Requires:  library -- configuration library data type
  !            molecname -- name(s) of species
  !            withVels -- True means that velocities will be stored
  !            centered -- True => shift to COM and net momentum zero
  !            opt_spcname -- optional species name(s) in datafile
  !------------------------------------------------------------------------
  Subroutine conflib_basicInit(library,molecnames,withVels,centered,opt_spcnames)
    Type(conflibrary), Intent(InOut)                 :: library
    Character(len=strLen), Dimension(:), Intent(In)  :: molecnames
    Logical, Intent(In)                              :: withVels,centered
    Character(len=strLen), Dimension(:), Intent(In), Optional :: opt_spcnames

    Integer     :: i, spcnum, spc, natoms, error, nspc

    !** Set default flags, store intranrg, posns, and if needed velocities
    library%stores_posns = .True.

    ! this vels will be changed later based on whether datadile stores 
    ! it or not
    library%stores_vels = withVels
    library%stores_intranrg = .True.
    library%centered = centered

    !** Allocate space for the origin information for each species
    nspc = molecules_getnsorbs()
    Allocate(library%sourcefiles(nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    
    Allocate(library%spc_orignames(nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    
    library%sourcefiles = 'NULL'
    Do spc = 1,nspc
      library%sourcefiles(spc) = 'TOREAD'
      library%spc_orignames(spc) = molecules_name(spc)
    End Do

    !** Set the origin species name(s) if they are to be different
    If (Present(opt_spcnames)) Then
      Do i = 1,Size(molecnames)
        spc = molecules_gettype(molecnames(i))
        If (molecnames(i) /= library%spc_orignames(spc)) Then
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              ' input molecule name does not match that in molecules ',&
              Trim(molecnames(i))
          Stop
        End If
        library%spc_orignames(spc) = opt_spcnames(i)
      End Do
    End If

    !** Nullify storage for configurations
    library%nconfs = 0
    Nullify(library%confs)

  End Subroutine conflib_basicInit

  !--------------------------------------------------------------------
  ! Allocates memory and reads in data from a datafile for a single 
  ! species type.
  ! Requires:  library -- configuration library data type
  !            molecname -- name of species
  !            datafile_name -- name of datafile to be read
  !            optnconfs -- opt number of configurations from datafile
  !--------------------------------------------------------------------
  Subroutine conflib_dataInit(library, molecname, datafile_name, opt_nconfs)
    Type(conflibrary), Intent(InOut)   :: library
    Character(len=strLen), Intent(In)  :: molecname
    Character(len=strLen), Intent(In)  :: datafile_name
    Integer, Intent(In), Optional      :: opt_nconfs

    Integer             :: datfile_spcnum, nspcs_datafile, nconfs, spc, configno
    Integer             :: i, j, count, ndats, natoms, nmolecs, init_nmolecs
    Integer             :: oldnconfs
    Logical             :: energies
    Real(kind=RDbl)     :: intra_nrg,nonintra_nrg
    Type(CONFILE)       :: configfile
    Type(AtMolCoords)   :: tempConf
    Type(VecType), Dimension(:,:), Pointer :: rvecs, vvecs

    !** Nullify the pointers
    Nullify(rvecs)
    Nullify(vvecs)

    !** Store the source file name
    spc = molecules_gettype(molecname)
    library%sourcefiles(spc) = datafile_name

    !** Get the number of atoms in this species type
    natoms = molecules_getnatoms(spc)

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,'(a,i3,2a)') "Initializing library for spc number: ", spc, &
        " from data file : ", Trim(datafile_name)

    !** Allocate space for the temporary configuration of this species type
    Call config_allocfields(tempConf, spc, 1)

    !** Open the data file and begin reading
    Call datafile_initin(configfile, datafile_name)
    If (datafile_hasvelocity(configfile)) Then
      Write(*,*)  "datafile has velocity stored"
      Write(*,*)  "this can be used for both H-GCMC and MD-CGCC"

    Else
      Write(*,*)  "datafile has no velocity"
      Write(*,*)  "this should be used only with H-GCMC, not MD-CGCC"
    Endif

    nspcs_datafile = configfile%xtrainfo%nspcs  !** HACK

    !** Set and check the number of configurations to be read
    nconfs = datafile_getnconfs(configfile)
    If (Present(opt_nconfs)) Then
      If (nconfs < opt_nconfs) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            ' data file does not contain enough configurations'
        Stop
      Else
        nconfs = opt_nconfs
      End If
    End If

    !** Datafile might have different name for species
    !** Find the index number of the species in the datafile
    Write(*,*) "going to check whether datafile contains "&
        //Trim(library%spc_orignames(spc))
    datfile_spcnum = -1
    Do i = 1,nspcs_datafile
      If (cleanstring(configfile%xtrainfo%spcnamelist(i)) == &
          cleanstring(library%spc_orignames(spc))) Then
        datfile_spcnum = i
        Exit
      End If
    End Do
    If (datfile_spcnum == -1) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' datafile does not contain specified species ',&
          Trim(library%spc_orignames(spc))
      Stop
    End If

    !** Read thru the whole datafile
    count = 0
    Do i = 1,nconfs
      energies = .True.
      Call datafile_readOneSpc(configfile, rvecs, vvecs, datfile_spcnum, &
          intra_nrg, nonintra_nrg, nmolecs)

      !** If more than one molecule, cannot know individual energies, flag this
      If (nmolecs > 1) Then
        energies = .False.
      End If

      If (i == 1) Then
        Write(*,'(a,i8)') "Estimated size of library : ", nconfs*nmolecs
        If (.Not.datafile_hasvelocity(configfile)) Then
          Write(*,'(a)') "library does not have velocities, assumed=0"
          library%stores_vels=.False.
        Endif
        !*** We can size the library only after we read the first configuration
        !*** because number of molecules will be known only then
        ! somewhat of a hack but no choice-- Shaji
        oldnconfs = library%nconfs
        Call conflib_resize(library,(library%nconfs + nconfs*nmolecs), &
            spc,energies)

        !** init_nmolecs will be used as the number of molecules of this spc
        init_nmolecs = nmolecs
      End If

      !** Make sure there are enough molecules to read
      If (nmolecs < init_nmolecs) Then
        Write(0,'(1x,2a,i4,2a,2i4)') __FILE__," : ",__LINE__, &
            ' While reading datafile: number of molecules in configuration', &
            ' is too small ',nmolecs,init_nmolecs
        Stop
      End If

      !** Place the coordinates in the library
      Do j = 1,init_nmolecs
        count = count + 1
        configno = oldnconfs + count
        ! this is necessary to make sure that junk data from 
        ! datafile are not enterd into the library
        If (library%stores_vels) Then
          tempConf%coords(1:natoms,1)%v = vvecs(1:natoms,j)
        Else
          tempConf%coords(1:natoms,1)%v=vector_zerovec() 
        Endif
        tempConf%coords(1:natoms,1)%rp = rvecs(1:natoms,j)
        If (energies) Then
          Call conflib_update(library%confs(configno),tempConf, &
              intra_nrg,nonintra_nrg)
        Else
          Call conflib_update(library%confs(configno),tempConf)
        End If
      End Do
    End Do

    !** Deallocate the temporary configuration
    Call config_clean(tempConf)

    !** Close the data file
    Call datafile_close(configfile)
    Call datafile_clean(configfile)

  End Subroutine conflib_dataInit

  !------------------------------------------------------------------------
  ! Update a particular single configuration
  ! Requires:  library -- configuration library data type
  !            coords -- configuration data to insert into library config
  !            intra_nrg -- intramolecular energy for this configuration
  !            nonintra_nrg -- nonintramolecular energy for this config
  !------------------------------------------------------------------------
  Subroutine conflib_update(conf,coords,intra_nrg,nonintra_nrg)
    Type(SingleConfig), Intent(InOut)      :: conf
    Type(AtMolCoords), Intent(In)          :: coords
    Real(kind=RDbl), Intent(In), Optional  :: intra_nrg, nonintra_nrg

    Integer                   :: a,natoms
    Real(kind=RDbl)           :: m_sum, mass
    Type(VecType)             :: m_r_sum, m_v_sum
    Real(kind=RDbl), Dimension(MAX_ATOMS*10) :: masses

    !** Get number of atoms
    natoms = molecules_getnatoms(conf%spc)

    !** Center the positions and velocities if necessary
    If (conf%centered) Then
      !** Zero the sums
      m_sum = zero
      m_r_sum = zero
      m_v_sum = zero
  
      !** Get the species number, number of atoms and masses
      natoms = molecules_getnatoms(conf%spc)
      Call molecules_getmasses(conf%spc,masses)
  
      Do a = 1,natoms
        mass    = masses(a)
        m_sum   = mass + m_sum
        m_r_sum = (mass * coords%coords(a,1)%rp) + m_r_sum
        m_v_sum = (mass * coords%coords(a,1)%v)  + m_v_sum
      End Do

      conf%r0 = m_r_sum/m_sum
      conf%v0 = m_v_sum/m_sum

    Else
      conf%r0 = 0.0_RDbl
      conf%v0 = 0.0_RDbl
    End If

    !** Save the coordinates, recentering if necessary
    Do a = 1,natoms
      conf%r(a) = coords%coords(a,1)%rp - conf%r0
    End Do
    If (conf%vels) Then
      Do a = 1,natoms
        conf%v(a) = coords%coords(a,1)%v  - conf%v0
      End Do
    End If

    !** Save the energies
    If (conf%energies) Then
      If ((.Not. Present(intra_nrg)).Or.(.Not. Present(intra_nrg))) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            ' library configuration energy flag on, must pass energies'
        Stop
      End If

      conf%intra_nrg(TOTAL_INDEX) = intra_nrg
      conf%nonintra_nrg = nonintra_nrg
    End If

  End Subroutine conflib_update
  
  !-------------------------------------------------------------------
  ! Resize the storage for the set of configurations
  ! Requires:  library -- configuration library data type
  !            nconfs -- new desired number of configurations
  !            newspc -- spc number for new configurations
  !            energies -- True if energies to be stored
  !-------------------------------------------------------------------
  Subroutine conflib_resize(library,nconfs,newspc,energies)
    Type(conflibrary), Intent(InOut)   :: library
    Integer, Intent(In)                :: nconfs,newspc
    Logical, Intent(In)                :: energies

    Integer         :: i, error, spc, noldconfs
    Type(SingleConfig), Dimension(:), Allocatable  :: tempconfs

    !** Return immediately if no resizing is needed
    If (nconfs == library%nconfs) Return

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'wanted: ',nconfs,'  current: ',library%nconfs

    !** Save existing configurations
    noldconfs = library%nconfs
    If (noldconfs > 0) Then
      Allocate(tempconfs(noldconfs), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Do i = 1,noldconfs
        Call conflib_initcopyconfig(tempconfs(i),library%confs(i))
      End Do
    End If

    !** Clean any existing configurations
    If (Associated(library%confs)) Then
      Do i = 1,library%nconfs
        Call conflib_cleanconfig(library%confs(i))
      End Do

      Deallocate(library%confs, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End If

    !** Allocate space for new configuration set
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) nconfs
    Allocate(library%confs(nconfs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    library%nconfs = nconfs

    !** Copy existing configurations back into set
    If (noldconfs > 0) Then
      Do i = 1,noldconfs
        Call conflib_initcopyconfig(library%confs(i),tempconfs(i))
      End Do
    End If

    !** Create space for new configurations
    Do i = (noldconfs+1),library%nconfs
      Call conflib_sizeconfig(library%confs(i),newspc,energies, &
          library%centered,library%stores_vels)
    End Do

    !** Deallocate space for temporary configurations
    If (noldconfs > 0) Then
      Do i = 1,noldconfs
        Call conflib_cleanconfig(tempconfs(i))
      End Do
      Deallocate(tempconfs, STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)
    End If

  End Subroutine conflib_resize

  !--------------------------------------------------------------------
  ! Size a single configuration
  ! Requires:  conf -- single configuration data type
  !            spc -- species number
  !            energies -- True if energies to be stored
  !            centered -- True => shift to COM and net momentum zero
  !            velocities -- True => initialize velocities also
  !--------------------------------------------------------------------
  Subroutine conflib_sizeconfig(conf,spc,energies,centered,velocities)
    Type(SingleConfig), Intent(InOut)     :: conf
    Integer, Intent(In)                   :: spc
    Logical, Intent(In)                   :: energies,centered,velocities

    Integer         :: natoms, error

    !** Store the species number and get number of atoms
    conf%energies = energies
    conf%centered = centered
    conf%vels = velocities
    conf%spc = spc
    natoms = molecules_getnatoms(spc)

    !** Allocate space for the configuration's postitions
    Allocate(conf%r(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Allocate space for the configuration's velocities if necessary
    If (conf%vels) Then
      Allocate(conf%v(natoms), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Else
      Nullify(conf%v)
    End If

    !** Allocate space for the configuration's intramolecular energies
    Allocate(conf%intra_nrg(NO_OF_INTRA_POTS), STAT=error)          
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

  End Subroutine conflib_sizeconfig

  !-------------------------------------------------------------
  ! Copy a single configuration into a new one.  The destination
  ! configuration type must already be initialized.
  ! Requires:  dest -- destination configuration data type 
  !            orig -- original configuration data type 
  !-------------------------------------------------------------
  Subroutine conflib_copyconfig(dest,orig)
    Type(SingleConfig), Intent(InOut)     :: dest
    Type(SingleConfig), Intent(InOut)     :: orig

    !** Copy the simple stuff
    dest%spc = orig%spc
    dest%energies = orig%energies
    dest%centered = orig%centered
    dest%r0 = orig%r0
    If (orig%vels) dest%v0 = orig%v0
    dest%weight = orig%weight

    !** Copy the energies
    If (orig%energies) Then
      dest%nonintra_nrg = orig%nonintra_nrg
      dest%intra_nrg = orig%intra_nrg
    Else
      dest%nonintra_nrg = zero
      dest%intra_nrg = zero
    End If

    !** Copy the positions and velocities
    dest%r = orig%r
    If (orig%vels) dest%v = orig%v

  End Subroutine conflib_copyconfig

  !-------------------------------------------------------------
  ! Initialize and Copy a single configuration into a new one.  
  ! Requires:  dest -- destination configuration data type 
  !            orig -- original configuration data type 
  !-------------------------------------------------------------
  Subroutine conflib_initcopyconfig(dest,orig)
    Type(SingleConfig), Intent(Out)       :: dest
    Type(SingleConfig), Intent(InOut)     :: orig

    !** Size the destination configuration
    Call conflib_sizeconfig(dest,orig%spc,orig%energies,orig%centered,orig%vels)

    !** Copy to destination
    Call conflib_copyconfig(dest,orig)

  End Subroutine conflib_initcopyconfig

  !------------------------------------------------------------------
  ! Picks a random index and returns the corresponding molecule's
  ! coordinates from the library, returns the "intra" energy of that 
  ! configuration also.
  ! Requires:  library -- configuration library data type
  !            comPos -- positions of al atoms wrt Center of Mass
  !            comVel -- velocities of al atoms wrt Center of Mass
  !            r0, v0 -- position,  velocity of Center of Mass
  !            intra -- total intramolecular energy stored here
  !------------------------------------------------------------------  
  Subroutine conflib_getrandomconfig(library,comPos,comVel,r0,v0,intra)
    Type(conflibrary), Intent(In)            :: library
    Type(VecType), Dimension(:), Intent(Out) :: comPos, comVel 
    Type(VecType), Intent(Out)               :: r0, v0 
    Real(kind=RDbl), Intent(Out)             :: intra  

    Integer     :: index,a,natoms
    
    !** Pick the configuration, Int() truncates
    index  = Int(library%nconfs * rranf()) + 1

    If (.Not. library%confs(index)%centered) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' WARNING: conflib_getrandomconfig is returning uncentered config'
    End If

    !** Return the promised quantities
    natoms = molecules_getnatoms(library%confs(index)%spc)
    comPos(1:natoms) = library%confs(index)%r
    r0               = library%confs(index)%r0
    v0               = library%confs(index)%v0
    intra            = library%confs(index)%intra_nrg(TOTAL_INDEX)

    !** Return the velocities if the are stored
    If (library%confs(index)%vels) Then
      comVel(1:natoms) = library%confs(index)%v
    Else
      Do a = 1,natoms
        comVel(a) = 0.0_RDbl
      End Do
    End If

  End Subroutine conflib_getrandomconfig

  !--------------------------------------------------------
  ! Returns the number of configurations stored 
  ! Requires:  library -- configuration library data type
  !--------------------------------------------------------
  Integer Function conflib_getnconfs(library)
    Type(conflibrary), Intent(In) :: library

    conflib_getnconfs = library%nconfs

  End Function conflib_getnconfs

  !------------------------------------------------------------------
  ! Cleans a configuration library data type
  ! Requires:  library -- configuration library data type
  !------------------------------------------------------------------  
  Subroutine conflib_clean(library)
    Type(conflibrary), Intent(InOut)         :: library

    Integer     :: i,error

    !** Clean the individual configurations
    Do i = 1,library%nconfs
      Call conflib_cleanconfig(library%confs(i))
    End Do
    
  End Subroutine conflib_clean

  !------------------------------------------------------------------
  ! Cleans a single configuration
  ! Requires:  conf -- single configuration data type
  !------------------------------------------------------------------  
  Subroutine conflib_cleanconfig(conf)
    Type(SingleConfig), Intent(InOut)         :: conf

    Integer     :: error

    !** Deallocate position storage
    Deallocate(conf%r, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)

    !** Deallocate velocity storage
    If (conf%vels) Then
      Deallocate(conf%v, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    End If

    !** Deallocate interactions information
    Deallocate(conf%intra_nrg, STAT=error )
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)

  End Subroutine conflib_cleanconfig

End Module conflib




