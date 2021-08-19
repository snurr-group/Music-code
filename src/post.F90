!------------------------------------------------------------------------------
! This is the module that deals with everything u need to do after the 
! music production run is completed. [ It can also be used while the run 
! is in progress, y simply specifying in the controlfile what 
! percentaage of data you are interested in ]
! Using one of the "datafile" the orginal control file is regenerated 
! and Then used to, initialize the various modules (atoms, molecules, simcell)
! and calculate statistics from the configuration file. For very specialized
! things, (like Diffusivity) It calls separate modules.
!------------------------------------------------------------------------------
Module post

  Use atom, Only: atom_getmass, atom_init
  Use config, Only: AtMolCoords, config_init, default_config_tag, &
      config_initdefaultsize, config_initmolecule, config_initfixed, &
      config_isfixed, config_checkandincr, config_getnmoles, &
      config_getdistance, config_getangle, config_gettorangle, &
      config_getMolecCOM, config_getCOMVel, config_setnmoles, config_applypbcs
  Use datafile, Only: datafile_initin, datafile_gencontrolfile, CONFILE, &
      datafile_readconfig, datafile_hasposition, datafile_hastime, &
      datafile_getextrainfo, datafile_setContentTag
  Use defaults, Only: RDbl, RSgl,strLen, lstrLen, d_ctrl_file, dashedline, &
      NO_OF_INTRA_POTS, TOTAL_INDEX, zero, one, pi, radToDeg, MAX_ATOMS,  &
      MAX_SORBS, calToJ,xlstrLen
  Use file, Only: file_open, file_close, file_settag, file_getname
  Use general, Only: general_init, general_getnoofconfigs, genParams
  Use histogram, Only: Histogram_Params, histogram_update, histogram_display, &
      histogram_reinit, histogram_init
  Use interact, Only: Interaction_Model, interact_init, interact_display, &
      interact_resize, interact_int, interact_changeAllNmoles
  Use mdpc, Only: Mdpostparams, MDScratch, mdpc_init, mdpc_display, &
      mdpc_createscratchfiles, mdpc_calcselfdiffusivity, mdpc_samplecf, &
      mdpc_calcRadProf
  Use mmath, Only : mmath_LSFit
  Use molecules, Only: molecules_getnthatom, molecules_getnatoms, &
      molecules_getnsorbs, molecules_name, molecules_gettype, molecs, &
      molecules_init, molecules_exists
  Use multihist, Only: MultiHist_Params, multihist_init, multihist_update
  Use simcell, Only: Simcell_Params, simcell_init, simcell_getzeoell, &
      simcell_getnuc, simcell_maptouc, simcell_getell, simcell_minimage, &
      simcell_getvolume
  Use smap, Only: Smap_Params, smap_init, smap_display, smap_getPosSiteType, &
      smap_getallsitetypes, smap_showxyz
  Use stats, Only: Statistics, stats_init, stats_update, stats_getcavg
  Use storebase, Only:  storebase_init,EnergyPlus,storebase_disp
  Use storestats, Only: Species_Stats, storestats_init, storestats_getcoul, &
      storestats_getnoncoul, storestats_getke, storestats_getintranrg
  Use storetop, Only: storetop_fastsum, storetop_initcopy, &
      Forcefield_Results, storetop_fillsub, storetop_setmap, storetop_extract
  Use utils, Only: filesrchstr, stripcmnt, genfilename, split, toint, toreal, &
      toupper, allocErrDisplay, int2str, real2str, cleanstring,str2seq
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag

  Implicit None
  Save 

  Private 
  Public :: post_init, post_getnfiles, post_getConfigFilename, &
      post_reInitStats, post_initconfileReading, post_samplecf, post_display, &
      post_writeAvgstoFile, post_readandavg, post_calcMSD, PostInfo, &
      nrgAvg, loadingAvg, sitingAvg, PoreMapParams, &
      structureAvg, DiffusParams, velocityParams, interactParams, radDistAvg, &
      post_calcRadProf

  Type PostInfo
    ! Basic info 
    Type(AtMolCoords), Dimension(:), Pointer :: species
    Type(Simcell_Params)  :: scell ! Simulation cell information
    Character(len=strLen) :: baseFilename ! configuration base filename
    Integer               :: first, last  ! the first and last config file nos 
    Character(len=strLen) :: newCtrlFile  ! new control filename
    Character(len=strLen) :: outFilename  ! base output filename
    Character(len=strLen) :: simType      ! GCMC or MD etc...
    Type(CONFILE), Dimension(:), Pointer :: configFiles
    Real(kind=RDbl)                      :: startskip, endskip ! in %
    Integer                              :: firstconf, lastconf

    ! Analysis Types
    Type(nrgAvg), Pointer                :: nrgs
    Type(loadingAvg), Pointer            :: loading
    Type(sitingAvg), Pointer             :: siting 
    Type(poremapParams), Pointer         :: poremap
    Type(structureAvg), Pointer          :: structure
    Type(DiffusParams), Pointer          :: MSD
    Type(velocityParams), Pointer        :: vels 
    Type(interactParams), Pointer        :: iact 
    Type(radDistAvg), Pointer            :: radDist

    ! array containing number of moles of each species
    Integer, Dimension(:), Pointer       :: nmole_list

    ! Energies
    Type(Species_Stats) :: spcstats ! in KJ/mol

    Real(kind=RDbl)   :: time 

    ! To be removed later
    Logical :: datfiles_are_old

  End Type PostInfo

  Type nrgAvg
    Integer                                    :: nblocks,statblocksize
    Real(kind=RDbl),Dimension(:,:), Pointer    :: output_array
    Type(Statistics)                           :: total
    Type(Statistics)                           :: totnrg_blockavgs 
    Type(Statistics) , Dimension(:,:), Pointer :: coulNrg
    Type(Statistics) , Dimension(:,:), Pointer :: noncoulNrg
    Type(Statistics) , Dimension(:), Pointer   :: ke 
    Type(Statistics) , Dimension(:,:), Pointer :: intra 
  End Type nrgAvg

  Type loadingAvg
    Integer                                   :: nblocks, statblocksize
    Real(kind=RDbl), Dimension(:,:,:), Pointer  :: output_array
    Type(Statistics) , Dimension(: ), Pointer :: totalavgs 
    Type(Statistics) , Dimension(: ), Pointer :: blockavgs ! not used now
    Real(kind=RDbl), Dimension(:,:,:), Pointer :: isotherm
  End Type loadingAvg

  Type sitingAvg
    Integer                                     :: nblocks, statblocksize
    Type(Smap_Params),Pointer                   :: sitemap
    Type(Simcell_Params),Pointer                :: scell
    Character(len=strlen)                       :: sitemapfile
    Type(Statistics) , Dimension(:,: ), Pointer :: siteavgs 
    Integer, Dimension(:,:), Pointer            :: site_inst
    Integer, Dimension(:),Pointer               :: sites
    Integer                                     :: nsites

    !for outputting Stats in  VERRBOSE mode
    Real(kind=RDbl), Dimension(:,:,:,:), Pointer  :: output_array

    Logical :: com_based            ! avgs based on COM
    Integer :: atom_num             ! or atom_num

    ! to plot isotherm at end
    Real(kind=RDbl), Dimension(:,:,:), Pointer :: isotherm

  End Type sitingAvg

  Type poremapParams
    Integer                                     :: nmaps
    Character(len=strlen),Dimension(:), Pointer :: sorbnames, mapname
    Real(kind=RDbl), Dimension(:,:,:,:), Pointer  :: fillArray 
    Real(kind=RDbl), Dimension(3)               :: uc_edge_l
    Integer, Dimension(3)                       :: nbins
    Integer, Dimension(3)                       :: uc_repeat_info
    Integer :: ncont
    Logical :: onlycom, atom_list_given, cont_div_exponential, dump_dmap
    Integer,Dimension(:),Pointer :: atom_list ! used only when nmaps==1
    Integer                      :: natoms    ! size of atom_list 
    Real(kind=RDbl) :: scaling_factor
    Real(kind=RDbl)                             :: dx, dy, dz
  End Type poremapParams

  !** Each of the members of stat corresponds to the average for one 
  !** particular intra molecular characteristc of particular sorb Types.
  Type structureAvg
    Integer                                     :: nparams
    Real(kind=RDbl), Dimension(:), Pointer      :: minval, maxval
    Integer, Dimension(:), Pointer              :: ndivs 
    Type(Histogram_Params),Dimension(:), Pointer      :: hist 
    Character(len=strlen), Dimension(:), Pointer:: intra_type, sorbnames
    Integer, Dimension(:,:), Pointer            :: atom_list
  End Type structureAvg

  !** creates histogram for velocity of each of the atoms of the molecule
  !** histograms correspond to x,y,z compoenents of velocity of each 
  !** atom(atom_comp_hist), or totalvelocity of atom (atom_hist), 
  !** or com velocity of molecule
  Type velocityParams
    Character(len=strLen)        :: spcname
    Integer                               :: spctype, natoms
    Real(kind=RDbl)                       :: minVel, maxVel, vol
    Integer                               :: ndivs 
    Type(Histogram_Params),Dimension(:,:), Pointer      :: atomcomp_hist
    Type(Histogram_Params),Dimension(:), Pointer      :: atom_hist
    Type(Histogram_Params),Pointer      :: com_hist
    Integer, Dimension(:), Pointer      :: display_list
  End Type velocityParams


  !** anything to do with nrg calcs
  !** for now we will have only one imodel
  Type interactParams
    Integer                               :: nimodels
    Type(Interaction_Model)               :: imodel
    Real(kind=RDbl), Dimension(:,:,:), Pointer :: coul, ncoul
    Real(kind=RDbl), Dimension(:,:), Pointer :: intra
    Type(Forcefield_Results),Dimension(:), Pointer  :: molsys_nrg
  End Type interactParams

  Type DiffusParams
    Type(Mdpostparams)                          :: mdpp
    Type(MDScratch), Dimension(:), Pointer   :: scratch
    Real(Kind=RDbl), Dimension(:,:), Pointer :: dself
  End Type DiffusParams

  !** For verious radial distribution functions
  !** Only one pair can be processed at a time
  Type radDistAvg
    Real(kind=RDbl)                             :: mindist, maxdist, vol
    Integer                                     :: nbins, ncalls, n_loops
    Type(Histogram_Params)                      :: hist
    Character(len=strlen),Dimension(2)          :: sorbnames
    Integer,Dimension(2)                        :: sorbs, atoms
    Logical,Dimension(2)                        :: com_based
    Type(Simcell_Params),Pointer                :: scell
  End Type radDistAvg



  Character(strLen), Parameter :: postIdString = "Post Processor Information"
  Character(strLen), Parameter :: postNrgAvgIdString = &
      "Post : Energy Average Info "
  Character(strLen), Parameter :: postLoadingIdString =&
      "Post : Loading Average Info "
  Character(strLen), Parameter :: postSitingIdString =&
      "Post : Siting Average Info "
  Character(strLen), Parameter :: postPoreMapIdString =&
      "Post : Pore Map Info"
  Character(strLen), Parameter :: postStructureIdString =&
      "Post : Structure Averages Info"
  Character(strLen), Parameter :: postMSDIdString =&
      "MD Post Code Information"
  Character(strLen), Parameter :: postVelHistIdString =&
      "Post : Velocity Histogram"
  Character(strLen), Parameter :: postInteractString =&
      "Post : Interaction Info"
  Character(strLen), Parameter :: postRadDistIdString =&
      "Post : Radial Distribution Info"

  ! to be removed later, as of now helps to deal with *.con files 
  ! made before feb-2002
  Character(strLen), Parameter :: postConFileAgeString =&
      "Age of datafile"

  ! general format for writing columns, need a utils for generating this
  Character(strLen), Dimension(7) :: REALCOLFORMAT = (/ &
      '(1f12.7)', '(2f12.7)','(3f12.7)','(4f12.7)','(5f12.7)','(6f12.7)',&
      '(7f12.7)' /)

  Integer, Parameter :: post_MAX_MOLEC=200

#ifdef SITE_HACK
  ! some hack combines IS and SC

  ! 3 sorbs, 3sites =1,2,3 (IS, SC, ZZ) head, and tail
  Real(kind=RDbl) , Dimension(3,3,3) :: HTSTATS=zero

  ! 3 sorbs, 3 sorbs, 3sites =1,2,3 (IS, SC, ZZ)
  Type(Statistics), Dimension(3,3,3) :: S_NC_STATS
  Type(Statistics), Dimension(3,3) :: S_INTRA_STATS

  Integer :: n_single_site=0, n_mixed=0, sunit1, sunit2 
#endif

Contains

  !----------------------------------------------------------------------------
  ! Initializes the post routine by reading the information from the
  ! control file. Pass off anything we don't know how to initialize.
  !----------------------------------------------------------------------------
  Subroutine post_init(postobj, controlFile, optTag)
    Type(PostInfo), Intent(InOut)           :: postobj
    Character(*), Intent(In)                :: controlFile
    Character(strLen), Optional, Intent(In) :: optTag
    Character(len=strLen),Dimension(strLen) :: fields
    Integer           :: unitno, lineno, nfields, nsorbates, i, error, ndats
    Character(strLen) :: tag, text, config_init_source, simType

    !** Decide which tag to use, the optional tag or default
    If (Present(optTag)) Then 
      tag = Trim(optTag)
    Else
      tag = Trim(postIdString)
    End If

    !** Nullify the pointers in the post variable
    Call post_nullify(postobj)
    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
        " Opening the post control file and reading it"
    !** Open the file if it isn't already open
    unitno = file_open(controlFile)

    !** Find the post control file tag
    lineno = filesrchstr(unitno,tag,text)

    !** Read in the type of simulation and decide how to initialize config
    Read(unitno,'(a)') text
    simType = Trim(stripcmnt(text))
    simType = toupper(simType)
    postobj%simType=simType
    !** Read in the base filename for the configuration files
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    postobj%baseFilename = trim(text)

    !** Read in the numbers of the config files to process
    !** Remember, these configuration files should share a
    !** common control file. For example, these could be 
    !** the config files from a GCMC simulation at different
    !** pressures.
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    !** split the line to get the start and end
    nfields = split(text,fields,',')
    postobj%first = toint(fields(1))
    If (nfields /= 1) Then
      postobj%last = toint(fields(2))
    Else
      postobj%last = toint(fields(1))
    End If

    !** Read in the name of the new control file to generate
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    postobj%newCtrlFile = trim(text)

    !** Read in the base name for the output files generated
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    postobj%outFilename = trim(text)

    !** Read in the initial percentage and final percentage of data to 
    !** be skipped. Example if startskip= 15 and endskip=33 then only data from
    !** 15% to 67% will be analysed. 
    Read(unitno,'(a)') text
    text = stripcmnt(text)

    !** split the line to get the start and end
    nfields = split(text,fields,',')
    postobj%startSkip = toreal(fields(1))
    postobj%endSkip = toreal(fields(2))

    If ( (postobj%startSkip>100).Or.(postobj%startSkip<0).Or. &
        (postobj%endSkip>100).Or.(postobj%endSkip<0).Or.      &
        ((postobj%startSkip+postobj%endSkip)>100) )         Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,&
          "      Wrong value of endskip/startskip specified in control file"
      Stop
    End If

    !** decide how to initialize species
    Select Case (simType)
    Case("GCMC", "CBGCMC", "HGCMC" )
      config_init_source= "DEFAULT"
    Case("MD")
      config_init_source= "CONTROLFILE"
    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Wrong simType specified in ctrl file"
      Stop
    End Select

    !** Make the configuration file list, generate and store their names
    Call post_makeConfigList(postobj,postobj%baseFilename,postobj%first, &
        postobj%last)

    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
        " Looking for section specifying old or new datafile "

    ! to be removed later, as of now helps to deal with *.con files 
    ! which were made before feb-2002
    ! This routine can be called only after post_makeconfiglist
    Call post_initConFilesAge(postobj,controlFile)


    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
        " Re-Generating the ctrlfile used for original simulation "
    !** Reconstruct the original control file from the first
    !** configuration file
    Call post_extractControlFile(postobj%configFiles(1), &
        postobj%newCtrlFile)

    Write(*,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
        " Reading the Re-Generated ctrl file to init  : ", &
        "genral/atom/molecules/simcell/config"
    !** Read in the control file to initialize ONLY : atoms,
    !** molecules, simcell, configuration, and (MAYBE?) forcefield
    Call post_readControlFile(postobj,postobj%newCtrlFile,config_init_source)

    !-------------------------------------------------------------------
    !** CAUTION : All the below operations can be done only after atms, *** 
    !** molecules, general; etc has been initialized.                   ***
    !-------------------------------------------------------------------

    ! Assumes all datafiles have same number of configs
    ! should be changed later
    nsorbates=molecules_getnsorbs()
    ndats=general_getnoofconfigs()
    postobj%firstconf= Int(ndats*(postobj%startskip/100))+1
    postobj%lastconf=ndats-Int(ndats*(postobj%endskip/100))

    !** array with number of moles 
    Allocate(postobj%nmole_list(nsorbates), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"loading-inst")
    postobj%nmole_list=0

    !** Initialize the species energies statistics variable. This collects
    !** information about the configuration energy values (kinetic, 
    !** potential coulombic, potential noncoulombic, intramolecular)
    Call storestats_init(postobj%spcstats, .False. )


    Do i=1,Size(postobj%configFiles)
      If (postobj%configFiles(i)%datafile_is_old) &
          Call datafile_setContentTag(postobj%configFiles(i))
    End Do

    !** Initialize parameters used for averaging, making maps, MSD etc.
    Call post_initAnalTypes(postobj,controlFile)

    !** set content-tag from general
  End Subroutine post_init

  !--------------------------------------------------------------------
  ! Checks whether the *.con.* file is old or not
  ! This routine can be called only after post_makeconfiglist
  !--------------------------------------------------------------------
  Subroutine post_initConFilesAge(postobj,ctrlfile)
    Type(PostInfo), Intent(InOut)     :: postobj 
    Character(len=strLen), Intent(In) :: ctrlfile
    Integer :: lineno, i, unitno
    Character(len=strLen) :: tag, text
    Logical :: rewind_flag, oldflag

    rewind_flag=.true.
    tag=postConFileAgeString
    unitno=file_open(ctrlfile)
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno==0) Then
      ! the post control file does contain any information regarding age 
      ! of datafile
      oldflag=.False.
    Else
      Read(unitno,*) text
      text=Adjustl(stripcmnt(text))
      Select Case (text)
      Case("OLD_DATAFILE")
        oldflag=.True.
      Case("NEW_DATAFILE")
        oldflag=.False.
      Case default
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Not able to understand the option :"//Trim(text)
        stop
      End Select
    Endif

    postobj%datfiles_are_old=oldflag
    If (oldflag) Then
      Write(*,*) " Make sure that this post code is being run on old  &
          & *.con files ( made before 02-01-2002) "
      Do i=1,Size(postobj%configFiles)
        postobj%configFiles(i)%datafile_is_old=.True.
      End Do
    Else
      Write(*,*) " Make sure that this post code is being run on new  &
          & *.con files ( made after 02-01-2002) "
      Do i=1,Size(postobj%configFiles)
        postobj%configFiles(i)%datafile_is_old=.False.
      End Do
    Endif
  End Subroutine post_initConFilesAge

  !----------------------------------------------------------------------------
  ! Nullifies the pointers in the PostInfo variable
  !----------------------------------------------------------------------------
  Subroutine post_nullify(postobj)
    Type(PostInfo), Intent(InOut) :: postobj

    Nullify(postobj%nrgs)
    Nullify(postobj%loading)
    Nullify(postobj%siting)
    Nullify(postobj%poremap)
    Nullify(postobj%structure)
    Nullify(postobj%MSD)
    Nullify(postobj%vels)
    Nullify(postobj%iact)
    Nullify(postobj%radDist)

  End Subroutine post_nullify

  !----------------------------------------------------------------------------
  ! Generates the names of the configuration files given the base name and
  ! the first and last configuration file numbers
  !----------------------------------------------------------------------------
  Subroutine post_makeConfigList(pparams,baseFile,first,last)
    Type(PostInfo), Intent(InOut)     :: pparams
    Character(len=strLen), Intent(In) :: baseFile
    Integer, Intent(In) :: first,last
    Integer :: error, i

    !** Allocate space to save the configuration filenames
    Allocate(pparams%configFiles(last-first+1), stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate pparams%configFiles array of size ", &
          last-first+1
      Stop
    End If

    !** Generate the configuration file names given the base name 
    !** and the first number. Store them for later.
    Do i = 1, Size(pparams%configFiles)
      pparams%configFiles(i)%name = genfilename(baseFile,first+i-1)
    End Do
  End Subroutine post_makeConfigList

  !----------------------------------------------------------------------------
  ! Extracts the original control file from the given configuration file
  !----------------------------------------------------------------------------
  Subroutine post_extractControlFile(configFile,newCtrlFile)
    Type(CONFILE), Intent(InOut) :: configFile
    Character(len=strLen), Intent(In) :: newCtrlFile
    Logical :: oldflag

    ! determine whther music-2-2 *.con.* file
    oldflag= configFile%datafile_is_old

    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
        " Reading the header of datafile"
    !** First, open the configuration file for reading    
    Call datafile_initin(configFile,configFile%name, oldflag )

    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
        " writing re-generated ctrlfile"
    !** Generate the new control file from the config file
    Call datafile_gencontrolfile(newCtrlFile, configFile)

    !** Set the default control file tag, which may be useful elsewhere
    Call file_settag(newCtrlFile,d_ctrl_file)

    !** Close the file
    Close(configFile%unitno)

  End Subroutine post_extractControlFile

  !----------------------------------------------------------------------------
  ! Reads in the original control file to initialize the atoms, molecules, 
  ! simulation cell, configuration, and forcefield information
  ! - config_init_source tells how to initialize the configuration
  !   "DEFAULT" Implies a gcmc type initialization
  !   "CONTROLFILE" Implies use of control file 
  !----------------------------------------------------------------------------
  Subroutine post_readControlFile(postobj,newCtrlFile,config_init_source)
    Type(PostInfo), Intent(InOut) :: postobj
    Character(*), Intent(In) :: newCtrlFile
    Character(*), Intent(In) :: config_init_source
    Integer                  :: nmtypes, error

    !** Initialize the general parameters
    Call general_init(newCtrlfile)
    Write(*,*) "Initialized general section from new ctrlfile "

    !** Initialize atoms
    Call atom_init()
    Write(*,*) "Initialized atoms section from new ctrlfile "

    !** Initialize the molecules
    nmtypes = molecules_init()
    Write(*,*) "Initialized molecules section from new ctrlfile "

    !** Allocate the species array
    Allocate(postobj%species(nmtypes),stat=error)
    If (error /= 0) Then
      Write(0,'(a,1x,i5,a,i5)') __FILE__,__LINE__, &
          " : Could not allocate species array of size ",nmtypes
      Stop
    End If

    !** Initialize the simulation cell
    Call simcell_init(postobj%scell, newCtrlFile)
    Write(*,*) "Initialized simcell section from new ctrlfile "

    !** Initialize the configuration    
    If (Trim(config_init_source)=="DEFAULT") Then
      Call post_dfltConfigInit(postobj%species, postobj%scell, newCtrlfile)
      Write(*,*) "Initialized config section from new ctrlfile, made space &
          & for  default number of molecules "

      !** If there are restartfiles used for config init , then we have to 
      !** make sure that they are still kept there, Otherwise below 
      !** initializations could lead to wrong calculations
    Elseif(Trim(config_init_source)=="CONTROLFILE") Then
      Call config_init(postobj%species, postobj%scell, newCtrlFile)
      Write(*,*) "Initialized config section from new ctrlfile "
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Wrong config_init_source specified"
    End If

  End Subroutine Post_readControlFile

  !--------------------------------------------------------------------
  ! Initializes the sorbates data structure 
  ! This routine is for doing post code of gcmc type of moves
  ! If you use the usual config_init routine there is a danger that 
  ! you will not have enough space to read the first data point in 
  ! the config file
  !--------------------------------------------------------------------
  Subroutine post_dfltConfigInit(sorbates,scell, ctrlfile )
    Type(AtMolCoords), Dimension(:), Pointer :: sorbates
    Type(Simcell_Params),  Intent(InOut)     :: scell 
    Character(len=strLen), Intent(in)        :: ctrlfile
    Character(len=strLen)  :: srchstr, sourcetype, filename, tag, sorbname
    Integer    :: nsorbates
    Integer    :: sorb
    Integer    :: unitno, lineno, sorbsread, error, ios, i, j

    !** Open the ctrl_file if it is not opened
    unitno = file_open(ctrlfile)

    tag = default_config_tag
    !** Find the configuration intialization section
    lineno = filesrchstr(unitno, tag, srchstr)
    If (lineno == 0) Then
      Write(0,'(1x, 4a)') 'Could not find the tag "', Trim(tag), &
          '" in the control file ', d_ctrl_file
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

    !** Allocate the sorbates structure
    nsorbates = molecules_getnsorbs()  ! number of molecule types
    Allocate(sorbates(nsorbates), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__," sorbates")

    Do i=1, nsorbates
      sorbates(i)%isinit = .False.
    End Do

    !** Read the sorbate type and the source information
    Do i = 1, nsorbates
      sorbsread = 0
      Read(unitno, *, IOSTAT=ios) sorbname

      !** If ios /= 0 the we have not found at least one of the sorbates 
      !** though we don't know which one.  So we loop through to find 
      !** which one

      If (ios /= 0) Then
        Do j=1, nsorbates
          sorbname = molecules_name(j)
          If (.Not. sorbates(j)%isinit) Then
            Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
                " Could not find initial configuration information for ", &
                Trim(sorbname)
            Stop
          End If
        End Do
      End If

      sorb = molecules_gettype(sorbname)
      If (sorb==0) Then
        Write(*,*) " Wrong sorbname specified in the config section of the &
            &control file, sorbname -",Trim(sorbname)
      End If

      !** Set the initialization flag
      sorbates(sorb)%isinit = .True.

      !** Set the fixed coordinate flag to false
      sorbates(sorb)%fixed = .False.

      !** Get the information from one of many sources      
      Read(unitno, *, IOSTAT=ios) sourcetype, filename
      sorbates(sorb)%sourcetype = sourcetype
      sorbates(sorb)%filename = filename

      !** Now go and get the initial configuration
      !** This is a dummy of config_init in config.f90
      !** we would prefer not to use restart files etc for post code
      Select Case (toupper(Trim(sourcetype)))
      Case ("EQUILFILE")
        Call config_initdefaultsize(sorbates, sorb)
      Case ("RESTARTFILE")
        Call config_initdefaultsize(sorbates, sorb)
      Case ("CRASHFILE")
        Call config_initdefaultsize(sorbates, sorb)
      Case ("MEPRESTARTFILE")
        Call config_initdefaultsize(sorbates, sorb)
      Case ("CUSTOM")
        Call config_initdefaultsize(sorbates, sorb)
      Case ("GCMC")
        Call config_initdefaultsize(sorbates, sorb)
      Case ("DFLT_SIZE")
        Call config_initdefaultsize(sorbates, sorb)
      Case ("MOLECULE")
        !Just read the coordinates from molecule
        Call config_initmolecule(sorbates, sorb)
      Case ("FIXED")
        !** In this case, the coordinates are fixed and the molecule
        !** is never moved
        Call config_initFixed(sorbates,sorb)

        Cycle

      Case ("NOCONFIG")
        !There is no configuration for this molecule. Skip it.
        !This will probably only be used for the zeolite.
        sorbates(sorb)%natoms = 0
        sorbates(sorb)%nmoles = 0
        Nullify(sorbates(sorb)%aslow)
        Nullify(sorbates(sorb)%afast)
        Nullify(sorbates(sorb)%coords)
        Nullify(sorbates(sorb)%gcoords)
        !** Indicate that this sorbate is fixed in space, just
        !** in case any old control files are used
        sorbates(sorb)%fixed = .True.
        Cycle
      Case Default
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(1x,2a)')  &
            "Don't know how to initialize configurations for ", sourcetype
        Write(0,'(1x,2a)') &
            "Options are: EQUILFILE RESTARTFILE GCMC MOLECULE CCBGCMC NOCONFIG" 
        Stop
      End Select
    End Do
    Close(unitno)


  End Subroutine post_dfltConfigInit


  !----------------------------------------------------------------------------
  ! Returns the number of different configuration files whose names are
  ! stored in pparams%configFiles
  !----------------------------------------------------------------------------
  Integer Function post_getnfiles(pparams)
    Type(PostInfo), Intent(In) :: pparams

    post_getnfiles = Size(pparams%configFiles,1)
  End Function post_getnfiles

  !----------------------------------------------------------------------------
  ! Returns the configuration filename corresponding to the given number
  ! in the pparams%configFiles list
  !----------------------------------------------------------------------------
  Function post_getConfigFilename(pparams,fileno)
    Type(PostInfo), Intent(In) :: pparams
    Integer, Intent(In) :: fileno
    Character(len=strLen) :: post_getConfigFilename

    post_getConfigFilename = pparams%configFiles(fileno)%name
  End Function post_getConfigFilename

  !----------------------------------------------------------------------------
  ! Displays the parameters stored in pparams
  !----------------------------------------------------------------------------
  Subroutine post_display(pparams,unitno)
    Type(PostInfo), Intent(In) :: pparams
    Integer, Intent(In) :: unitno
    Integer :: i,nsites
    Integer, Dimension(:), Pointer :: allsites
    Write(unitno,'(a)') dashedline
    Write(unitno,'(a)') "Post Processor General Information"
    Write(unitno,'(a)') dashedline
    Write(unitno,'(2a)') "Configuration files base filename   : ", &
        trim(pparams%baseFilename)
    Write(unitno,'(2a)') "Configuration files to analyze      : ", &
        trim(pparams%configFiles(1)%name)
    Do i = 2, Size(pparams%configFiles)
      Write(unitno,'(38x,a)') Trim(pparams%configFiles(i)%name)
    End Do
    Write(unitno,'(2a)') "Reconstructed control filename      : ", &
        trim(pparams%newCtrlFile)
    Write(unitno,'(2a)') "Post processor output base filename : ", &
        trim(pparams%outFilename)
    Write(unitno,'(a)')
    Write(unitno,'(a)') "The following analysis will be done,"
    Write(unitno,'(a)') "modify post code ctrl file , &
        &if you dont need some of them"
    If (associated(pparams%nrgs)) &
        Write(unitno,'(a)') "ENERGY STATISTICS : "

    If (Associated(pparams%loading)) &
        Write(unitno,'(a)') "LOADING STATISTICS (ISOTHERMS) : "
    If (associated(pparams%siting)) Then
      Write(unitno,'(a)') "SITING STATISTICS : "
      Call smap_getallsitetypes(pparams%siting%sitemap, allsites, nsites)
      Write(*,'(1x,a,i4,a)') "original sitemap contains: ", nsites, "sites"
      Write(unitno,*) "they are : ", allsites
      Write(unitno,*) "check that you have correct sites in postcode ctrlfile"
      Call smap_display(pparams%siting%sitemap, 6, 6)
    Endif
    If (Associated(pparams%poremap)) &
        Write(unitno,'(a)') "MAKE POREMAP BASED ON PROBABILITIES : "
    If (Associated(pparams%structure)) &
        Write(unitno,'(a)') "STRETCHING, ANGLE, TORSION DISTRIBUTIONS"
    If (Associated(pparams%MSD)) Then
      Write(unitno,'(a)') "DIFFUSIVITY MSD"
      ! Next flag added by simon 1/31/04
      If (Associated(pparams%MSD%mdpp%radial)) &
          Write(unitno,'(a)') "RADIAL PROFILE"
    End If
    If (Associated(pparams%radDist))Then
      Write(unitno,'(a)') "RADIAL DISTRIBUTION CALCULATIONS"
      Write(*,'(1x,a,2f12.2)') "Minimum  and Maximum Distance, Ang. : ", &
          pparams%radDist%mindist, pparams%radDist%maxdist
      Write(*,'(1x,a,i5)') "Number of bins : ",pparams%radDist%nbins
      Write(*,'(1x,a,f12.2)') "Simcell Volume (Ang^3) : ",pparams%radDist%vol
      If (pparams%radDist%com_based(1)) Then
        Write(*,'(1x,a)') "Center of Mass of Species-1 willbe used &
            &for radial distribution functions"
      Else
        Write(*,'(1x,a,i4,a)') "Atom number :",pparams%radDist%atoms(1), &
            " Species-1 willbe used for radial distribution functions"
      Endif
      If (pparams%radDist%com_based(2)) Then
        Write(*,'(1x,a)') "Center of Mass of Species-2 willbe used &
            &for radial distribution functions"
      Else
        Write(*,'(1x,a,i4,a)') "Atom number :",pparams%radDist%atoms(2), &
            " of Mass of Species-2 willbe used &
            &for radial distribution functions"
      Endif
      Write(*,'(1x,3a,i3)') "SPC(1) : ",&
          Trim(pparams%radDist%sorbnames(1)), "Type :",pparams%radDist%sorbs(1)
      Write(*,'(1x,3a,i3)') "SPC(2) : ",&
          Trim(pparams%radDist%sorbnames(2)), "Type :",pparams%radDist%sorbs(2)
    Endif

    If (Associated(pparams%vels)) &
        Write(unitno,'(a)') "VELOCITY DISTRIBUTIONS"

    If (Associated(pparams%iact)) Then
      Write(unitno,'(a)') "ENERGY CALCULATIONS (params given below) "
      Call interact_display(pparams%iact%imodel,2,6)
    Endif

    Write(unitno,'(a)') "---------- End of post code params --------- "    
    Write(unitno,'(a)')
  End Subroutine post_display


  !----------------------------------------------------------------------------
  ! Writes information about what it should read in from a control file
  ! during post_init
  !----------------------------------------------------------------------------
  Subroutine post_sampleCF(unitno)
    Integer, Intent(In) :: unitno


    !** General section
    Write(unitno, '(a30)') '### Required section ######'
    Write(unitno,'(a30)') "-- "//Trim(postIdString)//" ---"
    Write(unitno,'(a,t30,a)') 'Character','# Type of simulation, GCMC, MD, &
        & NVTMC'
    Write(unitno,'(a,t30,a)') 'Character','# basename for config files'
    Write(unitno,'(a,t30,a)') 'Integer, Integer',&
        '# first and last file numbers'
    Write(unitno,'(a,t30,a)') 'Character','# new ctrlfile regenerated'
    Write(unitno,'(a,t30,a)') 'Character','# Base name for output files '
    Write(unitno,'(a,t30,a)') 'Real, Real',&
        '# Percentages of data to skipped at start and end '

    !** Sections that will work with all type of simulations
    Write(unitno, '(a)') '### Optional sections ######'
    Write(unitno, '(a)') '### Delete unwanted/unsuitable sections ######'

    Write(unitno,'(a)') "-- "//Trim(postNrgAvgIdString)//" ---"
    Write(unitno,'(a,t30,a)') 'Integer','# Number of blocks for stats' 

    Write(unitno,'(a)') "-- "//Trim(postLoadingIdString)//" ---"
    Write(unitno,'(a,t30,a)') 'Integer','# Number of blocks for stats' 

    Write(unitno,'(a)') "-- "//Trim(postSitingIdString)//" ---"
    Write(unitno,'(a,t30,a)') 'Integer','# Number of blocks for stats' 
    Write(unitno,'(a,t30,a)') 'Character','# Name of sitemap'
    Write(unitno,'(a,t30,a)') 'Integers','# Indices of sites, comma seperated'

    Write(unitno,'(a)') "-- "//Trim(postPoreMapIdString)//" ---"
    Write(unitno,'(a,t30,a)') 'Integer','# Number of poremaps ' 
    Write(unitno,'(a,t30,a)') 'Integer(3)','# Number of divisions in x,y,z dir'
    Write(unitno,'(a,t30,a)') 'Character','# Name of sorbate for 1st map'
    Write(unitno,'(a,t30,a)') 'Character','# Name of 1st map'
    Write(unitno,'(a)') '..Add more sorbate names, and poremap names if &
        &you have more than one sorbate of interest..'
    Write(unitno,*) 

    Write(unitno,'(a)') "-- "//Trim(postStructureIdString)//" ---"

    !** sections that will work only with specific simulation types
    Write(unitno, '(a)') '### Below sections might not work for all &
        &type of sims ####'
    Write(unitno, '(a)') '### Delete unwanted/unsuitable sections ######'
    Call mdpc_sampleCF(unitno)

  End Subroutine post_sampleCF



  !----------------------------------------------------------------------------
  ! Some of the analysis types will be initialized here, These are the ones 
  ! which are used for ensemble averaging, making maps etc..
  ! 
  !----------------------------------------------------------------------------
  Subroutine post_initAnalTypes(pparams,controlFile)
    Type(PostInfo), Intent(InOut)     :: pparams
    Character(len=strLen), Intent(In) :: controlFile
    Character(len=strLen)             :: tag,text 
    Integer :: unitno, lineno

    Logical :: rewind_flag = .True.
    unitno=file_open(controlFile)

    !** Find the analysis types' flag in the control file , 
    !** If found initialize them

    tag=postNrgAvgIdString
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno/=0) Call post_initNrgAvg(pparams%nrgs,unitno)

    tag=postLoadingIdString
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno/=0) Call post_initLoading(pparams, pparams%loading,unitno)

    tag=postSitingIdString
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno/=0) Call post_initSiting(pparams,pparams%siting, &
        pparams%scell,unitno)

    tag=postPoreMapIdString
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno/=0) Call post_initPoreMap(pparams,pparams%poremap,unitno)

    tag=postStructureIdString
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno/=0) Call post_initStructure(pparams,pparams%structure, &
        unitno)

    tag=postRadDistIdString
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno/=0) Call post_initRadDist(pparams,pparams%RadDist, &
        unitno)

    tag=postMSDIdString
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno/=0) Call post_initMSD(pparams,controlFile)

    tag=postVelHistIdString
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno/=0) Call post_initVelHist(pparams,pparams%vels,unitno)

    tag=postInteractString
    lineno = filesrchstr(unitno,tag,text,rewind_flag)
    If (lineno/=0) Call post_initInteract(pparams,pparams%iact,unitno)



  End Subroutine post_initAnalTypes




  !----------------------------------------------------------------------------
  ! Read the required info for the nrg averaging 
  !----------------------------------------------------------------------------
  Subroutine post_initNrgAvg(nrgs,unitno)
    Type(nrgAvg), Pointer :: nrgs
    Integer, Intent(In) :: unitno
    Integer           :: error, nsorbates, i

    Allocate(nrgs,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"nrg-params")

    ! block size used for statistics
    Read(unitno,*) nrgs%nblocks
    nsorbates=molecules_getnsorbs()

    Allocate(nrgs%coulNrg(nsorbates, nsorbates), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"coulNrg")

    Allocate(nrgs%noncoulNrg(nsorbates, nsorbates), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"noncoulNrg")

    Allocate(nrgs%intra(nsorbates,NO_OF_INTRA_POTS), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"intraNrg")

    Allocate(nrgs%ke(nsorbates), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"KinNrg")


    Allocate(nrgs%output_array(nrgs%nblocks, 4), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"output array")

  End Subroutine post_initNrgAvg

  !----------------------------------------------------------------------------
  ! Read the required info for the loading averaging and allocate memory
  !----------------------------------------------------------------------------
  Subroutine post_initLoading(pparams,loading,unitno)
    Type(PostInfo), Intent(InOut)    :: pparams
    Type(loadingAvg), Pointer :: loading
    Integer, Intent(In) :: unitno
    Integer           :: error, nsorbates, nsims

    Allocate(loading,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"loading-params")

    ! block size used for statistics
    Read(unitno,*) loading%nblocks

    nsorbates=molecules_getnsorbs()
    Allocate(loading%totalavgs(nsorbates ), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"loading-avg")

    Allocate(loading%output_array(nsorbates, loading%nblocks, 4), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"output array")

    ! the dimesnions are sims, nsorbates, 2(for pressure and loading)
    nsims=pparams%last - pparams%first + 1
    Allocate(loading%isotherm(nsims,nsorbates, 2), STAT=error)



  End Subroutine post_initLoading



  !----------------------------------------------------------------------------
  ! Read the required info for the velocity histograms 
  !----------------------------------------------------------------------------
  Subroutine post_initVelHist(pparams,vels,unitno)
    Type(PostInfo), Intent(InOut)    :: pparams
    Type(velocityParams), Pointer         :: vels 
    Integer, Intent(In) :: unitno

    Integer           :: error, nsorbates, spctype, natoms, nfields, ndivs, i
    Real(kind=RDbl)   :: minVel, maxVel
    Character(len=strLen) :: text, spcname, tempstring
    Character(len=strLen),Dimension(strLen) :: fields 

    Allocate(vels,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"velhist-params")

    ! block size used for statistics
    Read(unitno,'(a)') text
    spcname=cleanstring(text)
    spctype=molecules_gettype(spcname)
    If (spctype==0)  Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "wrong molecule in ctrlf file for vel histogram"
      Stop
    Endif
    natoms=molecules_getnatoms(spctype)
    vels%natoms=natoms
    vels%spctype=spctype
    vels%spcname=spcname

    Allocate(vels%atomcomp_hist(natoms,3), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"velocity-components")


    Allocate(vels%atom_hist(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"atom-velocity")

    Allocate(vels%com_hist, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"com-velocity")


    ! this atom list is read for display purpose all hsitograms will 
    ! not be displayed
    Read(unitno,'(a)') text
    nfields=split(cleanstring(text), fields, ",")
    Allocate(vels%display_list(nfields), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"com-velocity")
    Do i=1,nfields
      vels%display_list(i)=toint(fields(i))
    End Do

    Read(unitno,*) minVel, maxVel
    Read(unitno,*) ndivs
    vels%minVel=minVel
    vels%maxVel=maxVel
    vels%ndivs=ndivs

    Do i=1,natoms
      tempstring="xvelocity-histogram"
      Call histogram_init(vels%atomcomp_hist(i,1),ndivs,minVel,maxVel,&
          tempstring)
      tempstring="yvelocity-histogram"
      Call histogram_init(vels%atomcomp_hist(i,2),ndivs,minVel,maxVel,&
          tempstring)
      tempstring="zvelocity-histogram"
      Call histogram_init(vels%atomcomp_hist(i,3),ndivs,minVel,maxVel,&
          tempstring)
      tempstring="velocity-histogram"
      Call histogram_init(vels%atom_hist(i),ndivs,minVel,maxVel,&
          tempstring)
    End Do
    tempstring="com velocity-histogram"
    Call histogram_init(vels%com_hist,ndivs,minVel,maxVel,&
        tempstring)
  End Subroutine post_initVelHist



  !----------------------------------------------------------------------------
  ! Read the required info for interaction calculation
  !----------------------------------------------------------------------------
  Subroutine post_initInteract(pparams,iact,unitno)
    Type(PostInfo), Intent(InOut)    :: pparams
    Type(interactParams), Pointer         :: iact 
    Integer, Intent(In) :: unitno

    Integer           :: error, nsorbates, spctype, natoms, nfields, ndivs, i
    Real(kind=RDbl)   :: minVel, maxVel
    Character(len=strLen) :: c_file, text, spcname, tempstring
    Character(len=strLen),Dimension(strLen) :: fields 

    Allocate(iact,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"interact-params")

    ! ok for now we are just going to interact as for a normal program
    ! later we might add more fields here
    iact%nimodels=1

    nsorbates=molecules_getnsorbs()

    ! imodel needs ctrlfile_name
    c_file=file_getname(unitno)

    Write(*,'(1x,2a,i4,a)') __FILE__," : ",__LINE__," initializing interact" 
    Call interact_init(iact%imodel, c_file, pparams%scell, &
        pparams%species, 'MC' , postInteractString)

    ! arrays for storing mol-sys nrgs
    Allocate(iact%coul(nsorbates,post_MAX_MOLEC,nsorbates), &
        iact%ncoul(nsorbates,post_MAX_MOLEC,nsorbates),&
        iact%intra(nsorbates,post_MAX_MOLEC),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"interact-params")
    iact%coul=zero
    iact%ncoul=zero
    iact%intra=zero
    ! arrays for getting mol-sys nrgs from interact
    Allocate(iact%molsys_nrg(nsorbates),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"interact-params")

    Do i=1,nsorbates
      If (config_isfixed(pparams%species(i))) Cycle
      ! now we need to shape these molsys_nrgs properly, copy it from storage
      ! (/i,1,0/)  - implies for one molecule of species i
      Call storetop_initcopy(iact%molsys_nrg(i) ,iact%imodel%results(1), &
          (/i,1,0/) )
      ! do some copying stuff -may be not reqd ?
      Call storetop_fillsub(iact%molsys_nrg(i) , (/i,1,0/) )
      Call storetop_setmap(iact%molsys_nrg(i),.True., (/i,1,0/) )
    End Do


  End Subroutine post_initInteract

  !----------------------------------------------------------------------------
  ! Read the required info for the siting averaging and allocate memory
  !----------------------------------------------------------------------------
  Subroutine post_initSiting(pparams,siting,scell,unitno)
    Type(PostInfo), Intent(InOut)    :: pparams
    Type(sitingAvg), Pointer         :: siting
    Type(Simcell_Params), Intent(In) :: scell

    Integer, Intent(In) :: unitno
    Integer           :: error, nsorbates, s, sorb, nsites, nsims, atom_num
    Integer, Dimension(:), Pointer :: allsites
    Logical :: hasposition

    Character(len=strLen) :: line, sitemapname,filename
    Character(len=strLen), Dimension(strLen) :: fields 

    Allocate(siting,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"siting-params")

    ! block size used for statistics
    Read(unitno,*) siting%nblocks 
    Read(unitno,*) sitemapname
    sitemapname=Trim((stripcmnt(sitemapname)))
    nsorbates=molecules_getnsorbs()
    hasposition=datafile_hasposition(pparams%configFiles(1))
    If (.Not.hasposition) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "config file does not contain positions &
          &so no sitemap info can be generated -Sorry"
      Stop
    End If

    !** Read the index number of sites to watch for
    Read(unitno,'(a)') line
    line=cleanstring(line)
    siting%nsites = split(line, fields, ",")
    Allocate(siting%sites(siting%nsites),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Allocate(siting%siteavgs(siting%nsites,nsorbates),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Allocate(siting%output_array(siting%nsites, nsorbates, &
        siting%nblocks, 4), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"output array")

    Allocate(siting%site_inst(siting%nsites,nsorbates),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Do s=1,siting%nsites
      siting%sites(s)=toint(fields(s))
    End Do

    siting%sitemapfile=sitemapname
    Call smap_init(siting%sitemap,pparams%scell,sitemapname)
    Call smap_getallsitetypes(siting%sitemap, allsites, nsites)
    Write(*,'(a,i4,a)') "original sitemap contains : ", nsites, " sites"
    Write(*,*) "they are ", allsites
    Write(*,*) "check whether you have correct sites in post code ctrlfile"

    Call smap_display(siting%sitemap, 6, 6)

    !
    Read(unitno,'(a)') line
    line=cleanstring(line)
    If (Trim(line)=="COM") Then
      siting%com_based=.True.
    Else
      siting%com_based=.False.
      atom_num=toint(line,error)
      siting%atom_num=atom_num
      If (error/=0) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) " Error in smap section, atom_num is specified worng"
        Write(*,*) "Line used is : ",Trim(line)
        Stop
      Endif
    Endif

    Read(unitno,'(a)') line
    line=cleanstring(line)
    If (Trim(line)=="DISPLAY_SMAP") Then
      Write(*,*) " Displaying the smap as a  rasmol plottable file "
      Write(*,*) " Use rasmol on file smap.xyz to look at various sites "
      filename='smap.xyz'
      Call smap_showxyz(siting%sitemap,filename,allsites(1:nsites))
    Endif

    !allocate isotherm array
    ! the dimesnions are sims, nsorbates, 2+nsites 
    !                                     (pressure, total, site-loading)
    nsims=pparams%last - pparams%first + 1
    Allocate(siting%isotherm(nsims,nsorbates, 2+siting%nsites), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__) 

  End Subroutine post_initSiting

  !----------------------------------------------------------------------------
  ! Read the required info for making a pore map and allocate memory
  ! If you dont know what is a pore-map see section marked "POREMAP_DETAILS" 
  ! at the End of this file
  !----------------------------------------------------------------------------
  Subroutine post_initPoremap(pparams,poremap, unitno)
    Type(PostInfo), Intent(InOut) :: pparams
    Type(poremapParams), Pointer :: poremap
    Integer, Intent(In) :: unitno
    Integer           :: error, i, nfields, natoms, sorb, tempatomno
    Integer,Dimension(MAX_ATOMS)  :: int_array

    Character(len=strlen) :: text
    Character(len=strlen),dimension(strLen) :: fields 
    Logical :: hasposition

    hasposition=datafile_hasposition(pparams%configFiles(1))
    If (.Not.hasposition) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "config file does not contain positions &
          &so no pore map info can be generated -Sorry"
      Stop
    End If

    Allocate(poremap,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"poremap-params")
    Write(*,*) "Initializing porempa/dmap generation section"

    ! number of dmaps to be made , one for each sorbate
    Read(unitno,*) poremap%nmaps

    Allocate(poremap%sorbnames(poremap%nmaps), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"poremap")

    Allocate(poremap%mapname(poremap%nmaps), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"poremap")

    Read(unitno,*) poremap%nbins(1),poremap%nbins(2),poremap%nbins(3)
    Allocate(poremap%fillArray(poremap%nmaps,poremap%nbins(1), &
        poremap%nbins(2), poremap%nbins(3)), STAT=error)
    poremap%fillArray=0

    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__," fill array")
    ! number of unitcells in each direction
    Read(unitno,*) poremap%uc_repeat_info 

    Do i=1,poremap%nmaps
      Read(unitno,*) text
      poremap%sorbnames(i)=Trim(stripcmnt(text))
      sorb=molecules_gettype(Trim(poremap%sorbnames(i)))
      If (sorb<1) Then
        Write(*,*) "Wrong name for sorbate : "//&
            Trim(poremap%sorbnames(i))//" - specified in poremap section"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End If
      Read(unitno,*) text
      poremap%mapname(i)=Trim(stripcmnt(text))
    End Do

    Read(unitno,*) poremap%ncont
    Read(unitno,*) poremap%scaling_factor
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    poremap%atom_list_given=.False.
    Read(unitno,'(a)') text
    text=cleanstring(text)
    Select Case(Trim(text))
    Case("COM")
      poremap%onlycom=.True.
    Case ("ALLATOMS")
      poremap%onlycom=.False.
    Case default
      ! see whether we have an atom list
      natoms=str2seq(text,int_array)

      If (natoms==0) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Error in poremap-post code section, expecting strings :"
        Write(*,*) "--COM-- OR --ALLATOMS--, but found : -"//Trim(text)//"--"
        Stop
      Else
        poremap%natoms=natoms
        poremap%onlycom=.False.
        poremap%atom_list_given=.True.
      Endif
    End Select

    If (poremap%atom_list_given) Then
      If ((poremap%nmaps)>1) Then
        Write(*,*) "lazy to code... if u want individual atom poremap/dmap"
        Write(*,*) "then can have only one poremap"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      Allocate(poremap%atom_list(natoms),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      poremap%atom_list(1:natoms)=int_array(1:natoms)
    Endif

    Read(unitno,*) text
    text=cleanstring(text)
    Select Case(Trim(text))
    Case("LOGARITHMIC")
      poremap%cont_div_exponential=.True.
    Case ("ARITHMETIC") 
      poremap%cont_div_exponential=.False.
    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Error initing poremap-post code section, expecting strings :"
      Write(*,*) "--LOGARITHMIC-- OR --ARITHMETC--, but found : --"&
          //Trim(text)//"--"
      Stop
    End Select

    Read(unitno,*) text
    text=cleanstring(text)
    Select Case(Trim(text))
    Case("DUMP_DMAP_YES")
      poremap%dump_dmap=.True.
    Case ("DUMP_DMAP_NO") 
      poremap%dump_dmap=.False.
    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Error initing poremap-post code section, expecting strings :"
      Write(*,*) "--DUMP_DMAP_YES-- OR --DUMP_DMAP_NO--, but found : --"&
          //Trim(text)//"--"
      Stop
    End Select

    poremap%uc_edge_l=simcell_getell(pparams%scell, .True.)
    poremap%dx=poremap%uc_edge_l(1)/poremap%nbins(1)
    poremap%dy=poremap%uc_edge_l(2)/poremap%nbins(2)
    poremap%dz=poremap%uc_edge_l(3)/poremap%nbins(3)

    Write(*,*) "Finshed initializing poremap/dmap generation section"
  End Subroutine post_initPoremap

  !----------------------------------------------------------------------------
  ! Read the required info for doing averages of molecular
  !structure. These include Bond angle averaging, Bond length averaging,
  !Torsion and any other intramolecular variable that you can think of
  !----------------------------------------------------------------------------
  Subroutine post_initStructure(pparams, structure, unitno)
    Type(PostInfo), Intent(InOut) :: pparams
    Type(structureAvg), Pointer :: structure
    Integer, Intent(In) :: unitno
    Integer           :: error, i, j, nfields, natoms
    Character(len=strlen) :: text, sorbname, temptag
    Character(len=strlen),Dimension(strLen) :: fields
    Logical :: hasposition

    hasposition=datafile_hasposition(pparams%configFiles(1))

    If (.Not.hasposition) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "config file does not contain positions &
          & so no structure info can be generated -Sorry"
      Stop
    End If



    Allocate(structure,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"structure-params")

    ! number of structure analysis to be made 
    Read(unitno,*) structure%nparams

    Allocate(structure%hist(structure%nparams), &
        structure%minval(structure%nparams), &
        structure%maxval(structure%nparams), &
        structure%ndivs(structure%nparams), &
        STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"structure")

    Allocate(structure%sorbnames(structure%nparams), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"structure")

    Allocate(structure%intra_type(structure%nparams), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"structure")

    ! assume that maximum number of atoms required for intra-molecular 
    ! characteristic is 10. For example bond-angles requires 3 atoms, 
    ! torsion requires 4 atoms ------HACK-----
    Allocate(structure%atom_list(structure%nparams,10), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"structure")

    !** Read thru individual structure averages
    Do i=1,structure%nparams
      Read(unitno,'(a)') text
      text=cleanstring(text)
      nfields=split(text, fields, ",")
      Write(*,*) nfields, text 
      If (nfields<5) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      If (fields(1)=="ANGLES") Then
        sorbname=Trim(fields(2))
        structure%intra_type(i)=fields(1)
        If (molecules_exists(sorbname)) Then
          structure%sorbnames(i)=sorbname
        Else
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "Wrongsorbate name in this angle avg section"
          Stop
        Endif
        structure%minval(i)=toreal(fields(3))
        structure%maxval(i)=toreal(fields(4))
        structure%ndivs(i)=toint(fields(5))
        temptag=fields(1)
        Call histogram_init( structure%hist(i), structure%ndivs(i), &
            structure%minval(i), structure%maxval(i), temptag)
        !** Read the line containing angles-list
        Read(unitno,*) text
        text=stripcmnt(text)
        natoms=split(text,fields,"-")
        If (natoms/=3) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) 'natoms', natoms
          Write(*,*) "something wrong in ctrlfile ?"
          stop
        End If
        Do j=1,3
          structure%atom_list(i,j)=toint(fields(j))
        End Do

      Elseif(fields(1)=="BONDLENGTHS") Then

        sorbname=Trim(fields(2))
        structure%intra_type(i)=fields(1)
        If (molecules_exists(sorbname)) Then
          structure%sorbnames(i)=sorbname
        Else
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "sorbate name in this length avg section"
          Stop
        Endif

        structure%minval(i)=toreal(fields(3))
        structure%maxval(i)=toreal(fields(4))
        structure%ndivs(i)=toint(fields(5))
        temptag=fields(1)
        Call histogram_init( structure%hist(i), structure%ndivs(i), &
            structure%minval(i), structure%maxval(i), temptag)

        !** Read the line containing bondlength-list
        Read(unitno,*) text
        text=stripcmnt(text)
        natoms=split(text,fields,"-")
        If (natoms/=2) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) 'natoms', natoms
          Write(*,*) "something wrong in ctrlfile ?"
          Stop
        End If
        Do j=1,2
          structure%atom_list(i,j)=toint(fields(j))
        End Do

      Elseif(fields(1)=="TORSIONANGLE") Then

        sorbname=Trim(fields(2))
        structure%intra_type(i)=fields(1)
        tEmptag=fields(1)
        If (molecules_exists(sorbname)) Then
          structure%sorbnames(i)=sorbname
        Else
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "sorbate name in this torsion-angle avg section"
          Stop
        Endif

        structure%minval(i)=toreal(fields(3))
        structure%maxval(i)=toreal(fields(4))
        structure%ndivs(i)=toint(fields(5))

        Call histogram_init( structure%hist(i), structure%ndivs(i), &
            structure%minval(i), structure%maxval(i), temptag)

        !** Read the line containing bondlength-list
        Read(unitno,*) text
        text=stripcmnt(text)
        natoms=split(text,fields,"-")
        If (natoms/=4) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) 'natoms', natoms
          Write(*,*) "something wrong in ctrlfile ?"
          Stop
        End If
        Do j=1,4
          structure%atom_list(i,j)=toint(fields(j))
        End Do

      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End If
    End Do

  End Subroutine post_initStructure


  !----------------------------------------------------------------------------
  ! Read the required info for doing averages of radial distribtuion functions 
  !----------------------------------------------------------------------------
  Subroutine post_initRadDist(pparams, radDist, unitno)
    Type(PostInfo), Intent(InOut) :: pparams
    Type(radDistAvg), Pointer :: radDist
    Integer, Intent(In) :: unitno
    Integer           :: error, i, nfields
    Character(len=strlen) :: text, temptag
    Character(len=strlen),Dimension(strLen) :: fields
    Logical :: hasposition

    hasposition=datafile_hasposition(pparams%configFiles(1))

    If (.Not.hasposition) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "config file does not contain positions &
          & so rad. dist. functions can not be generated -Sorry"
      Stop
    End If

    Allocate(radDist,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"rad.dist.-params")

    ! min, max and no ofbins 
    Read(unitno,*) radDist%mindist,  radDist%maxdist
    Read(unitno,*) radDist%nbins

    Do i=1,2
      Read(unitno,'(a)') text
      text=cleanstring(text)
      nfields=split(text, fields, ",")

      If (nfields/=2) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      radDist%sorbnames(i)=cleanstring(fields(1))
      fields(2)=cleanstring(fields(2))
      radDist%sorbs(i)=molecules_gettype(fields(1))
      If (  (radDist%sorbs(i)<1).Or.&
          (radDist%sorbs(i)>molecules_getnsorbs()) ) Then
        Write(*,*) "wrong sorbname in rad. dist. func. section"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      If (Trim(fields(2))=="COM") Then
        radDist%com_based(i)=.True.
        radDist%atoms(i)=0
      Else
        radDist%com_based(i)=.False.
        radDist%atoms(i)=toint(fields(2))
      Endif

      If (radDist%atoms(i)>molecules_getnatoms( radDist%sorbs(i))) Then
        Write(*,*) "wrong atom-number in rad. dist. func. section"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif


    end do


    temptag= "rad. dist. func."
    Call histogram_init( radDist%hist, radDist%nbins, &
        radDist%mindist, radDist%maxdist, temptag)

    raddist%vol=simcell_getvolume(pparams%scell)

  End Subroutine post_initRadDist


  !----------------------------------------------------------------------------
  ! Read the required info for the nrg averaging 
  !----------------------------------------------------------------------------
  Subroutine post_initMSD(pparams,postFile)
    Type(PostInfo), Intent(InOut) :: pparams
    Character(len=strLen), Intent(in) :: postFile

    Integer :: nfile, error, simno

    Write(*,*) " Initializing for MD Diffusivities calculation"

    Allocate(pparams%MSD, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__," diffusivity" )
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Call mdpc_init(pparams%MSD%mdpp, postFile, pparams%newCtrlFile, &
        pparams%species, pparams%outFileName )
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Call mdpc_display(pparams%MSD%mdpp,6)
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    nfile = post_getnfiles(pparams)

    !** We have the total number of config files to analyze and the
    !** number of sorbates in our system. Allocate space for the 
    !** scratch array, which stores the unit numbers and other useful
    !** information for the scratch files we use in the analysis
    Allocate(pparams%MSD%scratch(nfile), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"scratch-mdpc")
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    !** Loop through the config files we are to analyze and create a
    !** Fortran scratch file for each. Fortran scratch files are
    !** temporary files that only last as long as they are kept open
    !** by the running program. We use these scratch files to temporarily
    !** store information for processing.
    Do simno=1,nfile
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      !** Create the scratch files we need for futher analysis. The last
      !** option specifies what to write to the scratch files. We want
      !** diffusivites, so we write positions and time information
      Call mdpc_createScratchFiles(pparams%MSD%mdpp, &
          pparams%MSD%scratch(simno),&
          pparams%species, pparams%spcstats, &
          pparams%configFiles(simno), 101,0)
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      !** Update the user
      Write(0,'(a,i4,2a,i4)') __FILE__,__LINE__,": ", &
          "Created scratch files for sim ",simno
    End Do
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    !** Allocate space for the diffusivities. We have 4 diffusivities for
    !** each molecule type (at least).
    Allocate(pparams%MSD%dself(molecules_getnsorbs(),4),stat=error)
    If (error /= 0) Then
      Write(0,*) __FILE__,__LINE__,": Could not allocate dself array of &
          & size 2x4"
      Stop
    End If
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
  End Subroutine post_initMSD

  !------------------------------------------------------------------
  ! Re initialize all the statistics variables stored in all the fields
  ! odf pparams
  !------------------------------------------------------------------
  Subroutine post_calcMSD(postobj)
    Type(PostInfo), Intent(InOut) :: postobj

    !** Update the user
    Write(0,*) "Starting to Calculating diffusivities"

    !** Call mdpc to calculate the diffusivites. If we have more than
    !** one configuration file, it will average the diffusivity over
    !** all the configuration files if requested in order to improve
    !** the statistics. 
    Call mdpc_calcSelfDiffusivity(postobj%MSD%mdpp,postobj%MSD%scratch,&
        postobj%MSD%dself,.True.,.True.)
    Write(0,*) __FILE__,__LINE__,postobj%MSD%dself(:,4)

  End Subroutine post_calcMSD


  !---------------------------------------------------------------------
  ! Calculate the radial profiles
  ! Added by Simon 01/30/04
  !-------------------------------------------------------------------
  Subroutine post_calcRadProf(postobj)
    Type(PostInfo), Intent(InOut) :: postobj

    !currently not doing own files,
    ! so the scratch files have to be created when the diffusion 
    !calculations are done

    ! First we need to create the scratch files using the radial flag
    ! Call mdpc_createScratchFiles(postobj%MSD%mdpp,postobj%MSD%scratch,&
    !    postobj%species,postobj%spcstats, postobj%configFiles,101,1)

    !** Call mdpc to calculate the radial profiles
    Call mdpc_calcRadProf(postobj%MSD%mdpp,postobj%MSD%scratch(1))

  End Subroutine post_calcRadProf


  !------------------------------------------------------------------
  ! Re initialize all the statistics variables stored in all the fields
  ! odf pparams
  !------------------------------------------------------------------
  Subroutine post_reInitStats(pparams)
    Type(PostInfo), Intent(InOut) :: pparams
    Integer :: blocksize, nsorbs, ndats,  i, j

    !** Get the number of molecules
    nsorbs=molecules_getnsorbs()

    !** Initialize nrg statistics
    If (Associated(pparams%nrgs)) Then
      blocksize=Int((one*(pparams%lastconf-pparams%firstconf))/&
          pparams%nrgs%nblocks)+1
      If (blocksize==0) Then
        Write(*,*) " Too many number of blocks specified for nrg avgs"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End If
      pparams%nrgs%statblocksize=blocksize
      Call stats_init(pparams%nrgs%total, ' ', blocksize,.False., 'f15.3')

      Do i=1,nsorbs
        Do j=1, nsorbs
          Call stats_init(pparams%nrgs%coulNrg(i,j),' ',blocksize, &
              .False.,'f15.3')
          Call stats_init(pparams%nrgs%noncoulNrg(i,j),' ',blocksize,&
              .False.,'f15.3')
        End Do
        Call stats_init(pparams%nrgs%ke(i),' ',blocksize,.False.,'f15.3')
        Do j=1, NO_OF_INTRA_POTS
          Call stats_init(pparams%nrgs%intra(i,j),' ',blocksize, &
              .False., 'f15.3')
        End Do
      End Do
    End If

    !** Initialize loading statistics
    If (Associated(pparams%loading)) Then
      blocksize=Int((one*(pparams%lastconf-pparams%firstconf))/ &
          pparams%loading%nblocks)+1
      If (blocksize==0) Then
        Write(*,*) " Too many number of blocks specified for loading avgs"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End If
      pparams%loading%statblocksize=blocksize
      Do i=1,nsorbs
        Call stats_init(pparams%loading%totalAvgs(i),' ',blocksize, &
            .False.,'f15.3')
      End Do
    End If

    !** Initialize siting statistics
    If (Associated(pparams%siting)) Then
      blocksize=Int((one*(pparams%lastconf-pparams%firstconf))/&
          pparams%siting%nblocks)+1
      If (blocksize==0) Then
        Write(*,*) " Too many number of blocks specified for siting avgs"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End If
      pparams%siting%statblocksize=blocksize
      Do i=1,pparams%siting%nsites
        Do j=1,nsorbs
          Call stats_init(pparams%siting%siteavgs(i,j),' ',blocksize, &
              .False.,'f15.3')
        End Do
      End Do
    End If

    !** Initialize intramolecular structure histogram statistics
    If (Associated(pparams%structure)) Then
      Do i=1,pparams%structure%nparams
        Call histogram_reinit(pparams%structure%hist(i))
      End Do
    End If

    !** Initialize radial distribution functions
    If (Associated(pparams%radDist)) Then
      Call histogram_reinit(pparams%radDist%hist)
      pparams%radDist%n_loops=0
      pparams%radDist%ncalls=0


    End If

    !** Initialize intramolecular structure histogram statistics
    If (Associated(pparams%vels)) Then
      Do i=1,pparams%vels%natoms

        Call histogram_reinit(pparams%vels%atomcomp_hist(i,1))
        Call histogram_reinit(pparams%vels%atomcomp_hist(i,2))
        Call histogram_reinit(pparams%vels%atomcomp_hist(i,3))
        Call histogram_reinit(pparams%vels%atom_hist(i))
      End Do
      Call histogram_reinit(pparams%vels%com_hist)
    End If

  End Subroutine post_reInitStats



  !----------------------------------------------------------------------
  ! makes the datafile corresponding to simno ready to read from
  ! Reads the headre etc and brings the file pointer to the first data
  !----------------------------------------------------------------------
  Subroutine post_initconfilereading(pparams,simno)
    Type(PostInfo), Intent(InOut) :: pparams
    Integer, Intent(in)           :: simno

    Integer               :: tempunitnum
    Character(len=strLen) :: configfilename
    Logical :: is_old
    !** First make sure that datafile is closed, then open it, 
    !** Read its header etc
    configfilename=pparams%configFiles(simno)%name
    tempunitnum=file_close(configfilename)

    is_old = pparams%configFiles(simno)%datafile_is_old
    Call datafile_initin(pparams%configFiles(simno),configfilename, is_old)

  End Subroutine post_initconfilereading

  !----------------------------------------------------------------------
  ! Reads the configfile, update the averages
  ! configfile was generated during the run. The run might still be in 
  ! progress
  !----------------------------------------------------------------------
  Subroutine post_readandavg(pparams,simno)
    Type(PostInfo), Intent(InOut) :: pparams
    Integer, Intent(in)           :: simno

    Integer :: ndats, nmoles_now, i, l, spc, nsorbs
    Real(kind=RDbl) :: totnrg
    Logical :: hastime
    Logical :: fast, recalc, nointra, accels, success
    Real(kind=RDbl) :: max_nrg, dummydt, dummyT



#ifdef SITE_HACK
    Integer :: j,k, sum
    Logical , Save :: first_time=.True.
    Real(kind=RDbl) :: nrg1, nrg2, in1, in2, pressure, mix_ratio
    If (first_time) Then
      sunit1=file_open('hack-low.dat')
      sunit2=file_open('hack-high.dat')

      Write(sunit1,'(a)')"# ------ some site hack stuff for low alkane -----"
      Write(sunit1,'(a)')&
          "# press     | mix  | SC-SC | ZZ-ZZ | SC-ZZ | SC-sili&
          &| ZZ-sili|SC-intra|ZZ-intra"
      Write(sunit1,'(a)')&
          "# kpa       |ratio | %     | %     | %     | KJ/mol &
          &| KJ/mol |KJ/mol  | KJ/mol"


      Write(sunit2,'(a)')"# ------ some site hack stuff for high alkane -----"
      Write(sunit2,'(a)')&
          "# press     | mix  | SC-SC | ZZ-ZZ | SC-ZZ | SC-sili&
          &| ZZ-sili|SC-intra|ZZ-intra"
      Write(sunit2,'(a)')&
          "# kpa       |ratio | %     | %     | %     | KJ/mol &
          &| KJ/mol |KJ/mol  | KJ/mol"

      first_time=.False.
    Endif

    Do i=1,3
      Do j=1,3
        Do k=1,3
          Call stats_init(S_NC_STATS(i,j,k),' ',10,.False., 'f15.3')
        end do
        Call stats_init(S_INTRA_STATS(i,j),' ',10,.False., 'f15.3')
      end do
    end do

    n_single_site=0
    n_mixed=0
#endif

    ndats = post_noofconfigs()
    nsorbs = molecules_getnsorbs()
    hastime=datafile_hastime(pparams%configFiles(simno))
    !** Inform the user
    Write(*,'(a,i6,a)') " Expecting total of ",ndats, &
        " configurations in datafile :"//&
        Trim(pparams%configFiles(simno)%name)
    Write(*,'(a,i7,a,i7,a)') " This post code run will read from &
        & configuration - ", pparams%firstconf, " to configuration - ",&
        pparams%lastconf, " including both"


    !** read thru the data file, calculate necessary ensemble averages
    !** calculate various maps to be generated
    Do i=1, ndats
#ifdef DEBUG
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "reading data from control file, datapoint : ", i
#endif

      !** Check and increment the size of sorbates, If required 
      !** needed for simulations where number of molecs change
      Do l=1, nsorbs
        If (.Not.config_isfixed(pparams%species(l))) Then
          !use previous value of nmoles
          nmoles_now=Int(pparams%nmole_list(l))
          Call config_checkandincr(pparams%species, l, nmoles_now)
        End If
      End Do

      ! Reading one configuration from datafile
      ! Only the values stored will be returned, rest will be returned as zero
      ! content_tag stored in general.F90 can be used to check , what values 
      ! are scored
      ! spcstats gets nrg in KJ/mol
      If (hastime) Then
        Call datafile_readconfig(pparams%configFiles(simno), &
            pparams%species, &
            pparams%nmole_list,pparams%spcstats, totnrg, pparams%time)
      Else
        Call datafile_readconfig(pparams%configFiles(simno), &
            pparams%species, &
            pparams%nmole_list,pparams%spcstats, totnrg)
      End If

      ! Do the analysis only if we are in the asked range
      If ( (i>=pparams%firstconf) .And. (i<=pparams%lastconf) ) Then

        !** 1) set number of moles in species 2) do pbcs
        Do spc=1,nsorbs
          Call config_setnmoles(pparams%species(spc), &
              pparams%nmole_list(spc) )

          Call config_applypbcs(pparams%species,pparams%scell,(/spc,0,0/))
        End Do


        If (Associated(pparams%nrgs)) Then
          ! spcstats has nrg in KJ/mol
          Call post_updateNrgAvgs(pparams%nrgs, pparams%spcstats, &
              totnrg)
          ! imodel checks start here
          If (Associated(pparams%iact)) Then
            Call post_calcNewNrgs(pparams%iact,pparams)
          Endif


        End if

        !** Loading information
        If (Associated(pparams%loading)) Then
          Do spc=1,molecules_getnsorbs()
            Call stats_update(pparams%loading%totalavgs(spc), &
                one*pparams%nmole_list(spc))
          End Do
        End If

        !** Siting information
        If (Associated(pparams%siting)) Then
          Call post_updateSitingAvgs(pparams, pparams%siting, &
              pparams%species)
        End If

        !** velocity histograms information
        If (Associated(pparams%vels)) Then
          Call post_updateVelHists(pparams, pparams%vels, &
              pparams%species, pparams%nmole_list)
        End If

        !** Pore map.
        If (Associated(pparams%poremap)) Then
          Call post_updatePoreMapVals(pparams, pparams%poremap,  &
              pparams%species)
        End If

        !** Structure averages : Related to averaging and making histograms of
        ! bond lengths, angles, torsions etc. Useful for checking
        ! whether the configurational phase space is sampled adequately
        If (Associated(pparams%structure)) Then
          Call post_updateStructureAvgs(pparams%structure, pparams%species)
        End If

        !** Radial Distribution functions
        If (Associated(pparams%radDist)) Then
          Call post_updateRadDistAvgs(pparams%radDist, pparams%scell, &
              pparams%species)
        End If

        !** write some of the running averages stored in the stats variables
        Call post_writeoutputs(pparams,i)

      Else
        !** avoid overshooting the no of data points to be read
        !** Simulations might be in progress
        If (i>pparams%lastconf) Exit
      End If

    End Do

#ifdef SITE_HACK
    ! normalise
    Do i=1,3
      sum=0
      Do j=1,3
        Do k=1,3
          sum=sum+Int(HTSTATS(i,j,k))
        end do
      end do
      If (sum/=0) HTSTATS(i,1:3,1:3)=HTSTATS(i,1:3,1:3)*(100.00/sum)
    end do
    pressure=datafile_getextrainfo(pparams%configfiles(simno), &
        'pressure')
    mix_ratio=n_mixed*one/(n_mixed+n_single_site)
    nrg1=stats_getcavg(S_NC_STATS(1,3,2)) ! sc
    nrg2=stats_getcavg(S_NC_STATS(1,3,3)) ! zz
    in1=stats_getcavg(S_INTRA_STATS(1,2))
    in2=stats_getcavg(S_INTRA_STATS(1,3))
    Write(sunit1,'(1e12.3,f7.2,3f7.1,2f10.2,2f8.2)') pressure, mix_ratio,  &
        HTSTATS(1, 2, 2), HTSTATS(1, 3, 3), &
        HTSTATS(1, 2, 3)+ HTSTATS(1, 3, 2), nrg1, nrg2, in1, in2

    nrg1=stats_getcavg(S_NC_STATS(2,3,2)) ! sc
    nrg2=stats_getcavg(S_NC_STATS(2,3,3)) ! zz
    in1=stats_getcavg(S_INTRA_STATS(2,2))
    in2=stats_getcavg(S_INTRA_STATS(2,3))
    Write(sunit2,'(1e12.3,f7.2,3f7.1,2f10.2,2f8.2)') pressure, mix_ratio,  &
        HTSTATS(2, 2, 2), HTSTATS(2, 3, 3), &
        HTSTATS(2, 2, 3)+ HTSTATS(2, 3, 2), nrg1, nrg2, in1, in2


    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) " the matrix for sorb 1"
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,'(3f10.3)') HTSTATS(1, 1, 1),  HTSTATS(1, 1, 2),  HTSTATS(1, 1, 3)
    Write(*,'(3f10.3)') HTSTATS(1, 2, 1),  HTSTATS(1, 2, 2),  HTSTATS(1, 2, 3)
    Write(*,'(3f10.3)') HTSTATS(1, 3, 1),  HTSTATS(1, 3, 2),  HTSTATS(1, 3, 3)
    Write(*,*)
    Write(*,*) " the matrix for sorb 2"
    Write(*,'(3f10.3)') HTSTATS(2, 1, 1),  HTSTATS(2, 1, 2),  HTSTATS(2, 1, 3)
    Write(*,'(3f10.3)') HTSTATS(2, 2, 1),  HTSTATS(2, 2, 2),  HTSTATS(2, 2, 3)
    Write(*,'(3f10.3)') HTSTATS(2, 3, 1),  HTSTATS(2, 3, 2),  HTSTATS(2, 3, 3)
    Write(*,*)

    ! zero again
    HTSTATS=zero

#endif


  End Subroutine post_readandavg

  !---------------------------------------------------------------------
  ! calculate fresh nrgs with imodel
  !---------------------------------------------------------------------
  Subroutine post_calcNewNrgs(iact, pparams)
    Type(interactParams), Pointer        :: iact 
    Type(PostInfo), Intent(inout)  :: pparams

    Integer, Dimension(MAX_SORBS), Save :: new_nmoles, prev_nmoles, max_nmoles
    Logical , Save :: first_time=.True.

    Integer ::  i,j,m,nsorbs,calc_type,nmoles,nderivs,spc2
    Integer,Parameter:: FULL_SYS=0,MOL_SYS=1

    Logical :: fast, recalc, nointra, accels, success, incr_size,decr_size    
    Logical :: want_intra, want_ncoul, want_coul
    Real(kind=RDbl) :: max_nrg, dummydt, dummyT, intra, coul, ncoul, tot
    Real(kind=RDbl) :: i_tot,c_tot,nc_tot
    Real(kind=RDbl),Dimension(3) :: aux_params
    Character(len=xlstrLen)         :: string
    Type(EnergyPlus),Save :: ffout

    If (first_time) Then
      max_nmoles=0
      first_time=.False.
      nderivs = iact%imodel%results(1)%nderivs
      Call storebase_init(ffout,nderivs,.True.)  ! .True. for intra molecular
    Endif

    nsorbs=molecules_getnsorbs()
    calc_type=MOL_SYS


    fast=.True.
    recalc=.True.
    nointra=.False.
    accels=.False.
    max_nrg=100000.00_RDbl
    dummydt= 0.1
    dummyT=300.00
    aux_params=(/ max_nrg, dummydt, dummyT/)



    ! do pbcs REQUIRED OR NOT?? added by simon's changes in feb/mar-2004
    Do i=1,nsorbs
      Call config_applypbcs(pparams%species,pparams%scell,(/i,0,0/))
    Enddo

    Select Case(calc_type)
    Case(FULL_SYS)



      !
      ! Increase memmory if required
      ! ** Elaborate attempt to reduce the number of interact_resize calls
      ! ** Does not help. This routine leaks memory like crazy
      prev_nmoles=iact%imodel%results(1)%nmoles_system(1:nsorbs)
      Do i=1,nsorbs
        If (prev_nmoles(i)>max_nmoles(i)) max_nmoles(i)=prev_nmoles(i)
      End Do

      ! imodel can deal with old_nmoles, let's whether there is space for
      ! new_nmoles
      incr_size=.False.
      decr_size=.False.
      Do i=1,nsorbs
        new_nmoles(i) = pparams%nmole_list(i)
        If ( new_nmoles(i)> max_nmoles(i)) incr_size=.True.
        prev_nmoles(i) = pparams%nmole_list(i)
      End Do
      If (incr_size) Then
        !** resize imodel, this creates all kinds of memory problems
        Call interact_resize(iact%imodel, pparams%scell, pparams%species)
        max_nmoles=new_nmoles
      Else
        Call interact_changeAllNmoles(iact%imodel, new_nmoles(1:nsorbs))
      Endif

      !------------ end of increase memory --------------

      ! (/0,0,0/) means whole system
      success = interact_int(pparams%iact%imodel,&
          pparams%iact%imodel%results(1),pparams%species,&
          pparams%scell,fast,recalc,nointra,accels, &
          aux_params,(/0,0,0/))

      If (.Not.success) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        stop
      Endif

      ! compare total nrg vs stored energy

      !make sure the sums in the storage are internally consistent
      Call storetop_fastsum( pparams%iact%imodel%results(1))  
      ! assumes accurate mol-mol sums

#ifdef DEBUG
      Do i=1,nsorbs
        If (config_isfixed(pparams%species(i))) Cycle
        Do j=i,nsorbs
          Write(*,*) "SPC1, SPC2", i, j,  &
              pparams%iact%imodel%results(1)%ncoul%ab(i,j)%total%nrg*calToJ, &
              pparams%iact%imodel%results(1)%coul%ab(i,j)%total%nrg*calToJ
          ncoul=storestats_getnoncoul(pparams%spcstats, i, j, 'INST' )
          coul=storestats_getcoul(pparams%spcstats, i, j, 'INST' )
          Write(*,*) "SPC1, SPC2", i, j, ncoul,coul
        End Do
      End Do


      Do i=1,nsorbs
        Write(*,*) "SPC1", i,  &
            Sum(pparams%iact%imodel%results(1)%intra(i)%total%intranrg)*calToJ
        intra=storestats_getintranrg(pparams%spcstats, i, 'INST' )
        Write(*,*) "SPC1,", i,intra
      End Do
#endif




      ! our favorite case
    Case(MOL_SYS)

      iact%coul=zero
      iact%ncoul=zero
      iact%intra=zero


      Do i=1,nsorbs
        !        tot=zero
        If (config_isfixed(pparams%species(i))) Cycle        
        nmoles=config_getnmoles(pparams%species,i)
        Do m=1,nmoles
          ! make call to get nrg of this mol vs sys
          ! molsys_nrg is actually an spc_sys storage. 
          ! It is basically storage for an spc With one molecule in it
          success= interact_int(iact%imodel, iact%molsys_nrg(i), &
              pparams%species,&
              pparams%scell,fast,recalc,nointra,accels, &
              aux_params,(/i,m,0/),(/0,0,0/))

          If (.Not. success) Then
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            stop
          Endif

          !          tot=tot+iact%molsys_nrg(i)%total%nrg 
          ! +sum(iact%molsys_nrg(i)%total%intranrg)

          ! intra nrgies first
          want_intra=.True.
          want_ncoul=.False.
          want_coul=.False.
          success= storetop_extract(iact%molsys_nrg(i),  want_intra, &
              want_ncoul, want_coul,ffout, (/i,0,0/) )

          iact%intra(i,m)= sum(ffout%intranrg)*calToJ

          ! NCOUL : molsys_nrg is actually a spc-sys nrg 
          ! so use  (/i,0,0/) instead of  (/i,m,0/)
          Do spc2=1,nsorbs
            want_intra=.False.
            want_ncoul=.True.
            want_coul=.False.
            success= storetop_extract(iact%molsys_nrg(i), want_intra, &
                want_ncoul, want_coul,ffout, (/i,0,0/),(/spc2,0,0/))
            iact%ncoul(i,m,spc2)=  ffout%nrg*calToJ
          End Do

          ! COUL : molsys_nrg is actually a spc-sys nrg 
          ! so use  (/i,0,0/) instead of  (/i,m,0/)
          Do spc2=1,nsorbs
            want_intra=.False.
            want_coul=.True.
            want_ncoul=.False.
            success= storetop_extract(iact%molsys_nrg(i), want_intra, &
                want_ncoul, want_coul, ffout, (/i,0,0/),(/spc2,0,0/))
            iact%coul(i,m,spc2)=  ffout%nrg*calToJ
          end do

        End Do
      End Do
#ifdef DEBUG
      ! add up nrgs
      Do i=1,nsorbs-1
        Do spc2=i,nsorbs
          c_tot=zero
          nc_tot=zero
          Do m=1,config_getnmoles(pparams%species,i)
            c_tot=c_tot+iact%coul(i,m,spc2)
            nc_tot=nc_tot+iact%ncoul(i,m,spc2)
          end do
          If(i==spc2)  c_tot= c_tot/2
          If(i==spc2)  nc_tot= nc_tot/2
          !          Write(*,*) " ** SPC1, SPC2", i, spc2, nc_tot*calToj, c_tot*calToj
        End Do
      end do

      Do i=1,nsorbs
        i_tot=zero
        Do m=1,config_getnmoles(pparams%species,i)
          i_tot=i_tot+iact%intra(i,m)
        end do
        !        Write(*,*) " ** SPC1 ", i, i_tot*calToj
      End Do
#endif

    Case  Default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select



  End Subroutine post_calcNewNrgs

  !---------------------------------------------------------------------
  ! update the values in nrg average structure
  ! Note : 1) This basically copies one stats structure into another.
  !           spcstats has most of what we want. But we would like post.F90
  !           to have its own update routine
  !        2) change to call functions from storestats instead of 
  !           directly accessing the fields
  !        3) spcstats has nrg in KJ/mol
  !---------------------------------------------------------------------
  Subroutine post_updateNrgAvgs(nrgstat, spcstats, totalinst)
    Type(nrgAvg), Intent(inout)                   :: nrgstat
    Type(Species_Stats)   :: spcstats 
    Real(kind=RDbl), Intent(in)                  :: totalinst

    Integer :: i, j
    Real(kind=RDbl)                 :: value

    Call stats_update(nrgstat%total, totalinst )
    Do i=1,molecules_getnsorbs()
      Do j=i,molecules_getnsorbs()

        value=storestats_getnoncoul(spcstats, i, j, 'INST' )
        Call stats_update(nrgstat%noncoulNrg(i,j), value)

        value=storestats_getcoul(spcstats, i, j, 'INST' )
        Call stats_update(nrgstat%coulNrg(i,j), value)
      End Do

      value=storestats_getke(spcstats, i, 'INST' )
      Call stats_update(nrgstat%ke(i), value)

      value=storestats_getintranrg(spcstats, i, 'INST', 'TOTAL' )
      Call stats_update(nrgstat%intra(i, TOTAL_INDEX), value)

    End Do

  End Subroutine post_updateNrgAvgs


  !---------------------------------------------------------------------
  ! update the values in related to structure of molecules 
  !---------------------------------------------------------------------
  Subroutine post_updateStructureAvgs(structure, species)
    Type(structureAvg)                            :: structure
    Type(AtMolCoords), Dimension(:), Pointer      :: species

    Integer :: analtype, nmoles, m, sorb, a1, a2, a3, a4
    Real(kind=RDbl)                 :: value

    !** Loop through each type of analysis, each type refers to a 
    !** prticular bond or angle of a particular sorbate Type
    Do analtype=1,structure%nparams
      sorb=molecules_gettype(cleanstring(structure%sorbnames(analtype)))
      If ( (sorb<1) .or. (sorb> size(species)) )  Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      nmoles=config_getnmoles(species, sorb)
      Select Case(cleanstring(structure%intra_type(analtype)))
      Case("BONDLENGTHS") 
        a1 = structure%atom_list(analtype, 1)
        a2 = structure%atom_list(analtype, 2)
        Do m=1,nmoles
          value = config_getdistance(species, sorb, m, &
              a1, sorb, m, a2)
          Call histogram_update(structure%hist(analtype),value)
        End Do
      Case("ANGLES") 
        a1 = structure%atom_list(analtype, 1)
        a2 = structure%atom_list(analtype, 2)
        a3 = structure%atom_list(analtype, 3)
        Do m=1,nmoles
          value = config_getangle(species, sorb, m, a1, a2, a3)
          Call histogram_update(structure%hist(analtype), value*radToDeg)
        End Do

      Case("TORSIONANGLE") 
        a1 = structure%atom_list(analtype, 1)
        a2 = structure%atom_list(analtype, 2)
        a3 = structure%atom_list(analtype, 3)
        a4 = structure%atom_list(analtype, 4)
        Do m=1,nmoles
          value = config_gettorangle(species, sorb, m, a1, a2, a3, a4)
          Call histogram_update(structure%hist(analtype),value*radToDeg)
        End Do
      Case Default 
        Write(*,*) "can recognise option : ", &
            Trim(Trim(structure%intra_type(analtype)))
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End Select
    End Do

  End Subroutine post_updateStructureAvgs


  !---------------------------------------------------------------------
  ! update the values in related to structure of molecules 
  !---------------------------------------------------------------------
  Subroutine post_updateRadDistAvgs(param, scell, species )
    Type(radDistAvg)                            :: param
    Type(AtMolCoords), Dimension(:), Pointer      :: species
    Type(Simcell_Params), Intent(In) :: scell

    Integer :: spc1, spc2, m1, m2, mols1, mols2, a1, a2
    Real(kind=RDbl) :: dist, dist2
    Type(VecType) :: vec1, vec2, diff
    Logical :: silentFlag


    spc1=param%sorbs(1)
    spc2=param%sorbs(2)
    mols1=config_getnmoles(species,spc1)
    mols2=config_getnmoles(species,spc2)
    a1=param%atoms(1)
    a2=param%atoms(2)

    silentFlag=.True.

    ! Loop over first molecule and atom type
    Do m1=1,mols1
      param%n_loops=param%n_loops+1
      ! decide whether COM based or atom-based
      If (param%com_based(1)) Then
        vec1=config_getMolecCOM(species,spc1,m1)
      Else
        vec1=species(spc1)%coords(a1,m1)%r
      Endif

      ! Loop over each molecule of other species
      Do m2=1,mols2


        ! decide whether COM based or atom-based, for this species
        If (param%com_based(2)) Then
          vec2=config_getMolecCOM(species,spc2,m2)
        Else
          vec2=species(spc2)%coords(a2,m2)%r
        Endif

        ! calulate distance
        Call simcell_minimage(scell, vec1, vec2, diff, dist2)

        dist=Sqrt(dist2)

        If ((spc1==spc2).And.(m1==m2)) Cycle

        ! note this histogram will have to be reweighted later
        Call histogram_update(param%hist,dist,silentFlag)
        param%ncalls=param%ncalls+1
      End Do
    End Do
  End Subroutine post_updateRadDistAvgs


  !---------------------------------------------------------------------
  ! update the values in siting structure 
  !---------------------------------------------------------------------
  Subroutine post_updatePoreMapVals(pparams, poremap, species)
    Type(PostInfo), Intent(InOut)      :: pparams
    Type(PoreMapParams), Intent(inout)                   :: poremap
    Type(AtMolCoords), Dimension(:), Intent(inout)  :: species 

    Integer :: mp, sorb, natoms, a, al, m, x_index, y_index, z_index, nx
    Type(VecType) :: com, uc_com, vec, uc_vec

    ! loop over map sorbs
    Do mp=1,poremap%nmaps
      sorb=molecules_gettype(Trim(poremap%sorbnames(mp)))
      natoms=molecules_getnatoms(sorb)
      ! loop over moleucles
      Do m=1,config_getnmoles(species(sorb))

        If (poremap%onlycom) Then

          com=config_getMolecCOM(species, sorb, m)
          uc_com=simcell_maptouc(pparams%scell,com)
          x_index = Int(uc_com%comp(1)/poremap%dx)+1
          y_index = Int(uc_com%comp(2)/poremap%dy)+1
          z_index = Int(uc_com%comp(3)/poremap%dz)+1
          poremap%fillArray(mp,x_index,y_index,z_index)= &
              poremap%fillArray(mp,x_index,y_index,z_index)+1
        else

          Do a=1,natoms

            ! actual atom number
            al=a

            ! works only for one map
            If (poremap%atom_list_given) Then
              If (poremap%nmaps>1) Then
                Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
                Stop
              Endif
              If (a>poremap%natoms) Exit ! out of the natoms do loop
              al=poremap%atom_list(a)
              If (al>natoms) Then
                Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
                stop
              Endif
            Endif

            vec=species(sorb)%coords(al,m)%rp
            uc_vec=simcell_maptouc(pparams%scell,vec)
            x_index = Int(uc_vec%comp(1)/poremap%dx)+1
            y_index = Int(uc_vec%comp(2)/poremap%dy)+1
            z_index = Int(uc_vec%comp(3)/poremap%dz)+1
            poremap%fillArray(mp,x_index,y_index,z_index)= &
                poremap%fillArray(mp,x_index,y_index,z_index)+1
          end do
        endif
      End Do

    End Do

  End Subroutine post_updatePoreMapVals


  !---------------------------------------------------------------------
  ! Write poremap to files/ dmap too 
  !---------------------------------------------------------------------
  Subroutine post_writePoreMap(poremap)
    Type(PoreMapParams), Intent(in)                   :: poremap

    Integer, Dimension(:), Pointer :: contlevels
    Integer :: error, max, mp, c, j, k, l, m, ncells, i_ucx, i_ucy, i_ucz
    Integer :: x_ind,y_ind, z_ind, unitno, incr
    Character(len=4) :: atom_name, num_strg
    Character(len=lstrLen) :: newfilename
    Real(kind=RDbl) :: fac, scale , sum

    Real(kind=RSgl), Dimension(:,:,:), Pointer :: dmap_array
    Character(len=25) :: dmap_molec_name
    Character(len=strLen) :: dmapname
    Integer           :: dmap_nx, dmap_ny, dmap_nz
    Real(kind=RDbl)   :: dmap_sum

    !** Allocate the contour level array
    Allocate(contlevels(poremap%ncont),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Do mp=1,poremap%nmaps

      max=0
      !*** Find maximum
      Do k=1,poremap%nbins(1)
        Do l=1,poremap%nbins(2)
          Do m=1,poremap%nbins(3)
            If (poremap%fillArray(mp,k,l,m)>max) Then
              max=poremap%fillArray(mp,k,l,m)
            End If
          End Do
        End Do
      End Do

      ! assign contour levels
      If (poremap%cont_div_exponential) Then
        fac=Exp(Log(max*one)/poremap%ncont)
        contlevels(1)=1
        If (poremap%ncont>1) Then
          Do j=2,poremap%ncont
            contlevels(j)= (contlevels(j-1)*fac)+1
          End Do
        End If
      ELse
        contlevels(1)=1
        If (poremap%ncont>1) Then
          incr=(max-one)/poremap%ncont-1
          Do j=2,poremap%ncont
            contlevels(j)= contlevels(j-1) + incr
          End Do
          contlevels(poremap%ncont)= max 
        End If
      Endif

      Write(*,*) "contour levels are " , contlevels, "max was", max


      Do c=1,poremap%ncont

        !** Open the disk file to write to
        If (c<1000) Then
          num_strg=int2str(c)
          num_strg=Adjustl(num_strg)
        Else
          Write(*,*) "Too many contors specified, exiting"
          Stop
        End If

        newfilename=Trim(poremap%mapname(mp))//"."//Trim(num_strg)
        unitno=file_open(newfilename)
        atom_name= "A"//Trim(num_strg)

        ! default atom_name is given above, 
        ! the below names are more convenient
        Select Case(c)
        Case(1)
          atom_name="He"
        Case(2)
          atom_name="Ne"
        Case(3)
          atom_name="Ar"
        Case(4)
          atom_name="Kr" 
        Case(5)
          atom_name="Xe" 
        Case(6)
          atom_name="Rn" 
        Case default
          ! keep the default name only
        End Select

        ncells=0
        If (c==1) dmap_sum=zero
        !** check how many cells belong to this pore
        Do k=1,poremap%nbins(1)
          Do l=1,poremap%nbins(2)
            Do m=1,poremap%nbins(3)
              If (c==1) dmap_sum=dmap_sum+one*poremap%fillArray(mp,k,l,m)
              If (poremap%fillArray(mp,k,l,m) >=contlevels(c)) Then
                !** this is one of the members of the pore map
                ncells=ncells+1
              End If
            End Do
          End Do
        End Do

        Write(*,*) "Total fillings is ",ncells

        !Write number of points in rasmol file, ncells* no of unitcells
        Write(unitno,*) ncells * ( poremap%uc_repeat_info(1) * &
            poremap%uc_repeat_info(2) * poremap%uc_repeat_info(3) )
        Write(unitno,*) "some description"

        !** if "X" is the position in angstroms, then the *.xyz will have 
        !** "scale * X" as the position values
        scale = poremap%scaling_factor

        !** go through it again, this time write the pore map too
        Do k=1,poremap%nbins(1)
          Do l=1,poremap%nbins(2)
            Do m=1,poremap%nbins(3)
              If (poremap%fillArray(mp,k,l,m) >= contlevels(c)) Then
                !** this is one of the members of the pore map
                !** They have to be written to all unitcells
                Do i_ucx=1,poremap%uc_repeat_info(1)
                  Do i_ucy=1,poremap%uc_repeat_info(2)
                    Do i_ucz=1,poremap%uc_repeat_info(3)
                      x_ind=(i_ucx-1)*poremap%nbins(1) + k
                      y_ind=(i_ucy-1)*poremap%nbins(2) + l
                      z_ind=(i_ucz-1)*poremap%nbins(3) + m

                      Write(unitno,'(a,3f16.6)') Trim(atom_name), &
                          scale*x_ind*poremap%dx, &
                          scale*y_ind*poremap%dy, &
                          scale*z_ind*poremap%dz 

                    End Do
                  End Do
                End Do
                ! end of unit cells loop

              End If


            End Do
          End Do
        End Do
        !** end of rasmol%fillArray(sorb,:,:,:) loop

        Close(unitno)
      End Do
      !** end of conours loop

      !********* WRITING DMAP ****************************
      If (poremap%dump_dmap) Then
        dmap_molec_name=Trim(poremap%sorbnames(mp))
        dmap_nx=poremap%nbins(1)
        dmap_ny=poremap%nbins(2)
        dmap_nz=poremap%nbins(3)
        dmapname=Trim(poremap%mapname(mp))//".dmap"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "writing dmap (oepning as a new file) : ", Trim(dmapname)

        Allocate(dmap_array(dmap_nx, dmap_ny, dmap_nz),STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        dmap_array(1:dmap_nx, 1:dmap_ny,1:dmap_nz) = &
            poremap%fillArray(mp, 1:dmap_nx, 1:dmap_ny,1:dmap_nz) /dmap_sum
        unitno=file_open(dmapname,101)
        ! 2= probability density
        Write(unitno) dmap_molec_name, 2
        Write(unitno) dmap_nx, dmap_ny, dmap_nz
        Write(unitno) dmap_array
        close(unitno)
        Deallocate(dmap_array,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Endif

    End Do

  End Subroutine post_writePoreMap


  !---------------------------------------------------------------------
  ! Write to files 
  !---------------------------------------------------------------------
  Subroutine post_writestructureAvgs(pparams, structure, simno, unitno )
    Type(PostInfo), Intent(InOut)      :: pparams
    Type(structureAvg), Intent(inout)                   :: structure
    Integer,Intent(in) :: simno, unitno

    Integer :: analtype ,a1,a2,a3,a4

    Do analtype=1,structure%nparams 
      Write(unitno,'(a)') "## Histogram of "//&
          Trim(structure%intra_type(analtype))//" for the ensemble"
      Write(unitno,'(a)') "## configfile : "//&
          Trim(pparams%configFiles(simno)%name)
      Write(unitno, '(a)') "## Molecule : "//&
          cleanstring(structure%sorbnames(analtype))

      Select Case(cleanstring(structure%intra_type(analtype)))
      Case ("BONDLENGTHS")
        a1=structure%atom_list(analtype,1)
        a2=structure%atom_list(analtype,2)
        Write(unitno,'(a,2i6)') "## atoms ", a1,a2
      Case ("ANGLES")
        a1=structure%atom_list(analtype,1)
        a2=structure%atom_list(analtype,2)
        a3=structure%atom_list(analtype,3)
        Write(unitno,'(a,3i6)') "## atoms ", a1,a2,a3
      Case ("TORSIONANGLE")
        a1=structure%atom_list(analtype,1)
        a2=structure%atom_list(analtype,2)
        a3=structure%atom_list(analtype,3)
        a4=structure%atom_list(analtype,4)
        Write(unitno,'(a,4i6)') "## atoms ", a1,a2,a3,a4
      Case Default
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End Select

      Call histogram_display(structure%hist(analtype), unitno)
      Write(unitno,'(a)') " "
    End Do

  End Subroutine post_writestructureAvgs

  !---------------------------------------------------------------------
  ! Write to files 
  !---------------------------------------------------------------------
  Subroutine post_writeRadDistAvgs(pparams, raddist, simno, unitno, &
      num_density)
    Type(PostInfo), Intent(InOut)      :: pparams
    Type(radDistAvg), Intent(inout)                   :: raddist
    Integer,Intent(in) :: simno, unitno
    Real(kind=Rdbl), Dimension(4), Optional, Intent(in) :: num_density

    Integer :: a1, a2, error, nbins, i
    Real(kind=Rdbl)       :: density, r, low_r, high_r, dr, dV
    Real(kind=Rdbl),Dimension(:), Pointer :: histg

    nbins=raddist%nbins

    Write(*,*) "## Total histogram calls : ", raddist%ncalls
    Write(*,*) "## Total Number of mol1 loops : ", raddist%n_loops
    Write(*,*) "## Histogram Elements    : ", raddist%hist%nelements
    Write(*,*) "## % of data outside range ", &
        100*( one - (Real(raddist%hist%nelements)/raddist%ncalls))


    Write(unitno,'(a)') "## Unweghted Histogram of  radial &
        & distribution for the ensemble"
    Write(unitno,'(a)') "## configfile : "//&
        Trim(pparams%configFiles(simno)%name)
    Write(unitno, '(a)') "## Molecule-1 : "//&
        cleanstring(raddist%sorbnames(1))
    Write(unitno, '(a)') "## Molecule-2 : "//&
        cleanstring(raddist%sorbnames(2))
    a1=raddist%atoms(1)
    a2=raddist%atoms(2)
    Write(unitno,'(a,2i6)') "## atoms (If zero  implies COM) :", a1,a2
    Call histogram_display(raddist%hist, unitno)
    Write(unitno,'(a)') " "

    ! copy the values to temp_array
    Allocate(histg(nbins),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    histg=one*(raddist%hist%histg)

    If (Present(num_density)) Then
      density= num_density(raddist%sorbs(2))
      Write(unitno,'(a)') "## Properly Weighted Radial distribution for the ensemble" 
      Write(unitno,'(a,f16.4)') "## Density (spc 2) : ", density
    Else
      Write(unitno,'(a)') "## Weighted Radial distribution for the ensemble" 
      Write(unitno,'(a)') "## Multiply by number density to get proper g(r)"
      density=one
    Endif

    dr=(raddist%maxdist-raddist%mindist)/nbins
    Do i=1,nbins
      high_r=i*dr+raddist%mindist
      low_r=high_r-dr
      dV=(4.0*pi/3)*( ((high_r)**3) - ((low_r)**3) )
      r=(low_r+high_r)/2
      histg(i)=histg(i)/(density * dV * raddist%n_loops)
      Write(unitno,'(2f16.4)') r, histg(i)
    End Do
    Write(unitno,'(a)') " "

  End Subroutine post_writeRadDistAvgs




  !---------------------------------------------------------------------
  ! Write to files 
  !---------------------------------------------------------------------
  Subroutine post_writeVelHists(pparams, vels, simno, unitno )
    Type(PostInfo), Intent(InOut)      :: pparams
    Type(velocityParams), Intent(inout)                   :: vels
    Integer,Intent(in) :: simno, unitno

    Integer :: i,a 
    Write(unitno,'(a)') "## Velocity Histogram for Sorbate "//&
        Trim(vels%spcname)
    Do i=1,Size(vels%display_list)
      a=vels%display_list(i)
      Write(unitno,'(a,i3,a)') "## Histogram of x velocity for atom :",a,&
          " for the ensemble"
      Call histogram_display(vels%atomcomp_hist(a,1), unitno)

      Write(unitno,'(a,i3,a)') "## Histogram of y velocity for atom :",a,&
          " for the ensemble"
      Call histogram_display(vels%atomcomp_hist(a,2), unitno)

      Write(unitno,'(a,i3,a)') "## Histogram of y velocity for atom :",a,&
          " for the ensemble"
      Call histogram_display(vels%atomcomp_hist(a,3), unitno)


      Write(unitno,'(a,i3,a)') "## Histogram of TOTAL velocity for atom :",a,&
          " for the ensemble"
      Call histogram_display(vels%atom_hist(a), unitno)

    End Do

    Write(unitno,'(a,i3,a)') "## Histogram of com velocity for atom &
        & for the ensemble"
    Call histogram_display(vels%com_hist, unitno)

  End Subroutine post_writeVelHists


  !---------------------------------------------------------------------
  ! update the values in siting structure (based on com or one atom)
  !---------------------------------------------------------------------
  Subroutine post_updateSitingAvgs(pparams, siting, species)
    Type(PostInfo), Intent(In) :: pparams
    Type(sitingAvg), Intent(inout)                   :: siting 
    Type(AtMolCoords), Dimension(:), Intent(inout)  :: species 

    Integer :: nsorbs, sorb, s, a, m, temp_loading, sitetype
    Logical :: com_based
#ifdef SITE_HACK
    Integer :: headsite, tailsite, natoms, new_site, k
    Real(kind=RDbl) :: nc_nrg, i_nrg
    Logical :: single_site
#endif

    nsorbs=molecules_getnsorbs()

    !** Loop over all types of molecules
    Do sorb=1,nsorbs
      com_based=siting%com_based
      ! avoid zeolite type of species
      If (config_isfixed(species(sorb))) Cycle

      ! decide whether COM or atom-based
      If (.Not.(com_based)) Then
        a=siting%atom_num
        If (a>molecules_getnatoms(sorb)) Then
          ! we have less atoms for this sorb, so go bak to atom_based
          com_based=.True.
        Endif
      Endif
      ! loop over all site
      Do s=1,siting%nsites

        temp_loading =0

        ! assumes that config_setnmoles was called in datafile_readconfig
        Do m=1,config_getnmoles(species, sorb)
          If (com_based) Then
            sitetype=post_getsite(siting%sitemap,species,pparams%scell, &
                sorb,m)
          Else
            sitetype=post_getsite(siting%sitemap,species, pparams%scell, &
                sorb,m,a)
          Endif

#ifdef SITE_HACK
          ! we are going to combine site-1 and site-2
          If (sitetype==1) sitetype=2
#endif
          If (sitetype==siting%sites(s)) Then
            temp_loading =temp_loading +1
          End If



        End Do ! ** molecs loop

        Call stats_update(siting%siteavgs(s,sorb), temp_loading*one )
        siting%site_inst(s,sorb)=temp_loading
      End Do ! ** site loop
    End Do ! ** sorb loop

#ifdef SITE_HACK
    Do sorb=1,nsorbs
      If (config_isfixed(species(sorb))) Cycle
      natoms=molecules_getnatoms(sorb)
      Do m=1,config_getnmoles(species, sorb) 
        headsite=post_getsite(siting%sitemap,species,pparams%scell, &
            sorb,m,1)
        If (headsite==1) headsite=2             !MORE HACK

        tailsite=post_getsite(siting%sitemap,species,pparams%scell, &
            sorb,m,natoms)
        If (tailsite==1) tailsite=2             !MORE HACK

        !          If (headsite/=tailsite) Then 
        !          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        !          Stop
        !          Endif

        HTSTATS(sorb,headsite,tailsite)=HTSTATS(sorb,headsite,tailsite)+1

        ! find how many are mixed

        single_site=.True.
        Do a=2,natoms
          new_site=post_getsite(siting%sitemap,species,pparams%scell, &
              sorb,m,a)
          If (new_site==1) new_site=2
          If(new_site/=headsite) Then
            single_site=.False.
            n_mixed=n_mixed+1
            Exit
          Endif
          If (a==natoms) n_single_site= n_single_site+1
        End Do


        IF (associated(pparams%iact)) Then

          ! If single_site Then find energy and update
          If (single_site) Then
            Do k=1,nsorbs
              nc_nrg=pparams%iact%ncoul(sorb,m,k)
              Call stats_update(S_NC_STATS(sorb,k,headsite),nc_nrg) 
            end do
            i_nrg=pparams%iact%intra(sorb,m)
            Call stats_update(S_INTRA_STATS(sorb,headsite),I_nrg) 
          Endif
        Endif
      End Do !-- mol --!
    End Do !-- sorb --!
#endif



  End Subroutine post_updateSitingAvgs
  !------------------------------------------------------------------
  ! gives the site of an atom
  ! spc- type of molecule, m=molec number, a=atom number
  ! if a is not given uses COM
  !------------------------------------------------------------------
  Integer Function post_getsite(sitemap,species,scell,spc,m,a)
    Type(Smap_Params),Pointer                   :: sitemap
    Type(AtMolCoords), Dimension(:), Intent(inout)  :: species 
    Type(Simcell_Params),Intent(in)  :: scell 
    Integer,Intent(in) :: spc, m
    Integer,Intent(in),Optional :: a
    Type(VecType) :: pos, uc_vec
    If (present(a)) Then
      pos=species(spc)%coords(a,m)%rp
      uc_vec=simcell_maptouc(scell,pos)
    Else
      pos=config_getMolecCOM(species, spc, m)
      uc_vec=simcell_maptouc(scell,pos)
    Endif

    !** Update the siting averages based on com
#ifdef SMAPNOZERO
    post_getsite=smap_getPosSiteType(sitemap,uc_vec,.True.)
#endif
#ifndef SMAPNOZERO
    post_getsite=smap_getPosSiteType(sitemap,uc_vec)
#endif

  End Function post_getsite



  !---------------------------------------------------------------------
  ! update the velocity histograms 
  !---------------------------------------------------------------------
  Subroutine post_updateVelHists(pparams, vels, species, nmol_list)
    Type(PostInfo), Intent(In) :: pparams
    Type(velocityParams), Intent(inout)                   :: vels 
    Type(AtMolCoords), Dimension(:), Intent(inout)  :: species 
    Integer, Dimension(:), Intent(in)  :: nmol_list


    Integer :: nsorbs, nmols, natoms, a,m, temp_loading, sitetype,spc
    Type(VecType) :: com, uc_com
    Real(kind=RDbl) :: v
    spc=vels%spctype
    nmols=nmol_list(spc)
    natoms=molecules_getnatoms(spc)
    Do m=1,nmols
      Do a=1,natoms
        Call histogram_update(vels%atomcomp_hist(a,1),&
            species(spc)%coords(a,m)%v%comp(1))
        Call histogram_update(vels%atomcomp_hist(a,2),&
            species(spc)%coords(a,m)%v%comp(2))
        Call histogram_update(vels%atomcomp_hist(a,3),&
            species(spc)%coords(a,m)%v%comp(3))
        v=mag(species(spc)%coords(a,m)%v)
        ! just for symmetry in plotting we update +v, -v equally
        Call histogram_update(vels%atom_hist(a),v)
        Call histogram_update(vels%atom_hist(a),-v)


      End Do
      v=mag(config_getcomvel(species,spc,m))
      ! just for symmetry in plotting we update +v, -v equally
      Call histogram_update(vels%com_hist,v)
      Call histogram_update(vels%com_hist,-v)
    End Do

  End Subroutine post_updateVelHists


  !----------------------------------------------------------------------
  ! Write the stat variables at the end of each block
  !----------------------------------------------------------------------
  Subroutine post_writeoutputs(pparams,confignum)
    Type(PostInfo), Intent(InOut) :: pparams
    Integer , Intent(in) :: confignum
    Integer :: outputnum , itnum, blocksize, i, site, spc

    If (Associated(pparams%nrgs)) Then
      blocksize=pparams%nrgs%statblocksize
      If (Mod(confignum-pparams%firstconf+1,blocksize)==0)     Then

        outputnum=(confignum-pparams%firstconf+1)/blocksize

        ! original itn number
        itnum=confignum*genParams%iconfig

        ! ** output array is Real only
        pparams%nrgs%output_array(outputnum,1)=itnum*one
        pparams%nrgs%output_array(outputnum,2)=pparams%nrgs%total%inst
        pparams%nrgs%output_array(outputnum,3)=&
            pparams%nrgs%total%blockavg
        pparams%nrgs%output_array(outputnum,4)=pparams%nrgs%total%cumavg
      Elseif(confignum==pparams%lastconf) Then 
        ! some times the block size wont divide exactly into the total 
        ! number of configs so the last block will only be partially filled
        ! original itn number
        outputnum=Size(pparams%nrgs%output_array,1)
        itnum=confignum*genParams%iconfig

        ! ** output array is Real only
        pparams%nrgs%output_array(outputnum,1)=itnum*one
        pparams%nrgs%output_array(outputnum,2)=pparams%nrgs%total%inst
        pparams%nrgs%output_array(outputnum,3)=&
            pparams%nrgs%total%blockavg
        pparams%nrgs%output_array(outputnum,4)=pparams%nrgs%total%cumavg
      End If
    End If


    If (Associated(pparams%loading)) Then
      blocksize=pparams%loading%statblocksize
      If (Mod(confignum-pparams%firstconf+1,blocksize)==0)       Then
        outputnum=((confignum-pparams%firstconf+1)/blocksize) 
        itnum=confignum*genParams%iconfig

        Do spc=1,molecules_getnsorbs()
          ! original itn number
          pparams%loading%output_array(spc,outputnum,1)=(itnum)*one

          pparams%loading%output_array(spc,outputnum,2)= &
              pparams%loading%totalavgs(spc)%inst
          pparams%loading%output_array(spc,outputnum,3)= &
              pparams%loading%totalavgs(spc)%blockavg
          pparams%loading%output_array(spc,outputnum,4)= &
              pparams%loading%totalavgs(spc)%cumavg
        End Do
      End If
    End If



    If (Associated(pparams%siting)) Then
      blocksize=pparams%siting%statblocksize
      If (Mod(confignum-pparams%firstconf+1,blocksize)==0)       Then
        outputnum=((confignum-pparams%firstconf+1)/blocksize) 
        itnum=confignum*genParams%iconfig

        Do site=1,pparams%siting%nsites
          Do spc=1,molecules_getnsorbs()
            ! original itn number
            pparams%siting%output_array(site,spc,outputnum,1)=(itnum)*one
            pparams%siting%output_array(site,spc, outputnum,2)= &
                pparams%siting%siteavgs(site,spc)%inst
            pparams%siting%output_array(site,spc, outputnum,3)= &
                pparams%siting%siteavgs(site,spc)%blockavg
            pparams%siting%output_array(site,spc, outputnum,4)= &
                pparams%siting%siteavgs(site,spc)%cumavg
          Enddo
        End Do
      End If
    End If

  End Subroutine post_writeoutputs


  !-----------------------------------------------------------------------
  ! returns the number of configs used by post code, It includes the
  ! skipped parts of the run also.
  !-----------------------------------------------------------------------
  Integer Function post_noofconfigs()
    post_noofconfigs=general_getnoofconfigs()
  End Function post_noofconfigs

  !------------------------------------------------------------------------
  ! isoTdat - 1st column is pressure kPa,          (fug) 
  !           2nd column is loading molec/uc       (z)
  ! Fits to fug= K * z / (z_max -z)
  ! err= is error in z recalcualted from new fit in percetage of max loading
  !------------------------------------------------------------------------
  Subroutine post_langmuirFit(isoTdat, z_max, K, err)
    Real(kind=RDbl),Dimension(:,:), Pointer :: isoTdat
    Real(kind=RDbl),Intent(out) :: z_max, K, err
    Real(kind=RDbl),Dimension(:,:), Allocatable :: temp_array
    Integer :: nsims, sim, error
    Real(kind=RDbl) :: z, m, c, new_z, max
    nsims=Size(isoTdat,2)

    Allocate(temp_array(2,nsims),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do sim=1,nsims
      If (isoTdat(2,sim)<1.0e-8) Then
        Write(*,*) "error in Langmuir fit calculations"
        Write(*,*) "adjusting loading value to 1.0e-8 to avoid floating error"
        isoTdat(2,sim)=1.0e-8
        Write(*,*) "Everything is OK, except the langmuir fit"
      Endif
      If (isoTdat(1,sim)<1.0e-8) Then
        Write(*,*) "error in Langmuir fit calculations"
        Write(*,*) "adjusting some low pressure vlaues to 1.0e-7 to &
            &avoid floating error"
        isoTdat(1,sim)=1.0e-8
        Write(*,*) "Everything is OK, except the langmuir fit"
      Endif
      temp_array(1,sim)=one/isoTdat(1,sim)
      temp_array(2,sim)=one/isoTdat(2,sim)
    End Do

    ! do least square 
    Call mmath_LSFit(temp_array(1,1:nsims), temp_array(2,1:nsims), &
        nsims, m, c, err)

    z_max=1/c
    K=m/c
    err=zero
    max=1.0e-10
    Do sim=1,nsims
      z=isoTdat(2,sim)
      new_z=z_max*(isoTdat(1,sim))/(isoTdat(1,sim)+K )
      err=err+(new_z-z)**2
      If (z>max) max=z
    End Do
    ! MSQD in %
    err=100*Sqrt(err/nsims)/max
  End Subroutine post_langmuirFit
  !-----------------------------------------------------------------------
  ! We have already readthru the configfile, Now write all the results
  ! All simple averages are written to just one file( outfilename), 
  ! with appropriate comments
  ! Special things( like diffusivity data) are written into separate files 
  ! - Isotherm values are copied from stats to isothemr here ( inside 
  ! loading)
  !------------------------------------------------------------------------
  Subroutine post_writeAvgstoFile(pparams,simno)
    Type(PostInfo), Intent(Inout) :: pparams
    Integer, Intent(in) :: simno
    Integer :: unitno, i, sorb, sorb1, sorb2, nuc, actual_blocks, nsims, sim
    Integer :: site, isotunit, nsorbs, error
    Real(kind=RDbl) :: sorb1n, sorb2n, total_n, cumval, pressure, uc_loading
    Real(kind=RDbl) :: simvol, z_max, err, K
    Real(kind=RDbl),Dimension(4) :: num_density
    Real(kind=RDbl),Dimension(:,:), Pointer :: isoTdat
    Character(len=lstrLen) :: names, name, fmtstring, molecname, isothermfile
    Character(len=9) :: ITstring="isotherm."
    Character(len=14) :: siteITstring="site-isotherm."
    Character(len=15) :: str1, str2
    Logical :: gcmc
    unitno=file_open(pparams%outFilename)
    nuc=simcell_getnuc(pparams%scell) 

    nsorbs=molecules_getnsorbs()
    simvol=simcell_getvolume(pparams%scell)

    !** Write Nrg avgs
    If (Associated(pparams%nrgs)) Then
      !** This data is useful for looking at convergence of the run
      Write(unitno,'(a)') "## Total Energy averages , for the ensemble"
      Write(unitno,'(a)') "## configfile : "//&
          Trim(pparams%configFiles(simno)%name)
      Write(unitno,*) " Iter. No       Inst.(KJ)      Block (KJ)      &
          &Cumul(KJ) "
      ! ** actual number of blocks might be different from that 
      ! ** specified in ctrlfile because of truncation
      actual_blocks=((pparams%lastconf-pparams%firstconf+1)/ &
          pparams%nrgs%statblocksize)
      Do i=1,actual_blocks 
        Write(unitno,'(i10, 3f16.6)') Int(pparams%nrgs%output_array(i,1) ),&
            pparams%nrgs%output_array(i,2),  &
            pparams%nrgs%output_array(i,3), &
            pparams%nrgs%output_array(i,4) 
      End Do

      !** The individual sorb-sorb coul/noncoul/ke/intra nrgs
      Write(unitno,*) 
      If (Associated(pparams%loading)) Then
        Write(unitno,*) "## All the averages below are based on cum avg"
        Write(unitno,*) "## Unstored/not used nrgs might appear as zero"
        Write(unitno,'(2a)') "Sorbs:                                 ",&
            "Nrg-Total        mol/uc    Nrgs(KJ/mol)" 
        Do sorb1=1,nsorbs
          Do sorb2=sorb1,nsorbs
            names=Trim(molecules_name(sorb1))//"--"//&
                Trim(molecules_name(sorb2))

            If (config_isfixed(pparams%species(sorb1))) Then 
              sorb1n=zero
            Else
              sorb1n=stats_getcavg(pparams%loading%totalavgs(sorb1))
            End If

            If (config_isfixed(pparams%species(sorb2))) Then 
              sorb2n=zero
            Else
              sorb2n=stats_getcavg(pparams%loading%totalavgs(sorb2))
            End If

            !** toital number of molecules of sorb1 and sorb2
            total_n=sorb1n+sorb2n
            If (sorb1==sorb2)  total_n=sorb1n

            If (total_n<0.01e-1) Then
              Write(unitno,'(a )') (names(1:20))//" -- No Molecs --"
              Cycle
            Else
              cumval=stats_getcavg(pparams%nrgs%coulNrg(sorb1,sorb2))
              Write(unitno,'(a,3f15.5)') (names(1:20))//" Coulombic : ", &
                  cumval, total_n/nuc, cumval/total_n 
              cumval=stats_getcavg(pparams%nrgs%nonCoulNrg(sorb1,sorb2))
              Write(unitno,'(a,3f15.5)') (names(1:20))//" NonCoulom : ", &
                  cumval, total_n/nuc, cumval/total_n 
            End If
          End Do
        End Do

        Do sorb=1,nsorbs
          name=Trim(molecules_name(sorb))

          If (config_isfixed(pparams%species(sorb))) Then 
            total_n = zero
          Else
            total_n = stats_getcavg(pparams%loading%totalavgs(sorb))
          End If
          If (total_n<0.01e-1) Then
            Write(unitno,'(a )') (name(1:20))//" -- No Molecs --"
            Cycle              
          Else
            cumval=stats_getcavg(pparams%nrgs%ke(sorb))
            Write(unitno,'(a,3f15.5)') (name(1:20))//"  Kintetic : ", &
                cumval, total_n/nuc, cumval/total_n 

            cumval=stats_getcavg(pparams%nrgs%intra(sorb,TOTAL_INDEX))
            Write(unitno,'(a,3f15.5)') (name(1:20))//" Intra(Tot): ", &
                cumval, total_n/nuc, cumval/total_n 
          End If
        End Do

      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "If you want molar energies, &
            & specify loading avg tag in the control file."
      End If
    End If

    !** Write the loading averages
    If (Associated(pparams%loading)) Then
      Write(unitno, *) 


      !** get total number of unitcells
      Do sorb=1,nsorbs

        Write(unitno,'(a)') "## Per Unitcell Loading averages , &
            & for the ensemble"
        Write(unitno,'(a)') "## configfile : "//&
            Trim(pparams%configFiles(simno)%name)
        Write(unitno,*) " Iter. No    Inst.(molec/uc)   Block (molec/uc)   &
            &Cumul(molec/uc "
        Write(unitno, '(a)') " Molecule : "//Trim(molecules_name(sorb))
        If (config_isfixed(pparams%species(sorb))) Then
          Write(unitno, *) "FIXED CONFIG, SO NO LOADING AVG REPORTED FOR &
              &THIS SORB"
          Write(unitno, *) 
        Else

          ! ** actual number of blocks might be different from that 
          ! ** specified in ctrlfile because of truncation
          actual_blocks=((pparams%lastconf-pparams%firstconf+1)/ &
              pparams%loading%statblocksize)
          Do i=1,actual_blocks 
            Write(unitno,'(i10, 3f16.6)') &
                Int(pparams%loading%output_array(sorb,i,1) ),&
                pparams%loading%output_array(sorb,i,2)/nuc,  &
                pparams%loading%output_array(sorb,i,3)/nuc, &
                pparams%loading%output_array(sorb,i,4)/nuc 
          End Do !** blocks loop

          ! copy averages to isotherm
          ! pressure is obtained from datafile. Its the pressure of first element
          ! uc_loaing is the last cum value stored.
          ! these values are writte at the end of this routine
          pressure=datafile_getextrainfo(pparams%configfiles(simno), &
              'pressure')
          uc_loading= pparams%loading%output_array(sorb,actual_blocks,4)/nuc 
          pparams%loading%isotherm(simno, sorb, 1) = pressure
          pparams%loading%isotherm(simno, sorb, 2) = uc_loading
          num_density(1:nsorbs)=&
              pparams%loading%output_array(1:nsorbs,actual_blocks,4)/simvol
        End If
      End Do !** sorb loop
      nsims=(pparams%last-pparams%first+1)
      !** Output the isotherm values during the last output
      If (simno==nsims) Then
        ! allocate for copying isotherm
        Allocate(isoTdat(2,nsims),STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

        Do sorb=1,nsorbs
          molecname=Trim(molecules_name(sorb))
          isothermfile=ITstring//Trim(molecname)
          isotunit=file_open(isothermfile)
          Write(isotunit,*) "#------------OVERALL ISOTHERM-------------"
          Write(isotunit, '(a)') "# Molecule : "//Trim(molecules_name(sorb))
          If (config_isfixed(pparams%species(sorb))) Then
            Write(isotunit, *) "FIXED CONFIG, SO NO LOADING AVG REPORTED FOR &
                &THIS SORB"
            Write(isotunit, *) 
          Else
            Write(isotunit,*) "#  Pressure                  loading   "
            Write(isotunit,*) "#  kpa(1 st gcmc species)    (molec/uc)"

            Do sim=1,nsims
              isoTdat(1,sim)= pparams%loading%isotherm(sim,sorb,1)
              isoTdat(2,sim)= pparams%loading%isotherm(sim,sorb,2)
              str1=real2str( isoTdat(1,sim), 15)
              Write(isotunit,'(a,10x,1f12.7)') str1, isoTdat(2,sim)
            End Do
            ! do a langmuir fit always
            gcmc=((Trim(pparams%simType)=="GCMC").or.&
                (Trim(pparams%simType)=="HGCMC"))
            If ((nsims>1).and.(gcmc)) Then
              Call post_langmuirFit(isoTdat, z_max, K, err)
              Write(isotunit, '(a)')"# Results of a Langmuir Fit, &
                  &fug= K * z / (z_max -z)"
              Write(isotunit, '(a)')"# Where K is is kPa, z , and z_max &
                  & are zeo-loadings in molec/uc"
              str1=real2str( K, 15)
              Write(isotunit, '(3a,f16.5)')"# K =",str1, " ; z_max", z_max
              Write(isotunit, '(a,f16.5)')"# % Error in loading ", err
              If (err>5.0) Then
                Write(isotunit, '(a)')"# This high error could be from &
                    &statistically wrong valuesat low pressure values"
                Write(isotunit, '(a)')"# skip them and try again" 
              Endif
            Endif



          Endif
          Write(isotunit, *) 
          close(isotunit)
        End Do

      Endif



    End If

    !** Write the Siting averages
    If (Associated(pparams%siting)) Then
      Write(unitno, *) "       -------------------          "
      Do sorb=1,nsorbs
        If (config_isfixed(pparams%species(sorb))) Cycle

        Write(unitno,'(a)') "## Per Unitcell SITING averages for the ensemble"
        Write(unitno,'(a)') "## configfile : "//&
            Trim(pparams%configFiles(simno)%name)
        Write(unitno, '(a)') " Molecule : "//Trim(molecules_name(sorb))
        Write(unitno,*) " ## Site Type    Cumul(molec/uc) "
        Do i=1,pparams%siting%nsites
          Write(unitno,'(i10, f16.6)') pparams%siting%sites(i), &
              stats_getcavg(pparams%siting%siteavgs(i,sorb))
        End Do !** sites loop

        ! ** actual number of blocks might be different from that 
        ! ** specified in ctrlfile because of truncation
        actual_blocks=((pparams%lastconf-pparams%firstconf+1)/ &
            pparams%siting%statblocksize)

        If (Trim(genparams%displaymode)=="VERBOSE") Then
          !** now write the running averages
          Write(unitno,'(a)') " the running averages, per unit cell, per site"
          Do site=1,pparams%siting%nsites
            Write(unitno,'(a,i10)') "Loading in site",  &
                pparams%siting%sites(site)
            Do i=1,actual_blocks 
              Write(unitno,'(i10, 3f16.6)') &
                  Int(pparams%siting%output_array(site,sorb,i,1) ),&
                  pparams%siting%output_array(site,sorb,i,2)/nuc,  &
                  pparams%siting%output_array(site,sorb,i,3)/nuc, &
                  pparams%siting%output_array(site,sorb,i,4)/nuc 
            End Do !** blocks loop
            Write(unitno,*)

          End Do !** sites loop
        ENdif

        Write(unitno,*) 

        ! copy averages to isotherm, FOR later output of isotherm
        ! pressure is obtained from datafile. Its the pressure of first element
        ! uc_loaing is the last cum value stored.
        ! these values are writte at the end of this routine
        pressure=datafile_getextrainfo(pparams%configfiles(simno), &
            'pressure')
        pparams%siting%isotherm(simno, sorb, 1) = pressure
        pparams%siting%isotherm(simno, sorb, 2) = zero ! tot loading 
        Do site=1, pparams%siting%nsites
          uc_loading= pparams%siting%output_array(site, &
              sorb, actual_blocks, 4)/nuc 

          pparams%siting%isotherm(simno, sorb, 2+site) = uc_loading
          pparams%siting%isotherm(simno, sorb, 2) = &
              pparams%siting%isotherm(simno, sorb, 2) + uc_loading
        End do
      End Do !** sorb loop 

      nsims=(pparams%last-pparams%first+1)
      !** Output the isotherm values during the last output
      If (simno==nsims) Then
        Do sorb=1,nsorbs
          molecname=Trim(molecules_name(sorb))
          isothermfile=siteITstring//Trim(molecname)
          isotunit=file_open(isothermfile)
          Write(isotunit,*) "#------------OVERALL SITING ISOTHERM-------------"

          If (config_isfixed(pparams%species(sorb))) Then
            Cycle
          Else

            Write(isotunit, '(a)') &
                "#  Molecule : "//Trim(molecules_name(sorb))
            Write(isotunit,*) &
                "#  NOTE : Get accurate pressure values from isotherm.* files"

            Write(isotunit,*) "#  Pres  Totalloading  Siteloading(:)  "
            Write(isotunit,*) "#  kpa(1 st gcmc species)    (molec/uc)"
            Write(isotunit,*) "#  sites are : ",pparams%siting%sites
            fmtstring=cleanstring(REALCOLFORMAT(pparams%siting%nsites+2))
            Do sim=1,nsims
              Write(isotunit,fmtstring) &
                  pparams%siting%isotherm(sim,sorb,:)
            End Do
          endif
          Write(isotunit, *) 
          close(isotunit)
        End Do

      Endif




    End If

    !** Write the poremap
    If (Associated(pparams%poremap)) Then
      Call post_writePoreMap(pparams%poremap)
    End If

    !** Write the poremap
    If (Associated(pparams%structure)) Then
      Call post_writeStructureAvgs(pparams, pparams%structure, simno, &
          unitno)
    End If

    If (Associated(pparams%radDist)) Then
      If (Associated(pparams%loading)) Then
        Call post_writeRadDistAvgs(pparams, pparams%radDist, simno, &
            unitno, num_density)
      Else
        Write(*,*) " ---- > add loading section to calculate number &
            &density automatically for radial distrib calcualtions"
        Call post_writeRadDistAvgs(pparams, pparams%radDist, simno, &
            unitno)
      Endif
    Endif
    !** Write the poremap
    If (Associated(pparams%vels)) Then
      Call post_writeVelHists(pparams, pparams%vels, simno, &
          unitno)
    End If

  End Subroutine post_writeAvgstoFile

End Module post

! DMAP_DETAILS
!------------------------------------------------------------------
! Dmaps contain value of a probability density function in a rectangular grid.
! They are usefule for visualizing probability densities, pores, potentials 
! A code written by louis clark can be used for converting them to rasmol 
! readable forme. 
! In this post code the geneerated dmaps contain probability of 
! finding the given sorb at points inside one unit cell of the simcell.
! The format is similar to that of Louis
!**-----------------------------------------------------------------
!  The density map file format is 
!     (unformatted file ):
!     molecule name (Char*25), map type (Integer)
!     grid sizes (nbrx,nbry,nbrz) (Integer,Integer,Integer)
!     the map (Real*4 array(nbrx,nbry,nbrz))
!
!  The map types are:
!     1) potential map (arbitrary units)
!     2) probability map 
!     3) Log(prob) map (used for storing no's that would overrun real*4)
!  Written by Louis Clark (1997-1998)
!**-------------------------------------------------------------------------
! we need type= 2 for this postcode
!------------------------------------------------------------------

!--------------------------------------------------------------------------
! POREMAP_DETAILS
!------------------------------------------------------------------
! The simcell is divided into a cubic grid. All points in the grid that &
! were atleast visited once by a given sorbate is written to the pore map. 
! This basically shows all accessible points inside a zeolite. 
! Can be directly visualized using rasmol
!**-----------------------------------------------------------------
!  The pore map file format is 
!     (formatted file [text file , rasmol *.xyz format] ) :-
!     number of points{Integer}
!     molecule name+system description (Character(len=strLen)
!     x,y,z (Real(3))
!------------------------------------------------------------------


