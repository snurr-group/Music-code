!-------------------------------------------------------------
! This module contains the different data structures and 
! routines for doing the analysis of results of a HYBRID Run,
! There are many general calculations here that can be used by any 
! music postcode
!-------------------------------------------------------------
Module hybridpc

  Use defaults,    Only : strLen, RDbl, RSgl, MAX_SORBS, MAX_ATOMS, &
      MAX_SITES, MAX_DAT_FILES, d_ctrl_file, zero, one, INITIAL_SIZE
  Use utils,       Only : filesrchstr, genfilename, stripcmnt, split, toint, &
      maxn, toupper, isfileopen, allocErrDisplay
  Use vector,      Only : VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag, vector_getcom
  Use file,        Only :  file_open, file_getunit, file_settag
  Use histogram,   Only :  Histogram_Params, histogram_init, histogram_display
  Use multihist,   Only :  MultiHist_Params, multihist_init, &
      multihist_normalize, multihist_getnbins, multihist_update
  Use stats,       Only :  Statistics,stats_init,stats_update, &
      stats_getblock, stats_getcavg, stats_getstd
      
  Use general,     Only :  general_getcontenttag, genparams, &
      general_getnoofconfigs, general_init, general_display
  Use config,      Only :  AtMolCoords, config_defaultinit, config_display, &
      config_getnatoms, config_initdefaultsize, config_allocfields, &
      config_checkandincr, config_isfixed,default_config_tag,&
      config_initFixed,config_initmolecule
  Use datafile,    Only :  CONFILE, datafile_initin, datafile_gencontrolfile, &
      datafile_readconfig
  Use forcefield,  Only :  forcefield_init, forcefield_display  
  Use smap,        Only :  Smap_Params, smap_init, smap_getPosSiteType
  Use molecules,   Only : molecules_getnsorbs, molecules_gettype, &
      molecules_getnatoms, molecules_name, molecules_init, molecules_display, &
      molecs 
  Use molecule,    Only : MolecularParams, molecule_getname 
  Use atom,        Only : atom_init, atom_display
  Use simcell,     Only :  SimCell_Params, simcell_init, simcell_display, &
      simcell_maptouc, simcell_getell, simcell_getzeoell, simcell_geteff
  Use post,        Only : post_getbondlengths, post_getbondangles

  Implicit None
  Save

  Public 
 
  !** The tags for marking different sections of the control file
  Character(len=strLen), Parameter    :: default_hybridpc_tag = &
      "HYBRID Post Code Information"
  Character(len=strLen), Parameter    :: default_loading_tag = &
      "Loading/Energy Information"
  Character(len=strLen), Parameter    :: default_smap_tag = &
      "Site Map Information"
  Character(len=strLen), Parameter    :: default_rasmol_tag = &
      "Rasmol Pore Map Information"
  Character(len=strLen), Parameter    :: default_blength_tag = &
      "Bond Length Processing info"
  Character(len=strLen), Parameter    :: default_bangle_tag = &
      "Bond Angle Processing info"

  !** variables related to intramolecular-stuff
  Type HYBRIDPCIntraVars
    ! Num of sorbs for which length shud be analysed, angl shud be  anlzd
    Integer                                 :: lengthsorbs,anglesorbs
    
    ! type of each sorb
    Integer,Dimension(:),Pointer            :: len_sorbtype,angl_sorbtype
    
    ! name of each sorb
    Character(len=strLen),Dimension(:),Pointer :: len_sorbnames,angl_sorbnames
    Integer,Dimension(MAX_SORBS)           :: noofbonds,noofangles
    
    ! type of the "bond"/"angl" for each sorb
    Integer,Dimension(MAX_SORBS,MAX_ATOMS) :: angl_types
    
    Integer,Dimension(MAX_SORBS,MAX_ATOMS,2) :: bond_types
    
    
    ! Histogram variables for lengths angles
    Type(Histogram_Params),Dimension(:,:),Pointer    :: avgLengths,avgAngles
    
    ! min and max value. Histigrams are made only between these values
    Real(kind=RDbl) :: minbondlength,maxbondlength,minbondangle,maxbondangle
    
    ! Number of bins into which the interval (max-min) is divided
    Integer :: noofBondBins,noofAngleBins
    Character(len=strLen)         :: anglOutFile,lenOutFile
    
  End Type HYBRIDPCIntraVars
  



  !** stuff to be read from postcode-control file
  Type HYBRIDPCgeneral
    Integer                       :: N0   !** serial number of first datafile
    Integer                       :: NN   !** serial number of last datafile
    Character(len=strLen)         :: basename
    
    ! Name of new ctrlfile that will be regenerated drom the datafile
    Character(len=strLen)         :: new_ctrlfile
    Real(kind=RDbl)       :: startskip !** % of initial data to be skipped
    Real(kind=RDbl)       :: endskip   !** % of final data to be skipped
    
    ! Hacked up to take care of writes to datafile during "moves"
    Integer                       :: iters,configStep
    
    ! Number of blocks to which the datafile should be divided during
    ! loading/energy averages
    Integer                       :: nblocks
  End Type HYBRIDPCgeneral
  


  !** Variables for map-related stuff
  Type HYBRIDPCMapStuff
    Type(Smap_Params), Pointer     :: smap
    Character(len=strLen)          :: smapname
    Integer                        :: n_sites !** number of sites in smap
    Integer, Dimension(MAX_SITES )  :: sitetype !** Types of sites
    Real(kind=RDbl), Dimension(MAX_SITES,MAX_SORBS)  :: site_loading 
    Integer :: n_dmaps  ! number of dmaps to be generated
    Type(Multihist_Params),Dimension(:),Pointer          :: dmap
    
    ! Name of dmap, name of sorb for which dmap should be made
    Character(len=strLen),Dimension(:),Pointer :: dmapName,dmapSorb
    
  End Type HYBRIDPCMapStuff


  !** Variables for generating a rasmol file for all points visited by 
  !** a sorbate 
  Type HYBRIDPCrasmol
    Character(len=strLen)          :: basename !basename
    Integer                        :: nx,ny,nz  !** number of divisions
    Integer                        :: nucx,nucy,nucz  !** number of unitcells
    Real(kind=RDbl)                :: dx,dy,dz  !** unitcell divisions
    Integer, Dimension(:,:,:,:), Pointer  :: fillArray !** filling data
    Character(len=strLen),Dimension(:),Pointer :: rasSorbs !** sorbates
    Integer                        :: nsorbs, ncont  !** ncont= no of contours
    Real(kind=RDbl)                :: scaling_factor
  End Type HYBRIDPCrasmol
  

  !** temp. variables which pertain to averaging of one datafile
  Type HYBRIDPCdat
    !variables for averaging among all the data in the datafile
    Type(Statistics),Dimension(MAX_SORBS)::nmole_list
    Type(Statistics),Dimension(MAX_SORBS,MAX_SORBS,2)::nrg_list
    Type(Statistics) ::energy,p_energy
    
    !** variables for finding standard devn and averages among the blocks
    !**these are for checking convergence
    Type(Statistics),Dimension(MAX_SORBS)::nmoleblock
    Type(Statistics) ::nrgblock
  End Type HYBRIDPCdat
  

  !** The combination of all of the above types
  Type HYBRIDPcParams
    Type(HYBRIDPCgeneral)           :: gen
    Type(HYBRIDPCdat)               :: dat
    Type(HYBRIDPCMapStuff),Pointer  :: map
    Type(HYBRIDPCRasmol),Pointer    :: rasmol
    Type(HYBRIDPCIntraVars)  :: intra
    Type(AtMolCoords), Dimension(:), Pointer :: sorbates
    Type(SimCell_Params)            :: scell 
    
    Character(len=strLen)         :: hybrid_tag
    Character(len=strlen):: commentString
    Integer :: nsims,content_tag
    
    Integer,Dimension(MAX_SORBS):: molec_natoms
    Character(len=strlen),Dimension(MAX_SORBS):: molecnames
    Integer:: nads,nsorbates
    
    !list of datafiles to be processed
    Type(CONFILE),Dimension(MAX_DAT_FILES) :: dflist
    
    ! contains name of all output files for loading averages
    Character(len=strlen),Dimension(MAX_DAT_FILES,MAX_SORBS)   :: outfile 
    
    ! contains name of all output files for energy averages
    Character(len=strlen),Dimension(MAX_DAT_FILES)             :: en_file
    Character(len=2*strlen),Dimension(MAX_SORBS)             :: IsoThermLoad
    Character(len=2*strlen)             :: IsoThermNrg
    Integer :: iters,configStep
  End Type HYBRIDPcParams
  
  Type(HYBRIDPcParams):: hpost
  
Contains


  !----------------------------------------------------------------
  !This subroutine does the necessary initialisations
  !for reading what is written in the file - config_filename
  !It reads the header, generates the original control file that 
  !was used for making the given  config file,and initialises all structures
  !---------------------------------------------------------------
  Subroutine hybridpc_init(postfile)
    Character(len=strlen),Intent(in)::postfile
    Integer :: postUnit
    
    hpost%hybrid_tag= default_hybridpc_tag
    !** read info from post code-control file
    Call hybridpc_readPostFile(postfile,postUnit,hpost%hybrid_tag)
    
    !** Tells what is the format of datafile
    hpost%content_tag=general_getContentTag()
    
    !** get general details so that they can be used anywhere within module
    Call hybridpc_getMolecDetails(hpost%nsorbates,hpost%nads &
        ,hpost%molecnames,hpost%molec_natoms)
    
  End Subroutine hybridpc_init

  
  !------------------------------------------
  ! Initialises atoms,molecules,simcell,zeolite et
  !------------------------------------------
  Subroutine hybridpc_initGeneral()
    Integer :: dUnit=6
    
    Call file_settag(hpost%gen%new_ctrlfile, d_ctrl_file)
    Write(*,*) 
    Write(*,*) "Reading from the re-generated control file"
    Write(*,*) 
    Call general_init(hpost%gen%new_ctrlfile)
    Write(*,*) "##################################################"
    Write(*,*) "THE GENERAL PARAMS USED FOR THE SIMULATION"
    Call general_display(dUnit)
    Write(*,*)
    
    Call atom_init()
    Write(*,*) "##################################################"
    Write(*,*) "THE ATOMS USED FOR THE SIMULATION"
    Call atom_display(dUnit)
    
    Write(*,*) 
    Call molecules_init()
    Write(*,*) "##################################################"
    Write(*,*) "THE MOLECULES USED FOR THE SIMULATION"
    Call molecules_display(dUnit)
    
    Write(*,*) 
!    Call zeolite_init()
    Call simcell_init(hpost%scell, hpost%gen%new_ctrlfile)

    Write(*,*) "##################################################"
    Write(*,*) "THE SIMULATION-CELL VARIABLES USED FOR THE SIMULATION"
    Call simcell_display(hpost%scell,2,dUnit)

    !** We dont want to depend on the original run's initialisation here
    !** Probably the restart file might not exist, So allocate some default 
    !** using config_initdefaultsize
     Call hybridpc_configInit(hpost%sorbates)


    Write(*,*) "##################################################"
    Write(*,*) "THE INITIAL CONFIGURATION USED FOR THE SIMULATION"
    Call config_display(hpost%sorbates)
    
    Call forcefield_init(hpost%gen%new_ctrlfile, hpost%scell)
    Call forcefield_display(hpost%sorbates, hpost%scell)
    Write(*,*) 
    
  End Subroutine hybridpc_initGeneral


  !--------------------------------------------------------------------
  ! Initializes the "sorbates" data structure for "nsorbates" given
  ! the "ctrl_file" and the optional control file tag "optconfig_tag"
  !--------------------------------------------------------------------
  Subroutine hybridpc_configInit(sorbates)
    Type(AtMolCoords), Dimension(:), Pointer :: sorbates
    Character(len=strLen)  :: srchstr, sourcetype, filename, tag, sorbname
    Integer    :: nsorbates
    Integer    :: sorb
    Integer    :: unitno, lineno, sorbsread, error, ios, i, j
    Real(kind=RDbl), Dimension(Size(sorbates,1)) :: sznrg

      tag = default_config_tag

    !** Open the ctrl_file if it is not opened
    unitno = file_open(hpost%gen%new_ctrlfile)

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
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, & 
          " Could not allocate 'sorbates'"
      Stop
    End If

    !** Initialize the sourcetype field
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
      Endif
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
      !** we would prefer not to use reastart files etc
      Select Case (toupper(Trim(sourcetype)))
      Case ("EQUILFILE")
        Call config_initdefaultsize(hpost%sorbates, sorb)
      Case ("RESTARTFILE")
        Call config_initdefaultsize(hpost%sorbates, sorb)
      Case ("CRASHFILE")
        Call config_initdefaultsize(hpost%sorbates, sorb)
      Case ("MEPRESTARTFILE")
        Call config_initdefaultsize(hpost%sorbates, sorb)
      Case ("CUSTOM")
        Call config_initdefaultsize(hpost%sorbates, sorb)
      Case ("GCMC")
        Call config_initdefaultsize(sorbates, sorb)
      Case ("DFLT_SIZE")
        Call config_initdefaultsize(hpost%sorbates, sorb)
      Case ("MOLECULE")
        !Just read the coordinates from molecule
        Call config_initmolecule(sorbates, sorb)
      Case ("CCBGCMC")
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(1x,2a)')  &
            "Amit is a hacker."
        Stop
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

  End Subroutine hybridpc_configInit


  !-----------------------------------------------------------------
  ! Reads the section containing specifics of the HYBRID Post code,
  !and initialises the "var" part of post
  !Also sets the array dflist, which contains the name of files to be analysed
  !-----------------------------------------------------------------
  Subroutine hybridpc_readPostFile(postfile,postUnit,opt_tag)
    Character(*), Intent(in)  :: postfile
    Integer ,Intent(out) :: postUnit
    Character(*), Optional, Intent(in)  :: opt_tag
    Character(len=strLen),Dimension(20)          :: fields
    Type(CONFILE)                 :: firstconfig
    
    Character(len=strlen)::basename,angleSorbName
    Character(len=strLen)   :: tag, line
    Integer                 :: unitno, lineno,i,nfiles,k,j
    
    !** open postfile if not already opened
    unitno=file_open(postfile)
    
    If (Present(opt_tag)) Then
      tag = opt_tag
    Else
      tag = default_hybridpc_tag
    End If
    
    !** Find the HYBRID post code tag
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif
    
    !** Read the general hybrid post code stuff
    Write(*,*)
    Write(*,*) "Reading the post code control file"
    Read(unitno,*) hpost%gen%N0
    Read(unitno,*) hpost%gen%NN
    Read(unitno,*) hpost%gen%new_ctrlfile
    Read(unitno,*) hpost%gen%basename

    !*** generates the New Control File from the given config file
    firstconfig%name=genfilename(hpost%gen%basename, hpost%gen%N0)
    Write(*,*) 
    Write(*,*) "Reading the first config file("//Trim(firstconfig%name)//")"
    Write(*,*) " -and re-generating the control-file used for the simulation"
    Call datafile_initin(firstconfig,0)
    Call datafile_gencontrolfile(hpost%gen%new_ctrlfile,firstconfig)
    Close(firstconfig%unitno)
    
    Read(unitno,*) hpost%gen%startskip
    Read(unitno,*) hpost%gen%endskip
    If ( (hpost%gen%startskip>100).Or.(hpost%gen%startskip<0).Or. &
        (hpost%gen%endskip>100).Or.(hpost%gen%endskip<0).Or.      &
        ((hpost%gen%endskip+hpost%gen%startskip)>100) )         Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,&
          "      Wrong value of endskip/startskip specified in control file"
      Stop
    Endif

    Read(unitno,*) hpost%commentString
    hpost%commentString=Trim(stripcmnt(hpost%commentString))//"#"    

    
    !** stuff should be got from gcmc params, ****** HACK **********
    !** but because of compiler problems I cant do it
    !** I might have to make a pseudogcmc_init and copy it int his module
    !** so that gcmc-details can be got from the control file 
    Read(unitno,*) hpost%gen%iters
    Read(unitno,*) hpost%gen%configStep
    
    Read(unitno,*)    
    
    !****  loading/energy stuff
    Read(unitno,*) hpost%gen%nblocks
    
    !** Use the molecules /atoms etc which the original run had used
    Call hybridpc_initGeneral()
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    !**** Initialise Site Map Information
    Call hybridpc_initSitemap(unitno)
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    !**** Initialise rasmol file information
    !**** rasmol file shows all points visited by a particular sorbate
    !**** will be usefule for visualising cavities, pores etc
    Call hybridpc_initRasmol(unitno)

    !**** Initialise Bond length averaging section
    Call hybridpc_initBLength(unitno)
    Write(*,*) "BONDS IS", hpost%intra%lengthsorbs
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__    
    !**** Initialise Bond Angle averaging section
    Call hybridpc_initBangle(unitno)
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__    
    basename=hpost%gen%basename
    nfiles=hybridpc_GetNoOfFiles()
    !*** Fixes the filename of config files to be processed
    Do i=1,nfiles
      If (nfiles > MAX_DAT_FILES) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Maximum no. of config files exceeded "
        Stop
      End If
      hpost%dflist(i)%name=genfilename(basename, hpost%gen%N0+i-1)
    End Do
    postUnit=unitno
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
  End Subroutine hybridpc_readPostFile
  
  !-----------------------------------------------------------------
  ! Reads the section of post-code control fil
  ! containing specifics of sitemaps
  !-----------------------------------------------------------------  
  Subroutine hybridpc_initBLength(unitno)
    Integer ,Intent(in) :: unitno
    Character(len=2*strLen)   :: tag, line
    Character(len=strLen)   :: lengthtag,tempstr
    Character(len=strLen),Dimension(10)   :: fields
    Character(len=strLen),Dimension(2)    :: atomnum
    Integer :: lineno,i,j,k,error,nsorbs,nbonds,natoms,sorbtype
    Integer :: maxBonds


    Rewind(unitno)

    !** Locate the bond-length section in control file
    tag=default_blength_tag
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(*,*) " NO Bond Length Processing Specified "
      hpost%intra%lengthsorbs=0
      Return
    Endif
    
    !** number of sorbates-for which blength analysis is rqd.
    Read(unitno,*) nsorbs
    hpost%intra%lengthsorbs = nsorbs
    Allocate(hpost%intra%len_sorbtype(nsorbs),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(hpost%intra%len_sorbnames(nsorbs),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
    Do i=1,nsorbs
      Read(unitno,'(a)') line
      hpost%intra%len_sorbnames(i)=Trim(stripcmnt(line))
      sorbtype=molecules_gettype( hpost%intra%len_sorbnames(i) )
      If(sorbtype<=0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__ ,&
            "Wrong Molecule Name"//hpost%intra%len_sorbnames(i) 
        Stop
      Else
        hpost%intra%len_sorbtype(i)=sorbtype
      Endif
      
      Read(unitno,'(a)') line
      line=Trim(stripcmnt(line))
      nbonds= split(line, fields, ",")
      Do j=1,nbonds
        natoms= split( fields(j), atomnum, "-" )
        hpost%intra%bond_types(i,j,1)=toint( atomnum(1) )
        hpost%intra%bond_types(i,j,2)=toint( atomnum(2) )
        Write(*,*) " Atom numbers :",atomnum(1:2)
      End Do
      hpost%intra%noofbonds(i)=nbonds
    End Do

    Write(*,*) "BONDS IS",nbonds
    
    maxBonds=maxn( hpost%intra%noofbonds(1:nsorbs) )
    Write(*,*) "Maximum DEBUG no of bonds is ",maxBonds
    Allocate(hpost%intra%avgLengths(nsorbs,maxBonds),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    !** histogram stuff
    Read(unitno,*) hpost%intra%noofBondBins
    Read(unitno,*) hpost%intra%minbondlength,hpost%intra%maxbondlength
    Do i=1,nsorbs
      Do j=1,hpost%intra%noofbonds(i)
    tempstr=" "
    Write(tempstr,'(a,i4,a,i4)') "sorb-",i," bond-",j
    Write(*,*) "tempstr",tempstr

    lengthtag=" Lengths : "
    lengthtag=Trim(lengthtag)//Trim(tempstr)
    Write(*,*) "ltag",lengthtag
        Call histogram_init( hpost%intra%avgLengths(i,j), &
            hpost%intra%noofBondBins, hpost%intra%minbondlength, &
            hpost%intra%maxbondlength,lengthtag )
        End do
      End  Do
      
      Read(unitno,*) line
      hpost%intra%lenOutFile=Trim(stripcmnt(line))
      
    End Subroutine hybridpc_initBLength
    
    
    !-----------------------------------------------------------------
    ! Reads the section of post-code control fil
    ! containing specifics of sitemaps
    !-----------------------------------------------------------------  
    Subroutine hybridpc_initBAngle(unitno)
      Integer ,Intent(in) :: unitno
      Character(len=2*strLen)   :: tag, line
      Character(len=strLen),Dimension(10)   :: fields
      Integer :: lineno,i,j,k,error,nsorbs,nanglesPL_1,sorbtype,maxAngles
      
      !** Wasted 45 minutes figuring that this rewind was missing
      ! --poor grad student, 2001 AD.
      Rewind(unitno)
 
      tag=default_bangle_tag
      lineno = filesrchstr(unitno, tag, line)
      If (lineno == 0) Then
        Write(*,*) " NO Bond Angle  Processing Specified "
      hpost%intra%anglesorbs=0
      Return
    Endif
    
    !** number of sorbates
    Read(unitno,*) nsorbs
    hpost%intra%anglesorbs = nsorbs
    Allocate(hpost%intra%angl_sorbtype(nsorbs),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(hpost%intra%angl_sorbnames(nsorbs),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Do i=1,nsorbs
    Read(unitno,'(a)') line
    line=Trim(stripcmnt(line))
    ! Note : nanglesPL_1 = no of angles + 1
    nanglesPl_1 = split(line, fields, ",")
    hpost%intra%angl_sorbnames(i)=fields(1)    
    sorbtype=molecules_gettype(fields(1))
      If(sorbtype<=0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__ ,&
            "Wrong Molecule Name",fields(1)
        Stop
      Else
        hpost%intra%angl_sorbtype(i)=sorbtype
      Endif
      Do j=2,nanglesPL_1
        hpost%intra%angl_types(i,j-1)=toint(fields(j)    )
      End Do
      hpost%intra%noofangles(i)=nanglesPL_1-1
    End Do

    maxAngles=maxn( hpost%intra%noofangles(1:nsorbs) )
    Write(*,*) "Maximum DEBUG no of angles is ",maxAngles
    Allocate(hpost%intra%avgAngles(nsorbs,maxAngles),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    !** histogram stuff
    Read(unitno,*) hpost%intra%noofAngleBins
    Read(unitno,*) hpost%intra%minbondAngle,hpost%intra%maxbondAngle
    Do i=1,nsorbs
      Do j=1,hpost%intra%noofangles(i)
        Call histogram_init( hpost%intra%avgAngles(i,j), &
            hpost%intra%noofAngleBins, hpost%intra%minbondAngle, &
            hpost%intra%maxbondAngle )
        End do
    End  Do

    Read(unitno,*) line
    hpost%intra%anglOutFile=Trim(stripcmnt(line))

  End Subroutine hybridpc_initBAngle


  !-----------------------------------------------------------------
  ! Reads the section of post-code control fil
  ! containing specifics of sitemaps
  !-----------------------------------------------------------------  
  Subroutine hybridpc_initSiteMap(unitno)
    Integer ,Intent(in) :: unitno
    Character(len=strLen)   :: tag, line, sorbname,filename
    Character(len=strLen),Dimension(10)   :: fields
    Integer :: lineno,i,j,k,n_dmaps,error,dim
    Integer, Dimension(3) :: nbins
    Real(kind=RDbl),Dimension(3) :: highvalue,lowvalue
    Logical:: map_yes
    
    If (general_getcontenttag()<3) Then
      Write(*,*) "config file does not contain positions &
          &so no sitemap info can be generated -Sorry"
    Nullify(hpost%map)
      Return
      Stop
    Endif
    
    map_yes=.False.
    Nullify(hpost%map)

    !** Find the SiteMap Tag
    tag=default_smap_tag
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the smap tag ", tag, " in the control file"
      map_yes=.False.
      Return
    Else
      map_yes=.True.
      Allocate(hpost%map,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Read(unitno,*) line
      hpost%map%smapname=Trim(line)
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      If (Trim(toupper(line)) == 'NULL') Then
        hpost%map%n_sites= 0
      Else
        hpost%map%n_sites = split(line, fields, ",")
        Do j=1, hpost%map%n_sites
          hpost%map%sitetype(j) =toint(fields(j))
          Do k=1,MAX_SORBS
            hpost%map%site_loading(j,k)= zero
          End Do
        End Do
      End If
    Endif

    If (toupper(hpost%map%smapname) /= "NULL") Then
      Call smap_init(hpost%map%smap,hpost%scell,hpost%map%smapname)
    Else
      Nullify(hpost%map%smap)
    End If

    If (map_yes) Then

    Read(unitno,*) n_dmaps
    hpost%map%n_dmaps=n_dmaps
    Allocate(hpost%map%dmap(n_dmaps),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(hpost%map%dmapSorb(n_dmaps),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(hpost%map%dmapName(n_dmaps),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i=1,n_dmaps
      Read(unitno,*) sorbname

      hpost%map%dmapSorb(i)=Trim(stripcmnt(sorbname))
      Read(unitno,*) filename
      Read(unitno,*) nbins(1),nbins(2),nbins(3)
      hpost%map%dmapName(i)=Trim(stripcmnt(filename))
      dim=3
      Do j=1,3
        lowvalue(j)=zero
        Write(*,*) hpost%map%smap%ell%comp
        !** The edge length of the unitcell
        highvalue(j)=hpost%map%smap%ell%comp(j)
      End Do
      Call multihist_init(hpost%map%dmap(i),nbins,lowvalue,highvalue,dim)
      Write(*,*) "High Values for dmaps",highvalue(1:3)
      Write(*,*) "low Values for dmaps",lowvalue(1:3)
    End Do
    
    Endif
    
  End Subroutine hybridpc_initSiteMap
  


  !-----------------------------------------------------------------    
  ! Initialise rasmol file information
  ! rasmol file shows all points visited by a particular sorbate
  ! will be usefule for visualising cavities, pores etc
  !-----------------------------------------------------------------
  Subroutine hybridpc_initRasmol(unitno)
    Integer ,Intent(in) :: unitno
    Character(len=strLen)   :: tag, line, sorbname,filename
    Character(len=strLen),Dimension(10)   :: fields
    Integer :: lineno,i,j,k,n_dmaps,error,dim,nx,ny,nz
    Integer :: nsorbs
    Integer, Dimension(3) :: nbins
    Real(kind=RDbl),Dimension(3) :: highvalue,lowvalue,ell
    Logical:: rasmol_yes

    If (general_getcontenttag()<3) Then
      Write(*,*) "config file does not contain positions &
          &so no rasmol poremap info can be generated -Sorry"
         Nullify(hpost%rasmol)
         Return
      Stop
    Endif
    
    
    rasmol_yes=.False.
    Nullify(hpost%rasmol)
    Rewind(unitno)
    
    !** Find the Rasmol-Pore Map Tag
    tag=default_rasmol_tag
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(*,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the rasmol pore map tag: ", tag, &
          " in the control file"
      Write(*,*) " ---So Rasmol Pore Maps will not be generated--- "
      rasmol_yes=.False.
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      rasmol_yes=.True.
      Allocate(hpost%rasmol,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Read(unitno,*) line
      hpost%rasmol%basename =Trim(stripcmnt(line))
      Read(unitno,*) nx,ny,nz
      Read(unitno,*) hpost%rasmol%nucx, hpost%rasmol%nucy, hpost%rasmol%nucz
      hpost%rasmol%nx=nx
      hpost%rasmol%ny=ny
      hpost%rasmol%nz=nz
      !** edge lengths of the unitcell
      eff=simcell_geteff(hpost%scell, .True.)
      Write(*,'(a,3f15.5)') " edge lengths of unit cell ",eff(1),eff(2),eff(3)
      hpost%rasmol%dx=eff(1)/nx
      hpost%rasmol%dy=eff(2)/ny
      hpost%rasmol%dz=eff(3)/nz
      Read(unitno,*) nsorbs
      hpost%rasmol%nsorbs=nsorbs
      
      Allocate(hpost%rasmol%rasSorbs(nsorbs),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"hpost%rasmol")
      
      Allocate(hpost%rasmol%fillArray(nsorbs,nx,ny,nz),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"hpost%rasmol")
      
      !** assign all elements to zero
      hpost%rasmol%fillArray=0
      
      Do i=1,nsorbs
        Read(unitno,'(a)') line
        hpost%rasmol%rasSorbs(i)=Trim(stripcmnt(line))
      End Do
      
      Read(unitno,*) hpost%rasmol%ncont
      !** the factor by which values are scaled while writing to *.xyz
      Read(unitno,*) hpost%rasmol%scaling_factor
      
    Endif

  End Subroutine hybridpc_initRasmol




  !------------------------------------------
  ! Initialises the arrays and filenames for averaging varying quantities
  ! of the run, outfile contains the name of output files for each 
  ! config file and each molecule
  !------------------------------------------
  Subroutine hybridpc_initfilenames()
    Character(len=strlen)                     ::tempstrg
    Integer :: i,j,k,nfiles
    
    nfiles=hybridpc_GetNoOfFiles()
    Do i=1,nfiles
      Do j=1,hpost%nsorbates
        !only molecules with more than 0 atoms are counted. ie things like
        ! zeolites are excluded
        If (hpost%molec_natoms(j)>0) Then
          tempstrg=hpost%dflist(i)%name
          hpost%outfile(i,j)= &
              Trim(tempstrg)//"."//Trim(hpost%molecnames(j))
        Endif
      End Do

    !*** names of the energy file for the i th config file
    tempstrg=genfilename(hpost%gen%basename,hpost%gen%N0+i-1)
    tempstrg=Trim(tempstrg)//".energy"
    hpost%en_file(i)=tempstrg

    End Do

    !*** names of the  file for the Overall Isotherms
      Do j=1,hpost%nsorbates
        !only molecules with more than 0 atoms are counted. ie things like
        ! zeolites are excluded
        If (hpost%molec_natoms(j)>0) Then
          tempstrg=Trim(hpost%gen%basename)//".IsoThermLoad"
          hpost%IsoThermLoad(j)= &
              Trim(tempstrg)//"."//Trim(hpost%molecnames(j))
        Endif
      End Do
    hpost%IsoThermNrg=Trim(hpost%gen%basename)//".IsoThermNrg"
    

  End Subroutine hybridpc_initfilenames
  
  !------------------------------------------------------------------
  ! Reinitialises or Initialises all the stats variable, before
  ! starting to read a new config file
  !------------------------------------------------------------------  
  Subroutine hybridpc_reInitStats()
    Integer::j,blocksize,k
    
    Call hybridpc_getblocksize(blocksize)
    Do j=1,hpost%nsorbates
      Call stats_init(hpost%dat%nmole_list(j),' ',blocksize,'f15.3',.False.)
      Call stats_init(hpost%dat%nmoleblock(j),' ',hpost%gen%nblocks,&
          .False.,'f15.3')
      Do k=1,hpost%nsorbates
        Call stats_init(hpost%dat%nrg_list(j,k,1),' ',blocksize,.False.,'f15.3')
        Call stats_init(hpost%dat%nrg_list(j,k,2),' ',blocksize,.False.,'f15.3')
      End Do
      
    End Do
    Call stats_init(hpost%dat%energy,' ',blocksize,.False.,'e15.4')
    Call stats_init(hpost%dat%p_energy,' ',blocksize,.False.,'e15.4')
    Call stats_init(hpost%dat%nrgblock,' ',hpost%gen%nblocks,.False.,&
        'e15.4')
  End Subroutine hybridpc_reInitStats
  
  
  !------------------------------------------------------------------
  ! Gets the names of molecules, and number of atoms in each of them
  !------------------------------------------------------------------  
  Subroutine hybridpc_getMolecDetails(nsorbates,nads,molecnames,molec_natoms)
    Integer,Intent(out):: nads
    Integer,Intent(out):: nsorbates
    Character(len=strlen),Dimension(:),Intent(inout)  :: molecnames
    Integer,Dimension(:),Intent(inout)  :: molec_natoms
    Integer :: i
    
    nsorbates = molecules_getnsorbs()
    nads=nsorbates
    Do i=1,nsorbates
      molecnames(i)=molecule_getname(molecs%mo(i))
      molec_natoms(i)=config_getnatoms(hpost%sorbates,i)
      If (molec_natoms(i)==0) nads=nads-1
    End Do
  End Subroutine hybridpc_getMolecDetails
  
  !------------------------------------------
  ! Gives the number of files to be processed
  !------------------------------------------
  Integer Function hybridpc_GetNoOfFiles()
    hybridpc_GetNoOfFiles=(hpost%gen%NN-hpost%gen%N0)+1
  End Function hybridpc_GetNoOfFiles

  !---------------------------------------------------------
  !gets the name of configfile to be opened
  !------------------------------------------------------------
  Subroutine hybridpc_getconfigfile(sim_n,configfile)
    Type(confile), Intent(inout)       :: configfile
    Integer,Intent(in)::sim_n
    configfile%name=hpost%dflist(sim_n)%name
  End Subroutine hybridpc_getconfigfile

  !---------------------------------------------------------
  !gets the fugacity used for generating the sim_n th config file
  ! for the given sorbate : sorb
  !------------------------------------------------------------
  Real(kind=RDbl) Function hybridpc_getFugacity(sim_n,sorb)
    Integer,Intent(in)::sim_n,sorb
    Integer :: i,actualIndex,sorbtype
    Character(len=strLen) :: molecname
    
    hybridpc_getFugacity=zero
!    molecname=hpost%molecnames(sorb)
!
!    Do i=1,hpost%nads
!      sorbtype=hpost%hybridparams%gcmcsorbs(i)%sorbtype
!      If (sorbtype==molecules_gettype(molecname)) Then
!        hybridpc_getFugacity= &
!            hpost%hybridparams%gcmcsorbs(i)%fuglist(sim_n)%pressure
!      Endif
!    End Do
  End Function hybridpc_getFugacity

  !---------------------------------------------------------
  ! calculates the blocksize, using given starting point , ending point
  ! and number of blocks.
  !------------------------------------------------------------
  Subroutine hybridpc_getblocksize(blocksize)
    Integer,Intent(out)::blocksize
    Integer            ::ndats,nstart,nend,nblocks

    ndats=hybridpc_ndats()

    !** calculates starting pt and endinf pt, from the given "percentage vals"
    nstart=Int(ndats*(hpost%gen%startskip/100))
    nend=ndats-Int(ndats*(hpost%gen%endskip/100))
    nblocks=hpost%gen%nblocks

    !** truncates the blocksize to an integer
    blocksize=Int((one*(nend-nstart))/nblocks)
    !** prevent division by zero
    If (blocksize<1) blocksize=1
  End Subroutine hybridpc_getblocksize


  !--------------------------------------------------------------------
  !Opens the output files and returns the unit numbers
  ! each sim_n denotes a config file
  !--------------------------------------------------------------------
  Subroutine hybridpc_openoutfiles(sim_n,newunits,en_unit,prefx)  
    Integer,Intent(in) :: sim_n
    Integer,Dimension(:),Intent(out) :: newunits
    Integer,Intent(out) :: en_unit
    Character(*),Intent(in),Optional:: prefx
    Real (kind=RDbl)    :: pressure

    Character(len=strlen):: comnt1
    Character(len=strlen):: molecname
    Character(len=5):: temp_str
    Integer:: i
    
    If (Present(prefx)) Then 
      comnt1=Trim(prefx)
    Else 
      comnt1="#"
    Endif

    !** opens files for molecules with more than 0 atoms.
    Do i=1,hpost%nsorbates
    pressure = hybridpc_getFugacity(sim_n,i)
      If (hpost%molec_natoms(i)/=0) Then
        newunits(i)=file_getunit(hpost%outfile(sim_n,i))
        Open(file=hpost%outfile(sim_n,i),unit=newunits(i))
     Write(newunits(i),'(a)') Trim(comnt1)//"###############################"
        Write(newunits(i),'(2a)') Trim(comnt1)//"HYBRID Loading: &
            &Number of Molecules vs iteration number For ",&
            Trim(hpost%molecnames(i))
        Write(newunits(i),'(a,f15.5)') Trim(comnt1)//"Pressure  : ",pressure
        Write(newunits(i),'(a)') Trim(comnt1)//&
     "   Iter.No          Inst.  Block.Average  Cumul.Average         ST.Dev"
      Endif
    End Do
    
    !*** opens the energy file for the sim_n th config file
    en_unit=file_getunit(hpost%en_file(sim_n))
    Open(file=hpost%en_file(sim_n),unit=en_unit)
    Write(en_unit,'(a)') Trim(comnt1)//"##################################"
    Write(en_unit,'(a)') Trim(comnt1)//"HYBRID data :Energy vs &
        &Iteration Number "
    Write(en_unit,'(a)') Trim(comnt1)//&
    "   Iter.No          Inst.  Block.Average  Cumul.Average         ST.Dev"

    
  End Subroutine hybridpc_openoutfiles
  
  !-----------------------------------------
  ! Writes the averaged data( no of molecs) into the files, so that 
  ! they can be plotted and viewed
  !-----------------------------------------
  Subroutine  hybridpc_writeoutfiles(units,en_unit,itn)
    Integer,Dimension(MAX_SORBS),Intent(in):: units
    Integer,Intent(in):: en_unit
    Integer,Intent(in):: itn

    Integer :: i

    Do i=1,hpost%nsorbates
      If (hpost%molec_natoms(i)>0) Then
        Write(units(i),'(i10,4f15.3)') itn,hpost%dat%nmole_list(i)%inst, &
            hpost%dat%nmole_list(i)%blockavg, &
            hpost%dat%nmole_list(i)%cumavg !, &
            !Sqrt(hpost%dat%nmole_list(i)%cumdevsq)
        
      Endif
    End Do
    Write(en_unit,'(i10,4f15.3)') itn,hpost%dat%energy%inst, &
        hpost%dat%energy%blockavg, &
        hpost%dat%energy%cumavg, &
        Sqrt(hpost%dat%energy%cumdevsq)
    
  End Subroutine hybridpc_writeoutfiles
  
  !-----------------------------------------
  !Closes the output files corresponding to sim_n
  !-----------------------------------------
  Subroutine hybridpc_closeoutfiles(newunits,en_unit)  
    Integer,Intent(in) :: en_unit
    Integer,Dimension(:),Intent(in) :: newunits
    
    Integer :: i
    Do i=1,hpost%nsorbates
      If (hpost%molec_natoms(i)>0) Close(unit=newunits(i))
    End Do
    Close(en_unit)
  End Subroutine hybridpc_closeoutfiles
  

  !-----------------------------------------
  !Closes the output files corresponding to sim_n
  !-----------------------------------------
  Subroutine hybridpc_statsupdate(nmole_list,potlist,energy)  
    Integer,Dimension(:),Intent(in) :: nmole_list
    Real(kind=RDbl),Dimension(:,:,:),Intent(in) ::potlist
    Real(kind=RDbl),Intent(in) :: energy
    Real(kind=RDbl) :: p_energy
    Integer :: i,j
    p_energy=zero
    Do i=1,hpost%nsorbates

      !** avoid zeolites 
      If (hpost%molec_natoms(i)>0) Then
        Call stats_update(hpost%dat%nmole_list(i),nmole_list(i)*1.0_RDbl)
      Endif

      Do j=i,hpost%nsorbates
        !** non coul interaction energy
        p_energy=potlist(i,j,1)+p_energy
        Call stats_update(hpost%dat%nrg_list(i,j,1),potlist(i,j,1)*1.0_RDbl)

        !** coul coul interaction energy
        p_energy=potlist(i,j,2)+p_energy
        Call stats_update(hpost%dat%nrg_list(i,j,2),potlist(i,j,2)*1.0_RDbl)
      End Do

    End Do

    !** Update total energies 
    Call stats_update(hpost%dat%energy,energy)
    Call stats_update(hpost%dat%p_energy,p_energy)
    
  End Subroutine hybridpc_statsupdate
  
  !---------------------------------------------------------
  !Opens the Isotherm files, writes their headers
  !------------------------------------------------------------
  Subroutine hybridpc_openIsothermFiles()  
    Character(len=strlen):: comnt1
    Character(len=5):: temp_str
    Integer :: i,en_unit
    Integer,Dimension(MAX_SORBS)::newunits

      comnt1=hpost%commentString

    !** opens files for molecules with more than 0 atoms.
    Do i=1,hpost%nsorbates
      !only molecules with more than 0 atoms are counted. ie things like
      ! zeolites are excluded
      If (hpost%molec_natoms(i)>0) Then
        newunits(i)=file_getunit(hpost%IsoThermLoad(i))
        Open(file=hpost%IsoThermLoad(i),unit=newunits(i))
        Write(newunits(i),'(a)') Trim(comnt1)//"##############################"
        Write(newunits(i),'(2a)') Trim(comnt1)//"Isotherm From HYBRID: &
            &Number of Molecules vs Fugacity ",Trim(hpost%molecnames(i))
        Write(newunits(i),'(a)') Trim(comnt1)//&
            "    Sim.No  Pressure        Loading         Loading/uc     STD(Loading) "
      Endif
    End Do
    !** opens file for writing energies into
    en_unit=file_getunit(hpost%IsoThermNrg)
    Open(file=hpost%IsoThermNrg,unit=en_unit)
    Write(en_unit,'(a)') Trim(comnt1)//"###############################"
    Write(en_unit,'(a)') Trim(comnt1)//"Energies From HYBRID: &
        &Energy vs Fugacity "
    Write(en_unit,'(a)') Trim(comnt1)//&
            "    Sim.No  Total Energy    Energy/Molec    STD(Total) "

  End Subroutine hybridpc_openIsothermFiles


  !-------------------------------------------------------------------------
  ! Initialises all the histogram variables like bond-angles
  ! and bond-lengths
  !-------------------------------------------------------------------------
  Subroutine hybridpc_histogramInit(displacements, initcoords, sorb )
    Real(kind=RDbl),Dimension(:,:,:),Pointer :: displacements
    Type(AtMolCoords),Dimension(:),Pointer  :: initcoords    
    Integer ,Intent(in) :: sorb
    Integer :: natoms,nAngles,maxnatoms, maxnAngles,i,j,k,error
    Logical :: is_bending
    Type(MolecularParams) :: molecptr
    
    maxnatoms=0
    Do i=1,hpost%nsorbates
      natoms = molecules_getnatoms(i)    
      If ( ( .not.config_isfixed(hpost%sorbates(i)) ) .and. (natoms>maxnatoms) ) &
          maxnatoms=natoms
      !** even if one sorbate has bond angles is_bending should be true
    End Do
    
    Allocate ( displacements(hpost%nsorbates,maxnatoms, &
        INITIAL_SIZE),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
    !** initial co-ordiantes array
    Allocate(initcoords(hpost%nsorbates),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
    i=sorb
    !** displace ment -histogram params
    Do j=1,maxnatoms
      Do k=1,INITIAL_SIZE
        displacements(i,j,k)=zero
      End Do
    End Do

    Do i=1,hpost%nsorbates
    Call config_allocfields (initcoords(i),i,INITIAL_SIZE) 
    End do

  End Subroutine hybridpc_histogramInit

  !---------------------------------------------------------
  !Calculates the no-of data to be read from the config-file
  !------------------------------------------------------------
  Integer Function hybridpc_ndats()
    Integer :: iters,configStep
    !** till now we assumed that general takes care of "ndats",
    !** now its going to be done form moves also ---So need some ReCODING Here
    !    since configs are written during integration too
    iters=hpost%gen %iters          ! no of iterations during one integration
    configStep=hpost%gen%configStep ! no. itns between writes to config file
    hybridpc_ndats=( iters / configStep ) * genparams%niterations +  &
        general_getnoofconfigs()
  End Function hybridpc_ndats

  !---------------------------------------------------------
  !Calculates the no-of data to be read from the config-file
  !------------------------------------------------------------
  Subroutine hybridpc_intraAvg(nmole_list)
    Integer ,Dimension(:),Intent(in)::nmole_list
    Integer :: i,j,sorbtype,no_of_bonds,no_of_angles
    Integer,Dimension(10,2) :: bond_types
    Integer,Dimension(10)   :: angle_types

    !** Get Bond-Length average data
    If (hpost%intra%lengthsorbs>0) Then
      Do i=1,hpost%intra%lengthsorbs
        sorbtype=hpost%intra%len_sorbtype(i)
        no_of_bonds=hpost%intra%noofbonds(i)
        If (no_of_bonds>0) Then
          !** hacking here
          If (no_of_bonds>10) Then
            Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__        
            Stop
          Endif
          
          Do j=1,no_of_bonds
            bond_types(j,1)=hpost%intra%bond_types(i,j,1)
            bond_types(j,2)=hpost%intra%bond_types(i,j,2)
          End Do
          
          !** Update avgLengths() histograms 
          Call post_getbondlengths(hpost%sorbates,nmole_list, &
              hpost%intra%avgLengths(i,1:no_of_bonds),&
              hpost%intra%len_sorbtype(i),bond_types,no_of_bonds)
        Endif
        
      End Do
    Endif
    
    !** Get Bond-Angle average data
    If (hpost%intra%anglesorbs>0) Then
      Do i=1,hpost%intra%anglesorbs
        sorbtype=hpost%intra%angl_sorbtype(i)
        no_of_angles=hpost%intra%noofangles(i)
        If (no_of_angles>0) Then
          !** More hacking
          If (no_of_angles>10) Then
            Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__        
            Stop
          Endif
          Do j=1,no_of_angles      
            angle_types=hpost%intra%angl_types(i,1:no_of_angles)
          End Do
          
          !** Update avgLengths() histograms  
          Call post_getbondAngles(hpost%sorbates,nmole_list, &
              hpost%intra%avgAngles(i,1:no_of_angles),&
              hpost%intra%angl_sorbtype(i),angle_types,no_of_angles)
        Endif
      End Do
    Endif
  End Subroutine Hybridpc_intraAvg
  
  
  !---------------------------------------------------------
  !Calculates the no-of data to be read from the config-file
  !------------------------------------------------------------
  Subroutine hybridpc_writeNrgLoad(curr,nstart,newunits,en_unit)
    Integer ,Intent(in) :: nstart,curr,en_unit
    Integer,Dimension(:),Pointer :: newunits
    Integer :: j,blocksize, k,itnScale
    Real(kind=RDbl) :: blockavg
    
    Call hybridpc_getblocksize(blocksize)      
    
    !** write energy/loading to output files
    If (Mod(curr-nstart, blocksize) == 0) Then
      
      itnScale=genparams%iconfig
      
      Call hybridpc_writeoutfiles(newunits,en_unit,curr*itnScale)
      
      !** Update the deviation stats between the block avg values
      Do k=1,hpost%nsorbates
        If (hpost%molec_natoms(k)/=0) Then
          blockavg=stats_getblock(hpost%dat%nmole_list(k))
          Call stats_update(hpost%dat%nmoleblock(k),blockavg)
        Endif
      End Do

      !** block average of pot-rng
      blockavg=stats_getblock(hpost%dat%p_energy)
      Call stats_update(hpost%dat%nrgblock,blockavg)
      
    Endif
  End Subroutine Hybridpc_writeNrgLoad
  
  !---------------------------------------------------------
  !At the end, writes the values to Isotherm Files
  !------------------------------------------------------------
  Subroutine hybridpc_writeIsothermFiles(sim_n)
    Integer,Intent(in) :: sim_n
    Integer :: cavgsum,k,unitno
    Real(kind=RDbl) ::cavg, cavg_perunitcell,std_amongblocks,pressure
    Real(kind=RDbl) ::cumnrg,cumnrg_permolec
    
    cavgsum=0
    Do k=1,hpost%nsorbates
      !** avoid zeolite type of molecules
      If (hpost%molec_natoms(k)/=0) Then
        unitno=file_getunit(hpost%IsoThermLoad(k))
        cavg=stats_getcavg(hpost%dat%nmoleblock(k))
        cavg_perunitcell=cavg/hpost%scell%ncells
        
        ! std of loading average between the block values
        std_amongblocks=    stats_getstd(hpost%dat%nmoleblock(k)) 
        
        !total number of molecules/unit cell for all sorbates
        cavgsum=cavgsum+cavg
        pressure= hybridpc_getFugacity(sim_n,k)
        Write(unitno,'(a,i4,a,f12.5,a,f12.5,a,f12.5,a,f12.5)')  &
            "   ",sim_n,"    ",pressure,"    ",cavg,"    ",&
            cavg_perunitcell,"    ",std_amongblocks
      Endif
    End Do


    !** enrgy values to Isotherm file
    unitno=file_getunit(hpost%IsoThermNrg)
    cumnrg=stats_getcavg(hpost%dat%nrgblock)
    If (cavgsum>0) Then 
      cumnrg_permolec=cumnrg/cavgsum
    Else
      Write(unitno,*) "Number of molecs ZERO"
    Endif
    
    std_amongblocks= stats_getstd(hpost%dat%nrgblock)
    Write(unitno,'(a,i4,a,f12.5,a,f12.5,a,f12.5)')  &
        "   ",sim_n,"    ",cumnrg,"    ",cumnrg_permolec,"    ",&
        std_amongblocks 
  End Subroutine hybridpc_writeIsothermFiles

  !---------------------------------------------------------
  !Writes bond length statistics Bond Angle Statistics etc
  !------------------------------------------------------------
  Subroutine hybridpc_writeIntra()
    Real(kind=RDbl):: temp_mean,temp_mean_sq,std_sq
    Integer :: i,j,unitno,sorbno

    If (hpost%intra%lengthsorbs>0) Then 

      unitno=file_open(hpost%intra%lenOutFile)
      
      Write(unitno,*) " # Bond Length Statistics "
      Write(unitno,*) "##########################"
      Do i=1,hpost%intra%lengthsorbs
        Write(unitno,*) "# Molecule Name : ",Trim(hpost%intra%len_sorbnames(i))
        sorbno=hpost%intra%len_sorbtype(i)
        If ( .Not.config_isfixed( hpost%sorbates(sorbno) )  ) Then
          Do j=1,hpost%intra%noofbonds(i)
            Write(unitno,'(a,2i4)') "# Bond atoms ",&
                hpost%intra%bond_types(i,j,1:2)
            Call histogram_display(hpost%intra%avgLengths(i,j),unitno)
            temp_mean=hpost%intra%avgLengths(i,j)%stat%cumavg
            temp_mean_sq=hpost%intra%avgLengths(i,j)%stat%cumsqavg
            Write(unitno,*) "# avg bond statistics are ----------"
            Write(unitno,*) "# Min Value ",hpost%intra%avglengths(i,j)%min,&
                "Max Value ",  hpost%intra%avglengths(i,j)%max
            Write(unitno,*) "# Mean ",temp_mean
            std_sq = temp_mean_sq - temp_mean * temp_mean 
            If( std_sq>zero) Then
              Write(unitno,*) "# STD ",Sqrt(std_sq)
            Else
              Write(unitno,*) "# STD-Sq/ Wong Value- should be +ve ",std_sq
            Endif
            Write(unitno,*) 
          End Do
        Endif
      End Do
      
      Close(unitno)
    Endif

    
    If (hpost%intra%anglesorbs>0) Then
      unitno=file_open(hpost%intra%anglOutFile)
      Write(unitno,*) "# Bond Angle Statistics "
      Write(unitno,*) "##########################"
      Do i=1,hpost%intra%anglesorbs
        Write(unitno,*) "# Molecule Name :: ",&
            Trim(hpost%intra%angl_sorbnames(i))
        sorbno=hpost%intra%angl_sorbtype(i)
        
        If ( .Not.config_isfixed( hpost%sorbates(sorbno) )  ) Then
          Do j=1,hpost%intra%noofangles(i)
            Write(unitno,*) "# Angle Number :: ",hpost%intra%angl_types(i,j)
            Call histogram_display(hpost%intra%avgAngles(i,j),unitno)
            temp_mean=hpost%intra%avgAngles(i,j)%stat%cumavg
            temp_mean_sq=hpost%intra%avgAngles(i,j)%stat%cumsqavg
            Write(unitno,*) "# avg angle statistics are ----------"
            Write(unitno,*) "# Min Value ",hpost%intra%avgAngles(i,j)%min,&
                "Max Value ",  hpost%intra%avgAngles(i,j)%max
            Write(unitno,*) "# Mean ",temp_mean
            std_sq = temp_mean_sq - temp_mean * temp_mean 
            If( std_sq>zero) Then
              Write(unitno,*) "# STD ",Sqrt(std_sq)
            Else
              Write(unitno,*) "# STD-Sq/ Wong Value- should be +ve ",std_sq
            Endif
            Write(unitno,*) 
          End Do
        Endif
      End Do
      Close(unitno)
    Endif

  End Subroutine Hybridpc_writeIntra
  
  !-----------------------------------------------------------------
  !Writes an xyz file in rasmol format for each sorbate specified 
  ! in the poremap section. If hpost%rasmol%cont > 1 then more than one 
  ! contours will be generated. The xyz file corresponds to all the points
  ! visited by the sorbate during the GCMC simulation
  ! If a particular grid is visted 3-times then its value 
  ! hpost%rasmol%fillArray()=3
  !-----------------------------------------------------------------
  Subroutine hybridpc_writePoreMap()
    Integer ::i,j,k,l,m,i_ucx,i_ucy,i_ucz,x_ind,y_ind,z_ind
    Integer ::ncells,error,unitno
    Real(kind=RDbl) :: max,fac,scale
    Character(len=3) :: num_strg
    Character(len=5) :: atom_name !the symbol to be written to rasmol file
    Character(len=strLen) :: newfilename
    Integer,Dimension(:),Pointer :: contlevels

    !** Allocate the contour level array
    Allocate(contlevels(hpost%rasmol%ncont),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Do i=1,hpost%rasmol%nsorbs
      max=0
      !*** Find maximum
      Do k=1,hpost%rasmol%nx
        Do l=1,hpost%rasmol%ny
          Do m=1,hpost%rasmol%nz
            If (hpost%rasmol%fillArray(i,k,l,m)>max) Then
              max=hpost%rasmol%fillArray(i,k,l,m)
            Endif
          End Do
        End Do
      End Do

      fac=Exp(Log(max*one)/hpost%rasmol%ncont)
      contlevels(1)=1
      If (hpost%rasmol%ncont>1) Then
        Do j=2,hpost%rasmol%ncont
          contlevels(j)= (contlevels(j-1)*fac)+1
        End Do
        
      Endif
      Write(*,*) "contlevels are " , contlevels, "max was", max
      Do j=1,hpost%rasmol%ncont
        
        
        !** Open the disk file to write to
        If (j<1000) Then
          Write(num_strg,'(i3)') j
          num_strg=Adjustl(num_strg)
        Else
          Write(*,*) "Too many contors specified, exiting"
          Stop
        Endif

        newfilename=Trim(hpost%rasmol%basename)//"."//&
            Trim(hpost%rasmol%rasSorbs(i))//"."//Trim(num_strg)
        unitno=file_open(newfilename)
        
        atom_name= "A"//Trim(num_strg)

        ncells=0
        !** check how many cells belong to this pore
        Do k=1,hpost%rasmol%nx
          Do l=1,hpost%rasmol%ny
            Do m=1,hpost%rasmol%nz
              If (hpost%rasmol%fillArray(i,k,l,m) >=contlevels(j)) Then
                !** this is one of the members of the pore map
                ncells=ncells+1
              Endif
            End Do
          End Do
        End Do
        
        Write(*,*) "Total fillings is ",ncells

        !Write number of points in rasmol file, ncells* no of unitcells
        Write(unitno,*) ncells * ( hpost%rasmol%nucx * hpost%rasmol%nucy * &
            hpost%rasmol%nucz )
        
        Write(unitno,*) "some description"
        

        !** if "X" is the position in angstroms, then the *.xyz will have 
        !** "scale * X" as the position values
        scale = hpost%rasmol%scaling_factor

        !** go through it again, this time write the pore map too
        Do k=1,hpost%rasmol%nx
          Do l=1,hpost%rasmol%ny
            Do m=1,hpost%rasmol%nz
 !SDEBUG
 !SDEBUG
!   CORRECT one
!              If (hpost%rasmol%fillArray(i,k,l,m) >= contlevels(j)) Then

!   WRONG one
              If (hpost%rasmol%fillArray(i,k,l,m) < contlevels(j)) Then




                !** this is one of the members of the pore map
                !** They have to be written to all unitcells
                Do i_ucx=1,hpost%rasmol%nucx
                  Do i_ucy=1,hpost%rasmol%nucy
                    Do i_ucz=1,hpost%rasmol%nucz

                      x_ind=(i_ucx-1)*hpost%rasmol%nx + k
                      y_ind=(i_ucy-1)*hpost%rasmol%ny + l
                      z_ind=(i_ucz-1)*hpost%rasmol%nz + m

                      !** factor of 10 used here to make sure rasmol file 
                      !** is scaled well
                      Write(unitno,'(a,3f16.6)') Trim(atom_name), &
                          scale*x_ind*hpost%rasmol%dx, &
                          scale*y_ind*hpost%rasmol%dy, &
                          scale*z_ind*hpost%rasmol%dz 

                    End Do
                  End Do
                End Do
                ! end of unit cells loop

              Endif

              
            End Do
          End Do
        End Do
        !** end of rasmol%fillArray(sorb,:,:,:) loop

        Close(unitno)
      End Do
      !** end of conours loop

    End Do
    !** end of nsorbs loop

  End Subroutine hybridpc_writePoreMap



  !----------------------------------------------------
  ! A subroutine for updating sitemap/dmap stuff using the 
  ! positions in current co-ordinates in sorbates 
  !----------------------------------------------------
  Subroutine hybridpc_mapAvg(oldwt,newwt)
    Real(kind=RDbl),Intent(in) :: oldwt,newwt
    Type(VecType),Dimension(MAX_ATOMS)::xyzs
    Type(VecType)::com,uc_com
    Integer:: sorb,molec,maxMolec,i,site_type,natoms,k,dmapNum
    Logical:: dmap_yes
    Integer,Dimension(hpost%map%n_sites,hpost%nsorbates) :: temp_Loading

    !** Loop over all types of molecules
    Do sorb=1,hpost%nsorbates
      
      !** avoid zeolite type of stuff
      maxMolec=Int(hpost%dat%nmole_list(sorb)%inst)
      natoms=hpost%molec_natoms(sorb)
      If (maxMolec<=0) Cycle
      
      !** flag to see whether dmap analysis is reqd for this sorbate
      dmap_yes=.False.
      Do i=1,hpost%map%n_dmaps
        If (   Trim( toupper(hpost%map%dmapSorb(i)) ) == &
            toupper(molecules_name(sorb))  ) Then
          dmap_yes=.True.
          dmapNum=i
        Endif
     End Do

      !** Loop over all molecules of this sorb
      Do i=1,hpost%map%n_sites
        temp_loading(i,sorb)=0

        Do molec=1,maxMolec
          
          !** positions  of all atoms in xyz co-ords
          xyzs=  hpost%sorbates(sorb)%coords(1:natoms,molec)%rp
          
          !** arithmetic com; not mass weighted
          com=vector_getcom(xyzs(1:natoms))
          
          !** com in unitcell
          uc_com=simcell_maptouc(hpost%scell,com)

          !** write dmap during the first loop in sites inly
          If (dmap_yes.and.(i==1)) &
              Call multihist_update(hpost%map%dmap(dmapNum),uc_com%comp)  
          !** Update the siting averages based on com
          site_type=smap_getPosSiteType(hpost%map%smap,com)
          If (site_type==hpost%map%sitetype(i)) Then
            temp_loading(i,sorb)=temp_loading(i,sorb)+1
          Endif
        End Do  !end of molec-loop

      hpost%map%site_loading(i,sorb) =  ( hpost%map%site_loading(i,sorb)  &
          * oldwt + newwt * temp_loading(i,sorb) )  &
          /(oldwt + newwt)

      End Do  ! end of site-loop

      
    End Do    ! end of sorb-loop

  End Subroutine hybridpc_mapAvg




  !----------------------------------------------------
  ! A subroutine for updating rasmol-pore map info
  ! If the COM of the sorbate lies in a particular tube its 
  ! fallArray value is increased by 1
  !----------------------------------------------------
  Subroutine hybridpc_rasmolAvg()
    Type(VecType)::com,uc_com
    Logical :: poremap_yes
    Integer:: sorb,molec,maxMolec,i,site_type,natoms,k,poreMapNum
    Integer :: x_index, y_index, z_index
    ! for storing positions of atoms in one molecule
    Type(VecType),Dimension(MAX_ATOMS)::xyzs    


    !** Loop over all types of molecules
    Do sorb=1,hpost%nsorbates
      
      !** avoid zeolite type of stuff
      maxMolec=Int(hpost%dat%nmole_list(sorb)%inst)
      natoms=hpost%molec_natoms(sorb)
      If (maxMolec<=0) Cycle
      
      !** flag to see whether dmap analysis is reqd for this sorbate
      poremap_yes=.False.
      Do i=1,hpost%rasmol%nsorbs
        If (   Trim( toupper(hpost%rasmol%rasSorbs(i)) ) == &
            toupper(molecules_name(sorb))  ) Then
          poremap_yes=.True.
          poreMapNum=i
        Endif
      End Do

      !** Loop over all molecules of this sorb
      Do molec=1,maxMolec
          
          !** positions  of all atoms in xyz co-ords
          xyzs=  hpost%sorbates(sorb)%coords(1:natoms,molec)%rp
          
          !** arithmetic com; not mass weighted
          com=vector_getcom(xyzs(1:natoms))
          
          !** com in unitcell
          uc_com=simcell_maptouc(hpost%scell,com)
          x_index = Int(uc_com%comp(1)/hpost%rasmol%dx)+1
          y_index = Int(uc_com%comp(2)/hpost%rasmol%dy)+1
          z_index = Int(uc_com%comp(3)/hpost%rasmol%dz)+1
          hpost%rasmol%fillArray(poreMapNum,x_index,y_index,z_index)= &
              hpost%rasmol%fillArray(poreMapNum,x_index,y_index,z_index)+1
        End Do  !end of molec-loop

    End Do    ! end of sorb-loop

  End Subroutine hybridpc_rasmolAvg





  
  !------------------------------------------------------------------
  ! Write the dmaps to the required file
  ! dmaps contain probability of finding the given sorb in 
  ! one unit cell of the simcell.
  ! The format is sililar to that of Louis. For processing dmaps use 
  ! louies code , which converts dmaps to contours
  !**-----------------------------------------------------------------
  !  The density map file format is (unformatted):
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
  ! we need type= 2 
  !------------------------------------------------------------------
  Subroutine hybridpc_writedmaps()
    Real(kind=RSgl),Dimension(:,:,:),Pointer :: prob
    Integer :: map_type=2,bin_index,nbx,nby,nbz
    Integer :: i,ix,jy,kz,unitno,error
    Real(kind=RDbl) :: sum
    Character(len=25)::molecName
    Character(len=strLen) :: dmapfilename


    Do i=1,hpost%map%n_dmaps
      Call multihist_normalize(hpost%map%dmap(i))
      molecName=Trim(hpost%map%dmapSorb(i))  
      nbx=multihist_getnbins(hpost%map%dmap(i),1)
      nby=multihist_getnbins(hpost%map%dmap(i),2)
      nbz=multihist_getnbins(hpost%map%dmap(i),3)
      Allocate(prob(nbx,nby,nbz),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)  
      bin_index=0
      sum=zero
      Do ix=1,nbx
        Do jy=1,nby
          Do kz=1,nbz
            bin_index=bin_index+1
            prob(ix,jy,kz)=hpost%map%dmap(i)%probability(bin_index)
            sum=sum+prob(ix,jy,kz)
          End Do
        End Do
      End Do

      Write(*,*) "THE Sum in dmaps is ",sum,nbx,nby,nbz,bin_index
      dmapfilename=hpost%map%dmapName(i)
      
      unitno = isfileopen(dmapfilename)
      !** If already not opened, Open it
      If (unitno < 0) Then
        unitno = file_getunit(dmapfilename)
        Open(file=dmapfilename, unit=unitno,form='unformatted')
      Endif
      
      
      Write(unitno) molecname,map_type
      Write(unitno) nbx,nby,nbz
      Write(unitno) prob
      Close(unitno)
      Deallocate(prob,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End Do
  End Subroutine hybridpc_writedmaps
  
  !------------------------------------------------------------------
  !Write the dmaps to the required file
  !------------------------------------------------------------------
  Subroutine hybridpc_writeSiting(sim_n)
    Integer,Intent(in) :: sim_n
    Character(len=strLen) :: sitefile,tempstrg,cmnt
    Integer :: unitno,i,j
    
    cmnt=hpost%commentString
    tempstrg=hpost%dflist(sim_n)%name
    sitefile=Trim( Trim(tempstrg) // ".siting" )
    unitno=file_open(sitefile)
    Write(unitno,*) "Site Type  ,  Loading"
    Do i=1,hpost%nsorbates
      If (.Not.config_isfixed(hpost%sorbates(i))) Then
        Write(unitno,*) Trim(cmnt)//"Molecule_name : "//Trim(molecules_name(i))
        Do j=1,hpost%map%n_sites
          Write(unitno,*) hpost%map%sitetype(j),hpost%map%site_loading(j,i)
        End Do
        Write(unitno,*) "----------------------"
      Endif
    End Do
  End Subroutine hybridpc_writeSiting
  !----------------------------------------------------------------------------
  ! Writes a sample of the post code control file information to unit =unitno
  !----------------------------------------------------------------------------
  Subroutine hybridpc_sampleCF(unitno)
    Integer, Intent(In) :: unitno
    
    Write(unitno,'(a)')" ---  "//Trim(default_hybridpc_tag)//" --- "
    Write(unitno,'(a,t30,a)') 'Integer','# Index Number of first config file'
    Write(unitno,'(a,t30,a)') 'Integer','# Index Number of last config file' 
    Write(unitno,'(a,t30,a)') 'Character', &
        '# Name for new regenerated ctrlfile '
    Write(unitno,'(a,t30,a)') 'Character', '# base name of config files '
    Write(unitno,'(a,t30,a)') 'Real',&
        '# Percentage of initial configs to be skipped '
    Write(unitno,'(a,t30,a)') 'Real',&
        '# Percentage of configs at the end of the run to be skipped '
    Write(unitno,'(a,t30,a)') 'Character', &
        '# Comment character to be used in output files '
    Write(unitno,'(a,t30,a)') 'Integer','# Junk-Hack=100' 
    Write(unitno,'(a,t30,a)') 'Integer','# Junk-Hack=1000'
    Write(unitno,*) ' more optional stuff to follow '
  End Subroutine hybridpc_sampleCF


  !---------------------------------------------------------
  !Calculates the averages in a datafile
  !sim_n= serial number of the run among the runs to be processed
  !------------------------------------------------------------
  Subroutine hybridpc_getaverages(sim_n)
    Integer,Intent(in)           :: sim_n
    
    Type(confile)                :: configfile
    Integer                      :: en_unit,ndats,nstart,nend
    Integer                      :: i,j,k,l,natoms,maxnatoms
    Integer                      :: unitno,error,curr,nmoles_now
    Integer,Dimension(:),Pointer :: newunits
    Integer,Dimension(:),Pointer :: nmole_list
    Real(kind=RDbl)              :: energy,blockavg,cumnrg,cumnrg_permolec
    Real(kind=RDbl),Dimension(MAX_SORBS,MAX_SORBS,2) ::potlist
    Real(kind=RDbl),Dimension(MAX_SORBS) ::intra,ke
    Real(kind=RDbl),Dimension(:,:,:),Pointer :: displacements
    Type(AtMolCoords),Dimension(:),Pointer  :: initcoords    
    Integer :: hack_unit,hack_n
 !SDEBUG
 !SDEBUG
 !   hack_unit=file_open('hack_nrg')
 !SDEBUG 
    !** Total number of data points in the config file( after end of run) 
    ndats=hybridpc_ndats()    
    Write(*,*) "Expected number of data : ", ndats 

    !** Re-Initialise all stat variable for loading/energy
    Call hybridpc_reInitStats()

    !** itn.numbers from which averaging of data should start and end
    !** the values outside this range are neglected
    nstart = ndats * ( hpost%gen%startskip / 100 )
    nend= ndats- Int ( ndats * ( hpost%gen%endskip/ 100) )

    !** from the list of datafiles, get the datafile corresponding to sim_n
    Call hybridpc_getconfigfile(sim_n,configfile)
    
    Write(*,*)
    Write(*,*)"Reading and doing the average for data in configfile - "//&
        Trim(configfile%name)
    Write(*,*)

    !** Temperrory variable for reading nmole_list
    Allocate ( nmole_list(hpost%nsorbates),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** unit numbers for  writing loading averages
    Allocate ( newunits(hpost%nsorbates),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Call hybridpc_openoutfiles(sim_n,newunits,en_unit,hpost%commentString)

    !** initialise the datafile for reading from it
    Call datafile_initin(configfile,1)
    Write(*,'(a,i8)') "Data points expected in config file : ", ndats
    Write(*,'(a,i7,a,i7,a)') "This run will read from point-",nstart-1,&
        " to point-",nend,"including both"

    !*** Main Loop for reading and averaging
    !***
    Do j=1,ndats
      
      !** Check and increment the size of sorbates, If required 
      !** needed for GCMC
      Do l=1,hpost%nsorbates
        If (.Not.config_isfixed(hpost%sorbates(l))) Then
          !use previous value of nmoles
          nmoles_now=Int(hpost%dat%nmole_list(l)%inst)
          Call config_checkandincr(hpost%sorbates, l, nmoles_now)
        Endif
      End Do

      !      Write(*,*) "Reading config -",j
      Call datafile_readconfig(hpost%sorbates, configfile, nmole_list, &
          energy,potlist,intra,ke)
  
      If ( (j>=nstart) .And. (j<=nend) ) Then
          
        !** nmolelist may not be useful if nvt-mc
        Call hybridpc_statsupdate(nmole_list,potlist,energy)
  
        !SDEBUG
        !SDEBUG
  !      Write(hack_unit,*) " CONFIG_NO : ",j
  !      Do hack_n=1,18
  !        Write(hack_unit,'(3f15.5)')hpost%sorbates(1)%coords(hack_n,1)%rp
  !      End Do
  !      Write(hack_unit,*) "INTRA_NRG", intra(1)
        !SDEBUG
        !SDEBUG
   
        !** Do the averaging for intra-molecular variables(angles/lengths)
        Call hybridpc_intraAvg(nmole_list)
        !** Do sitemap analysis/ Dmap analysis
        If (Associated(hpost%map)) Call hybridpc_mapAvg((j-nstart)*one,one)
   
        !** Updaate rasmol-dmap file details
        If (Associated(hpost%rasmol)) Call hybridpc_rasmolAvg()
   
        !** displacment from initial position FOR only hybrid-nvt-mc
        !        If (j==1) Then
        !          Call post_copycoords(hpost%sorbates,initcoords,nmole_list)
        !        Else If (j==ndats) Then
        !        Call post_displacements(hpost%sorbates,initcoords, &
        !            displacements,nmole_list)
        !        Endif
        !** write NRg and Loading to output files/If necessary
        curr=j
   
        Call hybridpc_writeNrgLoad(curr,nstart,newunits,en_unit)
   
      Else
   
        !** avoid overshooting the no of data points to be read
        If (j>nend) Exit
   
      Endif
      
    End Do   !** ndats-loop, over all configurations in the datafile

    !** Close files containing loading/energy averages
    Call hybridpc_closeoutfiles(newunits,en_unit)

    !** writes the loading values for current fugacity into IsothermFile 
    Call hybridpc_writeIsothermFiles(sim_n)

    !** Out put intra analysis results
    Call hybridpc_writeIntra()

    !** Out put siting details
    If (Associated(hpost%map)) Call hybridpc_writeSiting(sim_n)
    
    !** output dmaps
    If (Associated(hpost%map)) Call hybridpc_writedmaps()

    !** Write rasmol file fore pore maps
    If (Associated(hpost%rasmol)) Call hybridpc_writeporemap()

  End Subroutine hybridpc_getaverages
  
  
End Module hybridpc







