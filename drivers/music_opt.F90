!-------------------------------------------------------------------
! This is the driver reading a bunch of data from an xyz file
! and optimizing the force field to match the energies in the xyz file
! xyz file should contain single molecule configs
! 2nd line of each xyz entry should contain energy there (hartree units)
! input  : a control file, and music.xyz
! Calculates \sum(nrg-nrg0)^2 error and optimizes to minimise the error
!-------------------------------------------------------------------

Program music

  Use defaults, Only: strLen, RDbl, d_ctrl_file, d_res_file, pcalls, &
      dashedline, MAX_ROUTINES, one, zero, calToj, STRETCH_INDEX, &
      BENDING_INDEX, TORSION_INDEX, TOTAL_INDEX
  Use utils, Only: genfilename, int2str, real2str, allocerrdisplay,&
      split, toreal, filesrchstr
  Use file, Only: file_getunit, file_settag, file_open, file_close
  Use vector, Only: VecType
  Use commandline, Only : commandline_init
  Use general, Only: genparams, general_init, general_sampleCF
  Use wtime, Only: wtime_init, wtime_display, wtime_starttime, wtime_stoptime
  Use atom, Only: atom_display, atom_init, atom_sampleCF
  Use molecules, Only:molecules_display, molecules_init, &
      molecules_sampleCF, molecules_getnatoms
  Use datafile, Only: CONFILE, datafile_writeconfig, datafile_initout
  Use simcell, Only: SimCell_Params, simcell_display, simcell_init, &
      simcell_sampleCF
  Use config, Only: config_writerestartfile, config_init, config_sampleCF, &
      AtMolCoords, config_display, config_dump, config_displaycoords, &
      config_changerp, config_xyz2config, config_allocfields, &
      config_setnmoles, config_xyz2molec, config_copymolec
  Use interact, Only: Interaction_Model, interact_init, interact_display
  !  Use movie, Only: MovieInfo, movie_init, movie_display, movie_makemovie
  Use stopcriteria, Only: STOP_CHECK, stopcriteria_init, stopcriteria_check
  Use ffcall, Only: ffcall_nrg
  Use forcefield, Only: Forcefield_info
  Use intramolecular, Only: IntramolecularInfo
  Use bsmodel, Only : StretchInfo, bsmodel_setHarK
  Use bbmodel, Only : BendingInfo, bbmodel_setHarK
  Use tormodel, Only : TorsionInfo, tormodel_setA1A2
  Use random, Only: rranf
  Implicit None

  Type(StretchInfo),Pointer :: bspntr
  Type(BendingInfo),Pointer :: bbpntr
  Type(TorsionInfo),Pointer :: torpntr

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  Character(len=strLen) :: xyzfile="music.xyz"
  Integer,Parameter :: MAXCONFS=10000
  Real(Kind=RDbl)            :: minerr=100000_RDbl;
  ! search window width
  Real(kind=RDbl) , Parameter :: INIT_SEARCH_WIDTH=1.00_RDbl 
  Real(kind=RDbl) , Parameter :: FINAL_WIDTH=0.10_RDbl 
  Integer, Parameter :: MAX_OPT_ITERS =5000

  Character(len=2*strLen)    :: longline
  Character(len=strLen),Dimension(10)    :: fields
  Integer :: sarray(2), barray(3), tarray(4)
  Real(kind=RDbl) , Dimension(MAXCONFS) :: nrg
  Real(kind=RDbl) , Dimension(50,50,2) :: stretchlist
  Real(kind=RDbl) , Dimension(50,50,3) :: bendinglist
  Real(kind=RDbl) , Dimension(50,50,4) :: torsionlist
  Integer , Dimension(50) :: stretchlistnpairs
  Integer , Dimension(50) :: bendinglistntriplets
  Integer , Dimension(50) :: torsionlistnquartets
  Real(kind=RDbl) , Dimension(50) :: stretch
  Real(kind=RDbl) , Dimension(50) :: bending
  Real(kind=RDbl) , Dimension(50) :: torsion


  Character(len=strLen), Parameter :: stretch_tag = "STRETCH GROUPS"
  Character(len=strLen), Parameter :: bending_tag = "BENDING GROUPS"
  Character(len=strLen), Parameter :: torsion_tag = "TORSION GROUPS"

  Real(Kind=RDbl)            :: intranrg, nrgerr

  !** Tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action
  !** Simulation Cell Variable(s)
  Character(len=strLen),Parameter :: simcell_tag="Simulation Cell Information"
  Type(SimCell_Params)            :: scell
  !** Configuration Variable(s)
  Type(AtMolCoords), Dimension(:), Pointer :: species
  Type(AtMolCoords),Dimension(:), Pointer  :: storedcoords
  Character(len=strLen), Parameter         :: config_tag = & 
      "Configuration Initialization"
  !** Interaction/Forcefield Parameters and Statistics
  Type(Interaction_Model)                  :: imodel

  !** Configuration file
  Type(CONFILE) :: configfile

  !** MD Simulation time
  Real(Kind=RDbl)            :: mtime



  !** set the default display unit
  Integer, Parameter         :: dUnit = 6

  Integer                    :: nsims, simno, iter, error, nconfs, nspc
  Integer                    :: configunitno, restartunitno   
  Character(len=strLen)      :: comment, filename, restartfile

  Character(len=strLen)    :: line
  Integer :: subset(3), spc, a1, a2, a3, a4
  Integer :: nunit, i, j, k, nfields, natoms, unitno, nparams
  Integer :: nstretch, nbending, ntorsion, npairs, ntriplets, nquartets
  Integer :: optiter, maxoptiter, itype, ninter, inter, typeinters(3)
  Real(kind=RDbl) :: width, reducefac, currentval, newval, ran
  Logical :: success=.False., found


  !-------------------------------------------------------------
  ! Do basic initialization, based on command line input
  !-------------------------------------------------------------
  Call commandline_init(ctrl_filename,action)

  If (action=="WriteSampleCtrlfile") Then
    Write(0,'(a)') dashedline
    Call general_sampleCF(0)
    Write(0,'(a)') dashedline
    Call atom_sampleCF(0)
    Write(0,'(a)') dashedline
    Call molecules_sampleCF(0)
    Write(0,'(a)') dashedline
    Call simcell_sampleCF(0)
    Write(0,'(a)') dashedline
    Call config_sampleCF(0)
    Write(0,'(a)') dashedline
    Stop
  Else If (action /= "DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  End If

  Call file_settag(ctrl_filename, d_ctrl_file)

  !----------------------------------------------
  ! Welcome the user
  !----------------------------------------------
  Write(dUnit,'(a)') "Welcome to Music"
  Write(dUnit,'(2a)') "Reading from control file ",Trim(ctrl_filename)

  !----------------------------------------------
  ! Initialize the general parameters
  !----------------------------------------------
  Call general_init(ctrl_filename)

  !----------------------------------------------------------------------------
  ! Initialize atoms, molecules, simulation cell, configuration, forcefield
  !----------------------------------------------------------------------------
  Call atom_init()
  Call atom_display(dUnit)

  nspc = molecules_init()

  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,2,dUnit)
  Call molecules_display(dUnit)

  !** Allocate the species structure
  Allocate(species(nspc), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'species structure')    

  Call config_init(species, scell, ctrl_filename, config_tag)
  Call config_display(species,dUnit)
  Write(*,*) "-Right now this code is written to read one molecule "
  Write(*,*) "-from an xyz file and calculate its intra molecular energy"
  Write(*,*) "-this can be later generalized to read maultiple molecules from"
  Write(*,*) "-single xyz file"

  Call interact_init(imodel,ctrl_filename,scell,species,'MD')
  Call interact_display(imodel,2,dUnit)

  !----------------------------------------------------------
  ! Initialize the move type information
  !----------------------------------------------------------
  Call wtime_init()
  Call wtime_starttime(MAX_ROUTINES)
  Call wtime_init(1)
  Call wtime_init(2)
  Call wtime_init(3)

  Call stopcriteria_init(ctrl_filename)


  !***********************************************************************
  ! give space for string cords
  ! asuumes 1 spc, max 10000 structures
  !***********************************************************************
  subset=(/1,1,0/)
  Allocate(storedcoords(1), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
  spc=1
  Call config_allocfields(storedcoords(1),spc, MAXCONFS)
  Call config_setnmoles(storedcoords(1),MAXCONFS)


  !***********************************************************************
  ! read file and store coord, nrg etc
  !***********************************************************************
  nconfs=0
  Do 
    nconfs=nconfs+1

    If (nconfs>MAXCONFS) Then
      nconfs=nconfs-1
      Write(*,*) "not enough space for reading all molecules"
      Exit
    Endif


    !read position form xyz file
    Call config_xyz2molec(storedcoords,scell,spc,nconfs,xyzfile,success)
    If (.Not.success) Then
      nconfs=nconfs-1
      Exit
    Endif

  End Do
  Write(*,'(a,i3,2a)') "Found ", nconfs, " configuartaion in file :", xyzfile



  !***********************************************************************
  ! read nrg for nrg.dat file
  !***********************************************************************
  Write(*,*) "Reading nrg in hartrees from xyzfile comment lines"
  call file_close(xyzfile)
  nunit=file_open(xyzfile,110)
  natoms=molecules_getnatoms(1)
  Do i=1,nconfs
    Read(nunit,*)  ! natoms
    Read(nunit,'(a)') longline ! comment line
    nfields=split(longline,fields,":")
    nrg(i)=toreal(fields(2),error)
    If (error/=0) Then
      Write(*,*) "problem reading nrgs from xyz file"
      Write(*,*) longline
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    nrg(i)=nrg(i)*627.5 ! kcal/mol
    ! skip all other lines
    Do j=1,natoms
      Read(nunit,*)
    End Do

  End Do
  Close(nunit)

  !********************************************************************
  ! read stretch section from ctrlfile
  !********************************************************************
  unitno=file_open(ctrl_filename)
  Rewind(unitno)
  !read ctrlfile for stretch groups: can be generated from an MD
  found = filesrchstr(unitno,stretch_tag,line,.True.)
  If (.Not.found) Then
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Endif

  Read(unitno,*) nstretch !number of stretch type groups
  typeinters(STRETCH_INDEX)=nstretch

  Do i=1,nstretch
    ! number of stretch  interaction in each group and the 
    ! value of spring contant
    Read(unitno,*) nparams, stretch(i) 


    stretchlistnpairs(i)=nparams
    ! atom list for each of those stretch interactions
    Do j=1, nparams
      Read(unitno,*) a1, a2
      stretchlist(i,j,1)=a1
      stretchlist(i,j,2)=a2
    End Do

  End Do


  !********************************************************************
  ! read bending section from ctrlfile
  !********************************************************************
  unitno=file_open(ctrl_filename)
  Rewind(unitno)
  !read ctrlfile for bending groups: can be generated from an MD
  found = filesrchstr(unitno,bending_tag,line,.True.)
  If (.Not.found) Then
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Endif

  Read(unitno,*) nbending !number of stretch type groups
  typeinters(BENDING_INDEX)=nbending

  Do i=1,nbending
    ! number of stretch  interaction in each group and the 
    ! value of spring contant
    Read(unitno,*) nparams, bending(i) 


    bendinglistntriplets(i)=nparams
    ! atom list for each of those stretch interactions
    Do j=1, nparams
      Read(unitno,*) a1, a2, a3
      bendinglist(i,j,1)=a1
      bendinglist(i,j,2)=a2
      bendinglist(i,j,3)=a3
    End Do

  End Do



  !********************************************************************
  ! read torsion section from ctrlfile
  !********************************************************************
  unitno=file_open(ctrl_filename)
  Rewind(unitno)
  !read ctrlfile for bending groups: can be generated from an MD
  found = filesrchstr(unitno,torsion_tag,line,.True.)
  If (.Not.found) Then
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Endif

  Read(unitno,*) ntorsion !number of stretch type groups
  typeinters(TORSION_INDEX)=ntorsion

  Do i=1,ntorsion
    ! number of stretch  interaction in each group and the 
    ! value of spring contant
    Read(unitno,*) nparams, torsion(i) 


   torsionlistnquartets(i)=nparams
    ! atom list for each of those stretch interactions
    Do j=1, nparams
      Read(unitno,*) a1, a2, a3, a4
      torsionlist(i,j,1)=a1
      torsionlist(i,j,2)=a2
      torsionlist(i,j,3)=a3
      torsionlist(i,j,4)=a4
    End Do

  End Do





  
  !**********************************************************
  ! optimize parameters
  !**********************************************************
  simno = 1
  width=INIT_SEARCH_WIDTH ! search window width
  maxoptiter=MAX_OPT_ITERS
  !after each iteration serach window is reduced by this factor
  reducefac=(FINAL_WIDTH)**(one/maxoptiter)

  Do optiter=1, maxoptiter
    width=width*reducefac
    itype=Int(3*rranf())+1 !stretch or bend or torsion
!    itype=1
    Select Case(itype)

      ! optimize stretch
    Case (STRETCH_INDEX)
      ninter=typeinters(itype) ! total stretch groups
      ran=rranf()
      inter=Int(ninter*ran)+1 !which atom types?
      currentval=stretch(inter) !current value of this spring constant
      newval=currentval*(1+width*(rranf()-0.5)) !suggested value
      stretch(inter)=newval

      ! set all values from stretch array
      Do k=1,ninter
      !get all atom lists of this type
      npairs=stretchlistnpairs(k)
      Do i=1,npairs
        sarray(1)=   stretchlist(k,i,1)
        sarray(2)=   stretchlist(k,i,2)
        bspntr=>imodel%ff(1)%iparams(1)%stretch
        Call bsmodel_setHarK(bspntr,sarray,stretch(k))
      End Do
    End Do

      ! optimize bending parameters
    Case (BENDING_INDEX)
      ninter=typeinters(itype) ! total stretch groups
      ran=rranf()
      inter=Int(ninter*ran)+1 !which atom types?
      currentval=bending(inter) !current value of this spring constant
      newval=currentval*(1+width*(rranf()-0.5)) !suggested value
      bending(inter)=newval

      ! set all values from stretch array
      Do k=1,ninter
      !get all atom lists of this type
      ntriplets=bendinglistntriplets(k)
      Do i=1,ntriplets
        barray(1)=   bendinglist(k,i,1)
        barray(2)=   bendinglist(k,i,2)
        barray(3)=   bendinglist(k,i,3)
        bbpntr=>imodel%ff(1)%iparams(1)%bending
        Call bbmodel_setHarK(bbpntr,barray,bending(k))
      End Do
    End Do




      ! optimize bending parameters
    Case (TORSION_INDEX)
      ninter=typeinters(itype) ! total stretch groups
      ran=rranf()
      inter=Int(ninter*ran)+1 !which atom types?
      currentval=torsion(inter) !current value of this spring constant
      newval=currentval*(1+width*(rranf()-0.5)) !suggested value
      torsion(inter)=newval

      ! set all values from stretch array
      Do k=1,ninter
      !get all atom lists of this type
      nquartets=torsionlistnquartets(k)
      Do i=1,nquartets
        tarray(1)=  torsionlist(k,i,1)
        tarray(2)=  torsionlist(k,i,2)
        tarray(3)=  torsionlist(k,i,3)
        tarray(4)=  torsionlist(k,i,4)
        torpntr=>imodel%ff(1)%iparams(1)%torsion
        Call tormodel_setA1A2(torpntr,tarray,torsion(k))
      End Do
    End Do




    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select

    !***********************************************************************
    ! Calculate nrg of each structure and add up error with new forcefield
    !***********************************************************************
    nrgerr=zero
    Do iter = 1, nconfs

      Call config_copymolec(storedcoords(1), iter, species(1), 1, .False.)

      ! main force field call, in kJ/mol
      Call ffcall_nrg(imodel,scell,species,"INTRA",intranrg)

      intranrg=intranrg/calToj ! convert to kcal

      nrgerr=nrgerr+((nrg(iter)-intranrg)**2)

    End Do

    nrgerr=Sqrt(nrgerr/nconfs)

    If (nrgerr<minerr) Then
      minerr=nrgerr

    Select Case(itype)
    Case (STRETCH_INDEX)      ! accept new ff
      stretch(inter)=newval
    Case (BENDING_INDEX)      ! accept new ff
      bending(inter)=newval
    Case (TORSION_INDEX)      ! accept new ff
      torsion(inter)=newval
    Case default
     Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
     Stop
    End Select
      Write(*,*) "Total:", nrgerr, " Min ", minerr, "kcal/mol"
      Write(*,'(a,3f12.5)') "constants : ", stretch(1:3)
      Write(*,'(a,4f12.5)') "constants : ", bending(1:4)
      Write(*,'(a,1f12.5)') "constants : ", torsion(1)
    Else
      ! reject changes
    Select Case(itype)
    Case (STRETCH_INDEX)    
      stretch(inter)=currentval
    Case (BENDING_INDEX)    
      bending(inter)=currentval
    Case (TORSION_INDEX)    
      torsion(inter)=currentval
    Case default
     Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
     Stop
    End Select
    Endif


  End Do
  !---------------------------------------
  ! Get the elapsed time
  !---------------------------------------
  Call wtime_stoptime(MAX_ROUTINES)
  Call wtime_display("Main Program", MAX_ROUTINES)

End Program music

