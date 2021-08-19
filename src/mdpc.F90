!-------------------------------------------------------------
! This module contains the different data structures and 
!routines for doing the analysis of results of a MD Run,
! ie. The post code stuff.
!-------------------------------------------------------------
Module mdpc

  Use defaults,    Only : zero, one, RDbl, strlen, MAX_SORBS, MAX_DAT_FILES, &
      d_ctrl_file, INITIAL_SIZE, dashedline, zeroTolerance, lstrLen, &
      dashedline2
  Use utils,       Only : genfilename, isfileopen, filesrchstr, stripcmnt, &
      split, toint, allocerrdisplay
  Use vector,      Only : VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag, vector_getcomp
  Use general,     Only : general_init, general_getnoofconfigs, &
      GenSim_Params, genparams, general_getWriteStep
  Use file,        Only : file_getunit, file_getname , file_settag, &
      file_open, file_close
  Use stats,       Only : Statistics, stats_update
  Use atom,        Only : atom_init, atom_display, atom_getmass
  Use mmath,       Only : mmath_LSFit
  Use molecules,   Only : molecules_init, molecules_gettype, &
      molecules_name, molecules_getatype
  Use config,      Only : config_init, config_getnatoms, config_getnmoles, &
      AtMolcoords
  Use datafile,    Only : CONFILE, datafile_initin, datafile_gencontrolfile, &
      datafile_getSimNum, datafile_hastime, datafile_readconfig
  Use simcell,     Only : Simcell_Params, simcell_init, simcell_display
  Use histogram,   Only : Histogram_Params, histogram_init, histogram_update, &
      histogram_display
  Use interact, Only: Interaction_Model
  Use storestats, Only: Species_Stats
!!$  Use md, Only: MDInfo, md_init, md_getdt

  Implicit None
  Save

  Private
  Public :: MDPostParams, MSDParams, MDScratch, RadProfparams,  mdpc_init, &
      mdpc_createScratchFiles, mdpc_calcSelfDiffusivity, mdpc_display, &
      mdpc_sampleCF, mdpc_calcRadProf

  !** Information about the self diffusivity analysis, which includes
  !** the number (nsorbs) and types (molecList) of molecules to calc
  !** Dself for, the list of points (timeList) at which the MSD is to 
  !** be calculated when constructing a plot of MSD vs. time, the 
  !** amount of time to skip between sucessive MSD calculations in the
  !** ensemble average (blockSkip), and any times to skip at the 
  !** beginning (loSkip) and end (hiSkip) of the configuration file.
  Type MSDparams
    Integer         :: nsorbs
    Integer         :: nTimes
    Real(Kind=RDbl) :: blockSkip
    Real(Kind=RDbl) :: loSkip, hiSkip
    Integer, Dimension(:), Pointer         :: molecList   
    Real(Kind=RDbl), Dimension(:), Pointer :: timeList
  End Type MSDparams

  !** Data type with information for this module. Contains info
  !** about the MD simulation, the base filename for output files,
  !** and pointers to each of the analysis types.
  Type MDPostParams
!!$    Type(MDInfo)               :: mdParams
    Character(len=strLen)         :: mdTag
    Type(MSDparams), Pointer      :: dself
    Type (RadProfparams), Pointer :: radial !@@ added 01/28/04
    Character(len=strlen)         :: outname
  End Type MDPostParams

  !** Data type that contains information about the scratch files
  !** created from the configuration files
  Type MDScratch
    Integer, Dimension(:), Pointer :: unitno
    Integer, Dimension(:), Pointer :: nmoles
    Real(Kind=RDbl) :: totTime
    Integer :: fType
  End Type MDScratch

  !@@ RadProfparams added By Simon on 1/28/04
  ! Information about the Radial profile analysis. Includes the number
  ! of sorbates (nsorbs) and types (molecList) of molecules to calculate
  ! the radial profile, and the time values that will be used in the 
  !calculation. Also have some information about the pore (rpore, 
  ! Xc, Yc) and the number of bins in which the radius will be divided.
  Type RadProfparams
    Integer                                :: nsorbs
    Integer, Dimension(:), Pointer         :: molecList
    Real(Kind=RDbl), Dimension(:), Pointer :: timeList
    Real(Kind=RDBl)                        :: rpore, Xc, Yc,dZ    
    Integer                                :: Nbins, ntimes
  End Type RadProfparams

  !** These variables are used for the MD params hack
  Real(kind=RDbl) :: timeStep
  Character(len=strLen) :: ensemble

  !** Define an interface for the mdpc_calcMSD routines
  Interface mdpc_calcMSD
    Module Procedure mdpc_calcMSDFixed
    Module Procedure mdpc_calcMSDVariableMult
  End Interface

  !** Default id tag for the post code information section.
  !** This section should contain one or more id strings for
  !** post processing tasks.
  Character(len=strLen), Parameter :: mdpc_idstring="MD Post Code Information"

  !** This is the list of valid post processing tasks. 
  Character(len=strLen), Parameter :: mdpc_dselfString="MD Post Diffusivities"

  !@@ Post processing tag for calculating the radial profile. Simon 1/28/04
  Character(len=strLen), Parameter :: mdpc_radprofString = "MD Radial Profile"

  !** This is the tag used in the original control file for the MD section
  Character(len=strLen), Parameter    :: default_md_tag = "MD Information"



Contains

  !----------------------------------------------------------------------------
  ! This subroutine initializes the md post code by reading (1) the MD post
  ! code section from the post code control file, and (2) reading in the 
  ! MD section from the original control file stored in the configuration 
  ! file.
  !----------------------------------------------------------------------------
  Subroutine mdpc_init(mdpp,postfile,newCtrlFile,species,outname,optTag)
    Type(MDPostParams), Intent(InOut) :: mdpp
    Character(len=strlen),Intent(in)::postfile, newCtrlFile
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Character(len=strLen), Intent(In) :: outname
    Character(len=strLen), Intent(In), Optional :: optTag

    !** Check for the optional optTag used to search the control file.
    !** If it is not present, then use the default mdpc_idstring. Call
    !** the subroutine to read in information from the control file.
    If (Present(optTag)) Then
      Call mdpc_readPostFile(mdpp,postfile,optTag)
    Else
      Call mdpc_readPostFile(mdpp,postfile,mdpc_idstring)
    End If
    mdpp%mdTag= "MD Post Code Information"

    !** Set the base output filename passed to us.
    mdpp%outname = outname

!!$    !** With the information about the post code initialized, we need
!!$    !** to initialize the information about the MD run(s) stored in 
!!$    !** the configuration file (time step, ensemble, and so on). Call
!!$    !** md_init to read this information in from the reconstructed
!!$    !** control file.
!!$    Call md_init(mdpp%mdParams,species,newCtrlFile,default_md_tag)


    !** Moves can be a bit nasty to compile on some platforms, so we
    !** can cheat and just use this hack to read the time step in
    !** from the regenerated control file.
    Call mdpc_readCtrF(newCtrlFile,default_md_tag)

  End Subroutine mdpc_init


  !----------------------------------------------------------------------------
  ! Reads the control file that was generated from the configuration file
  ! to get the necessary oinformation regarding the md RUN . I dont want to 
  ! use "md_init" routine because it has lot of extra stuff, not required here
  !----------------------------------------------------------------------------
  Subroutine mdpc_readCtrF(ctrlfile,tag)
    Character(*), Intent(In) :: ctrlfile
    Character(*), Intent(In) :: tag

    Character(len=255) :: text
    Character(len=strLen), Dimension(strLen) :: params

    Integer :: unitno, useMD, j, error, i, k, dUnit 
    Integer :: mdTypes

    !** Open the control file if it is not yet open
    unitno = file_open(ctrlfile)

    !** Find the MD section
    useMD = filesrchstr(unitno,tag,text,.True.)
    If (useMD /= 0) Then

      Write(*,'(a)')
      Write(*,'(a)') dashedline
      Write(*,'(a)') "MD Simulation Information"
      Write(*,'(a)') dashedline

      Read(unitno,'(a)') text
      text = stripcmnt(text)
      j = split(text,params)
      !** Number of MD moves listed
      mdTypes = toint(params(1))
      Read(unitno,'(a)') text
      !** Read in all the MD move types
      ! I am sort of cheating , at present just assuming only one movetype
      ! is there lets proceed, When you have Multiple timesteps rewrite this

      !** Reading only the necessary info out of the MD section
      Do i = 1, mdTypes
        !Read the blank line separating the move types
        Read(unitno,*)
        Read(unitno,'(a)') text
        Read(unitno,'(a)') text
        Read(unitno,*) timeStep
        Write(*,'(a,f10.4)') "Time step used : " ,timeStep
        !        !* time step
        !        text = stripcmnt(text)
        !        j = split(text,params)
        !        postmd%Tstep = toreal(params(1))
        Read(unitno,'(a)') text
        Read(unitno,'(a)') text
        Read(unitno,'(a)') text
        Read(unitno,*) ensemble
        Write(*,'(2a)')      "Ensemble used  : ", ensemble

        Read(unitno,'(a)') text
        Read(unitno,'(a)') text
        Read(unitno,'(a)') text
        Read(unitno,'(a)') text
!!$        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$        Write(*,*) text
      End Do
    End If
  End Subroutine mdpc_readCtrF

  !----------------------------------------------------------------------------
  ! Reads the actions that need to be performed by the post code and 
  ! initializes these actions (such as calculating self diffusivities).
  ! Valid actions must have an idstring defined in this module.
  !----------------------------------------------------------------------------
  Subroutine mdpc_readPostFile(mdpp,postfile,opt_tag)
    Type(MDPostParams), Intent(InOut) :: mdpp
    Character(len=strLen), Intent(in)          :: postfile
    Character(len=strLen), Optional, Intent(In):: opt_tag

    Character(len=strLen)   :: tag, line
    Integer :: unitno, nactions, error
    Integer                 :: lineno,i

    !** Nullify the pointers to each of the possible post processing tasks
    Nullify(mdpp%dself)
    Nullify(mdpp%radial)

    !** Open the postfile if it is not opened
    unitno = file_open(postfile)

    !** Check for the presence of the optional data tag to search for
    tag = mdpc_idstring
    If (Present(opt_tag)) tag = opt_tag

    !** Find the MD post code tag in the file. If it is not there,
    !** report the error.
    lineno = filesrchstr(unitno, tag, line, .True.)
    If (lineno == 0) Then
      Write(0,'(1x,a,i4,3a)') __FILE__,__LINE__,": ", &
          " Post control file does not contain the tag ",trim(tag)
      Stop
    End If

    !** Read the number of post processing actions to take
    Read(unitno,*) nactions

    !** Loop through and read in the actions. Note that the
    !** Actions should correspond to the tags for each action.
    Do i = 1, nactions

      !** Read in the line and strip the comment from it.
      Read(unitno,'(a)') line
      line = stripcmnt(line)

      !** Find the analysis this corresponds to
      Select Case (Trim(line))

      Case (mdpc_dselfString)
        Allocate(mdpp%dself,stat=error)
        If (error /= 0) Then
          Write(0,'(a,i4,2a)') __FILE__,__LINE__,": ", &
              "Error allocating pointer mdpp%dself"
          Stop
        End If


        !@@ added 01/28/04
      Case (mdpc_radprofString)
        Allocate(mdpp%radial,stat=error)
        If (error /= 0) Then
          Write(0,'(a,i4,2a)') __FILE__,__LINE__,": ", &
              "Error allocating pointer mdpp%radial"
          Stop
        End If

      Case Default
        Write(0,'(a,i4,3a)') __FILE__,__LINE__,": ", &
            "Could not identify post processor action string ",trim(line)
        Stop
      End Select

    End Do ! Do i = 1, nactions

    !** Now initialize each of the tasks by checking the allocation
    !** status of the pointers and calling the proper subroutines.

    If (Associated(mdpp%dself)) Then
      !** Initialize the information needed to calculate the self
      !** diffusivities.
      Call mdpc_initDself(mdpp%dself,unitno)
    End If

    !@@ added 01/28/04
    If (Associated(mdpp%radial)) Then
      !** Initialize the information needed to calculate the radial
      !** profile
      Call mdpc_initRadial(mdpp%radial,unitno)
    End If

  End Subroutine mdpc_readPostFile

  !----------------------------------------------------------------------------
  ! Initializes the parameters for performing self diffusivity calculations.
  ! This includes the number and names of sorbates to calculate the self 
  ! diffusivity for, the number of configurations to skip at the beginning
  ! and end of the configuration file, the time points at which to calculate 
  ! the MSDs, and the amount of time to skip between calculations.
  !----------------------------------------------------------------------------
  Subroutine mdpc_initDself(dparams,unitno,optTag)
    Type(MSDparams), Intent(InOut) :: dparams
    Integer, Intent(In)   :: unitno
    Character(*), Intent(In), Optional :: optTag

    Character(len=lstrLen):: line
    Character(len=strLen) :: tag
    Character(len=strLen), Dimension(strLen) :: fields
    Real(kind=RDbl) :: t0, tt, dt
    Integer :: i, nfields, lineno, error

    !** Search for the location of the diffusivity information in the
    !** control file
    tag = mdpc_dselfString
    If (Present(optTag)) tag = optTag

    !** Find the line. The first instance is the one listed in the 
    !** general MD post code section of the code. The second should 
    !** be the section containing the parameters.
    lineno = filesrchstr(unitno, tag, line, .True.)
    lineno = filesrchstr(unitno, tag, line)

    If (lineno == 0) Then
      Write(0,'(a,i4,3a)') __FILE__,__LINE__,": ", &
          "Control file does not contain tag ",trim(tag)
      Write(0,'(2a)') "Make sure the tag is located AFTER the general MD ", &
          "post code information section."
      Stop
    End If

    !** Read in the number of species to do calculations on.
    Read(unitno,*) dparams%nsorbs

    !** Allocate space for the list of molecules we will be calculating
    !** self diffusivities for.
    Allocate(dparams%molecList(dparams%nsorbs),stat=error)
    If (error /= 0) Then
      Write(0,'(2(a,i4))') __FILE__, __LINE__, &
          ": Could not allocate dparams%molecList of size ",dparams%nsorbs
      Stop
    End If

    !** Read in each of the molecule names we wish to calculate the 
    !** self diffusivity for and store their molecule type in the 
    !** molecList.
    Do i = 1, dparams%nsorbs
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields)
      dparams%molecList(i) = molecules_gettype(fields(1))
    End Do

    !** Read in the number of points we would like to plot on the 
    !** MSD vs time plot for calculating the self diffusivity.
    Read(unitno,*) dparams%nTimes

    !** Warn if the number is very small, stop if it is one.
    If (dparams%nTimes == 1) Then
      Write(0,'(2(a,i4))') __FILE__, __LINE__, &
          ": Must have at least 2 points for a linear fit. nTimes = ", &
          dparams%nTimes
      Stop
    End If
    If (dparams%nTimes < 5) Then
      Write(0,'(2(a,i4))') __FILE__, __LINE__, &
          ": WARNING: Requested number of points for linear fit is small! &
          & nTimes = ", dparams%nTimes
    End If

    !** Allocate space for the time list
    Allocate(dparams%timeList(dparams%ntimes),stat=error)
    If (error /= 0) Then
      Write(0,'(2(a,i4))') __FILE__, __LINE__, &
          ": Could not allocate dparams%timeList of size ",dparams%nTimes
      Stop
    End If

    !** Read in the first and last points at which to calculate the 
    !** MSD for the MSD vs time plot.
    Read(unitno,*) t0,tt

    !** Read in the amount to time to skip at the beginning and end of 
    !** the configuration files.
    Read(unitno,*) dparams%loSkip, dparams%hiSkip

    !** Create the list of times at which to calculate the MSD.
    dt=(tt-t0)/Real(dparams%nTimes-1,Kind=RDbl)
    Do i = 1, dparams%nTimes
      dparams%timeList(i) = t0 + Real(i-1,kind=RDbl)*dt
    End Do

    !** Read in the amount of time to skip between ensemble average
    !** calculations
    Read(unitno,*) dparams%blockSkip

  End Subroutine mdpc_initDself


  !----------------------------------------------------------------------------
  ! Initializes the parameters for performing radial profile calculations.
  ! This includes the number and names of sorbates to calculate the radial 
  ! profile for, tha data of the pore and the time points at which to 
  ! do the calculations (mainly initializes a RadProfParams data type
  ! structure).
  !@@ Added by Simon on 1/30/04 
  !----------------------------------------------------------------------------
  Subroutine mdpc_initRadial(RPparams,unitno)
    Type(RadProfparams), Intent(InOut) :: RPparams
    Integer, Intent(In)   :: unitno
    Character(len=lstrLen):: line
    Character(len=strLen) :: tag
    Character(len=strLen), Dimension(strLen) :: fields
    Real(kind=RDbl) :: t0, tt, dt
    Integer :: i, nfields, lineno, error, ntimes


    !** Search for the location of the radial profile information in the
    !** control file
    tag = mdpc_radprofString

    !** Find the line. The first instance is the one listed in the 
    !** general MD post code section of the code. The second should 
    !** be the section containing the parameters.
    lineno = filesrchstr(unitno, tag, line, .True.)
    lineno = filesrchstr(unitno, tag, line)

    If (lineno == 0) Then
      Write(0,'(a,i4,3a)') __FILE__,__LINE__,": ", &
          "Control file does not contain tag ",trim(tag)
      Write(0,'(2a)') "Make sure the tag is located AFTER the general MD ", &
          "post code information section."
      Stop
    End If

    !** Read in the number of species to do calculations on.
    Read(unitno,*) RPparams%nsorbs

    !** Allocate space for the list of molecules we will be calculating
    !** self diffusivities for.
    Allocate(RPparams%molecList(RPparams%nsorbs),stat=error)
    If (error /= 0) Then
      Write(0,'(2(a,i4))') __FILE__, __LINE__, &
          ": Could not allocate RPparams%molecList of size ",RPparams%nsorbs
      Stop
    End If

    !** Read in each of the molecule names we wish to calculate the 
    !** radial profile for and store their molecule type in the 
    !** molecList.
    Do i = 1, RPparams%nsorbs
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields)
      RPparams%molecList(i) = molecules_gettype(fields(1))
    End Do

    ! Read in the pore radius
    Read(unitno,*) RPparams%rpore

    ! Read in the X,Y of the center of the pore
    Read(unitno,*) RPparams%Xc,RPparams%Yc


    !** Read in the first and last time values at which to calculate the 
    !** radial profile.
    Read(unitno,*) t0,tt

    ! Read in the Dt to use
    Read(unitno,*) dt

    !** Read in the number of bins to divide the radius in
    Read(unitno,*) RPparams%Nbins

    !** Read in the dZ to divide the pore in Z 
    Read(unitno,*) RPparams%dZ

    ntimes = Aint((tt -t0)/dt) + 1

    ! store the size of the time vector 
    RPparams%ntimes = ntimes

    !** Allocate space for the time list
    Allocate(RPparams%timeList(RPparams%ntimes),stat=error)
    If (error /= 0) Then
      Write(0,'(2(a,i4))') __FILE__, __LINE__, &
          ": Could not allocate RPparams%timeList of size ",RPparams%nTimes
      Stop
    End If

    !** Create the list of times at which to calculate the Radial Profile
    Do i = 1, nTimes
      RPparams%timeList(i) = t0 + Real(i-1,kind=RDbl)*dt
    End Do

    Write(*,'(a)')'salio de mdpc_initradial'

  End Subroutine mdpc_initRadial


  !----------------------------------------------------------------------------
  ! Creates the COM scratch file by reading configurations from the 
  ! specified config file and writing the COM coordinates and (optionally)
  ! the time and velocities to the scratch file. It returns the unit
  ! numbers of the scratch files stored by molecule type.
  !
  ! fileOpt is a 3 digit binary (IJK) that has the following values:
  !
  ! I = 1, write positions     I = 0, no positions written
  ! J = 1, write velocities    J = 0, no velocities written
  ! K = 1, write time          K = 0, no time written
  !
  ! Its default value is 100, positions only.
  !
  ! @@ Added by Simon 1/30/04
  ! radialOpt is used to determine if the scratch file is being created to 
  ! caclulate radial profile, its values are:
  ! 0 false (default)
  ! 1 true 
  !----------------------------------------------------------------------------
  Subroutine mdpc_createScratchFiles(mdpp, mdpScratch, species, spcnrgs, &
      configFile, fileOpt,radialOpt)
    Type(MDPostParams), Intent(In) :: mdpp
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(CONFILE), Intent(InOut)  :: configFile
    Type(MDScratch), Intent(InOut):: mdpScratch 
    Type(Species_Stats), Intent(InOut) :: spcnrgs
    Integer, Intent(In), Optional :: fileOpt
    Integer, Intent(In), Optional :: radialOpt  

    !** Temperory energy storage variable for reading and passing species
    !** system energies
    !LC    Type (Species_System_Energies), Dimension(:), Pointer :: spcsysnrgs

!!$    Real(Kind=RDbl), Dimension(Size(species,1)) :: intra, ke
!!$    Real(Kind=RDbl), Dimension(Size(species,1),Size(species,1),2) :: pot
    Integer, Dimension(Size(species,1)) :: nmolesList, checkMoles
    Integer :: nconfig, molecType, recLen, fType, sn, nspec,nmoles, error
    Integer :: i, j, n, tempunitno
    Real(Kind=RDbl) :: time, totenergy, dt
    Character(len=strLen) :: fname
    Integer :: nsorbs
    Integer :: radial

    !** Zero the nmoles lists
    checkMoles = 0
    nmolesList = 0

    !** Check on the file options
    fType = 100
    If (Present(fileOpt)) fType = fileOpt

    radial = 0
    If (Present(radialOpt)) radial = radialOpt

    !** Initialize the config file and request the content tag to be assigned
    !SDEBUG    Call datafile_initin(configFile,1)
    Call datafile_initin(configFile, configFile%name)

    !** Get the number of species
    nspec = Size(species, 1)
    !LC    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"spcsysnrgs")
    !LC    Do i=1,nspec
    !LC      !** Initialize the energy reading variable for each sorbate
    !LC      Call interact_initspcsys(spcsysnrgs(i))
    !LC    End Do

    !** Allocate the arrays of scratch
    Allocate(mdpScratch%unitno(nspec),stat=error)
    If (error /= 0) Then
      Write(0,'(a,i5,2a,i4)') __FILE__,__LINE__,": ", &
          "Could not allocate scratch%unitno of size ",nspec
      Stop
    End If

    Allocate(mdpScratch%nmoles(nspec),stat=error)
    If (error /= 0) Then
      Write(0,'(a,i5,2a,i4)') __FILE__,__LINE__,": ", &
          "Could not allocate scratch%nmoles of size ",nspec
      Stop
    End If

    !** Zero the number of molecules and the unit numbers
    mdpScratch%nmoles = 0
    mdpScratch%unitno = 0

    !** Open the scratch files for writing

    !@@ Decide between radial and self. Added on 1/30/04 by Simon    
    If (radial==0) Then 
      nsorbs = mdpp%dself%nsorbs
    Else 
      nsorbs = mdpp%radial%nsorbs
    End If

    Do n = 1, nsorbs

      !** Get the species type
      If (radial==0) Then 
        i = mdpp%dself%molecList(n)
      Else 
        i = mdpp%radial%molecList(n)
      End If

      !** Get the number of molecules of the species
      nmoles = config_getnmoles(species(i))

      !** Get the record length for the file we want to write
#ifndef NAGCOMPILER
      Write(*,*) "this feature (IOLENGTH) has not been checked with", &
          "non nag compilers"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
#endif
#ifdef NAGCOMPILER
      Select Case (fType)

      Case (1)
        Inquire(IOLENGTH=recLen) time
      Case (10)
        Inquire(IOLENGTH=recLen) species(i)%coords(1,1:nmoles)%v
      Case (11)
        Inquire(IOLENGTH=recLen) species(i)%coords(1,1:nmoles)%v, time
      Case (100)
        Inquire(IOLENGTH=recLen) species(i)%coords(1,1:nmoles)%rp
      Case (101)
        Inquire(IOLENGTH=recLen) species(i)%coords(1,1:nmoles)%rp, time
      Case (111)
        Inquire(IOLENGTH=recLen) species(i)%coords(1,1:nmoles)%rp, &
            species(i)%coords(1,1:nmoles)%v, time
      End Select
#endif
      !** We need a unique scratch file number. This will allow us to 
      !** have 99 scratch files for each config file, i.e., 99 molec types
      sn = datafile_getSimNum(configFile)*100+i
      fname = genfilename(configFile%name,sn)

      !** Open the scratch file and get the unit number
      mdpScratch%unitno(i) = file_open(recLen,fname)

!!$      !MDEBUG
!!$      Write(0,'(a,i5,a,2i5,a)') __FILE__,__LINE__," Scratch file number ", &
!!$          mdpScratch%unitno(i),sn,Trim(fname)

    End Do

    !** We need to know how many records there are in the config file
    nconfig = general_getnoofconfigs()

!!$    !** This is the time step used in the simulation, but the time between
!!$    !** writes to the configuration file is different. Get the steps
!!$    !** between writes to find the total time.
!!$    dt = md_getdt(mdpp%mdParams) * general_getWriteStep()

    !** Here we use the hack to prevent compiling of the moves module.
    !** Get the size of the time inbetween 2 successive records in the 
    !** configuration file. It is the simulation time step * number of
    !** steps between writes to the configuration file.
    dt = timeStep * general_getWriteStep()

    !** Loop through the records and write
    Do i = 1, nconfig

      !** Read a single configuration from the config file
      If (datafile_hastime(configFile)) Then
        Call datafile_readconfig(configfile, species, &
            nmolesList,spcnrgs,totenergy,time)
      Else
        !** Calculate the time
        time = dt * Real(i,Kind=RDbl) 
        Call datafile_readconfig(configfile, species, nmolesList, &
            spcnrgs, totenergy)
      End If

      !** Check to make sure the number of molecules is the same as the
      !** first configuration read in. If not, we should stop since we
      !** do MD in the NVE or NVT ensemble. Ok, this is not a thorough
      !** check since we compare the sums, but it's a start.
      If (i == 1) Then
        checkMoles = nmolesList
      Else
        If (Sum(checkMoles) /= Sum(nmolesList)) Then
          Write(0,'(a,i4,2a)') __FILE__,__LINE__,": ", &
              "Error, configurations from config file have different number &
              & of molecules"
          Stop
        End If
      End If

      !** Loop over all the species types
      Do j = 1, nsorbs

        !** Get the molecule type 
        !@@ Decision between radial and self added on 1/30/04 Simon.
        If (radial==0) Then 
          molecType = mdpp%dself%molecList(j)
        Else 
          molecType = mdpp%radial%molecList(j)
        End If

        !** Call writeScratch to create the COM coordinates and
        !** write them to the scratch file
        Call mdpc_WriteScratchRecord(species(molecType), &
            mdpScratch%unitno(molecType),molecType,i,time,fType)
      End Do
    End Do

    !** Record the total time for the simulation
    mdpScratch%totTime = time

    !** Record the number of molecules of each type for the simulation
    mdpScratch%nmoles = nmolesList

    !** Record the file type options
    mdpScratch%fType = fType

    !** Close the config file
    tempunitno=file_close(configFile%name)
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
  End Subroutine mdpc_createScratchFiles


  !----------------------------------------------------------------------------
  ! Writes the center of mass of all molecules to the scratch file
  !so that they can be accessed directly
  !----------------------------------------------------------------------------
  Subroutine mdpc_WriteScratchRecord(species,scratchUnit,mType,recn,time,fType)
    Type(AtMolCoords), Intent(In) :: species
    Integer,Intent(in) :: scratchUnit, mType
    Integer,Intent(in) :: recn
    Real(Kind=RDbl), Optional, Intent(In) :: time
    Integer, Intent(In), Optional :: fType

    Integer :: err,natoms,nmoles
    !** assuming that nobody wants to do 4-Dimensional simulation
    !** the dimension of COMarray is assigned to be 3
    Real(kind=RDbl),Dimension(Size(species%coords,2),3):: COMrp
    Real(kind=RDbl),Dimension(Size(species%coords,2),3):: COMv

    !** Get the number of atoms, molecules, and the time (if required)
    natoms = config_getnatoms(species)
    nmoles = config_getnmoles(species)

    !** Calculate the center of mass coordinates for the species
    !** We only pass Comrp(:,1:3) since there are only 3 dimensions
    Call mdpc_com(species,mType,COMrp(:,1:3),COMv(:,1:3))

    !** Write the requested data to the scratch file
    Select Case (fType)

    Case (1)
      Write(unit=scratchUnit,IOSTAT=err,REC=recn) time
    Case (10)
      Write(unit=scratchUnit,IOSTAT=err,REC=recn) COMv(1:nmoles,1:3)
    Case (11)
      Write(unit=scratchUnit,IOSTAT=err,REC=recn) COMv(1:nmoles,1:3), time
    Case (100)
      Write(unit=scratchUnit,IOSTAT=err,REC=recn) COMrp(1:nmoles,1:3)
    Case (101)
      Write(unit=scratchUnit,IOSTAT=err,REC=recn) COMrp(1:nmoles,1:3),time
    Case (111)
      Write(unit=scratchUnit,IOSTAT=err,REC=recn) COMrp(1:nmoles,1:3), &
          COMv(1:nmoles,1:3),time
    End Select

    !** Report any errors
    If (err/=0) Then
      Write(0,'(2a,i4,a,i3,a,i8,a,i3)') __FILE__," : ",__LINE__, &
          " Error in opening the direct access file unit ",scratchUnit, &
          " for record ",recn,", error number ",err
      Stop
    End If

  End Subroutine mdpc_WriteScratchRecord

  !----------------------------------------------------------------------------
  ! This routine removes the scratch files by closing their unit numbers.
  !----------------------------------------------------------------------------
  Subroutine mdpc_removeScratch(mdpScratch)
    Type(MDScratch), Dimension(:), Intent(InOut) :: mdpScratch
    Integer :: i, j

    Do i = 1, Size(mdpScratch,1)
      Do j = 1, Size(mdpScratch(i)%unitno,1)
        !** Close the file
        Close(mdpScratch(i)%unitno(j))
      End Do
      !** Zero the data
      mdpScratch(i)%nmoles = 0
      mdpScratch(i)%unitno = 0
      mdpScratch(i)%totTime= 0
      mdpScratch(i)%fType  = 0
    End Do
  End Subroutine mdpc_removeScratch

  !----------------------------------------------------------------------------
  ! Writes the center of mass of all molecules to a file so you may
  ! manually review the data
  !----------------------------------------------------------------------------
  Subroutine  mdpc_WriteCOMFile(species,COM,tempunit,molecType,time)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Integer,Intent(in):: tempunit
    Integer,Intent(in):: molecType
    Real(Kind=RDbl), Dimension(:,:), Intent(In) :: COM
    Real(Kind=RDbl), Optional, Intent(In) :: time
    Integer :: natoms,nmoles
    !** assuming that nobody wants to do 4-Dimensional simulation
    !** the dimension of COMarray is assigned to be 3
    Real(kind=RDbl),Dimension(Size(species(molecType)%coords,2),3):: COMarray
    Real(kind=RDbl),Dimension(3):: dist
    Real(Kind=RDbl) :: dist2

    !** Get the number of atoms, molecules, and the time (if required)
    natoms = config_getnatoms(species,molecType)
    nmoles = config_getnmoles(species,molecType)

    !** Calculate the center of mass coordinates for the species
    !** We only pass COMarray(:,1:3) since there are only 3 dimensions
    Call mdpc_com(species(molecType),molecType,COMarray(:,1:3))

    dist(1) = COMarray(1,1)-COM(1,1)
    dist(2) = COMarray(1,2)-COM(1,2)
    dist(3) = COMarray(1,3)-COM(1,3)
    dist2 = dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3)

    !** Write the requested data to the scratch file
    If (Present(time)) Then
      Write(tempunit,'(3f20.3)') time,Sqrt(dist2), &
          mag(species(molecType)%coords(1,1)%v)
    Else
      Write(tempunit,'(3f20.3)') dist2
    End If

  End Subroutine mdpc_WriteCOMFile

  !----------------------------------------------------------------------------
  ! Calculates the center of mass, of all molecules of type nsorb,
  ! and returns in COM array. If the optional real 2D array COMv
  ! is passed, the COM velocities are calculated as well.
  !----------------------------------------------------------------------------
  Subroutine mdpc_com(species,molecType,COM,COMv)  
    Type(AtMolCoords), Intent(In) :: species
    Integer, Intent(In) :: molecType
    Real(kind=RDbl),Dimension(:,:),Intent(out):: COM
    Real(kind=RDbl),Dimension(:,:), Intent(out), Optional :: COMv
    Logical :: vels

    Real(kind=RDbl):: tmass,atmass
    Integer:: i,j,nmoles,natoms

    !** Get the number of molecules of type molecType
    nmoles=config_getnmoles(species)

    !** Zero the center of mass values
    COM = 0.0_RDbl
    vels = .False.
    If (Present(COMv)) Then
      COMv = 0.0_RDbl
      vels = .True.
    End If

    !**** the second index 1,2,3 represents x,y, and z respectively
    Do i=1,nmoles
      natoms=config_getnatoms(species)
      tmass=0.0000
      If (natoms == 1) Then
        COM(i,:) = species%coords(1,i)%rp
        If (vels) COMv(i,:) = species%coords(1,i)%v
      Else
        Do j=1,natoms

          !** Get the atomic mass for averaging
          atmass=atom_getmass(molecules_getatype(molecType,j))
          tmass=tmass+atmass

          !** Center of mass positions
          COM(i,1)=COM(i,1)+atmass*species%coords(j,i)%rp%comp(1)
          COM(i,2)=COM(i,2)+atmass*species%coords(j,i)%rp%comp(2)
          COM(i,3)=COM(i,3)+atmass*species%coords(j,i)%rp%comp(3)

          !** Center of mass velocities
          If (vels) Then
            COMv(i,1)=COMv(i,1)+atmass*species%coords(j,i)%v%comp(1)
            COMv(i,2)=COMv(i,2)+atmass*species%coords(j,i)%v%comp(2)
            COMv(i,3)=COMv(i,3)+atmass*species%coords(j,i)%v%comp(3)
          End If

        End Do
        COM(i,1)=COM(i,1)/tmass
        COM(i,2)=COM(i,2)/tmass
        COM(i,3)=COM(i,3)/tmass

        If (vels) Then
          COMv(i,1) = COMv(i,1)/tmass
          COMv(i,2) = COMv(i,2)/tmass
          COMv(i,3) = COMv(i,3)/tmass
        End If

      End If
    End Do
  End Subroutine mdpc_com

  !----------------------------------------------------------------------------
  ! Calculates the self diffusivities given the scratch file unit numbers. The
  ! scratch files should contain at least positions, but may also contain
  ! velocities and time information as well. Once again, this is specified
  ! with the fileOpt integer (see mdpc_createScratchFiles for its
  ! definition). 
  !
  ! The self diffusivities are returned in the 2D array dself, whose 
  ! dimensions are dself(Sim. Num, Species Type). 
  !
  ! There are two additional options. In order to improve statistics, the
  ! mean squared displacements may be averaged over several configuration
  ! files. Pass avgConfig=.True. for this. If you wish to write the mean
  ! squared displacements to a file for plotting later, pass writeMSD=.True.
  !----------------------------------------------------------------------------
  Subroutine mdpc_calcSelfDiffusivity(mdpp,mdpScratch,dself,avgConf,writeMSD)
    Type(MDPostParams), Intent(In):: mdpp
    Type(MDScratch), Dimension(:), Intent(In)   :: mdpScratch 
    Real(Kind=RDbl), Dimension(:,:), Intent(Out) :: dself
    Logical, Intent(In), Optional :: avgConf
    Logical, Intent(In), Optional :: writeMSD

    Logical :: avgMSD, dumpMSD
    Integer :: molecType, i, j
    Real(Kind=RDbl), Dimension(mdpp%dself%nTimes,4) :: MSD
    Real(Kind=RDbl), Dimension(mdpp%dself%nTimes,4) :: fit
    Real(Kind=RDbl), Dimension(4)   :: slope, intercept, dev
    Character(len=strLen) :: filename

    !** Initialize the self diffusivity values to zero
    dself = 0.0_RDbl

    !** In order to calculate the diffusivities, we must:
    !** (1) compute the mean squared displacements (MSDs)
    !** (2) fit a line to the MSD data
    !** (3) calculate the slope

    !** Check to see if we are averaging over many config files
    avgMSD = .False.
    If (Present(avgConf)) avgMSD = avgConf

    !** Check to see if we are dumping the MSD data to a file
    dumpMSD = .False.
    If (Present(writeMSD)) dumpMSD = writeMSD

    !** Do each sorbate type separately
    Do i = 1, mdpp%dself%nsorbs

      !** Get the molecule type
      molecType = mdpp%dself%molecList(i)

      !** First, get the MSDs. Pass the avgConfig flag if it
      !** is present.
      If (avgMSD) Then
        Call mdpc_CalcMSDVariableMult(mdpp%dself,mdpScratch,molecType,MSD)
      Else
!!$        Do j = 1, Size(mdpScratch,1)
!!$          Call mdpc_CalcMSDVariable(mdpp%dself,(/mdpScratch(j)/),molecType,MSD)
!!$        End Do
      End If

      !** Now that we have the MSDs at a number of points, calculate
      !** the slope. Calculating the slope of a MSD vs time plot yields
      !** the self diffusivity.
      Do j = 1, 4
        Call mmath_LSFit(mdpp%dself%timeList, MSD(:,j), &
            mdpp%dself%nTimes, slope(j), intercept(j), dev(j))

        !** Store the self diffusivity
        dself(i,j) = intercept(j)

        If (dumpMSD) Then 
          !** Check to see if we are dumping the data to a file. If so,
          !** fit a line to the data, then write the MSDs and the fitted
          !** line to the file
          Call mdpc_GetLine(mdpp%dself%timeList,mdpp%dself%nTimes,slope(j), &
              intercept(j),fit(:,j))
        End If
      End Do

      !** Now dump the file if requested
      If (dumpMSD) Then
        If (avgMSD) Then
          filename = Trim(mdpp%outname)//"."//Trim(molecules_name(molecType)) &
              //".MSD"
        End If
        Call mdpc_writeMSD(mdpp,filename,MSD,fit)
      End If

    End Do

    !** Display the results of the calculation.
    Call mdpc_displaySlopes(slope,dev,intercept,6)

  End Subroutine mdpc_calcSelfDiffusivity

  !----------------------------------------------------------------------
  !Calculates the MSD, using the data given in the scratch file,
  !It takes the unit number of the scratch file(unitno), number of the 
  !molecule(nsorb).  
  !----------------------------------------------------------------------
  Subroutine mdpc_CalcMSDFixed(unitno,nsorb,MSD)  
    Integer,Intent(in) :: unitno
    Integer,Intent(in) :: nsorb
    Real(kind=RDbl),Dimension(:,:),Pointer::MSD
!!$    
!!$    Real(kind=RDbl),Dimension(:,:),Pointer::COM_0,COM_t
!!$    Real(kind=RDbl) ::Tot_t,tskip,time_diff,t,t0
!!$    Real(Kind=RDbl) :: timeInterval, totalTime
!!$    Real(kind=RDbl),Dimension(4) ::dist
!!$    Real(kind=RDbl),Dimension(INITIAL_SIZE,4) ::distlist
!!$
!!$    Integer:: i,j,k,maxj,nconf,tot_configs,nskip,nct,nc0,err
!!$    Type(Statistics) :: MSDstats
!!$    
!!$    !** Allocate the size of the arrays for storing MSD values
!!$    Call mdpc_MSDinit(MSD,COM_0,COM_t)
!!$
!!$    !** Loop among the different times at which MSD should be calculated    
!!$    !** Each time is the ensemble average mean squared displacement for
!!$    !** that given time interval. We accomplish this by averaging the 
!!$    !** MSD for a given time interval over a number of configurations
!!$    !** from the configuration file.
!!$    Do i=1,postmd%MSD%nTimes
!!$      
!!$      !** Zero the MSD values
!!$      MSD(i,:)=0
!!$      
!!$      !** Get the time interval for which we are calculating an MSD
!!$      timeInterval = postmd%MSD%TimeList(i)
!!$      
!!$      !** Get the number of configurations to skip for calculating the
!!$      !** ensemble average MSD
!!$      tskip = postmd%MSD%dat_skip
!!$      
!!$      !** Get the total number of stored configurations
!!$      tot_configs = general_getnoofconfigs()
!!$      
!!$      !** Get the total simulation time
!!$      totalTime = mdpc_GetTotTime()
!!$
!!$      !** Calculate the number of configurations that add up to 
!!$      !** the given time interval
!!$      nconf = Nint(tot_configs*(timeInterval/totalTime))
!!$      
!!$      !** Calculate the number of configurations that add up to 
!!$      !** the desired skip time
!!$      nskip=Nint(tot_configs*(tskip/totalTime))
!!$
!!$      !** Report the information to the user
!!$      Write(*,'(a)') " "
!!$      Write(*,'(a)') "*************************************************"
!!$      Write(*,'(a,f12.4,a)') "Calculating Mean Square Displacement  For Time &
!!$          &= ", timeInterval," ps"
!!$      Write(*,'(a,i6,a)') "Uses ",nconf, " configurations during each interval"
!!$      Write(*,'(a,i6,a)') "Skips ",nskip, " configurations between &
!!$          & each interval"
!!$
!!$      !** Warn if there are not enough configurations to calculate the average
!!$      If (nconf==0) Then
!!$        nconf=1
!!$        Write(*,'(a,f12.4,2a)') "Time specified too small (",timeInterval, &
!!$            " ps) ; Less than the difference between two consecutive & 
!!$            &configurations in the scratch file : ",Trim(file_getname(unitno))
!!$      End If
!!$
!!$      !** Warn is there is no value for skipping between MSD points
!!$      !** or if the calculation requires too many configurations
!!$      If (nskip<1) Then
!!$        Write(0,'(1x,a,2(a,i4))') __FILE__," : ",__LINE__, &
!!$            " nskip error in mdpc, nskip = "
!!$        Stop
!!$      Else If (nconf>tot_configs) Then
!!$        Write(0,'(1x,a,3(a,i4))') __FILE__," : ",__LINE__, &
!!$            " nconf > tot_configs in mdpc, nconf = ",nconf, &
!!$            "; tot_configs = ",tot_configs
!!$        Stop
!!$      End If
!!$
!!$      !** Calculate the actual number of values to be averages to get 
!!$      !** the ensemble average
!!$      maxj=Int(Real(tot_configs-nconf)/nskip)
!!$      If (maxj<1) Then
!!$        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$        Stop
!!$      End If
!!$
!!$      !**** Calculate MSD's at various regions of the data, for the 
!!$      !**** given time difference
!!$
!!$      Write(*,'(a,i6,a)') "Will do the calculation in ",maxj, " different &
!!$          & intervals"
!!$
!!$      !** Loop through each individual calculation and add it to the sum
!!$      Do j=1,maxj
!!$
!!$        !** Starting configuration for this value
!!$        nc0=1+(j-1)*nskip
!!$        Read(UNIT=unitno,IOSTAT=err,REC=nc0) COM_0
!!$        If (err/=0) Write(0,'(1x,2a,i4)') __FILE__," Read error at ",__LINE__
!!$
!!$        !** Ending configuration for this value
!!$        nct=1+nconf+(j-1)*nskip
!!$        Read(UNIT=unitno,IOSTAT=err,REC=nct) COM_t
!!$        If (err/=0) Write(0,'(1x,2a,i4)') __FILE__," Read error at ",__LINE__
!!$
!!$        !** Calculate the actual MSD for these two points
!!$        Call mdpc_distsum(COM_0,COM_t,dist,distlist)
!!$
!!$        Do k=1,postmd%nmoles
!!$          Call histogram_update(MSD_hist(i),distlist(k,4))
!!$       End Do
!!$
!!$        !** Add the values to the MSD sum
!!$        Do k=1,4
!!$          MSD(i,k)=MSD(i,k)+dist(k)/maxj/postmd%nmoles
!!$        End Do
!!$
!!$      End Do
!!$
!!$    End Do
!!$    
  End Subroutine mdpc_CalcMSDFixed

  !----------------------------------------------------------------------------
  ! Calculates the MSD using the data given in the scratch file.
  ! It takes the unit number of the scratch file(unitno) and the
  ! molecule type. It deals with varying time intervals between 
  ! configurations in the file by interpolating between stored 
  ! points. This subroutine averages over several configuration files
  ! in order to improve statistics. The scratch file information (unit nos,
  ! total time, number of molecules) are passed in with scratch.
  !----------------------------------------------------------------------------
  Subroutine mdpc_CalcMSDVariableMult(dsParams,scratch,mtype,MSD)  
    Type(MSDparams), Intent(In)              :: dsParams 
    Type(MDScratch), Intent(In), Dimension(:):: scratch
    Integer, Intent(In)                      :: mtype
    Real(kind=RDbl), Dimension(:,:)          :: MSD

    Real(kind=RDbl), Dimension(:,:), Pointer :: COM_0,COM_t
    Real(kind=RDbl) :: tskip,dt,startT
    Real(kind=RDbl), Dimension(4) ::dist
    Real(kind=RDbl), Dimension(:,:), Allocatable ::distlist
    Logical, Dimension(Size(scratch,1))      :: skipme

    Integer:: i,j,k,maxj,nconf,nskip,err

    Real(kind=RDbl) :: t, t0, rnmoles, interval
    Integer :: totConfigs, n, filen

    !** Allocates the arrays for storing MSD values. WARNING: This will
    !** break if you attempt to use a different number of molecules in
    !** each of the config files
    Allocate(COM_t(scratch(1)%nmoles(mtype),3), &
        COM_0(scratch(1)%nmoles(mtype),3), &
        distlist(scratch(1)%nmoles(mtype),4),STAT=err)
    If (err/= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, & 
          " Could not allocate 'MSD'"
      Stop
    End If

    !*** Loop among the different configs at which MSD should be calculated    
    !** Get the amount of time to skip between calculations
    tskip = dsParams%blockSkip

    !** Get the total number of configurations in the file. It must be the
    !** same for each configuration file for the time being. This could
    !** be coded so it would take different numbers of tot configurations.
    totConfigs = general_getnoofconfigs()


    !** We go from longest to shortest time in order to weed out 
    !** data files that aren't long enough
    skipMe = .False.
    Do i = dsParams%nTimes,1,-1

      !** Get the time interval we are interested in.
      interval = dsParams%TimeList(i)

      !** Set the counter to zero. This is the number of averages
      !** we have
      n = 0

      !** Zero the MSD values for this time point
      MSD(i,:) = 0.0_RDbl

      !** Loop through the different scratch files and average over them
      !** as well
      Do filen = 1, Size(scratch,1)

        !** If we weeded out this file, skip it
        If (skipMe(filen)) Cycle

        !** Calculate the average time between configs in the simulation
        dt = scratch(filen)%totTime/Real(totConfigs,Kind=RDbl)

        !** Calculate the number of configurations that add up to 
        !** the desired skip time
        nskip=Nint(totConfigs*(tskip/scratch(filen)%totTime))

        !** Calculate the approximate number of configurations
        !** that we need for the given time interval
        nconf=Nint(interval/scratch(filen)%totTime*Real(totConfigs,Kind=RDbl))

        !** Update the user about our diabolical plans
        Write(*,'(a)') " "
        Write(*,'(a)') dashedline
        Write(*,'(a,f18.2,a)') "Calculating Mean Square Displacement for time &
            & = ",interval," ps"
        Write(*,'(a,i12,a)') "Using ",nconf, &
            " configurations during each interval"
        Write(*,'(a,i6,a)') "Skiping ",nskip, " configurations between &
            & each interval"

        !** Check to make sure we actually have enough configurations to 
        !** calculate a MSD
        If (nconf==0) Then
          nconf=1
          Write(*,'(a,f12.4,2a)') "Time specified too small (",interval," ps);&
              & Less than the difference between two consecutive & 
              &configurations in the scratch file : "
        End If

        !** Check to make sure the we aren't using too many configurations
        !** and that we are actually skipping some configurations
        If (nskip<1) Then
          Write(0,'(a,i4,2a,i6)') __FILE__,__LINE__,":", &
              " Error: Specified time to skip between configurations is &
              & too small."
          Write(0,'(a,f8.1,a,5x,a,i10)') "Time : ",tskip," ps","Configs : ", &
              nskip
          Stop
        Else If (nconf > totConfigs) Then
          Write(0,'(a,i5,2a,i6)') __FILE__,__LINE__,":", &
              " ERROR, number of configurations required for the given time &
              &interval is too large."
          Write(0,'(a,f8.1,a,2(5x,a,i10))') "Time : ",interval," ps", &
              "Req. Configs : ",nconf,"Total Configs : ",totConfigs
          Write(0,'(2(a,f8.1,5x))') "dt : ",dt,"Total Time : ", &
              scratch(filen)%totTime
          !NEW STUFF FOR BACKWARDS STEP
          skipMe(filen) = .True.
          Write(0,'(a,i5,a,i5)') __FILE__,__LINE__, &
              " I'm skipping this scratch file from now on : ",filen
          cycle
        End If

        !** Figure out how many values we will average together to get
        !** the ensemble average MSD for the given time
        !** This is not the actual number we will use! With a flexible
        !** time scale, only some will fall in the appropriate range.
        maxj=Int(Real(totConfigs-nconf)/nskip)
        If (maxj<1) Then
          Write(0,'(1x,a,i5,2a,i6)') __FILE__,__LINE__," : ", &
              "Number of values to average for MSD is too small. maxj = ",maxj
          Write(0,'(1x,3(a,i6))') "totConfigs=",totConfigs," nconf=",nconf, &
              " nskip=",nskip
          Stop
        End If

        Write(*,'(a,i6,a)') "Will do the calculation in ",maxj, " different &
            & intervals"

        !** Number of molecules of the given species, as a real. We
        !** can get this from the scratch information. Remember, this
        !** number had better be the same for each config file unless
        !** you have run your systems at infinite dillution!
        rnmoles = Real(scratch(filen)%nmoles(mtype),Kind=RDbl)

        !** Set the start time
        Read(Unit=scratch(filen)%unitno(mtype),IOSTAT=err,REC=1) COM_0, startT

        !** Loop through the different configurations
        Do j=1,maxj

          !** Calculate the correct starting time for this interval
          !** Add in any time to skip at the beginning of the configuration
          !** file (loSkip) at this point.
          t0 = Real(j-1,kind=RDbl)*tskip+startT+dsParams%loSkip

          !** Calculate the correct ending time for this interval
          !          t = interval + Real(j-1,kind=RDbl)*tskip+startT
          ! From what I understand this is correct above is wrong 
          !      -shaji.02/28/2003
          t = interval + Real(j-1,kind=RDbl)*tskip+startT+dsParams%loSkip

          !** Check to make sure we haven't gone out of range.
          !** Add in any time to skip at the end of the configuration
          !** file (hiSkip) here.
          If (t > (scratch(filen)%totTime - dsParams%hiSkip)) Exit

          !** Search the records and find the best approximation
          !** by linear interpolation between stored points
          Call mdpc_findCOM(scratch(filen)%unitno(mtype),t0,dt, &
              totConfigs,COM_0)

          !** Search the records and find the best approximation
          !** by linear interpolation between stored points
          Call mdpc_findCOM(scratch(filen)%unitno(mtype),t,dt,totConfigs,COM_t)

          !** Increment the counter
          n = n + 1

          !** Calculate the MSD for the given points
          Call mdpc_distsum(COM_0,COM_t,dist,distlist)

          !** Update the four values (x, y, z, and total MSD)
          Do k=1,4
            MSD(i,k)=MSD(i,k)+dist(k)/rnmoles
          End Do

        End Do

      End Do

      !** Report the number used to calculate the ensemble average
      Write(0,'(2a,i5,a,i4)') __FILE__,":",__LINE__, &
          " Number of values used to calculate ensemble average : ",n

      If (n==0) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) " No values used, so setting MSD to zero"
        MSD(i,:) = zero
      Else
        !** Divide by the number of accepted MSD values to get the
        !** ensemble average MSD
        MSD(i,:) = MSD(i,:)/Real(n,kind=RDbl)
      Endif
    End Do

    !** Report any files that were skipped
    Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
        " Files skipped because there weren't enough time points : "
    Do filen = 1, Size(skipMe,1)
      If (skipMe(filen)) Write(0,'(a,i5)') "File number skipped : ",filen
    End Do

  End Subroutine mdpc_CalcMSDVariableMult


  !----------------------------------------------------------------------------
  ! Calculates the radial profile given the scratch file unit number.
  ! The scratch file should contain positions, but may also contain
  ! velocities and time information. This is specified with the fileOpt 
  ! integer (see mdpc_createScratchFiles for its definition). 
  !
  ! The radial profile (fraction of particles on each bin) and the "axial 
  ! flux profile" (number of particles crossing that bin flowing) are 
  ! returned in arrays named profile and flux
  ! profile and flux have dimensions of Nbins x (Species+1)  
  !
  ! So far it uses only one configuration file 
  ! Simon 1/30/04
  !----------------------------------------------------------------------------
  Subroutine mdpc_calcRadProf(mdpp,mdpScratch)
    Type(MDPostParams), Intent(In):: mdpp
    Type(MDScratch), Intent(In)   :: mdpScratch 

    Real(Kind=RDbl), Dimension(:,:),Pointer  :: profile,flux
    Integer                          :: molecType, i, j, k
    Character(len=strLen)            :: filename

    Integer          :: totConfigs, nmoles
    Real(kind=RDbl), Dimension(:,:), Pointer :: COM_0,COM_1

    Integer          :: ig, zmove, Zplace0, Zplace1, err,q,nc
    Integer          :: bi, bf, b
    Real (Kind=RDbl) :: rpore2inv, dist2, thickinv, dtConf, xp, yp, zp


    ! Calculate this quantities that will be used in the bin placing
    rpore2inv = 1.0/(mdpp%radial%rpore*mdpp%radial%rpore)
    thickinv = 1.0/mdpp%radial%dZ

    ! Allocates the array to store the results
    Allocate(profile(mdpp%radial%Nbins,mdpp%radial%nsorbs +1), STAT=err)
    If (err/= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, & 
          " Could not allocate profile"
      Stop
    End If

    ! Allocates the array to store the results
    Allocate(flux(mdpp%radial%Nbins, mdpp%radial%nsorbs +1 ),STAT=err)
    If (err/= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, & 
          " Could not allocate profile"
      Stop
    End If

    !** Initialize the profile array to zero
    profile = 0.0
    flux = 0.0


    ! First we need to determine the number of config and the dt 
    ! used to dump configurations some information about the conf file   
    ! Number of configurations:
    totConfigs = general_getnoofconfigs()
    ! determine dtConf used to dump the configurations 
    dtConf = mdpScratch%totTime/Real(totConfigs,Kind=RDbl)

    ! In order to calculate the radial profile, we need to calculate
    ! the distance from each particle to the center of the pore.
    ! and then with that distance calculate the bin where it's located
    ! At the end we need to normalize with the total number of 
    ! configurations used in the calculations.  

    ! Do each sorbate type separately
    Do i = 1, mdpp%radial%nsorbs

      ! Get the molecule type
      molecType = mdpp%radial%molecList(i)

      ! Allocates the arrays for storing COM values
      Allocate(COM_0(mdpScratch%nmoles(molectype),3), &
          COM_1(mdpScratch%nmoles(molectype),3),STAT=err)
      If (err/= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, & 
            " Could not allocate COM "
        Stop
      End If


      ! Loop through the different times
      Do j = 1, mdpp%radial%ntimes

        ! Get the COM for the given time
        ! So far the calculations only uses 1 configuration file, 
        ! thats why mdpSratch(1)%..

        Call mdpc_findCOM(mdpScratch%unitno(molecType),&
            mdpp%radial%timeList(j),dtConf,totConfigs,COM_0)

        ! Get the size of the COM array (should be the number of molecules)
        nmoles = Size(COM_0,1)

        ! Loop through the particles
        Do k = 1, nmoles
          ! Calculate the square of the distance between the particle and the 
          ! center of the pore
          dist2 = (COM_0(k,1) - mdpp%radial%Xc)**2 + (COM_0(k,2) - &
              mdpp%radial%Yc)**2

          ! Place in the bin
          ig = Aint(dist2 * rpore2inv * mdpp%radial%Nbins + 1)

          ! If is "inside" the wall, place in the last bin
          If (ig > mdpp%radial%Nbins ) Then 
            ig = mdpp%radial%Nbins
          End If

          ! Add one to the counter of that bin for all molecules
          ! and one to the species counter
          profile(ig,i) = profile(ig,i)+ 1.0         
          profile(ig,mdpp%radial%nsorbs + 1) = &
              profile(ig,mdpp%radial%nsorbs+1)+ 1.0

          ! Now go on to determine the flux in each bin
          If (j > 1) Then
            ! COM_0 is the positions at this step, and COM_1 is the position 
            ! at the previous step 
            Zplace0 = Aint(COM_0(k,3)*thickinv)
            Zplace1 = Aint(COM_1(k,3)*thickinv)

            ! The difference will tell if the particle had a net movement and
            ! how many planes in Z it crossed
            zmove = abs(Zplace0 - Zplace1)

            ! If it moved, we register the radial bins where it crossed
            If (zmove /= 0 ) Then

              ! We need to determine where it crossed and add 1 to the flux 
              ! profile of that bin

              Do nc = 1, zmove

                ! determine if the particle is moving foward or backwards,
                ! this is used to determine the crossing plane in case of 
                !multiple crossings
                If (Zplace1>Zplace0) Then
                  ! The particle is moving backwards
                  zp= mdpp%radial%dZ*(Zplace1 - nc + 1)
                Else
                  ! the particle is moving foward
                  zp= mdpp%radial%dZ*(Zplace1 + nc)
                End If

                ! Now determine the crossing coordinates  
                xp= COM_1(k,1) + ( COM_0(k,1)-COM_1(k,1) )*&
                    ( zp-COM_1(k,3))/(COM_0(k,3)-COM_1(k,3))

                yp=  COM_1(k,2) + ( COM_0(k,2)-COM_1(k,2) )*&
                    ( zp-COM_1(k,3))/(COM_0(k,3)-COM_1(k,3))

                ! Place the crossing point on a radial bin
                dist2 = (xp - mdpp%radial%Xc)**2 + (yp - &
                    mdpp%radial%Yc)**2
                bi = Aint(dist2 * rpore2inv * mdpp%radial%Nbins + 1)

                ! Add the crossing to the flux profile
                flux(bi,mdpp%radial%nsorbs + 1) = flux(bi, &
                    mdpp%radial%nsorbs + 1) + 1.0 
                flux(bi,i) = flux(bi,i) + 1.0

              End Do

              !calculate the bin in the previous step
              ! dist2 = (COM_1(k,1) - mdpp%radial%Xc)**2 + (COM_1(k,2) - &
              !     mdpp%radial%Yc)**2
              ! bi = Aint(dist2 * rpore2inv * mdpp%radial%Nbins + 1)

              ! We know the present bin: ig              
              ! Now, to do the "DO" from the smaller value:
              ! If (bi<ig) Then
              !   bf = ig
              ! Else
              !   bf = bi
              !   bi = ig
              ! End If
              ! 
              ! Do b = bi,bf
              ! Add a fraction 1/(crossed bins) to each bin crossed 
              ! in the counter of all molecules and in each species counter 

              !   flux(b, (mdpp%radial%nsorbs + 1) ) = flux(b, &
              !      ( mdpp%radial%nsorbs + 1) ) + 1.0/(bf-bi+1.0)
              !   flux(b,i) = flux(b,i) + 1.0/(bf-bi+1.0)

              ! End Do

            End If

          End If

        End Do

        ! Store the present positions temporary in COM_1 (used in flux)
        COM_1 = COM_0

      End Do

      ! Normalize the results with the number of configurations used
      profile = profile/Real(mdpp%radial%ntimes)
      flux = flux/Real(mdpp%radial%ntimes)

    End Do

    ! write the results
    filename = Trim(mdpp%outname)//".RADIAL"
    Write(*,'(a)')'creating file : ',filename
    Call mdpc_writeRadial(mdpp,profile,flux,filename)

  End Subroutine mdpc_calcRadProf


  !----------------------------------------------------------------------------
  ! Gets the COM coordinates from the scratch file for a given time.
  ! If no exact match for the time is found, it will provide a linear 
  ! estimate of the value.
  !----------------------------------------------------------------------------
  Subroutine mdpc_findCOM(unitno,time,dt,lastRec,COM)
    Integer, Intent(In) :: unitno, lastRec
    Real(Kind=RDbl), Intent(In) :: time, dt
    Real(kind=RDbl), Dimension(:,:), Intent(Out) :: COM
    Real(Kind=RDbl), Dimension(Size(COM,1),Size(COM,2)) :: COMhi, COMlo
    Real(kind=RDbl) :: t, thi, tlo
    Integer :: guessRec, err, i, hiRec, loRec

    !** Zero
    COMhi = 0.0_RDbl
    COMlo = 0.0_RDbl

    !** Make a guess as to the record 
    guessRec = Nint(time/dt)
    If (guessRec == 0) guessRec = 1
    !** Check that record
    Read(UNIT=unitno,IOSTAT=err,REC=guessRec) COM, t
    If (err /= 0) Then
      Write(0,'(a,i6,a,i3)') __FILE__,__LINE__, &
          " : File error ",err
      Stop
    End If

    !** See if we should return now
    If (Abs(time-t) < zeroTolerance) Return

    !** Now see if we need to interpolate
    If (t < time) Then
      !** Loop through the records until we find the next greatest time
      Do i = guessRec+1, lastRec
        Read(Unit=unitno,IOSTAT=err,REC=i) COMhi, thi
        If (err /= 0) Then
          Write(0,'(a,i6,a,i3)') __FILE__,__LINE__, &
              " : File error ",err
          Stop
        End If
        If (thi < time) Cycle
        hiRec = i
        If (thi > time) Exit
      End Do
      loRec = hiRec-1
      Read(Unit=unitno,IOSTAT=err,REC=loRec) COMlo, tlo
      If (err /= 0) Then
        Write(0,'(a,i6,a,i3)') __FILE__,__LINE__, &
            " : File error ",err
        Stop
      End If
    Else If (t > time) Then
      !** Loop through the records until we find the next smallest time
      Do i = guessRec-1, 1, -1
        Read(Unit=unitno,IOSTAT=err,REC=i) COMlo, tlo
        If (err /= 0) Then
          Write(0,'(a,i6,a,i3)') __FILE__,__LINE__, &
              " : File error ",err
          Stop
        End If
        If (tlo > time) Cycle
        loRec = i
        If (tlo < time) Exit
      End Do
      hiRec = loRec+1
      Read(Unit=unitno,IOSTAT=err,REC=hiRec) COMhi, thi
      If (err /= 0) Then
        Write(0,'(a,i6,a,i3)') __FILE__,__LINE__, &
            " : File error ",err
        Stop
      End If
    End If

    !** Special case of loRec == 1, check to see if we are
    !** within the real tolerance
    If (Abs(tlo-time) < zeroTolerance) Then
      COM = COMlo
      Return
    Else If (Abs(thi - time) < zeroTolerance) Then
      COM = COMhi
      Return
    End If

    !** Calculate the COM using the high and low values using
    !** linear interpolation between points
    Do i = 1, Size(COM,1)
      COM(i,1) = (COMhi(i,1)-COMlo(i,1))/(thi-tlo)*(time-tlo)+COMlo(i,1)
      COM(i,2) = (COMhi(i,2)-COMlo(i,2))/(thi-tlo)*(time-tlo)+COMlo(i,2)
      COM(i,3) = (COMhi(i,3)-COMlo(i,3))/(thi-tlo)*(time-tlo)+COMlo(i,3)
    End Do

  End Subroutine mdpc_findCOM

  !----------------------------------------------------------------------------
  ! Calculates the displacements given two center of mass arrays COM1 and COM2.
  ! It returns the scalar displacements of each dimension (x,y,z in dimensions
  ! 1, 2, and 3 respectively), plus the total scalar displacement (in
  ! dimension 4) in tdist_sq. It returns the vectors in dist_list.
  !----------------------------------------------------------------------------
  Subroutine mdpc_distsum(COM1,COM2,tdist_sq,dist_list)
    Real(kind=RDbl), Dimension(:,:), Intent(In)  :: COM1
    Real(kind=RDbl), Dimension(:,:), Intent(In)  :: COM2
    Real(kind=RDbl), Dimension(:,:), Intent(out) :: dist_list 
    Real(kind=RDbl), Dimension(4), Intent(out)   :: tdist_sq

    Type(VecType):: vec1,vec2,dist    
    Integer:: nmoles,i,j,DIM

    !** Get the size of the COM array (should be the number of molecules)
    nmoles = Size(COM1,1)

    !** 3 dimensions
    DIM=3

    !** Zero the values
    tdist_sq = 0.0_RDbl
    dist_list= 0.0_RDbl

    !** Loop through the COM array and calculate the displacement.
    Do i=1,nmoles
      !**** using the vector module operator overloading stuff.
      vec1 = COM1(i,1:DIM)
      vec2 = COM2(i,1:DIM)
      dist = vec1 - vec2
      !first 3 componenets are x,y,z distances squared, the 4th component will
      ! be the distance between two molecules
      Do j=1,DIM
        dist_list(i,j)=(vector_getcomp(dist,j))**2
        tdist_sq(j) = tdist_sq(j) + dist_list(i,j)
      End Do
      dist_list(i,4)=dist_list(i,1)+dist_list(i,2)+dist_list(i,3)
    End Do

    tdist_sq(4)=tdist_sq(1)+tdist_sq(2)+tdist_sq(3)
  End Subroutine mdpc_distsum

  !---------------------------------------------------------
  !gives the number of the molecule " name" in the sorbates
  ! array 
  !---------------------------------------------------------
  Integer Function mdpc_GetSorb_n(name)
    Character(*), Intent(IN)         :: name
    mdpc_GetSorb_n=molecules_gettype(name)
  End Function mdpc_GetSorb_n

  !------------------------------------------
  ! Gives the y coordinates,given x values and the array length N, with 
  ! slope=m , y-intercept =c ; ie line  y=mx+c
  !------------------------------------------
  Subroutine mdpc_GetLine(x,N,m,c,y)
    Real(kind=RDbl), Dimension(:),Intent(in)::x
    Real(kind=RDbl), Dimension(:),Intent(out)::y
    Real(kind=RDbl),Intent(in)::m,c
    Integer,Intent(in):: N

    Integer::i
    If ( .Not.((Size(x)<N).Or.(Size(y)<N)) ) Then
      Do i=1,N
        y(i)=m*x(i)+c
      End Do
    Else
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "array sizes not correct"
    Endif
  End Subroutine  mdpc_GetLine

  !----------------------------------------------------------------------------
  ! Writes the MSD calculated , and the least square fitted MSD values 
  ! to the given file denoted by the unitno. Please note that conversion to
  ! m^2/s assumes that the internal units are Ang^2/ps!
  !----------------------------------------------------------------------------
  Subroutine mdpc_writeMSD(mdpp,filename,MSD,fit)
    Type(MDPostParams), Intent(In) :: mdpp
    Character(*), Intent(In) :: filename
    Real(kind=RDbl), Dimension(:,:), Intent(In) :: MSD, fit

    Integer :: j, unitno

    !** Open the file
    unitno = file_open(filename)

    !** Write some header information to the file so we can figure out
    !** what it is.
    Write(unitno,*) "% MSD from an MD Run, for one sorbate"
    Write(unitno,*) "% Column 1: Time in ps"
    Write(unitno,*) "% Column 2..4: MSD in x,y,z directions, Ang^2"
    Write(unitno,*) "% Column 5: MSD in 3D (sum of columns 2..4), Ang^2"
    Write(unitno,*) "% Column 6..8: (Fitted) MSD in x,y,z directions, Ang^2"
    Write(unitno,*) &
        "% Column 9: (Fitted) MSD in 3D (sum of columns 6..8), Ang^2"
    Write(unitno,*)"% MSD from an MD Run, for one sorbate"

    !** Dump the data
    Do j=1,mdpp%dself%nTimes
      Write(unitno,'(9e14.3)') mdpp%dself%timeList(j),MSD(j,1:4),fit(j,1:4)
    End Do

    !** Close it up
    Close(unitno)

  End Subroutine  mdpc_writeMSD

  !------------------------------------------
  ! Writes the MSD calculated , and the least square fitted MSD values 
  !to the given file denoted by the unitno
  !------------------------------------------
  Subroutine mdpc_displaySlopes(slope,dev,intc,optunit)
    Real (kind=RDbl),Dimension(4),Intent(in):: slope
    Real (kind=RDbl),Dimension(4),Intent(in):: dev
    Real (kind=RDbl),Dimension(4),Intent(in):: intc
    Integer,Intent(in),Optional::optunit

    Integer::unit
    If (Present(optunit)) Then
      unit=optunit
    Else
      unit=6
    Endif
    Write(unit,*)
    Write(unit,'(a)') dashedline
    Write(unit,'(a)') "The diffusivites Calculated by mean square diplacement"
    Write(unit,'(a,t32,3(a,12x),a)') "COMPONENTS             |","X","Y","Z", &
        "Total"
    Write(unit,'(a,4e13.4)') "Diffusivity (Ang^2/ps) |", slope(1)/2, &
        slope(2)/2,slope(3)/2, (slope(4))/6
    Write(unit,'(a,4e13.4)') "Diffusivity (m^2/s)    |", slope(1)/2*1.0E-8, &
        slope(2)/2*1.0E-8,slope(3)/2*1.0E-8, (slope(4))/6*1.0E-8
    Write(unit,'(a,4(f9.4,4x))') "% deviation in LS Fit  |",&
        dev(1),dev(2),dev(3),dev(4)
    Write(unit,'(a,4(f9.4,4x))') "Intercepts             |",&
        intc(1),intc(2),intc(3),intc(4)
    Write(unit,'(a)') dashedline
  End Subroutine  mdpc_displayslopes


  !------------------------------------------------------------------
  ! Writes the radial profile and the radial flux profile 
  ! simon 30/01/04
  !-----------------------------------------------------------------
  Subroutine mdpc_writeRadial(mdpp,profile,flux,filename)
    Type(MDPostParams), Intent(In):: mdpp
    Real(Kind=RDbl), Dimension(:,:), Intent (In)  :: profile,flux   
    Character(len=strLen), Intent (In)            :: filename 
    Integer :: j, unitno, i
    Character(len=strLen) :: frmt

    !** Open the file
    unitno = file_open(filename)

    ! Write some header information to the file so we can figure out
    ! what it is.
    Write(unitno,*) "% Radial analysis from an MD Run"
    Write(unitno,*) "% Column  1: # of bin"

    Do i=1,mdpp%radial%nsorbs
      Write(unitno,'(a,i2,a,a)') " % Column ",i+1,": ",&
          molecules_name( mdpp%radial%molecList(i))
    End Do

    Write (unitno,'(a,i2,a)')  " % Column ", i+1, ": Total particles" 

    Write(unitno,'(a)') dashedline  
    Write(unitno,*) "Radial Profile (average particles per bin)"
    Write(unitno,'(a)') dashedline
    ! Dump the data
    Write(frmt,'(a,i4,a)') "(a,",mdpp%radial%nsorbs+1,"(f8.4,1x))"
    Write(unitno,frmt) " ",(profile(j,:),j=1,mdpp%radial%Nbins) 

    Write(unitno,'(a)') dashedline  
    Write(unitno,*) "Radial Flux Profile (crossings per bin)"
    Write(unitno,'(a)') dashedline  
    Write(frmt,'(a,i4,a)') "(a,",mdpp%radial%nsorbs+1,"(f8.4,1x))"
    Write(unitno,frmt) " ",(flux(j,:),j=1,mdpp%radial%Nbins) 

    !** Close it up
    Close(unitno)

  End Subroutine  mdpc_writeRadial


  !----------------------------------------------------------------------------
  ! Displays parameters for the MD post code uses in its various analyses.
  !----------------------------------------------------------------------------
  Subroutine mdpc_display(mdpp,optUnit)
    Type(MDPostParams), Intent(In) :: mdpp
    Integer, Intent(In), Optional  :: optUnit

    Integer :: dunit, i
    Character(len=strLen) :: frmt

    !** Check for the optional display unit
    dunit = 6
    If (Present(optUnit)) dunit = optUnit

    !** Present a nice header to everyone
    Write(dunit,'(a)') dashedline
    Write(dunit,'(a)') "MD Post Processing Code Information"
    Write(dunit,'(a)') dashedline

    !** Display information about each of the analysis types
    If (Associated(mdpp%dself)) Then
      Write(dunit,'(a)') "Self Diffusivity Calculation"
      Write(dunit,'(a)') dashedline

      Write(dunit,'(a,2x)') "Species to calculate Dself for : ", &
          (Trim(molecules_name(mdpp%dself%molecList(i))),i=1,mdpp%dself%nsorbs)
      Write(frmt,'(a,i4,a)') "(a,",mdpp%dself%nTimes,"(f8.2,1x))"
      Write(dunit,frmt)     "Times to calculate MSDs at     : ", &
          (mdpp%dself%timeList(i),i=1,mdpp%dself%nTimes)
      Write(dunit,'(a,f8.2)') "Time to skip between block avg : ", &
          mdpp%dself%blockSkip
      Write(dunit,'(a,f8.2)') "Time to ignore at beginning    : ", &
          mdpp%dself%loSkip
      Write(dunit,'(a,f8.2)') "Time to ignore at end          : ", &
          mdpp%dself%hiSkip
    End If

    If (Associated(mdpp%radial)) Then
      Write(dunit,'(a)') dashedline
      Write(dunit,'(a)') "Radial Profile calculation"
      Write(dunit,'(a)') dashedline
      Write(dunit,'(a,2x)') "Species to calculate radial profile for : ", &
          (Trim(molecules_name(mdpp%radial%molecList(i))),&
          i=1,mdpp%radial%nsorbs)
      Write(dunit,'(a,f8.2)')     "Pore radius                 : ",&
          mdpp%radial%rpore
      Write(dunit,'(a,f8.2,f8.2)')"X and Y center coordinates  : ", &
          mdpp%radial%Xc,mdpp%radial%Yc
      Write(frmt,'(a,i4,a)') "(a,",mdpp%radial%nTimes,"(f8.2,1x))"
      Write(dunit,frmt)           "Times to calculate profiles : ", &
          (mdpp%radial%timeList(i),i=1,mdpp%radial%nTimes)
      Write(dunit,'(a,i4,a)')     "Bins to divide the radius   : ", &
          mdpp%radial%Nbins
      Write(dunit,'(a,f8.2)')     "dZ to divide the pore       : ", &
          mdpp%radial%dZ
      Write(dunit,'(a)') dashedline
    End If

  End Subroutine mdpc_display
  
  !----------------------------------------------------------------------------
  ! Writes information about what it should read in from a control file
  ! during mdpc_init
  !----------------------------------------------------------------------------
  Subroutine mdpc_sampleCF(unitno)
    Integer , Intent(in) :: unitno
#ifndef NAGCOMPILER
    Write(*,*) "mdpc_smaplecf has not been checke dwiht non-NAG compilers"
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
#endif

    Write(unitno,'(a)') "-- "//Trim(mdpc_idstring)//" ---"
    Write(unitno,'(a,t30,a)') 'Integer','# number of actions'
    Write(unitno,'(a,t30,a)') 'Character','# Tag for diffusivity section'


    Write(unitno,'(a)') "-- "//Trim(mdpc_dselfstring)//" ---"
    Write(unitno,'(a,t30,a)') 'Integer','# number of sorbates'
    Write(unitno,'(a,t30,a)') 'Character','# Name of sorbate-1'
    Write(unitno,'(a,t30,a)') 'Integer','# number of time values for &
        &diffusivity'
    Write(unitno,'(a,t30,a)') 'Real, Real','# higest and lowest time &
        &values, (picoseconds)'
    Write(unitno,'(a,t30,a)') 'Real, Real','# initial and final data to be &
        &skipped in % '
    Write(unitno,'(a,t30,a)') 'Real','# time skipped bewteen calculations &
        &of displ.(ps)'

    !@@ Added by simon 1/30/04
    Write(unitno,'(a)') "-- "//Trim(mdpc_radprofString)//" ---"
    Write(unitno,'(a,t30,a)') 'Integer','# number of sorbates'
    Write(unitno,'(a,t30,a)') 'Character','# Name of sorbates'
    Write(unitno,'(a,t30,a)') 'Real','# pore radius (A)'
    Write(unitno,'(a,t30,a)') 'Real, Real','# Xcenter, Ycenter of &
        & the pore'
    Write(unitno,'(a,t30,a)') 'Real, Real','# Highest and lowest &
        & time values (ps) '
    Write(unitno,'(a,t30,a)') 'Real','# time intervals (ps)'
    Write(unitno,'(a,t30,a)') 'Integer','# bins to divide the radius in' 
    Write(unitno, '(a,t30,a)') 'Real','# dZ to divide the pore in Z & 
        & and determine flux profile'

  End Subroutine mdpc_sampleCF


End Module mdpc





