!---------------------------------------------
! This module contains the overall simulation
! parameters	
!---------------------------------------------
Module general

  Use defaults, Only: RDbl, strLen, lstrLen, dashedline
  Use utils, Only: isfileopen, filesrchstr, stripcmnt
  Use file, Only: file_getunit
  Use random, Only: random_getiseed, random_init

  Implicit None

  Private
  Public :: GenSim_Params, general_getfiledesc, general_getcurrentiter, &
      general_getcurrenttime, general_getContentTag, genparams, &
      general_init, general_getnoofconfigs, general_display, &
      general_sampleCF, general_getWriteStep, general_setCurrentSimno, &
      general_getCurrentSimno

  !** The tag marking the beginning of the General section
  Character(len=strLen), Parameter    :: default_general_tag = &
      "General Information"
  
  Type GenSim_Params
    Character(len=lstrLen)    :: filedesc
    Integer                   :: currentiteration, currentsimno
    Real(kind=Rdbl)           :: currenttime
    Integer                   :: niterations
    Integer                   :: iprint
    Integer                   :: icrash
    Integer                   :: iconfig
    Integer                   :: simstart
    Integer                   :: content_tag = -1 ! default value
    Character(len=strLen)     :: restartfile
    Character(len=strLen)     :: configfile
    Character(len=10*strLen)  :: siminfo_string
    Character(len=strLen)     :: displaymode
    Character(len=strLen)     :: ctrlfile
  End Type GenSim_Params

  Type(GenSim_Params), Save   :: genparams

  Interface general_display
    Module Procedure general_displayUnfrmt
    Module Procedure general_displaySummary
  End Interface

Contains
 
  !------------------------------------------------------
  ! Initialises the general parameters by reading the general section
  ! of control file
  !------------------------------------------------------
  Subroutine general_init(ctrl_filename, opt_tag)
    Character(*), Intent(in)  :: ctrl_filename
    Character(*), Optional, Intent(in)  :: opt_tag
    
    Character(len=strLen)   :: tag, line
    Character(len=strlen)   :: Date,Time,Zone
    Integer,Dimension(8)    :: Values
    Integer                 :: unitno, lineno, iseed

    If (Present(opt_tag)) Then
      tag = opt_tag
    Else
      tag = default_general_tag
    End If

    !** Open the ctrl_file if it is not opened
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
      Open(file=ctrl_filename, unit=unitno)
    Endif
    
    !** Find the General section
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif

    Write(*,*) "Initializing General Parameters"
    !** Read the necessary stuff
    Read(unitno, '(a)') genparams%filedesc
    Read(unitno, *) genparams%niterations
    Read(unitno, *) genparams%iprint
    Read(unitno, *) genparams%icrash
    Read(unitno, *) genparams%iconfig
    Read(unitno, *) genparams%simstart
    Read(unitno, *) iseed
    Read(unitno, *) genparams%content_tag
    Read(unitno,'(a)') line
    genparams%restartfile = stripcmnt(line)
    Read(unitno,'(a)') line
    genparams%configfile = stripcmnt(line)
!    Read(unitno, *) genparams%siminfo_string
    Call Date_and_time(Date,Time,Zone,Values)
    genparams%siminfo_string="This Simulation was done in   Year : &
     &"//date(1:4)//",   Month &
     &: "//date(5:6)//",   on date : "//date(7:8)//",    at Time : &
     &"//time(1:2)//" hrs "//time(3:4)//" minutes. ###  The time specified &
     &(local time) varies from GMT by "//zone(1:3)//" hrs"//zone(4:5)//" &
     &minutes.  The control file used was "//Trim(ctrl_filename)//". "
    genparams%currentiteration = 0
    genparams%currenttime = 0.0_RDbl

    !** Initialize the random no. generator
    Call random_init(iseed)

  End Subroutine general_init

  !----------------------------------------------------------------------------
  ! Writes information about what it should read in from a control file to
  ! the unit unitno
  !----------------------------------------------------------------------------
  Subroutine general_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a30)') default_general_tag
    Write(unitno,'(a,t30,a)') 'Character','# File description'
    Write(unitno,'(a,t30,a)') 'Integer','# Number of iterations'
    Write(unitno,'(a,t30,a)') 'Integer','# Steps between writes to I/O'
    Write(unitno,'(a,t30,a)') 'Integer', &
        '# Steps between writes to restart file'
    Write(unitno,'(a,t30,a)') 'Integer','# Steps between writes to config file'
    Write(unitno,'(a,t30,a)') 'Integer','# Start numbering simulations from...'
    Write(unitno,'(a,t30,a)') 'Integer','# Random seed'
    Write(unitno,'(a,t30,a)') '[1234]','# Config file output type'
    Write(unitno,'(a,t30,a)') 'Character','# Restart file filename'
    Write(unitno,'(a,t30,a)') 'Character','# Configuration file filename'
  End Subroutine general_sampleCF

  !------------------------------------------------------
  ! Return the file description
  !------------------------------------------------------
!!$  Character(len=lstrLen) Function general_getfiledesc()
  Function general_getfiledesc()
    Character(len=lstrLen) :: general_getfiledesc
    general_getfiledesc = genparams%filedesc
  End Function general_getfiledesc

  !------------------------------------------------------
  ! Return the tag that specifies what is in the config file
  !------------------------------------------------------
  Integer Function general_getContentTag()
    general_getContentTag = genparams%content_tag
  End Function general_getContentTag

  !------------------------------------------------------
  ! Set current simno
  !------------------------------------------------------
  Subroutine general_setCurrentSimno(simno)
    Integer , Intent(in) :: simno
    genparams%currentsimno=simno
  End Subroutine  general_setCurrentSimno

  !------------------------------------------------------
  ! Return the current simno
  !------------------------------------------------------
  Integer Function general_getCurrentSimno()
    general_getcurrentsimno = genparams%currentsimno
  End Function general_getCurrentSimno

  !------------------------------------------------------
  ! Return the current no. of iterations
  !------------------------------------------------------
  Integer Function general_getcurrentiter()
    general_getcurrentiter = genparams%currentiteration
  End Function general_getcurrentiter

  !------------------------------------------------------
  ! Return the current simulation time
  !------------------------------------------------------
  Real(kind=RDbl) Function general_getcurrenttime()
    general_getcurrenttime = genparams%currenttime
  End Function general_getcurrenttime

  !------------------------------------------------------
  ! Returns the configuration file name
  !------------------------------------------------------
!!$  Character(len=strLen) Function general_getconfigfilename()
  Function general_getconfigfilename()
    Character(len=strLen) :: general_getconfigfilename
    general_getconfigfilename = genparams%configfile
  End Function general_getconfigfilename

  !------------------------------------------------------
  ! Return the simulation start index
  !------------------------------------------------------
  Integer Function general_getsimstart()
    general_getsimstart = genparams%simstart
  End Function general_getsimstart

  !------------------------------------------------------
  ! Return the number of configurations in the index file
  !------------------------------------------------------
  Integer Function general_getnoofconfigs()
    general_getnoofconfigs=genparams%niterations/genparams%iconfig
  End Function general_getnoofconfigs

  !----------------------------------------------------------------------------
  ! Returns the number of steps between writes to the config file.
  !----------------------------------------------------------------------------
  Integer Function general_getWriteStep()
    general_getWriteStep = genparams%iconfig
  End Function general_getWriteStep

  !--------------------------------------------------------------------
  ! Write out the fields of "genparams" without any formatting to the 
  ! output unit "unitno"
  !--------------------------------------------------------------------
  Subroutine general_displayUnfrmt(unitno)
    Integer, Intent(in)              :: unitno
    
    Write(unitno,*) genparams%filedesc
    Write(unitno,*) genparams%currentiteration
    Write(unitno,*) genparams%currenttime
    Write(unitno,*) genparams%niterations
    Write(unitno,*) genparams%iprint
    Write(unitno,*) genparams%icrash
    Write(unitno,*) genparams%iconfig
    Write(unitno,*) genparams%simstart
    Write(unitno,*) random_getiseed()
    Write(unitno,*) genparams%restartfile
    Write(unitno,*) genparams%configfile
  End Subroutine general_displayUnfrmt

  !--------------------------------------------------------------------
  ! Write out the fields of "genparams" without any formatting to the 
  ! output unit "unitno"
  !--------------------------------------------------------------------
  Subroutine general_displaySummary(unitno,nspc)
    Integer, Intent(in) :: unitno
    Integer, Intent(In) :: nspc
    Integer :: i
    Character(len=nspc) :: spaces

    !** Make the spaces
    spaces = ""
    Do i = 1, nspc
      spaces = spaces//" "
    End Do

    !** Write the formated parameters
    Write(unitno,'(a)') dashedline
    Write(unitno,'(a)') "General Simulation Parameters"
    Write(unitno,'(a)') dashedline
    Write(unitno,'(3a)')     spaces,"Description  : ",trim(genparams%filedesc)
    Write(unitno,'(2a,i10)') spaces,"Interations  : ",genparams%niterations
    Write(unitno,'(2a,i10,a)') spaces,"Screen output: ",genparams%iprint, &
        " steps"
    Write(unitno,'(2a,i10,a)') spaces,"Restart write: ",genparams%icrash, &
        " steps"
    Write(unitno,'(2a,i10,a)') spaces,"Config writes: ",genparams%iconfig, &
        " steps"
    Write(unitno,'(2a,f10.1)') spaces,"Random seed  : ",random_getiseed()
    Write(unitno,'(3a)') spaces,"Restart file : ",trim(genparams%restartfile)
    Write(unitno,'(3a)') spaces,"Config file  : ",trim(genparams%configfile)
    Write(unitno,'(2a,i2)') spaces,"Content Tag  : ",genparams%content_tag
  End Subroutine general_displaySummary

End Module general


