!----------------------------------------------------------------------------
! This module has the necessary datatypes and subroutines for working with 
! the binary (non-ascii) data file which stores all data during simulation 
! 
! The object of type CONFILE should be initialized during a simulation
! for writing Data (positions, energies, time) from a data file.  Similarly
! during post processing it should be initialized for reading the same
! data from the data file. The CONFILE object contains info about what are
! things written into it, how many datapoints are written into it etc.
!
! Things to be stored are 'nmoles', 'total nrg', 'time', 'pair_energy', 
! 'intra nrg', 'posns' and 'velocities' 
! The datatype ContentDetails contains logical flags for each of the 
! above fields to say whether they are to be stored or not
! Note that 'nmoles' is always stored
!
! Important routines : 
! datafile_initout -- initialize a CONFILE object for writing the output 
!                     during a simulation 
! datafile_initin  -- initialize a CONFILE object for reading the input 
!                     from there for use during post-processing
! datafile_readconfig -- read a configuration (a configuration means
!                        whatever was written into the datafile during 
!                        one iteration). It acts as an interface for two 
!                        routines. One for reading for old datafiles, and 
!                        one for new.
! datafile_getExtra_Info -- takes a string and returns some extra info 
!                           corresponding to that string (only if it was 
!                           originally stored)
! datafile_writeconfig -- write a configuration into the data file
! datafile_gencontrolfile -- regenerate the control file that was used 
!                            during the original simulation , by reading the 
!                            header of the datafile
! datafile_sampleCF -- datafile section of the control file.
!
! NOTES :
!  1) at some point the old datafile sections should be removed 
!  2) if any new fields need to be added to the datafile then add them in
!     the Extra_Info, so that reading of older datafiles will not have to be
!     modified everytime. (see datafile_readExtraInfo)
! 
!!$! FOR DEALING WITH OLD DATA FILES !!$!  CONTENT_TAG DEFINITIONS:
!!$!  content_tag=1 only nmoles, and total nrg written !!$!
!!    content_tag=2 only nmoles, nrg and pair nrg written !!$!
!!    content_tag=3 only nmoles, nrg ,pair nrg and coords written !!$!
!!    content_tag=4 nmoles, nrg, pairnrg, coords and velocity are written
!!$!  content_tag=5 nmoles, nrg, pairnrg, coords, velocity and time
! are written
!----------------------------------------------------------------------------

Module datafile

  Use config, Only: AtMolCoords, config_getnatoms, config_getnmoles, &
      config_isfixed
  Use defaults, Only: RDbl, strLen, lstrLen, MAX_CTRL_LINES, MAX_SORBS, &
      d_con_file, d_dat_file, TOTAL_INDEX, zero
  Use file, Only: file_getunit, file_settag, file_checkfileopen,file_open
  Use general, Only: general_getContentTag, GenSim_Params, genparams
  Use molecules, Only: molecules_getnsorbs, molecules_gettype, &
      molecules_name, molecules_getnatoms
  Use storestats, Only: Species_Stats, storestats_getnoncoul, &
      storestats_getcoul, storestats_getke, storestats_getintranrg, &
      storestats_updateEnergySS, storestats_updatenrg
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, &
      tolower, cleanstring, allocErrDisplay
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_zerovec

  Implicit None
  Save

  Private
  Public :: CONFILE, FileHeader_Info, Extra_Info, ContentDetails, &
      datafile_writeconfig, datafile_setExtraInfo, &
      datafile_initout, datafile_initin, datafile_readconfig, &
      datafile_gencontrolfile, datafile_close, datafile_getsimnum, &
      datafile_hastime, datafile_hasposition, datafile_hasvelocity, &
      datafile_getExtraInfo, datafile_sampleCF, datafile_getnconfs, &
      datafile_readOneMolec, datafile_readOneSpc, datafile_setContentTag, &
      datafile_clean

  !** Tells what should be written to a datafile
  Type ContentDetails
    Logical :: tot_nrg, pair_nrg, intra_nrg
    Logical :: posns, vels
    Logical :: time
  End Type ContentDetails

  ! Somethings which dont fit anywhere else. These are simulation specific and 
  ! need not be there for all simulations. More fields can be added here 
  ! without breaking postcode (hopefully)
  Type Extra_Info
    Integer                    :: simno
    Real(kind=Rdbl)            :: simPressure , simTemperature
    Character(len=strLen)      :: simType
    Character(len=10*strLen)   :: simInfoString
    Integer                    :: nspcs
    Character(len=strLen), Dimension(MAX_SORBS) :: spcnamelist
    Logical, Dimension(MAX_SORBS) :: is_fixed 
  End Type Extra_Info

  Type FileHeader_Info
    Character(len=strLen)   :: ctrl_filename
    Integer                 :: nlines 
    Character(len=2*lstrLen), Dimension(MAX_CTRL_LINES) :: ctrl_line
    Integer                 :: simno
    Character(len=10*strLen)   :: siminfo_string
  End Type FileHeader_Info

  Type CONFILE
    Character(len=strlen) :: name
    Integer               :: simno
    Integer               :: unitno
    Type(FileHeader_Info) :: header
    Type(ContentDetails)  :: is_storing

    !* Data during how many iterations will be write into this file
    Integer               :: niters       

    !* and how often will it be written
    Integer               :: iconfig

    !* Some other extra info the main driver program would like to stor
    Type(Extra_Info)      :: xtrainfo

    !-------------------------------------------------------------
    ! OBSOLETE, This is left here for working woth old *.con files
    Integer :: content_tag   
    Logical :: datafile_is_old
    !-------------------------------------------------------------
  End Type CONFILE

  Character(len=strLen), Parameter :: default_datafile_tag = &
      "Main Datafile Information"
  Integer, Parameter               :: MAX_EXTRA_INFO_LINES=100

Contains  

  !--------------------------------------------------------------------------
  ! Initializes the given file configfile for writing to it
  ! Requires:  ctrl_filename -- control file for the sim
  !            configfile -- the object to be intialized
  !            sim_n -- index of this simulation
  !            niter -- data during how many iters will be wriiten here 
  !            iconfig -- and how often
  !            opt_section_tag -- section which contains details of datafile
  !--------------------------------------------------------------------------
  Subroutine datafile_initout(configfile,ctrl_filename, species, &
      sim_n, niters, iconfig, simType, opt_section_tag)
    Character(*), Intent(In)                        :: ctrl_filename
    Type(CONFILE), Intent(InOut)                    :: configfile
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Integer,Intent(In)                              :: sim_n, niters, iconfig
    Character(*), Intent(In)                        :: simType
    Character(*), Intent(In), Optional              :: opt_section_tag 

    Integer                       :: funitno,error, nspcs, i
    Character(len=strLen)         :: section_tag=default_datafile_tag
    Character(len=strLen)         :: tempstring

    If (present(opt_section_tag)) section_tag=opt_section_tag

    !SDEBUG funitno = file_open(configfile%name,101)
    funitno = isfileopen(configfile%name)
    If (funitno < 0) Then
      funitno = file_getunit(configfile%name)
      Open(file=configfile%name,unit=funitno, status='unknown', &
          form='unformatted',IOSTAT=error)

      If (error /= 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            "  Could not open file ", Trim(configfile%name)
        Write(*,*) funitno,configfile%name
        Stop
      End If
    Else
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "datafile already initialised, Something Wrong ??!!" 
      Stop
    End If

    configfile%unitno=funitno
    configfile%simno=sim_n
    configfile%niters=niters 
    configfile%iconfig=iconfig

    Call file_settag(configfile%name, d_con_file)

    ! decide what are things to be stored in this binary file by reading 
    ! the datafile section from control file
    Call datafile_initContents(configfile, ctrl_filename, section_tag)

    ! write the original control file to the header
    Call datafile_writeheader(configfile,ctrl_filename)

    ! Write the storage info also to header so that while reading from 
    ! here we will know what to Read, ONLY FOR new *.con.* files
    Write(configfile%unitno) configfile%is_storing

    configfile%xtrainfo%simno=sim_n
    configfile%xtrainfo%simInfoString=configfile%header%siminfo_string

    nspcs=molecules_getnsorbs()
    configfile%xtrainfo%nspcs = nspcs 

    Do i=1,nspcs
      configfile%xtrainfo%spcnamelist(i)=molecules_name(i)
      configfile%xtrainfo%is_fixed(i)=config_isfixed(species(i))
    End Do

    tempstring=cleanstring(simType)
    Select Case(Trim(tempstring))
    Case("MC", "MD", "HMC", "GCMC", "NVTMC")
      configfile%xtrainfo%simType = tempstring
    Case default
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          'Unknown simtype "', Trim(simType),'" passed here'
      Stop
    End Select

    !** Other fields in xtrainfo should have been set by the main prgram 
    !** before calling initout otherwise they will be meaning less
    Call datafile_writeExtraInfo(configfile)

  End Subroutine datafile_initout

  !----------------------------------------------------------
  ! Initializes some fields of the given configfile 
  !----------------------------------------------------------
  Subroutine datafile_initContents(configfile,ctrlfile, tag)
    Type(CONFILE), Intent(inout)       :: configfile
    Character(*), Intent(in)           :: ctrlfile,tag

    Character(len=5*strlen)            :: line 
    Character(len=strlen)              :: temp_field 
    Character(len=strlen),Dimension(20):: fields 
    Integer                            :: ctrlunit, lineno, nfields, i

    ctrlunit=file_open(ctrlfile)

    !** Find the General section
    lineno = filesrchstr(ctrlunit, tag, line, .True.)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' Could not find the datafile tag "',Trim(tag),'"'
      Write(0,'(1x,a)') 'in the control file.  This control file section is required'
      Stop
    End If

    Read(ctrlunit,'(a)') line
    line=adjustl(Trim(stripcmnt(line)))
    nfields=split(line,fields,",")
    !    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    !    Write(*,*) line, fields
    ! by default we will not store any of these
    configfile%is_storing%tot_nrg   = .False.
    configfile%is_storing%pair_nrg  = .False.
    configfile%is_storing%intra_nrg = .False.

    configfile%is_storing%posns     = .False.
    configfile%is_storing%vels      = .False.
    configfile%is_storing%time      = .False.
    Do i=1,nfields
      temp_field=Tolower(Adjustl( Trim(fields(i)) ))
      Select Case(temp_field)
      Case("energy")
        configfile%is_storing%tot_nrg   = .True.
      Case("intra_energy")
        configfile%is_storing%intra_nrg   = .True.
      Case("pair_energy")
        configfile%is_storing%pair_nrg   = .True.
      Case("velocity")
        configfile%is_storing%vels   = .True.
      Case("position")
        configfile%is_storing%posns   = .True.
      Case("time")
        configfile%is_storing%time   = .True.
      Case default
        Write(*,*) "Wrong field : "//Trim(temp_field)//&
            " specified in datafile section"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End Select
    End Do
  End Subroutine datafile_initContents

  !-------------------------------------------------------------------------
  ! Subroutine acts an interface
  ! Requires:  configfile -- configuration file data type
  !            filename -- the filename
  !-------------------------------------------------------------------------
  Subroutine datafile_initin(configfile,filename, is_old)
    Type(CONFILE), Intent(InOut)       :: configfile
    Character(*), Intent(In)           :: filename
    Logical,Intent(in), Optional :: is_old

    configfile%datafile_is_old=.False.
    If (present(is_old)) configfile%datafile_is_old=is_old
    If (configfile%datafile_is_old) Then
      configfile%name = filename
      Call datafile_initinOld(configfile)
    Else
      Call datafile_initinNew(configfile,filename)
    End If

  End Subroutine datafile_initin


  !---------------------------------------------------------
  ! Acts as in interface 
  !---------------------------------------------------------
  Subroutine datafile_readconfig(configfile, species, nmole_list, &
      spcstats, totnrg, time)
    Type(CONFILE), Intent(in)       :: configfile
    Type(AtMolCoords), Dimension(:), Intent(Inout)              :: species
    Integer, Dimension(:)           ::  nmole_list
    Type(Species_Stats), Intent(InOut) :: spcstats
    Real(Kind=RDbl), Intent(Out)           :: totnrg
    Real(Kind=RDbl), Intent(Out), Optional :: time
    If (configfile%datafile_is_old) Then
      Call datafile_readconfigOld(configfile, species, nmole_list, &
          spcstats, totnrg, time)
    Else
      Call datafile_readconfigNew(configfile, species, nmole_list, &
          spcstats, totnrg, time)
    Endif
  End Subroutine datafile_readconfig

  !----------------------------------------------------------
  ! Initializes the given configfile for Reading from it
  ! the disk file should have name "filename"
  !----------------------------------------------------------
  Subroutine datafile_initinNew(configfile,filename,nooutput)
    Type(CONFILE), Intent(inout)       :: configfile
    Character(*), Intent(in)           :: filename
    Logical, Intent(In), Optional      :: nooutput

    Type(FileHeader_Info)  :: header
    Integer                         :: funitno, ios

    configfile%name=Adjustl(Trim(filename))
    funitno=isfileopen(configfile%name)
    If (funitno>0) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "The datafile is already open. &
          & Should not try to re-init it again"
      Stop
    Else
      funitno = file_open(configfile%name,100)
    Endif

    configfile%unitno=funitno
    Call file_settag(configfile%name, d_dat_file)

    !** Read the header, storage_info
    If (Present(nooutput)) Then
      Call datafile_readheader(configfile,header,nooutput)
    Else
      Call datafile_readheader(configfile,header)
    End If

    Read(configfile%unitno) configfile%is_storing
    configfile%simno=header%simno
    Call datafile_readExtraInfo(configfile, ios)
    If (ios/=0) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Looks like datafile's Extra_Info section was modified after &
          & this *.con.* file was produced"
      Write(*,*) "Extra Info could be not read correctly"
      Write(*,*) "Try to recompile this code with an earlier version &
          & of datafile.F90 "
      Write(*,*) " Also make sure the MAX_EXTRA_INFO_LINES and MAX_SORBS are",&
          "same as the one used previously"
    Endif

  End Subroutine datafile_initinNew



  !-------------------------------------------------------
  ! Tries to write xtrainfo for the given configuration file
  ! It writes all fields of extrainfo, except for siminfostring 
  !-------------------------------------------------------
  Subroutine datafile_writeExtraInfo(configfile)
    Type(CONFILE), Intent(InOut)       :: configfile

    Character(len=strLen), Dimension(MAX_EXTRA_INFO_LINES) :: tempstrings
    Character(len=3*strLen) :: XlongString
    Character(len=strLen) :: blank

    Integer :: i,nspcs

    If (MAX_EXTRA_INFO_LINES<10+2*MAX_SORBS) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    blank=Repeat(" ",strLen)

    !** make all blank
    tempstrings=blank

    Write(tempstrings(1),'(i5,5x,a)')       configfile%xtrainfo%simno, &
        " # simno"
    Write(tempstrings(2),'(f12.5,5x,a)')    configfile%xtrainfo%simPressure, &
        " # simPressure"
    Write(tempstrings(3),'(f12.5,5x,a)') configfile%xtrainfo%simTemperature, &
        " # simTemperature"
    XlongString=blank
    XlongString=Trim(configfile%xtrainfo%simType)//"     "//" # simType"
    Write(tempstrings(4),'(a)') XlongString(1:strLen)

    Write(tempstrings(5),'(i5,5x,a)')       configfile%xtrainfo%nspcs, &
        " # nspcs"

    nspcs=configfile%xtrainfo%nspcs
    !** write species names
    Do i=1,MAX_SORBS
      If (i<=nspcs) Then
        XlongString=blank
        XlongString=Trim(configfile%xtrainfo%spcnamelist(i))
        Write(tempstrings(5+i),'(a,i5)')  XlongString(1:20)//" # spc",i
      Else
        Write(tempstrings(5+i),'(a,5x,a)')  "NOSPECIES", " # spc"
      Endif
    End Do

    !** write whether they are fixed or not
    Do i=1,MAX_SORBS
      If (i<=nspcs) Then
        If (configfile%xtrainfo%is_fixed(i)) Then
          Write(tempstrings(5+MAX_SORBS+i),'(a,5x,a,i5)')  &
              "FIXED", " # spc", i
        Else
          Write(tempstrings(5+MAX_SORBS+i),'(a,5x,a,i5)')  &
              "NOT_FIXED", " # spc", i
        Endif
      Else
        Write(tempstrings(5+MAX_SORBS+i),'(a,5x,a,i5)')  &
            "NOSPECIES", " # spc", i
      Endif
    End Do

    ! ** adding more things here
    Write(tempstrings(6+2*MAX_SORBS),'(i12,5x,a)')  configfile%niters, &
        " # niters"
    Write(tempstrings(7+2*MAX_SORBS),'(i12,5x,a)')  configfile%iconfig, &
        " # iconfig"
    !** write these strings to actual ctrlfile
    Do i=1,MAX_EXTRA_INFO_LINES
      Write(configfile%unitno) tempstrings(i)
    End Do

  End Subroutine datafile_writeExtraInfo



  !-------------------------------------------------------
  ! Tries to tread xtrainfo for the given configuration file
  ! If successful ios=zero
  !-------------------------------------------------------
  Subroutine datafile_readExtraInfo(configfile, ios )
    Type(CONFILE), Intent(InOut)       :: configfile
    Integer :: ios
    Character(len=strLen), Dimension(MAX_EXTRA_INFO_LINES) :: tempstrings
    Character(len=strLen)                                  :: temp1
    Integer :: i, error, nspcs

    !** Read info into the temperory strings array
    Do i=1,MAX_EXTRA_INFO_LINES
      Read(configfile%unitno,IOSTAT=ios) tempstrings(i)
      If (ios/=0) Return
    End Do


    !** Read actual info
    Write(*,*) "Reading info from the extra-info string array"
    Read(tempstrings(1),*)   configfile%xtrainfo%simno
    Read(tempstrings(2),*,IOSTAT=error)   configfile%xtrainfo%simPressure
    If (error/=0) Then
      Write(*,*) "#### string : ", Trim(tempstrings(3))
      Write(*,*) "#### Pressure written into datafile is not clear ####"
      Write(*,*) "#### So assigning a value of zero                ####"
      configfile%xtrainfo%simPressure=zero
    Endif
    Read(tempstrings(3),*,IOSTAT=error)   configfile%xtrainfo%simTemperature
    If (error/=0) Then
      Write(*,*) "#### string : ", Trim(tempstrings(3))
      Write(*,*) "#### Temperature written into datafile is not clear ####"
      Write(*,*) "#### So assigning a value of zero                  ####"
      configfile%xtrainfo%simTemperature=zero
    Endif
    Read(tempstrings(4),*)   configfile%xtrainfo%simType
    Read(tempstrings(5),*)   configfile%xtrainfo%nspcs

    nspcs=configfile%xtrainfo%nspcs

    !** write species names
    Do i=1,MAX_SORBS
      If (i<=nspcs) Then
        Read(tempstrings(5+i),'(a)') temp1
        configfile%xtrainfo%spcnamelist(i)= cleanstring(temp1)
      Else
        configfile%xtrainfo%spcnamelist(i)="NOSPECIES"
      Endif
    End Do

    !** write whether they are fixed or not
    Do i=1,MAX_SORBS
      If (i<=nspcs) Then
        Read(tempstrings(5+MAX_SORBS+i),'(a)')   temp1
        If (cleanstring(temp1)=="FIXED") Then
          configfile%xtrainfo%is_fixed(i)=.true.
        else
          configfile%xtrainfo%is_fixed(i)=.False.
        Endif
      Else
        ! nothing here
      Endif
    End Do

    Read(tempstrings(6+2*MAX_SORBS),*)  configfile%niters
    Read(tempstrings(7+2*MAX_SORBS),*)  configfile%iconfig

  End Subroutine datafile_readExtraInfo


  !-------------------------------------------------------
  ! Close the given configuration file
  !-------------------------------------------------------
  Subroutine datafile_close(configfile)
    Type(CONFILE), Intent(InOut)       :: configfile

    Close(unit = configfile%unitno)
    configfile%unitno = -1

  End Subroutine datafile_close


  !----------------------------------------------------------------------------
  ! Returns the simulation number for the given config file
  !----------------------------------------------------------------------------
  Integer Function datafile_getSimNum(obj)
    Type(CONFILE), Intent(In) :: obj
    If (obj%datafile_is_old) Then
      datafile_getSimNum = obj%simno
    Else
      datafile_getSimNum = obj%xtrainfo%simno
    Endif
  End Function datafile_getSimNum

  !----------------------------------------------------------------------------
  ! Returns the simulation extra-info for the given config file
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function datafile_getExtraInfo(obj, typestring)
    Type(CONFILE), Intent(In) :: obj
    Character(*),Intent(in)     :: typestring

    Select Case(cleanstring(typestring))
    Case('pressure')
      datafile_getExtraInfo=obj%xtrainfo%simPressure
    Case('temperature')
      datafile_getExtraInfo=obj%xtrainfo%simTemperature 
    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select
  End Function datafile_getExtraInfo


  !----------------------------------------------------------------------------
  ! Sets the simulation extra-info for the current simulation 
  !----------------------------------------------------------------------------
  Subroutine datafile_setExtraInfo(obj, typestring, value)
    Type(CONFILE), Intent(Inout) :: obj
    Character(*)     :: typestring
    Real(kind=RDbl) :: value

    Select Case(cleanstring(typestring))
    Case('pressure')
      obj%xtrainfo%simPressure=value
    Case('temperature')
      obj%xtrainfo%simTemperature=value
    Case default
      Write(*,*) " Cant recognise the string : ", cleanstring(typestring)
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select
  End Subroutine datafile_setExtraInfo

  !------------------------------------------------------------
  ! Write the header of the data file linked to ctrl_filename
  ! The header contains the control file that was used to 
  !generate it and the current simulation number "simno"
  !------------------------------------------------------------
  Subroutine datafile_writeheader(configfile, ctrl_filename)
    Type(CONFILE), Intent(in)       :: configfile
    Character(*), Intent(in) :: ctrl_filename

    Type(FileHeader_Info)    :: headerinfo
    Character(len=lstrLen)   :: inputbuffer
    Integer                  :: nlines, error, ctrlunitno

    !** Open the control file if not already open
    ctrlunitno = isfileopen(ctrl_filename)
    If (ctrlunitno < 0) Then
      ctrlunitno = file_getunit(ctrl_filename)
      Open(unit = ctrlunitno, file=ctrl_filename, IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            "  Could not open file ", Trim(ctrl_filename)
        Stop
      End If
    End If

    Rewind(ctrlunitno)

    !** Get the lines of the control file and fill the headerinfo
    !** data structure
    nlines = 0
    Do 
      Read(ctrlunitno,'(a)',IOSTAT=error) inputbuffer
      If (error /= 0) Then
        Exit
      Endif
      nlines = nlines + 1
      headerinfo%ctrl_line(nlines) = inputbuffer
    End Do
    headerinfo%nlines = nlines
    headerinfo%ctrl_filename = ctrl_filename
    headerinfo%simno = configfile%simno
    headerinfo%siminfo_string = Trim(genparams%siminfo_string)

    !** Write the lines to the configuration file
    Rewind(configfile%unitno)
    Write(configfile%unitno) headerinfo

  End Subroutine datafile_writeheader


  !---------------------------------------------------------
  ! Read the configuration information from the data file
  ! configfile should be initialised
  !---------------------------------------------------------
  Subroutine datafile_readconfigNew(configfile, species, nmole_list, &
      spcstats, totnrg, time)
    Type(CONFILE), Intent(in)       :: configfile
    Type(AtMolCoords), Dimension(:), Intent(Inout)              :: species
    Integer, Dimension(:)           ::  nmole_list
    Type(Species_Stats), Intent(InOut) :: spcstats
    Real(Kind=RDbl), Intent(Out)           :: totnrg

    Real(Kind=RDbl), Intent(Out), Optional :: time

    Real(Kind=RDbl) :: temptime

    Integer :: unitno, i, natoms,nmoles

    unitno = isfileopen(configfile%name)
    If (configfile%unitno/=unitno) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(a)')"configfile shoud be initialised"
      Stop
    Endif

    If (configfile%is_storing%time) Then
      Read(unitno) temptime
      If (Present(time)) time=temptime
    Else If (Present(time)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Error, You want time returned, but this datafile does not store it"
      Write(0,*) "Data file (*.con.*) :"//Trim(configfile%name)
      Stop
    Endif

    !** Read the energy from the configuration file
    Call datafile_readnrgNew(configfile, spcstats, totnrg)

    Do i = 1, molecules_getnsorbs()
      natoms = config_getnatoms(species,i)    
      nmole_list(i)=0

      If (config_isfixed(species(i))) Cycle
      Read(unitno) nmoles

      nmole_list(i)=nmoles
      If (nmoles==0) Cycle
      If (nmoles>Size(species(i)%coords,2)) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Size of species : ",Size(species(i)%coords,2)
        Write(*,*) "spctype ", i
        Write(*,*) "nmoles", nmoles
        Write(*,*) " Species is not big enough to read the information &
            & in the datafile, The Calling program should allocate more &
            & memory for species structure"
        Stop
      Endif

      If (configfile%is_storing%posns) &
          Read(unitno) species(i)%coords(1:natoms,1:nmoles)%rp

      If (configfile%is_storing%vels) &
          Read(unitno) species(i)%coords(1:natoms,1:nmoles)%v

    End Do

  End Subroutine datafile_readconfigNew


  !---------------------------------------------------------
  ! Read the configuration of a single moleucle from control file
  ! returns its intramoleuclar energy back 
  ! spc -- the index in the species array
  ! datfilespc -- the index used in original simulation, could be different
  !---------------------------------------------------------
  Subroutine datafile_readOneMolec(configfile, species, datfilespc, intranrg)
    Type(CONFILE), Intent(in)       :: configfile
    Type(AtMolCoords),  Intent(Inout)              :: species
    Integer, Intent(in) :: datfilespc
    Real(Kind=RDbl), Intent(Out)           :: intranrg

    Real(Kind=RDbl) :: temptime,totalu

    Type(VecType) , Dimension(:,:) ,Pointer ,Save :: tempRVecArray 
    Type(VecType) , Dimension(:,:) ,Pointer ,Save :: tempVVecArray 
    Integer :: unitno, i, natoms,nmoles, spctype, maxnatoms, maxnmoles
    Integer :: nspcs, error

    Character(len=strlen) :: spcname
    Logical, Save :: first_time=.True.

    If (first_time) Then 
      first_time=.False.
      maxnatoms=1
      maxnmoles=1
      Allocate(tempRVecArray(maxnatoms, maxnmoles), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Allocate(tempVVecArray(maxnatoms, maxnmoles), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Endif

    unitno = isfileopen(configfile%name)
    If (configfile%unitno/=unitno) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(a)')"configfile shoud be initialised"
      Stop
    Endif

    If (configfile%is_storing%time) Then
      Read(unitno) temptime
    Endif

    !** Read the energy from the configuration file
    Call datafile_readnrgPartial(configfile, datfilespc, intranrg, totalu)

    nspcs=configfile%xtrainfo%nspcs 
    Do i = 1, nspcs
      spcname=configfile%xtrainfo%spcnamelist(i)
      spctype=molecules_gettype(spcname)
      natoms = molecules_getnatoms(spctype)

      If (configfile%xtrainfo%is_fixed(i)) Cycle
      Read(unitno) nmoles

      If ((natoms>maxnatoms) .Or. (nmoles>maxnmoles) ) Then
        maxnatoms=natoms
        maxnmoles=nmoles

        Deallocate(tempRVecArray, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        Allocate(tempRVecArray(maxnatoms, maxnmoles), STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        Deallocate(tempVVecArray, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        Allocate(tempVVecArray(maxnatoms, maxnmoles), STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Endif
!!$      If (datfilespc==i) Then
!!$        If (nmoles/=1) Then
!!$          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$          Write(*,*) "as of now only single molecule &
!!$              & configurations can be read in"
!!$          Stop
!!$        Endif
!!$      Endif

      If (nmoles==0) Cycle

      If (configfile%is_storing%posns) &
          Read(unitno) tempRVecArray(1:natoms,1:nmoles)

      If (configfile%is_storing%vels) &
          Read(unitno) tempVVecArray(1:natoms,1:nmoles)

      If (datfilespc==i) Then
        species%coords(1:natoms, 1)%rp=tempRVecArray(1:natoms,1)
        species%coords(1:natoms, 1)%v=tempVVecArray(1:natoms,1)
      Endif
    End Do

  End Subroutine datafile_readOneMolec

  !-----------------------------------------------------------------------------
  ! Read the configuration and energy of a single species from control file.
  ! Requires:  configfile -- configuration file data type
  !            rp -- returned atomic parent coordinates           
  !            v -- returned atomic velocities
  !            datfilespc -- species number in datafile (spc # in orig sim!)
  !            intranrg -- summed intramolecular energy for all molecules of spc
  !            totalnrg -- summed noncoul+coul energy for all molecules of spc
  !            num_molecs -- returned number of molecules
  ! NOTE: 'totalnrg' may not be available, in which case zero is returned
  !-----------------------------------------------------------------------------
  Subroutine datafile_readOneSpc(configfile, rp, v, datfilespc, &
      intranrg, totalnrg, num_molecs)
    Type(CONFILE), Intent(In)                :: configfile
    Integer, Intent(In)                      :: datfilespc
    Real(Kind=RDbl), Intent(Out)             :: intranrg,totalnrg
    Type(VecType), Dimension(:,:), Pointer   :: rp, v
    Integer, Intent(Out)                     :: num_molecs

    Integer               :: unitno, i, natoms,nmoles, spctype
    Integer               :: nspcs, error, current_asize, current_msize
    Logical               :: checked_sizes
    Real(kind=RDbl)       :: temptime
    Character(len=strlen) :: spcname

    !** Later used to make sure that size of rp and v are enough to read 
    !** the whole Data
    checked_sizes = .False.

    unitno = isfileopen(configfile%name)
    If (configfile%unitno /= unitno) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(a)') "configfile must be initialised"
      Stop
    End If

    If (configfile%is_storing%time) Then
      Read(unitno) temptime
    End If

    !** Read the energy from the configuration file
    Call datafile_readnrgPartial(configfile, datfilespc, intranrg, totalnrg)

    nspcs = configfile%xtrainfo%nspcs 
    Do i = 1,nspcs
      spcname=configfile%xtrainfo%spcnamelist(i)
      spctype=molecules_gettype(spcname)
      If (spctype==0) Then
        Write(*,*) "this molecules is not initialized : ", Trim(spcname)
        Write(*,*) " add it to moelcule list"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      End If
      natoms = molecules_getnatoms(spctype)

      If (configfile%xtrainfo%is_fixed(i)) Cycle
      Read(unitno) nmoles
      If (nmoles==0) Cycle
      If (datfilespc==i) Then
        num_molecs=nmoles

        !** make sure rp and v are big enough
        If (.Not.checked_sizes) Then

          If (Associated(rp)) Then
            current_asize=Size(rp,1)
            current_msize=Size(rp,2)

            If ( (num_molecs>current_msize).or. (natoms>current_asize) ) Then
              Nullify(rp,v)
              Allocate(rp(natoms,num_molecs), v(natoms,num_molecs),STAT=error) 
              If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
            Endif

          Else

            Allocate(rp(natoms,num_molecs), v(natoms,num_molecs),STAT=error) 
            If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
          Endif

          checked_sizes=.True.
        Endif

        If (datafile_hasposition(configfile)) Then
          Read(unitno) rp(1:natoms,1:nmoles)
        Else
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
        If (datafile_hasvelocity(configfile)) Then
          Read(unitno) v(1:natoms,1:nmoles)
        Else
          v(1:natoms,1:nmoles)=vector_zerovec()
        Endif


      Endif
    End Do

  End Subroutine datafile_readOneSpc

  !---------------------------------------------------------
  ! Write the configuration information to the data file
  ! pairwise energies and Total energies are written
  ! Total energy includes Kinetic Energy also.( For details see writenrg)
  !---------------------------------------------------------
  Subroutine datafile_writeconfig(species,spcstats,configfile,time)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(Species_Stats), Intent(In)             :: spcstats
    Type(confile), Intent(In)                   :: configfile
    Real(Kind=RDbl), Intent(In), Optional       :: time

    Integer ::  i, natoms, unitno,nmoles

    unitno = file_checkfileopen(configfile%name)
    If (configfile%unitno /= unitno) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Configuration file must be initialized"
      Write(*,*) configfile%unitno,unitno
      Write(*,*) configfile%name
      Stop
    End If

    If (present(time)) Then
      If (configfile%is_storing%time) Then
        Write(unitno) time
      Else 
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) " Passed time here, but not expected to store it. &
            & is this OK?"
      End If
    Else
      If (configfile%is_storing%time) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) " Did not pass time here, but expected to store it(?). &
            & This is not OK"
        Write(*,*) " Check your simulation type etc."
        Stop
      Endif
    Endif

    Call datafile_writenrg(configfile, spcstats)

    Do i = 1, molecules_getnsorbs()
      natoms = config_getnatoms(species,i)
      nmoles = config_getnmoles(species,i)
      If (config_isfixed(species(i))) Cycle
      If (natoms == 0) Cycle
      Write(unitno) nmoles
      If (nmoles == 0) Cycle
      If (configfile%is_storing%posns) &
          Write(unitno) species(i)%coords(1:natoms,1:nmoles)%rp
      If (configfile%is_storing%vels) &
          Write(unitno) species(i)%coords(1:natoms,1:nmoles)%v
    End Do
  End Subroutine datafile_writeconfig

  !--------------------------------------------------------
  ! Generates the control file that was used to create the
  ! the configuration file from the parameters. It takes
  ! unit no. of the configuration file
  !--------------------------------------------------------
  Subroutine datafile_gencontrolfile(optfilename,configfile)
    Character(len=strLen), Intent(in), Optional :: optfilename
    Type(confile), Intent(in)       :: configfile

    Integer      :: i, ctrlunitno
    Character(len=strLen)   :: filename
    Type(FileHeader_Info) :: header
    Logical :: nooutput


    !get the header info
    nooutput = .True.         ! suppress writing info to screen again
    Rewind(configfile%unitno)
    Call datafile_readheader(configfile,header, nooutput)

    !** Decide the filename
    If (Present(optfilename)) Then
      filename = optfilename
    Else
      filename = Trim(header%ctrl_filename)//".duplicate"
    End If

    ctrlunitno = file_getunit(filename)
    Open(unit=ctrlunitno, file=filename)
    Do i=1, header%nlines
      Write(ctrlunitno, '(a)') Trim(header%ctrl_line(i))
    End Do
    Close(ctrlunitno)
  End Subroutine datafile_gencontrolfile

  !--------------------------------------------------------
  ! Reads the header written at the top of the given control file
  ! the configfile has to be initialised before calling this
  !--------------------------------------------------------
  Subroutine datafile_readheader(configfile,header,nooutput)
    Type(CONFILE), Intent(in)       :: configfile
    Type(FileHeader_Info),Intent(out)  :: header
    Logical, Optional, Intent(In)   :: nooutput

    Logical :: suppressOutput 
    Integer :: unitno

    !** Check for optional nooutput
    suppressOutput = .False.
    If (Present(nooutput)) suppressOutput = nooutput

    unitno = isfileopen(configfile%name)
    If (configfile%unitno/=unitno) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*)"ctrlfile shoud be initialised"
      Stop
    End If

    Rewind(configfile%unitno)

    Read(configfile%unitno) header

    If (suppressOutput) Return

    Write(*,*) " Reading the Header of config file : ",Trim(configfile%name)
    Write(*,*) " INFO STRING, about this config file :--"
    Write(*,'(a)') Trim(header%siminfo_string)
    Write(*,*)

  End Subroutine datafile_readheader

  !----------------------------------------------------------------------
  ! Writes all the intra pair enrgies, and total enrg to the config file
  !----------------------------------------------------------------------
  Subroutine datafile_writenrg(configfile, spcstats)
    Type(CONFILE), Intent(in)       :: configfile
    Type(Species_Stats), Intent(In)   :: spcstats            

    Integer           :: i,j, unitno
    Real(kind = RDbl) :: total_u,unoncoul,ucoul,kin_nrg,intra_nrg

    unitno=configfile%unitno

    total_u=0.0
    Do i = 1, molecules_getnsorbs()
      Do j = i, molecules_getnsorbs()

        !** Write the potential energies
        unoncoul = storestats_getnoncoul(spcstats,i,j,'inst')
        ucoul = storestats_getcoul(spcstats,i,j,'inst')
        total_u = total_u+unoncoul+ucoul
        If (configfile%is_storing%pair_nrg) Write(unitno) unoncoul,ucoul
      End Do

      !** Write the kinetic and total intramolecular energies
      kin_nrg = storestats_getke(spcstats,i,'inst')
      intra_nrg = storestats_getintranrg(spcstats,i,'inst')
      total_u = total_u+kin_nrg+intra_nrg
      Write(unitno) kin_nrg,intra_nrg

      !** Write individual intra nrgs
      If (configfile%is_storing%intra_nrg) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            ' not prepared to write individual intramolecular energies'
        Stop
      End If
    End Do

    !** Write the total energy
    If (configfile%is_storing%tot_nrg) Then
      Write(unitno) total_u
    End If

  End Subroutine datafile_writenrg

  !-----------------------------------------------------------------
  ! Updates a species statistics data type with all the intra pair, 
  ! intra and kinetic and total energies from the config file.  
  ! Returns the total energy.
  ! Requires:  configfile -- configuration file data type
  !            spcstats -- species statistics to update
  !            total_u -- returned total energy
  !-----------------------------------------------------------------
  Subroutine datafile_readnrgNew(configfile, spcstats, total_u)
    Type(CONFILE), Intent(In)            :: configfile    
    Type(Species_Stats), Intent(InOut)   :: spcstats            
    Real(kind = RDbl),Intent(Out)        :: total_u

    Integer         :: i,j,unitno
    Real(Kind=RDbl) :: noncoulEnergy, coulEnergy, ke, intraEnergy

    unitno=configfile%unitno
    total_u=0.0

    Do i = 1, molecules_getnsorbs()
      Do j = i, molecules_getnsorbs()
        noncoulEnergy = zero
        coulEnergy = zero

        If (configfile%is_storing%pair_nrg) &
            Read(unitno) noncoulEnergy, coulEnergy

        !** Update the energy in the storage structure
        Call storestats_updateEnergySS(spcstats,i,j,'NONCOUL',noncoulEnergy)
        Call storestats_updateEnergySS(spcstats,i,j,'COUL',coulEnergy)
      End Do

      !** Read the kinetic and intramolecular energies
      Read(unitno) ke, intraEnergy

      !** Update the kinetic and intramolecular energies in spcstats
      Call storestats_updatenrg(spcstats,i,'KE',ke)
      Call storestats_updatenrg(spcstats,i,'INTRA',intraEnergy)
      ! individual intra energies
      If (configfile%is_storing%intra_nrg) Then 
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            ' not prepared to read individual intramolecular energies'
        Stop
      End If

    End Do

    !** Read the total energy
    If (configfile%is_storing%tot_nrg) Read(unitno) total_u

  End Subroutine datafile_readnrgNew

  !-----------------------------------------------------------------------
  ! Reads and returns intra energy and total energies stored 
  ! for one species.
  ! Requires:  configfile -- configuration file data type
  !            spc -- species number IN THE DATAFILE (matching orig sim)
  !            intranrg -- summed intramolecular energy for spc
  !            totalnrg -- summed noncoul+coul energy for spc
  ! NOTE: 'totalnrg' may not be available, in which case zero is returned
  !-----------------------------------------------------------------------
  Subroutine datafile_readnrgPartial(configfile, spc, intranrg, totalnrg)
    Type(CONFILE), Intent(In)          :: configfile    
    Integer, Intent(In)                :: spc
    Real(kind=RDbl), Intent(Out)       :: intranrg, totalnrg

    Integer            :: i,j,unitno,nspcs 
    Real(Kind=RDbl)    :: noncoulEnergy, coulEnergy, ke, temp_intra, total_u

    unitno = configfile%unitno

    !** nspcs of original simulation
    nspcs = configfile%xtrainfo%nspcs

    !** Make sure species number and total match
    If ((spc > nspcs) .Or. (spc < 1))Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "wrong value of spc", spc
      Stop
    End If

    !** Read the energies at the beginning of an entry
    totalnrg = zero
    Do i = 1,nspcs
      Do j = i,nspcs 
        noncoulEnergy = zero
        coulEnergy = zero

        If (configfile%is_storing%pair_nrg) &
            Read(unitno) noncoulEnergy, coulEnergy
      End Do

      !** Sum all noncoul + coul energies for the desired species number
      If (i == spc) totalnrg = totalnrg + noncoulEnergy + coulEnergy

      !** Read the kinetic and intramolecular energies
      Read(unitno) ke, temp_intra
      If (i == spc) intranrg = temp_intra

      !** Read the individual intra energies
      If (configfile%is_storing%intra_nrg) Then 
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            ' not prepared to read individual intramolecular energies'
        Stop
      End If
    End Do

    !** Read the total energy
    If (configfile%is_storing%tot_nrg) Read(unitno) total_u

  End Subroutine datafile_readnrgPartial

  !----------------------------------------------------------------------------
  ! Returnss the number of configurations stored 
  !----------------------------------------------------------------------------
  Integer Function datafile_getnconfs(obj)
    Type(CONFILE), Intent(In) :: obj
    datafile_getnconfs = obj%niters / obj%iconfig
  End Function datafile_getnconfs

  !----------------------------------------------------------------------------
  ! Returns TRUE if time is recorded in the config file.
  !----------------------------------------------------------------------------
  Logical Function datafile_hastime(obj)
    Type(CONFILE), Intent(In) :: obj
    datafile_hastime = obj%is_storing%time 
  End Function datafile_hastime

  !----------------------------------------------------------------------------
  ! Returns TRUE if positions are recorded in the config file.
  !----------------------------------------------------------------------------
  Logical Function datafile_hasposition(obj)
    Type(CONFILE), Intent(In) :: obj

    If (obj%datafile_is_old) Then
      datafile_hasposition = .False.
      If (obj%content_tag>2)         datafile_hasposition = .True.
    Else
      datafile_hasposition = obj%is_storing%posns 
    Endif
  End Function datafile_hasposition

  !----------------------------------------------------------------------------
  ! Returns TRUE if veocities are recorded in the config file.
  !----------------------------------------------------------------------------
  Logical Function datafile_hasvelocity(obj)
    Type(CONFILE), Intent(In) :: obj
    datafile_hasvelocity = obj%is_storing%vels 
  End Function datafile_hasvelocity





  !----------------------------------------------------------------------------
  ! Writes a sample of the required input for this module to unit unitno
  !----------------------------------------------------------------------------
  Subroutine datafile_sampleCF(unitno)
    Integer, Intent(In) :: unitno
    Write(unitno,'(a)') "--------  "//Trim(default_datafile_tag)//" --------"
    Write(unitno,'(a)') " Energy, intra_energy, position, Velocity, pair_energy, time "//" # contents of datafile"
    Write(unitno,'(2a)') '# Above fields are case-insensitive and order ', &
        'does not matter. '
    Write(unitno,'(a)') "# Delete words corresponding to variables &
        & you dont want stored."
    Write(unitno,*) 
  End Subroutine datafile_sampleCF

  !----------------------------------------------------------------------------
  ! ALL THE ROUTINES BELOW CAN BE DISCARDED, IF WE DONT WANT TO DO POST-CODE OF
  ! OLD DATAFILES ANYMORE
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Reads the time from the datafile OBSOLETE-should be deleted later
  !----------------------------------------------------------------------------
  Subroutine datafile_readtime(unitno,contentTag,time)
    Integer, Intent(In) :: unitno
    Integer, Intent(In) :: contentTag
    Real(Kind=RDbl), Intent(Out) :: time

    !** Check the content tag to make sure we should be doing this
    Select Case (contentTag)
    Case(1)
    Case(2)
    Case(3)
    Case(4)
    Case(5)
      Read(unitno) time
    Case Default
      Write(0,'(1x,2a,i4,a,i3,a)') __FILE__," : ",__LINE__, &
          " Specified datafile content tag ",contentTag," invalid"
      Stop
    End Select

  End Subroutine datafile_readtime

  !----------------------------------------------------------
  ! left here ability to read old *.con.* files
  ! Initializes the given file configfile for Reading from it
  !needs only configfile%name to be specified
  ! gen_init is the flag which shows genparams are initialised
  ! or not, configfile%content_tag can be got only if gen_init=1
  ! The header is read here; next "Read" will read configurations
  !----------------------------------------------------------
  Subroutine datafile_initinOLD(configfile)
    Type(CONFILE), Intent(inout)       :: configfile

    Type(FileHeader_Info)  :: header
    Integer                         :: funitno,error

    funitno = isfileopen(configfile%name)
    If (funitno < 0) Then
      funitno = file_getunit(configfile%name)
      Open(file=configfile%name,unit=funitno, status='unknown', &
          form='unformatted',IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            "  Could not Open file ", Trim(configfile%name)
        Stop
      End If
    Else
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "datafile already initialised" 
      Stop
    End If

    configfile%unitno=funitno
    Call file_settag(configfile%name, d_dat_file)
    Call datafile_readheader(configfile,header)
    configfile%simno=header%simno

  End Subroutine datafile_initinOLD
  
  !--------------------------------------------------------
  ! sets the content-tag by reading from general
  !--------------------------------------------------------
  Subroutine datafile_setContentTag(obj)
    Type(CONFILE), Intent(inout)       :: obj
    obj%content_tag=general_getContentTag()
  End Subroutine datafile_setContentTag

  !---------------------------------------------------------
  ! Read the configuration information from the data file
  ! configfile should be initialised
  !---------------------------------------------------------
  Subroutine datafile_readconfigOld(configfile, species, nmole_list, &
      spcstats, totnrg, time)
    Type(CONFILE), Intent(in)       :: configfile
    Type(AtMolCoords), Dimension(:), Intent(Inout)              :: species
    Integer, Dimension(:)           ::  nmole_list
    Type(Species_Stats), Intent(InOut) :: spcstats
    Real(Kind=RDbl), Intent(Out)           :: totnrg
    Real(Kind=RDbl), Intent(Out), Optional :: time
    Integer:: content
    Integer :: unitno, i, natoms,nmoles

    content=configfile%content_tag
    unitno = isfileopen(configfile%name)
    If (configfile%unitno/=unitno) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(a)')"configfile shoud be initialised"
      Stop
    Endif

    If ((.Not.(Present(time))).And.(content==5)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Error, this routine requires time to be passed"
      Stop
    Endif

    If (Present(time)) Then
      Call datafile_readtime(unitno, configfile%content_tag,time)
    End If

    !** Read the energy from the configuration file
    !LC    Call datafile_readnrg(unitno, configfile%content_tag, energy, totnrg )
    Call datafile_readnrgOld(spcstats, unitno, configfile%content_tag, totnrg)

!!$    totnrg = 0.0_RDbl

    Do i = 1, molecules_getnsorbs()
      natoms = config_getnatoms(species,i)    
      nmole_list(i)=0

      If (config_isfixed(species(i))) Cycle
      Read(unitno) nmoles
      nmole_list(i)=nmoles
      If (nmoles==0) Cycle
      If (nmoles>Size(species(i)%coords,2)) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Size of species : ",Size(species(i)%coords,2)
        Write(*,*) "spctype ", i
        Write(*,*) "nmoles", nmoles
        Write(*,*) " Species is not big enough to read the information &
            & in the datafile, The Calling program should allocate more &
            & memory for species structure"
        Stop
      Endif
      Select Case (configfile%content_tag)
      Case (1)
      Case (2)
      Case (3)
        Read(unitno) species(i)%coords(1:natoms,1:nmoles)%rp
      Case (4)
        Read(unitno) species(i)%coords(1:natoms,1:nmoles)%rp
        Read(unitno) species(i)%coords(1:natoms,1:nmoles)%v
      Case (5)
        Read(unitno) species(i)%coords(1:natoms,1:nmoles)%rp
        Read(unitno) species(i)%coords(1:natoms,1:nmoles)%v
      Case Default
        Write(0,'(1x,2a,i4,a,i3,a)') __FILE__," : ",__LINE__, &
            " Specified datafile content tag ", configfile%content_tag, &
            " invalid"
        Stop
      End Select

    End Do
  End Subroutine datafile_readconfigOld

  !----------------------------------------------------------------------------
  ! Writes the time to the datafile
  ! --------------OBSOLETE----------------
  !----------------------------------------------------------------------------
  Subroutine datafile_writetime(unitno,contentTag,time)
    Integer, Intent(In) :: unitno
    Integer, Intent(In) :: contentTag
    Real(Kind=RDbl), Intent(In) :: time

    !** Check the content tag to make sure we should be doing this
    Select Case (contentTag)
    Case(1)
    Case(2)
    Case(3)
    Case(4)
    Case(5)
      Write(unitno) time
    Case Default
      Write(0,'(1x,2a,i4,a,i3,a)') __FILE__," : ",__LINE__, &
          " Specified datafile content tag ",contentTag," invalid"
      Stop
    End Select

  End Subroutine datafile_writetime

  !----------------------------------------------------------------------
  ! Writes all the intra pair enrgies, and total enrg to the config file
  !----------------------------------------------------------------------
  Subroutine datafile_writenrgOld(spcstats,unitno,tag)
    Type(Species_Stats), Intent(In)   :: spcstats            
    Integer, Intent(In)               :: unitno
    Integer, Intent(In)               :: tag        

    Integer           :: i,j    
    Real(kind = RDbl) :: total_u,unoncoul,ucoul,kin_nrg,intra_nrg

    total_u=0.0
    Do i = 1, molecules_getnsorbs()
      Do j = i, molecules_getnsorbs()

        !** Write the potential energies
        unoncoul = storestats_getnoncoul(spcstats,i,j,'inst')
        ucoul = storestats_getcoul(spcstats,i,j,'inst')
        total_u = total_u+unoncoul+ucoul
        Select Case(tag)
        Case(1)
        Case(2)
          Write(unitno) unoncoul,ucoul
        Case(3)
          Write(unitno) unoncoul,ucoul
        Case(4)
          Write(unitno) unoncoul,ucoul
        Case(5)
          Write(unitno) unoncoul,ucoul
        End Select
      End Do

      !** Write the kinetic and intramolecular energies
      kin_nrg = storestats_getke(spcstats,i,'inst')
      intra_nrg = storestats_getintranrg(spcstats,i,'inst')
      total_u = total_u+kin_nrg+intra_nrg
      Write(unitno) kin_nrg,intra_nrg

    End Do

    !** Write the total energy
    Write(unitno) total_u

  End Subroutine datafile_writenrgOld

  !--------------------------------------------------------
  ! Reads all the intra pair, intra  and kinetic and total energies
  ! from the config file
  !--------------------------------------------------------
  Subroutine datafile_readnrgOld(spcstats,unitno,tag,total_u)
    Type(Species_Stats), Intent(InOut)   :: spcstats            
    Integer,Intent(in)::unitno
    Integer,Intent(in)::tag
    !LC    Type(Species_System_Energies), Dimension(:), Pointer :: energy
    Real(kind = RDbl),Intent(out)                     :: total_u
    Real(Kind=RDbl) :: noncoulEnergy, coulEnergy, ke, intraEnergy

    Integer:: i,j

    total_u=0.0

    Do i = 1, molecules_getnsorbs()
      Do j = i, molecules_getnsorbs()
        !** Read the potential energy
        Select Case(tag)
        Case(1)
          !LC          energy(i)%noncoulnrg(j)=zero
          !LC          energy(i)%noncoulnrg(j)=zero
          noncoulEnergy = zero
          coulEnergy = zero
        Case(2)
          !LC          Read(unitno) energy(i)%noncoulnrg(j), energy(i)%coulnrg(j)
          Read(unitno) noncoulEnergy, coulEnergy
        Case(3)
          !LC          Read(unitno) energy(i)%noncoulnrg(j), energy(i)%coulnrg(j)
          Read(unitno) noncoulEnergy, coulEnergy
        Case(4)
          !LC          Read(unitno) energy(i)%noncoulnrg(j), energy(i)%coulnrg(j)
          Read(unitno) noncoulEnergy, coulEnergy
        Case(5)
          !LC          Read(unitno) energy(i)%noncoulnrg(j), energy(i)%coulnrg(j)
          Read(unitno) noncoulEnergy, coulEnergy
        End Select

        !** Update the energy in the storage structure
        Call storestats_updateEnergySS(spcstats,i,j,'NONCOUL',noncoulEnergy)
        Call storestats_updateEnergySS(spcstats,i,j,'COUL',coulEnergy)

      End Do

      !** Read the kinetic and intramolecular energies
      !LC      Read(unitno) energy(i)%ke, energy(i)%intra(TOTAL_INDEX)
      Read(unitno) ke, intraEnergy

      !** Update the kinetic and intramolecular energies in spcstats
      Call storestats_updatenrg(spcstats,i,'KE',ke)
      Call storestats_updatenrg(spcstats,i,'INTRA',intraEnergy)

    End Do

    !** Read the total energy
    read(unitno) total_u

  End Subroutine datafile_readnrgOld

  !--------------------------------------------------------------------------
  ! Cleans the configuration file data type
  ! Requires:  configfile -- the object to be intialized
  !--------------------------------------------------------------------------
  Subroutine datafile_clean(configfile)
    Type(CONFILE), Intent(InOut)                    :: configfile

    !** no pointers or allocated types

  End Subroutine datafile_clean

End Module datafile
