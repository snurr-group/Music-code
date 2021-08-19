!------------------------------------------------------------------------------
! This module primarily handles the assignment of unit numbers to specific
! filenames and the opening of files.  The data structures here keep track
! of a file name and "tag" for each open file.  The idea is to prevent the
! use of already used unit numbers and allow easy access to already open files
! 
! The PROCEDURE FOR OPENING A FILE is as follows:
!    get the unit number and open the file simultaneously using:
!      unitno = file_open(filename,file_status,tag)
!    Where:
!      file_status -- a 3-bit integer that giving the expected/desired status
!      tag -- an optional tag to give the file, to simplify later access
!    The routine uses the data structure to keep track of used units and
!    will find the existing unit number if the file is already open
!
! Needed Improvements:
! 1) There is no provision for opening and getting the unit number of 
!    a file using just the tag.  This means that it's a pain to use tags
!    in comparison to simply using the filename and file_open
! 2) Consider removing the optional status of the fileopt in file_open,
!    this would make it safer.
!------------------------------------------------------------------------------

Module file

  Use defaults, Only: strLen, lstrLen, MAX_FILES, FIRST_UNIT

  Implicit None
  Save

  Private
  Public :: file_getunit, file_gettype, file_settag, file_open, file_getname, &
      file_checkfileopen, file_close, file_append

  Type Filename_Info
    Character(len=strLen)        :: name,tag
  End Type Filename_Info

  Type Filename_Set
    Type (Filename_Info), Dimension(MAX_FILES)  :: unit
    Integer                                     :: nunits
  End Type Filename_Set

  Type (Filename_Set)                           :: files

  Interface file_open
    Module Procedure file_openRegular
    Module Procedure file_openScratch
  End Interface


Contains

  !----------------------------------------------------------------------------
  ! Returns the filename, given the unit number
  !----------------------------------------------------------------------------
  Function file_getname(unit)
    Character(len=strLen) :: file_getname
    Integer, Intent(IN)                      :: unit

    file_getname = files%unit(unit)%name
  End Function file_getname

  !----------------------------------------------------------------------------
  ! Returns a unique unit number for a given filename
  !----------------------------------------------------------------------------
  Integer Function file_getunit(filename)
    Character(*), Intent(IN)                 :: filename
    Integer                                  :: unit
    Logical, SAVE                            :: flagstart = .TRUE.
    Logical                                  :: newname
    Integer                                  :: i

    If (flagstart) Then
      flagstart = .FALSE.
      Do i = 1,(FIRST_UNIT - 1)
        files%unit(i)%name = 'BLANK'
        files%unit(i)%tag = 'BLANK'
      End Do
      files%nunits = FIRST_UNIT
      files%unit(files%nunits)%name = trim(filename)
      file_getunit = files%nunits 
      Return
    End If

    newname = .TRUE.
    Do unit = 1,files%nunits
      If (trim(filename) == file_getname(unit)) Then
        file_getunit = unit
        newname = .FALSE.
        Exit
      End If
    End Do

    If (newname) Then
      files%nunits = files%nunits + 1
      If (files%nunits > MAX_FILES) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Maximum no. of files exceeded "
        Stop
      End If

      If (files%nunits == 6) files%nunits = files%nunits + 1  !skip unit 6
      files%unit(files%nunits)%name = filename
      file_getunit = files%nunits
    End If
  End Function file_getunit

  !----------------------------------------------------------------------------
  ! Returns the filename and unit number of the file with the specific tag type
  !----------------------------------------------------------------------------
  Subroutine file_gettype(tag,filename,unit)
    Character(*), Intent(IN)                 :: tag
    Character(*), Intent(OUT)                :: filename
    Integer, Intent(OUT), Optional           :: unit

    Integer                                  :: i,unitno

    unitno = 0
    filename = ''

    Do i = 1,files%nunits
      If (files%unit(i)%tag == trim(tag)) Then
        unitno = i
        filename = files%unit(i)%name
      End If
    End Do

    If (Present(unit)) unit = unitno

  End Subroutine file_gettype

  !----------------------------------------------------------------------------
  ! Sets the tag for a specific file, initializes the file if necessary
  !----------------------------------------------------------------------------
  Subroutine file_settag(filename,tag)
    Character(*), Intent(IN)                 :: tag,filename
    Integer                                  :: i
    Logical                                  :: found

    found = .FALSE.
    Do i = 1,files%nunits
      If (files%unit(i)%name == trim(filename)) Then
        files%unit(i)%tag = tag
        found = .TRUE.
      End If
    End Do

    If (.NOT. found) Then
      i = file_getunit(filename)
      files%unit(i)%tag = tag
    Endif

  End Subroutine file_settag

  !----------------------------------------------------------------------------
  ! Given a filename checks whether it is opened, if yes returns the 
  ! unit number, If not then opens it and returns the new unit number
  ! Requires: infile -- file name to open
  !           fileopt -- optional 3-bit number expecting/specifying file status
  !           tag -- optional file tag
  !           success -- optional flag to be passed back reporting success
  !   
  ! fileopt (IJK) numbers are:
  !   100 -- Unformatted, Old       J = 0/1 = unformatted/formatted file
  !   101 -- Unformatted, New       K = 0/1 = old/new file, 2=unknown
  !   110 -- Formatted, Old
  !   111 -- Formatted, New  
  !   210 -- Formatted, Old, APPEND position
  !   nothing or 0 -- no assumptions about previous file status DANGEROUS!
  !----------------------------------------------------------------------------
  Integer Function file_openRegular(infile,fileopt,tag,success)
    Character(*), Intent(IN)           :: infile
    Integer, Intent(In), Optional      :: fileopt  ! 3-bit number
    Character(*), Intent(IN), Optional :: tag
    Logical, Intent(OUT), Optional :: success

    Integer                                  :: unitno, ftype, ios
    Character(len=lstrLen)                   :: error_message,filename

    ! default value
    If (Present(success)) success=.true.

    !** Check for file options
    If (Present(fileopt)) Then
      ftype = fileopt
    Else
      ftype = 0
    End If

    filename = infile
    If (filename(1:2) == './') Then
      filename = filename(3:)
    End If

    !** Check whether already open
    unitno = file_checkfileopen(Trim(filename))

    !** If not open, Open it
    ios = 0
    If (unitno < 0) Then
      unitno = file_getunit(Trim(filename))
      Select Case (ftype)

      Case (0)
        Open(file=Trim(filename), unit=unitno,IOSTAT=ios)
        error_message = 'unidentified error type, nothing assumed'
      Case (100)
        Open(file=Trim(filename), unit=unitno, form="UNFORMATTED", &
            status="OLD", IOSTAT=ios)
        error_message = 'file does not already exist as binary file'
      Case (101)
        Open(file=Trim(filename), unit=unitno, form="UNFORMATTED", &
            status="NEW", IOSTAT=ios)
        error_message = 'file may already exist as binary file'
      Case (102)
        Open(file=Trim(filename), unit=unitno, form="UNFORMATTED", &
            status="UNKNOWN", IOSTAT=ios)
        error_message = 'file may already exist as binary file'
      Case (110)
        Open(file=Trim(filename), unit=unitno, form="FORMATTED", &
            status="OLD", IOSTAT=ios)
        error_message = 'file does not already exist as text file'
      Case (111)
        Open(file=Trim(filename), unit=unitno, form="FORMATTED", &
            status="NEW", IOSTAT=ios) 
        error_message = 'file may already exist as text file'
      Case (210)
        Open(file=Trim(filename), unit=unitno, form="FORMATTED", &
            status="OLD", IOSTAT=ios, POSITION='APPEND')
        error_message = 'file does not already exist as text file'
      Case Default
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(a)') 'what is going on'
        Stop
      End Select

    Else 
      ! if file is already open , then make sure that that is 
      ! consistent With the option fileopt
      Select Case (ftype)

      Case (101,111)
        ios=-1
        error_message="file is already open, cant open again as NEW"
      End Select

    End If

    If (ios/=0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(4a)') 'Error While opening file: ','"',Trim(filename),'"'
      Write(0,'(a)') Trim(error_message)
      Write(0,'(a)') '- consider checking for tabs in control file'
      Write(0,'(a)') '- check the PATH and environment variables'
      If (Present(tag)) Write(*,*) &
          "This file was supposed to be tagged as :",Trim(tag)
      If (Present(success)) Then
        ! let the calling module decide what to do
        success=.False.
        file_openRegular = -1
        Return
      Else
        Stop
      Endif
    End If

    If (Present(tag)) Call file_settag(Trim(filename),tag) 
    file_openRegular = unitno

  End Function file_openRegular
  
  !----------------------------------------------------------------------------
  ! Opens a new scratch file for writing given the record length. It returns
  ! the unit number of the scratch file. Multiple scratch files may be
  ! openned by specifying the optName option.
  !----------------------------------------------------------------------------
  Integer Function file_openScratch(recLen,optName)
    Integer, Intent(In) :: recLen
    Character(len=strLen), Optional, Intent(In) :: optName
    Integer :: error
    Character(len=strLen) :: filename

    !** Make up a fake name to get a unit number
    filename = "scratch"
    If (Present(optName)) filename = optName

    !** Get a unit number
    file_openScratch = file_getunit(trim(filename))

    !** Open the scratch file for writing.
    Open(unit=file_openScratch,STATUS="SCRATCH",iostat=error,RECL=recLen, &
        access="DIRECT",form="unformatted")

    If (error /= 0) Then
      Write(0,'(a,i5,a,i3)') __FILE__,__LINE__,": Could not open file, ", &
          "file error = ",error
      Stop
    End If

  End Function file_openScratch

  !---------------------------------------------------------------------------
  ! Checks to see if the file "filename" is open.  If it is then the 
  ! associated unit number is returned else a negative number is returned.
  ! *** THIS IS A COPY OF "isfileopen" from utils.f90
  ! Requires:  filename -- name of the file
  !---------------------------------------------------------------------------
  Integer Function file_checkfileopen(filename)
    Character(*), Intent(In)       :: filename

    Integer      :: unitno
    Logical      :: connected

    Inquire(File=filename, Opened=connected, Number=unitno)
    If (connected) Then
      file_checkfileopen = unitno
    Else
      file_checkfileopen = -1
    End If

  End Function file_checkfileopen  
  
  !----------------------------------------------------------------------
  ! Closes a file(filename), if its already open, Returns the old unitno
  ! Requires:  filename -- name of the file
  !----------------------------------------------------------------------
  Integer Function file_close(filename)
    Character(*), Intent(In)     :: filename

    Integer       :: unitno

    unitno = file_checkfileopen(filename)
    file_close = unitno

    !** Close the file if the unitno is reasonable
    If (unitno > 0) Close(unitno)

  End Function file_close

  !------------------------------------------------------------------------
  ! Appends the contents of one file to another file
  ! Requires:  base -- name of the base file to append to
  !            add -- name of the file whose contents to append
  !            kill -- if True, open the base file so it erases old stuff
  !------------------------------------------------------------------------
  Subroutine file_append(base,add,kill)
    Character(*), Intent(In)         :: base,add
    Logical, Intent(In), Optional    :: kill

    Integer            :: baseunit,addunit,ios
    Logical            :: erase
    Character(len=255) :: line

    !** Close the base file and determine how to reopen
    baseunit = file_close(base)
    erase = .False.
    If (Present(kill)) Then
      If (kill) erase = .True.
    End If

    !** Reopen base file, with file positioned to erase or append
    If (erase) Then
      baseunit = file_open(base)
    Else
      baseunit = file_open(base,210)
    End If

    !** Make sure the file to be added is open and rewound
    addunit = file_open(add,110)
    Rewind(unit=addunit)

    !** Read the add file line by line and write lines to base file
    Do
      Read(addunit,'(a)',IOSTAT=ios) line
      If (ios /= 0) Exit
      Write(baseunit,'(a)') Trim(line)
    End Do

    !** Close the file that was added
    addunit = file_close(add)

  End Subroutine file_append

End Module file
