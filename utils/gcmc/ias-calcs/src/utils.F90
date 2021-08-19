!-------------------------------------------------------------------------
! This module contains general utility subroutines and functions.
! Some of the general areas covered are:
! 1) File operations, such as 
! 2) String handling commands similar to some of those in Perl
! 3) Basic array operations such as inversions and searches
! 4) Routines for generating combinatorial quantities
! 5) Routines for determining machine specific precisions
!-------------------------------------------------------------------------

Module utils

  Use defaults, Only: strLen, lstrLen, RDbl, COMMENT_CHAR, twopi
#ifdef NAGCOMPILER
  Use f90_unix, Only: getenv
#endif

  Implicit None
  Save

  Private
  Public :: allocErrDisplay, isfileopen, split, getpath, filesrchstr, &
      stripcmnt, toreal, toint, toupper, tolower, findint, findstr, &
      genfilename, getlineswpair, realswap, getinvangle, makestr, sumlogical, &
      combine, erfc, filetoStringArr, maxn, findgt, getMachPrec, isdigit, &
      filesrchwocomment, deallocErrDisplay, comb, arrnorm, &
      multarrvec, multvecvec, dispvec, isblank, readblank, int2str, real2str, &
      str2seq,fileSrchStrAll

  Interface minn
    Module Procedure mini
    Module Procedure minr
  End Interface

  Interface maxn
    Module Procedure maxi
    Module Procedure maxr
  End Interface

  Interface inarray
    Module Procedure utils_inarrayi
  End Interface

  Interface filesrchstr
    Module Procedure fileSrchString
    Module Procedure fileSrchStrings
  End Interface

Contains
!----------------------------------------------------------
! String Utilities
!----------------------------------------------------------
  !-----------------------------------
  ! Converts the string "s1" to upper case
  !-----------------------------------
!!$  Character(len=strLen) Function toupper(s1)
  Function toupper(s1)
    Character(len=strLen) :: toupper
    Character(*), Intent(IN) :: s1
    Character(len=strLen)    :: temp
    Character           :: c
    Integer             :: length, LowerToUpper, i
 
    temp = s1
!    temp = Adjustl(s1)    ! Remove any leading blanks
    length = len(temp)
    LowerToUpper = iachar("A") - iachar("a") 

    toupper = ""
    Do i=1, length
       c = temp(i:i)
       If ( c >= "a" .AND. c <= "z") Then
          c = achar(iachar(c) + LowerToUpper)
       Endif
       toupper(i:i) = c
    Enddo
    Return
  End Function toupper

  !-----------------------------------
  ! Converts the string "s1" to lower case
  !-----------------------------------
!!$  Character(len=strLen) Function tolower(s1)
  Function tolower(s1)
    Character(len=strLen) :: tolower
    Character(*), Intent(IN) :: s1
    Character(len=strLen)    :: temp
    Character           :: c
    Integer             :: length, LowerToUpper, i
    
    temp = Adjustl(s1)    ! Remove any leading blanks
    length = len(temp)
    LowerToUpper = iachar("A") - iachar("a") 

    tolower = ""
    Do i=1, length
       c = temp(i:i)
       If ( c >= "A" .AND. c <= "Z") Then
          c = Achar(Iachar(c) - LowerToUpper)
       Endif
       tolower(i:i) = c
    Enddo
    Return
  End Function tolower


  !--------------------------------------------
  ! Returns the number of occurences of "delim"
  ! at the beginning of the string "line". 
  !--------------------------------------------
  Integer Function skip(line, delim)
    Character(*), Intent(in)  :: line
    Character, Intent(IN)     :: delim
    Integer    :: i, n_occur
    
    n_occur=0
    Do i=1,Len(line) 
      If (line(i:i) /= delim) Exit
     n_occur=i
    Enddo
    skip = n_occur
    Return
  End Function skip

  !-----------------------------------------------------
  ! Takes a string and returns the fields in an array
  ! where the fields are delimited using the optional 
  ! field "delim".  Default "delim" is space
  !-----------------------------------------------------
  Integer Function split(line, fields, delim)
    Character(*), Intent(in)                :: line
    Character(*), Dimension(:), Intent(out) :: fields
    Character, Intent(in), Optional         :: delim

    Character                 :: delimeter
    Character(len(line))      :: temp
    Integer     :: length, start, indx, nfields, maxfields, nskip

    ! make sure the fields are all blank
    fields(Lbound(fields,1):Ubound(fields,1)) = ""

    ! The default delimeter is a space
    If (Present(delim)) Then
      delimeter = delim
    Else
      delimeter = " "
    Endif

    ! Get the upper bound on number of fields
    maxfields = Ubound(fields,1)

    ! Remove any leading and trailing spaces
    temp = Trim(Adjustl(line))
    length = Len_trim(temp)
    If (length == 0) Then
      split = 1
      Return
    End If

    start = 1
    nfields = 0
    Do 
      ! Check if we've reached the end of the string
      If (start > length) Exit

      indx = Index(temp(start:length), delimeter)

      ! Make sure we are within bounds
      If (nfields == maxfields) Then
        Write(0,*) 'Size of the array is too small to hold all the fields'
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        split = nfields
        Return
      Endif
      nfields = nfields + 1

      If (indx /= 0) Then
        ! Get the field excluding the delimeter
        fields(nfields) = temp(start:start+indx-2)
        start = start + indx 
        If (start > length) Exit

        ! Skip over multiple delimeter occurrences
        nskip = skip(temp(start:length),delimeter)
        start = start + nskip
      Else
        fields(nfields) = temp(start:length)
        Exit
      End If
    End Do
    split = nfields
    Return
  End Function split

  !-----------------------------------------------------
  ! Takes a string and returns the fields in an array
  ! where the fields are delimited using the optional 
  ! field "delim".  Default "delim" is space
  !-----------------------------------------------------
!!$  Character(len=lstrLen) Function combine(fields, delim)
  Function combine(fields, delim)
    Character(len=lstrLen) :: combine
    Character, Intent(in), Optional         :: delim
    Character(*), Dimension(:), Intent(in) :: fields
    Character                 :: delimeter
    Integer     :: maxfields, i, leng
    Character(len=lstrLen) :: line

    ! make sure the line is blank
    !** Initialize the return line
    combine = ""
    line = ""

    !** The default delimeter is a space
    If (Present(delim)) Then
      delimeter = delim
    Else
      delimeter = " "
    Endif

    !** Get the upper bound on number of fields
    maxfields = Ubound(fields,1)

    !** Loop through the fields and add them to the line
    line = trim(fields(1))
    Do i = 2, maxfields
      leng = len_trim(line)
      line = line(1:leng)//delimeter//Trim(fields(i))
    End Do

    !** Trim any trailing spaces
    combine = Trim(line)

    Return
  End Function combine

  !------------------------------------------------------------------
  ! This function checks to see if a character string is "blank".  
  ! It will return a false value if there is anything but spaces, 
  ! tabs or '-' characters in the string.
  !------------------------------------------------------------------
  Logical Function isblank(line)
    Character(*), Intent(In)           :: line

    Integer                            :: i

    isblank = .True.

    Do i = 1,Len(Trim(line))
      If ((line(i:i) /= ' ').And.(line(i:i) /= '-').And. &
          (Ichar(line(i:i)) /= 9)) Then
        isblank = .False.
        Return
      End If
    End Do

  End Function isblank

  !------------------------------------------------------------------
  ! This subroutine reads a line that is expected to be blank.  It
  ! Returns an error message and stops if this is not the case.
  !------------------------------------------------------------------
  Subroutine readblank(unit,filename,lineno,errormsg)
    Integer, Intent(In)                :: unit
    Character(*), Intent(In)           :: filename
    Integer, Intent(In)                :: lineno
    Character(*), Intent(In), Optional :: errormsg

    Character(len=strLen) :: line

    Read(unit,*) line

    If (.Not.isblank(line)) Then
      Write(0,'(1x,2a,i4)') 'ERROR: expected to read a blank line at: ', &
          Trim(filename),lineno
      If (Present(errormsg)) Then
        Write(0,'(1x,a)') Trim(errormsg)
      End If
      Stop
    End If

  End Subroutine readblank

  !-------------------------------------------------------
  ! This function strips the "line" of the comment string
  ! which is defined by an optional argument "optcomment"
  ! The default is COMMENT_CHAR
  !-------------------------------------------------------
!!$  Character(len=lstrLen) Function stripcmnt(line, optcomment)
  Function stripcmnt(line, optcomment)
    Character(len=lstrLen)             :: stripcmnt
    Character(*), Intent(in)           :: line
    Character, Optional, Intent(in)    :: optcomment

    Character                             :: commentchar
    Character(len=lstrLen), Dimension(50) :: newline
    Integer                               :: nfields

    If (Present(optcomment)) Then
      Write(0,*) __FILE__,__LINE__,line
      Write(0,*) __FILE__,__LINE__,optcomment
    End If

    If (Present(optcomment)) Then
      commentchar = optcomment
    Else
      commentchar = COMMENT_CHAR
    End If

    newline = ""
    nfields = split(line, newline, commentchar)
    stripcmnt = newline(1)

  End Function stripcmnt

  !-----------------------------------------------------
  ! Takes an array of strings and converts the elements
  ! to real.  It returns the no. of fields converted
  !-----------------------------------------------------
  Real(kind=RDbl) Function toreal(strfield, opterror)
    Character(*), Intent(in) :: strfield
    Integer, Optional, Intent(out) :: opterror
    Integer                 :: error

    Read(strfield, *, IOSTAT=error) toreal
    If (Present(opterror)) Then
       opterror = error
    Else
       If (error /= 0) Then
          Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(0,'(3a,i4)') 'toreal error. String passed here is : " ',&
              Trim(strfield), '", error ',error 
          Stop
       End If
    End If
  End Function toreal

  !-----------------------------------------------------
  ! Takes an array of strings and converts the elements
  ! to integer.  It returns the no. of fields converted
  !-----------------------------------------------------
  Integer Function toint(strfield, opterror)
    Character(*), Intent(in) :: strfield
    Integer, Optional, Intent(out) :: opterror

    Integer         :: error
    Read(strfield, *, IOSTAT=error) toint
    If (Present(opterror)) Then
      opterror = error
    Else
       If (error /= 0) Then
         Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
             ' toint error'
         Write(0,*) "String used is ",trim(strfield)
          Stop
       End If
    End If
  End Function toint

  !-----------------------------------------------------
  ! Takes an integer and returns a string.  This is 
  ! useful because it can then be trimmed
  !-----------------------------------------------------
  Function int2str(int)
    Character(len=strLen)          :: int2str
    Integer, Intent(In)            :: int

    Write(int2str,'(i30)') int
    int2str = Adjustl(int2str)

  End Function int2str

  !-------------------------------------------------------------------------
  ! Takes a real and returns a string.  This is useful because it can then 
  ! be trimmed.  It also switches to scientific notation if practical.  If
  ! a suggested number of character is passed, it will either restrict the
  ! length to that number of characters or increase the number of characters
  ! as necessary to make the display reasonable.  
  ! Requires: number -- real number to convert to string
  !           sugnchar -- suggested number of characters in output, optional
  ! NOTE: Because it is a fairly complex routine, it should not be used near
  ! inner loops as it will slow down the code.  
  !-------------------------------------------------------------------------
  Function real2str(number,sugnchar)
    Character(len=strLen)          :: real2str
    Real(kind=RDbl), Intent(In)    :: number
    Integer, Intent(In), Optional  :: sugnchar

    Integer                        :: nchar,xtra,decimal_places,minchar
    Integer                        :: mindecimals
    Logical                        :: exp_form
    Character(len=10)              :: form
    Real(kind=RDbl)                :: logn,abslogn

    !** set defaults
    exp_form = .False.
    nchar = 12       !** default number of characters
    mindecimals = 1  !** minimum number of decimal places behind the point

    If (number == 0.0_RDbl) Then
      logn = 1
    Else
      logn = Log10(Abs(number))
    End If
    abslogn = Abs(logn)

    !** get the minimum required number of characters to avoid overflow
    xtra = 0
    If (number < 0) xtra = xtra + 1
    minchar = Max(0,Int(logn))
    minchar = minchar+xtra+1+mindecimals

    !** use the optional number of characters if present and reasonable
    If (Present(sugnchar)) Then
      nchar = Max(sugnchar,minchar)
      If (Abs(number) > 1) nchar = Max(sugnchar,minchar+1)
    End If

    !** decide if exponential notation is better
    decimal_places = Max((nchar-minchar),mindecimals)
    If ((logn+decimal_places) < 1) Then  !** switch to exp form
      nchar = Max(7,sugnchar)
      decimal_places = Max((nchar-minchar),mindecimals)
    End If

    If ((nchar > 6).And.((logn < -4.0_RDbl).Or.(logn > 6.0_RDbl))) Then
      exp_form = .True.
    End If

    !** set the format    
    If (exp_form) Then
      decimal_places = Max((nchar-6),mindecimals)
      Write(form,'(5a)') '(e',Trim(int2str(nchar)),'.', &
          Trim(int2str(decimal_places)),')'
    Else
      Write(form,'(5a)') '(f',Trim(int2str(nchar)),'.', &
          Trim(int2str(decimal_places)),')'
    End If
!LC    Write(*,*) nchar,'  ',Trim(form)

    Write(real2str,form) number

    real2str = Adjustl(real2str)

  End Function real2str

  !---------------------------------------------------------------------------
  ! Converts a string containing a sequence of integers into an array of 
  ! integers.  Format of the string is expected to be:  3,5-8,12 
  ! (for example).  This would be converted to: 3,5,6,7,8,12
  ! Requires: string -- string to process
  !           nums -- output array of integers
  !-----------------------------------------------------------------------------
  Integer Function str2seq(string,nums)
    Character(*), Intent(In)           :: string
    Integer, Dimension(:), Intent(Out) :: nums

    Integer                    :: i,j,n,n1,n2,nfields,nsubfields
    Character(len=strLen), Dimension(strLen) :: fields,subfields

    n = 0
    str2seq = -1

    !** split the string into comma separated chunks, then process
    nfields = split(string,fields,',')
    Do i = 1,nfields
      nsubfields = split(fields(i),subfields,'-')
      Select Case (nsubfields)
      Case(0) 
        Cycle
      Case(1)
        n = n+1
        If (n > Size(nums)) Then
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              " passed array too small"
          Stop
        End If
        nums(n) = toint(subfields(1))
      Case(2)
        n1 = toint(subfields(1))
        n2 = toint(subfields(2))
        Do j = n1,n2
          n = n + 1
          If (n > Size(nums)) Then
            Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
                " passed array too small"
            Stop
          End If
          nums(n) = j
        End Do
      Case Default
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            " Unprocessable string: ",Trim(string)
        Stop
      End Select
    End Do

    str2seq = n
    
  End Function str2seq

  !---------------------------------------------------
  ! Checks if a character is a digit
  !---------------------------------------------------
  Logical Function isdigit(char)
    Character, Intent(in) :: char
    
    If (char < '9' .And. char > '0') Then
      isdigit = .True.
    Else
      isdigit = .False.
    End If
  End Function isdigit

  !------------------------------------------------------
  ! Generates a compact string given two numbers "num1"
  ! and "num2" and a string in between.  Useful when
  ! you want to print out mathematical operations explicitly
  ! such as "num1/num2"
  !------------------------------------------------------
!!$  Character(len=strLen) Function makestr(num1, str, num2)
  Function makestr(num1, str, num2)
    Character(len=strLen) :: makestr
    Integer, Intent(in)   :: num1, num2
    Character(*), Intent(in)  :: str

    Character(len=strLen)     :: str1, str2
    
    Write(str1,'(i12)') num1
    Write(str2,'(i12)') num2
    makestr = Trim(Adjustl(str1))//Trim(str)//Trim(Adjustl(str2))
    Return
  End Function makestr
  
!----------------------------------------------------------
! File Utilities
!----------------------------------------------------------

  !-----------------------------------------------------------------
  ! Check if a file has "tabs" in it. The ascii value of "tab" is 9.
  ! The function returns the line number where the tab was found 
  ! otherwise it returns a zero.
  !-----------------------------------------------------------------
  Integer Function check_tabs(filename)
    Character(*), Intent(in) :: filename

    Character(len=lstrLen) :: inputline
    Integer                :: i, lineno, error
    
    check_tabs = 0
    lineno = 0

    Open(unit=17, file=filename)
    Do
      Read(17, '(a)', Iostat=error) inputline
      If (error /= 0) Exit
      
      lineno = lineno + 1
      Do i=1, Len(Trim(inputline))
        If (Ichar(inputline(i:i)) == 9) Then
          check_tabs = lineno
          Exit
        End If
      End Do
    End Do

    Close(unit=17)
  End Function check_tabs

  !-----------------------------------------------------
  ! This function skips "nlines" in the file connected to
  ! "unitno"
  !-----------------------------------------------------
  Subroutine skiplines(unitno, nlines)
    Integer, Intent(in) :: unitno, nlines

    Integer        :: i

    Rewind(unitno)
    Do i=1, nlines
      Read(unitno,*)
    End Do
  End Subroutine skiplines

  !----------------------------------------------------
  ! This function returns a string with the simulation
  ! no. "simno" appeded to the filename "filename"
  ! Requires: filename -- base string to add to
  !           simno -- integer to tack to the end
  !           secondno -- optional second number
  !----------------------------------------------------
!!$  Character(len=strLen) Function genfilename(filename, simno, secondno)
  Function genfilename(filename, simno, secondno)
    Character(len=strLen)          :: genfilename
    Character(*), Intent(In)       :: filename
    Integer, Intent(In)            :: simno
    Integer, Intent(In), Optional  :: secondno

    Character(len=strLen) :: i

    !** This is a hack to get it to work with the Lahey compiler
    i = int2str(simno)

    Write(genfilename, '(3a)') Trim(filename),'.',Trim(i)

    If (Present(secondno)) Then
      i = int2str(secondno)
      Write(genfilename, '(3a)') Trim(genfilename),'.',Trim(i)
    End If

  End Function genfilename

  !--------------------------------------------
  ! Get the pathname stored in the environment
  ! variable "s1"
  !--------------------------------------------
!!$  Character(len=lstrLen) Function getpath(s1)
  Function getpath(s1)
    Character(len=lstrLen) :: getpath
    Character(*), Intent(IN) :: s1
    Character(len=lstrLen)    :: path
    
    Call getenv(s1,path)
    If (path == ' ') Then
      Write(0,*)
      Write(0,'(3a)') 'Please set your ',s1,' environment variable'
      Write(0,'(2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    getpath = path
    If (path(len(path):len(path)) /= '/') Then
      getpath = trim(path) // '/'
    Endif
    Return
  End Function getpath

  !----------------------------------------------------------------------------
  ! Takes a string and a file unit number and returns the number of lines in
  ! the file that contain this string.  It also optionally returns the line
  ! numbers at which the strings occur.
  ! Requires: unitno -- unit number of an open file
  !           srchstr -- search string
  !           lines -- array of matching line numbers, optional
  !----------------------------------------------------------------------------
  Integer Function fileSrchStrAll(unitno,srchstr,lines)
    Integer, Intent(IN)                           :: unitno
    Character(*), Intent(IN)                      :: srchstr
    Integer, Dimension(:), Intent(OUT), Optional  :: lines

    Integer                         :: nmatch, ios, indx, length, line
    Character(len=4)                :: first_pass_str
    Character(len=255)              :: str

    length = Len(Trim(srchstr))
    If (length > 4) Then
      first_pass_str = srchstr(1:4)
    Else
      first_pass_str = srchstr
    End If

    !** rewind file
    Rewind(unitno)

    nmatch = 0 
    line = 0

    Do 
      Read(unitno,'(a)', IOSTAT=ios) str
!LC      Write(*,*) Trim(str)

      !** Check for end-of-file
      If (ios /= 0) Exit
      line = line + 1

      !** look for a shortened version first to speed (?) things up
      indx = Index(str,first_pass_str,.False.)
      If (indx /= 0) Then
        indx = Index(str,Trim(srchstr),.False.)
        If (indx /= 0) Then
          nmatch = nmatch + 1
          If (Present(lines)) lines(nmatch) = line
        End If
      End If

    End Do

    fileSrchStrAll = nmatch

  End Function fileSrchStrAll

  !---------------------------------------------------
  ! Takes a string and a file unit number and returns 
  ! the line in the associated text file with that 
  ! string.  The search has failed if return value is 0.
  ! If the string is found the return value is the line
  ! number where the string was found.
  !---------------------------------------------------
  Integer Function fileSrchString(unitno, srchstr, line, rewindfile)
    Integer, Intent(IN)             :: unitno
    Character(*), Intent(IN)        :: srchstr
    Character(*), Intent(OUT)       :: line
    Logical, Intent(IN), Optional   :: rewindfile
    Character(len=255)              :: str
    Integer                         :: ios, indx

    If (Present(rewindfile)) Then
      If (rewindfile) Rewind(unitno)
    End If

    fileSrchString = 0    
    Do 
      fileSrchString = fileSrchString + 1
      Read(unitno,'(a)', IOSTAT=ios) str

!LC      Write(*,*) Trim(str)
      !** Check for end-of-file
      If (ios /= 0) Then
        fileSrchString = 0
        line = ' '
        Exit
      End If

      indx = Index(toupper(str), Trim(toupper(srchstr)), .False.)
      If (indx /= 0) Then
        line = str
        Exit
      End If
    End Do
  End Function fileSrchString

  !---------------------------------------------------
  ! Same as fileSrchString, but first makes sure that 
  ! the line doesn't contain a comment
  !---------------------------------------------------
  Integer Function fileSrchwoComment(unitno, srchstr, line, rewindfile)
    Integer, Intent(IN)             :: unitno
    Character(*), Intent(IN)        :: srchstr
    Character(*), Intent(OUT)       :: line
    Logical, Intent(IN), Optional   :: rewindfile
    Character(len=255)              :: str
    Integer                         :: ios, indx

    If (Present(rewindfile)) Then
      If (rewindfile) Rewind(unitno)
    End If

    fileSrchwoComment = 0    
    Do 
      fileSrchwoComment = fileSrchwoComment + 1
      Read(unitno,'(a)', IOSTAT=ios) str
      str = stripcmnt(str)  !**the only difference compared to fileSrchString
      !** Check for end-of-file
      If (ios /= 0) Then
        fileSrchwoComment = 0
        line = ' '
        Exit
      End If

      indx = Index(toupper(str), Trim(toupper(srchstr)), .False.)
      If (indx /= 0) Then
        line = str
        Exit
      End If
    End Do
  End Function fileSrchwoComment

  !----------------------------------------------------------------------------
  ! Searches file unitno for a line that contains all the strings stored in
  ! the character array srchstrs. If found, it returns the line and line no.,
  ! otherwise, it returns 0. Optionally, it will rewind to the start of the 
  ! file if rewindfile is True.
  !----------------------------------------------------------------------------
  Integer Function fileSrchStrings(unitno, srchstrs, line, rewindfile)
    Integer, Intent(IN)             :: unitno
    Character(*), Intent(IN), Dimension(:) :: srchstrs
    Character(*), Intent(OUT)       :: line
    Logical, Intent(IN), Optional   :: rewindfile

    Character(len=255)              :: str
    Character(len=strLen), Dimension(strLen) :: fields
    Integer                         :: ios, nstr, i, j, k, nfields
    Logical :: found
    Integer, Dimension(Size(srchstrs,1)) :: foundpos

    nstr = Size(srchstrs,1)

    If (Present(rewindfile)) Then
      If (rewindfile) Rewind(unitno)
    End If

    fileSrchStrings = 0    
    Do 
      fileSrchStrings = fileSrchStrings + 1
      Read(unitno,'(a)', IOSTAT=ios) str
      ! Check for end-of-file
      If (ios /= 0) Then
        fileSrchStrings = 0
        line = ' '
        Exit
      End If

      nfields = split(str,fields)
      !** Check if all the strings exist
      Do i = 1, nstr
        found = .False.
        Do j = 1, nfields
          If (Trim(tolower(srchstrs(i))) == Trim(tolower(fields(j)))) Then
            ! We want a one-to-one mapping from elements of "srchstrs" to 
            ! element of "fields".  So, make sure that the field "j" does
            ! not already have a mapping
            found = .True.
            Do k=1, i-1
              If (j == foundpos(k)) Then
                found = .False.
                Exit
              End If
            End Do
            ! If found is still true then we have truly found a match
            ! and we can move on to the next field
            If (found) Then
              foundpos(i) = j
              Exit
            End If
          End If
        End Do
        If (.Not.found) Exit
      End Do
      If (found) Then 
        line = str
        Exit
      End If
    End Do
  End Function fileSrchStrings

  !---------------------------------------------------------------------------
  ! Checks to see if the file "filename" is open.  If it is then the 
  ! associated unit number is returned else a negative number is returned.
  !---------------------------------------------------------------------------
  Integer Function isfileopen(filename)
    Character(*), Intent(in)       :: filename
    Integer      :: unitno
    Logical      :: connected

    Inquire(File=trim(filename), Opened=connected, Number=unitno)
    If (connected) Then
      isfileopen = unitno
    Else
      isfileopen = -1
    End If
    Return
  End Function isfileopen

  !---------------------------------------------------
  ! Gets the no. of characters in a file
  !---------------------------------------------------
  !---------------------------------------------------
  ! Returns a array of integer with the index holding
  ! the no. of lines in the file and the second index
  ! the no. of characters in the file "unitno"
  !---------------------------------------------------
  Function getfilestats(unitno)
    Integer, Intent(in) :: unitno
    Integer, Dimension(2) :: getfilestats
    
    Character(150)   :: inputbuffer
    Integer          :: nlines, nchars, error
    
    Rewind(unitno)
    nchars = 0
    nlines = 0
    Do 
      nlines = nlines + 1
      Read(unitno,'(a)',IOSTAT=error) inputbuffer
      Write(*,'(i8,":",a)') nlines, Trim(inputbuffer)
      If (error /= 0) Exit
      nchars = nchars + Len(Trim(inputbuffer))
    End Do
    
    getfilestats = (/nlines, nchars/)
  End Function getfilestats

  !----------------------------------------------------------------------------
  ! Picks lines out of a file that have a certain pair of names in their
  ! first two words (order irrelevant).  File must already be open!
  !----------------------------------------------------------------------------
  Subroutine getlineswpair(unit,name1,name2,delim,nlines,lines)
    Integer, Intent(IN)                      :: unit
    Character(*), Intent(IN)                 :: name1,name2
    Character(len=1), Intent(IN)             :: delim
    Integer, Intent(OUT)                     :: nlines
    Character(*), Dimension(:), Intent(OUT)  :: lines

    Integer                                  :: lineno
    Integer                                  :: nfields
    Character(len=strLen), Dimension(10)     :: fields
    Integer                                  :: name1field,name2field
    Character(len=2*lstrLen)                 :: line

    nlines = 0
    Do 
      lineno = filesrchstr(unit,name1,line,.FALSE.)
      If (lineno == 0) Then
        Exit
!        Write(0,'(1x,2a,i4, 6a)') __FILE__," : ",__LINE__, &
!            " Could not find the pair ", Trim(name1)," ",Trim(name2), &
!            " in file ", Trim(filename)
      Else
        nfields = split(line,fields,delim)
        !** Check if the name is in field1 or field2
        If (Trim(fields(1)) == Trim(name1)) Then
          name1field = 1
          name2field = 2
        Else
          name1field = 2
          name2field = 1
        End If
        If (Trim(fields(name2field)) /= Trim(name2)) Cycle
      End If

      nlines = nlines + 1
      If (nlines > Size(lines)) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " passed string array too small"
        Stop
      End If
      lines(nlines) = line
      
    End Do
    Return
  End Subroutine getlineswpair

!---------------------------------------------------------------
! Array Utilities  
!---------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Searches an integer array to find the first instance of number. 
  ! Returns the position or 0 If not found. If the optional parameter 
  ! all is passed, Then the position of all occurances of number in 
  ! list are returned in all
  !----------------------------------------------------------------------------
  Integer Function findint(list,number,all)
    Integer, Dimension(:), Intent(In) :: list
    Integer, Intent(In) :: number
    Integer, Dimension(:), Intent(out), Optional :: all
    Integer :: i, j
    
    j = 0

    If (Size(list) < 1) Then
      findint = 0
      Return
    End If
    Do i = 1,Size(list)
      If (number == list(i)) Then
        findint = i
        If (Present(all)) Then
          j = j + 1
          all(j) = i
        Else
          Return
        End If
      End If
    End Do
    If (j == 0) findint = 0
  End Function findint

  !----------------------------------------------------------------------------
  ! Searches a STRING (len=strLen) array to find the first instance of 
  ! a given string.  Returns the position or 0 If not found. If the 
  ! optional parameter all is passed, then the position of all occurances 
  ! of this string in list are returned in all
  !----------------------------------------------------------------------------
  Integer Function findstr(list,srchstr,all)
    Character(len=strLen), Dimension(:), Intent(In) :: list
    Character(*), Intent(In)                        :: srchstr
    Integer, Dimension(:), Intent(Out), Optional    :: all

    Integer :: i, j
        
    j = 0

    If (Size(list) < 1) Then
      findstr = 0
      Return
    End If

    Do i = 1,Size(list)
!      Write(*,*) __FILE__,i,Trim(srchstr),Trim(list(i))
      If (Trim(srchstr) == Trim(list(i))) Then
        findstr = i
        If (Present(all)) Then
          j = j + 1
          all(j) = i
        Else
!          Write(*,*) 'found it',__FILE__,i,Trim(srchstr),Trim(list(i))
          Return
        End If
      End If
    End Do

    If (j == 0) findstr = 0

  End Function findstr

  !----------------------------------------------------------------------------
  ! Searches an integer array to find the first instance greater than number. 
  ! Returns the position or 0 If not found. If the optional parameter 
  ! all is passed, Then the position of all occurances of number in 
  ! list are returned in all and the number of instances found 
  !----------------------------------------------------------------------------
  Integer Function findgt(list,number,all)
    Integer, Dimension(:), Intent(In) :: list
    Integer, Intent(In) :: number
    Integer, Dimension(:), Intent(out), Optional :: all
    Integer :: i, j
    
    j = 0

    If (Size(list) < 1) Then
      findgt = 0
      Return
    End If
    Do i = 1,Size(list)
      If (number < list(i)) Then
        findgt = i
        If (Present(all)) Then
          j = j + 1
          all(j) = i
          findgt = j
        Else
          Return
        End If
      End If
    End Do
    If (j == 0) findgt = 0
  End Function findgt

  !---------------------------------------------
  ! find minimum in an integer or real array
  !---------------------------------------------
  Integer Function mini(list)
    Integer, Dimension(:) :: list
    Integer :: i

    mini = list(1)
    Do i = 2,Size(list)
      If (list(i) < mini) Then
        mini = list(i)
      End If
    End Do
  End Function mini

  Real(kind=RDbl) Function minr(list)
    Real(kind=RDbl), Dimension(:) :: list
    Integer :: i
    
    minr = list(1)
    Do i = 2,Size(list)
      If (list(i) < minr) Then
        minr = list(i)
      End If
    End Do
  End Function Minr

  !-----------------------------------------------
  ! finds the maximum in an array of ints or reals
  !-----------------------------------------------
  Integer Function maxi(list)
    Integer, Dimension(:) :: list
    Integer :: i

    maxi = list(1)
    Do i = 2,Size(list)
      If (list(i) > maxi) Then
        maxi = list(i)
      End If
    End Do
  End Function maxi

  Real(kind=RDbl) Function maxr(list)
    Real(kind=RDbl), Dimension(:) :: list
    Integer :: i
    
    maxr = list(1)
    Do i = 2,Size(list)
      If (list(i) > maxr) Then
        maxr = list(i)
      End If
    End Do
  End Function Maxr

  !----------------------------------------------
  ! Returns the index of the first instance of 
  ! n in the array, 0 if not found
  !----------------------------------------------
  Integer Function utils_inarrayi(n,list)

    Integer, Intent(In) :: n
    Integer, Dimension(:), Intent(In) :: list
    Integer :: i

    Do i = 1,Size(list)
      If (list(i) == n) Then
        utils_inarrayi = i
        Return
      End If
    End Do

    utils_inarrayi = 0

  End Function utils_inarrayi

  !-------------------------------------------------------------------
  ! This routine multiplies two arrays as matrices and returns the
  ! product in "prod"
  !-------------------------------------------------------------------
  Subroutine multarray(arr1, arr2, prod)
    Real(kind=RDbl), Dimension(:,:) :: arr1
    Real(kind=RDbl), Dimension(:,:) :: arr2
    Real(kind=RDbl), Dimension(:,:) :: prod

    Integer    :: dim1_1, dim1_2, dim2_1, dim2_2, dimp_1, dimp_2
    Integer    :: row, col, j
    Real(kind=RDbl) :: sum
    
    !** Make sure the dimensions of the array are consistent with matrix
    !** multiplication
    dim1_1 = Size(arr1, 1)   ! No. of rows
    dim1_2 = Size(arr1, 2)   ! No. of columns
    dim2_1 = Size(arr2, 1)
    dim2_2 = Size(arr2, 2)
    dimp_1 = Size(prod, 1)
    dimp_2 = Size(prod, 2)

    If (dim1_2 /= dim1_1) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

    If (dim1_1 /= dimp_1 .Or. dim2_2 /= dimp_2) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

    Do row=1, dim1_1
      Do col=1, dim2_2
        sum = 0.0_RDbl
        Do j=1, dim1_2
          sum = sum + arr1(row, j)*arr2(j, col)
        End Do
        prod(row, col) = sum
      End Do
    End Do
  End Subroutine multarray

  !--------------------------------------------------------------
  ! This routine multiplies an array("arr") and a vector("vec")
  ! and returns the result in the array "prod"
  !-------------------------------------------------------------
  Subroutine multarrvec(arr, vec, prod)
    Real(kind=RDbl), Dimension(:,:), Intent(in) :: arr
    Real(kind=RDbl), Dimension(:), Intent(in)   :: vec
    Real(kind=RDbl), Dimension(:), Intent(out)  :: prod

    Integer      :: i, j, arrdim1, arrdim2, vecdim, proddim
    Real(kind=RDbl)   :: vecel
    
    arrdim1 = Size(arr, 1)
    arrdim2 = Size(arr, 2)
    vecdim  = Size(vec)
    proddim = Size(prod)

    If (vecdim /= proddim .Or. arrdim2 /= vecdim) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Array dimensions don't match"
      Stop
    End If
    
    prod = 0.0_RDbl
    Do j = 1, arrdim2
      vecel = vec(j)
      Do i=1, arrdim1
        prod(i) = prod(i) + arr(i, j)*vecel
      End Do
    End Do
  End Subroutine multarrvec

  !---------------------------------------------------------------
  ! Multiplies two vectors "vec1" and "vec2" and returns the "prod"
  !---------------------------------------------------------------
  Subroutine multvecvec(vec1, vec2, prod)
    Real(kind=RDbl), Dimension(:), Intent(in) :: vec1, vec2
    Real(kind=RDbl), Intent(out)              :: prod

    Integer      :: i

    prod = 0.0_RDbl
    Do i=1, Size(vec1, 1)
      prod = prod + vec1(i)*vec2(i)
    End Do
  End Subroutine multvecvec

  !---------------------------------------------------------------
  ! Displays the contents of the array "arr" and outputs to the
  ! optional unitno "optunitno"
  !---------------------------------------------------------------
  Subroutine disparray(arr, label, fmt, optunitno)
    Real(kind=RDbl), Dimension(:, :), Intent(in) :: arr
    Character(*), Intent(in)        :: label
    Character(*), Intent(in)        :: fmt
    Integer, Optional, Intent(in)   :: optunitno

    Integer     :: i, j, unitno, dim1, dim2
    Character(len=strLen)     :: strfmt
    
    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6      ! Standard output
    End If

    strfmt = "("//Trim(fmt)//")"
    dim1 = Size(arr, 1)
    dim2 = Size(arr, 2)
    Write(unitno, '(a)') label
    Do i=1, dim1
      Do j=1, dim2
        Write(unitno, strfmt, Advance='No') arr(i,j)
      End Do
      Write(unitno, *)
    End Do
  End Subroutine disparray

  !---------------------------------------------------------------
  ! Displays the contents of the array "arr" and outputs to the
  ! optional unitno "optunitno"
  !---------------------------------------------------------------
  Subroutine dispvec(vec, fmt, optunitno)
    Real(kind=RDbl), Dimension(:), Intent(in) :: vec
    Character(*), Intent(in)        :: fmt
    Integer, Optional, Intent(in)   :: optunitno

    Integer     :: i, unitno, dim1
    Character(len=strLen)     :: strfmt
    
    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6      ! Standard output
    End If

    strfmt = "("//Trim(fmt)//")"
    dim1 = Size(vec, 1)
    Do i=1, dim1
      Write(unitno, strfmt, Advance='No') vec(i)
    End Do
  End Subroutine dispvec

  !-----------------------------------------------------------------
  ! This routine gets the norm of a vector "arr"
  !-----------------------------------------------------------------
  Real(kind=RDbl) Function arrnorm(arr)
    Real(kind=RDbl), Dimension(:), Intent(in) :: arr
    
    Integer       :: i
    
    arrnorm = 0.0_RDbl
    Do i=1, Size(arr)
      arrnorm = arrnorm + arr(i)*arr(i)
    End Do

    If (arrnorm > 0.0_RDbl) Then
      arrnorm = Sqrt(arrnorm)
    Else
      arrnorm = 0.0_Rdbl
    End If
    Return
  End Function arrnorm

!-----------------------------------------------------------
! Some other random stuff  
!-----------------------------------------------------------
  !-------------------------------------------------------------
  ! Calculates the factorial of a number
  !-------------------------------------------------------------
  Integer Function factorial(n)
    Integer, Intent(in) :: n
    
    Integer    :: i

    factorial = 1
    Do i=1, n
      factorial = factorial*i
    End Do
  End Function factorial

  !------------------------------------------------------------------
  ! Gets the no. of combinations of "n" objects taken "c" at a time
  !------------------------------------------------------------------
  Integer Function comb(n, c)
    Integer, Intent(in) :: n, c

    comb = factorial(n)/(factorial(n-c)*factorial(c))
  End Function comb

  !---------------------------------------------------------------
  ! This function returns the machine precision of the variable
  ! "eps".  This algorithm has been taken from "Numerical Methods"
  ! by Kahaner, Moler and Nash.
  !---------------------------------------------------------------
  Real(kind=RDbl) Function getMachPrec(eps)
    Real(kind=RDbl), Intent(inout) :: eps
    
    Real(kind=RDbl)  :: oldeps, eps1, lasteps
        
    !** Save the value of eps
    oldeps = eps

    !** Get the machine precision.
    eps = 1.0
    eps1 = eps + 1.0
    Do 
      If (eps1 <= 1.0) Exit
      lasteps = eps
      eps = eps/2.0
      eps1 = eps + 1.0
    Enddo

    eps = oldeps
    getMachPrec =  lasteps
    Return
  End Function getMachPrec

  !---------------------------------------------------------------
  ! This function returns the machine precision of the variable
  ! "eps".  This algorithm has been taken from "Numerical Methods"
  ! by Kahaner, Moler and Nash.
  !---------------------------------------------------------------
  Real(kind=RDbl) Function getMachRange(eps)
    Real(kind=RDbl), Intent(inout) :: eps
    
    Real(kind=RDbl)  :: oldeps, lasteps
        
    !** Save the value of eps
    oldeps = eps

    !** Get the machine precision.
    eps = 1.0
    Do 
      If (eps <= 0.0) Exit
      lasteps = eps
      eps = eps/2.0
    Enddo

    !** Restore the initial value
    eps = oldeps
    getMachRange =  lasteps
    Return
  End Function getMachRange

  !------------------------------------------------------------------
  ! This routine returns an angle "theta" in the correct quadrant given
  ! values of sin(theta), cos(theta)
  !------------------------------------------------------------------
  Subroutine getinvangle(sintheta, costheta, theta, optrange)
    Real(kind=RDbl), Intent(in) :: sintheta, costheta
    Real(kind=RDbl), Intent(out)  :: theta
    Logical, Optional, Intent(in) :: optrange

    Logical    :: zero_twopi
    
    If (Present(optrange)) Then
      zero_twopi = optrange
    Else
      zero_twopi = .True.
    End If

    theta = Acos(costheta)      ! theta is between 0 and pi
    If (sintheta < 0.0_RDbl) Then
      If (zero_twopi) Then
        theta = twopi - theta
      Else
        theta = -theta
      End If
    End If
  End Subroutine getinvangle

  !-----------------------------------------------------
  ! Swaps two (double) real numbers
  !-----------------------------------------------------

  Subroutine realswap(real1, real2)
  
    Real(Kind=RDbl), Intent(inout)  :: real1, real2
    
    Real(Kind=RDbl)                 :: temp

    temp = real1
    real1 = real2
    real2 = temp

  End Subroutine realswap


  !----------------------------------------------------------------------------
  ! Sums the logical array with .AND. and returns the result
  !----------------------------------------------------------------------------
  Logical Function sumlogical(values)
    Logical, Dimension(:), Intent(In) :: values
    Integer :: i

    sumlogical = values(1)
    Do i = 2, Size(values,1)
      sumlogical = (sumlogical .And. values(i))
    End Do
  End Function sumlogical
  
  !--------------------------------------------------------------
  ! displays pointer allocation errors, and stops the simulation
  !--------------------------------------------------------------
  Subroutine allocErrDisplay(filename, lineno,opt_variableName)
    Character(*) ,Intent(in) :: filename
    Integer,Intent(in)    :: lineno
    Character(*) ,Intent(in),Optional :: opt_variableName

    Character(len=strlen) :: varname
    
    If (Present(opt_variableName)) Then
      varname=opt_variableName
    Else
      varname=" "
    Endif
    
    Write(*,'(1x,a,i4)') "Error while allocating pointer variable " &
        //Trim(varname)//" in - "//Trim(filename)//" at line  no :",lineno

    Stop

  End Subroutine allocErrDisplay

  !----------------------------------------------------------------
  ! displays pointer deallocation errors, and stops the simulation
  !----------------------------------------------------------------
  Subroutine deallocErrDisplay(filename,lineno,opt_variableName)
    Character(*) ,Intent(in) :: filename
    Integer,Intent(in)    :: lineno
    Character(*) ,Intent(in),Optional :: opt_variableName

    Character(len=strlen) :: varname
    
    If (Present(opt_variableName)) Then
      varname=opt_variableName
    Else
      varname=" "
    Endif
    
    Write(*,'(1x,a,i4)') "Error while deallocating pointer variable " &
        //Trim(varname)//" in - "//Trim(filename)//" at line  no :",lineno

    Stop

  End Subroutine deallocErrDisplay
  
  !----------------------------------------------------------------------------
  !  This subroutine returns an approximation to the complementary error 
  !  Function
  !  Reference:  Abramowitz and Stegun, Handbook of Mathematical Functions,
  !              National Bureau of Standards, formula 7.1.26        
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function erfc(x)
    Real(Kind=RDbl), Intent(In) :: x
    Real(Kind=RDbl) :: t,xsq,tp
    Real(Kind=RDbl), Parameter :: a1 = 0.254829592, a2 = -0.284496736, &
        a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429, p  =  0.3275911

    t  = 1.0/(1.0 + p*x)
    xsq = x*x
    tp = t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))
    erfc = tp*Exp(-xsq)
    Return
  End Function erfc

  !---------------------------------------------------------
  !counts the no of lines in a text file( "unitno" is that of the text file )
  !Copies those lines into "stringArray"
  !this is useful while doing post-code things control files can be 
  ! copied easily and written to output files . 
  ! no of lines in the file is passed back if "nlines" is present
  !---------------------------------------------------------
  Subroutine FileToStringArr(unitno,stringArray,no_of_lines)
    Integer,Intent(in) :: unitno
    Character(len=2*strLen),Dimension(:),Pointer :: stringArray 
    Integer,Intent(out),Optional :: no_of_lines
    Integer :: i,error,nlines
    Character (len=2*strLen) :: inputbuffer

    !** rewind and count the no-of-lines from beginning
    Rewind(unitno)
    nlines = 0
    Do 
      Read(unitno,'(a)',IOSTAT=error) inputbuffer
      If (error /= 0) Then
        Exit
      Endif
      nlines = nlines + 1
    End Do
    !** read values into stringArray
    Rewind(unitno)
    Allocate(stringArray(nlines),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i=1,nlines
      Read(unitno,'(a)') inputbuffer
      stringArray(i)=inputbuffer
    End Do
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) stringArray
    If (present(no_of_lines)) no_of_lines=nlines
  End Subroutine FileToStringArr
  
End Module utils

