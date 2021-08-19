!-----------------------------------------------------------------------------
! deals with an isotherm file
! contains a header and isotherm in one of the formats : 1) P, n  
! 2) P, n, sigma(surface pressure)
!
!Example : 
! ----------------------------------------------------------------
! %Molecule Name : Butane  # forcefields =used
! %Description  : silicalite gcmc
! %Temperature   : 300      # K
! %Max Loading/uc : 4.00     # molecules/unitcell
! %Number of Points : 2
! %Type of file     : 2
! % Pressure  N       sigma
! % kpa       /uc     /uc
!  1.00       3.85    0.87
!  2.00       4.00    0.93
!-------------------------------------------------------------------
!Note 1.)the "%" symbols are necessary
!     2.) The "#" symbols and the text following are optional they r comments
!     3.) Text between "%" and ":" should be same as above
!     4.) first 6 lines should adhere to above format, lines-7&8 are for 
!         units not used anywhere. data should be on lines 9 onwards
!     5.) Type-1 -> data has first two coulmns only
!         Type-2 -> data has all 3 columns 
!-----------------------------------------------------------------------------

Module isotherm

  Use defaults, Only: strLen, lstrLen, RDbl
  Use utils, Only   : stripcmnt, split, allocErrDisplay, toint, toreal
  Use file, Only : file_open
  Implicit None
  Save

  Private
  Public :: isotherm_init, isotherm_write, isotherm_reallocate, &
      isotherm_copyheader, Isotherm_Info

  Type Isotherm_Info
    Character(len=strLen)        :: sorbname
    Character(len=lstrLen)       :: description 
    Real(kind=RDbl)              :: Temp
    Real(kind=RDbl)              :: maxUCLoading
    Integer                      :: ndat
    Integer                      :: fileType 
    Real(kind=RDbl), Dimension(:), Pointer :: N, P, sigma
  End Type Isotherm_Info


Contains

  !----------------------------------------------------------------------------
  ! Reads the isotherm and initializes
  !----------------------------------------------------------------------------
  Subroutine isotherm_init(isoTobj, isoTfile)
    Type(Isotherm_Info), Intent(inout) :: isoTobj
    Character(len=strLen), Intent(In)                       :: isoTfile 
    Integer                :: unitno, error, i, nfields
    Character(len=strlen) :: line, tag, str, fname
    Character(len=strlen),Dimension(strlen) :: fields 

    !** Open ctrlfile
    fname=adjustl(trim(isoTfile))
    unitno=file_open(fname)

    Read(unitno,'(a)') str
    str=stripcmnt(str)
    nfields=split(str,fields,":")
    If (Trim(Adjustl(fields(1)))=="%Molecule Name") Then
      isoTobj%sorbname=Adjustl(Trim(fields(2)))
    Else
      Write(*,*) "check header "
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Read(unitno,'(a)') str
    str=stripcmnt(str)
    nfields=split(str,fields,":")
    If (Trim(Adjustl(fields(1)))=="%Description") Then
      isoTobj%description=Trim(fields(2))
    Else
      Write(*,*) "check header "
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Read(unitno,'(a)') str
    str=stripcmnt(str)
    nfields=split(str,fields,":")
    If (Trim(Adjustl(fields(1)))=="%Temperature") Then
      isoTobj%Temp=toreal(Trim(fields(2)))
    Else
      Write(*,*) "check header "
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Read(unitno,'(a)') str
    str=stripcmnt(str)
    nfields=split(str,fields,":")
    If (Trim(Adjustl(fields(1)))=="%Max Loading/uc") Then
      isoTobj%maxUCLoading=toreal(Trim(fields(2)))
    Else
      Write(*,*) "check header "
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Read(unitno,'(a)') str
    str=stripcmnt(str)
    nfields=split(str,fields,":")
    If (Trim(Adjustl(fields(1)))=="%Number of Points") Then
      isoTobj%ndat=toint(Trim(fields(2)))
    Else
      Write(*,*) "check header "
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Read(unitno,'(a)') str
    str=stripcmnt(str)
    nfields=split(str,fields,":")
    If (Trim(Adjustl(fields(1)))=="%Type of file") Then
      isoTobj%fileType=toint(Trim(fields(2)))
      If ((isoTobj%fileType>2).Or.(isoTobj%fileType<1)) Then
        Write(*,*) "filetype should be 1 or 2"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        stop
      Endif
    Else
      Write(*,*) "check header "
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    !** read the two blank lines
    Read(unitno,'(a)') str
    Read(unitno,'(a)') str

    Allocate(isoTobj%N(isoTobj%ndat),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(isoTobj%P(isoTobj%ndat),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(isoTobj%sigma(isoTobj%ndat),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Select Case (isoTobj%filetype)
    Case(1)
      Do i=1,isoTobj%ndat
        Read(unitno,*) isoTobj%P(i), isoTobj%N(i)
      End Do
    Case(2)
      Do i=1,isoTobj%ndat
        Read(unitno,*) isoTobj%P(i), isoTobj%N(i), isoTobj%sigma(i)
      End Do
    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select
    Close(unitno)
  End Subroutine isotherm_init


  !-----------------------------------------------------------------------
  ! Copies the header of the isotherm, Dont copy ndat,filetype
  !-----------------------------------------------------------------------
  Subroutine isotherm_copyheader(fromIso , toIso)
    Type(Isotherm_Info), Intent(in) :: fromIso
    Type(Isotherm_Info), Intent(out) :: toIso
    toIso%sorbname=fromIso%sorbName
    toIso%description=fromIso%description
    toIso%Temp=fromIso%Temp
    toIso%maxUCloading=fromIso%maxUCloading
  End Subroutine isotherm_copyheader

  !-----------------------------------------------------------------------
  ! (Re)Allocates memory for isothewrm 
  !-----------------------------------------------------------------------
  Subroutine isotherm_reallocate(isoTobj, ndat)
    Type(Isotherm_Info), Intent(inout) :: isoTobj
    Integer , Intent(in) :: ndat
    Integer :: error
!!$
!!$    If (Size(isoTobj%N,1)>0) Then
!!$      Deallocate(isoTobj%N, isoTobj%P, isoTobj%sigma, STAT=error)
!!$      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
!!$    Endif
      Allocate(isoTobj%N(ndat), isoTobj%P(ndat), isoTobj%sigma(ndat), &
          STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End Subroutine isotherm_reallocate

  !-----------------------------------------------------------------------
  ! (Re)Allocates memory for isothewrm 
  !-----------------------------------------------------------------------
  Subroutine isotherm_write(isoTobj, filename)
    Type(Isotherm_Info), Intent(inout) :: isoTobj
    Character(*), Intent(In)                       :: filename
    Integer :: unitno, i
    unitno=file_open(filename)
    Write(unitno,'(a)')       "%Molecule Name    : "//Trim(isoTobj%sorbname)
    Write(unitno,'(a)')       "%Description      : "//Trim(isoTobj%description)
    Write(unitno,'(a,f10.2)') "%Temperature      : ", isoTobj%Temp
    Write(unitno,'(a,f10.2)') "%Max Loading/uc   : ", isoTobj%maxUCLoading
    Write(unitno,'(a,i8)')    "%Number of Points : ", isoTobj%ndat
    Write(unitno,'(a,i4)')    "%Type of file     : ", isoTobj%fileType
    Write(unitno,'(a)')       "%  Pressure     N         sigma      "
    Write(unitno,'(a)')       "%  kpa          mol/uc    mol/uc     "
    If (isoTobj%filetype==1) Then
      Do i=1,isoTobj%ndat
        Write(unitno, *) isoTobj%P(i), isoTobj%N(i)
      End Do
    Elseif (isoTobj%filetype==2) Then
Write(*,*)       isoTobj%ndat
      Do i=1,isoTobj%ndat
        Write(unitno, *) isoTobj%P(i), isoTobj%N(i), isoTobj%sigma(i)
      End Do
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

  End Subroutine isotherm_write



End Module isotherm
