!-----------------------------------------------------------------------------
! deals with an isotherm file
! an isotherm has 4 junk lines and 3 columns
! column-1 : partial pressure
! column-2  : loading molec.uc
! column-3 : Total pressure
!-------------------------------------------------------------------
Module isotherm

  Use defaults, Only: strLen, lstrLen, RDbl, zero, one
  Use file, Only : file_open
  Use random, only : rranf
  Use utils, Only   : stripcmnt, split, allocErrDisplay, toint, toreal

  Implicit None
  Save

  Private
  Public :: Isot, isotherm_init, isotherm_guess, isotherm_MSerror, &
      isotherm_fit

  Type Isot
    Character(len=strLen)        :: sorbname
    Character(len=lstrLen)       :: description 
    Real(kind=RDbl)              :: Temp
    Real(kind=RDbl)      :: startP, endP, startN, endN, maxN, minN, maxP
    Integer                      :: ndat
    Integer                      :: fileType 
    Real(kind=RDbl), Dimension(:), Pointer :: N, P
  End Type Isot


Contains

  !----------------------------------------------------------------------------
  ! Reads the isotherm and initializes
  !----------------------------------------------------------------------------
  Subroutine isotherm_init(iso, isoTfile)
    Type(Isot), Intent(inout) :: iso
    Character(len=strLen), Intent(In)                       :: isoTfile 
    Integer                :: unitno, error, i, nfields, npts
    Character(len=strlen) :: line, tag, str, fname
    Character(len=strlen),Dimension(strlen) :: fields 
    Real(kind=RDbl) :: partial, load, totalp
    !** Open ctrlfile
    fname=adjustl(trim(isoTfile))
    unitno=file_open(fname)

    Read(unitno,'(a,///)') str ! skip total of 4 lines

    ! count the number of points
    npts=0
    Do
      Read(unitno,*,IOSTAT=error) partial,load, totalp
      If (error/=0) Then
        Exit
      Else
        npts=npts+1
      Endif
    End Do
    If(npts<3) Then
      Write(*,*) "at least three points required"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Allocate(iso%N(npts),iso%P(npts),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    iso%ndat=npts

    Rewind(unitno)
    Read(unitno,'(a,///)') str ! skip total of 4 lines

    iso%minN=1.0e70_RDbl
    iso%maxN=-1.0e70_RDbl
    Do i=1,npts
      Read(unitno,*,IOSTAT=error) partial,load, totalp
      iso%N(i)=load
      iso%P(i)=totalp
      If (load>iso%maxN) Then
        iso%maxN=load
        iso%maxP=totalP
      Endif
      if (load<iso%minN) iso%minN=load
      If ((i>2).And.totalp<iso%P(i-1) ) Then
        Write(*,*) "the pressure should be in increasing order"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        stop
      endif
    End Do
    iso%startP=iso%P(1)
    iso%endP=iso%P(npts)
    iso%startN=iso%N(1)
    iso%endN=iso%N(npts)
    Write(*,*) "Finished reading this isotherm file"
    Write(*,*) "Total number of points : " ,npts
    Write(*,*) "Maximum Load:", iso%maxN,"Minimum Load:", iso%minN
    Write(*,*)
  End Subroutine isotherm_init

  !----------------------------------------------------------------------------
  ! Reads the isotherm and initializes
  !----------------------------------------------------------------------------
  Subroutine isotherm_guess(iso, K, Lmax, a , Pmax)
    Type(Isot), Intent(in) :: iso
    Real(kind=RDbl) , Intent(out) :: K, Lmax, a , Pmax

    Integer :: nlinear
    ! assume first 20% of points are in linear
    nlinear=Int(iso%ndat*0.2)
    If(nlinear<3) nlinear=3
    K= (iso%N(nlinear)-iso%N(1)) / (iso%P(nlinear)-iso%P(1))
    Lmax=iso%maxN
    Pmax=iso%maxP
    a=0.2*K
  End Subroutine isotherm_guess

  !----------------------------------------------------------------------------
  ! adjusts K, Lmax, and a to minimize MSE
  !----------------------------------------------------------------------------
  Subroutine isotherm_fit(iso, K, Lmax, a , Pmax, tol)
    Type(Isot), Intent(in) :: iso
    Real(kind=RDbl) , Intent(inout) :: K, Lmax, a 
    Real(kind=RDbl) , Intent(in) :: Pmax, tol

    Integer :: t , i, uno
    Integer, Parameter :: MAX_ITER=3000, KTRIALS=200, LTRIALS=200, ATRIALS=200
    Real(kind=RDbl), Parameter:: REDUCFAC=0.9995
    Real(kind=RDbl) :: newK, newLmax, newa, mserror, minerror, max, min, dX 
    Real(kind=RDbl) :: kwidth, lwidth, awidth
    mserror=isotherm_MSerror(iso, newK, newLmax, newa ,Pmax)
    If (mserror<tol) Return

    newK=K
    newLmax=Lmax
    newa=a
    minerror=1.0e100_RDbl

    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    a=zero
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__

    max=K*2
    min=K/2
    kwidth=max-min
    max=Lmax*2
    min=Lmax/2
    lwidth=max-min
    max=a*2
    min=-a*2
    awidth=max-min

    Do i=1,MAX_ITER
      Write (*,*) i
      kwidth=kwidth*REDUCFAC
      lwidth=lwidth*REDUCFAC
      awidth=awidth*REDUCFAC

      newa=a
      Do t=0,KTRIALS
        newk=K+kwidth*(rranf()-0.5)
        mserror=isotherm_MSerror(iso, newK, newLmax, newa , Pmax)
        If (mserror<minerror) Then
          K=newK
          minerror=mserror
        endif
      Enddo

      newK=K
      Do t=0,LTRIALS
        newLmax=Lmax+lwidth*(rranf()-0.5)
        mserror=isotherm_MSerror(iso, newK, newLmax, newa , Pmax)
        If (mserror<minerror) Then
          Lmax=newLmax
          minerror=mserror
        Endif
      Enddo

      newLmax=Lmax
      Do t=0,ATRIALS
        newa=a+awidth*(rranf()-0.5)
        mserror=isotherm_MSerror(iso, newK, newLmax, newa , Pmax)
        If (mserror<minerror) Then
          a=newa
          minerror=mserror
        endif
      Enddo
      If (minerror<tol) Exit
    End Do

    Write(*,*) "Found minimum in ",i, "trials. ", "Error :", minerror
    Write(*,*) "parametrs :", K, Lmax, a

    uno=file_open("test.dat")

    Do i=1,iso%ndat
      Write(uno,*) iso%P(i), iso%N(i), isotherm_function(iso%P(i), K, Lmax, a/Pmax/Pmax)
    End Do
    close(uno)

  End Subroutine isotherm_fit

  !----------------------------------------------------------------------------
  ! Calculates the Mean Square error
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function isotherm_MSerror(iso, K, Lmax, a , Pmax)
    Type(Isot), Intent(in) :: iso
    Real(kind=RDbl) , Intent(in) :: K, Lmax, a , Pmax

    Real(kind=RDbl) :: error,p,l    
    Integer :: i
    error=zero
    Do i=1,iso%ndat
      p=iso%P(i)
      ! take care of zero loading case also
      If (Lmax<1.0e-20) Then
        l=zero
      Else
        l=isotherm_function(p, K, Lmax, (a/Pmax/Pmax))
      Endif
      error=error+( (l-iso%N(i)) * (l-iso%N(i)) )
    End Do
    error=error/iso%ndat
    isotherm_MSerror=error
  End Function isotherm_MSerror


  !----------------------------------------------------------------------------
  ! This is the fitting function
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function isotherm_function(p,K, Lmax, c)
    Real(kind=RDbl) , Intent(in) :: p, K, Lmax, c
    isotherm_function=(K*p) / ( one+ ((K/Lmax)*p) + c*p*p )
  End Function isotherm_function

!!$  !-----------------------------------------------------------------------
!!$  ! (Re)Allocates memory for isothewrm 
!!$  !-----------------------------------------------------------------------
!!$  Subroutine isotherm_write(isoTobj, filename)
!!$    Type(Isotherm_Info), Intent(inout) :: isoTobj
!!$    Character(*), Intent(In)                       :: filename
!!$    Integer :: unitno, i
!!$    unitno=file_open(filename)
!!$    Write(unitno,'(a)')       "%Molecule Name    : "//Trim(isoTobj%sorbname)
!!$    Write(unitno,'(a)')       "%Description      : "//Trim(isoTobj%description)
!!$    Write(unitno,'(a,f10.2)') "%Temperature      : ", isoTobj%Temp
!!$    Write(unitno,'(a,f10.2)') "%Max Loading/uc   : ", isoTobj%maxUCLoading
!!$    Write(unitno,'(a,i8)')    "%Number of Points : ", isoTobj%ndat
!!$    Write(unitno,'(a,i4)')    "%Type of file     : ", isoTobj%fileType
!!$    Write(unitno,'(a)')       "%  Pressure     N         sigma      "
!!$    Write(unitno,'(a)')       "%  kpa          mol/uc    mol/uc     "
!!$    If (isoTobj%filetype==1) Then
!!$      Do i=1,isoTobj%ndat
!!$        Write(unitno, *) isoTobj%P(i), isoTobj%N(i)
!!$      End Do
!!$    Elseif (isoTobj%filetype==2) Then
!!$Write(*,*)       isoTobj%ndat
!!$      Do i=1,isoTobj%ndat
!!$        Write(unitno, *) isoTobj%P(i), isoTobj%N(i), isoTobj%sigma(i)
!!$      End Do
!!$    Else
!!$      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$      Stop
!!$    Endif
!!$
!!$  End Subroutine isotherm_write
!!$


End Module isotherm
