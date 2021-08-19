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
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Read(unitno,'(a,///)') str ! skip total of 4 lines
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    ! count the number of points
    npts=0
    Do
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
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
  Subroutine isotherm_guess(iso, K, Lmax, a , b, Pmax)
    Type(Isot), Intent(in) :: iso
    Real(kind=RDbl) , Intent(out) :: K, Lmax, a , b, Pmax

    Integer :: nlinear
    ! assume first 20% of points are in linear
    nlinear=Int(iso%ndat*0.2)
    If(nlinear<3) nlinear=3
    K= (iso%N(nlinear)-iso%N(1)) / (iso%P(nlinear)-iso%P(1))
    Lmax=iso%maxN
    Pmax=iso%maxP
    a=one
    b=zero
  End Subroutine isotherm_guess

  !----------------------------------------------------------------------------
  ! adjusts K, Lmax, and a to minimize MSE
  !----------------------------------------------------------------------------
  Subroutine isotherm_fit(iso, K, Lmax, a , b, Pmax, tol)
    Type(Isot), Intent(in) :: iso
    Real(kind=RDbl) , Intent(inout) :: K, Lmax, a, b 
    Real(kind=RDbl) , Intent(in) :: Pmax, tol

    Integer :: t , i, uno
    Integer, Parameter :: MAX_ITER=1000, KTRIALS=50, LTRIALS=50, ATRIALS=50, BTRIALS=50
    Real(kind=RDbl), Parameter:: REDUCFAC=0.975
    Real(kind=RDbl) :: newK, newLmax, newa, newb, mserror, minerror, max, min, dX 
    Real(kind=RDbl) :: kwidth, lwidth, awidth, bwidth
    mserror=isotherm_MSerror(iso, K, Lmax, a ,b,Pmax)
    If (mserror<tol) Return

    newK=K
    newLmax=Lmax
    newa=a
    newb=b
    minerror=1.0e100_RDbl

    max=K*2
    min=K/2
    kwidth=max-min
    max=Lmax*2
    min=Lmax/2
    lwidth=max-min
    max=a*2
    min=-a*2
    awidth=max-min
    bwidth=(1/iso%endP)
    Do i=1,MAX_ITER

      kwidth=kwidth*REDUCFAC
      lwidth=lwidth*REDUCFAC
      awidth=awidth*REDUCFAC
      bwidth=bwidth*REDUCFAC

      newb=b
      newa=a
      Do t=0,KTRIALS
        newk=K+kwidth*(rranf()-0.5)
        If (newK<zero) newK=-newK ! -ve values will create havoc
        mserror=isotherm_MSerror(iso, newK, newLmax, newa , newb, Pmax)
        If (mserror<minerror) Then
          K=newK
          minerror=mserror
        endif
      Enddo

      newK=K
      Do t=0,LTRIALS
        newLmax=Lmax+lwidth*(rranf()-0.5)
        If (newLmax<zero) newLmax=-newLmax ! -ve values will create havoc
        mserror=isotherm_MSerror(iso, newK, newLmax, newa , newb, Pmax)

        If (mserror<minerror) Then
          Lmax=newLmax
          minerror=mserror
        Endif
      Enddo

      newLmax=Lmax
      Do t=0,ATRIALS
        newa=a+awidth*(rranf()-0.5)
        mserror=isotherm_MSerror(iso, newK, newLmax, newa , newb, Pmax)
        If (mserror<minerror) Then
          a=newa
          minerror=mserror
        endif
      Enddo


      newa=a
      Do t=0,BTRIALS
        newb=b+bwidth*(rranf()-0.5)
        If ((1+newb*iso%endP)<zero) newb=zero
        mserror=isotherm_MSerror(iso, newK, newLmax, newa , newb, Pmax)
        If (mserror<minerror) Then
          b=newb
          minerror=mserror
        endif
      Enddo

      If (minerror<tol) Exit
    End Do

    Write(*,*) "Found minimum in ",i, "trials. ", "Error :", minerror
    Write(*,*) "parametrs :", K, Lmax, a, b

    uno=file_open("test.dat")

    Do i=1,iso%ndat
      Write(uno,*) iso%P(i), iso%N(i), isotherm_function(iso%P(i), K, Lmax, a, b)
    End Do
    close(uno)

  End Subroutine isotherm_fit

  !----------------------------------------------------------------------------
  ! Calculates the Mean Square error
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function isotherm_MSerror(iso, K, Lmax, a , b, Pmax)
    Type(Isot), Intent(in) :: iso
    Real(kind=RDbl) , Intent(in) :: K, Lmax, a , b, Pmax

    Real(kind=RDbl) :: error,p,l    
    Integer :: i
    error=zero
    Do i=1,iso%ndat
      p=iso%P(i)
      ! take care of zero loading case also
      If (Lmax<1.0e-20) Then
        l=zero
      Else
!        l=isotherm_function(p, K, Lmax, (a/Pmax/Pmax))
        l=isotherm_function(p, K, Lmax, a, b)
      Endif
      error=error+( (l-iso%N(i)) * (l-iso%N(i)) )
    End Do
    error=error/iso%ndat
    isotherm_MSerror=error
  End Function isotherm_MSerror

!!$
!!$  !----------------------------------------------------------------------------
!!$  ! This is the fitting function
!!$  !----------------------------------------------------------------------------
!!$  Real(kind=RDbl) Function isotherm_function(p,K, Lmax, c)
!!$    Real(kind=RDbl) , Intent(in) :: p, K, Lmax, c
!!$    isotherm_function=(K*p) / ( one+ ((K/Lmax)*p) + c*p*p )
!!$  End Function isotherm_function
!!$
#ifndef langmuir
  !----------------------------------------------------------------------------
  ! This is the fitting function Toth Eqn, Heuchel , snurr and buss eqn-12
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function isotherm_function(p,K, Lmax, c, b)
    Real(kind=RDbl) , Intent(in) :: p, K, Lmax, c,b
    isotherm_function=(K*p) / ( one+ (((K/Lmax)*p/(1+b*P))**c) )**(one/c)
  End Function isotherm_function
#endif

#ifdef langmuir
  !----------------------------------------------------------------------------
  ! This is the fitting function for langmuir simple fit 
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function isotherm_function(p,K, Lmax, c, b)
    Real(kind=RDbl) , Intent(in) :: p, K, Lmax, c,b
    isotherm_function=(K*p) / ( one+ K*p/Lmax )
  End Function isotherm_function
#endif

End Module isotherm
