!-----------------------------------------------------------------------------
! Takes isotherms and calculates the sigma
! it divides the isothemr into two regions low pressure and high pressure 
! values. The simpson integartion is done with higher accuracy at lower  
! pressures
!-----------------------------------------------------------------------------

Module sigmacalc

  Use defaults, Only: strLen, lstrLen, RDbl, zero
  Use isotherm, Only: Isotherm_Info, isotherm_init, isotherm_copyheader, &
      isotherm_write, isotherm_reallocate
  Use utils, Only   : filesrchstr, stripcmnt, toupper, toreal, allocErrDisplay
  Use file, Only : file_open
  Implicit None
  Save

  Private
  Public :: sigmacalc_init, sigmacalc_write

  Type Sigmacalc_Info
    Character(len=strLen)        :: inputfile, outputfile , sorbname
    Real(kind=RDbl)              :: min_P, max_P, mid_P
    Integer                      :: ndat
    Type(Isotherm_Info)          :: isoT
  End Type Sigmacalc_Info
  
  Type(SigmaCalc_Info), Dimension(:), Pointer :: sigmaPars

  Character(len=strLen), PARAMETER :: sigmacalc_section_tag = &
      " Sigma Calculation Section"

Contains

  !----------------------------------------------------------------------------
  ! Reads the ctrlfile and initializes
  !----------------------------------------------------------------------------
  Subroutine sigmacalc_init(ctrlfile)
    Character(*), Intent(In)                       :: ctrlfile
    Integer                :: lineno, unitno, error, i ,j, nisotherms, max_j
    Character(len=strlen) :: line, tag, midPtext
    Real(kind=RDbl) :: tempReal ,max_n_p, min_n_p, n_p, dx

    logical :: midPisAuto

    !** Open ctrlfile
    unitno=file_open(ctrlfile)
    Write(*,'(1x,2a,i4)') "reading from ctrlfile "
    !** Find the sigma section
    tag=sigmacalc_section_tag
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", Trim(tag), " in the control file"
      Stop
    Endif

    !** Read the variables
    Read(unitno,*) nisotherms


    Allocate(sigmaPars(nisotherms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    ! Read individual isotherms params
      Write(*,*) "reading Info fro each isotherm"
      
    Do i=1,Size(sigmaPars,1)
      Read(unitno,*)            !  a separator line

      !** read the raw info
      Read(unitno,*) sigmaPars(i)%sorbname
      Read(unitno,*) sigmaPars(i)%inputfile
      Read(unitno,*) sigmaPars(i)%outputfile
      Read(unitno,*) sigmaPars(i)%min_p, sigmaPars(i)%max_p
      Read(unitno,*) midPtext
      Read(unitno,*) sigmaPars(i)%ndat
      Write(*,*) "finsihed reading for isotherm", i
      !** process the raw text
      sigmaPars(i)%sorbname=Adjustl(Trim(stripcmnt(sigmaPars(i)%sorbname)))
      sigmaPars(i)%inputfile=Trim(stripcmnt(sigmaPars(i)%inputfile))
      sigmaPars(i)%outputfile=Trim(stripcmnt(sigmaPars(i)%outputfile))
      !** make sure ndiv is odd (for simpson integration)
      If (Mod(sigmaPars(i)%ndat,2)==0) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Increasing number of divisions to", sigmaPars(i)%ndat+1
        sigmaPars(i)%ndat=sigmaPars(i)%ndat+1
      Endif
      midPtext=stripcmnt(midPtext)
      midPisAuto=.False.
      If (Trim(toupper(midPtext))=="UNIFORM") Then
        !** we want to treat the whole pressure range uniformally
        sigmaPars(i)%mid_P=  (sigmaPars(i)%min_P+sigmaPars(i)%max_P)/2
      Else if (Trim(toupper(midPtext))=="AUTO") Then
        !** we want the program to estimate the low and high pressure regions 
        midPisAuto=.True.
      Else
        !** its given by the user
        tempReal=toreal(midPText, error)
        If (error==0) Then
          sigmaPars(i)%mid_P=  tempReal
        Else
          Write(*,*) "Toreal error, string used is : ", Trim(midPText)
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
      Endif

      !** Now we have  the reqd info for this isotherm open it read it 
      !** and do some initializatons
      Call isotherm_init(sigmaPars(i)%isoT, sigmaPars(i)%inputfile)

      ! check whether the molecule names match
      If (Trim(sigmaPars(i)%isoT%sorbname)/=Trim(sigmaPars(i)%sorbname)) &
          Then
        Write(*,*) "The sorb names do not match in inputfile and ctrlfile"
        Write(*,*) "name in inputfile : ",Trim(sigmaPars(i)%isoT%sorbname)
        Write(*,*) "name in ctlfile   : ",Trim(sigmaPars(i)%sorbname)
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      max_n_p= -1.00
      min_n_p= 1.00e+20
      max_j=-1
      !** If you dont understand this dont use auto option
      If (midPisAuto)  Then
        Do j=1,sigmaPars(i)%isoT%ndat
          n_p=sigmaPars(i)%isoT%N(j)/sigmaPars(i)%isoT%P(j)
          If (n_p<min_n_p) min_n_p=n_p
          If (n_p>max_n_p) Then 
            max_n_p=n_p
            max_j=j
          Endif
        End Do
        If (max_j>0) Then
          dx= sigmaPars(i)%isoT%P(sigmaPars(i)%isoT%ndat)-&
              sigmaPars(i)%isoT%P(max_j)
          dx=dx/2
          sigmaPars(i)%mid_P= sigmaPars(i)%isoT%P(max_j)+dx
          Write(*,*) "mid_P generated by AUTO option", sigmaPars(i)%mid_P
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        else
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif

      Endif

    End Do ! ** end of isotherms loop

  End Subroutine sigmacalc_init

  !------------------------------------------------------------------
  ! writes the inter/extra polated n, p, and sigma values
  !------------------------------------------------------------------
  Subroutine sigmacalc_write()
    Integer :: i,j,ndat, error
    Real(kind=RDbl) :: min_P, max_P, mid_P
    Real(kind=Rdbl) :: lowPdx, highPdx, sigma_p
    Real(kind=RDbl) :: n_p_first, n_p_second, n_p_third
    Real(kind=RDbl),Dimension(:),Pointer :: lowPList, lowSigmaList, lowNList
    Real(kind=RDbl),Dimension(:),Pointer :: highPList, highSigmaList, highNList

    Type(Isotherm_Info) :: temp_iso
    !** siome default memory
    ndat=10
    Allocate(lowPList(ndat), lowSigmaList(ndat), &
        lowNList(ndat), highPList(ndat), highSigmaList(ndat), &
        highNList(ndat), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)


    !** Loop over isotherms list
    Do i=1,Size(sigmaPars,1)
      ndat=sigmaPars(i)%ndat
      !** First we need to interpolate and extra polate the given data

      ! check memory
      Deallocate(lowPList, lowSigmaList, lowNList, highPList, &
          highSigmaList, highNList, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Allocate(lowPList(ndat), lowSigmaList(ndat), lowNList(ndat), &
          highPList(ndat), highSigmaList(ndat), highNList(ndat), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      ! reinitialize for each isotherm

      min_P=sigmaPars(i)%min_P  ! should be zero
      max_P=sigmaPars(i)%max_P
      mid_P=sigmaPars(i)%mid_P

      If ((min_P<0).Or.(min_P>mid_P).Or.(min_P>max_P)) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      If (mid_P>max_P) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      !** generate list of pressures at which inter/extra polations is neede
      lowPdx=(mid_P-min_P)/(ndat-1)
      highPdx=(max_P-mid_P)/(ndat-1)

      Do j=1,ndat
        lowPlist(j)=min_P+(j-1)*lowPdx
        highPlist(j)=mid_P+(j-1)*highPdx
      End Do

      Call sigmacalc_interpolate(lowPlist(1:ndat), lowNList(1:ndat), &
          sigmaPars(i) )

      Call sigmacalc_interpolate(highPlist(1:ndat), highNList(1:ndat), &
          sigmaPars(i) )

      !** Now calculate the integral using simpsons integral and store 
      !** the sigma values

      sigma_p=zero
      ! trapezoidal rule for first value 0 + \integ n/p * p
      lowSigmaList(1)=sigma_p+lowNList(1)
      sigma_p=lowSigmaList(1)

      ! 1. first lowP
      Do j=1,(ndat-1)/2
        n_p_second=lowNList(2*j)/lowPList(2*j)
        n_p_third=lowNList(2*j+1)/lowPList(2*j+1)
        If (j==1) Then
          n_p_first=n_p_second
          else
        n_p_first=lowNList(2*j-1)/lowPList(2*j-1)
        endif
        lowSigmaList(2*j+1)=sigma_p+&
            (n_p_first+4*n_p_second+n_P_third)*lowPdx/3
        sigma_p=lowSigmaList(2*j+1)
        !use average for valus at intermediate point
        lowSigmaList(2*j)=(lowSigmaList(2*j+1)+lowSigmaList(2*j-1))/2
      End Do

      ! 1. then highP
      highSigmaList(1)=sigma_p
      Do j=1,(ndat-1)/2
        n_p_first=highNList(2*j-1)/highPList(2*j-1)
        n_p_second=highNList(2*j)/highPList(2*j)
        n_p_third=highNList(2*j+1)/highPList(2*j+1)
        highSigmaList(2*j+1)=sigma_p+&
            (n_p_first+4*n_p_second+n_P_third)*highPdx/3
        sigma_p=highSigmaList(2*j+1)
        !use average for valus at intermediate point
        highSigmaList(2*j)=(highSigmaList(2*j+1)+highSigmaList(2*j-1))/2
      End Do

      !** Now write this into a type-2 isotherm file
      Call isotherm_copyheader(sigmaPars(i)%isoT,temp_iso)

      Call isotherm_reallocate(temp_iso, 2*ndat-1)

      temp_iso%filetype=2
      temp_iso%ndat=2*ndat-1
      temp_iso%N(1:ndat)=lowNList(1:ndat)
      temp_iso%N(ndat+1:2*ndat-1)=highNList(2:ndat)
      temp_iso%P(1:ndat)=lowPList(1:ndat)
      temp_iso%P(ndat+1:2*ndat-1)=highPList(2:ndat)
      temp_iso%sigma(1:ndat)=lowSigmaList(1:ndat)
      temp_iso%sigma(ndat+1:2*ndat-1)=highSigmaList(2:ndat)

      Call isotherm_write(temp_iso, sigmaPars(i)%outputfile)

    End Do !** end of isotherms list

  End Subroutine sigmacalc_write

  !-------------------------------------------------------------
  ! Calculates the value of N at given values of P, using the given 
  ! isotherm as the base
  !-------------------------------------------------------------
  Subroutine sigmacalc_interpolate(plist, nlist, sigmapar)
    Real(kind=RDbl), Dimension(:), Intent(in) :: plist
    Real(kind=RDbl), Dimension(:), Intent(out) :: nlist
    Type(SigmaCalc_Info) , Intent(in) :: sigmapar

    Real(kind=RDbl) :: minpval, maxpval, min_n_p, max_n
    Real(kind=RDbl) :: high_n_p, low_n_p, high_p, low_p , slope , n_p
    Integer :: i, ndat, high, low

    ndat=Size(plist,1)

    minpval=sigmapar%isoT%P(1)
    maxpval=sigmapar%isoT%P(sigmapar%isoT%ndat)
    min_n_p=sigmapar%isoT%N(2)/sigmapar%isoT%P(2)         ! value at j=2
    max_n  =sigmapar%isoT%maxUCLoading                    ! value at max_p

    Do i=1,ndat
      If (plist(i)<minpval) Then
!!$        Write(*,*) "WARNING", plist(i), i, minpval
!!$        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        ! assume constant n/p and extrapolate
        nlist(i)=plist(i)*min_n_p
      Elseif(plist(i)>maxpval) Then
!!$        Write(*,*) "WARNING", plist(i), maxpval
!!$        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        ! assume constant n  and extrapolate
        nlist(i)=max_n 
      Else
        ! Assume log (n/P) vs log(P) is linear and interpolate
        ! first need to find the interval 
        high=sigmacalc_getinterval(plist(i), sigmapar%isoT%P)
        low=high-1
        If ((plist(i)>sigmapar%isoT%P(high)).or.(plist(i)<sigmapar%isoT%P(low))) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Interpolation check"
        Write(*,*) low, high
        Write(*,*) sigmapar%isoT%P(low), sigmapar%isoT%P(high), plist(i)
        Write(*,*)
        stop
        endif

        high_n_p= sigmapar%isoT%N(high) / sigmapar%isoT%P(high)
        high_p  = sigmapar%isoT%P(high)
        low_n_p = sigmapar%isoT%N(low) / sigmapar%isoT%P(low)
        low_p   = sigmapar%isoT%P(low)
        ! log interpolation
        ! log (y/y1) =log(x/x1) * [log(y2/y1) / log(x2/x1)]
        ! y = y1 * ( {x/x1}^[log(y2/y1)/log(x2/x1)] )
        slope=Log(high_n_p/low_n_p)/Log(high_p/low_p)
        n_p= low_n_p* ( (plist(i)/low_p)**(slope) )
        nlist(i)= n_p*plist(i)          ! returning n not n_p
        If (((high_n_p-n_p)*(low_n_p-n_p))>0) Then
          Write(*,*) n_p, low_n_p, high_n_p
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
      Endif
    End Do
  End Subroutine sigmacalc_interpolate



  !-----------------------------------------------
  ! given a value of P finds the index of the smallest pressure array member 
  ! With value higher  than P 
  ! Heavily hacked dont reuse anywhere else copy from bmap code instead
  !-----------------------------------------------
  Integer Function sigmacalc_getinterval(P, parray)
    Real(kind=RDbl), Dimension(:), Intent(in) :: parray
    Real(kind=RDbl), Intent(in) :: P
    
    Integer   :: low, high, mid

    Real(kind=RDbl) :: mid_P
    low = 1
    high = size(parray)

    !** Do a binary search to find the first index i so that parray(i) > P
    Do
!!$$$$      If (low >= high) Exit
      If (low > high) Exit
      mid = (low + high)/2
      mid_P = parray(mid) 

      If (P < mid_P) Then
        high = mid - 1
      Else If (P > mid_P ) Then
        low = mid + 1
        mid = low
      Else If (P == mid_P) Then
        mid=mid+1
        Exit
      Endif
    Enddo

    !** low is equal to high
    sigmacalc_getinterval=mid
    Return

  End Function sigmacalc_getinterval

End Module sigmacalc
