!-------------------------------------------------------------------------
! This module handles sampling from an arbitrary one-dimensional 
! discretized distribution.  Currently the supported distributions are:
!   Gaussian distribution:  1=sigma, 2=center, 3=half-width, 4=binsize
!
! All distributions are defined from 'lo' to 'hi' and have values (1..N)
! that are centered in the middle of the bin:
!   lo---1---|---2---|---3---|....|---N---hi
! 
! Important routines:
!   distrib_init -- initialize a distribution using idstring and params
!   distrib_pick -- pick from a distribution
!   distrib_prob -- get the probability of being within a given bin
!   distrib_dump -- dump a data file tabulating the important quantities
!
! Possibile Improvements:
! 1) Use analytical or numerical integration to get cumulative weight.
!    The current summation method is inaccurate, especially for large bins
!-------------------------------------------------------------------------

Module distrib

  Use defaults, Only: RDbl, strLen, lstrLen, Pi, zero
  Use utils, Only: toint, toupper, findstr, int2str, real2str, &
      allocErrDisplay, deallocErrDisplay
  Use file, Only: file_open
  Use random, Only: rranf

  Implicit None
  Save

  Private
  Public :: Distribution_Params, distrib_init, distrib_pick, distrib_prob, &
      distrib_dump, distrib_display, distrib_clean

  Character(len=strLen), Dimension(1) :: distrib_idstrings = &
      (/'GAUSSIAN'/)  

  Type Distribution_Params
    Character(len=strLen)           :: idstring
    Integer                         :: nbins,id
    Real(kind=RDbl)                 :: lo,hi,binsize,normfactor
    Real(kind=RDbl), Dimension(:), Pointer  :: dparams
    Real(kind=RDbl), Dimension(:), Pointer  :: fn,binctr,cum
  End Type Distribution_Params

Contains
  !---------------------------------------------------------------------------
  ! Initializes the distribution
  ! Requires:  params -- distribution parameters
  !            idstring -- string identifying distribution type
  !            dparams -- array of distribution parameters
  !---------------------------------------------------------------------------
  Subroutine distrib_init(params,idstring,dparams)
    Type(Distribution_Params), Intent(Out)      :: params
    Character(*), Intent(In)                    :: idstring
    Real(kind=RDbl), Dimension(:), Intent(In)   :: dparams

    Integer         :: bin,nparams,error
    Real(kind=RDbl) :: halfbin

    !** Find the numerical ID of the distribution
    params%idstring = idstring
    params%id = findstr(distrib_idstrings,idstring)

    !** Interpret the input distribution parameters
    Select Case(params%id)
    Case(1)  !** Gaussian: 1=sigma, 2=center, 3=half-width, 4=binsize
      nparams = 2
      params%lo = dparams(2) - dparams(3)
      params%hi = dparams(2) + dparams(3)

      params%nbins = Int((params%hi - params%lo)/dparams(4))
      params%binsize = (params%hi - params%lo)/params%nbins

    Case Default
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' Could not identify the specified function ',Trim(params%idstring)
      Stop
    End Select

    !** Store the parameters
    Allocate(params%dparams(nparams),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    params%dparams = dparams(1:nparams)

    !** Allocate memory for the distribution
    Allocate(params%fn(params%nbins),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    Allocate(params%binctr(params%nbins),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    Allocate(params%cum(params%nbins),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        

    !** Create the discrete distribution
    halfbin = params%binsize/2.0_RDbl
    params%normfactor = 0.0_RDbl
    Do bin = 1,params%nbins
      params%binctr(bin) = params%lo + halfbin + params%binsize*(bin - 1)
      params%fn(bin) = distrib_evalfn(params,params%binctr(bin))
      params%normfactor = params%normfactor + params%fn(bin)
      params%cum(bin) = params%normfactor
    End Do

    !** Normalize the cumulative distribution
    params%cum = params%cum/params%normfactor

  End Subroutine distrib_init

  !---------------------------------------------------------------------------
  ! Evaluate the specified function at a given value
  ! Requires:  params -- distribution parameters
  !            value -- input value to function
  !---------------------------------------------------------------------------
  Real(kind=RDbl) Function distrib_evalfn(params,value)
    Type(Distribution_Params), Intent(In)      :: params
    Real(kind=RDbl), Intent(In)                :: value

    Real(kind=RDbl)        :: prefactor,exponent

    Select Case(params%id)
    Case(1)  !** Gaussian distribution: 1=sigma, 2=center
      prefactor = 1.0_RDbl/(params%dparams(1)*Sqrt(2.0_RDbl*Pi))
      exponent = -1.0_RDbl*((value - params%dparams(2))**2)/ &
          (2.0_RDbl*params%dparams(1)**2)
      distrib_evalfn = prefactor*Exp(exponent)

    Case Default
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' Could not identify the specified function ',Trim(params%idstring)
      Stop
    End Select

  End Function distrib_evalfn

  !---------------------------------------------------------------------------
  ! Pick from a distribution.  This translates to picking a bin according
  ! to its cumulative weight, then returning the bin center value.
  ! Requires:  params -- distribution parameters
  !            prob -- optional returned probability of picking bin
  !---------------------------------------------------------------------------
  Real(kind=RDbl) Function distrib_pick(params,prob)
    Type(Distribution_Params), Intent(In)      :: params
    Real(kind=RDbl), Intent(Out), Optional     :: prob

    Integer               :: bin
    Real(kind=RDbl)       :: randno

    !** Get a uniformly distributed randum number
    randno = rranf()

    !** Match this to a bin using the cumulative distribution values
    Do bin = 1,params%nbins
      If (params%cum(bin) > randno) Exit
    End Do
    
    distrib_pick = params%binctr(bin)

    !** Also return the probability of picking this bin if requested
    If (Present(prob)) Then
      prob = params%fn(bin)/params%normfactor
    End If

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) bin,distrib_pick

  End Function distrib_pick

  !---------------------------------------------------------------------------
  ! Given a value, determine which bin it's in, then return the probability
  ! of being in that bin.
  ! Requires:  params -- distribution parameters
  !            value -- given input value 
  !---------------------------------------------------------------------------
  Real(kind=RDbl) Function distrib_prob(params,value)
    Type(Distribution_Params), Intent(In)      :: params
    Real(kind=RDbl), Intent(In)                :: value

    Integer               :: bin
    Real(kind=RDbl)       :: frac

    !** Returns with zero probability if outside range
    distrib_prob = zero
    If ((value <= params%lo).Or.(value >= params%hi)) Return

    !** fractional distance along binning axis
    frac = (value - params%lo)/(params%hi - params%lo)

    !** Convert to bin number
    bin = Int(frac*params%nbins) + 1

    distrib_prob = params%fn(bin)/params%normfactor

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) value,bin,distrib_prob

  End Function distrib_prob

  !---------------------------------------------------------------------------
  ! Display the distribution
  ! Requires:  params -- distribution parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine distrib_display(params,indent,unit)
    Type(Distribution_Params), Intent(In)    :: params
    Integer, Intent(In)                      :: indent,unit

    Integer                     :: i,nref
    Real(kind=RDbl)             :: midpt
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2,string3,string4

    blank = Repeat(' ',indent)

    Write(unit,'(3a)') blank,'Distribution type: ',Trim(params%idstring)
    string1 = real2str(params%lo,8)
    string2 = real2str(params%hi,8)
    Write(unit,'(1x,5a)') blank,'value range: ',Trim(string1), &
        ' -> ',Trim(string2)
    midpt = params%lo + (params%hi - params%lo)/2.0_RDbl
    string1 = real2str(midpt,8)
    Write(unit,'(1x,3a)') blank,'Center of range is ',Trim(string1)
    string1 = int2str(params%nbins)
    string2 = real2str(params%binsize,5)
    Write(unit,'(1x,3a,2x,a)') blank,'number of bins, size: ',Trim(string1), &
        Trim(string2)

    If (params%nbins > 10) Then
      Write(unit,'(1x,2a)') blank,'suppressing bin content output'
      Return
    End If

    Write(unit,'(1x,3a)') blank,'bin number, bin center, function, cumulative fn'
    Do i = 1,params%nbins
      string1 = int2str(i)
      string2 = real2str(params%binctr(i),8)
      string3 = real2str(params%fn(i),8)
      string4 = real2str(params%cum(i),8)
      Write(unit,'(1x,a,a4,2x,3(a,2x))') blank,Trim(string1),Trim(string2), &
          Trim(string3),Trim(string4)
    End Do

  End Subroutine distrib_display

  !---------------------------------------------------------------------------
  ! Dump the distribution to a file
  ! Requires:  params -- distribution parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine distrib_dump(params,filename)
    Type(Distribution_Params), Intent(In)    :: params
    Character(*)                             :: filename

    Integer                     :: i,unit
    Real(kind=RDbl)             :: midpt
    Character(len=strLen)       :: string1,string2,string3,string4

    !** Open the file
    unit = file_open(filename)

    Write(unit,'(2a)') '# Distribution type: ',Trim(params%idstring)
    string1 = real2str(params%lo,8)
    string2 = real2str(params%hi,8)
    Write(unit,'(4a)') '# value range: ',Trim(string1), &
        ' -> ',Trim(string2)
    midpt = params%lo + (params%hi - params%lo)/2.0_RDbl
    string1 = real2str(midpt,8)
    Write(unit,'(2a)') '# Center of range is ',Trim(string1)
    string1 = int2str(params%nbins)
    string2 = real2str(params%binsize,5)
    Write(unit,'(2a,2x,a)') '# number of bins, size: ',Trim(string1), &
        Trim(string2)

    Write(unit,'(3a)') '# bin number, bin center, function, cumulative fn'
    Do i = 1,params%nbins
      string1 = int2str(i)
      string2 = real2str(params%binctr(i),12)
      string3 = real2str(params%fn(i),12)
      string4 = real2str(params%cum(i),12)
      Write(unit,'(a,t10,3(a,2x))') Trim(string1),Trim(string2), &
          Trim(string3),Trim(string4)
    End Do

  End Subroutine distrib_dump

  !---------------------------------------------------------------------------
  ! Clean the distribution information
  ! Requires:  params -- distribution parameters
  !---------------------------------------------------------------------------
  Subroutine distrib_clean(params)
    Type(Distribution_Params), Intent(InOut)    :: params

    Integer           :: error

    Deallocate(params%dparams,STAT=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
    Deallocate(params%fn,STAT=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
    Deallocate(params%binctr,STAT=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
    Deallocate(params%cum,STAT=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        

  End Subroutine distrib_clean


End Module distrib
