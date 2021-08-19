!-------------------------------------------------------------------
! This module provides the necessary functionality for generating a
! multihist.
!-------------------------------------------------------------------
Module multihist

  Use defaults,    Only : RDbl, one, zero, STATS_BLOCKSIZE, strlen
  Use stats,    Only : Statistics, stats_init, stats_update

  Implicit None
  Save

  Private 
  Public :: MultiHist_Params, multihist_init, multihist_update, &
      multihist_display, multihist_normalize, multihist_getnbins

  !** the indices corresp. to x,y and z dimensions
  Integer,Parameter :: HIST_XDIM=1, HIST_YDIM=2, HIST_ZDIM=3

  Type Multihist_Params
    Integer                        :: dim,nelements,totBins
    Real(Kind=RDbl)                :: totVolume
    Integer,Dimension(:),Pointer   :: nbins,temp_bin_no,multi_fac
    !** For ex: nbins(x)=20, nbins(y)=13, nbins(z)=9
    !** temp_bin_no(x)=3(0..19),temp_bin_no(y)=7(0..12),temp_bin_no(z)=8 (0..8)
    !** THEN "bin number"= 8*multifac(z) + 7*multi_fac(y) +3*multi_fac(x)
    !** multi_fac(z)= 1  ,(least power)
    !** multi_fac(y)= 9    = nbins(z)
    !** multi_fac(x)= 9*13 = nbins(y)*nbins(z)
    Real(kind=RDbl),Dimension(:),Pointer  :: lowvalue,highvalue
    Real(kind=RDbl),Dimension(:),Pointer  :: min,max,binsize
    Type(Statistics),Dimension(:),Pointer :: stat
    !** In the array of x,y,z values z varies fast, then y then x

    Real(kind=RDbl), Dimension(:), Pointer :: probability
    Integer , Dimension(:), Pointer :: histg
    Logical :: is_init=.false.
  End Type Multihist_Params

Contains
  !----------------------------------------------------------
  ! Initializes the multihist object "histobj" by supplying 
  ! to it the no. of bins "nbins", the low value "lowvalue" and the
  ! high value "highvalue"
  !----------------------------------------------------------
  Subroutine multihist_init(histobj, nbins, lowvalue, highvalue,dim)
    Type(Multihist_Params), Intent(out) :: histobj
    Integer, Dimension(:),Intent(in):: nbins
    Integer, Intent(in)      :: dim
    Real(kind=RDbl), Dimension(:),Intent(in) :: lowvalue,highvalue
    Real(kind=RDbl)                          :: fac

    Integer         :: error,totnbins,i

    If (histobj%is_init) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) " attempt to re-initialize here"
      Stop
    Endif      

    histobj%dim=dim
    histobj%nelements = 0
    Call multihist_AllocDim(histobj,dim)
    If (.not.((Size(highvalue,1)==dim).and.(Size(lowvalue,1)==dim).and. &
        (Size(nbins,1)==dim) )  ) Then
            Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Write(*,*) "Wrong dimensions here in multihistogram"
            Stop
            Endif
    histobj%totVolume=one
    Do i=1,histobj%dim
      !** Assumes rectangular array in multi dimensions
      histobj%totVolume=histobj%totVolume*(highvalue(i) - lowvalue(i))
      histobj%lowvalue(i)=lowvalue(i)
      histobj%highvalue(i)=highvalue(i)
      histobj%min(i)=highvalue(i)+one
      histobj%max(i)=lowvalue(i)-one
      ! STATS_BLOCKSIZE is the default blocksize
      Call stats_init(histobj%stat(i),'Hist. Object',STATS_BLOCKSIZE,.False.)
      histobj%nbins(i)=nbins(i)
      histobj%binsize(i)   = (highvalue(i) - lowvalue(i))/nbins(i)
    End Do
    
    fac=1
    Do i=1,histobj%dim
      histobj%multi_fac(histobj%dim-i+1)=fac
      fac=fac*histobj%nbins(histobj%dim-i+1)
      Write(*,*) fac
    End Do
    

    !** total number of bins
    histobj%totBins=fac
    totnbins=fac

    Allocate(histobj%histg(totnbins), STAT = error)
    If (error/=0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, "Error while & 
          & allocating multihist variable"
      Stop
    Endif

    Allocate(histobj%probability(totnbins), STAT = error)
    If (error/=0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, "Error while & 
          & allocating multihist variable"
      Stop
    Endif
    
    !set all elemnts to zero
    histobj%histg = 0
    histobj%probability = zero
    histobj%is_init=.true.
  End Subroutine multihist_init

  !----------------------------------------------------------
  ! Allocates all fields according to the dimensions
  !----------------------------------------------------------
  Subroutine multihist_AllocDim(histobj,dim)
    Type(Multihist_Params), Intent(inout) :: histobj
    Integer ,Intent(in):: dim
    Integer:: error

    Allocate(histobj%nbins(dim), STAT=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"nbins")

    Allocate(histobj%temp_bin_no(histobj%dim),stat=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"temp_bin_no")

    Allocate(histobj%multi_fac(histobj%dim),stat=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"multifac")

    Allocate(histobj%lowvalue(dim), STAT=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"lowvalue")

    Allocate(histobj%highvalue(dim), STAT=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"highvalue")

    Allocate(histobj%min(dim), STAT=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"min")

    Allocate(histobj%max(dim), STAT=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"max")

    Allocate(histobj%binsize(dim), STAT=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"binsize")

    Allocate(histobj%stat(dim), STAT=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"stat")

  End Subroutine multihist_AllocDim
  
  !----------------------------------------------------------
  ! Initializes the multihist object "histobj" by supplying 
  ! to it the low value "low", the  high value "high" and
  ! the size of each bin "binsize"
  !----------------------------------------------------------
!  Subroutine multihist_init2(histobj, lowvalue, highvalue, binsize)
!    Type(Multihist_Params), Intent(inout) :: histobj
!    Real(kind=RDbl), Intent(in) :: lowvalue
!    Real(kind=RDbl), Intent(in) :: highvalue
!    Real(kind=RDbl), Intent(in) :: binsize
!    
!    Integer         :: error
!    Integer         :: nbins
!    
!    histobj%nelements = 0
!    nbins = (highvalue - lowvalue)/binsize
!    histobj%nbins     = nbins
!    histobj%lowvalue    = lowvalue
!    histobj%highvalue   = highvalue
!    histobj%binsize   = binsize
!    Allocate(histobj%histg(nbins), STAT = error)
!    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
!
!    Allocate(histobj%probability(nbins), STAT = error)
!    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
!
!    !** sets the min and max to aribtrarily high and low values respectidely
!    histobj%min=highvalue+one
!    histobj%max=lowvalue-one
!
!    ! STATS_BLOCKSIZE is the default blocksize
!    Call stats_init(histobj%stat,'Hist. Object',STATS_BLOCKSIZE,.False.)
!    !set all elemnts to zero
!    histobj%histg = 0
!    histobj%probability = zero
!  End Subroutine multihist_init2

  !-----------------------------------------------------------
  ! Returns the size of each bin in "histobj" in the given axis
  !-----------------------------------------------------------
  Real(kind=RDbl) Function multihist_getbinsize(histobj,direction)
    Type(Multihist_Params), Intent(in) :: histobj
    Integer ,Intent(in):: direction 
    Call multihist_checkinit(histobj)
    multihist_getbinsize = histobj%binsize(direction)
    Return
  End Function multihist_getbinsize
  
  !-----------------------------------------------------------
  ! Returns the no. of bins in "histobj" in the given axis
  !-----------------------------------------------------------
  Integer Function multihist_getnbins(histobj,direction)
    Type(Multihist_Params), Intent(in) :: histobj
    Integer ,Intent(in):: direction 
    Call multihist_checkinit(histobj)
    multihist_getnbins = histobj%nbins(direction)
    Return
  End Function multihist_getnbins

  !----------------------------------------------------------
  ! If the bounds of bin i are x_i, x_i+1 then the values which
  ! belong to the bin are given by [x_i, x_i+1) ie the lower value
  ! is inclusive.
  !----------------------------------------------------------
  Subroutine multihist_update(histobj, value)
    Type(Multihist_Params), Intent(inout) :: histobj
    Real(kind=RDbl), Dimension(:),Intent(in)  :: value

    Integer   :: binno, nelements, i
    !** Check dimension of input array
    Call multihist_checkinit(histobj)

    If (Size(value,1)/=histobj%dim) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "WRONG Values passed to histogram object"
      Stop
    Else
      Do i=1,histobj%dim
        !** Make sure the value is within bounds
        If ( (value(i) < histobj%lowvalue(i)) .Or. &
            (value(i) > histobj%highvalue(i)) ) Then
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              " The value is not within the initialized bounds"
          Write(0,'(a, 2e16.8)') 'Bounds: ', histobj%lowvalue(i),&
              histobj%highvalue(i)
          Write(0,'(a, e16.8)')  'Value : ', value(i)
          Stop
        End If
      End Do
    
    !** Update the multihist and the no. of elements
    nelements=histobj%nelements
    binno=0
    Do i=1,histobj%dim
      !** Note temp-bin numers start from 0..
      histobj%temp_bin_no(i) = ( value(i) - histobj%lowvalue(i) ) / &
          ( histobj%binsize(i) ) 
      binno = binno + (histobj%temp_bin_no(i) ) * histobj%multi_fac(i)
      Call stats_update(histobj%stat(i),value(i))
      If (value(i)>histobj%max(i)) Then
        histobj%max(i)=value(i)
      Else If (value(i)<histobj%min(i))  Then
        histobj%min(i)=value(i)
      Endif
    End Do

    !** but our arrays are indexed from 1..
    binno=binno+1
    
    histobj%histg(binno) = histobj%histg(binno) + 1
    histobj%nelements = nelements + 1
    Endif
  End Subroutine multihist_update

  !----------------------------------------------------------
  ! Make the area under the multihist unity
  !----------------------------------------------------------
  Subroutine multihist_normalize(histobj)
    Type(Multihist_Params), Intent(inout) :: histobj
    
    Integer   :: i
    Real(kind=RDbl) :: dV
    Call multihist_checkinit(histobj)
    dV=histobj%totVolume/histobj%totBins ! diff volume

    !** Integral of probability should be one.
    !** Divide each bin by total_no_of_elements*binsize
    If ( (histobj%nelements*dV) > zero ) &
    histobj%probability = histobj%histg/(histobj%nelements*dV)

  End Subroutine multihist_normalize

  !----------------------------------------------------------
  ! Do the necessary cleanup such as freeing up memory
  !----------------------------------------------------------
  Subroutine multihist_cleanup(histobj)
    Type(Multihist_Params), Intent(out) :: histobj

    Integer :: error
    
    Deallocate(histobj%histg, Stat=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"dealloc-histobj")
    Deallocate(histobj%probability, Stat=error)
    If (error/=0) Call multihist_allocErrDisplay(__LINE__,"dealloc-histobj")
  End Subroutine multihist_cleanup

  !----------------------------------------------------------
  ! Write the multihist to the optional unit "optunitno"
  !----------------------------------------------------------
  Subroutine multihist_display(histobj, optunitno)
    Type(Multihist_Params), Intent(inout) :: histobj
    Integer, Intent(in), Optional  :: optunitno
    Integer :: unitno
    Call multihist_checkinit(histobj)
    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6
    End If
    Select Case(histobj%dim)
    Case(1) 
      Call multihist_display1D(histobj, unitno)
!    Case(2) 
!      Call multihist_display2D(histobj, unitno)
    Case(3) 
      Call multihist_display3D(histobj, unitno)
    Case Default
      Write(*,*) "No display routine available currently"
      Stop
    End Select
    
  End Subroutine multihist_display
  
  !----------------------------------------------------------
  ! Write the multihist to the optional unit "optunitno"
  !----------------------------------------------------------
  Subroutine multihist_display1D(histobj, optunitno)
    Type(Multihist_Params), Intent(inout) :: histobj
    Integer, Intent(in), Optional  :: optunitno

    Integer   :: unitno, i
    Real(kind=RDbl) :: value
    Call multihist_checkinit(histobj)
    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6
    End If

    Call multihist_normalize(histobj)

    Do i=1, histobj%totBins
      value  = multihist_getvalue(histobj,1,i)
      Write(unitno, '(i10, e16.6, i10, f12.6)') i, value, histobj%histg(i), &
          histobj%probability(i)
    End Do

    Write(unitno,*) "cumAvg,  cumAvSq, Min, Max "
    Write(unitno, '(4e16.6)') histobj%stat(1)%cumavg,histobj%stat(1)%cumsqavg,&
        histobj%min(1),histobj%max(1)

  End Subroutine multihist_display1D

  !----------------------------------------------------------
  ! Write the multihist to the optional unit "optunitno"
  !----------------------------------------------------------
  Subroutine multihist_display3D(histobj, optunitno)
    Type(Multihist_Params), Intent(inout) :: histobj
    Integer, Intent(in), Optional  :: optunitno
    
    Integer   :: unitno, i,j,k,bin_index
    Real(kind=RDbl) :: value,x,y,z
        Call multihist_checkinit(histobj)
    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6
    End If
    

    Write(*,*) "DISPLAY STARTED"
    Write(unitno,*) "dim=",histobj%dim,"No of bins, x-",histobj%nbins(1) &
        ," y-",histobj%nbins(2)," z-",histobj%nbins(3)
    Write(unitno,*) histobj%multi_fac
    Call multihist_normalize(histobj)
    bin_index=0
    Do i=1,histobj%nbins(HIST_XDIM)
      x=multihist_getValue(histobj,HIST_XDIM,i)
      Do j=1,histobj%nbins(HIST_YDIM)
        y=multihist_getValue(histobj,HIST_YDIM,j)
        Do k=1,histobj%nbins(HIST_ZDIM)
          bin_index=bin_index+1
          z=multihist_getValue(histobj,HIST_ZDIM,k)
          If (histobj%histg(bin_index)>zero) Then
            Write(unitno,'(3f10.5,i9,e14.5,i9)') x,y,z,&
                bin_index,histobj%probability(bin_index)&
                ,histobj%histg(bin_index)
          Endif
        End Do
      End Do
    End Do
    
    Write(unitno,*) "Dim , cumAvg,  cumAvSq, Min, Max "
    Do i=1,histobj%dim
      Write(unitno, '(i4,4e16.6)') i,histobj%stat(i)%cumavg,&
          histobj%stat(i)%cumsqavg,&
          histobj%min(i),histobj%max(i)
    End Do
    
  End Subroutine multihist_display3D
  
  
  !------------------------------------------------------------------
  !Returns the value in the given dimension corresponding to "val_index"
  !------------------------------------------------------------------
  Real(kind=RDbl) Function multihist_getValue(histobj,dim_index,val_index)
  Type(Multihist_Params), Intent(inout) :: histobj
  Integer, Intent(in)  :: dim_index,val_index
    Call multihist_checkinit(histobj)
    multihist_getValue= histobj%lowvalue(dim_index) +  &
          (val_index-0.50_RDbl)*histobj%binsize(dim_index)
  End Function multihist_getValue

  !--------------------------------
  !Checks whether initlzd, Stops if not
  !--------------------------------
  Subroutine multihist_checkinit(histobj)
    Type(Multihist_Params), Intent(in) :: histobj  
    If (.Not.(histobj%is_init)) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Without initialisation a histogram is being accessed"
      Stop
    Endif
    
  End Subroutine multihist_checkinit

  !--------------------------------------------------------------
  ! displays pointer allocation errors, and stops the simulation
  !--------------------------------------------------------------
  Subroutine multihist_allocErrDisplay(lineno,opt_variableName)
    Integer,Intent(in)    :: lineno
    Character(*) ,Intent(in),Optional :: opt_variableName

    Character(len=strlen) :: varname, filename
    
    If (Present(opt_variableName)) Then
      varname=opt_variableName
    Else
      varname=" "
    Endif
    
    filename= __FILE__
    
    Write(*,'(1x,a,i4)') "Error while allocating pointer variable " &
        //Trim(varname)//" in - "//Trim(filename)//" at line  no :",lineno

    Stop

  End Subroutine multihist_allocErrDisplay

End Module multihist

