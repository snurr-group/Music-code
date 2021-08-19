!-------------------------------------------------------------------
! This module provides the necessary functionality for generating a
! histogram.
!-------------------------------------------------------------------
Module histogram

  Use defaults, Only : RDbl, strLen, one, zero, STATS_BLOCKSIZE
  Use stats,    Only : Statistics, stats_init, stats_update

  Implicit None
  Save

  Private
  Public :: Histogram_Params, histogram_init, histogram_reinit, &
      histogram_update, histogram_display

  Type Histogram_Params
    Integer      :: nelements
    Integer      :: nbins
    Real(kind=RDbl) :: lowvalue
    Real(kind=RDbl) :: highvalue
    Real(kind=RDbl) :: min,max
    Type(Statistics) :: stat
    Real(kind=RDbl) :: binsize
    Real(kind=RDbl), Dimension(:), Pointer :: probability
    Integer , Dimension(:), Pointer :: histg
    Logical    :: is_init=.false.
    Character(len=strlen) :: tag
  End Type Histogram_Params

  Character(len=strlen) :: default_HIST_tag="untagged histogram"

  Interface histogram_init
    Module Procedure histogram_init1
    Module Procedure histogram_init2
  End Interface

Contains
  !----------------------------------------------------------
  ! Initializes the histogram object "histobj" by supplying 
  ! to it the no. of bins "nbins", the low value "lowvalue" and the
  ! high value "highvalue"
  !----------------------------------------------------------
  Subroutine histogram_init1(histobj, nbins, lowvalue, highvalue,opttag)
    Type(Histogram_Params), Intent(out) :: histobj
    Integer, Intent(in)    :: nbins
    Real(kind=RDbl), Intent(in) :: lowvalue
    Real(kind=RDbl), Intent(in) :: highvalue
    Character(len=strLen) ,Intent(in),Optional :: opttag

    Real(kind=RDbl) :: binsize
    Integer         :: error
    Character(len=strLen)  :: tag

    If (histobj%is_init) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) " attempt to re-initialize here"
      Stop
    Endif

    If (Present(opttag)) Then
      histobj%tag=opttag
    Else
      histobj%tag=default_HIST_Tag
    Endif

    histobj%nelements = 0
    histobj%nbins     = nbins
    histobj%lowvalue    = lowvalue
    histobj%highvalue   = highvalue
    !** sets the min and max to aribtrarily high and low values respectidely
    histobj%min=highvalue+one
    histobj%max=lowvalue-one
    ! STATS_BLOCKSIZE is the default blocksize
    Call stats_init(histobj%stat,'Hist. Object',STATS_BLOCKSIZE,.False.)
    

    histobj%binsize   = (highvalue - lowvalue)/nbins

    Allocate(histobj%histg(nbins), STAT = error)
    If (error/=0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, "Error while & 
          & allocating histogram variable. Tag : "//Trim(histobj%tag)
      Stop
    Endif

    Allocate(histobj%probability(nbins), STAT = error)
    If (error/=0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, "Error while & 
          & allocating histogram variable. Tag : "//Trim(histobj%tag)
      Stop
    Endif
    
    !set all elemnts to zero
    histobj%histg = 0
    histobj%probability = zero
    histobj%is_init=.true.
  End Subroutine histogram_init1

  !----------------------------------------------------------
  ! Initializes the histogram object "histobj" by supplying 
  ! to it the low value "low", the  high value "high" and
  ! the size of each bin "binsize"
  !----------------------------------------------------------
  Subroutine histogram_init2(histobj, lowvalue, highvalue, binsize)
    Type(Histogram_Params), Intent(out) :: histobj
    Real(kind=RDbl), Intent(in) :: lowvalue
    Real(kind=RDbl), Intent(in) :: highvalue
    Real(kind=RDbl), Intent(in) :: binsize

    Integer         :: error
    Integer         :: nbins

    histobj%nelements = 0
    nbins = (highvalue - lowvalue)/binsize
    histobj%nbins     = nbins
    histobj%lowvalue    = lowvalue
    histobj%highvalue   = highvalue
    histobj%binsize   = binsize

    Allocate(histobj%histg(nbins), STAT = error)
    If (error/=0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, "Error while & 
          & allocating histogram variable. Tag : "//Trim(histobj%tag)
      Stop
    Endif

    Allocate(histobj%probability(nbins), STAT = error)
    If (error/=0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, "Error while & 
          & allocating histogram variable. Tag : "//Trim(histobj%tag)
      Stop
    Endif

    !** sets the min and max to aribtrarily high and low values respectidely
    histobj%min=highvalue+one
    histobj%max=lowvalue-one

    ! STATS_BLOCKSIZE is the default blocksize
    Call stats_init(histobj%stat,'Hist. Object',STATS_BLOCKSIZE,.False.)
    !set all elemnts to zero
    histobj%histg = 0
    histobj%probability = zero
  End Subroutine histogram_init2

  !----------------------------------------------------------
  ! sets all values to zero, as during init
  !----------------------------------------------------------
  Subroutine histogram_reinit(histobj)
    Type(Histogram_Params), Intent(inout) :: histobj

    Real(kind=RDbl) :: binsize
    Integer         :: error
    Character(len=strLen)  :: tag

    If (.Not.(histobj%is_init)) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) " attempt to re-initialize here, but never initialized?"
      Stop
    Endif

    histobj%nelements = 0

    !** sets the min and max to aribtrarily high and low values respectidely
    histobj%min=histobj%highvalue+one
    histobj%max=histobj%lowvalue-one

    ! STATS_BLOCKSIZE is the default blocksize
    Call stats_init(histobj%stat,'Hist. Object',STATS_BLOCKSIZE,.False.)
    
    !set all elemnts to zero
    histobj%histg = 0
    histobj%probability = zero

  End Subroutine histogram_reinit

  !-----------------------------------------------------------
  ! Returns the size of each bin in "histobj"
  !-----------------------------------------------------------
  Real(kind=RDbl) Function histogram_getbinsize(histobj)
    Type(Histogram_Params), Intent(in) :: histobj

    Call histogram_checkinit(histobj)
    histogram_getbinsize = histobj%binsize
    Return
  End Function histogram_getbinsize
  
  !-----------------------------------------------------------
  ! Returns the no. of bins in "histobj"
  !-----------------------------------------------------------
  Integer Function histogram_getnbins(histobj)
    Type(Histogram_Params), Intent(in) :: histobj

    Call histogram_checkinit(histobj)

    histogram_getnbins = histobj%nbins
    Return
  End Function histogram_getnbins

  !----------------------------------------------------------
  ! If the bounds of bin i are x_i, x_i+1 then the values which
  ! belong to the bin are given by [x_i, x_i+1) ie the lower value
  ! is inclusive.
  !----------------------------------------------------------
  Subroutine histogram_update(histobj, value,silentFlag)
    Type(Histogram_Params), Intent(inout) :: histobj
    Real(kind=RDbl), Intent(in)  :: value
    Logical , Optional, Intent(in) :: silentFlag

    Integer   :: binno
    Integer   :: nelements

    Call histogram_checkinit(histobj)

    !** Update max and min irrespecive of whether within bounds or not
    if (value>histobj%max) histobj%max=value
    if (value<histobj%min) histobj%min=value

    !** Make sure the value is within bounds
    If (value < histobj%lowvalue .Or. (value > histobj%highvalue)) Then
      If (Present(silentFlag)) Then
        If (silentFlag) Return
      Endif
        
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " The value is not within the initialized bounds"
      Write(0,'(a, 2e16.8)') 'Bounds: ', histobj%lowvalue, histobj%highvalue
      Write(0,'(a, e16.8)')  'Value : ', value
      Return
    End If
    
    !** Update the histogram and the no. of elements
    nelements=histobj%nelements
    binno = (value - histobj%lowvalue)/(histobj%binsize) + 1
!    Write(*,'(a,i4,3f10.3)')"binno", binno,value, histobj%lowvalue,histobj%highvalue
!    Write(*,*)"binno", size(histobj%histg)
!    Write(*,*)"binno", binno,histobj%histg(binno)
!    Write(*,*) histobj%tag,histobj%lowvalue,histobj%highvalue
    histobj%histg(binno) = histobj%histg(binno) + 1


    Call stats_update(histobj%stat,value)

    histobj%nelements = nelements + 1

  End Subroutine histogram_update

  !----------------------------------------------------------
  ! Make the area under the histogram unity
  !----------------------------------------------------------
  Subroutine histogram_normalize(histobj)
    Type(Histogram_Params), Intent(inout) :: histobj
    
    Integer   :: i
    Real(kind=RDbl) :: factor

    Call histogram_checkinit(histobj)
    If (histobj%nelements>0) Then
      !** Divide each bin by total_no_of_elements*binsize
      histobj%probability = histobj%histg/(histobj%nelements*histobj%binsize)
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
  End Subroutine histogram_normalize

  !----------------------------------------------------------
  ! Do the necessary cleanup such as freeing up memory
  !----------------------------------------------------------
  Subroutine histogram_cleanup(histobj)
    Type(Histogram_Params), Intent(out) :: histobj

    Integer :: error
    
    Deallocate(histobj%histg, Stat=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " Could not DEallocate 'histobj%histg'"
      Stop
    End If
  End Subroutine histogram_cleanup

  !----------------------------------------------------------
  ! Write the histogram to the optional unit "optunitno"
  !----------------------------------------------------------
  Subroutine histogram_display(histobj, optunitno)
    Type(Histogram_Params), Intent(inout) :: histobj
    Integer, Intent(in), Optional  :: optunitno

    Integer   :: unitno, i
    Real(kind=RDbl) :: lbound, ubound, value, binsize
    Call histogram_checkinit(histobj)
    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6
    End If
    If (histobj%nelements>0) Then
      Call histogram_normalize(histobj)
    Else
      Write(unitno,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(unitno,*) " Warning : This histogram does not contain any data"
      Return
    Endif
    binsize = histobj%binsize
    Do i=1, histobj%nbins
      lbound = histobj%lowvalue +  (i-1)*binsize
      ubound = lbound + binsize
      value  = (lbound + ubound)/2.0_RDbl
      Write(unitno, '(i10, e16.6, i10, f12.6)') i, value, histobj%histg(i), &
          histobj%probability(i)
    End Do
    Write(unitno,*) "##   mean      mean(sq)     min       max "
    Write(unitno, '(a,4e16.6)') "## ", histobj%stat%cumavg,&
        histobj%stat%cumsqavg, histobj%min,histobj%max
    Write(unitno,*) " "

  End Subroutine histogram_display
  
  
  !--------------------------------
  !Checks whether initlzd, Stops if not
  !--------------------------------
  
  Subroutine histogram_checkinit(histobj)
    Type(Histogram_Params), Intent(in) :: histobj  
    If (.Not.(histobj%is_init)) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Without initialisation a histogram is being accessed"
      Stop
    Endif
    
  End Subroutine histogram_checkinit
    
End Module histogram





