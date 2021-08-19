!--------------------------------------------------------------------
! This module contains a data type for storing a value and 
! statistics for this value.
! 
! It can also accomodate Umbrella sampling style averages:
!             1/N sum_i A_i*weight_i     1/N sum_i A_i*weight_i  
!   <A>_o =  ------------------------ = ------------------------ 
!               1/N sum_i weight_i        normalization factor   
! Note that the difference in this case is that the input quantity 
! is the quantity being averaged times a weighting factor.  A 
! normalization constant is then required to get a proper average.
! See the header of the umbrella module for more notes.  I'm not sure
! how to apply the standard deviation to the umbrella average 
! quantities, so the std calculation ignores the normalization factor.
! 
! Important routines:
!   stats_init -- initilize the data type
!   stats_update -- update the instant value and recalculate stats
!   stats_getvalue -- get the instantaneous value
!   stats_getblock -- get the current block average
!   stats_getcavg -- get the cumulative average
!   stats_getstd -- get the standard deviation
!   stats_reset -- reset the counters in the data type
!   stats_copy -- copies all fields
!
! Needed Improvements:
! 1) Recalculating the statistics at every update is SLOW, should 
!    be reworked so that it's smarter.
!--------------------------------------------------------------------

Module stats

  Use defaults, Only: RDbl, strLen, lstrLen

  Implicit None
  Save

  Private
  Public :: Statistics, stats_getvalue, stats_getblock, stats_getcavg, &
      stats_getstd, stats_update, stats_init, stats_setvalue,stats_display, &
      stats_reset, stats_getfield, Assignment(=)

  Type Statistics
    Integer                   :: blocksize,blockiters,cumiters
    Logical                   :: usampling
    Real(kind=RDbl)           :: inst,weight
    Real(kind=RDbl)           :: blockavg
    Real(kind=RDbl)           :: cumavg
    Real(kind=RDbl)           :: cumsqavg
    Real(kind=RDbl)           :: cumdevsq
    Real(kind=RDbl)           :: blknormsum,normsum
    Real(kind=RDbl)           :: blknormfactor,normfactor
    Character(len=lstrLen)    :: label
    Character(len=strLen)     :: format
  End Type Statistics

  Interface Assignment(=)
    Module Procedure stats_copy
  End Interface

Contains

  !----------------------------------------------------------------------
  ! Initialize the data and various counters
  ! Requires:  obj -- data structure to initialize
  !            label -- a string label
  !            blocksize -- block size for analysis
  !            usample -- flag indicating if Umbrella sampling is used
  !            format -- optional format for the displays
  !----------------------------------------------------------------------
  Subroutine stats_init(obj, label, blocksize, usample, format)
    Type(Statistics), Intent(InOut)    :: obj
    Character(*), Intent(In)           :: label
    Integer, Intent(In)                :: blocksize
    Logical, Intent(In)                :: usample
    Character(*), Intent(In), Optional :: format
    
    obj%label = label
    If (Present(Format)) Then
!       Write(obj%format,*) "(a, ' : ', 4"//Adjustl(Trim(Format))//")"
      obj%format= "(a, 4"//Adjustl(Trim(Format))//")"
    Else
      obj%format = "(a, 4e13.6)"
    End If

    obj%usampling  = usample

    obj%blocksize  = blocksize
    obj%blockiters = 0
    obj%cumiters   = 0
    obj%inst       = 0.0_RDbl
    obj%weight     = 1.0_RDbl
    obj%blockavg   = 0.0_RDbl
    obj%cumavg     = 0.0_RDbl
    obj%cumsqavg   = 0.0_RDbl
    obj%cumdevsq   = 0.0_RDbl

    obj%blknormsum    = 0.0_RDbl
    obj%normsum       = 0.0_RDbl
    obj%blknormfactor = 1.0_RDbl
    obj%normfactor    = 1.0_RDbl

  End Subroutine stats_init

  !--------------------------------------------------------------
  ! Reset the data and various counters
  ! Requires:  obj -- statistics data structure 
  !--------------------------------------------------------------
  Subroutine stats_reset(obj)
    Type(Statistics), Intent(InOut)  :: obj

    obj%blockiters = 0
    obj%cumiters   = 0
    obj%inst       = 0.0_RDbl
    obj%weight     = 1.0_RDbl
    obj%blockavg   = 0.0_RDbl
    obj%cumavg     = 0.0_RDbl
    obj%cumsqavg   = 0.0_RDbl
    obj%cumdevsq   = 0.0_RDbl

    obj%blknormsum    = 0.0_RDbl
    obj%normsum       = 0.0_RDbl
    obj%blknormfactor = 1.0_RDbl
    obj%normfactor    = 1.0_RDbl

  End Subroutine stats_reset

  !-----------------------------------------------------------------
  ! Set the initial value of the fields to "value"
  ! Requires:  obj -- statistics data structure 
  !            value -- intitial value
  !-----------------------------------------------------------------
  Subroutine stats_setvalue(obj, value)
    Type(Statistics), Intent(InOut)  :: obj
    Real(kind=RDbl), Intent(In)      :: value

    obj%inst      = value
    obj%blockavg  = value
    obj%cumavg    = value
    obj%cumsqavg  = value*value
    obj%cumdevsq  = 0.0_RDbl

  End Subroutine stats_setvalue

  !----------------------------------------------------------------------
  ! Update the various means and calculate the standard deviation 
  ! over the "step" steps
  ! Requires:  obj -- statistics data structure 
  !            value -- new, current value
  !            weight -- weighting factor, required for Umbrella sampling
  !----------------------------------------------------------------------
  Subroutine stats_update(obj, value, weight)
    Type(Statistics), Intent(InOut)          :: obj
    Real(kind=RDbl), Intent(In)              :: value
    Real(kind=RDbl), Intent(In), Optional    :: weight

    Real(kind=RDbl)          :: val

    !** Weight the value if we're doing Umbrella sampling
    If (obj%usampling) Then
      val = weight*value
      obj%normsum = obj%normsum + weight
    Else
      val = value
    End If
    obj%inst = val
    
    !** Calculate the block average
    obj%blockiters = obj%blockiters + 1
    obj%blockavg = (obj%blockavg*(obj%blockiters-1.0_RDbl) + value)/ &
        (obj%blockiters*1.0_RDbl)
    If (obj%blockiters == obj%blocksize) Then
      obj%blockiters = 0
    End If

    !** Calculate the cumulative average and the variance = sigma**2
    obj%cumiters = obj%cumiters + 1
    obj%cumavg   = (obj%cumavg  *(obj%cumiters-1.0_RDbl) + value)/ &
        (obj%cumiters*1.0_RDbl)
    obj%cumsqavg = (obj%cumsqavg*(obj%cumiters-1.0_RDbl) + value*value)/ &
        (obj%cumiters*1.0_RDbl)
    obj%cumdevsq = (obj%cumsqavg - obj%cumavg**2)

    !** Normalize the block and cumulative averages for Umbrella sampling
    If (obj%usampling) Then
      obj%normsum = obj%normsum + weight
      obj%blknormsum = obj%blknormsum + weight
      obj%normfactor = obj%normsum/obj%cumiters
      obj%blknormfactor = obj%blknormsum/obj%cumiters

      obj%blockavg = obj%blockavg/obj%blknormfactor
      obj%cumavg = obj%cumavg/obj%normfactor
    End If

  End Subroutine stats_update

  !--------------------------------------------------------------
  ! Returns the current value as stored in the object "obj"
  ! Requires:  obj -- statistics data structure 
  !--------------------------------------------------------------
  Real(kind=RDbl) Function stats_getvalue(obj)
    Type(statistics), Intent(In)   :: obj    

    stats_getvalue = obj%inst

  End Function stats_getvalue


  !--------------------------------------------------------------
  ! Returns the value as stored in the object "obj"
  ! field is a string which decides what to return
  ! Requires:  obj -- statistics data structure 
  !            field -- string; tells what is needed
  !--------------------------------------------------------------
  Real(kind=RDbl) Function stats_getfield(obj,field)
    Type(statistics), Intent(In)   :: obj  
    Character(len=strLen)  , Intent(In)     :: field  
    Select Case(Trim(field))
    Case ('INST')
      stats_getfield = obj%inst
    Case ('CAVG')
      stats_getfield = obj%cumavg
    Case ('BLOCK')
      stats_getfield = obj%blockavg
    Case ('STD')
      stats_getfield = stats_getstd(obj)
    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "stats did not understand the field passed here",&
          Trim(field)
      Write(*,*) " legitimate field values : CAVG or BLOCK or INST or STD" 
      Stop
    End Select
  End Function stats_getfield

  !--------------------------------------------------------------
  ! Returns the block average as stored in the object "obj"
  ! Requires:  obj -- statistics data structure 
  !--------------------------------------------------------------
  Real(kind=RDbl) Function stats_getblock(obj)
    Type(statistics), Intent(In)   :: obj    

    stats_getblock = obj%blockavg

  End Function stats_getblock

  !--------------------------------------------------------------
  ! Returns the cumulative average value stored in "obj"
  ! Requires:  obj -- statistics data structure 
  !--------------------------------------------------------------
  Real(kind=RDbl) Function stats_getcavg(obj)
    Type(statistics), Intent(In)   :: obj    

    stats_getcavg = obj%cumavg

  End Function stats_getcavg

  !--------------------------------------------------------------
  ! Returns the std value stored in "obj"
  ! Requires:  obj -- statistics data structure 
  !--------------------------------------------------------------
  Real(kind=RDbl) Function stats_getstd(obj)
    Type(statistics), Intent(In)   :: obj    

    If (obj%cumdevsq < 0.0_RDbl) Then
      stats_getstd = 0.0_RDbl
    Else
      stats_getstd = Sqrt(obj%cumdevsq)
    End If

  End Function stats_getstd

  !--------------------------------------------------------------
  ! Returns an array of values such that index 1) inst. value
  ! 2) block avg. 3) cum avg. 4) stddev
  ! Requires:  obj -- statistics data structure 
  !--------------------------------------------------------------
  Function stats_getvalues(obj)
    Type(Statistics), Intent(In)   :: obj
    Real(kind=RDbl), Dimension(4)  :: stats_getvalues

    stats_getvalues(1) = obj%inst
    stats_getvalues(2) = obj%blockavg
    stats_getvalues(3) = obj%cumavg
    stats_getvalues(4) = Sqrt(obj%cumdevsq)

  End Function stats_getvalues

  !----------------------------------------------------------------
  ! Copies all the fields from fromStat to ToStat
  ! Requires:  toStat -- destination statistics data structure 
  !            fromStat -- origin statistics data structure 
  !----------------------------------------------------------------
  Subroutine stats_copy(toStat, fromStat)
    Type(Statistics), Intent(In)   :: fromStat
    Type(Statistics), Intent(Out)  :: toStat

    ToStat%label      = fromStat%label
    ToStat%usampling  = fromStat%usampling
    ToStat%format     = fromStat%format       
    ToStat%blocksize  = fromStat%blocksize    
    ToStat%blockiters = fromStat%blockiters   
    ToStat%cumiters   = fromStat%cumiters     
    ToStat%inst       = fromStat%inst         
    ToStat%weight     = fromStat%weight
    ToStat%blockavg   = fromStat%blockavg     
    ToStat%cumavg     = fromStat%cumavg       
    ToStat%cumsqavg   = fromStat%cumsqavg     
    ToStat%cumdevsq   = fromStat%cumdevsq      

    ToStat%blknormsum    = fromStat%blknormsum   
    ToStat%normsum       = fromStat%normsum      
    ToStat%blknormfactor = fromStat%blknormfactor
    ToStat%normfactor    = fromStat%normfactor   

  End Subroutine stats_copy

  !------------------------------------------------------------
  ! Print out the various statistics
  ! Requires:  obj -- statistics object
  !            indent -- number of spaces from left margin
  !            optunit -- optional output unit
  !            optflag -- optional flag ?
  !------------------------------------------------------------
  Subroutine stats_display(obj, indent, optunit,optflag)
    Type(statistics), Intent(In)  :: obj
    Integer, Intent(In)           :: indent
    Integer, Optional, Intent(In) :: optunit,optflag

    Integer                :: unitno
    Real(kind=RDbl)        :: stddev
    Character(len=lstrLen) :: label
    Character(len=indent)  :: blank

    blank = Repeat(' ',indent)
    unitno = 6
    If (Present(optunit)) unitno = optunit
    
    If (Present(optflag)) Then
      If (optflag==0) Then
        label = blank//Trim(obj%label)
      Else
        label = ''
      Endif
    Else
      label = blank//Trim(obj%label)
    Endif
    
    !** Prevent round-off error.  If the variance is zero
    !** then obj%cumdevsq can become negative due to the
    !** finite numerical precision.
    If (obj%cumdevsq < 0.0_RDbl) Then
      stddev = 0.0_RDbl
    Else
      stddev = Sqrt(obj%cumdevsq)     
    End If

    Write(unitno,obj%format) Trim(label), obj%inst, &
        obj%blockavg, obj%cumavg, stddev

  End Subroutine stats_display

End Module stats  








