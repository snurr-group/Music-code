!-------------------------------------------------------------------------
! This module handles the timing of routines.  It contains a declared data
! structure that stores the start-time, total-time and number of calls
! to wtime for each tracked routine.
!  - wtime_initall and wtime_initone initialize the clock and zero the rest
!  - wtime_starttime resets the clock
!  - wtime_stoptime increments the total-time and the number of calls
!  - wtime_display displays the elapsed time and time per call
! NOTE : uses external subroutine cpu_time for evaluating cpu time spent
!-------------------------------------------------------------------------

Module wtime

  Use defaults, Only: strLen, lstrLen, MAX_ROUTINES, RDbl
  Use utils, Only: int2str, real2str

  Implicit None
  Save

  Private
  Public :: wtime_initall, wtime_initone, wtime_init, wtime_starttime, &
      wtime_stoptime, wtime_display

  Type Timing_Info
    Real :: cpustarttime
    Real :: totaltime
    Integer         :: ncalls
  End Type Timing_Info

  Type(Timing_Info), Dimension(MAX_ROUTINES)      :: simtime

  Interface wtime_init
    Module Procedure wtime_initall
    Module Procedure wtime_initone
  End Interface
  
Contains        
  !-----------------------------------------------------------
  ! This initializes the array that holds the times for each
  ! routine that we are timing
  !-----------------------------------------------------------
  Subroutine wtime_initall

    Real      :: cputime
    Integer   :: i

    !** cpu_time - external subroutine part of fortran90
    Call cpu_time(cputime)

    Do i=1, MAX_ROUTINES
      simtime(i)%cpustarttime = cputime
      simtime(i)%totaltime = 0.0
      simtime(i)%ncalls = 0
    End Do

  End Subroutine wtime_initall

  !-----------------------------------------------------------
  ! This initializes the time value for the routine 
  ! indexed by "routineno"
  !-----------------------------------------------------------
  Subroutine wtime_initone(routineno)
    Integer, Intent(in) :: routineno
    Real      :: cputime

    !** cpu_time - external subroutine part of fortran90
    Call cpu_time(cputime)

    simtime(routineno)%cpustarttime = cputime
    simtime(routineno)%totaltime = 0.0
    simtime(routineno)%ncalls = 0
  End Subroutine wtime_initone
  
  !-----------------------------------------------------------
  ! Reset the starttime of routine indexed by "routineno"
  !-----------------------------------------------------------
  Subroutine wtime_starttime(routineno)
    Integer, Intent(in)  :: routineno
    Real               :: cputime
    Real, External     :: etime

    !** Get the current CPU time
    !** cpu_time - external subroutine part of fortran90
    Call cpu_time(cputime)
    
    simtime(routineno)%cpustarttime = (cputime*1.0)

  End Subroutine wtime_starttime

  !----------------------------------------------------------------
  ! This routine updates the time spent in the routine "routineno"
  !----------------------------------------------------------------
  Subroutine wtime_stoptime(routineno)
    Integer, Intent(in)  :: routineno
    Real               :: cputime, etime
    Real    :: deltatime
    
    !** Get the current CPU time
    !** cpu_time - external subroutine part of fortran90
    Call cpu_time(cputime)

    !** store the elapsed time
    deltatime = cputime - simtime(routineno)%cpustarttime
    simtime(routineno)%totaltime = simtime(routineno)%totaltime + deltatime
    simtime(routineno)%ncalls = simtime(routineno)%ncalls + 1
  End Subroutine wtime_stoptime
  
  !---------------------------------------------------------
  ! Print out the time for the routine routineno
  ! Requires: label -- string label for output beginning
  !           routineno -- routine number
  !           optunitno -- optional output unit number
  !---------------------------------------------------------
  Subroutine wtime_display(label,routineno,optunitno)
    Character(*),  Intent(In)      :: label
    Integer, Intent(In)            :: routineno
    Integer, Optional, Intent(In)  :: optunitno
    
    Integer                :: unitno
    Real(kind=RDbl)        :: elapsedtime          !* CPU time elapsed
    Real(kind=RDbl)        :: tpercall,totaltime   !* change of units
    Character(len=strLen)  :: tpc_units,total_units
    Character(len=strLen)  :: string1,string2,string3
    Character(len=lstrLen) :: outstr 

    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6
    End If

    !** set the total time and associated units
    elapsedtime = simtime(routineno)%totaltime
    If (elapsedtime > 86400.0e0) Then
      totaltime = elapsedtime/86400.0e0
      total_units = "days"
    Else If (elapsedtime > 3600.0e0) Then
      totaltime = elapsedtime/3600.0e0
      total_units = "hrs"
    Else If (elapsedtime > 60.0e0) Then
      totaltime = elapsedtime/60.0e0
      total_units = "min"
    Else If (elapsedtime < 0.1e0) Then
      totaltime = elapsedtime*1000.0e0
      total_units = "msec"
    Else
      totaltime = elapsedtime
      total_units = "sec"
    End If

    !** set the time per call and associated units
    tpercall = elapsedtime/Real(simtime(routineno)%ncalls)
    If (tpercall > 86400.0e0) Then
      tpercall = tpercall/86400.0e0
      tpc_units = "days"
    Else If (tpercall > 3600.0e0) Then
      tpercall = tpercall/3600.0e0
      tpc_units = "hrs"
    Else If (tpercall > 60.0e0) Then
      tpercall = tpercall/60.0e0
      tpc_units = "min"
    Else If (tpercall < 0.1e0) Then
      tpercall = tpercall*1000.0e0
      tpc_units = "msec"
    Else
      tpc_units = "sec"
    End If

    string1 = real2str(totaltime,6)
    string2 = int2str(simtime(routineno)%ncalls)
    string3 = real2str(tpercall,6)
    Write(outstr,'(11a)') Trim(label),'   Time: ',Trim(string1),' ', &
        Trim(total_units),'   Number of calls: ',Trim(string2), &
        '   Time per call: ',Trim(string3),' ',Trim(tpc_units)

    Write(unitno, '(a)') Trim(Adjustl(outstr))

  End Subroutine wtime_display

End Module wtime


