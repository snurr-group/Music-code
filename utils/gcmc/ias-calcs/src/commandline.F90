!----------------------------------------------------
! Deals with command line input
! Important routine : commandline_init
! Idea is take some values from command line input 
!----------------------------------------------------
Module commandline

  Use defaults, Only : strLen, d_ctrl_file
  Use utils, Only : allocErrDisplay
  Use file, Only : file_open
#ifdef NAGCOMPILER
  Use f90_unix, Only: getarg,iargc
#endif

  Implicit None
  Save

  Private 
  Public :: commandline_init


  Character(len=strLen),Dimension(:), Pointer :: CmndLineInputs
  Character(len=strLen) :: ProgName
  Logical :: CmndLineInit=.False.

Contains
  !--------------------------------------------------------
  !- Stores the command line inputs, and name of the executable
  !- checks the command line options typed
  !- returns the ctrlfile name
  !- returns a recommended action to the calling program
  !--------------------------------------------------------
  Subroutine commandline_init(ctrl_filename,action)
    Character(len=strLen),Intent(out) :: ctrl_filename,action

    Integer :: no_of_args,i,error,argno,ctrl_unit
#ifndef NAGCOMPILER
    Integer :: iargc !** iargc is a fortran function
#endif
    Logical :: ctrl_file_found

    ! ** iargc returns the number of command line arguments, 
    no_of_args=iargc()
#ifdef NAGCOMPILER
    no_of_args = 1
#endif

    !** getarg(i,arg) returns the i'th command line argument
    !** zero 'th argument is the command itself

    !** name of executable
    Call getarg(0, ProgName)

    If (no_of_args<1) Then
      Call commandline_writehelp()
      stop
    Endif

    Allocate(CmndLineInputs(no_of_args),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Do i=1,no_of_args
      Call getarg(i, CmndLineInputs(i))
    End Do

    !** Take appropriate action for each option entered
    ctrl_file_found=.False.
    Do argno=1,no_of_args
      
      If (Trim(CmndLineInputs(argno))=='--samplecf') Then
        action='WriteSampleCtrlfile'
        Return
      Elseif  (Trim(CmndLineInputs(argno))=='--help') Then
        Call commandline_writehelp()
        Stop
      Elseif (argno==no_of_args) Then
        !** Now this has to be control file name, and that should be 
        !** the last argument entered at the command line

        !** Check whether ctrlfile exists, and open it if exists
        ctrl_filename=CmndLineInputs(no_of_args)
        Write(*,*) "Looking for Control file : ",Trim(ctrl_filename)
        ctrl_unit= file_open(ctrl_filename,110,d_ctrl_file)
        ctrl_file_found=.True.
        action="DoSimulation"
        Return
      Else 
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Call commandline_writehelp()
        Write(*,*) "Looks like something wrong with way you typed"
        !** or something wrong With shaji
        Write(*,*) "Wrong Option : ", Trim(CmndLineInputs(argno))
        Stop
      Endif

    End Do  ! end of options loop

    If (.Not.ctrl_file_found) Then
      Call commandline_writehelp()
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    
  End Subroutine commandline_init

  !---------------------------------------------------------------------
  ! Help for the clueless person trying to use the music executable
  !---------------------------------------------------------------------
  Subroutine commandline_writehelp()
      Write(0,*) 'Usage: ',Trim(progname),' [options] run_filename'
      Write(0,*) '   available options are:'
      Write(0,*) '   --help      : repeats this message'
      Write(0,*) '   --samplecf  : Writes a sample control file to output'
      Write(0,*) '   -VERBOSE, or -SILENT, or -NORMAL'
      Write(0,*) '               : display mode '
      Write(0,*) 
      Write(0,*) "Please Enter run_filename as the last argument"


  End Subroutine commandline_writehelp

End Module commandline
