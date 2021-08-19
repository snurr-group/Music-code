!--------------------------------------------------------------------------
! Deals with command line input
!
! Important routine: commandline_init
!
! Idea is take some values from command line input and set some fields
! (of module general.F90) that are common to the whole simulation
!--------------------------------------------------------------------------

Module commandline

  Use defaults, Only : strLen, d_ctrl_file
  Use utils, Only : allocErrDisplay
  Use file, Only : file_open
  Use general, Only : genparams
#ifdef NAGCOMPILER
  Use f90_unix, Only: getarg,iargc
#endif

  Implicit None
  Save

  Private 
  Public :: commandline_init

  Logical                   :: CmndLineInit=.False.
  Character(len=strLen)     :: ProgName
  Character(len=strLen), Dimension(:), Pointer :: CmndLineInputs

Contains
  !-------------------------------------------------------------------------
  ! Reads the command line and does various important tasks
  ! 1) Stores the command line inputs, and name of the executable
  ! 2) checks the command line options typed
  ! 3) opens the file and returns the ctrlfile name
  ! 4) returns a recommended action to the calling program
  ! 5) returns uninterpretable command line pieces if optional space passed
  ! Requires:  ctrl_filename -- output control file filename
  !            action -- recommand action
  !            xtra_args -- optional string array space for xtra arguments
  ! NOTE: it will not open the control filename (last argument) if the
  !       the optional 'xtra_args' is passed
  !-------------------------------------------------------------------------
  Subroutine commandline_init(ctrl_filename,action,xtra_args)
    Character(len=strLen), Intent(Out)                :: ctrl_filename,action
    Character(*), Dimension(:), Intent(Out), Optional :: xtra_args

    Integer          :: i,error
    Integer          :: no_of_args,argno,ctrl_unit,nxtra
    Logical          :: ctrl_file_found,error_flag
    Character(len=1) :: firstchar,secondchar
#ifndef NAGCOMPILER
    Integer          :: iargc !** iargc is a fortran function
#endif

    !** iargc returns the number of command line arguments, 
    no_of_args = iargc()

    !** Use getarg(i,arg) to return the i'th command line argument
    !** name of executable is zero'th argument
    Call getarg(0, ProgName)

    !** Dump help to screen if only the program name was input
    If (no_of_args < 1) Then
      Call commandline_writehelp()
      Stop
    End If

    !** Allocate space for the rest of the command line chunks and get 'em
    Allocate(CmndLineInputs(no_of_args),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1,no_of_args
      Call getarg(i, CmndLineInputs(i))
    End Do

    !** Take appropriate action for each argument on command line
    action = ''
    ctrl_file_found = .False.
    Do argno = 1,no_of_args
      If (argno == no_of_args) Then
        !** Now this should be the control file name, and that should 
        !** be the last argument entered at the command line
        !** its still possible that this is not ctrlfile. It could be --samplecf
        !** or --help options. So we will weight for other checks to finish
        !** befor opening the ctrlfile
        ctrl_filename = CmndLineInputs(no_of_args)
        ctrl_file_found=.True.

        If (Index(ctrl_filename,'@') /= 0) Then  !** it's not a filename
          Return   !** it's some sort of command, return it 
        End If

        !** Check whether ctrlfile exists, and open it if exists
        If (.Not. Present(xtra_args)) Then
          Write(*,*) "Looking for Control file : ",Trim(ctrl_filename)
          ctrl_unit = file_open(ctrl_filename,110,d_ctrl_file)
          ctrl_file_found = .True.
          action = "DoSimulation"
        End If

        Return
      End If

      !** Interpret or store the other types of command line arguments
      nxtra = 0
      error_flag = .False.
      firstchar = CmndLineInputs(argno)(1:1)

      If (firstchar == '-') Then
        secondchar = CmndLineInputs(argno)(2:2)

        !** Handle the processing based on the number of preceeding '-' chars
        If (secondchar == '-') Then  !** Handle double '--' command line chunks

          ctrl_file_found = .False. ! ctrlfile should be last
          Select Case(Trim(CmndLineInputs(argno)))
          Case ('--samplecf') 
            action = 'WriteSampleCtrlfile'
            Return
          Case ('--help') 
            Call commandline_writehelp()
            Stop
          Case Default
            error_flag = .True.
          End Select
        Else   !** Handle single '-' command line chunks

          ctrl_file_found = .False.       ! ctrlfile should be last
          Select Case(Trim(CmndLineInputs(argno)))
          Case ('-VERBOSE') 
            genparams%displaymode="VERBOSE"
          Case ('-SILENT') 
            genparams%displaymode="SILENT"
          Case ('-NORMAL') 
            genparams%displaymode="NORMAL"
          Case Default
            !** Either store this uninterpretable chunk or set error flag
            If (Present(xtra_args)) Then
              nxtra = nxtra + 1
              xtra_args(nxtra) = CmndLineInputs(argno)(2:)
            Else
              error_flag = .True.
            End If
          End Select
        End If

      Else
        If (argno == no_of_args) Then
          !  this is the last argument, most probbaly ctrlfile
        Else
          error_flag = .True.       
        End If
      End If

      !** Give the standard error message if something is uninterpretable
      If (error_flag) Then
        Call commandline_writehelp()
        Write(0,'(2a)') 'Could not interpret command line option: ', &
            Trim(CmndLineInputs(argno))
        Stop
      End If

    End Do     !** end of commandline chunk interpretation loop

    !** Check whether ctrlfile exists, and open it if exists
    If (ctrl_file_found) Then
      Write(*,*) "Looking for Control file : ",Trim(ctrl_filename)
      ctrl_unit = file_open(ctrl_filename,110,d_ctrl_file)
      action = "DoSimulation"
      Return
    Else   
      !** Complain if the control file name was not found
      If (.Not. ctrl_file_found) Then
        Call commandline_writehelp()
        Write(0,'(2a)') 'Could not find control file name on command line'
        Stop
      End If
    End If
    
  End Subroutine commandline_init

  !---------------------------------------------------------------------
  ! Help for the clueless person trying to use the music executable
  ! Requires:  nothing
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
