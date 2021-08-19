!-----------------------------------------------------------------------
! This contains the main program for extrapolating/interpolating
! adsorption Isotherm data and calculating surface pressure using Myer's
! method for binary IAS calculations
!-----------------------------------------------------------------------
Program IAS

  Use commandline, Only: commandline_init
  Use file, Only: file_open
  Use utils, Only: genfilename,allocErrDisplay,int2str
  Use defaults, Only: 

  Implicit None

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action

  !-------------------------------------------------------------
  ! Do basic initialization, based on command line input
  !-------------------------------------------------------------
  Call commandline_init(ctrl_filename,action)

  If (action=="WriteSampleCtrlfile") Then
    Write(0,'(a)') dashedline
    Call sigmacalc_sampleCF(0)
    Write(0,'(a)') dashedline
    Stop
  Elseif(action/="DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Endif
      
  Call file_settag(ctrl_filename, d_ctrl_file)

  !* Read the ctrlfile
  Call sigmacalc_init(ctrl_filename)

  !** calculate sigma and write it
  Call sigmacalc_write()

End Program IAS



