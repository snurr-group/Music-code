!-----------------------------------------------------------------------
! This contains the main program for extrapolating/interpolating
! adsorption Isotherm data and calculating surface pressure using Myer's
! method for binary IAS calculations
!-----------------------------------------------------------------------
Program IAS

  Use commandline, Only: commandline_init
  Use file, Only: file_open, file_settag
  Use utils, Only: genfilename,allocErrDisplay,int2str
  Use defaults, Only: dashedline, strLen, d_ctrl_file
  Use sigmacalc, Only: sigmacalc_init, sigmacalc_write
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
!!$    Call sigmacalc_sampleCF(0)
    Write(0,'(a)') dashedline
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Elseif(action/="DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Endif
      
  Call file_settag(ctrl_filename, d_ctrl_file)

  !* Read the ctrlfile
  Call sigmacalc_init(ctrl_filename)
  Write(*,*) "Analyzing and writing the results"
  !** calculate sigma and write it
  Call sigmacalc_write()
  Write(*,*) "Finished"
End Program IAS



