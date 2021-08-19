!-----------------------------------------------------------------------
! Driver for displaying bmap/pmap
!-----------------------------------------------------------------------
Program music
  Use bmap, Only: BiasMap_Params, bmap_init, bmap_makeXYZ, DISPLAY_ONLY_FLAG
  Use commandline, Only: commandline_init
  Use defaults, Only: strlen, RDbl, dashedline, d_res_file, d_ctrl_file, &
      ROUT_MAIN,TOTAL_INDEX
  Use file, Only: file_settag, file_open
  Use simcell, Only: SimCell_Params, simcell_init, simcell_display, &
      simcell_samplecf, simcell_getzeoell, simcell_geteff
  Use utils, Only: genfilename,allocErrDisplay,int2str,filesrchstr,toupper, &
      stripcmnt

  Implicit None


  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action

  !** Simulation Cell Variable(s)
  Character(len=strLen),Parameter :: simcell_tag='Simulation Cell Information'
  Type(SimCell_Params)            :: scell

  Type(BiasMap_Params),Pointer            :: biasmap
  Character(len=strLen) :: bmap_TOXYZ_TAG="bmap to XYZ info"
  Character(len=strLen)           :: line,scaling,xyzfile,descr,biasfilename
  Character(len=3)                :: atsym
  
  Real(kind=RDbl),Dimension(3)    :: eff
  Logical                         :: logflag
  Integer, Dimension(3)           :: gridsize
  Integer                         :: unitno, lineno
  Real(kind=RDbl)                 :: tk, minpc, maxpc
  
  !-------------------------------------------------------------
  ! Do basic initialization, based on command line input
  !-------------------------------------------------------------
  Call commandline_init(ctrl_filename,action)

  If (action=='WriteSampleCtrlfile') Then
    Call simcell_sampleCF(0)
    Write(0,'(a)') dashedline

    Write(0,'(a)') '--- '//Trim(bmap_TOXYZ_TAG)//' --- '
    Write(0,'(a,t30,a)') 'Integer(3)','# nx, ny and nz for the new grid'
    Write(0,'(a,t30,a)') 'Character', '# LOGARITHMIC or ARITHMETIC - &
        & type of scaling'
    Write(0,'(a,t30,a)') 'Integer(2)', '# min % and max % to be dispalyed, &
        &0-100'
    Write(0,'(a,t30,a)') 'Character', '# name of new xyz file'
    Write(0,'(a,t30,a)') 'Character', '# Description to be written to xyz file'
    Write(0,'(a,t30,a)') 'Character', '# Atom symbol, ex: C, Si, O, CH3'
    Write(0,'(a,t30,a)') 'Character', '# pmap file which is used as the bmap '
    Write(0,'(a,t30,a)') 'Real', '# Temperature used for biasing, Deg.K'
    
    Stop
  Elseif(action/='DoSimulation') Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__,' : ',__LINE__
    Stop
  Endif
      
      
  Call file_settag(ctrl_filename, d_ctrl_file)

  !----------------------------------------------------------------------------
  ! Initialize atoms, molecules
  !----------------------------------------------------------------------------
  Call simcell_init(scell, ctrl_filename, simcell_tag)
  Call simcell_display(scell,6,6)

  eff=simcell_geteff(scell, .True.)
  unitno=file_open(ctrl_filename)
  lineno=filesrchstr(unitno, bmap_TOXYZ_TAG, line)
  If (lineno==0)    Then
    Write(0,'(1x,2a,i4,a)') __FILE__,' : ',__LINE__, &
        ' Could not find the tag '//Trim(bmap_TOXYZ_TAG) &
        //' in the control file'
    Stop
  Endif
  !** read bmap display section from ctrlfile
  Read(unitno,*) gridsize(1),gridsize(2),gridsize(3)
  Read(unitno,'(a)') scaling
  If (Trim(toupper(stripcmnt(scaling)))=='LOGARITHMIC') Then
    logflag=.True.  
  Else If(Trim(toupper(stripcmnt(scaling)))=='ARITHMETIC') Then
    logflag=.false.
  Else
    Write(unitno,*) 'This line should be saying - LOGARITHMIC or - ARITHMETIC'
    Stop
  Endif
  Read(unitno,*) minpc, maxpc
  Read(unitno,*) xyzfile
  Read(unitno,*) descr
  Read(unitno, *) atsym
  Read(unitno,*) biasfilename
  Read(unitno, *) tk

  ! set the bmap  DISPLAY_ONLY_FLAG
  DISPLAY_ONLY_FLAG=.True.

  Call bmap_init(biasmap, eff, tk, biasfilename, scell)

  Call bmap_makeXYZ(biasmap, gridsize, logflag, minpc, maxpc, &
      xyzfile, descr, atsym )

End Program music




