!------------------------------------------------------------------------------
! This module initializes the variables required by EMD for calcuation of 
! L coefficients. The COM velocities are periodically written to a binary file
! which can be later post-processed
!------------------------------------------------------------------------------
Module emd
  Use config, Only: AtMolCoords, config_getnmoles, config_isfixed, &
      config_getSpcCOMVel
  Use defaults, Only: RDbl, strLen, STATS_BLOCKSIZE, one, zero, scalef, &
      MAX_SORBS
  Use file, Only: file_open
  Use molecules, Only: molecules_getnatoms, molecules_AtomMass, &
      molecules_getatype, molecules_getnsorbs, molecules_getmass, &
      molecules_gettype, molecules_name
  Use simcell, Only: SimCell_Params, simcell_getvolume
  Use utils, Only: stripcmnt, split, tolower, toupper, filesrchstr, &
      allocErrDisplay, DeallocErrDisplay, int2str, real2str, toreal, &
      cleanstring
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag, vector_zerovec


  Implicit None
  Save

  Private
  Public :: emd_init, emd_display, emd_output, EMD_IS_ON, EMDINFO, EMDParams

  Type EMDInfo
    Character(len=strLen)           :: filename
    Integer                         :: unitno
    Real(kind=RDbl)                 :: tstep
    Type(VecType) ,Dimension(:), Pointer                 :: vels_sum
  End Type EMDInfo

  Type(EMDInfo) , Pointer :: EMDParams
  Logical                 :: EMD_IS_ON=.False.

  Character(len=strLen) :: EMD_tag= " EMD Details "

Contains

  !----------------------------------------------------------------------------
  ! checks whether EMD is required, Initializes EMD
  ! Requires:  ctrlfile -- main control file 
  ! This works only with constant number of molecules
  !----------------------------------------------------------------------------
  Subroutine emd_init(ctrlfile,species,scell,Tk)
    Character(*), Intent(in) :: ctrlfile
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(SimCell_Params), Intent(In)               :: scell
    Real(kind=RDbl),Intent(in) :: Tk
    Integer               :: spc, error, nsorbs, useEMD, nmoles, noutsorbs
    Integer               :: unitno, outunit
    Character(len=strLen) :: text, name
    Real(kind=RDbl) :: volume
    unitno=file_open(ctrlfile)
    useEMD= filesrchstr(unitno,EMD_tag,text,.True.)

    If (useEMD==0) Then
      EMD_IS_ON=.False.
      Nullify(EMDParams)
      Return
    else
      EMD_IS_ON=.True.
      Allocate(EMDParams,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Write(*,*) "Initializing the EMD details by reading from the ctrlfile"
      Read(unitno,*) EMDParams%filename
      EMDParams%filename=cleanstring(EMDParams%filename)
      Read(unitno,*) EMDParams%tstep
    Endif

    ! Write the header of the outout file
    EMDParams%unitno=file_open(EMDParams%filename,102)
    outunit=EMDParams%unitno
    volume=simcell_getvolume(scell)
    Write(outunit) volume, Tk
    nsorbs=molecules_getnsorbs()
    ! ** count non fixed sorbs
    noutsorbs=0
    Do spc=1,nsorbs
      If (config_isfixed(species(spc))) Cycle
      noutsorbs=    noutsorbs+1
    End Do
    Write(outunit) noutsorbs
    Do spc=1,nsorbs
      If (config_isfixed(species(spc))) Cycle
      nmoles=config_getnmoles(species(spc))
      name=molecules_name(spc)
      Write(outunit) spc,nmoles,name
      noutsorbs=    noutsorbs+1
    End Do
    Allocate(EMDParams%vels_sum(noutsorbs),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    EMDParams%vels_sum=vector_zerovec()
  End Subroutine emd_init

  !------------------------------------------------
  ! periodically write velocities to a file
  !------------------------------------------------
  Subroutine emd_output(species,time)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Real(kind=RDbl) , Intent(in) :: time

    Type(VecType) :: comvel
    Integer :: nmoles, outspc, spc, nspc, outunit

    Logical, Save :: first_time = .True.
    Real(kind=RDbl) , Save :: last_time = zero
    ! ** 
    If (first_time) Then
      first_time=.False.
      last_time=time
    Endif
    !** we use 0.9999 to avoid comparison of two reals using ==
    !** 
    If ((time-last_time)>(0.999999*EMDParams%tstep)) Then
      last_time=time
      outspc=0
      nspc=molecules_getnsorbs()
      Do spc=1,nspc
        If (config_isfixed(species(spc))) cycle
        outspc=outspc+1
        nmoles=config_getnmoles(species,spc)
        comvel=config_getSpcCOMVel(species,spc)
        EMDParams%vels_sum(outspc)=comvel*nmoles
      End Do
      outunit=EMDParams%unitno
      Write(outunit) time, EMDParams%vels_sum(1:outspc)
    Endif

  End Subroutine emd_output

  !----------------------------------------------------------------------------
  ! Displays information about EMD parameters
  !----------------------------------------------------------------------------
  Subroutine emd_display(unitno,nskip)
    Integer, Intent(In)              :: unitno
    Integer, Intent(In)              :: nskip

    Integer               :: dUnit
    Character(len=nskip)  :: blank

    dUnit = unitno
    blank=Repeat(" ", nskip)
    Write(dUnit,'(2a)')blank, "EMD Details : "
    Write(dUnit,'(2a,t40,2a)')blank, &
        "Ouput Filename ",":  ",Trim(EMDParams%filename)
    Write(dUnit,'(2a,t40,a,f10.4,a)')blank, &
        "Ouput Timestep ",":  ",EMDParams%tstep, " pico second"
    Write(dUnit,*) 
  End Subroutine emd_display

End Module emd

