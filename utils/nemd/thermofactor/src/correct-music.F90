!-----------------------------------------------------------------------
! This contains the main program AND SUBROUTINES 
! for fitting an adsorption Isotherm
! Each isotherm is for a particular gas composition
! Isotherms contain : partial pressure, loading,total pressure
!-----------------------------------------------------------------------
Program fitting

  Use commandline, Only: commandline_init
  Use file, Only: file_open, file_settag
  Use utils, Only: genfilename,allocErrDisplay,int2str, stripcmnt, &
      split, toreal, toint
  Use defaults, Only: dashedline, strLen, d_ctrl_file, RDbl
  Use isotherm, Only: Isot, isotherm_init , isotherm_guess , &
      isotherm_MSerror, isotherm_fit

  Implicit None

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action

  Type(Isot), Dimension(:),Pointer :: isos

  Character(len=strLen)   :: sorbname
  Character(len=strLen), Dimension(strlen)     :: fields 
  Character(len=2*strLen) :: text ,outfile
  Character(len=strLen), Dimension(:),Pointer     :: isofile

  Integer                 :: i, j, ctrlunit, nfields, nisot, error
  Integer                 :: outunit
  Real(kind=RDbl), Dimension(:), Pointer :: gas_y

  Real(kind=RDbl) :: P, tolerance, mserror
  Real(kind=RDbl) :: K, Lmax, a, Pmax

  !-------------------------------------------------------------
  !
  ! Do basic initialization, based on command line input
  !
  !-------------------------------------------------------------
  Call commandline_init(ctrl_filename,action)

  If (action=="WriteSampleCtrlfile") Then
    Write(0,'(a)') dashedline
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Elseif(action/="DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Endif

  Call file_settag(ctrl_filename, d_ctrl_file)
  ctrlunit=file_open(ctrl_filename,110)

  !------------------------------------------------------------------
  !
  !                      Read the ctrlfile
  !
  !------------------------------------------------------------------
  Write(*,*) "Read ctrlfile :"//Trim(ctrl_filename)
  ! get the name of the molecule 
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  sorbname=Trim(Adjustl(fields(2)))

  ! get the name of the ouput file
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  outfile=Trim(Adjustl(fields(2)))

  ! least square error tolerance
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  tolerance=toreal(fields(2))

  ! no of isotherms
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  nisoT=toint(fields(2))
  Allocate(isofile(nisoT), gas_y(nisoT),isos(nisoT),STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

  ! blankline
  Read(ctrlunit,*) text

  Write(*,*) "Reading list of isotherms from  ctrlfile :"

  ! get the name of the file and corresponding gas phase composition 
  ! of this compound
  Do i=1,nisoT
    Read(ctrlunit,'(a)') text
    nfields=split(text,fields,":")
    isofile(i)=Trim(Adjustl(fields(1)))
    gas_y(i)=toreal(fields(2))
  Enddo

  Write(*,*) "Finished reading from  ctrlfile :"

  ! write the header of the output file
  outunit=file_open(outfile)
  Write(outunit,*) "%Molecule Name        : ",Trim(sorbname)  
  Write(outunit,*) "%Filename   y     K     Lmax,    c/P/P, error"



  Do i=1,nisoT
    Write(*,*) "Reading from  isotherm file : "//Trim(isofile(i))
    Call isotherm_init(isos(i),isofile(i))

    Write(*,*) "Guessing the parameters"
    Call isotherm_guess(isos(i),K,Lmax,a,PMax)
    Write(*,*) "Guesses are : ", K, Lmax, a, PMax
    mserror=isotherm_MSerror(isos(i),K, Lmax, a, PMax)
    Write(*,*) "Error with this Guess : ", mserror
    ! Eqn is N = K*P / ( 1+ (K/Lmax)*P + a * (P/PMax)**2)
    ! where K, Lmax and a are adjustable
    ! at the end K, Lmax and a/Pmax^2 will be reported

    Call isotherm_fit(isos(i), K, Lmax, a, PMax, tolerance)
    mserror=isotherm_MSerror(isos(i),K, Lmax, a, PMax)
    ! write to outpu tfile
    Write(outunit,'(a,t20,f10.4,1x,e12.3,1x,f10.5,1x,2e12.3)') &
        Trim(isofile(i)), gas_y(i), K, Lmax, a/Pmax/Pmax, mserror 

  Enddo


!!$
!!$
!!$
!!$
!!$
!!$  !------------------------------------------------------------------
!!$  !  Write the outputs to output file 
!!$  !  
!!$  !                      
!!$  !------------------------------------------------------------------
!!$  Write(*,*) outputfile
!!$  outunit=file_open(outputfile)
!!$
!!$  Write(outunit,'(a)')       "% Binary adsorption from IAS calcs"
!!$  Write(outunit,'(a)')       "% Molecule-1  : "//Trim(sorbname(1))
!!$  Write(outunit,'(a)')       "% Molecule-2  : "//Trim(sorbname(2))
!!$  Write(outunit,'(a,f10.5)') "% Temperature : ",isoT(1)%Temp
!!$  Write(outunit,'(a)')       "% Details     : "//Trim(isoT(1)%description) 
!!$  Write(outunit,'(a)')       "% Pressure-1 Pressure-2 Loading-1 Loading-2 "
!!$  Write(outunit,'(a)')       "% kpa        kpa        mol/uc    mol/uc    "
!!$  Do i=1,npts
!!$    Write(outunit,'(3x,4f10.5)') p_i(1,i), p_i(2,i), n(i)*x(1,i), n(i)*x(2,i)
!!$  End Do
!!$  close(outunit)
!!$  Write(*,*) "Finished writing output to : "//Trim(outputfile)
!!$
!!$

End Program fitting



