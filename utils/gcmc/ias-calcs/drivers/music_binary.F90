!-----------------------------------------------------------------------
! This contains the main program AND SUBROUTINES 
! for estimating the binary adsorption Isotherm from single component 
! two isotherms for IAS theory
! ctrlfile format :
! Toluene, iso.tol.dat     # sorb-1, isotherm file(see format in isotherm.F90
! Xylene , iso.xyl.dat     # sorb-1, isotherm file(see format in isotherm.F90
! binary.dat               # output file
! 3                        # number of pressure points
! 0.1     1.9              # partial pressure toluene first the n xylene
! 0.3     1.6
! 0.9     1.8 
!
!-----------------------------------------------------------------------
Program IAS_binary

  Use commandline, Only: commandline_init
  Use file, Only: file_open, file_settag
  Use utils, Only: genfilename,allocErrDisplay,int2str, stripcmnt, split
  Use defaults, Only: dashedline, strLen, d_ctrl_file, RDbl
  Use isotherm, Only: isotherm_init, Isotherm_Info
  Use iasbinary, Only: iasbinary_interpolate,iasbinary_SolveGivenSigP
  Implicit None

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action

  Type(Isotherm_Info), Dimension(2) :: isoT

  Character(len=strLen), Dimension(2)          :: sorbname
  Character(len=strLen), Dimension(strlen)     :: fields 
  Character(len=2*strLen) :: inputfile, text ,outputfile

  Integer                 :: i, j, ctrlunit, nfields, npts, error,top, bottom
  Integer                 :: outunit
  Real(kind=RDbl), Dimension(:,:), Pointer :: p_i, y, x
  Real(kind=RDbl), Dimension(:), Pointer   :: Pmix, n
  Real(kind=RDbl) :: P, sigma1, sigma2, err, low_sig, high_sig, trial_sig
  Real(kind=RDbl) :: x1, y1, x2, y2, n_tot
  Logical outOfRange
  Integer , Parameter :: MAX_SEARCH = 100000000
  Real(kind=RDbl), Parameter :: MOLE_FRAC_TOLERANCE =1.00e-6


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

  ! get the name of the molecule and its single component isotherm
  Do i=1,2
    Read(ctrlunit,'(a)') text
    text=stripcmnt(text)
    nfields=split(text,fields,",")
    If (nfields/=2) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif

    sorbname(i)=Trim(adjustl(fields(1)))
    inputfile=Trim(fields(2))
    Write(*,*) "Reading isotherm from "//Trim(inputfile)//&
        " for sorbate "//Trim(sorbname(i))
    Call isotherm_init(isoT(i),inputfile)
    If (Trim(sorbname(i))/=Trim(isoT(i)%sorbname)) Then
      Write(*,*) "Wrong sorbname"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
  End Do
  Read(ctrlunit,*) text
  outputfile=Trim(Adjustl(stripcmnt(text)))
  Write(*,*) " Reading the pressure list from ctrlfile"
  Read(ctrlunit,*) npts
  Allocate(p_i(2,npts), y(2,npts), Pmix(npts), x(2,npts), n(npts), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
!!$  Allocate(p1(npts), p2(npts), y1(npts), y2(npts), P(npts), x1(npts), &
!!$      x2(npts), n(npts), STAT=error)
!!$  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
  !** Given variables : p_i : partial pressures in gas phase
  !** derived from above : y, Pmix : gas mole fractions, Total pressure
  !** To be calculated : x, n : mole fractions in ads. phase, totalloading
  Do i=1,npts
    Read(ctrlunit,*) p_i(1,i), p_i(2,i)
    Pmix(i)=p_i(1,i)+p_i(2,i)
    y(1,i)=p_i(1,i)/Pmix(i)
    y(2,i)=1-y(1,i)
  End Do
  Write(*,*) " Finished reading the pressure list from ctrlfile "






  !------------------------------------------------------------------
  !  for each point (given T, P, y_i) calculate the binary adsn isotherm
  !  Refer to Myer's IAS paper AIChe 1965 Fig.1 for details
  !                      
  !------------------------------------------------------------------
  Do i=1,npts

    P= Pmix(i) ! total pressure

    ! value at which vertical line cuts the sigma vs Pressure 
    ! curves at pressure=P. Note sigma=surface pressure eqvt to pi/RT
    sigma1=iasbinary_interpolate(P, isoT(1)%P, isoT(1)%sigma ,outOfRange)
    sigma2=iasbinary_interpolate(P, isoT(2)%P, isoT(2)%sigma ,outOfRange)


    !** decide which is the top curve and which is the bottom curve
    If (sigma1>sigma2) Then
      top=1
      bottom=2
    Else
      top=2
      bottom=1
    Endif

    err=Max(1.00_RDbl,1/MOLE_FRAC_TOLERANCE)      ! some large value
    high_sig= iasbinary_interpolate(P, isoT(top)%P, isoT(top)%sigma ,&
        outOfRange)
    low_sig= iasbinary_interpolate(P, isoT(bottom)%P, isoT(bottom)%sigma, outOfRange )

    ! do binary search for the correct sigma
    Do j=1,MAX_SEARCH
      If (j==MAX_SEARCH) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*)"You have a convergence problem"
        Stop
      Endif
      trial_sig=(high_sig+low_sig)/2
      Call iasbinary_SolveGivenSigP( trial_sig, P, isoT(top), isoT(bottom), &
          x1, x2, y1, y2, n_tot)
      err=(Abs(y1-y(top,i)))+(Abs(y2-y(bottom,i)))
      If (err< MOLE_FRAC_TOLERANCE) Then
        x(top,i)=x1
        x(bottom,i)=x2
        n(i)= n_tot
        Exit
      Endif
      ! fix  interval for next trial
      If (y1>y(top,i)) Then
        high_sig=trial_sig
      Else
        low_sig=trial_sig
      Endif

    End Do ! end of binary seach loop, j

  End Do ! end of data points loop, i





  !------------------------------------------------------------------
  !  Write the outputs to output file 
  !  
  !                      
  !------------------------------------------------------------------
  Write(*,*) outputfile
  outunit=file_open(outputfile)

  Write(outunit,'(a)')       "% Binary adsorption from IAS calcs"
  Write(outunit,'(a)')       "% Molecule-1  : "//Trim(sorbname(1))
  Write(outunit,'(a)')       "% Molecule-2  : "//Trim(sorbname(2))
  Write(outunit,'(a,f10.5)') "% Temperature : ",isoT(1)%Temp
  Write(outunit,'(a)')       "% Details     : "//Trim(isoT(1)%description) 
  Write(outunit,'(a)')       "% Pressure-1 Pressure-2 Loading-1 Loading-2 "
  Write(outunit,'(a)')       "% kpa        kpa        mol/uc    mol/uc    "
  Do i=1,npts
    Write(outunit,'(3x,4f10.5)') p_i(1,i), p_i(2,i), n(i)*x(1,i), n(i)*x(2,i)
  End Do
  close(outunit)
  Write(*,*) "Finished writing output to : "//Trim(outputfile)



End Program IAS_BINARY



