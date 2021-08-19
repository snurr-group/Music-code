!-----------------------------------------------------------------------
! This contains the main program AND SUBROUTINES 
! for interpolating P,y values to match a given N1, N2
! Ultimate aim is to generate dlnf/dlnn as a function of loading
!-----------------------------------------------------------------------
Program interpol

  Use commandline, Only: commandline_init
  Use file, Only: file_open, file_settag
  Use utils, Only: genfilename,allocErrDisplay,int2str, stripcmnt, &
      split, toreal, toint
  Use defaults, Only: dashedline, strLen, d_ctrl_file, RDbl, one, zero
  Use params, Only: ParamsSet, params_init , params_guess , &
      params_solve, params_eval, params_copy, params_interpol

  Implicit None

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action

  ! params for both compounds
  Type(ParamsSet), Dimension(:),Pointer :: p1,p2
  Type(ParamsSet) :: trial_iso1, trial_iso2

  Character(len=strLen)   :: sorb1,sorb2
  Character(len=strLen), Dimension(strlen)     :: fields 
  Character(len=strLen), Dimension(10)     :: ffields 
  Character(len=2*strLen) :: text ,outfile
  Character(len=strLen)      :: infile1, infile2, string

  Integer                 :: nfields, niso1, niso2, error
  Integer                 :: i, iy, i5, iN2, nN2, j, k, ctrlunit, isounit 
  Integer                 :: outunit, ydivN, minN2, maxN2, incrN2
  Real(kind=RDbl), Dimension(:), Pointer :: gas_y1, gas_y2
  Real(kind=RDbl) :: f1_01, f1_10, f1_11, f1_12, f1_21
  Real(kind=RDbl) :: f2_01, f2_10, f2_11, f2_12, f2_21
  Integer                 :: smally_index, largey_index
  Real(kind=RDbl) :: P, ntolerance, currentN2, currentN1, deltaN, dN1, dN2
  Real(kind=RDbl) :: N1, N2, y, dY, low, high, tempP,Nactual,Nerror
  Logical :: found_soln, found_fullsoln



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

  ! get the name of the output file
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  outfile=Trim(Adjustl(fields(2)))
  Write(*,*) "outfile : ",Trim(outfile)

  ! get the name of the molecules
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  sorb1=Trim(Adjustl(fields(2)))
  Write(*,*) "First sorbate : ",Trim(sorb1)
  Write(*,*) "******** Make sure this sorbate does not show a maximum *****"
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  infile1=Trim(Adjustl(fields(2)))
  Write(*,*) "infile for sorbate 1 : ",Trim(infile1)

  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  nfields=split(fields(2),ffields,",")
  currentN1=toint(ffields(1))
  Write(*,*) "N1 for sorbate 1 : ",currentN1


  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  sorb2=Trim(Adjustl(fields(2)))
  Write(*,*) "Second sorbate : ",Trim(sorb1)

  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  infile2=Trim(Adjustl(fields(2)))
  Write(*,*) "infile for sorbate 2 : ",Trim(infile1)

  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  nfields=split(fields(2),ffields,",")
  minN2=toint(ffields(1))
  maxN2=toint(ffields(2))
  incrN2=toint(ffields(3))
  Write(*,*) "minN1 and maxN1 and increment for sorbate 2 : ",&
      minN2, maxN2, incrN2

  ! deltaN for numerical differentiation wrt N
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  deltaN=toreal(fields(2))
  Write(*,*) "deltaN used for numer. differentiation : ", deltaN

  ! divisions in y while checking for the correct y
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  ydivN=toint(fields(2))
  Write(*,*) "divisions in y while checking for the correct y", ydivN

  ! no of isotherms
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  ntolerance=toreal(fields(2))
  Write(*,*) "Allowed error in N1, N2", ntolerance


  ! open first sorbates file and read isotherm-parameters
  isounit=file_open(infile1)
  Write(*,*) "Reading isotherm params file : ", infile1
  Read(isounit,'(//a)') string
  nfields=split(string,fields,":")
  niso1=toint(fields(2))
  Allocate(p1(niso1),gas_y1(niso1), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
  Do i=1,niso1
    Call params_init(p1(i),isounit)
    gas_y1(i)=p1(i)%y
  Enddo
  close(isounit)

  ! open first sorbates file and read isotherm-parameters
  isounit=file_open(infile2)
  Write(*,*) "Reading isotherm params file : ", infile2
  Read(isounit,'(//a)') string
  nfields=split(string,fields,":")
  niso2=toint(fields(2))
  If (niso2/=niso1) Then
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Endif
  Allocate(p2(niso2),gas_y2(niso2), STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
  Do i=1,niso2
    Call params_init(p2(i),isounit)
    gas_y2(i)=1-p2(i)%y
    If( Abs(p1(i)%y-p2(i)%y)> 0.00001) Then
      Write(*,*) "Both files should have same value of y(gas comp. &
          &of first sorbate in binary-gcmc)"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
  Enddo


  !open outfile
  outunit=file_open(outfile)
  Write(outunit,'(3a)') "# Sorbates : " , Trim(sorb1), Trim(sorb2)
  Write(outunit,'(a)')  "# f1/c2 actually refers to dlnf1/ dlnc1 "
  Write(outunit,'(a)')  "# N1   N2    P       y    f1       f2        f1/c1   f1/c2   f2/c2   f2/c1"



  ! series of loops

  ! over different N2 values
  nN2=(maxN2-minN2)/incrN2
  Do iN2=0,nN2
    currentN2=minN2+iN2*incrN2
    If (currentN2>maxN2) Exit

    ! evaluate at five different sets of N for numerical differenciation
    !    01             N2 changes along  rows, N1 along columns
    ! 10 11 12 
    !    21
    Do i5=1,5
      dN1=0
      dN2=0
      Select Case(i5)
      Case(1) ! 11 
      Case(2) ! 01
        dN2=-1
      Case(3) ! 10
        dN1=-1
      Case(4) ! 12
        dN1=1
      Case(5) ! 21
        dN2=1
      End Select

      N1=currentN1+dN1*deltaN
      N2=currentN2+dN2*deltaN
      If ((currentN1<0).Or.(currentN2<0)) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      Write(*,'(a,t40,2f7.2)') "doing evaluations at N1, N2 :", N1, N2
      dY=(one/ydivN)
      found_fullsoln=.False.
      Do iy=1,ydivN-1

        y=dY*iy
        smally_index=1
        largey_index=0
        found_fullsoln=.False.
        ! find which interval we lie in
        Do j=1,niso1
          ! check whether it is too close
          If (Abs(y-gas_y1(j))<dY/100) Then
            smally_index=j
            largey_index=j
            Exit
          Endif
          If ((y>gas_y1(j)).And.(y<gas_y1(j+1))) Then
            smally_index=j
            largey_index=j+1
            Exit
          Endif
        End Do
        Write(*,*) "--------------------------------------"
        Write(*,*) "solving for : ",y, smally_index, largey_index

        !Interpolate
        If (smally_index==largey_index) Then
          Call params_copy(trial_iso1,p1(smally_index))
          Call params_copy(trial_iso2,p2(smally_index))
        Elseif (smally_index==(largey_index-1)) Then
          Call params_interpol(trial_iso1,p1(smally_index),p1(largey_index),y)
          Call params_interpol(trial_iso2,p2(smally_index),p2(largey_index),y)
        Else
          Write(*,*) "could not find interval"
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
!!$        Write(*,*) smally_index, largey_index
!!$        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$        Write(*,*) "using params : ", trial_iso1%K, trial_iso1%Lmax
!!$        Write(*,*) "using params : ", trial_iso1%c, trial_iso1%b
!!$        Write(*,*)  p1(smally_index)%K,p1(largey_index)%K
!!$        Write(*,*)  p1(smally_index)%Lmax,p1(largey_index)%Lmax
!!$        Write(*,*)  p1(smally_index)%c,p1(largey_index)%c
!!$        Write(*,*)  p1(smally_index)%b,p1(largey_index)%b
!!$        Write(*,*)  p1(smally_index)%y,p1(largey_index)%y

        low=1.0e-8
        high=-one
        ! Guess a solution using Langmuir
        tempP=params_guess(trial_iso1,N1)
        If (1+trial_iso1%b*tempP<0.0000001) Then
          found_fullsoln=.False.
          Cycle
        endif
        !!        Write(*,*) "Guess P : ", tempP
        Nactual=params_eval(trial_iso1,tempP)
        If (Nactual>N1) Then
          high=tempP
        Else
          Do k=1,30
            tempP=2*tempP
            If (1+trial_iso1%b*tempP<0.0000001) Then
              Exit
            Endif
            Nactual=params_eval(trial_iso1,tempP)
            If (Nactual>N1) Then
              high=tempP
              Exit
            Endif
          End Do
!!          Write(*,*) tempP, trial_iso1%b*tempP
          !            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          !            Write(*,*) "looks like no soltion for this y"
          !            Stop
          If ((tempP>1.0e8).or.((trial_iso1%b*tempP<(0.9999999)))) Then
            found_fullsoln=.False.
            Cycle
          Endif

          If (high<zero) then
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
          endif

        Endif
        !        Write(*,*) "high and low :", high, low

        !Use high and low and bisect till we reach a solution
        found_soln=.False.
        Do k=1,2000
          tempP=(low+high)/2
!!          Write(*,*) tempP, trial_iso1%b*tempP
          Nactual=params_eval(trial_iso1,tempP)
          Nerror=Abs(Nactual-N1)
          If (Nerror<ntolerance) Then
            found_soln=.True.
            exit
          Endif
          If (Nactual>N1) Then
            high=tempP
          else
            low=tempP
          endif
        End Do

        If (.Not.found_soln) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
        !    Write(*,*) "Nerror ", Nerror
        ! find N2
        Nactual=params_eval(trial_iso2,tempP)
        Nerror=Abs(Nactual-N2)
        !        Write(*,*) "Nerror 2", Nerror
        If (Nerror<ntolerance) Then
          found_fullsoln=.True.
          exit
        Endif

      End Do
      If (found_fullsoln) Then
        Write(*,'(a,t40,2e12.3)') "solution, f1, f2 :", tempP*y, tempP*(1-y)
        ! check soln
        Call params_interpol(trial_iso1,p1(smally_index),p1(largey_index),y)
        Call params_interpol(trial_iso2,p2(smally_index),p2(largey_index),y)
        !      Write (*,'(e12.3,2f12.5,e12.3)') trial_iso1%K,  trial_iso1%Lmax,  trial_iso1%c,  trial_iso1%b
        !      Write (*,'(e12.3,2f12.5,e12.3)') trial_iso2%K,  trial_iso2%Lmax,  trial_iso2%c,  trial_iso2%b
        Write(*,'(a,t40,2f17.5)') "Evaluated N1, N2", &
            params_eval(trial_iso1,tempP), params_eval(trial_iso2,tempP)
        Write(*,'(a,t40,i3,a,i3)')"Interpolation was done betwenn:" ,&
            smally_index, " and ", largey_index
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        stop
      endif

      Select Case(i5)
      Case(1) ! 11
        f1_11= tempP*y
        f2_11= tempP*(1-y)
      Case(2) ! 01
        f1_01= tempP*y
        f2_01= tempP*(1-y)
      Case(3) ! 10
        f1_10= tempP*y
        f2_10= tempP*(1-y)
      Case(4) ! 12
        f1_12= tempP*y
        f2_12= tempP*(1-y)
      Case(5) ! 21
        f1_21= tempP*y
        f2_21= tempP*(1-y)
      End Select


    End do

    Write(*,'(a,2e12.3,a)') "fugacity     : ", f1_11, f2_11, "kPa"
    Write(*,'(a,1e12.3,f12.5)') "P, y     : ", f1_11+f2_11, f1_11/(f1_11+f2_11)
    Write(*,'(a,2e12.3,a)') "dlnf1 / dlnc1: ", &
        currentN1*(log(f1_12/f1_11))/deltaN, currentN1*(log(f1_11/f1_10))/deltaN
    Write(*,'(a,2e12.3,a)') "dlnf1 / dlnc2: ", &
        currentN2*(log(f1_21/f1_11))/deltaN, currentN2*(log(f1_11/f1_01))/deltaN
    Write(*,'(a,2e12.3,a)') "dlnf2 / dlnc1: ", &
        currentN1*(log(f2_12/f2_11))/deltaN, currentN1*(log(f2_11/f2_10))/deltaN
    Write(*,'(a,2e12.3,a)') "dlnf2 / dlnc2: ", &
        currentN2*(log(f2_21/f2_11))/deltaN, currentN2*(log(f2_11/f2_01))/deltaN
    Write(*,*) "----------------------"
    Write(*,*) 

    Write(outunit,'(2f5.1,e9.2,f5.2,2e9.2,4f8.3)') currentN1, currentN2, &
        f1_11+f2_11, f1_11/(f1_11+f2_11), f1_11, f2_11,         &
        currentN1*(Log(f1_12/f1_10))/(2*deltaN), &
        currentN2*(Log(f1_21/f1_01))/(2*deltaN), &
        currentN2*(Log(f2_21/f2_01))/(2*deltaN), &
        currentN1*(Log(f2_12/f2_10))/(2*deltaN)
  End do
End Program interpol


