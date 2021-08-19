! format of in.dat for heptane  nextline sould be first line of in.dat 
! 450, 540.2, 2.736, 0.351, 1.0
!
!
!-----------------------------------------------
! does pr calculations of VLE for one compound
!-----------------------------------------------
Program pr
  Implicit None
  ! see sandle chapter 4.7 for meaning of these symbols
  Real :: Tc, Pc, omega, T, P, guessP, kappa, alpha, a, b, dZ, mwt
  Real :: AA, BB, C0, C1, C2, newf, prevf, Z, ZV, ZL, lnf_P, fV, fL, err
  Integer :: rcount, itercount, i
  Real,dimension(3) :: roots=0.00
  Real, Parameter  :: R = 8.314, R2=1.4142
  Integer, Parameter  :: NZdivs=1000000
  Write(*,*) "going to read inputfile in.dat"
  Write(*,*) "make sure you have the following there T, Tc, Pc, omega, guessP"
  Write(*,*) "units K, Mpa"
  Open(unit=11,file='in.dat')
  Read(11,*) T, Tc, Pc, omega, guessP

  ! convert to pascal
  Pc=Pc*1e6
  guessP=guessP*1e6

  ! calculate a,b
  kappa=  0.37464 + 1.54226*omega - 0.26992*omega*omega
  alpha= ( 1 + kappa * (1 - sqrt(T/Tc)) )**2
  a = 0.45724 * (R * Tc)**2 / Pc * alpha
  b = 0.0778 * R * Tc / Pc

  P=guessP
  itercount=0
  Do
    Write(*,*) "iteration trial pressure", itercount, P
    itercount=itercount+1
    ! AA=A in sandler , BB=B in sandler
    ! section 5.5
    AA=a*P/R/R/T/T
    BB=P*b/R/T

    ! solve PR eqn for compressibility
    !
    ! Searching for zero using Newton-Rhapson
    ! Z**3 + C2*Z**2 + C1*Z + C0 = 0
    ! Z: compressibility facor
    !

    C2 = BB - 1
    C1 = AA - 3*BB*BB - 2*BB
    C0 = -AA*BB + BB*BB + BB*BB*BB
    prevf=C0
    rcount=0
    dZ=1.00000000/NZdivs
    !    Open(unit=12,file="err")
    !      Write(12,*) 0.000,C0
    Do i=1,NZdivs
      Z=i*dZ
      newf=Z**3 + C2*Z**2 + C1*Z + C0
      !      Write(12,*) Z,newf
      If (Newf/prevf<0.000000000)Then
        rcount=rcount+1
        If (rcount>3) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
        roots(rcount)=Z
      Endif
      prevf=newf
    End Do

    If (rcount/=3) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    ! they should be in ascending order
    ZV=roots(3)
    ZL=roots(1)

    ! calculate gas fugacity, 5.4-14a
    lnf_P= (ZV-1) - Log(ZV-BB) -(AA/2/BB/R2) * &
        (  Log( (ZV+(1+R2)*BB)/ (ZV+(1-R2)*BB) )  )
    fV=Exp(lnf_P)*P

    ! calculate liq fugacity, 5.4-14b
    lnf_P= (ZL-1) - Log(ZL-BB) -(AA/2/BB/R2) * &
        (  Log( (ZL+(1+R2)*BB)/ (ZL+(1-R2)*BB) )  )
    fL=Exp(lnf_P)*P

    err=Abs(fL/fV-1)
    If (err>0.0001) Then
      P=P*fL/fV
    Else
      Write(*,*) "found solution"
      Exit
    Endif
  End Do
  Write(*,*) "equilibrium vapor presure", P
  Write(*,*) "volume of liquid", ZL*R*T/P
  Write(*,*) "volume of gas", ZV*R*T/P

  Write(*,*) "density liq mol/cm^3", 1e-6/(ZL*R*T/P)
  Write(*,*) "density gas mol/cm^3", 1e-6/(ZV*R*T/P)

  Write(*,*) "enter mol wt if you would like the density in g/cm^3, &
      &(-1 to exit)"
  Read(*,*) mwt
  If (mwt<0) Stop

  Write(*,*) "density liq g/cm^3", mwt*1e-6/(ZL*R*T/P)
  Write(*,*) "density gas g/cm^3", mwt*1e-6/(ZV*R*T/P)

End Program pr
