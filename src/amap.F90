!------------------------------------------------------------------------------
! This module deals with analyical form of a potential map 
!
!------------------------------------------------------------------------------

Module amap

  Use defaults, Only: RDbl, RSgl, dashedline, strlen,pi, one, zero, kcalmole_kb
  Use utils, Only: allocerrdisplay, deallocerrdisplay, cleanstring, toupper
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/),vector_getnormsq, mag
  Use simcell, Only: SimCell_Params, simcell_maptouc, simcell_isortho
  Use storebase, Only: EnergyPlus, storebase_nderivs
  Use file

  Implicit None
  Save

  Private
  Public :: AMap_Info, amap_init, amap_int, amap_display

  !** The core information of an analytical map
  Type AMap_Info
    Real(kind=RDbl) :: density, sig, eps, pore_radius 
    Real(kind=RDbl) :: V1,V2
    Real(kind=RDbl) :: F1,F2
    Real(kind=RDbl) :: X0,X1
    Integer         :: debugunit
  End Type AMap_Info

Contains
  !----------------------------------------------------------------------------
  ! Initializes the analytical map
  !----------------------------------------------------------------------------
  Subroutine amap_init(map,eff,unitno)
    Type(AMap_Info), Intent(Out)   :: map
    Real(kind=RDbl), Dimension(3), Intent(In) :: eff
    Integer,Intent(in) :: unitno

    Real(kind=RDbl) :: density, sigma, epsilon, pore_radius,X0,X1
    Integer         :: error,nbrx,nbry,nbrz,nfsize
    Integer         :: debug, ios
    Character(len=strLen) :: line,tag

    Read(unitno,'(a)') line
    tag=cleanstring(line)
    If (Trim(toupper(tag))/="CYLI-Z-AMAP") Then
      Write(*,*) "first line should be: CYLI-Z-AMAP"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif

    Read(unitno,*) density
    Read(unitno,*) sigma
    Read(unitno,*) epsilon
    Read(unitno,*) pore_radius
    Read(unitno,*) X0, X1

    map%density = density
    map%sig = sigma
    map%eps = epsilon*kcalmole_kb
    map%pore_radius = pore_radius
    map%V1 = 4*pi*density*map%eps*(sigma**12)/45
    map%V2 = -2*pi*density*map%eps*(sigma**6)/3
    map%F1 = 9*map%V1
    map%F2 = 3*map%V2
    map%X0 = X0
    map%X1 = X1
    
!    open (unit=150, file='trial', status='replace', action='write', &
!      iostat=ios)
!    If (ios/=0) Then
!      Write(*,*)'Error opening file trial'
!      Stop
!    End If
!   map%debugunit = 150

  End Subroutine amap_init

  !----------------------------------------------------------------------------
  ! Calculates the POTENTIAL and, optionally the GRADIENTS from the slit pore 
  ! parametrization formula.  Applys PBCs to the passed atomic position vector.
  ! Function returns True if calculation was successful, False otherwise.
  ! Requires:  map -- analytical map structure
  !            atvec -- position vector of atom
  !            simcell -- simulation cell data structure
  !----------------------------------------------------------------------------
  Logical Function amap_int(map,results,atvec,simcell)
    Type(AMap_Info), Intent(In)            :: map
    Type(EnergyPlus), Intent(InOut)        :: results
    Type(VecType), Intent(In)              :: atvec
    Type(SimCell_Params), Intent(In)       :: simcell  
    Logical                                :: orthoflag
    Real(kind=RDbl)                        :: norm, mul_fac, Rp
    Real(kind=RDbl)                        :: RpPr, RpMr, RpPr3, RpMr3, RpPr9, RpMr9

    Type(VecType)                          :: rvect

    !Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__

    amap_int=.False.
    results%nrg=zero

    Rp=map%pore_radius

    ! rvect is the component of atvec in xy plane using the center of the cylinder as the origin 
    rvect%comp(1) = atvec%comp(1) - map%X0
    rvect%comp(2) = atvec%comp(2) - map%X1
    rvect%comp(3) = zero
    norm = mag(rvect)

    !Write(*,*) NORM, RP, MAP%SIG, rvect
    If (norm>(Rp-(map%sig/10))) Return

    RpPr= Rp+norm
    RpMr= Rp-norm
    RpPr3= (RpPr)*(RpPr)*(RpPr)
    RpMr3= (RpMr)*(RpMr)*(RpMr)
    RpPr9= (RpPr3)*(RpPr3)*(RpPr3)
    RpMr9= (RpMr3)*(RpMr3)*(RpMr3)

    !Calculate the potential
    results%nrg = map%V1*( one/RpPr9 + one/RpMr9 ) &
        + map%V2*( one/RpPr3 +  one/RpMr3)


    !** Do the analytical calculation
    Select Case(storebase_nderivs(results))
    Case (0)
      amap_int=.True.

    Case (1)
      If (norm<1.0e-9)  Then
        mul_fac=zero
      Else
        !Calculate the forces in X and Y
        mul_fac = (map%F1*( one/RpPr9/RpPr - one/RpMr9/RpMr ) &
            + map%F2*(one/RpPr3/RpPr - one/RpMr3/RpMr))/norm
      Endif



!!$( map%F1*( 1/(map%pore_radius + norm)**10 - &
!!$          1/(map%pore_radius - norm)**10 ) + map%F2*( 1/(map%pore_radius &
!!$          + norm)**4 - 1/(map%pore_radius - norm)**4 ) )/norm

      results%force(1) = mul_fac*rvect%comp(1)
      results%force(2) = mul_fac*rvect%comp(2)
      results%force(3) = zero
      amap_int=.True.

      !store in a file the position, potential and force
!      Write(map%debugunit,*) rvect, norm, results%nrg, results%force

    Case default
      Write(*,*)' Only energies and forces are available from &
          & analytical map'
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select

  End Function amap_int

  !--------------------------------------------------
  ! Display the header information
  !--------------------------------------------------
  Subroutine amap_display(a_map, nspc, optunit)
    Type(AMap_Info), Intent(in) :: a_map
    Integer, Intent(in)  :: nspc ! No. of spaces from the left column
    Integer, Optional, Intent(in) :: optunit

    Integer   :: unitno, i
    Character(len=nspc) :: spc

    spc = Repeat(" ",nspc)
    unitno=6
    If (Present(optunit))       unitno = optunit

    Write(unitno,'(2a)') spc, dashedline
    Write(unitno,'(2a)') spc,"Amap Information Section:"
    Write(unitno,'(2a,f15.2)')  spc, "density",a_map%density 
    Write(unitno,'(2a,2f15.2)')  spc, "sigma, epsilon", a_map%sig, &
        a_map%eps
    Write(unitno,'(2a,3f15.2)')  spc, "pore_radius, axis",a_map%pore_radius, &
        a_map%X0, a_map%X1
  End Subroutine amap_display

End Module amap
