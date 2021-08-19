!------------------------------------------------------------------------------
! This module contains the interpolation routines that read potentials
! and, optionally, forces from a pre-tabulated map.  Using such a map
! greatly speeds simulation time if the simulation cell contains a
! significant number of fixed atoms.  The interpolation routine (hermite)
! was inherited from R.L. June and is based on other published routines.
! Interfaces are provided here to distinguish between interpolations that
! do and do not also require forces.
!------------------------------------------------------------------------------

Module interpolate

  Use defaults, Only: RDbl, RSgl
  Use utils, Only: allocerrdisplay, deallocerrdisplay
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use simcell, Only: SimCell_Params, simcell_maptouc, simcell_isortho
  Use mapheader, Only: Header_Info, mapheader_getIndices, mapheader_sizes
  Use storebase, Only: EnergyPlus, storebase_nderivs

  Implicit None
  Save

  Private
  Public :: Map_Core_Info, interpolate_initcore, interpolate_int, &
      interpolate_writecore, interpolate_readcore, interpolate_npts, &
      interpolate_clean

  !** The core information of any map
  Type Map_Core_Info
    Integer                                :: npts
    Real(kind=RSgl), Dimension(:), Pointer :: r, rx, ry, rz, rxy
    Real(kind=RSgl), Dimension(:), Pointer :: rxz, ryz, rxyz
    !** Note that the convention for storing the pointers is z,y,x
    Integer, Dimension(:,:,:), Pointer     :: ptr
  End Type Map_Core_Info

Contains
  !----------------------------------------------------------------------------
  ! Initializes the core of any map, only sizes arrays given information in
  ! the passed header
  ! Requires:  map -- core information to initialize
  !----------------------------------------------------------------------------
  Subroutine interpolate_initcore(map,header)
    Type(Map_Core_Info), Intent(Out)   :: map
    Type(Header_Info), Intent(In)      :: header

    Integer          :: error,nbrx,nbry,nbrz,nfsize

    !** Get the size information from the header
    Call mapheader_sizes(header,nbrx,nbry,nbrz,nfsize)
    map%npts = nfsize

    !** Allocate the pointer array
    Allocate(map%ptr(nbrz,nbry,nbrx), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%ptr')   

    !** Allocate potential 
    Allocate(map%r(nfsize),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%r')   

    !** Force, x direction
    Allocate(map%rx(nfsize),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rx')   

    !** Force, y direction
    Allocate(map%ry(nfsize),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%ry')   

    !** Force, z direction
    Allocate(map%rz(nfsize),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rz')   

    !** Second order derivative dxdy
    Allocate(map%rxy(nfsize),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rxy')   

    !** Second order derivative dxdz
    Allocate(map%rxz(nfsize),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rxz')   

    !** Second order derivative dydz
    Allocate(map%ryz(nfsize),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%ryz')   

    !** Third order derivative dxdydz
    Allocate(map%rxyz(nfsize),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rxyz')   

  End Subroutine interpolate_initcore

  !----------------------------------------------------------------------------
  ! Provides an interface to extract the POTENTIAL and, optionally the 
  ! GRADIENTS from the hermite interpolation procedure.  Applys PBCs to the 
  ! passed atomic position vector.  Function returns True if interpolation 
  ! was successful, False otherwise.
  ! Requires:  map -- potential map structure
  !            header -- header information about map
  !            pot -- output potential
  !            atvec -- position vector of atom
  !            simcell -- simulation cell data structure
  !            factor -- optional factor to multiply results by
  !----------------------------------------------------------------------------
  Logical Function interpolate_int(map,header,results,atvec,simcell,factor)
    Type(Map_Core_Info), Intent(In)        :: map
    Type(Header_Info), Intent(In)          :: header
    Type(EnergyPlus), Intent(InOut)        :: results
    Type(VecType), Intent(In)              :: atvec
    Type(SimCell_Params), Intent(In)       :: simcell  
    Real(kind=RDbl), Intent(In), Optional  :: factor

    Logical                     :: orthoflag
    Real(kind = RDbl)           :: p(4)
    Type(VecType)               :: newvec

    !** Map the coordinates to the unit cell
    newvec = simcell_maptouc(simcell,atvec)

    !** Check if the simcell is orthorhombic
    orthoflag = simcell_isortho(simcell)

    !** Do the interpolation with the appropriate routine
    Select Case(storebase_nderivs(results))
    Case (0)
      Call pothermite(header,orthoflag,newvec%comp,&
          map%ptr, map%r, map%rx, map%ry, map%rz, &
          map%rxy, map%rxz,map%ryz,map%rxyz, interpolate_int, p(1))

      If (Present(factor)) Then
        results%nrg = results%nrg + p(1)*factor
      Else
        results%nrg = results%nrg + p(1)        
      End If

    Case (1)
      Call forcehermite(header,orthoflag,newvec%comp,&
          map%ptr, map%r, map%rx, map%ry, map%rz, &
          map%rxy, map%rxz,map%ryz,map%rxyz, interpolate_int, p)

      If (Present(factor)) Then
        results%nrg = results%nrg + p(1)*factor
        results%force(1) = results%force(1) + p(2)*factor
        results%force(2) = results%force(2) + p(3)*factor
        results%force(3) = results%force(3) + p(4)*factor
      Else
        results%nrg = results%nrg + p(1)
        results%force(1) = results%force(1) + p(2)
        results%force(2) = results%force(2) + p(3)
        results%force(3) = results%force(3) + p(4)
      End If

    Case Default
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Only energies and gradients are available from map interpolation'
      Stop      
    End Select

  End Function interpolate_int

  !----------------------------------------------------------------------------
  ! Returns number of points tabulated in the map
  ! Requires:  map -- potential map structure
  !----------------------------------------------------------------------------
  Integer Function interpolate_npts(map)
    Type(Map_Core_Info), Intent(In)        :: map

    interpolate_npts = map%npts

  End Function interpolate_npts

  !----------------------------------------------------------------------------
  ! Hermite - a piecewise cubic hermite, 3D interpolation routine
  !
  ! Written by: R.L. June, 5/18/89
  !
  ! References: Shultz, M., Spline Analysis, p. 29
  !             Moler and others, Numerical Methods and Software, Ch 4
  ! 
  ! Output variables:
  !
  !        mapflag Logical - flags the fact that the tabulated data was
  !                not located and that the function value must be
  !                explicitly calculated.
  !
  !        p    Real*8, dimension(4) - interpolated function values p(1)
  !             is the potential, p(1,2,3) are forces in x,y,z directions
  !----------------------------------------------------------------------------

  Subroutine forcehermite(header,orthoflag,coord,&
      ptr,f,fx,fy,fz,fxy,fxz,fyz,fxyz, mapflag, p)
  
    Type(Header_Info), Intent(IN)             :: header
    Logical, Intent(IN)                       :: orthoflag
    Real(kind=RDbl), Dimension(3), Intent(IN) :: coord
    Integer, Dimension(:,:,:), Pointer        :: ptr
    Real(kind=RSgl), Dimension(:), Pointer    :: f,fx,fy,fz
    Real(kind=RSgl), Dimension(:), Pointer    :: fxy,fxz,fyz,fxyz 
    Real(kind=RDbl), Intent(OUT)              :: p(4)
    Logical, Intent(OUT)                      :: mapflag

    Integer                                   :: ia,ii,jj,kk,i,m
    Real(kind=RDbl)                           :: h00,h01,h10,h11
    Real(kind=RDbl)                           :: sum1,sum2,sum3,sum4,t,u,v
!!$    Real(kind=RDbl)                           :: h0,hd0,h1,hd1
!!$    Real(Kind=RDbl) :: xx,ds
!!$    Real(kind=RDbl)                           :: coord2(3),dist(3)
    Integer, Dimension(3) :: indices
                                              
    Parameter (ia=64)                         
                                              
    Real(kind=RDbl)                           :: hi(ia),hj(ia),hk(ia),ft(ia)
    Real(kind=RDbl)                           :: hdi(ia),hdj(ia),hdk(ia)
    Real(kind=RDbl), Dimension(8)             :: hx0,hx1,hy0,hy1,hz0,hz1
    Real(kind=RDbl), Dimension(8)             :: hdx0,hdx1,hdy0,hdy1,hdz0,hdz1
!!$    Integer                                   :: flag
    Integer                                   :: mi(8),iiplus,jjplus,kkplus

!!$    !** Define the hermite basis functions as FORTRAN statement functions
!!$    !** This defines what is like a macro in C.
!!$    h0(xx,ds)=xx*xx*(3.0d0-2.0d0*xx)
!!$    hd0(xx,ds)=6.0d0*xx*(1.0d0-xx)
!!$    h1(xx,ds)=ds*xx*xx*(xx-1.0d0)
!!$    hd1(xx,ds)=ds*xx*(3.0d0*xx-2.0d0)

    !**Setup for the interpolation
    !*** get the array indexes - assumes an evenly spaced grid.
!!$    ii = int(coord(1)/header%stepx)+1
!!$    jj = int(coord(2)/header%stepy)+1
!!$    kk = int(coord(3)/header%stepz)+1
    Call mapheader_getIndices(header,coord(1),coord(2),coord(3),indices)
    ii = indices(1)
    jj = indices(2)
    kk = indices(3)

#if CHECKPOT        
    !*** prevent array out of bounds error by checking the extremities of the
    !*** array indicies.
    flag = 0
    If (ii > nbrx) flag = 1
    If (ii < 1) flag = 1
    If (jj > nbry) flag = 1
    If (jj < 1) flag = 1
    If (kk > nbrz) flag = 1
    If (kk < 1) flag = 1
    If (flag == 1) Then
      Write(*,*) 'Hermite was not passed a principal unit cell point'
      Stop
    Endif
#endif

    !** calculate the scaled coordinates
    t = (coord(1) - header%xgrid(ii))/header%stepx
    u = (coord(2) - header%ygrid(jj))/header%stepy
    v = (coord(3) - header%zgrid(kk))/header%stepz

    !** get the nearest grid points
    iiplus = ii + 1
    jjplus = jj + 1
    kkplus = kk + 1

    !** use periodic boundary conditions if this is an orthorhombic unit cell
    !** NOTE: pbcs will NOT work with non-orthorhombic unit cells
    If (orthoflag) Then
      If (ii == header%nbrx) iiplus = 1 
      If (jj == header%nbry) jjplus = 1 
      If (kk == header%nbrz) kkplus = 1 
    Endif

!!$#if DEBUG
!!$    coord2 = coord
!!$    coord2 = coord2 - origin
!!$    Call slant(coord2)
!!$
!!$    dist(1) = min(coord2(1), (Abs(coord2(1) - alpha)))
!!$    dist(2) = min(coord2(2), (Abs(coord2(2) - beta)))
!!$    dist(3) = min(coord2(3), (Abs(coord2(3) - gamma)))
!!$
!!$    If (dist(1) < 0.05d0) Then
!!$      Write(10,*) 'Hermite: molecule at border x ',istep,debugflag,dist(1)
!!$    Endif
!!$
!!$    If (dist(2) < 0.05d0) Then
!!$      Write(10,*) 'Hermite: molecule at border y ',istep,debugflag,dist(2)
!!$    Endif
!!$
!!$    If (dist(3) < 0.05d0) Then
!!$      Write(10,*) 'Hermite: molecule at border z ',istep,debugflag,dist(3)
!!$    Endif
!!$#endif
!!$
!!$#if DEBUG
!!$    If (debugflag == 9) Then
!!$      Write(*,'(1x,a,8i8,4f8.3)') 'Hermite: ',ii,jj,kk,iiplus,jjplus,kkplus
!!$    Endif
!!$#endif

    !** get the appropriate indices in the interpolation arrays
    mi(1) = ptr(kk    ,jj    ,ii    )
    mi(2) = ptr(kk    ,jj    ,iiplus)
    mi(3) = ptr(kk    ,jjplus,ii    )
    mi(4) = ptr(kkplus,jj    ,ii    )
    mi(5) = ptr(kk    ,jjplus,iiplus)
    mi(6) = ptr(kkplus,jj    ,iiplus)
    mi(7) = ptr(kkplus,jjplus,ii    )
    mi(8) = ptr(kkplus,jjplus,iiplus)

    !** mapflag flags the fact that we have missed the potential summation
    !** mapflag is reset in the calling routine. note that we must avoid
    !** overflow of the integer dynamic range
    mapflag = .True.
    Do i = 1,8
      If (mi(i) == 0) Then
        mapflag = .False.
        p = 0.0_RDbl
        Return
      End if
    End Do

    !-----------------------------------------------------------------------
    !  Look up the required function values and derivatives
    !-----------------------------------------------------------------------

    Do i = 1,8
      m = (i-1)*8
      ft(m+1) = f(mi(i))
      ft(m+2) = fx(mi(i))
      ft(m+3) = fy(mi(i))
      ft(m+4) = fz(mi(i))
      ft(m+5) = fxy(mi(i))
      ft(m+6) = fxz(mi(i))
      ft(m+7) = fyz(mi(i))
      ft(m+8) = fxyz(mi(i))
    EndDo

    !--------------------------------------------------------------------------
    !  Do the interpolation
    !--------------------------------------------------------------------------

    !*** look up the basis function values at the input coordinates

    h00=h0(1-t, header%stepx)
    h01=h0(    t, header%stepx)
    h10=h1(1-t,-header%stepx)
    h11=h1(    t, header%stepx)

    hx0(1)=h00
    hx0(2)=h01
    hx0(3)=h00
    hx0(4)=h00
    hx0(5)=h01
    hx0(6)=h01
    hx0(7)=h00
    hx0(8)=h01

    hx1(1)=h10
    hx1(2)=h11
    hx1(3)=h10
    hx1(4)=h10
    hx1(5)=h11
    hx1(6)=h11
    hx1(7)=h10
    hx1(8)=h11

    Do i = 1,8
      m=(i-1)*8
      hi(m+1)=hx0(i)
      hi(m+2)=hx1(i)
      hi(m+3)=hx0(i)
      hi(m+4)=hx0(i)
      hi(m+5)=hx1(i)
      hi(m+6)=hx1(i)
      hi(m+7)=hx0(i)
      hi(m+8)=hx1(i)
    End Do

    !*** do the derivative function
    h00=-hd0(1-t, header%stepx)
    h01=hd0(    t, header%stepx)
    h10=-hd1(1-t,-header%stepx)
    h11=hd1(    t, header%stepx)

    hdx0(1)=h00
    hdx0(2)=h01
    hdx0(3)=h00
    hdx0(4)=h00
    hdx0(5)=h01
    hdx0(6)=h01
    hdx0(7)=h00
    hdx0(8)=h01

    hdx1(1)=h10
    hdx1(2)=h11
    hdx1(3)=h10
    hdx1(4)=h10
    hdx1(5)=h11
    hdx1(6)=h11
    hdx1(7)=h10
    hdx1(8)=h11

    Do i=1,8
      m=(i-1)*8
      hdi(m+1)=hdx0(i)
      hdi(m+2)=hdx1(i)
      hdi(m+3)=hdx0(i)
      hdi(m+4)=hdx0(i)
      hdi(m+5)=hdx1(i)
      hdi(m+6)=hdx1(i)
      hdi(m+7)=hdx0(i)
      hdi(m+8)=hdx1(i)
    End Do
    !-------------------------------------------------------------------
    !*** do the y coordinate
    h00=h0(1-u, header%stepy)
    h01=h0(    u, header%stepy)
    h10=h1(1-u,-header%stepy)
    h11=h1(    u, header%stepy)

    hy0(1)=h00
    hy0(2)=h00
    hy0(3)=h01
    hy0(4)=h00
    hy0(5)=h01
    hy0(6)=h00
    hy0(7)=h01
    hy0(8)=h01

    hy1(1)=h10
    hy1(2)=h10
    hy1(3)=h11
    hy1(4)=h10
    hy1(5)=h11
    hy1(6)=h10
    hy1(7)=h11
    hy1(8)=h11

    Do i=1,8
      m=(i-1)*8
      hj(m+1)=hy0(i)
      hj(m+2)=hy0(i)
      hj(m+3)=hy1(i)
      hj(m+4)=hy0(i)
      hj(m+5)=hy1(i)
      hj(m+6)=hy0(i)
      hj(m+7)=hy1(i)
      hj(m+8)=hy1(i)
    End Do

    !*** do the derivative

    h00=-hd0(1-u, header%stepy)
    h01=hd0(    u, header%stepy)
    h10=-hd1(1-u,-header%stepy)
    h11=hd1(    u, header%stepy)

    hdy0(1)=h00
    hdy0(2)=h00
    hdy0(3)=h01
    hdy0(4)=h00
    hdy0(5)=h01
    hdy0(6)=h00
    hdy0(7)=h01
    hdy0(8)=h01

    hdy1(1)=h10
    hdy1(2)=h10
    hdy1(3)=h11
    hdy1(4)=h10
    hdy1(5)=h11
    hdy1(6)=h10
    hdy1(7)=h11
    hdy1(8)=h11

    Do i=1,8
      m=(i-1)*8
      hdj(m+1)=hdy0(i)
      hdj(m+2)=hdy0(i)
      hdj(m+3)=hdy1(i)
      hdj(m+4)=hdy0(i)
      hdj(m+5)=hdy1(i)
      hdj(m+6)=hdy0(i)
      hdj(m+7)=hdy1(i)
      hdj(m+8)=hdy1(i)
    End Do
    !-------------------------------------------------------------------
    !*** do the z coordinate basis functions
    h00=h0(1-v, header%stepz)
    h01=h0(    v, header%stepz)
    h10=h1(1-v,-header%stepz)
    h11=h1(    v, header%stepz)

    hz0(1)=h00
    hz0(2)=h00
    hz0(3)=h00
    hz0(4)=h01
    hz0(5)=h00
    hz0(6)=h01
    hz0(7)=h01
    hz0(8)=h01

    hz1(1)=h10
    hz1(2)=h10
    hz1(3)=h10
    hz1(4)=h11
    hz1(5)=h10
    hz1(6)=h11
    hz1(7)=h11
    hz1(8)=h11

    Do i=1,8
      m=(i-1)*8
      hk(m+1)=hz0(i)
      hk(m+2)=hz0(i)
      hk(m+3)=hz0(i)
      hk(m+4)=hz1(i)
      hk(m+5)=hz0(i)
      hk(m+6)=hz1(i)
      hk(m+7)=hz1(i)
      hk(m+8)=hz1(i)
    End Do

    !*** evaluate the derivative

    h00=-hd0(1-v, header%stepz)
    h01=hd0(    v, header%stepz)
    h10=-hd1(1-v,-header%stepz)
    h11=hd1(    v, header%stepz)

    hdz0(1)=h00
    hdz0(2)=h00
    hdz0(3)=h00
    hdz0(4)=h01
    hdz0(5)=h00
    hdz0(6)=h01
    hdz0(7)=h01
    hdz0(8)=h01

    hdz1(1)=h10
    hdz1(2)=h10
    hdz1(3)=h10
    hdz1(4)=h11
    hdz1(5)=h10
    hdz1(6)=h11
    hdz1(7)=h11
    hdz1(8)=h11

    Do i=1,8
      m=(i-1)*8
      hdk(m+1)=hdz0(i)
      hdk(m+2)=hdz0(i)
      hdk(m+3)=hdz0(i)
      hdk(m+4)=hdz1(i)
      hdk(m+5)=hdz0(i)
      hdk(m+6)=hdz1(i)
      hdk(m+7)=hdz1(i)
      hdk(m+8)=hdz1(i)
    End Do
    !--------------------------------------------------------------------------
    !*** evaluate the summation for the expressions

    sum1=0.0_RDbl
    sum2=0.0_RDbl
    sum3=0.0_RDbl
    sum4=0.0_RDbl

    Do i=1,64
      sum1=sum1+ft(i)*hi(i)*hj(i)*hk(i)
      sum2=sum2+ft(i)*hdi(i)*hj(i)*hk(i)
      sum3=sum3+ft(i)*hi(i)*hdj(i)*hk(i)
      sum4=sum4+ft(i)*hi(i)*hj(i)*hdk(i)
    End Do

    p(1) = sum1
    p(2) = -sum2/header%stepx
    p(3) = -sum3/header%stepy
    p(4) = -sum4/header%stepz

#if DEBUG
    If (debugflag == 9) Then
      Write(*,'(1x,a,8i8,4f8.3)') 'Hermite: ',(mi(i),i=1,8),p
    Endif
#endif

    Return
  End Subroutine forcehermite

  !----------------------------------------------------------------------------
  ! Hermite - a piecewise cubic hermite, 3D interpolation routine
  !
  ! Written by: R.L. June, 5/18/89
  !
  ! References: Shultz, M., Spline Analysis, p. 29
  !             Moler and others, Numerical Methods and Software, Ch 4
  ! 
  ! Output variables:
  !
  !        mapflag Logical - flags the fact that the tabulated data was
  !                not located and that the function value must be
  !                explicitly calculated.
  !
  !        p    Real*8, dimension(4) - interpolated function values p(1)
  !             is the potential
  !----------------------------------------------------------------------------
  Subroutine pothermite(header,orthoflag,coord,&
      ptr,f,fx,fy,fz,fxy,fxz,fyz,fxyz, mapflag, p)
  
    Type(Header_Info), Intent(IN)               :: header
    Logical, Intent(IN)                         :: orthoflag
    Real(kind=RDbl), Dimension(3), Intent(IN)   :: coord
    Integer, Dimension(:,:,:), Pointer          :: ptr
    Real(kind=RSgl), Dimension(:), Pointer      :: f,fx,fy,fz
    Real(kind=RSgl), Dimension(:), Pointer      :: fxy,fxz,fyz,fxyz 
    Real(kind=RDbl), Intent(OUT)                :: p
    Logical, Intent(OUT)                        :: mapflag

    Integer                                     :: ia,ii,jj,kk,i,m
    Integer, Dimension(3) :: indices
    Real(kind=RDbl)                             :: h00,h01,h10,h11
    Real(kind=RDbl)                             :: sum1,t,u,v
!!$    Real(kind=RDbl)                             :: h0,hd0,h1,hd1
!!$    Real(kind=RDbl) :: xx,ds
!!$    Real(kind=RDbl)                             :: coord2(3),dist(3)

    Parameter (ia=64)

    Real(kind=RDbl)                             :: hi(ia),hj(ia),hk(ia),ft(ia)
!!$    Real(kind=RDbl)                             :: hdi(ia),hdj(ia),hdk(ia)
    Real(kind=RDbl), Dimension(8)               :: hx0,hx1,hy0,hy1,hz0,hz1
!!$    Integer                                     :: flag
    Integer                                     :: mi(8),iiplus,jjplus,kkplus

    !** Now functions in this module
!!$    !** Define the hermite basis functions as FORTRAN statement functions
!!$    !** This defines what is like a macro in C.
!!$    h0(xx,ds)=xx*xx*(3.0d0-2.0d0*xx)
!!$    h1(xx,ds)=ds*xx*xx*(xx-1.0d0)

    !**Setup for the interpolation
    !*** get the array indexes - assumes an evenly spaced grid.
!!$    ii = int(coord(1)/header%stepx)+1
!!$    jj = int(coord(2)/header%stepy)+1
!!$    kk = int(coord(3)/header%stepz)+1
    Call mapheader_getIndices(header,coord(1),coord(2),coord(3),indices)
    ii = indices(1)
    jj = indices(2)
    kk = indices(3)

#if CHECKPOT        
    !*** prevent array out of bounds error by checking the extremities of the
    !*** array indicies.
    flag = 0
    If (ii > nbrx) flag = 1
    If (ii < 1) flag = 1
    If (jj > nbry) flag = 1
    If (jj < 1) flag = 1
    If (kk > nbrz) flag = 1
    If (kk < 1) flag = 1
    If (flag == 1) Then
      Write(*,*) 'Hermite was not passed a principal unit cell point'
      Stop
    Endif
#endif

    !** calculate the scaled coordinates
    t = (coord(1) - header%xgrid(ii))/header%stepx
    u = (coord(2) - header%ygrid(jj))/header%stepy
    v = (coord(3) - header%zgrid(kk))/header%stepz

    !** get the nearest grid points
    iiplus = ii + 1
    jjplus = jj + 1
    kkplus = kk + 1

    !** use periodic boundary conditions if this is an orthorhombic unit cell
    !** NOTE: pbcs will NOT work with non-orthorhombic unit cells
    If (orthoflag) Then
      If (ii == header%nbrx) iiplus = 1 
      If (jj == header%nbry) jjplus = 1 
      If (kk == header%nbrz) kkplus = 1 
    Endif

    !** get the appropriate indices in the interpolation arrays
    mi(1) = ptr(kk    ,jj    ,ii    )
    mi(2) = ptr(kk    ,jj    ,iiplus)
    mi(3) = ptr(kk    ,jjplus,ii    )
    mi(4) = ptr(kkplus,jj    ,ii    )
    mi(5) = ptr(kk    ,jjplus,iiplus)
    mi(6) = ptr(kkplus,jj    ,iiplus)
    mi(7) = ptr(kkplus,jjplus,ii    )
    mi(8) = ptr(kkplus,jjplus,iiplus)

    !** mapflag flags the fact that we have missed the potential summation
    !** mapflag is reset in the calling routine. note that we must avoid
    !** overflow of the integer dynamic range
    mapflag = .True.
    Do i = 1,8
      If (mi(i) == 0) Then
        mapflag = .False.
        p = 0.0_RDbl
        Return
      Endif
    EndDo

    !----------------------------------------------------------------------------
    !  Look up the required function values and derivatives
    !----------------------------------------------------------------------------

    Do i = 1,8
      m = (i-1)*8
      ft(m+1) = f(mi(i))
      ft(m+2) = fx(mi(i))
      ft(m+3) = fy(mi(i))
      ft(m+4) = fz(mi(i))
      ft(m+5) = fxy(mi(i))
      ft(m+6) = fxz(mi(i))
      ft(m+7) = fyz(mi(i))
      ft(m+8) = fxyz(mi(i))
    EndDo

    !----------------------------------------------------------------------------
    !  Do the interpolation
    !----------------------------------------------------------------------------

    !*** look up the basis function values at the input coordinates

    h00=h0(1-t, header%stepx)
    h01=h0(    t, header%stepx)
    h10=h1(1-t,-header%stepx)
    h11=h1(    t, header%stepx)

    hx0(1)=h00
    hx0(2)=h01
    hx0(3)=h00
    hx0(4)=h00
    hx0(5)=h01
    hx0(6)=h01
    hx0(7)=h00
    hx0(8)=h01

    hx1(1)=h10
    hx1(2)=h11
    hx1(3)=h10
    hx1(4)=h10
    hx1(5)=h11
    hx1(6)=h11
    hx1(7)=h10
    hx1(8)=h11

    Do i = 1,8
      m=(i-1)*8
      hi(m+1)=hx0(i)
      hi(m+2)=hx1(i)
      hi(m+3)=hx0(i)
      hi(m+4)=hx0(i)
      hi(m+5)=hx1(i)
      hi(m+6)=hx1(i)
      hi(m+7)=hx0(i)
      hi(m+8)=hx1(i)
    End Do

    !-------------------------------------------------------------------
    !*** do the y coordinate

    h00=h0(1-u, header%stepy)
    h01=h0(    u, header%stepy)
    h10=h1(1-u,-header%stepy)
    h11=h1(    u, header%stepy)

    hy0(1)=h00
    hy0(2)=h00
    hy0(3)=h01
    hy0(4)=h00
    hy0(5)=h01
    hy0(6)=h00
    hy0(7)=h01
    hy0(8)=h01

    hy1(1)=h10
    hy1(2)=h10
    hy1(3)=h11
    hy1(4)=h10
    hy1(5)=h11
    hy1(6)=h10
    hy1(7)=h11
    hy1(8)=h11

    Do i=1,8
      m=(i-1)*8
      hj(m+1)=hy0(i)
      hj(m+2)=hy0(i)
      hj(m+3)=hy1(i)
      hj(m+4)=hy0(i)
      hj(m+5)=hy1(i)
      hj(m+6)=hy0(i)
      hj(m+7)=hy1(i)
      hj(m+8)=hy1(i)
    End Do

    !-------------------------------------------------------------------
    !*** do the z coordinate basis functions

    h00=h0(1-v, header%stepz)
    h01=h0(    v, header%stepz)
    h10=h1(1-v,-header%stepz)
    h11=h1(    v, header%stepz)

    hz0(1)=h00
    hz0(2)=h00
    hz0(3)=h00
    hz0(4)=h01
    hz0(5)=h00
    hz0(6)=h01
    hz0(7)=h01
    hz0(8)=h01

    hz1(1)=h10
    hz1(2)=h10
    hz1(3)=h10
    hz1(4)=h11
    hz1(5)=h10
    hz1(6)=h11
    hz1(7)=h11
    hz1(8)=h11

    Do i=1,8
      m=(i-1)*8
      hk(m+1)=hz0(i)
      hk(m+2)=hz0(i)
      hk(m+3)=hz0(i)
      hk(m+4)=hz1(i)
      hk(m+5)=hz0(i)
      hk(m+6)=hz1(i)
      hk(m+7)=hz1(i)
      hk(m+8)=hz1(i)
    End Do

    !--------------------------------------------------------------------------
    !*** evaluate the summation for the expressions

    sum1=0.0_RDbl

    Do i=1,64
      sum1=sum1+ft(i)*hi(i)*hj(i)*hk(i)
    End Do

    p = sum1

#if DEBUG
    If (debugflag == 9) Then
      Write(*,'(1x,a,8i8,4f8.3)') 'Hermite: ',(mi(i),i=1,8),p
    Endif
#endif

    Return
  End Subroutine pothermite

  !----------------------------------------------------------------------------
  ! Hermite basis function h0
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function h0(xx,ds)
    Real(Kind=RDbl), Intent(In) :: xx
    Real(Kind=RDbl), Intent(In) :: ds

    h0=xx*xx*(3.0d0-2.0d0*xx)
  End Function h0

  !----------------------------------------------------------------------------
  ! Hermite basis function h1
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function h1(xx,ds)
    Real(Kind=RDbl), Intent(In) :: xx
    Real(Kind=RDbl), Intent(In) :: ds

    h1 = ds*xx*xx*(xx-1.0d0)
  End Function h1

  !----------------------------------------------------------------------------
  ! Hermite basis function hd0
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function hd0(xx,ds)
    Real(Kind=RDbl), Intent(In) :: xx
    Real(Kind=RDbl), Intent(In) :: ds

    hd0 = 6.0d0*xx*(1.0d0-xx)
  End Function hd0

  !----------------------------------------------------------------------------
  ! Hermite basis function hd1
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function hd1(xx,ds)
    Real(Kind=RDbl), Intent(In) :: xx
    Real(Kind=RDbl), Intent(In) :: ds

    hd1 = ds*xx*(3.0d0*xx-2.0d0)
  End Function hd1

  !----------------------------------------------------------------------------
  ! Reads a map core from a given unit number
  ! Requires:  map -- map core to read from file
  !            unit -- filename to write from
  !----------------------------------------------------------------------------
  Subroutine interpolate_readcore(map,unit)
    Type(Map_Core_Info), Intent(InOut)  :: map
    Integer, Intent(In)                 :: unit

    Read(unit) map%ptr
    Read(unit) map%r
    Read(unit) map%rx
    Read(unit) map%ry
    Read(unit) map%rz
    Read(unit) map%rxy
    Read(unit) map%rxz
    Read(unit) map%ryz
    Read(unit) map%rxyz

  End Subroutine interpolate_readcore

  !----------------------------------------------------------------------------
  ! Writes a map core to a given unit number
  ! Requires:  map -- map core to write to file
  !            unit -- filename to write to
  !----------------------------------------------------------------------------
  Subroutine interpolate_writecore(map,unit)
    Type(Map_Core_Info), Intent(In)  :: map
    Integer, Intent(In)              :: unit

    !** Now write the data points to the file
    !** First is the 3D pointer array
    Write(unit) map%ptr
    Write(unit) map%r
    Write(unit) map%rx
    Write(unit) map%ry
    Write(unit) map%rz
    Write(unit) map%rxy
    Write(unit) map%rxz
    Write(unit) map%ryz
    Write(unit) map%rxyz

  End Subroutine interpolate_writecore

  !----------------------------------------------------------------
  ! Cleanup the core map pointers
  ! Requires:  map -- potential map to clean
  !----------------------------------------------------------------
  Subroutine interpolate_clean(map)
    Type(Map_Core_Info), Intent(InOut) :: map

    Integer    :: error

    If (Associated(map%ptr)) Then
      Deallocate(map%ptr, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'core%ptr')
    End If

    If (Associated(map%r)) Then
      Deallocate(map%r, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%r')
    End If

    If (Associated(map%rx)) Then
      Deallocate(map%rx, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rx')
    End If

    If (Associated(map%ry)) Then
      Deallocate(map%ry, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%ry')
    End If

    If (Associated(map%rz)) Then
      Deallocate(map%rz, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rz')
    End If

    If (Associated(map%rxy)) Then
      Deallocate(map%rxy, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rxy')
    End If

    If (Associated(map%rxz)) Then
      Deallocate(map%rxz, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rxz')
    End If

    If (Associated(map%ryz)) Then
      Deallocate(map%ryz, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%ryz')
    End If

    If (Associated(map%rxyz)) Then
      Deallocate(map%rxyz, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'core%rxyz')
    End If

  End Subroutine interpolate_clean

End Module interpolate
