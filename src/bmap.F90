!------------------------------------------------------
! This module contains the data structure and code
! to read in and work with the bias map for insertions 
! into a zeolite
!-----------------------------------------------------
Module bmap

  Use defaults, Only: strLen, RDbl, RSgl, MAX_MAPS, kcalmole_kb, dashedline,&
      zero, one, default_MAX_EXP
  Use utils, Only: split, getpath, isfileopen,allocErrDisplay, cleanstring
  Use file, Only: file_getunit, file_open
  Use mapheader, Only: mapheader_display, Header_Info, mapheader_cleanup, &
      mapheader_getheader
  Use random, Only: rranf
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use general, Only: genparams
  Use simcell, Only: SimCell_Params, simcell_isortho, simcell_pbcptFund, simcell_maptouc, simcell_getcellorigin
  Implicit None
  Save

  Private
  Public :: BiasMap_Params, CubeletInfo, bmap_display, bmap_getbiasindex, &
      bmap_getncubelets, bmap_getcellwt, bmap_getcellCOM, bmap_getbiaspt, &
      bmap_getbiaswt, bmap_init, bmap_makeXYZ, DISPLAY_ONLY_FLAG

  !*** Bias Potential Maps ***
  Type CubeletInfo
    Real(kind=Rdbl)             :: x, y, z
    Real(kind=Rdbl)             :: wt
    Real(kind=Rdbl)             :: cum_prob
  End Type CubeletInfo

  Type BiasMap_Params
    Character(len=strLen)       :: filename
    Real(kind=Rdbl), Dimension(:), Pointer :: rb
    Real(kind=Rsgl), Dimension(:), Pointer :: rb_temp
    Integer, Dimension(:,:,:), Pointer     :: ptrb
    Type(Header_Info)        :: header
    Type(CubeletInfo), Dimension(:), Pointer :: cube
    Integer, Dimension(:,:,:), Pointer  :: cube_index
    Integer                  :: ncubelets, numtab
    Real(kind=RDbl)          :: cum_wt
    Real(kind=RDbl)          :: ucutoff  ! incl. only cubelets w/ u <= ucutoff
    Real(kind=RDbl)          :: rti      ! 1.0/(Gas*tk), used for the biasing
    Logical                  :: reversebias
    Real(kind=RDbl) :: stepx, stepy, stepz
    Real(kind=RDbl), Dimension(3) :: eff
    Type(SimCell_Params) :: simcell
  End Type BiasMap_Params

  Type(BiasMap_Params), Allocatable, Dimension(:), Target  :: bmaps
  Integer                                    :: total_bmaps = 0

  ! any point where energy is higer than (in kcal) this is not considered 
  Real(kind=RDbl) , Parameter :: BMAP_DEFAULT_CUTOFF= 30.00_RDbl
  Logical :: DISPLAY_ONLY_FLAG = .False.
  Real(kind=RDbl) , Parameter :: TRESHOLD = 0.10_RDbl
  Logical :: gorthoflag = .True.
Contains
  !------------------------------------------------------------------------
  ! Initialize the bias map 'bmap'.  The bias can be used
  ! to preferentially sample the low energy regions of the zeolite or the
  ! high energy regions.  The optional flag 'optreverse' is used to decide
  ! how the biasing will be done.  If 'optreverse' is present and true then
  ! we bias towards the high energy regions otherwise if absent or false then
  ! we bias towards the low energy regions
  !------------------------------------------------------------------------
  Subroutine bmap_init(bmap, eff, tk, filename, simcell, optreverse)
    Type(BiasMap_Params), Pointer :: bmap
    Real(kind=RDbl), Dimension(:), Intent(in) :: eff
    Real(kind=RDbl), Intent(in):: tk   ! temperature in Kelvin
    Character(*), Intent(in)   :: filename
    Type(SimCell_Params), Intent(IN)             :: simcell
    Logical, Optional, Intent(in) :: optreverse

    Logical       :: reversebias

    Integer                        :: i, mapno
    Integer                        :: error,unitno
    Integer                        :: nbrx,nbry,nbrz,nfsize
    Character(len=150)             :: tempstring

    !** Initialize the list that keeps holds all the bias maps
    If (total_bmaps == 0) Then
      Allocate(bmaps(MAX_MAPS), Stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End If

    !** Choose how we bias
    If (Present(optreverse)) Then
      reversebias = optreverse
    Else
      reversebias = .False.
    End If

    !** Check if the map has already been read in
    Do i = 1,total_bmaps
      If ((filename == bmaps(i)%filename).And. &
          (reversebias .Eqv. bmaps(i)%reversebias)) Then
!LC        Write(*,*) "Pointing to the already initialized map from :",&
!LC            Trim(bmaps(i)%filename)
        bmap => bmaps(i)
        Return
      End If
    End Do

    !** Read the new bias map    
    If (MAX_MAPS == total_bmaps) Then
      Write(0,'(1x,2a,i4, a, i4)') __FILE__," : ",__LINE__, &
          "Total number of bmaps exceeded ", MAX_MAPS
      Stop
    Endif
    total_bmaps = total_bmaps + 1
    mapno = total_bmaps
    bmaps(mapno)%filename = filename
    bmaps(mapno)%reversebias = reversebias

    !** Get the path and generate the full filename
    tempstring = Trim(getpath('PMAPDIR'))//Trim(filename)

    !** Open the file if it is not opened
    unitno = file_open(tempstring,100)
    
    !** Get the bmap header information 
    Call mapheader_getheader(bmaps(mapno)%header, eff, unitno)

    !** Allocate the necessary memory
    nbrx = bmaps(mapno)%header%nbrx
    nbry = bmaps(mapno)%header%nbry
    nbrz = bmaps(mapno)%header%nbrz
    nfsize = bmaps(mapno)%header%nfsize


    ! copy some necessary info here. step length in each direction
    ! and the total lengths of the box in each direction
    bmaps(mapno)%stepx=bmaps(mapno)%header%stepx
    bmaps(mapno)%stepy=bmaps(mapno)%header%stepy
    bmaps(mapno)%stepz=bmaps(mapno)%header%stepz
    bmaps(mapno)%eff(1:3)=bmaps(mapno)%header%eff(1:3)
    Allocate(bmaps(mapno)%ptrb(nbrz,nbry,nbrx), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    
    Allocate(bmaps(mapno)%rb(nfsize), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(bmaps(mapno)%rb_temp(nfsize), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
    !** Now read in the indexes and potentials and close the file
    Read(unitno) bmaps(mapno)%ptrb
    Read(unitno) bmaps(mapno)%rb_temp  ! pmap contains single-precision-real 
    ! convert to double precision
    bmaps(mapno)%rb(1:nfsize) = bmaps(mapno)%rb_temp(1:nfsize)

    ! dont need this anymore
    Deallocate(bmaps(mapno)%rb_temp, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Close(unitno)

    !** Initialize the temperature field
    bmaps(mapno)%rti = (1.0_RDbl)/(tk*kcalmole_kb)

    !** Initialize the cubelet information
    bmaps(mapno)%ucutoff = BMAP_DEFAULT_CUTOFF
    Call bmap_cubeinit(bmaps(mapno), simcell)

    !** After initializing cubes
    !** We can deallocate ptrb and rb since they will be used no longer
    Deallocate(bmaps(mapno)%ptrb, bmaps(mapno)%rb, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Assign the pointer to the newly read bmap
    bmap => bmaps(mapno)


  End Subroutine bmap_init

  !-------------------------------------------------
  ! Initializes the cubelet information
  !-------------------------------------------------
  Subroutine bmap_cubeinit(bmap, simcell)
    Type(BiasMap_Params), Intent(inout) :: bmap 
    Type(SimCell_Params), Intent(IN)             :: simcell
    Integer    :: nfsizeb, nbrx, nbry, nbrz, numtab, ncubelets_nonortho
    Integer    :: error, i, j, k, ptr, count=0
    Real(kind=RDbl)  :: wt, u, rti, sum, expnnt, max, min, minu, maxu, &
        allowed_min, maxnrgval
    Real(kind=RDbl)  :: delta
    Type(VecType)                       :: position, positionpbc, deltavec
    Logical :: orthoflag


    !** Copy simcell data for local use
    bmap%simcell = simcell
    rti    = bmap%rti
    numtab = 0
    bmap%numtab = 0
    bmap%cum_wt = 0.0_RDbl
    nfsizeb = bmap%header%nfsize  ! number of tabulated points(nodes) in pmap
    nbrx    = bmap%header%nbrx    ! number of total nodes in pmap in xdirn
    nbry    = bmap%header%nbry    ! number of total nodes in pmap in ydirn
    nbrz    = bmap%header%nbrz    ! number of total nodes in pmap in zdirn

    ! holds the cubes (not nodes)
    Allocate(bmap%cube_index(nbrx-1,nbry-1, nbrz-1), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    bmap%cube_index=0       ! all untabualted cubes will have 0 as the index


    orthoflag=simcell_isortho(simcell)
    gorthoflag = orthoflag

    !** There is a problem with the values in rb. We are going to be adding up
    !** exp(-U/kT) to find the %cum_prob values. If U is very large This
    !** could lead to values of zeroes being added up, later giving problems
    !** in bmap_getbiasindex routine. So we adjust the rb values. in cube_init

    max = -1.00
    min = 1.00e20
    maxu = -1.00e20
    minu = 1.00e20
    numtab = 0
    sum=zero
    ncubelets_nonortho = 0
    !** First go through all values decide the required size 
    !** adjust rb values 
    Do i=1,nbrx-1
      Do j=1,nbry-1
        Do k=1,nbrz-1
          bmap%cube_index(i,j,k)=0
          ptr = bmap%ptrb(k,j,i) ! this was how it was indexed in pmap (k,j,i)
          !** Test of orthogonality and then exclusion of points outside  of the non-orthorombic unit cell
          If(.Not.orthoflag) Then
            position%comp(1)=(bmap%header%stepx) * real(i-1)+0.5_RDbl*(bmap%header%stepx)
            position%comp(2)=(bmap%header%stepy) * real(j-1)+0.5_RDbl*(bmap%header%stepy)
            position%comp(3)=(bmap%header%stepz) * real(k-1)+0.5_RDbl*(bmap%header%stepz)
            Call simcell_pbcptFund(simcell, position, positionpbc)
            deltavec=position-positionpbc
            delta=deltavec*deltavec
            If(delta>TRESHOLD) cycle ! this means the point lays outside non-orthorombic boundaries
            ncubelets_nonortho = ncubelets_nonortho + 1
          End If

          If(ptr /= 0) Then
            If (bmap%rb(ptr) <= bmap%ucutoff ) Then
              numtab = numtab + 1
              bmap%cube_index(i,j,k)=numtab
              If (u<minu) minu=u
              If (u>maxu) maxu=u

    u = bmap%rb(ptr)  ! energy from pmap
              expnnt=-rti*u

              If (expnnt>default_MAX_EXP) Then
                expnnt=default_MAX_EXP
              Elseif(expnnt<(-default_MAX_EXP)) Then
                expnnt=-default_MAX_EXP
              Endif

              ! I dont know why somebody wants reversbias 
              ! I guess it was used for debugging 
              If (.Not. bmap%reversebias) Then
                wt = Exp(expnnt)
              Else
                ! ** THIS IS NOT TESTED, MAY BE NOT NECESSARY
                Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
                Stop
                wt = Exp(-expnnt)
              End If
              If (wt> max) max=wt
              If (wt< min) min=wt
              sum=sum+wt
            End If
          End If
        End Do
      End Do
    End Do

    !** check max min and sum. min belongs to high energy, we have to 
    !** keep the nrg below a maxnrgval so that min/sum > 1.0e-6 (to
    !** avoid round off errors)
    ! just an approximate adjustment, might not be the best
    If (min/sum> 1.00e-6) Then
      ! its OK, no need to adjust %rb
    Elseif (DISPLAY_ONLY_FLAG) Then
      maxnrgval= bmap%ucutoff
    Else
      Write(*,*) "--------------- WARNING during bmap_cubeinit ---------------"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Automatically adjusting bmap energies here"
      ! Adjust the values in rb
      allowed_min = sum* 1.0e-6
      maxnrgval = -Log(allowed_min)/rti
      Write(*,*) "maxnrgval :", maxnrgval
      Write(*,*) "For biasing purposes all energies above this &
          & will be converted to this value"
      Write(*,*) "--------------- END WARNING --------------------------------"
    End If

    Allocate(bmap%cube(numtab), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    max = -1.00
    min =1.00e20
    numtab=0
    ! Note that we want the number of tabulated cubelets, not
    ! the number of tabulated nodes. Nodes run from 1 to nbrx
    Do i=1,nbrx-1
      Do j=1,nbry-1
        Do k=1,nbrz-1

          !** Test of orthogonality and then exclusion of points outside  of the non-orthorombic unit cell
          If(.Not.orthoflag) Then
            position%comp(1)=(bmap%header%stepx) * real(i-1)+0.5_RDbl*(bmap%header%stepx)
            position%comp(2)=(bmap%header%stepy) * real(j-1)+0.5_RDbl*(bmap%header%stepy)
            position%comp(3)=(bmap%header%stepz) * real(k-1)+0.5_RDbl*(bmap%header%stepz)
            call simcell_pbcptFund(simcell, position, positionpbc)
            deltavec=position-positionpbc
            delta=deltavec*deltavec
            If(delta>TRESHOLD) cycle ! this means the point lays outside non-orthorombic boundaries
          End If


          ptr = bmap%ptrb(k,j,i)

          If (ptr /= 0) Then  !** this point was in pmap
            If (bmap%rb(ptr) <= bmap%ucutoff ) Then

              !** check whether we need to adjust %rb values to avoid 
              !** roundoff errors
              If (bmap%rb(ptr)>maxnrgval) bmap%rb(ptr)=maxnrgval

#ifdef ALKANE_CBGCMC
              ! to make sure that alkane growth is possible from all parts of zeolite
              ! not just starting from the intersections
              If (bmap%rb(ptr)>zero) bmap%rb(ptr)=one
              If (bmap%rb(ptr)<zero) bmap%rb(ptr)=zero
#endif

              numtab = numtab + 1
              bmap%cube(numtab)%x = bmap%header%xgrid(i)
              bmap%cube(numtab)%y = bmap%header%ygrid(j)
              bmap%cube(numtab)%z = bmap%header%zgrid(k)
              u = bmap%rb(ptr)
              expnnt=-rti*u
              If (expnnt>default_MAX_EXP) Then
                expnnt=default_MAX_EXP
              Elseif(expnnt<(-default_MAX_EXP)) Then
                expnnt=-default_MAX_EXP
              Endif

              If (.Not. bmap%reversebias) Then
                wt = Exp(expnnt)
              Else
                wt = Exp(-expnnt)
              End If
              If (wt> max) max=wt
              If (wt< min) min=wt

              bmap%cube(numtab)%wt = wt
              bmap%cum_wt = bmap%cum_wt + wt 

            End If
          End If
        End Do
      End Do
    End Do

    bmap%numtab = numtab                    ! no of tabulated cubelets
    bmap%ncubelets = bmap%header%ncubelets  ! total cubelets
    If(.Not.orthoflag) Then
      bmap%ncubelets = ncubelets_nonortho
    End If

    If (numtab==0) Then
      Write(*,*) "No cubes in bmap somehting wrong"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

    !** Set up the cumulative probability which will used to
    !** sample from the distribution
    sum = zero
    Do i = 1,numtab
      sum = sum + bmap%cube(i)%wt
      bmap%cube(i)%cum_prob = sum/bmap%cum_wt
    End Do


#ifdef DEBUG
    Write(*,*) "------ some important bmap details for your perusal --------"
    Write(*,*) "NRg CUTOFF", BMAP_DEFAULT_CUTOFF
    Write(*,*) "x,y,z divisions/uc", nbrx-1, nbry-1, nbrz-1
    Write(*,*) "Total cubes", (nbrx-1)*(nbry-1)*(nbrz-1)
    Write(*,*) "ncubelets in pmap ", bmap%ncubelets
    Write(*,*) "Tabulated points in original pmap : ",  nfsizeb
    Write(*,*) "bmap size (no of tabulated cubes) : ",  numtab 
    Write(*,*) "% volume filled" ,numtab*100.00/((nbrx-1)*(nbry-1)*(nbrz-1))
    Write(*,*) "max and min, sum", max, min, sum
   Write(*,*) "maxu and minu    ", maxu, minu
#endif

  End Subroutine bmap_cubeinit

  !----------------------------------------------------
  ! Gets the total no. of cubelets in the unit cell
  !----------------------------------------------------
  Integer Function bmap_getncubelets(bmap)
    Type(BiasMap_Params), Intent(in) :: bmap
    bmap_getncubelets = bmap%ncubelets
  End Function bmap_getncubelets

  !---------------------------------------------------------------
  ! Get the bias-wt. of the cubelet which contains the point "com"
  ! opt_stop -> True  = stop the simulation, if there is error
  !             False = dont stop, retirn a -1.0 as  bias value
  !---------------------------------------------------------------
  Real(kind=RDbl) Function bmap_getcellwt(bmap, com, opt_stop)
    Type(BiasMap_Params), Intent(in) :: bmap
    Real(kind=RDbl), Dimension(3), Intent(in)   :: com
    Logical, Intent(in) , Optional :: opt_stop
    
    Integer         :: kindx, jindx, iindx, index
    Real(kind=RDbl) :: wt, newwt
    Logical         :: abortflag
    Real(kind=RDbl)  :: delta
    Type(VecType)                       :: position, positionpbc, deltavec

    abortflag=.True.
    If (Present(opt_stop)) Then
      abortflag= opt_stop
    Endif
    iindx = Int(com(1) / bmap%stepx) + 1
    jindx = Int(com(2) / bmap%stepy) + 1
    kindx = Int(com(3) / bmap%stepz) + 1

    If(.Not.gorthoflag) Then
          position%comp(1)=(bmap%stepx) * real(iindx-1)+0.5_RDbl*(bmap%header%stepx)
          position%comp(2)=(bmap%stepy) * real(jindx-1)+0.5_RDbl*(bmap%header%stepy)
          position%comp(3)=(bmap%stepz) * real(kindx-1)+0.5_RDbl*(bmap%header%stepz)
          call simcell_pbcptFund(bmap%simcell, position, positionpbc)
          deltavec=position-positionpbc
          delta=deltavec*deltavec
          If(delta>TRESHOLD) Then
           iindx = Int(positionpbc%comp(1) / bmap%stepx) + 1
           jindx = Int(positionpbc%comp(2) / bmap%stepy) + 1
           kindx = Int(positionpbc%comp(3) / bmap%stepz) + 1
           ! we need to take care of the situation where the non-ortho-cell 
           ! is not inside the bounding box always 
           If (iindx<=0) iindx=iindx+bmap%header%nbrx-1
           If (jindx<=0) jindx=jindx+bmap%header%nbry-1
           If (kindx<=0) kindx=kindx+bmap%header%nbrz-1
          End If
    End If
    index = bmap%cube_index(iindx, jindx, kindx)
    If (index == 0) Then
      If (abortflag) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__,&
          "  ERROR: Unable to get bias map index"
      Write(*,'(a,3f10.5)') "The vector used is ",com
      Write(*,*) "Indices are ",iindx, jindx, kindx
      Write(*,*) " Bmap used is ",bmap%filename
      wt = 0.0_RDbl
        Stop
      Else
        bmap_getcellwt=-1.00_RDbl
        Return
      Endif
    Else

!!$      energy = bmap%rb(index)
!!$      If (energy>BMAP_DEFAULT_CUTOFF) Then
!!$        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$        Write(*,*) "This cube was not used for insertion, no deletion also&
!!$            & should be allowed"
!!$        Write(*,*) "Try increasing default cutoff, or decrease magnitude &
!!$            & of moves like translation or integration"
!!$        Stop
!!$      Endif
!!$      If (.Not. bmap%reversebias) Then
!!$        expnnt = -energy*bmap%rti
!!$      Else
!!$        expnnt = energy*bmap%rti
!!$      End If
!!$
!!$      If (expnnt < -DEFAULT_MAX_EXP) Then
!!$        expnnt = -DEFAULT_MAX_EXP
!!$      Else If (expnnt > DEFAULT_MAX_EXP) Then
!!$        expnnt = DEFAULT_MAX_EXP
!!$      End If
!!$      ! old way of getting wts
!!$      oldwt  = Exp(expnnt)/bmap%cum_wt

      ! new way of getting wts
      newwt  = bmap%cube(index)%wt / bmap%cum_wt
!!$      If ((Abs((newwt-oldwt)/oldwt)) > 0.000001 ) Then
!!$        Write(*,*) "wt from new", newwt, oldwt
!!$        Stop
!!$      Endif
    End If
    
    bmap_getcellwt = newwt
    Return
  End Function bmap_getcellwt

  !----------------------------------------------------
  ! Get the bias weight of the cubelet indexed "index"
  !----------------------------------------------------
  Real(kind=RDbl) Function bmap_getbiaswt(bmap, index)
    Type(BiasMap_Params), Intent(in) :: bmap
    Integer, Intent(in)  :: index
    bmap_getbiaswt = bmap%cube(index)%wt/bmap%cum_wt

  End Function bmap_getbiaswt

  !------------------------------------------------
  ! Returns the center of mass of the bmap_cube 
  !------------------------------------------------
  Type(VecType) Function bmap_getcellCOM(bmap, index)
    Type(BiasMap_Params), Intent(in) :: bmap
    Integer, Intent(in)  :: index
    Real(kind=RDbl)   :: x, y, z

    x = bmap%cube(index)%x + bmap%stepx/2
    y = bmap%cube(index)%y + bmap%stepy/2
    z = bmap%cube(index)%z + bmap%stepz/2

    bmap_getcellCOM= (/ x, y, z /)
  
  End Function bmap_getcellCOM


  !------------------------------------------------
  ! Returns a point in the cube indexed by "index"
  !------------------------------------------------
  Type(VecType) Function bmap_getbiaspt(bmap, index)
    Type(BiasMap_Params), Intent(in) :: bmap
    Integer, Intent(in)  :: index
    Real(kind=RDbl)   :: x, y, z
    Type(VecType)                       :: position, positionpbc

    x = bmap%cube(index)%x
    y = bmap%cube(index)%y
    z = bmap%cube(index)%z

    x = x + rranf()*bmap%stepx
    y = y + rranf()*bmap%stepy
    z = z + rranf()*bmap%stepz

    ! we need to take care of the low probable event that rranf() returns 1.00
    ! and the x + stepx becomes just above the bmap edgelength because of 
    ! roundoff errors
    If(gorthoflag) Then
    If ( x >= bmap%eff(1) ) x = 0.9999990 *  bmap%eff(1)
    If ( y >= bmap%eff(2) ) y = 0.9999990 *  bmap%eff(2)
    If ( z >= bmap%eff(3) ) z = 0.9999990 *  bmap%eff(3)
    Else
    position%comp(1) = x
    position%comp(2) = y
    position%comp(3) = z
    Call simcell_pbcptFund(bmap%simcell, position, positionpbc)
    x = positionpbc%comp(1)
    y = positionpbc%comp(2)
    z = positionpbc%comp(3)
    End If
 
    bmap_getbiaspt = (/x, y, z/)
    Return
  End Function bmap_getbiaspt

  !-----------------------------------------------
  ! Randomly picks a cubelet according to its bias weight
  ! It returns the index of the cubelet.
  !-----------------------------------------------
  Integer Function bmap_getbiasindex(bmap)
    Type(BiasMap_Params), Intent(in) :: bmap
    
    Integer   :: low, high, mid, n
    Real(kind=Rdbl)  :: key, cum_prob
    Real(kind=Rdbl)  :: wn, wn_1

    low = 1
    high = bmap%numtab

    !** rranf returns RDbl, but key(Rdbl) has only first 1/2 of the sign.digits
    key  = rranf()

    !** Do a binary search to find the first index with a cumulative
    !** probability higher than key
    Do
      If (low > high) Exit
      mid = (low + high)/2
      cum_prob = bmap%cube(mid)%cum_prob
      If (key < cum_prob) Then
        high = mid - 1
      Else If (key > cum_prob) Then
        low = mid + 1
        mid = low
      Else If (key == cum_prob) Then
        Exit

      Endif
    Enddo
    

    If (mid>1) Then
    wn= bmap%cube(mid)%cum_prob
    wn_1=bmap%cube(mid-1)%cum_prob
    If ( (key>wn) .Or. (key<wn_1) ) Then
      Write(*,*) "KEY ",  key
      Write(*,*) "wt(n), wt(n-1)", wn, wn_1
      Write(*,*) low, mid, high
      Stop
    Endif
  Endif

    !** low is equal to high
    bmap_getbiasindex = mid
!    nrg= -Log(bmap%cube(mid)%wt)/bmap%rti
!!$    If (bmap%cube(mid)%wt<5) Then
!!$      Write(*,*) low, mid, high
!!$      Write(*,*) wn, wn_1
!!$      Write(*,*) bmap%cube(mid)%wt, bmap%cube(mid-1)%wt
!!$      Write(*,*) key
!!$      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$      Stop
!!$    Endif
!!$    freqArr(mid) = freqArr(mid)+1
!!$    If (genparams%currentiteration==genparams%niterations) Then
!!$      Do n=1,bmap%numtab
!!$!!        Write(*,*) freqArr(n), bmap%cube(n)%wt
!!$      End Do
!!$    Endif
!!$    Return

  End Function bmap_getbiasindex

  !------------------------------------------------
  ! Cleanup the pmap pointers
  !------------------------------------------------
  Subroutine bmap_cleanup(bmap)
    Type(BiasMap_Params), Intent(inout) :: bmap
    Integer    :: error(2), i

    Call mapheader_cleanup(bmap%header)
    If (Associated(bmap%ptrb)) Then
      Deallocate(bmap%ptrb, STAT=error(1))
    Endif
    If (Associated(bmap%rb)) Then
      Deallocate(bmap%rb, STAT=error(2))
    Endif
    Do i=1, 2
      If (error(i) /= 0) Then
        Write(0,'(1x,2a,i4, a, i4)') __FILE__," : ",__LINE__, &
            "Error in deallocation", i
        Stop
      End If
    End Do
  End Subroutine bmap_cleanup

  !-----------------------------------------------------------------
  ! Displays the bmap as an xyz file for display using rasmol
  ! A new grid is made and values ar evaluated at the center of the cubes 
  ! in the New grid.
  ! Requires - 
  !          bmap     - the map
  !          gridsize - gridsize( nx, ny, nz) of the new grid
  !          logscale - whether use a log scale or not
  !          minpc, maxpc - any point which has a value between minpc % or 
  !                       maxpc % will be, displayed. 0% refers to the 
  !                       minimum, 100% refers to the maximum
  !          xyzfilename - name of rasmol file 
  !          descr       - small description to be written to rasmol file
  !          atsym       - symbol of atom(ex: "C" for carbon, "CH3" for methyl
  !-----------------------------------------------------------------
  Subroutine bmap_makeXYZ( bmap, gridsize, logflag, minpc, maxpc, &
      xyzfilename, descr, atsym)
    Type(BiasMap_Params), Intent(in)  :: bmap
    Integer, Dimension(3), Intent(in) :: gridsize
    Logical, Intent(in)              :: logflag
    Real(kind=RDbl)                   :: minpc, maxpc
    Character(len=strLen)             :: xyzfilename, descr
    Character(len=3)                  :: atsym
    Type(VecType) :: point, ucpoint, outpoint
    Real(kind=RDbl),Dimension(:,:,:), Pointer :: xyz
    Logical,Dimension(:,:,:), Pointer :: untabulated
    Real(kind=RDbl)                           :: lx,ly,lz, x, y, z
    Real(kind=RDbl)                           :: gridDx, gridDy, gridDz, wt
    Real(kind=RDbl)                           :: minwt, maxwt, minval, maxval
    Integer                                   :: error,unitno, i, j, k, numpts
Integer :: li, lj, lk

    Write(*,*) "Displaying bmap as an xyz file into file : "//Trim(xyzfilename)
    unitno=file_open(xyzfilename)

    Allocate(xyz(gridsize(1),gridsize(2),gridsize(3)), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"bmap_xyz")
    Allocate(untabulated(gridsize(1),gridsize(2),gridsize(3)), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"bmap_untabulated")
    untabulated=.True.

    lx=(bmap%header%stepx) * (bmap%header%nbrx-1)
    ly=(bmap%header%stepy) * (bmap%header%nbry-1)
    lz=(bmap%header%stepz) * (bmap%header%nbrz-1)
    Write(*,'(a,3f12.5)') "Edge lengths of the cube ", lx, ly, lz

    !** dx , dy and dz for the grid
    gridDx=(lx/gridsize(1))
    gridDy=(ly/gridsize(2))
    gridDz=(lz/gridsize(3))

    !** Loop thru all points and find the minimum and maximum
    minwt=1.000E40_RDbl
    maxwt=-minwt
    Do i=1,gridsize(1)
      x=i*gridDx - gridDx/2
      Do j=1,gridsize(2)
        y=j*gridDy - gridDy/2
        Do k=1,gridsize(3)
          z=k*gridDz - gridDz/2
          ! wt= -1 if the cell is not tabulated
          wt=bmap_getcellwt(bmap,(/x,y,z/), .False.)
          If (wt< zero) Then
            untabulated(i,j,k)=.True.
          Else
            untabulated(i,j,k)=.False.
            If (wt>maxwt) maxwt=wt
            If (wt<minwt) minwt=wt
            xyz(i,j,k) = wt
          Endif
        End Do
      End Do
    End Do

    !** decide the bounds
    If (logflag) Then
      If (minwt<zero) Then
        Write(*,*) "can't use logscale, some how there are &
            &-ve values in bmap"
        Stop
      Endif
      minval=minwt* Exp((minpc/100)*Log(maxwt/minwt))
      maxval=minwt* Exp((maxpc/100)*Log(maxwt/minwt))
      Write(*,*) "using logscale"
      Write(*,*) "minwt :", minwt,  " maxwt :", maxwt
      Write(*,*) "minval:", minval, " maxval:", maxval
      Write(*,*) "wt is exp(-beta U )"
    Else
      minval=minwt+(minpc/100)*(maxwt-minwt)
      maxval=minwt+(maxpc/100)*(maxwt-minwt)
      Write(*,*) "using arithmetic-scale"
      Write(*,*) "minwt:", minwt, " maxwt:" ,maxwt
      Write(*,*) "minval:", minval, " maxval:", maxval
      Write(*,*) "wt is exp(-beta U )"
    Endif

    !** Loop thru again this time couting number of points, 
    !** and seeing how big the file is going to be 
    !** xyz file
    numpts=0

    Do i=1,gridsize(1)
      Do j=1,gridsize(2)
        Do k=1,gridsize(3)
          If (.Not.(untabulated(i,j,k))) Then
            If ((xyz(i,j,k) >= minval) .And. (xyz(i,j,k) <= maxval)) Then
              numpts=numpts+1
            Endif
          Endif
        End Do
      End Do
    End Do
    Write(*,*) "There will be ", numpts, " points written into the file; make &
        & sure you have enough disk space"

    !** Finall write it
    Write(unitno,'(i6)') numpts
    Write(unitno,'(a)') "from_bmap : "//Trim(descr)
    Do i=1,gridsize(1)
      Do j=1,gridsize(2)
        Do k=1,gridsize(3)
          If (.Not.(untabulated(i,j,k))) Then
            If ((xyz(i,j,k)>=minval) .And. (xyz(i,j,k) <= maxval)) Then
              point=(/ (i*gridDx-gridDx/2),  &
                  (j*gridDy-gridDy/2), (k*gridDz-gridDz/2) /)
              ucpoint=simcell_maptouc(bmap%simcell,point)

              Do li=0,0
                Do lj=0,0
                  Do lk=0,0
              outpoint=ucpoint+simcell_getcellorigin(bmap%simcell, li, lj, lk)
                    Write(unitno,'(a,3f15.6)') "Ne", outpoint
                  end do
                end do
              end do
!!$              Write(unitno,'(a,3f15.6)') "Ne", (i*gridDx-gridDx/2), &
!!$                  (j*gridDy-gridDy/2), (k*gridDz-gridDz/2)
            Endif
          Endif
        End Do
      End Do
    End Do

    Deallocate(xyz)
    Deallocate(untabulated)

    Close(unitno)

  End Subroutine bmap_makeXYZ
  
  !------------------------------------------------------------------
  ! Display the bias map information
  ! Requires:  bmap -- bias map parameters
  !            displayheader -- if True, display bmap header also
  !            indent -- indent
  !            optunit -- unit number to dump into
  !------------------------------------------------------------------
  Subroutine bmap_display(bmap, displayheader, indent, optunit)
    Type(BiasMap_Params), Intent(In) :: bmap
    Logical                          :: displayheader
    Integer, Intent(In)              :: indent 
    Integer, Optional, Intent(In)    :: optunit

    Integer                     :: i, unitno
    Character(len=indent)       :: blank

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    blank = Repeat(' ',indent)

    Write(unitno,'(a,2x,3a)') blank, "Bias map filename: ",Trim(bmap%filename)
    Write(unitno,'(a,2x,a,g14.8)') blank, "Cumulative Wt: ", bmap%cum_wt
    Write(unitno,'(a,2x,a,i8)') blank, "No. of tabulated cubelets: ", &
        bmap%numtab
    Write(unitno,'(a,2x,a,f12.5)') blank, "Cutoff Potential(kcal): ", &
        bmap%ucutoff
    If (displayheader) Call mapheader_display(bmap%header, indent+2, unitno)

  End Subroutine bmap_display

End Module bmap




