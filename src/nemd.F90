!------------------------------------------------------------------------------
! This module initializes the variables required by NEMD and runs NEMD 
! calculations
! Refer to Maginn et al for details
! UNITS : Color Field is read from control file in Kcal/mol/ang/color
!------------------------------------------------------------------------------
Module nemd
  Use config, Only: config_getconfig, AtMolCoords, config_getnatoms, &
      config_getnmoles, config_isfixed, config_getSpcCOMVel, &
      config_getSpcDensity, config_getMolecCOM
  Use defaults, Only: RDbl, strLen, STATS_BLOCKSIZE, one, zero, scalef, &
      MAX_SORBS
  Use file, Only: file_open
  Use molecules, Only: molecules_getnatoms, molecules_AtomMass, &
      molecules_getatype, molecules_getnsorbs, molecules_getmass, &
      molecules_gettype
  Use simcell, Only: SimCell_Params, simcell_getell
  Use stats, Only: Statistics, stats_getvalue, stats_getcavg, stats_getstd, &
      stats_init, stats_update
  Use utils, Only: stripcmnt, split, tolower, toupper, filesrchstr, &
      allocErrDisplay, DeallocErrDisplay, int2str, real2str, toreal, &
      cleanstring
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag


  Implicit None
  Save

  Private
  Public :: nemd_init, nemd_display, nemd_sampleCF, nemd_clean, NEMD_IS_ON, &
      NEMDParams, NEMDInfo, nemd_updateNrg, nemd_getpe, nemd_updateflux, &
      nemd_simdisplay
  Type NEMDInfo
    Real(kind=RDbl)                 :: T, totalcharge
    Real(kind=RDbl),Dimension(3)    :: field
    Real(kind=RDbl),Dimension(:), Pointer   :: color
    Type(VecType),Dimension(:), Pointer     :: force
    Real(kind=RDbl)                 :: current_nrg
    Type(Statistics)                :: pe
    Type(Statistics),Dimension(:,:),Pointer :: Vflux
    Integer                         :: Nflux
    Integer                         :: sorb,fluxdirection 
    ! fluxdirection : 1=x, 2=y, 3=z
    Character(len=strLen)           :: sorbname

    ! info required for calculating distance traveled
    Integer                         :: ntotmols
    Type(VecType), Dimension(:), Pointer :: com_0 ! a linear array for all 
    ! molecs of all species excluding the fixed species
  End Type NEMDInfo

  Type(NEMDInfo) , Pointer :: NEMDParams
  Logical                  :: NEMD_IS_ON=.False.

  Character(len=strLen) :: NEMD_tag= " NEMD Details "
  Integer, Parameter :: NSKIP_fluxupdate=4
  ! skips so many before updating fluxes
  Integer :: fluxdiff
Contains

  !----------------------------------------------------------------------------
  ! checks whether NEMD is required, Initializes NEMD
  ! Requires:  ctrlfile -- main control file 
  ! This works only with constant number of molecules
  ! other wise the flux calcualtions might lead to problems
  !----------------------------------------------------------------------------
  Subroutine nemd_init(ctrlfile,species)
    Character(*), Intent(in) :: ctrlfile
    Type(AtMolCoords), Dimension(:), Intent(In) :: species

    Integer               :: a, spc, m, i, error, nfields,unitno, useNEMD
    Integer               :: sorbno, natoms, nspc, ntotmols, m_index
    Real(kind=RDbl)       :: totalcharge, totalmass, amass,  mfac, color, max
    Character(len=strLen) :: text, sorbname
    Character(len=strLen), Dimension(strLen) :: fields

    unitno=file_open(ctrlfile)
    useNEMD= filesrchstr(unitno,NEMD_tag,text,.True.)

    If (useNEMD==0) Then
      NEMD_IS_ON=.False.
      Nullify(NEMDParams)
      Return
    else
      NEMD_IS_ON=.True.
      Allocate(NEMDParams,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Write(*,*) "Initializing the NEMD variables by reading from the ctrlfile"
      Read(unitno,*) NEMDParams%field
      !      Read(unitno,*) NEMDParams%color
      Read(unitno,'(a)') text
    Endif

    ! process the string and get all charges
    ! decide which is the nemd-sorb
    text=cleanstring(text)
    nfields=split(text,fields,",")
    sorbname=cleanstring(fields(1))
    sorbno=molecules_gettype(sorbname)
    natoms=molecules_getnatoms(sorbno)
    NEMDParams%sorbname=sorbname
    NEMDParams%sorb=sorbno
    If (config_getnmoles(species,sorbno)==0) Then
      Write(*,*) "You are applying the field on a sorb with zero moles"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    Allocate(NEMDParams%color(natoms),NEMDParams%force(natoms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    If (Trim(cleanstring(fields(2)))=="COM") Then
      totalcharge=toreal(fields(3))
      totalmass=molecules_getmass(sorbno)
      Write(*,*) totalcharge, totalmass
      Do a=1,natoms
        amass= molecules_AtomMass(sorbno, a)
        mfac=amass/totalmass
        color=mfac*totalcharge
        NEMDParams%color(a)=color
        Write(*,*) amass, mfac, color
      End Do
    Else
      If (nfields/=(natoms+1)) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Not enough charges on the NEMD line for sorb : "//&
            Trim(sorbname)
        Write(*,*) "Number of charges specified should be equal to number &
            &of atoms : ", natoms 
        Stop
      Endif
      totalcharge=zero
      Do a=1,natoms
        color=toreal(fields(a+1))
        NEMDParams%color(a)=color
        totalcharge=totalcharge+color
      End Do
    Endif
    NEMDParams%totalcharge=totalcharge


    ! F = -dV/dr
    ! Units field = Kcal/mol/ang./color
    !       force = Kcal/mol/ang.
    ! make force an natoms-array
    Do a=1,natoms
      Do i=1,3
        NEMDParams%force(a)%comp(i)= -NEMDParams%field(i)*NEMDParams%color(a)
      End Do
    End Do
    Call stats_init(NEMDParams%pe,"NEMD Total PE",STATS_BLOCKSIZE, &
        .False.,'f10.3')
    nspc=molecules_getnsorbs()
    Allocate(NEMDParams%VFlux(3,nspc),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    
    Do i=1,nspc
      Call stats_init(NEMDParams%VFlux(1,i),"NEMD x-VelocityFlux", &
          STATS_BLOCKSIZE, .False.,'f10.3')
      Call stats_init(NEMDParams%VFlux(2,i),"NEMD y-VelocityFlux", &
          STATS_BLOCKSIZE, .False.,'f10.3')
      Call stats_init(NEMDParams%VFlux(3,i),"NEMD z-VelocityFlux", &
          STATS_BLOCKSIZE, .False.,'f10.3')
    End Do
    NEMDParams%NFlux=0

    ! decide which is the main flux direction
    ! L_ij's will be reported only in this direction
    max=zero
    Do i=1,3
      If (Abs(NEMDParams%field(i))>max) Then
        max=Abs(NEMDParams%field(i))
        NEMDParams%fluxdirection=i
      Endif
    Enddo

    !*** initialize com arrays store initial number of molecules

    ! count number of molecules
    ntotmols=0
    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle
      ntotmols=ntotmols+config_getnmoles(species,spc)
    End Do
    NEMDParams%ntotmols=ntotmols
    Allocate(NEMDParams%com_0(ntotmols),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    ! fill up the array
    m_index=0
    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle
      Do m=1,config_getnmoles(species,spc)
        m_index=m_index+1
        NEMDParams%com_0(m_index)=config_getMolecCOM(species,spc,m)
      End Do
    End Do

  End Subroutine nemd_init


  !----------------------------------------------------------------------------
  ! Evaluate NEMD color field potentials nrg in Kcal/mol
  ! this value helps us to make sure that the energy is conserved
  !----------------------------------------------------------------------------
  Subroutine nemd_updateNrg(species)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Integer :: nspc, spc, natoms, a, nmoles, m
    Real(kind=RDbl) :: nrg, pos_nrg
    Type(VecType)          :: posvec

    nspc=molecules_getnsorbs()
    nrg=zero

    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle
      If (.not.(NEMDParams%sorb==spc) ) Cycle 
      natoms = molecules_getnatoms(spc)
      nmoles = config_getnmoles(species,spc)
      Do a=1,natoms
        Do m=1,nmoles
          ! current position
          posvec=species(spc)%coords(a,m)%rp

          !nrg at this position
          pos_nrg=nemd_vecNrg(posvec)

          ! multiply by color charge
          nrg=nrg+pos_nrg*NEMDParams%color(a)
        end do
      end do
    End Do
    ! NRG in units of Kcal/mol
    NEMDParams%current_nrg=nrg   
    Call stats_update( NEMDParams%pe,nrg)


  End Subroutine nemd_updateNrg

  !----------------------------------------------------------------------------
  ! Updates the flux once in every Nskip steps (usually once in 20 steps)
  ! Two types of Flux definitions are there. VFlux and NFlux
  ! VFlux : based on average velocity
  ! NFlux : this is based on number of molecules that cross the middle plane
  ! NFlux is just for debug-purposes,(ideally NFlux=VFlux)
  ! Requires :
  ! species - the configuration
  ! maxdist - estimate of maximum distance traveled
  !----------------------------------------------------------------------------
  Subroutine nemd_updateFlux(species,scell,maxdist)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: species
    Type(SimCell_Params), Intent(In) :: scell
    Real(kind=RDbl),Intent(in) :: maxdist
    Integer :: nspc, spc, nmols1, m, error, LRcount, RLcount, now
    Integer,Save :: count=1
    Integer, Parameter :: LEFT=0,RIGHT=1,FAR=2
    Logical,Save :: first_time=.True.
    Integer, Dimension(:),Pointer,Save:: pos_index
    Real(kind=RDbl),Dimension(3),Save :: ell
    Real(kind=RDbl),Save :: diagonal
    Real(kind=RDbl)      :: newpos
    Type(VecType) :: comvel


    ! decide whether to update or not
    If (count< NSKIP_fluxupdate) Then
      count=count+1
      Return
    Else
      count=1
    Endif

    If (first_time) Then
      first_time=.False.
      ell=simcell_getell(scell)
      diagonal=Sqrt(ell(1)*ell(1)+ell(2)*ell(2)+ell(3)*ell(3))
      nmols1=config_getnmoles(species,1)
      Allocate(pos_index(nmols1),STAT=error)
    Endif

    ! check how much the molecules have travled on an average
    If (maxdist>diagonal/5) Then
      Write(*,*) "maxdist, diagonal : ",maxdist, diagonal 
      Write(*,*) " Need to update flux more often"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif



    ! update VFLux, based on average velocities
    ! so it is not flux, it is the average velocity
    nspc=molecules_getnsorbs()
    Do spc=1,nspc

      If (config_isfixed(species(spc))) Cycle

      comvel=config_getSpcCOMVel(species,spc)

      Call stats_update( NEMDParams%VFlux(1,spc),comvel%comp(1)) ! x
      Call stats_update( NEMDParams%VFlux(2,spc),comvel%comp(2)) ! y
      Call stats_update( NEMDParams%VFlux(3,spc),comvel%comp(3)) ! z
    End Do

    ! update NFlux only in x direction, based on x=ell(1)/2 plane
    ! counts how many have %r < ell(1)/2
    spc=1 ! only for 1st species
    nmols1=config_getnmoles(species,spc)
    LRcount=0
    RLcount=0
    Do m=1,nmols1
      newpos=species(spc)%coords(1,m)%r%comp(1) ! x posn
      If (newpos<ell(1)/2) Then
        If (newpos>ell(1)/4) Then
          now=LEFT
        else
          now=FAR
        Endif
      else
        If (newpos<(3*ell(1)/4)) Then
          now=RIGHT
        else
          now=FAR
        Endif
      endif
      If ((now==LEFT).And.(pos_index(m)==RIGHT)) RLcount=RLcount+1
      If ((now==RIGHT).And.(pos_index(m)==LEFT)) LRcount=LRcount+1
      pos_index(m)=now

    End Do

    ! total molecules that ever crossed x=ell(1)/2
    NEMDParams%NFlux =  NEMDParams%NFlux+LRcount-RLcount 
    fluxdiff=LRcount-RLcount 
  End Subroutine nemd_updateFlux


  !----------------------------------------------------------------------------
  ! Returns the total exrernal field potential energy type requested by sType
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function nemd_getpe(sType)
    Character(*), Intent(In)         :: sType

    Select Case(tolower(sType))

    Case ('inst')
      nemd_getpe = stats_getvalue(NEMDParams%pe)
    Case ('cavg')
      nemd_getpe = stats_getcavg(NEMDParams%pe)
    Case ('std')
      nemd_getpe = stats_getstd(NEMDParams%pe)
    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select

  End Function nemd_getpe

  !----------------------------------------------------------
  ! Return the nrg given a position vector, in "kcal/mol/color"
  !----------------------------------------------------------
  Real(kind=RDbl) Function nemd_vecNrg(vec)
    Type(VecType),Intent(in)           :: vec
    nemd_vecNrg=vec%comp(1)*NEMDParams%field(1) + &
        vec%comp(2)*NEMDParams%field(2) + &
        vec%comp(3)*NEMDParams%field(3) 
  End Function nemd_vecNrg

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine nemd_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a,t30,a)') '------','Lazy Programmer. fill in pleeeese'
  End Subroutine nemd_sampleCF


  !----------------------------------------------------------------------------
  ! Displays information about NEMD details during simulation
  !----------------------------------------------------------------------------
  Subroutine nemd_simdisplay(species,scell,time,unitno,nskip)
    Integer, Intent(In)              :: unitno
    Integer, Intent(In)              :: nskip
    Real(kind=RDbl), Intent(In)  :: time 
    Type(SimCell_Params), Intent(In)            :: scell
    Type(AtMolCoords), Dimension(:), Intent(In)    :: species
    Type(VecType) :: new_com, diff_vec
    Integer               :: dUnit,nspc,spc,i,m,nmols,nflux,dir,m_index
    Character(len=nskip)  :: blank
    Character(len=20)  :: str1
    Character(len=100)  :: Lstring
    Real(kind=RDbl) :: density, area, flux, charge, field
    Real(kind=RDbl),Dimension(3) :: avgvels,ell
    Real(kind=RDbl),Dimension(MAX_SORBS) :: L

    dUnit = unitno
    blank=Repeat(" ", nskip)
    Write(dUnit,'(2a)')blank, "NEMD Sim Details : "
    nspc=molecules_getnsorbs()
    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle
      Do i=1,3
        ! avg. vel. : these have units of velocity ang/ps
        avgvels(i)=stats_getcavg(NEMDParams%VFlux(i,spc))
      end do
      str1=int2str(spc)
      Write(dUnit,'(2a,3e12.3,a)') blank,&
          "Vel. of compound :"//Trim(str1), avgvels, "  [ang./ps]"
    End Do

    ! calculate L_ij s
    dir=NEMDParams%fluxdirection ! main direction in which we do NEMD
    Lstring=Repeat(" ",100)
    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle
      Do i=1,3
        ! avg. vel. : these have units of velocity ang/ps
        avgvels(i)=stats_getcavg(NEMDParams%VFlux(i,spc))
      end do
      density=config_getSpcDensity(species,scell,spc) ! molec/ang^3
      flux= avgvels(dir)*density                        ! molec/ps/ang^2
      charge=NEMDParams%totalcharge ! this is expected to be equal to 1
      field=NEMDParams%field(dir)   ! kcal/mol/ang/colorcharge
      L(spc)=flux/(-one*charge*field)    ! (molec/ps/ang)/(kcal/mol)
      Lstring=Trim(Lstring)//"   "//Trim(real2str(L(spc)))
    End Do

    nflux= NEMDParams%NFlux
    ell=simcell_getell(scell)
    area=ell(1)*ell(3)
    Write(dUnit,'(2a,1e12.3,a)')blank,"Number based x-flux of spc-1", &
        nflux/area/time,           "  [molec/ang.^2/ps]"
    Write(dUnit,'(2a,i12)') blank, &
        " Number of molecules that crossed yz plane : ", nflux
    Write(dUnit,'(2a,i2,a)')blank, "L_ij coeff.s for different i, &
        &(j=NEMDsorb=",NEMDParams%sorb,") [molec/ang/ps]/[kcal/mol]"
    Write(dUnit,'(3a)')blank, "  :  ", Trim(Lstring)

    ! calculate fluxes/velocities based on total distance traveled
    m_index=0

    Write(dUnit,'(2a)')blank,"-- flux from total species displacements ---"
    Do spc=1,nspc
      If (config_isfixed(species(spc))) Cycle
      diff_vec=zero ! initialize total species displacement to zero
      nmols=config_getnmoles(species,spc)
      If (nmols==0) Cycle
      Do m=1,nmols
        m_index=m_index+1
        new_com=config_getMolecCOM(species,spc,m)
        diff_vec=diff_vec+(new_com-NEMDParams%com_0(m_index))
      End Do
      avgvels=(diff_vec)/(time*nmols)
      str1=int2str(spc)
      Write(dUnit,'(2a,3e12.3,a)') blank,&
          "Vel. of compound :"//Trim(str1), avgvels, "  [ang./ps]"
    End Do

  End Subroutine nemd_simdisplay


  !----------------------------------------------------------------------------
  ! Displays information about NEMD parameters
  !----------------------------------------------------------------------------
  Subroutine nemd_display(unitno,nskip)
    Integer, Intent(In)              :: unitno
    Integer, Intent(In)              :: nskip

    Integer               :: dUnit,a, natoms
    Character(len=nskip)  :: blank
    Real(kind=RDbl) :: mass
    Type(VecType) :: force, accel

    dUnit = unitno
    blank=Repeat(" ", nskip)
    Write(dUnit,'(2a)')blank, "NEMD Details : "
    Write(dUnit,'(2a,t40,a,3e10.3)')blank, &
        "Potential Coefficients (x,y,z) Kcal/mol/ang.",":",NEMDParams%field 
    Write(dUnit,'(2a,i2,a)') blank, &
        "fluxes and L_ij will be reported in direction-",&
        NEMDParams%fluxdirection, ". (x=1, y=2, z=3)"
    Write(dUnit,'(2a,t40,3a,i4)')blank, "Name and number of sorbate",&
        ":",Trim(NEMDParams%sorbname)," ", NEMDParams%sorb
    Write(dUnit,'(2a,t40,a)')blank, &
        "Charges/force/acceleration for each atom ",":"
    Write(dUnit,'(4a)') blank,"Atom Num.  Charge      Force         Accel."
    Write(dUnit,'(4a)') blank,"           Color       kcal/ang/mol  ang/ps^2"
    natoms=molecules_getnatoms(NEMDParams%sorb)
    Do a=1,natoms
      force=NEMDParams%force(a)
      mass=molecules_AtomMass(NEMDParams%sorb,a)
      ! onvert accel to ang/ps^2
      accel=(force/mass)*scalef
      Write(dUnit,'(a,i8,1f10.5,1e15.4,1f15.4)') blank,a,NEMDParams%color(a),&
          mag(force),mag(accel)
    End Do
    Write(dUnit,*) 
  End Subroutine nemd_display

  !----------------------------------------------------------------------------
  ! Cleans up the NEMD structure
  ! Requires:  nemd -- the nemd structure
  !----------------------------------------------------------------------------
  Subroutine nemd_clean()
    Integer :: error
    If (NEMD_IS_ON) Then
      !** Nothing
      If (Associated(NEMDParams)) Then
        Deallocate(NEMDParams%force, NEMDParams%color,STAT=error)
        If (error/=0) Call DeAllocErrDisplay(__FILE__,__LINE__)
        Nullify(NEMDParams)
      End If
    End If

  End Subroutine nemd_clean

End Module nemd

