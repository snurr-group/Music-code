!-------------------------------------------------------------------------------
! This module contains data types designed to store species-specific statistics
! 
! Needed Improvements:
! 1) finish requires comments
!-----------------------------------------------------------------------------

Module storestats

  Use defaults, Only: RDbl, strLen, lstrLen,STATS_BLOCKSIZE,dashedline,scalepe, &
      STRETCH_INDEX, BENDING_INDEX, TORSION_INDEX, INTRA_NRG_BLOCKSIZE, &
      CONSTRAINT_INDEX, INTRAPAIR_INDEX, INTRACOUL_INDEX, TOTAL_INDEX, &
      EXTERNAL_INDEX, NO_OF_INTRA_POTS
  Use utils, Only: toupper, int2str, real2str, &
      allocErrDisplay, deallocErrDisplay
  Use stats, Only: stats_getvalue, Statistics, stats_getblock, stats_getcavg, &
      stats_getstd, stats_init, stats_update,  stats_getfield
  Use molecules, Only: molecules_name,molecules_getnsorbs
  Use storebase, Only: EnergyPlus

  Implicit None
  Save

  Private
  Public :: IntraEnergies, Species_Stats, storestats_init, &
      storestats_getke, storestats_getcoul, storestats_getnoncoul, &
      storestats_getintranrg, storestats_initIntraEnergy, &
      storestats_updatetemp, storestats_updatenrg, storestats_gettemp, &
      storestats_updateEnergySS, storestats_displaynrg, storestats_clean, &
      storestats_incrcoulnrg, storestats_incrnoncoulnrg, &
      storestats_incrIntraNrg, storestats_updateIntraEnergy, &
      storestats_incrAllIntraNrg, storestats_incrkinnrg, &
      storestats_getallintra, storestats_displayavgs

  Type IntraEnergies
    Type(Statistics)     :: bs,bb,tor,ip,ic,con,ext,total
  End Type IntraEnergies
    
  Type Species_Stats
    Integer                                        :: nmolec_types 
    Logical                                        :: usampling
    Type(Statistics), Dimension(:,:), Pointer      :: noncoulnrg
    Type(Statistics), Dimension(:,:), Pointer      :: coulnrg
    Type(IntraEnergies), Dimension(:), Pointer     :: intranrg
    Type(Statistics), Dimension(:), Pointer        :: tmole
    Type(Statistics), Dimension(:), Pointer        :: tatom
    Type(Statistics), Dimension(:), Pointer        :: ke
  End Type Species_Stats                           

Contains
  !----------------------------------------------------------------------------
  ! Initialize the molecule information
  ! Requires:  spcstats -- storage structure for species statistics
  !            usample -- True indicates that umbrella sampling is being used
  !----------------------------------------------------------------------------
  Subroutine storestats_init(spcstats,usample)
    Type(Species_Stats), Intent(Out)      :: spcstats
    Logical, Intent(In)                   :: usample

    Integer                               :: i,j,ios,error, nspcs

    spcstats%usampling = usample
    spcstats%nmolec_types = molecules_getnsorbs()

    !** Initialize the coulombic and non-coulombic energy structures
    nspcs = spcstats%nmolec_types
    Allocate (spcstats%noncoulnrg(nspcs, nspcs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'spcstats%ncoulnrg')
    Allocate (spcstats%coulnrg(nspcs, nspcs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'spcstats%coulnrg')
    Allocate(spcstats%intranrg(nspcs),Stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'spcstats%intranrg')
    Allocate(spcstats%ke(nspcs),Stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'spcstats%ke')
    Allocate(spcstats%tatom(nspcs),Stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'spcstats%tatom')
    Allocate(spcstats%tmole(nspcs),Stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'spcstats%tmole')

    Do i = 1,nspcs
      Do j = 1,nspcs
        Call stats_init(spcstats%noncoulnrg(i, j),"Polarization Energy", &
            STATS_BLOCKSIZE,spcstats%usampling,'f10.3')
        Call stats_init(spcstats%coulnrg(i, j), "Coulombic Energy", &
            STATS_BLOCKSIZE,spcstats%usampling,'f10.3')
      End Do
      Call stats_init(spcstats%intranrg(i)%bs, &
          "Intramolecular Energy - Bond Stretching", STATS_BLOCKSIZE, &
          spcstats%usampling,'f10.3')
      Call stats_init(spcstats%intranrg(i)%bb, &
          "Intramolecular Energy - Bond Bending", STATS_BLOCKSIZE, &
          spcstats%usampling,'f10.3')
      Call stats_init(spcstats%intranrg(i)%tor, &
          "Intramolecular Energy - Torsion", STATS_BLOCKSIZE, &
          spcstats%usampling,'f10.3')
      Call stats_init(spcstats%intranrg(i)%ip, &
          "Intramolecular Energy - Intra Particle", STATS_BLOCKSIZE, &
          spcstats%usampling,'f10.3')
     Call stats_init(spcstats%intranrg(i)%ic, &
          "Intramolecular Energy - Intra Particle", STATS_BLOCKSIZE, &
          spcstats%usampling,'f10.3')
      Call stats_init(spcstats%intranrg(i)%con, &
          "Intramolecular Energy - Constraint", STATS_BLOCKSIZE, &
          spcstats%usampling,'f10.3')
      Call stats_init(spcstats%intranrg(i)%ext, &
          "Intramolecular Energy - External", STATS_BLOCKSIZE, &
          spcstats%usampling,'f10.3')
      Call stats_init(spcstats%intranrg(i)%total, &
          "Intramolecular Energy - Total", STATS_BLOCKSIZE, &
          spcstats%usampling,'f10.3')
          
      Call stats_init(spcstats%tatom(i),"Atomic Temperature", &
          STATS_BLOCKSIZE,spcstats%usampling,'f10.3')
      Call stats_init(spcstats%tmole(i),"Molecular Temperature", &
          STATS_BLOCKSIZE,spcstats%usampling,'f10.3')
      Call stats_init(spcstats%ke(i),"Kinetic Energy", &
          STATS_BLOCKSIZE,spcstats%usampling,'f12.4')
    End Do

  End Subroutine storestats_init

  !-----------------------------------------------------
  ! Increment the non-coulombic energy
  ! Requires:  spcstats -- storage structure for species statistics
  !-----------------------------------------------------
  Subroutine storestats_incrnoncoulnrg(spcstats,spc1, spc2, deltanrg)
    Type(Species_Stats), Intent(InOut)     :: spcstats
    Integer, Intent(in)         :: spc1, spc2
    Real(kind=RDbl), Intent(in) :: deltanrg

    Real(kind=RDbl)     :: oldnrg

    oldnrg = stats_getvalue(spcstats%noncoulnrg(spc1, spc2))
    Call stats_update(spcstats%noncoulnrg(spc1, spc2), oldnrg+deltanrg)
    Call stats_update(spcstats%noncoulnrg(spc2, spc1), oldnrg+deltanrg)
  End Subroutine storestats_incrnoncoulnrg

  !-----------------------------------------------------
  ! Increment the coulombic energy
  ! Requires:  spcstats -- storage structure for species statistics
  !-----------------------------------------------------
  Subroutine storestats_incrcoulnrg(spcstats,spc1, spc2, deltanrg)
    Type(Species_Stats), Intent(InOut)     :: spcstats
    Integer, Intent(in)         :: spc1, spc2
    Real(kind=RDbl), Intent(in) :: deltanrg

    Real(kind=RDbl)     :: oldnrg

    oldnrg = stats_getvalue(spcstats%coulnrg(spc1, spc2))
    Call stats_update(spcstats%coulnrg(spc1, spc2), oldnrg+deltanrg)
    Call stats_update(spcstats%coulnrg(spc2, spc1), oldnrg+deltanrg)
  End Subroutine storestats_incrcoulnrg

  !-----------------------------------------------------
  ! Increment the kinetic energy
  ! Requires:  spcstats -- storage structure for species statistics
  !-----------------------------------------------------
  Subroutine storestats_incrkinnrg(spcstats,spc, deltanrg)
    Type(Species_Stats), Intent(InOut)     :: spcstats
    Integer, Intent(in)         :: spc
    Real(kind=RDbl), Intent(in) :: deltanrg
    
    Real(kind=RDbl)     :: oldnrg
    
    oldnrg = stats_getvalue(spcstats%ke(spc ))
    Call stats_update(spcstats%ke(spc ), oldnrg+deltanrg)
  End Subroutine storestats_incrkinnrg
  
  !-----------------------------------------------------
  ! Increment the intramolecular energy
  ! WARNING: works only on total intra energy for now
  ! Requires:  spcstats -- storage structure for species statistics
  !-----------------------------------------------------
  Subroutine storestats_incrIntraNrg(spcstats,spc,deltanrg)
    Type(Species_Stats), Intent(InOut)     :: spcstats
    Integer, Intent(in)         :: spc
    Real(kind=RDbl), Intent(in) :: deltanrg
    Real(kind=RDbl)             :: oldnrg

    oldnrg = storestats_getintranrg(spcstats,spc,'inst')
    Call stats_update(spcstats%intranrg(spc)%total, oldnrg+deltanrg)

  End Subroutine storestats_incrIntraNrg

  !-----------------------------------------------------
  ! Increments all the intramolecular energy
  ! Requires:  spcstats -- storage structure for species statistics
  !-----------------------------------------------------
  Subroutine storestats_incrAllIntraNrg(spcstats,molec,deltanrg)
    Type(Species_Stats), Intent(InOut)     :: spcstats
    Integer, Intent(in)         :: molec
    Real(Kind=RDbl), Dimension(:), Intent(In) :: deltaNrg
    Real(kind=RDbl)             :: oldnrg

    oldnrg = storestats_getintranrg(spcstats,molec,'inst','STRETCH')
    Call stats_update(spcstats%intranrg(molec)%bs, &
        oldnrg+deltanrg(STRETCH_INDEX))

    oldnrg = storestats_getintranrg(spcstats,molec,'inst','BENDING')
    Call stats_update(spcstats%intranrg(molec)%bb, &
        oldnrg+deltanrg(BENDING_INDEX))

    oldnrg = storestats_getintranrg(spcstats,molec,'inst','TORSION')
    Call stats_update(spcstats%intranrg(molec)%tor, &
        oldnrg+deltanrg(TORSION_INDEX))

    oldnrg = storestats_getintranrg(spcstats,molec,'inst','INTRAPAIR')
    Call stats_update(spcstats%intranrg(molec)%ip, &
        oldnrg+deltanrg(INTRAPAIR_INDEX))

    oldnrg = storestats_getintranrg(spcstats,molec,'inst','INTRACOUL')
    Call stats_update(spcstats%intranrg(molec)%ic, &
        oldnrg+deltanrg(INTRACOUL_INDEX))

    oldnrg = storestats_getintranrg(spcstats,molec,'inst','CONSTRAINT')
    Call stats_update(spcstats%intranrg(molec)%con, &
        oldnrg+deltanrg(CONSTRAINT_INDEX))

    oldnrg = storestats_getintranrg(spcstats,molec,'inst','EXTERNAL')
    Call stats_update(spcstats%intranrg(molec)%ext, &
        oldnrg+deltanrg(EXTERNAL_INDEX))

    oldnrg = storestats_getintranrg(spcstats,molec,'inst','TOTAL')
    Call stats_update(spcstats%intranrg(molec)%total, &
        oldnrg+deltanrg(TOTAL_INDEX))

  End Subroutine storestats_incrAllIntraNrg

  !----------------------------------------------------------------------------
  ! Updates the ss energy
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Subroutine storestats_updateEnergySS(spcstats,spc1,spc2,eType,nrg)
    Type(Species_Stats), Intent(InOut)     :: spcstats
    Integer, Intent(In)         :: spc1, spc2
    Character(*), Intent(In)    :: eType
    Real(kind=RDbl)             :: nrg

    !SDEBUG
    ! this might be called a lot, using toupper not a good idea
    !    Select Case (ToUpper(eType))
    !SDEBUG

    Select Case (eType)
    Case ('NONCOUL','noncoul')
      Call stats_update(spcstats%noncoulnrg(spc1,spc2),nrg)
      If (spc1 /= spc2) Then
        Call stats_update(spcstats%noncoulnrg(spc2,spc1),nrg)
      End If

    Case ('COUL','coul')
      Call stats_update(spcstats%coulnrg(spc1,spc2),nrg)
      If (spc1 /= spc2) Then
        Call stats_update(spcstats%coulnrg(spc2,spc1),nrg)
      End If
    Case default
      Write(*,*) "could not find string : ", eType
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select

  End Subroutine storestats_updateEnergySS

  !----------------------------------------------------------------------------
  ! Updates the energies of a molecule type
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Subroutine storestats_updatenrg(spcstats,spc1,eType,nrg)
    Type(Species_Stats), Intent(InOut)     :: spcstats
    Integer, Intent(In)         :: spc1
    Character(*), Intent(In)    :: eType
    Real(kind=RDbl)             :: nrg

    Select Case (ToUpper(eType))
    Case ('KE')
      If (Associated(spcstats%ke)) Then
        Call stats_update(spcstats%ke(spc1),nrg)
      End If
    Case ('INTRA')
      Call stats_update(spcstats%intranrg(spc1)%total,nrg)
    Case ('CONST')
      Call stats_update(spcstats%intranrg(spc1)%con,nrg)
    End Select

  End Subroutine storestats_updatenrg

  !----------------------------------------------------------------------------
  ! Updates the intramolecular energy stats given molecule spc and intraNrg
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Subroutine storestats_updateIntraEnergy(spcstats,molec,intraNrg)
    Type(Species_Stats), Intent(InOut)        :: spcstats
    Integer, Intent(In)                       :: molec
    Real(Kind=RDbl), Dimension(:), Intent(In) :: intraNrg

    Call stats_update(spcstats%intranrg(molec)%bs, &
        intraNrg(STRETCH_INDEX))
    Call stats_update(spcstats%intranrg(molec)%bb, &
        intraNrg(BENDING_INDEX))
    Call stats_update(spcstats%intranrg(molec)%tor, &
        intraNrg(TORSION_INDEX))
    Call stats_update(spcstats%intranrg(molec)%ip, &
        intraNrg(INTRAPAIR_INDEX))
    Call stats_update(spcstats%intranrg(molec)%ic, &
        intraNrg(INTRACOUL_INDEX))
    Call stats_update(spcstats%intranrg(molec)%con, &
        intraNrg(CONSTRAINT_INDEX))
    Call stats_update(spcstats%intranrg(molec)%ext, &
        intraNrg(EXTERNAL_INDEX))
    Call stats_update(spcstats%intranrg(molec)%total, &
        intraNrg(TOTAL_INDEX))

  End Subroutine storestats_updateIntraEnergy

  !----------------------------------------------------------------------------
  ! Updates the temperature statistics
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Subroutine storestats_updateTemp(spcstats,spc,temp,tType)
    Type(Species_Stats), Intent(InOut)     :: spcstats
    Integer, Intent(In) :: spc
    Real(Kind=RDbl), Intent(In) :: temp
    Character(*), Intent(In) :: tType
    Select Case (tType)
    Case ('atomic')
      Call stats_update(spcstats%tatom(spc),temp)
    Case ('molecular')
      Call stats_update(spcstats%tmole(spc),temp)
    End Select
  End Subroutine storestats_updateTemp

  !----------------------------------------------------------------------------
  ! Initializes the statistics for the intramolecular energies
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Subroutine storestats_initIntraEnergy(spcstats,molec)
    Type(Species_Stats), Intent(InOut)     :: spcstats
    Integer, Intent(In) :: molec

    Call stats_init(spcstats%intranrg(molec)%bs,"Bond Stretch Energy",&
        INTRA_NRG_BLOCKSIZE,spcstats%usampling)
    Call stats_init(spcstats%intranrg(molec)%bb,"Bond Bending Energy",&
        INTRA_NRG_BLOCKSIZE,spcstats%usampling)
    Call stats_init(spcstats%intranrg(molec)%tor,"Bond Torsion Energy", &
        INTRA_NRG_BLOCKSIZE,spcstats%usampling)
    Call stats_init(spcstats%intranrg(molec)%ip,"Intra-pair Energy",&
        INTRA_NRG_BLOCKSIZE,spcstats%usampling)
    Call stats_init(spcstats%intranrg(molec)%ic,"Intra-coul Energy",&
        INTRA_NRG_BLOCKSIZE,spcstats%usampling)
    Call stats_init(spcstats%intranrg(molec)%con,"Constraint Energy",&
        INTRA_NRG_BLOCKSIZE,spcstats%usampling)
    Call stats_init(spcstats%intranrg(molec)%ext,"External Energy",&
        INTRA_NRG_BLOCKSIZE,spcstats%usampling)
    Call stats_init(spcstats%intranrg(molec)%total, &
        "Total Intramolecular Energy", INTRA_NRG_BLOCKSIZE,spcstats%usampling)

  End Subroutine storestats_initIntraEnergy

  !----------------------------------------------------------------------------
  ! Returns the values of the intramolecular energy
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function storestats_getintranrg(spcstats,spc,sType,eType)
    Type(Species_Stats), Intent(In)     :: spcstats
    Integer, Intent(In) :: spc
    Character(*), Intent(In) :: sType
    Character(*), Intent(In), Optional :: eType
    Type(Statistics) :: nrg

    If (Present(eType)) Then
      Select Case(toupper(eType))
      Case ('STRETCH')
        nrg = spcstats%intranrg(spc)%bs
      Case ('BENDING')
        nrg = spcstats%intranrg(spc)%bb
      Case ('TORSION')
        nrg = spcstats%intranrg(spc)%tor
      Case ('INTRAPAIR')
        nrg = spcstats%intranrg(spc)%ip
      Case ('INTRACOUL')
        nrg = spcstats%intranrg(spc)%ic
      Case ('CONSTRAINT')
        nrg = spcstats%intranrg(spc)%con
      Case ('EXTERNAL')
        nrg = spcstats%intranrg(spc)%ext
      Case ('TOTAL')
        nrg = spcstats%intranrg(spc)%total
      End Select
    Else
      nrg = spcstats%intranrg(spc)%total
    End If

    Select Case(ToUpper(sType))
    Case ('INST')
      storestats_getintranrg = stats_getvalue(nrg)
    Case ('CAVG')
      storestats_getintranrg = stats_getcavg(nrg)
    Case ('BLOCK')
      storestats_getintranrg = stats_getblock(nrg)
    Case ('STD')
      storestats_getintranrg = stats_getstd(nrg)
    End Select

  End Function storestats_getintranrg

  !----------------------------------------------------------------------------
  ! Returns All the values of the intramolecular energy
  ! Requires:  spcstats -- storage structure for species statistics
  !            spc -- species number
  !            sType -- statistics type
  !----------------------------------------------------------------------------
  Subroutine storestats_getallintra(spcstats,spc,sType,intra)
    Type(Species_Stats), Intent(In)  :: spcstats
    Integer, Intent(In)              :: spc
    Character(*), Intent(In)         :: sType
    Type(EnergyPlus), Intent(InOut)  :: intra          

    Type(Statistics) :: nrg

    Select Case(ToUpper(sType))
    Case ('INST')
      nrg = spcstats%intranrg(spc)%bs
      intra%intranrg(STRETCH_INDEX) = stats_getvalue(nrg)
      nrg = spcstats%intranrg(spc)%bb
      intra%intranrg(BENDING_INDEX) = stats_getvalue(nrg)
      nrg = spcstats%intranrg(spc)%tor
      intra%intranrg(TORSION_INDEX) = stats_getvalue(nrg)
      nrg = spcstats%intranrg(spc)%ip
      intra%intranrg(INTRAPAIR_INDEX) = stats_getvalue(nrg)
     nrg = spcstats%intranrg(spc)%ic
      intra%intranrg(INTRACOUL_INDEX) = stats_getvalue(nrg)
      nrg = spcstats%intranrg(spc)%con
      intra%intranrg(CONSTRAINT_INDEX) = stats_getvalue(nrg)
      nrg = spcstats%intranrg(spc)%ext
      intra%intranrg(EXTERNAL_INDEX) = stats_getvalue(nrg)
      nrg = spcstats%intranrg(spc)%total
      intra%intranrg(TOTAL_INDEX) = stats_getvalue(nrg)

    Case ('CAVG')
      nrg = spcstats%intranrg(spc)%bs
      intra%intranrg(STRETCH_INDEX) = stats_getcavg(nrg)
      nrg = spcstats%intranrg(spc)%bb
      intra%intranrg(BENDING_INDEX) = stats_getcavg(nrg)
      nrg = spcstats%intranrg(spc)%tor
      intra%intranrg(TORSION_INDEX) = stats_getcavg(nrg)
      nrg = spcstats%intranrg(spc)%ip
      intra%intranrg(INTRAPAIR_INDEX) = stats_getcavg(nrg)
      nrg = spcstats%intranrg(spc)%ic
      intra%intranrg(INTRACOUL_INDEX) = stats_getcavg(nrg)
      nrg = spcstats%intranrg(spc)%con
      intra%intranrg(CONSTRAINT_INDEX) = stats_getcavg(nrg)
      nrg = spcstats%intranrg(spc)%ext
      intra%intranrg(EXTERNAL_INDEX) = stats_getcavg(nrg)
      nrg = spcstats%intranrg(spc)%total
      intra%intranrg(TOTAL_INDEX) = stats_getcavg(nrg)

    Case ('BLOCK')
      nrg = spcstats%intranrg(spc)%bs
      intra%intranrg(STRETCH_INDEX) = stats_getblock(nrg)
      nrg = spcstats%intranrg(spc)%bb
      intra%intranrg(BENDING_INDEX) = stats_getblock(nrg)
      nrg = spcstats%intranrg(spc)%tor
      intra%intranrg(TORSION_INDEX) = stats_getblock(nrg)
      nrg = spcstats%intranrg(spc)%ip
      intra%intranrg(INTRAPAIR_INDEX) = stats_getblock(nrg)
      nrg = spcstats%intranrg(spc)%ic
      intra%intranrg(INTRACOUL_INDEX) = stats_getblock(nrg)
      nrg = spcstats%intranrg(spc)%con
      intra%intranrg(CONSTRAINT_INDEX) = stats_getblock(nrg)
      nrg = spcstats%intranrg(spc)%ext
      intra%intranrg(EXTERNAL_INDEX) = stats_getblock(nrg)
      nrg = spcstats%intranrg(spc)%total
      intra%intranrg(TOTAL_INDEX) = stats_getblock(nrg)

    Case ('STD')
      nrg = spcstats%intranrg(spc)%bs
      intra%intranrg(STRETCH_INDEX) = stats_getstd(nrg)
      nrg = spcstats%intranrg(spc)%bb
      intra%intranrg(BENDING_INDEX) = stats_getstd(nrg)
      nrg = spcstats%intranrg(spc)%tor
      intra%intranrg(TORSION_INDEX) = stats_getstd(nrg)
      nrg = spcstats%intranrg(spc)%ip
      intra%intranrg(INTRAPAIR_INDEX) = stats_getstd(nrg)
      nrg = spcstats%intranrg(spc)%ic
      intra%intranrg(INTRACOUL_INDEX) = stats_getstd(nrg)
      nrg = spcstats%intranrg(spc)%con
      intra%intranrg(CONSTRAINT_INDEX) = stats_getstd(nrg)
      nrg = spcstats%intranrg(spc)%ext
      intra%intranrg(EXTERNAL_INDEX) = stats_getstd(nrg)
      nrg = spcstats%intranrg(spc)%total
      intra%intranrg(TOTAL_INDEX) = stats_getstd(nrg)

    End Select

  End Subroutine storestats_getallintra

  !----------------------------------------------------------------------------
  ! Return the temperature of the requested molecule and type
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Real(Kind=Rdbl) Function storestats_gettemp(spcstats,spc,tType,sType)
    Type(Species_Stats), Intent(In)     :: spcstats
    Integer, Intent(In) :: spc
    Character(*), Intent(In) :: tType, sType
    Type(Statistics), Pointer :: temppntr

    Select Case(ToUpper(tType))

    Case ('TATOM')
      temppntr => spcstats%tatom(spc)
    Case ('TMOLE')
      temppntr => spcstats%tmole(spc)
    End Select

    Select Case(ToUpper(sType))
    Case ('INST')
      storestats_gettemp = stats_getvalue(temppntr)
    Case ('BLOCK')
      storestats_gettemp = stats_getblock(temppntr)
    Case ('CAVG')
      storestats_gettemp = stats_getcavg(temppntr)
    Case ('STD')
      storestats_gettemp = sqrt(stats_getstd(temppntr))
    End Select

  End Function storestats_gettemp

  !----------------------------------------------------------------------------
  ! Returns the noncoulombic interaction energy for the given pair
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function storestats_getnoncoul(spcstats,spc1,spc2,sType)
    Type(Species_Stats), Intent(In)     :: spcstats
    Integer, Intent(In) :: spc1, spc2
    Character(*), Intent(In) :: sType

    Select Case(ToUpper(sType))
    Case ('INST')
      storestats_getnoncoul = stats_getvalue(spcstats%noncoulnrg(spc1,spc2))
    Case ('CAVG')
      storestats_getnoncoul = stats_getcavg(spcstats%noncoulnrg(spc1,spc2))
    Case ('BLOCK')
      storestats_getnoncoul = stats_getblock(spcstats%noncoulnrg(spc1,spc2))
    Case ('STD')
      storestats_getnoncoul = stats_getstd(spcstats%noncoulnrg(spc1,spc2))
    End Select

  End Function storestats_getnoncoul

  !----------------------------------------------------------------------------
  ! Returns the coulombic interaction energy for the given pair
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function storestats_getcoul(spcstats,spc1,spc2,sType)
    Type(Species_Stats), Intent(In)     :: spcstats
    Integer, Intent(In)      :: spc1, spc2
    Character(*), Intent(In) :: sType

    Select Case(ToUpper(sType))
    Case ('INST')
      storestats_getcoul = stats_getvalue(spcstats%coulnrg(spc1,spc2))
    Case ('CAVG')
      storestats_getcoul = stats_getcavg(spcstats%coulnrg(spc1,spc2))
    Case ('BLOCK')
      storestats_getcoul = stats_getblock(spcstats%coulnrg(spc1,spc2))
    Case ('STD')
      storestats_getcoul = stats_getstd(spcstats%coulnrg(spc1,spc2))
    End Select

  End Function storestats_getcoul

  !----------------------------------------------------------------------------
  ! Returns the kinetic energy for the given molecule
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function storestats_getke(spcstats,spc,sType)
    Type(Species_Stats), Intent(In)     :: spcstats
    Integer, Intent(In)      :: spc
    Character(*), Intent(In) :: sType

    Select Case(ToUpper(sType))
    Case ('INST')
      storestats_getke = stats_getvalue(spcstats%ke(spc))
    Case ('CAVG')
      storestats_getke = stats_getcavg(spcstats%ke(spc))
    Case ('BLOCK')
      storestats_getke = stats_getblock(spcstats%ke(spc))
    Case ('STD')
      storestats_getke = stats_getstd(spcstats%ke(spc))
    End Select

  End Function storestats_getke

  !-------------------------------------------------------------------------
  ! Displays the different energies, 
  !  1) pairwise -coulombic, noncoulombic
  !  2) intramolecular energies
  ! Requires:  spcstats -- storage structure for species statistics
  !            indent -- number of spaces to indent
  !            optunit -- unit to write into, optional
  !-------------------------------------------------------------------------
  Subroutine storestats_displaynrg(spcstats,indent,optunit)
    Type(Species_Stats), Intent(In)   :: spcstats
    Integer, Intent(In)               :: indent
    Integer, Optional, Intent(In)     :: optunit

    Integer                  :: unitno, i, j
    Character(len=strLen)    :: molecname,string1,string2
    Character(len=indent)    :: blank
    Real(kind=RDbl)          :: coulnrg,noncoulnrg
    Real(kind=RDbl)          :: bsnrg,bbnrg,tornrg,ipnrg,icnrg,connrg,extnrg,totnrg
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    blank = Repeat(' ',indent)

    Do i = 1,spcstats%nmolec_types
      molecname = molecules_name(i)
      Write(unitno,'(2a,a)') blank,"Species Name : ",molecname
      string1 = int2str(i)
      Write(unitno,'(3a)') blank,"Species Type : ",Trim(string1)
      Write(unitno,'(2a)') blank,"Pairwise Energies(noncoul, coul) in kJ/mol"
      Do j = 1,spcstats%nmolec_types
        molecname = molecules_name(j)
        coulnrg = stats_getvalue(spcstats%coulnrg(i,j))
        noncoulnrg = stats_getvalue(spcstats%noncoulnrg(i,j))
        string1 = real2str(noncoulnrg,8)
        string2 = real2str(coulnrg,8)
        Write(unitno, '(a,3x,a15,3(2x,a))') blank,Trim(molecname), " : ", &
            Trim(string1), Trim(string2)
      End Do

      bsnrg = stats_getvalue(spcstats%intranrg(i)%bs)
      bbnrg = stats_getvalue(spcstats%intranrg(i)%bb)
      tornrg = stats_getvalue(spcstats%intranrg(i)%tor)
      ipnrg = stats_getvalue(spcstats%intranrg(i)%ip)
      icnrg = stats_getvalue(spcstats%intranrg(i)%ic)
      connrg = stats_getvalue(spcstats%intranrg(i)%con)
      extnrg = stats_getvalue(spcstats%intranrg(i)%ext)
      totnrg = stats_getvalue(spcstats%intranrg(i)%total)

      Write(unitno,'(2a)') blank,"Intramolecular Energies &
          &(bs,bb,tor,ip,con,total) in kJ/mol"
      Write(unitno,'(a,1x)',Advance='No') blank
      string1 = real2str(bsnrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(bbnrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(tornrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(ipnrg,8) 
     Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(icnrg,8)
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(connrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(totnrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      Write(unitno,*) 
      Write(unitno,*) 
    End Do

  End Subroutine storestats_displaynrg



  !-------------------------------------------------------------------------
  ! Displays the different energy averages; similar to _displaynrg 
  !  1) pairwise -coulombic, noncoulombic
  !  2) intramolecular energies
  ! Requires:  spcstats -- storage structure for species statistics
  !            atype -- type of average "BLOCK" or "CUM" or "STD" or "INST"
  !            indent -- number of spaces to indent
  !            unitno -- unit to write into
  !-------------------------------------------------------------------------
  Subroutine storestats_displayavgs(spcstats,avg_type,indent,unitno)
    Type(Species_Stats), Intent(In)     :: spcstats
    Character(len=strLen), Intent(In)   :: avg_type
    Integer, Intent(In)                 :: indent
    Integer, Intent(In)                 :: unitno

    Integer                  :: i, j
    Character(len=strLen)    :: molecname,string1,string2, atype
    Character(len=indent)    :: blank
    Real(kind=RDbl)          :: coulnrg,noncoulnrg
    Real(kind=RDbl)          :: bsnrg,bbnrg,tornrg,ipnrg,icnrg
    Real(kind=RDbl)          :: connrg,extnrg,totnrg

    blank = Repeat(' ',indent)

    atype=toupper(avg_type) ! atype = CAVG, BLOCK, INST, STD

    Do i = 1,spcstats%nmolec_types
      molecname = molecules_name(i)
      string1 = int2str(i)
      Write(unitno,'(6a)') blank,"Species Name (type): ",Trim(molecname),"  (", &
          Trim(string1),")"
      Write(unitno,'(4a)') blank," Pairwise ",Trim(aType), &
          " Energies(noncoul, coul) in kJ/mol"
      Do j = 1,spcstats%nmolec_types
        coulnrg = stats_getfield(spcstats%coulnrg(i,j), atype)
        noncoulnrg = stats_getfield(spcstats%noncoulnrg(i,j), atype)
        molecname = molecules_name(j)
        string1 = real2str(noncoulnrg,8)
        string2 = real2str(coulnrg,8)
        Write(unitno, '(a,3x,a15,3(2x,a))') blank,Trim(molecname), " : ", &
            Trim(string1), Trim(string2)
      End Do
      bsnrg  = stats_getfield(spcstats%intranrg(i)%bs,aType)
      bbnrg  = stats_getfield(spcstats%intranrg(i)%bb,aType)
      tornrg = stats_getfield(spcstats%intranrg(i)%tor,aType)
      ipnrg  = stats_getfield(spcstats%intranrg(i)%ip,aType)
      icnrg  = stats_getfield(spcstats%intranrg(i)%ic,aType)
      connrg = stats_getfield(spcstats%intranrg(i)%con,aType)
      extnrg = stats_getfield(spcstats%intranrg(i)%ext,aType)
      totnrg = stats_getfield(spcstats%intranrg(i)%total,aType)

      Write(unitno,'(2a)') blank," Intramolecular Energies &
          &(bs,bb,tor,ip,con,total) in kJ/mol"
      Write(unitno,'(a,1x)',Advance='No') blank
      string1 = real2str(bsnrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(bbnrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(tornrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(ipnrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(icnrg,8)
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(connrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      string1 = real2str(totnrg,8) 
      Write(unitno,'(2x,a)',Advance='No') Trim(string1)
      Write(unitno,*) 
    End Do

  End Subroutine storestats_displayavgs

  !----------------------------------------------------------------------------
  ! Cleanup the molecules allocations
  ! Requires:  spcstats -- storage structure for species statistics
  !----------------------------------------------------------------------------
  Subroutine storestats_clean(spcstats)
    Type(Species_Stats), Intent(InOut)     :: spcstats

    Integer :: error

    Deallocate (spcstats%noncoulnrg, STAT=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__,'noncoulnrg')

    Deallocate (spcstats%coulnrg, STAT=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__,'coulnrg')

    Deallocate(spcstats%intranrg,Stat=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__,'intranrg')

    Deallocate(spcstats%ke,Stat=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__,'ke')

    Deallocate(spcstats%tatom,Stat=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__,'tatom')

    Deallocate(spcstats%tmole,Stat=error)
    If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__,'tmole')

  End Subroutine storestats_clean

End Module storestats


