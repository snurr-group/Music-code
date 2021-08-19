!----------------------------------------
! The Virial equation of state
!----------------------------------------
Module virial

  Use defaults, Only: RDbl, strLen, Rgas
  Use thermoprop, Only: SysThermo_Params, CompThermo_Params, &
      thermoprop_getsorbindx, thermoprop_init, thermoprop_settemp, &
      thermoprop_setpressure
  Use molecules, Only: molecules_gettype

  Implicit None
  Save

  Private
  Public :: Virial_EOS, getfugacity, setpressure, settemperature, &
      getcompressibility, virial_init

  Type Virial_EOS
    Type(SysThermo_Params)      :: thparams
    Real(kind=RDbl), Dimension(:,:), Pointer :: Bij   ! 2nd Virial Coefficients
    ! Binary params required for the mixture calculations
    Type(CompThermo_Params), Dimension(:,:), Pointer :: binparams
    Integer   :: ncomps
  End Type Virial_EOS

  Interface setpressure
    Module Procedure virial_setcomppressure
  End Interface

  Interface settemperature
    Module Procedure virial_setsystemp
  End Interface
  
  Interface getcompressibility
    Module Procedure virial_getsysZ
  End Interface

  Interface getfugacity
    Module Procedure virial_getfugacity
  End Interface

Contains
  
  !---------------------------------------------------------
  ! Initialize the parameters needed for the virial equation
  ! of state.
  !---------------------------------------------------------
  Subroutine virial_init(vrparams, unitno)
    Type(Virial_EOS), Intent(inout)  :: vrparams
    Integer, Intent(in) :: unitno
    Integer             :: i, j, error, sorbtype, ncomps
    Character(len=strLen)   :: molec_name
    Real(kind=RDbl) :: tc1, tc2, pc1, pc2, vc1, vc2, zc1, zc2
    Real(kind=RDbl) :: acentric1, acentric2, acentric
    Real(kind=Rdbl) :: tc, pc, vc, zc

    !** Get the total no. of components and initialize the
    !** fields of each component and the overall system to zero
    Read(unitno, *) vrparams%ncomps
    Call thermoprop_init(vrparams%thparams, vrparams%ncomps)

    !** Initialize the Virial coefficient array
    ncomps = vrparams%ncomps
    Allocate(vrparams%Bij(ncomps, ncomps), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate 'Bij'"
      Stop
    End If

    !** Initialize the array for holding the binary critical parameters
    Allocate(vrparams%binparams(ncomps, ncomps), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate 'Bij'"
      Stop
    End If

    !** Read the necessary parameters for each molecule type
    Do i=1, vrparams%ncomps
      Read(unitno, *) molec_name
      sorbtype = molecules_gettype(molec_name)
      If (sorbtype == 0) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            " Molecule ", molec_name, " does not exist in molecules list"
        Stop
      End If
      vrparams%thparams%cmp(i)%sorbtype = sorbtype
      Read(unitno, *) tc, pc, vc
      Read(unitno, *) zc, acentric
      vrparams%thparams%cmp(i)%tcrit = tc
      vrparams%thparams%cmp(i)%pcrit = pc
      vrparams%thparams%cmp(i)%vcrit = vc  ! Only required if we have mixtures
      vrparams%thparams%cmp(i)%zcrit = zc  ! Only required if we have mixtures
      vrparams%thparams%cmp(i)%acentric = acentric
       
      ! Read the blank line
      If (i /= vrparams%ncomps) Then
        Read(unitno, *)
      End If
    End Do

    !** Generate the BINARY critical properties.
    !** This is based on correlations developed by Prausnitz
    !** See Smith and Van Ness, 3rd Ed., pg. 272
    Do i=1, vrparams%ncomps
      acentric1 = vrparams%thparams%cmp(i)%acentric
      vc1 = vrparams%thparams%cmp(i)%vcrit
      pc1 = vrparams%thparams%cmp(i)%pcrit
      tc1 = vrparams%thparams%cmp(i)%tcrit
      zc1 = vrparams%thparams%cmp(i)%zcrit
      Do j=1, vrparams%ncomps
        acentric2 = vrparams%thparams%cmp(j)%acentric
        vc2 = vrparams%thparams%cmp(j)%vcrit
        pc2 = vrparams%thparams%cmp(j)%pcrit
        tc2 = vrparams%thparams%cmp(j)%tcrit
        zc2 = vrparams%thparams%cmp(j)%zcrit

        !** Calculate the binary parameters
        acentric = (acentric1 + acentric2)/2.0_RDbl
        tc = Sqrt(tc1 * tc2)
        zc = (zc1 + zc2)/2.0_RDbl
        vc = ((vc1**(1.0/3.0) + vc2**(1.0/3.0))/2.0_RDbl)**3
        ! Calculate the value of critical pressure only for
        ! cross terms.
        If (i /= j) Then
          pc = zc*Rgas*tc*1.0e-3/vc   ! convert to kPa
        Else
          pc = pc1
        End If
        
        vrparams%binparams(i,j)%acentric = acentric
        vrparams%binparams(i,j)%tcrit = tc
        vrparams%binparams(i,j)%zcrit = zc
        vrparams%binparams(i,j)%vcrit = vc
        vrparams%binparams(i,j)%pcrit = pc
      End Do
    End Do
  End Subroutine virial_init

  !------------------------------------------------------
  ! Set the temperature of the mixture to "tk"
  !------------------------------------------------------
  Subroutine virial_setsystemp(vrparams, tk)
    Type(Virial_EOS), Intent(inout) :: vrparams
    Real(kind=RDbl), Intent(in)   :: tk
    Integer    :: i, j

    Call thermoprop_settemp(vrparams%thparams, tk)

    !** Generate the Virial coefficients now that we have
    !** the temperature
    Do i=1, vrparams%ncomps
      Do j=1, vrparams%ncomps
        Call virial_gencoeff(vrparams, i, j)
      End Do
    End Do
  End Subroutine virial_setsystemp

  !------------------------------------------------------
  ! Set the partial pressure of the component sorbtype
  !------------------------------------------------------
  Subroutine virial_setcomppressure(vrparams, sorbtype, pp)
    Type(Virial_EOS), Intent(inout) :: vrparams
    Integer, Intent(in)          :: sorbtype
    Real(kind=RDbl), Intent(in)  :: pp
    
    Call thermoprop_setpressure(vrparams%thparams, sorbtype, pp)
  End Subroutine virial_setcomppressure
  
  !-----------------------------------------------------------
  ! Generates the second virial coefficient for
  ! the "sorb1"-"sorb2" combination using the Pitzer
  ! and Abbott correlations.  See Perry's 6th Ed. pg. 3-268, 
  ! Smith and Van Ness 3rd Ed. pg. 86, Sandler 3rd Ed. pg. 286
  !-----------------------------------------------------------
  Subroutine virial_gencoeff(vrparams, indx1, indx2)
    Type(Virial_EOS), Intent(inout) :: vrparams
    Integer, Intent(in)          :: indx1, indx2 
    Real(kind=RDbl)  :: tk, B0, B1, tred, tcrit, pcrit, acentric, Bij
!!$    Real(kind=RDbl)  :: acentric_1, acentric_2, acentric_12
    
    !** Get the reduced temperature
    tk       = vrparams%thparams%sys%tk
    pcrit    = vrparams%binparams(indx1, indx2)%pcrit
    tcrit    = vrparams%binparams(indx1, indx2)%tcrit
    acentric = vrparams%binparams(indx1, indx2)%acentric
    tred     = tk/tcrit
    
    !** Use the Abbott correlation
    B0 = 0.083 - 0.422/(tred**1.6)
    B1 = 0.139 - 0.172/(tred**4.2)
    
    !** Use the Pitzer correlation
    Bij = (B0 + acentric*B1)*Rgas*tcrit/(pcrit*1.0e3)
    vrparams%Bij(indx1, indx2) = Bij
  End Subroutine virial_gencoeff

  !----------------------------------------------------
  ! Returns the overall virial coefficient for a system
  !----------------------------------------------------
  Real(kind=RDbl) Function virial_getbmix(vrparams)
    Type(Virial_EOS), Intent(in) :: vrparams

    Integer            :: i, j, ncomps
    Real(kind=RDbl)    :: bmix, bij, mfrac1, mfrac2

    bmix = 0.0_RDbl
    ncomps = vrparams%ncomps
    Do i=1, ncomps
      mfrac1 = vrparams%thparams%mfrac(i)
      Do j=1, ncomps
        mfrac2 = vrparams%thparams%mfrac(j)
        bij = vrparams%bij(i,j)
        bmix = bmix + mfrac1*mfrac2*bij
      End Do
    End Do
    virial_getbmix = bmix
    Return
  End Function virial_getbmix

  !------------------------------------------------------------
  ! Returns the fugacity of the sorbate type "sorbtype"
  !------------------------------------------------------------
  Real(kind=RDbl) Function virial_getfugacity(vrparams, sorbtype)
    Type(Virial_EOS), Intent(in) :: vrparams
    Integer, Intent(in)          :: sorbtype 
    
    Integer             :: indx, ncomps, j
    Real(kind=RDbl)     :: zmix, pp, psys, bij, bijsum, tk, fug, mfrac
    
    !** Get the sorbate index
    indx = thermoprop_getsorbindx(vrparams%thparams, sorbtype)

    !** Get the partial pressure
    pp = vrparams%thparams%cmp(indx)%p

    !** Get the compressibility and pressure of the system
    psys = vrparams%thparams%sys%p*1.0e3
    tk   = vrparams%thparams%sys%tk
    zmix = virial_getsysZ(vrparams)
    
    !** Get the summation of bij's
    ncomps = vrparams%ncomps
    bijsum = 0.0_RDbl
    Do j=1, ncomps
      mfrac = vrparams%thparams%mfrac(j)
      bij = vrparams%bij(indx, j)
      bijsum = bijsum + mfrac*bij
    End Do

    !** Now calculate the fugacity
    fug = 2.0_RDbl*psys*bijsum/(zmix*Rgas*tk) - Log(zmix)
    virial_getfugacity = Exp(fug)*pp  
    Return
  End Function virial_getfugacity

  !-----------------------------------------------------------
  ! Returns the compressibility of the sorbate type "sorbtype"
  ! Use the formulae given in Sandler, 3rd Ed. Pg. 286, 408
  ! to get the compressibility.  Note that this is different
  ! from the one used by Smith and Van Ness.
  !-----------------------------------------------------------  
  Real(kind=RDbl) Function virial_getsysZ(vrparams)
    Type(Virial_EOS), Intent(in) :: vrparams

    Real(kind=RDbl)    :: tk, pr, bmix, discrim
    
    !** Compressibility
    pr    = vrparams%thparams%sys%p
    tk    = vrparams%thparams%sys%tk
    bmix  = virial_getbmix(vrparams)

    If (pr == 0.0_RDbl .OR. tk == 0.0_RDbl) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(1x,a)') "Please set the system temperature and pressure"
      Stop
    End If
    
    !** We need to solve a quadratic eqn. to get z.
    !** Check to see if the discriminant >= 0
    discrim = 1.0_RDBl + (4.0_RDbl*bmix*pr*1.0e3/(Rgas*tk))

    If (discrim < 0.0_RDbl) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          "  The system is not physical"
      Stop
    End If
    virial_getsysZ = 0.5_RDbl*(1.0_RDbl + Sqrt(discrim))
    
  End Function virial_getsysZ

End Module virial

 


