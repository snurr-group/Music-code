!---------------------------------------------------------
! The IDEAL GAS equation of state.  This module is 
! only for the sake of completeness.  Most of the stuff
! is trivial
!---------------------------------------------------------
Module ideal

  Use defaults, Only: RDbl, strLen,one
  Use thermoprop, Only: SysThermo_Params, thermoprop_getsorbindx, &
      thermoprop_sysinit, thermoprop_settemp, thermoprop_setpressure
  Use molecules, Only: molecules_gettype
  Use utils, Only : toreal,split,stripcmnt

  Implicit None
  Save

  Private
  Public :: Ideal_EOS, setpressure, ideal_init, settemperature, &
      getcompressibility, getfugacity, ideal_getConfInteg

  Type Ideal_EOS
    Type(SysThermo_Params)   :: thparams
    Integer                  :: ncomps
  End Type Ideal_EOS

  Interface setpressure
    Module Procedure ideal_setcomppressure
  End Interface

  Interface settemperature
    Module Procedure ideal_setsystemp
  End Interface
  
  Interface getcompressibility
    Module Procedure ideal_getsysZ
  End Interface

  Interface getfugacity
    Module Procedure ideal_getfugacity
  End Interface

Contains

  !----------------------------------------------
  ! Initialize the fields of "idparams"
  !----------------------------------------------
  Subroutine ideal_init(idparams, unitno)
    Type(Ideal_EOS), Intent(inout) :: idparams
    Integer, Intent(in) :: unitno
    Character(len=strLen)   :: molec_name
    Character(len=5*strLen)   :: text
    Character(len=strLen),Dimension(10) :: fields

    Integer             :: i, sorbtype, nfields
    Real(kind=RDbl)     :: Zig
    !** Get the total no. of components and initialize the
    !** fields of each component and the overall system to zero
    Read(unitno, *) idparams%ncomps
    Call thermoprop_sysinit(idparams%thparams, idparams%ncomps)
    
    Do i = 1,idparams%ncomps
      Read(unitno, '(a)') text 
      text = Trim(stripcmnt(text))
      nfields = split(text,fields,",")

      If (nfields>1) Then
        !** That means we have molecule name and Zig listed on one line
        sorbtype = molecules_gettype(fields(1))
        Zig = toreal(fields(2))
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "config integral for molecule "//Trim(fields(1))//" is",Zig
      Else
        sorbtype = molecules_gettype(fields(1))
        Zig = one
      Endif

      If (sorbtype == 0) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            "Molecule ", Trim(molec_name), " does not exist in molecs list"
        Stop
      End If

      !** Only things of interest to us for ideal gas params
      idparams%thparams%cmp(i)%sorbtype = sorbtype
      idparams%thparams%cmp(i)%Zig = Zig 

    End Do

  End Subroutine ideal_init

  !------------------------------------------------------
  ! Set the temperature of the mixture to "tk"
  !------------------------------------------------------
  Subroutine ideal_setsystemp(idparams, tk)
    Type(Ideal_EOS), Intent(inout) :: idparams
    Real(kind=RDbl), Intent(in)   :: tk
    Call thermoprop_settemp(idparams%thparams, tk)
  End Subroutine ideal_setsystemp

  !------------------------------------------------------
  ! Set the partial pressure of the component sorbtype
  !------------------------------------------------------
  Subroutine ideal_setcomppressure(idparams, sorbtype, pp)
    Type(Ideal_EOS), Intent(inout) :: idparams
    Integer, Intent(in)          :: sorbtype
    Real(kind=RDbl), Intent(in)  :: pp
    
    Call thermoprop_setpressure(idparams%thparams, sorbtype, pp)
  End Subroutine ideal_setcomppressure
  
  !------------------------------------------------------------
  ! Returns the fugacity of the sorbate type "sorbtype"
  !------------------------------------------------------------
  Real(kind=RDbl) Function ideal_getfugacity(idparams, sorbtype)
    Type(Ideal_EOS), Intent(in)  :: idparams
    Integer, Intent(in)          :: sorbtype 
    
    Integer             :: indx
    Real(kind=RDbl)     :: pp
    
    !** Get the sorbate index
    indx = thermoprop_getsorbindx(idparams%thparams, sorbtype)

    !** Get the partial pressure
    pp = idparams%thparams%cmp(indx)%p

    !** Now calculate the fugacity
    ideal_getfugacity = pp  
    Return
  End Function ideal_getfugacity


  !-----------------------------------------------------------
  ! Returns the compressibility of the sorbate type "sorbtype"
  ! Use the formulae given in Sandler, 3rd Ed. Pg. 286, 408
  ! to get the compressibility.  Note that this is different
  ! from the one used by Smith and Van Ness.
  !-----------------------------------------------------------  
  Real(kind=RDbl) Function ideal_getsysZ(idparams)
    Type(Ideal_EOS), Intent(in) :: idparams
    ideal_getsysZ = 1.0_RDbl
    Return
  End Function ideal_getsysZ


  !-----------------------------------------------------------
  ! Returns the configurational integral (Z/Omega) for Ideal chains 
  !-----------------------------------------------------------  
  Real(kind=RDbl) Function ideal_getConfInteg(idparams,sorbtype)
    Type(Ideal_EOS), Intent(in) :: idparams
    Integer, Intent(in) :: sorbtype
    Integer :: indx

    !** index of molecule sorbtype in thermoprop array
    indx = thermoprop_getsorbindx(idparams%thparams, sorbtype)
    ideal_getConfInteg = idparams%thparams%cmp(indx)%Zig
    Return
  End Function ideal_getConfInteg
  
End Module ideal
