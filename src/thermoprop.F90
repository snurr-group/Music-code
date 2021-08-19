!--------------------------------------------------
! This module contains the thermodynamic parameters 
! and values for a single component
!--------------------------------------------------
Module thermoprop

  Use utils, Only:
  Use defaults, Only: RDbl, strLen

  Implicit None
  Save

  Private
  Public :: SysThermo_Params, CompThermo_Params, thermoprop_getsorbindx, &
      thermoprop_sysinit, thermoprop_settemp, thermoprop_setpressure, &
      thermoprop_init

  !** Thermodynamic parameters for a single component
  Type CompThermo_Params
    Integer               :: sorbtype
    Real(kind=RDbl)       :: tcrit, pcrit, vcrit, zcrit, acentric, tboil
    Real(kind=RDbl)       :: tref, pref

    Real(kind=RDbl)       :: tk, p, z, molv, rho, entr, enth
    Real(kind=RDbl)       :: tred, pred, vred
    Real(kind=RDbl)       :: Zig  !** conf integral of Ideal gas chain molecs
  End Type CompThermo_Params


  !** Thermodynamic parameters for the entire system
  Type SysThermo_Params
    ! Parameters of each component
    Type(CompThermo_Params), Dimension(:), Pointer  :: cmp

    ! Parameters for the overall system 
    ! Some fields don't necessarily apply to multicomponent systems
    Type(CompThermo_Params)                  :: sys
    Real(kind=RDbl), Dimension(:), Pointer   :: mfrac ! Composition
    Integer           :: ncomps          ! No. of components
  End Type SysThermo_Params

  Interface thermoprop_init
    Module Procedure thermoprop_sysinit
  End Interface

  Interface display
    Module Procedure thermoprop_sysdisplay
    Module Procedure thermoprop_compdisplay
  End Interface

  Interface thermoprop_setpressure
    Module Procedure thermoprop_setcomppressure
  End Interface
  
  Interface thermoprop_settemp
    Module Procedure thermoprop_setsystemp
  End Interface

Contains
  
  !----------------------------------------------------------------------------
  ! Sets all variables to zero
  !----------------------------------------------------------------------------
  Subroutine thermoprop_compinit(thparams)
    Type(CompThermo_Params), Intent(inout)  ::thparams
    
    thparams%tcrit = 0.0_RDbl
    thparams%pcrit = 0.0_RDbl
    thparams%vcrit = 0.0_RDbl
    thparams%zcrit = 0.0_RDbl
    thparams%acentric = 0.0_RDbl
    thparams%tboil = 0.0_RDbl
    thparams%tref  = 0.0_RDbl
    thparams%pref  = 0.0_RDbl

    thparams%tk    = 0.0_RDbl
    thparams%p     = 0.0_RDbl
    thparams%z     = 0.0_RDbl
    thparams%molv  = 0.0_RDbl
    thparams%rho   = 0.0_RDbl
    thparams%enth  = 0.0_RDbl
    thparams%entr  = 0.0_RDbl
    thparams%tred  = 0.0_RDbl
    thparams%pred  = 0.0_RDbl
    thparams%vred  = 0.0_RDbl
  End Subroutine thermoprop_compinit
 
  !----------------------------------------------------------
  ! Set the system size and initialize all the fields to zero
  !----------------------------------------------------------
  Subroutine thermoprop_sysinit(thparams, ncomps)
    Type(SysThermo_Params), Intent(inout) :: thparams
    Integer, Intent(in)   :: ncomps
    Integer    :: i, error
    
    !** Allocate the necessary memory    
    thparams%ncomps = ncomps
    Allocate(thparams%cmp(ncomps), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Could not allocate 'cmp'"
      Stop
    End If
    Allocate(thparams%mfrac(ncomps), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Could not allocate 'mfrac'"
      Stop
    End If

    !** Initialize the system fields to zero
    Call thermoprop_compinit(thparams%sys)
    thparams%mfrac = 0.0_RDbl

    !** Initialize the individual fields
    Do i=1, ncomps
      Call thermoprop_compinit(thparams%cmp(i))
    End Do
  End Subroutine thermoprop_sysinit

  !-------------------------------------------------------
  ! Get the index of the thparams that hold the parameters
  ! for the sorbate type
  !-------------------------------------------------------
  Integer Function thermoprop_getsorbindx(thparams, sorbtype)
    Type(SysThermo_Params), Intent(in)   :: thparams
    Integer, Intent(in)  :: sorbtype
    Integer    :: i
    
    thermoprop_getsorbindx = 0
    Do i=1, thparams%ncomps
      If (thparams%cmp(i)%sorbtype == sorbtype) Then
        thermoprop_getsorbindx = i
      End If
    End Do
  End Function thermoprop_getsorbindx

  !----------------------------------------------
  ! Sets the temperature of the system and every
  ! component in the system
  !---------------------------------------------
  Subroutine thermoprop_setsystemp(thparams, tk)
    Type(SysThermo_Params), Intent(inout) :: thparams
    Real(kind=RDbl), Intent(in)   :: tk
    Real(kind=RDbl)    :: tcrit
    Integer   :: i

    !** Set the temperature of all the components as well
    !** as the reduced temperature of each one
    thparams%sys%tk = tk
    Do i=1, thparams%ncomps
      thparams%cmp(i)%tk = tk
      tcrit = thparams%cmp(i)%tcrit
      ! It is possible that tcrit may not be set e.g. for ideal gas
      If (tcrit /= 0.0_RDbl) Then
        thparams%cmp(i)%tred = tk/tcrit
      End If
    End Do
  End Subroutine thermoprop_setsystemp

  !------------------------------------------------------
  ! Set the temperature of each component
  !------------------------------------------------------
  Subroutine thermoprop_setcomppressure(thparams, sorbtype, pp)
    Type(SysThermo_Params), Intent(inout) :: thparams
    Integer, Intent(in)          :: sorbtype 
    Real(kind=RDbl), Intent(in)  :: pp
    Real(kind=RDbl)              :: oldpp, pcrit, p
    Integer   :: sorbindx, i

    !** Find the sorbate type
    sorbindx = thermoprop_getsorbindx(thparams, sorbtype)
    If (sorbindx == 0) Then
      Write(0,'(1x,2a,i4, a, i4)') __FILE__," : ",__LINE__, &
          " Could not find parameters for sorbate type ", sorbtype
      Stop
    End If

    !** Set the partial pressure
    oldpp = thparams%cmp(sorbindx)%p
    thparams%cmp(sorbindx)%p = pp
    
    !** Update the total pressure
    thparams%sys%p = thparams%sys%p + pp - oldpp

    !** Update all the mole fraction
    Do i=1, thparams%ncomps
      p = thparams%cmp(i)%p
      thparams%mfrac(i) = pp/thparams%sys%p
    End Do

    !** Set the reduced pressure
    pcrit = thparams%cmp(sorbindx)%p
    ! It is possible that pcrit may not be set e.g. for ideal gas
    If (pcrit /= 0.0_RDbl) Then
      thparams%cmp(sorbindx)%pred = thparams%cmp(sorbindx)%p/pcrit
    End If

  End Subroutine thermoprop_setcomppressure

  !-----------------------------------------------
  ! Display the values of the parameters and the 
  ! variables for the component "comp"
  !-----------------------------------------------
  Subroutine thermoprop_compdisplay(thparams, comp)
    Type(CompThermo_Params), Intent(in) :: thparams
    Integer, Intent(in)    :: comp
  End Subroutine thermoprop_compdisplay

  !-----------------------------------------------
  ! Display the values of the parameters and the 
  ! variables for all the components
  !-----------------------------------------------
  Subroutine thermoprop_sysdisplay(thparams)
    Type(SysThermo_Params), Intent(in) :: thparams
  End Subroutine thermoprop_sysdisplay

  !-----------------------------------------------------------
  ! Deallocate the memory and do whatever cleanup is required
  !-----------------------------------------------------------
  Subroutine thermoprop_cleanup(thparams)
    Type(SysThermo_Params), Intent(inout) :: thparams
    Integer    :: error
    
    Deallocate(thparams%cmp, thparams%mfrac, STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not deallocate 'cmp' or 'mfrac'"
      Stop
    End If
    
  End Subroutine thermoprop_cleanup
End Module thermoprop
