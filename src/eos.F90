!-----------------------------------------------------------------
! This module is the gateway to the different equations of state.
!-----------------------------------------------------------------

Module eos

  Use defaults, Only: RDbl, strLen
  Use utils, Only: isfileopen, filesrchstr,toupper
  Use file, Only: file_getunit
  Use ideal, Only: Ideal_EOS,getcompressibility,getfugacity,&
      setpressure,settemperature,ideal_getConfInteg,ideal_init
  Use virial, Only: Virial_EOS,getcompressibility,getfugacity,&
      setpressure,settemperature,virial_init

  Implicit None
  Save

  Private
  Public :: EOS_Models,eos_init,eos_setpressure,eos_settemperature,&
      eos_getfugacity,eos_getConfInteg

  Type EOS_Models
    Type(Ideal_EOS), Pointer   :: ideal
    Type(Virial_EOS), Pointer  :: virial
    Character(len=strLen)      :: eosmodel
  End Type EOS_Models

Contains
  
  !----------------------------------------------------------------------
  ! The optional argument "optlineno" is passed for the case when
  ! the "ctrl_tag" is listed as a parameter in another section of the
  ! control file.  For e.g. the GCMC section of the control file has
  ! the tag of the section where the eqn. of state parameters can be
  ! found.  Now, when searching for the "ctrl_tag" we don't want to pick
  ! the tag string in the GCMC section.  So, we pass the line number of
  ! GCMC section and make sure that the line number found for the searched
  ! tag is less than "optlineno".  Which means that the eqn. of state
  ! parameter section used by GCMC has to be defined BEFORE the 
  ! corresponding GCMC section.
  !-----------------------------------------------------------------------
  Subroutine eos_init(eosparams, ctrl_filename, ctrl_tag, optlineno)
    Type(EOS_Models), Intent(inout) :: eosparams
    Character(*), Intent(in)  :: ctrl_filename
    Character(*), Intent(in)  :: ctrl_tag
    Integer, Optional, Intent(in)  :: optlineno
    
    Character(len=strLen)     :: line
    Integer     :: error, unitno, lineno

    !** Check if the file is open
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
      Open(unit=unitno, file=ctrl_filename, status='old')
    End If

    !** Find the ctrl_tag
    lineno = filesrchstr(unitno, ctrl_tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", Trim(ctrl_tag), " in the control file"
      Stop
    End If
    ! If the optional lineno is present make sure it is less than
    ! "lineno".
    If (Present(optlineno)) Then
      If (optlineno <= lineno) Then
        Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
            " Could not find the section '", Trim(ctrl_tag), "'"
        Write(0,'(1x,2a)') "Either it is not defined in the control file or", &
            " it is defined AFTER the section referring to it"
        Stop
      End If
    End If

    !** Nullify the pointers to the models
    Nullify(eosparams%virial)
    Nullify(eosparams%ideal)

    !** Read the rest of the stuff
    Read(unitno,*) eosparams%eosmodel

    !** Initialize the pointer depending on the equation of state
    Select Case (Trim(toupper(eosparams%eosmodel)))
    Case ("VIRIAL")
      Allocate(eosparams%virial, STAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Could not allocate 'virial'"
        Stop
      End If
      Call virial_init(eosparams%virial, unitno)

    Case ("IDEAL")
      Allocate(eosparams%ideal, STAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Could not allocate 'ideal'"
        Stop
      End If
      Call ideal_init(eosparams%ideal, unitno)

    Case Default
      Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
          " Equation of state ", Trim(eosparams%eosmodel), " not supported"
      Stop
    End Select
    
  End Subroutine eos_init

  !--------------------------------------------------------
  ! Set the pressure of the sorbate type "sorbtype"
  !--------------------------------------------------------
  Subroutine eos_setpressure(eosparams, sorbtype, pp)
    Type(EOS_Models), Intent(inout)  :: eosparams
    Real(kind=RDbl), Intent(in)      :: pp
    Integer, Intent(in)              :: sorbtype
    
    If (Associated (eosparams%virial)) Then
      Call setpressure(eosparams%virial, sorbtype, pp)
    Else If (Associated (eosparams%ideal)) Then
      Call setpressure(eosparams%ideal, sorbtype, pp)
    End If
  End Subroutine eos_setpressure

  !--------------------------------------------------------
  ! Set the temperature of the system
  !--------------------------------------------------------
  Subroutine eos_settemperature(eosparams, tk)
    Type(EOS_Models), Intent(inout)  :: eosparams
    Real(kind=RDbl), Intent(in)      :: tk
    
    If (Associated (eosparams%virial)) Then
      Call settemperature(eosparams%virial, tk)
    Else If (Associated (eosparams%ideal)) Then
      Call settemperature(eosparams%ideal, tk)
    End If
  End Subroutine eos_settemperature

  !--------------------------------------------------------
  ! Get the fugacity of the sorbate type "sorbtype"
  !--------------------------------------------------------
  Real(kind=RDbl) Function eos_getfugacity(eosparams, sorbtype)
    Type(EOS_Models), Intent(inout)  :: eosparams
    Integer, Intent(in)      :: sorbtype
    
    If (Associated (eosparams%virial)) Then
      eos_getfugacity = getfugacity(eosparams%virial, sorbtype)
    Else If (Associated (eosparams%ideal)) Then
      eos_getfugacity = getfugacity(eosparams%ideal, sorbtype)
    End If
  End Function eos_getfugacity

  !---------------------------------------------------------
  ! Get the compressibility of the system
  !---------------------------------------------------------
  Real(kind=RDbl) Function eos_getcompressibility(eosparams)
    Type(EOS_Models), Intent(inout)  :: eosparams
    If (Associated (eosparams%virial)) Then
      eos_getcompressibility = getcompressibility(eosparams%virial)
    Else If (Associated (eosparams%ideal)) Then
      eos_getcompressibility = getcompressibility(eosparams%ideal)
    End If
  End Function eos_getcompressibility


  !---------------------------------------------------------
  ! Gets the Configurational Integral of Ideal gas chains 
  !---------------------------------------------------------
  Real(kind=RDbl) Function eos_getConfInteg(eosparams,sorbtype)
    Type(EOS_Models), Intent(inout)  :: eosparams
    Integer,Intent(in) :: sorbtype
    If (Associated (eosparams%ideal)) Then
      eos_getConfInteg = ideal_getConfInteg(eosparams%ideal,sorbtype)
    Else 
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "EOS Ideal should be used while trying to access &
          & confinteg "
      Stop
    End If
  End Function eos_getConfInteg
  
End Module eos

