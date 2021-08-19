Module mcparams

#ifdef USEONLY
  Use defaults, Only: RDbl, strLen
  Use utils, Only: toupper
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use matrix, Only: MatrixType
  Use rigidmcparams, Only: RigidMolecGCMCParams, RigidMolecNVTParams, &
      rigidmcparams_gettransmat, &
      rigidmcparams_init, rigidmcparams_adjustdeltatrans, &
      rigidmcparams_adjustdeltarot, rigidmcparams_setconstrnormal, &
      rigidmcparams_displaynvtstats, rigidmcparams_display
  Use rigidgc, Only: rigidgc_displaycoords
  Use simcell, Only: SimCell_Params
  Use branchedmcparams, Only: brmcparams_initgcmcparams, &
      brmcparams_adjustdeltatrans, brmcparams_adjustdeltarot, &
      brmcparams_setconstrnormal, brmcparams_gettransmat, &
      BranchedMolecGCMCParams, BranchedMolecNVTParams, init, &
      display, brmcparams_displaynvtstats
#else
  Use defaults
  Use utils
  Use vector
  Use matrix
  Use rigidmcparams
  Use rigidgc
  Use simcell
  Use branchedmcparams
#endif

  Implicit None
  Save

#ifdef PRIVATE
  Private
  Public :: mcparams_init, GenCoordsGCMCParams, GenCoordsNVTParams
#endif

  Type Move_Stats
    Integer       :: att
    Integer       :: succ
  End Type Move_Stats
  
  Type GenCoordsGCMCParams
    Type(RigidMolecGCMCParams), Pointer     :: rigidgcmcparams
    Type(BranchedMolecGCMCParams), Pointer  :: branchedgcmcparams
  End Type GenCoordsGCMCParams

  Type GenCoordsNVTParams
    Type(RigidMolecNVTParams), Pointer      :: rigidnvtparams
    Type(BranchedMolecNVTParams), Pointer   :: branchednvtparams
 End Type GenCoordsNVTParams

  Interface mcparams_init
    Module Procedure mcparams_initgcmcparams
    Module Procedure mcparams_initnvtparams
  End Interface

  Interface Assignment(=)
    Module Procedure mcparams_copygcmcparams
    Module Procedure mcparams_copynvtparams
  End Interface

Contains
  !--------------------------------------------------------------
  ! Initializing the parameters required for GCMC simulations
  !--------------------------------------------------------------
  Subroutine mcparams_initgcmcparams(gcparams, modeltype, simcell, &
      sorbtype, tk, filename)
    Type(GenCoordsGCMCParams), Intent(inout)  :: gcparams
    Character(*), Intent(in)  :: modeltype
    Type(SimCell_Params), Intent(in) :: simcell
    Integer, Intent(in)       :: sorbtype
    Real(kind=RDbl), Intent(in) :: tk
    Character(*), Intent(in)    :: filename
    
    Integer      :: error

    !** Nullify all the pointer instances in gcmodel
    Nullify(gcparams%rigidgcmcparams)
    Nullify(gcparams%branchedgcmcparams)

    Select Case (Trim(toupper(modeltype)))
    Case('RIGID')
      Allocate(gcparams%rigidgcmcparams, STAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Could not allocate 'rigidparams'"
        Stop
      End If
      Call rigidmcparams_init(gcparams%rigidgcmcparams, &
          simcell, sorbtype, tk, filename)
    Case('BRANCHED')
      Allocate(gcparams%branchedgcmcparams, STAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Could not allocate 'branchedparams'"
        Stop
      End If
      Call brmcparams_initgcmcparams(gcparams%branchedgcmcparams, &
          simcell, sorbtype, tk, filename)
    Case Default
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not get parameters for modeltype ", Trim(modeltype)
      Stop
    End Select
  End Subroutine mcparams_initgcmcparams

  !--------------------------------------------------------------
  ! Initializing the parameters required for NVT MC simulations
  !--------------------------------------------------------------
  Subroutine mcparams_initnvtparams(gcparams, modeltype, sorbtype, filename)
    Type(GenCoordsNVTParams), Intent(inout)  :: gcparams
    Character(*), Intent(in)  :: modeltype
    Integer, Intent(in)       :: sorbtype
    Character(*), Intent(in)  :: filename
    
    Integer      :: error

    !** Nullify all the pointer instances in gcmodel
    Nullify(gcparams%rigidnvtparams)
    Nullify(gcparams%branchednvtparams)

    Select Case (Trim(toupper(modeltype)))
    Case('RIGID')
      Allocate(gcparams%rigidnvtparams, STAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Could not allocate 'rigidparams'"
        Stop
      End If
      Call rigidmcparams_init(gcparams%rigidnvtparams, sorbtype, filename)
    Case('BRANCHED')
      Allocate(gcparams%branchednvtparams, STAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Could not allocate 'branchedparams'"
        Stop
      End If
      Call init(gcparams%branchednvtparams, sorbtype, filename)
    Case Default
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not get parameters for modeltype ", Trim(modeltype)
      Stop
    End Select
  End Subroutine mcparams_initnvtparams

  !---------------------------------------------------------------
  ! This subroutine COPIES the CONTENTS of gcoords2 into gcoords1
  ! The default assignment operator just copies POINTERS
  !---------------------------------------------------------------
  Subroutine mcparams_copygcmcparams(gcoords1, gcoords2)
    Type(GenCoordsGCMCParams), Intent(inout) :: gcoords1
    Type(GenCoordsGCMCParams), Intent(in) :: gcoords2
 
    If (Associated(gcoords2%rigidgcmcparams)) Then
      gcoords1%branchedgcmcparams = gcoords2%branchedgcmcparams
      gcoords1%rigidgcmcparams = gcoords2%rigidgcmcparams
    End If
  End Subroutine mcparams_copygcmcparams
    

  !---------------------------------------------------------------
  ! This subroutine COPIES the CONTENTS of gcoords2 into gcoords1
  ! The default assignment operator just copies POINTERS
  !---------------------------------------------------------------
  Subroutine mcparams_copynvtparams(gcoords1, gcoords2)
    Type(GenCoordsNVTParams), Intent(inout) :: gcoords1
    Type(GenCoordsNVTParams), Intent(in) :: gcoords2
 
    If (Associated(gcoords2%rigidnvtparams)) Then
      gcoords1%rigidnvtparams = gcoords2%rigidnvtparams
      gcoords1%rigidnvtparams = gcoords2%rigidnvtparams
    End If
  End Subroutine mcparams_copynvtparams

  !----------------------------------------------------
  ! Adjust the rotation jump length for NVT parameters
  !----------------------------------------------------
  Subroutine mcparams_adjustdeltatrans(nvtparams, movstats)
    Type(GenCoordsNVTParams), Intent(inout)   :: nvtparams
    Type(Move_Stats), Intent(in)           :: movstats

    Real(kind=Rdbl)           :: ratio

    !** Calculate the acceptance ratio
    If (movstats%att /= 0) Then
      ratio = movstats%succ/(movstats%att*1.0_RDbl)
    Else
      Return
    End If

    If (Associated(nvtparams%rigidnvtparams)) Then
      Call rigidmcparams_adjustdeltatrans(nvtparams%rigidnvtparams, ratio)
    Else If (Associated(nvtparams%branchednvtparams)) Then
      Call brmcparams_adjustdeltatrans(nvtparams%branchednvtparams, &
          ratio)
    Endif
  End Subroutine mcparams_adjustdeltatrans

  !-----------------------------------------------------------
  ! Adjust the rotation jump length for NVT parameters
  !-----------------------------------------------------------
  Subroutine mcparams_adjustdeltarot(nvtparams, movstats)
    Type(GenCoordsNVTParams), Intent(inout)  :: nvtparams
    Type(Move_Stats), Intent(in)          :: movstats

    Real(kind=Rdbl)           :: ratio

    !** Calculate the acceptance ratio
    If (movstats%att /= 0) Then
      ratio = movstats%succ/(movstats%att*1.0_RDbl)
    Else
      Return
    End If

    If (Associated(nvtparams%rigidnvtparams)) Then
      Call rigidmcparams_adjustdeltarot(nvtparams%rigidnvtparams, ratio)
    Else If (Associated(nvtparams%branchednvtparams)) Then
      Call brmcparams_adjustdeltarot(nvtparams%branchednvtparams, ratio)
    Endif
  End Subroutine mcparams_adjustdeltarot
  
  !----------------------------------------------------------
  ! Set the constraint plane for the nvt parameters
  !----------------------------------------------------------
  Subroutine mcparams_setconstrnormal(nvtparams, normal)
    Type(GenCoordsNVTParams), Intent(inout)  :: nvtparams
    Type(VecType), Intent(in)             :: normal
    
    If (Associated(nvtparams%rigidnvtparams)) Then
      Call rigidmcparams_setconstrnormal(nvtparams%rigidnvtparams, normal)
    Else If (Associated(nvtparams%branchednvtparams)) Then
      Call brmcparams_setconstrnormal(nvtparams%branchednvtparams, &
          normal)
    Endif
  End Subroutine mcparams_setconstrnormal

  !--------------------------------------------------------------------
  ! Returns the transformation matrix that goes from the fixed system
  ! to the system where the normal is aligned with the 'x-axis'
  !-------------------------------------------------------------------
  Type(MatrixType) Function mcparams_gettransmat(nvtparams)
    Type(GenCoordsNVTParams), Intent(in)  :: nvtparams

    If (Associated(nvtparams%rigidnvtparams)) Then
      mcparams_gettransmat = &
          rigidmcparams_gettransmat(nvtparams%rigidnvtparams)
    Else If (Associated(nvtparams%branchednvtparams)) Then
      mcparams_gettransmat = &
          brmcparams_gettransmat(nvtparams%branchednvtparams)
    Endif
  End Function mcparams_gettransmat

  !-----------------------------------------------------------------
  ! Calls the appropriate display routine for the GCMC parameters
  !-----------------------------------------------------------------
  Subroutine mcparams_displaygcmcparams(gcparams, nspc, optunit)
    Type(GenCoordsGCMCParams), Intent(in) :: gcparams
    Integer, Intent(in)   :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer    :: unitno

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    If (Associated(gcparams%rigidgcmcparams)) Then
      Call rigidmcparams_display(gcparams%rigidgcmcparams, nspc, optunit)
    Else If (Associated(gcparams%branchedgcmcparams)) Then
      Call display(gcparams%branchedgcmcparams, nspc, optunit)
    End If

  End Subroutine mcparams_displaygcmcparams

  !-----------------------------------------------------------
  ! Calls the appropriate display routine for NVTMC parameters
  !-----------------------------------------------------------
  Subroutine mcparams_displaynvtstats(nvtparams, nspc, optunit)
    Type(GenCoordsNVTParams), Intent(in) :: nvtparams
    Integer, Intent(in)   :: nspc
    Integer, Optional, Intent(in) :: optunit
    
    Character(len=nspc) :: spc
    Integer    :: unitno
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    If (Associated(nvtparams%rigidnvtparams)) Then
      Call rigidmcparams_displaynvtstats(nvtparams%rigidnvtparams, &
          nspc, optunit)
    Else If (Associated(nvtparams%branchednvtparams)) Then
      Call brmcparams_displaynvtstats(nvtparams%branchednvtparams, &
          nspc, optunit)
    Endif
  End Subroutine mcparams_displaynvtstats

  !-----------------------------------------------------------
  ! Calls the appropriate display routine for NVTMC parameters
  !-----------------------------------------------------------
  Subroutine mcparams_displaynvtparams(nvtparams, nspc, optunit)
    Type(GenCoordsNVTParams), Intent(in) :: nvtparams
    Integer, Intent(in)   :: nspc
    Integer, Optional, Intent(in) :: optunit
    
    Character(len=nspc) :: spc
    Integer    :: unitno
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If
    
    If (Associated(nvtparams%rigidnvtparams)) Then
      Call rigidmcparams_display(nvtparams%rigidnvtparams, nspc, optunit)
    Else If (Associated(nvtparams%branchednvtparams)) Then
      Call display(nvtparams%branchednvtparams, nspc, optunit)
    Endif
  End Subroutine mcparams_displaynvtparams

End Module mcparams
