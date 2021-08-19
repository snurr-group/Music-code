Module rigidmcparams

  Use defaults, Only: RDbl, strLen, Pi, dashedline
  Use file, Only: file_getunit
  Use vector, Only: VecType, vector_getnorm, vector_display, Assignment(=), &
      Operator(+), Operator(-), Operator(*), Operator(/)
  Use matrix, Only: MatrixType, matrix_genrotcz, matrix_getinv, &
      matrix_genrotcy, Operator(*), Assignment(=)
  Use bmap, Only: bmap_display, BiasMap_params, bmap_init
  Use simcell, Only: SimCell_Params, simcell_getzeoell, simcell_geteff
  Use molecules, Only: molecules_getnatoms, molecules_getcomdefn

  Implicit None
  Save

  Private
  Public :: RigidMolecGCMCParams, RigidMolecNVTParams,  &
      rigidmcparams_gettransmat, &
      rigidmcparams_init, rigidmcparams_adjustdeltatrans, &
      rigidmcparams_adjustdeltarot, rigidmcparams_setconstrnormal, &
      rigidmcparams_displaynvtstats, rigidmcparams_display
  
  Type RigidMolecGCMCParams
    Type(BiasMap_Params), Pointer :: bmap
    Type(SimCell_Params), Pointer :: simcell
    Character(len=strLen)         :: biasfilename
    Type(VecType), Dimension(:), Pointer   :: com_defn
  End Type RigidMolecGCMCParams

  Type RigidMolecNVTParams
    Logical           :: constrained
    Type(MatrixType)  :: constrtrans ! Transforming to the constrained system
    Type(MatrixType)  :: constrinvtrans
    Type(VecType)     :: normal
    Logical           :: scale_trans, scale_rot
    Real(kind=RDbl)   :: deltatrans
    Real(kind=RDbl)   :: deltarot
    Type(VecType), Dimension(:), Pointer   :: com_defn
  End Type RigidMolecNVTParams

  
  Interface rigidmcparams_init
    Module Procedure rigidmcparams_initnvtparams
    Module Procedure rigidmcparams_initgcmcparams
  End Interface

  Interface rigidmcparams_display
    Module Procedure rigidmcparams_displaynvtparams
    Module Procedure rigidmcparams_displaygcmcparams
  End Interface

Contains
  !------------------------------------------------------
  ! Read the information that will be required to do GCMC
  ! type of moves.
  !------------------------------------------------------
  Subroutine rigidmcparams_initgcmcparams(gcmodelparams, simcell, sorbtype, &
      tk, filename)
    Type(RigidMolecGCMCParams), Intent(inout) :: gcmodelparams
    Type(SimCell_Params), Target, Intent(in)  :: simcell
    Integer, Intent(in)         :: sorbtype
    Real(kind=RDbl), Intent(in) :: tk
    Character(*), Intent(in)    :: filename

    Character(len=strLen)   :: biasfilename
    Real(kind=RDbl), Dimension(3) :: eff
    Integer    :: unitno, natoms, error
    
    !** Set the simcell pointer
    gcmodelparams%simcell => simcell

    !** Read the bias filename
    unitno = file_getunit(filename)
    Read(unitno,*) biasfilename
    gcmodelparams%biasfilename = biasfilename

    !** Initialize the bias map
    eff = simcell_geteff(simcell, .True.)
    Call bmap_init(gcmodelparams%bmap, eff, tk, biasfilename, simcell)

    !** Initialize the molecular definition
    natoms = molecules_getnatoms(sorbtype)
    Allocate(gcmodelparams%com_defn(natoms), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Could not allocate memory for 'com_defn'"
      Stop
    End If
    gcmodelparams%com_defn = molecules_getcomdefn(sorbtype)

  End Subroutine rigidmcparams_initgcmcparams

  !------------------------------------------------------
  ! Read the information that will be required to do NVT
  ! Monte Carlo type of moves.
  !------------------------------------------------------
  Subroutine rigidmcparams_initnvtparams(gcmodelparams, sorbtype, filename)
    Type(RigidMolecNVTParams), Intent(inout) :: gcmodelparams
    Integer, Intent(in)        :: sorbtype
    Character(*), Intent(in)   :: filename

    Integer    :: unitno, natoms, error, scale_trans, scale_rot

    !** Get the delta displacements
    unitno = file_getunit(filename)
    Read(unitno,*) gcmodelparams%deltatrans, scale_trans
    Read(unitno,*) gcmodelparams%deltarot, scale_rot

    !** Set the constrained flag to .False.
    gcmodelparams%constrained = .False.

    !** Set the flag to decide whether we scale the jump lengths or not
    gcmodelparams%scale_trans = .False.
    gcmodelparams%scale_rot   = .False.
    If (scale_trans == 1) Then
      gcmodelparams%scale_trans = .True.
    End If
    If (scale_rot == 1) Then
      gcmodelparams%scale_rot = .True.
    End If

    !** Initialize the molecular definition
    natoms = molecules_getnatoms(sorbtype)
    Allocate(gcmodelparams%com_defn(natoms), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Could not allocate memory for 'com_defn'"
      Stop
    End If
    gcmodelparams%com_defn = molecules_getcomdefn(sorbtype)
  End Subroutine rigidmcparams_initnvtparams
  
  !------------------------------------------------------------------
  ! Read the information that will be required to do NVT
  ! Monte Carlo constrained to a plane type of moves.
  !------------------------------------------------------------------
  Subroutine rigidmcparams_setconstrnormal(gcmodelparams, normal)
    Type(RigidMolecNVTParams), Intent(inout) :: gcmodelparams
    Type(VecType), Intent(in)                :: normal

    Real(kind=RDbl)     :: norm

    gcmodelparams%normal = normal

    !** Normalize the normal and set the constrained flag
    norm = vector_getnorm(gcmodelparams%normal)
    ! Check if the norm is greater than 0.00001 else it is zero
    If (norm > 1.0e-5_RDbl) Then
      gcmodelparams%constrained = .True.
      gcmodelparams%normal = gcmodelparams%normal/norm
      Call mcparams_genconstrmatrices(gcmodelparams)
    Else
      gcmodelparams%constrained = .False.
    End If
  End Subroutine rigidmcparams_setconstrnormal

  !----------------------------------------------------------------
  ! Generates the transformation matrices that go between the fixed
  ! coordinate system and the coordinate system where the normal to
  ! constrained plane is along the x-axis.
  !----------------------------------------------------------------
  Subroutine mcparams_genconstrmatrices(gcmodelparams)
    Type(RigidMolecNVTParams), Intent(inout) :: gcmodelparams
    
    Type(VecType)     :: normal, normal_z, normal_zy
    Type(MatrixType)  :: mrotz, mroty
    Real(kind=RDbl)   :: sintheta, costheta
    Real(kind=RDbl), Dimension(3)  :: comp

    normal = gcmodelparams%normal
    comp   = normal

    !** First do a rotation to put the normal in the x-z plane
    costheta = comp(1)/Sqrt(comp(1)**2 + comp(2)**2)
    sintheta = comp(2)/Sqrt(comp(1)**2 + comp(2)**2)
    mrotz = matrix_genrotcz(sintheta, costheta)
    normal_z = mrotz*normal

    !** Now generate the matrix for rotation about the y-axis to align
    !** the x-axis with the transformed vector
    comp = normal_z
    costheta = comp(1)/Sqrt(comp(1)**2 + comp(3)**2)
    sintheta = comp(3)/Sqrt(comp(1)**2 + comp(3)**2)
    mroty = matrix_genrotcy(-sintheta, costheta)

    !** Now multiply the two transformation matrices to get the overal
    !** transformation matrix.  Also get its inverse
    gcmodelparams%constrtrans = mroty*mrotz
    gcmodelparams%constrinvtrans = matrix_getinv(gcmodelparams%constrtrans)
  End Subroutine mcparams_genconstrmatrices

  !--------------------------------------------------------------------
  ! Returns the transformation matrix that goes from the fixed system
  ! to the system where the normal is aligned with the 'x-axis'
  !-------------------------------------------------------------------
  Type(MatrixType) Function rigidmcparams_gettransmat(nvtparams)
    Type(RigidMolecNVTParams), Intent(in) :: nvtparams
    
    rigidmcparams_gettransmat = nvtparams%constrtrans
  End Function rigidmcparams_gettransmat

  !------------------------------------------------------------
  ! Adjust the jump length to try and keep the acceptance ratio
  ! close to 50%
  !------------------------------------------------------------
  Subroutine rigidmcparams_adjustdeltatrans(nvtparams, ratio)
    Type(RigidMolecNVTParams), Intent(inout) :: nvtparams
    Real(kind=RDbl), Intent(in)       :: ratio

    If (.Not. nvtparams%scale_trans) Return

    If (ratio < 0.49) Then
      nvtparams%deltatrans = Max(nvtparams%deltatrans*0.95_RDbl, 0.01_RDbl)
    Else If (ratio > 0.51) Then
      nvtparams%deltatrans = Min(nvtparams%deltatrans*1.05_RDbl, 1.0_RDbl)
    Endif
  End Subroutine rigidmcparams_adjustdeltatrans

  !------------------------------------------------------------
  ! Adjust the jump length to try and keep the acceptance ratio
  ! close to 50%
  !------------------------------------------------------------
  Subroutine rigidmcparams_adjustdeltarot(nvtparams, ratio)
    Type(RigidMolecNVTParams), Intent(inout) :: nvtparams
    Real(kind=RDbl), Intent(in)         :: ratio

    If (.Not. nvtparams%scale_rot) Return
    
    If (ratio < 0.49) Then
      nvtparams%deltarot = Max(nvtparams%deltarot*0.95_RDbl, 0.01_RDbl)
    Else If (ratio > 0.51) Then
      nvtparams%deltarot = Min(nvtparams%deltarot*1.05_RDbl, Pi)
    Endif
  End Subroutine rigidmcparams_adjustdeltarot

  !------------------------------------------------------------
  ! Display the GCMC parameters
  !------------------------------------------------------------
  Subroutine rigidmcparams_displaygcmcparams(gcparams, nspc, optunit)
    Type(RigidMolecGCMCParams), Intent(in) :: gcparams
    Integer, Intent(in)   :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer     :: i, j, unitno, nsorbs, sorbtype 

    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If  
    Write(unitno,'(2a)') spc, dashedline
    Write(unitno,'(2a)') spc, "Generalized Coordinates GCMC Parameters:"
    Write(unitno, '(a,2x,2a)') spc,"Bias Filename       : ",  &
        gcparams%biasfilename
    Call bmap_display(gcparams%bmap, nspc+2, unitno)
    
  End Subroutine rigidmcparams_displaygcmcparams

  !--------------------------------------------------------
  ! Just displays the generalized coordinate parameters
  ! for NVT MC simulations
  !--------------------------------------------------------
  Subroutine rigidmcparams_displaynvtstats(nvtparams, nspc, optunit)
    Type(RigidMolecNVTParams), Intent(in) :: nvtparams
    Integer, Intent(in)   :: nspc
    Integer, Optional, Intent(in) :: optunit

    Integer      :: i, unitno
    Character(len=nspc) :: spc

    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If  

    Write(unitno,'(a,2x,a,f8.3)') spc, "Delta Translate      : ", &
        nvtparams%deltatrans
    Write(unitno,'(a,2x,a,f8.3)') spc, "Delta Rotate         : ", &
        nvtparams%deltarot
  End Subroutine rigidmcparams_displaynvtstats

  !------------------------------------------------------------
  ! Display the NVTMC parameters
  !------------------------------------------------------------
  Subroutine rigidmcparams_displaynvtparams(nvtparams, nspc, optunit)
    Type(RigidMolecNVTParams), Intent(in) :: nvtparams
    Integer, Intent(in)   :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer     :: i, j, unitno, nsorbs, sorbtype 

    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do
    
    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If  
    Write(unitno,'(2a)') spc, dashedline
    Write(unitno,'(2a)') spc, "Generalized Coordinates NVTMC Parameters:"
    If (nvtparams%constrained) Then
      Write(unitno,'(a,2x,2a)') spc, "Constrained ?   : ", "True"
      Write(unitno,'(a,2x,2a)') spc,  "Normal          : ", &
          vector_display(nvtparams%normal, "f8.3")
    Else
      Write(unitno,'(a,2x,2a)') spc, "Constrained ?   : ", "False"
    End If
    Write(unitno,'(a,2x,a,f8.3)') spc, "Delta Translate : ", &
        nvtparams%deltatrans
    Write(unitno,'(a,2x,a,f8.3)') spc, "Delta Rotate    : ", &
        nvtparams%deltarot
  End Subroutine rigidmcparams_displaynvtparams


End Module rigidmcparams


