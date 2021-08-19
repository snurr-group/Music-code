Module branchedmcparams

  Use defaults, Only: strLen, RDbl, dashedline, Pi
  Use file, Only: file_getunit
  Use utils, Only:
  Use vector, Only: VecType, mag, vector_display, vector_getnormsq, &
      vector_getcomp, vector_iscollinear, vector_crossprod, &
      Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnorm
  Use matrix, Only: MatrixType, matrix_getinv, matrix_genrotcy, &
      matrix_genrotcz, Operator(*), Assignment(=)
  Use bmap, Only: BiasMap_Params, bmap_init, bmap_display
  Use simcell, Only: SimCell_Params, simcell_getzeoell, simcell_geteff
  Use molecules, Only: molecules_getcomdefn, molecules_getnatoms

  Implicit None
  Save

  Private
  Public :: BranchedMolecGCMCParams, BranchedMolecNVTParams, &
      brmcparams_gettransmat, brmcparams_initgcmcparams, &
      brmcparams_adjustdeltatrans, brmcparams_adjustdeltarot, &
      brmcparams_setconstrnormal, init, display, &
      brmcparams_displaynvtstats

  Type BranchedMolecGCMCParams
    Type(BiasMap_Params), Pointer :: bmap
    Type(SimCell_Params), Pointer :: simcell
    Character(len=strLen)         :: biasfilename
    Type(VecType), Dimension(:), Pointer   :: com_defn
  End Type BranchedMolecGCMCParams

  Type BranchedMolecNVTParams
    Logical           :: constrained
    Type(MatrixType)  :: constrtrans ! Transforming to the constrained system
    Type(MatrixType)  :: constrinvtrans
    Type(VecType)     :: normal
    Real(kind=RDbl)   :: deltatrans
    Real(kind=RDbl)   :: deltarot
    Type(VecType), Dimension(:), Pointer   :: com_defn
  End Type BranchedMolecNVTParams

  
  Interface init
    Module Procedure brmcparams_initnvtparams
    Module Procedure brmcparams_initgcmcparams
  End Interface

  Interface display
    Module Procedure brmcparams_displaynvtparams
    Module Procedure brmcparams_displaygcmcparams
  End Interface

Contains
  !------------------------------------------------------
  ! Read the information that will be required to do GCMC
  ! type of moves.
  !------------------------------------------------------
  Subroutine brmcparams_initgcmcparams(gcmodelparams, simcell, sorbtype, &
      tk, filename)
    Type(BranchedMolecGCMCParams), Intent(inout) :: gcmodelparams
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

  End Subroutine brmcparams_initgcmcparams

  !------------------------------------------------------
  ! Read the information that will be required to do NVT
  ! Monte Carlo type of moves.
  !------------------------------------------------------
  Subroutine brmcparams_initnvtparams(gcmodelparams, sorbtype, filename)
    Type(BranchedMolecNVTParams), Intent(inout) :: gcmodelparams
    Integer, Intent(in)        :: sorbtype
    Character(*), Intent(in)   :: filename

    Integer    :: unitno, natoms, error

    !** Get the delta displacements
    unitno = file_getunit(filename)
    Read(unitno,*) gcmodelparams%deltatrans
    Read(unitno,*) gcmodelparams%deltarot
    
    !** Set the constrained flag to .False.
    gcmodelparams%constrained = .False.

    !** Initialize the molecular definition
    natoms = molecules_getnatoms(sorbtype)
    Allocate(gcmodelparams%com_defn(natoms), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Could not allocate memory for 'com_defn'"
      Stop
    End If
    gcmodelparams%com_defn = molecules_getcomdefn(sorbtype)
  End Subroutine brmcparams_initnvtparams
  
  !------------------------------------------------------------------
  ! Read the information that will be required to do NVT
  ! Monte Carlo constrained to a plane type of moves.
  !------------------------------------------------------------------
  Subroutine brmcparams_setconstrnormal(gcmodelparams, normal)
    Type(BranchedMolecNVTParams), Intent(inout) :: gcmodelparams
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
  End Subroutine brmcparams_setconstrnormal

  !----------------------------------------------------------------
  ! Generates the transformation matrices that go between the fixed
  ! coordinate system and the coordinate system where the normal to
  ! constrained plane is along the x-axis.
  !----------------------------------------------------------------
  Subroutine mcparams_genconstrmatrices(gcmodelparams)
    Type(BranchedMolecNVTParams), Intent(inout) :: gcmodelparams
    
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
  Type(MatrixType) Function brmcparams_gettransmat(nvtparams)
    Type(BranchedMolecNVTParams), Intent(in) :: nvtparams
    
    brmcparams_gettransmat = nvtparams%constrtrans
  End Function brmcparams_gettransmat

  !------------------------------------------------------------
  ! Adjust the jump length to try and keep the acceptance ratio
  ! close to 50%
  !------------------------------------------------------------
  Subroutine brmcparams_adjustdeltatrans(nvtparams, ratio)
    Type(BranchedMolecNVTParams), Intent(inout) :: nvtparams
    Real(kind=RDbl), Intent(in)       :: ratio
    
    If (ratio < 0.49) Then
      nvtparams%deltatrans = Max(nvtparams%deltatrans*0.95_RDbl, 0.01_RDbl)
    Else If (ratio > 0.51) Then
      nvtparams%deltatrans = Min(nvtparams%deltatrans*1.05_RDbl, 1.0_RDbl)
    Endif
  End Subroutine brmcparams_adjustdeltatrans

  !------------------------------------------------------------
  ! Adjust the jump length to try and keep the acceptance ratio
  ! close to 50%
  !------------------------------------------------------------
  Subroutine brmcparams_adjustdeltarot(nvtparams, ratio)
    Type(BranchedMolecNVTParams), Intent(inout) :: nvtparams
    Real(kind=RDbl), Intent(in)         :: ratio
    
    !** Calculate the acceptance ratio
    If (ratio < 0.49) Then
      nvtparams%deltarot = Max(nvtparams%deltarot*0.95_RDbl, 0.01_RDbl)
    Else If (ratio > 0.51) Then
      nvtparams%deltarot = Min(nvtparams%deltarot*1.05_RDbl, Pi)
    Endif
  End Subroutine brmcparams_adjustdeltarot

  !------------------------------------------------------------
  ! Display the GCMC parameters
  !------------------------------------------------------------
  Subroutine brmcparams_displaygcmcparams(gcparams, nspc, optunit)
    Type(BranchedMolecGCMCParams), Intent(in) :: gcparams
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
    
  End Subroutine brmcparams_displaygcmcparams

  !--------------------------------------------------------
  ! Just displays the generalized coordinate parameters
  ! for NVT MC simulations
  !--------------------------------------------------------
  Subroutine brmcparams_displaynvtstats(nvtparams, nspc, optunit)
    Type(BranchedMolecNVTParams), Intent(in) :: nvtparams
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
  End Subroutine brmcparams_displaynvtstats

  !------------------------------------------------------------
  ! Display the NVTMC parameters
  !------------------------------------------------------------
  Subroutine brmcparams_displaynvtparams(nvtparams, nspc, optunit)
    Type(BranchedMolecNVTParams), Intent(in) :: nvtparams
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
  End Subroutine brmcparams_displaynvtparams


End Module branchedmcparams


