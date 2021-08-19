!------------------------------------------------------------------------------
! This routine calculates the torsion angle potential and 
! force.  It works for the Jorgensen function type:
! V = SUM(k=1,3) [ V_k/2(1 +/- Cos(k phi) ]  (k=1,3 [+]; k=2 [-])
!   Where phi is the angle between the two plane normals
!   at phi = Pi, V = 0.0
! 
!         i+1       i-1
!          o         o      plane1: defined by bonds bi-1, bi
!         / \       /       plane2: defined by bonds bi, bi+1
!   bi+1 /   \ bi  / bi-1     
!       /     \   /          
!      /       \ /           
! i+2 o       i o     
!
! In the 'trans' conformation above, the angle is zero because
! the two plane normals are parallel.  In the 'gauche' conformation,
! it is Pi.
! 
! Note: since the Ryckaert/June parameters and the Jorgensen parameters
!       are inter-relatable, this module converts the Jorgensen parameters
!       to the other type and uses 'cosexpansion' to calculate the interactions
!------------------------------------------------------------------------------

Module torjorg

  Use defaults, Only: strLen, RDbl, lstrLen, one, zero
  Use utils, Only: toreal
  Use vector, Only: VecType
  Use cosexpansion, Only: CosExpansionModel, cosexpansion_getinteraction, &
      cosexpansion_cleanup

  Implicit None
  Save

  Private
  Public :: JorgTorsionModel, torjorg_idstring, torjorg_init, &
      torjorg_display, torjorg_clean, torjorg_getinteraction

  Character(len=strLen), Parameter  :: torjorg_idstring = 'JORGTORSION'
  
  Type JorgTorsionModel
    Real(kind=RDbl), Dimension(3)  :: vi
    Type(CosExpansionModel)        :: cosexp_params
  End Type JorgTorsionModel

Contains
  
  !--------------------------------------------------------------
  ! Initialize the cosexpansion parameters directly
  ! Requires:  params -- Jorgensen torsion parameters to init
  !            info -- strings containing the three parameters
  !--------------------------------------------------------------
  Subroutine torjorg_init(params,info)
    Type(JorgTorsionModel)           :: params
    Character(*), Dimension(:)       :: info

    Integer      :: i

    !** Convert strings to V_k's
    Do i = 1,3
      params%vi(i) = toreal(info(i))
    End Do

    !** Convert V_k's to cosexpansion parameters
    params%cosexp_params%cn(0) = 0.5_RDbl*(params%vi(1) + &
        3.0_RDbl*params%vi(2) + params%vi(3))
    params%cosexp_params%cn(1) = 0.5_RDbl*(3.0_RDbl*params%vi(3) - params%vi(1))
    params%cosexp_params%cn(2) = -1.0_RDbl*params%vi(2)
    params%cosexp_params%cn(3) = -2.0_RDbl*params%vi(3)
    params%cosexp_params%cn(4) = 0.0_RDbl
    params%cosexp_params%cn(5) = 0.0_RDbl
    
  End Subroutine torjorg_init
   
  !------------------------------------------------------------------------
  ! Call the cosexpansion routine to get the interactions
  ! Requires:  params -- Jorgensen torsion parameters
  !            r -- position vectors
  !            u -- potential energy
  !            vf -- force on each atom (optional)
  !------------------------------------------------------------------------
  Subroutine torjorg_getinteraction(params,r,ifvec,u,vf)
    Type(JorgTorsionModel), Intent(In)                 :: params
    Type(VecType), Intent(In), Dimension(:)            :: r
    Logical, Intent(In) :: ifvec
    Real(kind=RDbl), Intent(Out)                       :: u
    Type(VecType), Intent(Out), Dimension(:), Optional :: vf

    If (Present(vf)) Then
      Call cosexpansion_getinteraction(params%cosexp_params,r,ifvec,u,vf)
    Else
      Call cosexpansion_getinteraction(params%cosexp_params,r,ifvec,u)
    End If

  End Subroutine torjorg_getinteraction

  !-------------------------------------------------------------------
  ! Dumps the contents of the Jorgensen torsion parameters to string
  ! Requires:  params -- Jorgensen torsion parameters
  !-------------------------------------------------------------------
  Function torjorg_display(params)
    Character(lstrLen)          :: torjorg_display
    Type(JorgTorsionModel)      :: params

    Integer :: i
    
    Write(torjorg_display,'(1x,3f7.3,a)') (params%vi(i), i=1,3),' kcal/mol'

  End Function torjorg_display

  
  !---------------------------------------------------------------
  ! Clean-up the parameters
  ! Requires:  params -- Jorgensen torsion parameters to clean
  !---------------------------------------------------------------
  Subroutine torjorg_clean(params)
    Type(JorgTorsionModel), Intent(InOut)  :: params

    Integer :: error

    Call cosexpansion_cleanup(params%cosexp_params)
    
  End Subroutine torjorg_clean
  
End Module torjorg


