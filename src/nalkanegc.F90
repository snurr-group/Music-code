!------------------------------------------------------------
! This module defines the structure for holding generalized
! coordinates of united atom n-alkanes.
!------------------------------------------------------------
Module nalkanegc

  Use defaults, Only: RDbl
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)

  Implicit None
  Save

  Private
  Public :: nAlkane, nalkanegc_init

  Type nAlkane
    Type(VecType)     :: startpos
  End Type nAlkane

!!$  Interface init
!!$    Module Procedure nalkanegc_init
!!$  End Interface

Contains
  Subroutine nalkanegc_init(gcoords, sorbtype)
    Type(nAlkane), Pointer  :: gcoords
    Integer, Intent(in)     :: sorbtype

    Integer  :: error

    Allocate(gcoords, STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          "Could not allocate 'gcmodel'"
      Stop
    Endif
    gcoords%startpos = 0.0_RDbl
  End Subroutine nalkanegc_init

End Module nalkanegc

  
