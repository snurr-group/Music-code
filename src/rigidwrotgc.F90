!------------------------------------------------------------
! This module defines the structure for holding generalized
! coordinates of rigid species with one rotation angle flexible
! such as methanol.
!------------------------------------------------------------
Module rigidwrotgc

  Use defaults, Only: RDbl
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)

  Implicit None
  Save

  Private
  Public  :: RigidwRot, rigidwrotgc_init

  Type RigidwRot
    Type(VecType)     :: startpos
  End Type RigidwRot

!!$  Interface init
!!$    Module Procedure rigidwrot_init
!!$  End Interface

Contains
  Subroutine rigidwrotgc_init(gcoords, sorbtype)
    Type(RigidwRot), Pointer  :: gcoords
    Integer, Intent(in)       :: sorbtype

    Integer  :: error

    Allocate(gcoords, STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          "Could not allocate 'gcmodel'"
      Stop
    Endif
    gcoords%startpos = 0.0_RDbl
  End Subroutine rigidwrotgc_init
  
End Module rigidwrotgc

  
