!------------------------------------------------------------------------------
! This module returns a random vector
!------------------------------------------------------------------------------
Module randomvec

  Use defaults, Only: RDbl
  Use random, Only: random_gaussian
  Use vector, Only: VecType

  Implicit None
  Save

  Private
  Public :: randomvec_gaussian

Contains

  !----------------------------------------------------------------
  ! Returns a vector of random numbers from a gaussian distribution
  !----------------------------------------------------------------
  Type(VecType) Function randomvec_gaussian(ave,sigma)
    Real(Kind=RDbl), Intent(In) :: ave, sigma
    Integer :: i

    Do i = 1,3
      randomvec_gaussian%comp(i) = random_gaussian(ave,sigma)
    End Do

  End Function randomvec_gaussian

End Module randomvec
