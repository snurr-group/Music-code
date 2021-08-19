!------------------------------------------------------------------------
! This module handles matching two different atomic  structures.  At the
! moment it only includes a routine for determining if two different
! structures are identical.
!
! For now, structures are specified using cartesian coordinates and an
! elemental symbol on each atom.  Using readstruc's Structure data type.
!
! Need Improvments:
! 1) addition of routines for matching structures using translational and
!    rotational degrees of freedom.
!------------------------------------------------------------------------

Module match

  Use defaults, Only: RDbl,lstrLen
  Use utils, Only: allocErrDisplay,deallocErrdisplay,toupper
  Use vector, Only: VecType,vector_getdistsq,vector_display
  Use readstruc, Only: Structure,readstruc_xform
  Use visxyz, Only: visxyz_dump, XYZ_Entry, visxyz_make

  Implicit None
  Save

  Private
  Public :: match_compare,match_display

Contains
  !----------------------------------------------------------------------------
  ! Compares two different atomic coordinate structures.  Order of atoms in 
  ! passed arrays are irrelevant.  Atoms are matched based on their elemental
  ! symbols and a difference tolerance between position vectors.  Function
  ! returns the first encountered atom of the second structure that does
  ! not have a match in the first structure.
  ! Requires:  struc1 -- structure 1
  !            struc2 -- structure 2
  !            tolerance -- tolerance for matching positions (vector norm)
  !            display -- if True, display the match
  !----------------------------------------------------------------------------
  Integer Function match_compare(struc1,struc2,tolerance,display)
    Type(Structure), Intent(In)     :: struc1,struc2           
    Real(kind=RDbl), Intent(In)     :: tolerance
    Logical, Intent(In)             :: display

    Integer                            :: a1,a2,natoms
    Real(kind=RDbl)                    :: tolsq,distsq
    Integer, Dimension(struc1%natoms)  :: map1to2,map2to1

    match_compare = 0

    !** Check that all arrays have the same size
    natoms = struc1%natoms
    If (struc1%natoms /= struc2%natoms) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' passed structures are sized inconsistently'
      Stop            
    End If

    !** initialize atom maps
    map1to2 = 0
    map2to1 = 0

    tolsq = tolerance**2

    !** loop over atoms in structure 2 and find matches
    Do a2 = 1,natoms  
      Do a1 = 1,natoms  
        If (map1to2(a1) == 0) Then
          distsq = vector_getdistsq(struc1%coords(a1),struc2%coords(a2))
!LC          Write(*,*) a1,a2,Sqrt(distsq)
!LC          Write(*,*) '   ',Trim(vector_display(struc1%coords(a1),'f8.3')), &
!LC              Trim(vector_display(struc2%coords(a2),'f8.3'))
          If ((distsq <= tolsq).And. &
              (Toupper(struc1%elements(a1)) == Toupper(struc2%elements(a2))))Then
            map1to2(a1) = a2
            map2to1(a2) = a1
          End If
        End If
      End Do
      If (map2to1(a2) == 0) Then  !** couldn't find a match for a2
        match_compare = a2
        Exit
      End If
    End Do

    !** Display the attempted match if desired
    If (display) Then
      Call match_display(struc1,struc2,6,map2to1)
    End If

  End Function match_compare

  !----------------------------------------------------------------------------
  ! Just dumps the two structures to the screen.  The routine will reorder
  ! the second set of atoms, if a mapping array is passed
  ! Requires:  struc1 -- structure 1
  !            struc2 -- structure 2
  !            unit -- unit to dump into
  !            map2to1 -- optional array mapping atom from 2 to atom in 1
  !----------------------------------------------------------------------------
  Subroutine match_display(struc1,struc2,unit,map2to1)
    Type(Structure), Intent(In)                  :: struc1,struc2           
    Integer, Intent(In)                          :: unit
    Integer, Dimension(:), Intent(In), Optional  :: map2to1

    Integer                     :: a1,a2,natoms,match
    Character(len=lstrLen)      :: string1,string2,string

    !** Check that all arrays have the same size
    natoms = struc1%natoms
    If (struc1%natoms /= struc2%natoms) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' passed structures are sized inconsistently'
      Stop            
    End If

    !** dump the coordinates
    Write(unit,'(8x,a,i5)') 'Number of atoms: ',natoms
    Write(unit,'(8x,a,25x,a)') 'Structure 1','Structure 2'
    Do a1 = 1,natoms
      string = vector_display(struc1%coords(a1),'f8.3')
      Write(string1,'(i4,2x,a,2x,a)') a1,struc1%elements(a1),Trim(string)
      If (Present(map2to1)) Then
        Do a2 = natoms,1,-1
          match = map2to1(a2)
          If (map2to1(a2) == a1) Exit
        End Do
        If (match == 0) Then
          string2 = 'could not find match'
        Else
          string = vector_display(struc2%coords(match),'f8.3')
          Write(string2,'(i4,2x,a,2x,a)') a2,struc2%elements(match),Trim(string)
        End If
      Else
        string = vector_display(struc2%coords(a1),'f8.3')
        Write(string2,'(i4,2x,a,2x,a)') a1,struc2%elements(a1),Trim(string)
      End If
      Write(unit,'(a,3x,a)') Trim(string1),Trim(string2)
    End Do

  End Subroutine match_display

End Module match
