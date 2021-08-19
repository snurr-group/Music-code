!-----------------------------------------------------------
! This will most probably contain the data structure
! defined by Jeff Endelman for non-ring branched molecules
!----------------------------------------------------------
Module branchedcoords 

  Use defaults, Only: RDbl, zero, pi, twopi, one, degTorad, & 
      MAX_ATOMS , strLen
  Use file, Only: 
  Use random, Only: rranf
  Use utils, Only: swap 
  Use vector, Only: VecType, mag, vector_display, vector_getnormsq, &
      vector_getcomp, vector_iscollinear, vector_crossprod, &
      Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnorm       
  Use matrix, Only: MatrixType, Assignment(=), Operator(*) 
  Use atom, Only: atom_getname      
  Use molecule, Only: MolecularParams, molecule_getnatoms, &
      molecule_getnthatom, molecule_getnthstereo  
  Use molecules, Only: molecules_getmolecule, molecules_getnatoms  
  Use connects, Only: AtomList, connects_nconnected, connects_nthcon 

  Implicit None
  Save

  Private
  Public :: BranchedMolec, NodeType, branchedcoords_initcoords, &
      branchedcoords_toxyz, branchedcoords_place, &
      branchedcoords_unplace, branchedcoords_placenode, & 
      branchedcoords_coord,  & 
      branchedcoords_build_tree, branchedcoords_set_stereo,& 
      branchedcoords_find_min, branchedcoords_determine_stereo, &
      branchedcoords_equil_ba, branchedcoords_update_ta, &
      branchedcoords_display, branchedcoords_readrestartinfo, & 
      branchedcoords_dumprestartinfo, branchedcoords_writenodeinfo, &
      branchedcoords_copycoords  , branchedcoords_AssignNode, &
      branchedcoords_dumpNode, branchedcoords_findnode, &
      branchedcoords_markhigheratoms, branchedcoords_copyvalues

  ! The branched molecule is stored as a doubly linked tree.
  ! Each node has a predescessor (P) and three descendants:
  ! center (C), left (L), and right (R).

  ! The orientation of these children is important to stereochemistry.
  ! If the P atom is placed in the back, then a clockwise reading
  ! of the children would be C, L, R

  Type NodeType
    Type(NodeType), Pointer   :: P
    Type(NodeType), Pointer   :: C
    Type(NodeType), Pointer   :: L
    Type(NodeType), Pointer   :: R
    Integer                   :: atom_num
    Integer                   :: angle_index ! index in the angle library
    Integer                   :: atom_type
    Integer                   :: placed  ! 1 = has been placed, 0 = no
    Real(kind=RDbl)           :: bond_length      

    Character                 :: stereo
    ! 'S' and 'R' have same meaning as in chemistry, except that instead
    ! of atomic number, atom_num is used (higher number = higher priority).
    ! Thus, if we place lowest atom_num in the back, 'S' means we can read
    ! the atom numbers in descending order clockwise.  
    ! 'N' means no stereochemistry, or racemic mixture

    ! all angles in radians    
    ! ba = bond angle, ta = torsion angle

    Real(kind=RDbl)           :: baPC
    Real(kind=RDbl)           :: baPL
    Real(kind=RDbl)           :: baPR
    Real(kind=RDbl)           :: baCL
    Real(kind=RDbl)           :: baCR
    Real(kind=RDbl)           :: baLR
    Real(kind=RDbl)           :: taC
    Real(kind=RDbl)           :: taL
    Real(kind=RDbl)           :: taR

    ! ***Frame library for Node having more than two 
    ! connections  
    Integer                     ::  frame_index 

  End Type NodeType

  Type BranchedMolec
    Type(VecType)             :: startpos
    Type(NodeType), Pointer   :: root
  End Type BranchedMolec

  Interface init
    Module Procedure branchedcoords_initcoords
  End Interface


Contains

  !------------------------------------------------------------------------
  ! Initializes the generalized coordinates as a tree 
  !------------------------------------------------------------------------

  Subroutine branchedcoords_initcoords(gcoords, sorbtype)

    Type(BranchedMolec),Pointer   :: gcoords
    Integer, Intent(in)           :: sorbtype

    Integer                       :: error
    Integer                       :: num_atoms
    Type(MolecularParams)         :: molecule

    molecule = molecules_getmolecule(sorbtype)
    num_atoms = molecule_getnatoms(molecule)

    Allocate(gcoords, STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          "Could not allocate 'gcoords'"
      Stop
    Endif
    gcoords%startpos = 0.0_RDbl

    ! Begin allocating memory for the tree
    Allocate(gcoords%root, stat=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          "Could not allocate 'gcoords%root'"   
      Stop
    End If

    Nullify(gcoords%root%P, gcoords%root%C, gcoords%root%L, gcoords%root%R)
    gcoords%root%atom_num = 1  
    gcoords%root%bond_length = 0.0_RDbl
    Nullify(gcoords%root%P)

    ! This will build the entire tree by recursion

    Call branchedcoords_build_tree(molecule, gcoords%root)

    ! Build the proper stereochemistry into the tree
    Call branchedcoords_set_stereo(gcoords%root)

    ! Set bond angles to equil values
    Call branchedcoords_equil_ba(molecule, gcoords%root)

    ! Make torsion angles consistent
    Call branchedcoords_update_ta(gcoords%root)

  End Subroutine branchedcoords_initcoords

  !-------------------------------------------------
  !  Copies the branchedcoords from incoords to outcoords
  !-------------------------------------------------
  Subroutine branchedcoords_copycoords(outcoords,incoords)
    Type(BranchedMolec), Intent(inout)  :: outcoords
    Type(BranchedMolec), Intent(in)     :: incoords

    outcoords%startpos = incoords%startpos

    If (Associated(incoords%root)) Then
      Call branchedcoords_copynode(outcoords%root, incoords%root)
    End If

  End Subroutine branchedcoords_copycoords

  !-----------------------------------------------
  !  Copies the tree structure
  !----------------------------------------------
  Recursive Subroutine branchedcoords_copynode(outptr, inptr)

    Type(NodeType), Pointer  ::  inptr, outptr

    outptr%atom_num = inptr%atom_num
    outptr%atom_type = inptr%atom_type
    outptr%placed = inptr%placed
    outptr%bond_length = inptr%bond_length
    outptr%stereo = inptr%stereo
    outptr%baPC = inptr%baPC
    outptr%baPL = inptr%baPL
    outptr%baPR = inptr%baPR
    outptr%baCL = inptr%baCL
    outptr%baCR = inptr%baCR
    outptr%baLR = inptr%baLR
    outptr%taC = inptr%taC
    outptr%taL = inptr%taL
    outptr%taR = inptr%taR

    If (Associated(inptr%C)) Call branchedcoords_copynode(outptr%C,inptr%C)
    If (Associated(inptr%L)) Call branchedcoords_copynode(outptr%L,inptr%L)
    If (Associated(inptr%R)) Call branchedcoords_copynode(outptr%R,inptr%R) 

  End Subroutine branchedcoords_copynode


  !-------------------------------------------------------------------------
  !  Copies the tree structure from one molec type to another molec type
  !  lastptr is the last associated pointer for the shorter of these molecules
  !  does not copy : atom_num, atom_type, bondlength
  !-------------------------------------------------------------------------
  Recursive Subroutine branchedcoords_copyvalues(inptr, outptr, lastinptr, &
      lastoutptr)

    Type(NodeType), Pointer  ::  inptr, outptr, lastinptr, lastoutptr
    Nullify(lastinptr)
    Nullify(lastoutptr)
!    outptr%atom_num = inptr%atom_num
!    outptr%atom_type = inptr%atom_type
    outptr%placed = inptr%placed
!    outptr%bond_length = inptr%bond_length
    outptr%stereo = inptr%stereo
    outptr%baPC = inptr%baPC
    outptr%baPL = inptr%baPL
    outptr%baPR = inptr%baPR
    outptr%baCL = inptr%baCL
    outptr%baCR = inptr%baCR
    outptr%baLR = inptr%baLR
    outptr%taC = inptr%taC
    outptr%taL = inptr%taL
    outptr%taR = inptr%taR

    If ( (Associated(inptr%C)) .and.  (Associated(outptr%C)) ) Then
      Call branchedcoords_copyvalues(inptr%C, outptr%C, lastinptr, lastoutptr)
    Else
      lastinptr=>inptr
      lastoutptr=>outptr
    Endif
    If ( (Associated(inptr%L)) .Or. (Associated(inptr%R)) ) Then
      Write(*,*) " this routine works only for linear alkanes"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
  End Subroutine branchedcoords_copyvalues



  !-----------------------------------------------
  !  Makes atoms placed
  !-----------------------------------------------

  Recursive Subroutine branchedcoords_place(ptr)

    Type(NodeType), Pointer         :: ptr

    If (Associated(ptr)) Then
      ptr%placed = 1
      If (Associated(ptr%C)) Then
        Call branchedcoords_place(ptr%C)
        If (Associated(ptr%L)) Then
          Call branchedcoords_place(ptr%L)
          If (Associated(ptr%R)) Then
            Call branchedcoords_place(ptr%R)
          End If
        End If
      End If

    End If
  End Subroutine branchedcoords_place


  !-----------------------------------------------
  ! marks atoms higher than atomnum in the tree
  ! If higher then .True.
  ! higher means those which get "placed" earlier than this atom
  ! - atomlist should be .False. when this is called from outside
  ! - cpu intensive dont call too often
  !-----------------------------------------------
  Recursive Subroutine branchedcoords_markhigheratoms(ptr,atomlist,atomnum,&
      found)
    Logical,Dimension(:),Intent(inout) :: atomlist
    Type(NodeType), Pointer         :: ptr
    Integer,Intent(in) ::atomnum
    Logical,Intent(out) :: found

    found=.False.

    If (.Not.Associated(ptr%P)) Then
      ! root ptr
      If (ptr%atom_num==1) Then
        atomlist=.False.
      Else
        Write(*,*) "first atom should be root ptr"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Endif

    If (Associated(ptr)) Then

      If (ptr%atom_num == atomnum) Then
        found=.True.
        Return
      Endif
      If (atomlist(ptr%atom_num))  Then
        Write(*,*) "this should be false"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Else
        atomlist(ptr%atom_num)=.True.
      Endif
      If (Associated(ptr%C)) Then
        Call branchedcoords_markhigheratoms(ptr%C,atomlist,atomnum,found)
        If (found)   Return

        If (Associated(ptr%L)) Then
          Call branchedcoords_markhigheratoms(ptr%L,atomlist,atomnum,found)
          If (found)   Return

          If (Associated(ptr%R)) Then
            Call branchedcoords_markhigheratoms(ptr%R,atomlist,&
                atomnum,found)
            If (found)   Return
          End If
        End If
      End If
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop

    End If

  End Subroutine branchedcoords_markhigheratoms



  !------------------------------------------------
  !  Unplaces the atoms after a node 
  !------------------------------------------------

  Recursive Subroutine branchedcoords_unplace(ptr)

    Type(NodeType), Pointer    :: ptr

    If (Associated(ptr)) Then
      ptr%placed = 0
      Call branchedcoords_unplace(ptr%C)
      Call branchedcoords_unplace(ptr%L)
      Call branchedcoords_unplace(ptr%R)
    End If
  End Subroutine branchedcoords_unplace

  !------------------------------------------------
  !  Places a node only  
  !------------------------------------------------

  Recursive Subroutine branchedcoords_placenode(ptr)

    Type(NodeType), Pointer    :: ptr

    If (Associated(ptr)) Then
      ptr%placed = 1
    End If
  End Subroutine branchedcoords_placenode

  !----------------------------------------------------------------
  !  Rotates local coordinates onto global coordinate frame
  !  First call should pass identity matrix in global_rotation
  !  Right now this only works if vec_size = 3.  
  !  Terminates when reaches an atom that has not been placed.
  !----------------------------------------------------------------

  Recursive Subroutine branchedcoords_coord(ptr, xyzcoords, global_rotation)

    Type(NodeType), Pointer                        :: ptr
    Type(VecType), Dimension(:), Intent(InOut)     :: xyzcoords
    Type(MatrixType), Intent(in)                   :: global_rotation

    Type(MatrixType)                               :: local_rotation
    Type(VecType)                                  :: bond
    Real(Kind=RDbl)                                :: theta, phi

    If (Associated(ptr) .AND. ptr%placed == 1) Then
      If (Associated(ptr%C)) Then
        theta = pi - ptr%baPC
        phi = ptr%taC
        local_rotation = Reshape((/cos(theta),sin(theta)*cos(phi), &
            Sin(theta)*Sin(phi),sin(theta),-cos(theta)*cos(phi), &
            -Cos(theta)*Sin(phi),zero,sin(phi),-cos(phi)/), (/3, 3/)) 

        !     local_rotation = matrix_m_mult_m(global_rotation, local_rotation)
        local_rotation = global_rotation*local_rotation 

        bond = (/ptr%C%bond_length,zero,zero/)       

        xyzcoords(ptr%C%atom_num) = xyzcoords(ptr%atom_num) + &
            local_rotation*bond   

        If (Associated(ptr%C%C)) &    
            Call branchedcoords_coord(ptr%C, xyzcoords, local_rotation)

      End If

      If (Associated(ptr%L)) Then
        theta = pi - ptr%baPL
        phi = ptr%taL

        local_rotation = Reshape((/cos(theta),sin(theta)*cos(phi), &
            Sin(theta)*Sin(phi),sin(theta),-cos(theta)*cos(phi), &
            -Cos(theta)*Sin(phi),zero,sin(phi),-cos(phi)/), (/3, 3/)) 

        !     local_rotation = matrix_m_mult_m(global_rotation, local_rotation)
        local_rotation = global_rotation*local_rotation 

        bond = (/ptr%L%bond_length,zero,zero/)    

        xyzcoords(ptr%L%atom_num) = xyzcoords(ptr%atom_num) + &
            local_rotation*bond   

        Call branchedcoords_coord(ptr%L, xyzcoords, local_rotation)
      End If

      If (Associated(ptr%R)) Then
        theta = pi - ptr%baPR
        phi = ptr%taR

        local_rotation = Reshape((/cos(theta),sin(theta)*cos(phi), &
            Sin(theta)*Sin(phi),sin(theta),-cos(theta)*cos(phi), &
            -Cos(theta)*Sin(phi),zero,Sin(phi),-Cos(phi)/), (/3, 3/)) 

        !     local_rotation = matrix_m_mult_m(global_rotation, local_rotation)
        local_rotation = global_rotation*local_rotation     

        bond = (/ptr%R%bond_length,zero,zero/)   

        xyzcoords(ptr%R%atom_num) = xyzcoords(ptr%atom_num) + &
            local_rotation*bond   

        Call branchedcoords_coord(ptr%R, xyzcoords, local_rotation)
      End If
    End If
  End Subroutine branchedcoords_coord


  !-----------------------------------------------------------------------
  ! Generates xyz coordinates from generalized coordinates
  ! Only works for vec_size = 3
  !-----------------------------------------------------------------------

  Subroutine branchedcoords_toxyz(gcoords, xyzcoords)

    Type(BranchedMolec), Pointer                    :: gcoords
    Type(VecType), Dimension(:), Intent(out)        :: xyzcoords

    Type(MatrixType)                                :: identity

    ! Atom1 must be the root of the tree
    xyzcoords(1) = gcoords%startpos

    !   Call matrix_m_eq_2darray(identity, &
    !        Reshape((/one,zero,zero,zero,one,zero,zero,zero,one/),(/3,3/)))
    identity = Reshape((/one,zero,zero,zero,one,zero,zero,zero,& 
        one/),(/3,3/)) 

    Call branchedcoords_coord(gcoords%root,xyzcoords,identity)  

  End Subroutine branchedcoords_toxyz


  !------------------------------------------------------
  ! Returns a pointer (node) to the node with atom_num.
  ! Can only search below root.
  !------------------------------------------------------
  Recursive Subroutine branchedcoords_findnode(root, atom_num, node, found)

    Type(NodeType), Pointer    :: root, node
    Integer, Intent(In)        :: atom_num
    Logical, Intent(Out)       :: found  

    found = .False.

    If (Associated(root)) Then
      If (root%atom_num == atom_num) Then
        node => root
        found = .True.
      Else
        Call branchedcoords_findnode(root%C, atom_num, node, found)
        If (.Not.found) Then
          Call branchedcoords_findnode(root%L, atom_num, node, found)
        End If
        If (.Not.found) Then
          Call branchedcoords_findnode(root%R, atom_num, node, found)
        End If
      End If
    End If

  End Subroutine branchedcoords_findnode

  !--------------------------------------------------------------------
  ! Gets passed a pointer to a node, zeroes angles, connects 
  ! children, and enters in their information
  !--------------------------------------------------------------------
  Recursive Subroutine branchedcoords_build_tree(molecule, this_ptr)
    Type(MolecularParams)         :: molecule
    Type(NodeType), Pointer       :: this_ptr, new_ptr
    Integer                       :: i, error
    Integer                       :: num_connect
    Type(AtomList)                :: connection

    If (Associated(this_ptr)) Then
      this_ptr%stereo = molecule_getnthstereo(molecule,this_ptr%atom_num)
      this_ptr%atom_type = molecule_getnthatom(molecule,this_ptr%atom_num)
      this_ptr%placed = 0
      this_ptr%baPC = 0.0_RDbl
      this_ptr%baPL = 0.0_RDbl
      this_ptr%baPR = 0.0_RDbl
      this_ptr%baCL = 0.0_RDbl
      this_ptr%baCR = 0.0_RDbl
      this_ptr%baLR = 0.0_RDbl
      this_ptr%taL = 0.0_RDbl
      this_ptr%taC = 0.0_RDbl
      this_ptr%taR = 0.0_RDbl

      !      num_connect = molecule%atoms(this_ptr%atom_num)%nbonds
      !** for the special case of no connections
      If (molecule%natoms==1) Then
        num_connect=0
      Else
        num_connect = &
            connects_nconnected(molecule%connect(this_ptr%atom_num))  
      Endif

      Do i = 1, num_connect

        !** Get the ith connection
        Call connects_nthcon(molecule%connect(this_ptr%atom_num), &
            i, connection)

        !** Don't connect to previous atom
        If (Associated(this_ptr%P)) Then
          If (connection%atom == this_ptr%P%atom_num) Then
            Cycle  ! skip through this iteration of Do loop
          End If
        End If

        Allocate(new_ptr, STAT=error)
        If (error /= 0) Then
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              "Could not allocate 'new_ptr'"
          Stop
        Endif
        Nullify(new_ptr%P, new_ptr%C, new_ptr%L, new_ptr%R)

        new_ptr%bond_length = connection%length
        new_ptr%atom_num = connection%atom
        new_ptr%P => this_ptr

        If (.NOT.(Associated(this_ptr%C))) Then
          this_ptr%C => new_ptr
          Call branchedcoords_build_tree(molecule, new_ptr)
        Else If (.NOT.(Associated(this_ptr%L))) Then
          this_ptr%L => new_ptr
          Call branchedcoords_build_tree(molecule, new_ptr)
        Else
          this_ptr%R => new_ptr
          Call branchedcoords_build_tree(molecule, new_ptr)
        End If
      End Do
    End If

  End Subroutine branchedcoords_build_tree

  !---------------------------------------------------------------------------
  !  Switches left and right children and angles if necessary to effect the
  !  proper stereochemistry.  If there is no stereochemistry,
  !  orientation is chosen randomly (won't matter if not a stereocenter,
  !  otherwise leads to racemic mixture).
  !---------------------------------------------------------------------------

  Recursive Subroutine branchedcoords_set_stereo(this_ptr)

    Type(NodeType), Pointer  :: this_ptr
    Type(NodeType), Pointer  :: temp_ptr
    Character                :: current, lowest

    Nullify(temp_ptr)

    If (Associated(this_ptr)) Then

      ! Can't have stereochemistry at first atom.  Should be terminal.
      If (Associated(this_ptr%P)) Then   

        ! Check whether valence is at least 3 and not the last atom 
        !       If (Associated(this_ptr%L) .AND. Associated(this_ptr%C%C) ) Then
        If (Associated(this_ptr%L)) Then
          If( .NOT. Associated(this_ptr%C%C) ) Return 
          ! If there is no stereochemistry, randomly make one descendant 
          ! left and the other right (thus it's randomly S or R)
          If (this_ptr%stereo == "N" ) Then

            If (rranf() < 0.5) Then
              ! Swap pointers
              temp_ptr => this_ptr%L
              this_ptr%L => this_ptr%R
              this_ptr%R => temp_ptr

              ! Swap angles
              Call swap(this_ptr%baPL, this_ptr%baPR)
              Call swap(this_ptr%baCL, this_ptr%baCR)
              Call swap(this_ptr%taL, this_ptr%taR)
            End If
          Else

            If (Associated(this_ptr%R)) Then
              ! Valence is 4
              ! Find lowest priority atom
              lowest = branchedcoords_find_min(this_ptr)
            Else
              ! Valence is 3. So nonexistent right child is lowest priority
              lowest = 'R'
            End If

            ! Determine current stereochemistry
            Select Case (lowest)
              ! All calls to branchedcoords_determine_stereo must pass
              ! pointers in clockwise orientation with unnamed pointer
              ! in the back.  See NodeType definition.

            Case ('R')

              current = branchedcoords_determine_stereo(this_ptr%P%atom_num, &
                  this_ptr%L%atom_num, this_ptr%C%atom_num)

            Case ('P')
              current = branchedcoords_determine_stereo(this_ptr%C%atom_num, & 
                  this_ptr%L%atom_num, this_ptr%R%atom_num)

            Case ('C')
              current = branchedcoords_determine_stereo(this_ptr%P%atom_num, &
                  this_ptr%R%atom_num, this_ptr%L%atom_num)

            Case ('L')
              current = branchedcoords_determine_stereo(this_ptr%P%atom_num, &
                  this_ptr%C%atom_num, this_ptr%R%atom_num)

            Case Default
              Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
                  "Unassociated pointer failure in set_stereo"   
              Stop

            End Select

            ! If current stereochemistry is not the desired one, swap
            If (current /= this_ptr%stereo) Then
              ! Swap pointers
              temp_ptr => this_ptr%L
              this_ptr%L => this_ptr%R
              this_ptr%R => temp_ptr

              ! Swap angles
              Call swap(this_ptr%baPL, this_ptr%baPR)
              Call swap(this_ptr%baCL, this_ptr%baCR)
              Call swap(this_ptr%taL, this_ptr%taR)
            End If

          End If
        End If
      End If

      If (Associated(this_ptr%C)) Call branchedcoords_set_stereo(this_ptr%C)
      If (Associated(this_ptr%L)) Call branchedcoords_set_stereo(this_ptr%L)
      If (Associated(this_ptr%R)) Call branchedcoords_set_stereo(this_ptr%R)
    End If

  End Subroutine branchedcoords_set_stereo

  !--------------------------------------------------------------------------
  !  Returns a letter indicating which connection P,C,L,or R has
  !  the smallest atom_num.  this_ptr and this_ptr%P MUST be associated, 
  !  otherwise returns 'F'
  !--------------------------------------------------------------------------

  Character Function branchedcoords_find_min(this_ptr)

    Type(NodeType), Pointer    :: this_ptr

    Integer                    :: min
    Character                  :: letter

    If (Associated(this_ptr) .AND. Associated(this_ptr%P)) Then
      letter = 'P'
      min = this_ptr%P%atom_num

      If (Associated(this_ptr%C)) Then
        If (this_ptr%C%atom_num < min) Then
          letter = 'C'
          min = this_ptr%C%atom_num
        End If
      End If

      If (Associated(this_ptr%L)) Then
        If (this_ptr%L%atom_num < min) Then
          letter = 'L'
          min = this_ptr%L%atom_num
        End If
      End If

      If (Associated(this_ptr%R)) Then
        If (this_ptr%R%atom_num < min) Then
          letter = 'R'
          min = this_ptr%R%atom_num
        End If
      End If

    Else
      letter = 'F'
    End If

    branchedcoords_find_min = letter

  End Function branchedcoords_find_min



  !-------------------------------------------------------------------------
  ! Determines the stereochemistry for a sequence of three numbers, which
  ! are correspond to atom numbers in clockwise fashion.  If these numbers
  ! are descending clockwise, 'R' is returned.  Otherwise 'S'.
  !----------------------------------------------------------------------
  Character Function branchedcoords_determine_stereo(num1, num2, num3)

    ! Define a type to be used only in this subroutine.  
    ! It allows for construction of circular linked list.
    Type LoopType
      Integer                 :: value
      Type(LoopType), Pointer :: next
    End Type LoopType

    Integer, Intent(In)              :: num1, num2, num3
    Type(LoopType), Target           :: one, two, three
    Type(LoopType), Pointer          :: max

    ! Construct the circular list
    one%value = num1
    two%value = num2
    three%value = num3
    one%next => two
    two%next => three
    three%next => one

    ! Now find largest number
    max => one
    If (max%value < two%value) max => two
    If (max%value < three%value) max => three

    ! Are the numbers decreasing or increasing?
    If (max%next%value > max%next%next%value) Then
      branchedcoords_determine_stereo = 'R'
    Else 
      branchedcoords_determine_stereo = 'S'
    End If

  End Function branchedcoords_determine_stereo

  !----------------------------------------------------------
  !  Initialize bond angles to their equilibrium values.
  !----------------------------------------------------------

  Recursive Subroutine branchedcoords_equil_ba(molecule, ptr)

    Type(MolecularParams), Intent(IN)  :: molecule
    Type(NodeType), Pointer            :: ptr

    Integer                            :: atom,atomC,atomP,atomL,atomR


    If (Associated(ptr)) Then
      atom = molecule_getnthatom(molecule,ptr%atom_num)

      If (Associated(ptr%P)) Then
        atomP = molecule_getnthatom(molecule,ptr%P%atom_num)
      End If

      If (Associated(ptr%C)) Then
        atomC = molecule_getnthatom(molecule,ptr%C%atom_num)
      End If

      If (Associated(ptr%L)) Then
        atomL = molecule_getnthatom(molecule,ptr%L%atom_num)
      End If

      If (Associated(ptr%R)) Then
        atomR = molecule_getnthatom(molecule,ptr%R%atom_num)
      End If

      If (Associated(ptr%P) .AND. Associated(ptr%C)) Then
        !LC            ptr%baPC = degTorad* &
        !LC            molecule%bending%bbparams(atomP,atom,atomC)%hara%thetaeq
      End If

      If (Associated(ptr%P) .AND. Associated(ptr%L)) Then
        !LC            ptr%baPL = degTorad* &
        !LC            molecule%bending%bbparams(atomP,atom,atomL)%hara%thetaeq
      End If

      If (Associated(ptr%P) .AND. Associated(ptr%R)) Then
        !LC            ptr%baPR = degTorad* &
        !LC            molecule%bending%bbparams(atomP,atom,atomR)%hara%thetaeq
      End If

      If (Associated(ptr%C) .AND. Associated(ptr%L)) Then
        !LC            ptr%baCL = degTorad* &
        !LC            molecule%bending%bbparams(atomC,atom,atomL)%hara%thetaeq
      End If

      If (Associated(ptr%C) .AND. Associated(ptr%R)) Then
        !LC            ptr%baCR = degTorad* &
        !LC            molecule%bending%bbparams(atomC,atom,atomR)%hara%thetaeq
      End If

      If (Associated(ptr%L) .AND. Associated(ptr%R)) Then
        !LC            ptr%baLR = degTorad* &
        !LC            molecule%bending%bbparams(atomL,atom,atomR)%hara%thetaeq
      End If

      If (Associated(ptr%C)) Call branchedcoords_equil_ba(molecule, ptr%C)
      If (Associated(ptr%R)) Call branchedcoords_equil_ba(molecule, ptr%R)
      If (Associated(ptr%L)) Call branchedcoords_equil_ba(molecule, ptr%L)

    End If

  End Subroutine branchedcoords_equil_ba


  !-------------------------------------------------------------------
  !  Updates taL and taR to be geometrically consistent with baPC,
  !  baPL, baPR, baCL, and baCR.  baLR is ignored.
  !------------------------------------------------------------------

  Recursive Subroutine branchedcoords_update_ta(ptr)

    Type(NodeType), Pointer          :: ptr
    Real(Kind=RDbl)                  :: cosine_phi, sine_phi
    Real(Kind=RDbl)                  :: piece1, piece2, piece3

    If (Associated(ptr)) Then

      If (Associated(ptr%L)) Then

        piece1 = cos(ptr%baCL)-cos(ptr%baPC)*cos(ptr%baPL)
        piece2 = sqrt((sin(ptr%baPC)*sin(ptr%baPL))**2-(cos(ptr%baCL)- &
            cos(ptr%baPC)*cos(ptr%baPL))**2)
        piece3 = sin(ptr%baPC)*sin(ptr%baPL)
        cosine_phi = (piece1*cos(ptr%taC)+sin(ptr%taC)*piece2)/piece3

        If (cosine_phi > 1) cosine_phi = 1
        If (cosine_phi < -1) cosine_phi = -1

        sine_phi = (piece1*sin(ptr%taC)-cos(ptr%taC)*piece2)/piece3
        ptr%taL = acos(cosine_phi)

        ! Fortan uses principal branch of arccosine function.  Use sine to
        ! determine whether other sign is needed.
        If (sine_phi < 0) ptr%taL = -ptr%taL

        Call branchedcoords_update_ta(ptr%L)
      End If

      If (Associated(ptr%R)) Then

        piece1 = cos(ptr%baCR)-cos(ptr%baPC)*cos(ptr%baPR)
        piece2 = sqrt((sin(ptr%baPC)*sin(ptr%baPR))**2-(cos(ptr%baCR)- &
            cos(ptr%baPC)*cos(ptr%baPR))**2)
        piece3 = sin(ptr%baPC)*sin(ptr%baPR)

        cosine_phi = (piece1*cos(ptr%taC)-sin(ptr%taC)*piece2)/piece3

        If (cosine_phi > 1) cosine_phi = 1
        If (cosine_phi < -1) cosine_phi = -1

        sine_phi = (piece1*sin(ptr%taC)+cos(ptr%taC)*piece2)/piece3

        ptr%taR = acos(cosine_phi)

        If (sine_phi < 0) ptr%taR = -ptr%taR

        Call branchedcoords_update_ta(ptr%R)
      End If

      If (Associated(ptr%C)) Call branchedcoords_update_ta(ptr%C)
    End If

  End Subroutine branchedcoords_update_ta

  !-------------------------------------------------
  !  Displays generalized coordinates
  !-------------------------------------------------

  Subroutine branchedcoords_display(gcoords, molecno, unitno)

    Type(BranchedMolec), Pointer      :: gcoords
    Integer, Intent(In)               :: molecno, unitno

    If (Associated(gcoords)) Then
      Write(unitno,*) "Generalized Coordinates for Molecule #", molecno
      Write(unitno,*) "First Atom Position: ", &
          Trim(vector_display(gcoords%startpos,"f13.5"))
      Write(unitno,*)
      Call branchedcoords_display_loop(gcoords%root, unitno)
    End If

  End Subroutine branchedcoords_display

  !-----------------------------------------------------
  ! Displays data at every node below the node passed
  ! in the initial call.
  !-----------------------------------------------------

  Recursive Subroutine branchedcoords_display_loop(ptr, unitno)

    Type(NodeType), Pointer :: ptr
    Integer, Intent(in)     :: unitno

    If (Associated(ptr)) Then
75    Format(A,T25,F8.2)
      Write (unitno,'(A,T25,I8)') "Atom Number: ", ptr%atom_num
      Write (unitno,'(A,T18,A15)') "Atom Type: ", &
          atom_getname(ptr%atom_type)
      Write (unitno,'(A,T25,A8)') "Stereochemistry: ", ptr%stereo
      Write (unitno,fmt=75) "Bond Length: ", ptr%bond_length
      Write (unitno,fmt=75) "baPC: ", ptr%baPC
      Write (unitno,fmt=75) "baPL: ", ptr%baPL
      Write (unitno,fmt=75) "baPR: ", ptr%baPR
      Write (unitno,fmt=75) "baCL: ", ptr%baCL
      Write (unitno,fmt=75) "baCR: ", ptr%baCR
      Write (unitno,fmt=75) "baLR: ", ptr%baLR
      Write (unitno,fmt=75) "taC: ", ptr%taC
      Write (unitno,fmt=75) "taL: ", ptr%taL
      Write (unitno,fmt=75) "taR: ", ptr%taR
      Write (unitno,fmt=75) 
      Call branchedcoords_display_loop(ptr%C, unitno)
      Call branchedcoords_display_loop(ptr%L, unitno)
      Call branchedcoords_display_loop(ptr%R, unitno)
    End If

  End Subroutine branchedcoords_display_loop


  !----------------------------------------------------------------
  ! Read the generalized coordinates from the restart file
  !----------------------------------------------------------------
  Subroutine branchedcoords_readrestartinfo(gcoords, sorbtype,unitno)
    Type(BranchedMolec), Pointer   :: gcoords
    Integer, Intent(in)            :: sorbtype,unitno

    Real(kind=Rdbl)                :: x, y, z
    Real(kind=RDbl), Dimension(MAX_ATOMS)     :: angle1   
    Real(kind=RDbl), Dimension(MAX_ATOMS)     :: angle2   
    Integer, Dimension(MAX_ATOMS)  :: num      
    Integer                        :: i, natoms  
    If (.NOT. Associated(gcoords%root%C)) Then
      Read(unitno,*) x, y, z
      gcoords%startpos = (/x, y, z/)
      gcoords%root%placed=1
    ELseIf (.NOT. Associated(gcoords%root%C%C)) Then
      Read(unitno,*) x, y, z, gcoords%root%baPC, gcoords%root%taC, &
          gcoords%root%C%taC
      gcoords%startpos = (/x, y, z/)
      gcoords%root%placed=1
      gcoords%root%C%placed=1
    Else
      natoms = molecules_getnatoms(sorbtype)  
      Read(unitno,*) x, y, z
      Read(unitno,*) angle1(1), angle2(1), num(1)
      Read(unitno,*) (angle1(i),angle2(i),num(i),i=2,natoms-1)  
      gcoords%startpos = (/x, y, z/)
      !** Copying the values from angles arrays into Torsion and bond angles
      Call  branchedcoords_AssignNode(gcoords%root,angle1,angle2,num) 

    End If

  End Subroutine branchedcoords_readrestartinfo


  !-----------------------------------------------------------
  ! Assign the angles from the restartfile to the nodes  
  !-----------------------------------------------------------
  Recursive Subroutine branchedcoords_AssignNode(ptr,angle1,angle2,num) 
    Type(NodeType), Pointer        :: ptr
    Real(kind=RDbl), Dimension(:), Intent(in) :: angle1 
    Real(kind=RDbl), Dimension(:), Intent(in) :: angle2
    Integer, Dimension(:), Intent(in) :: num   

    Integer                        :: i   

    i  = ptr%atom_num  
    ptr%baPC= angle1(i)          
    ptr%taC= angle2(i) 

    !SDEBUG
    !!OBSOLETE
    ptr%angle_index= num(i)     
    !SDEBUG

    ptr%placed=1
    If( Associated(ptr%C)) Then 
      Call branchedcoords_AssignNode(ptr%C,angle1,angle2,num) 
    End If

  End Subroutine branchedcoords_AssignNode

  !-----------------------------------------------------------
  ! Write out the generalized branched coordinates to the restart file
  ! for one molecule
  !-----------------------------------------------------------
  Subroutine branchedcoords_dumprestartinfo(gcoords, sorbtype,unitno)
    Type(BranchedMolec), Pointer   :: gcoords
    Integer, Intent(in)            :: sorbtype,unitno

    Real(kind=RDbl), Dimension(MAX_ATOMS)     :: angle1   
    Real(kind=RDbl), Dimension(MAX_ATOMS)     :: angle2   
    Integer, Dimension(MAX_ATOMS)  :: num      
    Integer  ::  i,natoms  
    Character(len=strlen) :: temp_string

    !SDEBUG
    !*** This num() is nonsense, should not be here, 
    !*** and removed from restartfile
    !SDEBUG

    temp_string= Trim(vector_display(gcoords%startpos, "f13.5"))
    If (.NOT. Associated(gcoords%root%C)) Then
      Write(unitno,'(a, 2f13.5)') temp_string, & 
          gcoords%root%baPC, gcoords%root%taC

    ElseIf (.NOT. Associated(gcoords%root%C%C)) Then
      Write(unitno,'(a, 3f13.5)')  temp_string, &
          gcoords%root%baPC, gcoords%root%taC, gcoords%root%C%taC

    Else
      Call  branchedcoords_DumpNode(gcoords%root%C,angle1,angle2,num)
      natoms  = molecules_getnatoms(sorbtype)   

      Write(unitno,'(a)')  temp_string
      Write(unitno,'(2f13.5,i4)') gcoords%root%baPC,gcoords%root%taC, 0
      Write(unitno,'(2f13.5,i4)') ( angle1(i), angle2(i), 0, i=2,natoms-1 )
      !** Note  this zero is the arbitrarily assigned value for the 
      !** obsolete num() variable

    End If

  End Subroutine branchedcoords_dumprestartinfo


  !-----------------------------------------------------------
  ! Assign the angles from the restartfile to the nodes
  !-----------------------------------------------------------
  Recursive Subroutine branchedcoords_DumpNode(ptr,angle1,angle2,num)
    Type(NodeType), Pointer        :: ptr
    Real(kind=RDbl), Dimension(:), Intent(out) :: angle1 
    Real(kind=RDbl), Dimension(:), Intent(out) :: angle2 
    Integer, Dimension(:), Intent(out) :: num    

    Integer                        :: i 

    i = ptr%atom_num  
    angle1(i) = ptr%baPC  
    angle2(i) = ptr%taC 
    !SDEBUG
    !! OBSOLETE
    num(i)  = ptr%angle_index
    !SDEBUG

    If( Associated(ptr%C)) Then
      Call branchedcoords_DumpNode (ptr%C, angle1, angle2, num)
    End If

  End Subroutine branchedcoords_DumpNode

  !-----------------------------------------------------------
  ! Write out the node information to the restart file
  !-----------------------------------------------------------
  Recursive Subroutine branchedcoords_writenodeinfo(ptr, unitno)
    Type(NodeType), Pointer        :: ptr   
    Integer, Intent(in)            :: unitno

    If (Associated(ptr%C)) Then 
      Write(unitno,'(2f13.5)') ptr%baPC, ptr%taC
      ptr = ptr%C
    End If

  End Subroutine branchedcoords_writenodeinfo

End Module branchedcoords


