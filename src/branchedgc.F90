!-----------------------------------------------------------
! This will most probably contain the data structure
! defined by Jeff Endelman for non-ring branched molecules
!----------------------------------------------------------
Module branchedgc

  Use defaults, Only: RDbl, zero, pi, twopi, one, degTorad
  Use file, Only:
  Use random, Only: rranf
  Use utils, Only: realswap
  Use vector, Only: VecType, mag, vector_display, vector_getnormsq, &
      vector_getcomp, vector_iscollinear, vector_crossprod, &
      Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnorm
  Use matrix, Only: MatrixType, Assignment(=), Operator(*)
  Use bmap, Only: bmap_getbiasindex, bmap_getncubelets, bmap_getbiaswt, &
      bmap_getbiaspt, bmap_getcellwt
  Use simcell, Only: simcell_getcellorigin, simcell_maptouc
  Use mcparams, Only:
  Use molecules, Only: molecules_getmolecule
  Use connects, Only: AtomList, nconnected, connects_nth_connection
  Use molecule, Only: molecule_getnatoms, molecule_getnthstereo, &
      molecule_getnthatom, MolecularParams
  Use branchedmcparams, Only: BranchedMolecNVTParams, BranchedMolecGCMCParams
  Use atom, Only: atom_getname

  Implicit None
  Save

  Private
  Public :: BranchedMolec, NodeType, branchedgc_initcoords, &
      branchedgc_toxyz, insert, delete, translate, rotate

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
     
  End Type NodeType

    
  Type BranchedMolec
     Type(VecType)             :: startpos
     Type(NodeType), Pointer   :: root
  End Type BranchedMolec
  
  Interface Assignment(=)
    Module Procedure branchedgc_copycoords
  End Interface
  
  Interface init
     Module Procedure branchedgc_initcoords
  End Interface
  
  Interface find_node
     Module Procedure branchedgc_find_node
  End Interface

  Interface delete
    Module Procedure branchedgc_delete
  End Interface

  Interface insert
    Module Procedure branchedgc_insert
  End Interface

  Interface translate
    Module Procedure branchedgc_randomtranslate
  End Interface

  Interface rotate
    Module Procedure branchedgc_randomrotate
  End Interface

Contains
  
  !------------------------------------------------------------------------
  ! Initializes the generalized coordinates as a tree 
  !------------------------------------------------------------------------
  
  Subroutine branchedgc_initcoords(gcoords, sorbtype)
    
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
    Call branchedgc_build_tree(molecule, gcoords%root)

    ! Build the proper stereochemistry into the tree
    Call branchedgc_set_stereo(gcoords%root)

    ! Set bond angles to equil values
    Call branchedgc_equil_ba(molecule, gcoords%root)

    ! Make torsion angles consistent
    Call branchedgc_update_ta(gcoords%root)

  End Subroutine branchedgc_initcoords

  !-----------------------------------------------------------
  ! Inserts first atom of the molecule according to a
  ! bias map and randomly chooses angular degrees of freedom.
  ! It also returns the bias of the insertion
  ! in "biasfactor"
  !-----------------------------------------------------------
  Subroutine branchedgc_insert(gcoords, xyzcoords, gcparams, biasfactor)
    Type(BranchedMolec), Pointer                   :: gcoords
    Type(VecType), Dimension(:), Intent(inout)  :: xyzcoords
    Type(BranchedMolecGCMCParams), Intent(in) :: gcparams
    Real(kind=RDbl), Intent(out)   :: biasfactor

    Integer              :: natoms, i, index, nx, ny, nz
    Integer              :: icellx,icelly,icellz, ncubelets
    Real(kind=RDbl)      :: costheta
    Type(VecType)        :: startpos

    !** Get the position of the first atom using the bias map
    !** 1.) Pick a point in the unit cell according to the biasmap
    index = bmap_getbiasindex(gcparams%bmap)
    ncubelets  = bmap_getncubelets(gcparams%bmap)
    biasfactor = 1.0_RDBl/bmap_getbiaswt(gcparams%bmap,index)/ncubelets
    startpos   = bmap_getbiaspt(gcparams%bmap, index)

    !** 2.) Select one of the unit cells in the simulation cell
    icellx = Int(rranf()*gcparams%simcell%nx)
    icelly = Int(rranf()*gcparams%simcell%ny)
    icellz = Int(rranf()*gcparams%simcell%nz)

    !** 3.) Translate the pt. to the appropriate unit cell
    startpos = startpos + &
        simcell_getcellorigin(gcparams%simcell, icellx, icelly, icellz)
    gcoords%startpos = startpos
    
    !** If no. of atoms is 1 we are done
    If (Size(xyzcoords, 1) == 1) Then
      xyzcoords(1) = startpos
      Return
    End If

    !** Generate the three rotational dof at random
    ! Sample theta from a cosine distribution.
    costheta = 1 - rranf() * 2.0_RDbl
    gcoords%root%baPC = Acos(costheta)
    gcoords%root%taC = pi - rranf()*twopi
    gcoords%root%C%taC = pi - rranf()*twopi

    !** Update ta
    Call branchedgc_update_ta(gcoords%root)

    !** Generate the xyz coordinates
    !  How should placing be handled???
    Call branchedgc_place(gcoords%root)
    Call branchedgc_toxyz(gcoords, xyzcoords)
    Call branchedgc_unplace(gcoords%root)


  End Subroutine branchedgc_insert

  !-----------------------------------------------
  !  Makes atoms placed
  !-----------------------------------------------
  
  Recursive Subroutine branchedgc_place(ptr)

    Type(NodeType), Pointer         :: ptr

    If (Associated(ptr)) Then
      ptr%placed = 1
      Call branchedgc_place(ptr%C)
      Call branchedgc_place(ptr%L)
      Call branchedgc_place(ptr%R)
    End If
  End Subroutine branchedgc_place

  !------------------------------------------------
  !  Unplaces atoms
  !------------------------------------------------

  Recursive Subroutine branchedgc_unplace(ptr)
    
    Type(NodeType), Pointer    :: ptr

    If (Associated(ptr)) Then
      ptr%placed = 0
      Call branchedgc_unplace(ptr%C)
      Call branchedgc_unplace(ptr%L)
      Call branchedgc_unplace(ptr%R)
    End If
  End Subroutine branchedgc_unplace

  !-----------------------------------------------------------------
  ! This routine gets the bias-weight for deleting a molecule
  !-----------------------------------------------------------------
  Subroutine branchedgc_delete(gcoords, gcparams, biasfactor)
    Type(BranchedMolec), Pointer              :: gcoords
    Type(BranchedMolecGCMCParams), Intent(in) :: gcparams
    Real(kind=RDbl), Intent(out)   :: biasfactor    

    Integer         :: index, ncubelets
    Real(kind=RDbl) :: wt
    Real(kind=RDbl), Dimension(3)  :: ucstartpos
    
    !** Map the position to the unit cell. Also convert the 
    !** vector type to an array
    ucstartpos = simcell_maptouc(gcparams%simcell, gcoords%startpos)

    !** Get the weigth of cubelet where the COM lies
    wt = bmap_getcellwt(gcparams%bmap, ucstartpos)

    !** Calculate the biasfactor
    ncubelets  = bmap_getncubelets(gcparams%bmap)
    biasfactor = ncubelets*wt

    Return
  End Subroutine branchedgc_delete

  !-------------------------------------------------------------
  ! Translates molecule and generates
  ! xyzcoords.
  !-------------------------------------------------------------
  Subroutine branchedgc_randomtranslate(gcoords, xyzcoords, gcparams)
    Type(BranchedMolec), Pointer                   :: gcoords
    Type(VecType), Dimension(:), Intent(inout)     :: xyzcoords
    Type(BranchedMolecNVTParams), Intent(in)       :: gcparams

    Type(VecType)       :: deltadisp, startpos
    Real(kind=RDbl)     :: deltatrans, xdelta, ydelta, zdelta
    
    deltatrans = gcparams%deltatrans
    If (gcparams%constrained) Then
      !** Do a coordinate transformation to the constrained plane
      !** system where the normal is along the x-axis
      startpos = gcoords%startpos
      startpos = gcparams%constrtrans*startpos

      !** Get the displacement in the y-z plane
      ydelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      zdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      deltadisp = (/0.0_RDbl, ydelta, zdelta/)
      startpos = startpos + deltadisp

      !** Tranform back to the original coordinate system
      gcoords%startpos = gcparams%constrinvtrans*startpos
    Else
      xdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      ydelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      zdelta =  ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltatrans
      deltadisp = (/xdelta, ydelta, zdelta/)
      gcoords%startpos = gcoords%startpos + deltadisp
    End If
    
    !** Generate the new coordinates
    Call branchedgc_place(gcoords%root)
    Call branchedgc_toxyz(gcoords, xyzcoords)
    Call branchedgc_unplace(gcoords%root)
  End Subroutine branchedgc_randomtranslate

  !-------------------------------------------------------------
  ! Rotate the molecule by generating new "eulerian" angles and 
  ! generates xyzcoords.
  !-------------------------------------------------------------
  Subroutine branchedgc_randomrotate(gcoords, xyzcoords, gcparams)
    Type(BranchedMolec), Pointer                   :: gcoords
    Type(VecType), Dimension(:), Intent(inout)     :: xyzcoords
    Type(BranchedMolecNVTParams), Intent(in)       :: gcparams

    Integer             :: natoms
    Type(VecType)       :: deltadisp, startpos
    Real(kind=RDbl)     :: deltarot, oldphi, phi, oldpsi, psi
    Real(kind=RDbl)     :: oldcostheta, costheta

    !** Check the no. of atoms and return if no. of atoms is 1
    natoms = Size(xyzcoords, 1)
    If (natoms == 1) Return

    deltarot = gcparams%deltarot

    !** Sample phi and psi uniformly between -pi,pi but theta
    !** needs to be sampled from a cosine distribution
    !** See Allen and Tildesley, Pg. 132-133, "Molecular Liquids"
    phi = gcoords%root%taC + ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
    gcoords%root%taC = phi - Anint(phi/twopi)*twopi

    psi = gcoords%root%C%taC + ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
    gcoords%root%C%taC = psi - Anint(psi/twopi)*twopi

    costheta = Cos(gcoords%root%baPC) + &
        ((rranf()*2.0_RDbl) - 1.0_RDbl)*deltarot
    costheta = costheta - Anint(costheta/2.0_RDbl)*2.0_RDbl
    gcoords%root%baPC = Acos(costheta)

    ! Update Torsion angles
    Call branchedgc_update_ta(gcoords%root)

    !** Generate the new coordinates
    Call branchedgc_place(gcoords%root)
    Call branchedgc_toxyz(gcoords, xyzcoords)
    Call branchedgc_unplace(gcoords%root)
  End Subroutine branchedgc_randomrotate
  
  
  !----------------------------------------------------------------
  !  Rotates local coordinates onto global coordinate frame
  !  First call should pass identity matrix in global_rotation
  !  Right now this only works if vec_size = 3.  
  !  Terminates when reaches an atom that has not been placed.
  !----------------------------------------------------------------
  
  Recursive Subroutine branchedgc_coord(ptr, xyzcoords, global_rotation)
    
    Type(NodeType), Pointer                        :: ptr
    Type(VecType), Dimension(:), Intent(InOut)     :: xyzcoords
    Type(MatrixType), Intent(in)                   :: global_rotation
    
    Type(MatrixType)                               :: local_rotation
    Type(VecType)                                  :: bond
    Real(Kind=RDbl)                                :: theta, phi
    
    If (Associated(ptr)) Then
      If (ptr%placed == 1) Then
       If (Associated(ptr%C)) Then
          theta = pi - ptr%baPC
          phi = ptr%taC
          
          !Call matrix_m_eq_2darray(local_rotation, &
          !     Reshape((/cos(theta),sin(theta)*cos(phi),sin(theta)*sin(phi), &
          !     sin(theta),-cos(theta)*cos(phi),-cos(theta)*sin(phi), &
          !     zero,sin(phi),-cos(phi)/), (/3, 3/)))
          local_rotation = Reshape((/cos(theta),sin(theta)*cos(phi), &
              Sin(theta)*Sin(phi),sin(theta),-cos(theta)*cos(phi), &
              -Cos(theta)*Sin(phi),zero,sin(phi),-cos(phi)/), (/3, 3/))

          ! local_rotation = matrix_m_mult_m(global_rotation, local_rotation)
          local_rotation = global_rotation*local_rotation

          bond = (/ptr%C%bond_length,zero,zero/)
          ! Call vector_vectype_eq_rarray(bond,(/ptr%C%bond_length,zero,zero/))
          
          xyzcoords(ptr%C%atom_num) = xyzcoords(ptr%atom_num) + &
              local_rotation*bond
          !local_rotation*bond = matrix_m_mult_v(local_rotation, bond)
         
          Call branchedgc_coord(ptr%C, xyzcoords, local_rotation)
       End If
       
       If (Associated(ptr%L)) Then
          theta = pi - ptr%baPL
          phi = ptr%taL
          
          !Call matrix_m_eq_2darray(local_rotation, &
          !     Reshape((/cos(theta),sin(theta)*cos(phi),sin(theta)*sin(phi), &
          !     sin(theta),-cos(theta)*cos(phi),-cos(theta)*sin(phi), &
          !     zero,sin(phi),-cos(phi)/), (/3, 3/)))
          local_rotation = Reshape((/cos(theta),sin(theta)*cos(phi), &
              Sin(theta)*Sin(phi),sin(theta),-cos(theta)*cos(phi), &
              -Cos(theta)*Sin(phi),zero,sin(phi),-cos(phi)/), (/3, 3/))
          
          ! local_rotation = matrix_m_mult_m(global_rotation, local_rotation)
          local_rotation = global_rotation*local_rotation
          
          ! Call vector_vectype_eq_rarray(bond,(/ptr%L%bond_length,zero,zero/))
          bond = (/ptr%L%bond_length,zero,zero/)

          ! local_rotation*bond = matrix_m_mult_v(local_rotation, bond)
          xyzcoords(ptr%L%atom_num) = xyzcoords(ptr%atom_num) + &
              local_rotation*bond
                
          Call branchedgc_coord(ptr%L, xyzcoords, local_rotation)
       End If
       
       If (Associated(ptr%R)) Then
          theta = pi - ptr%baPR
          phi = ptr%taR
          
          !Call matrix_m_eq_2darray(local_rotation, &
          !     Reshape((/cos(theta),sin(theta)*cos(phi),sin(theta)*sin(phi), &
          !     sin(theta),-cos(theta)*cos(phi),-cos(theta)*sin(phi), &
          !     zero,sin(phi),-cos(phi)/), (/3, 3/)))
          local_rotation = Reshape((/cos(theta),sin(theta)*cos(phi), &
              Sin(theta)*Sin(phi),sin(theta),-cos(theta)*cos(phi), &
              -Cos(theta)*Sin(phi),zero,Sin(phi),-Cos(phi)/), (/3, 3/))
          
          ! local_rotation = matrix_m_mult_m(global_rotation, local_rotation)
          local_rotation = global_rotation*local_rotation

          ! Call vector_vectype_eq_rarray(bond,(/ptr%R%bond_length,zero,zero/))
          bond = (/ptr%R%bond_length,zero,zero/)

          xyzcoords(ptr%R%atom_num) = xyzcoords(ptr%atom_num) + local_rotation*bond
          ! local_rotation*bond = matrix_m_mult_v(local_rotation, bond)

          Call branchedgc_coord(ptr%R, xyzcoords, local_rotation)
        End If
      End If
    End If
  End Subroutine branchedgc_coord
  
  
  !-----------------------------------------------------------------------
  ! Generates xyz coordinates from generalized coordinates
  ! Only works for vec_size = 3
  !-----------------------------------------------------------------------
  
  Subroutine branchedgc_toxyz(gcoords, xyzcoords)
    
    Type(BranchedMolec), Pointer                    :: gcoords
    Type(VecType), Dimension(:), Intent(out)        :: xyzcoords
    
    Type(MatrixType)                                :: identity

    ! Atom1 must be the root of the tree
    xyzcoords(1) = gcoords%startpos
  
  
    !Call matrix_m_eq_2darray(identity, &
    !     Reshape((/one,zero,zero,zero,one,zero,zero,zero,one/),(/3,3/)))
    identity = Reshape((/one,zero,zero,zero,one,zero,zero,zero,one/),(/3,3/))


    Call branchedgc_coord(gcoords%root,xyzcoords,identity)  

  End Subroutine branchedgc_toxyz
  
  
  !------------------------------------------------------
  ! Returns a pointer (node) to the node with atom_num.
  ! Can only search below root.
  !------------------------------------------------------
  Recursive Subroutine branchedgc_find_node(root, atom_num, node, found)
    
    Type(NodeType), Pointer    :: root, node
    Integer, Intent(In)        :: atom_num
    Integer, Intent(Out)       :: found  ! 1 if successful, 0 if not
    
    found = 0
    
    If (Associated(root)) Then
       If (root%atom_num == atom_num) Then
          node => root
          found = 1
       Else
          Call branchedgc_find_node(root%C, atom_num, node, found)
          If (found /= 1) Then
             Call branchedgc_find_node(root%L, atom_num, node, found)
          End If
          If (found /= 1) Then
             Call branchedgc_find_node(root%R, atom_num, node, found)
          End If
       End If
    End If
    
  End Subroutine branchedgc_find_node
  
  !--------------------------------------------------------------------
  ! Gets passed a pointer to a node, zeroes angles, connects 
  ! children, and enters in their information
  !--------------------------------------------------------------------
  Recursive Subroutine branchedgc_build_tree(molecule, this_ptr)
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
       
       num_connect = nconnected(molecule%connect(this_ptr%atom_num))
       
       Do i = 1, num_connect
          
          ! Get the ith connection
          Call connects_nth_connection(molecule%connect(this_ptr%atom_num), &
               i, connection)
          
          ! Don't connect to previous atom
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
             Call branchedgc_build_tree(molecule, new_ptr)
          Else If (.NOT.(Associated(this_ptr%L))) Then
             this_ptr%L => new_ptr
             Call branchedgc_build_tree(molecule, new_ptr)
          Else
             this_ptr%R => new_ptr
             Call branchedgc_build_tree(molecule, new_ptr)
          End If
       End Do
    End If

  End Subroutine branchedgc_build_tree

  !---------------------------------------------------------------------------
  !  Switches left and right children and angles if necessary to effect the
  !  proper stereochemistry.  If there is no stereochemistry,
  !  orientation is chosen randomly (won't matter if not a stereocenter,
  !  otherwise leads to racemic mixture).
  !---------------------------------------------------------------------------

  Recursive Subroutine branchedgc_set_stereo(this_ptr)

    Type(NodeType), Pointer  :: this_ptr
    
    Type(NodeType), Pointer  :: temp_ptr
    Character                :: current, lowest
    Real(Kind=RDbl)          :: temp_angle
    
    Nullify(temp_ptr)
    
    If (Associated(this_ptr)) Then
      
      ! Can't have stereochemistry at first atom.  Should be terminal.
      If (Associated(this_ptr%P)) Then   
        
        ! Check whether valence is at least 3
        If (Associated(this_ptr%L)) Then
          
          ! If there is no stereochemistry, randomly make one descendant 
          ! left and the other right (thus it's randomly S or R)
          If (this_ptr%stereo == "N") Then
            If (rranf() < 0.5) Then
              ! Swap pointers
              temp_ptr => this_ptr%L
              this_ptr%L => this_ptr%R
              this_ptr%R => temp_ptr
              
              ! Swap angles
              Call realswap(this_ptr%baPL, this_ptr%baPR)
              Call realswap(this_ptr%baCL, this_ptr%baCR)
              Call realswap(this_ptr%taL, this_ptr%taR)
            End If
          Else
            
            If (Associated(this_ptr%R)) Then
              ! Valence is 4
              ! Find lowest priority atom
              lowest = branchedgc_find_min(this_ptr)
            Else
              ! Valence is 3. So nonexistent right child is lowest priority
              lowest = 'R'
            End If
            
            ! Determine current stereochemistry
            Select Case (lowest)
              ! All calls to branchedgc_determine_stereo must pass
              ! pointers in clockwise orientation with unnamed pointer
              ! in the back.  See NodeType definition.
              
            Case ('R')
              current = branchedgc_determine_stereo(this_ptr%P%atom_num, &
                  this_ptr%L%atom_num, this_ptr%C%atom_num)
              
            Case ('P')
              current = branchedgc_determine_stereo(this_ptr%C%atom_num, & 
                  this_ptr%L%atom_num, this_ptr%R%atom_num)
              
            Case ('C')
              current = branchedgc_determine_stereo(this_ptr%P%atom_num, &
                  this_ptr%R%atom_num, this_ptr%L%atom_num)
              
            Case ('L')
              current = branchedgc_determine_stereo(this_ptr%P%atom_num, &
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
              Call realswap(this_ptr%baPL, this_ptr%baPR)
              Call realswap(this_ptr%baCL, this_ptr%baCR)
              Call realswap(this_ptr%taL, this_ptr%taR)
            End If
            
          End If
        End If
      End If
      
      If (Associated(this_ptr%C)) Call branchedgc_set_stereo(this_ptr%C)
      If (Associated(this_ptr%L)) Call branchedgc_set_stereo(this_ptr%L)
      If (Associated(this_ptr%R)) Call branchedgc_set_stereo(this_ptr%R)
    End If
    
  End Subroutine branchedgc_set_stereo

  !--------------------------------------------------------------------------
  !  Returns a letter indicating which connection P,C,L,or R has
  !  the smallest atom_num.  this_ptr and this_ptr%P MUST be associated, 
  !  otherwise returns 'F'
  !--------------------------------------------------------------------------
  
  Character Function branchedgc_find_min(this_ptr)
    
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
    
    branchedgc_find_min = letter

  End Function branchedgc_find_min



  !-------------------------------------------------------------------------
  ! Determines the stereochemistry for a sequence of three numbers, which
  ! are correspond to atom numbers in clockwise fashion.  If these numbers
  ! are descending clockwise, 'R' is returned.  Otherwise 'S'.
  !----------------------------------------------------------------------
  Character Function branchedgc_determine_stereo(num1, num2, num3)

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
         branchedgc_determine_stereo = 'R'
      Else 
         branchedgc_determine_stereo = 'S'
      End If
         
    End Function branchedgc_determine_stereo

    !----------------------------------------------------------
    !  Initialize bond angles to their equilibrium values.
    !----------------------------------------------------------
    
    Recursive Subroutine branchedgc_equil_ba(molecule, ptr)
      
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
            ptr%baPC = degTorad* &
            molecule%bending%bbparams(atomP,atom,atomC)%hara%thetaeq
         End If
     
         If (Associated(ptr%P) .AND. Associated(ptr%L)) Then
            ptr%baPL = degTorad* &
            molecule%bending%bbparams(atomP,atom,atomL)%hara%thetaeq
         End If
  
         If (Associated(ptr%P) .AND. Associated(ptr%R)) Then
            ptr%baPR = degTorad* &
            molecule%bending%bbparams(atomP,atom,atomR)%hara%thetaeq
         End If
        
         If (Associated(ptr%C) .AND. Associated(ptr%L)) Then
            ptr%baCL = degTorad* &
            molecule%bending%bbparams(atomC,atom,atomL)%hara%thetaeq
         End If

         If (Associated(ptr%C) .AND. Associated(ptr%R)) Then
            ptr%baCR = degTorad* &
            molecule%bending%bbparams(atomC,atom,atomR)%hara%thetaeq
         End If

         If (Associated(ptr%L) .AND. Associated(ptr%R)) Then
            ptr%baLR = degTorad* &
            molecule%bending%bbparams(atomL,atom,atomR)%hara%thetaeq
         End If

         If (Associated(ptr%C)) Call branchedgc_equil_ba(molecule, ptr%C)
         If (Associated(ptr%R)) Call branchedgc_equil_ba(molecule, ptr%R)
         If (Associated(ptr%L)) Call branchedgc_equil_ba(molecule, ptr%L)

      End If

    End Subroutine branchedgc_equil_ba

    
    !-------------------------------------------------------------------
    !  Updates taL and taR to be geometrically consistent with baPC,
    !  baPL, baPR, baCL, and baCR.  baLR is ignored.
    !------------------------------------------------------------------

    Recursive Subroutine branchedgc_update_ta(ptr)
        
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

            Call branchedgc_update_ta(ptr%L)
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

            Call branchedgc_update_ta(ptr%R)
         End If            

         If (Associated(ptr%C)) Call branchedgc_update_ta(ptr%C)
      End If
      
    End Subroutine branchedgc_update_ta

    !-------------------------------------------------
    !  Displays generalized coordinates
    !-------------------------------------------------

    Subroutine branchedgc_display(gcoords, molecno, unitno)
      
      Type(BranchedMolec), Pointer      :: gcoords
      Integer, Intent(In)               :: molecno, unitno
      
      If (Associated(gcoords)) Then
        Write(unitno,*) "Generalized Coordinates for Molecule #", molecno
        Write(unitno,*) "First Atom Position: ", &
            Trim(vector_display(gcoords%startpos,"f13.5"))
        Write(unitno,*)
        Call branchedgc_display_loop(gcoords%root, unitno)
      End If

    End Subroutine branchedgc_display

    !-----------------------------------------------------
    ! Displays data at every node below the node passed
    ! in the initial call.
    !-----------------------------------------------------
      
    Recursive Subroutine branchedgc_display_loop(ptr, unitno)
      
      Type(NodeType), Pointer :: ptr
      Integer, Intent(in)     :: unitno
      
      If (Associated(ptr)) Then
 75      Format(A,T25,F8.2)
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
         Call branchedgc_display_loop(ptr%C, unitno)
         Call branchedgc_display_loop(ptr%L, unitno)
         Call branchedgc_display_loop(ptr%R, unitno)
       End If
      
     End Subroutine branchedgc_display_loop

     !-------------------------------------------------
     !  Copies the branchedcoords
     !-------------------------------------------------

     Subroutine branchedgc_copycoords(outcoords,incoords)
       Type(BranchedMolec), Intent(inout)  :: outcoords
       Type(BranchedMolec), Intent(in)     :: incoords

       outcoords%startpos = incoords%startpos

       If (Associated(incoords%root)) Then
         Call branchedgc_copynode(outcoords%root, incoords%root)
       End If
       
     End Subroutine branchedgc_copycoords
       
     !-----------------------------------------------
     !  Copies the tree structure
     !----------------------------------------------

     Recursive Subroutine branchedgc_copynode(outptr, inptr)

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
       
       If (Associated(inptr%C)) Call branchedgc_copynode(outptr%C,inptr%C)
       If (Associated(inptr%L)) Call branchedgc_copynode(outptr%L,inptr%L)
       If (Associated(inptr%R)) Call branchedgc_copynode(outptr%R,inptr%R)

     End Subroutine branchedgc_copynode

End Module branchedgc









