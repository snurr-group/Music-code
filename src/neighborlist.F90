!------------------------------------------------------------------------------
! This module handles the neighborlists for a given sorbate-sorbate pair
! The call is: 
!    getneighbors(nlists,sorb1,molec1,atom1,sorb2,molec2,list,size,error) 
! for the Atom_Atom neighborlist.  Remove atom1 & molec2 for the Molec_Molec
! neighborlist.  Note that, for a given sorb1,sorb2 pair, BOTH neighborlist
! types can be defined and used.  However, for a given pair, only a maximum
! of one type of neighbor list can currently be used.
!
! Written by Louis Clark (summer 1999)
!------------------------------------------------------------------------------

Module neighborlist

  Use defaults, Only: RDbl, strLen, MAX_SORBS, MAX_MOLECULES
  Use config, Only: AtMolCoords, config_getnmoles, config_getrp, &
      config_getMolecCOM
  Use vector, Only: VecType, vector_getdistsq, Assignment(=), Operator(+), &
      Operator(-), Operator(*), Operator(/)
  Use molecules, Only: molecules_getnatoms

  Implicit None
  Save

  Private

  Interface neighborlist_getneighbors
    Module Procedure neighborlist_getatoms
    Module Procedure neighborlist_getmolecules
  End Interface

  Interface neighborlist_display
    Module Procedure neighborlist_display_aa
    Module Procedure neighborlist_display_mm
  End Interface

  !----------------------------------------------------------------------------
  ! How to use the individual neighbor lists:
  ! NOTE1: There is a separate neighbor list for each sorbate-sorbate pair,
  !        it is accessed using nlist(sorb1,sorb2)
  ! NOTE2: There are two types of neighbor lists, Atom-Atom or Molec-Molec,
  ! - start(molec2,atom1,molec1) gives the starting index in %list() for atom2
  ! - end(molec2,atom1,molec1) gives the ending index in %list() for atom2
  ! - start(molec2,molec1) gives the starting index in %list() for molec2
  ! - end(molec2,molec1) gives the ending index in %list() for molec2
  ! - size is the number of atom2's stored in %list()
  ! - rad is the cutoff radius
  ! - rad2 is the cutoff radius squared
  ! - skin is the skin thickness
  ! - type stores the neighbor list type, i.e. 'Atom-Atom' or 'Molec-Molec'
  !----------------------------------------------------------------------------

  Type Last_Coords
    Type(VecType), Dimension(:,:), Pointer     :: r0
  End Type Last_Coords

  Type Last_COM_Coords
    Type(VecType), Dimension(:), Pointer       :: com0
  End Type Last_COM_Coords

  Type Atom_Atom_Neighbor_List
    Character(len=strLen)                               :: type
    Integer                                    :: sorb1,sorb2,size
    Real                                       :: rad,rad2,skin
    Integer, Pointer, Dimension(:,:,:)         :: start,end
    Type(Last_Coords), Dimension(2)            :: lastr          !(sorb1,sorb2)
    Integer, Pointer, Dimension(:)             :: list
  End Type

  Type Molec_Molec_Neighbor_List
    Character(len=strLen)                               :: type
    Integer                                    :: sorb1,sorb2,size
    Real                                       :: rad,rad2,skin
    Integer, Pointer, Dimension(:)             :: start,end
    Type(Last_COM_Coords), Dimension(2)        :: lastcom        !(sorb1,sorb2)
    Integer, Pointer, Dimension(:)             :: list
  End Type

  Type Neighbor_Lists
    Type(Atom_Atom_Neighbor_List), Pointer     :: aa
    Type(Molec_Molec_Neighbor_List), Pointer   :: mm
  End Type Neighbor_Lists
  
  Type(Neighbor_Lists), Dimension(MAX_SORBS,MAX_SORBS) :: nlists


  Contains

  !----------------------------------------------------------------------------
  ! Get the molec2 atoms that are in the neighborlist for sorb1,molec1,atom1
  !----------------------------------------------------------------------------
  Subroutine neighborlist_getatoms(sorbates,nlists,sorb1,molec1,atom1,&
      sorb2,molec2,list,size,error)

    Implicit None

    Type(AtMolCoords), Dimension(:), Intent(in)    :: sorbates
    Type(Neighbor_Lists), Intent(INOUT), Dimension(:,:) :: nlists
    Integer, Intent(IN)                   :: sorb1,molec1,atom1,sorb2,molec2
    Integer, Intent(OUT), Dimension(:)    :: list
    Integer, Intent(OUT)                  :: size
    Logical, Intent(OUT)                  :: error

    Integer                               :: list_start,list_end

    !--------------------------------------------------------------------------
    ! Make sure the list is initialized then get the list
    !--------------------------------------------------------------------------
    error = Associated(nlists(sorb1,sorb2)%aa)
    If (error) Return
    
    Call neighborlist_update_aa(sorbates,nlists,sorb1,sorb2)

    list_start = nlists(sorb1,sorb2)%aa%start(molec2,atom1,molec1)
    list_end = nlists(sorb1,sorb2)%aa%end(molec2,atom1,molec1)
    size = list_end - list_start + 1

    list = nlists(sorb1,sorb2)%aa%list(list_start:list_end)

  End Subroutine neighborlist_getatoms


  !----------------------------------------------------------------------------
  ! Get the sorb2 molecules that are in the neighborlist for sorb1,molec1
  !----------------------------------------------------------------------------
  Subroutine neighborlist_getmolecules(sorbates,nlists,sorb1,molec1,sorb2,  &
      list,size,error)

    Implicit None

    Type(AtMolCoords), Dimension(:), Intent(in)         :: sorbates
    Type(Neighbor_Lists), Intent(INOUT), Dimension(:,:) :: nlists
    Integer, Intent(IN)                   :: sorb1,molec1,sorb2
    Integer, Intent(OUT), Dimension(:)    :: list
    Integer, Intent(OUT)                  :: size
    Logical, Intent(OUT)                  :: error

    Integer                               :: list_start,list_end

    !--------------------------------------------------------------------------
    ! Make sure the list is initialized then get the list
    !--------------------------------------------------------------------------
    error = Associated(nlists(sorb1,sorb2)%mm)
    If (error) Return
    
    Call neighborlist_update_mm(sorbates,nlists,sorb1,sorb2)

    list_start = nlists(sorb1,sorb2)%mm%start(molec1)
    list_end = nlists(sorb1,sorb2)%mm%end(molec1)
    size = list_end - list_start + 1
    list = nlists(sorb1,sorb2)%mm%list(list_start:list_end)

  End Subroutine neighborlist_getmolecules


  !----------------------------------------------------------------------------
  ! Update the neighborlists as needed
  !----------------------------------------------------------------------------
  Subroutine neighborlist_update(sorbates, nlists,sorb1,sorb2)

    Implicit None

    Type(AtMolCoords), Dimension(:), Intent(in)         :: sorbates
    Type(Neighbor_lists), Dimension(:,:), Intent(INOUT) :: nlists
    Integer, Intent(IN)                                 :: sorb1,sorb2

    Integer                            :: flag

    !--------------------------------------------------------------------------

    flag = 0
    If (Associated(nlists(sorb1,sorb2)%aa)) Then
      Call neighborlist_update_aa(sorbates,nlists,sorb1,sorb2)
      flag = 1
    End If 

    If (Associated(nlists(sorb1,sorb2)%mm)) Then
      Call neighborlist_update_mm(sorbates,nlists,sorb1,sorb2)
      flag = 1
    End If 
    
    If (flag == 0) Then
      Write(0,'(2a,i4,a,2i4)') __FILE__,": ",__LINE__, &
           'Neither neighborlist type is associated ',sorb1,sorb2
      Stop
    Endif

  End Subroutine neighborlist_update


  !----------------------------------------------------------------------------
  ! Update the Atom-Atom neighborlist as needed
  !----------------------------------------------------------------------------
  Subroutine neighborlist_update_aa(sorbates,nlists,sorb1,sorb2)

    Implicit None

    Type(AtMolCoords), Dimension(:), Intent(in)  :: sorbates
    Type(Neighbor_lists), Dimension(:,:), Intent(INOUT) :: nlists
    Integer, Intent(IN)                :: sorb1,sorb2
                                    
    Logical                            :: update
    Integer                            :: i,j,k,error,size
    Integer                            :: molec1,molec2,atom1,atom2
    Integer                            :: nmolec1,natom1,nmolec2,natom2
    Integer                            :: m1_low,m1_high,m2_low,m2_high
    Real(kind = RDbl)                  :: maxdisp2,maxdisp,dist2
    Character(len=strLen)                       :: list_type
    Type(VecType)                      :: nowvec,lastvec,atvec1,atvec2

    !--------------------------------------------------------------------------

    nmolec1 = config_getnmoles(sorbates,sorb1)
    nmolec2 = config_getnmoles(sorbates,sorb2)
    natom1 = molecules_getnatoms(sorb1)
    natom2 = molecules_getnatoms(sorb2)
    update = .FALSE.
    
    !--------------------------------------------------------------------------
    ! Calculate the maximum displacement for the first sorbate
    !--------------------------------------------------------------------------
    maxdisp2 = 0.0_RDbl
    Do j = 1,nmolec1
      Do k = 1,natom1
        lastvec = nlists(sorb1,sorb2)%aa%lastr(1)%r0(atom1,molec1)
        nowvec = config_getrp(sorbates,(/sorb1,molec1,atom1/))
        nlists(sorb1,sorb2)%aa%lastr(1)%r0(atom1,molec1) = nowvec
        maxdisp2 = Max(vector_getdistsq(nowvec,lastvec),maxdisp2)
      End Do
    End Do

    !--------------------------------------------------------------------------
    ! Calculate the maximum displacement for the second sorbate
    !--------------------------------------------------------------------------
    Do j = 1,nmolec2
      Do k = 1,natom2
        lastvec = nlists(sorb1,sorb2)%aa%lastr(2)%r0(atom2,molec2)
        nowvec = config_getrp(sorbates,(/sorb2,molec2,atom2/))
        nlists(sorb1,sorb2)%aa%lastr(2)%r0(atom2,molec2) = nowvec
        maxdisp2 = Max(vector_getdistsq(nowvec,lastvec),maxdisp2)
      End Do
    End Do

    !--------------------------------------------------------------------------
    ! Update this neighborlist if the displacement is greater than 1/2*skin
    ! NOTE: the factor of 1/2 accounts for the possibility that two maximum
    !       displacement molecules are moving directly towards each other
    !--------------------------------------------------------------------------

    maxdisp = Sqrt(maxdisp2)
    update = (maxdisp > 0.5_RDbl*nlists(sorb1,sorb2)%aa%skin) 

    If (update) Then
      nlists(sorb1,sorb2)%aa%list = 0
      nlists(sorb1,sorb2)%aa%start = 0
      nlists(sorb1,sorb2)%aa%end = 0

      If (sorb1 == sorb2) Then
        m1_low = 1
        m1_high = nmolec1 - 1
      Else
        m1_low = 1
        m1_high = nmolec1
        m2_low = 1
        m2_high = nmolec2
      Endif

      size = 0
      Do molec1 = m1_low,m1_high
        If (sorb1 == sorb2) Then
          m2_low = 1
          m2_high = nmolec1 - 1
        End If

        Do molec2 = m2_low,m2_high
          Do atom1 = 1,natom1
      
            atvec1 = config_getrp(sorbates,(/sorb1,molec1,atom1/))
            nlists(sorb1,sorb2)%aa%start(molec2,atom1,molec1) = size + 1
      
            Do atom2 = 1,natom2
              atvec2 = config_getrp(sorbates,(/sorb2,molec2,atom2/))
              dist2 = vector_getdistsq(atvec1,atvec2)
              If (dist2 <= nlists(sorb1,sorb2)%aa%rad2) Then
                size = size + 1
                nlists(sorb1,sorb2)%aa%list(size) = atom2
              Endif
            End Do     
      
            nlists(sorb1,sorb2)%aa%end(molec2,atom1,molec1) = size
      
          End Do
        End Do
      End Do

      nlists(sorb1,sorb2)%aa%size = size

    Endif
          
  End Subroutine neighborlist_update_aa


  !----------------------------------------------------------------------------
  ! Update the Molec_Molec neighborlist as needed
  !----------------------------------------------------------------------------
  Subroutine neighborlist_update_mm(sorbates,nlists,sorb1,sorb2)

    Implicit None

    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Type(Neighbor_lists), Dimension(:,:), Intent(INOUT) :: nlists
    Integer, Intent(IN)                :: sorb1,sorb2
                                    
    Logical                            :: update
    Integer                            :: i,j,k,error,size
    Integer                            :: molec1,molec2,nmolec1,nmolec2
    Integer                            :: m1_low,m1_high,m2_low,m2_high
    Real(kind = RDbl)                  :: maxdisp2,maxdisp,dist2
    Type(VecType)                      :: nowvec,lastvec,mvec1,mvec2

    !--------------------------------------------------------------------------

    nmolec1 = config_getnmoles(sorbates,sorb1)
    nmolec2 = config_getnmoles(sorbates,sorb2)
    update = .FALSE.
    
    !--------------------------------------------------------------------------
    ! Calculate the maximum displacement for the first sorbate
    !--------------------------------------------------------------------------
    maxdisp2 = 0.0_RDbl
    Do i = 1,nmolec1
      lastvec = nlists(sorb1,sorb2)%mm%lastcom(1)%com0(molec1)
      nowvec = config_getMolecCOM(sorbates,sorb1,molec1)
      nlists(sorb1,sorb2)%mm%lastcom(1)%com0(molec1) = nowvec
      maxdisp2 = Max(vector_getdistsq(nowvec,lastvec),maxdisp2)
    End Do

    !--------------------------------------------------------------------------
    ! Calculate the maximum displacement for the second sorbate
    !--------------------------------------------------------------------------
    Do i = 1,nmolec2
      lastvec = nlists(sorb1,sorb2)%mm%lastcom(2)%com0(molec2)
      nowvec = config_getMolecCOM(sorbates,sorb2,molec2)
      nlists(sorb1,sorb2)%mm%lastcom(2)%com0(molec2) = nowvec
      maxdisp2 = Max(vector_getdistsq(nowvec,lastvec),maxdisp2)
    End Do

    !--------------------------------------------------------------------------
    ! Update this neighborlist if the displacement is greater than 1/2*skin
    ! NOTE: the factor of 1/2 accounts for the possibility that two maximum
    !       displacement molecules are moving directly towards each other
    !--------------------------------------------------------------------------

    maxdisp = Sqrt(maxdisp2)
    update = (maxdisp > 0.5_RDbl*nlists(sorb1,sorb2)%mm%skin) 

    If (update) Then
      nlists(sorb1,sorb2)%mm%list = 0
      nlists(sorb1,sorb2)%mm%start = 0
      nlists(sorb1,sorb2)%mm%end = 0

      If (sorb1 == sorb2) Then
        m1_low = 1
        m1_high = nmolec1 - 1
      Else
        m1_low = 1
        m1_high = nmolec1
        m2_low = 1
        m2_high = nmolec2
      Endif

      size = 0
      Do molec1 = 1,nmolec1
        mvec1 = nlists(sorb1,sorb2)%mm%lastcom(1)%com0(molec1)
        nlists(sorb1,sorb2)%mm%start(molec1) = size + 1

        If (sorb1 == sorb2) Then
          m2_low = 1
          m2_high = nmolec1 - 1
        End If

        Do molec2 = 1,nmolec2
          mvec2 = nlists(sorb1,sorb2)%mm%lastcom(2)%com0(molec2)
          dist2 = vector_getdistsq(mvec1,mvec2)
          If (dist2 <= nlists(sorb1,sorb2)%mm%rad2) Then
            size = size + 1
            nlists(sorb1,sorb2)%mm%list(size) = molec2
          End If
          nlists(sorb1,sorb2)%mm%end(molec2) = size
        End Do
      End Do
      nlists(sorb1,sorb2)%mm%size = size
    End If 
          
  End Subroutine neighborlist_update_mm


  !----------------------------------------------------------------------------
  ! Initialize the neighborlist
  !----------------------------------------------------------------------------
  Subroutine neighborlist_init(nlists,sorb1,sorb2,list_type,rad,skin)

    Implicit None

    Integer, Intent(IN)                 :: sorb1,sorb2
    Character(len=strLen), Intent(IN)   :: list_type
    Real(kind = RDbl)                   :: rad,skin
    Type(Neighbor_Lists), Intent(INOUT), Dimension(:,:) :: nlists
 
    Integer                             :: error,i,j
    Integer                             :: natom1,natom2
    Logical                             :: flagstart = .TRUE.

    !--------------------------------------------------------------------------
    ! Initialize the pointers
    !--------------------------------------------------------------------------
    If (flagstart) Then
      flagstart = .FALSE.
      
      Do i = 1,MAX_SORBS
        Do j = 1,MAX_SORBS
          If (i <= j) Then
            nullify(nlists(i,j)%aa)
            nullify(nlists(i,j)%mm)
          Else
            nlists(i,j)%aa => nlists(j,i)%aa
            nlists(i,j)%mm => nlists(j,i)%mm
          Endif          
        EndDo
      EndDo
    Endif

    !--------------------------------------------------------------------------
    ! Select the correct neighborlist type and allocate it
    ! Put the passed parameters into the neighborlist
    !--------------------------------------------------------------------------
    natom2 = molecules_getnatoms(sorb2)
    natom1 = molecules_getnatoms(sorb1)

    Select Case (list_type)

    Case ('Atom-Atom') 
      Allocate(nlists(sorb1,sorb2)%aa, STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a,2i4)') __FILE__,": ",__LINE__, &
             'unable to allocate nlist%aa ',sorb1,sorb2
        Stop
      Endif
      nlists(sorb2,sorb1)%aa => nlists(sorb1,sorb2)%aa

      nlists(sorb1,sorb2)%aa%type = list_type
      nlists(sorb1,sorb2)%aa%rad = rad
      nlists(sorb1,sorb2)%aa%rad2 = rad**2
      nlists(sorb1,sorb2)%aa%skin = skin
      nlists(sorb1,sorb2)%aa%sorb1 = sorb1
      nlists(sorb1,sorb2)%aa%sorb2 = sorb2

      Allocate(nlists(sorb1,sorb2)%aa%lastr(1)%r0(natom1,MAX_MOLECULES),&
          STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%aa%lastr%r0'
        Stop
      Endif
      Allocate(nlists(sorb1,sorb2)%aa%lastr(2)%r0(natom2,MAX_MOLECULES),&
          STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%aa%lastr%r0'
        Stop
      Endif
      Allocate(nlists(sorb1,sorb2)%aa%start &
          (MAX_MOLECULES,natom1,MAX_MOLECULES),STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%aa%start'
        Stop
      Endif
      Allocate(nlists(sorb1,sorb2)%aa%end &
          (MAX_MOLECULES,natom2,MAX_MOLECULES),STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%aa%end'
        Stop
      Endif
      Allocate(nlists(sorb1,sorb2)%aa%list(MAX_MOLECULES*natom2), STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%aa%list'
        Stop
      Endif

    Case ('Molec-Molec')
      Allocate(nlists(sorb1,sorb2)%mm, STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a,2i4)') __FILE__,": ",__LINE__, &
             'unable to allocate nlist%mm ',sorb1,sorb2
        Stop
      Endif
      nlists(sorb2,sorb1)%mm => nlists(sorb1,sorb2)%mm

      nlists(sorb1,sorb2)%mm%type = list_type
      nlists(sorb1,sorb2)%mm%rad = rad
      nlists(sorb1,sorb2)%mm%rad2 = rad**2
      nlists(sorb1,sorb2)%mm%skin = skin
      nlists(sorb1,sorb2)%mm%sorb1 = sorb1
      nlists(sorb1,sorb2)%mm%sorb2 = sorb2

      Allocate(nlists(sorb1,sorb2)%mm%lastcom(1)%com0(MAX_MOLECULES),&
          STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%mm%lastcom%com0'
        Stop
      Endif
      Allocate(nlists(sorb1,sorb2)%mm%lastcom(2)%com0(MAX_MOLECULES),&
          STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%mm%lastcom%com0'
        Stop
      Endif
      Allocate(nlists(sorb1,sorb2)%mm%start(MAX_MOLECULES),STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%mm%start'
        Stop
      Endif
      Allocate(nlists(sorb1,sorb2)%mm%end(MAX_MOLECULES),STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%mm%end'
        Stop
      Endif
      Allocate(nlists(sorb1,sorb2)%mm%list(MAX_MOLECULES), STAT=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'unable to allocate nlists()%mm%list'
        Stop
      Endif
      
    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
           'Could not recognize neighborlist type ',list_type
      Stop      
    End Select

  End Subroutine neighborlist_init


  !----------------------------------------------------------------------------
  ! Displays a given sorb1,sorb2 Atom-Atom neighborlist
  !----------------------------------------------------------------------------
  Subroutine neighborlist_display_aa(sorbates,nblist)

    Implicit None

    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Type(Atom_Atom_Neighbor_List), Intent(IN)   :: nblist

    Integer                         :: i,sorb1,sorb2
    Integer                         :: molec1,molec2,atom1,atom2
    Integer                         :: nmolec1,nmolec2,natom1,natom2
    Integer                         :: list_start,list_end,size

    !--------------------------------------------------------------------------
    ! Dump the list to screen
    !--------------------------------------------------------------------------

    nmolec1 = config_getnmoles(sorbates,nblist%sorb1)
    natom1 = molecules_getnatoms(nblist%sorb1)
    nmolec2 = config_getnmoles(sorbates,nblist%sorb2)

    Write(*,*)
    Do molec1 = 1,nmolec1
      Write(*,'(1x,a,i4)') 'Sorbate1: Molecule1 ',molec1
      Do atom1 = 1,natom1
        Write(*,'(4x,a,i4)') 'Sorbate1: Molecule1: Atom1 ',atom1
        Do molec2 = 1,nmolec2
          Write(*,'(7x,a,i4)') 'Sorbate1: Molecule1: Atom1: Molecule2 ',&
              molec2
          list_start = nblist%start(molec2,atom1,molec1)
          list_end = nblist%end(molec2,atom1,molec1)
          size = list_end - list_start + 1
          Write(*,'(10x,a,40i4)') 'Atom2 neighbors: ', &
               (nblist%list((list_start+i)),i=1,size)
        End Do
      End Do
    End Do

  End Subroutine neighborlist_display_aa


  !----------------------------------------------------------------------------
  ! Displays a given sorb1,sorb2 Molec-Molec neighborlist
  !----------------------------------------------------------------------------
  Subroutine neighborlist_display_mm(sorbates,nblist)

    Implicit None

    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Type(Molec_Molec_Neighbor_List), Intent(IN) :: nblist

    Integer                         :: i,sorb1,sorb2
    Integer                         :: molec1,molec2
    Integer                         :: nmolec1,nmolec2
    Integer                         :: list_start,list_end,size

    !--------------------------------------------------------------------------
    ! Dump the list to screen
    !--------------------------------------------------------------------------

    nmolec1 = config_getnmoles(sorbates,nblist%sorb1)
    nmolec2 = config_getnmoles(sorbates,nblist%sorb2)

    Write(*,*)
    Do molec1 = 1,nmolec1
      Write(*,'(1x,a,i4)') 'Sorbate1: Molecule1 ',molec1
      Do molec2 = 1,nmolec2
        list_start = nblist%start(molec1)
        list_end = nblist%end(molec1)
        size = list_end - list_start + 1
        Write(*,'(4x,a,i4)') 'Molecule2 neighbors: ', &
            (nblist%list(list_start+i),i=1,size)
      End Do
    End Do

  End Subroutine neighborlist_display_mm
  
End Module neighborlist
