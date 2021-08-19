!---------------------------------------------------------------------------- 
! This module contains the data structure and code 
! for working with cavities in the simulation cell. It can be 
! used as a verlet neighbor list with some additional programming work
! For now this will be used for cavity+nrg bias rigid body insertions.
!    -----------------------------------------
!
! Detalis of cubelist :
! - the field cube_array contains the actual used cavities, it is a 1-D array
! - the field cubeindex maps the acutal 3-D positions into the cube_array
!    If cube_index(1,273, 15)=1176, That means the cube at x=1*dx, y=273*dy 
!    and z=15*dz is stored at cub_array(1176)
! More Details:
!  This Module is very very complicated, But a  simplified picture is here:
!  Assume 2 1-D arrays Cubes(c=1..ncubes) and Molecs(m=1:nmoles)
!  Assume single component mono atomic
! The occupancy can be shown by connections between these arrays
! Suppose c=10 is connected to m=(11, 15); that means cube 10 contains 
! molecules  11, and 15
! Similarly if m=11 is connected to c=10 and c=12 Then that means 
! molecule 11 occupies both cubes 10 and 12. 
! When we insert we make connections, when we delete we break connection
! when we translate/rotate we change connections 
!
! Complimetary 1-D arrays : For faster access of information there are 
! many complimentary arrays stored here.
! example : Cubes <-> Molecs
!           sorblist <-> csorblist
! If cubes(c) contains molecule=m, then molecs(m) will contains cube=c info
! Other points:
! - all sorbs are not used for cavity analysis. 
!   "sorblist" in indexed from 1:"number of sorbs used in cavity analysis"
!   "csorblist(some-sorb)" will give the index of "some-sorb" 
!   in the "sorblist" array [some-sorb is the molecules_type()]
! WARNING : watch the following defaults:
!  Integer, Parameter :: MAX_MOLECS_PER_CUBE=5
!  Integer, Parameter :: INITIAL_NUM_OF_MOLECS=50
!  Integer, Parameter :: MEMORY_STEP_SIZE=10! 
!----------------------------------------------------- 

Module cavitylist
  Use config, Only: AtMolCoords, config_getnmoles
  Use defaults, Only: strLen, lstrLen, RDbl, RSgl, MAX_MAPS, dashedline,&
      one, zero, d_ctrl_file, MAX_ATOMS
  Use file, Only: file_open
  Use molecules, Only: molecules_gettype, molecules_getnsorbs, &
      molecules_getnatoms
  Use mapheader, Only: mapheader_display, Header_Info, mapheader_cleanup, &
      mapheader_getheader
  Use random, Only: rranf
  Use simcell, Only: Simcell_Params, simcell_getell, simcell_getnuc
  Use utils, Only: split, getpath, allocErrDisplay, filesrchstr, int2str, &
      stripcmnt, findint
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)

  Implicit None
  Save

  Private
  Public :: Cavity_Params, Cube_Info, Sorb_Info, Molec_Info, &
      cavitylist_init, cavitylist_makeXYZ, &
      cavitylist_checkcavity, cavitylist_displayCubeInfo, &
      cavitylist_displaySorbInfo,       cavitylist_display, &
      cavitylist_checkandIncrMemory, cavitylist_updateMolecInfo, &
      cavitylist_getCubeIndices


  Type Cavity_Params

    Integer                            :: nsorbs
    Real(kind=RDbl)                    :: dx, dy, dz ! cube size
    Real(kind=RDbl)                    :: lx, ly, lz ! simcell edge Lengths
    Integer                            :: nx, ny, nz ! nx*dx = lx
    Integer                            :: nucx, nucy, nucz  
    Integer                            :: ncubes


    ! dimension= nsorbs(used for cavity analysis)
    Integer, Dimension(:), Pointer         :: sorblist
    Character(len=strLen), Dimension(:), Pointer  :: sorbnames
    Type (Sorb_Info), Dimension(:),Pointer :: sorbs 

    ! dimension = total number of sorbate types in the simulation
    Integer, Dimension(:), Pointer         :: csorblist
    Logical, Dimension(:), Pointer         :: is_cavity_sorb

    ! dimension = nx, ny, nz
    Integer, Dimension(:,:,:), Pointer     :: cube_index

    ! dimension = ncubes
    Type(Cube_Info), Dimension(:), Pointer :: cube_array 

    !** Variables related to the energy map
    Character(len=lstrLen)    :: pmapName
    Type(Header_Info)        :: header
    Real(kind=RDbl)          :: ucutoff

  End Type Cavity_Params

  !** Information stored in each cube
  Type Cube_Info
    Integer :: ix, iy, iz            !** indices in the full 3-D array
    Integer                         :: nsorbs
    Integer, Dimension(:), Pointer  :: sorblist
    Integer, Dimension(:), Pointer  :: nmolecs
    Integer, Dimension(:,:),Pointer :: moleclist
    Logical :: is_init=.False.
  End Type Cube_Info

  !** Information about each sorbate, this can be considered as very
  !** rough, lattice like information about the position of sorbates, these
  !** need not be initialized for unwanted sorbtypes
  Type Sorb_Info
    Integer                                  :: nmolecs
    Type (Molec_Info), Dimension(:), Pointer :: molecs
    Integer                                  :: currentsize
  End Type Sorb_Info

  !** Information about each molecule's position, in the cavities
  Type Molec_Info
    Integer                                  :: ncubes
    Integer, Dimension(:), Pointer           :: cube_list 
    !** cubelist will usually be assigned a dimension = 
    !**                                 molecules_getnatoms(sorbtype)
    ! ncubes = natoms, there will be duplicacies
  End Type Molec_Info

  Interface cavitylist_checkcavity
    Module Procedure cavitylist_checkcavity1
    Module Procedure cavitylist_checkcavity2
  End Interface




  Integer, Parameter :: MAX_MOLECS_PER_CUBE=5
  Integer, Parameter :: INITIAL_NUM_OF_MOLECS=200
  Integer, Parameter :: MEMORY_STEP_SIZE=10

  Character(len=strLen) :: default_cavitylist_tag=  "Cavity Bias Info"
  Integer            :: number_of_cavity_maps=0
  Type(Cavity_Params), Pointer ::  cavTMap
  Logical :: cavitylist_initialized=.False.
Contains
  !------------------------------------------------------------------------
  ! Initialize the 'cavitylist'.  If one is already initilized point to &
  ! the already initialized one
  ! Requires: cavitylist -- the empty cubes  object
  !           sorbates, scell -- the coordinates, simulation cell
  !           ctrl_filename -- the file with init info
  !           opt_tag -- tag marking the cavity bias section
  !           opt_found -- If no tag is found, then no cavity, simply return
  !------------------------------------------------------------------------
  Subroutine cavitylist_init(cubelist, sorbates, scell, ctrl_filename, &
opt_found,      opt_tag) 
    Type(Cavity_Params), Pointer              :: cubelist
    Type(AtMolCoords), Dimension(:), Intent(in)  :: sorbates
    Type( Simcell_Params), Intent(in)            :: scell
    

    Character(*),                  Intent(in) :: ctrl_filename
    Character(*), Optional, Intent(in)        :: opt_tag
    Logical     , Optional, Intent(out)        :: opt_found
    Real(kind=RDbl), Dimension(3)             :: ell

    Character(len=strLen) :: tag, srchstr
    Integer               :: unitno, error, lineno
    If (present(opt_found)) opt_found=.True.
    If (Present(opt_tag)) Then
      tag=opt_tag
    Else
      tag=default_cavitylist_tag
    End If

    !** Check whether CavTMap is already initialised
    If (number_of_cavity_maps==0) Then
      Allocate(cavTMap,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Nullify(cubelist)
      cubelist=>cavTMap
      ! and then continue rest of initialization
      number_of_cavity_maps=number_of_cavity_maps+1

    Else
      Nullify(cubelist)
      cubelist=>cavTMap
      Return

    Endif

    ! only one list containing cavityies should be initialized
    If (cavitylist_initialized) Then
      Write(*,*) "Only one list containing cavities should be initialized"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif
    cavitylist_initialized=.True.

    !** Open the ctrl_file if it is not opened and find the reqd section
    unitno = file_open(ctrl_filename)
    Rewind(unitno)
    lineno = filesrchstr(unitno, tag, srchstr)
    If (lineno == 0) Then
      If (present(opt_found)) then
        opt_found=.False.
        Return
      Else
        Write(0,'(1x, 4a)') 'Could not find the tag "', Trim(tag), &
            '" in the control file ', d_ctrl_file
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    End If
    If (.Not.Associated(cubelist)) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "CubeList should be associated before calling this &
          &routine"
      stop
    Endif

    !** Read the important fields from ctrlfile and initialize
    Call cavitylist_ReadCtrlfile(cubelist, unitno)

    Allocate(cubelist%cube_index(cubelist%nx,cubelist%ny,cubelist%nz), &
        STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__, "3D-Array of cubes")

    cubelist%cube_index=0

    !** using vector assignment overloading here
    ell=simcell_getell(scell)

    cubelist%lx=ell(1)
    cubelist%ly=ell(2)
    cubelist%lz=ell(3)
    cubelist%dx= cubelist%lx / cubelist%nx
    cubelist%dy= cubelist%ly / cubelist%ny
    cubelist%dz= cubelist%lz / cubelist%nz

    ! number of unit cells in each direction
    cubelist%nucx=simcell_getnuc(scell,'x')
    cubelist%nucy=simcell_getnuc(scell,'y')
    cubelist%nucz=simcell_getnuc(scell,'z')

    ! we need to check that nx , ny and nz are exactly divisible by nucx,y,z
    If ( (Mod(cubelist%nx,cubelist%nucx)/=0).Or. &
        (Mod(cubelist%ny,cubelist%nucy)/=0).Or.&
        (Mod(cubelist%nz,cubelist%nucz)/=0) ) Then
      Write(*,*) "The nx, ny, nz specified in ctrlfile should be exactly &
          & divisible by number of unitcells in thgose respective direction. &
          & i.e. : simcell section should be consistant with the cavitylist &
          &section"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    !** First Initialize the Cubes part
    ! This is where information about which moleucles are contained in each
    ! cube will be stored
    If ((cubelist%pmapName)=="NULL") Then 
      Call cavitylist_initArrUniform(cubelist)
    Else 
      Call cavitylist_initArrPmap(cubelist)
    Endif

    !** Then initalise the individual sorbate info part
    ! This is where information about which cubes are occupied by each 
    ! molecule will be stored
    Nullify(cubelist%sorbs)
    Call cavitylist_initSorbInfo(cubelist,sorbates)



  End Subroutine cavitylist_init


  !-------------------------------------------------
  ! Initializes information about, in  which cell each of the molecules 
  ! are lying With the current configuration, by looping thru whole sorbate
  !-------------------------------------------------
  Subroutine cavitylist_initSorbInfo(cubelist,sorbates)
    Type(Cavity_Params), Intent(inout)              :: cubelist
    Type(AtMolCoords), Dimension(:), Intent(in)  :: sorbates

    Integer :: sorb, i_sorb, m, error, nsorbs, nmoles, newsize
    Integer :: natoms 

    !** Make sure this was already not initialized
    If (Associated(cubelist%sorbs)) Then
      Write(*,*) " Nullify cubelist%sorbs before calling this routine"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    nsorbs=molecules_getnsorbs()
    Allocate(cubelist%sorbs(cubelist%nsorbs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Do i_sorb=1,nsorbs
      If (cavitylist_checkCavitySorb(cubelist,i_sorb)) Then
        sorb=cubelist%csorblist(i_sorb)
        cubelist%sorbs(sorb)%nmolecs=0

        nmoles=config_getnmoles(sorbates, i_sorb)

        !** How much memory do we want to give initially
        If (nmoles<INITIAL_NUM_OF_MOLECS) Then
          newsize=INITIAL_NUM_OF_MOLECS
        Else
          newsize=nmoles + INITIAL_NUM_OF_MOLECS
        Endif

        Allocate(cubelist%sorbs(sorb)%molecs(newsize), STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        cubelist%sorbs(sorb)%currentsize=newsize

        natoms=molecules_getnatoms(i_sorb)

        Do m=1,newsize
          ! cubelist tells which cube is ocupied by the atoms of the molecule
          ! in the bin molecs(m).
          cubelist%sorbs(sorb)%molecs(m)%ncubes=0
          Allocate(cubelist%sorbs(sorb)%molecs(m)%cube_list(natoms), &
              STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
          cubelist%sorbs(sorb)%molecs(m)%cube_list=0

        End Do

        Do m=1,nmoles
          ! 101 - implies new molecule
          Call cavitylist_updateMolecInfo(cubelist,sorbates(i_sorb), &
              i_sorb , m, 101)
        End Do

      Endif
    End Do

  End Subroutine cavitylist_initSorbInfo

  !-----------------------------------------------------------------
  ! Changes the cubelist structure based on the changes in a particular 
  ! molecule
  !      action= 101       - Molecules was inserted new
  !            = 111       - Molecules was transfomed from old to new position
  !            = 110       - Molecules has to be deleted
  ! NOTE : - Should be used only with GCMC insert, delete, and translate
  !        - WIth any other driver be extremely careful
  !-----------------------------------------------------------------
  Subroutine cavitylist_updateMolecInfo(cubelist, sorbate, sorb, molec, &
      action)
    Type(Cavity_Params), Intent(inout)              :: cubelist
    Type(AtMolCoords), Intent(in)                   :: sorbate
    Integer, Intent(in) :: sorb,molec, action

    Integer, Dimension(MAX_ATOMS), Save :: CurrentCubes, new_list, old_list
    Integer :: cArrInd, i, new_c, old_c, natoms, lastMolec, last_c
    Integer :: posn_in_array

    ! index of this sorb in the cubelist
    cArrInd=cubelist%csorblist(sorb)    

    If (cArrInd==0) Then
      Write(*,*) "probably wrong sorbname in cavitylist section ?"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    natoms=molecules_getnatoms(sorb)

!!$    If (startdebug) Then
!!$    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$    Write(*,*) "Updating molecule InfO", "sorb", sorb, "cArrInd" , cArrInd
!!$    Write(*,*) "natoms", natoms, "molec", molec, "action",action
!!$    Write(*,*) 
!!$    Endif

    Select Case (action)

    Case(101) 
      ! WE ARE INSERTING
      !SDEBUG
      If (molec>=cubelist%sorbs(cArrInd)%currentsize) Then
        Write(*,*) "size problems here, insering molecule: ", molec, &
            "current size only : ", cubelist%sorbs(cArrInd)%currentsize
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif


      If (molec/=(cubelist%sorbs(cArrInd)%nmolecs+1)) Then
        Write(*,*) "inserted molecules index and number of current &
            & molecules dont match"
        Write(*,*) molec, cubelist%sorbs(cArrInd)%nmolecs,&
            cubelist%sorbs(cArrInd)%currentsize, cArrInd
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      !SDEBUG

      !** get index of cube for each atom
      Call cavitylist_getCubeIndices(cubelist, sorbate, molec, &
          natoms,CurrentCubes(1:natoms))

      Do i=1,natoms
        new_c=CurrentCubes(i)
        If (new_c>0) Then
          Call cavitylist_addMolec(cubelist, cArrInd, molec, new_c )
        Else
          Write(*,*) "The new suggested point is lying in an untabulated &
              & area, check whether the pmaps/bmaps are consistent"
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Stop
        Endif
      End Do

      !** Update relevant info store in sorbinfo
      cubelist%sorbs(cArrInd)%nmolecs = cubelist%sorbs(cArrInd)%nmolecs + 1
      cubelist%sorbs(cArrInd)%molecs(molec)%cube_list(1:natoms)=&
          CurrentCubes(1:natoms)


    Case(110)
      ! WE ARE DELETING THIS MOLECULE
      ! The following changes are done to cube_array
      ! - remove old ('molec') from all cubes previously occupied by it
      ! - change the index of last molecule('lastMolec') in array to 'molec'
      ! - remove 'lastMolec' from all arrays occupied by it, and add 
      !   'molec' to the same arrays

      !** 1.) Remove old molecule Info from cubelist
      Do i=1,natoms
        old_c=cubelist%sorbs(cArrInd)%molecs(molec)%cube_list(i)
        Call cavitylist_deleteMolec(cubelist, cArrInd, molec, old_c )
      End Do

      !** 2.) change the index of last molecule from lastMolec to molec
      lastMolec = cubelist%sorbs(cArrInd)%nmolecs 
      Do i=1,natoms
        last_c=cubelist%sorbs(cArrInd)%molecs(lastMolec)%cube_list(i)
        Call cavitylist_addMolec(cubelist, cArrInd, molec, last_c )
        Call cavitylist_deleteMolec(cubelist, cArrInd, lastMolec, last_c )
      End Do

      cubelist%sorbs(cArrInd)%nmolecs = cubelist%sorbs(cArrInd)%nmolecs -1
      !** finally copy the cubelist in sorbinfo
      cubelist%sorbs(cArrInd)%molecs(molec)%cube_list(1:natoms)=  &
          cubelist%sorbs(cArrInd)%molecs(lastMolec)%cube_list(1:natoms)
    Case(111)

      !** get index of cube for each atom
      Call cavitylist_getCubeIndices(cubelist, sorbate, molec, &
          natoms,CurrentCubes(1:natoms))
      old_list(1:natoms)=&
          cubelist%sorbs(cArrInd)%molecs(molec)%cube_list(1:natoms)

      !** Transformed molecules
      Do i=1,natoms
        new_c=CurrentCubes(i)
        old_c=old_list(i)
        If (old_c==0) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "The molecule was not inducted into the cavity before ?"
          Stop
        Else
          If (new_c>0) Then
            ! ** change new cube
            Call cavitylist_addMolec(cubelist, cArrInd, molec, new_c )
          Else
            Write(*,*) "The new suggested point is lying in an untabulated &
                & area, check whether the pmaps/bmaps are consistent"
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
          Endif

          !** change basic info in sorbs also
          cubelist%sorbs(cArrInd)%molecs(molec)%cube_list(i)=new_c
        Endif

      End Do

      !** Transformed molecules
      Do i=1,natoms
        new_c=CurrentCubes(i)
        old_c=old_list(i)
        If (new_c==old_c) Then
          ! do nothing, already added it
          Cycle
        ElseIf (old_c==0) Then
          Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
          Write(*,*) "The molecule was not inducted into the cavity before ?"
          Stop
        Else
          ! check whethter this should be really deleted 
          posn_in_array=findint(Currentcubes(1:natoms),old_c)
          If (posn_in_array==0) Then
            !** delete from old cube
            Call cavitylist_deleteMolec(cubelist, cArrInd, molec, old_c )
          Endif
        Endif

      End Do



    Case default
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop

    End Select

  End Subroutine Cavitylist_updateMolecInfo

  !-----------------------------------------------------------------
  ! Finds all cubes occupied by each atom of the molecule
  !-----------------------------------------------------------------
  Subroutine cavitylist_getCubeIndices(cubelist, sorbate, molec, &
      natoms,cubes)
    Type(Cavity_Params), Intent(in)              :: cubelist
    Type(AtMolCoords), Intent(in)                   :: sorbate
    Integer, Intent(in)                             :: molec, natoms 
    Integer, Dimension(:), Intent(inout) :: cubes
    Integer :: i

    Do i=1,natoms
      cubes(i)=cavitylist_getindexVec(cubelist, sorbate%coords(i,molec)%r)
    End Do

  End Subroutine cavitylist_getCubeIndices


  !----------------------------------------------------------------
  ! Returns the index of the cube in cube_array from a vector.
  !----------------------------------------------------------------
  Integer Function cavitylist_getIndexVec(cubelist, vec)
    Type(Cavity_Params), Intent(in)              :: cubelist
    Type(VecType), Intent(in)                    :: vec
    Integer :: ix, iy, iz

    ix= Int(vec%comp(1)/cubelist%dx) +1
    iy= Int(vec%comp(2)/cubelist%dy) +1
    iz= Int(vec%comp(3)/cubelist%dz) +1
    cavitylist_getIndexVec=cubelist%cube_index(ix, iy, iz)

  End Function cavitylist_getIndexVec

  !-----------------------------------------------------------------
  ! Make sure that all fields have enough memory.
  !-----------------------------------------------------------------
  Subroutine cavitylist_checkandIncrMemory(cubelist, sorbate, sorb )
    Type(Cavity_Params), Intent(inout)              :: cubelist
    Type(AtMolCoords), Intent(in)  :: sorbate
    Integer, Intent(in) :: sorb ! , molec
    Integer :: cArrInd, oldsize, newsize, error, m, natoms

    Type(Molec_Info), Dimension(:), Pointer :: TempMolecArray

    cArrInd=cubelist%csorblist(sorb)
    natoms=molecules_getnatoms(sorb)

    oldsize=cubelist%sorbs(cArrInd)%currentsize 
    ! we give a padding of 2
    If (cubelist%sorbs(cArrInd)%nmolecs>= (oldsize-2)) Then
      ! we are reaching the limits so increase the size of molecPtr 
      ! and molecPtrIndex
      newsize=oldsize+MEMORY_STEP_SIZE
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Increasing memory in in cavitylist"
      Write(*,*) "Current number of molecules", config_getnmoles(sorbate)
      Write(*,*) " Current capacity for molecs in cavitylist ", size(cubelist%sorbs(cArrInd)%molecs)
      Write(*,*) "New SIze", newsize


      !** Pointing temperory pointers to the old memory
      TempMolecArray=>cubelist%sorbs(cArrInd)%molecs

      !** Release our permanent arrays from old memory
      Nullify(cubelist%sorbs(cArrInd)%molecs)

      !** Allocate new memory for our permanent arrays
      Allocate(cubelist%sorbs(cArrInd)%molecs(newsize), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      cubelist%sorbs(sorb)%currentsize=newsize

      !** Copy fields into the new permanent array
      ! fill default values here
      Do m=1,newsize
        If (m<=oldsize) Then
          Call cavitylist_copyMolecInfo(cubelist%sorbs(cArrInd)%molecs(m), &
              TempMolecArray(m), sorb )
        Else
          ! give default values
          cubelist%sorbs(cArrInd)%molecs(m)%ncubes=0
          Allocate(cubelist%sorbs(cArrInd)%molecs(m)%cube_list(natoms), &
              STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
          cubelist%sorbs(cArrInd)%molecs(m)%cube_list=0
        Endif
      End Do

      !** Deallocate all extra memory
      Deallocate(TempMolecArray , STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)


    Endif

  End Subroutine cavitylist_checkandIncrMemory



  !-----------------------------------------------------------------
  ! Copies the values from oldinfo to newinfo 
  ! assumes that oldinfo was initialized, and newinfo was not
  !-----------------------------------------------------------------
  Subroutine cavitylist_copyMolecInfo(newinfo, oldinfo, sorb)
    Type(Molec_Info), Intent(inout) :: newinfo, oldinfo
    Integer , Intent(in) :: sorb
    Integer :: natoms, error

    natoms=molecules_getnatoms(sorb)
    If (Size(oldinfo%cube_list)/=natoms) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    newinfo%ncubes=oldinfo%ncubes

    Allocate(newinfo%cube_list(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    newinfo%cube_list(1:natoms)=  oldinfo%cube_list(1:natoms)

  End Subroutine cavitylist_copyMolecInfo


  !-----------------------------------------------------------------
  ! Deletes a molecules from the given cavity (cavity index = c )
  ! sorb is the index in cavitylist not in molecules or config
  !-----------------------------------------------------------------
  Subroutine cavitylist_deleteMolec(cubelist, sorb, molec, c )
    Type(Cavity_Params), Intent(inout)              :: cubelist
    Integer, Intent(in) :: sorb, molec, c

    Integer :: nmolecs, j, k
    Logical :: removed
    removed=.False.

    nmolecs=cubelist%cube_array(c)%nmolecs(sorb)

    If (nmolecs==0) Then
      ! maybe it was already removed while updating for a previous atom
      Return
    Elseif (nmolecs==1) Then
      !** Only one remove It
      cubelist%cube_array(c)%nmolecs(sorb)=0
      removed=.True.
    Elseif(nmolecs<0) Then
      !** if nmoelcs==0 delete should not have been attempted
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop

    Elseif(cubelist%cube_array(c)%moleclist(sorb,nmolecs)==molec)Then
      !** The last one is to be removed, no shifting needed
      cubelist%cube_array(c)%nmolecs(sorb)=nmolecs-1
      removed=.True.
    Else
      !** some body in middle has to removed, shifting of array reqd
      Do j=1,nmolecs-1

        ! find where is this molecule placed
        If (cubelist%cube_array(c)%moleclist(sorb,j)==molec)Then
          cubelist%cube_array(c)%nmolecs(sorb)=nmolecs-1

          ! shift all indices above this particular one  and 
          ! compress the moleclist
          Do k=j,nmolecs-1
            cubelist%cube_array(c)%moleclist(sorb,k)= &
                cubelist%cube_array(c)%moleclist(sorb,k+1)
          End Do

          removed=.True.
          Exit        !out of the outer j loop
        Endif
      End Do

    Endif  !** end of nmolecs check

    If (.Not.removed) Then
      Write(*,*) nmolecs
      Write(*,*) molec, cubelist%cube_array(c)%moleclist(sorb,1:nmolecs)
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

  End Subroutine cavitylist_deleteMolec

  !-----------------------------------------------------------------
  ! Adds a molecules to the given cavity (cavity index = c )
  ! sorb is the index in cavitylist not in molecules or config
  !-----------------------------------------------------------------
  Subroutine cavitylist_addMolec(cubelist, sorb, molec, c )
    Type(Cavity_Params), Intent(inout)              :: cubelist
    Integer, Intent(in) :: sorb, molec, c

    Integer :: nmolecs, j
    Logical :: addthis

    nmolecs=cubelist%cube_array(c)%nmolecs(sorb)
    addthis=.True.
    If (nmolecs>0) Then
      Do j=1,nmolecs
        If (cubelist%cube_array(c)%moleclist(sorb,j)==molec) Then
          addthis=.False.
          Exit
        Endif
      End Do
    Endif

    If (addthis) Then
      If ((nmolecs+1)>MAX_MOLECS_PER_CUBE) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
      cubelist%cube_array(c)%nmolecs(sorb)=nmolecs+1
      cubelist%cube_array(c)%moleclist(sorb,nmolecs+1)=molec
    Endif

  End Subroutine cavitylist_addMolec



  !-----------------------------------------------------------------
  ! Check whether the given sorb is used for  cavity anaylsis
  !-----------------------------------------------------------------
  Logical Function cavitylist_checkCavitySorb(cubelist, sorb)
    Type(Cavity_Params), Intent(in)              :: cubelist
    Integer, Intent(in) :: sorb
    cavitylist_checkCavitySorb=cubelist%is_cavity_sorb(sorb)
  End Function cavitylist_checkCavitySorb

  !-------------------------------------------------
  ! Initializes all the individual cubes when there is no pmap present 
  !-------------------------------------------------
  Subroutine cavitylist_initArrUniform(cubelist)
    Type(Cavity_Params), Intent(inout)              :: cubelist
    Integer :: ncubes, cube_no, i, j, k,error

    ! No pmap is given so we have to take all cubes into account
    ncubes= cubelist%nx * cubelist%ny * cubelist%nz
    cubelist%ncubes=ncubes
    Allocate(cubelist%cube_array(ncubes),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__, &
        Trim("cube_array ncubes= "//Trim(int2str(ncubes))))

    !** Loop thru all cubes and assign cube_index 
    cube_no=0
    Do i=1,cubelist%nx
      Do j=1,cubelist%ny
        Do k=1,cubelist%nz
          cube_no=cube_no+1
          cubelist%cube_index(i,j,k) = cube_no
        End Do
      End Do
    End Do

    !** Loop thru again this time initializing cube_array
    Allocate(cubelist%cube_array(cube_no), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    cube_no=0
    Do i=1,cubelist%nx
      Do j=1,cubelist%ny
        Do k=1,cubelist%nz
          If (cubelist%cube_index(i,j,k)/=0) Then
            cube_no=cube_no+1
            Call cavitylist_initcubeArray( cubelist, cube_no, (/i,j,k /) )
          Endif
        End Do
      End Do
    End Do

  End Subroutine cavitylist_initArrUniform

  !-------------------------------------------------
  ! Initializes all the individual cubes when there is a pmap specified 
  ! It avoids all cubes( of cubelist) in which all points( from pmap) 
  ! have energy above the specified cutoff
  ! Note that a pmap is usually specified for a fundcell
  ! but the simcell might contain more thatn one fundcell.
  !-------------------------------------------------
  Subroutine cavitylist_initArrPmap(cubelist)
    Type(Cavity_Params), Intent(inout)              :: cubelist

    Integer :: ncubes, cube_no, i, j, k, nbrx, nbry, nbrz, nfsize, error
    Integer :: i_cube, j_cube, k_cube, i_index, j_index, k_index, mapunit
    Integer :: nrgArrIndex, cube_error, i_ucx, j_ucy, k_ucz, max_cubes

    Real(kind=RDbl) :: fcx, fcy, fcz, x, y, z
    Character(len=lstrLen) :: mapfilename

    !** Temperory variables for reading the pmap
    Real(kind=RSgl), Dimension(:), Pointer :: map_cubearray
    Integer, Dimension(:,:,:), Pointer     :: map_cubeindex

    ! we are given a pmap, first open , 100 -> old unformatted file
    mapfilename=Trim(Trim(getpath('PMAPDIR'))//Trim(cubelist%pmapName))
    mapunit = file_open(mapfilename,100)

    ! fundcell dimensions
    fcx = cubelist%lx/cubelist%nucx
    fcy = cubelist%ly/cubelist%nucy
    fcz = cubelist%lz/cubelist%nucz

    !** Get the bmap header information, We could check whether  
    Call mapheader_getheader(cubelist%header, &
        (/ fcx, fcy, fcz /), mapunit)

    !** Allocate the necessary memory
    nbrx = cubelist%header%nbrx
    nbry = cubelist%header%nbry
    nbrz = cubelist%header%nbrz
    nfsize = cubelist%header%nfsize

    !** Now we are going to allocate the temperory arrays and read 
    !** energies into them
    Allocate(map_cubeindex(nbrz,nbry,nbrx), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"map_cubes")
    Allocate(map_cubearray(nfsize), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"map_cubes")


    !** Now read in the potential/force etc.. and close the file
    Read(mapunit) map_cubeindex
    Read(mapunit) map_cubearray
    Close(mapunit)
    max_cubes = cubelist%nx * cubelist%ny * cubelist%nz
    Write(*,*) " Maximum expected number of cubes",max_cubes
    cube_no=0


    !** We are going to scan pmap and add only those cubes from cavitylist 
    !** which includes the pmap-cube within some cutoff. 
    !** For this to work the pmap cubes should be smaller than cavitylist-cubes
    If (cubelist%nx>nbrx) Then
      Write(*,*) "Too many cavity divisions in x-axis, it should be < ", &
          nbrx-1
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    If (cubelist%ny>nbry) Then
      Write(*,*) "Too many cavity divisions in y-axis, it should be < ", &
          nbry-1
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    If (cubelist%nz>nbrz) Then
      Write(*,*) "Too many cavity divisions in z-axis, it should be < ", &
          nbrz-1
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif



    !** nbrx is number of nodes, nbrx-1 is number of cubes, In the Pmap
    Do i=1,nbrx-1
      Do j=1,nbry-1
        Do k=1,nbrz-1
          nrgArrIndex=map_cubeindex(k,j,i)
          If (nrgArrIndex/=0) Then
            ! that means the cube is tabulated in the pmap
            If (map_cubearray(nrgArrIndex)<cubelist%ucutoff) Then
              ! OK, we have found a point in the pmap that should
              ! be considered for cavity analysis. Now we have to  
              ! find the cube in which this point lies

              !** By defn in mapheader these are not midpoints of the 
              !** cube, rather they at the corner of the cube closest 
              !** to the origin, hence the addn of cubelength/2
              x=cubelist%header%xgrid(i)+cubelist%header%stepx/2
              y=cubelist%header%ygrid(j)+cubelist%header%stepy/2
              z=cubelist%header%zgrid(k)+cubelist%header%stepz/2

              !** Now we need to find the corresponding indices in cubelist
              ! Note that there are more than one unitcell in the simcell
              !  all points below are in the first unit cel
              i_cube=Int(x/cubelist%dx)+1
              j_cube=Int(y/cubelist%dy)+1
              k_cube=Int(z/cubelist%dz)+1

              cube_error=0
              If ((i_cube<0).Or.(i_cube>cubelist%nx)) cube_error=1
              If ((j_cube<0).Or.(j_cube>cubelist%ny)) cube_error=1
              If ((k_cube<0).Or.(k_cube>cubelist%nz)) cube_error=1
              If (cube_error==1) Then
                Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
                Write(*,*) i_cube, j_cube, k_cube
                Stop
              Endif


              !** Check whether we have already realised that this 
              !** cube is good
              If (cubelist%cube_index(i_cube,j_cube,k_cube)/=0) Then
                Cycle
              Else

                ! we found a new cube
                ! Do the initialization for cube here and all other unitcells
                Do i_ucx=0,cubelist%nucx-1
                  Do j_ucy=0,cubelist%nucy-1
                    Do k_ucz=0,cubelist%nucz-1

                      cube_no=cube_no+1
                      If (cube_no> max_cubes) Then
                        Write(*,*) "something terribly wrong with cubes"
                        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
                        Stop
                      Endif

                      !** Set the fields in individual cubes
                      i_index=i_cube+(cubelist%nx/cubelist%nucx)*i_ucx
                      j_index=j_cube+(cubelist%ny/cubelist%nucy)*j_ucy
                      k_index=k_cube+(cubelist%nz/cubelist%nucz)*k_ucz

                      !** This assignment is temperory (Just to mark 
                      ! occupied and unoccupied cubes. Will be changed during
                      !   cube initialization loop
                      cubelist%cube_index(i_index,j_index,k_index)=1

                    End Do
                  End Do
                End Do  !** end of unicell loop

              Endif
            Endif

          Endif !** end of -- If (nrgArrIndex/=0 )--

        End Do
      End Do
    End Do !** End of mapcubes loop

    Deallocate(map_cubeindex)
    Deallocate(map_cubearray)

    ncubes= cube_no
    cubelist%ncubes=ncubes

    !** Now we have successfully identified all the cube_indices that are 
    !** accessible, let us allocate and initialize the cube_array 

    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,'(a,i10,a)') " Allocating ", ncubes, &
        " cubes for stroring cavity info"
    Write(*,'(a,f12.6,a)') "cube array shrank to ", &
        (ncubes*one*100)/max_cubes,  "% because of the use of pmap"

    Allocate(cubelist%cube_array(cube_no), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    cube_no=0
    Do i=1,cubelist%nx
      Do j=1,cubelist%ny
        Do k=1,cubelist%nz
          If (cubelist%cube_index(i,j,k)/=0) Then
            cube_no=cube_no+1
            Call cavitylist_initcubeArray( cubelist, cube_no, (/i,j,k /) )

            !** Note : Previously we were just counting through the number of 
            !          cubes. Since we cycled in a different order there, 
            !           we need to correct the cube_index array here
            cubelist%cube_index(i, j, k)=cube_no

          Endif
        End Do
      End Do
    End Do

  End Subroutine cavitylist_initArrPmap


  !-------------------------------------------------
  ! Initializes the cube- arrIndex in the cube_array
  ! This is done individually for all cubes. So that unvisited (from bmap) 
  ! cubes are never initialized
  !-------------------------------------------------
  Subroutine cavitylist_initcubeArray(cubelist, arrIndex, cubeIndices  )
    Type(Cavity_Params), Intent(inout)              :: cubelist
    Integer, Intent(in)                             :: arrIndex
    Integer, Dimension(3), Intent(in)               :: cubeIndices

    Integer :: error
    !** Check whether already inited
    If (cubelist%cube_array(arrIndex)%is_init) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Else
      cubelist%cube_array(arrIndex)%is_init=.True.
    Endif

    cubelist%cube_array(arrIndex)%ix=cubeIndices(1)
    cubelist%cube_array(arrIndex)%iy=cubeIndices(2)
    cubelist%cube_array(arrIndex)%iz=cubeIndices(3)

    ! ** Just the list of sorbs for keeping track 
    cubelist%cube_array(arrIndex)%sorblist=> cubelist%sorblist

    !** array for storing number of molecules of each sorb
    Allocate(cubelist%cube_array(arrIndex)%nmolecs(cubelist%nsorbs),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    cubelist%cube_array(arrIndex)%nmolecs=0

    !** Index of each of those molecules, big memory use here, WATCH MEMORY
    Allocate(cubelist%cube_array(arrIndex)%moleclist(cubelist%nsorbs, &
        MAX_MOLECS_PER_CUBE),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    cubelist%cube_array(arrIndex)%moleclist=0
  End Subroutine cavitylist_initcubeArray

  !-------------------------------------------------
  ! Read the necessary fields from ctrlfile
  !-------------------------------------------------
  Subroutine cavitylist_readCtrlFile(cubelist, unitno )
    Type(Cavity_Params), Intent(inout)              :: cubelist
    Integer, Intent(in)                             :: unitno

    Integer :: i,error,sorbtype, totsorbs

    !** Read the names of all sorbates for cavity considerations 
    Read(unitno, *) cubelist%nsorbs

    Allocate(cubelist%sorbnames(cubelist%nsorbs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__, "Cube list")
    Allocate(cubelist%sorblist(cubelist%nsorbs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__, "Cube list")

    totsorbs=molecules_getnsorbs()
    Allocate(cubelist%is_cavity_sorb(totsorbs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__, "Cube list")
    Allocate(cubelist%csorblist(totsorbs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__, "Cube list")
    cubelist%is_cavity_sorb=.False.
    cubelist%csorblist=0
    Do i=1,cubelist%nsorbs
      ! Read the name of the sorbate
      Read(unitno,'(a)') cubelist%sorbnames(i)
      cubelist%sorbnames(i)=Trim(stripcmnt(cubelist%sorbnames(i)))
      sorbtype=molecules_gettype(cubelist%sorbnames(i))
      If (sorbtype==0) Then
        Write(*,*) "wrong sorb name in cavitylist section"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Else
        ! store the type of the sorbate
        cubelist%sorblist(i)=sorbtype
      Endif

      !** keep the revers eindices also
      cubelist%csorblist(sorbtype) =i

      If (.Not.(cubelist%is_cavity_sorb(sorbtype))) Then
        cubelist%is_cavity_sorb(sorbtype)=.True.
      Else
        Write(*,*) "Duplication ?"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

    End Do

    !** Initialize details about the geometry of the unit cell
    Read(unitno, *) cubelist%nx, cubelist%ny, cubelist%nz

    !** Now we need to calculate the required dimension for cube_array
    !** First let us read the ctrlfile and see whether any pmap is specified
    Read(unitno, *) cubelist%pmapName
    cubelist%pmapName=Trim(stripcmnt(cubelist%pmapName))

    If (Trim((cubelist%pmapName))/="NULL") Then 
      Read(unitno,*) cubelist%ucutoff
    Endif

  End Subroutine cavitylist_readCtrlFile

  !----------------------------------------------------
  ! Gets the total no. of cubelets in the unit cell
  !----------------------------------------------------
  Integer Function cavitylist_getncubes(cubelist)
    Type(Cavity_Params), Intent(in) :: cubelist
    cavitylist_getncubes = cubelist%ncubes
  End Function cavitylist_getncubes


  !--------------------------------------------------------------------
  ! Check whether the given position lies in  a cavity
  ! point - position to check
  !--------------------------------------------------------------------
  Logical Function cavitylist_checkcavity1(cubelist, point)
    Type(Cavity_Params), Pointer :: cubelist
    Type(VecType),Intent(in)     :: point 
    Integer :: index
    index=cavitylist_getindexvec(cubelist, point)
    cavitylist_checkcavity1=cavitylist_checkcavityFromIndex(cubelist,index)

  End Function cavitylist_checkcavity1

  !--------------------------------------------------------------------
  ! Check whether the given position lies in  a cavity
  ! If point is in one of the cubes2avoid then it is sid 
  ! to be not a cavity
  !  point - position to check
  !  cubes2avoid - these are empty sites, parts of cavity
  !--------------------------------------------------------------------
  Logical Function cavitylist_checkcavity2(cubelist, point, cubes2avoid )
    Type(Cavity_Params), Pointer :: cubelist
    Type(VecType),Intent(in)     :: point 
    Integer, Dimension(:), Intent(in) ::  cubes2avoid
    Integer :: index, arr_index
    Logical :: is_cavity

    index=cavitylist_getindexvec(cubelist, point)
!    Write(*,*) " to avoid :: ", cubes2avoid
!    Write(*,*) " index    :: ", index 

    is_cavity=.True.              
    arr_index=findint(cubes2avoid,index)     
    If (arr_index==0) is_cavity=cavitylist_checkcavityFromIndex(&
        cubelist,index)
!    Write(*,*) "is_cavity :: ", is_cavity
    cavitylist_checkcavity2=is_cavity

  End Function cavitylist_checkcavity2


  !--------------------------------------------------------------------
  ! Check whether the given position ( corrsponding to index) is a cavity
  !--------------------------------------------------------------------
  Logical Function cavitylist_checkcavityfromIndex(cubelist, index)
    Type(Cavity_Params), Pointer :: cubelist
    Integer, Intent(inout)       :: index
    Integer :: i

    cavitylist_checkcavityfromIndex=.True.

    If (index<1) Then
      ! If its lying in an untabulated place then that is not a cavity
      cavitylist_checkcavityfromIndex=.False.
      Return
    Endif

    Do i=1,cubelist%nsorbs
      If (cubelist%cube_array(index)%nmolecs(i)>0) Then
        cavitylist_checkcavityfromIndex=.False.
        Return
      Endif
    End Do

  End Function cavitylist_checkCavityFromIndex



  !-------------------------------------------------------
  ! Write sample control file
  !-------------------------------------------------------
  Subroutine cavitylist_samplecf(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a )') " --- "//Trim(default_cavitylist_tag)//" --- "
    Write(unitno,'(a,t30,a)') 'Integer',&
        '# Number of sorbates for cavity analysis'
    Write(unitno,'(a,t30,a)') 'Character',&
        '# Names of sorbates [next few lines too....]'
    Write(unitno,*) 
    Write(unitno,'(a,t30,a)') 'Integer, Integer, Integer',&
        '# Number of divisions along each axes, nx,ny, nz'
    Write(unitno,'(a,t30,a)') 'Character',&
        '# Name of a pmap(ex=sili.Methyl.pmap) or  the keyword- NULL '
    Write(unitno,'(a,t30,a)') 'Real',&
        '# Potential cut off in KJ/molName , not relevent if NULL'
  End Subroutine cavitylist_samplecf


  !-------------------------------------------------------
  ! Display the bias map information
  !-------------------------------------------------------
  Subroutine cavitylist_display(cubelist, nspc, optunit)
    Type(Cavity_Params), Intent(in) :: cubelist 
    Integer, Intent(in)  :: nspc ! No. of spaces from the left column
    Integer, Optional, Intent(in)    :: optunit

    Integer   :: unitno
    Character(len=nspc) :: spc

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) "Fill in the cavity bias stuff here"
!    stop

  End Subroutine cavitylist_display


!---------------- Some routines rsed for visualization during debugging -----!
  !-----------------------------------------------------------------
  ! Displays the cavitylist as an xyz file for display using rasmol
  ! An atom is displayed at the center of each cube. This is just for
  ! testing purposes not to be used by any other modules as of now
  ! Requires - 
  !          cavitylist     - the map
  !          xyzfilename - name of rasmol file 
  !          descr       - small description to be written to rasmol file
  !          atsym       - symbol of atom(ex: "C" for carbon, "CH3" for methyl
  !-----------------------------------------------------------------
  Subroutine cavitylist_makeXYZ( cubelist, xyzfilename, descr, atsym)
    Type(Cavity_Params), Intent(in)  :: cubelist
    Character(*)             :: xyzfilename, descr
    Character(len=3)                  :: atsym

    Integer         :: unitno, i, j, k
    Real(kind=RDbl) :: dx, dy, dz

    Write(*,*) "Displaying cavitylist as an xyz file into file : "//Trim(xyzfilename)
    unitno=file_open(xyzfilename)
    Write(unitno,*) cubelist%ncubes
    Write(unitno,*) descr
    dx=cubelist%dx
    dy=cubelist%dy
    dz=cubelist%dz
    Do i=1, cubelist%nx
      Do j=1, cubelist%ny
        Do k=1, cubelist%nz
          If (cubelist%cube_index(i,j,k)/=0) Then
            Write(unitno,'(a,t10,3f14.5)') atsym, i*dx-dx/2, &
                j*dy-dy/2, k*dz-dz/2  
          Endif
        End Do
      End Do
    End Do
    close(unitno)

  End Subroutine cavitylist_makeXYZ

  !-----------------------------------------------------------------
  ! Displays the occupied cubes  as an xyz file for display using rasmol
  ! An atom is displayed at the center of each cube. This is just for
  ! testing purposes not to be used by any other modules as of now
  ! Requires - 
  !          cavitylist     - the map
  !          xyzfilename - name of rasmol file 
  !          descr       - small description to be written to rasmol file
  !          atsym       - symbol of atom(ex: "C" for carbon, "CH3" for methyl
  !-----------------------------------------------------------------
  Subroutine cavitylist_displayCubeInfo( cubelist, xyzfilename, descr, atsym)
    Type(Cavity_Params), Pointer  :: cubelist
    Character(*), Intent(in)             :: xyzfilename, descr
    Character(len=3), Intent(in)                  :: atsym

    Integer         :: unitno, i, j, k, count
    Real(kind=RDbl) :: dx, dy, dz
    Type(VecType) :: rvec

    !** 1. Open Output file
    Write(*,*) "Displaying Occupied as an xyz file into file : "&
        //Trim(xyzfilename)
    unitno=file_open(xyzfilename)

    !** 2. Count cubes
    count=0
    dx=cubelist%dx
    dy=cubelist%dy
    dz=cubelist%dz
    Do i=1, cubelist%nx
      Do j=1, cubelist%ny
        Do k=1, cubelist%nz
          If (cubelist%cube_index(i,j,k)/=0) Then
            rvec=(/ i*dx-dx/2, j*dy-dy/2, k*dz-dz/2 /)
            If (.Not.(cavitylist_checkcavity(cubelist, rvec))) Then
              count=count+1
            Endif
          Endif
        End Do
      End Do
    End Do

    !** 3. Write the output file
    Write(unitno,*) count
    Write(unitno,*) descr
    Do i=1, cubelist%nx
      Do j=1, cubelist%ny
        Do k=1, cubelist%nz
          If (cubelist%cube_index(i,j,k)/=0) Then
            rvec=(/ i*dx-dx/2, j*dy-dy/2, k*dz-dz/2 /)
            If (.Not.(cavitylist_checkcavity(cubelist, rvec))) Then
              Write(unitno,'(a,t10,3f14.5)') atsym, rvec%comp(1:3) 
              Write(*,*) "occupied cube : ", i,j,k,cubelist%cube_index(i,j,k) 
            Endif
          Endif
        End Do
      End Do
    End Do
    close(unitno)

  End Subroutine cavitylist_displayCubeInfo


  !-----------------------------------------------------------------
  ! Displays the occupied cubes  as an xyz file for display using rasmol.
  ! The diffrence between this and previous is that, now it is based 
  ! on the info stored in SorbInfo field of the cubelist
  ! An atom is displayed at the center of each cube. This is just for
  ! testing purposes not to be used by any other modules as of now
  ! Requires - 
  !          cavitylist     - the map
  !          xyzfilename - name of rasmol file 
  !          descr       - small description to be written to rasmol file
  !          atsym       - symbol of atom(ex: "C" for carbon, "CH3" for methyl
  !-----------------------------------------------------------------
  Subroutine cavitylist_displaySorbInfo( cubelist, sorbates, &
      xyzfilename, descr,  atsym)
    Type(Cavity_Params), Intent(in)  :: cubelist
    Type(AtMolCoords), Dimension(:),Intent(in):: sorbates
    Character(*),Intent(in)             :: xyzfilename, descr
    Character(len=3),Intent(in)                  :: atsym

    Integer         :: unitno, i, ix, iy, iz
    Integer         :: m, a, count, natoms
    Real(kind=RDbl) :: dx, dy, dz
    Type(VecType)   :: rvec
    Integer, Dimension(MAX_ATOMS), Save :: CurrentCubes

    !** 1. Open Output file
    Write(*,*) "Displaying Molecules as an xyz file into file : "&
        //Trim(xyzfilename)
    unitno=file_open(xyzfilename)
    dx=cubelist%dx
    dy=cubelist%dy
    dz=cubelist%dz

    !** 2. Count cubes
    count=0
    Do i=1, molecules_getnsorbs()
      If (.Not.(cavitylist_checkCavitySorb(cubelist, i))) Cycle
      natoms=molecules_getnatoms(i)
      count=count+config_getnmoles(sorbates, i)*natoms 
    End Do

    !** 3. Write the output file
    Write(unitno,*) count
    Write(unitno,*) descr

!!$    startdebug=.True.

    Do i=1, molecules_getnsorbs()
      If (.Not.(cavitylist_checkCavitySorb(cubelist, i))) Cycle
      natoms=molecules_getnatoms(i)
      Do m=1,config_getnmoles(sorbates, i)
        Call cavitylist_getCubeIndices(cubelist, sorbates(i), m, &
            natoms, CurrentCubes(1:natoms) )
        Do a=1,natoms
          ix=cubelist%cube_array(CurrentCubes(a))%ix
          iy=cubelist%cube_array(CurrentCubes(a))%iy
          iz=cubelist%cube_array(CurrentCubes(a))%iz
          rvec=(/ ix*dx-dx/2, iy*dy-dy/2, iz*dz-dz/2 /)
          Write(unitno,'(a,t10,3f14.5)') atsym, rvec%comp(1:3)         
              Write(*,*) "sorb at : ", i,m,a,ix,iy,iz,currentcubes(a)
        End Do
      End Do
    End Do

    Close(unitno)

  End Subroutine cavitylist_displaySorbInfo

End Module cavitylist
