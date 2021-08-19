!--------------------------------------------------------
! This modules is a front end to the various generalized
! coordinate models supported.
!--------------------------------------------------------
Module gcmodels

  Use defaults, Only: RDbl, strLen, MAX_DOF, lstrLen
  Use utils, Only: toupper
  Use vector, Only: VecType
  Use molecules, Only: molecules_getgcmodeltype
  Use simcell, Only: SimCell_Params
  Use rigidcoords, Only: RigidMolec,rigidcoords_displaystr,rigidcoords_getcom,&
      rigidcoords_getgencoords, rigidcoords_initcoords, rigidcoords_clean, &
      rigidcoords_gcoords_eq_vectype, rigidcoords_toxyz, rigidcoords_fromxyz, &
      rigidcoords_setgencoords, rigidcoords_getmetrictensor, rigidcoords_pbc, &
      rigidcoords_readrestartfile, rigidcoords_dumprestartinfo, &
      rigidcoords_convstr, rigidcoords_getgenforce, rigidcoords_copy, &
      rigidcoords_changeflex, rigidcoords_changedof, rigidcoords_setfrom_rp
  Use branchedcoords, Only: BranchedMolec,NodeType,branchedcoords_initcoords, &
      branchedcoords_toxyz, branchedcoords_place, &
      branchedcoords_unplace, branchedcoords_placenode, &
      branchedcoords_readrestartinfo, branchedcoords_dumprestartinfo, & 
      branchedcoords_coord, branchedcoords_copycoords 

  Implicit None
  Save

  Private
  Public :: GeneralizedCoords, gcmodels_initcoords, gcmodels_readrestartfile, &
      gcmodels_dumprestartinfo, gcmodels_getcom, Assignment(=), &
      gcmodels_getgencoords, gcmodels_display, gcmodels_getgenforce, &
      gcmodels_getmetrictensor, gcmodels_setgencoords, gcmodels_changeflex, &
      gcmodels_toxyz, gcmodels_pbc, gcmodels_fromxyz, gcmodels_changedof, &
      gcmodels_setfrom_rp, gcmodels_clean

  Type GeneralizedCoords
    Type(RigidMolec), Pointer    :: rigidcoords
    Type(BranchedMolec), Pointer :: branchedcoords
    Logical :: nogcoords 
  End Type GeneralizedCoords

  Interface Assignment(=)
    Module Procedure gcmodels_copycoords
    Module Procedure gcmodels_gcoords_eq_vectype
  End Interface

  Interface gcmodels_display 
    Module Procedure gcmodels_displaystr 
  End Interface

Contains
  !-----------------------------------------------------------------
  ! Initialize the gcmodel coordinate models based on the modeltype 
  ! string.
  ! Requires: gcoords -- generalized coordinates
  !           spc -- species number
  !           ref_struc -- optional reference structure
  !----------------------------------------------------------------
  Subroutine gcmodels_initcoords(gcoords,spc,ref_struc)
    Type(GeneralizedCoords), Intent(InOut)  :: gcoords
    Integer, Intent(In)                     :: spc
    Type(VecType), Dimension(:), Optional   :: ref_struc

    Character(len=strLen)    :: modeltype

    !** Nullify all the pointer instances in gcmodel
    Nullify(gcoords%rigidcoords)
    Nullify(gcoords%branchedcoords)
    gcoords%nogcoords=.False.
    modeltype = molecules_getgcmodeltype(spc)

    !** Call the appropriate routine to complete initialization
    Select Case (Trim(toupper(modeltype)))
    Case('RIGID')
      If (Present(ref_struc)) Then
        Call rigidcoords_initcoords(gcoords%rigidcoords,spc,ref_struc)
      Else
        Call rigidcoords_initcoords(gcoords%rigidcoords,spc)
      End If
    Case('BRANCHED')
      Call branchedcoords_initcoords(gcoords%branchedcoords,spc)
    Case('NONE')
      gcoords%nogcoords=.True.
    Case Default
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          "could not understand gcmodel type : "//Trim(modeltype)
      Stop
    End Select

  End Subroutine gcmodels_initcoords

  !---------------------------------------------------------------
  ! This subroutine COPIES the CONTENTS of gcoords2 into gcoords1
  ! The default assignment operator just copies POINTERS
  !---------------------------------------------------------------
  Subroutine gcmodels_copycoords(gcoords1, gcoords2)
    Type(GeneralizedCoords), Intent(inout) :: gcoords1
    Type(GeneralizedCoords), Intent(in)    :: gcoords2

    If (Associated(gcoords2%rigidcoords)) Then 
      Call rigidcoords_copy(gcoords1%rigidcoords, gcoords2%rigidcoords)
    Else If (Associated(gcoords2%branchedcoords)) Then
      Call branchedcoords_copycoords(gcoords1%branchedcoords,& 
          gcoords2%branchedcoords ) 
    End If
  End Subroutine gcmodels_copycoords

  !----------------------------------------------------------
  ! This subroutine assigns a vector type to the generalized
  ! coordinate type.  Its functionality is dependent on the
  ! generalized coordinate type.  For the "rigid" molecule it
  ! assigns the center-of-mass of the molecule to "vec"
  !----------------------------------------------------------
  Subroutine gcmodels_gcoords_eq_vectype(gcoords, vec)
    Type(GeneralizedCoords), Intent(inout) :: gcoords
    Type(VecType), Intent(in)  :: vec

    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_gcoords_eq_vectype(gcoords%rigidcoords, vec)
    End If
  End Subroutine gcmodels_gcoords_eq_vectype

  !-------------------------------------------------------------------
  ! Generates the xyz-coordinates in "xyzcoords" from the generalized
  ! coordinates stores in "gcoords".  
  ! Requires: gcoords -- generalized coordinates
  !           xyzcoords -- output xyz coordinate array
  !           optref_defn -- optional reference coordinate definition
  !-------------------------------------------------------------------
  Subroutine gcmodels_toxyz(gcoords, xyzcoords)
    Type(GeneralizedCoords), Intent(In)               :: gcoords
    Type(VecType), Dimension(:), Intent(Out)          :: xyzcoords
    
    !** Call the  routine based on the kind of generalized coordinate model
    If (Associated(gcoords%rigidcoords)) Then
       Call rigidcoords_toxyz(gcoords%rigidcoords, xyzcoords)
     Else If (Associated(gcoords%branchedcoords)) Then
        Call branchedcoords_toxyz(gcoords%branchedcoords, xyzcoords)
     Else
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Could not find assigned generalized coords pointer"
      Stop
    End If
    
  End Subroutine gcmodels_toxyz

  !---------------------------------------------------------------------
  ! Takes the xyz coordinates and generates the generalized coordinates
  ! The vector of center-of-mass coordinates is an optional parameter
  !---------------------------------------------------------------------
  Subroutine gcmodels_fromxyz(gcoords, xyzcoords, spc)
    Type(GeneralizedCoords), Intent(InOut)   :: gcoords
    Type(VecType), Dimension(:), Intent(In)  :: xyzcoords
    Integer, Intent(In)                      :: spc
    
    !** Call the  routine based on the kind of generalized coordinate model
    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_fromxyz(gcoords%rigidcoords, xyzcoords, spc)
    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' fromxyz routine for branched model does not exist yet, sorry'
      Stop
    End If

  End Subroutine gcmodels_fromxyz

  !------------------------------------------------------------------
  ! Gets the Center-of-Mass of the molecule described by generalized
  ! coords "gcoords"
  !-----------------------------------------------------------------
  Type(VecType) Function gcmodels_getcom(gcoords)
    Type(GeneralizedCoords), Intent(inout) :: gcoords
    
    If (Associated(gcoords%rigidcoords)) Then
      gcmodels_getcom = rigidcoords_getcom(gcoords%rigidcoords)
    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' getcom routine for branched model does not exist yet, sorry'
      Stop
    End If
  End Function gcmodels_getcom

  !-----------------------------------------------------------------
  ! This function gets generalized coordinates stored in "gcoords"
  ! and returns them as an array
  !-----------------------------------------------------------------
  Function gcmodels_getgencoords(gcoords)
    Real(kind=RDbl), Dimension(MAX_DOF)  :: gcmodels_getgencoords
    Type(GeneralizedCoords), Intent(In)  :: gcoords

    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_getgencoords(gcoords%rigidcoords,&
          gcmodels_getgencoords(1:6))
    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' getgencoords routine for branched model does not exist yet, sorry'
      Stop
    End If

  End Function gcmodels_getgencoords

  !-----------------------------------------------------------------
  ! This subroutine sets generalized coordinates stored in "gcoords"
  ! from the values stored in the array
  ! Requires: gcoords -- generalized coordinates for one molecule
  !           genarray -- flatted list of reals to use in setting
  !-----------------------------------------------------------------
  Subroutine gcmodels_setgencoords(gcoords,genarray)
    Type(GeneralizedCoords), Intent(InOut)     :: gcoords
    Real(kind=RDbl), Dimension(:), Intent(In)  :: genarray
    If (gcoords%nogcoords) Return
    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_setgencoords(gcoords%rigidcoords,genarray)
    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' set state routine for branched model does not exist yet, sorry'
      Stop
    Else
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " could not find assigned generalized coords pointer"
      Stop
    End If

  End Subroutine gcmodels_setgencoords

  !-------------------------------------------------------------------------
  ! Changes the internal flexibility flag, only works on rigid coordinates
  ! Requires: gcoords -- generalized rigid coordinates for one molecule
  !           flexflag -- new value for internal flexibility
  !-------------------------------------------------------------------------
  Subroutine gcmodels_changeflex(gcoords,flexflag)
    Type(GeneralizedCoords), Intent(InOut)     :: gcoords
    Logical, Intent(In)                        :: flexflag
    If (gcoords%nogcoords) Return
    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_changeflex(gcoords%rigidcoords,flexflag)
    Else If (Associated(gcoords%branchedcoords)) Then
      Return
    Else
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " could not find assigned generalized coords pointer"
      Stop
    End If
    
  End Subroutine gcmodels_changeflex


  !-------------------------------------------------------------------------
  ! Changes the dof, only works on rigid coordinates
  ! Requires: gcoords -- generalized rigid coordinates for one molecule
  !           dof -- new value for internal flexibility
  !-------------------------------------------------------------------------
  Subroutine gcmodels_changedof(gcoords,dof)
    Type(GeneralizedCoords), Intent(InOut)     :: gcoords
    Integer, Intent(In)                        :: dof 

    If (gcoords%nogcoords) Return
    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_changedof(gcoords%rigidcoords,dof)
    Else If (Associated(gcoords%branchedcoords)) Then
      Return
    Else
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " could not find assigned generalized coords pointer"
      Stop
    End If
    
  End Subroutine gcmodels_changedof

  !--------------------------------------------------------
  ! Takes an array of forces in xyz coordinates and returns
  ! the forces in generalized coordinates
  !--------------------------------------------------------
  Subroutine gcmodels_getGenForce(gcoords, xyzforce, genforce)
    Type(GeneralizedCoords), Intent(In)          :: gcoords
    Type(VecType), Dimension(:), Intent(In)      :: xyzforce
    Real(kind=RDbl), Dimension(:), Intent(InOut) :: genforce
    
    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_getGenForce(gcoords%rigidcoords,xyzforce,genforce)
    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' getgenForce routine for branched model does not exist yet, sorry'
      Stop
    End If
  End Subroutine gcmodels_getGenForce

  !-----------------------------------------------------------------------
  ! Sets the generalized coordinates given a set of xyz coordinates
  ! Requires:  gcoords -- generalized rigid coordinates for one molecule
  !            spc -- species number
  !            rp -- set of xyz coordinate vectors for molecule
  !            com -- molecular center of mass
  !-----------------------------------------------------------------------
  Subroutine gcmodels_setfrom_rp(gcoords, spc, rp, com)
    Type(GeneralizedCoords), Intent(In)          :: gcoords
    Integer, Intent(In)                          :: spc
    Type(VecType), Dimension(:), Intent(In)      :: rp 
    Type(VecType), Intent(In)                    :: com

    !** Return immediately if there are no generalized coordinates
    If (gcoords%nogcoords) Return

    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_setfrom_rp(gcoords%rigidcoords,spc,rp,com)

    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          'not prepared to generated branched coordinates from xyz coordinates'
      Stop
    End If

  End Subroutine gcmodels_setfrom_rp

  !-------------------------------------------------------------------
  ! Gets the metric tensor for the transformation from the cartesian
  ! to the generalized coordinates.  
  !------------------------------------------------------------------
  Subroutine gcmodels_getMetricTensor(gcoords, sorb, G, ierr)
    Type(GeneralizedCoords), Intent(inout) :: gcoords
    Integer, Intent(in)      :: sorb
    Real(kind=RDbl), Dimension(:,:), Intent(out) :: G
    Integer, Intent(out)     :: ierr
    
    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_getMetricTensor(gcoords%rigidcoords, sorb, G, ierr)
    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' getMetricTensor routine for branched model does not exist yet, sorry'
      Stop
    End If
  End Subroutine gcmodels_getMetricTensor
  
  !-----------------------------------------------------------------------
  ! This routine does "periodic boundary conditions" on the generalized
  ! coordinates.  E.g. center-of-mass is translated to the simulation cell
  ! and all angles are between 0 and 2pi.
  !-----------------------------------------------------------------------
  Subroutine gcmodels_pbc(gcoords, simcell)
    Type(GeneralizedCoords), Intent(inout) :: gcoords
    Type(Simcell_Params), Intent(in)       :: simcell

    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_pbc(gcoords%rigidcoords, simcell)
    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' pbc routine for branched model does not exist yet, sorry'
      Stop
    End If
  End Subroutine gcmodels_pbc

  !----------------------------------------------------------
  ! Read the generalized coordinates to a restart file
  ! - need_update : a flag telling that the (xyz) coords need to be manually
  !                 updated by the program calling this routine
  !----------------------------------------------------------
  Subroutine gcmodels_readrestartfile(gcoords, sorbtype,unitno,need_update)
    Type(GeneralizedCoords), Intent(inout) :: gcoords
    Integer, Intent(in) :: sorbtype,unitno
    Logical, Intent(out) , Optional :: need_update
    If (Present(need_update)) need_update=.False.
    If (gcoords%nogcoords) Return
    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_readrestartfile(gcoords%rigidcoords, unitno)
    Else If (Associated(gcoords%branchedcoords)) Then
      ! If angles are read from restartfile, then the xyz coords could 
      ! be inacuurate, so they might have to be updated
      If (Present(need_update)) need_update=.True.
     Call branchedcoords_readrestartinfo(gcoords%branchedcoords, sorbtype, &
          unitno)
    End If

  End Subroutine gcmodels_readrestartfile

  !----------------------------------------------------------
  ! Dump the generalized coordinates of one molecule to a restart file
  !----------------------------------------------------------
  Subroutine gcmodels_dumprestartinfo(gcoords, sorbtype, unitno)
    Type(GeneralizedCoords), Intent(in) :: gcoords
    Integer, Intent(in)                 :: sorbtype,unitno
    If (gcoords%nogcoords) Return
    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_dumprestartinfo(gcoords%rigidcoords, unitno)
    Else If (Associated(gcoords%branchedcoords)) Then
      Call branchedcoords_dumprestartinfo(gcoords%branchedcoords, sorbtype, &
          unitno)
    End If

  End Subroutine gcmodels_dumprestartinfo

!!$  !------------------------------------------------------
!!$  !  Displays the generalized coordinates
!!$  !------------------------------------------------------
!!$  Subroutine gcmodels_display(sorbates, unitno)
!!$
!!$    Type(AtMolCoords), Dimension(:), Pointer :: sorbates
!!$    Integer, Intent(In)                      :: unitno
!!$
!!$    Integer                                  :: i
!!$
!!$    Do i = 1, molecules_getnsorbs()
!!$      If (Associated(sorbates(i)%gcoords(1)%branchedcoords)) Then
!!$        Call branchedgc_display(gcoords(1)%branchedcoords, i, unitno)
!!$    End If
!!$ 
!!$  End Subroutine gcmodels_display


  !-------------------------------------------------------------------
  ! Returns the generalized parameters as a string
  !-------------------------------------------------------------------
!!$  Character(len=lstrLen) Function gcmodels_displaystr(gcoords, optformat)
  Function gcmodels_displaystr(gcoords, optformat)
    Character(len=lstrLen) :: gcmodels_displaystr
    Type(GeneralizedCoords), Intent(in) :: gcoords
    Character(*), Optional, Intent(in)  :: optformat
    
    Character(len=strLen)        :: strformat

    If (Present(optformat)) Then
      strformat = optformat
    Else
      strformat = "f8.3"
    End If

    gcmodels_displaystr = " "
    If (Associated(gcoords%rigidcoords)) Then
      gcmodels_displaystr = rigidcoords_displaystr(gcoords%rigidcoords, strformat)
    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' displaystr routine for branched model does not exist yet, sorry'
      Stop
    Endif
  End Function gcmodels_displaystr

  !-------------------------------------------------------------------
  ! Returns the generalized parameters as a string
  !-------------------------------------------------------------------
  Subroutine gcmodels_convstr(gcoords, str)
    Type(GeneralizedCoords), Intent(in) :: gcoords
    Character(*), Intent(in)  :: str

    If (Associated(gcoords%rigidcoords)) Then
      Call rigidcoords_convstr(gcoords%rigidcoords, str)
    Else If (Associated(gcoords%branchedcoords)) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' convstr routine for branched model does not exist yet, sorry'
      Stop
    End If
  End Subroutine gcmodels_convstr

  !----------------------------------------------------------------------
  ! Clean the associated generalized coordinate model
  ! Requires:  gcoords -- generalized rigid coordinates for one molecule
  !----------------------------------------------------------------------
  Subroutine gcmodels_clean(gcoords)
    Type(GeneralizedCoords), Intent(InOut) :: gcoords

    If (Associated(gcoords%rigidcoords)) Then 
      Call rigidcoords_clean(gcoords%rigidcoords)
    Else If (Associated(gcoords%branchedcoords)) Then
      !**doesn't exist yet Call branchedcoords_clean(gcoords%branchedcoords)
    End If

  End Subroutine gcmodels_clean

End Module gcmodels




