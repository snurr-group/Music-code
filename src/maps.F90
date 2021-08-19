!------------------------------------------------------------------------------
! This module a pointer set with pointers to all of the different map types,
! it is the wrapper that selects a map type to initialize, evaluate etc.
!------------------------------------------------------------------------------

Module maps

  Use defaults, Only: RDbl, strLen, xlstrLen, lstrlen
  Use utils, Only: toupper,combine
  Use simcell, Only: SimCell_Params, simcell_getell, simcell_geteff
  Use vector, Only: VecType
  Use pmap, Only: Pmap_Params, pmap_idstring, emap_idstring, pmap_display, &
      pmap_clean, pmap_int, pmap_eint, pmap_init, pmap_disp, pmap_boxinfo
  Use gmap, Only: Generic_Map_Params, gmap_idstring, gmap_int, gmap_display, &
      gmap_clean, gmap_init, gmap_disp, gmap_boxinfo
  Use storebase, Only: EnergyPlus,storebase_nderivs

  Implicit None
  Save

  Private
  Public :: Maps_Pointer_Set, maps_init, maps_int, maps_display, maps_clean, &
      maps_nullify, maps_maptype, maps_disp, maps_boxinfo

  !** The pointer set containing all the map types
  Type Maps_Pointer_Set
    Character(len=xlstrLen)            :: init_string
    Type(Pmap_Params), Pointer         :: pmap
    Type(Pmap_Params), Pointer         :: emap
    Type(Generic_Map_Params), Pointer  :: gmap
  End Type Maps_Pointer_Set

Contains
  !------------------------------------------------------------------
  ! Initializes a particular map using input string.  The string
  ! must contain the idstring for the particular map in the first
  ! position.
  ! Requires:  ptrset -- pointer set of map types to initialize
  !            params -- strings containing initialization information
  !            simcell -- simulation cell information
  !------------------------------------------------------------------
  Subroutine maps_init(ptrset,params,simcell)
    Type(Maps_Pointer_Set), Intent(InOut)           :: ptrset
    Character(len=strLen), Dimension(:), Intent(In) :: params
    Type(SimCell_Params), Intent(In)                :: simcell

    Integer                               :: nfields
    Real(Kind=Rdbl), Dimension(3)         :: eff

    !** Nullify the pointers
    Call maps_nullify(ptrset)

    !** store the initialization string
    ptrset%init_string = combine(params)

    !** Decide which map to initialize and do it
    nfields = Size(params)
    Select Case(Trim(ToUpper(params(1))))
    Case (pmap_idstring)
      eff = simcell_geteff(simcell,.True.)
      Call pmap_init(ptrset%pmap,params(1:nfields),eff)

    Case (emap_idstring)
      eff = simcell_geteff(simcell,.True.)
      Call pmap_init(ptrset%emap,params(1:nfields),eff)

    Case (gmap_idstring)
      Call gmap_init(ptrset%gmap,params(2:nfields),simcell)

    Case Default
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Could not interpret map type designator ",Trim(params(1))
      Stop

    End Select

  End Subroutine maps_init

  !------------------------------------------------------------------
  ! Evaluates an interaction given the pointer set and input info
  ! Requires:  ptrset -- pointer set of map types
  !            simcell -- simulation cell data structure
  !            pt -- position vector for particle
  !            output -- resultant energy plus possible derivatives
  !            q -- optional charge for electrostatic maps
  !------------------------------------------------------------------
  Logical Function maps_int(ptrset,simcell,pt,output,q)
    Type(Maps_Pointer_Set), Intent(In)      :: ptrset
    Type(SimCell_Params), Intent(In)        :: simcell  
    Type(VecType), Intent(In)               :: pt
    Type(EnergyPlus), Intent(InOut)         :: output
    Real(kind=RDbl), Intent(In), Optional   :: q

    !** Call the correct routine to get the interactions
    If (Associated(ptrset%pmap)) Then
      maps_int = pmap_int(ptrset%pmap,simcell,pt,output)

    Else If (Associated(ptrset%emap)) Then
      maps_int = pmap_eint(ptrset%emap,simcell,pt,q,output)

    Else If (Associated(ptrset%gmap)) Then
      maps_int = gmap_int(ptrset%gmap,simcell,pt,output)  !** will need more info

    Else 
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " no map pointer associated, unexpected"
      Stop
    End If

  End Function maps_int

  !------------------------------------------------------------------
  ! Cleanup the pointer set
  ! Requires:  ptrset -- pointer set of map types to nullify
  !------------------------------------------------------------------
  Subroutine maps_nullify(ptrset)
    Type(Maps_Pointer_Set), Intent(InOut)    :: ptrset

    ptrset%init_string = ''
    Nullify(ptrset%pmap)
    Nullify(ptrset%emap)
    Nullify(ptrset%gmap)

  End Subroutine maps_nullify

  !------------------------------------------------------------------
  ! Return a string designating which map type is initialized
  ! Requires:  ptrset -- pointer set of map types 
  !------------------------------------------------------------------
  Function maps_maptype(ptrset)
    Character(len=strLen)                 :: maps_maptype
    Type(Maps_Pointer_Set), Intent(In)    :: ptrset

    If (Associated(ptrset%pmap)) Then
      maps_maptype = 'PMAP'      
    Else If (Associated(ptrset%emap)) Then
      maps_maptype = 'EMAP'      
    Else If (Associated(ptrset%gmap)) Then
      maps_maptype = 'GMAP'      
    Else
      maps_maptype = 'NONE'
    End If

  End Function maps_maptype

  !----------------------------------------------------------------------------
  ! Returns the box increments for a potential map if they are available.
  ! Requires:  ptrset -- pointer set of map types 
  !            boxincrements -- x,y,z increments for potential map
  !            boxsteps -- number of boxes in x,y,z directions
  !----------------------------------------------------------------------------
  Logical Function maps_boxinfo(ptrset,boxincrements,boxsteps)
    Type(Maps_Pointer_Set), Intent(In)          :: ptrset
    Real(Kind=RDbl), Dimension(3), Intent(Out)  :: boxincrements
    Integer, Dimension(3), Intent(Out)          :: boxsteps

    maps_boxinfo = .False.

    If (Associated(ptrset%pmap)) Then
      Call pmap_boxinfo(ptrset%pmap,boxincrements,boxsteps)
      maps_boxinfo = .True.
      Return
    Else If (Associated(ptrset%emap)) Then
      Call pmap_boxinfo(ptrset%emap,boxincrements,boxsteps)
      maps_boxinfo = .True.
      Return
    Else If (Associated(ptrset%gmap)) Then
      Call gmap_boxinfo(ptrset%gmap,boxincrements,boxsteps)
      maps_boxinfo = .True.
      Return
    End If

  End Function maps_boxinfo

  !----------------------------------------------------------------------------
  ! Display the Maps pointer set
  ! Requires:  ptrset -- pointer set of map types to display
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional display unit number
  !----------------------------------------------------------------------------
  Subroutine maps_display(ptrset,indent,unit)
    Type(Maps_Pointer_Set), Intent(In)    :: ptrset
    Integer, Intent(In)                   :: indent
    Integer, Intent(In)                   :: unit

    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)    

    Write(unit,'(3a)') blank,'Map pointer set init info: ', &
        Trim(ptrset%init_string)

    If (Associated(ptrset%pmap)) Then
      Call pmap_display(ptrset%pmap,indent,unit)
    End If

    If (Associated(ptrset%emap)) Then
      Call pmap_display(ptrset%emap,indent,unit)
    End If

    If (Associated(ptrset%gmap)) Then
      Call gmap_display(ptrset%gmap,indent,unit)
    End If

  End Subroutine maps_display

  !----------------------------------------------------------------
  ! Returns a short string giving pertinent information about maps
  ! Requires:  ptrset -- pointer set of map types 
  !----------------------------------------------------------------
  Function maps_disp(ptrset)
    Type(Maps_Pointer_Set), Intent(In)   :: ptrset
    Character(len=xlstrLen)              :: maps_disp

    Character(len=lstrLen)      :: head,tail

    head = 'No map associated'
    tail = ''

    If (Associated(ptrset%pmap)) Then
      head = Trim(pmap_idstring)
      tail = pmap_disp(ptrset%pmap)
    End If

    If (Associated(ptrset%emap)) Then
      head = Trim(emap_idstring)
      tail = pmap_disp(ptrset%emap)
    End If

    If (Associated(ptrset%gmap)) Then
      head = Trim(gmap_idstring)
      tail = gmap_disp(ptrset%gmap)
    End If
    
    Write(maps_disp,'(a,1x,a)') Trim(head),Trim(tail)

  End Function maps_disp

  !------------------------------------------------------------------
  ! Cleanup the pointer set
  ! Requires:  ptrset -- pointer set of map types to deallocate
  !------------------------------------------------------------------
  Subroutine maps_clean(ptrset)
    Type(Maps_Pointer_Set), Intent(InOut)    :: ptrset
    
    !** Clean the pmap if allocated
    If (Associated(ptrset%pmap)) Then
      Call pmap_clean(ptrset%pmap,.True.)
    End If

    !** Clean the emap if allocated
    If (Associated(ptrset%emap)) Then
      Call pmap_clean(ptrset%emap,.True.)
    End If

    !** Clean the gmap if allocated
    If (Associated(ptrset%gmap)) Then
      Call gmap_clean(ptrset%gmap,.True.)
    End If

    Call maps_nullify(ptrset)

  End Subroutine maps_clean

End Module maps

