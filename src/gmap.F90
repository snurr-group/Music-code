!------------------------------------------------------------------------------
! This module handles the generic potential maps, i.e. those that are non-
! electrostatic but still not fixed to a specific atom.  This flexibility
! can be gotten by pretabulating components of a specific potential type,
! such as the Lennard Jones 12-6 potential.  These maps allow us to 
! pre-tabulated the interactions (potential energy and force) between a 
! point particle and an entire FIXED species.  They are in the same module 
! because they are similar enough to use some of the same routines.  Note 
! that the storage for these maps is also in this module to avoid duplication; 
! pointers to this storage can be found in other data structures.
!
! Needed Improvements
! 1) probably should move storage to this module, as indicated above
! 2) need to swtich to using Map_Core_Info
!------------------------------------------------------------------------------

Module gmap

  Use defaults, Only: RDbl, strLen, RSgl, MAX_MAPS, lstrlen
  Use utils, Only: split, getpath, isfileopen, allocerrdisplay, &
      deallocerrdisplay, int2str
  Use file, Only: file_getunit, file_open
  Use vector, Only: VecType
  Use simcell, Only: SimCell_Params, simcell_getell
  Use mapheader, Only: mapheader_getheader, mapheader_cleanup, Header_Info, &
      mapheader_getIndices, mapheader_writeheader, mapheader_boxinfo
  Use interpolate, Only: Map_Core_Info, interpolate_int, interpolate_clean, &
      interpolate_initcore
  Use storebase, Only: EnergyPlus,storebase_nderivs

  Implicit None
  Save

  Private
  Public :: Generic_Map_Params, gmap_init, gmap_int, gmap_clean, & 
      gmap_display, gmap_idstring, gmap_disp, gmap_boxinfo

  !** Generic Maps, specifically for LJ potential type
  Type Generic_Map_Params
    Character(len=strLen)                  :: filename
    Type(Header_Info)                      :: header
    Integer                                :: npoints
    Type(Map_Core_Info)                    :: core12,core6
  End Type Generic_Map_Params

  Character(len=strLen), Parameter         :: gmap_idstring = 'GMAP'

Contains

  !----------------------------------------------------------------------------
  ! Initialize a generic map
  ! Requires:  map -- Generic map type to initialize
  !            fields -- array of strings containing initialization info
  !            simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Subroutine gmap_init(map,fields,simcell)
    Type(Generic_Map_Params), Intent(InOut)         :: map
    Character(len=strLen), Dimension(:), Intent(In) :: fields    
    Type(SimCell_Params), Intent(In)                :: simcell

    Integer                  :: i,nfields

    nfields = Size(fields)

    Write(0,'(2a,i6,a,3i4)') __FILE__,":",__LINE__, &
        ' Generic maps not yet available, sorry'
    Stop

    !** Interpret the fields and act accordingly
    Do i = 1,nfields
      Select Case(Trim(fields(i)))
      Case ('GENERATE')
        Write(0,'(2a,i6,a,3i4)') __FILE__,":",__LINE__, &
            ' GENERATE option not yet available for Generic maps'
        Stop
        
      Case Default
        Write(0,'(2a,i6,a,3i4)') __FILE__,":",__LINE__, &
            ' Could not understand field: ',Trim(fields(i))
        Stop

      End Select

    End Do

  End Subroutine gmap_init

  !----------------------------------------------------------------------------
  ! Get the interaction from a potential map at a specified point
  ! Requires:  pmap -- Potential map 
  !            simcell -- simulation cell data structure
  !            pt -- position vector for particle
  !            output -- resultant energy plus possible derivatives
  !----------------------------------------------------------------------------
  Logical Function gmap_int(map,simcell,pt,output)
    Type(Generic_Map_Params), Intent(In)  :: map
    Type(SimCell_Params), Intent(In)      :: simcell  
    Type(VecType), Intent(In)             :: pt
    Type(EnergyPlus), Intent(InOut)       :: output

    !** Get the interactions by interpolating from map
!    gmap_int = interpolate_int(map%core,map%header,output,pt,simcell)
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Generic Map evaluations not yet available'
    output%nrg = output%nrg + 0.0_Rdbl
    gmap_int = .False.

  End Function gmap_int

  !----------------------------------------------------------------------------
  ! Display the generic map structure
  ! Requires:  map -- Generic map type to display
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional display unit number
  !----------------------------------------------------------------------------
  Subroutine gmap_display(map,indent,unit)
    Type(Generic_Map_Params), Intent(In)  :: map
    Integer, Intent(In)                   :: indent
    Integer, Intent(In)                   :: unit

    Character(len=indent)   :: blank

    blank = Repeat(' ',indent)

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'gmap display routine not yet available'
    Stop

  End Subroutine gmap_display

  !----------------------------------------------------------------------------
  ! Returns the box increments for a potential map if they are available.
  ! Requires:  map -- potential map 
  !            boxincrements -- x,y,z increments for potential map
  !            boxsteps -- number of boxes in x,y,z directions
  !----------------------------------------------------------------------------
  Subroutine gmap_boxinfo(map,boxincrements,boxsteps)
    Type(Generic_Map_Params), Intent(In)        :: map
    Real(Kind=RDbl), Dimension(3), Intent(Out)  :: boxincrements
    Integer, Dimension(3), Intent(Out)          :: boxsteps

    Call mapheader_boxinfo(map%header,boxincrements,boxsteps)

  End Subroutine gmap_boxinfo

  !----------------------------------------------------------------
  ! Returns a short string giving pertinent information about map
  ! Requires:  map -- generic map for which to return info
  !----------------------------------------------------------------
  Function gmap_disp(map)
    Type(Generic_Map_Params), Intent(In) :: map
    Character(len=lstrLen)               :: gmap_disp

    Character(len=strLen)      :: string
    
    string = int2str(map%npoints)
    Write(gmap_disp,'(2a,3x,2a)') 'map from: ',Trim(map%filename), &
        'tabulated points: ',Trim(string)

  End Function gmap_disp

  !----------------------------------------------------------------------------
  ! Clean the generic map structure
  ! Requires:  map -- Generic map type to deallocate
  !            headertoo -- flag indicating if it should clean the header also
  !----------------------------------------------------------------------------
  Subroutine gmap_clean(map,headertoo)
    Type(Generic_Map_Params), Intent(InOut)   :: map
    Logical, Intent(In)                       :: headertoo

    Integer      :: error

    !** Clean the header if desired
    If (headertoo) Then
      Call mapheader_cleanup(map%header)
    End If

    !** Clean the cores
    Call interpolate_clean(map%core12)
    Call interpolate_clean(map%core6)

  End Subroutine gmap_clean

End Module gmap
