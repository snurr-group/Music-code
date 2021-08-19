!---------------------------------------------------------------------------
! A module for the creation, maintenance and dumping of densitymaps 
!---------------------------------------------------------------------------

Module dmap

  Use defaults, Only: RDbl, RSgl, strLen, lstrLen, calToJ, MAX_ATOMS
  Use file, Only: file_open
  Use utils, Only: deallocerrdisplay, allocerrdisplay
  Use vector, Only: VecType, Assignment(=)
  Use molecules, Only: molecules_getnatoms, molecules_name
  Use config, Only: AtMolCoords, config_getnmoles, config_getnatoms, &
      config_getnmoleslist, config_getMolecCOM, config_getr
  Use simcell, Only: SimCell_Params, simcell_maptouc
  Use interact, Only: Interaction_Model,interact_boxinfo

  Implicit None
  Save

  Private
  Public :: DensityMap, dmap_idstring, dmap_init, dmap_dump, dmap_update, &
      dmap_clean, dmap_display

  Type DensityMap
    Logical                               :: comflag
    Integer                               :: natoms,write_freq
    Integer, Dimension(MAX_ATOMS)         :: atoms
    Integer, Dimension(3)                 :: subset
    Real(kind=RDbl)                       :: norm
    Integer, Dimension(3)                 :: boxstep
    Real(kind=RDbl), Dimension(3)         :: boxincr
    Character(len=lstrLen)                :: filename
    Real(kind=RSgl), Dimension(:,:,:), Pointer  :: map
  End Type DensityMap

  Character(len=lstrLen), Parameter :: dmap_idstring = 'DENSITY_MAP'

Contains

  !----------------------------------------------------------------------------
  ! Initializes the density map
  ! Requires:  dmap -- density map data structure
  !            write_freq -- frequency of dumping to file
  !            subset -- subset to place in dmap (2nd index need not be defined)
  !            filename -- file name to write into
  !----------------------------------------------------------------------------
  Subroutine dmap_init(dmap,write_freq,subset,filename)
    Type(DensityMap), Intent(Out)       :: dmap
    Integer, Intent(In)                 :: write_freq
    Integer, Dimension(3), Intent(In)   :: subset
    Character(*), Intent(In)            :: filename

    Integer            :: a,m,natoms

    !** Initialize the basics
    dmap%write_freq = write_freq
    dmap%subset = subset
    dmap%filename = filename
    dmap%atoms = 0
    dmap%norm = 0.0_RDbl
    dmap%boxstep = 0
    dmap%boxincr = 0.0_RDbl
    dmap%comflag = .False.
    Nullify(dmap%map)

    !** Determine which atoms should be placed in map
    If (dmap%subset(1) == 0) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' dmap not currently setup for full-system evaluations'
      Stop
    End If
    natoms = molecules_getnatoms(subset(1))
    If (dmap%subset(3) <= 0) Then
      dmap%natoms = natoms
      dmap%comflag = .True.
      Do a = 1,natoms
        dmap%atoms(a) = a
      End Do      

    Else If (dmap%subset(3) > natoms) Then
      dmap%natoms = natoms
      Do a = 1,natoms
        dmap%atoms(a) = a
      End Do

    Else
      dmap%natoms = 1
      dmap%atoms(1) = dmap%subset(3)
    End If

  End Subroutine dmap_init

  !----------------------------------------------------------------------------
  ! Updates the density map, returns true.
  ! Requires:  dmap -- density map data structure
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            imodel -- interaction model information
  !            iter -- iteration number
  !            writetofile -- flag indicating if file should be written to
  !----------------------------------------------------------------------------
  Logical Function dmap_update(dmap,species,simcell,imodel,iter,writetofile)
    Type(DensityMap), Intent(InOut)             :: dmap
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(SimCell_Params), Intent(In)            :: simcell  
    Type(Interaction_Model), Intent(InOut)      :: imodel
    Integer, Intent(In)                         :: iter
    Logical, Intent(Out)                        :: writetofile

    Integer          :: a,m,nmoles,error
    Type(VecType)    :: coord

    !** Set defaults
    writetofile = .False.  !** never want to write to tracking file on return
    dmap_update = .True.

    !** Size the map if this is the first time
    If (.Not. Associated(dmap%map)) Then
      Call interact_boxinfo(imodel,species,dmap%boxincr,dmap%boxstep)

      Allocate(dmap%map(dmap%boxstep(1),dmap%boxstep(2),dmap%boxstep(3)), &
          STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'density map')      
    End If

    !** Update the map
    nmoles = config_getnmoles(species,dmap%subset(1))
    If (dmap%comflag) Then
      Do m = 1,nmoles
        coord = config_getMolecCOM(species,dmap%subset(1),m)
        Call dmap_add(dmap,simcell,coord)
      End Do

    Else 
      Do m = 1,nmoles
        Do a = 1,dmap%natoms
          coord = config_getr(species,(/dmap%subset(1),m,dmap%atoms(a)/))
          Call dmap_add(dmap,simcell,coord)
        End Do
      End Do      
    End If

    !** Write to file if necessary
    If (Mod(iter,dmap%write_freq) == 0) Then
      Call dmap_dump(dmap)
    End If

  End Function dmap_update

  !----------------------------------------------------------------------------
  ! Adds a single new point to the density map
  ! Requires:  dmap -- density map data structure
  !            simcell -- simulation cell information
  !            coord -- coordinate
  !----------------------------------------------------------------------------
  Subroutine dmap_add(dmap,simcell,coord)
    Type(DensityMap), Intent(InOut)      :: dmap
    Type(SimCell_Params), Intent(In)     :: simcell  
    Type(VecType), Intent(InOut)         :: coord

    Integer          :: ipt,jpt,kpt
    Type(VecType)    :: ucvec

    !** Map the coordinates to the unit cell
    ucvec = simcell_maptouc(simcell,coord)

    !** Get the grid point indices
    ipt = int(ucvec%comp(1)/dmap%boxincr(1)) + 1
    jpt = int(ucvec%comp(2)/dmap%boxincr(2)) + 1
    kpt = int(ucvec%comp(3)/dmap%boxincr(3)) + 1
    
    !** put this point in the position density map
    dmap%map(ipt,jpt,kpt) = dmap%map(ipt,jpt,kpt) + 1.0_RSgl

    !** keep track of normalization factor
    dmap%norm = dmap%norm + 1.0_RDbl

  End Subroutine dmap_add

  !----------------------------------------------------------------------------
  ! Normalize and write the density map to a file
  ! Requires:  dmap -- density map data structure
  !----------------------------------------------------------------------------
  Subroutine dmap_dump(dmap)
    Type(DensityMap), Intent(InOut)             :: dmap

    Integer                   :: maptype,unit
    Character(len=25)         :: molecname

    molecname = molecules_name(dmap%subset(1))

    unit = file_open(dmap%filename,102)

    maptype = 2
    Write(unit) molecname,maptype
    Write(unit) dmap%boxstep(1),dmap%boxstep(2),dmap%boxstep(3)
    Write(unit) dmap%map

    Close(unit)

  End Subroutine dmap_dump

  !----------------------------------------------------------------------------
  ! Display the density map (short for now)
  ! Requires:  dmap -- density map data structure
  !            indent -- no. of spaces from the left margin
  !            optunit -- unit number for display
  !----------------------------------------------------------------------------
  Subroutine dmap_display(dmap,indent,unit)
    Type(DensityMap), Intent(In)      :: dmap
    Integer, Intent(In)               :: indent, unit

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'displaying dmap, nothing yet'

  End Subroutine dmap_display

  !----------------------------------------------------------------------------
  ! Display the density map (short for now)
  ! Requires:  dmap -- density map data structure
  !----------------------------------------------------------------------------
  Subroutine dmap_clean(dmap)
    Type(DensityMap), Intent(InOut)             :: dmap

    Integer               :: error

    Deallocate(dmap%map, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)

  End Subroutine dmap_clean

End Module dmap

