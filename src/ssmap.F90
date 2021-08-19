!------------------------------------------------------------------------------
! This module handles interactions between atom, molecules or species
! with an entire FIXED species.  To do this, it uses the pre-tabulated
! map for the fixed species.  The full call sequence for species-species
! interactions is:
!   ssmap_ssint -> ssmap_msint -> ssmap_asint -> maps_int ->-> interpolate
!
! At the moment, three different types of maps are accounted for:
!   old pmap -- The old potential map with a header, specific to a single
!               atom type interacting with a fixed species
!       gmap -- A generic map specific only to the fixed species and the
!               potential type, such as Lennard-Jones 6-12.
!       emap -- The electrostatic map, gives electrostatic potential specific
!               to the fixed species.  This is multiplied by charge to get
!               the potential energy and the forces if necessary
! As usual, these different map types are handled using a set of pointers,
! one of which can be allocated and determines which type will be used.
! 
! Needed Improvements:
! 1) Generic map usage is not complete, need to write more routine in gmap
! 2) Also need a more general method of passing generic map variable parameters
!    and charge if needed
!------------------------------------------------------------------------------

Module ssmap

  Use defaults, Only: RDbl, strLen, lstrLen, MAX_SORBS, scalef, MAX_ATOMTYPES, &
      xlstrLen
  Use utils, Only: split, toupper, allocerrdisplay, deallocerrdisplay, int2str
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_display
  Use simcell, Only: SimCell_Params, simcell_getell
  Use atom, Only: atom_getname, atom_gettypename, atom_getntypes
  Use molecules, Only: molecules_name, molecules_getnatoms, molecules_getatype, &
      molecules_getnatomtypes, molecules_getcharge
  Use config, Only: AtMolCoords, config_getatype, config_getnmoles
  Use maps, Only: Maps_Pointer_Set, maps_init, maps_int, maps_display, &
      maps_clean, maps_nullify, maps_maptype, maps_disp, maps_boxinfo
  Use storebase, Only: EnergyPlus, storebase_disp, storebase_inc
  Use store, Only: Store_Level_Pair, store_idbranch

  Implicit None
  Save

  Private
  Public :: SSMapParams, ssmap_eval, ssmap_ssint, ssmap_msint, ssmap_asint, &
      ssmap_display, ssmap_clean, ssmap_idstring, ssmap_init, &
      ssmap_Oldidstring, ssmap_boxinfo

  Type SSMapParams
    Integer                                        :: spc1,spc2
    Logical                                        :: fast
    Character(len=lstrLen)                         :: line
!    Integer, Dimension(:), Pointer                 :: atype2map
    Type(Maps_Pointer_Set), Dimension(:), Pointer  :: maps
  End Type SSMapParams

  Character(len=lstrLen), Parameter          :: ssmap_idstring = 'MAP'
  Character(len=lstrLen), Parameter          :: ssmap_Oldidstring = 'BASICMAP'

Contains
  !----------------------------------------------------------------------------
  ! Handles the initialization of a map-type interaction calculation between
  ! atoms of one species (spc1) with the other entire species (spc2).  The 
  ! map type and filename for each atom type in the first species is specified
  ! using this notation in the spc-spc interaction file:
  !     atom_name@map_type@filename
  !  for example:
  !     hydrogen@PMAP@sili.hydrogen.pmap
  ! It should also be possible to substitute 'GENERATE' for 'filename'
  ! 
  ! Requires:  params -- SS_Pmap parameters
  !            spc1 -- 1st species number
  !            spc2 -- 2nd species number
  !            simcell -- simulation cell data structure
  !            line -- input initialization info line
  !----------------------------------------------------------------------------
  Subroutine ssmap_init(params,spc1,spc2,simcell,line)
    Type(SSMapParams), Intent(Out)      :: params
    Integer, Intent(In)                 :: spc1,spc2
    Character(*), Intent(In)            :: line
    Type(SimCell_Params), Intent(In)    :: simcell

    Integer                           :: nfields,i,j,n,error,atype
    Integer                           :: natypes,total_atypes
    Integer, Dimension(MAX_ATOMTYPES) :: atypes
    Character(len=lstrLen)            :: newline
    Character(len=strLen), Dimension(100)   :: fields,chunks

    !** Set defaults
    params%spc1 = spc1
    params%spc2 = spc2
    params%fast = .FALSE.
    params%line = line

    !** Allocate the set of maps, so that index corresponds to atom type
    total_atypes = atom_getntypes()
    Allocate(params%maps(total_atypes), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'params%maps')

    !** Get a list of the atom types in this species
    natypes = molecules_getnatomtypes(spc1,atypes)    

    !** Nullify the pointer sets
    Do atype = 1,total_atypes
      Call maps_nullify(params%maps(atype))
    End Do

    !** Split the line into pieces
    nfields = split(line,fields)

    !** Run through the fields and interpret here
    Do i = 5,nfields
      n = split(fields(i),chunks,'@')
      Select Case(Toupper(chunks(1)))
      Case('FAST')
        params%fast = .TRUE.

      Case('ALL')
        !** Initialize the map information for each atom type in species
        Do j = 1,natypes
          Call maps_init(params%maps(atypes(j)),chunks(2:n),simcell)
        End Do

      Case Default     !** assume that it's an atom name
        atype = atom_gettypename(chunks(1))
        If (atype == 0) Then
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              ' Could not find definition for atom name: ',Trim(chunks(1))
          Stop
        End If
!        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!        Write(*,*) 'atom name,type: ',Trim(chunks(1)),'  ',atype

        !** Initialize the actual map information
        Call maps_init(params%maps(atype),chunks(2:n),simcell)

      End Select
    End Do
    IF(.NOT. params%fast) THEN
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) 'Simulation will not run without the fast flag in the sorb_sorb_file!'
        STOP
    END IF

    !** Check to make sure that all the necessary maps were specified
    Do i = 1,natypes
      atype = atypes(i)
      If (maps_maptype(params%maps(atype)) == 'NONE') Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            ' missing map specification for atom type ', &
            Trim(atom_getname(atype))
        Stop      
      End If
    End Do

  End Subroutine ssmap_init

  !----------------------------------------------------------------------------
  ! Handles a single SPECIES-SPECIES interaction.  Assumes that the second
  ! species is the one for while the interpolation map is to be used.
  ! Requires:  ffparams -- spc-spc map forcefield parameters
  !            ffout -- forcefield output 
  !            spc1 -- 1st species number
  !            spc2 -- 2nd species number
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !            pcalcflag -- input flag requesting explicit calculation
  !            coul -- flags desire to calculate coulombic interactions
  !----------------------------------------------------------------------------
  Logical Function ssmap_ssint(ffparams,ffout,spc1,spc2,species,simcell, &
      pcalcflag,coul)
    Type(SSMapParams), Intent(In)               :: ffparams
    Type(Store_Level_Pair), Intent(InOut)       :: ffout
    Integer, Intent(In)                         :: spc1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(SimCell_Params), Intent(In)            :: simcell  
    Logical, Intent(In)                         :: pcalcflag,coul

    Integer                      :: m,level

    ssmap_ssint = .True.
    level = store_idbranch(ffout)

    If (level < 0) Then
      !** no molecule-detail, pass ffout as it is
      Do m = 1,config_getnmoles(species,spc1)
        ssmap_ssint = ssmap_msint(ffparams,ffout,spc1,m,spc2, &
            species,simcell,pcalcflag,coul)
        If (.Not. ssmap_ssint) Return
      End Do

    Else
      !** Pass storage specific to molecule
      Do m = 1,config_getnmoles(species,spc1)
        ssmap_ssint = ssmap_msint(ffparams,ffout%mi(m),spc1,m,spc2, &
            species,simcell,pcalcflag,coul)
        If (.Not. ssmap_ssint) Return
      End Do
    End If

  End Function ssmap_ssint

  !----------------------------------------------------------------------------
  ! Handles a single MOLECULE-SPECIES interaction.  Assumes that the second
  ! species is the one for while the interpolation map is to be used.
  ! Requires:  ffparams -- spc-spc map forcefield parameters
  !            ffout -- forcefield output 
  !            spc1 -- 1st species number
  !            molec1 -- 1st molecule number
  !            spc2 -- 2nd species number
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !            pcalcflag -- input flag requesting explicit calculation
  !            coul -- flags desire to calculate coulombic interactions
  !----------------------------------------------------------------------------
  Logical Function ssmap_msint(ffparams,ffout,spc1,molec1,spc2,species, &
      simcell,pcalcflag,coul)
    Type(SSMapParams), Intent(In)               :: ffparams
    Type(Store_Level_Pair), Intent(InOut)       :: ffout
    Integer, Intent(In)                         :: spc1,molec1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(SimCell_Params), Intent(In)            :: simcell  
    Logical, Intent(In)                         :: pcalcflag,coul

    Integer                      :: a,level
    Real(kind=RDbl)              :: q

    ssmap_msint = .True.

    !** Get the type of the passed storage
    level = store_idbranch(ffout)

    !** Call ssmap_asint with the appropriate part of the storage
    Select Case(level)
    Case(-1)   !** molecule-storage
      If (coul) Then
        Do a = 1,molecules_getnatoms(spc1)
          q = molecules_getcharge(spc1,a)
          ssmap_msint = ssmap_asint(ffparams,ffout%total,spc1,molec1,a, &
              spc2,species,simcell,pcalcflag,q)
          If (.Not. ssmap_msint) Return
        End Do
      Else 
        Do a = 1,molecules_getnatoms(spc1)
          ssmap_msint = ssmap_asint(ffparams,ffout%total,spc1,molec1,a, &
              spc2,species,simcell,pcalcflag)
          If (.Not. ssmap_msint) Return
        End Do
      End If

    Case(0)   !** "compressed" atom-storage
      If (coul) Then
        Do a = 1,molecules_getnatoms(spc1)
          q = molecules_getcharge(spc1,a)
          ssmap_msint = ssmap_asint(ffparams,ffout%binfo(a),spc1,molec1,a, &
              spc2,species,simcell,pcalcflag,q)
          If (.Not. ssmap_msint) Return
          Call storebase_inc(ffout%total,ffout%binfo(a)%nrg)
        End Do
      Else 
        Do a = 1,molecules_getnatoms(spc1)
          ssmap_msint = ssmap_asint(ffparams,ffout%binfo(a),spc1,molec1,a, &
              spc2,species,simcell,pcalcflag)
          If (.Not. ssmap_msint) Return
          Call storebase_inc(ffout%total,ffout%binfo(a)%nrg)
        End Do
      End If

    Case(10)   !** true atom-storage
      If (coul) Then
        Do a = 1,molecules_getnatoms(spc1)
          q = molecules_getcharge(spc1,a)
          ssmap_msint = ssmap_asint(ffparams,ffout%mi(a)%total,spc1,molec1,a, &
              spc2,species,simcell,pcalcflag,q)
          Call storebase_inc(ffout%total,ffout%mi(a)%total%nrg)
          If (.Not. ssmap_msint) Return
        End Do
      Else 
        Do a = 1,molecules_getnatoms(spc1)
          ssmap_msint = ssmap_asint(ffparams,ffout%mi(a)%total,spc1,molec1,a, &
              spc2,species,simcell,pcalcflag)
          Call storebase_inc(ffout%total,ffout%mi(a)%total%nrg)
          If (.Not. ssmap_msint) Return
        End Do
      End If

    Case Default
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' unexpected storage type passed'
      Stop      

    End Select

  End Function ssmap_msint

  !----------------------------------------------------------------------------
  ! Handles a single ATOM-SPECIES interaction.  Assumes that the second
  ! species is the one for while the interpolation map is to be used.
  ! Requires:  ffparams -- spc-spc map forcefield parameters
  !            ffout -- forcefield output 
  !            spc1 -- 1st species number
  !            molec1 -- 1st molecule number
  !            atom1 -- 1st atom number
  !            spc2 -- 2nd species number
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !            fast -- fast/slow interaction flag
  !            pcalcflag -- input flag requesting explicit calculation
  !            q -- optional charge, needed for emap calculations
  !----------------------------------------------------------------------------
  Logical Function ssmap_asint(ffparams,ffout,spc1,molec1,atom1,spc2,species, &
      simcell,pcalcflag,q)
    Type(SSMapParams), Intent(In)               :: ffparams
    Type(EnergyPlus), Intent(InOut)             :: ffout
    Integer, Intent(In)                         :: spc1,molec1,atom1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(SimCell_Params), Intent(In)            :: simcell  
    Logical, Intent(In)                         :: pcalcflag
    Real(kind=RDbl), Intent(In), Optional       :: q

    Integer                            :: atype
    Logical                            :: mapflag
    Real(kind=RDbl)                    :: pot
    Type(VecType)                      :: grad,atvec

    !** Get the atom coordinates and type
    atvec = species(spc1)%coords(atom1,molec1)%r
    atype = molecules_getatype(spc1,atom1)

    !** Interpolate to get interaction(s) at the requested point from the map
    If (Present(q)) Then
      mapflag = maps_int(ffparams%maps(atype),simcell,atvec,ffout,q)
    Else
      mapflag = maps_int(ffparams%maps(atype),simcell,atvec,ffout)
    End If

    !** Do the brute force calculation if allowed and necessary
    If (pcalcflag .And. (.Not. mapflag)) Then
      Write(0,'(3a,i4,a,i4)') "Need to perform brute force pot calc on ", &
          Trim(molecules_name(spc1))," molecule ",molec1," atom ",atom1
      Write(0,'(2a)') "This is not complete yet. Exiting.."
      Stop
    End If

    ssmap_asint = mapflag

  End Function ssmap_asint

  !----------------------------------------------------------------------------
  ! Returns the box increments for a potential map if they are available.
  ! Requires:  ffparams -- spc-spc map forcefield parameters
  !            boxincrements -- x,y,z increments for potential map
  !            boxsteps -- number of boxes in x,y,z directions
  !----------------------------------------------------------------------------
  Logical Function ssmap_boxinfo(ffparams,boxincrements,boxsteps)
    Type(SSMapParams), Intent(In)               :: ffparams
    Real(Kind=RDbl), Dimension(3), Intent(Out)  :: boxincrements
    Integer, Dimension(3), Intent(Out)          :: boxsteps

    Integer    :: i

    ssmap_boxinfo = .False.

    Do i = 1,Size(ffparams%maps)
      ssmap_boxinfo = maps_boxinfo(ffparams%maps(i),boxincrements,boxsteps)
      If (ssmap_boxinfo) Return
    End Do

  End Function ssmap_boxinfo

  !--------------------------------------------------------------------
  ! Determines if the interaction should be evaluated
  ! Requires:  ptr -- pointer to interaction parameters structure
  !            fast -- proposed evaluation speed
  !--------------------------------------------------------------------
  Logical Function ssmap_eval(ptr,fast)
    Type(SSMapParams), Pointer        :: ptr
    Logical, Intent(In)               :: fast

    ssmap_eval = .False.
    If (.Not. Associated(ptr)) Return
    ssmap_eval = (ptr%fast .And. fast)

  End Function ssmap_eval

  !----------------------------------------------------------------------------
  ! Displays the parameters for this interaction type
  ! Requires:  params -- SS_Pmap parameters
  !            indent -- no. of spaces from the left margin
  !            unit -- unit number to dump to
  !----------------------------------------------------------------------------
  Subroutine ssmap_display(params,indent,unit)
    Type(SSMapParams), Intent(In)      :: params
    Integer, Intent(In)                :: indent
    Integer, Intent(In)                :: unit

    Integer                              :: atype
    Character(len=indent)                :: blank
    Character(len=strLen)                :: string
    Character(len=xlstrLen)              :: output
    Character(len=strLen), Dimension(10) :: fields

    blank = Repeat(' ',indent)

    Write(unit,'(2a)') blank,Trim(params%line)

    !** Display information about the maps for each atomic species
    Do atype = 1,Size(params%maps)
      output = maps_disp(params%maps(atype))
      string = int2str(atype)
      Write(unit,'(5a,t24,a)') blank,Trim(atom_getname(atype)), &
          '(',Trim(string),'): ',Trim(output)
    End Do

  End Subroutine ssmap_display

  !----------------------------------------------------------------------------
  ! Cleans up after usage
  ! Requires:  params -- SS_Pmap parameters
  !----------------------------------------------------------------------------
  Subroutine ssmap_clean(params)
    Type(SSMapParams), Intent(InOut)   :: params

    Integer          :: i

    Do i = 1,Size(params%maps)
      Call maps_clean(params%maps(i))
    End Do

  End Subroutine ssmap_clean
  
End Module ssmap
