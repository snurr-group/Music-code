!------------------------------------------------------------------------------
! Used to make any type of map that you need. It needs to be higher than
! forcefield in order to make the maps. Otherwise, this would fit neatly into
! maps or ssmap or something like that.
!------------------------------------------------------------------------------
Module mapmaker

  Use defaults, Only: RDbl, strLen, dashedline, lstrLen, scalepe, Rsgl, &
      MAX_ATOMTYPES
  Use simcell, Only: SimCell_Params, simcell_geteff, simcell_getmolectype
  Use config, Only: AtMolCoords, config_makegrid, config_getnmoles, &
      config_getr, config_isfixed
  Use vector, Only: VecType, vector_display, Assignment(=)
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toreal, &
      toupper, findint
  Use file, Only: file_getunit, file_open
  Use molecules, Only: molecules_gettype, molecules_name, molecules_getnatoms,&
      molecules_getnthatom, molecules_getcharge, molecules_getnatomtypes, &
      molecules_getatype
  Use subinteract, Only: Subset_Interactions, subinteract_init, &
      subinteract_newnrg, subinteract_int
  Use interact, Only: Interaction_Model, interact_simpleint, &
      interact_getpotparameters
  Use gmap, Only: Generic_Map_Params
  Use general, Only: genparams
  Use pmap, Only: Pmap_Params, pmap_addpoint, pmap_newmap, pmap_write, &
      pmap_compact
  Use atom, Only: atom_getname, atom_getsymbol, atom_getntypes
  Use wtime, Only: wtime_init, wtime_starttime, wtime_stoptime, wtime_display
  Use storebase, Only: EnergyPlus, storebase_init

  Implicit None

  Private
  Public :: MapInfo, MapMakerParams, mapmaker_idstring, mapmaker_init, &
      mapmaker_display, mapmaker_makemap

  Type MapInfo
    Integer :: nmaps
    Type(MapMakerParams), Dimension(:), Pointer :: maps
    !** Subset interaction storage, species-dependent, needed for nrg calcs
  End Type MapInfo

  Type MapMakerParams
    Type(AtMolCoords)         :: grid
    Real(Kind=RDbl)           :: gridspace
    Integer                   :: mapsorb
    Integer                   :: probe
    Character(len=strLen)     :: filename, int, model
    Integer, Dimension(3)     :: npnts
    Real(Kind=RDbl)           :: hipotcut  ! Stored in kJ/mol
    Type(Pmap_Params)         :: map
    Type(Subset_Interactions) :: subint
  End Type MapMakerParams

  Character(len=strLen), Parameter :: mapmaker_idstring="Mapmaker Information"
  Real(kind=RDbl) :: VERY_LRG_NRG=1.00E10_RDbl

!!$L
    Type(EnergyPlus) :: ffout
  
    Real(kind=RDBl),Dimension(:,:,:),Allocatable,SAVE :: map_pot_params
Contains

  !----------------------------------------------------------------------------
  ! Reads from the control file the parameters for making a new map. 
  ! Also makes up the probe sorbate grid and stores it in sorbates
  !----------------------------------------------------------------------------
  Subroutine mapmaker_init(mmparams,imodel, simcell,sorbates,ctrlFile,optTag)
    Type(MapInfo), Intent(InOut)     :: mmparams
    Type(Interaction_Model), Intent(In), Target  :: imodel
    Type(SimCell_Params), Intent(In) :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: sorbates
    Character(*), Intent(In) :: ctrlFile
    Character(*), Optional, Intent(In) :: optTag
    Integer :: nomaps, lineno, unitno, error, nmaps, i, nfields, tempunit
    Character(len=strLen), Dimension(strLen) :: fields
    Real(Kind=RDbl), Dimension(3) :: eff
    Character(len=strLen) :: tag
    Character(len=lstrLen) :: line
    Integer :: j, unitn, probesorb, probe, probeAtomType, nprobeatoms
    Integer, Dimension(3) :: indices

    Integer                                 :: sorb1, sorb2, natypes1, natypes2
    Integer                                 :: a1, a2, natype
    Integer, Dimension(max_atomtypes)       :: atype1, atype2
    Real(kind=RDbl),Dimension(6)            :: pot_params


    !** Check to see if there is a custom tag we should search for
    If (Present(optTag)) Then
      tag = optTag
    Else
      tag = mapmaker_idstring
    End If

    !** Open the ctrlFile if it is not opened
    unitno = isfileopen(ctrlFile)
    If (unitno < 0) Then
      unitno = file_getunit(ctrlFile)
      Open(file=ctrlFile, unit=unitno)
    Endif

    !** Find the tag in the control file and begin reading at that
    !** point
    lineno = filesrchstr(unitno,tag,line,.True.)
    If (lineno == 0) Then
      Write(0,'(2a,i3,4a)') __FILE__,":",__LINE__, &
          " Could not find string ",Trim(mapmaker_idstring), &
          " in the control file ",trim(ctrlFile)
      Write(0,'(a)') " Expected format for the section in the control file: "
      Write(0,'(a)') trim(mapmaker_idstring)
      Write(0,'(a)') "(Int)    # no. of maps"
      Write(0,'(a)') "(blank line)"
      Write(0,'(a)') "(Char)   # Map sorbate type (must be a listed sorbate)"
      Write(0,'(a)') "(Char)   # Probe sorbate type (must be a listed sorbate)"
      Write(0,'(a)') "(Char)   # Interaction type (COUL, NCOUL, LJ, etc)"
      Write(0,'(a)') "(Real)   # Grid spacing, Ang"
      Write(0,'(a)') "(Real)   # High potential cutoff, kcal"
      Write(0,'(a)') "(Char)   # Map base filename"
      Stop
    End If

    !** Read how many maps we have to make
    Read(unitno,*) mmparams%nmaps

    !** Allocate mmparams
    Allocate(mmparams%maps(mmparams%nmaps),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i3,a)') __FILE__,":",__LINE__, &
          " Could not allocate mmparams%map."
      Stop
    End If

    !** Loop over the number of maps and read in the parameters.
    Do i = 1, mmparams%nmaps

      !MDEBUG 
      Write(*,*) "Reading map ",i

      !** Read the blank line
      Read(unitno,'(a)') line

      !** Read in the sorbate to map
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields)
      mmparams%maps(i)%mapsorb = molecules_gettype(trim(fields(1)))

      !** Check to make sure the molecule is both fixed
      !** and associated with the simulation cell
      If (.Not.(config_isfixed(sorbates(mmparams%maps(i)%mapsorb)))) Then
        Write(0,'(2a,i5,3a)') __FILE__,":",__LINE__, &
            " Molecule to map ", &
            Trim(molecules_name(mmparams%maps(i)%mapsorb)),&
            " must have a FIXED configuration"
        Stop
      End If

      If (.Not.(mmparams%maps(i)%mapsorb==simcell_getmolectype(simcell))) Then
        Write(0,'(2a,i5,3a)') __FILE__,":",__LINE__, &
            " Molecule to map ", &
            Trim(molecules_name(mmparams%maps(i)%mapsorb)),&
            " must be associated with the simulation cell"
        Stop
      End If

      !** Read in the probe molecule
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields)
      mmparams%maps(i)%probe = molecules_gettype(trim(fields(1)))

      ! These checks are good so that later there will not be any problem
      probe=mmparams%maps(i)%probe
      If (.Not.(probe==1)) Then
        Write(0,'(2a,i5,3a)') __FILE__,":",__LINE__, &
            " Probe to map : ", &
            Trim(fields(1))," must be the first sorbate; "//&
            "otherwise map generation will not work"
        Stop
      End If
      nprobeatoms=molecules_getnatoms(probe)
      If (nprobeatoms/=1) Then
        Write(0,'(2a,i5,3a)') __FILE__,":",__LINE__, &
            " Probe to map : ", &
            Trim(fields(1))," Should have only one atom. Now it has ",&
            nprobeatoms
        Stop
      Endif
      probeAtomType=molecules_getatype(probe,1)
      If (probeAtomType/=1) Then
        Write(0,'(2a,i5,2a)') __FILE__,":",__LINE__, &
            " Probe-atom to map  Should be the first in the atom list "
        Stop
      Endif

      !** Read in the interaction type to map
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields)
      mmparams%maps(i)%int = trim(fields(1))
      mmparams%maps(i)%model = trim(fields(2))

      !** Read in the grid spacing in Ang
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields)
      mmparams%maps(i)%gridspace = toreal(fields(1))

      !** Read in the high potential cutoff
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields)
      mmparams%maps(i)%hipotcut = toreal(fields(1))

      !** Read in the base filename
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields)
      mmparams%maps(i)%filename = trim(fields(1))

      !** Make the filename
      If (toupper(Trim(mmparams%maps(i)%filename)) == "AUTO") Then
        mmparams%maps(i)%filename = mapmaker_makeName(mmparams%maps(i))
      End If

      !check whether the file already exists, its a pain to rewrite 
      !over an existing mapfile by mistake
      Write(*,*) " pmap generated will have name : "//&
          Trim(mmparams%maps(i)%filename)
      Write(*,*) "Checking for existing pmap file with same name :- "
      Write(*,*) "It already exists, Please change the pmap-file name &
          & in ctrlfile"
      tempunit=file_open(mmparams%maps(i)%filename,101)

    End Do
    !
    ! Tina added
    !
    natype = atom_getntypes()
    Allocate(map_pot_params(1,natype,6))

    Do i = 1, mmparams%nmaps
      !** Keep the user up to speed
      Write(0,'(5a)') "Creating grid for ", &
          Trim(molecules_name(mmparams%maps(i)%mapsorb)), &
          "-",Trim(molecules_name(mmparams%maps(i)%probe)), &
          " map interaction"
      !** Create the grid that will be used. 
      !** Call config to construct the grid points
      Call config_makegrid(simcell,mmparams%maps(i)%probe, &
          mmparams%maps(i)%gridspace, &
          mmparams%maps(i)%grid,mmparams%maps(i)%npnts,.True.)
!!$      !MDEBUG
!!$      Write(0,*) __FILE__,__LINE__," dumping points to grid.xyz"
!!$      unitn = file_open('grid.xyz',111)
!!$      Do j = 1, mmparams%maps(i)%npnts(1)*mmparams%maps(i)%npnts(2)*mmparams%maps(i)%npnts(3)
!!$        Write(unitn,'(i6,2x,3f15.6)') j, &
!!$            mmparams%maps(i)%grid%coords(1,j)%r%comp(1), &
!!$            mmparams%maps(i)%grid%coords(1,j)%r%comp(2), &
!!$            mmparams%maps(i)%grid%coords(1,j)%r%comp(3)
!!$        Call mapheader_getIndices(mmparams%maps(i)%grid%coords(1,j)%r%comp(1),&
!!$            mmparams%maps(i)%grid%coords(1,j)%r%comp(2), &
!!$            mmparams%maps(i)%grid%coords(1,j)%r%comp(3),indices)
!!$        Write(unitn,'(i6,2x,3i5)') j,indices(3),indices(2),indices(1)
!!$      End Do
!!$      Close(unitn)
!!$      stop
      !
      ! Tina added
      ! Initialise some of the parameters to be written to mapheader
      !
      sorb1 = mmparams%maps(i)%probe
      sorb2 = mmparams%maps(i)%mapsorb
      atype1 = 0
      atype2 = 0
      natypes1 = molecules_getnatomtypes(sorb1,atype1)
      natypes2 = molecules_getnatomtypes(sorb2,atype2)
      natype = atom_getntypes()

      Do a1 = 1, natypes1
        If (findint(atype1(1:natypes1),a1) == 0) Cycle
        Do a2 = natypes1+1, natype
          If (findint(atype2(1:natypes2),a2) == 0) Cycle
          Call interact_getpotparameters(imodel,sorb1,sorb2,a1,a2,pot_params)
          map_pot_params(1,a2,:) = pot_params
        End do
      End Do
    End Do

    !** Check to make sure that we have enough unit cells
    !** to do the job. We need enough such that any edge length
    !** is twice as long as the potential cutoff
    !    cutoff = forcefield_cutoff(probe,map,.True.)
    !** Well, we'll do this once forcefield is worked out such
    !** that it can return the cutoff for the potential used.
    Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
        " WARNING: the number of UCs required has not been checked"

    If (size(mmparams%maps)>1) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
      !SDEBUG
    endif
    probesorb=mmparams%maps(1)%probe 
    !** initialize the subset interaction structure
    mmparams%maps(1)%subint%imodel=>imodel
!!$L   Call subinteract_init(mmparams%maps(1)%subint,imodel,'Molec_System', &
!!$L        'MC',(/0,0,0/),(/0,0,0/))

  End Subroutine mapmaker_init

  !----------------------------------------------------------------------------
  ! Hacks the information about the map into the map header. Here, we directly
  ! access the fields of the header variable, hence the hack name.
  !----------------------------------------------------------------------------
  Subroutine mapmaker_hackheader(mmparams,simcell,npoints)
    Type(MapMakerParams), Intent(InOut) :: mmparams
    Type(SimCell_Params), Intent(In)    :: simcell
    Real(Kind=RDbl), Dimension(3) :: steps, npnts
    Real(Kind=RDbl), Dimension(3) :: eff
    Integer, Intent(In), Optional :: npoints
    Integer, Dimension(MAX_ATOMTYPES) :: alist

    Integer :: natoms, atype, atom, i, error

    !** The number of points
     mmparams%map%header%nbrx = mmparams%npnts(1)
     mmparams%map%header%nbry = mmparams%npnts(2)
     mmparams%map%header%nbrz = mmparams%npnts(3)

    !** The number of cublets
     mmparams%map%header%ncubelets=(mmparams%npnts(1)-1)* &
        (mmparams%npnts(2)-1)*(mmparams%npnts(3)-1)

    !** This is the number of TABULATED points. If passed,
    !** set it. Otherwise it is 0 for now.
     If (Present(npoints)) Then
      mmparams%map%header%nfsize = npoints
     Else
      mmparams%map%header%nfsize = 0
     End If

    !** Get the effective edge lengths of the simulation cell
     eff = simcell_geteff(simcell,.True.)
     mmparams%map%header%eff = eff

    !** The step size
     npnts = Int(eff/mmparams%gridspace) + 1
     steps = eff/Real(npnts-1,Kind=RDbl)

     mmparams%map%header%stepx = steps(1)
     mmparams%map%header%stepy = steps(2)
     mmparams%map%header%stepz = steps(3)
 
    !** High potential cutoff
     mmparams%map%header%pcut_off = mmparams%hipotcut

    !** Map and probe molecules
     mmparams%map%header%mapMolecule = molecules_name(mmparams%mapsorb)
     mmparams%map%header%probeMolecule = molecules_name(mmparams%probe)

    !** Store the map molecule atoms, charges, names, and symbols
     mmparams%map%header%nzatoms = &
        molecules_getnatomtypes(mmparams%mapsorb,alist)
     Allocate(mmparams%map%header%zatom(mmparams%map%header%nzatoms),stat=error)
     If (error /= 0) Then
       Write(0,'(2a,i5,a,i5)') __FILE__,":",__LINE__, &
          " Could not allocate header%zatoms of size ", &
          mmparams%map%header%nzatoms
      Stop
     End If
#ifdef DEBUG
Write(*,*) "SOME PMAP INFO"
    Write(*,*) Size(map_pot_params,1)
    Write(*,*) Size(map_pot_params,2)
    Write(*,*) Size(map_pot_params,3)
    Write(*,*) map_pot_params(1,1,:)
    Write(*,*) map_pot_params(1,2,:)
    Write(*,*) map_pot_params(1,3,:)
    Write(*,*) map_pot_params(1,4,:)
    Write(*,*) map_pot_params(1,5,:)
    Write(*,*) map_pot_params(1,6,:)
#endif
     Do i = 1, mmparams%map%header%nzatoms
      atype = alist(i)
      mmparams%map%header%zatom(i)%atom_name = atom_getname(atype)
      Write(*,*) atype, mmparams%map%header%zatom(i)%atom_name 
      mmparams%map%header%zatom(i)%symbol = atom_getsymbol(atype)
      mmparams%map%header%zatom(i)%charge = &
          molecules_getcharge(mmparams%mapSorb,i)      
      mmparams%map%header%zatom(i)%A12 = map_pot_params(1,atype,1)
      mmparams%map%header%zatom(i)%B6 = map_pot_params(1,atype,2)
      mmparams%map%header%zatom(i)%C2 = map_pot_params(1,atype,3)
      mmparams%map%header%zatom(i)%D0 = map_pot_params(1,atype,4)
      mmparams%map%header%zatom(i)%rcuthi = map_pot_params(1,atype,5)
      mmparams%map%header%zatom(i)%rcutlo = map_pot_params(1,atype,6)
     End Do
  
  End Subroutine mapmaker_hackheader

  !----------------------------------------------------------------------------
  ! Creates a copy of mmparams in cpparams EXECPT for the map itself
  !----------------------------------------------------------------------------
  Subroutine mapmaker_copy(mmparams,cpparams)
    Type(MapMakerParams), Intent(In) :: mmparams
    Type(MapMakerParams), Intent(Out):: cpparams

    cpparams%npnts = mmparams%npnts
    cpparams%hipotcut = mmparams%hipotcut
    cpparams%mapsorb = mmparams%mapsorb
    cpparams%probe = mmparams%probe
    cpparams%gridspace = mmparams%gridspace

  End Subroutine mapmaker_copy

  !----------------------------------------------------------------------------
  ! Calls the correct routines to make the map. The sorbate we are mapping
  ! is mapsorb, and it is probed using probesorb.
  !----------------------------------------------------------------------------
  Subroutine mapmaker_makemap(mmparams,simcell,sorbates)
    Type(MapMakerParams), Intent(InOut) :: mmparams
    Type(SimCell_Params), Intent(In)    :: simcell
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: sorbates

    Integer         :: npmoles, nmmoles, npatoms, nmatoms, m, natoms
    Integer         :: mapsorb, probesorb,ap, molec, npts, i, index
    Logical         :: mapflag, fast, accept, skip_intra_flag, succ_flag, addpoint
    Type(VecType)   :: pcoord, r
    Real(Kind=RDbl) :: pot,cpot, max_nrg
    Character(len=strLen)       :: filename
    Real(Kind=RDbl), Dimension(1:8) :: hot
    Type(MapMakerParams)        :: tempMap

!!$L
    !**Initialize a structure to store energy for each point
    Call storebase_init(ffout,2)

    !** Start the timer
    Call wtime_init(1)
    Call wtime_starttime(1)
    Write(0,'(5a)') "Making map for ",Trim(molecules_name(mmparams%mapsorb)), &
        "-",Trim(molecules_name(mmparams%probe))," map interaction"

    !** WARNING this variable should be set elsewhere
    Write(0,'(a)') "WARNING: fast flag is not currently passed to makemap"

    !** The map will be constructed on a coord by coord basis. In
    !** sorbates, the probesorb is a collection of grid points that 
    !** represent the space surrounding the sorbate to be mapped,
    !** mapsorb.

    !** Get the mapsorb and probesorb
    mapsorb = mmparams%mapsorb
    probesorb = mmparams%probe

    !** get the number of atoms and molecules
    npatoms = molecules_getnatoms(probesorb)

    !** Remember, npmoles is the number of grid points for the map
    npmoles = config_getnmoles(mmparams%grid)
    nmatoms = molecules_getnatoms(mapsorb)
    nmmoles = config_getnmoles(sorbates(mapsorb))

    !** Copy the parameters
    Call mapmaker_copy(mmparams,tempMap)

    !** Initialize the map variables
    Call pmap_newmap(tempMap%map,npmoles,mmparams%npnts)
    Call pmap_newmap(mmparams%map,npmoles,mmparams%npnts)

    !** Hack in the header
    Call mapmaker_hackheader(tempMap,simcell)

!!$    !MDEBUG
!!$    npmoles = 10
!!$    natoms = 1

    max_nrg=VERY_LRG_NRG


    Do m = 1, npmoles

      !** Update the user as to the progess made
      If (Mod(m,1000) == 0) Then
        Write(0,'(a,i9,a,f6.2,a)') "Map progress: ",m," points calculated", &
            (m*100.0)/npmoles,"% of total"
      End If

      pot = 0.0_RDbl
      hot = 0.0_RDbl

      !** We do 1 point at a time. 
      sorbates(probesorb)%coords(1,1)%r = mmparams%grid%coords(1,m)%r
      sorbates(probesorb)%coords(1,1)%rp = mmparams%grid%coords(1,m)%r

      !** coord to probe at; probe should be single atom
      pcoord = config_getr(sorbates,(/probesorb,1,1/))

!!$      !** Write out the coord for debugging
!!$      Write(*,*) "pcoord : ",m,pcoord
      skip_intra_flag=.True.
      fast=.True.
      molec=1
      addpoint=.True.

      succ_flag = interact_simpleint(mmparams%subint%imodel,ffout,(/probesorb,1,1/), &
          (/mapsorb,0,0/),sorbates,simcell,fast,hot)
      pot=ffout%nrg

      !	print*, "Storebase: ", ffout%nrg, ffout%force(1), ffout%force(2), ffout%force(3), &
      !                                          ffout%hess(1), ffout%hess(2), ffout%hess(3), ffout%hess(4)
      !
      !  	print*, "Hots:       ", hot(1), hot(2), hot(3), hot(4), hot(5), hot(6), hot(7), hot(8)


!!$L      succ_flag = subinteract_int(mmparams%subint,sorbates, simcell,fast, &
!!$L          .True.,skip_intra_flag,(/max_nrg/),(/0,0,0/))
!!$L      pot=subinteract_newnrg(mmparams%subint,.True.)
!!$
!!$      !** Get the potential and forces at that point
!!$      Call forcefield_getmsint(sorbates,simcell,probesorb,1,.True.,pot, &
!!$          mapflag,hot)

!!$      Write(*,'(a,f15.3,a)') "Map potential : ",pot," kcal"
!!$
!!$      !MDEBUG
!!$      Write(*,*) "hot1 : ",vector_display(hot(1),'e14.3')
!!$      Write(*,*) "hot2 : ",vector_display(hot(2),'e14.3')
!!$      Write(*,*) "hot3 : ",vector_display(hot(3),'e14.3')
!!$

      !** no need to use this value ( most probably energy went very high
      !** so no need to tabulate
      If (.Not.(succ_flag)) Then
!!$L       Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$L       Write(*,*) " most probably the probe atoms is getting too close to "
!!$L       Write(*,*) " the map-sorb atoms. Check cutoffs in atom-atom file. "
!!$L        Write(*,*) "Use lower cutoff If necessary"
        pot = max_nrg
        hot(1) = max_nrg
      Endif
      !SDEBUG
!!$      If (pot==0.0_RDbl) Then
!!$        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$        Write(*,*) pot, tempMap%map%npoints+1, mapflag
!!$        Write(*,*) "Ar",pcoord
!!$      Endif
      !SDEBUG      

      !** Store them in the map
      Call pmap_addpoint(tempMap%map,pcoord,pot,hot,.True.)

    End Do

    !MDEBUG
    Write(0,*) __FILE__,__LINE__,"total points : ",tempMap%map%npoints, &
        mmparams%hipotcut

    !** Now we need to compact the map into the final map
    Call pmap_compact(tempMap%map,mmparams%map,mmparams%hipotcut/scalepe)

    !MDEBUG
    Write(0,*) __FILE__,__LINE__,"compact points : ",mmparams%map%npoints

    !** Hack in the header
    Call mapmaker_hackheader(mmparams,simcell,mmparams%map%npoints)

    If (Trim(genparams%displaymode) == "VERBOSE") Then
      npts=mmparams%map%core%npts
      Write(*,*) " Displaying the 5 points in the map NRG and Forces"
      Write(*,*) "    Index         NRG        F(x)        F(x)        F(z)"
      Do i=1,5
        index=Int((i*npts*1.00))/5
        Write(*,'(i10,4e12.4)') index,mmparams%map%core%r(index), &
            mmparams%map%core%rx(index), mmparams%map%core%ry(index), &
            mmparams%map%core%rz(index)
      End Do
      Write(*,*)
      Write(*,*) " Displaying the 5 points in the map 2nd and 3rd derivs"
      Write(*,*) "    Index         Vxy         Vxz         Vyz        Vxyz"
      Do i=1,5
        index=Int((i*npts*1.00))/5
        Write(*,'(i10,4e12.4)') index,mmparams%map%core%rxy(index), &
            mmparams%map%core%rxz(index), mmparams%map%core%ryz(index), &
            mmparams%map%core%rxyz(index)
      End Do
    Endif

    !** Write the new map to the file
    Call pmap_write(mmparams%map,mmparams%filename)

    !** Stop the time and report the time
    Call wtime_stoptime(1)
    Call wtime_display('Map creation time',1)

  End Subroutine mapmaker_makemap

  !----------------------------------------------------------------------------
  ! Creates the map filename based on the mapSorb name, probeSorb name,
  ! and interaction type
  !----------------------------------------------------------------------------
!!$  Character(len=strLen) Function mapmaker_makeName(mmparams)
  Function mapmaker_makeName(mmparams)
    Character(len=strLen) :: mapmaker_makeName
    Type(MapMakerParams), Intent(In) :: mmparams
    Character(len=strLen) :: filename

    !** Come up with the new filename
    filename = Trim(molecules_name(mmparams%mapSorb))

    !** If this is not a probe molecule, add the atom name
    If (toupper(Trim(molecules_name(mmparams%probe))) /= "PROBE") Then
      filename=trim(filename)//"."//Trim(molecules_name(mmparams%probe))
    End If

    !** Add the interaction type we are mapping to the name
    filename = trim(filename)//"."//Trim(mmparams%model)//".map"

    !** Set the return variable
    mapmaker_makeName = filename

  End Function mapmaker_makeName

  !----------------------------------------------------------------------------
  ! Compacts the map by saving only those points below the energy cutoff
  ! and those that are needed for the hermite interpolation
  !----------------------------------------------------------------------------
  Subroutine mapmaker_compact(bigMap,finalMap)
    Type(MapMakerParams), Intent(In) :: bigMap
    Type(MapMakerParams), Intent(In) :: finalMap

  End Subroutine mapmaker_compact

  !----------------------------------------------------------------------------
  ! Displays information stored in the Generic Map variable
  !----------------------------------------------------------------------------
  Subroutine mapmaker_display(mmparams,unit)
    Type(MapInfo), Intent(In) :: mmparams
    Integer, Intent(In)              :: unit

    Integer :: i

    !** Introduce yourself
    Write(unit,'(a)') dashedline
    Write(unit,'(a)') "Map Maker Information"
    Write(unit,'(a,i4)') "Number of maps : ",mmparams%nmaps

    !** Report the information about the individual maps
    Do i = 1, mmparams%nmaps
      Write(unit,'(a)') dashedline
      Write(unit,*) 
      Write(unit,'(a,i4)') "Map ",i
      Write(unit,'(a,i3,2a)') "Map sorbate     : (",mmparams%maps(i)%mapsorb, &
          ") ",trim(molecules_name(mmparams%maps(i)%mapsorb))
      Write(unit,'(a,i3,2a)') "Probe sorbate   : (",mmparams%maps(i)%probe, &
          ") ",trim(molecules_name(mmparams%maps(i)%probe))
      Write(unit,'(2a)')      "Interaction     :  ",mmparams%maps(i)%int
      Write(unit,'(2a)')      "Forcefield Model: ",mmparams%maps(i)%model
      Write(unit,'(a,3i4)')  "No. pnts x,y,z  :  ",mmparams%maps(i)%npnts(1), &
          mmparams%maps(i)%npnts(2),mmparams%maps(i)%npnts(3)
      Write(unit,'(a,f10.3)') "Grid spacing    :  ",mmparams%maps(i)%gridspace
      Write(unit,'(a,f10.3,a)') "Potential cutoff:  ", &
          mmparams%maps(i)%hipotcut, " kJ/mol"
      Write(unit,'(2a)') "Filename        :  ",Trim(mmparams%maps(i)%filename)
    End Do

    !** Make it look nice and write another dashed line
    Write(unit,'(a)') dashedline

  End Subroutine mapmaker_display

End Module mapmaker







