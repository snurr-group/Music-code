!----------------------------------------------------------------------------
! This module handles pre-specified mappings between two molecules.
! For general matching tools see 'match'.  The molecule file entry should
! look like this:
!
! MoleculeMap_Info: 
!    2  #number of maps
!   [species_name] 1-3@31-33;5-7@34-36 #new species, atom mapping  (old->new)
!
! See the str2intmap util function for mapping syntax specifications.
!----------------------------------------------------------------------------

Module molmap

  Use defaults, Only: RDbl, strLen, lstrLen
  Use utils, Only: filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, checkandstop, int2str, findint, str2intmap
  Use file, Only: file_open
  Use random, Only: rranf
  Use vector, Only: VecType
  Use molecules, Only: molecules_name,molecules_gettype,molecules_getfilename, &
      molecules_getlinewtag,molecules_getnatoms, molecules_getneighbors, &
      molecules_getmaxnatoms
  Use config, Only: AtMolCoords, config_getrp, config_getfixed

  Implicit None
  Save

  Private
  Public :: MolecMolec_Map,molmap_init,molmap_findbrokenbonds,molmap_map, &
      molmap_doublemap,molmap_brokenbond,molmap_mapatm,molmap_whichspc, &
      molmap_destspcs,molmap_display,molmap_clean

  !** The molecule->molecule mapping specifications
  !** individual atom number maps are atom_newspc = amap(atom_oldspc)
  Type MolecMolec_Map
    Integer                          :: spc        !** origin species number
    Integer                          :: nmaps,natoms,nbroken
    Integer, Dimension(:), Pointer   :: tospc    
    Integer, Dimension(:,:), Pointer :: brokenbond
    Integer, Dimension(:,:), Pointer :: amap       !** (spc,atom_oldspc)
    Character(len=lstrLen), Dimension(:), Pointer  :: mapstr
  End Type MolecMolec_Map

Contains
  !---------------------------------------------------------------------------
  ! Initializes a molecule map from a file.  If no filename is passed, it
  ! will be assumed that the information is in the molecule file.
  ! Requires:  mmmap -- molecule map to initialize
  !            spc -- species number
  !            filename -- filename where initialization data can be found
  !---------------------------------------------------------------------------
  Subroutine molmap_init(mmmap,spc,filename)
    Type(MolecMolec_Map), Intent(Out)        :: mmmap
    Integer, Intent(In)                      :: spc
    Character(*), Intent(In), Optional       :: filename

    Integer                              :: unit,i,n,error
    Integer                              :: sum,atm,mapno,newspc,natoms
    Character(len=strLen)                :: string1,string2
    Character(len=lstrLen)               :: tag,line
    Character(len=strLen), Dimension(20) :: fields

    !** Set the easy stuff
    Call molmap_nullify(mmmap)
    mmmap%spc = spc
    mmmap%nbroken = 0
    mmmap%natoms = molecules_getnatoms(spc)

    !** Find the molecule map specifications
    tag = 'MoleculeMap_Info'
    If (Present(filename)) Then     !** get from arbitrary file
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' not prepared to read molecule map from arbitrary file'
      Stop
    Else    !** get from the molecule file
      If (.Not. molecules_getlinewtag(spc,tag,line,unit)) Then
        string1 = molecules_getfilename(spc)
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Could not find molecule map in molecule file: ',Trim(string1)
        Write(0,'(2a)') 'Desired Tag: ',Trim(tag)
        Stop
      End If
    End If

    !** Read the number of maps and allocate memory
    Read(unit,'(i5)') mmmap%nmaps
    Allocate(mmmap%tospc(mmmap%nmaps),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'tospc')    
    Allocate(mmmap%amap(mmmap%nmaps,mmmap%natoms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'amap')    
    Allocate(mmmap%mapstr(mmmap%nmaps),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'mapstr')    
    mmmap%amap = 0

    !** Read the specifications
    sum = 0
    i = 0
    Do 
      i = i + 1
      Read(unit,'(a)') line
      line = stripcmnt(line)
      n = split(line,fields)

      !** Interpret the species name that this species is being mapped into
      mmmap%tospc(i) = molecules_gettype(fields(1))
      If (mmmap%tospc(i) == 0) Then
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Could not identify molecule name: "',Trim(fields(1)),'"'
        Stop
      End If

      !** Get the map itself
      mmmap%mapstr(i) = fields(2)
      sum = sum + str2intmap(fields(2),mmmap%amap(i,:))

      If (i == mmmap%nmaps) Exit
    End Do

    !** Make sure there aren't more mappings than atoms
    If (sum > mmmap%natoms) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' More atom mappings than atoms in molecule'
      Call molmap_display(mmmap,2,0)
      Stop
    End If

    !** Make sure none of the mapped atom numbers exceed number of atoms
    Do atm = 1,mmmap%natoms
      newspc = molmap_whichspc(mmmap,atm,mapno)
      If (newspc == 0) Cycle
      natoms = molecules_getnatoms(newspc)
      If (mmmap%amap(mapno,atm) > natoms) Then
        string1 = int2str(atm)
        string2 = molecules_name(newspc)
        Write(0,'(1x,2a,i4,4a)') __FILE__," : ",__LINE__, &
            ' new atom number for atom ',Trim(string1), &
            ' exceeds natoms of species ',Trim(string2)
        Call molmap_display(mmmap,2,0)
        Stop
      End If
    End Do

    !** Make sure all atoms are mapped
    Do atm = 1,mmmap%natoms
      If (molmap_whichspc(mmmap,atm) == 0) Then
        string1 = int2str(atm)
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Mapping is incomplete, atom number ',Trim(string1),' unmapped'
        Call molmap_display(mmmap,2,0)
        Stop
      End If
    End Do

  End Subroutine molmap_init

  !---------------------------------------------------------------------------
  ! Nullify the pointers in a molecule map
  ! Requires:  mmmap -- molecule map to initialize
  !---------------------------------------------------------------------------
  Subroutine molmap_nullify(mmmap)
    Type(MolecMolec_Map), Intent(InOut)        :: mmmap

    Nullify(mmmap%tospc)
    Nullify(mmmap%brokenbond)
    Nullify(mmmap%amap)
    Nullify(mmmap%mapstr)

  End Subroutine molmap_nullify

  !---------------------------------------------------------------------------
  ! Find the bonds that are broken during a given molecule-molecule mapping.
  ! Stores the resulting information in the map type itself for later use.
  ! Requires:  mmmap -- molecule map for which to initialize broken bond info
  !---------------------------------------------------------------------------
  Subroutine molmap_findbrokenbonds(mmmap)
    Type(MolecMolec_Map), Intent(InOut)        :: mmmap

    Integer                   :: i,error,npairs,nneighbors
    Integer                   :: atm,basespc,connectspc
    Integer, Dimension(100)   :: neighbors
    Integer, Dimension(100,2) :: pairs

    npairs = 0
    pairs = 0

    !** Loop through the atoms in the molecule and check connections
    Do atm = 1,mmmap%natoms
      !** Get the neighbors
      nneighbors = molecules_getneighbors(mmmap%spc,atm,neighbors)

      !** Get the species number that this atom would map into
      basespc = molmap_whichspc(mmmap,atm)

      !** Determine if any of the connecting atoms are in a different species
      Do i = 1,nneighbors
        
        !** Skip those neighboring atoms that have already been checked as 'atm'
        If (neighbors(i) < atm) Cycle

        !** Get the species number the neighbor would map into
        connectspc = molmap_whichspc(mmmap,neighbors(i))
        If (connectspc == 0) Cycle

        If (connectspc /= basespc) Then
          npairs = npairs + 1
          pairs(npairs,1) = atm
          pairs(npairs,2) = neighbors(i)
        End If
      End Do
    End Do

    !** Allocate space and save the broken bond atom numbers
    mmmap%nbroken = npairs
    Allocate(mmmap%brokenbond(mmmap%nbroken,2),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    
    mmmap%brokenbond(1:npairs,:) = pairs(1:npairs,:)

  End Subroutine molmap_findbrokenbonds

  !---------------------------------------------------------------------------
  ! Returns one of the atom number pairs that specify a broken bond in the
  ! mapping.  Also returns the number of broken bonds optionally
  ! Requires:  mmmap -- molecule map 
  !            bondno -- broken bond number in list
  !            spcpair -- returned species numbers for bondpair atoms
  !            bondpair -- returned atom numbers in broken bond
  !            nbroken -- optionally returned number of broken bonds
  !---------------------------------------------------------------------------
  Subroutine molmap_brokenbond(mmmap,bondno,spcpair,bondpair,nbroken)
    Type(MolecMolec_Map), Intent(In)        :: mmmap
    Integer, Intent(In)                     :: bondno
    Integer, Dimension(2), Intent(Out)      :: spcpair,bondpair
    Integer, Intent(Out), Optional          :: nbroken

    If (Present(nbroken)) nbroken = mmmap%nbroken

    If (.Not. Associated(mmmap%brokenbond)) Then
      Write(0,'(1x,2a,i4,a,i3)') __FILE__," : ",__LINE__, &
          ' Broken bond information has not been initialized'
      Stop
    End If

    If (bondno > mmmap%nbroken) Then
      Write(0,'(1x,2a,i4,a,i3)') __FILE__," : ",__LINE__, &
          ' Requested broken bond number does not exist ',bondno
      Stop
    End If

    bondpair = mmmap%brokenbond(bondno,:)
    spcpair(1) = molmap_whichspc(mmmap,bondpair(1))
    spcpair(2) = molmap_whichspc(mmmap,bondpair(2))

  End Subroutine molmap_brokenbond

  !---------------------------------------------------------------------------
  ! Find the species number that the molecule map would map a particular 
  ! atom number into.  Returns zero if it couldn't find a match.
  ! Requires:  mmmap -- molecule map 
  !            atm -- atom number
  !            mapno -- optional returned map number
  !---------------------------------------------------------------------------
  Integer Function molmap_whichspc(mmmap,atm,mapno)
    Type(MolecMolec_Map), Intent(In)         :: mmmap
    Integer, Intent(In)                      :: atm
    Integer, Intent(Out), Optional           :: mapno
    
    Integer        :: i,newatm

    molmap_whichspc = 0
    Do i = 1,mmmap%nmaps
      newatm = mmmap%amap(i,atm)
      If (newatm <= 0) Cycle
      molmap_whichspc = mmmap%tospc(i)
      If (Present(mapno)) mapno = i
      If (molmap_whichspc /= 0) Return
    End Do

  End Function molmap_whichspc

  !---------------------------------------------------------------------------
  ! Return the destination species numbers that the molecule map contains.
  ! Requires:  mmmap -- molecule map 
  !            spclist -- list of species numbers
  !---------------------------------------------------------------------------
  Integer Function molmap_destspcs(mmmap,spclist)
    Type(MolecMolec_Map), Intent(In)         :: mmmap
    Integer, Dimension(:), Intent(Out)       :: spclist
    
    molmap_destspcs = mmmap%nmaps
    spclist(1:mmmap%nmaps) = mmmap%tospc

  End Function molmap_destspcs

  !---------------------------------------------------------------------------
  ! Returns the atom number in the new species
  ! Requires:  mmmap -- molecule map 
  !            atm -- old atom number
  !            spc -- optional returned new species number
  !---------------------------------------------------------------------------
  Integer Function molmap_mapatm(mmmap,atm,spc)
    Type(MolecMolec_Map), Intent(In)         :: mmmap
    Integer, Intent(In)                      :: atm
    Integer, Intent(Out), Optional           :: spc
    
    Integer        :: i

    Do i = 1,mmmap%nmaps
      molmap_mapatm = mmmap%amap(i,atm)
      If (molmap_mapatm <= 0) Cycle
      If (Present(spc)) spc = mmmap%tospc(i)
      If (molmap_mapatm > 0) Return
    End Do

  End Function molmap_mapatm

  !---------------------------------------------------------------------------
  ! Use a molecule map to produce xyz coordinates for one of the designated
  ! species types.  This routine can either map into or out of the origin
  ! species type specified in the molecule map type.  Essentially it takes
  ! the principal xyz coordinates from config for a single molecule and 
  ! reorders some of all of them into a different species type.
  ! Requires:  mmmap -- molecule map to utilize
  !            species -- species data structure
  !            ospc -- origin species number
  !            omol -- origin molecule number
  !            dspc -- destination species number to map into
  !            coords -- new coordinates for specified spc number
  !---------------------------------------------------------------------------
  Subroutine molmap_map(mmmap,species,ospc,omol,dspc,coords)
    Type(MolecMolec_Map), Intent(In)             :: mmmap
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: ospc,omol,dspc
    Type(VecType), Dimension(:), Intent(InOut)   :: coords

    Integer       :: index,toatm,atm,natoms

    !** Act based on what the origin species is in the molecule map
    If (ospc == mmmap%spc) Then   !** map from map origin species to aux species
      !** Get the index of the species in the maps
      index = findint(mmmap%tospc,dspc)
      If (index == 0) Then
        Write(0,'(1x,2a,i4,a,i3)') __FILE__," : ",__LINE__, &
            ' Could not find destination species in map ',dspc
        Stop
      End If
  
      !** Map the atoms of the origin species individually
      Do atm = 1,mmmap%natoms
        toatm = mmmap%amap(index,atm)
!        Write(*,*) atm,' -> ',toatm
        If (toatm <= 0) Cycle
        coords(toatm) = config_getrp(species,(/mmmap%spc,omol,atm/))
      End Do

    Else If (findint(mmmap%tospc,ospc) /= 0) Then  !** aux species to map origin
      !** Get the index of the species in the maps
      index = findint(mmmap%tospc,ospc)
      natoms = molecules_getnatoms(ospc)

      !** Map the atoms of the origin species individually
      Do atm = 1,mmmap%natoms
        toatm = mmmap%amap(index,atm)
        If (toatm <= 0) Cycle
        coords(atm) = config_getrp(species,(/ospc,omol,toatm/))
      End Do      

    Else
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' Origin species is not represented in given molecule map'
      Stop
    End If

  End Subroutine molmap_map

  !---------------------------------------------------------------------------
  ! Make two map transformations.  The first map takes the coordinates to
  ! species that are common to the two maps.  The second mapping takes the
  ! coordinates to a species that is not present in the first map.
  ! For now, it is assumed that the mapping goes like this:
  !  ospc -> spc1,...,spcN -> dspc 
  ! Where ospc and dspc are the origin species of mmmap1 and mmmap2
  ! Requires:  mmmap1 -- 1st molecule map
  !            mmmap2 -- 2nd molecule map
  !            species -- species data structure
  !            ospc -- origin species number
  !            omol -- origin molecule number
  !            dspc -- destination species number to map into
  !            coords -- new coordinates for specified dspc number
  !---------------------------------------------------------------------------
  Subroutine molmap_doublemap(mmmap1,mmmap2,species,ospc,omol,dspc,coords)
    Type(MolecMolec_Map), Intent(In)             :: mmmap1,mmmap2
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Integer, Intent(In)                          :: ospc,omol,dspc
    Type(VecType), Dimension(:), Intent(InOut)   :: coords

    Integer                            :: i,j,index,toatm,atm,natoms
    Integer                            :: nfixed,max,error,putspc
    Logical, Save                      :: firsttime = .True.
    Character(len=strLen)              :: string
    Integer, Dimension(mmmap1%nmaps)   :: dspcmap,spclist
    Type(VecType), Dimension(:), Allocatable, Save  :: tempcoords

    !** Make sure the number of destination species is the same
    If (mmmap1%nmaps /= mmmap2%nmaps) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' number of destination species is not the same'
      Stop
    End If

    !** Make sure the two maps have common destination species
    !** Map index in mmmap1%tospc to mmmap2%tospc
    dspcmap = 0   
    Do i = 1,mmmap1%nmaps
      dspcmap(i) = findint(mmmap2%tospc,mmmap1%tospc(i))
      If (dspcmap(i) == 0) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            ' Could not match mmmap1 species ',mmmap1%tospc(i)
        string = int2str(mmmap1%tospc)
        Write(0,'(1x,2a)') 'mmmap1 species list: ',Trim(string)
        string = int2str(mmmap2%tospc)
        Write(0,'(1x,2a)') 'mmmap2 species list: ',Trim(string)
        Stop
      End If
    End Do

    !** Allocate temporary space for coordinates
    If (firsttime) Then
      firsttime = .False.
      nfixed = config_getfixed(species,spclist)
      max = molecules_getmaxnatoms(spclist(1:nfixed))
      Allocate(tempcoords(max),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End If 

    Do i = 1,mmmap1%nmaps
      !** Conduct first mapping:  ospc -> spc1,...,spcN
      putspc = mmmap2%tospc(dspcmap(i))
      natoms = molecules_getnatoms(putspc)
      Call molmap_map(mmmap1,species,ospc,omol,putspc,tempcoords(1:natoms))

      !** Conduct second mapping:  spc1,...,spcN -> dspc 
      !** Get the index of the species in the maps
      index = findint(mmmap2%tospc,putspc)
      natoms = molecules_getnatoms(putspc)

      !** Map the atoms of the final destination species individually
      Do atm = 1,mmmap2%natoms
        toatm = mmmap2%amap(index,atm)
        If (toatm <= 0) Cycle  !** skip those that aren't in map
        coords(atm) = tempcoords(toatm)
      End Do      
    End Do

  End Subroutine molmap_doublemap

  !----------------------------------------------------------------------
  ! Display the molecule map data type
  ! Requires:  mmmap -- molecule map data type
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine molmap_display(mmmap,indent,unit)
    Type(MolecMolec_Map), Intent(In)        :: mmmap
    Integer, Intent(In)                     :: indent,unit

    Integer                     :: i,nnew,atm
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2
    Integer, Dimension(20)      :: newatm

    blank = Repeat(' ',indent)

    string1 = molecules_name(mmmap%spc)

    Write(unit,'(4a)') blank,'Molecule map for ',Trim(string1),' species'
    Do i = 1,mmmap%nmaps
      string1 = molecules_name(mmmap%tospc(i))
      Write(unit,'(2x,5a)') blank,'Map to ',Trim(string1),': ', &
          Trim(mmmap%mapstr(i))
    End Do

    If (Associated(mmmap%brokenbond)) Then
      Do i = 1,mmmap%nbroken
        Write(unit,'(2x,2a,i3,3x,2i4)') blank,'broken bond pair ',i, &
            mmmap%brokenbond(i,1),mmmap%brokenbond(i,2)
      End Do
    End If

    Return

    Write(unit,'(2a)') blank,'Raw map: original atom, new species, new atom(s)'
    Do atm = 1,mmmap%natoms
      string1 = ''
      nnew = 0
      Do i = 1,mmmap%nmaps
        If (mmmap%amap(i,atm) > 0) Then
          nnew = nnew + 1
          newatm(nnew) = mmmap%amap(i,atm)
          string2 = molecules_name(mmmap%tospc(i))
          Write(string1,'(a,1x,a)') Trim(string1),Trim(string2)
        End If
      End Do
      If (string1 == '') Then
        string1 = 'NONE!'
        Write(unit,'(i4,2x,a)') atm,Trim(string1)
      Else
        Write(unit,'(i4,2x,a,2x,5i4)') atm,Trim(string1),newatm(1:nnew)
      End If
    End Do

  End Subroutine molmap_display

  !----------------------------------------------------------------------
  ! Clean the molecule map data type
  ! Requires:  mmmap -- molecule map data type
  !----------------------------------------------------------------------
  Subroutine molmap_clean(mmmap)
    Type(MolecMolec_Map), Intent(InOut)        :: mmmap

    Integer          :: error

    Deallocate(mmmap%tospc,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'tospc')    
    Deallocate(mmmap%amap,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'amap')    
    Deallocate(mmmap%mapstr,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'mapstr')    

  End Subroutine molmap_clean

End Module molmap
