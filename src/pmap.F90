!------------------------------------------------------------------------------
! This module handles the non-generic potential maps, i.e. the old pmap and
! the electrostatic map types.  These maps allow us to pre-tabulated the 
! interactions (potential energy and force) between a point particle and an
! entire FIXED species.  They are in the same module because they are
! similar enough to use some of the same routines.  Note that the storage for
! these maps is also in this module to avoid duplication; pointers to this
! storage can be found in other data structures.
!
! Needed Improvements
! 1) need to activate the GENERATE option and integrate with mapmaker.F90
!------------------------------------------------------------------------------
Module pmap

  Use defaults, Only: RDbl, strLen, RSgl, MAX_MAPS, lstrLen
  Use utils, Only: split, getpath, allocErrDisplay, int2str, cleanstring
  Use file, Only: file_open
  Use vector, Only: VecType
  Use simcell, Only: SimCell_Params, simcell_getell
  Use amap, Only:  AMap_Info, amap_init, amap_int, amap_display
  Use mapheader, Only: mapheader_getheader, mapheader_cleanup, Header_Info, &
      mapheader_getIndices, mapheader_writeheader, mapheader_display, &
      mapheader_boxinfo
  Use interpolate, Only: Map_Core_Info, interpolate_int, interpolate_writecore, &
      interpolate_initcore, interpolate_readcore, interpolate_clean, &
      interpolate_npts
  Use storebase, Only: EnergyPlus, storebase_nderivs

  Implicit None
  Save

  Private
  Public :: Pmap_Params, pmap_pmapinit, pmap_emapinit, pmap_int, pmap_eint, &
      pmap_compare, pmap_display, pmap_clean, pmap_idstring, emap_idstring, &
      pmap_newmap, pmap_compact, pmap_addpoint, pmap_write, pmap_init, &
      pmap_disp, pmap_boxinfo

  !** Potential Map data type, used for both pmaps and emaps
  Type Pmap_Params
    Logical :: analytical
    Character(len=strLen)       :: filename
    Type(Header_Info)           :: header
    Integer                     :: npoints
    Type(Map_Core_Info)         :: core
    Type(AMap_Info)         :: a_map
  End Type Pmap_Params

  Character(len=strLen), Parameter         :: pmap_idstring = 'PMAP'
  Character(len=strLen), Parameter         :: emap_idstring = 'EMAP'

  !** The old potential map storage
  Integer            :: total_pmaps = 0
  Type(Pmap_Params), Dimension(:), Pointer     :: pmaps

  !** The Electrostatic map storage
  Integer            :: total_emaps =  0
  Type(Pmap_Params), Dimension(:), Pointer     :: emaps

Contains

  !----------------------------------------------------------------------------
  ! Initialize a potential map.  This routine interprets the strings and forks
  ! to _pmapinit or _emapinit
  ! Requires:  map -- Potential map pointer to set
  !            params -- parameters to initialize map with, 1st should be type
  !            eff -- unit cell edge lengths
  !----------------------------------------------------------------------------
  Subroutine pmap_init(map,params,eff)
    Type(Pmap_Params), Pointer                      :: map
    Character(len=strLen), Dimension(:), Intent(In) :: params
    Real(Kind=Rdbl), Dimension(3), Intent(In)       :: eff

    If (Size(params) > 2) Then
      Write(0,'(2a,i6,a,i3)') __FILE__,":",__LINE__, &
          ' to many parameters passed ',Size(params)
      Stop
    End If

    If (Trim(params(2)) == 'GENERATE') Then
      Write(0,'(2a,i6,a,3i4)') __FILE__,":",__LINE__, &
          ' GENERATE option not yet available for Potential maps'
      Stop      
    End If

    !** Call the appropriate pmap or emap init routine
    Select Case(Trim(params(1)))
    Case(Trim(pmap_idstring))
      Call pmap_pmapinit(map,eff,params(2))

    Case(Trim(emap_idstring))
      Call pmap_emapinit(map,eff,params(2))

    Case Default
      Write(0,'(2a,i6,2a)') __FILE__,":",__LINE__, &
          ' Could not interpret string ',Trim(params(1))
      Stop      
    End Select

  End Subroutine pmap_init

  !---------------------------------------------------------------------------
  ! Initialize a pmap using provided information
  ! Requires:  pmap -- pmap to initialize (pointer)
  !            eff -- unit cell edge lengths
  !            filename -- file containing map, to read from
  !---------------------------------------------------------------------------
  Subroutine pmap_pmapinit(pmap, eff, filename)
    Type(Pmap_Params), Pointer                :: pmap
    Real(kind=RDbl), Dimension(3), Intent(In) :: eff
    Character(*), Intent(In)                  :: filename

    Integer                        :: i,mapno
    Integer                        :: error,unitno
    Logical :: analytical
    Character(len=strLen) :: a_filename
    !** Initialize the list the holds all the sitemaps read in
    If (total_pmaps == 0) Then
      Allocate(pmaps(MAX_MAPS), Stat=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not allocate 'pmaps'"
        Stop
      End If
    End If

    analytical=.False.
    analytical = pmap_checkanalytical(filename,a_filename)

    !** Check if the map has already been read in
    Do i = 1,total_pmaps
      If (filename == pmaps(i)%filename) Then
        pmap => pmaps(i)
        Return
      End If
    End Do

    !** Check that we haven't exceeded the number of maps
    If (MAX_MAPS == total_pmaps) Then
      Write(0,'(1x,2a,i4, a, i4)') __FILE__," : ",__LINE__, &
          "Total number of pmaps exceeded ", MAX_MAPS
      Stop
    End if

    total_pmaps = total_pmaps + 1
    mapno = total_pmaps
    pmaps(mapno)%filename = filename
    pmaps(mapno)%npoints = 0
    If (analytical) Then
      Call pmap_initanalytical(pmap,pmaps,mapno,eff,a_filename)
    Else
      !** Actually read the potential map
      Call pmap_readmap(pmap,pmaps,mapno,eff,filename)
    Endif

    pmap%analytical=analytical

  End Subroutine pmap_pmapinit
  
  !---------------------------------------------------------------------------
  ! Redirects the call to amap to initialize analytical map
  ! Requires:  pmap -- potential map pointer to initialize
  !            pmaps -- array of maps
  !            mapno -- index in maps array to read map into
  !            eff -- unit cell edge lengths
  !            filename -- text file containing analytical function info 
  !---------------------------------------------------------------------------
  Subroutine pmap_initanalytical(pmap,pmaps,mapno,eff,filename)
    Type(Pmap_Params), Pointer                :: pmap
    Type(Pmap_Params), Dimension(:), Pointer  :: pmaps
    Integer, Intent(In)                       :: mapno
    Real(kind=RDbl), Dimension(:), Intent(In) :: eff
    Character(*), Intent(In)                  :: filename
    
    Character(len=lstrLen)     :: filewpath
    Integer                    :: unitno

    !** Assign the pointer to the newly read map
    pmap => pmaps(mapno)

    !** Get the path and generate the full filename
    Write(filewpath,'(2a)') Trim(getpath('PMAPDIR')),Trim(filename)

    !** Open the file if it is not already open
    unitno = file_open(filewpath,110)

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Reading analytical pmap: ',Trim(filewpath)

    !** Get the amap info 
    Call amap_init(pmap%a_map, eff, unitno)
    Call amap_display(pmap%a_map, 3, 6)


  End Subroutine pmap_initanalytical


  !-----------------------------------------------
  ! checks the fname string to decide analytical or not
  !-----------------------------------------------
  Logical Function pmap_checkanalytical(fname, a_fname)
    Character(len=strLen),intent(in) :: fname
    Character(len=strLen),Intent(out) :: a_fname
    Integer :: nfields
    Character(len=strLen),dimension(3) :: fields
    nfields=split(fname,fields,"-")

    pmap_checkanalytical=.False.

    If (nfields/=2) Then
      pmap_checkanalytical=.False.
    Else
      If ((fields(1)=="ANALYTICAL")) Then
        pmap_checkanalytical=.True.
        a_fname=cleanstring(fields(2))
      else
        Write(*,*)"Could not understand string: ", Trim(fname)
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      endif
    endif
  End Function pmap_checkanalytical

  !---------------------------------------------------------------------------
  ! Initialize the pmap "map".  We need to pass the zeolite edge lengths "ell"
  ! to set up the tabulated grid.
  !---------------------------------------------------------------------------
  Subroutine pmap_emapinit(emap, eff, filename)
    Type(Pmap_Params), Pointer                :: emap
    Real(kind=RDbl), Dimension(:), Intent(In) :: eff
    Character(*), Intent(In)                  :: filename
                                   
    Integer                        :: i,mapno
    Integer                        :: error,unitno

    !** Initialize the list that holds all the emaps
    If (total_emaps == 0) Then
      Allocate(emaps(MAX_MAPS), Stat=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not allocate 'emaps'"
        Stop
      End If
    End If

    !** Check if the map has already been read in
    Do i = 1,total_pmaps
      If (filename == emaps(i)%filename) Then
        emap => emaps(i)
        Return
      End If
    End Do

    !** Check that we haven't exceeded the number of maps
    If (MAX_MAPS == total_emaps) Then
      Write(0,'(1x,2a,i4, a, i4)') __FILE__," : ",__LINE__, &
          "Total number of pmaps exceeded ", MAX_MAPS
      Stop
    Endif

    total_emaps = total_emaps + 1
    mapno = total_emaps
    emaps(mapno)%filename = filename
    emaps(mapno)%npoints = 0

    !** Actually read the electrostatic map
    Call pmap_readmap(emap,emaps,mapno,eff,filename)

  End Subroutine pmap_emapinit

  !----------------------------------------------------------------------------
  ! Get the interaction from a potential map at a specified point
  ! Requires:  pmap -- Potential map 
  !            simcell -- simulation cell data structure
  !            pt -- position vector for particle
  !            output -- resultant energy plus possible derivatives
  !----------------------------------------------------------------------------
  Logical Function pmap_int(map,simcell,pt,output)
    Type(Pmap_Params), Intent(In)        :: map
    Type(SimCell_Params), Intent(In)     :: simcell  
    Type(VecType), Intent(In)            :: pt
    Type(EnergyPlus), Intent(InOut)        :: output

    If (map%analytical) Then

      pmap_int = amap_int(map%a_map, output, pt, simcell)

    Else
      !** Get the interactions by interpolating from map
      pmap_int = interpolate_int(map%core,map%header,output,pt,simcell)
    Endif
  End Function pmap_int

  !----------------------------------------------------------------------------
  ! Get the interaction from an electrostatic potential map at a specified point
  ! Requires:  pmap -- Potential map 
  !            simcell -- simulation cell data structure
  !            pt -- position vector for particle
  !            output -- resultant energy plus possible derivatives
  !            q -- charge of point particle
  !----------------------------------------------------------------------------
  Logical Function pmap_eint(map,simcell,pt,q,output)
    Type(Pmap_Params), Intent(In)        :: map
    Type(SimCell_Params), Intent(In)     :: simcell  
    Type(VecType), Intent(In)            :: pt
    Real(kind=RDbl), Intent(In)          :: q
    Type(EnergyPlus), Intent(InOut)      :: output

    !** Get the interactions by interpolating from map
    pmap_eint = interpolate_int(map%core,map%header,output,pt,simcell,q)

  End Function pmap_eint

  !---------------------------------------------------------------------------
  ! Read a potential map from file
  ! Requires:  map -- potential map pointer to initialize
  !            maplist -- array of maps
  !            newidx -- index in maps array to read map into
  !            eff -- unit cell edge lengths
  !            filename -- file containing map, to read from
  !---------------------------------------------------------------------------
  Subroutine pmap_readmap(map,maplist,newidx,eff,filename)
    Type(Pmap_Params), Pointer                :: map
    Type(Pmap_Params), Dimension(:), Pointer  :: maplist
    Integer, Intent(In)                       :: newidx
    Real(kind=RDbl), Dimension(:), Intent(In) :: eff
    Character(*), Intent(In)                  :: filename

    Integer                    :: unitno,error
    Character(len=lstrLen)     :: filewpath

    !** Get the path and generate the full filename
    Write(filewpath,'(2a)') Trim(getpath('PMAPDIR')),Trim(filename)

    !** Open the file if it is not already open
    unitno = file_open(filewpath,100)

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Reading potential map: ',Trim(filewpath)

    !** Get the pmap header information 
    Call mapheader_getheader(maplist(newidx)%header, eff, unitno)
    Call mapheader_display(maplist(newidx)%header, 3, 6)

    !** Allocate the necessary memory
    Call interpolate_initcore(maplist(newidx)%core,maplist(newidx)%header)

    !** Now read in the potential/force etc.. 
    Call interpolate_readcore(maplist(newidx)%core,unitno)

    !** Store the number of points in the map
    maplist(newidx)%npoints = interpolate_npts(maplist(newidx)%core)

    Close(unit=unitno)

    !** Assign the pointer to the newly read map
    map => maplist(newidx)

  End Subroutine pmap_readmap

  !----------------------------------------------------------------------------
  ! Writes a map to the given filename
  ! Requires:  map -- Generic map type to write to file
  !            filename -- filename to write to
  !----------------------------------------------------------------------------
  Subroutine pmap_write(map,filename)
    Type(Pmap_Params), Intent(InOut)  :: map
    Character(*), Intent(In)          :: filename

    Integer :: unitno

    !** Open the file if it is not already open
    unitno = file_open(filename,100,'map file')

    !** Write the header to the file
    Call mapheader_writeHeader(map%header,unitno)

    !** Now write the data points to the file
    Call interpolate_writecore(map%core,unitno)

    Close(unit=unitno)

  End Subroutine pmap_write

  !----------------------------------------------------------------------------
  ! Sets up the variables for a new map
  ! Requires:  map -- Generic map type to initialize
  !            nptnts -- number of points 
  !            pnts -- number of grid points along each axis, an array
  !----------------------------------------------------------------------------
  Subroutine pmap_newmap(map,npnts,pnts)
    Type(Pmap_Params), Intent(InOut)    :: map
    Integer, Intent(In)                 :: npnts
    Integer, Intent(In), Dimension(:)   :: pnts
    
    Integer :: error

    !** obviously this is HACK
    map%header%nfsize = npnts
    map%header%nbrx = pnts(1)
    map%header%nbry = pnts(2)
    map%header%nbrz = pnts(3)

    !** Allocate space for these points that need to be tabulated
    Call pmap_initTabulated(map,npnts)

    !** Set the pointer index to zero
    map%npoints = 0

  End Subroutine pmap_newmap

  !----------------------------------------------------------------------------
  ! Sizes the tabulated points arrays (potentials and HOTs)
  ! NOTE : Also destroys values stored in map%core%ptr
  ! Requires:  map -- Generic map type to initialize
  !            nptnts -- number of points 
  !            oldMap -- flag indicating the the map needs to be deallocated
  !----------------------------------------------------------------------------
  Subroutine pmap_initTabulated(map,npnts,oldMap)
    Type(Pmap_Params), Intent(InOut)        :: map
    Integer, Intent(In)                     :: npnts
    Logical, Intent(In), Optional           :: oldMap

    Integer :: error

    !** Check to see if the oldMap flag is passed. This means we 
    !** need to deallocate first.
    If (Present(oldMap)) Then
      If (oldMap) Then
        Call pmap_clean(map,.False.)
      End If
    End If

    !** obviously this is HACK
    map%header%nfsize = npnts

    !** Reallocate the core
    Call interpolate_initcore(map%core,map%header)

  End Subroutine pmap_initTabulated

  !---------------------------------------------------------------------------
  ! Stores a point in the map given the positions, potentials, and HOTs
  ! Requires:  map -- Generic map type to operate on
  !              r -- position vector?
  !            pot -- potential ?
  !            hot -- vector of ?
  !        storeIt -- flag?
  !---------------------------------------------------------------------------
  Subroutine pmap_addpoint(map,r,pot,hot,storeIt)
    Type(Pmap_Params), Intent(InOut) :: map
    Type(VecType), Intent(In)               :: r
    Real(Kind=RDbl), Intent(In)             :: pot
    Real(Kind=RDbl), Dimension(:), Intent(In) :: hot
    Logical, Intent(In)                     :: storeIt

    Integer               :: i,j,k,n
    Integer, Dimension(3) :: ind

    !** Get the pointer indices. We want the nearest integer, so pass .True.
    Call mapheader_getIndices(map%header,r%comp(1),r%comp(2), &
        r%comp(3),ind,.True.)

    !** Put them into i j k for easier coding
    i = ind(1)
    j = ind(2)
    k = ind(3)

    !** If we're not storing it, set the pointer to zero and
    !** exit
    If (.Not.storeIt) Then
      map%core%ptr(k,j,i) = 0
      Return
    Else
      !** Increment the number of points
      map%npoints = map%npoints + 1
      !** Add it to the "pointer" array
      map%core%ptr(k,j,i) = map%npoints
      !** Set i to the counter value, just to make the references
      !** below cleaner
      n = map%npoints
    End If

    !** All the quantites are stored as real singles

    !** The potential
    map%core%r(n) = Real(pot,Kind=RSgl)

    !** The higher order terms are stored as individual x,y,z components
    map%core%rx(n) = hot(2)
    map%core%ry(n) = hot(3)
    map%core%rz(n) = hot(4)

    map%core%rxy(n) = hot(5)
    map%core%rxz(n) = hot(6)
    map%core%ryz(n) = hot(7)

    !** all three components of hot(8) should be the same
    map%core%rxyz(n) = hot(8)

  End Subroutine pmap_addpoint

  !----------------------------------------------------------------------------
  ! Compacts a map by excluding those points above the cutoff, but
  ! keeps those below the cutoff and those needed for hermite interpolation
  ! Requires:  bigMap -- Generic map type to operate on
  !            lilMap -- littler map to create
  !            cutoff -- cut-off distance for keeping points
  !----------------------------------------------------------------------------
  Subroutine pmap_compact(bigMap,lilMap,cutoff)
    Type(Pmap_Params), Intent(In)     :: bigMap
    Type(Pmap_Params), Intent(InOut)  :: lilMap
    Real(Kind=RDbl), Intent(In)              :: cutoff

    Integer         :: i, j, k, ii, jj, kk, indx, jndx, kndx
    Integer         :: naccept, reject, acceptCut, acceptNeighbor
    Integer         :: nx, ny, nz
    Integer         :: thisOne, lilOne, bigOne, error
    Integer, Dimension(:,:,:), Pointer :: lilOneTempPtr
    !** Initialize the acccepted counter to zero
    naccept = 0
    reject = 0
    acceptCut = 0
    acceptNeighbor = 0

    !** Set the lilMap ptrs to zero
    lilMap%core%ptr(:,:,:) = 0

    !** Number of points in each direction
    nx = Size(bigMap%core%ptr,3)
    ny = Size(bigMap%core%ptr,2)
    nz = Size(bigMap%core%ptr,1)

    !** Loop over all the points
    Do i = 1, nx
      Do j = 1, ny
        Do k = 1, nz

          !** Get the pointer number
          thisOne = bigMap%core%ptr(k,j,i)

          !** Make sure this is not a skipped point
          If (thisOne == 0) Then 
            Write(0,'(2a,i6,a,3i4)') __FILE__,":",__LINE__, &
                " WARNING: Found a pointer of 0 in ptr arrary at ",k,j,i
          End If

          !** Warn if the point is exactly equal to zero
          If (bigMap%core%r(thisOne) == 0.0_RDbl) Then
            Write(0,'(2a,i6,a,3i6,a,f8.3)') __FILE__,":",__LINE__, &
                " Warning point ",k,j,i," is exactly 0 : ",bigMap%core%r(thisOne)
          End If

          !** Is it below the cutoff?
          If ((bigMap%core%r(thisOne) < Real(cutoff,kind=RSgl))) Then

            !** Add it in there if it isn't already there 
            If (lilMap%core%ptr(k,j,i) == 0) Then
              acceptCut = acceptCut + 1
              naccept = naccept + 1
              lilMap%core%ptr(k,j,i) = naccept
            End If

            !** Check the surrounding points, make sure
            !** those needed for hermite are included
            Do ii = -1, 1
              Do jj = -1, 1
                Do kk = -1, 1
                  If ((ii == 0).And.(jj == 0).And.(kk == 0)) Cycle
                  indx = i + ii
                  jndx = j + jj
                  kndx = k + kk

                  !** Invoke PBCs
                  If (indx > nx) indx = 1
                  If (jndx > ny) jndx = 1
                  If (kndx > nz) kndx = 1
                  If (indx < 1)  indx = nx
                  If (jndx < 1)  jndx = ny
                  If (kndx < 1)  kndx = nz

                  !** Add the points if it isn't there but is greater than
                  !** the cutoff
                  If ((lilMap%core%ptr(kndx,jndx,indx) == 0).And. &
                      (bigMap%core%r(bigMap%core%ptr(kndx,jndx,indx)) &
                      >= Real(cutoff,RSgl))) Then
                    acceptNeighbor = acceptNeighbor + 1
                    naccept = naccept+1
                    lilMap%core%ptr(kndx,jndx,indx) = naccept
                  End If
                End Do
              End Do
            End Do
          Else
            !** Check to make sure that we aren't excluding a
            !** neighboring point needed for interpolation
            If (lilMap%core%ptr(k,j,i) /= 0) Cycle
            lilMap%core%ptr(k,j,i) = 0
            reject = reject + 1
          End If
        End Do
      End Do
    End Do

    !** Set the number of points in the map
    lilMap%npoints = naccept

    !** We need to store the values in lilMap%core%ptr somewhere 
    !** before we re-Allocate it at pmap_inittabulated
    lilOneTempPtr=>lilMap%core%ptr
    Nullify(lilMap%core%ptr)

    !SDEBUG
!!$    Allocate(lilOneTempPtr(nz,ny,nk), STAT=error)
!!$    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
!!$    lilOneTempPtr=lilMap%core%ptr
    !SDEBUG


    !** After all that, we now have the number of tabulated points,
    !** naccept, and the complete pointer array. Now, size the 
    !** compact array of potentials and derivates and copy them in
    !** for tabulated points

    Call pmap_initTabulated(lilMap,naccept,.True.)

    !** Copy the ptr array that we stored earlier and free the tempptr memory
    lilMap%core%ptr=lilOneTempPtr
    Deallocate(lilOneTempPtr,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Do i = 1, nx
      Do j = 1, ny
        Do k = 1, nz
          lilOne = lilMap%core%ptr(k,j,i)
          bigOne = bigMap%core%ptr(k,j,i)
          If (lilOne /= 0) Then
            lilMap%core%r(lilOne)   = bigMap%core%r(bigOne)
            lilMap%core%rx(lilOne)  = bigMap%core%rx(bigOne)
            lilMap%core%ry(lilOne)  = bigMap%core%ry(bigOne)
            lilMap%core%rz(lilOne)  = bigMap%core%rz(bigOne)
            lilMap%core%rxy(lilOne) = bigMap%core%rxy(bigOne)
            lilMap%core%rxz(lilOne) = bigMap%core%rxz(bigOne)
            lilMap%core%ryz(lilOne) = bigMap%core%ryz(bigOne)
            lilMap%core%rxyz(lilOne)= bigMap%core%rxyz(bigOne)
          End If
        End Do
      End Do
    End Do

    !MDEBUG
    Write(0,*) __FILE__,__LINE__," Accepted ( < cutoff) : ",acceptCut
    Write(0,*) __FILE__,__LINE__," Accepted ( neighbor) : ",acceptNeighbor
    Write(0,*) __FILE__,__LINE__," Accepted (  total  ) : ",naccept
    Write(0,*) __FILE__,__LINE__," Rejected ( > cutoff) : ",reject

  End Subroutine pmap_compact

  !----------------------------------------------------------------------------
  ! Compares two maps and reports any differences to the screen
  ! Requires:  thismap -- potential map 1
  !            thatmap -- potential map 1
  !            tolerance -- tolerance above which to display
  !            unitno -- unit number to send display to
  !----------------------------------------------------------------------------
  Subroutine pmap_compare(thisMap,thatMap,tolerance,unitno)
    Type(Pmap_Params), Pointer    :: thisMap
    Type(Pmap_Params), Pointer    :: thatMap
    Integer, Intent(In)           :: unitno
    Real(Kind=RSgl), Intent(In)   :: tolerance

    Integer         :: i, j, k, pot
    Real(kind=RSgl) :: diff, thisPot, thatPot

    !** Compare the pointer arrays, make sure they are the same
    Do i = 1, thisMap%header%nbrx
      Do j = 1, thisMap%header%nbry
        Do k = 1, thisMap%header%nbrz
          If (((thatMap%core%ptr(k,j,i) /= 0).And. &
              (thisMap%core%ptr(k,j,i) /= 0)).Or. &
              ((thatMap%core%ptr(k,j,i) == 0).And. &
              (thisMap%core%ptr(k,j,i) == 0))) Then
            Cycle
          Else
            Write(unitno,'(a,3i4,2(2x,a,i6),2x,2f8.3)') &
                " Mismatch : ",k,j,i,"thisMap : ",thisMap%core%ptr(k,j,i), &
                "thatMap : ",thatMap%core%ptr(k,j,i),&
                thisMap%core%r(thisMap%core%ptr(k,j,i)), &
                thatMap%core%r(thatMap%core%ptr(k,j,i))
          End If
        End Do
      End Do
    End Do

    !** Mismatched potential counter
    pot = 0

    !** Compare the potentials, see if they are within tol
    Do i = 1, thisMap%header%nbrx
      Do j = 1, thisMap%header%nbry
        Do k = 1, thisMap%header%nbrz
          thisPot = thisMap%core%r(thisMap%core%ptr(k,j,i))
          thatPot = thatMap%core%r(thisMap%core%ptr(k,j,i))
          diff = Abs(thisPot - thatPot)
          If (diff > tolerance) Then
            pot = pot + 1
            Write(unitno,'(a,3i4,2(2x,a,f10.4))') "Mismatch potential : ", &
                k,j,i,"thisMap : ",thisPot,"thatMap : ",thatPot
          End If
        End Do
      End Do
    End Do

    !** Report the number of mismatches
    Write(unitno,'(a)') ''
    Write(unitno,'(a,i7)') "Total number of mismatched potentials : ",pot

  End Subroutine pmap_compare

  !----------------------------------------------------------------------------
  ! Returns the box increments for a potential map if they are available.
  ! Requires:  map -- potential map 
  !            boxincrements -- x,y,z increments for potential map
  !            boxsteps -- number of boxes in x,y,z directions
  !----------------------------------------------------------------------------
  Subroutine pmap_boxinfo(map,boxincrements,boxsteps)
    Type(Pmap_Params), Intent(In)               :: map
    Real(Kind=RDbl), Dimension(3), Intent(Out)  :: boxincrements
    Integer, Dimension(3), Intent(Out)          :: boxsteps

    Call mapheader_boxinfo(map%header,boxincrements,boxsteps)

  End Subroutine pmap_boxinfo

  !----------------------------------------------------------------
  ! Display the map information
  ! Requires:  map -- potential map to display
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional display unit number
  !----------------------------------------------------------------
  Subroutine pmap_display(map,indent,unit)
    Type(Pmap_Params), Intent(In) :: map
    Integer, Intent(In)           :: indent
    Integer, Intent(In)           :: unit
    
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'display routine not ready, what to display?'
    Stop

  End Subroutine pmap_display

  !----------------------------------------------------------------
  ! Returns a short string giving pertinent information about map
  ! Requires:  map -- potential map for which to return info
  !----------------------------------------------------------------
  Function pmap_disp(map)
    Type(Pmap_Params), Intent(In) :: map
    Character(len=lstrLen)        :: pmap_disp

    Character(len=strLen)      :: string
    
    string = int2str(map%npoints)
    Write(pmap_disp,'(2a,3x,2a)') 'map from: ',Trim(map%filename), &
        'tabulated points: ',Trim(string)

  End Function pmap_disp

  !----------------------------------------------------------------------------
  ! Cleanup the pmap pointers
  ! Requires:  map -- potential map to clean
  !            headertoo -- flag indicating if it should clean the header also
  !----------------------------------------------------------------------------
  Subroutine pmap_clean(map,headertoo)
    Type(Pmap_Params), Intent(InOut)   :: map
    Logical, Intent(In)                :: headertoo
    
    !** Clean the header if desired
    If (headertoo) Then
      Call mapheader_cleanup(map%header)
    End If

    Call interpolate_clean(map%core)

  End Subroutine pmap_clean

End Module pmap



