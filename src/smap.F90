!------------------------------------------------------------------------
! This module provides methods for working with sitemaps
!------------------------------------------------------------------------
Module smap

  Use defaults, Only: RDbl, strLen,dashedline,MAX_MAPS
  Use utils, Only: isfileopen, filesrchstr, toupper, split, stripcmnt, toint, &
      getpath, findint, allocErrdisplay
  Use file, Only: file_getunit
  Use simcell, Only: SimCell_Params, simcell_maptouc,simcell_getzeoell
  Use config,only:AtMolcoords
  Use vector,Only: VecType,vector_display,vector_getcomp,Assignment(=)
  Use stats, Only: Statistics
  Use file, Only: file_open
  Implicit None
  Save

  Private
  Public :: Smap_Params,smap_init,Smap_getPosSiteType, smap_display, &
      smap_getAllSiteTypes, smap_showxyz

  !** Site Maps
  Type Smap_Params
    Character(len=strLen)     :: filename
    Integer                   :: nsitetypes
    Integer                   :: nbrx, nbry, nbrz
    Integer, Dimension(:,:,:), Pointer   :: sitetype
    Type(VecType)              :: ell    ! zeolite edge length
    Type(VecType)              :: cubell ! Edge length of each cubelet
    Type(Simcell_Params), Pointer :: scell
  End Type Smap_Params

  Type(Smap_Params), Allocatable, Dimension(:), Target  :: smaps
  Integer                                         :: total_smaps = 0

Contains
  !--------------------------------------------------------------------
  ! This routine initializes and reads in a sitemap
  !--------------------------------------------------------------------
  Subroutine smap_init(smap, simcell, filename)
    Type(Smap_Params), Pointer               :: smap
    Type(Simcell_Params), Target, Intent(in) :: simcell
    Character(*), Intent(in)         :: filename

    Integer     :: error, unitno, smapno, i, nbrx, nbry, nbrz, nsitetypes
    Character(len=strLen) :: smapfile
    Real(kind=RDbl)       :: stepx, stepy, stepz, ellx, elly, ellz

    !** Initialize the list the holds all the sitemaps read in
    If (total_smaps == 0) Then
      Allocate(smaps(MAX_MAPS), Stat=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Could not allocate smaps"
        Stop
      End If
    End If

    !** Make sure that the sitemap has not already been read in
    smapno = 0
    Do i=1, total_smaps
      If (Trim(filename) == Trim(smaps(i)%filename)) Then
        smapno = i
        Exit
      End If
    End Do
    If (smapno /= 0) Then
      smap => smaps(smapno)
      Return
    End If
    total_smaps = total_smaps + 1
    smapno = total_smaps

    !** Set the filename and the simulation cell/zeolite this sitemap
    !** corresponds to
    smaps(smapno)%filename = filename
    smaps(smapno)%scell    => simcell

    !** Open the smap file
    smapfile = Trim(getpath('SMAPDIR'))//Trim(filename)
    unitno = file_getunit(smapfile)
    Open(file=smapfile, unit=unitno, status='old', &
        form='unformatted', IOSTAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not open the sitemap file ", Trim(smapfile)
      Stop
    End If
    Rewind(unitno)

    !** Read the number of sites and the no. of cubelets in the 3 directions
    Read(unitno) nsitetypes, nbrx, nbry, nbrz
    smaps(smapno)%nsitetypes = nsitetypes
    smaps(smapno)%nbrx   = nbrx
    smaps(smapno)%nbry   = nbry
    smaps(smapno)%nbrz   = nbrz

    !** Allocate the space for the sitemap
    Allocate(smaps(smapno)%sitetype(nbrx, nbry, nbrz), Stat=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not allocate memory for 'smaps()%sitetype'"
      Stop
    End If

    !** Read the sitemap
    Read(unitno) smaps(smapno)%sitetype

    !** Get the unit cell edge length
    smaps(smapno)%ell = simcell_getzeoell(simcell)

    !** Get the size of each cubelet
    ellx = vector_getcomp(smaps(smapno)%ell, 1)
    elly = vector_getcomp(smaps(smapno)%ell, 2)
    ellz = vector_getcomp(smaps(smapno)%ell, 3)
    stepx = ellx/(nbrx-1.0_RDbl)
    stepy = elly/(nbry-1.0_RDbl)
    stepz = ellz/(nbrz-1.0_RDbl)
    smaps(smapno)%cubell = (/stepx, stepy, stepz/)

    !** Set the smap pointer to the map just read in
    smap => smaps(smapno)
    Return
  End Subroutine smap_init

  !-----------------------------------------------------------------------
  ! This function returns the site type of the point "pos" in the sitemap
  ! "smap".  It returns a negative number if an error occurs
  !-----------------------------------------------------------------------
  Integer Function smap_getPosSiteType(smap, pos, nozero)
    Type(Smap_Params), Intent(in) :: smap
    Type(VecType), Intent(in) :: pos
    Type(VecType) :: ucPos
    Logical , Optional , Intent(in) :: nozero
    Integer    :: iindex, jindex, kindex
    Real(kind=RDbl) :: stepx, stepy, stepz
    Integer :: stype
#ifdef SMAPNOZERO
    Integer :: i,ix, iy, iz, iix, iiy, iiz
    Integer :: nbrx, nbry, nbrz, unitno
#endif
    stype=-1
    ucPos = simcell_maptouc(smap%scell,pos)

    !** Edge length of the cubelets
    stepx = vector_getcomp(smap%cubell, 1)
    stepy = vector_getcomp(smap%cubell, 2)
    stepz = vector_getcomp(smap%cubell, 3)

    !** Get the cubelet size
    !** Get the indices
    iindex = Int(vector_getcomp(ucPos, 1)/stepx) + 1
    jindex = Int(vector_getcomp(ucPos, 2)/stepy) + 1
    kindex = Int(vector_getcomp(ucPos, 3)/stepz) + 1
    stype = smap%sitetype(iindex, jindex, kindex)
#ifdef SMAPNOZERO
    If (stype==0) Then
      If (Present(nozero)) Then
        if (.not.nozero) then 
          smap_getPosSiteType = stype
          Return
        endif
      Else
        smap_getPosSiteType = stype
        return
      Endif
!unitno=file_open("temp.xyz")
!Write(unitno,*) "C  ",ucpos%comp(1), ucpos%comp(2), ucpos%comp(3)
      stype=-1
      nbrx=smap%nbrx
      nbry=smap%nbry
      nbrz=smap%nbrz

      ! search up to 6 layers for the nearest nonzero site   
      Do i=1,6
        Do ix=-i,i
          iix = iindex+ix
          If (iix>nbrx) iix=iix-nbrx
          Do iy=-i,i
            iiy = jindex+iy
            If (iiy>nbry) iiy=iiy-nbry
            Do iz=-i,i
              iiz = kindex+iz
              If (iiz>nbrz) iiz=iiz-nbrz

              stype = smap%sitetype(iix, iiy, iiz)
              If (stype>0) Exit

            end do
            If (stype>0) Exit
          end do
          If (stype>0) Exit
        end do
        If (stype>0) Exit
      End Do
    Endif
    If (stype<1)  Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
#endif
    smap_getPosSiteType = stype
  End Function smap_getPosSiteType

  !---------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------
  Subroutine smap_UpdateSiting(smap, sorbates, smap_stats)
    Type(Smap_Params), Intent(in) :: smap
    Type(AtMolCoords), Dimension(:), Intent(in)   :: sorbates
    Type(Statistics), Dimension(:), Intent(inout) :: smap_stats
  End Subroutine smap_updatesiting

  !---------------------------------------------------------------------
  ! returns all site types calibrated
  !---------------------------------------------------------------------
  Subroutine smap_getAllSiteTypes(smap, sitetypes, nsites)
    Type(Smap_Params), Intent(in) :: smap
    Integer, Dimension(:), Pointer :: sitetypes
    Integer, Intent(out):: nsites

    Integer :: i,j,k,error,stype, sites_found
    Integer,Parameter ::  NON_SITE_INTEGER=-11111111

    !** i guess this is how smap code was originaly written
    nsites =     smap%nsitetypes + 1

    Nullify(sitetypes)
    Allocate(sitetypes(nsites),STAT=error)
    sitetypes=(NON_SITE_INTEGER)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    sites_found=0
    Do i=1,smap%nbrx
      Do j=1,smap%nbry
        Do k=1,smap%nbrz

          stype=smap%sitetype(i,j,k)

          If (stype==NON_SITE_INTEGER) then
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
          Endif
          !          Write(*,*) sitetypes,stype,findint(sitetypes,stype)
          If (findint(sitetypes,stype)==0) Then
            sites_found=sites_found+1
            If (sites_found>nsites) Then
              Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
              Stop
            Endif
            sitetypes(sites_found)=stype
          Endif

        End Do
      End Do
    End Do

    If (sites_found<nsites) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

  End Subroutine smap_getAllSiteTypes


  !---------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------
  !  Subroutine smap_getAllSites(smap, sitetype, cubeletlist)
  !    Type(Smap_Params), Intent(in) :: smap
  !    Integer, Intent(in) :: sitetype
  !    Type(IntVecType), Dimension(:), Intent(out) :: cubeletlist
  !  End Subroutine smap_getAllSites


  !-----------------------------------------------------------------------
  ! This routine takes as input a sitemap "smap" and a position "center"
  ! and find all the cubelets around this center which have the same
  ! site type.  The list of cubelets is returned in "cubeletlist"
  !-----------------------------------------------------------------------
  !  Subroutine smap_GrowSite(smap, center, cubeletlist)
  !    Type(Smap_Params), Intent(in) :: smap
  !    Type(VecType), Intent(in) :: center
  !    Type(IntVecType), Dimension(:), Intent(out) :: cubeletlist
  !  End Subroutine smap_GrowSite


  !----------------------------------------------------------------------
  ! This routine displays the sitemap as an xyz file
  ! excludes any site==0
  ! skips 2 points between every point for size reduction
  !----------------------------------------------------------------------
  Subroutine smap_showxyz(smap, filename,sites)
    Type(Smap_Params), Intent(in) :: smap
    Character(len=strLen) :: filename
    Integer, Dimension(:), Intent(in) :: sites
    Integer   :: i, j, k,s, passed_sites, totpoints, unitno, skipped
    Real(kind=RDbl) :: dx, dy, dz 
    Logical :: foundone
    Integer, Parameter :: SYM_SIZE=2,MAX_SYMS=6,SKIP=20
    ! SYM_SIZE : size of symbols, MAX_SYMS = maxinimum number of symbols, 
    ! SKIP : skips so many points
    Character(len=SYM_SIZE),Dimension(MAX_SYMS) :: asymbol
    Character(len=SYM_SIZE) :: symbol
    Character(len=MAX_SYMS*(SYM_SIZE+1)) :: atoms

    asymbol= (/"He", "Ne", "Ar", "Kr", "Xe", "Rn"/)
    atoms=Repeat(" ",MAX_SYMS*(SYM_SIZE+1))

    passed_sites=Size(sites,1)
    If (passed_sites>MAX_SYMS) Then

      Write(*,*) "sitetypes cant be more than :", MAX_SYMS
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    ! count total sites
    totpoints=0
    skipped=0
    Do s=1,passed_sites
      If (sites(s)==0) Cycle! for now we avoid zero
      foundone=.False.
      Do  i=1,smap%nbrx
        Do  j=1,smap%nbry
          Do  k=1,smap%nbrz
            If (smap%sitetype(i,j,k)==sites(s)) Then
              skipped=skipped+1
              If (skipped>SKIP) Then 
                skipped=0
                totpoints=totpoints+1
                foundone=.True.
              Endif
            Endif
          End Do
        End Do
      End Do
      If (foundone) atoms=Trim(atoms)//asymbol(s)//"-"
    End Do

    dx = vector_getcomp(smap%cubell, 1)
    dy = vector_getcomp(smap%cubell, 2)
    dz = vector_getcomp(smap%cubell, 3)


    unitno = file_open(filename)
    Write(unitno,'(i15)')  totpoints
    Write(unitno, '(a)')  "from smap :"//Trim( smap%filename)//Trim(atoms)

    skipped=0
    ! retrace the same loop, but this rime write the map
    Do s=1,passed_sites
      If (sites(s)==0) Cycle         ! for now we avoid zero
      symbol=asymbol(s)

      Do  i=1,smap%nbrx
        Do  j=1,smap%nbry
          Do  k=1,smap%nbrz
            If (smap%sitetype(i,j,k)==sites(s)) Then
              skipped=skipped+1
              If (skipped>SKIP) Then 
                skipped=0

                Write(unitno,'(a,5x,3f15.3)') symbol, (i-1)*dx+dx/2, &
                    (j-1)*dy+dy/2, (k-1)*dz+dz/2
              endif

            Endif
          End Do
        End Do
      End Do
    end do
    close(unitno)
  End Subroutine smap_showxyz

  !----------------------------------------------------------------------
  ! This routine displays the various fields of the sitemap
  !----------------------------------------------------------------------
  Subroutine smap_display(smap, nspaces, optunit)
    Type(Smap_Params), Intent(in) :: smap
    Integer, Intent(in)  :: nspaces ! No. of spaces from the left column
    Integer, Optional, Intent(in)    :: optunit

    Integer   :: i, unitno
    Character(len=nspaces) :: spaces, string

    spaces = Repeat(' ',nspaces)
    unitno = 6
    If (Present(optunit)) unitno = optunit


    Write(unitno, '(2a)') spaces, dashedline
    Write(unitno, '(2a)') spaces, "The Site Map Section:"
    Write(unitno, '(a,2x,2a)') spaces, "Filename       : ", smap%filename
    Write(unitno, '(a,2x,a,i14)') spaces, "No. of Sites   : ", smap%nsitetypes
    Write(unitno, '(a,2x,a,3i5)') spaces, "No. of Cubelets(x, y, z)  : ", &
        smap%nbrx, smap%nbry, smap%nbrz
    string=        Trim(vector_display(smap%cubell, "f8.3"))
    Write(unitno, '(a,2x,a,2x,a)') spaces, "Cubelet Edge Lengths      : ",&
        string

  End Subroutine smap_display
End Module smap
