Module postmep
  Use defaults
  Use utils
  Use file
  Use stats
  Use simcell
  
  Implicit None
  Save

  Integer, Parameter :: MAX_PTS = 100     ! Maximum points in each plane
  Integer, Parameter :: MAX_CENTERS = 50  ! Maximum clusters in each plane
  Integer, Parameter :: MAX_NPLANES = 300 ! Maximum no. of planes

  !** This structure holds the result of a single annealing run
  Type AnnealResult_Info
    Real(kind=RDbl), Dimension(7)   :: energies ! Assuming 4 sorbs (7 pairs)
    Type(VecType), Dimension(3)     :: lastpos  ! 3 Moveable Sorbates
    Type(VecType), Dimension(3)     :: lasteuler
    Type(VecType)                   :: normal
  End Type AnnealResult_Info

  !**  This structure stores the data for each cluster in a plane.  It also
  !** stores information regarding which cluster it is connected to backwards
  !** and forward
  Type Cluster_Info
    Integer          :: prevplane, nextplane
    Integer          :: prevcent,  nextcent
    Integer          :: npts              ! No. of points in the cluster
    Integer, Dimension(MAX_PTS)    :: pts !indices of the points in cluster
    Type(VecType), Dimension(3)    :: com !centers-of-mass of the 3 molecules
    Type(Statistics), Dimension(7) :: avgnrg !avgnrg of the pts in cluster
  End Type Cluster_Info

  !**
  Type Plane_Info
    Integer                                     :: nanneals
    Type(AnnealResult_Info), Dimension(MAX_PTS) :: annealno
    Type(Statistics), Dimension(7)              :: avgPlaneEnergy
    Integer                                     :: ncenters !No. of clusters
    Type(Cluster_Info), Dimension(MAX_CENTERS)  :: clusters
  End Type Plane_Info

  Type POSTMEP_Params
    Character(len=strLen)  :: inbasefilename
    Character(len=strLen)  :: outbasefilename
    Integer, Dimension(10) :: startext, endext
    Character(len=strLen)  :: pentname
    Integer         :: pentsorbtype ! Penetrant sorbate type
    Integer         :: nplanes
    Integer         :: nsorbs
    Integer         :: nsets
    Integer         :: lpno       ! last plane number
    Integer         :: ncontpaths ! No. of continuous paths
    Integer,Dimension(MAX_CENTERS) :: pathhead ! The head of a path in plane 1
    Type(VecType)   :: normal
    Character(len=strLen), Dimension(MAX_SORBS) :: sorbnames
    Type(Plane_Info), Dimension(MAX_NPLANES)    :: planeno
    Real(kind=RDbl) :: clustertol
    Integer         :: planetol
  End Type POSTMEP_Params

  Character(len=strLen), Parameter :: default_postmep_tag = &
      "Post MEP Section"

  Character(len=strLen), Dimension(7) :: nrglabel = &
      (/"Meth-Block1  ", "Meth-Block2  ", "Meth-Sili    ", "Block1-Block2", &
        "Block1-Sili  ", "Block2-Sili  ", "Total        "/)

Contains
  !--------------------------------------------------------------
  ! Initialize the postmep structure
  !--------------------------------------------------------------
  Subroutine postmep_init(pmepparams, ctrl_filename, optctrl_tag)
    Type(POSTMEP_Params), Intent(inout) :: pmepparams
    Character(*), Intent(in)         :: ctrl_filename
    Character(*), Optional, Intent(in) :: optctrl_tag

    Character(len=strLen)     :: tag, line
    Integer                   :: lineno, unitno, i, nplanes
    Real(kind=RDbl)           :: x1, x2, x3

    If (Present(optctrl_tag)) Then
      tag = optctrl_tag
    Else
      tag = default_postmep_tag
    End If

    !** Open the ctrl file
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
      Open(unit=unitno, file=ctrl_filename)
    End If
    
    !** Find the starting of the mep section
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Could not find the tag : ", tag
      Stop
    End If

    Read(unitno,*) pmepparams%inbasefilename
    Read(unitno,*) x1, x2, x3
    pmepparams%normal = (/x1, x2, x3/)
    Read(unitno,*) pmepparams%pentname
    pmepparams%pentsorbtype = molecules_gettype(pmepparams%pentname)
    Read(unitno,*) pmepparams%clustertol
    Read(unitno,*) pmepparams%planetol
    Read(unitno,*) pmepparams%nsets
    Read(unitno,*) 
    nplanes = 0
    Do i=1, pmepparams%nsets
      Read(unitno,*) pmepparams%startext(i), pmepparams%endext(i)
      nplanes = nplanes + pmepparams%endext(i) - pmepparams%startext(i) 
    End Do
    pmepparams%nplanes = nplanes
    Read(unitno,*)
    Read(unitno,*) pmepparams%outbasefilename
    Read(unitno,*) pmepparams%nsorbs
    Read(unitno,*)
    Do i=1, pmepparams%nsorbs
      Read(unitno,*) pmepparams%sorbnames(i)
    End Do

    !** Make sure that if the no. of sets is more than one then
    !** the first set has negative numbers for extensions.  This
    !** ensures that the backward runs are analyzed before the forward
    !** runs and hence that we are going forward in our coordinate axis
    If (pmepparams%nsets > 1) Then
      If (pmepparams%startext(1) > 0 .Or. pmepparams%endext(1) > 0) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
            " The 'backward' file extensions should be specified first"
        Stop
      End If
    End If
    
    !** Read the different files
    pmepparams%lpno = 0
    Do i=1, pmepparams%nsets
      Call postmep_readdata(pmepparams, i)
    End Do
  End Subroutine postmep_init

  !--------------------------------------------------------------
  ! Read the data from the different files
  !--------------------------------------------------------------
  Subroutine postmep_readdata(pmepparams, setno)
    Type(POSTMEP_Params), Intent(inout) :: pmepparams
    Integer, Intent(in)   :: setno

    Integer      :: startext, endext, i, nrgunitno, sorbno, error, nrecs, j
    Integer      :: junk, nsorbs, lpno
    Integer, Dimension(MAX_SORBS)  :: avgunitno
    Character(len=strLen)  :: nrgfilename, avgfilename, sorbname, dir
    Real(kind=RDbl)        :: x, y, z, theta, phi, psi
    Real(kind=RDbl), Dimension(7)  :: pairnrg

    !** Cycle through all the input energy files and read the stats
    startext = pmepparams%startext(setno)
    endext   = pmepparams%endext(setno)
    lpno     = pmepparams%lpno
    nsorbs   = pmepparams%nsorbs
    Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) 'Reading from set ', setno, startext, endext
    If (startext < 0 .And. endext < 0) Then
      dir = "backward/"
    Else
      dir = "forward/"
    End If

    Do i=startext, endext
      
      lpno = lpno + 1

      !** Initialize the statistics that get the average energy
      Do j=1, 7
        Call stats_init(pmepparams%planeno(lpno)%avgPlaneEnergy(j), &
            Trim(nrglabel(j)), 5,.False., "f8.3")
      End Do

      !** Open the energy and the average position files
      ! Open the energyfile
      nrgfilename = &
         genfilename(Trim(dir)//Trim(pmepparams%inbasefilename)//".nrg",abs(i))
      nrgunitno = file_getunit(nrgfilename)
      Open(unit=nrgunitno, file=nrgfilename, status='old', IOSTAT=error)
      If (error /= 0) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
            " Could not open file : ", nrgfilename
        Stop
      End If
      ! Open the average position files
      Do sorbno=1, nsorbs
        sorbname = pmepparams%sorbnames(sorbno)
        avgfilename = Trim(dir)//Trim(pmepparams%inbasefilename)// &
            ".lastpos."//Trim(sorbname)
        avgfilename = genfilename(Trim(avgfilename), abs(i))
        avgunitno(sorbno) = file_getunit(avgfilename)
        Open(unit=avgunitno(sorbno), file=avgfilename, status='old', &
            IOSTAT=error)
        If (error /= 0) Then
          Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
              " Could not open file : ", avgfilename
          Stop
        End If
      End Do
      
      !** Read the points
      nrecs = 0
      Read(nrgunitno, *) ! Header line
      Do 
        ! Read the energy.  The 7 fields are as follows
        ! mb1, mb2, mz, b1b2, b1z, b2z, totalnrg
        Read(nrgunitno, *, IOSTAT=error) junk, junk, (pairnrg(j), j=1, 7)
        
        If (error /= 0) Then
          Write(*,*) 'Done processing file : ', Trim(nrgfilename)
          Write(*,*) 'No. of records read  : ', nrecs
          Close(nrgunitno)
          Exit
        End If
        nrecs = nrecs + 1
        If (nrecs > MAX_PTS) Then
          Write(0,'(1x,2a,i4, a, i5)') __FILE__," : ",__LINE__, &
              " Maximum no. of points exceeded in plane ", i
          Stop
        End If
        pmepparams%planeno(lpno)%annealno(nrecs)%energies = pairnrg
        ! Update the stats
        Do j=1, 7
          Call stats_update(pmepparams%planeno(lpno)%avgPlaneEnergy(j), &
              pairnrg(j))
        End Do

        ! Read the average position
        Do sorbno=1, nsorbs
          Read(avgunitno(sorbno),*) junk, x, y, z, theta, phi, psi
          pmepparams%planeno(lpno)%annealno(nrecs)%lastpos(sorbno)=(/x, y, z/)
          pmepparams%planeno(lpno)%annealno(nrecs)%lasteuler(sorbno) = &
              (/theta, phi, psi/)
        End Do
      End Do
      pmepparams%planeno(lpno)%nanneals = nrecs
    End Do

    pmepparams%lpno = lpno
  End Subroutine postmep_readdata

  !-------------------------------------------------------------------------
  ! Generates the a continuous path by finding cluster centers in each plane
  ! and joining them
  !-------------------------------------------------------------------------
  Integer Function postmep_genPath(pmepparams, simcell)
    Type(POSTMEP_Params), Intent(inout) :: pmepparams
    Type(Simcell_Params), Intent(in) :: simcell

    Integer, Dimension(MAX_CENTERS)  :: paths
    Integer    :: i, j, error, npaths
    Integer    :: npts, nsorbs, sorbno
    Character(len=strLen)   :: outxyzfile, avgposfilename
    Real(kind=RDbl)  :: objvalue, totalnrg, clustertol

    !** Group the different paths within a plane into clusters/centers
    Do i=1, pmepparams%lpno
      npts = pmepparams%planeno(i)%nanneals
      nsorbs = pmepparams%nsorbs
      clustertol = pmepparams%clustertol
      Write(*,*) 'Plane No. : ', i
      Call postmep_gencenters(pmepparams%planeno(i), npts, clustertol, nsorbs) 
    End Do

    !** Generate paths by connecting centers
    Do j=1, pmepparams%lpno
      Do i=1, pmepparams%planeno(j)%ncenters
        objvalue = postmep_joincenters(pmepparams, j, i, 1)
      End Do
    End Do
    
    !** Find those paths which go from the first to the last plane
    npaths = postmep_findcontpaths(pmepparams, paths)
    Write(*,'(a,i4)') 'No. of paths ', npaths
    Write(*,*) paths(1:npaths)
    pmepparams%pathhead = paths
    postmep_genPath = npaths
  End Function postmep_genPath

  !------------------------------------------------------
  ! Generates the cluster centers in a plane
  !------------------------------------------------------
  Subroutine postmep_genCenters(planeparams, totalanneals, clustertol, nsorbs)
    Type(Plane_Info), Intent(inout) :: planeparams
    Integer, Intent(in)          :: totalanneals
    Real(kind=RDbl), Intent(in)  :: clustertol
    Integer, Intent(in)          :: nsorbs

    Type(VecType), Dimension(3) :: vec1
    Type(VecType)               :: vec2
    Integer         :: i, j, center, pts, sorbno, npts 
    Integer         :: ncenters
    Type(Statistics), Dimension(3, MAX_CENTERS) :: xcom, ycom, zcom
    Type(Statistics), Dimension(7, MAX_CENTERS) :: nrg
    Real(kind=RDbl), Dimension(7) :: pairnrg
    Real(kind=RDbl) :: leastdist, dist, xyz(3), totalnrg, stddev

    !** Initialize the centers
    ncenters = 0
    Do i=1, MAX_CENTERS
      Do sorbno = 1, nsorbs
        planeparams%clusters(i)%com(sorbno)  = 0.0_RDbl
        Call stats_init(xcom(sorbno, i), "x com", 5,.False., "f8.3")
        Call stats_init(ycom(sorbno, i), "y com", 5,.False., "f8.3")
        Call stats_init(zcom(sorbno, i), "z com", 5,.False., "f8.3")
      End Do
      planeparams%clusters(i)%npts    = 0
      Do j=1, 7
        Call stats_init(nrg(j, i), nrglabel(j), 5,.False., "f8.3")
      End Do
    End Do
    
    !** Start a center with the first point
    ncenters = 1
    i = 1
    ! Assign the center of mass to the first center of each sorbate type
    Do sorbno=1, nsorbs
      vec1(sorbno) = planeparams%annealno(i)%lastpos(sorbno)
      planeparams%clusters(ncenters)%npts = 1
      planeparams%clusters(ncenters)%pts(1) = 1
      xyz = vec1(sorbno)
      Call stats_update(xcom(sorbno, ncenters), xyz(1))
      Call stats_update(ycom(sorbno, ncenters), xyz(2))
      Call stats_update(zcom(sorbno, ncenters), xyz(3))
    End Do
    ! Assign the energies for that annealing to each center
    pairnrg = planeparams%annealno(1)%energies
    Do j=1, 7
      Call stats_update(nrg(j, ncenters), pairnrg(j)) 
    End Do

    !** Now assign the other points to the different centers
    Do i=2, totalanneals
      pairnrg = planeparams%annealno(i)%energies
      Do sorbno =1, nsorbs
        vec1(sorbno) = planeparams%annealno(i)%lastpos(sorbno)
      End Do

      !** Cycle through all the centers and find the one that is closest
      leastdist = 100.0_RDbl     ! Some large number
      Do j=1, ncenters
        dist = 0.0_RDbl
        Do sorbno = 1, nsorbs
          vec2 = (/stats_getcavg(xcom(sorbno, j)), &
              stats_getcavg(ycom(sorbno, j)),stats_getcavg(zcom(sorbno, j))/)
          dist = dist + mag(vec1(sorbno) - vec2)
        End Do
        If (dist < leastdist) Then
          center = j
          leastdist = dist
        End If
      End Do
      
      !** If the minimum distance is too big we want to start
      !** a new center
      If (leastdist > clustertol) Then
        ! Start a new center
        ncenters = ncenters + 1
        center = ncenters
      End If

      !** Now add the new point to the center "center"
      planeparams%clusters(center)%npts = planeparams%clusters(center)%npts + 1
      npts = planeparams%clusters(center)%npts 
      ! Add the index of the annealno which lies in this cluster
      planeparams%clusters(center)%pts(npts) = i
      
      !** Update energy and center-of-mass of the center
      Do sorbno = 1, nsorbs
        xyz = vec1(sorbno)
        Call stats_update(xcom(sorbno, center), xyz(1))
        Call stats_update(ycom(sorbno, center), xyz(2))
        Call stats_update(zcom(sorbno, center), xyz(3))      
      End Do
      Do j=1, 7
        Call stats_update(nrg(j, center), pairnrg(j))
      End Do
    End Do

    !** Now we are done with all the points. 
    !** Assign the center of mass to the cluster center
    Do i=1, ncenters
      Do sorbno = 1, nsorbs
        planeparams%clusters(i)%com(sorbno) = &
            (/stats_getcavg(xcom(sorbno, i)), &
            stats_getcavg(ycom(sorbno, i)),stats_getcavg(zcom(sorbno, i))/)
      End Do
      Do j= 1, 7
        planeparams%clusters(i)%avgnrg(j) = nrg(j, i)
      End Do
    End Do
    planeparams%ncenters = ncenters
    
    Write(*,*) 'No. of center ', ncenters
    Do i=1, ncenters
      pts = planeparams%clusters(i)%npts
      totalnrg = stats_getcavg(planeparams%clusters(i)%avgnrg(7))
      stddev   = stats_getstd(planeparams%clusters(i)%avgnrg(7))
      Write(*,'(2i3, 3a19, f9.2, f6.3)') i, pts, &
          (Trim(vector_display(planeparams%clusters(i)%com(j),"f6.2")),j=1,3),&
          totalnrg, stddev
    End Do
  End Subroutine postmep_genCenters

  !------------------------------------------------------------------------
  ! This routine looks ahead "depth" number of planes to find the optimal
  ! plane and center "cplane" and "cent1" should connnect to.  As it cycles
  ! through "1 - depth" number of planes it picks the optimal path at each
  ! level
  !------------------------------------------------------------------------
  Recursive Function postmep_joinCenters(pmepparams, cplane, cent1, depth) &
      Result (objvalue)
    Type(POSTMEP_Params), Intent(inout):: pmepparams
    Integer, Intent(in)                :: cplane, cent1, depth
    Real(kind=RDbl)                    :: objvalue

    Integer         :: npotcenters, i, j, nplane, nextplane, nextcent
    Integer         :: nsorbs
    Integer, Dimension(MAX_CENTERS) :: centers
    Real(kind=RDbl) :: leastobjvalue, pathvalue, dist, nrg1
    Type(VecType), Dimension(3)     :: com1
    Type(VecType)   :: com2, dispvect
    Character(len=strLen) :: temp, diststr

    !** Some sundrys
    nsorbs = pmepparams%nsorbs

    If (cplane == pmepparams%lpno .Or. depth == 0) Then
      objvalue = 0.0_RDbl
      Return
    Else
      Do i=1, nsorbs
        com1(i) = pmepparams%planeno(cplane)%clusters(cent1)%com(i)
      End Do

      nrg1 = &
          stats_getcavg(pmepparams%planeno(cplane)%clusters(cent1)%avgnrg(7))

      !** Get the potential centers in a subsequent plane
      npotcenters = &
          postmep_findcenters(pmepparams, com1(1),nrg1,cplane+1,nplane,centers)

      !** Pick a center that has the lowest objective value from a subsequent
      !** plane
      leastobjvalue = 1000.0_RDbl   ! Some large number
      Do i=1, npotcenters
        pathvalue = postmep_joincenters(pmepparams,nplane,centers(i),depth-1)
        dist = 0.0_Rdbl
        Do j=1, nsorbs
          com2 = pmepparams%planeno(nplane)%clusters(centers(i))%com(j)
          dist = dist + mag(com1(j) - com2)
        End Do
        pathvalue = pathvalue + dist
        If (leastobjvalue > pathvalue) Then
          leastobjvalue = pathvalue
          nextplane = nplane
          nextcent  = centers(i)
        End If
      End Do

      !** Set up the path to the center with the lowest objective value
      If (nplane /= 0) Then
        objvalue = leastobjvalue
        pmepparams%planeno(cplane)%clusters(cent1)%nextplane = nextplane
        pmepparams%planeno(cplane)%clusters(cent1)%nextcent  = nextcent
        pmepparams%planeno(nextplane)%clusters(nextcent)%prevplane = cplane
        pmepparams%planeno(nextplane)%clusters(nextcent)%prevcent = cent1
      Else
        pmepparams%planeno(cplane)%clusters(cent1)%nextplane = 0
        pmepparams%planeno(cplane)%clusters(cent1)%nextcent  = 0
        objvalue = 0.0_RDbl  
      Endif

!!$      Write(*,'(a,4i5,f10.3)') 'cplane, cent1, nextplane, nextcent', &
!!$          cplane, cent1, nextplane, nextcent, objvalue
    End If
  End Function postmep_joinCenters

  !------------------------------------------------------------------------
  ! Find all the centers which lie in the positive direction of "com1"
  ! The routine returns the plane no in "nplane" and the center number
  ! in "centno" where the centers were found.  The no. of centers found
  ! is returned in the function value.  "com1" is the center-of-mass of
  ! the penetrant.
  !------------------------------------------------------------------------
  Integer Function postmep_findCenters(pmepparams, com1, nrg1, cplane, &
      nplane, centers)
    Type(POSTMEP_Params), Intent(inout) :: pmepparams
    Type(VecType), Intent(in)           :: com1
    Real(kind=RDbl), Intent(in)         :: nrg1
    Integer, Intent(in)                 :: cplane
    Integer, Intent(out)                :: nplane
    Integer, Dimension(:), Intent(out)  :: centers

    Integer         :: cent2, next, i, j, planetol, npotcenters
    Type(VecType)   :: com2, normal, dispvect
    Real(kind=RDbl) :: leastdist, dist, direc, nrg2
    Integer, Dimension(MAX_CENTERS) :: nextplanecenters

    !** Some sundry information
    normal = pmepparams%normal
    nplane = 0
    planetol = pmepparams%planetol
    centers  = 0

    !** Find the closest centers in a positive direction.  Only use
    !** the position of penetrant molecule, which we assume is sorbate 1,
    !** to get the potential centers
    Do i=cplane, Min(pmepparams%lpno, cplane+planetol)
      npotcenters = 0 ! No. of potential centers in the next plane
      Do cent2 = 1, pmepparams%planeno(i)%ncenters
        com2 = pmepparams%planeno(i)%clusters(cent2)%com(1)
        nrg2 = stats_getcavg(pmepparams%planeno(i)%clusters(cent2)%avgnrg(7))
        dispvect = com2 - com1
        direc = dispvect*normal
        ! Find the potential centers in the next based on "direc"
        If (direc >= 0.0_RDbl .And. (Abs(nrg1-nrg2) < 2.0_RDbl)) Then
!!$        If ((Abs(nrg1-nrg2) < 5.0_RDbl)) Then
          npotcenters = npotcenters + 1
          centers(npotcenters) = cent2
        End If
      End Do
      If (npotcenters /= 0) Then
        nplane = i
        postmep_findcenters = npotcenters
        Return
      End If
    End Do
    nplane = 0
    postmep_findcenters = 0
    Return
  End Function postmep_findCenters

  !----------------------------------------------------------------
  ! Find continuous paths from the beginning to the end where the end
  ! is defined to within (pmepparams%lpno - planetol)
  !-----------------------------------------------------------------
  Integer Function postmep_findcontpaths(pmepparams, path)
    Type(POSTMEP_Params), Intent(in)  :: pmepparams
    Integer, Dimension(:), Intent(out):: path

    Integer     :: ncenters, j, lastplane, nextplane, nextcent
    Integer     :: planetol, npaths
    
    !** We define the end as any plane number within "planetol" of 
    !** the last plane
    planetol = pmepparams%planetol
    
    !** Find the continuous paths
    ncenters = pmepparams%planeno(1)%ncenters
    npaths = 0
    path   = 0
    Do j=1, ncenters
      lastplane = 1
      nextplane = pmepparams%planeno(1)%clusters(j)%nextplane      
      nextcent  = pmepparams%planeno(1)%clusters(j)%nextcent
      Do 
        If (nextplane == 0) Exit
        lastplane = nextplane
        nextplane =  &
            pmepparams%planeno(lastplane)%clusters(nextcent)%nextplane
        nextcent = pmepparams%planeno(lastplane)%clusters(nextcent)%nextcent
      End Do
      ! See if we made through to the last plane +/- planetol
      If (Abs(lastplane - pmepparams%lpno) < planetol) Then
        npaths = npaths + 1
        path(npaths) = j
      End If
    End Do
    postmep_findcontpaths = npaths
  End Function postmep_findcontpaths

  !---------------------------------------------------------------------
  ! This routine dumps the positions of all the centers to RASMOL file
  !---------------------------------------------------------------------
  Subroutine postmep_dumpCenters(pmepparams)
    Type(POSTMEP_Params), Intent(in) :: pmepparams

    Integer    :: i, j, error, outunitno, inunitno
    Integer    :: npts, nsorbs, sorbno, maxcenters
    Integer    :: lastplane, npaths, startcent
    Type(VecType) :: rp, r, com
    Type(IntVecType) :: cindx
    Character(len=strLen)   :: outxyzfile

    !** Open the xyz file
    outxyzfile = Trim(pmepparams%outbasefilename)//".centers.xyz"
    outunitno  = file_getunit(outxyzfile)
    Open(unit=outunitno, file=outxyzfile, IOSTAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not open file : ", outxyzfile
      Stop
    End If
    
    !** Dump the penetrant molecule
    Do i=1, pmepparams%lpno
      Do j=1, pmepparams%planeno(i)%ncenters
        Write(outunitno,'(a5, a30, 2i4)') "Ar", &
            Trim(vector_display( &
            pmepparams%planeno(i)%clusters(j)%com(1),"f8.3")), i, j
      End Do
      
      !** Dump the other molecules
      Do sorbno=2, pmepparams%nsorbs
        Do j=1, pmepparams%planeno(i)%ncenters
          Write(outunitno,'(a5, a30, 2i4)') "He", &
            Trim(vector_display( &
            pmepparams%planeno(i)%clusters(j)%com(sorbno),"f8.3")), i, j
        End Do
      End Do
    End Do
    Close(outunitno)
  End Subroutine postmep_dumpCenters

  !------------------------------------------------------------------
  ! Generate the energy statistics from each file.  Dump the energy
  ! of each plane by averaging over all the iterations.
  !------------------------------------------------------------------
  Subroutine postmep_dumpAverageEnergy(pmepparams)
    Type(POSTMEP_Params), Intent(in) :: pmepparams

    Integer      :: i, j, startext, endext, outunitno, inunitno, error 
    Integer      :: nrecs, junk
    Character(len=strLen)   :: nrgfilename, outnrgfile

    Real(kind=RDbl)   :: mb1, mb2, mz, b1b2, b1z, b2z, totalnrg, coord
    Type(Statistics)  :: statmb1, statmb2, statmz, statb1b2, statb1z, &
                         statb2z, stattotalnrg
    Character(len=lstrLen) :: junkstr
    Type(VecType)     :: normal

    !** Open the output energy file
    outnrgfile = Trim(pmepparams%outbasefilename)//".cumnrg"
    outunitno  = file_getunit(outnrgfile)
    Open(unit=outunitno, file=outnrgfile, IOSTAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not open file : ", outnrgfile
      Stop
    End If

    !** Cycle through all the input energy files and gather the stats
    Do i=1, pmepparams%lpno
      coord = pmepparams%normal*pmepparams%planeno(i)%annealno(1)%lastpos(1)
      Write(outunitno,'(i4, f8.3, 8f10.3)') i, coord, &
          (stats_getcavg(pmepparams%planeno(i)%avgPlaneEnergy(j)), j=1, 7), &
          stats_getstd(pmepparams%planeno(i)%avgPlaneEnergy(7))
    End Do
    Close(outunitno)
  End Subroutine postmep_dumpAverageEnergy

  !-------------------------------------------------------------------------
  ! Generate the energy statistics from each file.  This routine 
  ! should be called after the continuous paths have been generated.
  ! It then dumps the energy of each center in the path.  The center
  ! energy is obtained by averaging over the different simulated annealings
  ! comprising the center
  !------------------------------------------------------------------------
  Subroutine postmep_dumpPathEnergy(pmepparams, pathno)
    Type(POSTMEP_Params), Intent(in) :: pmepparams
    Integer, Intent(in)              :: pathno

    Integer      :: i, j, outunitno, inunitno, error, cent, plane, lastplane 
    Integer      :: ptno
    Character(len=strLen)   :: nrgfilename, outnrgfile
    Real(kind=RDbl)   :: totalstddev, coord
    Real(kind=RDbl), Dimension(7) :: pairnrg
    Character(len=lstrLen) :: junkstr
    Type(VecType)     :: normal, com

    !** Get the normal to obtain the reaction coordinate
    normal = pmepparams%normal

    !** Open the output energy file
    outnrgfile = Trim(pmepparams%outbasefilename)//".pathnrg"
    outnrgfile = genfilename(outnrgfile, pathno)
    outunitno  = file_getunit(outnrgfile)
    Open(unit=outunitno, file=outnrgfile, IOSTAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not open file : ", outnrgfile
      Stop
    End If

    !** Cycle through the path and dump the energy
    cent  = pmepparams%pathhead(pathno)
    plane = 1
    ptno = 0
    Do 
      If (cent == 0) Exit
      com = pmepparams%planeno(plane)%clusters(cent)%com(1)
      coord = com*normal
      Do j=1, 7
        pairnrg(j) = &
            stats_getcavg(pmepparams%planeno(plane)%clusters(cent)%avgnrg(j))
      End Do
      totalstddev = stats_getstd(&
          pmepparams%planeno(plane)%clusters(cent)%avgnrg(7))
      ptno = ptno + 1
      Write(outunitno, '(2i5, 8f9.3, f7.3)') ptno, plane, coord, &
          pairnrg(1:7), totalstddev
      lastplane = plane
      plane = pmepparams%planeno(lastplane)%clusters(cent)%nextplane
      cent  = pmepparams%planeno(lastplane)%clusters(cent)%nextcent
    End Do
    Close(outunitno)
  End Subroutine postmep_dumpPathEnergy

  !-------------------------------------------------------------
  ! This routine dumps the xyz coordinates of the path "pathno"
  ! which has as its first center "paths(pathno)"
  !-------------------------------------------------------------
  Subroutine postmep_dumpPath(pmepparams, pathno)
    Type(POSTMEP_Params), Intent(in)  :: pmepparams
    Integer, Intent(in)               :: pathno

    Integer   :: i, j, lastplane, plane, cent, outunitno, error 
    Character(len=strLen) :: outpathfile
    Real(kind=RDbl) :: totalnrg, stddev
    Type(VecType)   :: com

    !** Open the output energy file
    outpathfile = Trim(pmepparams%outbasefilename)//".pathstats"
    outpathfile = genfilename(outpathfile, pathno)
    outunitno  = file_getunit(outpathfile)
    Open(unit=outunitno, file=outpathfile, IOSTAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not open file : ", outpathfile
      Stop
    End If

    !** Cycle through the path and dump the coordinates
    cent  = pmepparams%pathhead(pathno)
    plane = 1
    Write(outunitno, '(a25, i5)') "Path No. :", pathno
    Do 
      If (cent == 0) Exit
      com = pmepparams%planeno(plane)%clusters(cent)%com(1)
      totalnrg = stats_getcavg(&
          pmepparams%planeno(plane)%clusters(cent)%avgnrg(7))
      stddev = stats_getstd(&
          pmepparams%planeno(plane)%clusters(cent)%avgnrg(7))
      Write(outunitno, '(i5, a40,2f12.3)') plane, &
          Trim(vector_display(com, "f12.3")), totalnrg, stddev
      lastplane = plane
      plane = pmepparams%planeno(lastplane)%clusters(cent)%nextplane
      cent  = pmepparams%planeno(lastplane)%clusters(cent)%nextcent
    End Do
    Close(outunitno)
  End Subroutine postmep_dumpPath


  !-----------------------------------------------------------------
  ! Generates a file with all the paths after mapping back to the
  ! unit cell coordinates.
  !-----------------------------------------------------------------
  Subroutine postmep_dumpAllPaths(pmepparams, simcell)
    Type(POSTMEP_Params), Intent(in) :: pmepparams
    Type(Simcell_Params), Intent(in) :: simcell

    Integer       :: startext, endext, i, j, sorbno, error, nsorbs 
    Integer       :: outunitno, inunitno, nrecs, junk
    Type(VecType) :: rp, r
    Type(IntVecType) :: cindx
    Character(len=strLen)   :: outxyzfile, avgposfilename
    Real(kind=RDbl)  :: x, y, z

    !** Open the xyz file
    outxyzfile = Trim(pmepparams%outbasefilename)//".allpaths.xyz"
    outunitno  = file_getunit(outxyzfile)
    Open(unit=outunitno, file=outxyzfile, IOSTAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " Could not open file : ", outxyzfile
      Stop
    End If
    
    nsorbs = pmepparams%nsorbs
    !** Cycle through all the input energy files and gather the stats
    Do i=1, pmepparams%lpno
      !**Get the x y z positions
      Do j=1, pmepparams%planeno(i)%nanneals
        !**Generate unit cell coordinates
        Do sorbno=1, nsorbs
          rp = pmepparams%planeno(i)%annealno(j)%lastpos(sorbno)
          Call simcell_pbc(simcell, rp, r, cindx)
        
          !** Write the positions to the rasmol file.  All the Methane atoms
          !** are written as Argon
          If (sorbno == 1) Then
            Write(outunitno,'(a,a35,i10)') "Ar", &
                 Trim(vector_display(r,"f9.3")), i
          Else
            Write(outunitno,'(a,a35,i10)') "Ne", &
                 Trim(vector_display(r,"f9.3")), i
          End If
        End Do
      End Do
    End Do
    Close(outunitno)
  End Subroutine postmep_dumpAllPaths

End Module postmep

