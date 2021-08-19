!----------------------------------------------------------------------------
! This module contains the data and the routines relevant to the simulation 
! cell parameters.  The simulation cells constructed using a multiple, in
! each direction of "fundamental" cells from the fundcell module.  Its 
! overall dimension and shape are also specified using this other data type.
! 
! Needed Improvements:
! 1) could still use more requires statements
!----------------------------------------------------------------------------

Module simcell

  Use defaults, Only: strLen, RDbl, dashedline, zero, radTodeg
  Use utils, Only: split, isfileopen, stripcmnt, filesrchstr, tolower, &
      toupper, getpath
  Use molecules, Only: molecules_gettype, molecules_getfilename, &
      molecules_expand, molecules_reconnect, molecules_dump2xyz
  Use atom, Only: atom_nbonds
  Use file, Only: file_getunit,file_open
  Use fundcell, Only: fundcell_construct, fundcell_getvolume, fundcell_angles,&
      fundcell_getell, fundcell_slant, fundcell_unslant, Fundamental_Cell, &
      fundcell_isortho, fundcell_geteff, fundcell_idstring, fundcell_init, &
      fundcell_display, fundcell_latticevec, fundcell_minwidth, &
      fundcell_snglMinImage, fundcell_multMinImage
  Use vector, Only: VecType, vector_getnormsq, Assignment(=), Operator(+), &
      Operator(-), Operator(*), Operator(/), IntVecType
  Use random, Only: rranf

  Implicit None
  Save

  Private
  Public  :: default_simcell_tag, SimCell_Params, simcell_minimage, &
      simcell_pbc, simcell_getzeoell, simcell_getcellorigin, simcell_maptouc, &
      simcell_getnuc, simcell_checkinit, simcell_getmolectype, simcell_init, &
      simcell_display, simcell_getvolume, simcell_maptosimcell, &
      simcell_uniformRandPt, simcell_getell, simcell_isortho, simcell_minwidth, &
      simcell_getfillsorb, simcell_geteff, simcell_sampleCF, simcell_angles, simcell_pbcpt, simcell_pbcptFund

  Character(len=strLen),Parameter :: default_simcell_tag = &
      "Simulation Cell Information"

  Type SimCell_Params
    Type(IntVecType)                :: size
    Integer                         :: nx,ny,nz
    Integer                         :: ncells
    Type(Fundamental_Cell)          :: fundC, fullC
    !** Type of boundary condition
    Integer :: bcx, bcy, bcz  ! bc = 0 (closed), 1 (periodic), 2 (open)
    Logical                         :: pbcx, pbcy, pbcz
    Integer                         :: molecType
  End Type SimCell_Params

  Interface simcell_pbc
    Module Procedure simcell_pbcpt
    Module Procedure simcell_pbcarray
  End Interface

  Interface simcell_minimage
    Module Procedure simcell_snglMinImage
    Module Procedure simcell_multMinImage
  End Interface

Contains

  !----------------------------------------------------------------------------
  ! Returns the sorbate type number of the simcell
  !----------------------------------------------------------------------------
  Integer Function simcell_getmolectype(simcell)
    Type(SimCell_Params), Intent(In) :: simcell
    simcell_getmolectype = simcell%molecType
  End Function simcell_getmolectype

  !----------------------------------------------------------------------------
  ! Initialize the fields of the simulation cell "simcell"
  !----------------------------------------------------------------------------
  Subroutine simcell_init(simcell, ctrl_filename, optsimcell_tag)
    Implicit None
    Type(SimCell_Params), Intent(OUT)           :: simcell
    Character(*), Intent(IN)                    :: ctrl_filename
    Character(len=strLen), Optional, Intent(IN) :: optsimcell_tag

    Integer                             :: unitno, lineno, i, nfields
    Integer                             :: sorb, unit
    Character(len=strLen)               :: tag, line, name
    Character(len=5*strLen)             :: filewpath
    Type(VecType)                       :: offset
    Real(kind=RDbl), Dimension(3)       :: eff
    Character(len=strLen), Dimension(strLen) :: fields

    If (Present(optsimcell_tag)) Then
      tag = optsimcell_tag
    Else
      tag = default_simcell_tag
    End If

    !** Find the tag in the control file that corresponds to the
    !** id tag for this routine
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
      Open(unit=unitno, file=ctrl_filename, status='old')
    Endif

    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          "Could not find the tag ", tag, "in the control file"
      Stop
    End If

    !** Read information from the control file
    !** Retrieve the fundamental cell information
    !** This may be of two forms: (1) a name of a sorbate from
    !** whose file the fundcell info will be read, (2) the name
    !** of an external file from which the fundcell info will
    !** be read, or (3) the tag NULL, which indicates the fundcell
    !** info will immediately follow.
    Read(unitno,'(a)') line
    line = stripcmnt(line)
    nfields = split(line,fields)

    !** Figure out what to fill the simulation cell with
    Select Case (toupper(Trim(fields(1))))

    Case ("NULL")
      !** There is no fill sorbate. It is an empty unit cell.
      !** Read in the parameters from the control file.
      unit = unitno
      simcell%molecType= -1
    Case Default
      !** Check to see if the name listed is a sorbate type
      name = fields(1)
      
      sorb = molecules_gettype(trim(name))

      If (sorb /= 0) Then
        !** It's a molecule. Set the name and type in the
        !** simcell variable. Get the molecule filename
        simcell%molecType = sorb
        name = molecules_getfilename(sorb)
      Else
        !** It's not a sorbate type. Set the appropriate variables
        !** to null values
        simcell%molecType = -1
      End If

      !** If it wasn't a molecule, we have the filename from 
      !** the control file. If it was a molecule, we now have 
      !** the molecule filename. Open the file and get the
      !** funcell information.

      unit = isfileopen(name)
      If (unit < 0) Then
        !** File isn't open. Open it.
        Write(filewpath,'(2a)') Trim(getpath('MOLSDIR')),Trim(name)
        unit = file_open(filewpath,110)
        !** Find the proper line to read from
        lineno = filesrchstr(unit,fundcell_idstring,line,.True.)
      End If

    End Select

    !** Pass the unit number and data structure to fundcell
    !** so it can initialize itself.
    Call fundcell_init(simcell%fundC,unit)

    !** Retrieve information that will allow us to build the simulation
    !** Read in the number of unit cells
    Read(unitno,*) (simcell%size%comp(i),i=1,3)
    simcell%nx = simcell%size%comp(1)
    simcell%ny = simcell%size%comp(2)
    simcell%nz = simcell%size%comp(3)
    simcell%ncells = simcell%size%comp(1)*simcell%size%comp(2)*&
        simcell%size%comp(3)

    !** Determine the periodicity
    simcell%pbcx = .False.
    simcell%pbcy = .False.
    simcell%pbcz = .False.
    Read(unitno,*) simcell%bcx, simcell%bcy, simcell%bcz

    !** If bc is 1, then it is periodic. If bc is 0 or 2, it is
    !** either closed (0) or open (2). Code must be added to handle
    !** these cases.
    If (simcell%bcx == 1) simcell%pbcx = .True.
    If (simcell%bcy == 1) simcell%pbcy = .True.
    If (simcell%bcz == 1) simcell%pbcz = .True.

    !** Now that the fundcell is initialized, we must
    !** construct the simcell from the fundamental cell
    Call fundcell_construct(simcell%fundC,simcell%size,simcell%fullC)

    !** With the fundamental cell initialized, we must 
    !** check to see if it is a molecule type and initialize any special
    !** things associated with it.
    If (simcell%molecType > 0) Then
      !** expand the contents of the fundcell to fill the entire simcell
      !** Get the size of the fundamental cell
      eff = simcell_geteff(simcell,.True.)

      !** Set the offset, which is the length of the fundcell in the
      !** x direction.
      offset = fundcell_latticevec(simcell%fundC,1)
      Call molecules_expand(simcell%molecType,offset,simcell%nx-1,.False.)

      !** Set the offset, which is the length of the fundcell in the
      !** y direction.
      offset = fundcell_latticevec(simcell%fundC,2)
      Call molecules_expand(simcell%molecType,offset,simcell%ny-1,.False.)

      !** Set the offset, which is the length of the fundcell in the
      !** z direction.
      offset = fundcell_latticevec(simcell%fundC,3)
      Call molecules_expand(simcell%molecType,offset,simcell%nz-1,.False.)

      !** Regenerate the connections list taking into account the 
      !** periodic boundary conditions. This must be done manually
      !** since the PBC routines are contained in this module.
      Call molecules_reconnect(simcell%molecType,simcell%fullC)

!      Call molecules_dump2xyz(simcell%molecType,'full.xyz','pretty please')
    End If
    Close(unitno)
  End Subroutine simcell_init

  !----------------------------------------------------------------------------
  ! Writes a sample of the control file to the unit unitno
  !----------------------------------------------------------------------------
  Subroutine simcell_sampleCF(unitno)
    Integer, Intent(in) :: unitno

    Write(unitno,'(a)') default_simcell_tag
    Write(unitno,'(a,t30,a)') 'Character', &
        '# Fundamental cell type info (filename or zeolite)'
    Write(unitno,'(a,t30,a)') 'Integer,Integer,Integer', &
        '# Number of unit cells in x,y,z'
    Write(unitno,'(a,t30,a)') '[012],[012],[012]', &
        '# Simcell Boundaries, 0=Closed, 1=Periodic, 2=Open'
  End Subroutine simcell_sampleCF
  
  !----------------------------------------------------------------
  ! Returns the number of unit cells in the specified direction.
  ! If no direction is specified, returns the total number
  ! Requires: simcell -- simulation cell parameters
  !           direct -- optional direction
  !----------------------------------------------------------------
  Integer Function simcell_getnuc(simcell,direct)
    Type(SimCell_Params), Intent(In) :: simcell
    Character, Intent(In), Optional  :: direct

    If (Present(direct)) Then
      If (tolower(direct) == 'x') simcell_getnuc = simcell%nx
      If (tolower(direct) == 'y') simcell_getnuc = simcell%ny
      If (tolower(direct) == 'z') simcell_getnuc = simcell%nz
    Else
      simcell_getnuc = simcell%ncells
    End If

  End Function simcell_getnuc


  !----------------------------------------------------------------------
  ! Returns a unifrom random vector within the simcell, in the 
  ! fundamental unitcell. Assumes that the space is cartesian/orthogonal.
  ! Requires:  simcell -- simulation cell data structure
  !----------------------------------------------------------------------
    Type(VecType)Function simcell_uniformRandPt(simcell)
    Type(SimCell_Params), Intent(In) :: simcell

    Real(kind=RDbl), Dimension(3)    :: ell
    Type(VecType)  rvec, v1, v2, v3

    If(simcell%fundC%orthoflag) Then
    ell = fundcell_getell(simcell%fullC)
    
    !** assumes orthogonality
    simcell_uniformRandPt%comp(1) = ell(1) * rranf()
    simcell_uniformRandPt%comp(2) = ell(2) * rranf()
    simcell_uniformRandPt%comp(3) = ell(3) * rranf()

    Else

    v1=simcell%fullC%lvec(1)*rranf()
    v2=simcell%fullC%lvec(2)*rranf()
    v3=simcell%fullC%lvec(3)*rranf()
    rvec=(v1+v2)+v3

    simcell_uniformRandPt%comp(1) = rvec%comp(1) 
    simcell_uniformRandPt%comp(2) = rvec%comp(2) 
    simcell_uniformRandPt%comp(3) = rvec%comp(3) 

     End IF

  End Function simcell_uniformRandPt

  !----------------------------------------------------------------
  ! Returns the volume of the simulation cell "simcell"
  !----------------------------------------------------------------
  Real(kind=Rdbl) Function simcell_getvolume(simcell)
    Type(SimCell_Params), Intent(in) :: simcell
    
    simcell_getvolume = fundcell_getvolume(simcell%fullC)
  End Function simcell_getvolume

  !----------------------------------------------------------------------------
  ! Returns true if the simcell is orthorhombic
  !----------------------------------------------------------------------------
  Logical Function simcell_isortho(simcell)
    Type(SimCell_Params), Intent(In) :: simcell
    simcell_isortho = fundcell_isortho(simcell%fundC)
  End Function simcell_isortho

  !----------------------------------------------------------------------------
  ! Get the molecule type in the simcell
  !----------------------------------------------------------------------------
  Integer Function simcell_getfillsorb(simcell)
    Type(SimCell_Params), Intent(In) :: simcell
    simcell_getfillsorb = simcell%molecType
  End Function simcell_getfillsorb

  !----------------------------------------------------------------
  ! This function returns the origin of the unit cell (cx, cy, cz)
  ! making up the simulation cell.  The first unit cell is (0, 0, 0)
  !----------------------------------------------------------------
  Type(VecType) Function simcell_getcellorigin(simcell, cx, cy, cz)
    Type(SimCell_Params), Intent(in) :: simcell
    Integer, Intent(in)   :: cx, cy, cz
    Type(VecType)  rvec, v1, v2, v3

    Real(kind=RDbl)      :: x, y, z
    Real(kind=RDbl), Dimension(3)    :: ell

    ell = fundcell_getell(simcell%fundC)
    If(simcell%fundC%orthoflag) Then
    x = ell(1)*cx
    y = ell(2)*cy
    z = ell(3)*cz
    Else
    v1 = simcell%fundC%lvec(1)*cx
    v2 = simcell%fundC%lvec(2)*cy
    v3 = simcell%fundC%lvec(3)*cz
    rvec = (v1+v2)+v3
    x = rvec%comp(1)
    y = rvec%comp(2)
    z = rvec%comp(3)
    End If

    simcell_getcellorigin = (/x, y, z/)
    Return
  End Function simcell_getcellorigin

  !----------------------------------------------------------------------------
  ! Get the edge length of the simulation cell. If fund is passed as
  ! TRUE, return the fundamental cell edge lengths.
  !----------------------------------------------------------------------------
  Function simcell_getell(scell,fund)
    Type(SimCell_Params), Intent(in) :: scell
    Logical, Intent(in), Optional    :: fund
    Logical :: getfund
    Real(kind=RDbl), Dimension(3) :: simcell_getell 
    
    !** Check to see if the optional parameter has passed
    getfund = .False.
    If (Present(fund)) getfund = fund

    !** Get the fundamental cell length
    If (getfund) Then
      !** The user has requested the fundamental cell info
      simcell_getell = fundcell_getell(scell%fundC)
    Else
      !** Return the full simulation cell info
      simcell_getell = fundcell_getell(scell%fullC)
    End If

  End Function simcell_getell

  !----------------------------------------------------------------------------
  ! Returns the cell angles
  ! Requires: simcell -- simulation cell parameters
  !           deg_flag -- flags that angles if should be in degrees
  !----------------------------------------------------------------------------
  Function simcell_angles(simcell,deg_flag)
    Type(SimCell_Params), Intent(In) :: simcell
    Logical, Intent(In)              :: deg_flag
    Real(kind=RDbl), Dimension(3)    :: simcell_angles

    simcell_angles = fundcell_angles(simcell%fundC,deg_flag)

  End Function simcell_angles

  !----------------------------------------------------------------------------
  ! Converts a passed coordinate in the simcell to fractional coordinates
  ! Requires: simcell -- simulation cell parameters
  !           coord -- xyz coordinates of a point in the simcell
  !----------------------------------------------------------------------------
  Type(VecType) Function simcell_getfrac(simcell,coord)
    Type(SimCell_Params), Intent(In) :: simcell
    Type(VecType), Intent(In)        :: coord

    Real(kind=RDbl), Dimension(3)  :: ell

    ell = simcell_getell(simcell)

    simcell_getfrac%comp(1) = coord%comp(1)/ell(1)
    simcell_getfrac%comp(2) = coord%comp(2)/ell(2)
    simcell_getfrac%comp(3) = coord%comp(3)/ell(3)

  End Function simcell_getfrac

  !----------------------------------------------------------------------------
  ! Get the effective edge length of the simcell. If the optional
  ! argument FUND is passed as .True., it returns the edge lengths of the
  ! fundamental cell.
  !----------------------------------------------------------------------------
  Function simcell_geteff(simcell,fund)
    Type(SimCell_Params), Intent(In) :: simcell
    Logical, Optional, Intent(In) :: fund
    Real(Kind=RDbl), Dimension(3) :: simcell_geteff
    Logical :: getfund

    getfund = .False.
    If (Present(fund)) getfund = fund

    If (getfund) Then
      simcell_geteff = fundcell_geteff(simcell%fundC)
    Else
      simcell_geteff = fundcell_geteff(simcell%fullC)
    End If
  End Function simcell_geteff

  !--------------------------------------------------
  ! Get the edge length of the zeolite unit cell 
  ! making up the simulation cell
  !--------------------------------------------------
  Function simcell_getzeoell(simcell)
    Type(SimCell_Params), Intent(in) :: simcell
    Real(kind=RDbl), Dimension(3) :: simcell_getzeoell 
    
    simcell_getzeoell = fundcell_getell(simcell%fundC)
  End Function simcell_getzeoell

  !-----------------------------------------------------------------
  ! Get the minimum width of the simulation cell.  This is needed
  ! for checking cutoffs.  
  ! Requires:  simcell -- simulation cell data type
  !-----------------------------------------------------------------
  Real(kind=RDbl) Function simcell_minwidth(simcell)
    Type(SimCell_Params), Intent(In) :: simcell
    
    simcell_minwidth = fundcell_minwidth(simcell%fullC)
  End Function simcell_minwidth

  !----------------------------------------------------------------------------
  ! Does periodic boundary conditions on the coordinate "coord".  Returns
  ! the vector "pbccoord" which is within the unit cell. The optional
  ! argument "optcellindx" is the number of simulation cells the initial
  ! point was away from the principal simulation cell.
  !----------------------------------------------------------------------------
  Subroutine simcell_pbcpt(simcell, coord, pbccoord, optcellindx)
    Implicit None
    Type(SimCell_Params), Intent(in)   :: simcell
    Type(VecType), Intent(in)          :: coord
    Type(VecType), Intent(out)         :: pbccoord
    Type(IntVecType), Optional, Intent(out):: optcellindx
    Type(VecType)                      :: temp
    Real(kind=RDbl), Dimension(3)      :: coordarray
    Integer, Dimension(3)              :: indx
    Integer                            :: i

    coordarray = fundcell_slant(simcell%fullC,coord)
    Do i=1, 3
      indx(i) = Floor(coordarray(i)/simcell%fullC%ell(i))
      coordarray(i) = coordarray(i) -  indx(i)*simcell%fullC%ell(i)
    Enddo
    ! Convert coordarray to a vector and then unslant it
    temp = coordarray
    pbccoord = fundcell_unslant(simcell%fullC,temp)

    If (Present(optcellindx)) Then
      optcellindx = indx
    End If
    Return
  End Subroutine simcell_pbcpt

  !----------------------------------------------------------------------------
  ! Does periodic boundary conditions on the coordinate "coord".  Returns
  ! the vector "pbccoord" which is within the FUNDAMENTAL unit cell. The optional
  ! argument "optcellindx" is the number of simulation cells the initial
  ! point was away from the principal simulation cell.
  !----------------------------------------------------------------------------
  Subroutine simcell_pbcptFund(simcell, coord, pbccoord, optcellindx)
    Implicit None
    Type(SimCell_Params), Intent(in)   :: simcell
    Type(VecType), Intent(in)          :: coord
    Type(VecType), Intent(out)         :: pbccoord
    Type(IntVecType), Optional, Intent(out):: optcellindx
    Type(VecType)                      :: temp
    Real(kind=RDbl), Dimension(3)      :: coordarray
    Integer, Dimension(3)              :: indx
    Integer                            :: i

    coordarray = fundcell_slant(simcell%fundC,coord)
    Do i=1, 3
      indx(i) = Floor(coordarray(i)/simcell%fundC%ell(i))
      coordarray(i) = coordarray(i) -  indx(i)*simcell%fundC%ell(i)
    Enddo
    ! Convert coordarray to a vector and then unslant it
    temp = coordarray
    pbccoord = fundcell_unslant(simcell%fundC,temp)

    If (Present(optcellindx)) Then
      optcellindx = indx
    End If
    Return
  End Subroutine simcell_pbcptFund

  !----------------------------------------------------------------------------
  ! Does periodic boundary conditions on an array of coordinates "coords"
  !----------------------------------------------------------------------------
  Subroutine simcell_pbcarray(simcell, coords, pbccoords, cellindx)
    Type(SimCell_Params), Intent(in)         :: simcell
    Type(VecType), Dimension(:), Intent(in)  :: coords
    Type(VecType), Dimension(:), Intent(out) :: pbccoords
    Type(IntVecType), Dimension(:), Optional, Intent(out):: cellindx
    Type(VecType)    :: temp
    Real(kind=RDbl), Dimension(3)      :: coordarray
    Integer, Dimension(3)              :: indx
    Integer   :: nelem, elem, i
    
    nelem = Size(coords,1)
    Do elem = 1,nelem
      
      coordarray = fundcell_slant(simcell%fullC,coords(elem))
      Do i = 1,3
        indx(i) = Floor(coordarray(i)/simcell%fullC%ell(i))
        coordarray(i) = coordarray(i) -  indx(i)*simcell%fullC%ell(i)
      End Do

      temp = coordarray
      pbccoords(elem) = fundcell_unslant(simcell%fullC,temp)
      If (Present(cellindx)) Then
        cellindx(elem) = indx
      End If
    End Do

  End Subroutine simcell_pbcarray

  !------------------------------------------------------------------
  ! Get the minimum image distance between two position vectors
  ! Requires:  simcell -- simulation cell information
  !            atvec1 -- 1st position vector
  !            atvec2 -- 2nd position vector
  !            sepvec -- minimum image separation vector
  !            dist2 -- square of minimum image distance 
  !------------------------------------------------------------------
  Subroutine simcell_snglMinImage(simcell,atvec1,atvec2,sepvec,dist2)
    Type(SimCell_Params), Intent(In)       :: simcell
    Type(VecType), Intent(In)              :: atvec1,atvec2
    Type(VecType), Intent(Out)             :: sepvec
    Real(kind=RDbl), Optional, Intent(Out) :: dist2

    If (Present(dist2)) Then
      Call fundcell_snglMinImage(simcell%fullC,atvec1,atvec2,sepvec,dist2)
    Else
      Call fundcell_snglMinImage(simcell%fullC,atvec1,atvec2,sepvec)
    End If

  End Subroutine simcell_snglMinImage

  !----------------------------------------------------------------------------
  ! Returns the minimum image distances between all coords in atvecs1 and 
  ! atvecs2.  Separation vectors are calculated as (atom1 - atom2).  
  ! sepvecs is indexed so (atom1,atom2) is the distance between point 
  ! atvecs1(atom1) and atvecs2(atom2).  Optionally returns the distance
  ! squared in dists2.  
  ! Requires:  simcell -- simulation cell information
  !            atvecs1 -- 1st position vectors
  !            atvecs2 -- 2nd position vectors
  !            sepvecs -- minimum image separation vectors
  !            dists2 -- square of minimum image distances 
  !----------------------------------------------------------------------------
  Subroutine simcell_multMinImage(simcell,atvecs1,atvecs2,sepvecs,dists2)
    Type(SimCell_Params), Intent(in)                       :: simcell
    Type(VecType), Intent(In), Dimension(:)                :: atvecs1
    Type(VecType), Intent(In), Dimension(:)                :: atvecs2
    Type(VecType), Intent(Out), Dimension(:,:)             :: sepvecs
    Real(kind=RDbl), Dimension(:,:), Intent(Out), Optional :: dists2    

    If (Present(dists2)) Then
      Call fundcell_multMinImage(simcell%fullC,atvecs1,atvecs2,sepvecs,dists2)
    Else
      Call fundcell_multMinImage(simcell%fullC,atvecs1,atvecs2,sepvecs)
    End If

  End Subroutine simcell_multMinImage

  !----------------------------------------------------------------------------
  ! Translate a vector to a specific unit cell
  !----------------------------------------------------------------------------
  Type(VecType) Function simcell_translate(simcell,vec,uc)
    Implicit None
    Type(SimCell_Params), Intent(IN)   :: simcell
    Type(VecType), Intent(IN)          :: vec
    Type(IntVecType), Intent(IN)       :: uc

    Type(VecType)                      :: dummy
    Real(kind=RDbl), Dimension(3)      :: ell

!    dummy = vec - simcell%zeolite%fcell%origin
!    dummy = fundcell_slant(simcell%zeolite%fcell,dummy)

    ell = fundcell_getell(simcell%fullC)

    dummy = vec - simcell%fundC%origin
    dummy = fundcell_slant(simcell%fundC,dummy)

    dummy%comp(1) = dummy%comp(1) + uc%comp(1)*ell(1)
    dummy%comp(2) = dummy%comp(2) + uc%comp(2)*ell(2)
    dummy%comp(3) = dummy%comp(3) + uc%comp(3)*ell(3)
    
    dummy = fundcell_unslant(simcell%fundC,dummy)
    simcell_translate = dummy + simcell%fundC%origin    
    
  End Function simcell_translate

  !----------------------------------------------------------------------------
  ! Translate a vector to the primary unit cell
  ! Requires:  simcell -- the simulation cell
  !            vec -- the vector to be transformed
  !            cvec -- optional (output only) continuation vector
  ! NOTE: strongly consider just calling fundcell_maptocell
  !----------------------------------------------------------------------------
  Type(VecType) Function simcell_maptouc(simcell,vec,cvec)
    Type(SimCell_Params), Intent(IN)         :: simcell
    Type(VecType), Intent(IN)                :: vec
    Type(IntVecType), Intent(Out), Optional  :: cvec

    Integer                            :: i
    Type(VecType)                      :: dummy

    dummy = vec - simcell%fundC%origin
    dummy = fundcell_slant(simcell%fundC,dummy)
    dummy%comp = dummy%comp*simcell%fundC%ell_inv

    If (Present(cvec)) Then
      !** also get the continuation vectors
      Do i = 1,3
        cvec%comp(i) = floor(dummy%comp(i))
        dummy%comp(i) = dummy%comp(i) - cvec%comp(i)
      End Do
    Else
      !** simply translate back to the primary unit cell
      Do i = 1,3
        If (dummy%comp(i) >= 0.0d0) Then
           dummy%comp(i) = mod(dummy%comp(i), 1.0d0)
        Else
           dummy%comp(i) = mod(dummy%comp(i), 1.0d0) + 1.0d0
        End If
      End Do
    End If

    dummy%comp = dummy%comp*simcell%fundC%ell
    dummy = fundcell_unslant(simcell%fundC,dummy)
    dummy = dummy + simcell%fundC%origin    

    simcell_maptouc = dummy
  End Function simcell_maptouc

  !----------------------------------------------------------------------------
  ! Translate a vector to the primary simulation unit cell
  ! Requires:  simcell -- the simulation cell
  !            vec -- the vector to be transformed
  !            cvec -- optional (output only) continuation vector
  !----------------------------------------------------------------------------
  Type(VecType) Function simcell_maptosimcell(simcell,vec,cvec)
    Type(SimCell_Params), Intent(In)         :: simcell
    Type(VecType), Intent(In)                :: vec
    Type(IntVecType), Intent(Out), Optional  :: cvec

    Integer                            :: i
    Type(VecType)                      :: dummy

    dummy = vec - simcell%fundC%origin
    dummy = fundcell_slant(simcell%fundC,dummy)
    dummy%comp = dummy%comp*simcell%fullC%ell_inv

    If (Present(cvec)) Then
      !**also get the continuation vectors
      Do i = 1,3
        cvec%comp(i) = floor(dummy%comp(i))
        dummy%comp(i) = dummy%comp(i) - cvec%comp(i)
      End Do
    Else 
      !**simply translate back to the primary unit cell
      Do i = 1,3
        If (dummy%comp(i) >= 0.0d0) Then
          dummy%comp(i) = Mod(dummy%comp(i), 1.0d0)
        Else
          dummy%comp(i) = Mod(dummy%comp(i), 1.0d0) + 1.0d0
        End If
      End Do
    End If

    dummy%comp = dummy%comp*simcell%fullC%ell
    dummy = fundcell_unslant(simcell%fundC,dummy)

    simcell_maptosimcell = dummy
    
  End Function simcell_maptosimcell

  !----------------------------------------------------------------------------
  ! Checks to see if the atom structure has been initialized.  
  ! It does this by checking to see that the zeolite has been initialized
  !----------------------------------------------------------------------------
  Logical Function simcell_checkinit(simcell)
    Type(SimCell_Params), Intent(IN)   :: simcell

!!$    !MDEBUG
!!$    !** What do we check to see if everything is set?
!!$    If (Associated(simcell%zeolite)) Then
!!$      simcell_checkinit = .TRUE.
!!$    Else
!!$      simcell_checkinit = .FALSE.
!!$    End If
    simcell_checkinit = .TRUE.
    Return
  End Function simcell_checkinit

  !----------------------------------------------------------------------------
  ! Display the fields of the simulation cell.
  ! Requires:  simcell -- the simulation cell
  !            indent -- no. of spaces from the left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------------
  Subroutine simcell_display(simcell,indent,unit)
    Type(SimCell_Params), Intent(In)   :: simcell
    Integer, Intent(In)                :: indent
    Integer, Intent(In)                :: unit

    Character(len=indent)      :: blank

    blank = Repeat(' ',indent)
    
    Write(unit,'(2a)') blank,'The SIMCELL structure:'
    Write(unit,'(2a,3i4)') blank,'nx, ny, nz            : ', &
        simcell%nx,simcell%ny,simcell%nz
    Write(unit,'(2a,i4)') blank, 'ncells                : ', & 
        simcell%ncells
    Call fundcell_display(simcell%fullC,indent,unit)
    
  End Subroutine simcell_display

End Module simcell




