!---------------------------------------------------------------------------
! This module contains the data structure and the code to read in the 
! pmap/bias map headers.  The header contains all the information needed
! to read the map format and interpret/use its information.
!
! Needed Improvements:
! 1) the current header structure is out-dated, should improve, see start
!---------------------------------------------------------------------------

Module mapheader

  Use defaults, Only: RDbl, RSgl, dashedline
  Use utils, Only: allocerrdisplay, deallocerrdisplay, cleanstring

  Implicit None
  Save

  Private
  Public :: Header_Info, Map_Zeolite_Params, mapheader_display, &
      mapheader_getheader, mapheader_cleanup, mapheader_getIndices, &
      mapheader_writeHeader, mapheader_sizes, mapheader_boxinfo

  Type Map_Zeolite_Params
    Character(len=25)                      :: atom_name
    Character(len=2)                       :: symbol
    Real                                   :: A12, B6, C2, D0  !pot params
    Real                                   :: charge
    Real                                   :: rcuthi, rcutlo
  End Type Map_Zeolite_Params
  
  Type Header_Info
     ! nbrx, nbry, nbrz are the number of nodes NOT cubelets
    Integer                                :: nbrx,nbry,nbrz, nfsize
    Integer                                :: ncubelets
    Real(kind=RDbl), Dimension(:), Pointer :: xgrid, ygrid, zgrid
    Real(kind=RDbl)                        :: stepx, stepy, stepz
    Real(kind=RDbl), Dimension(3)          :: eff
    Type(Map_Zeolite_Params), Dimension(:), Pointer :: zatom    
    Real(kind=RSgl)                        :: pcut_off
    Integer                                :: nzatoms
    Character                              :: mapMolecule
    Character                              :: probeMolecule
    Character(len=60)           :: mapheader_identifier
  End Type Header_Info

  !** This is the new map header type. It needs to be completed.
!!$  Type MapHeader
!!$    Integer    :: nx, ny, nz    ! Number of points in x, y, z directions
!!$    Integer    :: ncubelets     ! Total number of cublets
!!$    Integer    :: ntab          ! Number of TABULATED points
!!$    ! xgrid, ygrid, zgrid,  aren't stored in the file
!!$    Real(kind=RDbl), Dimension(:), Pointer :: xgrid, ygrid, zgrid
!!$  End Type MapHeader

  Interface mapheader_writeHeader
    Module Procedure mapheader_makeOldHeader
  End Interface

Contains
  !--------------------------------------------------------
  ! Gets the pmap/bmap file header information
  ! Requires:  header -- the header information type
  !            eff -- unit cell edge lengths
  !--------------------------------------------------------
  Subroutine mapheader_getheader(header, eff, unitno)
    Type(Header_Info), Intent(InOut)          :: header
    Real(kind=RDbl), Dimension(:), Intent(In) :: eff
    Integer, Intent(In)                       :: unitno

    Type Zeolite_Atom_Types
      Integer                   :: o,si,al,h,p,b,ti,na,ca,k,other
    End Type Zeolite_Atom_types
    Type(Zeolite_Atom_types)    :: zeoatom_types
    Integer                     :: error, intjunk, i
    Integer                     :: nbrx,nbry,nbrz,nfsize
    Real(kind=RSgl)             :: realjunk
    Real(kind=RDbl)             :: stepx, stepy, stepz
    Character(len=60)           :: string60
    Character(len=40)           :: string40
    Character(len=25)           :: string25

    !** Now lets start the reading process
    Rewind(unitno)
    Read(unitno, IOSTAT=error) string40
    If (error /= 0) Then
      Write(0,'(1x,2a,i4,a,i3)') __FILE__," : ",__LINE__, &
          " does not seem like a pmap file, iostat error ",error
      Stop
    End If

    !** make sure the map has a header
    If (Trim(string40) /= 'file contains header') Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
           ' This file does not contain a header'
      Stop
    End If

    Read(unitno) string60
    header%mapheader_identifier= Trim(string60)
    Read(unitno) string25
    Read(unitno) realjunk
    Read(unitno) header%nzatoms
    Read(unitno) zeoatom_types

    !** Allocate the memory for zatoms
    Allocate(header%zatom(header%nzatoms), STAT=error)
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          "Could not allocate memory for 'zatoms'"
      Stop
    End If


    !** read the zeolite atoms
    Do i = 1,header%nzatoms
      Read(unitno) header%zatom(i)%atom_name, header%zatom(i)%symbol
      Read(unitno) header%zatom(i)%charge
      Read(unitno) header%zatom(i)%A12, header%zatom(i)%B6, &
          header%zatom(i)%C2, header%zatom(i)%D0
      ! header was change silightly in 2003-Jan to be consistent with 
      ! New pmap code
      IF(cleanstring(string60) == 'Old type header 2003') THEN 
         Read(unitno) header%zatom(i)%rcuthi, header%zatom(i)%rcutlo 
      ELSE
         header%zatom(i)%rcuthi = 0.
         header%zatom(i)%rcutlo = 0.
      END IF
    End Do
         
    Read(unitno) intjunk
    Read(unitno) header%pcut_off
    Read(unitno) realjunk
    Read(unitno) realjunk
    Read(unitno) intjunk
    Read(unitno) realjunk
    Read(unitno) nbrx,nbry,nbrz
    Read(unitno) nfsize

    header%nbrx = nbrx
    header%nbry = nbry
    header%nbrz = nbrz
    header%ncubelets = (nbrx-1)*(nbry-1)*(nbrz-1)
    header%nfsize = nfsize
    header%eff    = eff

    !** Set up the tabulation grid
    stepx = eff(1)/(nbrx - 1.0_RDbl)
    stepy = eff(2)/(nbry - 1.0_RDbl)
    stepz = eff(3)/(nbrz - 1.0_RDbl)
    header%stepx = stepx
    header%stepy = stepy
    header%stepz = stepz

    !** Allocate xgrid, ygrid and zgrid
    Allocate(header%xgrid(nbrx), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'header%xgrid')   
    Allocate(header%ygrid(nbry), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'header%ygrid')   
    Allocate(header%zgrid(nbrz), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'header%zgrid')   

    !** Set up the grids
    Do i = 1,nbrx
      header%xgrid(i) = (i-1) * stepx
    End Do
    Do i = 1,nbry
      header%ygrid(i) = (i-1) * stepy
    End Do
    Do i = 1,nbrz
      header%zgrid(i) = (i-1) * stepz
    End Do

  End Subroutine mapheader_getheader

  !----------------------------------------------------------------------------
  ! Writes the old style header to the beginning of the given unitno
  !----------------------------------------------------------------------------
  Subroutine mapheader_makeOldHeader(header, unitno)
    Type(Header_Info), Intent(inout) :: header
    Integer, Intent(in)   :: unitno

    Type Zeolite_Atom_Types
      Integer                   :: o,si,al,h,p,b,ti,na,ca,k,other
    End Type Zeolite_Atom_types
    Type(Zeolite_Atom_types)    :: zeoatom_types
    Integer                     :: intjunk, i
    Real(kind=RSgl)             :: realjunk
    Character(len=60)           :: string60
    Character(len=40)           :: string40
    Character(len=25)           :: string25
!!$    Character(len=150)          :: file

    !** Rewind the file to ensure we are at the beginning since the
    !** header MUST go there
    Rewind(unitno)

    intjunk = 0
    realjunk = 0.0

    !** Write the stupid line that tells you yes, there is a header
    Write(string40,'(a20,20x)') 'file contains header'
    Write(unitno) string40

    !** Write the description of the file
    Write(string60,'(a20,40x)') 'Old type header 2003'
    Write(unitno) string60

    Write(string25,'(a5,20x)') 'probe'
    Write(unitno) string25
    Write(unitno) realjunk

    Write(unitno) header%nzatoms
    Write(unitno) zeoatom_types

    Do i = 1, header%nzatoms
      Write(unitno) header%zatom(i)%atom_name, header%zatom(i)%symbol
      Write(unitno) header%zatom(i)%charge
      Write(unitno) header%zatom(i)%A12, header%zatom(i)%B6, &
          header%zatom(i)%C2, header%zatom(i)%D0
      write(unitno) header%zatom(i)%rcuthi, header%zatom(i)%rcutlo 
   End Do
         
    Write(unitno) intjunk
    Write(unitno) header%pcut_off
    Write(unitno) realjunk
    Write(unitno) realjunk       
    Write(unitno) intjunk
    Write(unitno) realjunk
    Write(unitno) header%nbrx,header%nbry,header%nbrz
    Write(unitno) header%nfsize

  End Subroutine mapheader_makeOldHeader

  !----------------------------------------------------------------------------
  ! Gets the critical sizes from the header
  ! Requires:  header -- header to take information from
  !            nbrx -- number of grid points in x-direction
  !            nbry -- number of grid points in y-direction
  !            nbrz -- number of grid points in z-direction
  !            nfsize -- number of tabulated points
  !----------------------------------------------------------------------------
  Subroutine mapheader_sizes(header,nbrx,nbry,nbrz,nfsize)
    Type(Header_Info), Intent(In) :: header
    Integer, Intent(Out)          :: nbrx,nbry,nbrz,nfsize

    nbrx   = header%nbrx
    nbry   = header%nbry
    nbrz   = header%nbrz
    nfsize = header%nfsize

  End Subroutine mapheader_sizes

  !----------------------------------------------------------------------------
  ! Returns the box increments for a potential map if they are available.
  ! Requires:  header -- header to take information from
  !            boxincrements -- x,y,z increments for potential map
  !            boxsteps -- number of boxes in x,y,z directions
  !----------------------------------------------------------------------------
  Subroutine mapheader_boxinfo(header,boxincrements,boxsteps)
    Type(Header_Info), Intent(In)               :: header
    Real(Kind=RDbl), Dimension(3), Intent(Out)  :: boxincrements
    Integer, Dimension(3), Intent(Out)          :: boxsteps

    boxincrements(1) = header%stepx
    boxincrements(2) = header%stepy
    boxincrements(3) = header%stepz
    boxsteps(1) = header%nbrx
    boxsteps(2) = header%nbry
    boxsteps(3) = header%nbrz

  End Subroutine mapheader_boxinfo

  !----------------------------------------------------------------------------
  ! Calculates the indices based on the passed coordinate and the
  ! step size stored in the header
  !----------------------------------------------------------------------------
  Subroutine mapheader_getIndices(header,x,y,z,indices,nearest)
    Type(Header_Info), Intent(In) :: header
    Real(Kind=RDbl), Intent(In)   :: x,y,z
    Integer, Dimension(:), Intent(Out) :: indices
    Logical, Optional, Intent(In) :: nearest
    Logical :: useNint

    !** Check to see if we want the INT or NINT (nearest integer) function
    useNint = .False.
    If (Present(nearest)) useNint = nearest

    !** Calculate the indices using the step size and xyz coordinates
    If (useNint) Then
      indices(1) = Int(Nint(x/header%stepx+1.0_RDbl))
      indices(2) = Int(Nint(y/header%stepy+1.0_RDbl))
      indices(3) = Int(Nint(z/header%stepz+1.0_RDbl))
    Else
      indices(1) = Int(x/header%stepx+1.0_RDbl)
      indices(2) = Int(y/header%stepy+1.0_RDbl)
      indices(3) = Int(z/header%stepz+1.0_RDbl)
    End If

  End Subroutine mapheader_getIndices

  !--------------------------------------------------
  ! Display the header information
  !--------------------------------------------------
  Subroutine mapheader_display(header, nspc, optunit)
    Type(Header_Info), Intent(in) :: header
    Integer, Intent(in)  :: nspc ! No. of spaces from the left column
    Integer, Optional, Intent(in) :: optunit

    Integer   :: unitno, i
    Character(len=nspc) :: spc

    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    Write(unitno,'(2a)') spc, dashedline
    Write(unitno,'(2a)') spc,"Pmap Header Section:"
    Write(unitno,'(a,2x,a,i8)') spc, 'No. of tabulated pts.: ', header%nfsize
    Write(unitno,'(a,2x,a,3f8.3)') spc, 'Edge Lengths         : ', header%eff
    Write(unitno,'(a,2x,a,3i5)') spc, 'Grid size (x, y, z)  : ', &
        header%nbrx, header%nbry, header%nbrz
    Write(unitno,'(a,2x,a,3f8.3)') spc, 'Grid spacing         : ', & 
        header%stepx, header%stepy, header%stepz
    Write(unitno,'(a,2x,a,f14.3)') spc, 'Potential Cutoff(kJ) : ', &
        header%pcut_off
    Write(unitno,'(a,2x,a,i5)') spc, 'No. of types of atoms: ', header%nzatoms

    Do i=1, header%nzatoms
      Write(unitno,'(a,4x,2a)') spc, 'Atom: ', Trim(header%zatom(i)%atom_name)
      Write(unitno,'(a,4x,2a)') spc, 'Symbol    : ', &
          Trim(header%zatom(i)%symbol)
      IF(Trim(cleanstring(header%mapheader_identifier)) == &
          'Old type header 2003') Then
      ELse
        ! this means premusic or music-2-2 maps
        Write(unitno,'(a,4x,a,f8.4)') spc, 'Charge    : ', &
            header%zatom(i)%charge
      Endif
      Write(unitno,'(a,4x,a,4f14.5)') spc, 'Parameters: ', &
          header%zatom(i)%A12, header%zatom(i)%B6, header%zatom(i)%C2, &
          header%zatom(i)%D0
      IF(header%zatom(i)%rcuthi > 0.) THEN
        Write(unitno,'(a,4x,a,2f14.3)') spc, 'Cutoffs (High, Low)  : ', & 
            header%zatom(i)%rcuthi, header%zatom(i)%rcutlo
      END IF
    End Do
    Write(unitno,*) 'Warning: charges are taken from molecule file and are not displayed here!'

  End Subroutine mapheader_display

  !---------------------------------------------------
  ! Deallocate the allocated memory
  !---------------------------------------------------
  Subroutine mapheader_cleanup(header)
    Type(Header_Info), Intent(inout) :: header
    Integer     :: error

    If (Associated(header%zatom)) Then
      Deallocate(header%zatom, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'header%zatom')   
    End If

  End Subroutine mapheader_cleanup

End Module mapheader




