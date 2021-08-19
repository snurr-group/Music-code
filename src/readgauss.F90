!-----------------------------------------------------------------------
! This module handles the reading and storage of formatted checkpoint
! files from Gaussian98.  These module is structured on the assumed that
! the majority of the G98 checkpoint file will be first read then 
! processed as a whole.  "readgauss_getsystem" currently reads the 
! coordinates, atomic numbers, gradients and the Hessian.  Additional
! routines are then available for processing this collection of 
! information.
!-----------------------------------------------------------------------
Module readgauss

  Use defaults, Only: strLen,lstrLen,RDbl,HTokJmol,bohrToAng
  Use file, Only: file_open
  Use utils, Only: split,filesrchstr,findstr,toreal,toint,allocerrdisplay, &
      deallocerrdisplay
  Use vector, Only: VecType,vector_display,Assignment(=),Operator(+), &
      Operator(-),Operator(*),Operator(/)
  Use matrix, Only: MatrixType,matrix_display,Assignment(=),Operator(*)
  Use hessian, Only: HessianInfo,Assignment(=),hessian_init, &
      hessian_scalehessian,hessian_dispmode,hessian_rotate,hessian_removert, &
      hessian_dispevalues,hessian_eigeninfo,hessian_clean,hessian_display
  Use readstruc, Only: Structure,readstruc_set
  Use visxyz, Only: XYZ_Entry, visxyz_dump, visxyz_make
  Use atom, Only: atom_atomicnos

  Implicit None
  Save

  Private
  Public :: GaussianSystem,readgauss_getsystem,readgauss_idstring, &
      readgauss_vissys,readgauss_dispmode,readgauss_getcoords, &
      readgauss_hessian,readgauss_coords,readgauss_natoms, readgauss_display, &
      readgauss_rotatesys,readgauss_getmodes,readgauss_clean, readgauss_getstruc

  Character(len=strLen), Parameter       :: readgauss_idstring = 'G98CHKPTFILE'

  Type GaussianSystem
    Integer               :: natoms,charge
    Real(kind=RDbl)       :: total_nrg
    Character(len=strLen) :: origin,length_units,energy_units
    Type(HessianInfo)     :: hessian
    Integer, Dimension(:), Pointer        :: atomnumbers
    Type(VecType), Dimension(:), Pointer  :: atomcoords,gradients
  End Type GaussianSystem

Contains

  !-----------------------------------------------------------------------
  ! Reads a Gaussian98 formatted checkpoint file and enters the important
  ! system parameters into the provided data type.  Returns information
  ! with length units of Angstroms and energy units of kJ/mol.
  ! Requires: filename -- file to read from
  !           system -- Gaussian checkpoint file system data structure
  !-----------------------------------------------------------------------
  Subroutine readgauss_getsystem(filename,system)
    Character(*), Intent(In)                   :: filename
    Type(GaussianSystem), Intent(InOut)        :: system

    Integer                        :: error,nentries,dim,atom,i,j,index
    Logical                        :: rflag
    Character(strLen)              :: header_string
    Integer, Dimension(10)         :: intarray  !dummy
    Real(kind=RDbl), Dimension(:), Allocatable   :: realarray
    Real(kind=RDbl), Dimension(:,:), Allocatable :: dummyhessian

    system%length_units = 'Bohr'
    system%energy_units = 'Hartrees'
    system%origin = filename
    header_string = 'Number of atoms'
    system%natoms = readgauss_getint(filename,header_string)
    header_string = 'Charge'
    system%charge = readgauss_getint(filename,header_string)
    header_string = 'Total Energy'
    system%total_nrg = readgauss_getreal(filename,header_string)

    !** get the atomic numbers
    header_string = 'Atomic numbers'
    Allocate(system%atomnumbers(system%natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(realarray(system%natoms*3), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    If (.Not.(readgauss_getsection(filename,header_string,nentries,rflag,&
        system%atomnumbers,realarray)).Or.(rflag)) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Unable to get: ',Trim(header_string)
      Stop  
    End If

    !** get the atomic coordinates
    header_string = 'Current cartesian coordinates'
    Allocate(system%atomcoords(system%natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    If (.Not.(readgauss_getsection(filename,header_string, &
        nentries,rflag,intarray,realarray)).Or.(.Not.rflag)) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Unable to get: ',Trim(header_string)
      Stop  
    End If
    atom = 0
    Do i = 1,nentries,3
      atom = atom + 1
      system%atomcoords(atom) = realarray(i:(i+2))
    End Do

    !** get the cartesian gradients
    header_string = 'Cartesian Gradient'
    Allocate(system%gradients(system%natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    If (.Not.(readgauss_getsection(filename,header_string, &
        nentries,rflag,intarray,realarray)).Or.(.Not.rflag)) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Unable to get: ',Trim(header_string)
      Stop  
    End If
    atom = 0
    Do i = 1,nentries,3
      atom = atom + 1
      system%gradients(atom) = realarray(i:(i+2))
    End Do

    !** Get the force constants, ie the Hessian
    header_string = 'Cartesian Force Constants'

    dim = ((system%natoms*3)**2 - system%natoms*3)/2 + system%natoms*3

    Deallocate(realarray)
    Allocate(realarray(dim), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    If (.Not.(readgauss_getsection(filename,header_string, &
        nentries,rflag,intarray,realarray)).Or.(.Not.rflag)) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Unable to get: ',Trim(header_string)
      Stop  
    End If

    dim = (system%natoms*3)
    Allocate(dummyhessian(dim,dim), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    !** place the force constants in a matrix
    index = 0
    Do i = 1,dim
      Do j = 1,i
        index = index + 1
        dummyhessian(i,j) = realarray(index)
        dummyhessian(j,i) = dummyhessian(i,j)
      End Do
    End Do

    Call hessian_init(system%hessian,dummyhessian,.True.)

    Call readgauss_changeunits(system)

!    Call readgauss_display(system,2,6)

  End Subroutine readgauss_getsystem

  !-----------------------------------------------------------------------
  ! Just get the Hessian from the system data structure
  ! Requires:  system -- Hessian data structure
  !            units -- optinonal returned units
  !-----------------------------------------------------------------------
  Type(HessianInfo) Function readgauss_hessian(system,units)
    Type(GaussianSystem), Intent(In)     :: system
    Character(*), Intent(Out), Optional  :: units

    readgauss_hessian = system%hessian
    If (Present(units)) Then
      Write(units,'(4a)') Trim(system%energy_units),'/', &
          Trim(system%length_units),'^2'
    End If

  End Function readgauss_hessian

  !-----------------------------------------------------------------------
  ! Just get the atomic coordinates and elements from the 
  ! system data structure
  ! Requires: system -- Gaussian checkpoint file system data structure
  !           atomcoords -- array of coordinate vectors
  !           elements -- array of elemental symbols
  !-----------------------------------------------------------------------
  Subroutine readgauss_coords(system,atomcoords,elements)
    Type(GaussianSystem), Intent(In)            :: system
    Type(VecType), Dimension(:), Intent(Out)    :: atomcoords
    Character(len=2), Dimension(:), Intent(Out) :: elements

    Integer            :: atom

    If (Size(atomcoords,1) < system%natoms) Then
      Write(0,'(2a,i4,a,2i4)') __FILE__,": ",__LINE__, &
          ' Supplied atomcoords array too small to contain all atoms', &
          Size(atomcoords,1),system%natoms
      Stop      
    End If
    If (Size(elements,1) < system%natoms) Then
      Write(0,'(2a,i4,a,2i4)') __FILE__,": ",__LINE__, &
          ' Supplied elements array too small to contain all atoms', &
          Size(elements,1),system%natoms
      Stop      
    End If

    Do atom = 1,system%natoms
      atomcoords(atom) = system%atomcoords(atom)
      elements(atom) = atom_atomicnos(system%atomnumbers(atom))
    End Do

  End Subroutine readgauss_coords

  !-----------------------------------------------------------------------
  ! Just put the atomic coordinates and elements into a data structure
  ! Requires:  system -- Gaussian checkpoint file system data structure
  !            struc -- new structure to create
  !-----------------------------------------------------------------------
  Subroutine readgauss_getstruc(system,struc)
    Type(GaussianSystem), Intent(In)            :: system
    Type(Structure), Intent(Out)                :: struc

    Integer            :: atom
    Character(len=2), Dimension(system%natoms)  :: elements

    Do atom = 1,system%natoms
      elements(atom) = atom_atomicnos(system%atomnumbers(atom))
    End Do

    Call readstruc_set(struc,system%atomcoords,elements)

  End Subroutine readgauss_getstruc

  !-----------------------------------------------------------------------
  ! Just get the number of atoms from the system data structure
  ! Requires: system -- Gaussian checkpoint file system data structure
  !-----------------------------------------------------------------------
  Integer Function readgauss_natoms(system)
    Type(GaussianSystem), Intent(In)            :: system

    readgauss_natoms = system%natoms
  End Function readgauss_natoms

  !-----------------------------------------------------------------------
  ! Reads a Gaussian98 formatted checkpoint file and extracts just the 
  ! coordinates and elements from the file.  Returns vectors in units of
  ! Angstroms.
  ! Requires: atomcoords -- array of coordinate vectors
  !           elements -- array of elemental symbols
  !           filename -- checkpoint file filename
  ! NOTE: HAS NOT BEEN TESTED!!
  !-----------------------------------------------------------------------
  Subroutine readgauss_getcoords(atomcoords,elements,filename)
    Type(VecType), Dimension(:), Pointer       :: atomcoords
    Character(len=2), Dimension(:), Pointer    :: elements
    Character(*), Intent(In)                   :: filename

    Integer                        :: error,atom,i
    Integer                        :: nentries,natoms
    Logical                        :: rflag
    Character(strLen)              :: header_string
    Integer, Dimension(:), Allocatable           :: intarray 
    Real(kind=RDbl), Dimension(:), Allocatable   :: realarray

    header_string = 'Number of atoms'
    natoms = readgauss_getint(filename,header_string)

    !** size storage arrays
    Allocate(atomcoords(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(elements(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** size temporary arrays
    Allocate(realarray(natoms*3), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(intarray(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** get the atomic coordinates
    header_string = 'Current cartesian coordinates'
    If (.Not.(readgauss_getsection(filename,header_string, &
        nentries,rflag,intarray,realarray)).Or.(.Not.rflag)) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Unable to get: ',Trim(header_string)
      Stop  
    End If
    atom = 0
    Do i = 1,nentries,3
      atom = atom + 1
      atomcoords(atom) = bohrToAng*realarray(i:(i+2))
    End Do

    !** get the atomic numbers
    header_string = 'Atomic numbers'
    If (.Not.(readgauss_getsection(filename,header_string,nentries,rflag,&
        intarray,realarray)).Or.(rflag)) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Unable to get: ',Trim(header_string)
      Stop  
    End If

    !** convert to element symbols
    Do i = 1,natoms
      elements = atom_atomicnos(intarray(i))
    End Do

    !** deallocate temporary variables
    Deallocate(realarray, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    Deallocate(intarray, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)    

  End Subroutine readgauss_getcoords

  !-------------------------------------------------------------------------
  ! Changes the system units from Hartrees and Bohr to kJ/mol and Angstroms
  !-------------------------------------------------------------------------
  Subroutine readgauss_changeunits(system)
    Type(GaussianSystem), Intent(InOut)    :: system

    Integer                        :: i
    Real(kind=RDbl)                :: factor

    If ((system%length_units /= 'Bohr').Or. &
        (system%energy_units /= 'Hartrees')) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' incoming units must be Hartrees and Bohr'
      Stop         
    End If

    system%length_units = 'Angstroms'
    system%energy_units = 'kJ/mol'

    system%total_nrg = system%total_nrg*HTokJmol

    !** want to convert hessian from Hartree/Bohr^2 to kJ/mol/Ang^2
    factor = HTokJmol / (bohrToAng**2)
    Call hessian_scalehessian(system%hessian,factor)

    !** convert the coordinates from Bohr to Angstrom
    Do i = 1,system%natoms
      system%atomcoords(i) = bohrToAng*system%atomcoords(i)
    End Do

    !** convert the gradients from Hartree/Bohr to kJ/mol/Angstrom
    Do i = 1,system%natoms
      system%gradients(i) = (HTokJmol/bohrToAng)*system%gradients(i)
    End Do

  End Subroutine readgauss_changeunits

  !-------------------------------------------------------------------------
  ! Rotate the system using a supplied rotation matrix.  Operates on the
  ! atomic coordinates, the gradients and the Hessian, but not on the 
  ! information derived from the Hessian.
  ! Requires: system -- Gaussian system data structure
  !           trans_vec -- translation vector
  !           rotmtx -- 3x3 rotation matrix
  !-------------------------------------------------------------------------
  Subroutine readgauss_rotatesys(system,trans_vec,rotmtx)
    Type(GaussianSystem), Intent(InOut)    :: system
    Type(VecType), Intent(In)              :: trans_vec
    Type(MatrixType), Intent(In)           :: rotmtx

    Integer                        :: a,natoms
    Type(VecType)                  :: vector

    natoms = system%natoms

    !** transform the atomic coordinates
    Do a = 1,natoms
      vector = system%atomcoords(a)
      vector = vector + trans_vec
      vector = rotmtx*vector
      system%atomcoords(a) = vector
    End Do

    !** transform the gradients
    Do a = 1,natoms
      system%gradients(a) = rotmtx*system%gradients(a)
    End Do

    !** transform the Hessian
    Call hessian_rotate(system%hessian,rotmtx)

  End Subroutine readgauss_rotatesys

  !-------------------------------------------------------------------------
  ! Processes the Hessian to get the eigenmodes.  Note that the rotational
  ! and translational eigenmodes are first removed.  Produces no output,
  ! just readies the Hessian information for later use.
  ! Requires: system -- Gaussian system data structure
  !-------------------------------------------------------------------------
  Subroutine readgauss_getmodes(system)
    Type(GaussianSystem), Intent(InOut)    :: system

    Integer                        :: dUnit

    dUnit = 6

    Write(*,*) 'Working on HESSIAN'
    Call hessian_eigeninfo(system%hessian)
  
    Write(*,*) 'After Eigenvector determination:'
!    Call hessian_display(system%hessian,2,dUnit)
    Call hessian_dispevalues(system%hessian,2,dUnit)

    Write(*,*) 'Removing rotational and translational components'
    Call hessian_removert(system%hessian,system%atomcoords)

    Write(*,*) 'Working on HESSIAN'
    Call hessian_eigeninfo(system%hessian)

    Call hessian_dispevalues(system%hessian,2,dUnit)

  End Subroutine readgauss_getmodes

  !--------------------------------------------------------------------
  ! Reads a INTEGER number corresponding to some header string 
  !--------------------------------------------------------------------
  Real(kind=RDbl) Function readgauss_getint(filename,header_string)
    Character(strLen), Intent(In)  :: filename,header_string

    Logical                        :: rflag
    Integer                        :: nentries
    Integer, Dimension(1)          :: intarray
    Real(kind=RDbl), Dimension(1)  :: realarray

    If ((readgauss_getsection(filename,header_string,nentries, &
        rflag,intarray,realarray)).And.(.Not.rflag)) Then
       readgauss_getint = intarray(1)
    Else 
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Unable to get integer number for: ',Trim(header_string)
      Stop  
    End If

  End Function readgauss_getint


  !--------------------------------------------------------------------
  ! Reads a REAL number corresponding to some header string 
  !--------------------------------------------------------------------
  Real(kind=RDbl) Function readgauss_getreal(filename,header_string)
    Character(strLen), Intent(In)              :: filename,header_string

    Logical                        :: rflag
    Integer                        :: nentries
    Integer, Dimension(1)          :: intarray
    Real(kind=RDbl), Dimension(1)  :: realarray

    If ((readgauss_getsection(filename,header_string,nentries, &
        rflag,intarray,realarray)).And.(rflag)) Then
       readgauss_getreal = realarray(1)
    Else 
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Unable to get real number for: ',Trim(header_string)
      Stop  
    End If

  End Function readgauss_getreal


  !--------------------------------------------------------------------
  ! Reads a portion of the file specified by the header string into
  ! either a integer or a real one-dimensional array.  Returns "False"
  ! if some error occurs.  For example, if it can't find the string
  ! Requires: filename -- name of file
  !           header_string -- header string
  !           n -- returned number of entries
  !           rflag -- logical flag true if it's a real
  !           intarray -- returned array of integers
  !           realarray -- returned array of reals
  !--------------------------------------------------------------------
  Logical Function readgauss_getsection(filename,header_string,nentries, &
      rflag,intarray,realarray)
    Character(strLen), Intent(In)              :: filename,header_string
    Integer, Intent(Out)                       :: nentries
    Logical, Intent(Out)                       :: rflag
    Integer, Dimension(:), Intent(Out)         :: intarray
    Real(kind=RDbl), Dimension(:), Intent(Out) :: realarray

    Integer                        :: unitno,lineno,index,i,ios,n
    Integer                        :: lo,hi,nfields,dim,rcheck
    Character(len=lstrLen)         :: line
    Character(len=strLen), Dimension(strLen) :: fields

    readgauss_getsection = .False.

    !** Open the file 
    unitno = file_open(filename,110)
    !** find the header string in the file
    lineno = filesrchstr(unitno,header_string,line,.True.)
    If (lineno == 0) Then
      readgauss_getsection = .False.
      Return
    End If

    !** split the line and get the data type
    rflag = .False.
    nfields = split(line,fields)
    index = findstr(fields,'I')
    rcheck = findstr(fields,'R')
!    Write(*,*) 'Line: ',Trim(line)
    If ((index + rcheck) == 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' neither tag I nor tag R present on line:'
      Write(0,'(a)') line
      Stop  
    End If
    If (rcheck /= 0) Then 
      rflag = .True.
      index = rcheck
    End If

    !** get the number of entries and determine if it's an array or not
    i = findstr(fields,'N=')
    If (i == 0) Then
      nentries = 1
      If (rflag) Then
        realarray(1) = toreal(fields(index+1))
!        Write(*,*) 'Real: ',realarray(1)
      Else
        intarray(1) = toint(fields(index+1))
!        Write(*,*) 'Integer: ',intarray(1)
      End If
      readgauss_getsection = .True.
      Return
    Else
      nentries = toint(fields(i+1))
    End If

    !** get the array of numbers
    hi = 0
    If (rflag) Then
      dim = Size(realarray,1)
      If (nentries > dim) Then
        Write(0,'(2a,i4,a,i4,a,i4)') __FILE__,": ",__LINE__, &
            ' Passed real array not large enough ',nentries,' > ',dim
        Stop          
      End If
      n = 5  !** number of reals on the line
      Do While (hi < nentries)
        lo = hi + 1
        hi = Min((lo+n-1),nentries)
        Read(unitno,*,IOSTAT=ios) realarray(lo:hi)
        If (ios /= 0) Then
          Write(0,'(2a,i4,a,i3)') __FILE__,": ",__LINE__, &
              ' Problems reading real array, numbers on line assumption: ',n
          Stop          
        End If
      End Do
    Else
      dim = Size(intarray,1)
      If (nentries > dim) Then
        Write(0,'(2a,i4,a,i4,a,i4)') __FILE__,": ",__LINE__, &
            ' Passed integer array not large enough ',nentries,' > ',dim
        Stop          
      End If
      n = 6  !** number of integer on the line
      Do While (hi < nentries)
        lo = hi + 1
        hi = Min((lo+n-1),nentries)
        Read(unitno,*,IOSTAT=ios) intarray(lo:hi)
        If (ios /= 0) Then
          Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              ' Problems reading integer array, numbers on line assumption: ',n
          Stop          
        End If
      End Do
    End If

    Close(unit=unitno)    

    readgauss_getsection = .True.

  End Function readgauss_getsection

  !-----------------------------------------------------------------------
  ! Dumps an xyz file of the system's atom coordinates
  !-----------------------------------------------------------------------
  Subroutine readgauss_vissys(system,filename,fmt,comment)
    Type(GaussianSystem), Intent(In)    :: system
    Character(*), Intent(In)            :: filename
    Character(*), Optional              :: fmt,comment

    Integer                              :: i
    Character(len=lstrLen)               :: final_comment
    Type(XYZ_Entry), Dimension(system%natoms)   :: entries

    !** set the comment
    If (.NOT.(Present(comment))) Then
      Write(final_comment,'(3a)') 'Gaussian system from: ', &
          Trim(system%origin),'  Created by MUSIC code'
    Else
      final_comment = comment
    End If

    Do i = 1,system%natoms
      entries(i) = visxyz_make(system%atomcoords(i), & 
          atom_atomicnos(system%atomnumbers(i)))
    End Do

    If (Present(fmt)) Then
      Call visxyz_dump(entries,filename,fmt,final_comment)
    Else
      Call visxyz_dump(entries,filename,'NULL',final_comment)
    End If

  End Subroutine readgauss_vissys

  !-----------------------------------------------------------------------
  ! Dump a movie of the structure movement along a specified mode
  ! Requires: system -- Gaussian data structure
  !           mode -- the mode to display
  !           filename -- file to write into
  !-----------------------------------------------------------------------
  Subroutine readgauss_dispmode(system,mode,filename)
    Type(GaussianSystem), Intent(In)   :: system
    Integer, Intent(In)                :: mode
    Character(*), Intent(In)           :: filename

    Integer                                    :: i
    Type(XYZ_Entry), Dimension(system%natoms)  :: entries

    Do i = 1,system%natoms
      entries(i) = visxyz_make(system%atomcoords(i), & 
          atom_atomicnos(system%atomnumbers(i)))
    End Do

    Call hessian_dispmode(system%hessian,entries,mode,filename)

  End Subroutine readgauss_dispmode

  !-----------------------------------------------------------------------
  ! Displays the data in the Gaussian system data structure
  ! Requires: system -- Gaussian data structure
  !           indent -- indentation from left margin
  !           unit -- unit to dump into
  !-----------------------------------------------------------------------
  Subroutine readgauss_display(system,indent,unit)
    Type(GaussianSystem), Intent(In)        :: system
    Integer, Intent(In)                     :: indent,unit

    Integer                     :: i,atom,lo,hi
    Character(len=indent)       :: blank
    Character(len=lstrLen)      :: string

    blank = Repeat(' ',indent)

    Write(unit,'(2a,i2)') blank," number of atoms:  ",system%natoms
    Write(unit,'(2a,i2)') blank," charge on system: ",system%charge
    Write(unit,'(2a,e12.4)') blank," total energy:     ",system%total_nrg

    Call hessian_display(system%hessian,indent,unit)

    If (Associated(system%atomnumbers)) Then
      Write(unit,'(2a)') blank," Atom numbers:"
      Do atom = 1,system%natoms,10
        lo = atom
        hi = Min((atom+9),system%natoms)
        Write(unit,'(2x,a,10i3)') blank,(system%atomnumbers(i),i=lo,hi)
      End Do
    End If

    If (Associated(system%atomcoords)) Then
      Write(unit,'(2a)') blank," Atom coordinates:"
      Do atom = 1,system%natoms
        string = vector_display(system%atomcoords(atom),'f8.3')
        Write(unit,'(2x,a,i3,3x,a)') blank,atom,Trim(string)
      End Do
    End If

    If (Associated(system%atomcoords)) Then
      Write(unit,'(2a)') blank," Atom potential energy gradients:"
      Do atom = 1,system%natoms
        string = vector_display(system%gradients(atom),'e12.4')
        Write(unit,'(2x,a,i3,3x,a)') blank,atom,Trim(string)
      End Do
    End If

  End Subroutine readgauss_display

  !-----------------------------------------------------------------------
  ! Cleans the data in the Gaussian system data structure
  ! Requires: system -- Gaussian data structure
  !-----------------------------------------------------------------------
  Subroutine readgauss_clean(system)
    Type(GaussianSystem), Intent(InOut)        :: system

    Integer               :: error

    If (Associated(system%atomnumbers)) Then
      Deallocate(system%atomnumbers, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'atomnumbers')
    End If

    If (Associated(system%atomcoords)) Then
      Deallocate(system%atomcoords, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'atomcoords')
    End If

    If (Associated(system%gradients)) Then
      Deallocate(system%gradients, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'gradients')
    End If

    Call hessian_clean(system%hessian)

  End Subroutine readgauss_clean

End Module readgauss


