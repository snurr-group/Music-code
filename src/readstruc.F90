!------------------------------------------------------------------------
! This module handles the reads, writes and does some processing of 
! structures from files in common chemical structure formats.
! Currently support extensions are: .car
!
! Needed Improvements
! 1) may need to break this into several modules if it gets too large
!------------------------------------------------------------------------

Module readstruc

  Use defaults, Only: RDbl, strLen, lstrLen, MAX_ATOMS
  Use utils, Only:  filesrchstr, split, toupper, allocErrDisplay, real2str, &
      int2str, deallocErrDisplay, str2seq, largestofseq, getext, toreal, &
      skiplines, firstchars, findint, combine, isreal, tolower
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), &
      Operator(-),vector_display
  Use matrix, Only: MatrixType,Assignment(=),Operator(*),matrix_display
  Use molecules, Only: molecules_getnatoms,molecules_getcoords
  Use visxyz, Only: XYZ_Entry,visxyz_make,visxyz_dump
  Use fundcell, Only: Fundamental_Cell, fundcell_center, fundcell_maptocell, &
      fundcell_simpleinit, fundcell_cellstring, fundcell_shift

  Implicit None
  Save

  Private
  Public :: Structure,readstruc_init,Assignment(=),readstruc_xform, &
      readstruc_visxyz,readstruc_set,readstruc_setfromspc,readstruc_natoms, &
      readstruc_copy,readstruc_display,readstruc_clean,readstruc_subcopy, &
      readstruc_toxyz, readstruc_recenter, readstruc_initcopy, &
      readstruc_dumpcar, readstruc_expand

  Type Structure
    Logical                                 :: fftypes_on
    Integer                                 :: natoms
    Character(len=lstrLen)                  :: origin
    Type(VecType), Dimension(:), Pointer    :: coords
    Character(len=2), Dimension(:), Pointer :: elements
    Character(len=4), Dimension(:), Pointer :: fftype
  End Type Structure

  Interface Assignment(=)
    Module Procedure readstruc_copy
  End Interface 

Contains
  !---------------------------------------------------------------------------
  ! Reads a structure from an arbitrary chemical structure file type and
  ! stores the results in the structure data type.  The filename may include
  ! a specification to just extract some of the atoms in the structure.
  ! Example:  structure.car:10-20,56-98    reads just atoms 10-20,56-98
  ! Requires:  struc -- structure data type to initialize
  !            location -- filename of input file, with extension
  !            fcell -- optional fundamental cell to initialize from file
  !---------------------------------------------------------------------------
  Subroutine readstruc_init(struc,location,fcell)
    Type(Structure), Intent(Out)                   :: struc
    Character(*), Intent(In)                       :: location
    Type(Fundamental_Cell), Intent(Out), Optional  :: fcell

    Integer                  :: nfields
    Character(len=strLen)    :: extension
    Character(len=strLen)    :: filename
    Character(len=lstrLen), Dimension(10) :: fields

    struc%origin = location

    !** Determine if there's an atom subset specification at the end
    nfields = split(location,fields,':')
    filename = fields(1)
    If (nfields > 1) filename = combine(fields(1:2))

    !** Get the extension
    extension = getext(fields(1))

    Select Case(ToUpper(Trim(extension)))
    Case ('CAR')
      If (Present(fcell)) Then
        Call readstruc_readcar(struc,filename,fcell)
      Else
        Call readstruc_readcar(struc,filename)
      End If

    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Could not interpret structure filename extension: ',Trim(extension)
      Stop
    End Select

  End Subroutine readstruc_init

  !---------------------------------------------------------------------------
  ! Initialize and copy from one structure to another.
  ! Requires:  struc -- structure data type to initialize
  !            location -- filename of input file, with extension
  !---------------------------------------------------------------------------
  Subroutine readstruc_initcopy(image,old)
    Type(Structure), Intent(Out)      :: image
    Type(Structure), Intent(In)       :: old
    Call readstruc_size(image,old%natoms,old%fftypes_on)
    Call readstruc_copy(image,old)

  End Subroutine readstruc_initcopy

  !---------------------------------------------------------------------------
  ! Reads a structure specifically from a .car file.
  ! Requires:  struc -- structure data type to initialize
  !            filespc -- filename of input file, with extension + extra
  !            fcell -- optional fundamental cell to initialize from file
  !---------------------------------------------------------------------------
  Subroutine readstruc_readcar(struc,filespec,fcell)
    Type(Structure), Intent(InOut)                :: struc
    Character(*), Intent(In)                      :: filespec
    Type(Fundamental_Cell), Intent(Out), Optional :: fcell

    Integer                        :: unit,n,ios,error
    Integer                        :: nfields,natoms,aline
    Logical                        :: partial
    Character(len=200)             :: line
    Real(kind=RDbl), Dimension(3)  :: posvec
    Character(len=strLen)          :: symbol,filename
    Character(len=lstrLen), Dimension(20) :: fields
    Integer, Dimension(:), Allocatable    :: nums

    !** Determine if there's an atom subset specification at the end
    partial = .False.
    nfields = split(filespec,fields,':')
    If (nfields > 1) Then
      n = largestofseq(fields(2))
      Allocate(nums(n), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      natoms = str2seq(fields(2),nums)
      partial = .True.
    End If
    filename = fields(1)

    !** Open the file
    unit = file_open(filename,110)
    natoms = 0

    !** Get the number of atoms in final structure
    If (partial) Then
      natoms = Size(nums)
    Else 
      !** Read the .car file once through to get the number of atoms
      Call skiplines(unit,5)
      ios = 0
      Do While (ios == 0)
        Read(unit,'(a)',IOSTAT=ios) line
        If (Index(line,'end') /= 0) Exit
        natoms = natoms + 1
      End Do
      If (ios /= 0) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' Error during read of file: ',Trim(filename)
        Stop
      End If
    End If

    !** Read the PBC portion of the file and initialize fcell if present
    If (Present(fcell)) Then
      Rewind(unit)
      Call skiplines(unit,2)  !** skip two lines to skip the 'PBC=' line
      n = filesrchstr(unit,'PBC',line,.False.)
      If (n == 0) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' Could not find PBC specification in file: ',Trim(filename)
        Stop
      End If
      nfields = split(line,fields) 
      Call fundcell_simpleinit(fcell,combine(fields(2:7)))
    End If

    !** Allocate the number of atoms in the structure
    Call readstruc_size(struc,natoms,.True.)

    !** Now read the coordinates
    Rewind(unit=unit)
    Call skiplines(unit,5)
    ios = 0
    aline = 0
    natoms = 0
    Do While (ios == 0)
      Read(unit,'(a)',IOSTAT=ios) line
      If (Index(line,'end') /= 0) Exit
      aline = aline + 1

      If (partial) Then
        If (findint(nums,aline) == 0) Cycle
      End If

      nfields = split(line,fields)
      natoms = natoms + 1
      posvec(1) = toreal(fields(2))
      posvec(2) = toreal(fields(3))
      posvec(3) = toreal(fields(4))
      struc%fftype(natoms) = Trim(fields(7))
      struc%coords(natoms) = posvec
      symbol = firstchars(fields(1))
      If (Len(Trim(symbol)) > 2) Then
        Write(0,'(2a,i4,3a)') __FILE__,": ",__LINE__, &
            ' WARNING: truncating symbol "',Trim(symbol),'"'
      End If
      struc%elements(natoms) = symbol(1:2)
    End Do
    If (ios /= 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Error during read of file: ',Trim(filename)
      Stop
    End If

    !** Deallocate space if necessary
    If (partial) Then
      Deallocate(nums, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    End If

  End Subroutine readstruc_readcar

  !---------------------------------------------------------------------------
  ! Dump a structure into a .car file.  The .car file format is very strict
  ! and I don't understand it.  This is just an attempt to get it right.
  ! Requires:  struc -- structure data type to initialize
  !            filename -- filename to dump into
  !            fcell -- optional fundamental cell to make PBC section from
  !---------------------------------------------------------------------------
  Subroutine readstruc_dumpcar(struc,filename,fcell)
    Type(Structure), Intent(InOut)                :: struc
    Character(*), Intent(In)                      :: filename
    Type(Fundamental_Cell), Intent(In), Optional  :: fcell

    Integer                  :: a,unit
    Real(kind=RDbl)          :: charge
    Character(len=2)         :: element
    Character(len=5)         :: symbol,idstring,fftype
    Character(len=strLen)    :: string1
    Character(len=lstrLen)   :: string

    !** Open the file
    unit = file_open(filename)    

    !** Write a generic header
    Write(unit,'(a)') '!BIOSYM archive 3'   !* no idea what this means
    If (Present(fcell)) Then
      Write(unit,'(a)') 'PBC=ON'
    Else
      Write(unit,'(a)') 'PBC=OFF'
    End If
    Write(unit,'(a)') 'Created by MUSIC'
    Write(unit,'(a)') '!DATE'
    If (Present(fcell)) Then
      string = fundcell_cellstring(fcell)
      Write(unit,'(3a)') 'PBC',Trim(string),' (P1)'
    Else
      Write(unit,'(a)') 'PBC nothing'
    End If

    !** Dump the coordinates to file
    Do a = 1,struc%natoms
      string1 = int2str(a)
      Write(symbol,'(2a)') Trim(ToUpper(struc%elements(a))),Trim(string1)
      Write(idstring,'(a)') 'XXX'
      If (struc%fftypes_on) Then
        Write(fftype,'(a)') Trim(tolower(struc%fftype(a)))
      Else
        Write(fftype,'(a)') Trim(tolower(struc%elements(a)))
      End If
      element = ToUpper(struc%elements(a))
      charge = 0.0_RDbl
      Write(unit,'(a,3f15.9,t52,a,i5,t64,a,t72,a,f7.3)') symbol, &
          struc%coords(a),idstring,a,fftype,element,charge
    End Do

    Write(unit,'(a)') 'end '
    Write(unit,'(a)') 'end '

  End Subroutine readstruc_dumpcar

  !---------------------------------------------------------------------------
  ! Set the structure using coordinates and elements from arrays
  ! Requires:  struc -- structure data type to initialize
  !            atomcoords -- atomic coordinate vectors
  !            elements -- elements for each atom
  !---------------------------------------------------------------------------
  Subroutine readstruc_set(struc,atomcoords,elements)
    Type(Structure), Intent(Out)                 :: struc
    Type(VecType), Dimension(:), Intent(In)      :: atomcoords
    Character(len=2), Dimension(:), Intent(In)   :: elements

    Integer        :: a,natoms

    !** Check size consistency
    natoms = Size(atomcoords)
    If (Size(atomcoords) /= Size(elements)) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Passed coordinate and elements arrays incorrectly sized'
      Stop
    End If

    Call readstruc_size(struc,natoms,.False.)

    Do a = 1,natoms
      struc%coords(a) = atomcoords(a)
      struc%elements(a) = elements(a)
    End Do

  End Subroutine readstruc_set

  !---------------------------------------------------------------------------
  ! Set the structure using coordinates and elements from molecules module
  ! Requires:  struc -- structure data type to initialize
  !            spc -- species number
  !---------------------------------------------------------------------------
  Subroutine readstruc_setfromspc(struc,spc)
    Type(Structure), Intent(Out)     :: struc
    Integer, Intent(In)              :: spc

    Integer        :: natoms

    !** Get number of atoms in species type
    natoms = molecules_getnatoms(spc)

    !** Allocate the structure
    Call readstruc_size(struc,natoms,.False.)

    !** Use molecules to fill it in
    Call molecules_getcoords(spc,struc%coords,struc%elements)

  End Subroutine readstruc_setfromspc

  !---------------------------------------------------------------------------
  ! Size the structure 
  ! Requires:  struc -- structure data type to initialize
  !            natoms -- number of atoms in the structure
  !            fftype_flag -- if True, size the fftype arrays also
  !---------------------------------------------------------------------------
  Subroutine readstruc_size(struc,natoms,fftype_flag)
    Type(Structure), Intent(Out)     :: struc
    Integer, Intent(In)              :: natoms
    Logical, Intent(In)              :: fftype_flag

    Integer        :: error

    struc%fftypes_on = .False.

    struc%natoms = natoms
    Allocate(struc%coords(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(struc%elements(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    If (fftype_flag) Then
      struc%fftypes_on = .True.
      Allocate(struc%fftype(natoms), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End If

  End Subroutine readstruc_size

  !---------------------------------------------------------------------------
  ! Copy one structure into another
  ! Requires:  struc1 -- structure data type to copy into
  !            struc2 -- structure data type to get data from
  !---------------------------------------------------------------------------
  Subroutine readstruc_copy(struc1,struc2)
    Type(Structure), Intent(InOut)    :: struc1
    Type(Structure), Intent(In)       :: struc2

    struc1%origin = struc2%origin
    struc1%coords = struc2%coords
    struc1%elements = struc2%elements
    If (struc2%fftypes_on) Then
      struc1%fftype = struc2%fftype
    End If

  End Subroutine readstruc_copy

  !---------------------------------------------------------------------------
  ! Copy a subset of one structure into another.  Note that atomnums must be
  ! correctly sized.
  ! Requires:  struc -- structure data type to initialize
  !            atomnums -- number of atoms to copy into new structure
  !---------------------------------------------------------------------------
  Type(Structure) Function readstruc_subcopy(struc,atomnums)
    Type(Structure), Intent(In)        :: struc
    Integer, Dimension(:), Intent(In)  :: atomnums

    Integer                    :: a,natoms,indx

    !** size the new structure
    natoms = Size(atomnums)
    Call readstruc_size(readstruc_subcopy,natoms,struc%fftypes_on)

    !** Make a new origin
    Write(readstruc_subcopy%origin,'(2a)') 'Subset of: ',Trim(struc%origin)

    !** Do the selective copying
    indx = 0
    Do a = 1,struc%natoms
      If (findint(atomnums,a) == 0) Cycle
      indx = indx + 1
      readstruc_subcopy%coords(indx) = struc%coords(a)
      readstruc_subcopy%elements(indx) = struc%elements(a)
      If (struc%fftypes_on) Then
        readstruc_subcopy%fftype(indx) = struc%fftype(a)
      End If
    End Do

    If (natoms /= indx) Then
      Write(0,'(2a,i4,a,i4,a,i4)') __FILE__,": ",__LINE__, &
          ' Error: Did not copy all of desired atoms, found only: ', &
          indx,' of ',natoms
      Write(0,'(2x,a)') 'There may be repeats or large numbers in the list'
      Stop
    End If

  End Function readstruc_subcopy

  !---------------------------------------------------------------------------
  ! Return the number of atoms in structure
  ! Requires:  struc -- structure data type to initialize
  !---------------------------------------------------------------------------
  Integer Function readstruc_natoms(struc)
    Type(Structure), Intent(In)     :: struc

    readstruc_natoms = struc%natoms

  End Function readstruc_natoms

  !-------------------------------------------------------------------------
  ! Translate the structure and rotate using a supplied rotation matrix.  
  ! Requires:  struc -- structure data type to initialize
  !            trans_vec -- translation vector
  !            rotmtx -- 3x3 rotation matrix
  !-------------------------------------------------------------------------
  Subroutine readstruc_xform(struc,trans_vec,rotmtx)
    Type(Structure), Intent(InOut)           :: struc
    Type(VecType), Intent(In)                :: trans_vec
    Type(MatrixType), Intent(In), Optional   :: rotmtx

    Integer                        :: a
    Type(VecType)                  :: vector

    !** transform the atomic coordinates
    Do a = 1,struc%natoms
      vector = struc%coords(a)
      vector = vector + trans_vec
      If (Present(rotmtx)) Then
        vector = rotmtx*vector
      End if
      struc%coords(a) = vector
    End Do

  End Subroutine readstruc_xform

  !-------------------------------------------------------------------------
  ! Recenter the structure by shifting all coordinates so that one atom is
  ! is in the center of the passed fundamental cell, then apply periodic 
  ! boundary conditions to make sure all atoms are within the cell.  Returns
  ! the translation vector used.
  ! Requires:  struc -- structure data type to recenter
  !            atm -- atom number to recenter structure around
  !            fcell -- fundamental cell data structure defining cell
  !-------------------------------------------------------------------------
  Subroutine readstruc_recenter(struc,atm,fcell,shift)
    Type(Structure), Intent(InOut)       :: struc
    Integer, Intent(In)                  :: atm
    Type(Fundamental_Cell), Intent(In)   :: fcell
    Type(VecType), Intent(Out), Optional :: shift

    Integer                        :: a
    Type(VecType)                  :: shift_vector,center

    !** Transform the atomic coordinates if atom number is not zero
    If (atm > 0) Then
      center = fundcell_center(fcell)
      shift_vector = center - struc%coords(atm)
      Call readstruc_xform(struc,shift_vector)
    Else
      shift_vector = (/0.0_RDbl,0.0_RDbl,0.0_RDbl/)
    End If

    !** Shift all coordinates to their images within the cell
    Do a = 1,struc%natoms
      struc%coords(a) = fundcell_maptocell(fcell,struc%coords(a))
    End Do

    If (Present(shift)) shift = shift_vector

  End Subroutine readstruc_recenter

  !-------------------------------------------------------------------------
  ! Expand the structure into its neighboring images
  ! Requires:  struc -- structure data type to expand
  !            fcell -- fundamental cell data structure defining cell
  !            images -- number of new images in each (a,b,c) direction
  !            newstruc -- resultant structure
  !            newfcell -- new fundamental cell data structure
  !-------------------------------------------------------------------------
  Subroutine readstruc_expand(struc,fcell,images,newstruc,newfcell)
    Type(Structure), Intent(In)                    :: struc
    Type(Fundamental_Cell), Intent(In)             :: fcell
    Integer, Dimension(3), Intent(In)              :: images
    Type(Structure), Intent(Out)                   :: newstruc
    Type(Fundamental_Cell), Intent(Out), Optional  :: newfcell

    Integer                        :: i,j,k,a,nimages,indx

    If (Present(newfcell)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' not prepared to expand fcell also'
      Stop      
    End If

    !** Calculate the number of total cells
    nimages = (images(1) + 1)*(images(2) + 1)*(images(3) + 1)

    !** Initialize the new structure
    Call readstruc_size(newstruc,struc%natoms*nimages,struc%fftypes_on)    

    !** Copy the basics
    newstruc%fftypes_on = struc%fftypes_on
    Write(newstruc%origin,'(2a)') 'Expanded from: ',Trim(struc%origin)
    newstruc%natoms = struc%natoms*nimages

    !** Create the new images by translation of the base
    indx = 0
    Do i = 0,images(1)
      Do j = 0,images(2)
        Do k = 0,images(3)

          !** Create the coordinates by shifting along lattice vectors
          Do a = 1,struc%natoms
            indx = indx + 1
            newstruc%elements(indx) = struc%elements(a)
            If (struc%fftypes_on) newstruc%fftype(indx) = struc%fftype(a)
            newstruc%coords(indx) = fundcell_shift(fcell, &
                struc%coords(a),(/i,j,k/))
          End Do

        End Do
      End Do
    End Do


  End Subroutine readstruc_expand

  !-----------------------------------------------------------------------
  ! Returns a set of xyz entries for the structure.  These, for example,
  ! can be used in the Hessian module to display eigenmodes.
  ! Requires:  struc -- structure data type to display
  !            entries -- .xyz file entry data type for each atom
  !-----------------------------------------------------------------------
  Subroutine readstruc_toxyz(struc,entries)
    Type(Structure), Intent(In)      :: struc
    Type(XYZ_Entry), Dimension(:)    :: entries

    Integer                 :: i,dim

    !** Make sure the array is large enough
    dim = Size(entries)
    If (dim < struc%natoms) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' passed xyz entries array not large enough'
      Stop
    End If

    Do i = 1,struc%natoms
      entries(i) = visxyz_make(struc%coords(i),struc%elements(i))
    End Do

  End Subroutine readstruc_toxyz

  !-----------------------------------------------------------------------
  ! Dumps an xyz file of the structure
  ! Requires:  struc -- structure data type to display
  !            filename -- file to dump to
  !            fmt -- optional string format for number display
  !            comment -- optional comment for .xyz file
  !-----------------------------------------------------------------------
  Subroutine readstruc_visxyz(struc,filename,fmt,comment)
    Type(Structure), Intent(In)      :: struc
    Character(*), Intent(In)         :: filename
    Character(*), Optional           :: fmt,comment

    Character(len=lstrLen)               :: final_comment
    Type(XYZ_Entry), Dimension(struc%natoms)   :: entries

    !** Set the comment
    If (.NOT.(Present(comment))) Then
      Write(final_comment,'(3a)') 'Structure from "', &
          Trim(struc%origin),'" Created by MUSIC code'
    Else
      final_comment = comment
    End If

    Call readstruc_toxyz(struc,entries)

    If (Present(fmt)) Then
      Call visxyz_dump(entries,filename,fmt,final_comment)
    Else
      Call visxyz_dump(entries,filename,'NULL',final_comment)
    End If

  End Subroutine readstruc_visxyz

  !---------------------------------------------------------------------------
  ! Display the structure
  ! Requires:  struc -- structure data type to display
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine readstruc_display(struc,indent,unit)
    Type(Structure), Intent(In)      :: struc
    Integer, Intent(In)              :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1
    Character(len=lstrLen)      :: string2

    blank = Repeat(' ',indent)

    Write(unit,'(2a,i5)') blank,'Number of atoms in structure: ',&
        struc%natoms
    Do i = 1,struc%natoms
      string2 = vector_display(struc%coords(i),'f12.6')
      string1 = ''
      If (struc%fftypes_on) Then
        Write(string1,'(a)') struc%fftype(i)
      End If
      Write(unit,'(2a,3x,a,i5)') blank,struc%elements(i),Trim(string2), &
          Trim(string1),i
    End Do

  End Subroutine readstruc_display

  !---------------------------------------------------------------------------
  ! Clean the structure
  ! Requires:  struc -- structure data type 
  !---------------------------------------------------------------------------
  Subroutine readstruc_clean(struc)
    Type(Structure), Intent(InOut)   :: struc

    Integer        :: error

    Deallocate(struc%coords, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    Deallocate(struc%elements, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)

    If (struc%fftypes_on) Then
      Deallocate(struc%fftype, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    End If

  End Subroutine readstruc_clean

End Module readstruc
