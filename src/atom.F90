!------------------------------------------------------------------------------
! This module contains the basic definition of each atom type.
! These are building blocks from which the molecule types are built.
! 
! Currently the atom type information consists of: 
!   Simulation info:  name, definition file name, symbol
!   Physical info:    equilibrium charge, mass
! 
! It is important to note that this information must be static 
! because it deals with a TYPE of atom, not necessarily unique,
! individual atoms in the simulation.
!------------------------------------------------------------------------------

Module atom

  Use defaults, Only: RDbl, strLen, MAX_ATOMTYPES, lstrLen, d_ctrl_file, &
      d_aa_file, dashedline, one, zero
  Use utils, Only: split, filesrchstr, getpath, stripcmnt, toreal, toint, &
      toupper, isblank
  Use file, Only: file_settag, file_gettype, file_getunit, file_open

  Implicit None
  Save

  Private
  Public :: Atomic_Params, atom_getmass, atom_gettypename, atom_invmass, &
      atom_nbonds, atom_getntypes, atom_getsymbol, atom_getname, &
      atom_checkinit, Atom_List, atoms, atom_init, atom_display, &
      atom_sampleCF, atom_atomicnos, atom_getcharge, atom_invmassMult, &
      atom_isIon, atom_iontype

  Type Atomic_Params
    Character(len=strLen)   :: atom_name
    Character(len=strLen)   :: atom_file
    Character(len=2)        :: symbol
    Real(kind=RDbl)         :: sscharge,szcharge
    Real(kind=RDbl)         :: atom_mass
    Real(kind=RDbl)         :: atom_massi
    Integer                 :: nbonds
    Character(len=5)        :: ionType  ! Values: Core, Shell, None
    Logical                 :: is_ion
  End Type Atomic_Params

  Type Atom_List
    Type(Atomic_Params), Dimension(MAX_ATOMTYPES) :: at
    Integer                                       :: natom_types
  End Type Atom_List

  Type(Atom_List)                                 :: atoms

  Interface getmass
    Module procedure atom_getmass
  End Interface

  Character(strLen), Parameter :: atom_idstring = "Atomic Types"
  Character(strLen), Parameter :: d_atom_file = "Atom File"

  !** provide a translation from atomic number to atomic symbol (1-18 so far)
  Character(len=2), Dimension(18), Parameter :: atom_atomicnos = (/'H ','He',&
      'Li','Be','B ','C ','N ','O ','F ','Ne', &
      'Na','Mg','Al','Si','P ','S ','Cl','Ar'/)

Contains

  !----------------------------------------------------------------------------
  ! Read the control file and find the atom-atom information
  !----------------------------------------------------------------------------
  Subroutine atom_init()
    Integer                               :: unit,lineno,nfields,atom,ios
    Character(len=strLen)                 :: field, path
    Character(len=lstrLen)             :: line,srchstr,filename,filewpath,line2
    Character(len=strLen), Dimension(10)  :: fields
    Logical :: success
    !** Read the atomic types section of the control file
    Call file_gettype(d_ctrl_file,filename,unit)

    Open(unit=unit, file=filename, status='old', IOSTAT=ios)
    If (ios /= 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          '  Could not open file ',trim(filename)
      Stop
    Endif

    srchstr = atom_idstring  !header line
    lineno = filesrchstr(unit, srchstr, line)
    If (lineno == 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          '  Could not find search string ',trim(srchstr)
      Stop
    End If

    !** Read in the number of atom types listed
    Read(unit,*) atoms%natom_types
    If (atoms%natom_types > MAX_ATOMTYPES) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " The maximum no. of atoms exceeded. "
      Write(0,'(a)') "Change MAX_ATOMTYPES in defaults.f90"
      Stop
    End If

    !** Read the blank line separating the number and types
    Read(unit,'(a)') line 
    If (.Not.isblank(stripcmnt(line))) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " The line after the number of atom types should be blank"
      Stop
    End If

    !** Read each type of atom, the name and filename
    Do atom = 1,atoms%natom_types
      Read(unit,'(a)') line
      line2 = stripcmnt(line)
      nfields = split(line2,fields)
      atoms%at(atom)%atom_name = fields(1)
      Read(unit,'(a)') line
      line2 = stripcmnt(line)
      nfields = split(line2,fields)
      atoms%at(atom)%atom_file = fields(1)
      If (atom /= atoms%natom_types) Read(unit,'(a)') line
    End Do

    Close(unit = unit)

    !** Read in the basic atom information from the given filename
    Do atom = 1,atoms%natom_types
      filename = atoms%at(atom)%atom_file
      path = getpath('ATOMSDIR')
      !!      Write(filewpath,'(3a)') Trim(path),'/',Trim(filename)
      !! If getpath has done a good job then you dont need to add the '/'
      Write(filewpath,'(2a)') Trim(path),Trim(filename)
      unit=file_open(filewpath, 110,d_atom_file,success)
      If (.Not.(success)) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Total expected atoms types         : ", atoms%natom_types
        Write(*,*) "Error while initializing atom type : ", atom
        Write(*,*) "Atom name                          : ", &
            atoms%at(atom)%atom_name
        Write(*,*) "Expected File                      : ", filewpath
        Stop
      Endif
!!!      Open(unit=unit, file=filewpath, status='old')

      !** Find the atom name
      srchstr = 'Atom_Name:'
      lineno = filesrchstr(unit, srchstr, line)
      If (lineno == 0) Then
        Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
            '  Could not find search string ',Trim(srchstr), &
            ' in file ',trim(filename)
        Stop
      End If
      !Verify the atom file contains the correct atom
      !      Write(0,*) __FILE__,__LINE__, line
      line2 = stripcmnt(line)
      nfields = split(line2,fields, ":")
      field = fields(2)  ! Compiler Bug: Can't pass an array to adjustl
      If (Trim(Adjustl(field)) /= Trim(atoms%at(atom)%atom_name)) Then 
        Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
            '  Could not find atom ',Trim(atoms%at(atom)%atom_name), &
            ' in file ',Trim(filename)
        Stop
      End If

      !** Get the other parameters
      srchstr = "Atom_Symbol"
      Rewind(unit)
      lineno = filesrchstr(unit, srchstr, line)
      If (lineno == 0) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            '  Could not find search string ',Trim(srchstr)
        Stop
      End If
      line = stripcmnt(line)
      nfields = split(line,fields, ":")
      field = fields(2)   ! Compiler Bug: Can't pass an array to adjustl
      atoms%at(atom)%symbol = Trim(Adjustl(field))
      srchstr = "Atom_SS_Charge"
      Rewind(unit)
      lineno = filesrchstr(unit, srchstr, line)
      If (lineno == 0) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            '  Could not find search string ',Trim(srchstr)
        Write(0,*) "assigning zero"
        atoms%at(atom)%sscharge = zero
      Else
        line = stripcmnt(line)
        nfields = split(line,fields,":")
        atoms%at(atom)%sscharge = toreal(fields(2))
      End If


      srchstr = "Atom_SZ_Charge"
      Rewind(unit)
      lineno = filesrchstr(unit, srchstr, line)
      If (lineno == 0) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            '  Could not find search string ',Trim(srchstr)
        Write(0,*) "assigning zero"
        atoms%at(atom)%szcharge = zero
      Else
        line = stripcmnt(line)
        nfields = split(line,fields,":")
        atoms%at(atom)%szcharge = toreal(fields(2))
      End If


      srchstr = "Atom_Mass"
      Rewind(unit)
      lineno = filesrchstr(unit, srchstr, line)
      If (lineno == 0) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            '  Could not find search string ',Trim(srchstr)
        Stop
      End If
      line = stripcmnt(line)
      nfields = split(line,fields, ":")
      atoms%at(atom)%atom_mass = toreal(fields(2))
      atoms%at(atom)%atom_massi = 1.0_RDbl/atoms%at(atom)%atom_mass

      srchstr = "Atom_Valency"
      Rewind(unit)
      lineno = filesrchstr(unit, srchstr, line)
      If (lineno == 0) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            '  Could not find search string ',Trim(srchstr)
        Stop
      End If
      line = stripcmnt(line)
      nfields = split(line,fields, ":")
      atoms%at(atom)%nbonds = toint(fields(2))

      !** Get info about the atom ion state
      srchstr = "Atom_Ion_Type"
      Rewind(unit)
      lineno = filesrchstr(unit,srchstr,line,.True.)
      If (lineno == 0) Then
        atoms%at(atom)%ionType = "NONE"
        atoms%at(atom)%is_ion = .False.
      Else
        line = stripcmnt(line)
        nfields = split(line,fields,":")
        Select Case (toupper(trim(adjustl(fields(2)))))

        Case ("CORE")
          atoms%at(atom)%ionType = "CORE"
          atoms%at(atom)%is_ion = .True.
        Case ("SHELL")
          atoms%at(atom)%ionType = "SHELL"
          atoms%at(atom)%is_ion = .True.
        Case ("NONE")
          atoms%at(atom)%ionType = "NONE"
          atoms%at(atom)%is_ion = .False.
        Case Default
          Write(0,'(a,i6,3a)') __FILE__,__LINE__, &
              ": Unrecognized ion type ",Trim(fields(2)), &
              ". Allowed options : SHELL, CORE, NONE"
          Stop 
        End Select
      End If

      Close(unit = unit)
    End Do

  End Subroutine atom_init

  !----------------------------------------------------------------------------
  ! Writes a sample of the required input for this routine to unit unitno
  !----------------------------------------------------------------------------
  Subroutine atom_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a)') atom_idstring
    Write(unitno,'(a,t30,a)') 'Integer','# Number of atom types listed'
    Write(unitno,'(a,t30,a)') 'Character', &
        '# atom-atom interaction file filename'
    Write(unitno,'(t30,a)') '# Blank line (required), then for each atom:'
    Write(unitno,'(a,t30,a)') 'Character','# Name of atom'
    Write(unitno,'(a,t30,a)') 'Character','# Atom file filename'
    Write(unitno,'(t30,a)') '# Blank line (required)'
  End Subroutine atom_sampleCF


  !----------------------------------------------------------------------------
  ! Checks to see if the atom structure has been initialized.  
  ! It does this by looking at natoms
  !----------------------------------------------------------------------------
  Logical Function atom_checkinit(atoms)
    Type(Atom_List), Intent(in)   :: atoms

    If (atoms%natom_types == 0) Then
      atom_checkinit = .FALSE.
    Else
      atom_checkinit = .TRUE.
    End If
    Return
  End Function atom_checkinit

  !----------------------------------------------------------------------------
  ! Returns the number of bonds the atom requires
  !----------------------------------------------------------------------------
  Integer Function atom_nbonds(atype)
    Integer, Intent(In) :: atype
    atom_nbonds = atoms%at(atype)%nbonds
  End Function atom_nbonds

  !----------------------------------------------------------------------------
  ! Returns the inverse mass (1/mass) of the atom type specified
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function atom_invmass(atype)
    Integer, Intent(In)           :: atype

    atom_invmass = atoms%at(atype)%atom_massi

  End Function atom_invmass

  !----------------------------------------------------------------------------
  ! Returns an array of inverse mass values corresponding to the atom types
  ! array passed.
  !----------------------------------------------------------------------------
  Subroutine atom_invmassMult(atype,invmass)
    Integer, Intent(IN), Dimension(:)          :: atype
    Real(Kind=RDbl), Dimension(:), Intent(Out) :: invmass
    Integer :: i

    !** Check the size. Die if invmass isn't at least as large at atype.
    If (Size(atype,1) > Size(invmass,1)) Then
      Write(0,'(a,i5,a)') __FILE__,__LINE__," : Error, atype is larger than &
          &invmass array"
      Stop
    End If

    Do i = 1, Size(atype,1)
      invmass(i) = atoms%at(atype(i))%atom_massi
    End Do
  End Subroutine atom_invmassMult

  !----------------------------------------------------------------------------
  ! Gets the atom type i.e. the index number given the SYMBOL
  !----------------------------------------------------------------------------
  Integer Function atom_gettypesym(element)
    Implicit None
    Character(2), Intent(IN)      :: element
    Integer                       :: i

    Do i=1, atoms%natom_types
      If (atoms%at(i)%symbol == toupper(element)) Then
        atom_gettypesym = i
        Return
      End If
    End Do
    atom_gettypesym = 0

  End Function atom_gettypesym

  !----------------------------------------------------------------------------
  ! Gets the atom type i.e. the index number given the NAME
  ! Returns 0 if the name is not in the list
  !----------------------------------------------------------------------------
  Integer Function atom_gettypename(name)
    Character(*), Intent(IN)      :: name
    Integer                       :: i

    Do i=1, atoms%natom_types
      If (toupper(atoms%at(i)%atom_name) &
          == toupper(name)) Then
        atom_gettypename = i
        Return
      End If
    End Do
    atom_gettypename = 0

  End Function atom_gettypename

  !----------------------------------------------------------------------------
  ! Get the symbol of atom type "type"
  !----------------------------------------------------------------------------
  Character(len=2) Function atom_getsymbol(Type)
    Implicit None
    Integer                     :: type

    atom_getsymbol = atoms%at(type)%symbol
    Return
  End Function atom_getsymbol

  !----------------------------------------------------------------------------
  ! Returns the name of the atom given the type
  !----------------------------------------------------------------------------
!!$  Character(len=strLen) Function atom_getname(atomno)
  Function atom_getname(atomno)
    Character(len=strLen) :: atom_getname
    Integer, Intent(in)   :: atomno

    atom_getname = atoms%at(atomno)%atom_name
  End Function atom_getname

  !----------------------------------------------------------------------------
  ! Returns the sscharge of an atom, given the atom type.  Will return 
  ! szcharge, if the optional flag is present and positive.
  ! Requires:  atomno -- atom type
  !            sz_flag -- flags desired return of szcharge
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function atom_getcharge(atomno,sz_flag)
    Integer, Intent(In)             :: atomno
    Logical, Intent(In), Optional   :: sz_flag

    atom_getcharge = atoms%at(atomno)%sscharge

    If (Present(sz_flag)) Then
      If (sz_flag) Then
        atom_getcharge = atoms%at(atomno)%szcharge
      End If
    End If

  End Function atom_getcharge

  !--------------------------------------------------
  ! Returns the number of atom types
  !--------------------------------------------------
  Integer Function atom_getntypes()
    atom_getntypes = atoms%natom_types
  End Function atom_getntypes

  !--------------------------------------------------
  ! Returns the mass of atom type atom
  !--------------------------------------------------
  Real(kind=RDbl) Function atom_getmass(atom)
    Integer, Intent(In) :: atom

    atom_getmass = atoms%at(atom)%atom_mass
  End Function atom_getmass

  !----------------------------------------------------------------------------
  ! Returns the value of the iontype of a given atom type
  !----------------------------------------------------------------------------
  Function atom_iontype(atom)
    Character(len=5)        :: atom_iontype
    Integer, Intent(In)     :: atom

    atom_iontype = atoms%at(atom)%ionType

  End Function atom_iontype

  !----------------------------------------------------------------------------
  ! Returns TRUE if the atom type is a core or shell, FALSE otherwise. If the
  ! optional argument isCore is passed, it will return TRUE if it is a core,
  ! FALSE otherwise.
  !----------------------------------------------------------------------------
  Logical Function atom_isIon(atom,isCore)
    Integer, Intent(In) :: atom
    Logical, Intent(Out), Optional :: isCore

    !** Default value of isCore
    If (Present(isCore)) isCore = .False.
    atom_ision= atoms%at(atom)%is_ion

!!$
!!$    !** Check if the atom is a core or shell
!!$    If (trim(atoms%at(atom)%ionType) == "NONE") Then
!!$      atom_isIon = .False.
!!$    Else
!!$      atom_isIon = .True.
!!$    End If

    !** Check if the atom is a core
    If (atom_isIon) Then 
      If ( ( Trim(atoms%at(atom)%ionType) == "CORE") &
          .And.Present(isCore) ) isCore = .True.
    Endif

  End Function atom_isIon

  !----------------------------------------------------------------------------
  ! Dump the fields of the "atoms" structure
  !----------------------------------------------------------------------------
  Subroutine atom_display(unit)
    Implicit None
    Integer              :: i,unit

    Write(unit,'(a)') dashedline
    Write(unit,'(a)') 'The ATOMS structure:'
    Write(unit,'(2x,a,i5)') 'No. of atom types: ', atoms%natom_types
    Do i=1, atoms%natom_types
      Write(unit,'(4x,2a)')     'Atom Name           : ', &
          Trim(atoms%at(i)%atom_name)
      Write(unit,'(4x,2a)')     'Atom Symbol         : ', & 
          Trim(atoms%at(i)%symbol)
      Write(unit,'(4x,2a)')     'Atom File           : ', & 
          Trim(atoms%at(i)%atom_file)
      Write(unit,'(4x,a, f8.2)')'Atom Mass           : ', & 
          atoms%at(i)%atom_mass
      Write(unit,'(4x,a,2f8.4)')'sscharge, szcharge  : ', & 
          atoms%at(i)%sscharge,atoms%at(i)%szcharge
      Write(unit,'(4x,a,i3)')   'Required # of bonds : ', &
          atoms%at(i)%nbonds
      Write(unit,'(4x,2a)')     'Ion Type            : ', &
          atoms%at(i)%ionType
      If (i /= atoms%natom_types) Write(unit,*)
    End Do
    Write(unit,'(a)') dashedline
  End Subroutine atom_display

End Module atom
