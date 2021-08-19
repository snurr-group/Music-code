!---------------------------------------------------------------------------
! This module handles the reading of GULP output files.  It provides a
! data structure to contain all the information in the file as well as
! routines that can be used without this data structure.  
!
! Currently it is set-up to read the initial configuration and, for
! optimizations, also the final configuration.  
!
! Important Routines:
!     _getsystem -- reads the whole output file and stores information
!                  in the system data structure.  This data structure
!                  is returned and can later be extracted from.
!    _getsnglnrg -- returns the total energy for a single point calculation
!
! Needed Improvements: 
! 1) incorporation of fundcell to generalize cell info
! 2) reading of gradients, hessians etc.
!---------------------------------------------------------------------------

Module readgulp

  Use defaults, Only: strLen,lstrLen,xlstrLen,RDbl,calToJ
  Use file, Only: file_open
  Use vector, Only: VecType,vector_display,Assignment(=),Operator(+), &
      Operator(-),Operator(*),Operator(/)
  Use utils, Only: toint,filesrchstr,split,toupper,toreal,allocerrdisplay
  Use visxyz, Only: XYZ_Entry, visxyz_dump, visxyz_make

  Implicit None
  Save

  Private
  Public :: GULPsystem,GULPconfig,readgulp_getsystem,readgulp_config2xyz, &
      readgulp_display,readgulp_getsnglnrg,readgulp_nrg,readgulp_forces

  Type GULPconfig
    Integer                                      :: natoms
    Real(kind=RDbl)                              :: total_nrg,total_charge
    Type(VecType), Dimension(:), Pointer         :: coord
    Real(kind=RDbl), Dimension(:), Pointer       :: charge
    Character(len=2), Dimension(:), Pointer      :: element
    Character(len=strLen), Dimension(:), Pointer :: label
  End Type GULPconfig

  Type GULPsystem
    Integer               :: nconfigs,natoms
    Character(len=strLen) :: origin,length_units,energy_units
    Real(kind=RDbl), Dimension(3)            :: celldim,angles
    Type(GULPconfig), Dimension(:), Pointer  :: config
  End Type GULPsystem

Contains

  !-----------------------------------------------------------------------
  ! Reads a GULP output file and puts the parameters into the system type
  ! Requires: filename -- GULP output filename
  !           system -- data structure to contain data
  !-----------------------------------------------------------------------
  Subroutine readgulp_getsystem(filename,system)
    Character(*), Intent(In)               :: filename
    Type(GULPsystem), Intent(InOut)        :: system

    Integer                        :: unit,error,i,j,nspacers,nfields,lineno
    Logical                        :: rflag,cluster
    Real(kind=RDbl)                :: number
    Character(len=xlstrLen)        :: line
    Character(len=lstrLen)         :: key_string
    Character(len=lstrLen), Dimension(2)  :: key_strings
    Character(len=strLen), Dimension(20)  :: fields

    system%length_units = 'Angstroms'
    system%energy_units = 'eV'
    system%origin = filename
    system%nconfigs = 0

    !** Open the file 
    unit = file_open(filename,110)
    Rewind(unit=unit)

    !** Determine the type of run from the output header
    lineno = 0
    nspacers = 0
    Do While ((nspacers < 3).And.(lineno < 20))
      lineno = lineno + 1
      Read(unit,'(a)') line
      If (line(2:6) == '*****') Then
        nspacers = nspacers + 1
        Cycle
      End If
      If (nspacers == 2) Then
        nfields = split(line,fields)
        Select Case (ToUpper(Trim(fields(2))))
        Case ('SINGLE')
          system%nconfigs = 1
        Case ('OPTIMISE')
          system%nconfigs = 2
        Case Default
          Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              ' Could not interpret keyword: ',Trim(Trim(fields(2)))
          Stop                    
        End Select
        Exit
      End If
      !NOTE could read other keywords here
    End Do

    !** allocate memory for the configurations
    Allocate(system%config(system%nconfigs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'system%config')

    !** get the number of atoms
    key_string = 'Number of irreducible atoms/shells'
    rflag = readgulp_getnumber(number,key_string,filename,.True.)
    If (.Not. rflag) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Unable to find number of atoms in file'
      Stop          
    End If
    system%natoms = Int(number)

    !** get the cell parameters
    cluster = .False.
    rflag = readgulp_getcell(system%celldim,system%angles,filename,.True.)
    If (.Not. rflag) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Unable to find cell parameters, assuming that it is not periodic'
      cluster = .True.
    End If

    !** get the initial configuration
    If (cluster) Then
      key_strings(1) = 'Cartesian coordinates of cluster'
    Else
      key_strings(1) = 'Fractional coordinates of asymmetric unit'    
    End If
    key_strings(2) = ''
    rflag = readgulp_getconfig(system%config(1),system%natoms, &
        key_strings,filename)
    If (.Not. rflag) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Unable to find initial configuration'
      Stop          
    End If

    !** scale the fractional coordinates  (HACK, accesses vector comp)
    !** should rewrite using fundcell, so it can handle non-ortho cells
    If (.Not. cluster) Then
      Do i = 1,system%natoms
        Do j = 1,3
          system%config(1)%coord(i)%comp(j) = &
              system%config(1)%coord(i)%comp(j)*system%celldim(j)
        End Do
      End Do
    End If

    Call readgulp_config2xyz(system%config(1),'initial.xyz','initial config')

    !** get the final configuration if it's there
    If (system%nconfigs == 2) Then
      key_strings(1) = 'Final'
      key_strings(2) = 'coordinates of atoms'    
      rflag = readgulp_getconfig(system%config(2),system%natoms, &
          key_strings,filename)
      If (.Not. rflag) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Unable to find final configuration'
        Stop          
      End If
  
      !** scale the fractional coordinates  (HACK, accesses vector comp)
      !** should rewrite using fundcell, so it can handle non-ortho cells
      If (.Not. cluster) Then
        Do i = 1,system%natoms
          Do j = 1,3
            system%config(2)%coord(i)%comp(j) = &
                system%config(2)%coord(i)%comp(j)*system%celldim(j)
          End Do
        End Do
      End If
  
      Call readgulp_config2xyz(system%config(2),'final.xyz','final config')
    End If 

    !** get the energies
    Call readgulp_getnrgs(filename,system)

  End Subroutine readgulp_getsystem

  !-----------------------------------------------------------------------
  ! Extract energies from a GULP output file, put in system structure
  ! Requires: filename -- GULP output filename
  !           system -- data structure to contain data
  !-----------------------------------------------------------------------
  Subroutine readgulp_getnrgs(filename,system)
    Character(*), Intent(In)               :: filename
    Type(GULPsystem), Intent(InOut)        :: system

    Integer                        :: unit,i,lineno,nspacers,nfields
    Real(kind=RDbl)                :: number
    Character(len=lstrLen)         :: key_string,line
    Character(len=strLen), Dimension(20)  :: fields

    unit = file_open(filename,110)
    Rewind(unit)

    key_string = 'Components of energy'

    system%energy_units = 'kcal/mol'

    Do i = 1,system%nconfigs
      lineno = filesrchstr(unit,key_string,line,.False.)
      If (lineno == 0) Then
        Write(0,'(2a,i4,a,i5)') __FILE__,": ",__LINE__, &
            ' Unable to find energies in file: ',Trim(filename)
        Stop         
      End If
  
      nspacers = 0
      Do While (nspacers < 4)
        Read(unit,'(a)') line
        If (line(1:5) == '-----') Then
          nspacers = nspacers + 1
          Cycle
        End If
        !NOTE: could get components of energy after nspacers = 1, convert units
        If (nspacers == 3) Then
          nfields = split(line,fields)
          system%config(i)%total_nrg = toreal(fields(5))
        End If
      End Do
    End Do

  End Subroutine readgulp_getnrgs

  !-----------------------------------------------------------------------
  ! Extract single-point energy from a GULP output file.  Returns total
  ! energy in units kcal/mol.
  ! Requires: filename -- GULP output filename
  !-----------------------------------------------------------------------
  Real(kind=RDbl) Function readgulp_getsnglnrg(filename)
    Character(*), Intent(In)               :: filename

    Integer                        :: unit,i,lineno,nspacers,nfields
    Real(kind=RDbl)                :: number
    Character(len=lstrLen)         :: key_string,line
    Character(len=strLen), Dimension(20)  :: fields

    unit = file_open(filename,110)
    Rewind(unit)

    key_string = 'Components of energy'

    lineno = filesrchstr(unit,key_string,line,.False.)
    If (lineno == 0) Then
      Write(0,'(2a,i4,a,i5)') __FILE__,": ",__LINE__, &
          ' Unable to find energies in file: ',Trim(filename)
      Stop         
    End If
    
    nspacers = 0
    Do While (nspacers < 4)
      Read(unit,'(a)') line
      If (line(1:5) == '-----') Then
        nspacers = nspacers + 1
        Cycle
      End If
      !NOTE: could get components of energy after nspacers = 1, convert units
      If (nspacers == 3) Then
        nfields = split(line,fields)
        readgulp_getsnglnrg = toreal(fields(5))
        readgulp_getsnglnrg = readgulp_getsnglnrg/calToJ !** convert to kcal/mol
      End If
    End Do

  End Function readgulp_getsnglnrg

  !-------------------------------------------------------------------------
  ! Reads the cell parameters from the GULP output. 
  ! Returns false if it can't find the number.
  ! Requires: lengths -- vector of cell lengths to get
  !           angles -- vector of cell angles to get
  !           filename -- GULP output file to get configuration from
  !           restart -- optional flag, used to prompt rewind of file
  !-------------------------------------------------------------------------
  Logical Function readgulp_getcell(lengths,angles,filename,restart)
    Real(kind=RDbl), Dimension(3), Intent(Out)  :: lengths,angles
    Character(*), Intent(In)                    :: filename
    Logical, Intent(In)                         :: restart

    Integer                     :: unit,lineno,nfields,idx,nfound,readcount
    Character(len=lstrLen)      :: line,key_string
    Character(len=strLen), Dimension(strLen) :: fields

    readgulp_getcell = .False.

    !** Open the file 
    unit = file_open(filename,110)

    key_string = 'Cell parameters'
    lineno = filesrchstr(unit,key_string,line,restart)
    If (lineno == 0) Return

    nfound = 0
    readcount = 0
    Do While (nfound < 6)
      Read(unit,'(a)') line
      readcount = readcount + 1
      If (Index(line,'=') == 0) Cycle
      nfields = split(line,fields)      

      idx = 0
      Do While (idx < nfields)
        idx = idx + 1
        If (fields(idx) == '=') Cycle
        Select Case (Trim(Toupper(Adjustl(fields(idx)))))
        Case ('A')
          lengths(1) = toreal(fields(idx+2))
          nfound = nfound + 1
        Case ('B')
          lengths(2) = toreal(fields(idx+2))
          nfound = nfound + 1
        Case ('C')
          lengths(3) = toreal(fields(idx+2))
          nfound = nfound + 1
        Case ('ALPHA')
          angles(1) = toreal(fields(idx+2))
          nfound = nfound + 1
        Case ('BETA')
          angles(2) = toreal(fields(idx+2))
          nfound = nfound + 1
        Case ('GAMMA')
          angles(3) = toreal(fields(idx+2))
          nfound = nfound + 1
        Case Default
          !nothing for now
        End Select
      End Do

      If (readcount > 10) Then
        Write(0,'(2a,i4,a,i5)') __FILE__,": ",__LINE__, &
            ' Read too many lines while search for cell parameters ',nfound
        Stop          
      End If

    End Do

    readgulp_getcell = .True.

  End Function readgulp_getcell


  !-------------------------------------------------------------------------
  ! Reads a single number from a GULP output.  This number is assumed
  ! to follow a key string.  For example: [key_string] = [number]
  ! Will rewind the file first, if true optional flag is passed.
  ! Returns false if it can't find the number.
  ! Requires: number -- number to get from file
  !           key_string -- string to find in file
  !           filename -- GULP output file to get configuration from
  !           restart -- optional flag, used to prompt rewind of file
  !-------------------------------------------------------------------------
  Logical Function readgulp_getnumber(number,key_string,filename,restart)
    Real(kind=RDbl), Intent(Out)      :: number
    Character(*), Intent(In)          :: key_string,filename
    Logical, Intent(In)               :: restart

    Integer                     :: unit,lineno,nfields,idx
    Character(len=lstrLen)      :: line
    Character(len=strLen), Dimension(strLen) :: fields

    readgulp_getnumber = .False.

    !** Open the file 
    unit = file_open(filename,110)

    lineno = filesrchstr(unit,key_string,line,restart)
    If (lineno == 0) Return
    If (Index(line,'=') == 0) Return

    nfields = split(line,fields,'=')

    Do idx = 1,nfields
!LC      Write(*,*) fields(idx)
      If (Index(Toupper(fields(idx)),Toupper(key_string)) /= 0) Exit
    End Do

    number = toreal(fields(idx+1))
!LC    Write(*,*) 'getnumber: ',number    

    readgulp_getnumber = .True.    

  End Function readgulp_getnumber

  !-------------------------------------------------------------------------
  ! Reads a configuration from the output file following a given set of 
  ! search strings.  The first search string is processed to determine 
  ! what type of information to expect in the configuration. 
  ! Requires: config -- data type to contain configuration (unsized)
  !           natoms -- number of atoms in configuration
  !           key_strings -- string(s) to find in file, flags configuration
  !           filename -- GULP output file to get configuration from
  !-------------------------------------------------------------------------
  Logical Function readgulp_getconfig(config,natoms,key_strings,filename)
    Type(GULPconfig), Intent(Out)          :: config
    Integer, Intent(In)                    :: natoms
    Character(*), Dimension(:), Intent(In) :: key_strings
    Character(*), Intent(In)               :: filename

    Integer                     :: i,unit,lineno,nsep,error,nstrings
    Logical                     :: frac_flag,charges
    Character(len=lstrLen)      :: line

    config%natoms = natoms
    readgulp_getconfig = .False.

    !** Count the number of strings
    nstrings = 1
    Do nstrings = 1,(Size(key_strings)-1)
      If (Trim(key_strings(nstrings+1)) == '') Exit
    End Do

    !** Open the file 
    unit = file_open(filename,110)

    !** find the key string(s) in the file
    If (nstrings > 1) Then
      lineno = filesrchstr(unit,key_strings,line,.True.)
    Else
      lineno = filesrchstr(unit,key_strings(1),line,.True.)
    End If
    If (lineno == 0) Then
      Return
    End If    

    !** allocate memory for the arrays
    Allocate(config%coord(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'config%coord')
    Allocate(config%element(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'config%element')
    Allocate(config%label(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'config%label')
    Allocate(config%charge(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'config%charge')

    !** interpret the key string
    charges = .False.
    If (Index(Toupper(key_strings(1)),'ASYMMETRIC UNIT') == 0) Then
      charges = .True.  
    End If 

    !** read from the file until the second separator ('---') is encountered
    nsep = 0
    Do While (nsep < 2)
      Read(unit,*) line
      If (line(1:3) == '---') nsep = nsep + 1      
    End Do

    !** actually do the reading
    Call readgulp_getlist(config%coord,config%element,config%label, &
        config%charge,unit)
    readgulp_getconfig = .True.

    !** dump the extra reals if they're not charges
    If (.Not. charges) Deallocate(config%charge)

  End Function readgulp_getconfig

  !-------------------------------------------------------------------------
  ! Reads a list from the output file at the current position.  It does 
  ! this in a generic fashion so it can be used to read initial and
  ! intemediate configurations, as well as gradients.
  ! Requires: vectors -- sized empty array of vectors 
  !           elements -- sized empty array of element symbols
  !           labels -- sized empty array of labels
  !           extras -- sized empty array of real numbers
  !           unit -- unit number, MUST be open and correctly positioned
  ! NOTE: all passed arrays must be the same size
  !-------------------------------------------------------------------------
  Subroutine readgulp_getlist(vectors,elements,labels,extras,unit)
    Type(VecType), Dimension(:), Intent(Out)         :: vectors
    Character(len=2), Dimension(:), Intent(Out)      :: elements
    Character(len=strLen), Dimension(:), Intent(Out) :: labels
    Real(kind=RDbl), Dimension(:), Intent(Out)       :: extras
    Integer, Intent(In)                              :: unit

    Integer                            :: i,ios,nfields,error
    Integer                            :: lineno,asize,direction,index,atomno
    Real(kind=RDbl), Dimension(3)      :: coord
    Character(len=lstrLen)             :: line
    Character(len=strLen), Dimension(strLen) :: fields
    
    asize = Size(vectors)
    If ((Size(elements) /= asize).Or.(Size(labels) /= asize).Or. &
        (Size(extras) /= asize)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed array sizes are incorrect'
      Write(0,'(4i4)') asize,Size(elements),Size(labels),Size(extras)
      Stop          
    End If

    Do lineno = 1,asize
      Read(unit,'(a)',IOSTAT=ios) line
!LC      Write(*,*) line
      If (ios /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Unexpected read error during configuration read'
        Stop          
      End If

      !** check for end of configuration
      If (line(1:3) == '---') Then
        Write(0,'(2a,i4,a,i5)') __FILE__,": ",__LINE__, &
            ' Unexpected end of configuration encountered ',lineno
        Stop
      End If

      nfields = split(line,fields)
      atomno = toint(fields(1))
      If (atomno /= lineno) Then
        Write(0,'(2a,i4,a,2i5)') __FILE__,": ",__LINE__, &
            ' mismatch between atom index and line number ',atomno,lineno
        Stop
      End If

      elements(atomno) = fields(2)(1:2)
      Select Case(Toupper(fields(3)))
      Case ('C')
        labels(atomno) = 'CORE'
      Case ('S')
        labels(atomno) = 'SHELL'
      Case Default
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' Could not understand label: ',Trim(fields(3))
        Stop
      End Select

      !** read the vector, skipping '*'s if they're there
      index = 3
      direction = 0
      Do While (direction < 3)
        index = index + 1
        If (Trim(Adjustl(fields(index))) /= '*') Then
          direction = direction + 1
          coord(direction) = toreal(fields(index))
        End If
      End Do
      vectors(atomno) = coord

      If (Trim(Adjustl(fields(index+1))) == '*') Then
        extras(atomno) = toreal(fields(index+2))
      Else
        extras(atomno) = toreal(fields(index+1))
      End If

    End Do

    Read(unit,*,IOSTAT=ios) line
    If (ios /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Unexpected read error during configuration read'
      Stop          
    End If
    If (line(1:3) /= '---') Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' WARNING: did not read all of the info during list read'
    End If

  End Subroutine readgulp_getlist

  !-------------------------------------------------------------------------
  ! Converts a configuration into an xyz file for visualization
  ! Requires: config -- configuration 
  !           filename -- output xyz filename
  !           comment -- optional comment 
  !-------------------------------------------------------------------------
  Subroutine readgulp_config2xyz(config,filename,comment)
    Type(GULPconfig), Intent(In)      :: config
    Character(*), Intent(In)          :: filename
    Character(*), Optional            :: comment

    Integer                              :: i
    Character(len=strLen)                :: display_format,newsymbol
    Character(len=lstrLen)               :: final_comment
    Type(XYZ_Entry), Dimension(config%natoms)   :: entries

    !** set the comment
    If (.NOT.(Present(comment))) Then
      Write(final_comment,'(a)') 'GULP configuration'
    Else
      final_comment = comment
    End If

    Do i = 1,config%natoms
      !** change the symbol, note that visxyz_make chops after two characters
      newsymbol = config%element(i)
      If (config%label(i) == 'CORE') Then
        Write(newsymbol,'(2a)') Trim(newsymbol),'c'
      Else If (config%label(i) == 'SHELL') Then
        Write(newsymbol,'(2a)') Trim(newsymbol),'s'
      End If

      entries(i) = visxyz_make(config%coord(i),newsymbol)
    End Do

    display_format = 'f12.5'
    Call visxyz_dump(entries,filename,display_format,final_comment)

  End Subroutine readgulp_config2xyz

  !-----------------------------------------------------------------------
  ! Return the total energy for the specified configuration
  ! Requires:  system -- GULP system data structure
  !            configno -- configuration number (1 or maybe 2)
  !            nrg_units -- units of energy
  !-----------------------------------------------------------------------
  Real(kind=RDbl) Function readgulp_nrg(system,configno,nrg_units)
    Type(GULPsystem), Intent(In)                  :: system
    Integer, Intent(In)                           :: configno
    Character(len=strLen), Intent(Out), Optional  :: nrg_units  
    
    readgulp_nrg = system%config(configno)%total_nrg
    If (Present(nrg_units)) Then
      nrg_units = system%energy_units
    End If

  End Function readgulp_nrg

  !-----------------------------------------------------------------------
  ! Return the forces on each for the specified configuration
  ! Requires:  system -- GULP system data structure
  !            configno -- configuration number (1 or maybe 2)
  !            forces -- returned forces vectors on each atom
  !            nrg_units -- units of energy
  !-----------------------------------------------------------------------
  Subroutine readgulp_forces(system,configno,forces,nrg_units)
    Type(GULPsystem), Intent(In)                  :: system
    Integer, Intent(In)                           :: configno
    Type(VecType), Dimension(:), Intent(Out)      :: forces
    Character(len=strLen), Intent(Out), Optional  :: nrg_units  

    Integer         :: a
    
    If (Present(nrg_units)) Then
      nrg_units = system%energy_units
    End If

    Do a = 1,1
      forces(a) = 0.0_RDbl
    End Do

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) ' need to write routines for force extraction'
    Stop

  End Subroutine readgulp_forces

  !-----------------------------------------------------------------------
  ! Display the contents of the system read from the GULP output file
  ! Requires: system -- data structure to contain data
  !           indent -- number of spaces to indent
  !           unit -- unit to write into 
  !-----------------------------------------------------------------------
  Subroutine readgulp_display(system,indent,unit)
    Type(GULPsystem), Intent(In)        :: system
    Integer, Intent(In)                 :: indent,unit

    Character(len=indent)       :: blank
    Character(len=strLen)       :: display

    blank = Repeat(' ',indent)

    Write(unit,'(3a)') blank,'Information from: ',Trim(system%origin)
    Write(unit,'(2x,3a)') blank,'length units: ',Trim(system%length_units)
    Write(unit,'(2x,3a)') blank,'energy units: ',Trim(system%energy_units)
    Write(unit,'(2x,2a,3f10.5)') blank,'initial cell lengths: ', &
        system%celldim
    Write(unit,'(2x,2a,3f10.5)') blank,'initial cell angles: ', &
        system%angles
    Write(unit,'(2x,2a,i4)') blank,'number of atoms: ',system%natoms
    Write(unit,'(2x,2a,i4)') blank,'number of configurations: ',system%nconfigs
    Write(unit,'(2x,2a)') blank,'(configurations not yet displayed)'

  End Subroutine readgulp_display

End Module readgulp


