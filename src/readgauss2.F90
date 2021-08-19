!-----------------------------------------------------------------------
! This module handles the reading of GAUSSIAN98 output files.  It 
! provides a data structure to contain all the information in the file 
! as well as routines for display, unit changes etc.  
!
! Currently only the total energy and the atomic coordinates are read
!
! Needed Improvements: 
! 1) transformation routines for rotating and translating the system
! 2) reading of the 'standard' orientation and last orientation
! 3) reading of other information like gradients and maybe charges
!-----------------------------------------------------------------------

Module readgauss2

  Use defaults, Only: strLen,lstrLen,xlstrLen,RDbl,HTokcalmol
  Use file, Only: file_open
  Use vector, Only: VecType,vector_display,Assignment(=),Operator(+), &
      Operator(-),Operator(*),Operator(/)
  Use utils, Only: toint,filesrchstr,split,toupper,toreal,allocerrdisplay, &
      fileSrchStrAll,real2str,deallocerrdisplay
  Use atom, Only: atom_atomicnos
  Use visxyz, Only: XYZ_Entry, visxyz_dump, visxyz_make, visxyz_entrydisplay

  Implicit None
  Save

  Private
  Public :: G98system,G98config,readgauss2_getsystem,readgauss2_config2xyz, &
      readgauss2_nrg,readgauss2_forces,readgauss2_charges, &
      readgauss2_display,readgauss2_clean

  Type G98config
    Integer                                      :: natoms
    Real(kind=RDbl)                              :: total_nrg,total_charge
    Type(VecType), Dimension(:), Pointer         :: coord
    Real(kind=RDbl), Dimension(:), Pointer       :: charge  !*not used
    Character(len=2), Dimension(:), Pointer      :: element
  End Type G98config

  Type G98system
    Integer                 :: nconfigs,natoms
    Character(len=xlstrLen) :: commandline
    Character(len=strLen)   :: origin,length_units,energy_units
    Type(G98config), Dimension(:), Pointer  :: config
  End Type G98system

Contains

  !-----------------------------------------------------------------------
  ! Reads a G98 output file and puts the parameters into the system type
  ! Requires:  filename -- Gaussian98 output filename
  !            system -- data structure to contain data
  !            config_types -- integer flag to select configuration types
  !            getfinal -- flag necessitating read of final configuration
  !    config_types = 0 -> only final configuration IF 'getfinal' is true
  !    config_types = 1 -> at least initial configuration
  !    config_types > 1 -> get configs seqeuentially until config_types 
  !-----------------------------------------------------------------------
  Subroutine readgauss2_getsystem(filename,system,config_types,getfinal)
    Character(*), Intent(In)              :: filename
    Type(G98system), Intent(InOut)        :: system
    Integer, Intent(In)                   :: config_types
    Logical, Intent(In)                   :: getfinal

    Integer                        :: unit,error,i,j,lineno
    Logical                        :: rflag,optimization
    Character(len=lstrLen)         :: key_string,line

    !** Set the units and defaults
    system%length_units = 'Angstroms'
    system%energy_units = 'Hartrees'
    system%origin = filename
    Nullify(system%config)

    !** Just a check for a reasonable call
    If (config_types < 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' passed config_types must be greater than or equal to 0'
      Stop          
    End If

    !** Open the file 
    unit = file_open(filename,110)

    !** Determine whether this was a single point calculation or an optimization
    optimization = .False.
    system%commandline = readgauss2_getcommandline(filename)
    If (Index(system%commandline,'opt') /= 0) optimization = .True.

    !** Get the number of atoms
    system%natoms = readgauss2_getnatoms(filename,optimization)
    If (system%natoms == 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Unable to get number of atoms in file'
      Stop          
    End If

    !** Read configurations if required
    Call readgauss2_readconfigs(filename,system,config_types,getfinal)

    !** Just read the final energy if it's not an optimization
    If (.Not. optimization) Then

      !** Make sure configurations aren't already allocated
      If (Associated(system%config)) Then
        Write(0,'(2a,i4,3a)') __FILE__,": ",__LINE__, &
            ' configurations have already been allocated, unexpected'
        Stop          
      End If

      !** Allocate memory for a single configuration
      system%nconfigs = 1
      Allocate(system%config(system%nconfigs), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'system%config')
      Nullify(system%config(1)%coord)
      Nullify(system%config(1)%charge)
      Nullify(system%config(1)%element)

      If ((Index(system%commandline,'MNDO') /= 0).Or. &
          (Index(system%commandline,'AM1') /= 0)) Then
        key_string = 'Energy'
      Else
        key_string = 'SCF Done:  E(RHF)'
      End If

      rflag = readgauss2_getnumber(system%config(1)%total_nrg,key_string, &
          filename,.True.)
      If (.Not. rflag) Then
        Write(0,'(2a,i4,3a)') __FILE__,": ",__LINE__, &
            ' Unable to read final SCF energy from "',Trim(filename),'"'
        Stop          
      End If
    End If

    !** Change energy units to kcal/mol
    Call readgauss2_changeunits(system)

    !** Make a movie of all the read configurations for debugging
!LC    Call readgauss2_configs2xyz(system%config,'movie.xyz')

  End Subroutine readgauss2_getsystem

  !-----------------------------------------------------------------------
  ! Reads the configurations from a G98 output file and puts the results
  ! into the system type
  ! Requires:  filename -- Gaussian98 output filename
  !            system -- data structure to contain data
  !            config_types -- integer flag to select configuration types
  !            getfinal -- flag necessitating read of final configuration
  !    config_types = 0 -> only final configuration IF 'getfinal' is true
  !    config_types = 1 -> at least initial configuration
  !    config_types > 1 -> get configs seqeuentially until config_types 
  !-----------------------------------------------------------------------
  Subroutine readgauss2_readconfigs(filename,system,config_types,getfinal)
    Character(*), Intent(In)              :: filename
    Type(G98system), Intent(InOut)        :: system
    Integer, Intent(In)                   :: config_types
    Logical, Intent(In)                   :: getfinal

    Integer                        :: unit,error,i,j,lineno
    Integer                        :: config_index,nconfigs_in_file,nconfigs
    Logical                        :: rflag,endflag
    Real(kind=RDbl)                :: nrg
    Character(len=lstrLen)         :: key_string,line

    !** Just a check for a reasonable call
    If (config_types < 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' passed config_types must be greater than or equal to 0'
      Stop          
    End If

    !** Open the file 
    unit = file_open(filename,110)

    !** If we must, get total number of configurations in file. I'm trying
    !** to keep this routine fast by avoiding reading as much as possible.
    !** The last "configuration" is, I *think*, a transformed structure
    !** for which no SCF calculation is done.  Hence the subtraction of 1
    system%nconfigs = config_types
    If ((getfinal).Or.(config_types > 1)) Then
      nconfigs_in_file = fileSrchStrAll(unit,'Input orientation',.True.)
      nconfigs_in_file = nconfigs_in_file - 1
      If (nconfigs_in_file < system%nconfigs) system%nconfigs = nconfigs_in_file 

      If ((getfinal).And.(nconfigs_in_file > system%nconfigs)) Then
        system%nconfigs = system%nconfigs + 1
      End If
    End If
    
    !** leave now if we're done
    If (system%nconfigs == 0) Return

    !** Allocate memory for the configurations
    Allocate(system%config(system%nconfigs), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'system%config')
    config_index = 0

    !** Set-up for reading
    Rewind(unit=unit)
    nconfigs = 0

    !** get configurations until we've reached number given by config_types
    If (config_types > 0) Then
      Do While (nconfigs < system%nconfigs)
        endflag = .False.
        config_index = config_index + 1

        !** get the configuration
        key_string = 'Input orientation'    
        rflag = readgauss2_getconfig(system%config(config_index),system%natoms, &
            key_string,filename)
        If (.Not. rflag) Then
          endflag = .True.
          system%nconfigs = config_index
          Exit
        End If
!LC        Write(*,*) 'got config ',config_index

        !** get the energy
        key_string = 'SCF Done:  E(RHF)'
        rflag = readgauss2_getnumber(nrg,key_string,filename,.False.)
        If (.Not. rflag) Then
          Write(0,'(2a,i4,a,i4)') __FILE__,": ",__LINE__, &
              ' Could not get energy for configuration ',config_index
          Stop          
        End If
        system%config(config_index)%total_nrg = nrg
!LC        Write(*,*) 'got energy ',nrg

        nconfigs = nconfigs + 1
      End Do

    End If

    If (config_index == system%nconfigs) endflag = .True.

    !** Get the final configuration if desired and it hasn't already been read
    If ((getfinal).And.(.Not. endflag)) Then

      !** advance the file configuration by configuration
      Do While (nconfigs < (nconfigs_in_file - 1))
        Write(*,*) 'skipping orientation'
        key_string = 'Input orientation'
        lineno = fileSrchStr(unit,key_string,line,.False.)
        nconfigs = nconfigs + 1
      End Do

      !** actually get the final configuration
      config_index = config_index + 1
      key_string = 'Input orientation'
      rflag = readgauss2_getconfig(system%config(config_index),system%natoms, &
          key_string,filename)
      If (.Not. rflag) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            ' Unable to get final configuration'
        Stop          
      End If
      Write(*,*) 'got config ',config_index

      !** get the energy
      key_string = 'SCF Done:  E(RHF)'
      rflag = readgauss2_getnumber(nrg,key_string,filename,.False.)
      If (.Not. rflag) Then
        Write(0,'(2a,i4,a,i4)') __FILE__,": ",__LINE__, &
            ' Could not get energy for configuration ',config_index
        Stop          
      End If
      system%config(config_index)%total_nrg = nrg

      nconfigs = nconfigs + 1
    End If

    If (system%nconfigs /= config_index) Then
      Write(0,'(2a,i4,a,i4)') __FILE__,": ",__LINE__, &
          ' nconfigs does not match ',config_index,system%nconfigs
      Stop          
    End If

  End Subroutine readgauss2_readconfigs

  !-------------------------------------------------------------------------
  ! Gets the command line from the beginning of a Gaussian98 output file.
  ! This 'line' is assumed to start with a '#P' and continue through one or
  ! more lines until a '-------' spacer line is found.
  ! Requires:  filename -- G98 output file 
  !-------------------------------------------------------------------------
  Function readgauss2_getcommandline(filename)
    Character(len=xlstrLen)           :: readgauss2_getcommandline
    Character(*), Intent(In)          :: filename

    Integer                    :: unit,lineno,nspacers
    Character(len=lstrLen)     :: line

    !** Open file and look for '#P' line
    unit = file_open(filename,110)
    Rewind(unit=unit)
    lineno = fileSrchStr(unit,'#P',line,.True.)

    If (lineno == 0) Then
      Write(0,'(2a,i4,3a)') __FILE__,": ",__LINE__, &
          ' Unable to get command line from "', &
          Trim(filename),'"'
      Stop          
    End If

    !** Read and store until a spacer line is encountered
    readgauss2_getcommandline = ''
    Do 
      !** See if it's a spacer line
      If (Index(line,'------------') /= 0) Exit

      !** Add to output
      Write(readgauss2_getcommandline,'(2a)') Trim(readgauss2_getcommandline), &
          Trim(line)
      
      !** Read the next line
      Read(unit,'(a)') line
    End Do

  End Function readgauss2_getcommandline

  !-----------------------------------------------------------------------
  ! Return the total energy for the configuration
  ! Requires:  system -- Gaussian98 system data structure
  !            configno -- configuration number (1 or maybe 2)
  !            nrg_units -- units of energy
  !-----------------------------------------------------------------------
  Real(kind=RDbl) Function readgauss2_nrg(system,configno,nrg_units)
    Type(G98system), Intent(In)                   :: system
    Integer, Intent(In)                           :: configno
    Character(len=strLen), Intent(Out), Optional  :: nrg_units  
    
    readgauss2_nrg = system%config(configno)%total_nrg
    If (Present(nrg_units)) Then
      nrg_units = system%energy_units
    End If

  End Function readgauss2_nrg

  !-----------------------------------------------------------------------
  ! Return the forces on each atom for the specified configuration
  ! Requires:  system -- Gaussian98 system data structure
  !            configno -- configuration number (1 or maybe 2)
  !            forces -- returned forces vectors on each atom
  !            nrg_units -- units of energy
  !-----------------------------------------------------------------------
  Subroutine readgauss2_forces(system,configno,forces,nrg_units)
    Type(G98system), Intent(In)                   :: system
    Integer, Intent(In)                           :: configno
    Type(VecType), Dimension(:), Intent(Out)      :: forces
    Character(len=strLen), Intent(Out), Optional  :: nrg_units  
    
    If (Present(nrg_units)) Then
      nrg_units = system%energy_units
    End If

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) ' need to write routines for force extraction'
    Stop

  End Subroutine readgauss2_forces

  !-----------------------------------------------------------------------
  ! Return the charge on each atom for the specified configuration
  ! Requires:  system -- Gaussian98 system data structure
  !            configno -- configuration number (1 or maybe 2)
  !            charges -- returned charges on each atom
  !-----------------------------------------------------------------------
  Subroutine readgauss2_charges(system,configno,charges)
    Type(G98system), Intent(In)                   :: system
    Integer, Intent(In)                           :: configno
    Real(kind=RDbl), Dimension(:), Intent(Out)    :: charges
    
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) ' need to write routines for charge extraction'
    Stop

  End Subroutine readgauss2_charges

  !-------------------------------------------------------------------------
  ! Gets the number of atoms in the system from the first configuration
  ! in the file labeled with "Input orientation".
  ! Requires:  filename -- G98 output file to get natoms from
  !            opt -- flag indicating if it's an optimization output
  !-------------------------------------------------------------------------
  Integer Function readgauss2_getnatoms(filename,opt)
    Character(*), Intent(In)          :: filename
    Logical, Intent(In)               :: opt

    Integer                           :: unit,lineno,nspacers
    Character(len=lstrLen)            :: line,key

    readgauss2_getnatoms = 0

    !** Set the key to look for in the output
    If (opt) Then
      key = 'Input orientation'
    Else
      key = 'Z-Matrix orientation'
    End If

    !** Open file and look for initial configuration
    unit = file_open(filename,110)
    lineno = fileSrchStr(unit,key,line,.True.)

    If (lineno == 0) Then
      Write(0,'(2a,i4,3a)') __FILE__,": ",__LINE__, &
          ' Unable to get natoms from initial configuration in "', &
          Trim(filename),'"'
      Stop          
    End If

    nspacers = 0
    Do While (nspacers < 3)
      Read(unit,'(a)') line

      !** see if it's a spacer line
      If (Index(line,'------------') /= 0) Then
        nspacers = nspacers + 1
        Cycle
      End If

      !** atom spec lines assumed to be within 2nd and 3rd spacer lines
      If (nspacers == 2) readgauss2_getnatoms = readgauss2_getnatoms + 1
    End Do

  End Function readgauss2_getnatoms

  !-------------------------------------------------------------------------
  ! Reads a single number from a G98 output.  This number is assumed
  ! to follow a key string.  For example: [key_string] = [number]
  ! Will rewind the file first, if true optional flag is passed.
  ! Returns false if it can't find the number.
  ! Requires: number -- number to get from file
  !           key_string -- string to find in file
  !           filename -- G98 output file to get configuration from
  !           restart -- optional flag, used to prompt rewind of file
  !-------------------------------------------------------------------------
  Logical Function readgauss2_getnumber(number,key_string,filename,restart)
    Real(kind=RDbl), Intent(Out)      :: number
    Character(*), Intent(In)          :: key_string,filename
    Logical, Intent(In)               :: restart

    Integer                     :: unit,lineno,nfields,idx
    Character(len=lstrLen)      :: line
    Character(len=strLen), Dimension(strLen) :: fields

    readgauss2_getnumber = .False.

    !** Open the file 
    unit = file_open(filename,110)

    !** find the key string and return if there isn't a '=' on the line
    lineno = filesrchstr(unit,key_string,.False.,line,restart)
    If (lineno == 0) Return
    If (Index(line,'=') == 0) Return

    !** Skip to the position in the line where the string is
    nfields = split(line,fields,'=')
    Do idx = 1,nfields
      If (Index(Toupper(fields(idx)),Toupper(key_string)) /= 0) Exit
    End Do

    !** Get the number after the string and the '='
    number = toreal(fields(idx+1))

    readgauss2_getnumber = .True.    

  End Function readgauss2_getnumber

  !-------------------------------------------------------------------------
  ! Reads a configuration from the output file following a given string
  ! Requires:  config -- data type to contain configuration (unsized)
  !            natoms -- number of atoms in configuration
  !            key_string -- string to find in file, flags configuration
  !            filename -- G98 output file to get configuration from
  ! In almost all cases, key_string should be "Input orientation", this
  ! is the configuration that G98 dumps before a calculation.  It also
  ! dumps this so-labeled configuration at the end of a run.
  !-------------------------------------------------------------------------
  Logical Function readgauss2_getconfig(config,natoms,key_string,filename)
    Type(G98config), Intent(Out)     :: config
    Integer, Intent(In)              :: natoms
    Character(*), Intent(In)         :: key_string,filename

    Integer                     :: i,unit,lineno,nsep,error
    Logical                     :: frac_flag,charges
    Character(len=lstrLen)      :: line

    config%natoms = natoms
    readgauss2_getconfig = .False.

    !** Open the file 
    unit = file_open(filename,110)

    !** find the key string in the file
    lineno = filesrchstr(unit,key_string,line,.False.)
    If (lineno == 0) Then
      Return
    End If    

    !** allocate memory for the arrays
    Allocate(config%coord(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'config%coord')
    Allocate(config%element(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'config%element')

    !** read from the file until the second separator ('---') is encountered
    nsep = 0
    Do While (nsep < 2)
      Read(unit,*) line
      If (line(1:3) == '---') nsep = nsep + 1      
    End Do

    !** actually do the reading
    Call readgauss2_getlist(config%coord,config%element,unit)
    readgauss2_getconfig = .True.
!LC    Write(*,*) 'done getting configuration'

  End Function readgauss2_getconfig

  !-------------------------------------------------------------------------
  ! Reads a list from the output file at the current position.  It does 
  ! this in a generic fashion so it can be used to read initial and
  ! intermediate configurations.
  ! Requires: vectors -- sized empty array of vectors 
  !           elements -- sized empty array of element symbols
  !           unit -- unit number, MUST be open and correctly positioned
  ! NOTE: all passed arrays must be the same size
  !-------------------------------------------------------------------------
  Subroutine readgauss2_getlist(vectors,elements,unit)
    Type(VecType), Dimension(:), Intent(Out)         :: vectors
    Character(len=2), Dimension(:), Intent(Out)      :: elements
    Integer, Intent(In)                              :: unit

    Integer                            :: i,ios,nfields,error
    Integer                            :: lineno,asize,direction,index,atomno
    Integer                            :: atomic_number
    Real(kind=RDbl), Dimension(3)      :: coord
    Character(len=lstrLen)             :: line
    Character(len=strLen), Dimension(strLen) :: fields
    
    asize = Size(vectors)
    If (Size(elements) /= asize) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed array sizes are incorrect'
      Write(0,'(2i4)') asize,Size(elements)
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

      atomic_number = toint(fields(2))
      elements(atomno) = atom_atomicnos(atomic_number)

      !** read the vector, could use this to skip junk 
      index = 3
      direction = 0
      Do While (direction < 3)
        index = index + 1
        direction = direction + 1
        coord(direction) = toreal(fields(index))
      End Do
      vectors(atomno) = coord

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

  End Subroutine readgauss2_getlist

  !-------------------------------------------------------------------------
  ! Converts a configuration into an xyz file for visualization
  ! Requires: config -- configuration 
  !           filename -- output xyz filename
  !           comment -- optional comment 
  !-------------------------------------------------------------------------
  Subroutine readgauss2_config2xyz(config,filename,comment)
    Type(G98config), Intent(In)       :: config
    Character(*), Intent(In)          :: filename
    Character(*), Optional            :: comment

    Integer                              :: i
    Character(len=strLen)                :: display_format
    Character(len=xlstrLen)              :: final_comment
    Type(XYZ_Entry), Dimension(config%natoms)   :: entries

    !** set the comment
    If (.NOT.(Present(comment))) Then
      Write(final_comment,'(a)') 'Gaussian98 configuration, '
      Write(final_comment,'(3a)') Trim(final_comment),' Energy: ', &
          Trim(real2str(config%total_nrg,12))
    Else
      final_comment = comment
    End If

    Do i = 1,config%natoms
      entries(i) = visxyz_make(config%coord(i),config%element(i))
    End Do

    display_format = 'f12.5'
    Call visxyz_dump(entries,filename,display_format,final_comment)

  End Subroutine readgauss2_config2xyz

  !-------------------------------------------------------------------------
  ! Converts a set of configurations into an xyz movie file for visualization
  ! Requires: configs -- configurations
  !           filename -- output xyz filename
  !           comment -- optional comment 
  !-------------------------------------------------------------------------
  Subroutine readgauss2_configs2xyz(configs,filename,comment)
    Type(G98config), Dimension(:), Intent(In) :: configs
    Character(*), Intent(In)                  :: filename
    Character(*), Optional                    :: comment

    Integer                              :: i,unit,n
    Character(len=strLen)                :: display_format
    Character(len=xlstrLen)              :: final_comment,base_comment
    Type(XYZ_Entry), Dimension(configs(1)%natoms)   :: entries

    display_format = 'f12.5'

    !** set the comment
    If (.NOT.(Present(comment))) Then
      Write(base_comment,'(a)') 'Gaussian98 configuration'
    Else
      base_comment = comment
    End If

    !** open the new file
    unit = file_open(filename)

    Do n = 1,Size(configs)
      Write(final_comment,'(3a)') Trim(base_comment),' Energy: ', &
          Trim(real2str(configs(n)%total_nrg,12))

      Write(unit,'(i10)') Size(entries)
      Write(unit,'(a)') Trim(final_comment)

      Do i = 1,configs(n)%natoms
        entries(i) = visxyz_make(configs(n)%coord(i),configs(n)%element(i))
        Write(unit,'(a,i5)') Trim(visxyz_entrydisplay(entries(i), &
            display_format)),i
      End Do

    End Do

  End Subroutine readgauss2_configs2xyz

  !-------------------------------------------------------------------------
  ! Changes the system units from Hartrees to kcal/mol
  ! Requires: system -- data structure to contain data
  !-------------------------------------------------------------------------
  Subroutine readgauss2_changeunits(system)
    Type(G98system), Intent(InOut)    :: system

    Integer                        :: n

    If (system%energy_units /= 'Hartrees') Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' incoming units must be Hartrees'
      Stop         
    End If

    system%energy_units = 'kcal/mol'

    Do n = 1,system%nconfigs
      system%config(n)%total_nrg = system%config(n)%total_nrg*HTokcalmol
    End Do

  End Subroutine readgauss2_changeunits

  !-----------------------------------------------------------------------
  ! Display the contents of the system read from the G98 output file
  ! Requires: system -- data structure to contain data
  !           indent -- number of spaces to indent
  !           unit -- unit to write into 
  !-----------------------------------------------------------------------
  Subroutine readgauss2_display(system,indent,unit)
    Type(G98system), Intent(In)        :: system
    Integer, Intent(In)                :: indent,unit

    Character(len=indent)       :: blank
    Character(len=strLen)       :: display

    blank = Repeat(' ',indent)

    Write(unit,'(3a)') blank,'Information from: ',Trim(system%origin)
    Write(unit,'(2x,3a)') blank,'length units: ',Trim(system%length_units)
    Write(unit,'(2x,3a)') blank,'energy units: ',Trim(system%energy_units)

    Write(unit,'(2x,2a,i4)') blank,'number of atoms: ',system%natoms
    Write(unit,'(2x,2a,i4)') blank,'number of configurations: ',system%nconfigs
    Write(unit,'(2x,2a)') blank,'(configurations not yet displayed)'

  End Subroutine readgauss2_display

  !-----------------------------------------------------------------------
  ! Clean the G98 system data structure
  ! Requires: system -- data structure containing data
  !-----------------------------------------------------------------------
  Subroutine readgauss2_clean(system)
    Type(G98system), Intent(InOut)        :: system
    
    Integer              :: n,error

    Do n = 1,system%nconfigs
      Deallocate(system%config(n)%coord, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'coord')
      Deallocate(system%config(n)%element, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'element')
    End Do

    Deallocate(system%config, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'config')

  End Subroutine readgauss2_clean

End Module readgauss2


