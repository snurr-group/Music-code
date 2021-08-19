!-----------------------------------------------------------------------
! This module handles the reading of Turbomole output files.  It 
! provides a data structure to contain all the information in the file 
! as well as routines for display, unit changes etc.  
!
! Needed Improvements: 
! 1) add support for charge analysis
!-----------------------------------------------------------------------

Module readturbo

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
  Public :: Turbomole_system,Turbomole_config,readturbo_getsystem, &
      readturbo_nrg,readturbo_forces,readturbo_charges, &
      readturbo_display,readturbo_clean

  Character(len=strLen)   :: readturbo_problemfile = 'dscf_problem'
  Character(len=strLen)   :: readturbo_controlfile = 'control'
  Character(len=strLen)   :: readturbo_nrgfile = 'energy'
  Character(len=strLen)   :: readturbo_gradfile = 'gradient'

  Type Turbomole_config
    Integer                                      :: natoms
    Real(kind=RDbl)                              :: total_nrg,total_charge
    Type(VecType), Dimension(:), Pointer         :: coord,force
    Real(kind=RDbl), Dimension(:), Pointer       :: charge  !*not used
    Character(len=2), Dimension(:), Pointer      :: element
  End Type Turbomole_config

  Type Turbomole_system
    Integer                 :: natoms
    Character(len=strLen)   :: length_units,energy_units
    Type(Turbomole_config)  :: config
  End Type Turbomole_system

Contains

  !-----------------------------------------------------------------------
  ! Reads a Turbomole output file and puts the parameters into the system 
  ! type data structure.  This is essentially the _init routine.
  ! Requires:  outfile -- Turbomole output filename
  !            rundir -- Turbomole run directory
  !            system -- data structure to contain data
  !-----------------------------------------------------------------------
  Subroutine readturbo_getsystem(outfile,rundir,system)
    Character(*), Intent(In)              :: outfile,rundir
    Type(Turbomole_system), Intent(InOut) :: system

    Integer                        :: unit,ios,i,j,lineno,nfields
    Logical                        :: rflag,endflag
    Real(kind=RDbl)                :: number,nrg
    Character(len=lstrLen)         :: key_string,filename
    Character(len=xlstrLen)        :: line,lastline
    Character(len=strLen), Dimension(strLen) :: fields

    !** Set defaults and do nullification
    system%length_units = 'Bohr'
    system%energy_units = 'Hartrees'
    system%config%natoms = 0
    system%config%total_charge = 0
    Nullify(system%config%coord)
    Nullify(system%config%force)
    Nullify(system%config%charge)
    Nullify(system%config%element)

    !** Make sure there isn't a problem file
    Write(filename,'(3a)') Trim(rundir),'/',Trim(readturbo_problemfile)
    Open(file=filename, unit=999, form="FORMATTED",  &
        status="OLD", IOSTAT=ios)
    If (ios == 0) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' Turbomole problem file exists, run stopped'
      Close(unit=999)
      Stop  
    End If

    !** Make sure there isn't an 'actual step' line in the control file
    Write(filename,'(3a)') Trim(rundir),'/',Trim(readturbo_controlfile)
    unit = file_open(filename,110)
    Rewind(unit)
    key_string = 'actual step'
    lineno = filesrchstr(unit,key_string,line)
    If (lineno /= 0) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' "actual step" detected in control file, something is wrong'
      Write(0,'(a)') Trim(line)
      Stop  
    End If
    Close(unit=unit)

    !** Open the energy file and read the last energy entry line
    Write(filename,'(3a)') Trim(rundir),'/',Trim(readturbo_nrgfile)
    unit = file_open(filename,110)
    Do 
      Read(unit,'(a)',IOSTAT=ios) line
      If (ios /= 0) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            ' energy file ended without an $end, please check'
        Write(0,'(a)') Trim(line)
        Stop  
      End If
      If (Index(line,'$end') /= 0) Exit
      lastline = line
    End Do

    !** Pull the energy out of the last line
    nfields = split(lastline,fields)
    system%config%total_nrg = toreal(fields(2))

    Call readturbo_changeunits(system)

  End Subroutine readturbo_getsystem

  !-----------------------------------------------------------------------
  ! Return the total energy for the configuration
  ! Requires:  system -- Turbomole system data structure
  !            nrg_units -- units of energy
  !-----------------------------------------------------------------------
  Real(kind=RDbl) Function readturbo_nrg(system,nrg_units)
    Type(Turbomole_system), Intent(In)            :: system
    Character(len=strLen), Intent(Out), Optional  :: nrg_units  
    
    readturbo_nrg = system%config%total_nrg
    If (Present(nrg_units)) Then
      nrg_units = system%energy_units
    End If

  End Function readturbo_nrg

  !-----------------------------------------------------------------------
  ! Return the forces on each atom for the specified configuration
  ! Requires:  system -- Turbomole system data structure
  !            forces -- returned forces vectors on each atom
  !            nrg_units -- units of energy
  !-----------------------------------------------------------------------
  Subroutine readturbo_forces(system,forces,nrg_units)
    Type(Turbomole_system), Intent(In)            :: system
    Type(VecType), Dimension(:), Intent(Out)      :: forces
    Character(len=strLen), Intent(Out), Optional  :: nrg_units  
    
    If (Present(nrg_units)) Then
      nrg_units = system%energy_units
    End If

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) ' need to write routines for force extraction'
    Stop

  End Subroutine readturbo_forces

  !-----------------------------------------------------------------------
  ! Return the charge on each atom for the specified configuration
  ! Requires:  system -- Turbomole system data structure
  !            charges -- returned forces vectors on each atom
  !-----------------------------------------------------------------------
  Subroutine readturbo_charges(system,charges)
    Type(Turbomole_system), Intent(In)            :: system
    Real(kind=RDbl), Dimension(:), Intent(Out)    :: charges
    
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) ' need to write routines for charge extraction'
    Stop

  End Subroutine readturbo_charges

  !---------------------------------------------------------------------------
  ! Changes the system units from Hartrees and Bohr to kcal/mol and Angstroms
  ! Requires:  system -- data structure to contain data
  !---------------------------------------------------------------------------
  Subroutine readturbo_changeunits(system)
    Type(Turbomole_system), Intent(InOut)            :: system

    Integer                        :: n

    If (system%energy_units /= 'Hartrees') Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' incoming units must be Hartrees'
      Stop         
    End If

    system%energy_units = 'kcal/mol'

    system%config%total_nrg = system%config%total_nrg*HTokcalmol

  End Subroutine readturbo_changeunits

  !-----------------------------------------------------------------------
  ! Display the contents of the system read from the Turbomol output file
  ! Requires: system -- data structure to contain data
  !           indent -- number of spaces to indent
  !           unit -- unit to write into 
  !-----------------------------------------------------------------------
  Subroutine readturbo_display(system,indent,unit)
    Type(Turbomole_system), Intent(In)    :: system
    Integer, Intent(In)                   :: indent,unit

    Character(len=indent)       :: blank
    Character(len=strLen)       :: display

    blank = Repeat(' ',indent)

    Write(unit,'(2x,3a)') blank,'length units: ',Trim(system%length_units)
    Write(unit,'(2x,3a)') blank,'energy units: ',Trim(system%energy_units)

    Write(unit,'(2x,2a,i4)') blank,'number of atoms: ',system%natoms
    Write(unit,'(2x,2a)') blank,'(configurations not yet displayed)'

  End Subroutine readturbo_display

  !-----------------------------------------------------------------------
  ! Clean the Turbomole system data structure
  ! Requires: system -- data structure containing data
  !-----------------------------------------------------------------------
  Subroutine readturbo_clean(system)
    Type(Turbomole_system), Intent(InOut)    :: system
    
    Integer              :: n,error

    If (Associated(system%config%coord)) Then
      Deallocate(system%config%coord, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'coord')
    End If

    If (Associated(system%config%force)) Then
      Deallocate(system%config%force, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'force')
    End If

    If (Associated(system%config%element)) Then
      Deallocate(system%config%element, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'element')
    End If

    If (Associated(system%config%charge)) Then
      Deallocate(system%config%charge, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'charge')
    End If

  End Subroutine readturbo_clean

End Module readturbo


