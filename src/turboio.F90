!-------------------------------------------------------------------------
! This module handles the interface of MUSIC with Turbomole.  All
! Input/Output (IO) calls go through this module.  It handles the
! creation of Turbomole input files, the calling of Turbomole and the
! translation of the output into something that MUSIC understands.  For
! now, it is setup only to do single-point calculations.  However, it is
! written to be easily generalized to optimizations.
!
! Needed Improvements: 
! 1) 
!-------------------------------------------------------------------------

Module turboio

  Use defaults, Only: strLen,lstrLen,RDbl,zeroTolerance,xlstrLen,bohrToAng
  Use file, Only: file_open, file_append
  Use general, Only: genparams
  Use utils, Only: toint,filesrchstr,split,toupper,toreal,allocerrdisplay, &
      int2str,firstchars,digits,deallocerrdisplay,syscall
  Use vector, Only: VecType,vector_display,Assignment(=),Operator(+), &
      Operator(-),Operator(*),Operator(/)
  Use readturbo, Only: Turbomole_system,readturbo_getsystem,readturbo_charges, &
      readturbo_display,readturbo_nrg,readturbo_forces
  Use atom, Only: atom_getsymbol,atom_getcharge,atom_getntypes,atom_getname, &
      atom_iontype
  Use molecules, Only: molecules_getnatoms,molecules_getcharge,molecules_name, &
      molecules_getnsorbs,molecules_gettype,molecules_getnatomtypes

  Implicit None

  Private
  Public :: Turbomole_Input,turboio_snglpt,turboio_init,turboio_display, &
      turboio_clean 

  Character(len=strLen)   :: turboio_outfile = 'turbomole_run.log'
  Character(len=strLen)   :: turboio_inputfile = 'control'

  !** Contains static information needed for every external call
  Type Turbomole_Input
    Logical                  :: calcq,file_input,ridft
    Character(len=lstrLen)   :: rundir,coords_file,path
  End Type Turbomole_Input

Contains

  !------------------------------------------------------------------------
  ! Initialize the Turobomole input data type.  The input data is taken 
  ! either in the form of a set of parameters or the name of a control file
  ! Requires:  specs -- data type to initialize
  !            calcq -- flags calculation of charges
  !            filenames -- filenames 
  !------------------------------------------------------------------------
  Subroutine turboio_init(specs,calcq,filenames)
    Type(Turbomole_Input), Intent(Out)                :: specs
    Logical, Intent(In)                               :: calcq    
    Character(len=lstrLen), Dimension(:), Intent(In)  :: filenames

    Integer                 :: error,i,unit
    Character(len=lstrLen)  :: line,filename

    If (calcq) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' Turbomole charge calculation not yet supported'
      Stop  
    End If

    !** Set defaults
    specs%ridft = .False.
    specs%calcq = calcq
    specs%file_input = .False.

    !** Set the run directory
    specs%rundir = filenames(1)
    Write(specs%path,'(2a)') Trim(specs%rundir),'/'
    specs%file_input = .True.

    specs%coords_file = 'coord'
    If (Size(filenames) > 1) Then
      specs%coords_file = filenames(2) 
    End If

    !** Open the control file, make sure it's there and read options
    Write(filename,'(3a)') Trim(specs%rundir),'/',Trim(turboio_inputfile)
    unit = file_open(filename,110)
    If (filesrchstr(unit,'$ridft',line,.True.) /= 0) specs%ridft = .True.
    Close(unit=unit)

  End Subroutine turboio_init

  !-----------------------------------------------------------------------
  ! Performs a single point energy calculation using Turbomole.  
  ! Requires:  specs -- specifications for input
  !            totalq -- total charge on system
  !            natoms -- number of atoms in arrays
  !            coords -- coordinates of each atom
  !            atypes -- atom type of each atom
  !            subsets -- subset identification of each atom
  !            charges -- output charge of each atom (not always filled)
  !            nrg -- resulting energy from Turbomole
  !            forces -- optional resulting forces on each atom
  !-----------------------------------------------------------------------
  Logical Function turboio_snglpt(specs,totalq,natoms,coords,atypes, &
      subsets,nrg,charges,forces)
    Type(Turbomole_Input), Intent(InOut)               :: specs
    Real(kind=RDbl), Intent(In)                        :: totalq
    Integer, Intent(In)                                :: natoms
    Type(VecType), Dimension(:), Intent(In)            :: coords
    Integer, Dimension(:), Intent(In)                  :: atypes
    Integer, Dimension(:,:), Intent(In)                :: subsets
    Real(kind=RDbl), Intent(Out)                       :: nrg  
    Real(kind=RDbl), Dimension(:), Intent(Out)         :: charges
    Type(VecType), Dimension(:), Intent(Out), Optional :: forces

    Integer                        :: unit,error,i,j
    Logical                        :: erflag
    Real(kind=RDbl)                :: number
    Character(len=lstrLen)         :: command
    Character(len=xlstrLen)        :: errormsg
    Type(Turbomole_system)         :: turbo_system

    !** set defaults
    turboio_snglpt = .False.
    nrg = 0.0_RDbl

    !** give feedback if -VERBOSE flag is used
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(*,'(1x,a)') 'Making external call to Turbomole:'
    End If

    !** Create the input file
    Call turboio_snglptfile(turboio_inputfile,specs,totalq,coords(1:natoms), &
        atypes(1:natoms),subsets(1:natoms,:))

    !** Set the command line
    If (specs%ridft) Then
      Write(command,'(4a)') 'cd ',Trim(specs%rundir),' && ridft >& ',Trim(turboio_outfile)
    Else
      Write(command,'(4a)') 'cd ',Trim(specs%rundir),' && dscf >& ',Trim(turboio_outfile)
    End If

    !** Run Turbomole
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(*,'(1x,2a,i4,6a)') __FILE__," : ",__LINE__,' Running Turbomole ',&
          '(input: ',Trim(turboio_inputfile),') (output: ', &
          Trim(turboio_outfile),')...'
    End If
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Running Turbomole'
    Write(errormsg,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Error executing Turbomole code, please check output file: ", &
          Trim(turboio_outfile)

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,*) 'skipping Tmole call for now'
    If (.Not. syscall(command,errormsg)) Stop

    !** Analyze the output
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(*,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' Analyzing Turbomole output (',Trim(turboio_outfile),')...'
    End If
    Call readturbo_getsystem(turboio_outfile,specs%rundir,turbo_system)
    nrg = readturbo_nrg(turbo_system)
    turboio_snglpt = .True.

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'energy: ',nrg

    !** Also extract the forces, if desired
    If (Present(forces)) Then
      Call readturbo_forces(turbo_system,forces)
    End If

    !** Also extract the charges, if desired
    If (specs%calcq) Then
      Call readturbo_charges(turbo_system,charges)
    End If

  End Function turboio_snglpt

  !---------------------------------------------------------------------------
  ! Creates a Turbomole control file for a single point calculation.  
  ! Requires:  filename -- file to create
  !            specs -- specifications for calculation input
  !            totalq -- total charge on system
  !            coords -- coordinates of each atom
  !            atypes -- atom type of each atom
  !            subsets -- subset identification of each atom
  !---------------------------------------------------------------------------
  Subroutine turboio_snglptfile(filename,specs,totalq,coords,atypes,subsets)
    Character(*), Intent(In)                     :: filename
    Type(Turbomole_Input), Intent(In)                :: specs
    Real(kind=RDbl), Intent(In)                  :: totalq
    Type(VecType), Dimension(:), Intent(In)      :: coords
    Integer, Dimension(:), Intent(In)            :: atypes
    Integer, Dimension(:,:), Intent(In)          :: subsets

    Integer                        :: unit,charge
    Character(len=xlstrLen)        :: keywords

    !** Create or read the header from a different file
    If (.Not. specs%file_input) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' Not ready to make custom Turbomole control file'
      Stop  
    End If

    !** Write coordinates to file
    Call turboio_wrtcoords(specs,coords,atypes)

  End Subroutine turboio_snglptfile

  !-----------------------------------------------------------------------
  ! Dumps the MUSIC coordinate structure into a Turbomole coordinate file
  ! Requires:  specs -- Turbomole input data
  !            coords -- coordinates of each atom
  !            atypes -- atom type of each atom
  !-----------------------------------------------------------------------
  Subroutine turboio_wrtcoords(specs,coords,atypes)
    Type(Turbomole_Input), Intent(In)                :: specs
    Type(VecType), Dimension(:), Intent(In)      :: coords
    Integer, Dimension(:), Intent(In)            :: atypes

    Integer                        :: unit,a,natoms
    Type(VecType)                  :: posvec
    character(len=2)               :: symbol
    Character(len=lstrLen)         :: string,filename

    Write(filename,'(3a)') Trim(specs%rundir),'/',Trim(specs%coords_file)
    unit = file_open(filename)

    !** Write first line
    Write(unit,'(a)') '$coord'

    !** Actually dump the coordinates to file in cartesian format
    natoms = Size(coords)
    Do a = 1,natoms
      !** Get element symbol
      symbol = atom_getsymbol(atypes(a))

      !** Convert coordinates to Bohr
      posvec = coords(a)/bohrToAng
    
      !** Write to file
      string = vector_display(posvec,'f15.8')
      Write(unit,'(5x,a,3x,a)') Trim(string),Trim(symbol)
    End Do

    !** Write last lines
    Write(unit,'(a)') '$user-defined bonds'
    Write(unit,'(a)') '$end'
    Write(unit,*) 

    !** Close the file
    close(unit=unit)

  End Subroutine turboio_wrtcoords

  !----------------------------------------------------------------------------
  ! Display the parameters
  ! Requires:  specs -- specifications for input
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine turboio_display(specs,indent,unitno)
    Type(Turbomole_Input), Intent(In)        :: specs
    Integer, Intent(In)                :: indent
    Integer, Intent(In), Optional      :: unitno

    Integer                     :: unit
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)
    unit = 6
    If (Present(unitno)) unit = unitno

    If (specs%file_input) Then
      Write(unit,'(5a)') blank,'Turbomole used in run directory: ', &
          Trim(specs%rundir)
    Else
      Write(unit,'(2a)') blank,'Turbomole used with user specified parameters'
    End If

  End Subroutine turboio_display

  !----------------------------------------------------------------------------
  ! Clean the parameters
  ! Requires:  specs -- specifications for input
  !----------------------------------------------------------------------------
  Subroutine turboio_clean(specs)
    Type(Turbomole_Input), Intent(InOut)        :: specs

  End Subroutine turboio_clean

End Module turboio


