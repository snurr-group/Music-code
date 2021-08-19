!-------------------------------------------------------------------------
! This module handles the interface of MUSIC with Gaussian98.  All
! Input/Output (IO) calls go through this module.  It handles the
! creation of Gaussian98 input files, the calling of Gaussian98 and the
! translation of the output into something that MUSIC understands.  For
! now, it is setup only to do single-point calculations.  However, it is
! written to be easily generalized to optimizations.
!
! Needed Improvements: 
! 1) 
!-------------------------------------------------------------------------

Module gaussio

  Use defaults, Only: strLen,lstrLen,RDbl,zeroTolerance,xlstrLen
  Use file, Only: file_open, file_append
  Use general, Only: genparams
  Use utils, Only: toint,filesrchstr,split,toupper,toreal,allocerrdisplay, &
      int2str,firstchars,digits,deallocerrdisplay,syscall
  Use vector, Only: VecType,vector_display,Assignment(=),Operator(+), &
      Operator(-),Operator(*),Operator(/)
  Use readgauss2, Only: G98system,readgauss2_getsystem,readgauss2_charges, &
      readgauss2_config2xyz,readgauss2_display,readgauss2_nrg,readgauss2_forces
  Use atom, Only: atom_getsymbol,atom_getcharge,atom_getntypes,atom_getname, &
      atom_iontype
  Use molecules, Only: molecules_getnatoms,molecules_getcharge,molecules_name, &
      molecules_getnsorbs,molecules_gettype,molecules_getnatomtypes

  Implicit None

  Private
  Public :: G98_Input,gaussio_snglpt,gaussio_init,gaussio_display,gaussio_clean

  Character(len=strLen)   :: gaussio_outfile = 'gaussian_run.log'
  Character(len=strLen)   :: gaussio_inputfile = 'gaussian_run.inp'

  !** Contains static information needed for every external call
  Type G98_Input
    Logical                  :: calcq,file_input
    Character(len=lstrLen)   :: head_file,tail_file
  End Type G98_Input

Contains

  !------------------------------------------------------------------------
  ! Initialize the G98 input data type.  The input data is taken either
  ! in the form of a set of parameters or generic header and tail files
  ! that will be used to construct the G98 input file.  If the generic
  ! files are present, they will be used to setup the control file.
  ! File contents:
  !   header file (#1) -- all the G98 commands before the comment, 
  !                       do not include the charge line
  !   tail file (#2) -- Any extra commands to include after the coordinates,
  !                     such as basis set definition, etc.
  ! Requires:  specs -- data type to initialize
  !            calcq -- flags calculation of charges
  !            filenames -- filenames for head (+tail?) of G98 control file
  !------------------------------------------------------------------------
  Subroutine gaussio_init(specs,calcq,filenames)
    Type(G98_Input), Intent(Out)          :: specs
    Logical, Intent(In)                   :: calcq
    Character(len=lstrLen), Dimension(:), Intent(In), Optional  :: filenames

    Integer         :: error,i

    !** Set defaults
    specs%calcq = calcq
    specs%file_input = .False.
    specs%head_file = 'NULL'
    specs%tail_file = 'NULL'

    !** Set the filenames
    If (Present(filenames)) Then
      specs%head_file = filenames(1)
      If (Size(filenames) > 1) specs%tail_file = filenames(2)
      specs%file_input = .True.
    End If

  End Subroutine gaussio_init

  !-----------------------------------------------------------------------
  ! Performs a single point energy calculation using Gaussian98.  
  ! Requires:  specs -- specifications for input
  !            totalq -- total charge on system
  !            natoms -- number of atoms in arrays
  !            coords -- coordinates of each atom
  !            atypes -- atom type of each atom
  !            subsets -- subset identification of each atom
  !            charges -- output charge of each atom (not always filled)
  !            nrg -- resulting energy from Gaussian98
  !            forces -- optional resulting forces on each atom
  !-----------------------------------------------------------------------
  Logical Function gaussio_snglpt(specs,totalq,natoms,coords,atypes, &
      subsets,nrg,charges,forces)
    Type(G98_Input), Intent(InOut)                     :: specs
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
    Type(G98system)                :: g98_system

    !** set defaults
    gaussio_snglpt = .False.
    nrg = 0.0_RDbl

    !** give feedback if -VERBOSE flag is used
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(*,'(1x,a)') 'Making external call to Gaussian98:'
    End If

    !** Create the input file
    Call gaussio_snglptfile(gaussio_inputfile,specs,totalq,coords(1:natoms), &
        atypes(1:natoms),subsets(1:natoms,:))

    !** Set the command line
    Write(command,'(4a)') 'g98 < ',Trim(gaussio_inputfile),' > ', &
        Trim(gaussio_outfile)

    !** Run Gaussian98
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(*,'(1x,2a,i4,6a)') __FILE__," : ",__LINE__,' Running Gaussian98 ',&
          '(input: ',Trim(gaussio_inputfile),') (output: ', &
          Trim(gaussio_outfile),')...'
    End If
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'Running Gaussian98'
    Write(errormsg,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " Error executing Gaussian98 code, please check output file: ", &
          Trim(gaussio_outfile)
    If (.Not. syscall(command,errormsg)) Stop

    !** Analyze the output
    If (Trim(genparams%displaymode) == "VERBOSE") Then
      Write(*,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          ' Analyzing Gaussian98 output (',Trim(gaussio_outfile),')...'
    End If
    Call readgauss2_getsystem(gaussio_outfile,g98_system,0,.False.)
    nrg = readgauss2_nrg(g98_system,1)
    gaussio_snglpt = .True.

    !** Also extract the forces, if desired
    If (Present(forces)) Then
      Call readgauss2_forces(g98_system,1,forces)
    End If

    !** Also extract the charges, if desired
    If (specs%calcq) Then
      Call readgauss2_charges(g98_system,1,charges)
    End If

  End Function gaussio_snglpt

  !---------------------------------------------------------------------------
  ! Creates a Gaussian98 control file for a single point calculation.  
  ! Requires:  filename -- file to create
  !            specs -- specifications for calculation input
  !            totalq -- total charge on system
  !            coords -- coordinates of each atom
  !            atypes -- atom type of each atom
  !            subsets -- subset identification of each atom
  !---------------------------------------------------------------------------
  Subroutine gaussio_snglptfile(filename,specs,totalq,coords,atypes,subsets)
    Character(*), Intent(In)                     :: filename
    Type(G98_Input), Intent(In)                  :: specs
    Real(kind=RDbl), Intent(In)                  :: totalq
    Type(VecType), Dimension(:), Intent(In)      :: coords
    Integer, Dimension(:), Intent(In)            :: atypes
    Integer, Dimension(:,:), Intent(In)          :: subsets

    Integer                        :: unit,charge
    Character(len=xlstrLen)        :: keywords

    !** Open the new control file
    unit = file_open(filename)

    !** Create or read the header from a different file
    If (specs%file_input) Then
      Call file_append(filename,specs%head_file,.True.)
      If (specs%tail_file /= 'NULL') Call file_append(filename,specs%tail_file)
    Else
      Call gaussio_makeheader(filename,specs)
    End If

    !** Write the space, comment and charge lines
    Write(unit,*)
    Write(unit,'(a)') 'Created by MUSIC'
    Write(unit,*)
    charge = NInt(totalq)
    If (Abs(totalq - Real(charge)) > 1.0e-4_RDbl) Then
      Write(0,'(1x,2a,i4,a,f8.3,e14.4)') __FILE__," : ",__LINE__, &
          " Could not round specified charge to whole number ",totalq, &
          Abs(totalq - Real(charge))
      Stop      
    End If
    Write(unit,'(i2,4x,i2)') charge,1

    !** Write coordinates to file
    Call gaussio_wrtcoords(filename,specs,coords,atypes)

    !** Write the tail of the control file if necessary
    If (specs%file_input) Then
      If (specs%tail_file /= 'NULL') Call file_append(filename,specs%tail_file)
    Else
      Call gaussio_maketail(filename,specs)
    End If

    !** Close the new control file
    Close(unit=unit)

  End Subroutine gaussio_snglptfile

  !---------------------------------------------------------------------------
  ! Creates and dumps the header of the Gaussian98 control file to the 
  ! specified file
  ! NOT FINISHED
  ! Requires:  filename -- file to create
  !            specs -- specifications for calculation input
  !---------------------------------------------------------------------------
  Subroutine gaussio_makeheader(filename,specs)
    Character(*), Intent(In)                     :: filename
    Type(G98_Input), Intent(In)                  :: specs

    Integer                        :: unit
    Character(len=xlstrLen)        :: keywords

    !** Open the input file or just get unit
    unit = file_open(filename)

    !** Set default keywords
    keywords = ''

    !** Add keywords if necessary
    If (specs%calcq) Then
      Write(keywords,'(2a)') Trim(keywords),' ?'
    End If

    !** Write keywords to file
    Write(unit,'(a)') Trim(keywords)
    Write(unit,*) 

  End Subroutine gaussio_makeheader

  !---------------------------------------------------------------------------
  ! Creates and dumps the tail of the Gaussian98 control file to the 
  ! specified file.
  ! NOT FINISHED
  ! Requires:  filename -- file to create
  !            specs -- specifications for calculation input
  !---------------------------------------------------------------------------
  Subroutine gaussio_maketail(filename,specs)
    Character(*), Intent(In)         :: filename
    Type(G98_Input), Intent(In)      :: specs

    Integer          :: unit

    !** Open the input file or just get unit
    unit = file_open(filename)

    Write(unit,*) 

  End Subroutine gaussio_maketail

  !-----------------------------------------------------------------------
  ! Dumps the MUSIC coordinate structure into a Gaussian98 control file
  ! Requires:  filename -- G98 output filename
  !            specs -- G98 forcefield input data
  !            coords -- coordinates of each atom
  !            atypes -- atom type of each atom
  !-----------------------------------------------------------------------
  Subroutine gaussio_wrtcoords(filename,specs,coords,atypes)
    Character(*), Intent(In)                     :: filename
    Type(G98_Input), Intent(In)                  :: specs
    Type(VecType), Dimension(:), Intent(In)      :: coords
    Integer, Dimension(:), Intent(In)            :: atypes

    Integer                        :: unit,a,natoms
    character(len=2)               :: symbol
    character(len=lstrLen)         :: string

    unit = file_open(filename,110)

    !** Actually dump the coordinates to file in cartesian format
    natoms = Size(coords)
    Do a = 1,natoms
      !** Get element symbol
      symbol = atom_getsymbol(atypes(a))
    
      !** Write to file
      string = vector_display(coords(a),'f12.8')
      Write(unit,'(a,3x,a,3x,a)') Trim(symbol),Trim(string)
    End Do

    Write(unit,*) 

  End Subroutine gaussio_wrtcoords

  !----------------------------------------------------------------------------
  ! Display the parameters
  ! Requires:  specs -- specifications for input
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine gaussio_display(specs,indent,unitno)
    Type(G98_Input), Intent(In)        :: specs
    Integer, Intent(In)                :: indent
    Integer, Intent(In), Optional      :: unitno

    Integer                     :: unit
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)
    unit = 6
    If (Present(unitno)) unit = unitno

    If (specs%file_input) Then
      Write(unit,'(5a)') blank,'Gaussian98 used with parameters from files: ', &
          Trim(specs%head_file),' and ',Trim(specs%tail_file)
    Else
      Write(unit,'(2a)') blank,'Gaussian98 used with user specified parameters'
    End If

  End Subroutine gaussio_display

  !----------------------------------------------------------------------------
  ! Clean the parameters
  ! Requires:  specs -- specifications for input
  !----------------------------------------------------------------------------
  Subroutine gaussio_clean(specs)
    Type(G98_Input), Intent(InOut)        :: specs

  End Subroutine gaussio_clean

End Module gaussio


