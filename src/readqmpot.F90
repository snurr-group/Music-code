!-----------------------------------------------------------------------
! This module handles the reading and storage of information from Marek
! Sierka's QM-Pot program.
!
! QM-Pot Hessians are stored in the file as the lower triangle in
! units of kcal/Ang
!
! Needed Improvements:
! 1) only deals with Hessians at the moment, should be expanded and
!    have a more general data type and init routine added.
!-----------------------------------------------------------------------
Module readqmpot

  Use defaults, Only: strLen,lstrLen,RDbl
  Use utils, Only: split,filesrchstr,findstr,toreal,toint,allocerrdisplay, &
      deallocerrdisplay
  Use file, Only: file_open,file_getunit
  Use vector, Only: VecType,vector_display,Assignment(=),Operator(+), &
      Operator(-),Operator(*),Operator(/)
  Use matrix, Only: MatrixType,matrix_display,Assignment(=),Operator(*)
  Use hessian, Only: HessianInfo,Assignment(=),hessian_init, &
      hessian_scalehessian,hessian_dispmode,hessian_rotate,hessian_removert, &
      hessian_dispevalues,hessian_eigeninfo,hessian_clean,hessian_display
  Use atom, Only: atom_atomicnos
  Use visxyz, Only: XYZ_Entry, visxyz_dump, visxyz_make

  Implicit None
  Save

  Private
  Public :: QMPot_Hessian,readqmpot_idstring,readqmpot_hessidstring, &
      readqmpot_readhess,readqmpot_hessian,readqmpot_natoms

  Character(len=strLen), Parameter       :: readqmpot_idstring = 'QMPOT'
  Character(len=strLen), Parameter       :: readqmpot_hessidstring = 'QMPOTHESS'

  Type QMPot_Hessian
    Character(len=lstrLen)  :: filename
    Integer                 :: natoms
    Type(HessianInfo)       :: hessian
  End Type QMPot_Hessian

Contains
  !-----------------------------------------------------------------------
  ! Reads a QM-Pot Hessian file and stores the results
  ! Requires:  system -- QM-Pot Hessian data structure
  !            filename -- file to read from
  !-----------------------------------------------------------------------
  Subroutine readqmpot_readhess(system,filename)
    Type(QMPot_Hessian), Intent(InOut)     :: system
    Character(*), Intent(In)               :: filename

    Integer                       :: i,j,error,unit,lineno,dim,ios
    Character(len=lstrLen)        :: line,srchstr
    Character(len=100000)         :: longstring
    Real(kind=RDbl), Dimension(:,:), Allocatable  :: array

    !** Open the file with extra-long record length
    unit = file_getunit(filename)
    Open(file=Trim(filename), unit=unit, form="FORMATTED", &
        status="OLD", IOSTAT=ios, RECL=50000) 
    If (ios /= 0) Then
      Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
           ' Could not open Hessian file: ',Trim(filename)
      Stop
    End If

    !** Find the Hessian section
    srchstr = '#matrix'
    lineno = filesrchstr(unit,srchstr,line,.True.)
    If (lineno == 0) Then
      Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
           ' Could not find search string "',Trim(srchstr), &
           '" in file ',Trim(filename)
      Write(0,'(2x,a)') ' consider erasing some lines before the keyword'
      Stop
    End If

    !** Get the Hessian dimension
    Read(unit,*) dim

    !** Record the number of atoms
    If (Mod(dim,3) /= 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Dimension of Hessian must be divisible by 3'
      Stop
    End If
    system%natoms = Int(Real(dim)/3)

    !** Allocate the temporary array
    Allocate(array(dim,dim), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'temporary Hessian')

    !** Read the Hessian in lower triangle form
    Do i = 1,dim
      Read(unit,*,IOSTAT=ios) array(i,1:i)
      If (ios /= 0) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' Error during Hessian read from ',Trim(filename)
        Stop
      End If
    End Do

    Close(unit=unit)

    !** Generate the upper triangle
    Do i = 1,dim
      Do j = i+1,dim
        array(i,j) = array(j,i)
      End Do
    End Do

    !** Initialize the hessian data type
    Call hessian_init(system%hessian,array)

    !** Deallocate the temporary array
    Deallocate(array, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'temporary Hessian')

  End Subroutine readqmpot_readhess

  !-----------------------------------------------------------------------
  ! Just get the Hessian from the system data structure
  ! Requires:  system -- QM-Pot Hessian data structure
  !            units -- optinonal returned units
  !-----------------------------------------------------------------------
  Type(HessianInfo) Function readqmpot_hessian(system,units)
    Type(QMPot_Hessian), Intent(In)      :: system
    Character(*), Intent(Out), Optional  :: units

    readqmpot_hessian = system%hessian
    If (Present(units)) Then
      units = 'kcal/mol/Ang^2'
!      Write(units,'(4a)') Trim(system%energy_units),'/', &
!          Trim(system%length_units),'^2'
    End If

  End Function readqmpot_hessian

  !-----------------------------------------------------------------------
  ! Just get the Hessian from the system data structure
  ! Requires:  system -- QM-Pot Hessian data structure
  !-----------------------------------------------------------------------
  Integer Function readqmpot_natoms(system)
    Type(QMPot_Hessian), Intent(In)     :: system

    readqmpot_natoms = system%natoms

  End Function readqmpot_natoms

  !-----------------------------------------------------------------------
  ! Display the structure
  ! Requires:  system -- QM-Pot Hessian data structure
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !-----------------------------------------------------------------------
  Subroutine readqmpot_display(system,indent,unit)
    Type(QMPot_Hessian), Intent(In)     :: system
    Integer, Intent(In)                 :: indent,unit

    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    Write(unit,'(3a)') blank," QM-Pot Hessian taken from: ", &
        Trim(system%filename)
    Write(unit,'(2a,i3)') blank," Number of atoms: ",system%natoms
    Write(unit,'(2a)') blank," The Hessian: "
    Call hessian_display(system%hessian,indent+2,unit)

  End Subroutine readqmpot_display

  !-----------------------------------------------------------------------
  ! Clean the structure
  ! Requires:  system -- QM-Pot Hessian data structure
  !-----------------------------------------------------------------------
  Subroutine readqmpot_clean(system)
    Type(QMPot_Hessian), Intent(InOut)     :: system

    Call hessian_clean(system%hessian)

  End Subroutine readqmpot_clean


End Module readqmpot
