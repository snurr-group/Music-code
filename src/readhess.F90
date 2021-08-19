!------------------------------------------------------------------------
! This module handles the reading and preliminary processing of Hessians
! from external sources.  It is basically a collection of routines that
! operate on the Hessian data type in conjunction with the external 
! read modules.
!------------------------------------------------------------------------

Module readhess

  Use defaults, Only: RDbl, strLen, lstrLen, MAX_ATOMS
  Use utils, Only:  filesrchstr, split, toupper, allocErrDisplay, real2str, &
      int2str, deallocErrDisplay, str2seq, stripcmnt
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), &
      Operator(-),vector_display
  Use matrix, Only: MatrixType,Assignment(=),Operator(*),matrix_display
  Use molecules, Only: molecules_getcomdefn,molecules_getlinewtag, &
      molecules_getfilename,molecules_getnatoms,molecules_getcoords, &
      molecules_origxform,molecules_atomsymbol,molecules_name
  Use hessian, Only: HessianInfo,Assignment(=),hessian_dispevalues, &
      hessian_dispmode,hessian_clean,hessian_display,hessian_submtx, &
      hessian_rotate
  Use readgauss, Only: GaussianSystem,readgauss_getsystem,readgauss_hessian, &
      readgauss_getstruc,readgauss_natoms,readgauss_vissys,readgauss_rotatesys, &
      readgauss_getmodes,readgauss_dispmode,readgauss_clean,readgauss_idstring, &
      readgauss_display
  Use readqmpot, Only: QMPot_Hessian,readqmpot_hessidstring, &
      readqmpot_readhess,readqmpot_hessian,readqmpot_natoms
  Use readstruc, Only: Structure,Assignment(=),readstruc_init,readstruc_xform, &
      readstruc_visxyz,readstruc_natoms,readstruc_display,readstruc_clean, &
      readstruc_subcopy,readstruc_setfromspc,readstruc_toxyz
  Use match, Only: match_compare
  Use visxyz, Only: XYZ_Entry,visxyz_make

  Implicit None
  Save

  Private
  Public :: Hessian_Specs,readhess_init,readhess_chkcoords, &
      readhess_units,readhess_display,readhess_clean

  Type Hessian_Specs
    Character(len=lstrLen)            :: filetype,filename,atomno_string
    Character(len=strLen)             :: units
    Integer                           :: nhessatoms
    Integer, Dimension(:), Pointer    :: atoms
    Type(HessianInfo)                 :: full_hessian
    Type(Structure)                   :: structure
    Type(GaussianSystem), Pointer     :: g98chkpt
    Type(QMPot_Hessian), Pointer      :: qmpot
  End Type Hessian_Specs

Contains
  !---------------------------------------------------------------------------
  ! Reads the hessian from a specific file type.  Also has the option to 
  ! compare the coordinates consistent with the read Hessian to a MuSiC
  ! stored molecule type.  This allows Hessian calculations to be done 
  ! externally and checked to make sure the structures match.
  ! Requires:  specs -- Hessian specifications and temporary storage
  !            hessian -- uninitialized hessian data structure
  !            filespec -- type and filename: FILETYPE@FILENAME
  !            spc -- compare coordinates against this stored species type
  !            hcoordloc -- location of coordinates matching Hessian
  !---------------------------------------------------------------------------
  Subroutine readhess_init(specs,hessian,filespec,spc,hcoordloc)
    Type(Hessian_Specs), Intent(Out)      :: specs
    Type(HessianInfo), Intent(Out)        :: hessian
    Character(len=lstrLen), Intent(In)    :: filespec
    Integer, Intent(In), Optional         :: spc 
    Character(*), Intent(In), Optional    :: hcoordloc

    Integer                                   :: unit,error
    Integer                                   :: nfields,nchunks
    Logical                                   :: compare
    Character(len=lstrLen)                    :: line,tag
    Character(len=lstrLen), Dimension(strLen) :: fields,chunks
    Integer, Dimension(MAX_ATOMS*10)          :: nums
    Type(HessianInfo)                         :: full_hessian    

    !** Nullify the pointers in the structure
    Nullify(specs%atoms)
    Nullify(specs%g98chkpt)
    Nullify(specs%qmpot)
    specs%nhessatoms = 0

    !** Decide if the generating structures need to be compared
    compare = .False.
    If (Present(spc)) compare = .True.

    !** Process the Hessian specifications
    nfields = split(filespec,fields,'@')
    If (nfields /= 2) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Hessian location specification must be ', &
          'FILETYPE@FILENAME[:atom_numbers] with optional'
      Write(0,'(3a)') 'atom numbers (ie 1-5,8,9), only: ', &
          Trim(filespec),' was passed in'
      Stop
    End If    
    specs%filetype = fields(1)
    
    !** Check to see if a Hessian sub-matrix is specified using atom numbers
    nchunks = split(fields(2),chunks,':')
    If (nchunks == 1) Then
      specs%filename = fields(2)
      specs%atomno_string = 'Full Hessian'
      If (compare) specs%nhessatoms = molecules_getnatoms(spc)

    Else If (nchunks == 2) Then
      specs%filename = chunks(1)
      specs%atomno_string = chunks(2)

      !** Interpret the atom numbers in the string
      specs%nhessatoms = str2seq(specs%atomno_string,nums)
      Allocate(specs%atoms(specs%nhessatoms), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      specs%atoms = nums(1:specs%nhessatoms)

    Else
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' could not interpret Hessian filename: ',Trim(fields(2))
      Stop
    End If

    !** Actually read the Hessian from the specified location and process
    Select Case(ToUpper(Trim(specs%filetype)))
    Case (readgauss_idstring)
      !** Allocate the pointer
      Allocate(specs%g98chkpt, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      !** Read the whole Gaussian 98 checkpoint file
      Call readgauss_getsystem(specs%filename,specs%g98chkpt)
!      Call readgauss_display(specs%g98chkpt,2,6)
      If (.Not. compare) specs%nhessatoms = readgauss_natoms(specs%g98chkpt)

      !** Make the structure comparison and adjustments, if necessary
      If (compare) Call readhess_chkcoords(specs,spc)

      !** Pull the Hessian out of the structure
      specs%full_hessian = readgauss_hessian(specs%g98chkpt,specs%units)

    Case (readqmpot_hessidstring)
      !** Allocate the pointer
      Allocate(specs%qmpot, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      !** Read the QM-Pot Hessian file and store
      Call readqmpot_readhess(specs%qmpot,specs%filename)
      If (.Not. compare) specs%nhessatoms = readqmpot_natoms(specs%qmpot)

      !** Make the structure comparison and adjustments, if necessary
      If (compare) Then
        If (.Not. Present(hcoordloc)) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              ' Must pass coordinate location if QM-Pot Hessian to be compared' 
          Stop
        End If
        Call readhess_chkcoords(specs,spc,hcoordloc)
      Else
        specs%full_hessian = readqmpot_hessian(specs%qmpot,specs%units)
      End If

    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Could not interpret Hessian file type specification: ', &
          Trim(specs%filetype)
      Stop
    End Select    


!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Call hessian_display(specs%full_hessian,2,6)
!    Stop

    !** Decide if we need to extract and return a sub-Hessian matrix
    If (specs%atomno_string == 'Full Hessian') Then
      hessian = specs%full_hessian
    Else
      hessian = hessian_submtx(specs%full_hessian,specs%atoms)
    End If

    !** Deallocate the temporary coordinate arrays
    If (Associated(specs%g98chkpt)) Then
      Call readgauss_clean(specs%g98chkpt)
    End If

!LC    Call hessian_dispevalues(params%hessian,2,6)

  End Subroutine readhess_init

  !---------------------------------------------------------------------------
  ! Check the coordinates corresponding to the newly read Hessian against
  ! an existing species structure in MuSiC.
  ! Requires:  specs -- Hessian specifications and temporary storage
  !            spc -- compare coordinates against this stored species type
  !            hcoordloc -- location of coordinates matching Hessian
  !---------------------------------------------------------------------------
  Subroutine readhess_chkcoords(specs,spc,hcoordloc)
    Type(Hessian_Specs), Intent(InOut)    :: specs
    Integer, Intent(In)                   :: spc 
    Character(*), Intent(In), Optional    :: hcoordloc

    Integer                     :: a,indx,error,nomatch_atom
    Integer                     :: natoms,natoms1,natoms2
    Real(kind = RDbl)           :: tol
    Character(len=strLen)       :: xyzfile,string
    Type(VecType)               :: trans_xform
    Type(MatrixType)            :: rot_xform
    Type(GaussianSystem)        :: orig
    Type(Structure)             :: refstruc,struc0,struc1,struc2

    !** Get the MuSiC species structure and number of atoms 
    Call readstruc_setfromspc(struc2,spc)
    natoms2 = readstruc_natoms(struc2)

    !** Check for consistency in number of atoms
    If (natoms2 /= specs%nhessatoms) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Number of atoms in MuSiC molecular species and Hessian do not match'
      Stop
    End If

    !** Get the transformation that was used to get the reference
    !** coordinates in the molecule structure from the molecule file
    Call molecules_origxform(spc,trans_xform,rot_xform)

    !** Read and transform the structure from which the Hessian came
    Select Case(ToUpper(Trim(specs%filetype)))
    Case (readgauss_idstring)  !** structure already in Gaussian system

      !** Pull the reference structure out of Gaussian system 
      Call readgauss_getstruc(specs%g98chkpt,refstruc)

      !** Pare down the atoms in the reference structure if need be
      If (specs%atomno_string == 'Full Hessian') Then
        struc0 = refstruc  
      Else
        struc0 = readstruc_subcopy(refstruc,specs%atoms)
      End If

      !** Rotate and transform the Gaussian system
      Call readgauss_rotatesys(specs%g98chkpt,trans_xform,rot_xform)      
#ifdef DEBUG
      Write(*,*) vector_display(trans_xform,'f8.3')
      Call matrix_display(rot_xform,'f12.5',6)
      Call readgauss_getmodes(specs%g98chkpt)
      Call readgauss_dispmode(specs%g98chkpt,9,'mode_xformed.xyz')
#endif

      !** Pull out the full Gaussian system coordinates for later comparison
      Call readgauss_getstruc(specs%g98chkpt,specs%structure)

      !** Pare down the atoms in the final structure if need be
      If (specs%atomno_string == 'Full Hessian') Then
        struc1 = specs%structure
      Else
        struc1 = readstruc_subcopy(specs%structure,specs%atoms)
      End If
      natoms1 = readstruc_natoms(struc1)

    Case (readqmpot_hessidstring)
      !** Read the structure from the coordinate location
      Call readstruc_init(specs%structure,hcoordloc)
      refstruc = specs%structure  
      struc0 = specs%structure    !* untransformed final structure

      !** Rotate and transform the structure
      Call readstruc_xform(specs%structure,trans_xform,rot_xform)      
      struc1 = specs%structure
      natoms1 = readstruc_natoms(struc1)

      !** Get the full Hessian
      specs%full_hessian = readqmpot_hessian(specs%qmpot)

      !** Rotate the Hessian
      Call hessian_rotate(specs%full_hessian,rot_xform)

    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Could not interpret Hessian file type specification: ', &
          Trim(specs%filetype)
      Stop
    End Select    

    !** Check the number of atoms
    natoms = readstruc_natoms(specs%structure)
    If (natoms /= specs%nhessatoms) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Number of atoms in structure does not match that in Hessian'
      Stop
    End If

    !** Compare the two sets of coordinates
    tol = 1.0e-6_RDbl
    nomatch_atom = match_compare(struc1,struc2,tol,.False.)
    If (nomatch_atom /= 0) Then
      string = real2str(tol,5)
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Molecule and Hessian structures do not match within: ', &
          Trim(string)
      string = int2str(nomatch_atom)
      Write(0,'(3a)') 'Could not match atom ',Trim(string), &
          ' of the second structure to an atom in the first'
      xyzfile = 'consistent.xyz'
      Write(0,'(2a)') 'Dumping coordinates that are consistent &
          &with Hessian to: ',Trim(xyzfile)
      Call readstruc_visxyz(struc0,xyzfile,'f12.5','Hessian consistent coords')
!      nomatch_atom = match_compare(struc1,struc2,tol,.True.)      
      Stop      
    End If

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call readgauss_getsystem(params%hessian_file,orig)
    Call readgauss_getmodes(orig)
    Call readgauss_dispmode(orig,9,'mode_orig.xyz')
    stop
#endif

    !** Deallocate the temporary structures
    Call readstruc_clean(refstruc)
    Call readstruc_clean(struc0)
    Call readstruc_clean(struc1)
    Call readstruc_clean(struc2)

  End Subroutine readhess_chkcoords

  !---------------------------------------------------------------------------
  ! Returns the units of the Hessian
  ! Requires:  specs -- Hessian specifications and temporary storage
  !---------------------------------------------------------------------------
  Function readhess_units(specs)
    Character(len=strLen)                 :: readhess_units
    Type(Hessian_Specs), Intent(In)       :: specs

    readhess_units = specs%units

  End Function readhess_units

  !---------------------------------------------------------------------------
  ! Display the specification information for the Hessian read method
  ! Requires:  specs -- Hessian specifications and temporary storage
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine readhess_display(specs,indent,unit)
    Type(Hessian_Specs), Intent(In)       :: specs
    Integer, Intent(In)                   :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2

    blank = Repeat(' ',indent)

    Write(unit,'(5a)') blank,'Hessian information taken from file: ', &
        Trim(specs%filename),' of type: ',Trim(specs%filetype)
    If (specs%atomno_string /= 'Full Hessian') Then
      Write(unit,'(4a)') blank,'only extracted sub-Hessian formed by ', &
          'atom numbers: ',Trim(specs%atomno_string)
    End If

  End Subroutine readhess_display

  !---------------------------------------------------------------------------
  ! Clean the specification information for the Hessian read method
  ! Requires:  specs -- Hessian specifications and temporary storage
  !---------------------------------------------------------------------------
  Subroutine readhess_clean(specs)
    Type(Hessian_Specs), Intent(InOut)       :: specs

    Integer            :: error

    Deallocate(specs%atoms, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'atoms')
    Call hessian_clean(specs%full_hessian)
    If (Associated(specs%g98chkpt)) Then
      Call readgauss_clean(specs%g98chkpt)
    End If

  End Subroutine readhess_clean

End Module readhess
