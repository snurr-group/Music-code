!-----------------------------------------------------------------------
! This module handles the data type and the processing of the second
! derivative matrix of the potential energy, the Hessian.  Currently,
! routines include those that remove the translational and rotational 
! components as well as those that find eigenvectors and eigenvalues.
!-----------------------------------------------------------------------
Module hessian

  Use defaults, Only: strLen,lstrLen,RDbl,HTokJmol,bohrToAng
  Use utils, Only: allocerrdisplay,deallocerrdisplay,int2str,real2str
  Use vector, Only: VecType,vector_crossprod,Assignment(=),Operator(+), &
      Operator(-),Operator(*),Operator(/)
  Use matrix, Only: MatrixType,Assignment(=),Operator(*), matrix_transpose
  Use gmatrix, Only: GenericMatrix,gmatrix_init,gmatrix_checkinit, &
      gmatrix_setmtx,gmatrix_eigeninfo,gmatrix_display,gmatrix_normcols, &
      gmatrix_simxform,gmatrix_projectout,gmatrix_makeorthonorm,gmatrix_copy, &
      gmatrix_clean,gmatrix_isorthonorm,gmatrix_issymm
  Use readstruc, Only: Structure,readstruc_toxyz,readstruc_natoms
  Use visxyz, Only: XYZ_Entry,visxyz_dispset

  Implicit None
  Save

  Private
  Public :: HessianInfo,Assignment(=),hessian_init,hessian_display, &
      hessian_eigeninfo,hessian_getneigen,hessian_scalehessian, &
      hessian_dispmode,hessian_dispevalues,hessian_removert,hessian_clean, &
      hessian_rotate,hessian_getmode,hessian_submtx,hessian_vismodes, &
      hessian_removetrans, hessian_vislastmodes

  Type HessianInfo
    Integer               :: dim
    Real(kind=RDbl), Dimension(:), Pointer  :: eigenvalues
    Type(GenericMatrix)                     :: hs            !the Hessian
    Type(GenericMatrix)                     :: eigenvectors  !in columns
  End Type HessianInfo

  Interface Assignment(=)
    Module Procedure hessian_equateinit
  End Interface

Contains

  !-----------------------------------------------------------------------
  ! Read the hessian from a 2D array into the data structure
  ! Requires: hessian -- the hessian data structure (uses only hessian)
  !           array -- the normal 2D array to take the hessian from
  !           orthonorm -- logical indication of orthonormality (optional)
  !-----------------------------------------------------------------------
  Subroutine hessian_init(hessian,array,orthonorm)
    Type(HessianInfo), Intent(Out)              :: hessian
    Real(kind=RDbl), Dimension(:,:), Intent(In) :: array
    Logical, Intent(In), Optional               :: orthonorm

    Integer          :: error

    If (Size(array,1) /= Size(array,2)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed array is not square, Hessian must be square'
      Stop
    End If

    hessian%dim = Size(array,1)
    If (Present(orthonorm)) Then
      hessian%hs = gmatrix_setmtx(array,orthonorm)
    Else
      hessian%hs = gmatrix_setmtx(array)
    End If

    !** allocate the eigenvalues
    Allocate(hessian%eigenvalues(hessian%dim), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    !** initialize the eigenvectors
    hessian%eigenvectors = gmatrix_init(hessian%dim,hessian%dim)

  End Subroutine hessian_init

  !-----------------------------------------------------------------------
  ! Initializes a hessian by copying from an old one
  ! Requires: old_structure -- the old hessian
  !           new_structure -- the new hessian
  !-----------------------------------------------------------------------
  Subroutine hessian_equateinit(new_structure,old_structure)
    Type(HessianInfo), Intent(Out)      :: new_structure
    Type(HessianInfo), Intent(In)       :: old_structure

    Integer            :: error

    new_structure%hs = gmatrix_copy(old_structure%hs)
    new_structure%dim = old_structure%dim
    new_structure%eigenvectors = gmatrix_copy(old_structure%eigenvectors)

    Allocate(new_structure%eigenvalues(new_structure%dim), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    new_structure%eigenvalues = old_structure%eigenvalues
    
  End Subroutine hessian_equateinit

  !-----------------------------------------------------------------------
  ! Scales each hessian element by a set factor.  This is used to change
  ! the units of the hessian.
  ! Requires: hessian -- the hessian data structure (uses only hessian)
  !           factor -- the scaling factor
  ! NOTE: this ONLY operates on the hessian NOT the eigen information
  !       HACK -- should be done in gmatrix
  !-----------------------------------------------------------------------
  Subroutine hessian_scalehessian(hessian,factor)
    Type(HessianInfo), Intent(InOut)    :: hessian
    Real(kind=RDbl), Intent(In)         :: factor

    Integer               :: i,j

    Do i = 1,hessian%dim 
      Do j = 1,hessian%dim 
        hessian%hs%mtx(i,j) = hessian%hs%mtx(i,j) * factor
      End Do
    End Do
    
  End Subroutine hessian_scalehessian

  !-----------------------------------------------------------------------
  ! Rotates the coordinate system that the Hessian is represented in
  ! Requires: hessian -- the hessian data structure (changes ONLY hessian)
  !           rotmtx -- the usual 3x3 rotation matrix
  ! This can be done on a per atom basis, i.e. take 3x3 chunks of the 
  ! Hessian and apply the transform:  new = transpose(rotmtx)*old*rotmtx
  ! NOTE: this ONLY operates on the hessian NOT the eigen information
  !       HACK -- should be done in gmatrix
  !-----------------------------------------------------------------------
  Subroutine hessian_rotate(hessian,rotmtx)
    Type(HessianInfo), Intent(InOut)    :: hessian
    Type(MatrixType), Intent(In)        :: rotmtx

    Integer               :: row,col
!!$    Logical               :: check
    Type(MatrixType)      :: old,xpose

    xpose = matrix_transpose(rotmtx)

    !** bite off 3x3 chunks of the hessian and rotate them
    !** could probably make this faster by exploiting symmetry
    Do row = 1,hessian%dim,3
      Do col = 1,hessian%dim,3
        old = hessian%hs%mtx(row:row+2,col:col+2)
        hessian%hs%mtx(row:row+2,col:col+2) = xpose*(old*rotmtx)
      End Do
    End Do

    !** check Hessian transformation
!LC    check = gmatrix_issymm(hessian%hs)
    
  End Subroutine hessian_rotate

  !-----------------------------------------------------------------------
  ! Simply returns the n smallest eigenvalues and their associated
  ! eigenvectors
  ! Requires: hessian -- the hessian data structure 
  !           num -- number of eigenvalues and eigenvectors
  !           array -- storage for the output 
  ! information is returned in array(vector_index,n). 
  !   n = 1   : eigenvalue
  !   n = 2-> : components of eigenvector
  !-----------------------------------------------------------------------
  Subroutine hessian_getneigen(hessian,num,array)
    Type(HessianInfo), Intent(In)                :: hessian
    Integer, Intent(In)                          :: num
    Real(kind=RDbl), Dimension(:,:), Intent(Out) :: array

    Integer                    :: i

    !** check the sizes of the passed array
    If (Size(array,1) < num) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed array index 1 is not large enough'
      Stop      
    End If
    If (Size(array,2) < (hessian%dim + 1)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed array index 2 is not large enough'
      Stop      
    End If

    If ((.Not.Associated(hessian%eigenvalues)).Or. &
        (.Not.gmatrix_checkinit(hessian%eigenvectors))) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Eigenvalues and eigenvectors must be first defined'
      Stop      
    End If

    Do i = 1,num
      array(i,1) = hessian%eigenvalues(i)
      array(i,2:(hessian%dim+1)) = hessian%eigenvectors%mtx(1:hessian%dim,i)
    End Do
    
  End Subroutine hessian_getneigen

  !-----------------------------------------------------------------------
  ! Returns the specified mode's eigenvalue and eigenvector
  ! Requires: hessian -- the hessian data structure 
  !           mode -- index number of mode (1...3N: lowest to highest)
  !           evalue -- returned eigenvalue
  !           evector -- returned eigenvector
  !-----------------------------------------------------------------------
  Subroutine hessian_getmode(hessian,mode,evalue,evector)
    Type(HessianInfo), Intent(In)              :: hessian
    Integer, Intent(In)                        :: mode
    Real(kind=RDbl), Intent(Out)               :: evalue
    Real(kind=RDbl), Dimension(:), Intent(Out) :: evector

    !** check the sizes of the passed array
    If (Size(evector) < hessian%dim) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed evector array is not large enough'
      Stop      
    End If

    evalue = hessian%eigenvalues(mode)
    evector(1:hessian%dim) = hessian%eigenvectors%mtx(1:hessian%dim,mode)
    
  End Subroutine hessian_getmode

  !-----------------------------------------------------------------------
  ! Remove the translational and rotational components from the hessian
  ! Requires: hessian -- the hessian data structure (uses only hessian)
  !           atomcoords -- the atomic coordinates (size must correct)
  !-----------------------------------------------------------------------
  Subroutine hessian_removert(hessian,atomcoords)
    Type(HessianInfo), Intent(InOut)         :: hessian
    Type(VecType), Dimension(:), Intent(In)  :: atomcoords

    Integer                    :: i,j,natoms,index
    Type(VecType)              :: axis
    Type(GenericMatrix)        :: pvectors
    Type(VecType), Dimension(Size(atomcoords)) :: cprods

    natoms = Size(atomcoords)
    pvectors = gmatrix_init(hessian%dim,6)
    
    !** set the unnormalized translational vectors
    Do j = 1,3
      Do i = 1,hessian%dim,3
        pvectors%mtx(i+j-1,j) = 1.0_RDbl
      End Do
    End Do

    !** set the unnormalized rotational vectors
    Do j = 1,3
      !** form appropriate axis vector
      axis%comp = 0.0_RDbl
      axis%comp(j) = 1.0_RDbl
      !** use cross-product to get normal instantaneous rotational vector
      Do i = 1,natoms
        cprods(i) = vector_crossprod(axis,atomcoords(i))
        index = 1 + (i-1)*3
        pvectors%mtx(index:index+2,j+3) = cprods(i)%comp
      End Do
    End Do

    Call gmatrix_makeorthonorm(pvectors)
!LC    Call gmatrix_display(pvectors,'f8.3',2,6)
    Call gmatrix_projectout(hessian%hs,pvectors)

  End Subroutine hessian_removert

  !-----------------------------------------------------------------------
  ! Remove the translational components from the hessian
  ! Requires: hessian -- the hessian data structure (uses only hessian)
  !           atomcoords -- the atomic coordinates (size must correct)
  !-----------------------------------------------------------------------
  Subroutine hessian_removetrans(hessian,atomcoords)
    Type(HessianInfo), Intent(InOut)         :: hessian
    Type(VecType), Dimension(:), Intent(In)  :: atomcoords

    Integer                    :: i,j,natoms,index
    Type(VecType)              :: axis
    Type(GenericMatrix)        :: pvectors
    Type(VecType), Dimension(Size(atomcoords)) :: cprods

    natoms = Size(atomcoords)
    pvectors = gmatrix_init(hessian%dim,3)
    
    !** set the unnormalized translational vectors
    Do j = 1,3
      Do i = 1,hessian%dim,3
        pvectors%mtx(i+j-1,j) = 1.0_RDbl
      End Do
    End Do

    Call gmatrix_makeorthonorm(pvectors)
!LC    Call gmatrix_display(pvectors,'f8.3',2,6)
    Call gmatrix_projectout(hessian%hs,pvectors)

  End Subroutine hessian_removetrans

  !-----------------------------------------------------------------------
  ! Find the eigenvalues and eigenvectors of the hessian
  ! Requires:  hessian -- the hessian data structure, must be initialized
  !-----------------------------------------------------------------------
  Subroutine hessian_eigeninfo(hessian)
    Type(HessianInfo), Intent(InOut)   :: hessian

    If (.Not.(Associated(hessian%eigenvalues))) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' eigenvalue array must already be associated'
      Stop
    End If

    If (.Not.(gmatrix_checkinit(hessian%eigenvectors))) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' eigenvector storage must already be associated'
      Stop
    End If

    Call gmatrix_eigeninfo(hessian%hs,hessian%eigenvectors,hessian%eigenvalues)

  End Subroutine hessian_eigeninfo

  !-----------------------------------------------------------------------
  ! Return a sub-matrix of the Hessian for the specified "atoms"
  ! Requires:  hessian -- the hessian data structure 
  !            atomnos -- an array of atom numbers, defining sub-Hessian
  !-----------------------------------------------------------------------
  Type(HessianInfo) Function hessian_submtx(hessian,atomnos)
    Type(HessianInfo), Intent(InOut)   :: hessian
    Integer, Dimension(:), Intent(In)  :: atomnos

    Integer                          :: i,j,a,lo,hi,oldlo,oldhi
    Integer                          :: natoms,oldnatoms
    Logical                          :: foundit
    Integer, Dimension(hessian%dim)  :: map
    Real(kind=RDbl), Dimension(Size(atomnos)*3,Size(atomnos)*3) :: array

    If (Mod(hessian%dim,3) /= 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Dimension of Hessian must be divisible by 3'
      Stop
    End If

    !** Get the number of old and new atoms
    oldnatoms = Int(hessian%dim/3)
    natoms = Size(atomnos)

    !** Make a mapping from each index of the old Hessian to the new
    lo = -2
    Do a = 1,oldnatoms
      !** Get the indices in the old hessian
      oldlo = 1 + (a-1)*3
      oldhi = oldlo + 2
      map(oldlo:oldhi) = 0

      !** Check to see if we should retain this atom's information
      foundit = .False.
      Do i = 1,natoms
        If (atomnos(i) == a) foundit = .True.
        If (foundit) Exit
      End Do

      !** Handle the retainment
      If (foundit) Then
        lo = lo + 3
        map(oldlo:oldhi) = (/lo,lo+1,lo+2/)
      End If
    End Do

    !** Create the dummy array containing the new Hessian
    Do i = 1,hessian%dim
      Do j = 1,hessian%dim
        If ((map(i) == 0).Or.(map(j) == 0)) Cycle
        array(map(i),map(j)) = hessian%hs%mtx(i,j)
      End Do
    End Do

    Call hessian_init(hessian_submtx,array)

  End Function hessian_submtx

  !-----------------------------------------------------------------------
  ! Dump a movie of the structure movement along a specified mode
  ! Requires:  hessian -- the hessian data structure
  !            atomcoords -- atomic coordinates in XYZ_Entry format
  !            mode -- eigenmode to visualize
  !            filename -- filename to dump to
  !-----------------------------------------------------------------------
  Subroutine hessian_dispmode(hessian,atomcoords,mode,filename)
    Type(HessianInfo), Intent(In)              :: hessian
    Type(XYZ_Entry), Dimension(:), Intent(In)  :: atomcoords
    Integer, Intent(In)                        :: mode
    Character(*), Intent(In)                   :: filename

    Integer                                 :: i,a,nframes
    Real(kind=RDbl)                         :: norm,step,scale
    Type(VecType), Dimension(Size(atomcoords))       :: dispvecs
    Real(kind=RDbl), Dimension(mode,(hessian%dim+1)) :: array

    !** get the 1->mode eigenmodes
    Call hessian_getneigen(hessian,mode,array)

    norm = 0.0_RDbl
    Do i = 2,hessian%dim
      norm = norm + array(mode,i)**2
    End Do
    norm = Sqrt(norm)

    !** set the parameters
    scale = 0.6_RDbl  
    nframes = 30
    step = scale/nframes

    a = 0
    Do i = 2,hessian%dim,3
      a = a + 1
      dispvecs(a) = array(mode,i:i+2)
      dispvecs(a) = dispvecs(a)/norm*step
    End Do

    Call visxyz_dispset(atomcoords,dispvecs,nframes,filename,'f12.5')
      
  End Subroutine hessian_dispmode

  !---------------------------------------------------------------------------
  ! Make an xyz movie of either a specific mode of the Hessian (if the mode
  ! is specified) or off all the negative eigenmodes.  Also displays the first
  ! 10 eigenvalues to the output string.
  ! Requires:  specs -- Hessian specifications and temporary storage
  !            struc -- structure corresponding to the Hessian
  !            basefile -- the base of the filename(s) to output
  !            outstring -- returned string containing smallest 10 mode info
  !            mode -- mode to visualize
  !---------------------------------------------------------------------------
  Subroutine hessian_vismodes(hessian,struc,basefile,outstring,mode)
    Type(HessianInfo), Intent(In)       :: hessian
    Type(Structure), Intent(In)         :: struc
    Character(len=lstrLen), Intent(In)  :: basefile
    Character(*), Intent(Out)           :: outstring
    Integer, Intent(In), Optional       :: mode

    Integer                    :: i,error,natoms,nmodes,negmodes
    Real(kind=RDbl)            :: evalue
    Character(len=lstrLen)     :: filename,string1,string2
    Type(XYZ_Entry), Dimension(:), Allocatable  :: atomcoords

    !** Size the atomic coordinates array
    natoms = readstruc_natoms(struc)
    Allocate(atomcoords(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    !** Make the structure into an xyz array
    Call readstruc_toxyz(struc,atomcoords)

    !** Write the smallest 10 or less eigenvalues to the output string
    nmodes = 10
    If (hessian%dim < 10) nmodes = hessian%dim
    string1 = int2str(nmodes)
    Write(outstring,'(3a)') 'Smallest ',Trim(string1),' eigenvalues: '
    Do i = 1,nmodes
      string1 = real2str(hessian%eigenvalues(i),8)
      string2 = int2str(i)
      Write(outstring,'(a,1x,4a)') Trim(outstring),Trim(string1),' (', &
          Trim(string2),')'
      If (i < nmodes) Write(outstring,'(2a)') Trim(outstring),','
    End Do

    !** Make the movie(s) of the mode(s)
    If (Present(mode)) Then
      string1 = int2str(mode)
      Write(filename,'(2a)') Trim(basefile),'.xyz'
      Write(*,'(4x,5a)') 'creating xyz movie of mode ',Trim(string1), &
          ' in "',Trim(filename),'"'
      Call hessian_dispmode(hessian,atomcoords,mode,filename)
    Else
      negmodes = 0
      Do i = 1,hessian%dim
        evalue = hessian%eigenvalues(i)
        If (i > 10) Exit  !** limit number of output files
        If (evalue < 0) Then
          negmodes = negmodes + 1          
          string1 = int2str(i)
          Write(filename,'(4a)') Trim(basefile),'_mode',Trim(string1),'.xyz'
          Write(*,'(4x,5a)') 'creating xyz movie of mode ',Trim(string1), &
              ' in "',Trim(filename),'"'
          Call hessian_dispmode(hessian,atomcoords,i,filename)
        End If
      End Do

      If (negmodes == 0) Then 
        Write(*,'(a)') 'No negative eigenmodes found, no movies made'
      End If

    End If

    !** Deallocate temporary storage
    Deallocate(atomcoords, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)    

  End Subroutine hessian_vismodes

  !---------------------------------------------------------------------------
  ! Make xyz movies of a specified number of the highest eigenmodes.
  ! Requires:  specs -- Hessian specifications and temporary storage
  !            struc -- structure corresponding to the Hessian
  !            basefile -- the base of the filename(s) to output
  !            outstring -- returned string containing mode info
  !            nmodes -- number of highest eigenvalue eigenmodes to visualize
  !---------------------------------------------------------------------------
  Subroutine hessian_vislastmodes(hessian,struc,basefile,outstring,nmodes)
    Type(HessianInfo), Intent(In)       :: hessian
    Type(Structure), Intent(In)         :: struc
    Character(len=lstrLen), Intent(In)  :: basefile
    Character(*), Intent(Out)           :: outstring
    Integer, Intent(In)                 :: nmodes

    Integer                    :: i,error,natoms,mode,negmodes,lastnmodes
    Real(kind=RDbl)            :: evalue
    Character(len=lstrLen)     :: filename,string1,string2
    Type(XYZ_Entry), Dimension(:), Allocatable  :: atomcoords

    !** Size the atomic coordinates array
    natoms = readstruc_natoms(struc)
    Allocate(atomcoords(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

    !** Make the structure into an xyz array
    Call readstruc_toxyz(struc,atomcoords)

    !** Write the smallest 10 or less eigenvalues to the output string
    lastnmodes = nmodes
    If (hessian%dim < 10) lastnmodes = hessian%dim
    string1 = int2str(lastnmodes)
    Write(outstring,'(3a)') 'Largest ',Trim(string1),' eigenvalues: '
    Do i = (hessian%dim - lastnmodes),hessian%dim
      string1 = real2str(hessian%eigenvalues(i),8)
      string2 = int2str(i)
      Write(outstring,'(a,1x,4a)') Trim(outstring),Trim(string1),' (', &
          Trim(string2),')'
      If (i < lastnmodes) Write(outstring,'(2a)') Trim(outstring),','
    End Do

    !** Make the movie(s) of the mode(s)
    Do i = (hessian%dim - lastnmodes),hessian%dim
      evalue = hessian%eigenvalues(i)
      string1 = int2str(i)
      Write(filename,'(4a)') Trim(basefile),'_mode',Trim(string1),'.xyz'
      Write(*,'(4x,5a)') 'creating xyz movie of mode ',Trim(string1), &
          ' in "',Trim(filename),'"'
      Call hessian_dispmode(hessian,atomcoords,i,filename)
    End Do

    !** Deallocate temporary storage
    Deallocate(atomcoords, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)    

  End Subroutine hessian_vislastmodes

  !-----------------------------------------------------------------------
  ! Displays the eigenvalues in the Hessian data structure
  ! Requires: hessian -- the Hessian data structure
  !           indent -- indentation from left margin
  !           unit -- unit to dump into
  !-----------------------------------------------------------------------
  Subroutine hessian_dispevalues(hessian,indent,unit)
    Type(HessianInfo), Intent(In)    :: hessian
    Integer, Intent(In)              :: indent,unit

    Integer                     :: i,j,lo,hi
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    If (Associated(hessian%eigenvalues)) Then
      Write(unit,'(2a,i3,a)') blank," Eigenvalues ",hessian%dim," :  "
      Do i = 1,hessian%dim,5
        lo = i
        hi = Min((i+4),hessian%dim)
        Write(unit,'(2x,a,10e14.5)') blank,(hessian%eigenvalues(j),j=lo,hi)
      End Do
    End If

  End Subroutine hessian_dispevalues

  !-----------------------------------------------------------------------
  ! Displays the data in the Hessian data structure
  ! Requires: hessian -- the Hessian data structure
  !           indent -- indentation from left margin
  !           unit -- unit to dump into
  !-----------------------------------------------------------------------
  Subroutine hessian_display(hessian,indent,unit)
    Type(HessianInfo), Intent(In)    :: hessian
    Integer, Intent(In)              :: indent,unit

    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    Write(unit,'(2a,i3)') blank," dimension of hessian:  ",hessian%dim

    Write(unit,'(2a)') blank," Hessian:  "
    Call gmatrix_display(hessian%hs,'f8.2',indent+2,unit)

    Call hessian_dispevalues(hessian,indent,unit)

    Write(unit,'(2a)') blank," Eigenvectors:  "
    Call gmatrix_display(hessian%eigenvectors,'f8.3',indent+2,unit)

  End Subroutine hessian_display

  !-----------------------------------------------------------------------
  ! Cleans the Hessian data structure
  ! Requires: hessian -- the Hessian data structure
  !-----------------------------------------------------------------------
  Subroutine hessian_clean(hessian)
    Type(HessianInfo), Intent(InOut)    :: hessian

    Integer               :: error

    If (Associated(hessian%eigenvalues)) Then
      Deallocate(hessian%eigenvalues, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'eigenvalues')
    End If

    Call gmatrix_clean(hessian%hs)
    Call gmatrix_clean(hessian%eigenvectors)

  End Subroutine hessian_clean

End Module hessian

