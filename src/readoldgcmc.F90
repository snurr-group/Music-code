!--------------------------------------------------------------------------
! This module handles the reading and storage of restart file from Louie's
! old GCMC code.  Maybe some other stuff too.  This is a separate module
! because copying some data type definitions from the old code makes things 
! easier.
!--------------------------------------------------------------------------

Module readoldgcmc

  Use defaults, Only: strLen,lstrLen,RDbl,MAX_SORBS,calToj
  Use utils, Only: allocErrDisplay,deallocErrdisplay,toupper,real2str,int2str
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), &
      Operator(-),Operator(*),vector_display
  Use matrix, Only: MatrixType, Assignment(=), Operator(*)
  Use molecules, Only: molecules_getnatoms,molecules_getcoords, &
      molecules_getcomdefn
  Use readstruc, Only: Structure,readstruc_set,readstruc_setfromspc, &
      readstruc_natoms,readstruc_visxyz,readstruc_clean
  Use match, Only: match_compare,match_display

  Implicit None
  Save

  Private
  Public :: OldGCMCSystem,Energies,NewStateInfo,readoldgcmc_init, &
      readoldgcmc_readfile,readoldgcmc_chkspc,readoldgcmc_nmoles, &
      readoldgcmc_molecstate

  !** Energy structure from the old code
  Type Energies
    Real(kind=RDbl)              :: total
    Real(kind=RDbl)              :: ss,sz
    Real(kind=RDbl)              :: sslj,ssel
    Real(kind=RDbl)              :: sstail
    Real(kind=RDbl)              :: szlj,szel
    Real(kind=RDbl)              :: tor
  End Type Energies

  !** Species energy structure from the old code
  Type Species_Energies
    Real(kind=RDbl)              :: sz
    Real(kind=RDbl)              :: szlj,szel
    Real(kind=RDbl)              :: tor
  End Type Species_Energies             

  !** Rigid coordinates from the old code
  Type StateInfo
    Integer            :: mo_type
    Real(kind=RDbl)    :: x, y, z, theta, phi, psi
    Real(kind=RDbl)    :: rot(5)
    Real(kind=RDbl)    :: szpot
  End Type StateInfo

  !** Temporary Rigid coordinates for this module
  Type NewStateInfo
    Real(kind=RDbl)    :: theta, phi, psi
    Type(VecType)      :: com
  End Type NewStateInfo

  !** Old reference atomic coordinate structure
  Type AtomInfo
    Character(len=24)                     :: name
    Integer                               :: at_type
    Real(kind=RDbl)                       :: x, y, z
    Real(kind=RDbl)                       :: wt
    Real(kind=RDbl)                       :: sscharge,szcharge
  End Type AtomInfo

  !** old statistics for move types
  Type MoveType
    Integer                               :: trans
    Integer                               :: rot
    Integer                               :: mrot
    Integer                               :: swap
    Integer                               :: insert
    Integer                               :: delete
  End Type MoveType

  !** Sripped-down version of the old molecules structure
  Type Molecular_Sorbate_Params
    Character(len=24)                      :: molec_name
    Integer                                :: natoms
    Type(AtomInfo), Dimension(:), Pointer  :: atoms

    !** stuff that appears in the restart file
    Integer                                :: nmoles     
    Type(Species_Energies)                 :: nrg
    Real                                   :: delta,deleulr
    Integer                                :: cumtotal
    Integer                                :: cumtotalsq
    Integer                                :: lastcumtotal
    Type(MoveType)                         :: attempts
    Type(MoveType)                         :: accepted
    Integer                                :: inwall
    Integer                                :: izero
    Real                                   :: fugacity,endfugacity
  End Type Molecular_Sorbate_Params

  Type OldGCMCSystem
    Integer                         :: total_nmoles,nspcs
    Type(Energies)                  :: oldenergies
    Type(NewStateInfo), Dimension(:,:), Pointer  :: state
  End Type OldGCMCSystem
  
  Type(Molecular_Sorbate_Params), Dimension(:), Pointer :: mo_params
  Type(StateInfo), Dimension(:), Allocatable    :: state

Contains

  !-------------------------------------------------------------------------
  ! Read the set of molecular reference configurations and names.  Then,
  ! read in the given restartfile.
  ! Requires: system -- the old GCMC system configuration data structure
  !           modefn_file -- molecular definition filename
  !           restart_file -- restart filename to read from
  !-------------------------------------------------------------------------
  Subroutine readoldgcmc_init(system,modefn_file,restart_file)
    Type(OldGCMCSystem), Intent(InOut)           :: system
    Character(*), Intent(In)                     :: modefn_file,restart_file

    Integer                     :: i,unit,error
    Integer                     :: nspc,spc
    Character(len=lstrLen)      :: spcnames

    !** make sure this hasn't already been done
    If (Associated(mo_params)) Then  
      Return
    End If

    unit = file_open(modefn_file,110)

    !** First line is the number of species
    Read(unit,*) nspc

    !** allocate memory
    Allocate(mo_params(nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'mo_params')

    !** Second line contains a list of species names (not currently used)
    Read(unit,*) spcnames

    Do spc = 1,nspc
      Read(unit,*) mo_params(spc)%molec_name

      Read(unit,*) mo_params(spc)%natoms
      Allocate(mo_params(spc)%atoms(mo_params(spc)%natoms), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      Do i = 1,mo_params(spc)%natoms
        Read(unit,*) mo_params(spc)%atoms(i)%name
        Read(unit,*) mo_params(spc)%atoms(i)%x,mo_params(spc)%atoms(i)%y, &
            mo_params(spc)%atoms(i)%z
        Read(unit,*) mo_params(spc)%atoms(i)%sscharge, &
            mo_params(spc)%atoms(i)%szcharge
      End Do
    End Do

    Close(unit=unit)

    !** read in the actual restart file
    Call readoldgcmc_readfile(system,restart_file)

    !** display the energies contained in the file
    Call readoldgcmc_displaynrg(system,2,6)

  End Subroutine readoldgcmc_init

  !-------------------------------------------------------------------------
  ! Read the restart file and store all the information for later access.
  ! The mo_params structure must first be initialized or this routine will
  ! not work.
  ! Requires: system -- the old GCMC system configuration data structure
  !           filename -- restart filename to read from
  !-------------------------------------------------------------------------
  Subroutine readoldgcmc_readfile(system,filename)
    Type(OldGCMCSystem), Intent(InOut)    :: system
    Character(*), Intent(In)              :: filename

    Integer              :: i,unit,error
    Integer              :: istep,iseed,mo_const,spc
    Integer, Dimension(MAX_SORBS)  :: counter

    unit = file_open(filename)

    Read(unit,*) istep, iseed
    Read(unit,*) system%oldenergies
    Read(unit,*) 
    Read(unit,*) mo_const

    Allocate(state(mo_const), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1,mo_const
       Read(unit,*) state(i)
    End Do

    Do i = 1,Size(mo_params)
       Read(unit,*) mo_params(i)%nmoles
       Read(unit,*) mo_params(i)%nrg
       Read(unit,*) mo_params(i)%delta
       Read(unit,*) mo_params(i)%deleulr
       Read(unit,*) mo_params(i)%cumtotal
       Read(unit,*) mo_params(i)%cumtotalsq
       Read(unit,*) mo_params(i)%lastcumtotal
       Read(unit,*) mo_params(i)%attempts
       Read(unit,*) mo_params(i)%accepted
       Read(unit,*) mo_params(i)%inwall
       Read(unit,*) mo_params(i)%izero
       Read(unit,*) mo_params(i)%fugacity
       Read(unit,*) mo_params(i)%endfugacity
    End Do    

    Close(unit=unit)

    !** Put basic information into the system configuration structure
    system%total_nmoles = mo_const
    system%nspcs = Size(mo_params)

    !** allocate the system's generalized coordinates for each species
    Do spc = 1,system%nspcs
      Allocate(system%state(spc,mo_params(spc)%nmoles), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    
    End Do
    
    !** Put the state information into the system configuration structure
    counter = 0
    Do i = 1,system%total_nmoles
      spc = state(i)%mo_type
      counter(spc) = counter(spc) + 1
      system%state(spc,counter(spc))%theta = state(i)%theta
      system%state(spc,counter(spc))%phi = state(i)%phi
      system%state(spc,counter(spc))%psi = state(i)%psi
      system%state(spc,counter(spc))%com = (/state(i)%x,state(i)%y,state(i)%z/)
    End Do

  End Subroutine readoldgcmc_readfile

  !-------------------------------------------------------------------------
  ! Returns the number of molecules of a particular species type
  ! Requires: system -- the old GCMC system configuration data structure
  !-------------------------------------------------------------------------
  Integer Function readoldgcmc_nmoles(spc,system)
    Type(OldGCMCSystem), Intent(In)  :: system
    Integer, Intent(In)              :: spc

    readoldgcmc_nmoles = mo_params(spc)%nmoles

  End Function readoldgcmc_nmoles

  !-------------------------------------------------------------------------
  ! Returns the rotation state of a specific molecule in the configuration
  ! Requires: system -- the old GCMC system configuration data structure
  !           molec -- the molecule number 
  !           spc -- the species number
  !           theta, phi, psi -- the euler angles from reference state
  !           com - the center of mass position vector
  !           xyzcoords -- xyzcoords consistent with generalized coords
  !-------------------------------------------------------------------------
  Subroutine readoldgcmc_molecstate(system,molec,spc,theta,phi,psi,com,xyzcoords)
    Type(OldGCMCSystem), Intent(In)          :: system
    Integer, Intent(In)                      :: spc,molec
    Real(kind=RDbl), Intent(Out)             :: theta,phi,psi
    Type(VecType), Intent(Out)               :: com   
    Type(VecType), Dimension(:), Intent(Out) :: xyzcoords

    Integer                               :: natoms,error
    Type(VecType), Dimension(:), Pointer  :: refcoords

    theta = system%state(spc,molec)%theta
    phi = system%state(spc,molec)%phi
    psi = system%state(spc,molec)%psi
    com = system%state(spc,molec)%com

    !** Check size of xyzcoords
    natoms = mo_params(spc)%natoms
    If (Size(xyzcoords) < natoms) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed xyzcoords not large enough'
      Stop
    End If

    !** get reference coordinates
    Allocate(refcoords(natoms), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Call molecules_getcomdefn(spc,refcoords)

    !** get xyzcoords
    Call readoldgcmc_state2xyz(system%state(spc,molec),refcoords,xyzcoords)

    Deallocate(refcoords, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)

  End Subroutine readoldgcmc_molecstate

  !-------------------------------------------------------------------------
  ! Returns the xyz coords given the OLD molecule state coords.  This has to
  ! be done here because a different convention was used in the old code.
  ! Requires: state -- old state coordinates
  !           refcoords -- reference coordinates of molecule's atoms
  !           xyzcoords -- new xyz coordinates
  !-------------------------------------------------------------------------
  Subroutine readoldgcmc_state2xyz(state,refcoords,xyzcoords)
    Type(NewStateInfo), Intent(In)             :: state
    Type(VecType), Dimension(:), Intent(In)    :: refcoords
    Type(VecType), Dimension(:), Intent(Out)   :: xyzcoords

    Integer                  :: a
    Real(kind=RDbl)          :: cospsi, cosphi1
    Real(kind=RDbl)          :: costhe, sinpsi, sinphi1, sinthe
    Type(MatrixType)         :: rotmtx
    Real(kind=RDbl), Dimension(3,3) :: mtx

    cospsi = cos(state%psi)
    cosphi1 = cos(state%phi)
    costhe = cos(state%theta)
    sinpsi = sin(state%psi)
    sinphi1 = sin(state%phi)
    sinthe = sin(state%theta)
    
    !*** June convention - opposite Goldstein's
    mtx(1,1) = cospsi*cosphi1-sinpsi*costhe*sinphi1
    mtx(1,2) = cospsi*sinphi1+sinpsi*costhe*cosphi1
    mtx(1,3) = sinpsi*sinthe
    mtx(2,1) = -sinpsi*cosphi1-cospsi*costhe*sinphi1
    mtx(2,2) = -sinpsi*sinphi1+cospsi*costhe*cosphi1
    mtx(2,3) = cospsi*sinthe
    mtx(3,1) = sinthe*sinphi1
    mtx(3,2) = -sinthe*cosphi1
    mtx(3,3) = costhe
    
    rotmtx = mtx

    Do a = 1,Size(refcoords)
      xyzcoords(a) = rotmtx*refcoords(a)
      xyzcoords(a) = xyzcoords(a) + state%com
    End Do

  End Subroutine readoldgcmc_state2xyz

  !-------------------------------------------------------------------------
  ! Compare the reference coordinates
  ! Requires: spc -- species to check
  !-------------------------------------------------------------------------
  Subroutine readoldgcmc_chkspc(spc)
    Integer, Intent(In)        :: spc

    Integer                    :: nomatch_atom,error,i
    Integer                    :: natoms1,natoms2
    Real(kind=RDbl)            :: tol
    Character(len=strLen)      :: atomname,xyzfile,string
    Type(Structure)            :: struc1,struc2
    Type(VecType), Dimension(:), Pointer     :: atomcoords1
    Character(len=2), Dimension(:), Pointer  :: elements1

    !** get the structure in mo_params
    natoms1 = mo_params(spc)%natoms
    Allocate(atomcoords1(natoms1), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(elements1(natoms1), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1,natoms1
      atomcoords1(i)%comp = (/mo_params(spc)%atoms(i)%x, &
          mo_params(spc)%atoms(i)%y,mo_params(spc)%atoms(i)%z/)
      atomname = mo_params(spc)%atoms(i)%name
      If (Index(ToUpper(Trim(atomname)),'CARBON') /= 0) Then
        elements1(i) = 'C '
      Else If (Index(ToUpper(Trim(atomname)),'HYDROGEN') /= 0) Then
        elements1(i) = 'H '
      Else
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' Unable to identify atom: ',Trim(atomname)
        Stop
      End If
    End Do

    !** place into a structure data type
    Call readstruc_set(struc1,atomcoords1,elements1)

    !** get the reference molecule coordinates, those in molecule structure
    Call readstruc_setfromspc(struc2,spc)
    natoms2 = readstruc_natoms(struc2)
    If (natoms2 /= natoms1) Then
      Write(0,'(2a,i4,a,2i5)') __FILE__,": ",__LINE__, &
          ' Number of atoms does not match: ',natoms1,natoms2
    End If

    !** compare the two sets of coordinates
    tol = 1.0e-6_RDbl
    nomatch_atom = match_compare(struc1,struc2,tol,.False.)
    If (nomatch_atom /= 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Molecule and Read Molecule structures do not match within: ', &
          Trim(real2str(tol,5))
      Write(0,'(a)') 'Structure 1 = read Molecule,  Structure 2 = internal'
      string = int2str(nomatch_atom)
      Write(0,'(3a)') 'Could not match atom ',Trim(string), &
          ' of the second structure to an atom in the first'
      Call match_display(struc1,struc2,6)
      xyzfile = 'consistent.xyz'
      Write(0,'(2a)') 'Dumping coordinates that are consistent &
          &with internal MUSIC structures to: ',Trim(xyzfile)
      Call readstruc_visxyz(struc2,xyzfile)
      Stop      
    End If    

    !** Deallocate coordinate arrays and structures
    Deallocate(atomcoords1, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    Deallocate(elements1, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    Call readstruc_clean(struc1)
    Call readstruc_clean(struc2)

  End Subroutine readoldgcmc_chkspc

  !-------------------------------------------------------------------------
  ! Display the energies that were read from the old restart file
  ! Requires: system -- the old GCMC system configuration data structure
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !-------------------------------------------------------------------------
  Subroutine readoldgcmc_displaynrg(system,indent,unit)
    Type(OldGCMCSystem), Intent(InOut)    :: system
    Integer, Intent(In)                   :: indent,unit

    Integer                     :: spc
    Real(kind=RDbl)             :: s
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)
    s = calToj
!LC    s = 4.1868_RDbl  !* this is what the old code used

    Write(unit,'(4a)') blank,"Old GCMC restart file total energy: ", &
        Trim(real2str(system%oldenergies%total*s,8)),'    all units are kJ/mol'
    Write(unit,'(2a,3(a12))') blank,'Overall Sorbate-Zeolite: ','SZ_Total', &
        'SZ_LJ','SZ_el'
    Write(unit,'(a,25x,3(4x,a))') blank, &
        Trim(real2str(system%oldenergies%sz*s,8)), &
        Trim(real2str(system%oldenergies%szlj*s,8)), &
        Trim(real2str(system%oldenergies%szel*s,8))

    Write(unit,'(2a,5(a12))') blank,'Overall Sorbate-Sorbate: ','SS_Total', &
        'SS_LJ','SS_el','SS_tail','Torsion'
    Write(unit,'(a,25x,5(4x,a))') blank, &
        Trim(real2str(system%oldenergies%ss*s,8)), &
        Trim(real2str(system%oldenergies%sslj*s,8)), &
        Trim(real2str(system%oldenergies%ssel*s,8)), &
        Trim(real2str(system%oldenergies%sstail*s,8)), &
        Trim(real2str(system%oldenergies%tor*s,8))

    Do spc = 1,Size(mo_params)
      Write(unit,'(a,a25,4(a12))') blank, &
          Trim(mo_params(spc)%molec_name)//' specific: ', &
          'SZ_Total','SZ_LJ','SZ_el','Torsion'
        
      Write(unit,'(a,25x,4(4x,a))') blank, &
          Trim(real2str(mo_params(spc)%nrg%sz*s,8)), &
          Trim(real2str(mo_params(spc)%nrg%szlj*s,8)), &
          Trim(real2str(mo_params(spc)%nrg%szel*s,8)), &
          Trim(real2str(mo_params(spc)%nrg%tor*s,8))
    End Do

  End Subroutine readoldgcmc_displaynrg

End Module readoldgcmc
