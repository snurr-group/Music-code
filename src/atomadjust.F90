!-----------------------------------------------------------------------------
! This module handles the adjustment of single atomic coordinates based
! on the position of some or all of the atoms in the remainder of the 
! molecule.  
!   AtomAdjust_Info 
!     [natoms_to_adjust]
!     [atom_number] RESETFROM@[reference_atom_number_sequence] [optional_keys]
!      OR
!     [atom_number] FORCEFIELDOPT@[intra_keyword_filename] [optional_keys]
! Example of one-atom control lines:
!
!   5 RESETFROM@7,8-9 PERPDUMMYATOM@7,8-9
! Atom 5 will be mapped to its original position relative to atoms 7,8,9 
! and a dummy atom placed perpendicular to to the plane of 8,9 at 7.  
!
!   5 PULLFROMPLANE@7,8-9@10@[distance]
! Atom 5 will be pulled the specified distance in the direction of atom 10
! from the plane formed by atoms 7,8-9.  This can be used to establish a
! directional preference that the forcefield optimization can then complete.
!
!   6 FORCEFIELDOPT@intramolecular_file  FTOL@10.0 STEPSIZE@0.1 RESETSTEPSIZE
! Atom 6 will be adjusted using a specified intramolecular forcefield
! AFTER the other adjustments are done.  The optimization will be done
! to a force tolerance of 10.0 kJ/mol/Ang with an initial stepsize of
! 0.1 Angstroms.  This stepsize will be reset every time a new
! optimization is performed.
!-----------------------------------------------------------------------------

Module atomadjust

  Use defaults, Only: RDbl, strLen, lstrLen, xlstrLen, calToJ
  Use utils, Only: filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, deallocErrDisplay, int2str, str2seq, findint, &
      findstr, combine, checkandstop, getpath, real2str, toreal
  Use file, Only: file_open
  Use random, Only: rranf
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getunitvec, vector_getnorm, mag, &
      vector_display,vector_crossprod,vector_getplanenorm,vector_ptonplane
  Use matrixops, Only: matrixops_ludcmp, matrixops_lubksb
  Use store, Only: Store_Level_Pair, store_zero, store_getforces, &
      store_sum, store_disp, store_clean
  Use storebase, Only: EnergyPlus,storebase_init,storebase_sumintra, &
      storebase_nrg
  Use intramolecular, Only: IntramolecularInfo, intramolecular_init, &
      intramolecular_initstore, intramolecular_molint, &
      intramolecular_display, intramolecular_clean
  Use optstep, Only: Optimization_Step,optstep_simpleinit,stpdesc_idstring, &
      optstep_stpdescmove,optstep_setstep, optstep_clean
  Use molecules, Only: molecules_getnatoms,molecules_getcoords, &
      molecules_name,molecules_getfilename,molecules_getlinewtag
  Use visxyz, Only: XYZ_Entry, visxyz_make, visxyz_dump

  Implicit None
  Save

  Private
  Public :: AtomAdjust_Params, RelativeReset, ForcefieldOptimize, &
      PullFromPlane,atomadjust_init, atomadjust_move, &
      atomadjust_display, atomadjust_clean

  Integer, Parameter       :: atomadjust_ntypes = 3
  Character(len=strLen), Dimension(atomadjust_ntypes) :: atomadjust_idstrings = &
      (/'RESETFROM    ','FORCEFIELDOPT','PULLFROMPLANE'/)  

  !** Parameters for making adjustments to individual atoms
  Type AtomAdjust_Params
    Integer               :: spc,totalnatms
    Type(RelativeReset), Pointer        :: reset
    Type(ForcefieldOptimize), Pointer   :: ffopt
    Type(PullFromPlane), Pointer        :: pull
  End Type AtomAdjust_Params

  !** Use a predefined mapping to adjust coordinates of specified atoms
  Type RelativeReset
    Integer                                  :: natms
    Integer, Dimension(:), Pointer           :: atms
    Integer, Dimension(:), Pointer           :: nrefatms
    Integer, Dimension(:,:), Pointer         :: refatms
    Real(kind=RDbl), Dimension(:,:), Pointer :: coeffs
  End Type RelativeReset

  !** Use a forcefield-based optimization to adjust specified atoms
  Type ForcefieldOptimize
    Integer                              :: natms,spc,maxsteps
    Logical                              :: resetstepsize
    Real(kind=RDbl)                      :: maxforce
    Real(kind=RDbl)                      :: refstepsize,stepsize
    Character(len=strLen)                :: intrafile
    Integer, Dimension(:), Pointer       :: atms,alist
    Type(Optimization_Step)              :: stepinfo
    Type(IntramolecularInfo)             :: iparams    !* forcefield params
    Type(Store_Level_Pair)               :: intra      !* intra nrgs+forces
    Type(VecType), Dimension(:), Pointer :: force
  End Type ForcefieldOptimize

  !** Pull a specified atom out of a plane
  Type PullFromPlane
    Integer                                  :: natms
    Logical, Dimension(:), Pointer           :: push   !** True => push instead
    Integer, Dimension(:), Pointer           :: atms,orientatms
    Integer, Dimension(:,:), Pointer         :: refatms
    Real(kind=RDbl), Dimension(:), Pointer   :: distance
  End Type PullFromPlane

Contains
  !---------------------------------------------------------------------------
  ! Initialize the atom position adjustment parameters.  Reads the information
  ! from the molecule file, or, optionally a passed filename.  See header for
  ! an example of expected input.
  ! Requires:  params -- atom adjustment parameters
  !            spc -- species number
  !            filename -- filename where initialization data can be found
  !---------------------------------------------------------------------------
  Subroutine atomadjust_init(params,spc,filename)
    Type(AtomAdjust_Params), Intent(Out)     :: params
    Integer, Intent(In)                      :: spc
    Character(*), Intent(In), Optional       :: filename

    Integer                    :: i,n,error,unit,natoms,idtype
    Integer                    :: nfields,nchunks
    Character(len=strLen)      :: string1,string2
    Character(len=xlstrLen)    :: tag,line
    Integer, Dimension(atomadjust_ntypes)               :: nofeach
    Character(len=strLen), Dimension(strLen)            :: fields,chunks
    Character(len=lstrLen), Dimension(:,:), Allocatable :: infostring

    !** Set defaults
    params%spc = spc
    Nullify(params%reset)
    Nullify(params%ffopt)
    Nullify(params%pull)

    !** Find the molecule map specifications
    tag = 'AtomAdjust_Info'
    If (Present(filename)) Then     !** get from arbitrary file
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' not prepared to read atom adjustment info from arbitrary file'
      Stop
    Else    !** get from the molecule file
      If (.Not. molecules_getlinewtag(spc,tag,line,unit)) Then
        string1 = molecules_getfilename(spc)
        Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
            ' Could not find molecule map in molecule file: ',Trim(string1)
        Write(0,'(2a)') 'Desired Tag: ',Trim(tag)
        Stop
      End If
    End If

    !** Read the number of adjustments and allocate memory
    Read(unit,'(a)') line
    line = stripcmnt(line)
    params%totalnatms = toint(line,'total number of atoms to adjust')

    Allocate(infostring(atomadjust_ntypes,params%totalnatms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        

    !** Read the adjustment lines themselves
    nofeach = 0
    Do i = 1,params%totalnatms
      Read(unit,'(a)') line
      nfields = split(line,fields)
      n = toint(fields(1),'atom number to adjust')   !** just checking

      !** Get the adjustment type and save the string
      nchunks = split(fields(2),chunks,'@')
      chunks(1) = ToUpper(Adjustl(chunks(1)))
      idtype = findstr(atomadjust_idstrings,chunks(1))
      If (idtype == 0) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            ' Could not interpret atom adjustment idstring: '
        Write(0,'(2x,5a)') '"',Trim(chunks(1)),'" in "',Trim(fields(2)),'"'
        Write(0,'(2x,2a)') 'input line was: ',Trim(line)
        Stop
      End If
      nofeach(idtype) = nofeach(idtype) + 1
      infostring(idtype,nofeach(idtype)) = combine(fields(1:nfields))
    End Do

    Do idtype = 1,atomadjust_ntypes
      If (nofeach(idtype) == 0) Cycle

      Select Case(idtype)
      Case(1)  !** reset relative position from others
        Allocate(params%reset,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)                
        Call atomadjust_initreset(params%reset,spc, &
            infostring(idtype,1:nofeach(idtype)),nofeach(idtype))

      Case(2)  !** forcefield-based optimization
        Allocate(params%ffopt,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)                
        Call atomadjust_initffopt(params%ffopt,spc, &
            infostring(idtype,1:nofeach(idtype)),nofeach(idtype))

      Case(3)  !** movement out of plane in specified direction
        Allocate(params%pull,STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)                
        Call atomadjust_initpull(params%pull,spc, &
            infostring(idtype,1:nofeach(idtype)),nofeach(idtype))

      End Select
    End Do

!    Call atomadjust_display(params,2,6)

    Deallocate(infostring,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        

  End Subroutine atomadjust_init

  !---------------------------------------------------------------------------
  ! Initialization for the reset parameters.  
  ! Requires:  params -- atom reset adjustment parameters
  !            spc -- species number
  !            infostrings -- initialization lines (correctly sized)
  !            natms -- number of atoms for which to initialize
  !---------------------------------------------------------------------------
  Subroutine atomadjust_initreset(params,spc,infostrings,natms)
    Type(RelativeReset), Intent(Out)         :: params
    Integer, Intent(In)                      :: spc
    Character(len=lstrLen), Dimension(:)     :: infostrings
    Integer, Intent(In)                      :: natms

    Integer                    :: i,j,k,n,error,unit,natoms
    Integer                    :: nfields,nchunks
    Character(len=strLen)      :: string1,string2
    Type(VecType)              :: dummy
    Integer, Dimension(20)     :: nums
    Type(VecType), Dimension(1000)           :: acoords
    Character(len=2), Dimension(1000)        :: elements
    Character(len=strLen), Dimension(strLen) :: fields,chunks

    !** Defaults
    params%natms = natms
    Nullify(params%coeffs)

    !** Allocate memory and initialize
    Allocate(params%atms(params%natms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    Allocate(params%nrefatms(params%natms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    Allocate(params%refatms(params%natms,10),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    params%refatms = 0
    Allocate(params%coeffs(params%natms,4),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    params%coeffs = 0.0_RDbl

    !** Interpret the adjustment lines themselves
    Do i = 1,params%natms
      nums = 0
      nfields = split(infostrings(i),fields)
      params%atms(i) = toint(fields(1),'atom number to adjust')
      nchunks = split(fields(2),chunks,'@')
      n = str2seq(chunks(2),nums)
      If (n > 4) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            ' Cannot use more than 4 reference atoms in AtomAdjust info'
        Stop
      End If
      params%refatms(i,1:n) = nums(1:n)
      
      !** Interpret the other keywords
      Do j = 3,nfields
        nchunks = split(fields(j),chunks,'@')
        Select Case(ToUpper(chunks(1)))
        Case('PERPDUMMYATOM')
          n = n + 1
          If (n > 4) Then
            Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
                ' Cannot use more than 4 reference atoms in AtomAdjust info'
            Write(0,'(2x,a)') 'only one many be a dummy atom'
            Stop
          End If

          !** Get the parameters for the dummy atom
          k = str2seq(chunks(2),nums)
          If (k /= 3) Then
            Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
                ' Expected to read 3 atoms numbers after PERPDUMMYATOM@'
          End If

          !** Store the atom numbers for the dummy atom as negative
          Do k = 1,3
            params%refatms(i,n+k-1) = -1*nums(k)
          End Do

        Case Default
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              ' Could not understand keyword: ',Trim(fields(j))
          Stop
        End Select

      End Do

      !** Finish up, set the number of reference atoms
      params%nrefatms(i) = n

    End Do

    !** Make sure that none of the atoms to be adjusted are used as refatms
    Do i = 1,params%natms
      Do j = 1,params%natms
        nums(1:10) = Abs(params%refatms(i,:))
        If (findint(nums,params%atms(j)) /= 0) Then
          Write(0,'(1x,2a,i4,a,i4)') __FILE__," : ",__LINE__, &
              ' Atom should not be both adjusted atom and reference atom ', &
              params%atms(j)
          Call atomadjust_dispreset(params,2,0)
!          Stop
        End If
      End Do
    End Do

    !** Create the adjustment coefficients
    natoms = molecules_getnatoms(spc)
    Call molecules_getcoords(spc,acoords(1:natoms),elements(1:natoms))
    Do i = 1,params%natms
      n = params%nrefatms(i)

      !** If there is a dummy atom, handle differently
      If (params%refatms(i,n) < 0) Then
        dummy = atomadjust_dummy(Abs(params%refatms(i,n:n+2)),acoords(1:natoms))
        Call atomadjust_initcoeffs(params%coeffs(i,1:n),params%atms(i), &
            (/params%refatms(i,1:n-1),natoms+1/),(/acoords(1:natoms),dummy/))
      Else
        Call atomadjust_initcoeffs(params%coeffs(i,1:n),params%atms(i), &
            params%refatms(i,1:n),acoords(1:natoms))
      End If
    End Do

#ifdef DEBUG
    params%coeffs(1,1:4) = 0.0_RDbl
    params%coeffs(1,1) = 1.0_RDbl
    params%nrefatms(1) = 1
    params%refatms(1,1:4) = 0
    params%refatms(1,1) = params%atms(1)
#endif

  End Subroutine atomadjust_initreset

  !---------------------------------------------------------------------------
  ! Initialization for the reset parameters.  
  ! Requires:  params -- forcefield optimization parameters
  !            spc -- species number
  !            infostrings -- initialization lines (correctly sized)
  !            natms -- number of atoms for which to initialize
  !---------------------------------------------------------------------------
  Subroutine atomadjust_initffopt(params,spc,infostrings,natms)
    Type(ForcefieldOptimize), Intent(Out)    :: params
    Integer, Intent(In)                      :: spc
    Character(len=lstrLen), Dimension(:)     :: infostrings
    Integer, Intent(In)                      :: natms

    Integer                    :: i,j,k,n,error,unit,natoms,atm
    Integer                    :: nfields,nchunks
    Character(len=strLen)      :: string1,string2
    Character(len=strLen), Dimension(strLen) :: fields,chunks

    !** Defaults
    params%maxsteps = 10000
    params%refstepsize = 5.0e-2_RDbl  !** Angstroms
    params%maxforce = 0.5_RDbl  !** kJ/mol/Ang
    params%resetstepsize = .False.  !** do not reset step for new optimization
    params%intrafile = ''
    params%spc = spc
    params%natms = natms

    !** Allocate memory and initialize
    Allocate(params%atms(params%natms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        

    !** Interpret the adjustment lines themselves
    Do i = 1,params%natms
      nfields = split(infostrings(i),fields)
      params%atms(i) = toint(fields(1),'atom number to adjust')

      !** Get the intramolecular keyword filename
      nchunks = split(fields(2),chunks,'@')
      If (params%intrafile == '') Then
        params%intrafile = chunks(2)
      Else
        If (chunks(2) /= params%intrafile) Then
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              ' multiple intramolecular files not yet supported for atomadjust'
          Stop
        End If
      End If
      
      !** Interpret the other keywords
      Do j = 3,nfields
        nchunks = split(fields(j),chunks,'@')
        Select Case(ToUpper(chunks(1)))
        Case('FORCETOL')
          params%maxforce = toreal(chunks(2),'maxforce in FORCETOL@[maxforce]')
        Case('STEPSIZE')
          params%refstepsize = toreal(chunks(2),'stepsize STEPSIZE@[size]')
        Case('RESETSTEPSIZE')
          params%resetstepsize = .True.
        Case Default
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              ' Could not understand keyword: ',Trim(fields(j))
          Stop
        End Select
      End Do
    End Do

    params%stepsize = params%refstepsize

    !** Use molecule for intramolecular filename if it's 'INMOLECULE'
    If (ToUpper(params%intrafile) == 'INMOLECULE') Then
      string1 = getpath('MOLSDIR')
      string2 = molecules_getfilename(spc)
      Write(params%intrafile,'(2a)') Trim(string1),Trim(string2)
    End If

    !** Initialize the forcefield for this species
    Call intramolecular_init(params%iparams,params%spc,params%intrafile)

    !** Initialize the intramolecular storage
    Call intramolecular_initstore(params%iparams,params%intra,spc,1,1,3)

    !** Initialize the steepest descent optimizer
    Call optstep_simpleinit(params%stepinfo,stpdesc_idstring,(/params%stepsize/))

    !** Make a mask to speed force zeroing
    natoms = molecules_getnatoms(params%spc)
    Allocate(params%alist(natoms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    params%alist = 0
    Do atm = 1,natoms
      If (findint(params%atms,atm) /= 0) params%alist(atm) = 1
    End Do

    !** Size the temporary forces array
    Allocate(params%force(natoms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        

    !** Convert maximum force to kcal/mol/Ang
    params%maxforce = params%maxforce/calToJ

  End Subroutine atomadjust_initffopt

  !---------------------------------------------------------------------------
  ! Initialization for the pull-from-plane parameters.  
  ! Requires:  params -- pull-from-plane parameters
  !            spc -- species number
  !            infostrings -- initialization lines (correctly sized)
  !            natms -- number of atoms for which to initialize
  !---------------------------------------------------------------------------
  Subroutine atomadjust_initpull(params,spc,infostrings,natms)
    Type(PullFromPlane), Intent(Out)         :: params
    Integer, Intent(In)                      :: spc
    Character(len=lstrLen), Dimension(:)     :: infostrings
    Integer, Intent(In)                      :: natms

    Integer                    :: i,j,k,n,error,unit,natoms,atm
    Integer                    :: nfields,nchunks
    Character(len=strLen)      :: string1,string2
    Character(len=strLen), Dimension(strLen) :: fields,chunks

    !** Defaults
    params%natms = natms

    !** Allocate memory and initialize
    Allocate(params%atms(params%natms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    Allocate(params%orientatms(params%natms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    Allocate(params%distance(params%natms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    Allocate(params%push(params%natms),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)        
    params%push = .False.
    Allocate(params%refatms(params%natms,3),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Interpret the adjustment lines themselves  
    !** PULLFROMPLANE@7,8-9@10@[distance]
    Do i = 1,params%natms
      nfields = split(infostrings(i),fields)
      params%atms(i) = toint(fields(1),'atom number to adjust')

      !** Split the command string and store chunks
      nchunks = split(fields(2),chunks,'@')
      n = str2seq(chunks(2),params%refatms(i,:))
      params%orientatms(i) = toint(chunks(3),'orientation atom')
      params%distance(i) = toreal(chunks(4),'pull from plane distance')
      
      !** Interpret the other keywords
      Do j = 3,nfields
        nchunks = split(fields(j),chunks,'@')
        Select Case(ToUpper(chunks(1)))
        Case('PUSH')
          params%push(i) = .True.
        Case Default
          Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
              ' Could not understand keyword: ',Trim(fields(j))
          Stop
        End Select
      End Do
    End Do

  End Subroutine atomadjust_initpull

  !---------------------------------------------------------------------------
  ! Given the order atomic coordinates of the molecule, make the adjustment
  ! to the specified atoms.
  ! Requires:  params -- atom adjustment parameters
  !            acoords -- list of all atomic position vectors in molecule
  !---------------------------------------------------------------------------
  Subroutine atomadjust_move(params,acoords)
    Type(AtomAdjust_Params), Intent(InOut)       :: params
    Type(VecType), Dimension(:), Intent(InOut)   :: acoords

    If (Associated(params%reset)) Then
      Call atomadjust_reset(params%reset,acoords)
    End If

    If (Associated(params%pull)) Then
      Call atomadjust_pull(params%pull,acoords)
    End If

    If (Associated(params%ffopt)) Then
      Call atomadjust_ffopt(params%ffopt,acoords)
    End If

  End Subroutine atomadjust_move

  !---------------------------------------------------------------------------
  ! Reset individual atom positions based on the their relative arrangement
  ! in the original molecule configuration.
  ! Requires:  params -- atom adjustment parameters
  !            acoords -- list of all atomic position vectors in molecule
  !---------------------------------------------------------------------------
  Subroutine atomadjust_reset(params,acoords)
    Type(RelativeReset), Intent(In)              :: params
    Type(VecType), Dimension(:), Intent(InOut)   :: acoords

    Integer          :: i,j,n
    Type(VecType)    :: coord,dummy

    Do i = 1,params%natms
!      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!      Write(*,*) i,'Coefficients: ',params%coeffs(i,:)

      coord = 0.0_RDbl
      n = params%nrefatms(i) 

      !** Handle the dummy atom if it's there
      If (params%refatms(i,n) < 0) Then 
        n = n - 1
        dummy = atomadjust_dummy(Abs(params%refatms(i,n+1:n+3)),acoords)
        coord = coord + params%coeffs(i,n+1)*dummy
      End If

      !** Handle the normal reference atoms
      Do j = 1,n
        coord = coord + params%coeffs(i,j)*acoords(params%refatms(i,j))
      End Do

      acoords(params%atms(i)) = coord
    End Do

  End Subroutine atomadjust_reset

  !---------------------------------------------------------------------------
  ! Do a forcefield optimization for individual atoms using a pre-specified
  ! intramolecular forcefield.
  ! Requires:  params -- atom adjustment parameters
  !            acoords -- list of all atomic position vectors in molecule
  !---------------------------------------------------------------------------
  Subroutine atomadjust_ffopt(params,acoords)
    Type(ForcefieldOptimize), Intent(InOut)      :: params
    Type(VecType), Dimension(:), Intent(InOut)   :: acoords

    Integer                  :: i,atm,natoms,lastdecrease
    Logical                  :: success,fast,converged
    Logical, SAVE            :: firsttime = .True.
    Real(kind=RDbl)          :: maxforce,unew,uold,udiff
    Character(len=strLen)    :: string1
    Type(EnergyPlus), SAVE   :: total
    Type(VecType), Dimension(Size(acoords))     :: savecoords
#ifdef DEBUG
    Character(len=2), Dimension(Size(acoords))  :: elements
    Type(XYZ_Entry), Dimension(Size(acoords))   :: vis

    Call molecules_getcoords(params%spc,savecoords,elements)a
#endif

    !** Initialize temporary total energy storage if this is the first time
    If (firsttime) Then
      Call storebase_init(total,1,.True.)
    End If

    !** Get number of atoms in molecule
    natoms = molecules_getnatoms(params%spc)

    !** Reset stepsize if desired
    If (params%resetstepsize) params%stepsize = params%refstepsize

    !** Perform the optimization steps
    converged = .False.
    lastdecrease = 0
    Do i = 1,params%maxsteps
      !** Zero the interactions storage structure
      Call store_zero(params%intra,.True.)

      !** Make the forcefield call
      fast = .True.
      success = intramolecular_molint(params%iparams, &
          params%intra%mi(1),acoords,fast)
      Call checkandstop(success,__FILE__,__LINE__, &
          ' intramolecular interactions evaluation problem in atomadjust')

      !** Get the new total energy for the molecule
      Call store_sum(params%intra,.True.,total)
      Call storebase_sumintra(total,.True.)
      unew = storebase_nrg(total)
      If (i /= 1) udiff = unew - uold

      !** Make sure the energy hasn't gone up in previous step
      If (i /= 1) Then
        If (udiff > 0.0_RDbl) Then  
          acoords = savecoords  !** reset coordinates
          params%stepsize = 0.75_RDbl*params%stepsize  !** decrease stepsize
          lastdecrease = i
!          Write(*,*) 'decreasing stepsize ',params%stepsize
        Else If ((i - lastdecrease) > 2) Then
          params%stepsize = 1.10_RDbl*params%stepsize  !** increase stepsize
!          Write(*,*) 'increasing stepsize ',params%stepsize
        End If
        Call optstep_setstep(params%stepinfo,params%stepsize)
      End If

      !** Extract forces from the storage structure
      params%force = store_getforces(params%intra,natoms,.False.,1)
      
      !** Zero forces on immobile atoms, look for maximum force
      maxforce = -1.0e10_RDbl
      Do atm = 1,natoms
        If (params%alist(atm) == 0) params%force(atm) = 0.0_RDbl
        If (mag(params%force(atm)) > maxforce) maxforce = mag(params%force(atm))
      End Do

      !** Check maximum force, exit if ok
      If (maxforce < params%maxforce) Then 
        converged = .True.
        Exit
      End If

      !** Save the coordinates, to reset if the move goes up in energy
      savecoords = acoords

#ifdef DEBUG
      !** Make an xyz movie file for debugging
      Do atm = 1,natoms
        vis(atm) = visxyz_make(acoords(atm),elements(atm))
      End Do
      Call visxyz_dump(vis,'optmovie.xyz')
#endif

      !** Call optimizer to take one step
      Call optstep_stpdescmove(params%stepinfo%stpdesc,acoords, &
          params%force,natoms)

      !** Save the previous energy so we can use in next step
      uold = unew
    End Do

    !** Make sure it converged
    If (.Not. converged) Then
      Write(0,'(2a,i4,a,i2)') __FILE__,":",__LINE__, &
          ' Unable to optimize selected atom positions within allowed nsteps'
      string1 = real2str(params%maxforce,8)
      Write(0,'(1x,3a)') 'max atomic force convergence criteria: ',&
          Trim(string1),' kcal/mol/Ang'
      Write(0,'(1x,a)') 'Force magnitudes (kcal/mol) on each atom:'
      Do atm = 1,natoms
        If (params%alist(atm) == 0) Cycle
        string1 = real2str(mag(params%force(atm)),8)
        Write(0,'(i3,2x,a)') atm,Trim(string1)
      End Do
      Write(0,'(1x,a)') 'Atom adjustment parameters:'
      Call atomadjust_dispffopt(params,2,0)
      Stop
    End If

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'required number of steps for atom optimization: ',i
    Write(*,*) 'last step size (Angstroms) ',params%stepsize
!    Stop

  End Subroutine atomadjust_ffopt

  !---------------------------------------------------------------------------
  ! Pull the chosen atoms from their planes.
  ! Requires:  params -- atom adjustment parameters
  !            acoords -- list of all atomic position vectors in molecule
  !---------------------------------------------------------------------------
  Subroutine atomadjust_pull(params,acoords)
    Type(PullFromPlane), Intent(In)              :: params
    Type(VecType), Dimension(:), Intent(InOut)   :: acoords

    Integer          :: i,j,n
    Type(VecType)    :: normal,planept,disp,avec

    Do i = 1,params%natms
      !** Get the normal to the plane
      normal = vector_getplanenorm((/acoords(params%refatms(i,1)), &
          acoords(params%refatms(i,2)),acoords(params%refatms(i,3))/))

      !** Project the atom's coordinates onto the plane using normal
      planept = vector_ptonplane(normal,acoords(params%refatms(i,1)), &
          acoords(params%atms(i)))
          
      !** Build a unit vector from planept to the orientation atom
      disp = vector_getunitvec(acoords(params%orientatms(i)) - planept)

      !** Reverse the displacement vector if we're pushing instead
      If (params%push(i)) disp = (-1.0_RDbl)*disp

      !** move the atom the specified distance along this vector from the plane
      acoords(params%atms(i)) = planept + params%distance(i)*disp
    End Do

  End Subroutine atomadjust_pull

  !---------------------------------------------------------------------------
  ! Create the coordinates for a dummy atom that is perpendicular to two
  ! given atoms at another atom.  The vector between the dummy atom and
  ! the first atom will be perpendicular to the plane formed by the last 
  ! two atoms and 10.0 Angstroms away.  The 10.0, rather than 1.0 seems to
  ! reduce the sensitivity to small changes in the relative reference atom
  ! positions.
  ! Requires:  atoms -- list of atom numbers
  !            acoords -- list of all atomic position vectors in molecule
  !---------------------------------------------------------------------------
  Type(VecType) Function atomadjust_dummy(atoms,acoords)
    Integer, Dimension(3), Intent(In)         :: atoms
    Type(VecType), Dimension(:), Intent(In)   :: acoords

    Type(VecType)       :: bvec1,bvec2,cprod

    bvec1 = acoords(atoms(2)) - acoords(atoms(1))
    bvec2 = acoords(atoms(3)) - acoords(atoms(1))
    cprod = vector_crossprod(bvec1,bvec2)
    atomadjust_dummy = 10.0_RDbl*vector_getunitvec(cprod)

  End Function atomadjust_dummy

  !--------------------------------------------------------------------------
  ! Construct a mapping that allows the coordinates of a select atom to
  ! be regenerated based on the positions of other specified atoms.
  !  atm = coeff(1)*ref(1) + coeff(2)*ref(2) + ... + coeff(4)*ref(4)
  ! where: ref() are the reference atom coordinate
  !        coeff() are the coefficients to be determined
  ! Requires:  coeffs -- coefficient to be created
  !            atm -- the atom number for which to build conversion coeffs
  !            anums -- list of reference atom numbers (only)
  !            acoords -- list of all atomic position vectors in molecule
  ! NOTE: make sure anums is sized correctly, four reference atoms for
  ! a 3D molecule, three for planar, two for linear.
  !--------------------------------------------------------------------------
  Subroutine atomadjust_initcoeffs(coeffs,atm,anums,acoords)
    Real(kind=RDbl), Dimension(:), Intent(Out) :: coeffs
    Integer, Intent(In)                        :: atm
    Integer, Dimension(:), Intent(In)          :: anums
    Type(VecType), Dimension(:), Intent(In)    :: acoords

    Integer                           :: np,dim,i,j
    Type(VecType)                     :: refcoord,alphacoord,jcoord,newcoord
    Real(kind=RDbl)                   :: nperm,norm
    Character(len=lstrLen)            :: tempstrg
    Integer, Dimension(Size(anums)-1)                       :: perm
    Real(kind=RDbl), Dimension(Size(anums)-1)               :: b
    Real(kind=RDbl), Dimension(Size(anums)-1,Size(anums)-1) :: A

    !** Basic information setup
    np = Size(anums)
    dim = np - 1  !dimension of molecule
    refcoord = acoords(anums(1))
    alphacoord = acoords(atm)

    !** Set-up system of equations (Ax=b) and solve
    Do i = 1,dim
      b(i) = alphacoord%comp(i) - refcoord%comp(i)
      Do j = 2,(dim+1)
        jcoord = acoords(anums(j)) 
        A(i,(j-1)) = jcoord%comp(i) - refcoord%comp(i)
      End Do
    End Do

    !** Solve the system of equations
    Call matrixops_ludcmp(A,dim,dim,perm,nperm)
    Call matrixops_lubksb(A,dim,dim,perm,b)

    !** Translate the information from the solution vector into Calphai's
    coeffs(1) = 1.0_RDbl

    Do i = 2,(dim + 1)
      coeffs(1) = coeffs(1) - b(i-1)
      coeffs(i) = b(i-1)
    End Do

    !** Check the normalization constraint
    norm = Sum(coeffs(1:np))
    If (Abs(norm - 1.0d0) > 1.0d-5) Then
      Write(0,'(2a,i4,a,i2)') __FILE__,":",__LINE__, &
          "problem with normalization of coefficients ",atm
      Stop
    End If

    !** Regenerate the reference atm atom coordinates as a check
    !** Problems with this regeneration can indicate that the constraints
    !** are poorly specified, for example, not divided well into the units
    newcoord = 0.0_RDbl
    Do i = 1,np
      newcoord = newcoord + &
          coeffs(i)*acoords(anums(i))
    End Do
    norm = mag(newcoord - alphacoord)

    If (Abs(norm) > 1.0d-3) Then  !NOTE: TOLERANCE HERE
      Write(0,'(2a,i4,a,2i3)') __FILE__,":",__LINE__, &
          " ERROR problem with regeneration of atm coord ",atm
      tempstrg = vector_display(alphacoord,'f12.5')
      Write(0,'(2a)') 'Original:    ',Trim(tempstrg)
      tempstrg = vector_display(newcoord,'f12.5')
      Write(0,'(2a)') 'Regenerated: ',Trim(tempstrg)
      Write(0,'(a,f12.5)') 'Norm(diff):  ',norm
      Write(0,'(a)') 'Consider increasing the number of primary atoms, up to 4'
      Write(0,'(a,i3)') 'Current number of reference atoms: ',np
!      Stop
    Else If (Abs(norm) > 1.0d-6) Then
      Write(0,'(2a,i4,a,i3)') __FILE__,":",__LINE__, &
          " WARNING potential problem with regeneration of atm coord ",atm
      Write(0,'(4x,a,i2,a,f8.5)') 'for atom ',atm, &
          ': abs(new_coord - old_coord) = ',norm
    End If

  End Subroutine atomadjust_initcoeffs

  !---------------------------------------------------------------------------
  ! Display the atom adjustment parameters
  ! Requires:  params -- atom adjustment parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine atomadjust_display(params,indent,unit)
    Type(AtomAdjust_Params), Intent(In)   :: params
    Integer, Intent(In)                   :: indent,unit

    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1

    blank = Repeat(' ',indent)

    string1 = molecules_name(params%spc)
    Write(unit,'(5a)') blank,"Atom adjustment parameters for species ", &
          Trim(string1)

    If (Associated(params%reset)) Then
      Call atomadjust_dispreset(params%reset,indent,unit)
    End If

    If (Associated(params%ffopt)) Then
      Call atomadjust_dispffopt(params%ffopt,indent,unit)
    End If    

    If (Associated(params%pull)) Then
      Call atomadjust_disppull(params%pull,indent,unit)
    End If    

  End Subroutine atomadjust_display

  !---------------------------------------------------------------------------
  ! Display the atom adjustment parameters
  ! Requires:  params -- atom adjustment parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine atomadjust_dispreset(params,indent,unit)
    Type(RelativeReset), Intent(In)   :: params
    Integer, Intent(In)               :: indent,unit

    Integer                     :: i,nref
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2,string3
    Integer, Dimension(3)       :: nums

    blank = Repeat(' ',indent)

    Do i = 1,params%natms
      nref = params%nrefatms(i)
      string1 = int2str(params%atms(i))
      Write(unit,'(3a)') blank,"Atom number: ",Trim(string1)
      If (params%refatms(i,params%nrefatms(i)) < 0) Then
        string1 = int2str(params%refatms(i,1:params%nrefatms(i)-1))
        nums = Abs(params%refatms(i,nref:nref+2))
        string2 = int2str(nums)
        Write(unit,'(1x,5a)') blank,"reference atoms: ", &
            Trim(string1),',DUMMY@',Trim(string2)
            
      Else
        Write(unit,'(1x,2a,4i3)') blank,"reference atoms: ", &
            params%refatms(i,1:params%nrefatms(i))
      End If
      If (Associated(params%coeffs)) Then
        Write(unit,'(1x,2a,4f8.3)') blank,"mapping coefficients : ", &
            params%coeffs(i,1:params%nrefatms(i))
      End If
    End Do
    
  End Subroutine atomadjust_dispreset

  !---------------------------------------------------------------------------
  ! Display the atom adjustment parameters
  ! Requires:  params -- atom adjustment parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine atomadjust_dispffopt(params,indent,unit)
    Type(ForcefieldOptimize), Intent(In)   :: params
    Integer, Intent(In)                    :: indent,unit

    Integer                     :: i,nref
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2,string3
    Integer, Dimension(3)       :: nums

    blank = Repeat(' ',indent)

    Do i = 1,params%natms
      string1 = int2str(params%atms(i))
      Write(unit,'(4a)') blank,"Atom number: ",Trim(string1), &
          '  -- adjust using forcefield optimization'
    End Do

    Write(unit,'(2a)') blank,'Parameters for optimizer:'
    string1 = int2str(params%maxsteps)
    Write(unit,'(2x,3a)') blank,'Maximum number of steps: ',Trim(string1)
    string1 = real2str(params%maxforce,8)
    Write(unit,'(2x,4a)') blank,'Maximum force on atom: ',Trim(string1), &
        ' kcal/mol/Ang'
    string1 = real2str(params%refstepsize,8)
    Write(unit,'(2x,4a)') blank,'Reference stepsize: ',Trim(string1),' Angstroms'
    string1 = real2str(params%refstepsize,8)
    Write(unit,'(2x,4a)') blank,'Current stepsize: ',Trim(string1),' Angstroms'
    If (params%resetstepsize) Then
      Write(unit,'(3x,3a)') blank,'Set size will be reset', &
          ' before each optimization'
    End If

    Write(unit,'(2a)') blank,'Intramolecular forcefield parameters:'
    Call intramolecular_display(params%iparams,params%spc,indent+2,unit)
    
  End Subroutine atomadjust_dispffopt

  !---------------------------------------------------------------------------
  ! Display the atom adjustment parameters
  ! Requires:  params -- atom adjustment parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !---------------------------------------------------------------------------
  Subroutine atomadjust_disppull(params,indent,unit)
    Type(PullFromPlane), Intent(In)   :: params
    Integer, Intent(In)               :: indent,unit

    Integer                     :: i,nref
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2,string3
    Integer, Dimension(3)       :: nums

    blank = Repeat(' ',indent)

    Do i = 1,params%natms
      string1 = int2str(params%atms(i))
      Write(unit,'(3a)') blank,"Atom number: ",Trim(string1)
      string1 = int2str(params%refatms(i,1:3))
      string2 = real2str(params%distance(i),4)
      string3 = int2str(params%orientatms(i))
        Write(unit,'(1x,3a,2x,3a)') blank,"pulled from plane formed by ", &
            Trim(string1),Trim(string2),' Angstroms towards atom ',Trim(string3)
    End Do
    
  End Subroutine atomadjust_disppull

  !---------------------------------------------------------------------------
  ! Clean the atom adjustment parameters
  ! Requires:  params -- atom adjustment parameters
  !---------------------------------------------------------------------------
  Subroutine atomadjust_clean(params)
    Type(AtomAdjust_Params), Intent(InOut)   :: params

    Integer                     :: error

    If (Associated(params%reset)) Then
      Deallocate(params%reset%atms,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
      Deallocate(params%reset%nrefatms,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
      Deallocate(params%reset%refatms,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
      Deallocate(params%reset%coeffs,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
    End If

    If (Associated(params%ffopt)) Then
      Deallocate(params%ffopt%atms,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
      Deallocate(params%ffopt%alist,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
      Call optstep_clean(params%ffopt%stepinfo)
      Call intramolecular_clean(params%ffopt%iparams)
      Call store_clean(params%ffopt%intra)
      Deallocate(params%ffopt%force,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
    End If

    If (Associated(params%pull)) Then
      Deallocate(params%pull%atms,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
      Deallocate(params%pull%distance,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
      Deallocate(params%pull%refatms,STAT=error)
      If (error/=0) Call DeallocErrDisplay(__FILE__,__LINE__)        
    End If    
    
  End Subroutine atomadjust_clean

End Module atomadjust
