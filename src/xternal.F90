!-------------------------------------------------------------------------------
! This module handles the evaluation of energy and other molecular interactions
! using external programs.  The basic idea is that either of:
!
!   xternal_intraint:  jump point for intramolecular interaction calculations
!   xternal_int:  jump point for species-species interaction calculations
!
! are called from another module, such as intramolecular or ssdriver.  These
! two routines then compile dynamic information such as atom coordinates,
! atom types, charges, etc. into temporary lists and pass these lists to
! the module that handles a specific external program:
!
!   Supported external programs:
!   1) GULP (Julian Gale) using gulpio module
!   2) Gaussian98 using gaussio module
!   3) Turbomole (Ahlrichs group) using turboio module
!
! Also passed are static chunks of information that are stored in data types
! maintained by these external program modules.  This information is 
! initialized during the forcefield initialization process.
!
! It is important to remember that one can only expect a certain amount of
! interactions detail from a generic external program.  For example, from 
! quantum chemistry programs, it is only possible to get the whole energy
! of the set of atoms passed, the forces or gradients on these atoms and
! maybe the second derivatives.  For this reason, all types of requests
! can go through a single routine, "xternal_send", that sends off the 
! dynamic information and receives an energy and maybe forces and second
! derivatives. 
!
! Evaluating partial system interactions with an external program is quite
! problematic.  For example, to get just the interactions between one species
! and another, without intra-species interactions, one needs to first put all
! the atoms in the simulation cell, then substract the interactions found by
! only putting individual species in the simulation cell.  Alternatively,
! one can turn off intramolecular interactions, if possible.
!
! Needed Improvements:
! 1) add force functionality
!-------------------------------------------------------------------------------

Module xternal

  Use defaults, Only: RDbl, strLen, lstrLen, xlstrlen, &
      EXTERNAL_INDEX, TOTAL_INDEX, NO_OF_INTRA_POTS, calToJ
  Use utils, Only: split, stripcmnt, findint, tolower, getpath, filesrchstr, &
      toint, toupper, toreal, combine, allocErrDisplay, deallocErrDisplay
  Use file, Only: file_getunit, file_getname, file_open
  Use vector, Only: VecType
  Use config, Only: AtMolCoords,config_getnatoms,config_nsubsetatoms, &
      config_getnmoles,config_getr,config_getatype
  Use molecules, Only: molecules_getcharge,molecules_getnatoms, &
      molecules_getatype,molecules_gettype,molecules_getcharges
  Use simcell, Only: SimCell_Params
  Use storetop, Only: storetop_level
  Use storebase, Only: storebase_disp,storebase_inc
  Use store, Only: Store_Level_Pair,store_display,store_nderivs,store_idbranch, &
      store_disp
  Use storesym, Only: Symmetric_Results, storesym_ptrs, storesym_nderivs, &
      storesym_display
  Use gulpio, Only: GULP_Forcefield_Input, gulpio_snglpt, gulpio_init, &
      gulpio_display, gulpio_clean
  Use gaussio, Only: G98_Input, gaussio_snglpt, gaussio_init, gaussio_display, &
      gaussio_clean
  Use turboio, Only: Turbomole_Input, turboio_snglpt, turboio_init, &
      turboio_display, turboio_clean

  Implicit None
  Save

  Private
  Public :: ExternalInfo, EXTERNAL_KEY, xternal_init, xternal_display, &
      xternal_intraint, xternal_int

  Type ExternalInfo
    Logical             :: shift
    Real(kind=RDbl)     :: nrgshift
    Type(GULP_Forcefield_Input), Pointer   :: gulp
    Type(G98_Input), Pointer               :: g98
    Type(Turbomole_Input), Pointer         :: turbo
  End Type ExternalInfo

  Character(len=strLen), Parameter  :: EXTERNAL_KEY = 'EXTERNAL'

  !** storage for temporary configuration to be passed to external programs
  Integer                                     :: temp_natoms
  Type(VecType), Dimension(:), Allocatable    :: temp_coords
  Integer, Dimension(:), Allocatable          :: temp_atypes
  Real(kind=RDbl), Dimension(:), Allocatable  :: temp_charges
  Integer, Dimension(:,:), Allocatable        :: temp_subsets
  Type(VecType), Dimension(:), Allocatable    :: temp_forces

Contains

  !----------------------------------------------------------------------------
  ! Initialize the external-evaluated interactions for a set of species.  
  ! Format of the instructions should be:
  !  'EXTERNAL@[s1]@[s2] [OPTION1]@[modifer1] [OPTION2]@[modifer1]...'
  ! The first modifer for the EXTERNAL (s1) defines the external program
  ! to be used and the second (s2) is an option for this program.
  ! Requires:  xparams -- external interaction evaluation parameters
  !            line -- instruction line from the control file
  !            spc_list -- list of species numbers
  !            simcell -- simulation cell information
  !            list -- list of atom numbers for interaction (set,1:Natoms)
  !            params -- array of strings containing type and parameters
  !----------------------------------------------------------------------------
  Subroutine xternal_init(xparams,line,spc_list,simcell,list,params)
    Type(ExternalInfo), Intent(Out)                             :: xparams
    Character(len=lstrLen), Intent(In)                          :: line
    Integer, Dimension(:), Intent(In)                           :: spc_list
    Type(SimCell_Params), Intent(In), Optional                  :: simcell
    Integer, Dimension(:,:), Intent(In), Optional               :: list
    Character(len=lstrLen), Dimension(:), Intent(Out), Optional :: params

    Integer                               :: i,j,error
    Integer                               :: nfields,nchunks,nsets,nspc
    Logical                               :: option_intraoff
    Character(len=lstrLen)                :: filename
    Character(len=lstrLen), Dimension(2)  :: filenames
    Character(len=strLen), Dimension(20)  :: fields,chunks
    Integer, Dimension(2)                 :: intraoff_spc

    !** Nullify the pointer set
    Call xternal_nullify(xparams)

    !** Set defaults
    xparams%shift = .False.
    xparams%nrgshift = 0.0_RDbl
    nspc = Size(spc_list)
    option_intraoff = .False.
    intraoff_spc = 0

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'incoming line: ',Trim(line)

    !** Interpret the line of commands, split first by spaces, then '@' chars
    nfields = split(line,fields)
    Do i = 1,nfields
      nchunks = split(fields(i),chunks,'@')

      !** Make sure the main EXTERNAL@[]@[] sequence is 1st chunk of 1st field
      If (i == 1) Then
        If (chunks(1) /= 'EXTERNAL') Then
          Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
              ' 1st chunk of 1st fields set must be "EXTERNAL"'
          Stop  
        End If 
      End If

      !** interpret the others and build-up modifier parameters
      If (i > 1) Then
        Select Case(ToUpper(chunks(1)))
        Case('SHIFT')   !** shifts energies by a fixed value (in kJ/mol)
          xparams%shift = .True.
          xparams%nrgshift = toreal(chunks(2), &
              'external energy shift value in kJ/mol')
          xparams%nrgshift = xparams%nrgshift/calToJ

        Case('INTRAOFF')   !** turns off intramolecular interaction 
          option_intraoff = .True.

          !** decide which species to turn it off for
          Select Case(chunks(2))
          Case('ALL')  !** all of the concerned species
            intraoff_spc(1) = spc_list(1)
            If (nspc > 1) intraoff_spc(2) = spc_list(2)

          Case Default !** just one, assumed to be a species name
            intraoff_spc(1) = molecules_gettype(chunks(2))
            intraoff_spc(2) = 0

          End Select
          
        Case Default
          Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
              ' Unable to interpret command chunk: ',Trim(chunks(1))
          Stop  
        End Select
      End If

    End Do

    !** Interpret the main EXTERNAL@[]@[] sequence and call initialization
    nchunks = split(fields(1),chunks,'@')
    Select Case (ToUpper(chunks(2)))  !** should be the identifier
    Case ('GULP')
      If (.Not. Present(simcell)) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            ' for GULP initializations, the simcell must be passed'
        Stop
      End If
      If (nchunks < 3) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            ' specifier for external GULP should be '
        Write(0,'(10x,a)') 'External@GULP@[filename] or External@GULP@USEFF'
        Stop  
      End If           
    
      !** Allocate the pointer
      Allocate(xparams%gulp, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'xparams%gulp') 
    
      If (chunks(3) == 'USEFF') Then  !** use internal forcefield params
        If (Present(list)) Then
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Call gulpio_init(xparams%gulp,filename,intraoff_spc, &
              simcell,list,params)
        Else
          Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
              ' paramter list must be passed for USEFF case'
          Stop
        End If

      Else    !** assume 'filename' contains the GULP forcefield specs
        filename = chunks(3)
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Call gulpio_init(xparams%gulp,filename,intraoff_spc,simcell)
      End If
    
    Case ('G98')
      !** Make sure enough information is present
      If (nchunks < 3) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            ' specifier for external Gaussian98 should be '
        Write(0,'(10x,a)') 'EXTERNAL@G98@[ctrl_file_header]@[ctrl_file_tail]'
        Write(0,'(10x,a)') 'where the @[ctrl_file_tail] portion is optional'
        Stop  
      End If           
    
      !** Allocate the pointer
      Allocate(xparams%g98, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'xparams%g98') 

      !** Extract the names for the header and maybe the tail of the control file
      filenames(1) = chunks(3)
      If (nchunks == 4) filenames(2) = chunks(4)

      !** Call the initialization routine
      Call gaussio_init(xparams%g98,.False.,filenames(1:nchunks-2))
        
    Case ('TURBO')
      !** Make sure enough information is present
      If (nchunks /= 3) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            ' specifier for external Turbomole should be '
        Write(0,'(10x,a)') 'EXTERNAL@TURBO@[run_directory]'
        Stop  
      End If           
    
      !** Allocate the pointer
      Allocate(xparams%turbo, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'xparams%turbo') 

      !** Extract the names for the header and maybe the tail of the control file
      filenames(1) = chunks(3)
      If (nchunks == 4) filenames(2) = chunks(4)

      !** Call the initialization routine
      Call turboio_init(xparams%turbo,.False.,filenames(1:1))

    Case Default
      Write(0,'(2a,i5,a,2i3)') __FILE__,":",__LINE__, &
          " Could not interpret external program ID string: ",Trim(chunks(2))
      Stop
          
    End Select

  End Subroutine xternal_init

  !----------------------------------------------------------------------------
  ! Nullify the external interaction pointer set.
  ! Requires:  xparams -- external interaction evaluation parameters
  !----------------------------------------------------------------------------
  Subroutine xternal_nullify(xparams)
    Type(ExternalInfo), Intent(InOut)   :: xparams

    Nullify(xparams%gulp)
    Nullify(xparams%g98)
    Nullify(xparams%turbo)

  End Subroutine xternal_nullify

  !----------------------------------------------------------------------------
  ! Prepares the local information needed for a calculation using input for
  ! intramolecular interactions.
  ! Requires:  xparams -- external interaction evaluation parameters
  !            storage -- intramolecular interaction level-pair for molecule
  !            coords -- coordinates for whole molecule indexed by atm num
  !            spc -- species number, so we can use molecules module
  !            simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Logical Function xternal_intraint(xparams,storage,coords,spc,simcell)
    Type(ExternalInfo), Intent(InOut)            :: xparams
    Type(Store_Level_Pair), Intent(InOut)        :: storage
    Type(VecType), Dimension(:), Intent(In)      :: coords
    Integer, Intent(In)                          :: spc
    Type(SimCell_Params), Intent(In), Optional   :: simcell  

    Integer                 :: a,natoms,nderivs
    Real(kind=RDbl)         :: nrg

    !** set default to False
    xternal_intraint = .False.

    !** Get the number of derivatives needed in output
    nderivs = store_nderivs(storage)
    If (nderivs > 1) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' xternal_intraint not yet setup to return 2nd derivatives'
      Stop
    End If

    !** Resize the temporary storage arrays if not they're not big enough
    natoms = Size(coords)
    Call xternal_inittemp(natoms)

    !** Fill in the temporary arrays
    temp_natoms = natoms
    Do a = 1,natoms
      temp_coords(a) = coords(a)
      temp_atypes(a) = molecules_getatype(spc,a)
      temp_subsets(a,:) = (/spc,0,0/)  !**HACK
      temp_forces(a) = VecType(0.0_RDbl)
    End Do
    Call molecules_getcharges(spc,temp_charges)

    !** Send the information to the external program
    If (Present(simcell)) Then
      xternal_intraint = xternal_send(xparams,nrg,natoms,nderivs,simcell)
    Else
      xternal_intraint = xternal_send(xparams,nrg,natoms,nderivs)
    End If

    !** Place the resultant interactions in the storage
    If (nderivs > 0) Then
      If (store_idbranch(storage) /= 0) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            ' wrong storage type for intramolecular external interactions'
        Stop
      End If

      Call storebase_inc(storage%total,nrg)
      Do a = 1,natoms
        Call storebase_inc(storage%binfo(a),temp_forces(a))
      End Do
    Else
      Call storebase_inc(storage%total,nrg,EXTERNAL_INDEX)
    End If

!LC    Call store_disp(storage,.False.,2,6)

  End Function xternal_intraint

  !-----------------------------------------------------------------------------
  ! Performs a generic (maybe species-species) subsystem interaction evaluation.
  ! Extracts coordinates (+charge) information from the whole system according 
  ! to passed subsets and passes this information to the driver that controls 
  ! the specific external program.  Returns True if evaluation was successful.
  ! Requires:  xparams -- external interaction evaluation parameters
  !            storage -- interaction output storage (spc-spc)
  !            species -- full system coordinates and charges
  !            simcell -- simulation cell information
  !            subset1 -- 1st subset 
  !            subset2 -- 2nd subset, optional
  !-----------------------------------------------------------------------------
  Logical Function xternal_int(xparams,storage,species,simcell,subset1,subset2)
    Type(ExternalInfo), Intent(InOut)                 :: xparams
    Type(Symmetric_Results), Intent(InOut)            :: storage
    Type(AtMolCoords), Dimension(:), Intent(In)       :: species
    Type(SimCell_Params), Intent(In)                  :: simcell  
    Integer, Dimension(:), Intent(In)                 :: subset1
    Integer, Dimension(:), Intent(In), Optional       :: subset2

    Integer                 :: natoms,nderivs,natoms1,natoms2
    Real(kind=RDbl)         :: nrg

    !** set default to False
    xternal_int = .False.

    !** Count the number of atoms to be sent externally
    natoms1 = config_nsubsetatoms(species,subset1)
    If (Present(subset2)) Then
      natoms2 = config_nsubsetatoms(species,subset2)
    End If

    !** Make sure the interacting species are present
    If (Present(subset2)) Then    
      If ((natoms1 == 0).Or.(natoms2 == 0)) Then
        xternal_int = .True.    !** no interactions, we're done
        Return   
      End If
      natoms = natoms1 + natoms2
    Else
      If (natoms1 == 0) Then
        xternal_int = .True.    !** no interactions, we're done
        Return   
      End If
      natoms = natoms1 
    End If

    !** Get the number of derivatives needed in output
    nderivs = storesym_nderivs(storage)
    If (nderivs > 1) Then
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' xternal_send not yet setup to return 2nd derivatives'
      Stop
    End If

    !** Resize the temporary storage arrays if not they're not big enough
    Call xternal_inittemp(natoms)

    !** Fill in the temporary arrays
    temp_natoms = 0
    Call xternal_copysubset(species,subset1)
    If (Present(subset2)) Then
      Call xternal_copysubset(species,subset2)      
    End If
    natoms = temp_natoms

    !** Send the information to the external program
    xternal_int = xternal_send(xparams,nrg,natoms,nderivs,simcell)

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'calculated energy is: ',nrg

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call storesym_display(storage,.False.,2,6)

    !** Place the resultant interactions in the storage
    If (Present(subset2)) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      !** This case is complicated because the output likely
      !** includes contributions from spc1 and spc2 with themselves
      !** How to place in storage?  
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' Not ready to place results of multiple subset evaluations'
      Stop

    Else
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      !** this is easier, get appropriate level pair set from
      !** storesym_ptrs, then insert energy
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' Not ready to place subset1-subset1 evaluations'
      Stop
    End If

  End Function xternal_int

  !----------------------------------------------------------------------------
  ! Perform a generic evaluation that is not specifically intramolecular.
  ! This routine takes coordinates and charges from the local temporary
  ! arrays and passes them to the appropriate external io module along with
  ! pre-determined information that may be needed to conduct the calculation.  
  ! It then returns interactions such as an energies and optionally forces.  
  ! Returned ENERGY AND LENGTH UNITS are kcal/mol and Angstroms
  ! Requires:  xparams -- external interaction evaluation parameters
  !            nrg -- resultant subsystem interaction energy
  !            natoms -- number of atoms in the temporary arrays
  !            nderivs -- number of derivatives expected in the results
  !            simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Logical Function xternal_send(xparams,nrg,natoms,nderivs,simcell)
    Type(ExternalInfo), Intent(InOut)          :: xparams
    Real(kind=RDbl), Intent(Out)               :: nrg
    Integer, Intent(In)                        :: natoms,nderivs
    Type(SimCell_Params), Intent(In), Optional :: simcell  

    Real(kind=RDbl)         :: totalq

    !** set default
    xternal_send = .False.

    !** Call the correct module for handling external programs
    If (Associated(xparams%gulp)) Then
      If (nderivs > 0) Then
        If (Present(simcell)) Then
          xternal_send = gulpio_snglpt(xparams%gulp,.True., &
              simcell,natoms,temp_coords(1:natoms), &
              temp_atypes(1:natoms),temp_charges(1:natoms), &
              temp_subsets(1:natoms,:),nrg,temp_forces(1:natoms))
        Else
          xternal_send = gulpio_snglpt(xparams%gulp,.True., &
              xparams%gulp%simcell,natoms,temp_coords(1:natoms), &
              temp_atypes(1:natoms),temp_charges(1:natoms), &
              temp_subsets(1:natoms,:),nrg,temp_forces(1:natoms))
        End If
      Else
        If (Present(simcell)) Then
          xternal_send = gulpio_snglpt(xparams%gulp,.True., &
              simcell,natoms,temp_coords(1:natoms), &
              temp_atypes(1:natoms),temp_charges(1:natoms), &
              temp_subsets(1:natoms,:),nrg)
        Else
          xternal_send = gulpio_snglpt(xparams%gulp,.True., &
              xparams%gulp%simcell,natoms,temp_coords(1:natoms), &
              temp_atypes(1:natoms),temp_charges(1:natoms), &
              temp_subsets(1:natoms,:),nrg)
        End If
      End If

    Else If (Associated(xparams%g98)) Then
      If (Present(simcell)) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            ' Cannot use simulation cell in a Gaussian98 calculation'
        Stop
      End If

      !** Get the total charge
      totalq = Sum(temp_charges(1:natoms))

      If (nderivs > 0) Then
        xternal_send = gaussio_snglpt(xparams%g98,totalq,natoms, &
            temp_coords(1:natoms),temp_atypes(1:natoms), &
            temp_subsets(1:natoms,:),nrg,temp_charges(1:natoms), &
            temp_forces(1:natoms))
      Else
        xternal_send = gaussio_snglpt(xparams%g98,totalq,natoms, &
            temp_coords(1:natoms),temp_atypes(1:natoms), &
            temp_subsets(1:natoms,:),nrg,temp_charges(1:natoms))
      End If

    Else If (Associated(xparams%turbo)) Then
      If (Present(simcell)) Then
        Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
            ' Cannot use simulation cell in a Turbomole calculation'
        Stop
      End If

      !** Get the total charge
      totalq = Sum(temp_charges(1:natoms))

      If (nderivs > 0) Then
        xternal_send = turboio_snglpt(xparams%turbo,totalq,natoms, &
            temp_coords(1:natoms),temp_atypes(1:natoms), &
            temp_subsets(1:natoms,:),nrg,temp_charges(1:natoms), &
            temp_forces(1:natoms))
      Else
        xternal_send = turboio_snglpt(xparams%turbo,totalq,natoms, &
            temp_coords(1:natoms),temp_atypes(1:natoms), &
            temp_subsets(1:natoms,:),nrg,temp_charges(1:natoms))
      End If

    Else 
      Write(0,'(2a,i6,a)') __FILE__,":",__LINE__, &
          ' No recognized external evaluation pointer associated'
      Stop
    End If

    !** Shift the resulting energy if previously specified (operating in kcal/mol)
    If (xparams%shift) Then
      nrg = nrg + xparams%nrgshift
    End If

  End Function xternal_send

  !----------------------------------------------------------------------------
  ! Allocates or resizes the temporary arrays
  ! Requires:  size -- desired size of temporary arrays
  !----------------------------------------------------------------------------
  Subroutine xternal_inittemp(size)
    Integer, Intent(In)     :: size

    Integer                 :: error
    Logical                 :: resize
    Logical, Save           :: firsttime = .True.

    !** Decide if the arrays should be allocated
    resize = .False.
    If (temp_natoms < size) resize = .True.

    !** Deallocate the arrays first, if this isn't the first time
    If (firsttime) Then
      firsttime = .False.
      resize = .True.
    Else If (resize) Then
      Deallocate(temp_coords,STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'temp_coords')
      Deallocate(temp_atypes,STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'temp_atypes')
      Deallocate(temp_charges,STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'temp_charges')
      Deallocate(temp_subsets,STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'temp_subsets')
      Deallocate(temp_forces,STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__,'temp_forces')
    End If

    !** Allocate the temporary storage arrays if they're too small
    If (resize) Then
      Allocate(temp_coords(size),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'temp_coords')
      Allocate(temp_atypes(size),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'temp_atypes')
      Allocate(temp_charges(size),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'temp_charges')
      Allocate(temp_subsets(size,3),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'temp_subsets')
      Allocate(temp_forces(size),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'temp_forces')
    End If

  End Subroutine xternal_inittemp

  !----------------------------------------------------------------------------
  ! Copy the coordinates, atom types, charges and subsets into the temporary
  ! storage arrays.
  ! Requires:  species -- full system coordinates and charges
  !            subset -- subset to extract from system
  !----------------------------------------------------------------------------
  Subroutine xternal_copysubset(species,subset)
    Type(AtMolCoords), Dimension(:), Intent(In)       :: species
    Integer, Dimension(:), Intent(In)                 :: subset

    Integer       :: spc,mol,atm,natoms,nmoles

    Do spc = 1,Size(species)
      If ((subset(1) /= 0).And.(subset(1) /= spc)) Cycle

      natoms = molecules_getnatoms(spc)
!      natoms = 10
      nmoles = config_getnmoles(species,spc)

      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'subset: ',subset,'      nmoles,natoms:',nmoles,natoms

      Do mol = 1,nmoles
        If ((subset(2) /= 0).And.(subset(2) /= mol)) Cycle

        Do atm = 1,natoms
          If ((subset(3) /= 0).And.(subset(3) /= atm)) Cycle
          temp_natoms = temp_natoms + 1
          temp_coords(temp_natoms) = config_getr(species,(/spc,mol,atm/))
          temp_atypes(temp_natoms) = config_getatype(species,spc,mol,atm)
          temp_charges(temp_natoms) = molecules_getcharge(spc,atm)
          temp_subsets(temp_natoms,:) = (/spc,mol,atm/)
        End Do
      End Do
    End Do

  End Subroutine xternal_copysubset

  !----------------------------------------------------------------------------
  ! Display the parameters
  ! Requires:  xparams -- external interaction evaluation parameters
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional output unit number, default is 6
  !----------------------------------------------------------------------------
  Subroutine xternal_display(xparams,indent,unitno)
    Type(ExternalInfo), Intent(In)    :: xparams
    Integer, Intent(In)               :: indent
    Integer, Intent(In), Optional     :: unitno

    Integer                     :: unit
    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)
    unit = 6
    If (Present(unitno)) unit = unitno

    Write(unit,'(2a)') blank,'External evaluation will be performed'

    If (Associated(xparams%gulp)) Then
      Call gulpio_display(xparams%gulp,indent+2,unit)
    End If
    If (Associated(xparams%g98)) Then
      Call gaussio_display(xparams%g98,indent+2,unit)
    End If
    If (Associated(xparams%turbo)) Then
      Call turboio_display(xparams%turbo,indent+2,unit)
    End If

  End Subroutine xternal_display

End Module xternal
