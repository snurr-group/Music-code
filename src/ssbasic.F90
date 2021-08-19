!------------------------------------------------------------------------------
! This module deals with simple interactions, both coulombic and non-coulombic
! Simple means pairwise-interactions where a cut-off is applied.
! The main purpose of this module is to convert interactions of multiple atom 
! structures to atom-atom interactions.  In doing this, it also prevents double
! counting when a flag is specified.
!
! A hierarchy of calls are provided for evaluating multi-atom interactions.  
!   ssbasic_ssint --------> ssbasic_msint -------> ssbasic_mmint
!        |                       |                      |
!  spc-spc interactions          |                      |
!                       molec-spc interactions          |
!                                               molec-molec interactions
! Or, alternatively, for just a single atom.
!                           ssbasic_asint -------> ssbasic_amint
!  NOTE: this sequence           |                      |
!  uses only EnergyPlus          |                      |
!                        atom-spc interactions          |
!                                               atom-molec interactions
!------------------------------------------------------------------------------

Module ssbasic

  Use defaults, Only: RDbl, strLen, lstrLen, xlstrLen, MAX_ATOMTYPES, scalef, &
      d_aa_file, d_ctrl_file, dashedline, SIZE_OF_HOT, zero
  Use file, Only: file_open, file_gettype
  Use utils, Only: toupper, split, toreal, filesrchstr, allocerrdisplay, &
      deallocerrdisplay, findint, checkandstop, int2str, swap
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getdistsq, vector_display
  Use atom, Only: atom_getntypes, atom_getname
  Use molecules, Only: molecules_getnatomtypes, molecules_getnatoms,  &
      molecules_getatype
  Use config, Only: config_getnmoles, config_getr, config_getMolecCTR, &
      AtMolCoords, config_getmolr, config_getmolq, config_getq, config_subset2xyz
  Use simcell, Only: SimCell_Params, simcell_minimage, simcell_minwidth
  Use nbrlist, Only: Nbrlist_Params, nbrlist_getlist
  Use pairmodel, Only: Pairmodel_Params, pairmodel_init, pairmodel_mmint, &
      pairmodel_mmint2, pairmodel_mmint3, pairmodel_mmint4, pairmodel_mmint5, &
      pairmodel_msint5, pairmodel_chkcutoffs, pairmodel_nullset, &
      pairmodel_countset, pairmodel_chkarray, pairmodel_makesym, &
      pairmodel_display, pairmodel_clean, pairmodel_getpotparameters
  Use smooth, Only: Smooth_Method,smooth_init,smooth_perform,smooth_perform1, &
      smooth_display,smooth_clean
  Use store, Only: Store_Level_Pair,store_display,store_idbranch,store_disp, &
      store_sum
  Use storebase, Only: EnergyPlus,storebase_nderivs,storebase_disp,storebase_init
  Use storesym, Only: Symmetric_Results, storesym_ptrs, storesym_display  

  Implicit None
  Save

  Private
  Public :: ssbasic_idstring, SSBasicParams, ssbasic_init, ssbasic_ssint, &
      ssbasic_asint, ssbasic_msint, ssbasic_mmint, ssbasic_display, &
      ssbasic_clean, ssbasic_getpotparameters

  Character(len=strLen), Parameter       :: ssbasic_idstring = 'BASIC'

  Type SSBasicParams
    Logical                                         :: smooth,comcut
    Logical                                         :: neighbor,fast
    Integer                                         :: spc1,spc2
    Real(kind = RDbl)                               :: smooth_rad,smooth_rad2
    Real(kind = RDbl)                               :: cutrad,cutrad2
    Character(len=lstrLen)                          :: line
    Character(len=strLen)                           :: potmodel
    Type(Smooth_Method), Pointer                    :: smooth_params
    Type(Pairmodel_Params), Dimension(:,:), Pointer :: pparams
    Type(Nbrlist_Params) , Pointer :: nlist
  End Type SSBasicParams

  Logical :: no_charges=.True.

Contains
  !----------------------------------------------------------------------------
  ! Handles the initialization for a species pair
  ! Requires:  params -- SSBasic parameters
  !            spc1 -- 1st species number
  !            spc2 -- 2nd species number
  !            simcell -- simulation cell information
  !            line -- line of commands to interpret
  !----------------------------------------------------------------------------
  Subroutine ssbasic_init(params,spc1,spc2,simcell,line)
    Type(SSBasicParams), Intent(InOut)   :: params
    Integer, Intent(In)                  :: spc1,spc2
    Type(SimCell_Params), Intent(In)     :: simcell
    Character(*), Intent(In)             :: line

    Integer                              :: nfields,n,i,unit,lineno
    Integer                              :: a1,a2,atype1,atype2
    Integer                              :: natom1,natom2,natypes,error
    Logical                              :: endit, has_charges
    Real(kind=RDbl)                      :: maxcutoff
    Character(len=xlstrLen)              :: typesline
    Character(len=lstrLen)               :: filename,srchstr,junk
    Character(len=strLen), Dimension(10) :: fields,chunks
    Integer, Dimension(MAX_ATOMTYPES)    :: atypes1,atypes2

    params%line = line
    params%spc1 = spc1
    params%spc2 = spc2
    nfields = split(line,fields)
    Nullify(params%nlist)

    !** Initialize the fields that govern the interaction between two
    !** molecules of a particular species-species pair 

    !** Is this a fast interaction?
    params%fast = .False.

    !** Should we use COM truncation?
    params%comcut = .False.

    !** Should we smooth the cutoff?
    params%smooth = .False.
    Nullify(params%smooth_params)

    !** Should we make a neighbor list?
    params%neighbor = .False.

    !** shift through fields 5->end and pass unknown strings to pairmodel
    !** it is assumed that these other strings are potential types
    typesline = ''
    Do i = 5,nfields
      endit = .False.
      n = split(fields(i),chunks,'@')

      Select Case(Toupper(chunks(1)))
      Case('FAST')
        params%fast = .TRUE.

      Case('NEIGHBOR')
        params%neighbor = .TRUE.
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
             ' Sorry, neighbor list not yet available'
        endit = .True.

      Case('COMTRUNC')
        params%comcut = .TRUE.
        params%cutrad = toreal(chunks(2))
        params%cutrad2 = params%cutrad**2
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
             ' WARNING: using approximation for COM during cut-off calculation'

      Case('COMCUT')
        params%comcut = .TRUE.
        params%cutrad = toreal(chunks(2))
        params%cutrad2 = params%cutrad**2
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
             ' WARNING: using approximation for COM during cut-off calculation'

        maxcutoff = simcell_minwidth(simcell)/2.0_RDbl
        If (params%cutrad > maxcutoff) Then
          Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
              ' Problem with mol-mol cut-off, must be less than half'
          Write(0,'(a,f8.4)') '  the minimum width of the simulation cell:', &
              maxcutoff
          Stop
        End If

      Case('SMOOTH')
        params%smooth = .TRUE.
        params%smooth_rad = toreal(chunks(2))
        params%smooth_rad2 = params%smooth_rad**2

      Case Default
        Write(typesline,'(a,2x,a)') Trim(typesline),Trim(fields(i))
      End Select

      !** Stop if necessary
      If (endit) Then
        Write(0,'(2a)') ' please turn option off, the command was: ',&
            Trim(fields(i))
        Write(0,'(2a)') 'command line: ',Trim(line)
        Stop        
      End If
    End Do
    IF(.NOT. params%fast) THEN
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) 'Simulation will not run without the fast flag in the sorb_sorb_file!'
        STOP
    END IF

    !** Initialize the smoothing method if necessary
    If (params%smooth) Then
      !** make sure there is a mol-mol cut-off if there is smoothing
      If (.Not. params%comcut) Then
        Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
            " A mol-mol cut-off must be specified if smoothing is specified"
        Stop    
      End If

      Allocate(params%smooth_params,stat=error)    
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)   
      Call smooth_init(params%smooth_params,'FNMULT',params%smooth_rad, &
          params%cutrad)
    End If

    !** Get the filename that contains the atom-atom interactions
    Call file_gettype(d_aa_file,filename,unit)
    If (filename == '') Then       
      !** Default, Read info from control file, not used much
      Call file_gettype(d_ctrl_file,filename,unit)
      unit = file_open(filename,110)
      
      srchstr = 'Atom-Atom Interactions'  !**header line
      lineno = filesrchstr(unit, srchstr, junk)
      If (lineno == 0) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
             'Could not find search string ',trim(srchstr)
        Stop
      End If 
    End If

    !** Allocate space for the atomtype-atomtype pairwise interactions
    natypes = atom_getntypes()
    Allocate(params%pparams(natypes,natypes),stat=error)    
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'pairparams')   

    !** Nullify the pointers
    Do a1 = 1,natypes
      Do a2 = 1,natypes
        Call pairmodel_nullset(params%pparams(a1,a2))
      End Do
    End Do

    !** Initialize the individual atomtype-atomtype pairwise interactions
    !** only those atom types present in the two species are initialized
    atypes1 = 0
    atypes2 = 0
    natom1 = molecules_getnatomtypes(spc1,atypes1)
    natom2 = molecules_getnatomtypes(spc2,atypes2) 
    Do atype1 = 1,natypes
      If (findint(atypes1(1:natom1),atype1) == 0) Cycle
      Do atype2 = 1,natypes
        If (findint(atypes2(1:natom2),atype2) == 0) Cycle
        Call pairmodel_init(params%pparams,atype1,atype2,typesline,filename,&
            has_charges)
        If (has_charges) no_charges=.False.
        Call pairmodel_chkcutoffs(params%pparams(atype1,atype2),simcell)
      End Do
    End Do

    !** Count initialized pair parameter types
    Do a1 = 1,natypes
      Do a2 = 1,natypes
        Call pairmodel_countset(params%pparams(a1,a2))
      End Do
    End Do

    !** Check for unreasonable specifications
    Do atype1 = 1,natypes
      If (findint(atypes1,atype1) == 0) Cycle
      Do atype2 = 1,natypes
        If (findint(atypes2,atype2) == 0) Cycle

        !** make sure there is a mol-mol cut-off for coulombic interactions
        If (Associated(params%pparams(atype1,atype2)%coul)) Then
          If (.Not. params%comcut) Then
            Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
                " We will use atome based cut off"
!            Stop
          End If
        End If
      End Do
    End Do

    !** Symmeterize the array, then check it
    Call pairmodel_makesym(params%pparams)
    Call pairmodel_chkarray(params%pparams)

  End Subroutine ssbasic_init

  !----------------------------------------------------------------------------
  ! Handles a single SPECIES-SPECIES interaction.  
  ! Requires:  ffparams -- spc-spc pairwise forcefield parameters
  !            grp1out -- spc1-spc2 forcefield output (mi-associated)
  !            grp2out -- spc2-spc1 forcefield output (im-associated)
  !            spc1 -- 1st species number
  !            spc2 -- 2nd species number
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !            allint -- True => evaluate ALL atom pairs, otherwise no doubles
  !----------------------------------------------------------------------------
  Logical Function ssbasic_ssint(ffparams,ffout,spc1,spc2,species, &
      simcell,fast,allint)
    Type(SSBasicParams), Intent(In)              :: ffparams
    Type(Symmetric_Results), Intent(InOut)       :: ffout
    Integer, Intent(In)                          :: spc1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell  
    Logical, Intent(In)                          :: fast,allint

    Integer                   :: m

    !** Set default return value
    ssbasic_ssint = .True.

    !** Check if the interaction "speeds" are complimentary
    If (fast .Neqv. (ffparams%fast)) Return

    Do m = 1,config_getnmoles(species,spc1)
      ssbasic_ssint = ssbasic_msint(ffparams,ffout,spc1,m,spc2, &
          species,simcell,fast,allint)
      If (.Not. ssbasic_ssint) Return
    End Do

  End Function ssbasic_ssint

#ifndef FAST_ASINT
  !----------------------------------------------------------------------------
  ! Handles a single ATOM-SPECIES interaction.  Skips interactions with the
  ! same molecule.  
  ! ********* OLD VERSION, NEAT CODE, BUT TOO SLOW **************
  !                  see the other faster version below 
  ! NOTE: This only operates on a single EnergyPlus structure, and not the 
  ! full storage structure because Shaji wanted it that way (and I'm lazy).
  ! Requires:  ffparams -- spc-spc pairwise forcefield parameters
  !            results -- single EnergyPlus structure, WILL NOT BE ZEROED
  !            spc1 -- 1st species number
  !            mol1 -- 1st molecule number
  !            atm1 -- 1st atom number
  !            spc2 -- 2nd species number
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !----------------------------------------------------------------------------
  Logical Function ssbasic_asint(ffparams,results,spc1,mol1,atm1,spc2,species, &
      simcell,fast,hot)
    Type(SSBasicParams), Intent(In)                   :: ffparams
    Type(EnergyPlus), Intent(InOut)                   :: results
    Integer, Intent(In)                               :: spc1,mol1,atm1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In)       :: species
    Type(SimCell_Params), Intent(In)                  :: simcell  
    Logical, Intent(In)                               :: fast
    Real(Kind=RDbl), Dimension(SIZE_OF_HOT), Optional :: hot

    Integer                                 :: m,natoms,skip
    Real(kind=RDbl), Dimension(SIZE_OF_HOT) :: hotmol

    !** Set default return value
    ssbasic_asint = .True.

    !** Check if the interaction "speeds" are complimentary
    If (fast .Neqv. (ffparams%fast)) Return

    !** Avoid intramolecular interaction
    skip = 0
    If (spc1 == spc2) skip = mol1

    natoms = molecules_getnatoms(spc2)
    If (Present(hot)) hot=zero

    Do m = 1,config_getnmoles(species,spc2)
      If (m == skip) Cycle
      If (.Not.Present(hot)) Then

        ssbasic_asint = ssbasic_amint(ffparams,results,spc1,mol1,atm1, &
            spc2,m,natoms,species,simcell,fast)
      Else

        ssbasic_asint = ssbasic_amint(ffparams,results,spc1,mol1,atm1, &
            spc2,m,natoms,species,simcell,fast,hotmol)
        hot = hot + hotmol
      End If
      If (.Not. ssbasic_asint) Return
    End Do

  End Function ssbasic_asint
#endif
#ifdef FAST_ASINT
  !----------------------------------------------------------------------------
  ! Handles a single ATOM-SPECIES interaction.  Skips interactions with the
  ! same molecule.  
  ! *** NEW VERSION , SIMILAR TO MUSIC-2-2, FASTER ****************************
  ! *needs more functions (hot, ctcdist etc), but as much as 3-6 times faster**
  ! NOTE: This only operates on a single EnergyPlus structure, and not the 
  ! full storage structure because Shaji wanted it that way (and I'm lazy).
  ! Requires:  ffparams -- spc-spc pairwise forcefield parameters
  !            results -- single EnergyPlus structure, WILL NOT BE ZEROED
  !            spc1 -- 1st species number
  !            mol1 -- 1st molecule number
  !            atm1 -- 1st atom number
  !            spc2 -- 2nd species number
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !----------------------------------------------------------------------------
  Logical Function ssbasic_asint(ffparams,results,spc1,mol1,atm1,spc2,&
      species, simcell, fast, hot )
    Type(SSBasicParams), Intent(In)              :: ffparams
    Type(EnergyPlus), Intent(InOut)              :: results
    Integer, Intent(In)                          :: spc1,mol1,atm1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell  
    Logical, Intent(In)                          :: fast
    Real(Kind=RDbl), Dimension(SIZE_OF_HOT), Optional :: hot

    Integer                   :: low, high, nmoles2
    If (Present(hot)) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif
    !** Set default return value
    ssbasic_asint = .True.
    nmoles2 = config_getnmoles(species,spc2)

    !** Check if the interaction "speeds" are complimentary
    If (fast .Neqv. (ffparams%fast)) Return
    If (spc1 == spc2) Then
      !** We don't want interactions of molec1 with itself.  So, we split
      !** it into two parts 1.) 1 to molec1-1 and 2.) molec1+1 to nmoles2
      low  = 1
      high = mol1 - 1
      If (low <= high) Then
        ssbasic_asint = ssbasic_atomgetas(ffparams, results, species, &
            simcell, spc1, &
            mol1, atm1, spc2, low, high)
        If ( .Not. ssbasic_asint) Return  
      End If
      low  = mol1 + 1
      high = nmoles2
      If (low <= high) Then

        ssbasic_asint = ssbasic_atomgetas(ffparams, results, species, &
            simcell, spc1, &
            mol1, atm1, spc2, low, high)
        If ( .Not. ssbasic_asint) Return  
      Endif
    Else
      low = 1
      high = nmoles2
      ssbasic_asint = ssbasic_atomgetas(ffparams, results, species, &
          simcell, spc1, &
          mol1, atm1, spc2, low, high)
      If ( .Not. ssbasic_asint) Return  

    Endif
  End Function ssbasic_asint
#endif

  !----------------------------------------------------------------------
  ! Handles a single ATOM-SPECIES interaction.  It decomposes the
  ! interactions into VECTOR-VECOTR interactions and passes the information
  ! Can not deal with charges/ctcdits etc, ( needs more coding for that )
  ! to the pairmodel module in one call.  
  ! Requires:  ffparams -- spc-spc pairwise forcefield parameters
  !            results -- atm1-spc2 Total energy structure output
  !            spc1 -- 1st species number
  !            mol1 -- 1st molecule number
  !            atm1 -- 1st atom number
  !            spc2 -- 2nd species number
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !----------------------------------------------------------------------
  Logical Function ssbasic_atomgetas(ffparams,results,species, simcell, &
      spc1,mol1,atm1, spc2,low,high)
    Type(SSBasicParams), Intent(In)              :: ffparams
    Type(EnergyPlus), Intent(InOut)              :: results
    Integer, Intent(In)                          :: spc1,mol1,atm1,low,high
    Integer, Intent(In)                          :: spc2
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell  

    Integer, Dimension(1)      :: atypes1
    Integer, Dimension(species(spc2)%natoms) :: atypes2
    Type(VecType), Dimension(1)::  atVecs1
    Type(VecType), Dimension(species(spc2)%natoms*(high-low+1))  :: atVecs2

    Integer  :: nvecs, natoms1, natoms2, natypes2, nmoles2, index, m, a

    Logical, Save                        :: firsttime = .True.
    Type(EnergyPlus), Save               :: opposite1,opposite0


    !** Allocate the dummy opposite interaction structure
    If (firsttime) Then
      firsttime=.False.
      Call storebase_init(opposite1,1)
      Call storebase_init(opposite0,0)
    End If

    natoms1=1
    atVecs1(1) = species(spc1)%coords(atm1,mol1)%r
    atypes1 = molecules_getatype(spc1,atm1)

    natoms2= molecules_getnatoms(spc2)
    natypes2 = molecules_getnatomtypes(spc2, atypes2(1:natoms2), .True.)
    nmoles2 = high - low + 1

    !** Flatten the species coords down to a 1D array
    !** atVecs2(1) = coords(1,1)%r
    !** atVecs2(2) = coords(2,1)%r
    !** atVecs2(n) = coords(a,m)%r,    n = natoms*(m-1)+a 
    !** Also, a = Mod(n,natoms) but if (a == 0), a = natoms
    !**       m = Ceiling(Real(n)/Real(natoms))

    nvecs = natoms2*nmoles2
    !** The loop below runs faster than this reshape statement
    !**   atVecs2(1:nvecs) = &
    !**   Reshape(species(spc2)%coords(1:natoms2,mbeg:mend)%r,(/nvecs/))
    index=0
    Do m=low,high
      Do a=1,natoms2
        index=index+1
        atVecs2(index)=species(spc2)%coords(a,m)%r
      End Do
    End Do
    !** Get the interactions
    ssbasic_atomgetas = pairmodel_msint5(ffparams%pparams,results,opposite0, &
        atypes1,atVecs1,atypes2,atVecs2,simcell) 

  End Function ssbasic_atomgetas

  !----------------------------------------------------------------------------
  ! Handles a single MOLECULE-SPECIES interaction
  ! Requires:  ffparams -- spc-spc pairwise forcefield parameters
  !            grp1out -- mol1-spc2 forcefield output (im-associated)
  !            grp2out -- spc2-mol1 forcefield output (mi-associated)
  !            spc1 -- 1st species number
  !            mol1 -- 1st molecule number
  !            spc2 -- 2nd species number
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !            allint -- True => evaluate ALL atom pairs, otherwise no doubles
  !----------------------------------------------------------------------------
  Logical Function ssbasic_msint(ffparams,ffout,spc1,mol1,spc2, &
      species,simcell,fast,allint)
    Type(SSBasicParams), Intent(In)              :: ffparams
    Type(Symmetric_Results), Intent(InOut)       :: ffout
    Integer, Intent(In)                          :: spc1,mol1,spc2
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell  
    Logical, Intent(In)                          :: fast,allint

    Integer         :: low,high,m,nmoles2,natoms1,natoms2,list_size,l
    Integer, Dimension(:), pointer :: list

    !** Set default return value
    ssbasic_msint = .True.

    !** Check if the interaction "speeds" are complimentary
    If (fast .Neqv. (ffparams%fast)) Return

    !** Get number of molecules in species2, return if zero
    nmoles2 = config_getnmoles(species,spc2)
    If (nmoles2 == 0) Return

    natoms1 = molecules_getnatoms(spc1)
    natoms2 = molecules_getnatoms(spc2)

    If (Associated(ffparams%nlist)) Then
      !***
      !***  Use Neghborlist

      ! nlist = neighborlist object
      ! list  = array containing number of neighbors 
      ! list_size =number of neighbors
      Call nbrlist_getlist(ffparams%nlist, spc1, mol1, spc2, list, list_size)

      ! case for MD, allint=.True.
      If (spc1 == spc2 .AND. allint) Then
!        Write(*,*) spc1, mol1, spc2, list_size
!        Write(*,*) list
        Do l = 1, list_size 
          m=list(l)
          If (m<=mol1) Cycle
            ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc1,mol1,natoms1, &
                spc2,m,natoms2,species,simcell,fast)
          If (.Not. ssbasic_msint) Return
        End Do

        ! this is not MD , for now just stop
      Else If (spc1 == spc2) Then
        Write(*,*) "This nbrlist is only for MD"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Else
        !** calculate full interaction because spc1 /= spc2
        Do l=1,list_size
          m=list(l)
          If (spc1 > spc2) Then  !** flip order
            ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc2,m,natoms2, &
                spc1,mol1,natoms1,species,simcell,fast)
          Else
            ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc1,mol1,natoms1, &
                spc2,m,natoms2,species,simcell,fast)
          End If
          If (.Not. ssbasic_msint) Return
        End Do
      End If

    Else
      !***
      !***  Without Neghborlist

      !** If the allint flag is set to .TRUE. then we can safely
      !** assume that interactions of mol1 with molecules [1 - (mol1-1)]
      !** of the same species type have already been calculated  
      If (spc1 == spc2 .AND. allint) Then
        low = mol1 + 1
        high = nmoles2
        If (low > high) Return

        Do m = low,high
          If (spc1 > spc2) Then  !** flip order
            ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc2,m,natoms2, &
                spc1,mol1,natoms1,species,simcell,fast)
          Else
            ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc1,mol1,natoms1, &
                spc2,m,natoms2,species,simcell,fast)
          End If
          If (.Not. ssbasic_msint) Return
        End Do

      Else If (spc1 == spc2) Then
        !** Don't want mol1-mol1 interactions. Split into 1st:  1 -> mol1-1
        low  = 1
        high = mol1 - 1
        If (low <= high) Then
          Do m = low,high
            If (spc1 > spc2) Then  !** flip order
              ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc2,m,natoms2, &
                  spc1,mol1,natoms1,species,simcell,fast)
            Else
              ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc1,mol1,natoms1, &
                  spc2,m,natoms2,species,simcell,fast)
            End If
            If (.Not. ssbasic_msint) Return
          End Do
        End If

        !** then into second part, 2nd:  mol1+1 -> nmoles2
        low  = mol1 + 1
        high = nmoles2
        If (low <= high) Then
          Do m = low,high
            If (spc1 > spc2) Then  !** flip order
              ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc2,m,natoms2, &
                  spc1,mol1,natoms1,species,simcell,fast)
            Else
              ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc1,mol1,natoms1, &
                  spc2,m,natoms2,species,simcell,fast)
            End If
            If (.Not. ssbasic_msint) Return
          End Do
        End If

      Else
        !** calculate full interaction because spc1 /= spc2
        low = 1
        high = nmoles2
        Do m = low,high
          If (spc1 > spc2) Then  !** flip order
            ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc2,m,natoms2, &
                spc1,mol1,natoms1,species,simcell,fast)
          Else
            ssbasic_msint = ssbasic_mmint(ffparams,ffout,spc1,mol1,natoms1, &
                spc2,m,natoms2,species,simcell,fast)
          End If
          If (.Not. ssbasic_msint) Return
        End Do
      End If


    Endif

  End Function ssbasic_msint

  !----------------------------------------------------------------------
  ! Handles a single MOLECULE-MOLECULE interaction.  It decomposes the
  ! interactions into atom-atom interactions and passes the information
  ! to the pairmodel module in one call.  The passed species numbers must
  ! obey the spc1 < spc2 convention.  Additionaly, this routine applies 
  ! one convention internally to the molecule pair numbering.  If the
  ! two species are the same, then mol1 < mol2.  These conventions are 
  ! necessary to maintain consistency in the interactions storage.
  ! Requires:  ffparams -- spc-spc pairwise forcefield parameters
  !            ffout -- symmetric results structure
  !            spc1 -- 1st species number
  !            mol1 -- 1st molecule number
  !            natoms1 -- number of atoms in 1st species
  !            spc2 -- 2nd species number
  !            mol2 -- 2st molecule number
  !            natoms2 -- number of atoms in 2nd species
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !----------------------------------------------------------------------
  Logical Function ssbasic_mmint(ffparams,ffout,spc1,mol1,natoms1, &
      spc2,mol2,natoms2,species,simcell,fast)
    Type(SSBasicParams), Intent(In)              :: ffparams
    Type(Symmetric_Results), Intent(InOut)       :: ffout
    Integer, Intent(In)                          :: spc1,mol1,natoms1
    Integer, Intent(In)                          :: spc2,mol2,natoms2
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell  
    Logical, Intent(In)                          :: fast

    Integer                              :: unit,nderivs,branch,m1,m2
    Integer                              :: natypes1,natypes2
    Real(kind=RDbl)                      :: ctcdist2,ctcdist
    Type(VecType)                        :: com1,com2,sepvec
    Character(len=lstrLen)               :: string1,string2
    Type(VecType), Dimension(natoms1)    :: atVecs1
    Type(VecType), Dimension(natoms2)    :: atVecs2
    Integer, Dimension(natoms1)          :: atypes1
    Integer, Dimension(natoms2)          :: atypes2
    Real(kind=RDbl), Dimension(natoms1)  :: charges1
    Real(kind=RDbl), Dimension(natoms2)  :: charges2
    Integer, Dimension(2)                :: indx
    Type(EnergyPlus), Pointer            :: nrgptr
    Type(Store_Level_Pair), Dimension(:), Pointer      :: abptr,baptr

    !** Set default
    ssbasic_mmint = .True.

    !** Get the storage pointers, accounts for m1<m2 convention internally
    !    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    !    Write(*,*) 'Into storesym_ptrs: ',spc1,mol1,'     ',spc2,mol2
    Call storesym_ptrs(ffout,(/spc1,mol1,0/),(/spc2,mol2,0/), &
        branch,indx,nrgptr,abptr,baptr)
    !    Call storesym_display(ffout,.False.,2,6)

    !** Copy the molecule and species numbers
    m1 = mol1
    m2 = mol2

    !** If it's same species interacting with itself, apply m1<m2 convention 
    If ((spc1 == spc2).And.(m1 > m2)) Call swap(m1,m2)

#ifdef DEBUG
    If (spc1 > spc2) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Stop
    End If
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'After ssbasic convention application: ',spc1,m1,'     ',spc2,m2
#endif

    !** Check the center-to-center cut-off if necessary
    If (ffparams%comcut) Then
      com1 = config_getMolecCTR(species,spc1,m1)
      com2 = config_getMolecCTR(species,spc2,m2)
      Call simcell_minimage(simcell,com1,com2,sepvec,ctcdist2)

      !      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      !      Write(*,*) m1,m2,Sqrt(ctcdist2)

      !** Determine if the molecule-molecule evaluation should not be done
      ssbasic_mmint=.True.
      If (ctcdist2 > ffparams%cutrad2) Return
    End If

    !** Extract the atomic coordinates for each molecule
    Call config_getmolr(species,spc1,m1,atVecs1)
    Call config_getmolr(species,spc2,m2,atVecs2)

    !** Get the number and types of atoms for each species type
    natypes1 = molecules_getnatomtypes(spc1,atypes1,.True.)
    natypes2 = molecules_getnatomtypes(spc2,atypes2,.True.)

    !** Get the charges on each atom
    If (no_charges) Then
    Else
      Call config_getmolq(species,spc1,m1,charges1)
      Call config_getmolq(species,spc2,m2,charges2)
    Endif
    !** Do the mol-mol evaluations, store and whatever detail and 
    !** make sure the mol-mol level interactions are correctly summed 
    Select Case(branch)
    Case(-2)  !** only total storage associated, no indices
      ssbasic_mmint = pairmodel_mmint4(ffparams%pparams,nrgptr,nrgptr, &
          atypes1,atVecs1,charges1,natoms1,atypes2,atVecs2, &
          charges2,natoms2,simcell) 

    Case(-1)  !** only total storage associated, with indices
      ssbasic_mmint = pairmodel_mmint3(ffparams%pparams,abptr(indx(1)), &
          baptr(indx(2)),atypes1,atVecs1,charges1,natoms1,atypes2,atVecs2, &
          charges2,natoms2,simcell) 
      Call store_sum(abptr(indx(1)),.True.)

    Case(0)   !** condensed storage (1D array) for each molecule (m1 in abptr)
      ssbasic_mmint = pairmodel_mmint2(ffparams%pparams,abptr(indx(1)), &
          baptr(indx(2)),atypes1,atVecs1,charges1,natoms1,atypes2,atVecs2, &
          charges2,natoms2,simcell) 
      Call store_sum(abptr(indx(1)),.True.)

    Case(11)  !** full atom-atom storage
      ssbasic_mmint = pairmodel_mmint(ffparams%pparams,abptr(indx(1)), &
          baptr(indx(2)),atypes1,atVecs1,charges1,natoms1,atypes2,atVecs2, &
          charges2,natoms2,simcell) 
      Call store_sum(abptr(indx(1)),.True.)

    Case Default
      Write(0,'(2a,i4,a,i10)') __FILE__,": ",__LINE__, &
          ' Problem with molecule-molecule storage pointer(s), branch ID: ', &
          branch
      Stop
    End Select

#ifdef DEBUG
    !** Give feedback if there was a problem
    Call checkandstop(ssbasic_mmint,__FILE__,__LINE__, &
        ' problem with BASIC forcefield calculation spc1,mol1,spc2,mol2,id:', &
        (/spc1,m1,spc2,m2,branch/))
#endif

    !** Smooth the mol-mol interactions if necessary
    If ((ffparams%smooth).And.(ffparams%comcut)) Then
      If (ctcdist2 > ffparams%smooth_rad2) Then
        ctcdist = Sqrt(ctcdist2)
        Call smooth_perform(ffparams%smooth_params,ctcdist,abptr(indx(1)))
        Call smooth_perform(ffparams%smooth_params,ctcdist,baptr(indx(2)))
      End If
    End If

  End Function ssbasic_mmint

  !----------------------------------------------------------------------
  ! Handles a single ATOM-MOLECULE interaction.  It decomposes the
  ! interactions into atom-atom interactions and passes the information
  ! to the pairmodel module in one call.  
  ! Requires:  ffparams -- spc-spc pairwise forcefield parameters
  !            results -- atm1-mol2 Total energy structure output
  !            spc1 -- 1st species number
  !            mol1 -- 1st molecule number
  !            atm1 -- 1st atom number
  !            spc2 -- 2nd species number
  !            mol2 -- 2st molecule number
  !            natoms2 -- number of atoms in 2nd species
  !            species -- species data structure
  !            simcell -- simulation cell data structure
  !            fast -- True (False) => evaluate "Fast" (Slow) interactions
  !----------------------------------------------------------------------
  Logical Function ssbasic_amint(ffparams,results,spc1,mol1,atm1, &
      spc2,mol2,natoms2,species,simcell,fast,hot)
    Type(SSBasicParams), Intent(In)              :: ffparams
    Type(EnergyPlus), Intent(InOut)              :: results
    Integer, Intent(In)                          :: spc1,mol1,atm1
    Integer, Intent(In)                          :: spc2,mol2,natoms2
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species
    Type(SimCell_Params), Intent(In)             :: simcell  
    Logical, Intent(In)                          :: fast
    Real(kind=RDbl), Dimension(:), Optional :: hot

    Integer                              :: unit,nderivs
    Integer                              :: natypes2
    Logical, Save                        :: firsttime = .True.
    Real(kind=RDbl)                      :: ctcdist2,ctcdist
    Type(VecType)                        :: com1,com2,sepvec
    Type(EnergyPlus), Save               :: opposite1,opposite0,opposite2
    Type(VecType), Dimension(1)          :: atVec1
    Type(VecType), Dimension(natoms2)    :: atVecs2
    Integer, Dimension(1)                :: atype1
    Integer, Dimension(natoms2)          :: atypes2
    Real(kind=RDbl), Dimension(1)        :: charge1
    Real(kind=RDbl), Dimension(natoms2)  :: charges2
    Real(kind=RDbl), Dimension(SIZE_OF_HOT) :: hotatm

    !** Allocate the dummy opposite interaction structure
    If (firsttime) Then
      firsttime=.False.
      Call storebase_init(opposite1,1)
      Call storebase_init(opposite0,0)
      Call storebase_init(opposite2,2)
    End If

    !** Check the center-to-center cut-off if necessary
    If (ffparams%comcut) Then
      com1 = config_getMolecCTR(species,spc1,mol1)
      com2 = config_getMolecCTR(species,spc2,mol2)
      Call simcell_minimage(simcell,com1,com2,sepvec,ctcdist2)

      !** Determine if the molecule-molecule evaluation should not be done
      If (ctcdist2 > ffparams%cutrad2) Return
    End If

    !** Extract the atomic coordinates for each atom
    atVec1 = config_getr(species,(/spc1,mol1,atm1/))
    Call config_getmolr(species,spc2,mol2,atVecs2)

    !** Get the number and types of atoms for 2nd species type
    atype1(1) = molecules_getatype(spc1,atm1)
    natypes2 = molecules_getnatomtypes(spc2,atypes2,.True.)

    !** Get the charges on each atom
    Call config_getq(species,spc1,mol1,atm1,charge1(1))
    Call config_getmolq(species,spc2,mol2,charges2)

    If (storebase_nderivs(results)==0) then
      !** Get the interactions
      If (.Not.Present(hot)) Then
        !** Get higher order terms if necessary: pairmodel_mmint5
        ssbasic_amint = pairmodel_mmint4(ffparams%pparams,results,opposite0, &
            atype1,atVec1,charge1,1,atypes2,atVecs2,charges2,natoms2,simcell) 
      Else
        ssbasic_amint = pairmodel_mmint5(ffparams%pparams,results,opposite0, &
            atype1,atVec1,charge1,1,atypes2,atVecs2,charges2,natoms2,simcell,hot) 
      End If
    ElseIf (storebase_nderivs(results)==1) then
      If (.Not.Present(hot)) Then
        !** Get higher order terms if necessary: pairmodel_mmint5
        ssbasic_amint = pairmodel_mmint4(ffparams%pparams,results,opposite1, &
            atype1,atVec1,charge1,1,atypes2,atVecs2,charges2,natoms2,simcell) 
      Else
        ssbasic_amint= pairmodel_mmint5(ffparams%pparams,results,opposite1, &
            atype1,atVec1,charge1,1,atypes2,atVecs2,charges2,natoms2,&
            simcell,hot) 
      End If
    Else
      !!** Assume that you call for higher order derivatives by default
      ssbasic_amint = pairmodel_mmint5(ffparams%pparams,results,opposite2, &
          atype1,atVec1,charge1,1,atypes2,atVecs2,charges2,natoms2,simcell,hot) 
    endif

    !** Smooth the mol-mol interactions if necessary
    If ((ffparams%smooth).And.(ffparams%comcut)) Then
      If (ctcdist2 > ffparams%smooth_rad2) Then
        ctcdist = Sqrt(ctcdist2)
        Call smooth_perform1(ffparams%smooth_params,ctcdist,results)
      End If
    End If

  End Function ssbasic_amint

  !----------------------------------------------------------------------------
  ! Displays the parameters for this interaction type
  ! Requires:  params -- SSBasic parameters
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional display unit number
  !----------------------------------------------------------------------------
  Subroutine ssbasic_display(params, indent, unit)
    Type(SSBasicParams), Intent(In)   :: params
    Integer, Intent(In)               :: indent
    Integer, Intent(In)               :: unit

    Integer                               :: nfields,i,j,aaunit
    Integer                               :: natypes1,natypes2,natoms1,natoms2
    Logical                               :: foundi,foundj
    Character(len=indent)                 :: blank
    Character(len=strLen)                 :: filename,name1,name2
    Character(len=strLen), Dimension(10)  :: fields
    Integer, Dimension(1000)              :: atypes1,atypes2

    blank = Repeat(' ',indent)

    nfields = split(params%line, fields)
    Write(unit,'(a,2a)',Advance='No') blank, "Input Line: "
    Do j = 3, nfields
      Write(unit,'(a,2x)',Advance='No') Trim(fields(j))
    End Do

    Write(unit,*)
    Write(unit,'(a,a,1x,l2)') blank, "Fast interaction?  ", params%fast
    Write(unit,'(a,a,1x,l2)') blank, "Use neighbor-list? ", params%neighbor
    Write(unit,'(a,a,1x,l2)') blank, "COM-based cut-off? ", params%comcut
    If (params%comcut) Then
      Write(unit,'(2x,2a,f8.3)') blank, "cut-off radius = ",params%cutrad
    End If
    Write(unit,'(a,a,1x,l2)') blank, "COM-COM smoothing? ", params%smooth
    If (params%smooth) Then
      Call smooth_display(params%smooth_params,indent+2,unit)
    End If

    Write(unit,'(2a)') blank, 'Atom-Atom PAIR PARAMETERS:'
    Call file_gettype(d_aa_file,filename,aaunit)
    Write(unit,'(2x,3a)') blank, 'Information taken from: ',Trim(filename)
    Write(unit,'(2a)') blank, 'Overview by atomic pairs:'

    !** Get lists of atom types in each species
    natypes1 = molecules_getnatomtypes(params%spc1,atypes1)
    natypes2 = molecules_getnatomtypes(params%spc2,atypes2)
    natoms1 = molecules_getnatoms(params%spc1)
    natoms2 = molecules_getnatoms(params%spc2)

    Do i = 1,Size(params%pparams,1)
      name1 = atom_getname(i)
      Do j = i,Size(params%pparams,2)
        name2 = atom_getname(j)

        !** make sure that at least one species contains these atom types
        foundi = .False.
        foundj = .False.
        If ((findint(atypes1(1:natypes1),i) /= 0) .Or. &
            (findint(atypes2(1:natypes2),i) /= 0)) foundi = .True.
        If ((findint(atypes1(1:natypes1),j) /= 0) .Or. &
            (findint(atypes2(1:natypes2),j) /= 0)) foundj = .True.
        If ((.Not. foundi).Or.(.Not. foundj)) Cycle

        Write(unit,'(a,2(a,i1),a,4a,t39)', Advance='No') blank, &
            '(',i,',',j,')  ',Trim(name1),', ',Trim(name2),': '
        Call pairmodel_display(params%pparams(i,j),indent+2,unit)
      End Do
    End Do

  End Subroutine ssbasic_display

  !----------------------------------------------------------------------------
  ! Cleans up after usage
  ! Requires:  params -- SSBasic parameters
  !----------------------------------------------------------------------------
  Subroutine ssbasic_clean(params)
    Type(SSBasicParams), Intent(InOut)   :: params

    Integer               :: i,j,error

    Do i = 1,Size(params%pparams,1)
      Do j = 1,Size(params%pparams,2)
        Call pairmodel_clean(params%pparams(i,j))
      End Do
    End Do

    If (Associated(params%smooth_params)) Then
      Call smooth_clean(params%smooth_params)
      Deallocate(params%smooth_params, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__)
    End If

    Deallocate(params%pparams, STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'pparams')

  End Subroutine ssbasic_clean

!------------------------------------------------------------------------------
! Tina added
! Subroutine to get potential parameters, the array pot_params contains
! A, B, C, D, hicut, locut and returns zero if the parameters are not
! defined for the potential 
!------------------------------------------------------------------------------
 Subroutine ssbasic_getpotparameters(params, a1, a2, pot_params)
   Type(ssbasicparams), Intent(in)             :: params
   Integer, Intent(IN)                         :: a1, a2
   Real(kind = Rdbl), Dimension(6),INTENT(OUT) :: pot_params
 
   Call pairmodel_getpotparameters(params%pparams(a1,a2),pot_params)

 End Subroutine ssbasic_getpotparameters 
  
End Module ssbasic


#ifdef OBSOLETE

!!$THIS SUBROUTINE SHOULD BE REMOVED
  !----------------------------------------------------------------------
  ! This routine calculates interactions between atom-atom_num and
  ! atoms-    1:(atom_num-4)
  ! **********THERE ARE LOT OF NOT_REQUIRED VARIABLES HERE CLEAN!!!
  ! Maybe obsolete, REMOVE
  !----------------------------------------------------------------------
  Subroutine ssbasic_IntraGetinteraction(species,params, molecPtr, &
      molec,atom_num, intrapot, ljflag)
    Type(AtMolCoords), Intent(In)               :: species
    Type(SSBasicParams), Intent(In)             :: params
    Type(MolecularParams), Intent(In) :: molecPtr
    Integer, Intent(In)                         :: molec,atom_num
    Real(kind=RDbl), Intent(Out)                :: intrapot
    Logical, Intent(Out)                        :: ljflag

    Integer          :: nvecs, natom_types, i
    Integer, Dimension(species%natoms) :: alists
    Integer, Dimension(2,species%natoms) :: apairs
    Type(VecType), Dimension(species%natoms-4) :: sepvec
    Type(VecType)         :: atVecs1,atVecs2

    nvecs = 0
    intrapot = 0.0_RDbl
    natom_types = molecule_getnatomtypes(molecPtr, alists, .True.)
  
    If( atom_num > 4) Then
      atVecs1 = species%coords(atom_num, molec)%rp

      Do i=1,atom_num-4
      !** The atoms that are separated by at least three chemical bonds
      !** are calculated
       nvecs = nvecs + 1
       atVecs2 = species%coords(i, molec)%rp
       sepvec(i) = atVecs1 - atVecs2
       apairs(1, i) = alists(i)
       apairs(2, i) = alists(atom_num)
      End Do

      Call lj_calc_interaction(params%ffparams%lj, apairs, &
          sepvec(1:nvecs), intrapot, ljflag)
     
    End If

  End Subroutine ssbasic_IntraGetinteraction
#endif




