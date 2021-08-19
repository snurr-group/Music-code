!------------------------------------------------------------------------------
! This module contains all the data structure and the routines necessary 
! to get the PAIRWISE INTERACTIONS.  
!
! The initialization routine is "pairmodel_init".  It recognizes the potential
! types and farms out the rest of the initilization to the modules containing
! the specific potentials.  It takes care of mixing for those potential types
! that support it.
!------------------------------------------------------------------------------

Module pairmodel

  Use defaults, Only: RDbl, strLen, lstrLen, xlstrLen, SIZE_OF_HOT, one, zero
  Use file, Only: file_open
  Use utils, Only: filesrchstr, toupper, tolower, split, deallocerrdisplay, &
      allocerrdisplay, cleanstring
  Use simcell, Only: SimCell_Params, simcell_minimage, simcell_minwidth
  Use vector, Only: VecType, vector_sumVecs, Assignment(=), Operator(+), &
      Operator(-), Operator(*), Operator(/), vector_display, mag
  Use atom, Only: atom_getntypes, atom_getname, atom_getsymbol
  Use molecules, Only: molecules_getnsorbs
  Use lj, Only: lj_idstring, lj_displaySigEps, lj_isinit, LJ_Params, lj_init, &
      lj_calc_mix, lj_calc_ab, lj_displayCutOffs, lj_display, lj_snglint, &
      lj_getcutoff, lj_chkequal,  lj_snglint2, lj_getABCD, lj_getcut
  Use lj206, Only: LJ206_Params, lj206_idstring, lj206_displaySigEps, &
      lj206_calc_mix, lj206_calc_ab, lj206_init, lj206_snglint, &
      lj206_displayCutOffs, lj206_display, lj206_isinit, lj206_getcutoff, &
      lj206_chkequal
  Use ljel, Only: ljel_idstring  
  Use coulomb, Only: coulomb_idstring, coulomb_isinit, Coulomb_Params, &
      coulomb_init, coulomb_disp, coulomb_displaycutoffs, coulomb_snglint, &
      coulomb_getcutoff, coul_chkequal
  Use buck, Only: Buckingham_Params, buck_idstring, buck_init, &
      buck_calc_interaction, buck_isinit, buck_displayCutoffs, &
      buck_displayParams, buck_snglint, buck_getcutoff, buck_chkequal
  Use benzpot, Only: BENZPOT_Params, benzpot_idstring, benzpot_init, &
      benzpot_snglint, benzpot_isinit, benzpot_display, benzpot_getcutoff, &
      benzpot_chkequal
  Use harmonic, Only: HarmonicParams, harmonic_idstring, harmonic_init, &
      harmonic_simpleinit, harmonic_snglint, harmonic_display
  Use storebase, Only: EnergyPlus, storebase_init, storebase_sumincarray
  Use store, Only: Store_Level_Pair, store_pack2D, store_unpack2Dinc, &
      store_display, store_sum, store_unpack2Dfast, store_idbranch

  Implicit None
  Save

  Private
  Public :: Pairmodel_Params, pairmodel_init, pairmodel_simpleinit, &
      pairmodel_mmint, pairmodel_mmint2, pairmodel_mmint3, pairmodel_mmint4, pairmodel_mmint5, &
      pairmodel_chkcutoffs, pairmodel_nullset, pairmodel_countset, &
      pairmodel_chkarray, pairmodel_makesym, pairmodel_display, &
      pairmodel_clean,  pairmodel_msint5, pairmodel_getpotparameters

  !** NAG compiler bug workaround
  Integer :: dummyPairmodel = strLen

  Type Pairmodel_Params
    Integer                          :: nterms,atype1,atype2
    Character(len=strLen)            :: identifiers
    Type(LJ206_Params), Pointer      :: lj206
    Type(LJ_Params), Pointer         :: lj
    Type(Coulomb_Params), Pointer    :: coul
    Type(Buckingham_Params), Pointer :: buck
    Type(Benzpot_Params), Pointer    :: cd    
    Type(HarmonicParams), Pointer    :: spr
  End Type Pairmodel_Params
      
Contains
  !-----------------------------------------------------------------------
  ! This subroutine initializes the parameters for the interaction of
  ! atom types "atype1" and "atype2".  The full 2-D array of pairmodel
  ! parameters must be passed in here because it must have access to 
  ! other parameters for possible parameter mixing.  If it needs a 
  ! parameter that isn't initialized yet, it will call itself recursively.
  ! Requires:  pparams -- full set of pairmodel parameter types
  !            atype1 -- 1st atom type
  !            atype2 -- 2nd atom type
  !            toinit -- string containing possible potential model IDs
  !            filename -- file to find atom-atom interactions in
  !-----------------------------------------------------------------------
  Recursive Subroutine pairmodel_init(pparams,atype1,atype2,toinit,&
      filename,has_charges)
    Type(Pairmodel_Params), Dimension(:,:), Intent(InOut) :: pparams
    Integer, Intent(In)                                   :: atype1,atype2
    Character(*), Intent(In)                              :: toinit,filename
    Logical, Intent(out)                                  :: has_charges
    Integer                              :: i,j,n,unit,error
    Integer                              :: lineno,npotmodels,nfields
    Character(len=strLen)                :: name1,name2
    Character(len=xlstrLen)              :: line
    Character(len=strLen), Dimension(30) :: id,fields

    has_charges=.False.
    !** save the atom types
    pparams(atype1,atype2)%atype1 = atype1
    pparams(atype1,atype2)%atype2 = atype2

    name1 = atom_getname(atype1)
    name2 = atom_getname(atype2)
    !** get the model IDs that need to be initialized this time
    npotmodels = pairmodel_ctrlinit(pparams(atype1,atype2),toinit,id)

    !** Open file
    unit = file_open(filename,110)

    !** Find the string in the file with atom names and the correct ID string
    Do i = 1,npotmodels
      lineno = filesrchstr(unit,(/name1,name2,id(i)/),line,.True.)

      !** If not found, die a horrible death
      If (lineno == 0) Then
        Write(0,'(2a,i4,3a)') __FILE__,': ',__LINE__,'  Failed to find ',&
            Trim(id(i)),' potential model ID string'
        Write(0,'(4x,6a)') ' for atom-atom pair ', &
            Trim(name1),'-',Trim(name2),' in file: ',Trim(filename)
        n = pairmodel_ctrlinit(pparams(atype1,atype2),toinit,id)
        Write(0,'(4x,30(a,2x))') ' Still need to find these IDs: ', &
            (Trim(id(j)),j=1,n)
        Stop

      Else
        !** Initialize the atype-atype pair
        Select Case (ToUpper(Trim(id(i))))
        Case (lj_idstring)
          !** Check if we need to calculate mixed parameters
          If (Index(ToUpper(line),'MIX') /= 0) Then

            !** make sure the pure terms are associated, if not get them
            If (.Not. Associated(pparams(atype1,atype1)%lj)) Then
              Call pairmodel_init(pparams,atype1,atype1,lj_idstring,&
                  filename,has_charges)
            End If
            If (.Not. Associated(pparams(atype2,atype2)%lj)) Then
              Call pairmodel_init(pparams,atype2,atype2,lj_idstring,filename,&
                  has_charges)
            End If

            !** perform the actual mixing and calculate the A&B terms
            Allocate(pparams(atype1,atype2)%lj,stat=error)    
            If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'lj pointer')
            Call lj_calc_mix(pparams(atype1,atype2)%lj,&
                pparams(atype1,atype1)%lj,pparams(atype2,atype2)%lj,line)
            Call lj_calc_AB(pparams(atype1,atype2)%lj)

          Else
            Allocate(pparams(atype1,atype2)%lj,stat=error)    
            If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'lj pointer')
            Call lj_init(pparams(atype1,atype2)%lj,line)
          End If

        Case (lj206_idstring)
          !** Check if we need to calculate mixed parameters
          If (Index(ToUpper(line),'MIX') /= 0) Then
            Write(0,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
            Write(0,*) 'WARNING: Mixing with LJ20-6 has not been checked'

            !** make sure the pure terms are associated, if not get them
            If (.Not. Associated(pparams(atype1,atype1)%lj206)) Then
              Call pairmodel_init(pparams,atype1,atype1,&
                  lj206_idstring,filename,has_charges)
            End If
            If (.Not. Associated(pparams(atype2,atype2)%lj206)) Then
              Call pairmodel_init(pparams,atype2,atype2,&
                  lj206_idstring,filename,has_charges)
            End If

            !** perform the actual mixing and calculate the A&B terms
            Allocate(pparams(atype1,atype2)%lj206,stat=error)    
            If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'lj206 pointer')
            Call lj206_calc_mix(pparams(atype1,atype2)%lj206,&
                pparams(atype1,atype1)%lj206,pparams(atype2,atype2)%lj206,line)
            Call lj206_calc_AB(pparams(atype1,atype2)%lj206)

          Else
            Allocate(pparams(atype1,atype2)%lj206,stat=error)    
            If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'lj206 Pointer')
            Call lj206_init(pparams(atype1,atype2)%lj206,line)
          End If

        Case (ljel_idstring)
          has_charges=.True.
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              'LJ-Coulomb combination depreciated, use separate specifications'
          Stop

        Case (coulomb_idstring)
          has_charges=.True.
          Allocate(pparams(atype1,atype2)%coul,stat=error)    
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'coul pointer')
          Call coulomb_init(pparams(atype1,atype2)%coul,line)

        Case (buck_idstring)
          Allocate(pparams(atype1,atype2)%buck,stat=error)    
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'buck pointer')
          Call buck_init(pparams(atype1,atype2)%buck,line)

        Case (benzpot_idstring)
          Allocate(pparams(atype1,atype2)%cd,stat=error)    
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'benzpot pointer')
          Call benzpot_init(pparams(atype1,atype2)%cd,line)

        Case (harmonic_idstring)
          Allocate(pparams(atype1,atype2)%spr,stat=error)    
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'spring pointer')
          nfields = split(line,fields)
          Call harmonic_init(pparams(atype1,atype2)%spr, &
              fields(nfields-1:nfields))
        
        Case ('OFF')
          !** Turn off pairwise potentials for this pair, nullify and return
          Call pairmodel_nullset(pparams(atype1,atype2))
          Return

        Case Default
          Write(0,'(2a,i4,2a)') __FILE__,': ',__LINE__, &
              ' Missing initialization code for: ',ToUpper(Trim(id(i)))
        End Select

      End If
    End Do

  End Subroutine pairmodel_init

  !-----------------------------------------------------------------------
  ! Determine which model types in a set have been initialized and which
  ! still need to be.  Takes a string containing keys as input and spits
  ! out an array of identifiers for the potential models that, according
  ! to the input string, still need to be initialized.  Returns the number
  ! of models that still need to be initialized.
  ! Requires:  pparams -- a single pairmodel parameter type
  !            line -- string containing possible potential model ids
  !            init_strings -- array of strings to initialize
  !-----------------------------------------------------------------------
  Integer Function pairmodel_ctrlinit(pparams,line,init_strings)
    Type(Pairmodel_Params), Intent(In)  :: pparams
    Character(*), Intent(In)            :: line
    Character(len=strLen), Dimension(:) :: init_strings

    Integer                               :: i,nstrings
    Character(len=strLen), Dimension(100) :: fields

    pairmodel_ctrlinit = 0
    init_strings = ''

    !** Get the potential model types and number
    nstrings = split(line,fields)
    If (nstrings == 0) Return

    Do i = 1,nstrings
      Select Case(Trim(ToUpper(fields(i))))
      Case (lj_idstring)
        If (.Not. Associated(pparams%lj)) Then
          pairmodel_ctrlinit = pairmodel_ctrlinit + 1
          init_strings(pairmodel_ctrlinit) = lj_idstring
        End If

      Case (lj206_idstring)
        If (.Not. Associated(pparams%lj206)) Then
          pairmodel_ctrlinit = pairmodel_ctrlinit + 1
          init_strings(pairmodel_ctrlinit) = lj206_idstring
        End If

      Case (coulomb_idstring)
        If (.Not. Associated(pparams%coul)) Then
          pairmodel_ctrlinit = pairmodel_ctrlinit + 1
          init_strings(pairmodel_ctrlinit) = coulomb_idstring
        End If

      Case (buck_idstring)
        If (.Not. Associated(pparams%buck)) Then
          pairmodel_ctrlinit = pairmodel_ctrlinit + 1
          init_strings(pairmodel_ctrlinit) = buck_idstring
        End If

      Case (benzpot_idstring)
        If (.Not. Associated(pparams%cd)) Then
          pairmodel_ctrlinit = pairmodel_ctrlinit + 1
          init_strings(pairmodel_ctrlinit) = benzpot_idstring
        End If

      Case Default
        Write(0,'(2a,i4,3a)') __FILE__,": ",__LINE__, &
            ' Unable to understand potential type string: "',Trim(fields(i)),'"'
        Stop

      End Select
    End Do

  End Function pairmodel_ctrlinit

  !--------------------------------------------------------------------
  ! Nullifys the pairwise potential model pointers in the set
  ! Requires:  pparams -- a single set of pairmodel parameters
  !--------------------------------------------------------------------
  Subroutine pairmodel_nullset(pparams)
    Type(Pairmodel_Params), Intent(Out)  :: pparams

    pparams%nterms = 0
    pparams%atype1 = 0
    pparams%atype2 = 0
    pparams%identifiers = ''

    !** nullify the pointers
    Nullify(pparams%lj)
    Nullify(pparams%lj206)
    Nullify(pparams%coul)
    Nullify(pparams%buck)
    Nullify(pparams%cd)
    Nullify(pparams%spr)

  End Subroutine pairmodel_nullset

  !-----------------------------------------------------------------------
  ! This subroutine initializes a single atype1-atype2 interaction without
  ! input from the atom-atom file
  ! Requires:  pparams -- single pair model parameter type
  !            atype1 -- 1st atom type
  !            atype2 -- 2nd atom type
  !            idstring -- string containing possible potential model IDs
  !            nums -- array of parameter constants for the models
  !-----------------------------------------------------------------------
  Subroutine pairmodel_simpleinit(pparams,atype1,atype2,idstring,nums)
    Type(Pairmodel_Params), Intent(Out)        :: pparams
    Integer, Intent(In)                        :: atype1,atype2
    Character(*), Intent(In)                   :: idstring
    Real(kind=RDbl), Dimension(:), Intent(In)  :: nums

    Integer                              :: i,j,n,unit,error
    Integer                              :: nlo,nhi,incr
    Character(len=strLen), Dimension(30) :: id

    !** Nullify the pointer set
    Call pairmodel_nullset(pparams)

    !** Save the atom types and idstring
    pparams%atype1 = atype1
    pparams%atype2 = atype2
    pparams%identifiers = idstring

    !** Split the idstring and initialize interactions in sequence
    pparams%nterms = split(idstring,id)

    nhi = 0
    Do i = 1,pparams%nterms
      !** determine the new start of the number list
      nlo = nhi + 1

      !** Initialize the atype-atype pair
      Select Case (ToUpper(Trim(id(i))))
      Case (lj_idstring)
        incr = 2   
        Allocate(pparams%lj,stat=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'lj pointer')
        nhi = nlo + incr - 1
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            ' LJ simple init routine not yet available'
        Stop
      
      Case (lj206_idstring)
        incr = 2   
        Allocate(pparams%lj206,stat=error)    
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'lj206 pointer')
        nhi = nlo + incr - 1
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            ' LJ 20-6 simple init routine not yet available'
        Stop
      
      Case (ljel_idstring)
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            'LJ-Coulomb combination depreciated, use separate specifications'
        Stop
      
      Case (coulomb_idstring)
        Allocate(pparams%coul,stat=error)    
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'coul pointer')
        nhi = nlo + incr - 1
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            ' Coulombic simple init routine not yet available'
        Stop
      
      Case (buck_idstring)
        Allocate(pparams%buck,stat=error)    
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'buck pointer')
        nhi = nlo + incr - 1
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            ' Buckingham simple init routine not yet available'
        Stop
      
      Case (benzpot_idstring)
        Allocate(pparams%cd,stat=error)    
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'benzpot pointer')
        nhi = nlo + incr - 1
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            ' Benzene Pot (C&D) simple init routine not yet available'
        Stop

      Case (harmonic_idstring)
        incr = 2
        Allocate(pparams%spr,stat=error)    
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'spring pointer')
        nhi = nlo + incr - 1
        Call harmonic_simpleinit(pparams%spr,nums(nlo:nhi))
      
      Case ('OFF')
        !** Turn off pairwise potentials for this pair, nullify and return
        Call pairmodel_nullset(pparams)
        Return
      
      Case Default
        Write(0,'(2a,i4,2a)') __FILE__,': ',__LINE__, &
            ' Missing initialization code for: ',ToUpper(Trim(id(i)))
      End Select
    End Do
    
  End Subroutine pairmodel_simpleinit

  !--------------------------------------------------------------------
  ! Link the pointers in one pairmodel parameter set to another.  Also
  ! copies the normal fields.
  ! Requires:  params1 -- 1st set of pairmodel parameters (new)
  !            params2 -- 2nd set of pairmodel parameters (original)
  !--------------------------------------------------------------------
  Subroutine pairmodel_link(params1,params2)
    Type(Pairmodel_Params), Intent(InOut)  :: params1,params2

    params1%nterms = params2%nterms
    params1%atype1 = params2%atype1
    params1%atype2 = params2%atype2
    params1%identifiers = params2%identifiers

    !** link the pointers
    params1%lj => params2%lj   
    params1%lj206 => params2%lj206
    params1%coul => params2%coul 
    params1%buck => params2%buck 
    params1%cd => params2%cd   

  End Subroutine pairmodel_link

  !--------------------------------------------------------------------
  ! Sets the number of terms and identifiers in the parameter set
  ! Requires:  pparams -- a single set of pairmodel parameters
  !--------------------------------------------------------------------
  Subroutine pairmodel_countset(pparams)
    Type(Pairmodel_Params), Intent(InOut)  :: pparams

    pparams%nterms = 0
    pparams%identifiers = ''

    If (Associated(pparams%lj)) Then
      pparams%nterms = pparams%nterms + 1
      Write(pparams%identifiers,'(a,2x,a)') Trim(pparams%identifiers), &
          Trim(lj_idstring)
    End If


    If (Associated(pparams%lj206)) Then
      pparams%nterms = pparams%nterms + 1
      Write(pparams%identifiers,'(a,2x,a)') Trim(pparams%identifiers), &
          Trim(lj206_idstring)
    End If

    If (Associated(pparams%coul)) Then
      pparams%nterms = pparams%nterms + 1
      Write(pparams%identifiers,'(a,2x,a)') Trim(pparams%identifiers), &
          Trim(coulomb_idstring)
    End If

    If (Associated(pparams%buck)) Then
      pparams%nterms = pparams%nterms + 1
      Write(pparams%identifiers,'(a,2x,a)') Trim(pparams%identifiers), &
          Trim(buck_idstring)
    End If

    If (Associated(pparams%cd)) Then
      pparams%nterms = pparams%nterms + 1
      Write(pparams%identifiers,'(a,2x,a)') Trim(pparams%identifiers), &
          Trim(benzpot_idstring)
    End If

  End Subroutine pairmodel_countset

  !----------------------------------------------------------------------------
  ! Handles the interaction between two "molecules".  These "molecules" are
  ! each defined by an array of coordinates and atom types.  The minimum image
  ! convention is employed here to get separation vectors.  
  ! Requires:  pparams -- pairmodel forcefield parameters indexed by atom-type
  !            grp1out -- group 1 forcefield output (mm field associated)
  !            grp2out -- group 2 forcefield output (mm field associated)
  !            atype1 -- array of atom types for molecule (group) 1
  !            atVecs1 -- array of coordinate vectors for molecule (group) 1
  !            charges1 -- charges on atoms of group 1
  !            natoms1 -- number of atoms in molecule (group) 1
  !            atype2 -- array of atom types for molecule (group) 2
  !            atVecs2 -- array of coordinate vectors for molecule (group) 2
  !            charges2 -- charges on atoms of group 2
  !            natoms2 -- number of atoms in molecule (group) 2
  !            simcell -- simulation cell data structure
  !----------------------------------------------------------------------------
  Logical Function pairmodel_mmint(pparams,grp1out,grp2out,atype1,atVecs1, &
        charges1,natoms1,atype2,atVecs2,charges2,natoms2,simcell) 
    Type(Pairmodel_Params), Dimension(:,:), Intent(In) :: pparams
    Type(Store_Level_Pair), Intent(InOut)              :: grp1out,grp2out
    Type(VecType), Dimension(:), Intent(In)            :: atVecs1,atVecs2
    Integer, Dimension(:), Intent(In)                  :: atype1,atype2
    Real(kind=RDbl), Dimension(:), Intent(In)          :: charges1,charges2
    Integer, Intent(In)                                :: natoms1,natoms2
    Type(SimCell_Params), Intent(In)                   :: simcell  

    Integer                                      :: a1,a2,idx
    Logical                                      :: flag
    Type(VecType), Dimension(natoms1,natoms2)    :: sepvec

    !** set defaults
    pairmodel_mmint = .True.

    !** Get the minimum image distance between all points in atVec1 and
    !** all points in atVec2, sepvec = atvec1 - atvec2
    Call simcell_minimage(simcell,atVecs1,atVecs2,sepvec)

    !** Do the evaluations, %mm(atom1,atom2) stores atom-atom interactions
    idx = 0
    Do a1 = 1,natoms1
      Do a2 = 1,natoms2
        If (Associated(pparams(atype1(a1),atype2(a2))%lj)) Then
          pairmodel_mmint = lj_snglint(pparams(atype1(a1),atype2(a2))%lj, &
              sepvec(a1,a2),grp1out%mm(a1,a2)%total,grp2out%mm(a1,a2)%total)
          If (.Not. pairmodel_mmint) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%lj206)) Then
          pairmodel_mmint = lj206_snglint(pparams(atype1(a1),atype2(a2))%lj206, &
              sepvec(a1,a2),grp1out%mm(a1,a2)%total,grp2out%mm(a1,a2)%total)
          If (.Not. pairmodel_mmint) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%coul)) Then
          pairmodel_mmint = coulomb_snglint(pparams(atype1(a1),atype2(a2))%coul,&
              sepvec(a1,a2),charges1(a1),charges2(a2), &
              grp1out%mm(a1,a2)%total,grp2out%mm(a1,a2)%total)
          If (.Not. pairmodel_mmint) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%buck)) Then
          pairmodel_mmint = buck_snglint(pparams(atype1(a1),atype2(a2))%buck, &
              sepvec(a1,a2),grp1out%mm(a1,a2)%total,grp2out%mm(a1,a2)%total)
          If (.Not. pairmodel_mmint) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%cd)) Then
          pairmodel_mmint = benzpot_snglint(pparams(atype1(a1),atype2(a2))%cd, &
              sepvec(a1,a2),grp1out%mm(a1,a2)%total,grp2out%mm(a1,a2)%total)
          If (.Not. pairmodel_mmint) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%spr)) Then
          pairmodel_mmint = harmonic_snglint(pparams(atype1(a1),atype2(a2))%spr,&
              sepvec(a1,a2),grp1out%mm(a1,a2)%total,grp2out%mm(a1,a2)%total)
          If (.Not. pairmodel_mmint) Return
        End If

      End Do
    End Do

  End Function pairmodel_mmint

  !----------------------------------------------------------------------------
  ! Same as mmint, but takes in binfo associated stuff
  ! Handles the interaction between two "molecules".  These "molecules" are
  ! each defined by an array of coordinates and atom types.  The minimum image
  ! convention is employed here to get separation vectors.  
  ! Requires:  pparams -- pairmodel forcefield parameters indexed by atom-type
  !            grp1out -- group 1 forcefield output (mm binfo associated)
  !            grp2out -- group 2 forcefield output (mm binfo associated)
  !            atype1 -- array of atom types for molecule (group) 1
  !            atVecs1 -- array of coordinate vectors for molecule (group) 1
  !            charges1 -- charges on atoms of group 1
  !            natoms1 -- number of atoms in molecule (group) 1
  !            atype2 -- array of atom types for molecule (group) 2
  !            atVecs2 -- array of coordinate vectors for molecule (group) 2
  !            charges2 -- charges on atoms of group 2
  !            natoms2 -- number of atoms in molecule (group) 2
  !            simcell -- simulation cell data structure
  !----------------------------------------------------------------------------
  Logical Function pairmodel_mmint2(pparams,grp1out,grp2out,atype1,atVecs1, &
        charges1,natoms1,atype2,atVecs2,charges2,natoms2,simcell) 
    Type(Pairmodel_Params), Dimension(:,:), Intent(In) :: pparams
    Type(Store_Level_Pair), Intent(InOut)              :: grp1out,grp2out
    Type(VecType), Dimension(:), Intent(In)            :: atVecs1,atVecs2
    Integer, Dimension(:), Intent(In)                  :: atype1,atype2
    Real(kind=RDbl), Dimension(:), Intent(In)          :: charges1,charges2
    Integer, Intent(In)                                :: natoms1,natoms2
    Type(SimCell_Params), Intent(In)                   :: simcell  

    Integer                                      :: a1,a2,idx
    Logical                                      :: flag
    Type(VecType), Dimension(natoms1,natoms2)    :: sepvec

    !** set defaults
    pairmodel_mmint2 = .True.

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,*) natoms1,natoms2

    !** Get the minimum image distance between all points in atVec1 and
    !** all points in atVec2, sepvec = atvec1 - atvec2
    Call simcell_minimage(simcell,atVecs1,atVecs2,sepvec)

    !** Do the evaluations, %binfo(atom) stores forces on atoms
    idx = 0

    Do a1 = 1,natoms1
      Do a2 = 1,natoms2
        If (Associated(pparams(atype1(a1),atype2(a2))%lj)) Then
!          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!          Write(*,*) a1,a2
!          Write(*,*) atype1(a1),atype2(a2)
!          Write(*,*) store_idbranch(grp1out),store_idbranch(grp2out)
!          Write(*,*) Size(grp1out%binfo),Size(grp2out%binfo)
          pairmodel_mmint2 = lj_snglint(pparams(atype1(a1),atype2(a2))%lj, &
              sepvec(a1,a2),grp1out%binfo(a1),grp2out%binfo(a2))
          If (.Not. pairmodel_mmint2) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%lj206)) Then
!          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          pairmodel_mmint2 = lj206_snglint(pparams(atype1(a1),atype2(a2))%lj206,&
              sepvec(a1,a2),grp1out%binfo(a1),grp2out%binfo(a2))
          If (.Not. pairmodel_mmint2) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%coul)) Then
          pairmodel_mmint2 = coulomb_snglint( &
              pparams(atype1(a1),atype2(a2))%coul,sepvec(a1,a2),&
              charges1(a1),charges2(a2),grp1out%binfo(a1),grp2out%binfo(a2))
          If (.Not. pairmodel_mmint2) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%buck)) Then
          pairmodel_mmint2 = buck_snglint(pparams(atype1(a1),atype2(a2))%buck, &
              sepvec(a1,a2),grp1out%binfo(a1),grp2out%binfo(a2))
          If (.Not. pairmodel_mmint2) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%cd)) Then
          pairmodel_mmint2 = benzpot_snglint(pparams(atype1(a1),atype2(a2))%cd, &
              sepvec(a1,a2),grp1out%binfo(a1),grp2out%binfo(a2))
          If (.Not. pairmodel_mmint2) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%spr)) Then
          pairmodel_mmint2 = harmonic_snglint( &
              pparams(atype1(a1),atype2(a2))%spr, &
              sepvec(a1,a2),grp1out%binfo(a1),grp2out%binfo(a2))
          If (.Not. pairmodel_mmint2) Return
        End If

      End Do
    End Do

  End Function pairmodel_mmint2

  !----------------------------------------------------------------------------
  ! Same as mmint, but takes in -1 associated stuff
  ! Handles the interaction between two "molecules".  These "molecules" are
  ! each defined by an array of coordinates and atom types.  The minimum image
  ! convention is employed here to get separation vectors.  
  ! Requires:  pparams -- pairmodel forcefield parameters indexed by atom-type
  !            grp1out -- group 1 forcefield output (total only associated)
  !            grp2out -- group 2 forcefield output (total only associated)
  !            atype1 -- array of atom types for molecule (group) 1
  !            atVecs1 -- array of coordinate vectors for molecule (group) 1
  !            charges1 -- charges on atoms of group 1
  !            natoms1 -- number of atoms in molecule (group) 1
  !            atype2 -- array of atom types for molecule (group) 2
  !            atVecs2 -- array of coordinate vectors for molecule (group) 2
  !            charges2 -- charges on atoms of group 2
  !            natoms2 -- number of atoms in molecule (group) 2
  !            simcell -- simulation cell data structure
  !----------------------------------------------------------------------------
  Logical Function pairmodel_mmint3(pparams,grp1out,grp2out,atype1,atVecs1, &
        charges1,natoms1,atype2,atVecs2,charges2,natoms2,simcell) 
    Type(Pairmodel_Params), Dimension(:,:), Intent(In) :: pparams
    Type(Store_Level_Pair), Intent(InOut)              :: grp1out,grp2out
    Type(VecType), Dimension(:), Intent(In)            :: atVecs1,atVecs2
    Integer, Dimension(:), Intent(In)                  :: atype1,atype2
    Real(kind=RDbl), Dimension(:), Intent(In)          :: charges1,charges2
    Integer, Intent(In)                                :: natoms1,natoms2
    Type(SimCell_Params), Intent(In)                   :: simcell  

    Integer                                      :: a1,a2,idx
    Logical                                      :: flag
    Type(VecType), Dimension(natoms1,natoms2)    :: sepvec
    Character(len=lstrlen)                       :: string

    !** set defaults
    pairmodel_mmint3 = .True.

    !** Get the minimum image distance between all points in atVec1 and
    !** all points in atVec2, sepvec = atvec1 - atvec2
    Call simcell_minimage(simcell,atVecs1,atVecs2,sepvec)

    !** Do the evaluations,  %total stores summed interactions
    idx = 0
    Do a1 = 1,natoms1
      Do a2 = 1,natoms2

        If (Associated(pparams(atype1(a1),atype2(a2))%lj)) Then
          pairmodel_mmint3 = lj_snglint(pparams(atype1(a1),atype2(a2))%lj, &
              sepvec(a1,a2),grp1out%total,grp2out%total)
          If (.Not. pairmodel_mmint3) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%lj206)) Then
          pairmodel_mmint3 = lj206_snglint(pparams(atype1(a1),atype2(a2))%lj206,&
              sepvec(a1,a2),grp1out%total,grp2out%total)
          If (.Not. pairmodel_mmint3) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%coul)) Then
          pairmodel_mmint3 = coulomb_snglint( &
              pparams(atype1(a1),atype2(a2))%coul,sepvec(a1,a2),&
              charges1(a1),charges2(a2),grp1out%total,grp2out%total)
          If (.Not. pairmodel_mmint3) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%buck)) Then
          pairmodel_mmint3 = buck_snglint(pparams(atype1(a1),atype2(a2))%buck, &
              sepvec(a1,a2),grp1out%total,grp2out%total)
          If (.Not. pairmodel_mmint3) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%cd)) Then
          pairmodel_mmint3 = benzpot_snglint( &
              pparams(atype1(a1),atype2(a2))%cd, &
              sepvec(a1,a2),grp1out%total,grp2out%total)
          If (.Not. pairmodel_mmint3) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%spr)) Then
          pairmodel_mmint3 = harmonic_snglint( &
              pparams(atype1(a1),atype2(a2))%spr, &
              sepvec(a1,a2),grp1out%total,grp2out%total)
          If (.Not. pairmodel_mmint3) Return
        End If

      End Do
    End Do

  End Function pairmodel_mmint3

  !----------------------------------------------------------------------------
  ! Same as mmint, mmint2, mmint3 but doesn't use level pair structure
  ! Handles the interaction between two "molecules".  These "molecules" are
  ! each defined by an array of coordinates and atom types.  The minimum image
  ! convention is employed here to get separation vectors.  
  ! Requires:  pparams -- pairmodel forcefield parameters indexed by atom-type
  !            grp1out -- group 1 forcefield output (total only associated)
  !            grp2out -- group 2 forcefield output (total only associated)
  !            atype1 -- array of atom types for molecule (group) 1
  !            atVecs1 -- array of coordinate vectors for molecule (group) 1
  !            charges1 -- charges on atoms of group 1
  !            natoms1 -- number of atoms in molecule (group) 1
  !            atype2 -- array of atom types for molecule (group) 2
  !            atVecs2 -- array of coordinate vectors for molecule (group) 2
  !            charges2 -- charges on atoms of group 2
  !            natoms2 -- number of atoms in molecule (group) 2
  !            simcell -- simulation cell data structure
  !----------------------------------------------------------------------------
  Logical Function pairmodel_mmint4(pparams,grp1out,grp2out,atype1,atVecs1, &
        charges1,natoms1,atype2,atVecs2,charges2,natoms2,simcell) 
    Type(Pairmodel_Params), Dimension(:,:), Intent(In) :: pparams
    Type(EnergyPlus), Intent(InOut)                    :: grp1out,grp2out
    Type(VecType), Dimension(:), Intent(In)            :: atVecs1,atVecs2
    Integer, Dimension(:), Intent(In)                  :: atype1,atype2
    Real(kind=RDbl), Dimension(:), Intent(In)          :: charges1,charges2
    Integer, Intent(In)                                :: natoms1,natoms2
    Type(SimCell_Params), Intent(In)                   :: simcell  

    Integer                                      :: a1,a2,idx
    Logical                                      :: flag
    Type(VecType), Dimension(natoms1,natoms2)    :: sepvec
    Character(len=lstrlen)                       :: string

    !** set defaults
    pairmodel_mmint4 = .True.

    !** Get the minimum image distance between all points in atVec1 and
    !** all points in atVec2, sepvec = atvec1 - atvec2
    Call simcell_minimage(simcell,atVecs1,atVecs2,sepvec)

    !** Do the evaluations 
    idx = 0
    Do a1 = 1,natoms1
      Do a2 = 1,natoms2

        If (Associated(pparams(atype1(a1),atype2(a2))%lj)) Then
          pairmodel_mmint4 = lj_snglint(pparams(atype1(a1),atype2(a2))%lj, &
              sepvec(a1,a2),grp1out,grp2out)
          If (.Not. pairmodel_mmint4) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%lj206)) Then
          pairmodel_mmint4 = lj206_snglint(pparams(atype1(a1),atype2(a2))%lj206,&
              sepvec(a1,a2),grp1out,grp2out)
          If (.Not. pairmodel_mmint4) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%coul)) Then
          pairmodel_mmint4 = coulomb_snglint( &
              pparams(atype1(a1),atype2(a2))%coul,sepvec(a1,a2),&
              charges1(a1),charges2(a2),grp1out,grp2out)
          If (.Not. pairmodel_mmint4) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%buck)) Then
          pairmodel_mmint4 = buck_snglint(pparams(atype1(a1),atype2(a2))%buck, &
              sepvec(a1,a2),grp1out,grp2out)
          If (.Not. pairmodel_mmint4) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%cd)) Then
          pairmodel_mmint4 = benzpot_snglint(pparams(atype1(a1),atype2(a2))%cd, &
              sepvec(a1,a2),grp1out,grp2out)
          If (.Not. pairmodel_mmint4) Return
        End If

        If (Associated(pparams(atype1(a1),atype2(a2))%spr)) Then
          pairmodel_mmint4 = harmonic_snglint( &
              pparams(atype1(a1),atype2(a2))%spr,sepvec(a1,a2),grp1out,grp2out)
          If (.Not. pairmodel_mmint4) Return
        End If
      End Do
    End Do

  End Function pairmodel_mmint4

  !----------------------------------------------------------------------------
  ! Same as mmint, mmint2, mmint3, mmint4 . But handles a specific case of 
  ! HOT calculations between an atom and a molecule   The minimum image
  ! convention is employed here to get separation vectors.  
  ! Requires:  pparams -- pairmodel forcefield parameters indexed by atom-type
  !            grp1out -- group 1 forcefield output (total only associated)
  !            grp2out -- group 2 forcefield output (total only associated)
  !            atype1 -- array of atom types for molecule (group) 1
  !            atVecs1 -- array of coordinate vectors for molecule (group) 1
  !            charges1 -- charges on atoms of group 1
  !            natoms1=1
  !            atype2 -- array of atom types for molecule (group) 2
  !            atVecs2 -- array of coordinate vectors for molecule (group) 2
  !            charges2 -- charges on atoms of group 2
  !            natoms2 -- number of atoms in molecule (group) 2
  !            simcell -- simulation cell data structure
  !----------------------------------------------------------------------------
  Logical Function pairmodel_mmint5(pparams,grp1out,grp2out,atype1,atVecs1, &
      charges1,natoms1,atype2,atVecs2,charges2,natoms2,simcell,hot) 
    Type(Pairmodel_Params), Dimension(:,:), Intent(In) :: pparams
    Type(EnergyPlus), Intent(InOut)                    :: grp1out,grp2out
    Type(VecType), Dimension(:), Intent(In)            :: atVecs1,atVecs2
    Integer, Dimension(:), Intent(In)                  :: atype1,atype2
    Real(kind=RDbl), Dimension(:), Intent(In)          :: charges1,charges2
    Integer, Intent(In)                                :: natoms1,natoms2
    Type(SimCell_Params), Intent(In)                   :: simcell 
    Real(kind=RDbl), Dimension(SIZE_OF_HOT) :: hot 

    Integer                                      :: a1,a2,idx
    Logical                                      :: flag, ljflag
    Type(VecType), Dimension(natoms1,natoms2)    :: sepvec
    Character(len=lstrlen)                       :: string
    Real(kind=RDbl), Dimension(SIZE_OF_HOT) :: hotatm

    hotatm = zero
    hot = zero

    !** set defaults
    pairmodel_mmint5 = .True.

    !** Get the minimum image distance between all points in atVec1 and
    !** all points in atVec2, sepvec = atvec1 - atvec2
    Call simcell_minimage(simcell,atVecs1,atVecs2,sepvec)

    !** Do the evaluations 
    idx = 0

    If (natoms1 /= 1) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'This routine calculates interactions between ', &
           '1 atom (natoms1) and 1 molecule'
      Write(*,*) 'You are doing something else. May be you should ', &
          'use pairmodel_mmint4 instead', natoms1
      Stop
    End If

    a1 = 1
    Do a2 = 1,natoms2
      If (Associated(pparams(atype1(a1),atype2(a2))%lj)) Then
        hotatm = zero
        pairmodel_mmint5 = lj_snglint(pparams(atype1(a1),atype2(a2))%lj, &
            sepvec(a1,a2),grp1out,grp2out,hotatm)
        hot = hot+hotatm
        If (.Not. pairmodel_mmint5) Return
      Else
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) 'Higher order terms can be calculated only &
            & for Lennard-Jones 12-6 interactions'
        Stop
      End If
    End Do

  End Function pairmodel_mmint5

  !----------------------------------------------------------------------------
  ! calculates molecule-sorbate interaction in one go, similar to _getmsnrgs
  ! from music-2-2. Not very neet here, this should have 
  ! functionalities similar to mmint4. But this is very fast and efficient
  ! Requires:  pparams -- pairmodel forcefield parameters indexed by atom-type
  !            grp1out -- group 1 forcefield output (total only associated)
  !            grp2out -- group 2 forcefield output (total only associated)
  !            atypes1 -- array of atom types for molecule (group) 1
  !            atVecs1 -- array of coordinate vectors for molecule (group) 1
  !            atype2 -- array of atom types for molecule (group) 2
  !            atVecs2 -- array of coordinate vectors for molecule (group) 2
  !            simcell -- simulation cell data structure
  !----------------------------------------------------------------------------
  Logical Function pairmodel_msint5(pparams,grp1out,grp2out,atypes1,atVecs1, &
      atypes2,atVecs2,simcell) 
    Type(Pairmodel_Params), Dimension(:,:), Intent(In) :: pparams
    Type(EnergyPlus), Intent(InOut)                    :: grp1out,grp2out
    Type(VecType), Dimension(:), Intent(In)            :: atVecs1,atVecs2
    Integer, Dimension(:), Intent(In)                  :: atypes1,atypes2
    Type(SimCell_Params), Intent(In)                   :: simcell  

    Integer, Dimension(2,Size(atVecs2,1)) :: apairs

    Integer    :: a1,a2,ncoords1, ncoords2, natypes, natoms1,  natoms2
    Type(VecType), Dimension(Size(atVecs1,1),Size(atVecs2,1)) :: sepvec

    !** set defaults
    pairmodel_msint5 = .True.

    !** Get some sizes
    natypes = atom_getntypes()
    ncoords1 = Size(atVecs1, 1)
    ncoords2 = Size(atVecs2, 1)
    natoms2  = Size(atypes2, 1)
    natoms1  = 1

    !** Generate the type information for each element
    apairs(2, 1:ncoords2) = Reshape(aTypes2, (/ncoords2/), atypes2)

    !** Get the minimum image distance between all points in atVec1 and
    !** all points in atVec2, sepvec = atvec1 - atvec2
    Call simcell_minimage(simcell,atVecs1,atVecs2,sepvec)

    !** Do the evaluations 
    Do a1 = 1,natoms1
      Do a2 = 1,ncoords2
        If (Associated(pparams(atypes1(a1),apairs(2,a2))%lj)) Then
          pairmodel_msint5 = &
              lj_snglint2( pparams(atypes1(a1),apairs(2,a2))%lj, &
              sepvec(a1,a2),grp1out )
          If (.Not. pairmodel_msint5) Return
        End If
      End Do
    End Do

  End Function pairmodel_msint5

  !------------------------------------------------------------------------
  ! Checks the cutoffs for each initialized potential type.  The cutoff
  ! must be less than half the minimum width of the simulation cell or
  ! the minimum image convention algorithm is not guaranteed to actually 
  ! give the minimum image.
  ! Requires:  pparams -- spc-spc pairwise forcefield parameters
  !            simcell -- simulation cell information
  !------------------------------------------------------------------------
  Subroutine pairmodel_chkcutoffs(pparams,simcell)
    Type(Pairmodel_Params), Intent(In)     :: pparams    
    Type(SimCell_Params), Intent(In)       :: simcell

    Logical                  :: problem
    Real(kind=RDbl)          :: maxcutoff
    Character(len=strLen)    :: name1,name2
    Character(len=lstrLen)   :: offenders

    !** Initialize
    problem = .False.
    offenders = ''
    maxcutoff = simcell_minwidth(simcell)/2.0_RDbl

    If (Associated(pparams%lj)) Then
      If (lj_getcutoff(pparams%lj) > maxcutoff) Then
        problem = .True.
        Write(offenders,'(1x,a)') lj_idstring
      End If
    End If

    If (Associated(pparams%lj206)) Then
      If (lj206_getcutoff(pparams%lj206) > maxcutoff) Then
        problem = .True.
        Write(offenders,'(1x,a)') lj206_idstring
      End If
    End If

    If (Associated(pparams%coul)) Then
      If (coulomb_getcutoff(pparams%coul) > maxcutoff) Then
        problem = .True.
        Write(offenders,'(1x,a)') coulomb_idstring
      End If
    End If

    If (Associated(pparams%buck)) Then
      If (buck_getcutoff(pparams%buck) > maxcutoff) Then
        problem = .True.
        Write(offenders,'(1x,a)') buck_idstring
      End If
    End If

    If (Associated(pparams%cd)) Then
      If (benzpot_getcutoff(pparams%cd) > maxcutoff) Then
        problem = .True.
        Write(offenders,'(1x,a)') benzpot_idstring
      End If
    End If

    If (problem) Then
      name1 = atom_getname(pparams%atype1)
      name2 = atom_getname(pparams%atype2)

      Write(0,'(2a,i4,4a)') __FILE__,": ",__LINE__, &
          ' Problem with pair potentials for: ',Trim(name1),' with ',Trim(name2)
      Write(0,'(a,f8.4)') '  Maximum cut-off distance must be less than half '
      Write(0,'(a,f8.4)') '  the minimum width of the simulation cell:',maxcutoff
      Write(0,'(2a)') '  Offending potential types are: ',Trim(offenders)
      Stop
    End If

  End Subroutine pairmodel_chkcutoffs

  !----------------------------------------------------------------------------
  ! Displays the pair_params structure
  ! Requires: pparams -- parameter set
  !           indent -- number of spaces to indent
  !           optunit -- optional unit to write into
  !----------------------------------------------------------------------------
  Subroutine pairmodel_display(pparams,indent,optunit)
    Type(Pairmodel_Params), Intent(In) :: pparams
    Integer, Intent(In)                :: indent
    Integer, Intent(In), Optional      :: optunit

    Integer                            :: unitno
    Character(len=indent)              :: blank
    Character(len=xlstrLen)            :: string

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    If (pparams%nterms == 0) Then
      Write(unitno,'(a)') 'OFF'
      Return
    Else
      Write(unitno,'(a)') Trim(pparams%identifiers)
    End If

    If (Associated(pparams%lj)) Then
      Call lj_display(pparams%lj,indent+2,unitno)
    End If

    If (Associated(pparams%lj206)) Then
      Call lj206_display(pparams%lj206,indent+2,unitno)
    End If

    If (Associated(pparams%coul)) Then
      string = coulomb_disp(pparams%coul)
      Write(unitno,'(3x,2a,1x,a)') blank,Trim(string)
      string = coulomb_displayCutoffs(pparams%coul)
      Write(unitno,'(3x,2a,1x,a)') blank,"Cutoffs",Trim(string)
    End If

    If (Associated(pparams%buck)) Then
      Write(unitno,'(3x,2a,1x,a)') blank, &
          Trim(buck_displayParams(pparams%buck))
      Write(unitno,'(3x,2a,1x,a)') blank,"Cutoffs", &
          Trim(buck_displayCutoffs(pparams%buck))
    End If

    If (Associated(pparams%cd)) Then
      Call benzpot_display(pparams%cd,indent+2,unitno)
    End If

  End Subroutine pairmodel_display

  !----------------------------------------------------------------------------
  ! Make a 2-D array of pairmodel parameters symmetric.
  ! Requires:  pparams -- parameter set array
  !----------------------------------------------------------------------------
  Subroutine pairmodel_makesym(params)
    Type(Pairmodel_Params), Dimension(:,:), Intent(InOut) :: params

    Integer                   :: i,j,dim
    Integer                   :: nterms1,nterms2
    Logical                   :: equal
    Character(len=lstrLen)    :: string,problem

    dim = Size(params,1)
    If (Size(params,2) /= dim) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed array must be square'
      Stop
    End If

    !** Loop over the upper triangle and make sure it's same as lower triangle
    Do i = 1,dim
      Do j = i+1,dim
        nterms1 = params(i,j)%nterms
        nterms2 = params(j,i)%nterms
        If ((nterms1 > 0).And.(nterms2 > 0)) Then
          If (.Not. pairmodel_chkequal(params(i,j),params(j,i),problem)) Then
            Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
                ' Both entries in symmetric pair initialized, what to do?'
            Stop          
          End If
        End If

        If ((nterms1 == 0).And.(nterms2 > 0)) Then
          Call pairmodel_link(params(i,j),params(j,i))
        End If

        If ((nterms1 > 0).And.(nterms2 == 0)) Then
          Call pairmodel_link(params(j,i),params(i,j))
        End If

        !** Just make sure that it worked
        If (.Not. pairmodel_chkequal(params(i,j),params(j,i),problem)) Then
          Write(0,'(2a,i4,a,2i3)') __FILE__,": ",__LINE__, &
              ' Problem with symmetrizing pair ',i,j
          Write(0,'(4x,a)') Trim(problem)
          Stop          
        End If

      End Do
    End Do

  End Subroutine pairmodel_makesym

  !----------------------------------------------------------------------------
  ! Checks to make sure a 2-D array of pairmodel parameters is symmetric.
  ! Requires:  pparams -- parameter set array
  !----------------------------------------------------------------------------
  Subroutine pairmodel_chkarray(pparams)
    Type(Pairmodel_Params), Dimension(:,:), Intent(In) :: pparams

    Integer                   :: i,j,dim
    Logical                   :: equal
    Character(len=lstrLen)    :: string,problem

    dim = Size(pparams,1)
    If (Size(pparams,2) /= dim) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Passed array must be square'
      Stop
    End If

    !** Loop over the upper triangle and make sure it's same as lower triangle
    Do i = 1,dim
      Do j = i+1,dim

        equal = pairmodel_chkequal(pparams(i,j),pparams(j,i),problem)
        If (.Not. equal) Then
          Write(0,'(2a,i4,a,2i3,2a)') __FILE__,": ",__LINE__, &
              ' Problem with atype pair ',i,j,': ',Trim(problem)
          Stop
        End If

      End Do
    End Do

  End Subroutine pairmodel_chkarray

  !----------------------------------------------------------------------------
  ! Checks a pair of pairmodel parameters sets to see if they are identical
  ! Requires:  params1 -- 1st parameter set
  !            params2 -- 2nd parameter set
  !            problem -- statement of problem
  !            tolerance -- optional specified tolerance
  !----------------------------------------------------------------------------
  Logical Function pairmodel_chkequal(params1,params2,problem,tolerance)
    Type(Pairmodel_Params), Intent(In)     :: params1,params2
    Character(*), Intent(Out)              :: problem
    Real(kind=RDbl), Intent(In), Optional  :: tolerance

    Logical          :: equal
    Real(kind=RDbl)  :: tol

    tol = 1.0e-6_RDbl
    If (Present(tolerance)) tol = tolerance

    pairmodel_chkequal = .True.
    problem = 'none'
    
    If (params1%nterms /= params1%nterms) Then
      pairmodel_chkequal = .False.
      problem = 'nterms not equal'
      Return
    End If

    !** Check LJ term
    If (Associated(params1%lj)) Then
      problem = 'second LJ not associated'
      If (.Not. Associated(params2%lj)) pairmodel_chkequal = .False.
      problem = 'LJ term not equal'
      equal = lj_chkequal(params1%lj,params2%lj,tol)
      If (.Not. equal) pairmodel_chkequal = .False.
    End If

    !** Check LJ 20-6 term    
    If (Associated(params1%lj206)) Then
      problem = 'second LJ206 not associated'      
      If (.Not. Associated(params2%lj206)) pairmodel_chkequal = .False.
      problem = 'LJ206 term not equal'
      equal = lj206_chkequal(params1%lj206,params2%lj206,tol)
      If (.Not. equal) pairmodel_chkequal = .False.
    End If

    !** Check Coulombic term    
    If (Associated(params1%coul)) Then
      problem = 'second coulombic term not associated'      
      If (.Not. Associated(params2%coul)) pairmodel_chkequal = .False.
      problem = 'Coulombic term not equal'
      equal = coul_chkequal(params1%coul,params2%coul,tol)
      If (.Not. equal) pairmodel_chkequal = .False.
    End If

    !** Check Buckingham term    
    If (Associated(params1%buck)) Then
      problem = 'second buckingham term not associated'      
      If (.Not. Associated(params2%buck)) pairmodel_chkequal = .False.
      problem = 'Buckingham term not equal'
      equal = buck_chkequal(params1%buck,params2%buck,tol)
      If (.Not. equal) pairmodel_chkequal = .False.
    End If

    !** Check CD term    
    If (Associated(params1%cd)) Then
      problem = 'second CD term not associated'      
      If (.Not. Associated(params2%cd)) pairmodel_chkequal = .False.
      problem = 'Buckingham term not equal'
      equal = benzpot_chkequal(params1%cd,params2%cd,tol)
      If (.Not. equal) pairmodel_chkequal = .False.
    End If

    problem = 'none'

  End Function pairmodel_chkequal

  !----------------------------------------------------------------------------
  ! Cleans the pair_params structure
  ! Requires:  pparams -- pairmodel parameters
  !----------------------------------------------------------------------------
  Subroutine pairmodel_clean(pparams)
    Type(Pairmodel_Params), Intent(InOut)     :: pparams

    Integer             :: error
    
    If (Associated(pparams%lj)) Then
      Deallocate(pparams%lj,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'lj')
    End If

    If (Associated(pparams%lj206)) Then
      Deallocate(pparams%lj206,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'lj206')
    End If

    If (Associated(pparams%coul)) Then
      Deallocate(pparams%coul,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'coul')
    End If

    If (Associated(pparams%buck)) Then
      Deallocate(pparams%buck,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'buck')
    End If

    If (Associated(pparams%cd)) Then
      Deallocate(pparams%cd,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'cd')
    End If
    
  End Subroutine pairmodel_clean

!------------------------------------------------------------------------
! Tina added
! Get pairpotential parameters
! pot_params(1) = A
! pot_params(2) = B
! pot_params(3) = C
! pot_params(4) = D
! pot_params(5) = hicut
! pot_params(6) = locut
!-------------------------------------------------------------------------

 Subroutine pairmodel_getpotparameters(params,pot_params)
 
  Type(pairmodel_params), INTENT(IN)           :: params 
  Real(kind = RDbl), Dimension(6), Intent(OUT) :: pot_params
  Real(kind = RDbl), Dimension(4)              :: params_abcd
  Real(kind = RDbl), Dimension(2)              :: params_cut
  Integer                                      :: i

  pot_params = 0._RDbl

  IF(cleanstring(params%identifiers) == 'LJ') THEN
    CALL lj_getABCD(params%lj,params_abcd)
    Do i = 1, 4
       pot_params(i) = params_abcd(i)
    END DO
    CALL lj_getcut(params%lj,params_cut)
    Do i = 5,6
       pot_params(i) = params_cut(i-4)
    END Do
  ELSE
    Write(0,*) __FILE__,': ',__LINE__
    write(0,*) 'Implement subroutines to get potential parameters for potential ',&
                cleanstring(params%identifiers),' before you go on'
    STOP
  END IF

 END Subroutine pairmodel_getpotparameters

End Module pairmodel



