!------------------------------------------------------------------------------
! This module contains the objects and methods for dealing with intra-pair
! interactions.  Like other potential model data structures, it contains a 
! list of pointers to specific intra-pair potential types.  The pointer
! that is initialized directs the code to do that type of calculation.
!
! Needed Improvements:
! 1) remove usage of molecule pointer, this is hack
! 2) reincorporate _getinteractionMult into call tree to speed things
!------------------------------------------------------------------------------

Module ipmodel

  Use defaults, Only: RDbl,strLen,lstrLen,INTRAPAIR_INDEX,degTorad,MAX_ATOMS
  Use utils, Only: toupper,combine,deallocerrdisplay,allocerrdisplay, &
      filesrchstr,split,toint,stripcmnt,findint
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/),vector_getnorm,mag
  Use lj, Only: LJ_Params, lj_displayparams, lj_init, lj_calc_interaction, &
      lj_idstring, lj_isinit, lj_displayCutoffs, lj_getcutoff, Assignment(=)
  Use buck, Only: Buckingham_Params, buck_init, buck_calc_interaction, &
      buck_idstring, buck_isinit, buck_getcutoff, buck_displayCutoffs, &
      buck_displayParams, Assignment(=)
  Use atom, Only: atom_gettypename,atom_getntypes,atom_getname
  Use molecule, Only: MolecularParams, & !HACK!
      molecule_getnatomtypes,molecule_getatomcoords
  Use molecules, Only: molecules_getnatomtypes, molecules_getnsorbs
  Use connects, Only: connects_makeTypeList,connects_allneighbors
  Use storebase, Only: storebase_inc
  Use store, Only: Store_Level_Pair,store_nderivs,store_idbranch
  Use simcell, Only: SimCell_Params, simcell_minimage

  Implicit None

  Private
  Public :: IntrapairModel, IntrapairInfo, ipmodel_ipinit, ipmodel_ipint, &
      ipmodel_ipcleanup, ipmodel_isinit, ipmodel_getcutoff, ipmodel_copy, &
      ipmodel_ipdisplay, ipmodel_ipinfo, ipmodel_displaycutoffs, isfast, &
      ipmodel_eval, INTRAPAIR_KEY

  Character(len=strLen), Parameter  :: INTRAPAIR_KEY = 'INTRAPAIR'

  Type IntrapairModel
    Type(LJ_Params), Pointer         :: lj
    Type(Buckingham_Params), Pointer :: buck
  End Type IntrapairModel

  Type IntrapairInfo
    Integer                                   :: nipairs, naway
    Logical                                   :: fast
    Character(len=strLen)                     :: model
    Integer, Dimension(:,:), Pointer          :: iplist
    Type(IntrapairModel), Dimension(:), Pointer   :: ipparams
    Type(IntrapairModel), Dimension(:,:), Pointer :: atype_list
  End Type IntrapairInfo

  Interface ipmodel_getinteraction
    Module Procedure ipmodel_getinteractionSngl
    Module Procedure ipmodel_getinteractionMult
  End Interface

  Interface isfast
    Module procedure ipmodel_isfastip
  End Interface

Contains
  !----------------------------------------------------------------------
  ! Initialize the intra pair parameters, assumes the input is 
  ! already allocated.
  ! Requires: params -- intra-pair pointer set, to be initialized
  !           modelName -- character string defining model type
  !           specs_line -- string containing model parameters
  ! If only the IntrapairModel variable is passed, it only 
  ! nullifies the two pointers and then exits.
  !----------------------------------------------------------------------
  Subroutine ipmodel_initparams(params,modelName,specs_line)
    Type(IntrapairModel), Intent(InOut)     :: params
    Character(*), Intent(In), Optional      :: modelName,specs_line

    Integer                      :: error

    !** Nullify all the pointers first
    Nullify(params%lj)
    Nullify(params%buck)

    !** Check if the modelName is present. If not, return.
    If (.Not.Present(modelName)) Return

    !** Figure out which model to initialize
    Select Case (Toupper(trim(modelName)))
    Case (lj_idstring)
      Allocate(params%lj,STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'lj')

      Call lj_init(params%lj,specs_line)

    Case (buck_idstring)
      Allocate(params%buck,STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'buck')

      Call buck_init(params%buck,specs_line)

    Case Default
      Write(0,*) 
      Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
          " : No such bond stretch model exists"

    End Select

  End Subroutine ipmodel_initparams

  !----------------------------------------------------------------------------
  ! Initialize intra-pair (1-4) interactions
  ! Requires: ptr -- pointer to intra-pair information 
  !           molecule -- pointer to molecule data structure
  !           filename -- file where intra-pair information can be found
  !----------------------------------------------------------------------------
  Recursive Subroutine ipmodel_ipinit(ptr,molec,filename,simcell)
    Type(IntrapairInfo), Pointer               :: ptr
    Type(MolecularParams), Pointer             :: molec
    Character(*), Intent(In)                   :: filename
    Type(Simcell_Params), Intent(In), Optional :: simcell

    Integer                          :: error,i,j,k,m
    Integer                          :: atom1,atom2,ninteract,unitno
    Integer                          :: atype1,atype2,nipairs
    Integer                          :: useIP,naway,natypes,n
    Real(Kind=RDbl)                  :: cut, dist, dist2
    Type(VecType)                    :: sepvec
    Character(len=strLen)            :: modelname,atom_name1,atom_name2
    Character(len=lstrLen)           :: key, specs_line
    Character(len=255)               :: text
    Character(len=strLen), Dimension(strLen) :: params
    Integer, Dimension(molec%natoms)         :: alist, atypes
    Integer, Dimension(Int(0.5*molec%natoms*(molec%natoms-1)),2) :: conlist

    If(.Not.Associated(ptr)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          " : passed pointer must associated"
      Stop
    End If

    !** nullify the atom-based parameters list, we may not use it
    Nullify(ptr%atype_list)

    !** Open the file containing the intramolecular interaction information
    unitno = file_open(filename,110)

    key = "Intrapair_Info"

    !** Read the type of interation for the molecule
    useip = filesrchstr(unitno,key,text,.True.)

    If (useip == 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "No intrapair potentials for this molecule- Something is wrong"
      Stop
    End If

    j = split(text,params)
    ptr%model = trim(params(3))
    modelname = trim(params(3))
    !** Fast or slow?
    If (Trim(toupper(params(4))) == "FAST") Then
      ptr%fast = .True.
    Else If (trim(toupper(params(4))) == "SLOW") Then
      ptr%fast = .False.
    Else
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          " : Must indicate 'Fast' or 'Slow' for ", key
      Stop
    End If
    
    !** Look for 'LISTED', 'GENERATED' or another filename    
    If (toupper(Trim(params(2))) == 'LISTED') Then
      !** intra-pair parameters are listed, get number
      Read(unitno,'(a)') text
      j = split(text,params)
      nipairs = toint(params(1))
      ptr%nipairs = nipairs

      !** initialize the model information pointer
      Allocate(ptr%ipparams(ptr%nipairs),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'ipparams')
    
      !** initialize the intra-pair list
      Allocate(ptr%iplist(nipairs,2),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'iplist')

      Do i = 1,nipairs
        Read(unitno,'(a)') text
        j = split(text,params)
        atom1 = toint(params(1))
        atom2 = toint(params(2))
        ptr%iplist(i,:) = (/atom1, atom2/)
        atom_name1 = atom_getname(molec%atoms(atom1)%atype)
        atom_name2 = atom_getname(molec%atoms(atom2)%atype)
        Write(specs_line,'(2(a,2x),a)') Trim(atom_name1),Trim(atom_name2), &
            Trim(combine(params(3:j)))
        Call ipmodel_initparams(ptr%ipparams(i),ptr%model,specs_line)
      End Do

    Else If (toupper(Trim(params(2))) == 'GENERATE') Then

      !** Zero the number of pairs
      ptr%nipairs = 0
    
      !** Read in the minimum separation distance between the atom pairs
      Read(unitno,'(a)') text
      j = split(text,params)
      naway = toint(params(1))
      ptr%naway = naway
    
      !** Warn for untested systems
      If ((naway /= 4).Or.(naway /= 1)) Then
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            " WARNING: Has not been verified for other than 1 or 4 away"
      End If
    
      !** Get the number of types to generate
      Read(unitno,'(a)') text
      j = split(text,params)
      ninteract = toint(params(1))

      !** Get the number of atom types and allocate atype-based storage
      
      ! PREVIOUS CODE - WRONG, this gives type of atoms in this molecule
      !      natypes = molecule_getnatomtypes(molec,atypes)

      !** rough estmate of total number of atom types
      natypes=0
      Do i=1,molecules_getnsorbs()
        natypes=natypes+molecules_getnatomtypes(i)
      End Do
      Allocate(ptr%atype_list(natypes,natypes),STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__)

      !** Loop through and nullify the individual parameters
      Do i = 1, natypes
        Do j = 1, natypes
          Call ipmodel_initparams(ptr%atype_list(i,j))
        End Do
      End Do

      ! for further calculations get natypes of this molecule
      natypes = molecule_getnatomtypes(molec,atypes)

      !** Read in each of the atype-based interaction types
      Do i = 1,ninteract
        
        !** Read in the line containing the atom pair and parameters
        Read(unitno,'(a)') text
        j = split(text,params)

        !** Get the atom types from the names
        atype1 = atom_gettypename(params(1))
        atype2 = atom_gettypename(params(2))

        !** Set up the parameters
        Call ipmodel_initparams(ptr%atype_list(atype1,atype2),ptr%model,text)
        Call ipmodel_initparams(ptr%atype_list(atype2,atype1),ptr%model,text)
      End Do

      !** Make the neighbor lists of interactions based on the types
      m = 0
      Do j = 1,(molec%natoms - 1)

        If (naway-1 > 0) Then 
          !** Get the list of all neighbors that are less than naway
          n = connects_allneighbors(molec%connect,j,alist,naway-1)
        Else
          n = 0
        End If

        Do k = j+1, molec%natoms
          atype1 = molec%atoms(j)%atype
          atype2 = molec%atoms(k)%atype
          !** Place in list if it's a valid pair of atoms
          If (ipmodel_isinit(ptr%atype_list(atype1,atype2))) Then

            !** Check to see if the atom is too close
            If ((findint(alist(1:n),k) == 0).Or.(n == 0)) Then

              !** Finally, check to see if the atom is too far away
              !** this looks for some special cutoff, not the normal one
              cut = ipmodel_getcutoff(ptr%atype_list(atype1,atype2),.True.)
              !** Get the minimum image distance if simcell is present
              If (Present(simcell)) Then
                Call simcell_minimage(simcell,molec%atoms(j)%r, &
                    molec%atoms(k)%r,sepvec,dist2)
                dist = Sqrt(dist2)
              Else
                dist = mag(molec%atoms(j)%r - molec%atoms(k)%r) 
              End If
              If (dist < cut) Then
                m = m+1
                conlist(m,:) = (/j,k/)
              End If

            End If
          End If
        End Do
      End Do

      !** Set the number of pairs
      ptr%nipairs = m

      !** Allocate the list
      Allocate(ptr%iplist(ptr%nipairs,2), STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'iplist')

      !** Copy the pairs into the structure
      ptr%iplist = conlist(1:ptr%nipairs,:)

      !** allocate memory for the individual intra-pair parameters
      Allocate(ptr%ipparams(ptr%nipairs), STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'ipparams')

      !** initialize the intra-pair parameters for each
      Do i = 1,ptr%nipairs
        atype1 = molec%atoms(ptr%iplist(i,1))%atype
        atype2 = molec%atoms(ptr%iplist(i,2))%atype
        Call ipmodel_copy(ptr%ipparams(i),ptr%atype_list(atype1,atype2))
      End Do
    
    Else
      !** assume it's a file we are reading in.
      Close(unitno)
      Call ipmodel_ipinit(ptr,molec,Trim(params(2)))
    End If

  End Subroutine ipmodel_ipinit

  !-------------------------------------------------------------------
  ! Initialize pointer set by copying from an existing set
  ! Requires: newset -- new intra-pair pointer set, to be initialized
  !           oldset -- old intra-pair pointer set
  !-------------------------------------------------------------------
  Subroutine ipmodel_copy(newset,oldset)
    Type(IntrapairModel), Intent(Out)     :: newset
    Type(IntrapairModel), Intent(In)      :: oldset

    Integer           :: error
    
    Nullify(newset%lj)
    Nullify(newset%buck)

    If (Associated(oldset%lj)) Then
      Allocate(newset%lj,STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'lj')
      newset%lj = oldset%lj
    Else If (Associated(oldset%buck)) Then
      Allocate(newset%buck,STAT=error)
      If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__,'buck')
      newset%buck = oldset%buck
    End If    

  End Subroutine ipmodel_copy

  !----------------------------------------------------------------------------
  ! Returns True if the interaction is initialized, else False
  ! Requires: params -- intra-pair pointer set
  !----------------------------------------------------------------------------
  Logical Function ipmodel_isinit(params)
    Type(IntrapairModel), Intent(In) :: params

    !** Set the flag to false
    ipmodel_isinit = .False.
    If (Associated(params%lj)) Then
      ipmodel_isinit = lj_isinit(params%lj)
    Else If (Associated(params%buck)) Then
      ipmodel_isinit = buck_isinit(params%buck)
    End If
  End Function ipmodel_isinit

  !----------------------------------------------------------------------------
  ! Returns the value of the cutoff. May return the neighbor cutoff if
  ! NEIGHBOR is present and True.
  ! Requires: params -- intra-pair pointer set
  !           neighbor -- flag to get neighbor cutoff
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function ipmodel_getcutoff(params,neighbor)
    Type(IntrapairModel), Intent(In) :: params
    Logical, Optional, Intent(In)    :: neighbor

    !** Pick the correct interaction
    If (Associated(params%buck)) Then
      If (Present(neighbor)) Then
        ipmodel_getcutoff = buck_getcutoff(params%buck,neighbor)
      Else
        ipmodel_getcutoff = buck_getcutoff(params%buck)
      End If
    End If

    If (Associated(params%lj)) Then
      If (Present(neighbor)) Then
        ipmodel_getcutoff = lj_getcutoff(params%lj,neighbor)
      Else
        ipmodel_getcutoff = lj_getcutoff(params%lj)
      End If
    End If

  End Function ipmodel_getcutoff

  !----------------------------------------------------------------------------
  ! Calculate the intrapair potential for atom pair atype1, atype2 and
  ! return the potential u and (optionally) the force vectors vf
  ! Requires: params -- intra-pair pointer set
  !           sepvec -- separation vector for a single pair
  !           u -- calculated potential energy
  !           vf -- calculated force, optional
  !----------------------------------------------------------------------------
  Subroutine ipmodel_getinteractionSngl(params,sepvec,u,vf)
    Type(IntrapairModel), Intent(In)      :: params
    Type(VecType), Intent(In)             :: sepvec
    Real(kind=RDbl), Intent(OUT)          :: u
    Type(VecType), Intent(Out), Optional  :: vf

    Logical       :: loFlag   

    !** check to see which model is associated for this given pair
    If (Associated(params%lj)) Then
      If (Present(vf)) Then
        Call lj_calc_interaction(params%lj,sepvec,u,loFlag,vf)
      Else
        Call lj_calc_interaction(params%lj,sepvec,u,loFlag)
      Endif
      
      If (.Not. loFlag) Then 
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' The distance is too close, LJ calculation failed '
        Stop
      End If
    End If

    If (Associated(params%buck)) Then
      loFlag = buck_calc_interaction(params%buck,sepvec,u,vf)
      If (.Not.loFlag) Then
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            ' The distance is too close, Buckingham calculation failed '
        Stop
      End If
    End If

  End Subroutine ipmodel_getinteractionSngl

  !----------------------------------------------------------------------------
  ! Choose the model and get the interaction. Works on a dimensioned
  ! list of separation vectors
  ! Requires: params -- array of intra-pair pointer set, one for each pair
  !           sepvec -- separation vectors
  !           u -- calculated potential energy
  !           vf -- calculated force, optional
  !----------------------------------------------------------------------------
  Subroutine ipmodel_getinteractionMult(params,sepvec,u,vf)
    Type(IntrapairModel), Dimension(:), Intent(In)     :: params
    Type(VecType), Dimension(:), Intent(In)            :: sepvec
    Real(kind=RDbl), Intent(Out)                       :: u
    Type(VecType), Dimension(:), Intent(Out), Optional :: vf

    Logical          :: loFlag   
    Integer          :: ip
    Real(kind=RDbl)  :: pot

    u = 0.0_RDbl

    !** Loop through the pairs and call the appropriate model
    Do ip = 1,Size(params)

      If (Associated(params(ip)%lj)) Then
        Call lj_calc_interaction(params(ip)%lj,sepvec(ip),pot,loFlag,vf(ip))
        If( .Not. loFlag) Then 
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              ' The distance is too close '
          Stop
        End If
      End If
      u = u + pot

      If (Associated(params(ip)%buck)) Then
        loFlag = buck_calc_interaction(params(ip)%buck,sepvec(ip),pot,vf(ip))
        If (.Not.loFlag) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              ' The distance is too close '
          Stop
        End If
      End If
      u = u + pot

    End Do

  End Subroutine ipmodel_getinteractionMult

  !----------------------------------------------------------------------------
  ! Calculate the intra-pair interaction using actual coordinates of atoms.
  ! Periodic boundary conditions will be applied to the atomic coordinates if
  ! the simcell is passed.
  ! Requires: params -- intra-pair forcefield parameters
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Subroutine ipmodel_ipint(params,ffout,coords,simcell)
    Type(IntrapairInfo), Intent(In)            :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell

    Integer                      :: i,nderivs,level
    Integer                      :: atom1, atom2  !atom numbers in molecule
    Real(kind=RDbl)              :: u,upart
    Type(VecType)                :: sepvec,force

    !** Get the type of the passed storage
    level = store_idbranch(ffout)

    !** Loop through the pairs
    nderivs = store_nderivs(ffout)
    Do i = 1, params%nipairs
      atom1 = params%iplist(i,1)
      atom2 = params%iplist(i,2)

      !** Get the separation vector
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,coords(atom1),coords(atom2),sepvec)
      Else
        sepvec = coords(atom1) - coords(atom2)
      End If



      !** Do the evaluation and store the results
      If (nderivs == 0) Then
        Call ipmodel_getinteraction(params%ipparams(i),sepvec,u)
      Else If (nderivs == 1) Then
        Call ipmodel_getinteraction(params%ipparams(i),sepvec,u,force)

        !** Store the results according to the storage form
        Select Case(level)
        Case(0)
          Call storebase_inc(ffout%binfo(atom1),force)
          Call storebase_inc(ffout%binfo(atom2),(force*(-1.0_RDbl)))
        Case(10)
          Call storebase_inc(ffout%mi(atom1)%total,force)
          Call storebase_inc(ffout%mi(atom2)%total,(force*(-1.0_RDbl)))
        Case Default
          Write(0,'(2a,i5,a,i3)') __FILE__,":",__LINE__, &
              ' Unable to handle storage type ',level
          Stop        
        End Select
            
      Else
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            ' Second order terms not presently available'
        Stop        
      End If

      !** Store the potential energy results according to the storage form
      Select Case(level)
      Case(-1)
        Call storebase_inc(ffout%total,u,INTRAPAIR_INDEX)
      Case(0)
        upart = 0.5_RDbl*u
        Call storebase_inc(ffout%binfo(atom1),upart,INTRAPAIR_INDEX)
        Call storebase_inc(ffout%binfo(atom2),upart,INTRAPAIR_INDEX)
      Case(10)
        upart = 0.5_RDbl*u
        Call storebase_inc(ffout%mi(atom1)%total,upart,INTRAPAIR_INDEX)
        Call storebase_inc(ffout%mi(atom2)%total,upart,INTRAPAIR_INDEX)
      Case Default
        Write(0,'(2a,i5,a,i3)') __FILE__,":",__LINE__, &
            ' Unable to handle storage type ',level
        Stop        
      End Select
    End Do

  End Subroutine ipmodel_ipint

#ifdef PREVIOUS
  !----------------------------------------------------------------------------
  ! Calculate the intra-pair interaction using actual coordinates of atoms.
  ! Periodic boundary conditions will be applied to the atomic coordinates if
  ! the simcell is passed.
  ! Requires: params -- intra-pair forcefield parameters
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Subroutine ipmodel_ipint(params,ffout,coords,simcell)
    Type(IntrapairInfo), Intent(In)            :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell

    Integer                      :: i,nderivs
    Integer                      :: atom1, atom2  !atom numbers in molecule
    Real(kind=RDbl)              :: u,upart
    Type(VecType)                :: sepvec,force

    !** Loop through the pairs
    nderivs = store_nderivs(ffout)
    Do i = 1, params%nipairs
      atom1 = params%iplist(i,1)
      atom2 = params%iplist(i,2)

      !** Get the separation vector
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,coords(atom1),coords(atom2),sepvec)
      Else
        sepvec = coords(atom1) - coords(atom2)
      End If

      !** Do the evaluation and store the results
      If (nderivs == 0) Then
        Call ipmodel_getinteraction(params%ipparams(i),sepvec,u)
      Else If (nderivs == 1) Then
        Call ipmodel_getinteraction(params%ipparams(i),sepvec,u,force)
            
        Call storebase_inc(ffout%mi(atom1)%total,force)
        Call storebase_inc(ffout%mi(atom2)%total,(force*(-1.0_RDbl)))
      Else
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            ' Second order terms not presently available'
        Stop        
      End If

      upart = 0.5_RDbl*u
      Call storebase_inc(ffout%mi(atom1)%total,upart,INTRAPAIR_INDEX)
      Call storebase_inc(ffout%mi(atom2)%total,upart,INTRAPAIR_INDEX)
    End Do

  End Subroutine ipmodel_ipint
#endif

  !----------------------------------------------
  ! Checks if the IP interaction is fast or slow
  !----------------------------------------------
  Logical Function ipmodel_isfastip(intrapair)
    Type(IntrapairInfo) :: intrapair

    ipmodel_isfastip = intrapair%fast
  End Function ipmodel_isfastip

  !--------------------------------------------------------------------
  ! Determines if the interaction should be evaluated
  ! Requires:  ptr -- pointer to interaction parameters structure
  !            fast -- proposed evaluation speed
  !--------------------------------------------------------------------
  Logical Function ipmodel_eval(ptr,fast)
    Type(IntrapairInfo), Pointer    :: ptr
    Logical, Intent(In)             :: fast

    ipmodel_eval = .False.
    If (.Not. Associated(ptr)) Return
    ipmodel_eval = (ptr%fast .And. fast)

  End Function ipmodel_eval

  !--------------------------------------------------------
  ! Display the model types and parameters
  ! Requires: params -- intra-pair pointer set
  !--------------------------------------------------------
!!$  Character(len=lstrLen) Function ipmodel_display(params)
  Function ipmodel_display(params)
    Character(len=lstrLen)              :: ipmodel_display
    Type(IntrapairModel), Intent(In)    :: params
    Character(len=lstrLen)              :: string

    !** LJ values
    If (Associated(params%lj)) Then
      string=lj_displayparams(params%lj)
      Write(ipmodel_display,'(a,2x,a)') "LJ", Trim(string)
    End If

    !** Buck values
    If (Associated(params%buck)) Then
      string=buck_displayParams(params%buck)
      Write(ipmodel_display,'(a,2x,a)') "Buck", Trim(string)
    End If

  End Function ipmodel_display

  !--------------------------------------------------------
  ! Display the model cutoffs
  ! Requires: params -- intra-pair pointer set
  !--------------------------------------------------------
!!$  Character(len=lstrLen) Function ipmodel_displayCutoffs(params)
  Function ipmodel_displayCutoffs(params)
    Character(len=lstrLen)              :: ipmodel_displayCutoffs
    Type(IntrapairModel), Intent(In)    :: params
    Character(len=lstrLen)              :: string

    !** LJ values
    If (Associated(params%lj)) Then
      string=lj_displayCutoffs(params%lj)
      Write(ipmodel_displayCutoffs,'(a)') Trim(string)
    End If

    !** Buck values
    If (Associated(params%buck)) Then
      string=buck_displayCutoffs(params%buck)
      Write(ipmodel_displayCutoffs,'(a)') Trim(string)
    End If

  End Function ipmodel_displayCutoffs

  !----------------------------------------------------------------------------
  ! Returns all tabulated bond-stretching information
  ! Requires:  info -- stretching data structure for a single species
  !            nsets -- number of potential sets returned
  !            list -- list of atom number pairs in molecule  (pair_index,1:2)
  !            params -- array of strings containing type and parameters
  !----------------------------------------------------------------------------
  Subroutine ipmodel_ipinfo(info,nsets,list,params)
    Type(IntraPairInfo), Intent(In)                     :: info
    Integer, Intent(Out)                                :: nsets
    Integer, Dimension(:,:), Intent(Out)                :: list
    Character(len=lstrLen), Dimension(:), Intent(Out)   :: params

    Integer          :: i

    !** copy the list
    nsets = info%nipairs 
!LC    Write(*,*) nsets, Size(list,1), Size(list,2)

    list(1:nsets,1:2) = info%iplist(1:nsets,1:2)

    !** get the parameters
    Do i = 1,nsets
      params(i) = ipmodel_display(info%ipparams(i))
    End Do

  End Subroutine ipmodel_ipinfo

  !------------------------------------------------------------------------
  ! Display the information about intra-molecular pairwise interactions
  ! Requires: ipinfo -- intra-pair information structure
  !           spc -- species number
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !------------------------------------------------------------------------
  Subroutine ipmodel_ipdisplay(ipinfo,spc,indent,unit)
    Type(IntrapairInfo)            :: ipinfo
    Integer, Intent(In)            :: spc,indent,unit

    Integer                        :: i,j,natypes
    Character(len=indent)          :: blank
    Character(len=strLen)          :: name1,name2, string
    Integer, Dimension(MAX_ATOMS)  :: atyp
    blank = Repeat(' ',indent)

    !** atype-based list if it exists    
    If (Associated(ipinfo%atype_list)) Then
      natypes = molecules_getnatomtypes(spc,atyp)
      Do i = 1, natypes
        Do j = i, natypes
          name1 = atom_getname(atyp(i))
          name2 = atom_getname(atyp(j))
          string=ipmodel_display(ipinfo%atype_list(atyp(i),atyp(j)))
          If (ipmodel_isinit(ipinfo%atype_list(atyp(i),atyp(j)))) Then
            Write(unit,'(a,3(a,2x))') blank,Trim(name1),Trim(name2),&
                Trim(string)
            string=ipmodel_displayCutoffs(ipinfo%atype_list(atyp(i),atyp(j)))
            Write(unit,'(4x,3a)') blank,"Cutoffs : ",&
                Trim(string)
          End If

        End Do
      End Do
    End If

    !** individual pairs and their parameters
    If (ipinfo%nipairs < 20) Then
      Write(unit,'(2a)') blank,'List of individual pairs (atom1,atom2,info)'
      Do i = 1,ipinfo%nipairs
        string=ipmodel_display(ipinfo%ipparams(i))
        Write(unit,'(2x,a,i2,a,2i4,4x,a)') blank,i,':',ipinfo%iplist(i,1), &
            ipinfo%iplist(i,2),Trim(string)
      End Do
    Else
      Write(unit,'(2x,3a)') blank,'number of intra-pairs > 20,', &
          ' not showing intra-pair doublets'
    End If
  End Subroutine Ipmodel_ipdisplay

  !--------------------------------------------------------
  ! Cleanup the intrapair pointer set
  ! Requires: params -- intra-pair pointer set
  !--------------------------------------------------------
  Subroutine ipmodel_cleanup(params)
    Type(IntrapairModel), Intent(InOut) :: params

    Integer          :: error
    
    If (Associated(params%lj)) Then
      Deallocate(params%lj,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'lj')
    End If

    If (Associated(params%buck)) Then
      Deallocate(params%buck,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'buck')
    End If

  End Subroutine ipmodel_cleanup

  !------------------------------------------
  ! Cleanup 1-4 info
  !------------------------------------------
  Subroutine ipmodel_ipcleanup(intrapair)
    Type(IntrapairInfo)            :: intrapair

    Integer       :: i,j,error

    Do i = 1,Size(intrapair%ipparams)
      Call ipmodel_cleanup(intrapair%ipparams(i))
    End Do

    Do i = 1,Size(intrapair%atype_list,1)
      Do j = 1,Size(intrapair%atype_list,2)
        Call ipmodel_cleanup(intrapair%atype_list(i,j))
      End Do
    End Do

    Deallocate(intrapair%iplist,stat=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'pairs')

  End Subroutine ipmodel_ipcleanup

End Module ipmodel














