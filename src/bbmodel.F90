!---------------------------------------------------------------------------
! The Bond Bending model module handles all bond bending interaction.  It
! serves to route calls to different bond bending potential types.
!---------------------------------------------------------------------------

Module bbmodel

  Use defaults, Only: RDbl,strLen,lstrLen,BENDING_INDEX,degTorad,xlstrLen
  Use utils, Only: toupper,combine,deallocerrdisplay,allocerrdisplay, &
      fileSrchwoComment,split,toint,stripcmnt,findint,str2seq,int2str, &
      condenseint
  Use file, Only: file_open
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/),vector_getnorm,vector_angle,vector_crossprod, &
      vector_getunitvec
  Use harmonica, Only: harmonica_init, harmonica_getinteraction, &
      harmonica_display, HarmonicAModel, harmonica_getparams,HarmonicaString
  Use atom, Only: atom_gettypename,atom_getntypes,atom_nbonds
  Use molecules, Only: molecules_name,molecules_getnatomtypes, &
      molecules_getnatoms,molecules_getconnects,molecules_getatype, &
      molecules_getcoords
  Use connects, Only: AtomConnections,connects_makeTypeList
  Use storebase, Only: storebase_inc
  Use store, Only: Store_Level_Pair,store_nderivs,store_idbranch
  Use simcell, Only: SimCell_Params, simcell_minimage

  Implicit None
  Save

  Private
  Public :: BondBendingModel, BendingInfo, bbmodel_bbinit, bbmodel_eval, &
      bbmodel_bbinfo, bbmodel_bbint, bbmodel_bbdisplay, isfast, &
      bbmodel_getparams, bbmodel_bbcleanup, BENDING_KEY, bbmodel_setHarK

  Character(len=strLen), Parameter  :: BENDING_KEY = 'BENDING'

  Type BondBendingModel
    Logical                       :: dummyatom
    Integer, Dimension(3)         :: baseatoms   !** used to construct dummyatom
    Type(HarmonicAModel), Pointer :: hara
  End Type BondBendingModel

  !** handles all 3-body interactions
  Type BendingInfo
    Integer                                :: nbends
    Logical                                :: fast
    Character(40)                          :: model
    Integer, Dimension(:,:), Pointer       :: bblist
    Logical, Dimension(:,:,:), Pointer     :: bbflag
    Type(BondBendingModel), Dimension(:), Pointer :: bbparams
  End Type BendingInfo

  Interface isfast
    Module procedure bbmodel_isfastbb
  End Interface

Contains

  !-------------------------------------------------------------------
  ! Initialize the bond-bending parameters, assumes the input is 
  ! already allocated.
  ! Requires:  bbparams -- bond bending type, to be initialized
  !            model -- character string defining model type
  !            angle -- bond angle
  !            params -- array of strings containing model parameters
  !            baseatoms -- optional atom numbers signalling dummy atom
  !-------------------------------------------------------------------
  Subroutine bbmodel_initparams(bbparams,model,angle,params,baseatoms)
    Type(BondBendingModel)                       :: bbparams
    Character(*), Intent(In)                     :: model
    Real(Kind=RDbl), Intent(In)                  :: angle
    Character(*), Dimension(:), Intent(In)       :: params
    Integer, Dimension(3), Intent(In), Optional  :: baseatoms

    Integer                   :: error

    !** Defaults
    bbparams%dummyatom = .False.
    bbparams%baseatoms = 0

    If (Present(baseatoms)) Then
      bbparams%dummyatom = .True.
      bbparams%baseatoms = baseatoms
    End If

    Select Case (toupper(Trim(model)))
    Case ('HARMONICA')
      Allocate(bbparams%hara,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

      If (Index(toupper(combine(params)),'CALC') /= 0) Then
        Call harmonica_init(bbparams%hara,params,angle)
      Else
        Call harmonica_init(bbparams%hara,params)
      End If

    Case Default
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          " : Unrecognized bb model type ",toupper(Trim(model))
      Stop
      
    End Select
    
  End Subroutine bbmodel_initparams

  !----------------------------------------------------------------------------
  ! Initialize the bond bending information
  ! Requires:  params -- bond-bending parameters
  !            spc -- species number
  !            filename -- file where bond bending information can be found
  !----------------------------------------------------------------------------
  Recursive Subroutine bbmodel_bbinit(params,spc,filename)
    Type(BendingInfo), Intent(Out)   :: params
    Integer, Intent(In)              :: spc
    Character(*), Intent(In)         :: filename

    Integer                          :: i,j,error,unitno,nchunks,n
    Integer                          :: atom1,atom2,atom3,nbends,pos,natoms
    Integer                          :: atype1,atype2,atype3
    Integer                          :: usebb,ntypes,natypes,maxv,nmatch
    Logical                          :: dummy
    Type(VecType), Dimension(3)      :: avecs
    Real(Kind=RDbl)                  :: angle
    Character(len=strLen)            :: key,string
    Character(len=255)               :: text
    Character(len=strLen), Dimension(strLen)     :: fields,chunks
    Integer, Dimension(3)                        :: baseatoms,nums
    Integer, Dimension(:), Allocatable           :: atypes
    Character(len=2), Dimension(:), Allocatable  :: elements
    Type(VecType), Dimension(:), Allocatable     :: coords
    Integer, Dimension(:,:), Allocatable         :: conlist
    Type(AtomConnections), Dimension(:), Pointer :: connections    
    Character(len=lstrLen), Dimension(:,:,:), Allocatable :: tripletparams

    !** Get number of atoms and make space for the atom types and elements
    natoms = molecules_getnatoms(spc)
    Allocate(atypes(natoms),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(coords(natoms),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(elements(natoms),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Get all the atomic coordinates
    Call molecules_getcoords(spc,coords,elements)

    !** Open the file containing the intramolecular interaction information
    unitno = file_open(filename)

    key = "Bond_Bending"

    !** Find the bond bending key word in the file
    usebb = fileSrchwoComment(unitno,key,text,.True.)

    If (usebb == 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,*) "No bond bending for this molecule- Something is wrong"
      Stop
    End If

    !** Get the number of atom types
    natypes = atom_getntypes()
    
    !** initialize the atype flag pointer
    Allocate(params%bbflag(natypes,natypes,natypes),stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    params%bbflag = .FALSE.
    
    !** Initialize the number of bends
    params%nbends = 0
    
    j = split(text,fields)
    params%model = trim(fields(3))
    !** Fast or slow?
    If (Trim(toupper(fields(4))) == "FAST") Then
      params%fast = .True.
    Else If (trim(toupper(fields(4))) == "SLOW") Then
      params%fast = .False.
    Else
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          " : Must indicate 'Fast' or 'Slow' for ", key
      Stop
    End If
    
    !** Look for 'LISTED', 'GENERATED' or another filename
    If (toupper(Trim(fields(2))) == 'LISTED') Then

      !** Bond bending parameters are listed, get number
      Read(unitno,'(a)') text
      j = split(text,fields)
      nbends = toint(fields(1))
      params%nbends = nbends
    
      !** initialize the model information pointer
      Allocate(params%bbparams(params%nbends),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bbparams')
    
      !** initialize the bblist
      Allocate(params%bblist(nbends,3),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
      Do i = 1,nbends
        dummy = .False.
        Read(unitno,'(a)') text
        j = split(text,fields)
        If (j < 5) Then
          Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
              ' A 3-body specification line must have 3 atoms and 2 parameters'
          Stop
        End If
        atom1 = toint(fields(1))
        atom2 = toint(fields(2))
        If (Index(ToUpper(fields(3)),'PERPDUMMYATOM') /= 0) Then
          nchunks = split(fields(3),chunks,'@')
          dummy = .True.
          n = str2seq(chunks(2),baseatoms)
          atom3 = -1
        Else
          atom3 = toint(fields(3))
        End If
        params%bblist(i,:) = (/atom1, atom2, atom3/)

        !** Make sure that one of the atom numbers isn't repeated
        nums = params%bblist(i,:)
        If (condenseint(nums) < 3) Then
          string = int2str(params%bblist(i,:))
          Write(0,'(2a,i5,3a)') __FILE__,":",__LINE__, &
              ' There is a repeat in the 3-atom angle set ',Trim(string)
          Stop          
        End If
    
        !** Indicate that that this atype triplet has interaction
        !** note, this is only used in a hack in bbmodel_reconnect
        atype1 = molecules_getatype(spc,atom1)
        atype2 = molecules_getatype(spc,atom2)
        If (dummy) Then
          atype3 = atype1  !** obviously HACK, but not very important
        Else
          atype3 = molecules_getatype(spc,atom3)
        End If
        params%bbflag(atype1,atype2,atype3) = .TRUE.
        params%bbflag(atype3,atype2,atype1) = .TRUE.
    
        !** Get the angle in case we need it
        avecs(1) = coords(atom1)
        avecs(2) = coords(atom2)
        If (dummy) Then
          avecs(3) = bbmodel_dummy(baseatoms,coords)
        Else
          avecs(3) = coords(atom3)
        End If
        angle = vector_angle(avecs(1),avecs(2),avecs(3))*degTorad
        If (dummy) Then
          Call bbmodel_initparams(params%bbparams(i),params%model, &
              angle,fields(4:j),baseatoms)
        Else
          Call bbmodel_initparams(params%bbparams(i),params%model, &
              angle,fields(4:j))
        End If
      End Do

    Else If (toupper(Trim(fields(2))) == 'GENERATE') Then

      !** We need to generate the bond bending list
      !** A list of atom triplets is needed for the list, just
      !** like we have a list of atom pairs for bond stretching
      Call molecules_getconnects(spc,connections)
    
      !** Read in the number of types list
      Read(unitno,'(a)') text
      j = split(text,fields)
      ntypes = toint(fields(1))

      !** Allocate temporary storage for atomtype-based parameters
      Allocate(tripletparams(natypes,natypes,natypes),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
      !** Get the number of atom types
      natypes = molecules_getnatomtypes(spc,atypes)
    
      !** Get the maximum valency of the atom types
      maxv = 0
      Do i = 1, natypes
        If (atom_nbonds(atypes(i)) > maxv) maxv = atom_nbonds(atypes(i))
      End Do
    
      !** Allocate temporary storage for the connectivity list.
      !** The total choices for the first atom are the number of atoms 
      !** in the molecule. The second is the maximum valency. The third, 
      !** then, has maxv*maxv atom possibilities
      Allocate(conlist(natoms*maxv**3,3),stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i5,a,i12)') __FILE__,":",__LINE__, &
            " Could not allocate conlist. Required size : ", &
            natoms*maxv**3
        Stop
      End If
    
      !** Get all the individual atom types
      natypes = molecules_getnatomtypes(spc,atypes,.True.)
    
      !** We can loop through the types of listed connections.
      Do i = 1, ntypes
        
        !** Read in the type and store it
        !** The type should consist of 3 atom type names and the
        !** parameters corresponding to the model that was specified
        Read(unitno,'(a)') text
        text = stripcmnt(text)
        j = split(text,fields)
        atype1 = atom_gettypename(fields(1))
        atype2 = atom_gettypename(fields(2))
        atype3 = atom_gettypename(fields(3))
        params%bbflag(atype1,atype2,atype3) = .TRUE.
        params%bbflag(atype3,atype2,atype1) = .TRUE.
        pos = findint((/atype1,atype2,atype3/),0)
        If (pos /= 0) Then
          Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
              " Bending could not recognized atom type ",fields(pos)
          Stop
        End If
    
        !** save the parameters
        tripletparams(atype1,atype2,atype3) = combine(fields(4:j))
        tripletparams(atype3,atype2,atype1) = combine(fields(4:j))
    
        !** With the specific type initialized, we need to build the 
        !** list of atom triplets that correspond to that type
        nmatch = connects_makeTypeList(connections, atypes, &
            (/atype1,atype2,atype3/),conlist(params%nbends+1:Size(conlist,1),:))
    
        !** Record the number of bond bending triplets
        params%nbends = nmatch + params%nbends
      End Do
    
      !** initialize the model information pointer
      Allocate(params%bbparams(params%nbends),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bbparams')
    
      !** Allocate the bblist. There are nbends, each with 3 atoms
      Allocate(params%bblist(params%nbends,3),stat=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,'bblist')
    
      !** Initialize the bond bending parameters
      Do i = 1,params%nbends



        atom1 = conlist(i,1)
        atom2 = conlist(i,2)
        atom3 = conlist(i,3)
        atype1=molecules_getatype(spc, atom1)
        atype2=molecules_getatype(spc, atom2)
        atype3=molecules_getatype(spc, atom3)
        params%bblist(i,:) = (/atom1, atom2, atom3/)

        avecs(1) = coords(atom1)
        avecs(2) = coords(atom2)
        avecs(3) = coords(atom3)
        angle = vector_angle(avecs(1),avecs(2),avecs(3))*degTorad
        j = split(tripletparams(atype1,atype2,atype3),fields)
        Call bbmodel_initparams(params%bbparams(i),params%model, &
            angle,fields(1:j))
      End Do
  
  
      !** Clean up temporary storage
      Deallocate(conlist,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'conlist')
      Deallocate(tripletparams,stat=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__,'tripletparams')
    
    Else
      !** assume it's a file we are reading in, recall this routine
      Close(unitno)
      Call bbmodel_bbinit(params,spc,Trim(fields(2)))
    End If

    Deallocate(atypes,stat=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    Deallocate(coords,stat=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    Deallocate(elements,stat=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)

  End Subroutine bbmodel_bbinit

  !----------------------------------------------------------------------------
  ! Calculating bond bending energy and (optionally) forces.  Either 3 atomic 
  ! coordinate vectors or a separation vector can be passed in.  
  ! Requires: bspntr -- bond stretching pointer set
  !           coords -- atomic coordinates (3), or separation vector
  !           sepvec -- flag indicating that separation vectors are passed
  !           u -- potential 
  !           vf -- optional forces on the two atoms
  !----------------------------------------------------------------------------
  Subroutine bbmodel_getinteraction(bbparams,coords,sepvec,u,vf)
    Type(BondBendingModel)                             :: bbparams
    Type(VecType), Dimension(:), Intent(In)            :: coords
    Logical, Intent(In)                                :: sepvec
    Real(kind=RDbl), Intent(Out)                       :: u       
    Type(VecType), Dimension(:), Intent(Out), Optional :: vf

    If (Present(vf)) Then
      If (Associated(bbparams%hara)) Then
        Call harmonica_getinteraction(bbparams%hara,coords,sepvec,u,vf)
      End If
    Else
      If (Associated(bbparams%hara)) Then
        Call harmonica_getinteraction(bbparams%hara,coords,sepvec,u)
      End If
    End If
    
  End Subroutine bbmodel_getinteraction

  !----------------------------------------------------------------------------
  ! Pick or calculate the three atomic coordinate vectors to be used in the 
  ! bond angle calculations.  
  ! Requires:  params -- bond bending forcefield parameters
  !            list -- list of atoms numbers
  !            coords -- atomic coordinates for a molecule
  !            avecs -- 3 output atomic vectors
  !----------------------------------------------------------------------------
  Subroutine bbmodel_atomvectors(bbparams,list,coords,avecs)
    Type(BondBendingModel), Intent(In)         :: bbparams
    Integer, Dimension(:), Intent(In)          :: list
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(VecType), Dimension(3), Intent(Out)   :: avecs

    !** Get the first two the easy way
    avecs(1) = coords(list(1))
    avecs(2) = coords(list(2))

    !** Get the third the easy way unless it's a dummy atom
    If (bbparams%dummyatom) Then
      avecs(3) = bbmodel_dummy(bbparams%baseatoms,coords)
    Else
      avecs(3) = coords(list(3))
    End If

  End Subroutine bbmodel_atomvectors

  !---------------------------------------------------------------------------
  ! Create the coordinates for a dummy atom that is perpendicular to two
  ! given atoms at another atom.  The vector between the dummy atom and
  ! the first atom will be perpendicular to the plane formed by the last 
  ! two atoms and 5.0 Angstroms away. 
  ! Requires:  atoms -- list of atom numbers
  !            acoords -- list of all atomic position vectors in molecule
  !---------------------------------------------------------------------------
  Type(VecType) Function bbmodel_dummy(atoms,acoords)
    Integer, Dimension(3), Intent(In)         :: atoms
    Type(VecType), Dimension(:), Intent(In)   :: acoords

    Type(VecType)       :: bvec1,bvec2,cprod

    bvec1 = acoords(atoms(2)) - acoords(atoms(1))
    bvec2 = acoords(atoms(3)) - acoords(atoms(1))
    cprod = vector_crossprod(bvec1,bvec2)
    bbmodel_dummy = 5.0_RDbl*vector_getunitvec(cprod)

  End Function bbmodel_dummy

  !----------------------------------------------------------------------------
  ! Calculate the bond bending interaction using actual coordinates of atoms.
  ! Periodic boundary conditions will be applied to the atomic coordinates if
  ! the simcell is passed.
  ! Requires: params -- bond bending forcefield parameters
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Subroutine bbmodel_bbint(params,ffout,coords,simcell)
    Type(BendingInfo), Intent(In)              :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell
    
    Integer                           :: i,nderivs,level
    Integer                           :: atom1, atom2, atom3
    Real(kind=RDbl)                   :: u,upart
    Type(VecType), Dimension(2)       :: sepvec
    Type(VecType), Dimension(3)       :: force,avecs

    !** Get the type of the passed storage
    level = store_idbranch(ffout)

    !** loop through the bending triplets and get interactions from bbmodel
    nderivs = store_nderivs(ffout)

    Do i = 1,params%nbends
      atom1 = params%bblist(i,1)
      atom2 = params%bblist(i,2)
      atom3 = params%bblist(i,3)
      Call bbmodel_atomvectors(params%bbparams(i),params%bblist(i,:), &
          coords,avecs)

      !** Get the separation vectors
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,avecs(1),avecs(2),sepvec(1))
        Call simcell_minimage(simcell,avecs(3),avecs(2),sepvec(2))
      Else
        sepvec(1) = avecs(1) - avecs(2)
        sepvec(2) = avecs(3) - avecs(2)
      End If

      !** Do the evaluation and store the results
      If (nderivs == 0) Then
        Call bbmodel_getinteraction(params%bbparams(i),sepvec,.True.,u)

      Else If (nderivs == 1) Then
        Call bbmodel_getinteraction(params%bbparams(i),sepvec,.True.,u,force)

        !** Store the results according to the storage form
        Select Case(level)
        Case(0)
          Call storebase_inc(ffout%binfo(atom1),force(1))
          Call storebase_inc(ffout%binfo(atom2),force(2))
          If (.Not. params%bbparams(i)%dummyatom) &
              Call storebase_inc(ffout%binfo(atom3),force(3))
        Case(10)
          Call storebase_inc(ffout%mi(atom1)%total,force(1))
          Call storebase_inc(ffout%mi(atom2)%total,force(2))
          If (.Not. params%bbparams(i)%dummyatom) &
              Call storebase_inc(ffout%mi(atom3)%total,force(3))
        Case Default
          Write(0,'(2a,i5,a,i3)') __FILE__,":",__LINE__, &
              ' Unable to handle storage type ',level
          Stop        
        End Select

      Else
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            ' Second order terms not available'
        Stop        
      End If

      !** Store the potential energy results according to the storage form
      Select Case(level)
      Case(-1)
        Call storebase_inc(ffout%total,u,BENDING_INDEX)
      Case(0)
        upart = u/3.0_RDbl
        Call storebase_inc(ffout%binfo(atom1),upart,BENDING_INDEX)
        Call storebase_inc(ffout%binfo(atom2),upart,BENDING_INDEX)
        If (.Not. params%bbparams(i)%dummyatom) &
            Call storebase_inc(ffout%binfo(atom3),upart,BENDING_INDEX)
      Case(10)
        upart = u/3.0_RDbl
        Call storebase_inc(ffout%mi(atom1)%total,upart,BENDING_INDEX)
        Call storebase_inc(ffout%mi(atom2)%total,upart,BENDING_INDEX)
        If (.Not. params%bbparams(i)%dummyatom) &
            Call storebase_inc(ffout%mi(atom3)%total,upart,BENDING_INDEX)
      Case Default
        Write(0,'(2a,i5,a,i3)') __FILE__,":",__LINE__, &
            ' Unable to handle storage type ',level
        Stop        
      End Select

    End Do

  End Subroutine bbmodel_bbint

#ifdef PREVIOUS
  !----------------------------------------------------------------------------
  ! Calculate the bond bending interaction using actual coordinates of atoms.
  ! Periodic boundary conditions will be applied to the atomic coordinates if
  ! the simcell is passed.
  ! Requires: params -- bond bending forcefield parameters
  !           ffout -- forcefield output for a single molecule
  !           coords -- atomic coordinates for a molecule
  !           simcell -- simulation cell information
  !----------------------------------------------------------------------------
  Subroutine bbmodel_bbint(params,ffout,coords,simcell)
    Type(BendingInfo), Intent(In)              :: params
    Type(Store_Level_Pair), Intent(InOut)      :: ffout
    Type(VecType), Dimension(:), Intent(In)    :: coords
    Type(SimCell_Params), Intent(In), Optional :: simcell
    
    Integer                           :: i,nderivs
    Integer                           :: atom1, atom2, atom3
    Real(kind=RDbl)                   :: u,upart
    Type(VecType), Dimension(2)       :: sepvec
    Type(VecType), Dimension(3)       :: grad

    !** loop through the bending triplets and get interactions from bbmodel
    nderivs = store_nderivs(ffout)
    Do i = 1,params%nbends
      atom1 = params%bblist(i,1)
      atom2 = params%bblist(i,2)
      atom3 = params%bblist(i,3)
      Call bbmodel_atomvectors(params%bbparams(i),params%bblist(i,:), &
          coords,avecs)

      !** get the separation vectors
      If (Present(simcell)) Then
        Call simcell_minimage(simcell,avecs(1),avecs(2),sepvec(1))
        Call simcell_minimage(simcell,avecs(3),avecs(2),sepvec(2))
      Else
        sepvec(1) = avecs(1) - avecs(2)
        sepvec(2) = avecs(3) - avecs(2)
      End If

      !** Do the evaluation and store the results
      If (nderivs == 0) Then
        Call bbmodel_getinteraction(params%bbparams(i),sepvec,.True.,u)
      Else If (nderivs == 1) Then
        Call bbmodel_getinteraction(params%bbparams(i),sepvec,.True.,u,grad)
        Call storebase_inc(ffout%mi(atom1)%total,grad(1))
        Call storebase_inc(ffout%mi(atom2)%total,grad(2))
        If (.Not. params%bbparams(i)%dummyatom) &
            Call storebase_inc(ffout%mi(atom3)%total,grad(3))
      Else
        Write(0,'(2a,i5,a)') __FILE__,":",__LINE__, &
            ' Second order terms not available'
        Stop        
      End If

      upart = u/3.0_RDbl
      Call storebase_inc(ffout%mi(atom1)%total,upart,BENDING_INDEX)
      Call storebase_inc(ffout%mi(atom2)%total,upart,BENDING_INDEX)
      If (.Not. params%bbparams(i)%dummyatom) &
          Call storebase_inc(ffout%mi(atom3)%total,upart,BENDING_INDEX)
    End Do

  End Subroutine bbmodel_bbint
#endif

  !--------------------------------------------------------------------------
  ! Returns the params used here
  ! As of now we have only one type of bending params with 2 parameters
  ! ktheta and theta_eq, They are respectively params(1) and params(2)
  !--------------------------------------------------------------------------
  Subroutine bbmodel_getparams(bbparams,model,params)
    Type(BondBendingModel),Intent(in)    :: bbparams
    Character(len=strLen),Intent(out)    :: model
    Real(kind=RDbl), Dimension(:), Intent(out)        :: params

    If (Associated(bbparams%hara)) Then
      model = HarmonicaString
      Call harmonica_getparams(bbparams%hara,params)
    Else
      Write(*,*) "Angle bending not initialized well ??"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    
  End Subroutine bbmodel_getparams

  !----------------------------------------------
  ! Checks if the BB interaction is fast or slow
  !----------------------------------------------
  Logical Function bbmodel_isfastbb(bending)
    Type(BendingInfo) :: bending

    bbmodel_isfastbb = bending%fast
  End Function bbmodel_isfastbb

  !--------------------------------------------------------------------
  ! Determines if the interaction should be evaluated
  ! Requires:  ptr -- pointer to interaction parameters structure
  !            fast -- proposed evaluation speed
  !--------------------------------------------------------------------
  Logical Function bbmodel_eval(ptr,fast)
    Type(BendingInfo), Pointer    :: ptr
    Logical, Intent(In)           :: fast

    bbmodel_eval = .False.
    If (.Not. Associated(ptr)) Return
    bbmodel_eval = (ptr%fast .And. fast)

  End Function bbmodel_eval

  !----------------------------------------------------------------------------
  ! Returns all tabulated bond-bending information
  ! Requires:  info -- bending data structure for a single species
  !            nsets -- number of potential sets returned
  !            list -- list of atom number pairs in molecule  (pair_index,1:2)
  !            params -- array of strings containing type and parameters
  !----------------------------------------------------------------------------
  Subroutine bbmodel_bbinfo(info,nsets,list,params)
    Type(BendingInfo), Intent(In)                       :: info
    Integer, Intent(Out)                                :: nsets
    Integer, Dimension(:,:), Intent(Out)                :: list
    Character(len=lstrLen), Dimension(:), Intent(Out)   :: params

    Integer          :: i

    !** copy the list
    nsets = info%nbends
!LC    Write(*,*) nsets, Size(list,1), Size(list,2)
    list(1:nsets,1:3) = info%bblist(1:nsets,1:3)

    !** get the parameters
    Do i = 1,nsets
      params(i) = bbmodel_display(info%bbparams(i))
    End Do

  End Subroutine bbmodel_bbinfo

  !-------------------------------------------------------
  ! Display the bb parameters
  !-------------------------------------------------------
!!$  Character(len=strLen) Function bbmodel_display(bbparams)
  Function bbmodel_display(bbparams)
    Character(len=strLen) :: bbmodel_display
    Character(strLen) :: mname
    Type(BondBendingModel) :: bbparams
    Character(strLen) :: vals

    vals = ""

    If (Associated(bbparams%hara)) Then
      mname="Harmonica"
      vals = Trim(harmonica_display(bbparams%hara))
    Else If (.Not.Associated(bbparams%hara)) Then
      mname = "unassociated"
    Else
      mname = "undetermined"
    End If

    Write(bbmodel_display,'(a,2x,a)') trim(mname),trim(vals)

  End Function bbmodel_display

  !-----------------------------------------------
  ! Display the information about bond bending
  !-----------------------------------------------
  Subroutine bbmodel_bbdisplay(ptr,spc,unit,indent)
    Type(BendingInfo), Pointer :: ptr
    Integer, Intent(In)        :: spc,unit,indent

    Integer                  :: i, atom1, atom2, atype1,atype2,atom3, atype3
    Character(len=strLen)    :: name
    Character(len=xlstrLen)  :: string1,string2
    Character(len=indent)    :: blank

    name = molecules_name(spc)

    blank = Repeat(' ',indent)

    If (Associated(ptr)) Then
      Write(unit,'(3a)') blank,'Bond Bending Information for ', &
          Trim(name)
      If (ptr%fast) Then
        Write(unit,'(2x,2a,l1)') blank,'Fast Interaction'
      Else
        Write(unit,'(2x,2a,l1)') blank,'Slow Interaction'
      End If
      Do i =1,ptr%nbends
        atom1 = ptr%bblist(i,1)
        atom2 = ptr%bblist(i,2)
        If (ptr%bbparams(i)%dummyatom) Then
          string2 = int2str(ptr%bbparams(i)%baseatoms)
          Write(string1,'(2a)') 'PERPDUMMYATOM@',Trim(string2)
        Else
          string1 = int2str(ptr%bblist(i,3))
        End If
        Write(unit,'(2x,2a,2i3,1x,a)') blank,'Bond bending groups :',atom1, &
            atom2,Trim(string1)
        string1 = bbmodel_display(ptr%bbparams(i))
        Write(unit,'(4x,3a)') blank,'Model and Params: ',Trim(string1)
      End Do
    Else
      Write(unit,'(3a)') blank,'No bond bending information for ',Trim(name)
    End If

  End Subroutine Bbmodel_bbdisplay

  !-------------------------------------------------------
  ! Cleanup 
  !-------------------------------------------------------
  Subroutine bbmodel_cleanup(bbparams)
    Type(BondBendingModel)    :: bbparams

    Integer :: error

    If (Associated(bbparams%hara)) Then
      Deallocate(bbparams%hara,STAT=error)
      Call deallocErrDisplay(__FILE__,__LINE__,'hara')
    End If

  End Subroutine bbmodel_cleanup

  !-------------------------------------------
  ! Cleanup bending info
  !-------------------------------------------
  Subroutine bbmodel_bbcleanup(bending)

    Type(BendingInfo) :: bending
    Integer :: i,error

    If (Associated(bending%bbparams)) Then
      Do i = 1,Size(bending%bbparams,1)
        Call bbmodel_cleanup(bending%bbparams(i))
      End Do
    End If
    Deallocate(bending%bblist,STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
          " : Error deallocating bending%blist"
      Stop
    End If
    Deallocate(bending%bbparams,STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
          " : Error deallocating bending"
      Stop
    End If
    
  End Subroutine bbmodel_bbcleanup

  !----------------------------------------------------------------------------
  ! set harmonica params to a set valu
  !----------------------------------------------------------------------------
  Subroutine bbmodel_setHarK(bbpntr, alist,harK)
    Type(BendingInfo), Pointer                 :: bbpntr
    Integer, Dimension(3), Intent(In)            :: alist
    Real(kind=RDbl), Intent(in)                       :: harK
    Integer :: index,atom1, atom2 ,atom3, i
    index=0
    Do i = 1,bbpntr%nbends
      index=i
      atom1 = bbpntr%bblist(i,1)
      atom2 = bbpntr%bblist(i,2)
      atom3 = bbpntr%bblist(i,3)
      If ( (atom1==alist(1)).And.(atom2==alist(2)).And.(atom3==alist(3)) ) Exit
    End Do
    If (index==0) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    If (Associated(bbpntr%bbparams(index)%hara)) Then
      bbpntr%bbparams(index)%hara%ktheta=harK
    Else 
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

  End Subroutine bbmodel_setHarK

End Module bbmodel
