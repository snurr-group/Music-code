!------------------------------------------------------------------------------
! Nudged Elastic Band (NEB) module.  This module handles the NEB method
! See these references:
!   - H. Jonsson and G. Mills and K. W. Jacobsen, "Nudged Elastic Band 
!     Method for Finding Minimum Energy Paths of Transitions", in 
!     book: Classical and Quantum Dynamics in Condensed Phase Simulations,
!     World Scientifc, pages 385-404.
!   - G. Henkelman and H. Jonsson, JCP 113 (2000) 9978-9985
!
! The strategy is to handle the different images each with their own 
! configuration and interaction storage structures.  Spring interactions
! between selected atoms in the images are dealt with manually using a 
! direct call to pairmodel.
!
! I find that it is best to run 100 cycles or so of the steepest descent
! type of step with 20 cycles of modified steepest descent "pre-optimization".
! After this, the modified velocity verlet method that Jonsson's group 
! recommends can be used.  The integration type steps are less stable, but
! probably more efficient since they use past information in terms of 
! accelerations.  They will certainly blow-up when used on super-high energy
! configuration produced by interpolation.
!
! Needed Improvements:
! 1) Spring interaction evaluations between images are done directly in
!    this module, this may well run counter to the philosophy in interact
! 2) The step size adjustments are done in this routine, this may take
!    too much control from optstep, especially since the step types are
!    so different.
!------------------------------------------------------------------------------

Module neb

  Use defaults, Only: RDbl, strLen, lstrLen, kjmole_kb, calToj, &
      zero, dashedline, shortdashedline, TUnit, velUnit, nrgUnit
  Use file, Only: file_getunit, file_open, file_close
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, &
      tolower, toupper, allocerrdisplay, int2str, real2str, str2seq, &
      checkandstop, deallocerrdisplay, findint
  Use general, Only: genparams
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag, vector_getunitvec, vector_display, &
      vector_getnorm
  Use atom, Only: atom_getsymbol, atom_getntypes
  Use molecules, Only: molecules_getnsorbs, molecules_name
  Use random, Only: random_gaussian
  Use config, Only: AtMolCoords, config_getnmoles, config_getnatoms, &
      config_isfixed, config_makeimages, config_config2xyz, &
      config_changefixed, config_getrp, config_getr, config_getatype, &
      config_checkinit, config_initcopy, config_copysubset, config_setxyz, &
      config_simpleinit
  Use simcell, Only: SimCell_Params, simcell_pbc, simcell_getfillsorb
  Use storebase, Only: EnergyPlus, storebase_init, storebase_project, &
      storebase_clean, storebase_inc, storebase_display
  Use storetop, Only: Forcefield_Results, storetop_totnrg, storetop_fastsum, &
      storetop_projectout, storetop_basicptrs, storetop_display
  Use interact, Only: Interaction_Model, interact_int
  Use forcefield, Only: forcefield_storelevel
  Use pairmodel, Only: Pairmodel_Params, pairmodel_simpleinit, &
      pairmodel_clean, pairmodel_mmint4, pairmodel_nullset
  Use optstep, Only: Optimization_Step, optstep_simpleinit, optstep_move, &
      optstep_setstep, optstep_stepsize, optstep_display, optstep_clean, &
      stpdesc_idstring, modstpdesc_idstring, modvv_idstring
  Use visxyz, Only: XYZ_ENTRY, visxyz_make, visxyz_dump
      
  Implicit None
  Save

  Private
  Public :: NEB_Info, neb_tag, neb_initbasic, neb_initrest, neb_movechain, &
      neb_configtag, neb_nimages, neb_images2xyz, neb_modimages2xyz, &
      neb_modimage2xyz, neb_sampleCF, neb_displaystats, neb_display, neb_clean

  Type NEB_Info
    Integer                     :: nimages,preopt,nspc
    Logical                     :: inter_on,project_intra
    Real(kind=RDbl)             :: spring_const, minstep, max_upnrg
    Character(len=lstrLen)      :: start_tag,end_tag,config_origin
    Character(len=lstrLen)      :: fixedspcs,steptype
    Type(Optimization_Step)     :: stepinfo
    Type(Pairmodel_Params), Dimension(:,:), Pointer :: sprconst

    !** Statistics for the images
    Integer, Dimension(:), Pointer           :: nresets,nmoves
    Real(kind=RDbl), Dimension(:), Pointer   :: totnrg,lastnrg  !* kJ/mol
    Real(kind=RDbl), Dimension(:), Pointer   :: stepsize  !* Angstroms or ps
    Real(kind=RDbl), Dimension(:,:), Pointer :: avgdist   !* image,spc (Angs)
  End Type NEB_Info

  Character(len=lstrLen), Parameter :: neb_tag = 'NEB Information'

Contains

  !----------------------------------------------------------------------------
  ! Initializes the NEB parameters from the control file. 
  ! Requires:  params -- NEB parameter data structure
  !            ctrl_filename -- control filename to find init info
  !----------------------------------------------------------------------------
  Subroutine neb_initbasic(params,ctrl_filename)
    Type(NEB_Info), Intent(InOut)    :: params
    Character(*), Intent(In)         :: ctrl_filename

    Integer                               :: i,j,unit,error,nfields
    Character(len=lstrLen)                :: line
    Character(len=255)                    :: text
    Character(len=strLen), Dimension(20)  :: fields

    !** Open the control file if it is not already open
    unit = file_open(ctrl_filename)

    !** Find the MD section
    If (filesrchstr(unit,neb_tag,text,.True.) == 0) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          " Could not find NEB parameters in control file"
      Write(0,'(2a)') 'Looking for Tag: ',Trim(neb_tag)
      Stop
    End If

    !** Read the basics
    Read(unit,'(i5)') params%nimages
    params%nimages = params%nimages + 2
    Read(unit,'(a)') line
    params%fixedspcs = stripcmnt(line)
    Read(unit,'(i5)') params%preopt
    Read(unit,'(f12.6)') params%spring_const
    Read(unit,'(f12.6)') params%max_upnrg
    If (params%max_upnrg < 0.0_RDbl) params%max_upnrg = 0.0_RDbl
    Read(unit,'(a)') line
    line = stripcmnt(line)
    Select Case (ToUpper(line))
    Case('NO')
      params%project_intra = .False.
    Case('YES')
      params%project_intra = .True.
    Case Default
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Must use YES or NO to specify intra-image projections ',Trim(line)
      Stop
    End Select
    Read(unit,'(a)') line
    params%steptype = stripcmnt(line)

    !** Set the inter-image interactions on or off
    params%inter_on = .True.
    If (params%spring_const < 1.0e-10) params%inter_on = .False.
      
    !** Read in the control file tag for the starting configuration
    Read(unit,'(a)') text
    text = stripcmnt(text)
    nfields = split(text,fields)
    params%start_tag = fields(1)

    !** Read in the command for generating in-between configurations
    Read(unit,'(a)') text
    text = stripcmnt(text)
    nfields = split(text,fields)
    params%config_origin = fields(1)

    !** Read in the control file tag for the ending configuration
    Read(unit,'(a)') text
    text = stripcmnt(text)
    nfields = split(text,fields)
    params%end_tag = fields(1)

  End Subroutine neb_initbasic

  !----------------------------------------------------------------------------
  ! Initializes the remainder of the NEB parameters once at least the end 
  ! points are established.  
  ! Requires:  params -- NEB parameter data structure
  !            image -- the set of configuration images
  !            imodel -- interaction model
  !            simcell -- simulation cell structure
  !----------------------------------------------------------------------------
  Subroutine neb_initrest(params,image,imodel,simcell)
    Type(NEB_Info), Intent(InOut)                     :: params
    Type(AtMolCoords), Dimension(:,:), Intent(InOut)  :: image
    Type(Interaction_Model), Intent(InOut)            :: imodel
    Type(SimCell_Params), Intent(In)                  :: simcell

    Integer           :: i,j,error,atype
    Integer           :: nimages,nfixed,natypes
    Logical           :: fast,success
    Real(kind=RDbl)   :: junk
    Character(len=3)  :: storelevel
    Character(len=strLen)    :: string
    Character(len=lstrLen)   :: filename
    Integer, Dimension(10)   :: nums

    nimages = params%nimages
    params%nspc = Size(image,2)

    !** Make sure the forcefield storage wasn't initialized at ATM-ATM level
    storelevel = forcefield_storelevel(imodel%ff)
    
    If (storelevel == 'ATM') Then
      Write(0,'(1x,2a,i4,3a)') __FILE__,' : ',__LINE__, &
          ' Cannot handle ATM-depth forcefield storage, (why?)'
      Stop      
    End If

    !** Fix the specified endpoint species if desired
    nfixed = str2seq(params%fixedspcs,nums)
    Do i = 1,nfixed
      If (nums(i) <= 0) Cycle
      Call config_changefixed(image(1,:),nums(i),.True.)
      Call config_changefixed(image(nimages,:),nums(i),.True.)
    End Do

    !** Initialize total energy storage for each image
    Allocate(params%totnrg(nimages), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(params%lastnrg(nimages), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    params%totnrg = 0.0_RDbl
    params%lastnrg = 0.0_RDbl

    !** Initialize the move statistics for each image
    Allocate(params%nmoves(nimages), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Allocate(params%nresets(nimages), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    params%nmoves = 0
    params%nresets = 0

    !** Fill in the images between the start and end points if required
    Select Case(ToUpper(params%config_origin))
    Case ('INTERPOLATE')
      Call config_makeimages(image(1,:),image(nimages,:), &
          image(2:nimages-1,:),simcell)

    Case ('READ')
      nfixed = str2seq(params%fixedspcs,nums)
      Do i = 2,params%nimages-1
        string = int2str(i)
        Write(filename,'(3a)') Trim(genparams%restartfile),Trim(string),'.res'
        Call config_simpleinit(image(i,:), simcell, 'RESTARTFILE', filename)

        !** Fix the species if necessary
        Do j = 1,nfixed
          If (nums(j) > 0) Call config_changefixed(image(i,:),nums(j),.True.)
        End Do
      End Do

    Case Default
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Could not interpret image initialization method: ', &
          Trim(params%config_origin)
      Stop
    End Select

    !** Get energies for at least the two end-point configurations
    Do i = 1,params%nimages
      If (.Not. config_checkinit(image(1,:))) Cycle
      
      !** Evaluate intra-image forces for whole system
      fast = .True.
      success = interact_int(imodel,image(i,:),junk,300.0_RDbl,simcell, &
          fast,.True.,.False.)
      Call checkandstop(success,__FILE__,__LINE__, &
          ' Forcefield evaluation problem')

      !** Sum the forcefield storage structure and extract total energy
      Call storetop_fastsum(imodel%ff%results)
      params%totnrg(i) = storetop_totnrg(imodel%ff%results)*calToj
    End Do

    !** Initialize spring constants between like-atype atoms
    If (params%inter_on) Then
      natypes = atom_getntypes()
      Allocate(params%sprconst(natypes,natypes), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Do atype = 1,natypes
        Do i = 1,natypes
          Call pairmodel_nullset(params%sprconst(atype,i))
        End Do
        Call pairmodel_simpleinit(params%sprconst(atype,atype),atype,atype, &
            'HARMONIC',(/params%spring_const,0.0_RDbl/))
      End Do
    Else
      Nullify(params%sprconst)
    End If

    params%minstep = 1.0e-4_RDbl
    !** Initialize stepsize storage for each image
    Allocate(params%stepsize(nimages), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Initialize modified steepest-descent method for first few steps
    If (params%preopt > 0) Then
      Call neb_initstep(params,modstpdesc_idstring)
    Else 
      Call neb_initstep(params,params%steptype)
    End If

    !** Initialize average distance storage for each image,species
    Allocate(params%avgdist(nimages,params%nspc), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

  End Subroutine neb_initrest

  !----------------------------------------------------------------------------
  ! Initializes the optimization step type
  ! Requires:  params -- NEB parameter data structure
  !            identifier -- a string identifying
  !            stepsize -- optional step size, will pick default otherwise
  ! NOTE: for now, always uses default with modified velocity verlet
  !----------------------------------------------------------------------------
  Subroutine neb_initstep(params,identifier,stepsize)
    Type(NEB_Info), Intent(InOut)         :: params
    Character(*), Intent(In)              :: identifier
    Real(Kind=RDbl), Intent(In), Optional :: stepsize

    Real(Kind=RDbl)       :: initial_step

    !** Identify method and pick parameters accordingly
    Select Case(ToUpper(identifier))
    Case(stpdesc_idstring)
      !** variable is the step size in Angstroms
      initial_step = 0.1_RDbl
      If (Present(stepsize)) initial_step = stepsize
      Call optstep_simpleinit(params%stepinfo,identifier,(/initial_step/))

    Case(modstpdesc_idstring)
      !** variables are the step size in Angstroms and smear factor
      initial_step = 0.1_RDbl
      If (Present(stepsize)) initial_step = stepsize
      Call optstep_simpleinit(params%stepinfo,identifier, &
          (/initial_step,0.90_RDbl/))

    Case(modvv_idstring)
      !** variables are the timestep in picoseconds and maximum step
      initial_step = 0.001_RDbl
      If (Present(stepsize)) initial_step = stepsize
!      Call optstep_simpleinit(params%stepinfo,identifier, &
!          (/initial_step,0.01_RDbl/))
      Call optstep_simpleinit(params%stepinfo,identifier, &
          (/initial_step/))

    Case Default
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          '  Could not identify optimization step method: ',Trim(identifier)
      Stop
    End Select

    !** Set the stepsize variables
    params%stepsize = initial_step 
    params%stepsize(1) = 0.0_RDbl
    params%stepsize(params%nimages) = 0.0_RDbl

  End Subroutine neb_initstep

  !----------------------------------------------------------------------------
  ! Moves the image chain through one iteration.
  ! Requires:  params -- NEB parameter data structure
  !            image -- the set of configuration images
  !            imodel -- interaction model
  !            simcell -- simulation cell structure
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------------
  Subroutine neb_movechain(params,image,imodel,simcell,indent,unit)
    Type(NEB_Info), Intent(InOut)                     :: params
    Type(AtMolCoords), Dimension(:,:), Intent(InOut)  :: image
    Type(Interaction_Model), Intent(InOut)            :: imodel
    Type(SimCell_Params), Intent(In)                  :: simcell
    Integer, Intent(In)                               :: indent,unit

    Integer                    :: i,j,error
    Integer                    :: nspc,nmols,natms,spc,mol,atm
    Logical                    :: success,fast,careful
    Logical, Save              :: firsttime = .True.
    Real(kind=RDbl)            :: junk,delnrg
    Character(len=indent)      :: blank
    Type(VecType), Dimension(:,:), Allocatable           :: tangent
    Type(AtMolCoords), Dimension(:,:), Allocatable, Save :: saveimage

    nspc = Size(image,2)

    !** Create storage for the saved configurations
    If (firsttime) Then
      firsttime = .False.
      Allocate(saveimage(params%nimages,nspc), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"saveimage")      
      Do i = 2,params%nimages-1
        Call config_initcopy(saveimage(i,:),image(i,:),.False.)
      End Do
    End If

    !** Switch to normal steepest descent after a few steps
    If ((params%nmoves(2) == params%preopt).And.(params%preopt /= 0)) Then
      blank = Repeat(' ',indent)
      Write(unit,'(2a)') blank,'Switching to normal optimization.  This means:'
      Write(unit,'(2a)') blank,'  1) intra-images force removed along tangent'
      Write(unit,'(2a)') blank,'  2) primary optimization method used'
      Call optstep_clean(params%stepinfo)
      Call neb_initstep(params,params%steptype)
    End If

    !** Save the previous energies
    params%lastnrg = params%totnrg

    !** Loop through the non-endpoint images, evaluate forces and move
    Do i = 2,params%nimages-1

      !** Evaluate intra-image forces for whole configuration
      fast = .True.
      success = interact_int(imodel,image(i,:),junk,300.0_RDbl,simcell, &
          fast,.True.,.False.)
      Call checkandstop(success,__FILE__,__LINE__, &
          ' Forcefield evaluation problem')

      !** Sum the forcefield storage structure and extract total energy
      Call storetop_fastsum(imodel%ff%results)
      params%totnrg(i) = storetop_totnrg(imodel%ff%results)*calToj

      !** Check the energy difference, reset to last config if necessary
      careful = .False.
      delnrg = params%totnrg(i) - params%lastnrg(i)
      If ((delnrg >= 0.0_RDbl).And.(params%nmoves(2) > 1)) Then

        !** Reset the previous configuration and energy if desired
        If (delnrg >= params%max_upnrg) Then
          Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
          Write(*,*) 'Resetting image ',i,delnrg
          Call config_copysubset(saveimage(i,:),image(i,:),(/0,0,0/),(/0,0,0/))
          params%totnrg(i) = params%lastnrg(i)
          params%nresets(i) = params%nresets(i) + 1
          careful = .True.
        End If

        !** Decrease the stepsize
        params%stepsize(i) = 0.75_RDbl*params%stepsize(i)
        If (params%stepsize(i) < params%minstep) Then
          params%stepsize(i) = params%minstep
        End If

      Else If (delnrg < 0.0_RDbl) Then
        params%stepsize(i) = params%stepsize(i)*1.10_RDbl
      End If
      
      !** Make sure that optstep has the correct step size for this image
      Call optstep_setstep(params%stepinfo,params%stepsize(i))

      !** Handle individual species separately
      Do spc = 1,nspc

        !** Skip this species if it's fixed
        If (config_isfixed(image(i,spc))) Cycle

        !** Get the number of atoms and molecules
        natms = config_getnatoms(image(i,:),spc)
        nmols = config_getnmoles(image(i,:),spc)

        !** Allocate space for the tangent vector(s)
        Allocate(tangent(natms,nmols), STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

        !** Evaluate tangent estimator
!        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!        Write(*,*) 'calculating tangents (image,spc): ',i,spc
        Call neb_calctangent(params,image,i,spc,tangent)

        !** Project out intra-image forces along tangent (don't do during preopt)
        If ((params%project_intra).And.(params%nmoves(2) >= params%preopt)) Then
          Do mol = 1,nmols
            Do atm = 1,natms
              Call storetop_projectout(imodel%ff%results, &
                  (/spc,mol,atm/),tangent(atm,mol))
            End Do
          End Do
        End If

        !** Evaluate inter-image forces, project and add to storage
        If (params%inter_on) Then
          Call neb_evalinter(params,image,i,spc,tangent,imodel,simcell)
        End If

        !** Deallocate space for the tangent vector(s)
        Deallocate(tangent, STAT=error)
        If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)    
      End Do

      !** Sum storage structure here?
!      Call storetop_fastsum(imodel%ff%results)

      !** Save the configuration, in case it later needs to be reset
      Call config_copysubset(image(i,:),saveimage(i,:),(/0,0,0/),(/0,0,0/))

      !** Move this image
      params%nmoves(i) = params%nmoves(i) + 1
      Call optstep_move(params%stepinfo,imodel,image(i,:),simcell,careful)
    End Do

  End Subroutine neb_movechain

  !----------------------------------------------------------------------------
  ! Evaluate, project and add the inter-image force between an image and
  ! its bounding images.  This routine operates only on a single species.
  ! Requires:  params -- NEB parameter data structure
  !            image -- the set of configuration images
  !            imageno -- image number
  !            spc -- species number
  !            tangent -- tangent vectors for each atom in species
  !            imodel -- interaction model
  !            simcell -- simulation cell structure
  !----------------------------------------------------------------------------
  Subroutine neb_evalinter(params,image,imageno,spc,tangent,imodel,simcell)
    Type(NEB_Info), Intent(InOut)                     :: params
    Type(AtMolCoords), Dimension(:,:), Intent(InOut)  :: image
    Integer, Intent(In)                               :: imageno,spc
    Type(VecType), Dimension(:,:), Intent(In)         :: tangent
    Type(Interaction_Model), Intent(InOut)            :: imodel
    Type(SimCell_Params), Intent(In)                  :: simcell

    Integer              :: nmols,natms,mol,atm,atype,error
    Logical              :: success
    Type(VecType)        :: atvec1
    Type(EnergyPlus)     :: results,opposite
    Type(EnergyPlus), Pointer   :: ptr
    Type(VecType), Dimension(2) :: atvecs2

    !** Initialize the temporary storage
    Call storebase_init(results,1,.False.)
    Call storebase_init(opposite,1,.False.)

    !** Get the number of atoms and molecules
    natms = config_getnatoms(image(imageno,:),spc)
    nmols = config_getnmoles(image(imageno,:),spc)

    !** Loop over atoms, handle interactions one at a time
    Do mol = 1,nmols
      Do atm = 1,natms

        !** Get atom type and coordinates 
        atype = config_getatype(image(imageno,:),spc,mol,atm)
        atvec1 = config_getr(image(imageno,:),(/spc,mol,atm/))
        atvecs2(1) = config_getrp(image(imageno-1,:),(/spc,mol,atm/))
        atvecs2(2) = config_getrp(image(imageno+1,:),(/spc,mol,atm/))

        !** Evaluate inter-image force
        success = pairmodel_mmint4(params%sprconst,results,opposite, &
            (/atype/),(/atvec1/),(/0.0_RDbl/),1,(/atype,atype/),atvecs2, &
            (/0.0_RDbl,0.0_RDbl/),2,simcell) 
        Call checkandstop(success,__FILE__,__LINE__, &
            ' inter-image spring force evaluation problem')
    
        !** Project out inter-image force perpendicular to tangent
        Call storebase_project(results,tangent(atm,mol))
    
        !** Add inter-image force to interactions storage (intra part)
        Call storetop_basicptrs(imodel%ff%results,.True.,.False., &
            (/spc,mol,atm/),(/0,0,0/),ptr)
        Call storebase_inc(ptr,results)

      End Do
    End Do

    !** Clean temporary storage
    Call storebase_clean(results)
    Call storebase_clean(opposite)

  End Subroutine neb_evalinter

  !----------------------------------------------------------------------------
  ! Calculate the tangent estimator for a single species of a given image.
  ! Returns unit vectors.
  ! Requires:  params -- NEB parameter data structure
  !            image -- the set of configuration images
  !            imageno -- image number
  !            spc -- species number
  !----------------------------------------------------------------------------
  Subroutine neb_calctangent(params,image,imageno,spc,tangent)
    Type(NEB_Info), Intent(InOut)                     :: params
    Type(AtMolCoords), Dimension(:,:), Intent(In)     :: image
    Integer, Intent(In)                               :: imageno,spc
    Type(VecType), Dimension(:,:), Intent(Out)        :: tangent

    Integer                :: i,j,atm,mol
    Integer                :: natms,nmols,curve
    Logical                :: upchain
    Real(kind=RDbl)        :: vforward,vback,vhere,delvmax,delvmin
    Type(VecType)          :: forward,back,vec
    Character(len=lstrLen) :: string

    !** Get the number of atoms and molecules
    natms = config_getnatoms(image(imageno,:),spc)
    nmols = config_getnmoles(image(imageno,:),spc)

    !** Get the energies 
    vback = params%totnrg(imageno-1)      !** V_i-1
    vhere = params%totnrg(imageno)        !** V_i
    vforward = params%totnrg(imageno+1)   !** V_i+1

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,'(a,3e14.4)') 'Energies (i-1, i, i+1): ',vback,vhere,vforward

    !** Establish the "curvature"
    If ((vforward > vhere) .And. (vhere > vback)) Then
      curve = 1   !** energy increasing in forward direction of chain
    Else If ((vforward < vhere) .And. (vhere < vback)) Then
      curve = -1  !** energy decreasing in forward direction of chain
    Else
      curve = 0
      upchain = .True.
      If (vforward < vback) upchain = .False.
      delvmax = Max(Abs(vforward - vhere), Abs(vback - vhere))
      delvmin = Min(Abs(vforward - vhere), Abs(vback - vhere))

      If ((Abs(delvmax) < 1.0e-6_RDbl) .And. (Abs(delvmin) < 1.0e-6_RDbl)) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            '  energy differences between images are too small'
        Stop
      End If
    End If

    !** Loop through the atoms and calculate the tangents
    Do atm = 1,natms
      Do mol = 1,nmols
        !** Calculate the appropriate unit separation vector(s)
        If (curve >= 0) Then    !** T+
          forward = config_getrp(image(imageno+1,:),(/spc,mol,atm/)) - &
              config_getrp(image(imageno,:),(/spc,mol,atm/))
          string = vector_display(forward,'e14.4')
        End If
        If (curve <= 0) Then    !** T-
          back = config_getrp(image(imageno,:),(/spc,mol,atm/)) - &
              config_getrp(image(imageno-1,:),(/spc,mol,atm/))
          string = vector_display(back,'e14.4')
        End If
    
        !** Calculate the tangent
        Select Case (curve)
        Case (-1)
          vec = back
        Case (1)
          vec = forward
        Case (0)
          If (upchain) Then
            vec = forward*delvmax + back*delvmin
          Else
            vec = forward*delvmin + back*delvmax
          End If
        End Select
        
        !** Normalize the tangent vector
        tangent(atm,mol) = vector_getunitvec(vec)

      End Do
    End Do

  End Subroutine neb_calctangent

  !----------------------------------------------------------------------------
  ! Calculate a distance estimator between the specified image and the previous
  ! image for statistical purposes.  The current estimator is the average 
  ! distance of each atom to corresponding atom in the previous image.
  ! Requires:  params -- NEB parameter data structure
  !            image -- the set of configuration images
  !            imageno -- image number
  !            spc -- species number
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function neb_calcdist(params,image,imageno,spc)
    Type(NEB_Info), Intent(InOut)                  :: params
    Type(AtMolCoords), Dimension(:,:), Intent(In)  :: image
    Integer, Intent(In)                            :: imageno,spc

    Integer                :: natms,nmols,atm,mol
    Real(kind=RDbl)        :: dist
    Type(VecType)          :: back

    !** Return zero if the first image is specified
    If (imageno == 1) Then
      neb_calcdist = 0.0_RDbl
      Return
    End If

    !** Get the number of atoms and molecules
    natms = config_getnatoms(image(imageno,:),spc)
    nmols = config_getnmoles(image(imageno,:),spc)

    !** Zero the distance estimator for this image,spc combination
    params%avgdist(imageno,spc) = 0.0_RDbl

    !** Loop through the atoms 
    Do atm = 1,natms
      Do mol = 1,nmols
        back = config_getrp(image(imageno,:),(/spc,mol,atm/)) - &
            config_getrp(image(imageno-1,:),(/spc,mol,atm/))
        dist = vector_getnorm(back)
        params%avgdist(imageno,spc) = params%avgdist(imageno,spc) + dist
      End Do
    End Do

    !** Average the distance estimator
    params%avgdist(imageno,spc) = params%avgdist(imageno,spc)/(natms*nmols)

    neb_calcdist = params%avgdist(imageno,spc)

  End Function neb_calcdist

  !----------------------------------------------------------------------------
  ! Returns the control file tags for the start and end configurations.
  ! Requires:  params -- NEB parameter data structure
  !            id -- 'START' or 'END'
  !----------------------------------------------------------------------------
  Function neb_configtag(params,id)
    Character(len=lstrLen)           :: neb_configtag
    Type(NEB_Info), Intent(InOut)    :: params
    Character(*), Intent(In)         :: id

    Select Case (ToUpper(id))
    Case ('START')
      neb_configtag = params%start_tag
    Case ('END')
      neb_configtag = params%end_tag
    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
          ' Could not interpret input string: ',Trim(id)
      Stop
    End Select 

  End Function neb_configtag

  !----------------------------------------------------------------------------
  ! Returns the number of images
  ! Requires:  params -- NEB parameter data structure
  !----------------------------------------------------------------------------
  Integer Function neb_nimages(params)
    Type(NEB_Info), Intent(InOut)     :: params

    neb_nimages = params%nimages

  End Function neb_nimages

  !-----------------------------------------------------------------------
  ! Dumps a modified image to an .xyz file.  It adds hydrogen atoms to
  ! each mobile species atom 1.0 Angstroms along the calculated tangent
  ! vector.  Useful for checking the tangent vectors.
  ! Requires:  params -- NEB parameter data structure
  !            image -- images data structure
  !            imageno -- image number to consider
  !            filename -- file to dump into
  !            comment -- comment for the xyz file
  !-----------------------------------------------------------------------
  Subroutine neb_modimage2xyz(params,image,imageno,filename,comment)
    Type(NEB_Info), Intent(InOut)                  :: params
    Type(AtMolCoords), Dimension(:,:), Intent(In)  :: image
    Integer, Intent(In)                            :: imageno
    Character(*), Intent(In)                       :: filename,comment

    Integer                  :: i,a,m,spc,n,lo,unitno,arraysize,error
    Integer                  :: nspc,nmols,natms,nfixed
    Type(VecType)            :: coord
    Integer, Dimension(20)   :: nums
    Type(XYZ_Entry), Dimension(:), Allocatable   :: entries
    Type(VecType), Dimension(:,:), Allocatable   :: tangent

    !** Cannot calculate tangents for 1st image, return
    If ((imageno == 1).Or.(imageno == params%nimages)) Then
      Write(0,'(1x,2a,i4,3a,i2)') __FILE__,": ",__LINE__, &
          ' WARNING: Cannot calculate tangents for 1st or last image'
      Write(0,'(1x,a)') '  Do not call neb_modimage2xyz with endpoint images'
      Return
    End If

    !** Get the number of species and the fixed species numbers
    nspc = Size(image,2)
    nfixed = str2seq(params%fixedspcs,nums)

    !** Count the number of atom spaces needed
    arraysize = 0
    Do spc = 1,nspc
      natms = config_getnatoms(image(imageno,:),spc)
      arraysize = arraysize + natms
      If (findint(nums(1:nfixed),spc) /= 0) Cycle
      arraysize = arraysize + natms
    End Do

    !** Size the entries array
    Allocate(entries(arraysize),STAT=error)
    If (error /= 0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Get the real coordinates and tangent vector "coordinates" for each spc
    n = 0
    Do spc = 1,nspc
      natms = config_getnatoms(image(imageno,:),spc)
      nmols = config_getnmoles(image(imageno,:),spc)

      !** Get the xyz entries for the real coordinates
      lo = n + 1
      n = n + config_setxyz(image(imageno,:),(/spc,0,0/),entries(lo:),.False.)

      If (findint(nums(1:nfixed),spc) /= 0) Cycle

      !** Allocate space for the tangent vector(s)
      Allocate(tangent(natms,nmols), STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)    

      !** Get the tangent vectors
      Call neb_calctangent(params,image,imageno,spc,tangent)

      !** Get the xyz entries for the tangent coordinates
      Do m = 1,nmols
        Do a = 1,natms
          coord = config_getrp(image(imageno,:),(/spc,m,a/)) + &
              0.5_RDbl*vector_getunitvec(tangent(a,m))
          n = n + 1
          If (n > arraysize) Then
            Write(0,'(1x,2a,i4,3a,i2)') __FILE__,": ",__LINE__, &
                ' entries array was not allocated large enough'
            Stop
          End If
          entries(n) = visxyz_make(coord,'H ')
        End Do
      End Do

      !** Deallocate the tangents
      Deallocate(tangent,STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__)
    End Do

    !** Dump the file
    Call visxyz_dump(entries,filename,'f8.3',comment)

    !** Deallocate the entries array
    Deallocate(entries,STAT=error)
    If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__)

  End Subroutine neb_modimage2xyz

  !-----------------------------------------------------------------------
  ! Dumps the current images to a single xyz movie file
  ! Requires:  params -- NEB parameter data structure
  !            images -- set of configurations making up the images
  !            filename -- file to dump into
  !-----------------------------------------------------------------------
  Subroutine neb_images2xyz(params,images,filename)
    Type(NEB_Info), Intent(In)                     :: params
    Type(AtMolCoords), Dimension(:,:), Intent(In)  :: images
    Character(*), Intent(In)                       :: filename

    Integer                  :: i,unitno
    Character(len=strLen)    :: comment

    Do i = 1,params%nimages
      Write(comment,'(a,i3)') ' NEB chain, image number: ',i
      Call config_config2xyz(images(i,:),filename,.False.,comment)
    End Do

    !** Close the movie file
    unitno = file_close(filename)

  End Subroutine neb_images2xyz

  !-----------------------------------------------------------------------
  ! Dumps the current images to a single xyz movie file WITH tangent 
  ! vectors display as hydrogen atoms.
  ! Requires:  params -- NEB parameter data structure
  !            images -- set of configurations making up the images
  !            filename -- file to dump into
  !-----------------------------------------------------------------------
  Subroutine neb_modimages2xyz(params,images,filename)
    Type(NEB_Info), Intent(InOut)                  :: params
    Type(AtMolCoords), Dimension(:,:), Intent(In)  :: images
    Character(*), Intent(In)                       :: filename

    Integer                  :: i,unitno
    Character(len=strLen)    :: comment

    Write(comment,'(a,i3)') ' NEB chain, image number: ',i
    Call config_config2xyz(images(1,:),filename,.False.,comment)

    Do i = 2,params%nimages-1
      Write(comment,'(a,i3)') ' NEB chain with tangents, image number: ',i
      Call neb_modimage2xyz(params,images,i,filename,comment)
    End Do

    Write(comment,'(a,i3)') ' NEB chain, image number: ',i
    Call config_config2xyz(images(params%nimages,:),filename,.False.,comment)

    !** Close the movie file
    unitno = file_close(filename)

  End Subroutine neb_modimages2xyz

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to specified unit
  ! Requires:  unitno -- unit number to dump into
  !----------------------------------------------------------------------------
  Subroutine neb_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(3a)') shortdashedline,Trim(neb_tag),shortdashedline
    Write(unitno,'(a,t30,a)') 'Integer','# Number of images between end-points'
    Write(unitno,'(a,t30,a)') 'Integer List','# Species numbers to freeze'
    Write(unitno,'(a,t30,a)') 'Integer','# Number of pre-optimization steps'
    Write(unitno,'(a,t30,2a)') 'Real','# Harmonic spring constant ', &
        'between image atoms'
    Write(unitno,'(a,t30,2a)') 'Real','# Maximum allowed one-step ', &
        'increase in energy (kJ/mol) '
    Write(unitno,'(a,t30,2a)') 'String','# Project out intra-image ', &
        'interactions along tangent (YES/NO)'
    Write(unitno,'(a,t30,2a)') 'String','# Primary step type identifier'
    Write(unitno,'(a,t30,2a)') 'String','# control file tag for ', &
        'starting configuration'
    Write(unitno,'(a,t30,2a)') 'String','# READ or INTERPOLATE  ', &
        'intermediate configurations'
    Write(unitno,'(a,t30,2a)') 'String','# control file tag for ', &
        'ending configuration'

  End Subroutine neb_sampleCF

  !--------------------------------------------------------------------------
  ! Display the NEB statistics for the current simulation
  ! Requires:  params -- GCMC simulation parameters
  !            image -- the set of configuration images
  !            indent -- no. of spaces from the left margin
  !            optunit -- optional unit number for display
  !--------------------------------------------------------------------------
  Subroutine neb_displaystats(params,image,indent,optunit)
    Type(NEB_Info), Intent(InOut)                  :: params
    Type(AtMolCoords), Dimension(:,:), Intent(In)  :: image
    Integer, Intent(In)                            :: indent
    Integer, Intent(In), Optional                  :: optunit

    Integer                      :: i,n,unitno,spc
    Real(kind=RDbl)              :: delnrg,fraction,avgdist
    Character(len=indent)        :: blank
    Character(len=strLen)        :: string1,string2,string3,string4,string5
    Integer, Dimension(10)       :: nums

    blank = Repeat(' ',indent)

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    nums = 0
    n = str2seq(params%fixedspcs,nums)

    Write(unitno,'(2a)') blank, dashedline
    Write(unitno,'(2a)') blank, "The NEB Statistics:"
    Write(unitno,'(2a,5(2x,a12))') blank, 'image','nrg (kJ/mol)', &
        'nrg change','step size  ','success rate','avg dist'
    Do i = 1,params%nimages
      !** Accumulate statistics for display
      delnrg = params%totnrg(i) - params%lastnrg(i)
      string1 = real2str(params%totnrg(i),8)
      string2 = real2str(delnrg,8)
      string3 = real2str(params%stepsize(i),8)
      fraction = 0.0_RDbl
      If (params%nmoves(i) /= 0) Then
        fraction = 1.0_RDbl - Real(params%nresets(i))/params%nmoves(i)
      End If
      string4 = real2str(fraction,5)

      !** Calculate average distance estimator from species average
      avgdist = 0.0_RDbl
      n = 0
      Do spc = 1,params%nspc
        If (findint(nums,spc) /= 0) Cycle
        n = n + 1
        avgdist = avgdist + neb_calcdist(params,image,i,spc)
      End Do
      If (n /= 0) avgdist = avgdist/Real(n)
      string5 = real2str(avgdist,6)

      !** Dump display for this image
      Write(unitno,'(a,i3,5(2x,a12))') blank,i,Trim(string1),Trim(string2), &
          Trim(string3),Trim(string4),Trim(string5)
    End Do

  End Subroutine neb_displaystats

  !-----------------------------------------------------------------------
  ! Displays the data in the NEB data structure
  ! Requires:  params -- NEB parameter data structure
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !-----------------------------------------------------------------------
  Subroutine neb_display(params,indent,unit)
    Type(NEB_Info), Intent(In)    :: params
    Integer, Intent(In)           :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2

    blank = Repeat(' ',indent)

    Write(unit,'(2a)') blank, dashedline
    Write(unit,'(2a)') blank,'Nudged Elastic Band Method Parameters:'

    string1 = int2str(params%nimages)
    Write(unit,'(3a)') blank,' Number of images: ',Trim(string1)
    If (params%fixedspcs .ne. '0') Then
      Write(unit,'(3a)') blank,' Fix species number(s): ',Trim(params%fixedspcs)
    End If
    string1 = int2str(params%preopt)
    Write(unit,'(3a)') blank,' Number of pre-optimization steps: ',Trim(string1)
    If (params%inter_on) Then
      string1 = real2str(params%spring_const,5)
      Write(unit,'(3a)') blank,' Spring constant (kcal/mol): ',Trim(string1)
    Else
      Write(unit,'(3a)') blank,' Inter-image interactions OFF'
    End If
    If (params%project_intra) Then
      Write(unit,'(3a)') blank,' Intra-image interaction projected out ',&
          'along tangent'
    Else
      Write(unit,'(3a)') blank,' Intra-image interaction NOT projected out ',&
          'along tangent'
    End If
    Write(unit,'(3a)') blank,' Intial chain images produced by: ', &
        Trim(params%config_origin)

    !** Display optimization step information
    Call optstep_display(params%stepinfo,indent+1,unit)

    Write(unit,'(2a)') blank,' Starting image statistics:'
    Write(unit,'(2x,2a,4(2x,a12))') blank, 'image','nrg (kJ/mol)', &
        'step size'
    Do i = 1,params%nimages
      Write(string1,'(f9.2)') params%totnrg(i)
      If (Abs(Log10(Abs(params%totnrg(i)))) > 5.0_RDbl) Then
        string1 = real2str(params%totnrg(i),9)
      End If
      string2 = real2str(params%stepsize(i),6)
      Write(unit,'(a,2x,i3,4(2x,a12))') blank,i,Trim(string1),Trim(string2)
    End Do

  End Subroutine neb_display 

  !-----------------------------------------------------------------------
  ! Cleans the NEB data structure
  ! Requires:  params -- NEB parameter data structure
  !-----------------------------------------------------------------------
  Subroutine neb_clean(params)
    Type(NEB_Info), Intent(InOut)    :: params

    Integer               :: i,j,error

    Deallocate(params%totnrg, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    Deallocate(params%lastnrg, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)

    Deallocate(params%nmoves, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    Deallocate(params%nresets, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    Deallocate(params%avgdist, STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)

    If (params%inter_on) Then
      Do i = 1,Size(params%sprconst,1)
        Do j = 1,Size(params%sprconst,2)
          Call pairmodel_clean(params%sprconst(i,j))
        End Do
      End Do

      Deallocate(params%sprconst, STAT=error)
      If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)
    End If

    Call optstep_clean(params%stepinfo)

  End Subroutine neb_clean
  
End Module neb



