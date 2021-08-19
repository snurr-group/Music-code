!------------------------------------------------------------------------------
! This module handles the Umbrella/Importance sampling.  Umbrella sampling is
! importance sampling intended to bridge two distributions.  This module is 
! a close extension to the interact module.
!
! The operating equation is the general importance sampling formulation:
!           <A/w exp[-beta V_o]>_w
!   <A>_o = ----------------------
!           <1/w exp[-beta V_o]>_w
! where A is the observable and <>_x denotes an ensemble average in the
! ensemble given by the 'x' distribution.  The subscript 'o' denotes the 
! real or reference ensemble, the one in which the quantity will be averaged.
! The 'w' denotes the weighting distribution, P(x) = w(x)dx, that produces 
! the ensemble in which the quantity will be sampled.
!
! In the special case that the 'w' ensemble is simply another forcefield,
! the above equation reduces to:
!           <A exp[-beta(V_o - V_w)]>_w       1/N sum_i A_i*weight_i
!   <A>_o = ---------------------------    = ------------------------
!            <exp[-beta(V_o - V_w)]>_w         normalization factor
! At the moment, this is the only formulation supported
!
! Important Routines:
!   umbrella_init -- initializes the umbrella sampling part of interaction model
!   umbrella_getwt -- returns the weighting and normalization factors
!
! Needed Improvements:
! 1) Add support for arbitrary weighting functions
!------------------------------------------------------------------------------

Module umbrella

  Use defaults, Only: RDbl,strLen,lstrLen,zero,one,dbgflag, &
      caltoj,kcalmole_kb,default_MAX_EXP   
  Use utils, Only: filesrchstr,stripcmnt,toupper,allocErrDisplay, &
      int2str,real2str,deallocErrDisplay,checkandstop,getdepth
  Use file, Only: file_open
  Use config, Only: AtMolCoords
  Use simcell, Only: SimCell_Params
  Use storestats, Only: Species_Stats
  Use forcefield, Only: Forcefield_Info, forcefield_display, forcefield_init, &
      forcefield_initstore
  Use storetop, Only: Forcefield_Results,storetop_display
  Use store, Only: Store_Level_Pair, store_idbranch, store_sum, store_disp
  Use storebase, Only: EnergyPlus, storebase_copy, storebase_init, &
      storebase_scalarmult, storebase_totnrg

  Implicit None
  Save

  Private 
  Public :: Umbrella_Params, umbrella_init, umbrella_idstring, &
      umbrella_getwt, umbrella_value, umbrella_display, umbrella_clean

  !** Stores all information about and from interaction model
  Type Umbrella_Params
    Integer                    :: nsamples  !** number of samples 
    Integer                    :: lastiter  !** last iteration number
    Real(kind=RDbl)            :: rti,TKelvin
    Real(kind=RDbl)            :: normfactor,normsum,weight
  End Type Umbrella_Params

  Character(len=strLen), Parameter :: umbrella_idstring = "UMBRELLA"

Contains
  !---------------------------------------------------------------------
  ! Initializes the Umbrella/Importance sampling parameters
  ! Requires:  params -- umbrella sampling parameters
  !            ctrl_filename -- the control filename to take info from
  !            species -- species data structure
  !            ffid -- forcefield IDs to initialize
  !---------------------------------------------------------------------
  Subroutine umbrella_init(params,ctrl_filename,species,ffid)
    Type(Umbrella_Params), Intent(Out)              :: params
    Character(*), Intent(In)                        :: ctrl_filename
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Character(*), Dimension(:), Intent(Out)         :: ffid

    Integer                   :: unit
    Character(len=lstrLen)    :: line

    !** Read the control file information
    unit = file_open(ctrl_filename,110)
    Read(unit,'(a)') line
    ffid(1) = stripcmnt(line)   !** first forcefield identifier
    Read(unit,'(a)') line
    ffid(2) = stripcmnt(line)   !** second forcefield identifier

    !** Set the defaults
    params%nsamples = 0
    params%lastiter = 0
    params%normfactor = one
    params%normsum = 0.0_RDbl
    params%weight = 0.0_RDbl

  End Subroutine umbrella_init

  !---------------------------------------------------------------------
  ! Calculate the weighting factor for the current configuration and
  ! the normalization factor needed to obtain averages.  It is
  ! assumed that the first value is the 'sampling' forcefield output and 
  ! the second is the 'real' output.  It requires the iteration number
  ! as input to be certain that points are not counted more than once
  ! in the calculation of the normalization factor.
  ! Requires:  params -- umbrella sampling parameters
  !            ffout -- energy only outputs from forcefield evaluations
  !            iter -- iteration number
  !            weight -- weighting factor for current configuration
  !            normfactor -- normalization factor for averaging
  !---------------------------------------------------------------------
  Subroutine umbrella_getwt(params,ffout,iter,weight,normfactor)
    Type(Umbrella_Params), Intent(InOut)       :: params
    Real(kind=RDbl), Dimension(2), Intent(In)  :: ffout
    Integer, Intent(In)                        :: iter
    Real(kind=RDbl), Intent(Out)               :: weight,normfactor

    Real(kind=RDbl)        :: deltaE,utemp

    !** Compare the iteration number to the last recorded number and act
    If (iter == params%lastiter) Then
      weight = params%weight
      normfactor = params%normfactor
    Else If (iter < params%lastiter) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          ' umbrella_getwt called for a past configuration'
      Stop      
    End If

    !** Increment and set counters
    params%lastiter = iter
    params%nsamples = params%nsamples + 1

    !** Calculate the energy difference
    deltaE = ffout(2) - ffout(1)

    !** Calculate the exponent for the weighting factor
    utemp = -1.0_RDbl*params%rti * deltaE

    !** Avoid overflow/underflow errors during exponentiation
    If (utemp > default_MAX_EXP) utemp = default_MAX_EXP   
    If (utemp < -default_MAX_EXP) utemp = -default_MAX_EXP

    !** Calculate the weighting factor
    weight = Exp(utemp)

    !** Calculate the normalization factor
    params%normsum = params%normsum + weight
    normfactor = params%normsum/params%nsamples

    !** Store the weighting and normalization factors
    params%weight = weight
    params%normfactor = normfactor

  End Subroutine umbrella_getwt

  !---------------------------------------------------------------------
  ! Transform a collection of outputs from the forcefield evaluations
  ! into an instantaneous value and a normalization factor.  The outputs
  ! from the forcefield evaluations are EnergyPlus structures.  It is
  ! assumed that the first is the 'sampling' forcefield output and the
  ! second is the 'real' output.  It requires the iteration number
  ! as input to be certain that points are not counted more than once
  ! in the calculation of the normalization factor.
  ! Requires:  params -- umbrella sampling parameters
  !            ffout -- outputs from forcefield evaluations
  !            iter -- iteration number
  !            info -- returned instanteous value(s) PRE-INITIALIZED!
  !            normfactor -- normalization factor for averaging
  ! NOTE: if called with iter==0, then it will return the value of the
  ! first (sampling) forcefield interactions
  !---------------------------------------------------------------------
  Subroutine umbrella_value(params,ffout,iter,info,normfactor)
    Type(Umbrella_Params), Intent(InOut)        :: params
    Type(EnergyPlus), Dimension(:), Intent(In)  :: ffout
    Integer, Intent(In)                         :: iter
    Type(EnergyPlus), Intent(Out)               :: info
    Real(kind=RDbl), Intent(Out)                :: normfactor

    Real(kind=RDbl)        :: weight,nrg1,nrg2

    !** Simply return the sampling forcefield interactions, if iter==0
    If (iter == 0) Then
      Call storebase_copy(info,ffout(1))
      Return
    End If

    !** Calculate the total energies (including intramolecular)
    nrg1 = storebase_totnrg(ffout(1),.True.)
    nrg2 = storebase_totnrg(ffout(2),.True.)

    !** Calculate the weighting and normalization factors
    Call umbrella_getwt(params,(/nrg1,nrg2/),iter,weight,normfactor)

    !** Copy the 'real' values into the output, then scale with weighting factor
    Call storebase_copy(info,ffout(2))
    Call storebase_scalarmult(info,weight)

  End Subroutine umbrella_value

  !---------------------------------------------------------------------
  ! Displays the parameters
  ! Requires:  params -- umbrella sampling parameters
  !            indent -- no. of spaces from the left margin
  !            unitno -- output unit number
  !---------------------------------------------------------------------
  Subroutine umbrella_display(params,indent,unitno)
    Type(Umbrella_Params), Intent(In) :: params
    Integer, Intent(In)               :: indent
    Integer, Intent(In)               :: unitno

    Integer                    :: unit
    Character(len=indent)      :: blank
    Character(len=strLen)      :: string1,string2

    blank = Repeat(' ',indent)

    string1 = real2str(params%weight,10)
    string2 = real2str(params%normfactor,10)
    Write(unitno,'(3a,2x,a)') blank,'Umbrella sampling wt,normfactor: ', &
        Trim(string1),Trim(string2)
    string1 = int2str(params%lastiter)
    string2 = int2str(params%nsamples)
    Write(unitno,'(a,1x,3a)') blank,'last iteration, nsamples: ', &
        Trim(string1),Trim(string2)

  End Subroutine umbrella_display
    
  !---------------------------------------------------------------------
  ! Clean the parameters
  ! Requires:  params -- umbrella sampling parameters
  !---------------------------------------------------------------------
  Subroutine umbrella_clean(params)
    Type(Umbrella_Params), Intent(InOut)      :: params

  End Subroutine umbrella_clean

End Module umbrella
