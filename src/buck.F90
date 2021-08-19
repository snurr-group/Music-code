!------------------------------------------------------------------------------
! This module contains the structure for defining Buckingham parameters
!------------------------------------------------------------------------------
Module buck

  Use defaults, Only: RDbl, strLen, zero, lstrlen
  Use utils, Only: split, toupper, toreal
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use storebase, Only: EnergyPlus, storebase_inc, storebase_nderivs, &
      storebase_disp

  Implicit None
  Save

  Private
  Public :: Buckingham_Params, buck_idstring, buck_calc_interaction, &
      buck_init, buck_isinit, buck_displayParams, buck_displayCutoffs, &
      buck_getcutoff, Assignment(=), buck_snglint, buck_chkequal

  Type Buckingham_Params
    Character(len=strLen)     :: atom_name1,atom_name2
    Real(kind=RDbl)           :: A, Rho, C, D
    Real(kind=RDbl)           :: hicut, locut, ncut ! High, low, and neighbor
    Real(kind=RDbl)           :: bar1, bar2, bar3
    Character(len=100)        :: line
    Character(len=80)         :: reference
    Logical                   :: off, initialized
  End Type Buckingham_Params

  Character(len=strLen), Parameter                 :: buck_idstring = 'BUCK'

  Interface Assignment(=)
    Module Procedure buck_copy
  End Interface

  Interface buck_calc_interaction
    Module Procedure buck_snglInteraction
    Module Procedure buck_multInteraction
  End Interface

Contains

  !----------------------------------------------------------------------------
  ! Initialize the interaction information from the file line
  !----------------------------------------------------------------------------
  Subroutine buck_init(params,line)
    Type(Buckingham_Params), Intent(INOUT)    :: params
    Character(*), Intent(IN)                  :: line

    Integer                                   :: nfields,n,i
    Character(len=strLen), Dimension(10)      :: fields,chunks

    !** Split the line into individual fields
    params%line = line
    nfields = split(line,fields)
    params%atom_name1 = Trim(fields(1))
    params%atom_name2 = Trim(fields(2))
    params%off = .False.

    !** Zero the coefficients
    params%A = 0.0_RDbl
    params%C = 0.0_RDbl
    params%D = 0.0_RDbl

    Do i = 4,nfields
      n = split(fields(i),chunks,'@')
      Select Case(Toupper(chunks(1)))

      Case('OFF')
        params%off = .True.
      Case('A')
        params%A = toreal(chunks(2))
      Case('RHO')
        params%rho = toreal(chunks(2))
      Case('C')
        params%C = toreal(chunks(2))
      Case('D')
        params%D = toreal(chunks(2))
      Case('HICUT')
        params%hicut = toreal(chunks(2))
      Case('LOCUT')
        params%locut = toreal(chunks(2))
      Case('NCUT')
        params%ncut = toreal(chunks(2))
      Case Default
        Write(0,'(2a,i4,a,2i4)') __FILE__,": ",__LINE__, &
             'Unable to identify BUCK interaction string ',trim(chunks(1))
        Stop
      End Select
    EndDo

    !** Set the initialized flag
    params%initialized = .True.

  End Subroutine buck_init

  !-------------------------------------------------------------------
  ! Initialize params by copying from an existing params
  ! Requires: newparams -- new params, to be initialized
  !           oldparams -- old params
  !-------------------------------------------------------------------
  Subroutine buck_copy(newparams,oldparams)
    Type(Buckingham_Params), Intent(Out)      :: newparams
    Type(Buckingham_Params), Intent(In)       :: oldparams

    newparams%atom_name1 = oldparams%atom_name1
    newparams%atom_name1 = oldparams%atom_name2
    newparams%A = oldparams%A
    newparams%Rho = oldparams%Rho
    newparams%C = oldparams%C
    newparams%D = oldparams%D
    newparams%hicut = oldparams%hicut
    newparams%locut = oldparams%locut
    newparams%ncut = oldparams%ncut
    newparams%bar1 = oldparams%bar1
    newparams%bar2 = oldparams%bar2
    newparams%bar3 = oldparams%bar3
    newparams%line = oldparams%line
    newparams%reference = oldparams%reference
    newparams%off = oldparams%off
    newparams%initialized = oldparams%initialized

  End Subroutine buck_copy

  !----------------------------------------------------------------------------
  ! Returns the value of the initialized flag
  !----------------------------------------------------------------------------
  Logical Function buck_isinit(params)
    Type(Buckingham_Params), Intent(In) :: params

    buck_isinit = params%initialized
  End Function buck_isinit

  !----------------------------------------------------------------------------
  ! Returns the value of the high end cutoff. If the NEIGHBOR flag is 
  ! present, check to see if the NEIGHBOR cutoff should be returned
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function buck_getcutoff(params,neighbor)
    Type(Buckingham_Params), Intent(In) :: params
    Logical, Optional, Intent(In)       :: neighbor

    !** Default value
    buck_getcutoff = params%hicut

    !** Check if we need to return the neighbor cutoff
    If (Present(neighbor)) Then
      If (neighbor) Then
        buck_getcutoff = params%ncut
      End If
    End If
  End Function buck_getcutoff

  !----------------------------------------------------------------------------
  ! Calculate the interaction 
  !----------------------------------------------------------------------------
  Logical Function buck_snglInteraction(params,sepvec,pot,fvec)
    Type(Buckingham_Params), Intent(IN)       :: params
    Type(VecType), Intent(IN)                 :: sepvec
    Real(kind = RDbl), Intent(OUT)            :: pot
    Type(VecType), Intent(OUT), Optional      :: fvec
    Real(kind=RDbl)                           :: r6, r8, r, A, C, rho, D
    Type(VecType)                             :: unitvec

    !** Zero the potential and the force
    pot = zero
    If (Present(fvec)) Then
      fvec = zero
    End If
    
    !** Set the return value to be true
    buck_snglInteraction = .True.

    !** Check to see if the interaction is on
    If (params%off) Return

    !** Get the scalar separation
    r  = mag(sepvec)
    r6 = r**6
    r8 = r6*r*r

    !** Make the variables neater
    A = params%A
    rho = params%rho
    C = params%C
    D = params%D

    !** Check for cutoff. If too low, set the return value to be false
    If (r > params%hicut) Return
    If (r < params%locut) Then
      buck_snglInteraction = .False.
      Return
    End If

    !** Get the potential
    pot = A*Exp(-r/rho) - C/r6 - D/r8

    !MDEBUG
!!$    Write(0,'(a,i5,a,5e15.5)') __FILE__,__LINE__,"buck params ", &
!!$        A*Exp(-r/rho), C/r6, 1/rho, r, 1.0_RDbl/(r**2)
!!$    Write(0,'(a,5e15.5)') "sepvec = ",sepvec
!!$    Write(0,'(a,5e15.5)') "buck force  ",A/rho*Exp(-r/rho)*r,6.0*C/r6

    !** Get the force if necessary
    If (Present(fvec)) Then
      !** Get the unit vector
      unitvec = sepvec/r

      !** Calculate the force. Force = -grad(potential)
      fvec = (-1.0_RDbl)*unitvec*((-A/rho)*Exp(-r/rho) + 6.0_RDbl*C/(r6*r) &
          + 8.0_RDbl*D/(r8*r))
    End If

  End Function buck_snglInteraction

  !----------------------------------------------------------------------------
  ! Calculate the interaction 
  !----------------------------------------------------------------------------
  Logical Function buck_multInteraction(params,apairs,sepvec,pot,fvec)
    Type(Buckingham_Params), Dimension(:,:), Intent(IN) :: params
    Integer, Dimension(:,:), Intent(In)       :: apairs
    Type(VecType), Intent(In), Dimension(:)   :: sepvec
    Real(kind = RDbl), Intent(OUT)            :: pot
    Type(VecType), Intent(OUT), Dimension(:), Optional :: fvec
    Real(Kind=Rdbl) :: r, r6, r8, A, rho, C, D
    Integer         :: i, a1,a2
    Type(VecType) :: unitvec
    
    !** Zero the potential and the force
    pot = zero
    If (Present(fvec)) fvec = VecType(0.0_RDbl)
    
    !** Set the return value to be true
    buck_multInteraction = .True.

    !** Loop over each of the separation vectors
    Do i = 1, Size(sepvec)

      !** Get the atom indicies
      a1 = apairs(i,1)
      a2 = apairs(i,2)

      !** Check to see if the nteraction is off
      If (params(a1,a2)%off) Cycle

      !** Get the scalar separation
      r  = mag(sepvec(i))
      r6 = r**6
      r8 = r6*r*r

      !** Check for cutoff. If too low, set the return value to be false
      If (r > params(a1,a2)%hicut) Cycle
      If (r < params(a1,a2)%locut) Then
        buck_multInteraction = .False.
        Return
      End If

      !** Let's simplify the following equations by renaming variables
      A = params(a1,a2)%A
      rho = params(a1,a2)%rho
      C = params(a1,a2)%C
      D = params(a1,a2)%D
      
      !** Get the potential
      pot = pot + A * Exp(-r/rho) - C/r6 - D/r8

      !** Get the force if necessary. Force = -grad(potential)
      !** Note that sepvec(i)/r = unitvec
      If (Present(fvec)) Then
        unitvec = sepvec(i) * (1.0_RDbl/r)
        fvec(i) = (-1.0_RDbl)*unitvec*((-A/rho)*Exp(-r/rho) &
            + 6.0_RDbl*C/(r6*r) + 8.0_RDbl*D/(r8*r))
      End If

    End Do

  End Function buck_multInteraction

  !---------------------------------------------------------------------------
  ! Calculate the interaction between a set of atomic points.  Returns False
  ! if any of the seperation vectors fall within the lower cutoff.  This 
  ! indicates atoms have gotten too close together and that the potential 
  ! would be extremely large.  The returned derivative components are those
  ! for the first atom in the separation vector defined by (atom1 - atom2).
  ! The optional second results structure must be passed if forces are to be
  ! calculated.
  ! Requires:  params -- Buckingham interaction parameters
  !            sepvec -- atom-atom separation vector
  !            results -- resultant energies and derivatives
  !            results2 -- resultant energies and opposite forces
  !---------------------------------------------------------------------------
  Logical Function buck_snglint(params,sepvec,results,results2)
    Type(Buckingham_Params), Intent(In)       :: params
    Type(VecType), Intent(In)                 :: sepvec
    Type(EnergyPlus), Intent(InOut)           :: results
    Type(EnergyPlus), Intent(InOut), Optional :: results2

    Real(kind=RDbl)           :: r8,r6,r,nrg
    Real(kind=RDbl)           :: A, C, rho, D
    Type(VecType)             :: force,unitvec

    !** Set the return value to be true
    buck_snglint = .True.

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'WARNING: not tested'
    Stop

    !** Check to see if the interaction is on
    If (params%off) Return

    !** Get the scalar separation
    r  = mag(sepvec)
    r6 = r**6
    r8 = r6*r*r

    !** Make the variables neater
    A = params%A
    rho = params%rho
    C = params%C
    D = params%D

    !MDEBUG
!!$    Write(0,'(a,i5,a,5e15.5)') __FILE__,__LINE__,"buck params ", &
!!$        A*Exp(-r/rho), C/r6, 1/rho, r, 1.0_RDbl/(r**2)
!!$    Write(0,'(a,5e15.5)') "sepvec = ",sepvec
!!$    Write(0,'(a,5e15.5)') "buck force  ",A/rho*Exp(-r/rho)*r,6.0*C/r6

    !** Check for cutoff. If too low, set the return value to be false
    If (r > params%hicut) Return
    If (r < params%locut) Then
      buck_snglint = .False.
      Return
    End If

    !** Get the potential
    nrg = A*Exp(-r/rho) - C/r6 - D/r8
    Call storebase_inc(results,nrg)
!    If (Present(results2)) Call storebase_inc(results2,nrg)

    !** Get the force if necessary
    If (storebase_nderivs(results) > 0) Then
      unitvec = sepvec/r
      force = (-1.0_RDbl)*unitvec* &
          ((-A/rho)*Exp(-r/rho) + 6.0_RDbl*C/(r6*r) + 8.0_RDbl*D/(r8*r))

      Call storebase_inc(results,force)
      Call storebase_inc(results2,(force*(-1.0_RDbl)))
    End If

    !** Calculate select second and third derivatives (HOT)
    If (storebase_nderivs(results) > 1) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Stop
      !PLEASE FILL THIS IN
    End If

  End Function buck_snglint

  !----------------------------------------------------------------------------
  ! Checks to see if two Buckingham data structures contain the same parameters
  ! Requires:  params1 -- 1st params
  !            params2 -- 2nd params
  !            tolerance -- optional tolerance for checking
  !----------------------------------------------------------------------------
  Logical Function buck_chkequal(params1,params2,tolerance)
    Type(Buckingham_Params), Intent(In)       :: params1,params2
    Real(kind=RDbl), Intent(In), Optional     :: tolerance

    Real(kind=RDbl)  :: tol,diff

    tol = 1.0e-6_RDbl
    If (Present(tolerance)) tol = tolerance

    buck_chkequal = .True.

    diff = Abs(params1%A - params2%A)
    If (diff > tol) buck_chkequal = .False.

    diff = Abs(params1%Rho - params2%Rho)
    If (diff > tol) buck_chkequal = .False.

    diff = Abs(params1%C - params2%C)
    If (diff > tol) buck_chkequal = .False.

    diff = Abs(params1%D - params2%D)
    If (diff > tol) buck_chkequal = .False.

  End Function buck_chkequal

  !----------------------------------------------------------------------------
  ! Display the pertinent information from the BUCK structure in string format
  ! Requires:  obj -- the parameters structure
  !----------------------------------------------------------------------------
  Function buck_displayParams(obj)
    Character(len=lstrLen)                :: buck_displayParams
    Type(Buckingham_Params), Intent(In)   :: obj
    
    If (obj%off) Then
      Write(buck_displayParams,'(a)') 'OFF'
    Else
      Write(buck_displayParams,'(3(a,f10.3,1x,a,2x))') "A: ",obj%A, &
          "kcal/mol","rho: ",obj%rho,"Ang","C: ",obj%C,"??"
    End If

  End Function buck_displayParams

  !----------------------------------------------------------------------------
  ! Display the pertinent information from the BUCK structure
  ! Requires:  obj -- the Buckingham parameters structure
  !            dunit -- display unit
  !            indent -- number of spaces to indent
  !----------------------------------------------------------------------------
  Subroutine buck_display(obj,dUnit,indent)
    Type(Buckingham_Params), Intent(IN)   :: obj
    Integer, Intent(In)                   :: dUnit
    Integer, Intent(In)                   :: indent

    Character(len=indent)        :: spaces
    
    !** Make the spaces
    spaces = Repeat(' ',indent)

    !** Diplay the information to the given display unit
    If (obj%off) Then
      Write(dUnit,'(2a)') spaces,'OFF'
    Else
      Write(dUnit,'(2a,f12.3,a)') spaces,'A   = ',obj%A,' kcal/mol'          
      Write(dUnit,'(2a,f12.3,a)') spaces,'rho = ',obj%rho,' Ang'
      Write(dUnit,'(2a,f12.3,a)') spaces,'C   = ',obj%C,' kcal/mol/Ang^6'
      Write(dUnit,'(2a,f12.3,a)') spaces,'D   = ',obj%D,' kcal/mol/Ang^8'
    End If
  End Subroutine buck_display

  !----------------------------------------------------------------------------
  ! Display the cutoffs for Buck
  ! Requires:  obj -- the Buckingham parameters structure
  !----------------------------------------------------------------------------
  Character(lstrLen) Function buck_displayCutoffs(obj)
    Type(Buckingham_Params), Intent(IN)      :: obj

    If (obj%off) Then
      buck_displayCutoffs = ''
    Else
      Write(buck_displayCutoffs,'(2(a,f8.3,1x,a,7x))') "High : ", &
          obj%hicut,"Ang","Low : ",obj%locut,"Ang" 
    End If

  End Function buck_displayCutoffs

End Module buck
