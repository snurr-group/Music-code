!------------------------------------------------------------------------------
! This module contains the routines specific to initializing and actually 
! calculating LENNARD-JONES INTERACTIONS.  Interactions can be calculated
! either between a single pair (lj_calc_snglInteraction) or for a set 
! of pairs (lj_calc_multInteraction).  Higher order terms 
! (HOT), meaning first and second derivatives, are also available.  The
! data structure (LJ_Params) contains all the potential parameterization
! information.  It is initialized using file information in "lj_init"
! 
! Input units from file MUST be:
! Sigma: Angstroms 
! Epsilon: Kelvin (implies that the quantity is epsilon/k_boltz)
! 
! Ouput units are in kcal/mol (for energy) and kcal/mol/Angstrom (for forces)
! 
! Needed Improvements: 
! 1) addition of smoothing features from the old code.
! 2) use of the "reference" parameter
!------------------------------------------------------------------------------

Module lj

  Use defaults, Only: RDbl, strLen, lstrLen, LOW_DIST_CUTOFF, kcalmole_kb, &
      HIGH_DIST_CUTOFF, xlstrlen, SIZE_OF_HOT
  Use utils, Only: toupper, split, toreal, getlineswpair, real2str, int2str
  Use file, Only: file_getunit
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnormsq, vector_display, mag
  Use atom, Only: atom_getmass
  Use storebase, Only: EnergyPlus, storebase_inc, storebase_nderivs, &
      storebase_disp

  Implicit None
  Save

  Private
  Public :: LJ_Params, lj_displayparams, lj_calc_interaction, lj_idstring, &
      lj_isinit, lj_displaySigEps, lj_calc_mix, lj_calc_ab, Assignment(=), &
      lj_init, lj_displayCutOffs, lj_calc_interactionHOT, lj_displayAB, &
      lj_getcutoff, lj_display, lj_grpint, lj_snglint, lj_grpint2, &
      lj_chkequal, lj_snglint2, lj_getABCD, LJ_getcut

  Interface Assignment(=)
    Module Procedure lj_copy
  End Interface

  Interface lj_calc_interaction
    Module procedure lj_calc_snglInteraction
    Module procedure lj_calc_multInteraction
  End Interface

  Interface lj_calc_interactionHOT
    Module procedure lj_calc_multInteractionHOT
  End Interface

  Type LJ_Params
    Character(len=strLen) :: atom_name1, atom_name2
    Character(len=strLen) :: source
    Logical             :: off
    Logical             :: initialized
    Real(kind=RDbl)     :: sig, eps
    Real(kind=RDbl)     :: A, B, C, D
    Character(len=100)  :: line
    Character(len=80)   :: reference   !** not being used, but good idea
    Real(kind=RDbl)     :: hicut,hicut2
    Real(kind=RDbl)     :: ncut,ncut2    ! neighbor cutoff
    Real(kind=RDbl)     :: locut, locut2 ! lower bound to prevent overflow

    !** stuff below has not yet been implemented!
    Real(kind=RDbl)     :: smoothrad,smoothrad2
    Real(kind=RDbl)     :: cor             !truncation correction
    Real(kind=RDbl)     :: ljer,eler       !errors at the truncation radius
    Real(kind=RDbl)     :: sradlj,sradel   !pots at smoothing radius
    Real(kind=RDbl)     :: sraddlj,sraddel !pot derivatives at smoothing radius
  End Type LJ_Params

  Character(len=strLen), Parameter         :: lj_idstring = 'LJ'

  Real(kind=RDbl), Parameter               :: LJ_NCUTOFF = 200.00_RDbl


Contains
  !----------------------------------------------------------------------------
  ! Initialize the interaction information from the file line
  ! Requires:  params -- LJ parameter structure to initialize
  !            line -- input command line to interpret
  !----------------------------------------------------------------------------
  Subroutine lj_init(params,line)
    Type(LJ_Params), Intent(InOut)     :: params
    Character(*), Intent(In)           :: line

    Integer                                   :: nfields,n,i,unit,ios
    Integer                                   :: nlines
    Character(len=150), Dimension(10)         :: lines
    Character(len=strLen)                     :: filename,name1,name2
    Character(len=strLen), Dimension(10)      :: fields,chunks
    Logical                                   :: calcAB

    !** allow for the possibility that the line only contains the filename
    !** where the actual potentials are stored
    If (Index(Toupper(line),'FILE') /= 0) Then
      nfields = split(line,fields)
      name1 = Trim(fields(1))
      name2 = Trim(fields(2))
      filename = ''

      Do i = 1,nfields
        n = split(fields(i),chunks,'@')
        If (n > 1) Then
          If (Trim(chunks(1)) == 'FILE') Then
            filename = trim(chunks(2))
          End If
        End If
      End Do

      !** Open the file and get the line
      If (filename /= '') Then    
        unit = file_getunit(filename)    
        Open(unit=unit, file=filename, status='old', IOSTAT=ios)
        If (ios /= 0) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
               '  Could not open file ',trim(filename)
          Stop
        End If 
        Call getlineswpair(unit,name1,name2,' ',nlines,lines)
        Close(unit=unit)

        If (nlines == 0) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
               '  Could not find matching line in file: ',Trim(filename)
          Stop        
        End If
        If (nlines > 1) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
               '  Too many matching lines in file: ',Trim(filename)
          Stop        
        End If
      Else
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
             '  Could not find filename in line: ',Trim(line)
        Stop        
      End If
    Else
      lines(1) = line
    End If

    nfields = split(lines(1),fields)
    params%line = lines(1)
    params%off = .FALSE.
    params%locut = LOW_DIST_CUTOFF
    params%locut2 = LOW_DIST_CUTOFF*LOW_DIST_CUTOFF
    params%hicut = HIGH_DIST_CUTOFF
    params%hicut2 = HIGH_DIST_CUTOFF*HIGH_DIST_CUTOFF
    params%ncut  = LJ_NCUTOFF
    params%ncut2 = LJ_NCUTOFF*LJ_NCUTOFF
    params%atom_name1 = Trim(fields(1))
    params%atom_name2 = Trim(fields(2))
    calcAB = .True.
    
    params%C = 0._RDbl
    params%D = 0._RDbl

    Do i = 4,nfields
      n = split(fields(i),chunks,'@')
      Select Case(toupper(chunks(1)))
        
      Case('OFF')
        params%off = .True.
        !** Set the remaining parameters to zero
        params%sig = 0.0_RDbl
        params%eps = 0.0_RDbl
        params%A = 0.0_RDbl
        params%B = 0.0_RDbl
        params%hicut = 0.0_RDbl
        params%hicut2 = 0.0_RDbl

      Case('SIG')
        If (Index(chunks(2),'MIX') == 0) Then
          params%sig = toreal(chunks(2))
        End If
        calcAB = .True.

      Case('EPS')
        If (Index(chunks(2),'MIX') == 0) Then
          params%eps = toreal(chunks(2))
        End If
        calcAB = .True.

      Case('TRUNC')
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            " WARNING: Use of the label TRUNC has been depreciated in favor", &
            " of HICUT"
        params%hicut = toreal(chunks(2))
        params%hicut2 = params%hicut**2

      Case('HICUT')
        params%hicut = toreal(chunks(2))
        params%hicut2 = params%hicut**2

      Case('LOWCUT') 
        params%locut = toreal(chunks(2))
        params%locut2 = params%locut**2
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            " WARNING: Use of the label LOWCUT has been depreciated in", &
            " favor of LOCUT"

      Case('LOCUT') 
        params%locut = toreal(chunks(2))
        params%locut2 = params%locut**2

      Case('SMOOTH')
        params%smoothrad = toreal(chunks(2))
        params%smoothrad2 = params%smoothrad**2
         Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              'LJ Smoothing not currently supported'
         Stop

      Case('A')
        !** This is the A value
        params%A = toreal(chunks(2))
        calcAB = .False.

      Case('B')
        !** This is the B value
        params%B = toreal(chunks(2))
        calcAB = .False.

      Case('C')
        !** This is the C value
        params%C = toreal(chunks(2))
        calcAB = .False.
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'LJ Calculation using C parameter not currently supported'
        Stop

      Case('D')
        !** This is the D value
        params%D = toreal(chunks(2))
        calcAB = .False.
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            'LJ Calculation using D parameter not currently supported'
        Stop

      Case Default
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            'Unable to identify LJ interaction string ',Trim(chunks(1))
        Stop
      End Select
    End Do

    !** Generate the AB or sig/eps parameters
    If (calcAB) Then
      Call lj_calc_AB(params)
      params%source = "SIGEPS"
    Else
      Call lj_calc_sigeps(params)
      params%source = "AB"
    End If

    params%initialized = .True.

   End Subroutine lj_init

  !-------------------------------------------------------------------
  ! Returns true if the fields of "params" have been initialized
  !-------------------------------------------------------------------
  Logical Function lj_isinit(params)
    Type(LJ_Params), Intent(in) :: params
  
    lj_isinit = params%initialized
  End Function lj_isinit

  !-------------------------------------------------------------------
  ! Initialize params by copying from an existing params
  ! Requires: newparams -- new params, to be initialized
  !           oldparams -- old params
  !-------------------------------------------------------------------
  Subroutine lj_copy(newparams,oldparams)
    Type(LJ_Params), Intent(Out)      :: newparams
    Type(LJ_Params), Intent(In)       :: oldparams

    newparams%atom_name1 = oldparams%atom_name1
    newparams%atom_name1 = oldparams%atom_name2
    newparams%source = oldparams%source
    newparams%off = oldparams%off
    newparams%initialized = oldparams%initialized
    newparams%sig = oldparams%sig
    newparams%eps = oldparams%eps
    newparams%A = oldparams%A
    newparams%B = oldparams%B
    newparams%C = oldparams%C
    newparams%D = oldparams%D
    newparams%line = oldparams%line
    newparams%reference = oldparams%reference
    newparams%hicut = oldparams%hicut
    newparams%hicut2 = oldparams%hicut2
    newparams%ncut = oldparams%ncut
    newparams%ncut2 = oldparams%ncut2
    newparams%locut = oldparams%locut
    newparams%locut2 = oldparams%locut2

  End Subroutine lj_copy

  !--------------------------------------------------------------------------
  ! Calculate the interaction 
  !--------------------------------------------------------------------------
  Subroutine lj_calc_snglInteraction(params,sepvec,pot,ljflag,fvec)
    Type(LJ_Params), Intent(IN)               :: params
    Type(VecType), Intent(IN)                 :: sepvec
    Real(kind = RDbl), Intent(OUT)            :: pot
    Type(VecType), Intent(OUT), Optional      :: fvec
    Logical, Intent(Out)                      :: ljflag
  
    Real(kind = RDbl)                         :: A,B,r2i,r6i,r12i,r2
  
    !** Initialize the flag
    ljflag = .True.
  
    !** Initialize the potential and (optionally) the force
    pot = 0.0_RDbl
    If (present(fvec)) fvec = VecType(0.0_RDbl)
  
    !** Check if the interaction is on
    If (params%off) Return
  
    !** Get the square of the distance
    r2 = vector_getnormsq(sepvec)
  
    !**Check if it is within the cut-off radius
    If (r2 > params%hicut2) Return
    If (r2 < params%locut2) Then
      ljflag = .False.
      Return
    End If
  
    A = params%A
    B = params%B
  
    r2i = 1.0_RDbl/r2
    r6i = r2i*r2i*r2i
    r12i= r6i*r6i
  
    If (Present(fvec)) Then
      fvec = sepvec * (12.0_RDbl*A*r12i - 6.0_RDbl*B*r6i)*r2i
    End If
  
    pot = A*r12i - B*r6i
  
  End Subroutine lj_calc_snglInteraction
  
!!$  remove when you see this next- shaji. tuesday dec 31 2002
!!$  !---------------------------------------------------------------------------
!!$  ! Calculate the interaction 
!!$  !---------------------------------------------------------------------------
!!$  Subroutine lj_calc_snglInteractionHOT(params,sepvec,pot,ljflag,hot)
!!$    Type(LJ_Params), Intent(IN)               :: params
!!$    Type(VecType), Intent(IN)                 :: sepvec
!!$    Real(kind = RDbl), Intent(OUT)            :: pot
!!$    Logical, Intent(Out)                      :: ljflag
!!$    Real(kind = RDbl), Dimension(:), Intent(Out)  :: hot
!!$    Real(kind = RDbl)                         :: A,B,r2i,r6i,r12i,r2,r16i
!!$  
!!$    !** Initialize the flag
!!$    ljflag = .True.
!!$    !** Initialize the potential and (optionally) the force
!!$    pot = 0.0_RDbl
!!$    hot = 0.0_RDbl
!!$  
!!$    !** Check if the interaction is on
!!$    If (params%off) Return
!!$  
!!$    !** Get the square of the distance
!!$    r2 = vector_getnormsq(sepvec)
!!$  
!!$    !**Check if it is within the cut-off radius
!!$    If (r2 > params%hicut2) Return
!!$    If (r2 < params%locut2) Then
!!$      ljflag = .False.
!!$      Return
!!$    End If
!!$  
!!$    A = params%A
!!$    B = params%B
!!$  
!!$    r2i = 1.0_RDbl/r2
!!$    r6i = r2i*r2i*r2i
!!$    r12i= r6i*r6i
!!$    r16i= r6i*r2i*r6i*r2i
!!$  
!!$    !** Calculate the potential
!!$    pot = A*r12i - B*r6i
!!$  
!!$    !** Calculate the higher order terms
!!$
!!$    hot(1)=pot
!!$    hot(2) = sepvec%comp(1) * (-12.0_RDbl*A*r12i + 6.0_RDbl*B*r6i)*r2i
!!$    hot(3) = sepvec%comp(2) * (-12.0_RDbl*A*r12i + 6.0_RDbl*B*r6i)*r2i
!!$    hot(4) = sepvec%comp(3) * (-12.0_RDbl*A*r12i + 6.0_RDbl*B*r6i)*r2i
!!$    hot(5) = (168.0_RDbl*A*r12i - 48.0_RDbl*B*r6i)*r2i*r2i &
!!$        *sepvec%comp(1)*sepvec%comp(2)
!!$    hot(6) = (168.0_RDbl*A*r12i - 48.0_RDbl*B*r6i)*r2i*r2i &
!!$        *sepvec%comp(1)*sepvec%comp(3)
!!$    hot(7) = (168.0_RDbl*A*r12i - 48.0_RDbl*B*r6i)*r2i*r2i &
!!$        *sepvec%comp(2)*sepvec%comp(3)
!!$    hot(8) = (-2688.0_RDbl*A*r16i + 480.0_RDbl*B*r12i)*r2i &
!!$        *sepvec%comp(1)*sepvec%comp(2)*sepvec%comp(3)
!!$
!!$    
!!$!!L    hot(1) = sepvec * (-12.0_RDbl*A*r12i + 6.0_RDbl*B*r6i)*r2i
!!$!!L    hot(2)%comp(1) = (168.0_RDbl*A*r12i - 48.0_RDbl*B*r6i)*r2i*r2i &
!!$!!L        *sepvec%comp(1)*sepvec%comp(2)
!!$!!L    hot(2)%comp(2) = (168.0_RDbl*A*r12i - 48.0_RDbl*B*r6i)*r2i*r2i &
!!$!!L        *sepvec%comp(1)*sepvec%comp(3)
!!$!!L    hot(2)%comp(3) = (168.0_RDbl*A*r12i - 48.0_RDbl*B*r6i)*r2i*r2i &
!!$!!L        *sepvec%comp(2)*sepvec%comp(3)
!!$!!L    hot(3)%comp(1) = (-2688.0_RDbl*A*r16i + 480.0_RDbl*B*r12i)*r2i &
!!$!!L        *sepvec%comp(1)*sepvec%comp(2)*sepvec%comp(3)
!!$!!L    hot(3)%comp(2) = hot(3)%comp(1)
!!$!!L    hot(3)%comp(3) = hot(3)%comp(1)
!!$  
!!$  End Subroutine lj_calc_snglInteractionHOT
  
  !---------------------------------------------------------------------------
  ! Calculate the interaction between a molecule and an entire sorbate 
  !---------------------------------------------------------------------------
  Subroutine lj_calc_multInteraction(params, apairs,sepvec,pot, ljflag, fvec)
    Type(LJ_Params), Dimension(:,:), Intent(In) :: params
    Integer, Dimension(:,:), Intent(In)       :: apairs
    Type(VecType), Dimension(:), Intent(IN)   :: sepvec
    Real(kind = RDbl), Intent(OUT) :: pot
    Logical, Intent(out)                      ::  ljflag   
    Type(VecType), Dimension(:), Intent(OUT), Optional :: fvec
  
    Real(kind = RDbl)                         :: A,B,r2i,r6i,r12i,r2
    Integer :: i, a1, a2
  
    !** Zero the forces
    pot = 0.0_RDbl
    If (present(fvec)) fvec = VecType(0.0_RDbl)
  
    !** Loop through all the sepvecs
    Do i = 1, Size(sepvec)
      !** Parameter indices
      a1 = apairs(1,i)
      a2 = apairs(2,i)
  
      !** Check to see if the interaction is on
      If (params(a1,a2)%off) Cycle
  
      !** Get the square of the distance
      r2 = vector_getnormsq(sepvec(i))
  
      !**Check if it is within the cut-off radius
      If (r2 > params(a1,a2)%hicut2) Cycle
      If (r2 < params(a1,a2)%locut2) Then 
        ljflag= .False.
        Return
      End If
  
      A = params(a1,a2)%A
      B = params(a1,a2)%B
  
      r2i = 1.0_RDbl/r2
      r6i = r2i*r2i*r2i
      r12i= r6i*r6i
  
      If (Present(fvec)) Then
        fvec(i) = sepvec(i) * (12.0_RDbl*A*r12i - 6.0_RDbl*B*r6i)*r2i
      End If
  
      pot = pot + A*r12i - B*r6i
      ljflag= .True. 
    End Do
  End Subroutine lj_calc_multInteraction

  !---------------------------------------------------------------------------
  ! Calculate the interaction between a set of atomic points.  Returns False
  ! if any of the seperation vectors fall within the lower cutoff.  This 
  ! indicates atoms have gotten too close together and that the potential 
  ! would be extremely large.  The returned derivative components are those
  ! for the first (?) atom in the separation vector defined by (atom1 - atom2)
  ! Requires:  params -- Lennard-Jones interaction parameters
  !            sepvec -- array of separation vectors
  !            ffout -- array of resulting energies and derivatives
  !---------------------------------------------------------------------------
  Logical Function lj_grpint(params,sepvec,results)
    Type(LJ_Params), Dimension(:), Intent(In)     :: params
    Type(VecType), Dimension(:), Intent(In)       :: sepvec
    Type(EnergyPlus), Dimension(:), Intent(InOut) :: results

    Integer                   :: pair
    Real(kind = RDbl)         :: A,B,r2i,r6i,r12i,r2,r16i,nrg
    Type(VecType)             :: grad

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) Size(sepvec),Size(results)

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Do pair = 1,Size(sepvec)
      Write(*,*) pair,Trim(storebase_disp(results(pair)))
    End Do
    Write(*,*) 
#endif

    Do pair = 1,Size(sepvec)

      lj_grpint = .False.

      !** Get the square of the distance
      r2 = vector_getnormsq(sepvec(pair))
!      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!      Write(*,*) Sqrt(r2),r2,params(pair)%hicut2
  
      !** Check the cutoffs
      If (r2 > params(pair)%hicut2) Cycle
      If (r2 < params(pair)%locut2) Return
  
      A = params(pair)%A
      B = params(pair)%B
  
      r2i = 1.0_RDbl/r2
      r6i = r2i*r2i*r2i
      r12i= r6i*r6i

      !** Calculate the potential energy
      nrg = A*r12i - B*r6i
      Call storebase_inc(results(pair),nrg)

      !** Calculate the gradients of the energy  
      If (storebase_nderivs(results(pair)) > 0) Then
        grad = sepvec(pair) * (12.0_RDbl*A*r12i - 6.0_RDbl*B*r6i)*r2i
        Call storebase_inc(results(pair),grad)
      End If

      !** Calculate select second and third derivatives (HOT)
      If (storebase_nderivs(results(pair)) > 1) Then
        r16i = r6i*r2i*r6i*r2i
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        stop
        !PLEASE FILL THIS IN
      End If
  
    End Do

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Do pair = 1,Size(sepvec)
      Write(*,'(i3,a,f10.4)') pair,Trim(storebase_disp(results(pair))), &
          mag(sepvec(pair))
    End Do
    Write(*,*) 
#endif

    lj_grpint = .True. 

  End Function lj_grpint

  !---------------------------------------------------------------------------
  ! attempt to speed by passing in nderivs
  !---------------------------------------------------------------------------
  Logical Function lj_grpint2(params,sepvec,results,nderivs)
    Type(LJ_Params), Dimension(:), Intent(In)     :: params
    Type(VecType), Dimension(:), Intent(In)       :: sepvec
    Type(EnergyPlus), Dimension(:), Intent(InOut) :: results
    Integer, Intent(In)                           :: nderivs

    Integer                   :: pair
    Real(kind = RDbl)         :: A,B,r2i,r6i,r12i,r2,r16i,nrg
    Type(VecType)             :: grad

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) Size(sepvec),Size(results)

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Do pair = 1,Size(sepvec)
      Write(*,*) pair,Trim(storebase_disp(results(pair)))
    End Do
    Write(*,*) 
#endif

    Do pair = 1,Size(sepvec)

      lj_grpint2 = .False.

      !** Get the square of the distance
      r2 = vector_getnormsq(sepvec(pair))
!      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!      Write(*,*) Sqrt(r2),r2,params(pair)%hicut2
  
      !** Check the cutoffs
      If (r2 > params(pair)%hicut2) Cycle
      If (r2 < params(pair)%locut2) Return
  
      A = params(pair)%A
      B = params(pair)%B
  
      r2i = 1.0_RDbl/r2
      r6i = r2i*r2i*r2i
      r12i= r6i*r6i

      !** Calculate the potential energy
      nrg = A*r12i - B*r6i
      Call storebase_inc(results(pair),nrg)

      !** Calculate the gradients of the energy  
      If (nderivs > 0) Then
        grad = sepvec(pair) * (12.0_RDbl*A*r12i - 6.0_RDbl*B*r6i)*r2i
        Call storebase_inc(results(pair),grad)
!        results(pair)%grad(1) = grad%comp(1)
!        results(pair)%grad(2) = grad%comp(2)
!        results(pair)%grad(3) = grad%comp(3)
      End If

      !** Calculate select second and third derivatives (HOT)
      If (nderivs > 1) Then
        r16i = r6i*r2i*r6i*r2i
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        stop
        !PLEASE FILL THIS IN
      End If
  
    End Do

#ifdef DEBUG
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Do pair = 1,Size(sepvec)
      Write(*,'(i3,a,f10.4)') pair,Trim(storebase_disp(results(pair))), &
          mag(sepvec(pair))
    End Do
    Write(*,*) 
#endif

    lj_grpint2 = .True. 

  End Function lj_grpint2

  !---------------------------------------------------------------------------
  ! Calculate the interaction between a set of atomic points.  ReturnsFalse
  ! if any of the seperation vectors fall within the lower cutoff.  This 
  ! indicates atoms have gotten too close together and that thepotential 
  ! would be extremely large.  The returned derivative components arethose
  ! for the first atom in the separation vector defined by (atom1 -atom2).
  ! The optional second results structure must be passed if forces areto be
  ! calculated.
  ! Requires:  params -- Lennard-Jones interaction parameters
  !            sepvec -- atom-atom separation vector
  !            results -- resultant energies and derivatives
  !            results2 -- resultant energies and opposite forces
  !---------------------------------------------------------------------------
  Logical Function lj_snglint(params,sepvec,results,results2,hot)
    Type(LJ_Params), Intent(In)               :: params
    Type(VecType), Intent(In)                 :: sepvec
    Type(EnergyPlus), Intent(InOut)           :: results
    Type(EnergyPlus), Intent(InOut), Optional :: results2
    Real(kind = RDbl), Dimension(SIZE_OF_HOT), Optional :: hot 

    Real(kind = RDbl)         :: A,B,r2i,r6i,r12i,r2,r16i,r10i,nrg,fac1
    Type(VecType)             :: force
    Real(kind = RDbl), Dimension(1:4) :: hess

    lj_snglint = .False.

    !** Get the square of the distance
    r2 = vector_getnormsq(sepvec)

    !** Check the cutoffs
    If (r2 > params%hicut2) Then
      lj_snglint = .True.      
      Return
    End If
    If (r2 < params%locut2) Return

    A = params%A
    B = params%B

    r2i = 1.0_RDbl/r2
    r6i = r2i*r2i*r2i
    r12i= r6i*r6i

    !** Calculate the potential energy
    nrg = A*r12i - B*r6i
    results%nrg = results%nrg + nrg    

    !** Calculate the gradients of the energy  
    If (storebase_nderivs(results) > 0) Then

      fac1=(12.0_RDbl*A*r12i - 6.0_RDbl*B*r6i)*r2i

      force%comp(1) = sepvec%comp(1) * fac1
      force%comp(2) = sepvec%comp(2) * fac1
      force%comp(3) = sepvec%comp(3) * fac1

      !** I dont like this here, not always necessary  SHAJI
      Call storebase_inc(results,force)
      Call storebase_inc(results2,(force*(-1.0_RDbl)))
    End If

    If (Present(hot)) Then    
      r10i = r12i*r2  !  note : r2= 1/r2i
      r16i = r10i*r6i

      !** Calculate select second and third derivatives (HOT)
      !** Calculate the higher order terms

      hot(1) = nrg

      hot(2) = -(force%comp(1))
      hot(3) = -(force%comp(2))
      hot(4) = -(force%comp(3))

      fac1=(168.0_RDbl*A*r16i - 48.0_RDbl*B*r10i)
      hot(5) = fac1 *sepvec%comp(1)*sepvec%comp(2)
      hot(6) = fac1 *sepvec%comp(1)*sepvec%comp(3)
      hot(7) = fac1 *sepvec%comp(2)*sepvec%comp(3)

      fac1=(-2688.0_RDbl*A*r16i + 480.0_RDbl*B*r10i)*r2i 
      hot(8) = fac1*sepvec%comp(1)*sepvec%comp(2)*sepvec%comp(3)

      hess(1) = hot(5)
      hess(2) = hot(6)
      hess(3) = hot(7)
      hess(4) = hot(8)

      !** Updating hessian in storage
      Call storebase_inc(results,hess)

    End If

    lj_snglint = .True. 

  End Function lj_snglint

  !---------------------------------------------------------------------------
  ! Calculate the interaction between tow atomic points.  Returns False
  ! if any of the seperation vectors fall within the lower cutoff. 
  ! FAST, useful for CBGCMC, no forces here
  ! Requires:  params -- Lennard-Jones interaction parameters
  !            sepvec -- atom-atom separation vector
  !            results -- resultant energies and derivatives
  !---------------------------------------------------------------------------
  Logical Function lj_snglint2(params,sepvec,results)
    Type(LJ_Params), Intent(In)               :: params
    Type(VecType), Intent(In)                 :: sepvec
    Type(EnergyPlus), Intent(InOut)           :: results

    Real(kind = RDbl)         :: A,B,r2i,r6i,r12i,r2,nrg

    lj_snglint2 = .False.

    !** Get the square of the distance
    r2 = vector_getnormsq(sepvec)

    !** Check the cutoffs
    If (r2 > params%hicut2) Then
      lj_snglint2 = .True.      
      Return
    End If
    If (r2 < params%locut2) Return

    A = params%A
    B = params%B

    r2i = 1.0_RDbl/r2
    r6i = r2i*r2i*r2i
    r12i= r6i*r6i

    !** Calculate the potential energy
    nrg = A*r12i - B*r6i
    results%nrg = results%nrg + nrg

    lj_snglint2 = .True. 

  End Function lj_snglint2

  
  
  !---------------------------------------------------------------------------
  ! Calculate the interaction between a molecule and an entire sorbate 
  !---------------------------------------------------------------------------
  Subroutine lj_calc_multInteractionHOT(params, apairs, sepvec, pot,ljflag,hot)
    Type(LJ_Params), Dimension(:,:), Intent(In) :: params
    Integer, Dimension(:,:), Intent(In)         :: apairs
    Type(VecType), Dimension(:), Intent(IN)     :: sepvec
    Real(kind = RDbl), Intent(OUT)              :: pot
    Logical, Intent(out)                        ::  ljflag   
    Real(kind = RDbl), Dimension(:), Intent(OUT)    :: hot
  
    Real(kind = RDbl) :: A,B,r2i,r6i,r12i,r2,r16i
    Integer :: i, a1, a2
  
    !** Zero the forces
    pot = 0.0_RDbl
  
    !** Zero the higher order terms
    hot = 0.0_RDbl
  
  
    !** Loop through all the sepvecs
    Do i = 1, Size(sepvec)
  
      !** Parameter indices
      a1 = apairs(1,i)
      a2 = apairs(2,i)
  
      !** Check to see if the interaction is on
      If (params(a1,a2)%off) Cycle
  
      !** Get the square of the distance
      r2 = vector_getnormsq(sepvec(i))
  
      !**Check if it is within the cut-off radius
      If (r2 > params(a1,a2)%hicut2) Then
  
        Cycle
      End If
      If (r2 < params(a1,a2)%locut2) Then 
        ljflag= .False.
        Return
      End If
  
      A = params(a1,a2)%A
      B = params(a1,a2)%B
  
      r2i = 1.0_RDbl/r2
      r6i = r2i*r2i*r2i
      r12i= r6i*r6i
      r16i= r6i*r2i*r6i*r2i
  
      !** Calculate the potential
      pot = pot + A*r12i - B*r6i


      !** Calculate the higher order terms
      hot(1)=pot
      hot(2) = hot(2) + sepvec(i)%comp(1) * (-12.0_RDbl*A*r12i+6.0_RDbl*B*r6i)*r2i
      hot(3) = hot(3) + sepvec(i)%comp(2) * (-12.0_RDbl*A*r12i+6.0_RDbl*B*r6i)*r2i
      hot(4) = hot(4) + sepvec(i)%comp(3) * (-12.0_RDbl*A*r12i+6.0_RDbl*B*r6i)*r2i
      hot(5) = hot(5)+ (168.0_RDbl*A*r12i-48.0_RDbl*B*r6i) &
          *r2i*r2i*sepvec(i)%comp(1)*sepvec(i)%comp(2)
      hot(6) = hot(6) + (168.0_RDbl*A*r12i-48.0_RDbl*B*r6i) &
          *r2i*r2i*sepvec(i)%comp(1)*sepvec(i)%comp(3)
      hot(7) = hot(7) + (168.0_RDbl*A*r12i-48.0_RDbl*B*r6i) &
          *r2i*r2i*sepvec(i)%comp(2)*sepvec(i)%comp(3)
      hot(8)= hot(8) + (-2688.0_RDbl*A*r16i + &
          480.0_RDbl*B*r12i)*r2i * sepvec(i)%comp(1)*sepvec(i)%comp(2)* &
          sepvec(i)%comp(3)
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "wrong hot(8)"
      Stop
      ljflag= .True. 
    End Do
  End Subroutine lj_calc_multInteractionHOT

  !---------------------------------------------------------------------------
  ! This routine calculates A, B parameters from sigma
  ! and epsilon.  A = 4*eps*(sig**12), B = 4*eps*(sig**6)
  !---------------------------------------------------------------------------
  Subroutine lj_calc_AB(obj)
    Implicit None
    Type(LJ_Params), Intent(INOUT)         :: obj
  
    obj%B = 4.0_Rdbl*obj%eps*(obj%sig**6)*kcalmole_kb
    obj%A = obj%b*obj%sig**6

 !   Write(*,*) 'LJ: sig, eps ',lj_displayparams(obj)
 !   Write(*,*) 'LJ: A12, B6 ',obj%A,obj%B
  End Subroutine lj_calc_AB


  !---------------------------------------------------------------------------
  ! This routine calculates the sigma and epsilon parameters from A and B
  ! Requires:  obj -- LJ parameter structure
  !---------------------------------------------------------------------------
  Subroutine lj_calc_sigeps(obj)
    Type(LJ_Params), Intent(InOut)   :: obj
  
    Real(kind=RDbl)           :: sig6

    !** Set Defaults
    obj%eps = 0.0_RDbl
    obj%sig = 0.0_RDbl

    !** Leave if one of the parameters is very small
    If ((Abs(obj%A) < 1.0e-10_RDbl) .Or. (Abs(obj%B) < 1.0e-10_RDbl)) Return

    !** Convert to sigma and epsilon
    sig6    = (obj%A/obj%B)
    obj%eps = (obj%B/(4.0_Rdbl*sig6*kcalmole_kb))
    obj%sig = sig6**(1/6.0)

  End Subroutine lj_calc_sigeps
  
  !---------------------------------------------------------------------------
  ! This routine mixes sigma and epsilon values according to the either
  ! the Lorentz-Bertholet mixing rules or the Jorgensen mixing rules
  ! Requires: params -- LJ params
  !---------------------------------------------------------------------------
  Subroutine lj_calc_mix(obj,pure1,pure2,line)
    Implicit None
    Type(LJ_Params), Intent(INOUT)   :: obj
    Type(LJ_Params), Intent(IN)      :: pure1,pure2
    Character(*), Intent(IN)         :: line
  
    Integer                               :: nfields,n,i
    Logical                               :: got_sig,got_eps
    Character(len=strLen), Dimension(10)  :: fields,chunks

    !** Assign the basic init fields first
    nfields = split(line,fields)
    obj%line = line
    obj%off  = .False.
    obj%locut = LOW_DIST_CUTOFF
    obj%hicut = HIGH_DIST_CUTOFF
    obj%atom_name1 = Trim(fields(1))
    obj%atom_name2 = Trim(fields(2))

    !** Do the mixing
    got_sig = .False.
    If (Index(Toupper(line),'SIG@LBMIX') /= 0) Then
      obj%sig = (pure1%sig + pure2%sig)/2.0_RDbl   !arithmetic
      got_sig = .True.
    Else If (Index(Toupper(line),'SIG@JORGMIX') /= 0) Then
      obj%sig = Sqrt(pure1%sig * pure2%sig)        !geometric
      got_sig = .True.
    End If

    got_eps = .False.  
    If (Index(Toupper(line),'EPS@LBMIX') /= 0) Then
      obj%eps = Sqrt((pure1%eps * pure2%eps))      !geometric
      got_eps = .True.
    Else If (Index(Toupper(line),'EPS@JORGMIX') /= 0) Then
      obj%eps = Sqrt((pure1%eps * pure2%eps))      !geometric
      got_eps = .True.
    End If
    
    If ((.Not. got_sig) .Or. (.Not. got_eps)) Then
      Write(0,'(2a,i4,a,2l2)') __FILE__,":",__LINE__, &
          " Could not find mixing rules specification for SIG and/or EPS ", &
          got_sig, got_eps
      Write(0,'(4x,a)') Trim(line)
      Stop
    End If

    !** Get the rest of the information from the line
    ! Skip the first 3, those are atom names and LJ keyword
    Do i = 4,nfields
      n = split(fields(i),chunks,'@')
      Select Case(toupper(chunks(1)))
      Case('HICUT')
        obj%hicut = toreal(chunks(2))
        obj%hicut2 = obj%hicut**2
      Case('TRUNC')
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            " WARNING: Use of the label TRUNC has been depreciated in favor", &
            " of HICUT"
        obj%hicut = toreal(chunks(2))
        obj%hicut2 = obj%hicut**2
      Case('SMOOTH')
        obj%smoothrad = toreal(chunks(2))
        obj%smoothrad2 = obj%smoothrad**2
      Case('LOCUT')
        obj%locut = toreal(chunks(2))
        obj%locut2 = obj%locut**2 
      Case('LOWCUT') 
        obj%locut = toreal(chunks(2))
        obj%locut2 = obj%locut**2
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            " WARNING: Use of the label LOWCUT has been depreciated in", &
            " favor of LOCUT"
      End Select
    End Do

    !** Generate the AB parameters
    Call lj_calc_AB(obj)

    obj%initialized = .True.
    obj%source = "MIXED"

  End Subroutine lj_calc_mix

  !----------------------------------------------------------------------------
  ! Returns the value of the high end cutoff. If the NEIGHBOR flag is 
  ! present, check to see if the NEIGHBOR cutoff should be returned
  ! Neighbour cut-off is larger value than  high cutoff. It is useful while
  ! making list of intra-pairs in the beginning during molecule 
  ! initialization.
  ! Requires: params -- LJ params
  !           neighbor -- flag to get neighbor cutoff
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function lj_getcutoff(params,neighbor)
    Type(LJ_Params), Intent(In)         :: params
    Logical, Optional, Intent(In)       :: neighbor

    !** Default value
    lj_getcutoff = params%hicut

    !** Check if we need to return the neighbor cutoff
    If (Present(neighbor)) Then
      If (neighbor) Then
        lj_getcutoff = params%ncut
      End If
    End If
  End Function lj_getcutoff

  !----------------------------------------------------------------------------
  ! Checks to see if two LJ type data structures contain the same parameters
  ! Requires:  params1 -- 1st LJ params
  !            params2 -- 2nd LJ params
  !            tolerance -- optional tolerance for checking
  !----------------------------------------------------------------------------
  Logical Function lj_chkequal(params1,params2,tolerance)
    Type(LJ_Params), Intent(In)            :: params1,params2
    Real(kind=RDbl), Intent(In), Optional  :: tolerance

    Real(kind=RDbl)  :: tol,diff

    tol = 1.0e-6_RDbl
    If (Present(tolerance)) tol = tolerance

    lj_chkequal = .True.

    diff = Abs(params1%A - params2%A)
    If (diff > tol) lj_chkequal = .False.

    diff = Abs(params1%B - params2%B)
    If (diff > tol) lj_chkequal = .False.

  End Function lj_chkequal
 
  !-------------------------------------------------------
  ! Returns a string with the properly formatted lj params
  ! Requires: params -- LJ params
  !-------------------------------------------------------
  Function lj_displayparams(ljparams)
    Character(len=strLen)       :: lj_displayparams
    Type(LJ_Params), Intent(In) :: ljparams
    Character(len=strLen)       :: streps, strsig
    strsig=real2str(ljparams%sig,5)
    streps=real2str(ljparams%eps,5)
    Write(lj_displayparams,'(1x,a,2x,a)') Trim(strsig), Trim(streps)
  End Function lj_displayparams

  !----------------------------------------------------------------------------
  ! Returns a string containing the sigma and epsilon values
  ! Requires: params -- LJ params
  !----------------------------------------------------------------------------
  Function lj_displaySigEps(ljparams)
    Character(len=lstrLen)      :: lj_displaySigEps
    Type(LJ_Params), Intent(In) :: ljparams

    Character(len=lstrLen)      :: string1,string2

    If (ljparams%off) Then
      Write(lj_displaySigEps,'(a)') "OFF"
    Else
      string1 = real2str(ljparams%sig,6)
      string2 = real2str(ljparams%eps,6)
      Write(lj_displaySigEps,'(2(2a,1x,a,3x))') "Sigma: ", &
          Trim(string1),"Ang","Epsilon/k_Boltz: ",Trim(string2),"K"
    End If

  End Function lj_displaySigEps

  !----------------------------------------------------------------------------
  ! Returns a string containing nicely formated A and B values
  ! Requires: params -- LJ params
  !----------------------------------------------------------------------------
  Function lj_displayAB(ljparams)
    Character(len=lstrLen)              :: lj_displayAB
    Type(LJ_Params), Intent(In)         :: ljparams

    Character(len=xlstrLen)             :: string1,string2

    If (ljparams%off) Then
      Write(lj_displayAB,'(a)') "OFF"
    Else
      string1 = real2str(ljparams%A,8) 
      string2 = real2str(ljparams%B,8)
      Write(lj_displayAB,'(2(2a,1x,a,3x))') "A: ", &
          Trim(string1),"kcal A^12/mol","B: ",Trim(string2),"kcal A^6/mol"
    End If

  End Function lj_displayAB

  !----------------------------------------------------------------------------
  ! Returns a string containing the cutoff Values
  ! Requires: params -- LJ params
  !----------------------------------------------------------------------------
  Function lj_displayCutOffs(ljparams)
    Character(len=lstrLen)          :: lj_displayCutOffs
    Type(LJ_Params), Intent(In)     :: ljparams
    
    Character(len=strLen)           :: string1,string2

    If (ljparams%off) Then
      Write(lj_displayCutoffs,'(a)') "OFF"
    Else
      string1 = real2str(ljparams%hicut,5)
      string2 = real2str(ljparams%locut,5)
      Write(lj_displayCutOffs,'(2(2a,1x,a,3x),f8.4)') &
          "high cutoff: ",Trim(string1),"Ang", &
          "low cutoff: ",Trim(string2),"Ang"
    End If

  End Function lj_displayCutOffs

  !-----------------------------------------------------------------------
  ! Displays the full set of LJ parameters for the structure
  ! Requires: params -- the LJ parameters structure
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !-----------------------------------------------------------------------
  Subroutine lj_display(params,indent,unit)
    Type(LJ_Params), Intent(In)         :: params
    Integer, Intent(In)                 :: indent,unit

    Character(len=indent)       :: blank
    Character(len=xlstrLen)     :: string

    blank = Repeat(' ',indent)

    If (.Not. params%initialized) Return

    If (params%off) Then
      Write(unit,'(5a)') blank,Trim(params%atom_name1),'-', &
          Trim(params%atom_name2),' Interaction OFF'
    Else
      Write(unit,'(6a)') blank,Trim(params%atom_name1),'-',&
          Trim(params%atom_name2),' LJ parameters from ',Trim(params%source)
      string = lj_displaySigEps(params)
      Write(unit,'(2x,2a)') blank,Trim(string)
      string = lj_displayAB(params)
      Write(unit,'(2x,2a)') blank,Trim(string)
      string = lj_displayCutOffs(params)
      Write(unit,'(2x,2a)') blank,Trim(string)
    End If

  End Subroutine lj_display

!-----------------------------------------------------------------------------
! Tina added
! Gets the A, B, C, D parameters of the Lennard-Jones potential
!-----------------------------------------------------------------------------

  Subroutine lj_getABCD(params_lj,pot_param)

   Type(LJ_params), Intent(IN)      :: params_lj
   Real(kind = RDbl), Dimension(4)  :: pot_param

   pot_param(1) = params_lj%A
   pot_param(2) = params_lj%B
   pot_param(3) = params_lj%C
   pot_param(4) = params_lj%D
 
  END Subroutine lj_getABCD 
!-----------------------------------------------------------------------------
! Tina added
! Gets the hicut and locut parameters of the Lennard-Jones potential
!-----------------------------------------------------------------------------

  Subroutine lj_getcut(params_lj,pot_param)

   Type(LJ_params), Intent(IN)      :: params_lj
   Real(kind = RDbl), Dimension(2)  :: pot_param

   pot_param(1) = params_lj%hicut
   pot_param(2) = params_lj%locut
 
  End Subroutine lj_getcut

End Module lj



