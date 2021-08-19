!------------------------------------------------------------------------------
! This module contains the routines specific to initializing and actually 
! calculating LENNARD-JONES INTERACTIONS.  Interactions can be calculated
! either between a single pair (lj206_calc_snglInteraction) or for a set 
! of pairs (lj206_calc_multInteraction).  Higher order terms 
! (HOT), meaning first and second derivatives, are also available.  The
! data structure (LJ206_Params) contains all the potential parameterization
! information.  It is initialized using file information in "lj206_init"
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

Module lj206

  Use defaults, Only: RDbl, strLen, lstrLen, LOW_DIST_CUTOFF, kcalmole_kb
  Use utils, Only: toupper, split, toreal, getlineswpair, real2str, int2str
  Use file, Only: file_getunit
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnormsq, vector_display
  Use atom, Only: atom_getmass
  Use storebase, Only: EnergyPlus, storebase_inc, storebase_nderivs, &
      storebase_disp

  Implicit None
  Save

  Private
  Public :: LJ206_Params, lj206_displayparams, lj206_snglint, &
      lj206_idstring, lj206_isinit, lj206_displaySigEps, lj206_calc_mix, &
      lj206_calc_ab, Assignment(=), lj206_init, lj206_displayCutOffs, &
      lj206_displayAB, lj206_getcutoff, lj206_display, lj206_chkequal

  Interface Assignment(=)
    Module Procedure lj206_copy
  End Interface

  Type LJ206_Params
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
    Real(kind=RDbl)     :: lj206er,eler       !errors at the truncation radius
    Real(kind=RDbl)     :: sradlj206,sradel   !pots at smoothing radius
    Real(kind=RDbl)     :: sraddlj206,sraddel !pot derivatives at smoothing radius
  End Type LJ206_Params

  Character(len=strLen), Parameter         :: lj206_idstring = 'LJ206'

  Real(kind=RDbl), Parameter               :: LJ206_NCUTOFF = 200.00_RDbl


Contains
  !----------------------------------------------------------------------------
  ! Initialize the interaction information from the file line
  !----------------------------------------------------------------------------
  Subroutine lj206_init(params,line)
    Type(LJ206_Params), Intent(INOUT)         :: params
    Character(*), Intent(IN)                  :: line

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
    params%ncut  = LJ206_NCUTOFF
    params%ncut2 = LJ206_NCUTOFF*LJ206_NCUTOFF
    params%atom_name1 = Trim(fields(1))
    params%atom_name2 = Trim(fields(2))
    calcAB = .True.

    Do i = 4,nfields
      n = split(fields(i),chunks,'@')
      Select Case(toupper(chunks(1)))
        
      Case('OFF')
        params%off = .TRUE.
        !** Set the remaining parameters to zero
        params%sig = 0.0_RDbl
        params%eps = 0.0_RDbl
        params%A = 0.0_RDbl
        params%B = 0.0_RDbl
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
              'LJ206 Calculation using C parameter not currently supported'
         Stop
       Case('D')
         !** This is the D value
         params%D = toreal(chunks(2))
         calcAB = .False.
         Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
              'LJ206 Calculation using D parameter not currently supported'
         Stop
       Case Default
         Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              'Unable to identify LJ206 interaction string ',trim(chunks(1))
         Stop
       End Select

     Enddo

     !** Generate the AB or sig/eps parameters
     If (calcAB) Then
       Call lj206_calc_AB(params)
       params%source = "SIGEPS"
     Else
       Call lj206_calc_sigeps(params)
       params%source = "AB"
     End If

     params%initialized = .True.

   End Subroutine lj206_init

  !-------------------------------------------------------------------
  ! Returns true if the fields of "params" have been initialized
  !-------------------------------------------------------------------
  Logical Function lj206_isinit(params)
    Type(LJ206_Params), Intent(in) :: params
  
    lj206_isinit = params%initialized
  End Function lj206_isinit

  !-------------------------------------------------------------------
  ! Initialize params by copying from an existing params
  ! Requires: newparams -- new params, to be initialized
  !           oldparams -- old params
  !-------------------------------------------------------------------
  Subroutine lj206_copy(newparams,oldparams)
    Type(LJ206_Params), Intent(Out)      :: newparams
    Type(LJ206_Params), Intent(In)       :: oldparams

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

  End Subroutine lj206_copy

  !---------------------------------------------------------------------------
  ! Calculate the interaction between a set of atomic points.  Returns False
  ! if any of the seperation vectors fall within the lower cutoff.  This 
  ! indicates atoms have gotten too close together and that the potential 
  ! would be extremely large.  The returned derivative components are those
  ! for the first atom in the separation vector defined by (atom1 - atom2).
  ! The optional second results structure must be passed if forces are to be
  ! calculated.
  ! Requires:  params -- Lennard-Jones 20-6 interaction parameters
  !            sepvec -- atom-atom separation vector
  !            results -- resultant energies and derivatives
  !            results2 -- resultant energies and opposite forces
  !---------------------------------------------------------------------------
  Logical Function lj206_snglint(params,sepvec,results,results2)
    Type(LJ206_Params), Intent(In)            :: params
    Type(VecType), Intent(In)                 :: sepvec
    Type(EnergyPlus), Intent(InOut)           :: results
    Type(EnergyPlus), Intent(InOut), Optional :: results2

    Real(kind = RDbl)         :: A,B,r2i,r6i,r20i,r2,nrg
    Type(VecType)             :: force

    lj206_snglint = .False.
    
    !** Get the square of the distance
    r2 = vector_getnormsq(sepvec)
    
    !** Check the cutoffs
    If (r2 > params%hicut2) Then
      lj206_snglint = .True.      
      Return
    End If
    If (r2 < params%locut2) Return
    
    A = params%A
    B = params%B

    r2i = 1.0_RDbl/r2
    r6i = r2i*r2i*r2i
    r20i = r6i*r6i*r6i*r2i
    
    !** Calculate the potential energy
    nrg = A*r20i - B*r6i
    Call storebase_inc(results,nrg)
!    If (Present(results2)) Call storebase_inc(results2,nrg)
!    results%total%nrg = results%total%nrg + nrg
    
    !** Calculate the gradients of the energy  
    If (storebase_nderivs(results) > 0) Then
      force = sepvec * (20.0_RDbl*A*r20i - 6.0_RDbl*B*r6i)*r2i
      Call storebase_inc(results,force)
      Call storebase_inc(results2,(force*(-1.0_RDbl)))
!    results%total%grad = results%total%grad + force
    End If
    
    !** Calculate select second and third derivatives (HOT)
    If (storebase_nderivs(results) > 1) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Stop
      !PLEASE FILL THIS IN
    End If
    
    lj206_snglint = .True. 

  End Function lj206_snglint

  !--------------------------------------------------------------------------
  ! Calculate the interaction 
  !--------------------------------------------------------------------------
  Subroutine lj206_calc_snglInteraction(params,sepvec,pot,lj206flag,fvec)
    Type(LJ206_Params), Intent(IN)               :: params
    Type(VecType), Intent(IN)                 :: sepvec
    Real(kind = RDbl), Intent(OUT)            :: pot
    Type(VecType), Intent(OUT), Optional      :: fvec
    Logical, Intent(Out)                      :: lj206flag
  
    Real(kind = RDbl)                         :: A,B,r2i,r6i,r20i,r2
  
    !** Initialize the flag
    lj206flag = .True.
  
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
      lj206flag = .False.
      Return
    End If
  
    A = params%A
    B = params%B
  
    r2i = 1.0_RDbl/r2
    r6i = r2i*r2i*r2i
    r20i = r6i*r6i*r6i*r2i
  
    If (Present(fvec)) Then
      fvec = sepvec * (20.0_RDbl*A*r20i - 6.0_RDbl*B*r6i)*r2i
    End If
  
    pot = A*r20i - B*r6i
  
  End Subroutine lj206_calc_snglInteraction
  
  
  !---------------------------------------------------------------------------
  ! Calculate the interaction 
  !---------------------------------------------------------------------------
  Subroutine lj206_calc_snglInteractionHOT(params,sepvec,pot,lj206flag,hot)
    Type(LJ206_Params), Intent(IN)               :: params
    Type(VecType), Intent(IN)                 :: sepvec
    Real(kind = RDbl), Intent(OUT)            :: pot
    Logical, Intent(Out)                      :: lj206flag
    Type(VecType), Dimension(3), Intent(Out)  :: hot
    Real(kind = RDbl)                         :: A,B,r2i,r6i,r20i,r2,r16i

    !MDEBUG
    Write(0,*) __FILE__,__LINE__," Sorry, I was too lazy to finish HOT &
        &for LJ 20-6 interactions."
    Stop
  
    !** Initialize the flag
    lj206flag = .True.
    !** Initialize the potential and (optionally) the force
    pot = 0.0_RDbl
  
    !** Check if the interaction is on
    If (params%off) Return
  
    !** Get the square of the distance
    r2 = vector_getnormsq(sepvec)
  
    !**Check if it is within the cut-off radius
    If (r2 > params%hicut2) Return
    If (r2 < params%locut2) Then
      lj206flag = .False.
      Return
    End If
  
    A = params%A
    B = params%B
  
    r2i = 1.0_RDbl/r2
    r6i = r2i*r2i*r2i
    r16i= r6i*r2i*r6i*r2i
    r20i = r6i*r6i*r6i*r2i
  
    !** Calculate the potential
    pot = A*r20i - B*r6i
  
    !** Calculate the higher order terms
    hot(1) = sepvec * (-20.0_RDbl*A*r20i + 6.0_RDbl*B*r6i)*r2i
!!$    hot(2)%comp(1) = (168.0_RDbl*A*r12i - 48.0_RDbl*B*r6i)*r2i*r2i &
!!$        *sepvec%comp(1)*sepvec%comp(2)
!!$    hot(2)%comp(2) = (168.0_RDbl*A*r12i - 48.0_RDbl*B*r6i)*r2i*r2i &
!!$        *sepvec%comp(1)*sepvec%comp(3)
!!$    hot(2)%comp(3) = (168.0_RDbl*A*r12i - 48.0_RDbl*B*r6i)*r2i*r2i &
!!$        *sepvec%comp(2)*sepvec%comp(3)
!!$    hot(3)%comp(1) = (-2688.0_RDbl*A*r16i + 480.0_RDbl*B*r12i)*r2i &
!!$        *sepvec%comp(1)*sepvec%comp(2)*sepvec%comp(3)
!!$    hot(3)%comp(2) = hot(3)%comp(1)
!!$    hot(3)%comp(3) = hot(3)%comp(1)
  
  End Subroutine lj206_calc_snglInteractionHOT
  
  !---------------------------------------------------------------------------
  ! Calculate the interaction between a molecule and an entire sorbate 
  !---------------------------------------------------------------------------
  Subroutine lj206_calc_multInteraction(params, apairs,sepvec,pot, lj206flag, fvec)
    Type(LJ206_Params), Dimension(:,:), Intent(In) :: params
    Integer, Dimension(:,:), Intent(In)       :: apairs
    Type(VecType), Dimension(:), Intent(IN)   :: sepvec
    Real(kind = RDbl), Intent(OUT) :: pot
    Logical, Intent(out)                      ::  lj206flag   
    Type(VecType), Dimension(:), Intent(OUT), Optional :: fvec
  
    Real(kind = RDbl)                         :: A,B,r2i,r6i,r20i,r2
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
        lj206flag= .False.
        Return
      End If
  
      A = params(a1,a2)%A
      B = params(a1,a2)%B
  
      r2i = 1.0_RDbl/r2
      r6i = r2i*r2i*r2i
      r20i = r6i*r6i*r6i*r2i
  
      If (Present(fvec)) Then
        fvec(i) = sepvec(i) * (20.0_RDbl*A*r20i - 6.0_RDbl*B*r6i)*r2i
      End If
  
      pot = pot + A*r20i - B*r6i
      lj206flag= .True. 
    End Do
  End Subroutine lj206_calc_multInteraction
  
  !---------------------------------------------------------------------------
  ! Calculate the interaction between a molecule and an entire sorbate 
  !---------------------------------------------------------------------------
  Subroutine lj206_calc_multInteractionHOT(params, apairs, sepvec, pot,lj206flag,hot)
    Type(LJ206_Params), Dimension(:,:), Intent(In) :: params
    Integer, Dimension(:,:), Intent(In)         :: apairs
    Type(VecType), Dimension(:), Intent(IN)     :: sepvec
    Real(kind = RDbl), Intent(OUT)              :: pot
    Logical, Intent(out)                        ::  lj206flag   
    Type(VecType), Dimension(3), Intent(OUT)    :: hot
  
    Real(kind = RDbl) :: A,B,r2i,r6i,r12i,r2,r16i
    Integer :: i, a1, a2

    !MDEBUG
    Write(0,*) __FILE__,__LINE__," Sorry, I was too lazy to finish HOT &
        &for LJ 20-6 interactions."
    Stop
  
    !** Zero the forces
    pot = 0.0_RDbl
  
    !** Zero the higher order terms
    hot(1) = 0.0_RDbl
    hot(2) = 0.0_RDbl
    hot(3) = 0.0_RDbl
  
  
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
        lj206flag= .False.
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

      !MDEBUG
 !!$      Write(*,*) "Atom accepted : ",i
 !!$      Write(*,'(2a)') "sepvec        : ",vector_display(sepvec(i),'f15.3')
 !!$      Write(*,*) 'Pcalc: r2, r ',r2,Sqrt(params(a1,a2)%locut2)
 !!$      Write(*,*) "pot           : ",A*r12i-B*r6i
 !!$      Write(*,*) "h1 : ",sepvec(i)*(-12.0_RDbl*A*r12i+6.0_RDbl*B*r6i)*r2i
 !!$        Write(*,*) "h2 : ",(168.0_RDbl*A*r12i-48.0_RDbl*B*r6i)*r2i*r2i
 !!$      Write(*,*) "h3 : ",(-2688.0_RDbl*A*r16i+480.0_RDbl*B*r12i)
 !!$      Write(*,*) "params ",A,B
 !!$      Write(*,*) "r18i, r14i: ",r16i*r2i,r12i*r2i

!!$      !** Calculate the higher order terms
!!$      hot(1) = hot(1) + sepvec(i) * (-12.0_RDbl*A*r12i+6.0_RDbl*B*r6i)*r2i
!!$      hot(2)%comp(1) = hot(2)%comp(1) + (168.0_RDbl*A*r12i-48.0_RDbl*B*r6i) &
!!$          *r2i*r2i*sepvec(i)%comp(1)*sepvec(i)%comp(2)
!!$      hot(2)%comp(2) = hot(2)%comp(2) + (168.0_RDbl*A*r12i-48.0_RDbl*B*r6i) &
!!$          *r2i*r2i*sepvec(i)%comp(1)*sepvec(i)%comp(3)
!!$      hot(2)%comp(3) = hot(2)%comp(3) + (168.0_RDbl*A*r12i-48.0_RDbl*B*r6i) &
!!$          *r2i*r2i*sepvec(i)%comp(2)*sepvec(i)%comp(3)
!!$      hot(3)%comp(1) = hot(3)%comp(1) + (-2688.0_RDbl*A*r16i + &
!!$          480.0_RDbl*B*r12i)*r2i * sepvec(i)%comp(1)*sepvec(i)%comp(2)* &
!!$          sepvec(i)%comp(3)
!!$      hot(3)%comp(2) = hot(3)%comp(1)
!!$      hot(3)%comp(3) = hot(3)%comp(1)

 !!$      !MDEBUG
 !!$      Write(*,*) "hot1 : ",hot(1)
 !!$      Write(*,*) "hot2 : ",hot(2)
      lj206flag= .True. 
    End Do

     !MDEBUG
 !!$    Write(*,*) "accepted : ",accept
 !!$    Write(*,*) "rejected : ",reject
 !!$    Write(*,*) "pot :      ",pot
 !!$    Write(*,*) "hot :      ",vector_display(hot(1),'e12.3')
 !!$    Write(*,*) "hot2:      ",vector_display(hot(2),'e12.3')
 !!$    Write(*,*) "hot3:      ",vector_display(hot(3),'e12.3')

  End Subroutine lj206_calc_multInteractionHOT

  !---------------------------------------------------------------------------
  ! This routine calculates A, B parameters from sigma
  ! and epsilon.  A = 4*eps*(sig**20), B = 4*eps*(sig**6)
  !---------------------------------------------------------------------------
  Subroutine lj206_calc_AB(obj)
    Implicit None
    Type(LJ206_Params), Intent(INOUT)         :: obj
  
    obj%B = 4.0_Rdbl*obj%eps*(obj%sig**6)*kcalmole_kb
    obj%A = obj%b*obj%sig**14
  
 !   Write(*,*) 'LJ206: sig, eps ',lj206_displayparams(obj)
 !   Write(*,*) 'LJ206: A12, B6 ',obj%A,obj%B
  End Subroutine lj206_calc_AB

  !---------------------------------------------------------------------------
  ! This routine calculates the sigma and epsilon parameters from A and B
  !---------------------------------------------------------------------------
  Subroutine lj206_calc_sigeps(obj)
    Implicit None
    Type(LJ206_Params), Intent(INOUT)   :: obj
  
    Real(kind=RDbl)           :: sig2
  
    sig2    = (obj%A/obj%B/obj%B/obj%B)
    obj%eps = (obj%B/(4.0_Rdbl*sig2*sig2*sig2*kcalmole_kb))
    obj%sig = Sqrt(sig2)
  End Subroutine lj206_calc_sigeps
  
  !---------------------------------------------------------------------------
  ! This routine mixes sigma and epsilon values according to the either
  ! the Lorentz-Bertholet mixing rules or the Jorgensen mixing rules
  ! Requires: params -- LJ206 params
  !---------------------------------------------------------------------------
  Subroutine lj206_calc_mix(obj,pure1,pure2,line)
    Type(LJ206_Params), Intent(INOUT)   :: obj
    Type(LJ206_Params), Intent(IN)      :: pure1,pure2
    Character(*), Intent(IN)         :: line
  
    Integer                                   :: nfields,n,i
    Character(len=strLen), Dimension(10)      :: fields,chunks

    !MDEBUG
    Write(0,*) __FILE__,__LINE__," Sorry, too lazy to finish mixing rules &
        &for LJ 20-6"
    Stop
  
    !** Assign the basic init fields first
    nfields = split(line,fields)
    obj%line = line
    obj%off  = .False.
    obj%locut = LOW_DIST_CUTOFF
    obj%atom_name1 = Trim(fields(1))
    obj%atom_name2 = Trim(fields(2))
  
    !** Do the mixing
    If (Index(Toupper(line),'SIG@LBMIX') /= 0) Then
      obj%sig = (pure1%sig + pure2%sig)/2.0_RDbl   !arithmetic
    Else If (Index(Toupper(line),'SIG@JORGMIX') /= 0) Then
      obj%sig = Sqrt(pure1%sig * pure2%sig)        !geometric
    End If
  
    If (Index(Toupper(line),'EPS@LBMIX') /= 0) Then
      obj%eps = Sqrt((pure1%eps * pure2%eps))      !geometric
    Else If (Index(Toupper(line),'EPS@JORGMIX') /= 0) Then
      obj%eps = Sqrt((pure1%eps * pure2%eps))      !geometric
    End If
  
    !** Get the rest of the information from the line
    ! Skip the first 3, those are atom names and LJ206 keyword
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
    Call lj206_calc_AB(obj)
  
    obj%initialized = .True.
    obj%source = "MIXED"
  End Subroutine lj206_calc_mix

  !----------------------------------------------------------------------------
  ! Returns the value of the high end cutoff. If the NEIGHBOR flag is 
  ! present, check to see if the NEIGHBOR cutoff should be returned
  ! Neighbour cut-off is larger value than  high cutoff. It is useful while
  ! making list of intra-pairs in the beginning during molecule 
  ! initialization.
  ! Requires: params -- LJ206 params
  !           neighbor -- flag to get neighbor cutoff
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function lj206_getcutoff(params,neighbor)
    Type(LJ206_Params), Intent(In)         :: params
    Logical, Optional, Intent(In)       :: neighbor

    !** Default value
    lj206_getcutoff = params%hicut

    !** Check if we need to return the neighbor cutoff
    If (Present(neighbor)) Then
      If (neighbor) Then
        lj206_getcutoff = params%ncut
      End If
    End If
  End Function lj206_getcutoff

  !----------------------------------------------------------------------------
  ! Checks to see if two LJ206 type data structures contain the same parameters
  ! Requires:  params1 -- 1st LJ206 params
  !            params2 -- 2nd LJ206 params
  !            tolerance -- optional tolerance for checking
  !----------------------------------------------------------------------------
  Logical Function lj206_chkequal(params1,params2,tolerance)
    Type(LJ206_Params), Intent(In)            :: params1,params2
    Real(kind=RDbl), Intent(In), Optional  :: tolerance

    Real(kind=RDbl)  :: tol,diff

    tol = 1.0e-6_RDbl
    If (Present(tolerance)) tol = tolerance

    lj206_chkequal = .True.

    diff = Abs(params1%A - params2%A)
    If (diff > tol) lj206_chkequal = .False.

    diff = Abs(params1%B - params2%B)
    If (diff > tol) lj206_chkequal = .False.

  End Function lj206_chkequal

  !-------------------------------------------------------
  ! Returns a string with the properly formatted lj206 params
  ! Requires: params -- LJ206 params
  !-------------------------------------------------------
   Function lj206_displayparams(lj206params)
     Character(len=strLen) :: lj206_displayparams
     Type(LJ206_Params), Intent(IN) :: lj206params

     Write(lj206_displayparams,'(1x,a,2x,a)') Trim(real2str(lj206params%sig,5)), &
         Trim(real2str(lj206params%eps,5))
   End Function lj206_displayparams

  !----------------------------------------------------------------------------
  ! Returns a string containing the sigma and epsilon values
  ! Requires: params -- LJ206 params
  !----------------------------------------------------------------------------
  Function lj206_displaySigEps(lj206params)
    Character(len=lstrLen)      :: lj206_displaySigEps
    Type(LJ206_Params), Intent(In) :: lj206params

    If (lj206params%off) Then
      Write(lj206_displaySigEps,'(a)') "OFF"
    Else
      Write(lj206_displaySigEps,'(2(2a,1x,a,3x))') "Sigma: ", &
          Trim(real2str(lj206params%sig,6)),"Ang", &
          "Epsilon/k_Boltz: ",Trim(real2str(lj206params%eps,6)),"K"
    End If

  End Function lj206_displaySigEps

  !----------------------------------------------------------------------------
  ! Returns a string containing nicely formated A and B values
  ! Requires: params -- LJ206 params
  !----------------------------------------------------------------------------
  Function lj206_displayAB(lj206params)
    Character(len=lstrLen)              :: lj206_displayAB
    Type(LJ206_Params), Intent(In)         :: lj206params

    If (lj206params%off) Then
      Write(lj206_displayAB,'(a)') "OFF"
    Else
      Write(lj206_displayAB,'(2(2a,1x,a,3x))') "A: ", &
          Trim(real2str(lj206params%A,8)),"kcal A^20/mol", &
          "B: ",Trim(real2str(lj206params%B,8)),"kcal A^6/mol"
    End If

  End Function lj206_displayAB

  !----------------------------------------------------------------------------
  ! Returns a string containing the cutoff Values
  ! Requires: params -- LJ206 params
  !----------------------------------------------------------------------------
  Function lj206_displayCutOffs(lj206params)
    Character(len=lstrLen)          :: lj206_displayCutOffs
    Type(LJ206_Params), Intent(In)     :: lj206params

    If (lj206params%off) Then
      Write(lj206_displayCutoffs,'(a)') "OFF"
    Else
      Write(lj206_displayCutOffs,'(2(2a,1x,a,3x),f8.4)') &
          "high cutoff: ",Trim(real2str(lj206params%hicut,5)),"Ang", &
          "low cutoff: ",Trim(real2str(lj206params%locut,5)),"Ang"
    End If

  End Function lj206_displayCutOffs

  !-----------------------------------------------------------------------
  ! Displays the full set of LJ206 parameters for the structure
  ! Requires: params -- the LJ206 parameters structure
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !-----------------------------------------------------------------------
  Subroutine lj206_display(params,indent,unit)
    Type(LJ206_Params), Intent(In)         :: params
    Integer, Intent(In)                    :: indent,unit

    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    If (.Not. params%initialized) Return

    If (params%off) Then
      Write(unit,'(5a)') blank,Trim(params%atom_name1),'-', &
          Trim(params%atom_name2),' Interaction OFF'
    Else
      Write(unit,'(6a)') blank,Trim(params%atom_name1),'-',&
          Trim(params%atom_name2),' LJ206 parameters from ',Trim(params%source)
      Write(unit,'(2x,2a)') blank,Trim(lj206_displaySigEps(params))
      Write(unit,'(2x,2a)') blank,Trim(lj206_displayAB(params))
      Write(unit,'(2x,2a)') blank,Trim(lj206_displayCutOffs(params))
    End If

  End Subroutine lj206_display

End Module lj206

