!------------------------------------------------------------------------------
! This module contains the data structure and routines for evaulating Lennard
! Jones and electrostatic interactions.  It is similar to the "lj" module,
! but it also includes electrostatic interactions.
!
! NOTE: the fact that the Lennard-Jones interactions are repeated in this
! module is unsatisfying.  If the data and call structure in "pairmodel"
! were better, this would not be necessary. 
!------------------------------------------------------------------------------

Module ljel

  Use defaults, Only: RDbl, strLen, lstrLen, LOW_DIST_CUTOFF, kcalmole_kb
  Use utils, Only: toupper, split, toreal, getlineswpair
  Use file, Only: file_getunit
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnormsq 
  Use atom, Only: atom_getmass

  Implicit None
  Save

  Private
  Public :: LJEL_Params, ljel_displayparams, ljel_calc_interaction, &
      ljel_idstring, ljel_isinit, ljel_displaySigEps, ljel_initparams, &
      ljel_calc_mix, ljel_calc_ab, ljel_init, ljel_displayCutOffs

  Interface ljel_calc_interaction
    Module procedure ljel_calc_multInteraction
  End Interface

  Type LJEL_Params
    Character(len=strLen) :: atom_name1, atom_name2
    Character(len=strLen) :: source
    Logical             :: off
    Logical             :: initialized
    Real(kind=RDbl)     :: sig, eps
    Real(kind=RDbl)     :: atom_charge1,atom_charge2   
    Real(kind=RDbl)     :: A, B, C, D
    Character(len=100)  :: line
    Character(len=80)   :: reference
    Real(kind=RDbl)     :: shield_const
    Real(kind=RDbl)     :: ljcutrad,ljcutrad2
    Real(kind=RDbl)     :: elcutrad,elcutrad2
    Real(kind=RDbl)     :: lowcut, lowcut2 ! lower bound to prevent overflow
    Real(kind=RDbl)     :: smoothrad,smoothrad2
    Real(kind=RDbl)     :: ljeler,eler       !errors at the truncation radius
    Real(kind=RDbl)     :: sradljel,sradel   !pots at smoothing radius
    Real(kind=RDbl)  :: sraddljel,sraddel !pot derivatives at smoothing radius
  End Type LJEL_Params

  Character(len=strLen), Parameter         :: ljel_idstring = 'LJEL'


Contains
  !----------------------------------------------------------------------------
  ! Initialize the interaction information from the file line
  !----------------------------------------------------------------------------
  Subroutine ljel_init(params,line)
    Type(LJEL_Params), Intent(INOUT)            :: params
    Character(*), Intent(IN)                  :: line

    Integer                                   :: nfields,n,i,unit,ios
    Integer                                   :: nlines
    Character(len=150), Dimension(20)         :: lines
    Character(len=strLen)                     :: filename,name1,name2
    Character(len=strLen), Dimension(20)      :: fields,chunks


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
    params%lowcut = LOW_DIST_CUTOFF
    params%lowcut2 = LOW_DIST_CUTOFF*LOW_DIST_CUTOFF
    params%atom_name1 = Trim(fields(1))
    params%atom_name2 = Trim(fields(2))

    Do i = 4,nfields
      n = split(fields(i),chunks,'@')

      Select Case(toupper(chunks(1)))
      Case('OFF')
        params%off = .TRUE.
      Case('SIG')
        If (Index(chunks(2),'MIX') == 0) Then
          params%sig = toreal(chunks(2))
        End If
      Case('EPS')
        If (Index(chunks(2),'MIX') == 0) Then
          params%eps = toreal(chunks(2))
        End If
      Case('LJTRUNC')
        params%ljcutrad = toreal(chunks(2))
        params%ljcutrad2 = params%ljcutrad**2
      Case('LJLOWCUT') 
        params%lowcut = toreal(chunks(2))
        params%lowcut2 = params%lowcut**2
      Case('Q1')
          params%atom_charge1 = toreal(chunks(2))
      Case('Q2')
          params%atom_charge2 = toreal(chunks(2))
      Case('ELTRUNC')
        params%elcutrad = toreal(chunks(2))
        params%elcutrad2 = params%elcutrad**2 
      Case('SHIELDCONST')
        params%shield_const = toreal(chunks(2))
      Case Default
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
             'Unable to identify LJEL interaction string ',trim(chunks(1))
        Stop
      End Select

    Enddo

    !** Generate the AB parameters
    Call ljel_calc_AB(params)

    !** Set the source
    params%initialized = .True.
    params%source = "SIGEPS"
  End Subroutine ljel_init

  !-------------------------------------------------------------------
  ! Returns true if the fields of "params" have been initialized
  !-------------------------------------------------------------------
  Logical Function ljel_isinit(params)
    Type(LJEL_Params), Intent(in) :: params
    
    ljel_isinit = params%initialized
  End Function ljel_isinit

  !-------------------------------------------------------------------
  ! This is the init routine that conforms to the intramolecular model
  !-------------------------------------------------------------------
  Subroutine ljel_initparams(ljelparams,params)

    Type(LJEL_Params), Intent(INOUT) :: ljelparams
    Character(*), Dimension(:) :: params

    ljelparams%sig = toreal(params(1))
    ljelparams%eps = toreal(params(2))
    !ljelparams%cutrad = toreal(params(3))
    !ljelparams%smoothrad = toreal(params(4))
    ljelparams%initialized = .True.
  End Subroutine ljel_initparams

  !----------------------------------------------------------------------------
  ! Calculate the interaction between a molecule and an entire sorbate 
  ! For some reason the interactions between atoms are calculated if the 
  ! molecules ( COM) are within the cutoff radius. 
  ! The cutoff is not based on an individual atom-atom distance
  ! This necessitates that all atom_atom cut offs (LJEL) should be same
  !----------------------------------------------------------------------------
  Subroutine ljel_calc_multInteraction(params, apairs, sepvec, pot, ljelflag, &
             ctcdis2, fvec)  
    Type(LJEL_Params), Dimension(:,:), Intent(In) :: params
    Integer, Dimension(:,:), Intent(In)       :: apairs
    Type(VecType), Dimension(:), Intent(IN)   :: sepvec
    Real(kind = RDbl), Dimension(:), Intent(IN)  :: ctcdis2
    Real(kind = RDbl), Intent(OUT) :: pot
    Real(kind = RDbl) :: shield_const  
    Logical, Intent(out)                      ::  ljelflag   
    Type(VecType), Dimension(:), Intent(OUT), Optional :: fvec

    Real(kind = RDbl)                :: A,B,r2i,r6i,r12i,r2,r 
    Integer :: i, a1, a2

    !** Zero the forces
    pot = 0.0_RDbl
    If (present(fvec)) fvec = VecType(0.0_RDbl)

    !** Loop through all the sepvecs
    Do i = 1, Size(sepvec)

      !** Parameter indices
      a1 = apairs(1,i)
      a2 = apairs(2,i)

      !** Assign the shielding constant for this interaction
      shield_const = params(a1,a2)%shield_const

      !** Get the square of the distance
      r2 = vector_getnormsq(sepvec(i))

      !**Check if it is within the cut-off radius
      If (r2 < params(a1,a2)%lowcut2) Then 
        ljelflag= .False.
        Return
      End If

      ! LJ interactions first      
      If (r2 < params(a1,a2)%ljcutrad2) Then 
        A = params(a1,a2)%A
        B = params(a1,a2)%B
        r2i = 1.0_RDbl/r2
        r6i = r2i*r2i*r2i
        r12i= r6i*r6i
      
        If (Present(fvec)) Then
          fvec(i) = sepvec(i) * (12.0_RDbl*A*r12i - 6.0_RDbl*B*r6i)*r2i
        End If
        pot = pot + A*r12i - B*r6i
      End If


      ! Electrostatic interactions second 
      ! the cutoff is judged by atom-atom distance, not by mol-mol dis.  
      If (ctcdis2(i) > params(a1,a2)%elcutrad2) Cycle
      !** Get the charges
      r = r2 ** 0.5 

      ! conversion constant(from atomic unit to kcal/mol)
      ! =9.0*10^9*(1.6*10^-19)^2/10^(-10)    *6.02*10^23 /4.184 =331 kcal/mol
      pot = 331.5*shield_const*params(a1,a2)%atom_charge1 * &
            params(a1,a2)%atom_charge2/r + pot
      If (Present(fvec)) Then
        fvec(i) = sepvec(i) * shield_const *331.5_RDbl* &
                  params(a1,a2)%atom_charge1*params(a1,a2)%atom_charge2/(r*r2) &
                  + fvec(i)
      End If 

      ljelflag= .True. 
    End Do
  End Subroutine ljel_calc_multInteraction

  !-------------------------------------------------------
  ! Returns a string with the properly formatted ljel params
  !-------------------------------------------------------
!!$  Character(len=strLen) Function ljel_displayparams(ljelparams)
  Function ljel_displayparams(ljelparams)
    Character(len=strLen) :: ljel_displayparams
    Type(LJEL_Params), Intent(IN) :: ljelparams

    Write(ljel_displayparams,'(1x,e10.4,2x,e10.4)') ljelparams%sig, &
        ljelparams%eps
  End Function ljel_displayparams

  !----------------------------------------------------------------------------
  ! This routine calculates A, B parameters from sigma
  ! and epsilon.  A = 4*eps*(sig**12), B = 4*eps*(sig**6)
  !----------------------------------------------------------------------------
  Subroutine ljel_calc_AB(obj)
    Implicit None
    Type(LJEL_Params), Intent(INOUT)         :: obj

    obj%B = 4.0_Rdbl*obj%eps*(obj%sig**6)*kcalmole_kb
    obj%A = obj%b*obj%sig**6

!    Write(*,*) 'LJEL: sig, eps ',ljel_displayparams(obj)
!    Write(*,*) 'LJEL: A12, B6 ',obj%A,obj%B
  End Subroutine ljel_calc_AB


  !----------------------------------------------------------------------------
  ! This routine calculates the sigma and epsilon parameters from A and B
  !----------------------------------------------------------------------------
  Subroutine ljel_calc_sigeps(obj)
    Implicit None
    Type(LJEL_Params), Intent(INOUT)   :: obj

    Real(kind=RDbl)           :: sig6

    sig6    = (obj%A/obj%B)
    obj%eps = (obj%B/(4.0_Rdbl*sig6*kcalmole_kb))
    obj%sig = sig6**(1/6.0)
  End Subroutine ljel_calc_sigeps


  !----------------------------------------------------------------------------
  ! This routine mixes sigma and epsilon values according to the either
  ! the Lorentz-Bertholet mixing rules or the Jorgensen mixing rules
  !----------------------------------------------------------------------------
  Subroutine ljel_calc_mix(obj,pure1,pure2,line)
    Implicit None
    Type(LJEL_Params), Intent(INOUT)   :: obj
    Type(LJEL_Params), Intent(IN)      :: pure1,pure2
    Character(*), Intent(IN)         :: line

    Integer                                   :: nfields,n,i
    Character(len=strLen), Dimension(10)      :: fields,chunks

    !** Assign the basic init fields first
    nfields = split(line,fields)
    obj%line = line
    obj%off  = .False.
    obj%lowcut = LOW_DIST_CUTOFF
    obj%atom_name1 = Trim(fields(1))
    obj%atom_name2 = Trim(fields(2))
    obj%atom_charge1 = pure1%atom_charge1 
    obj%atom_charge2 = pure2%atom_charge1  
    obj%elcutrad = pure1%elcutrad      
    obj%elcutrad2 = pure1%elcutrad2     

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
    ! Skip the first 3, those are atom names and LJEL keyword
    Do i = 4,nfields
      n = split(fields(i),chunks,'@')
      Select Case(toupper(chunks(1)))
      Case('LJTRUNC')
        obj%ljcutrad = toreal(chunks(2))
        obj%ljcutrad2 = obj%ljcutrad**2
      Case('LJSMOOTH')
        obj%smoothrad = toreal(chunks(2))
        obj%smoothrad2 = obj%smoothrad**2
      Case('LJLOWCUT')
        obj%lowcut = toreal(chunks(2))
        obj%lowcut2 = obj%lowcut**2 
      End Select
    End Do

    !** Generate the AB parameters
    Call ljel_calc_AB(obj)

    obj%initialized = .True.
    obj%source = "MIXED"
  End Subroutine ljel_calc_mix

  !----------------------------------------------------------------------------
  ! Display the pertinent information from the LJEL structure
  !----------------------------------------------------------------------------
  Subroutine ljel_display(obj,unit)
    Implicit None
    Type(LJEL_Params), Intent(IN)      :: obj
    Integer, Intent(IN)              :: unit
    
    Write(unit,'(4x,a)') Trim(obj%line)
  End Subroutine ljel_display

  !----------------------------------------------------------------------------
  ! Returns a string containing the sigma and epsilon values
  !----------------------------------------------------------------------------
!!$  Character(len=lstrLen) Function ljel_displaySigEps(ljelparams)
  Function ljel_displaySigEps(ljelparams)
    Character(len=lstrLen) :: ljel_displaySigEps
    Type(LJEL_Params), Intent(In) :: ljelparams

    Write(ljel_displaySigEps,'(2(a,f8.3,1x,a,3x))') "Sigma : ", &
        ljelparams%sig, &
        "Ang","Epsilon/k(Boltz) : ",ljelparams%eps,"K"
  End Function ljel_displaySigEps

  !----------------------------------------------------------------------------
  ! Returns a string containing the cutoff Values
  !----------------------------------------------------------------------------
!!$  Character(len=lstrLen) Function ljel_displayCutOffs(ljelparams)
  Function ljel_displayCutOffs(ljelparams)
    Character(len=lstrLen) :: ljel_displayCutOffs
    Type(LJEL_Params), Intent(In) :: ljelparams

    Write(ljel_displayCutOffs,'(2(a,f8.3,1x,a,3x),f8.4)') &
        "Trunc : ", ljelparams%ljcutrad,"Ang","LOWCUT : ",ljelparams%lowcut,"Ang"
  End Function ljel_displayCutOffs

End Module ljel

