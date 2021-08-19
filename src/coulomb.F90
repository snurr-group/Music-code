!------------------------------------------------------------------------------
! This module contains the data structure and routines for evaulating simple
! Coulomb electrostatic interactions.  It is similar to the "lj" module,
! but instead includes electrostatic interactions.
!
! Needed Improvements
! 1) add smoothing features
! 2) add possibility to use passed charges, currently too inflexible
! 3) NEEDS shifted coulombic potentials
!------------------------------------------------------------------------------

Module coulomb

  Use defaults, Only: RDbl, strLen, lstrLen, LOW_DIST_CUTOFF, kcalmole_kb, &
      xlstrLen, HIGH_DIST_CUTOFF, zero
  Use utils, Only: toupper, split, toreal, getlineswpair, real2str
  Use file, Only: file_getunit
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnormsq 
  Use storebase, Only: EnergyPlus, storebase_inc, storebase_nderivs, &
      storebase_disp

  Implicit None
  Save

  Private
  Public :: Coulomb_Params, coulomb_display, coulomb_snglint, coulomb_disp, &
      coulomb_idstring, coulomb_isinit, coulomb_init, coulomb_displayCutOffs, &
      coulomb_getcutoff, coul_chkequal

  Type Coulomb_Params
    Character(len=strLen) :: atom_name1, atom_name2
    Character(len=strLen) :: source
    Logical               :: off,initialized
    Character(len=100)    :: line
    Character(len=80)     :: reference
    Real(kind=RDbl)       :: shield_const
    Real(kind=RDbl)       :: locut, locut2   ! lower bound to prevent overflow
    Real(kind=RDbl)       :: hicut, hicut2   ! upper bound cutoff

    !** This stuff not yet in use
    Real(kind=RDbl)       :: smoothrad,smoothrad2
    Real(kind=RDbl)       :: eler     !errors at the truncation radius
    Real(kind=RDbl)       :: sradel   !pots at smoothing radius
    Real(kind=RDbl)       :: sraddel  !pot derivatives at smoothing radius

    !** for charges taken from aa_file
    Real(kind=RDbl)       :: Q1, Q2
    Logical               :: charges_fromAA_FILE
  End Type Coulomb_Params

  Character(len=strLen), Parameter         :: coulomb_idstring = 'COUL'

Contains
  !----------------------------------------------------------------------------
  ! Initialize the interaction information from the file line
  ! Requires:  params -- Coulombic parameters to initialize
  !            line -- string containing initialization information
  !----------------------------------------------------------------------------
  Subroutine coulomb_init(params,line)
    Type(Coulomb_Params), Intent(INOUT)       :: params
    Character(*), Intent(IN)                  :: line

    Integer                                   :: nfields,n,i,unit,ios
    Integer                                   :: nlines
    Character(len=150), Dimension(20)         :: lines
    Character(len=strLen)                     :: filename,name1,name2
    Character(len=strLen), Dimension(20)      :: fields,chunks
    Logical :: foundQ1, foundQ2
    Write(0,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(0,*) 'WARNING: coulombic interactions not checked against old code'

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

    !** Set the defaults
    params%line = lines(1)
    params%off = .FALSE.
    params%locut = LOW_DIST_CUTOFF
    params%locut2 = LOW_DIST_CUTOFF*LOW_DIST_CUTOFF
    params%hicut = zero
    params%hicut2 = zero
    params%atom_name1 = Trim(fields(1))
    params%atom_name2 = Trim(fields(2))
    params%shield_const = 1.0_RDbl

    params%charges_fromAA_FILE = .False.
    params%Q1 = zero 
    params%Q2 = zero
    foundQ1=.False.
    foundQ2=.False.

    !** Read and interpret the line
    nfields = split(lines(1),fields)
    Do i = 4,nfields
      n = split(fields(i),chunks,'@')

      Select Case(toupper(chunks(1)))
      Case('OFF')
        params%off = .TRUE.

      Case('Q1')
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            " WARNING: charges taken from atm-atm file, make sure you &
            & intended this"
        params%Q1 = toreal(chunks(2))
        foundQ1=.True.
      Case('Q2')
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            " WARNING: charges taken from atm-atm file, make sure you &
            & intended this"
        params%Q2 = toreal(chunks(2))
        foundQ2=.True.
      Case('TRUNC')
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            " WARNING: Use of the label TRUNC has been depreciated in favor", &
            " of HICUT"
        params%hicut = toreal(chunks(2))
        params%hicut2 = params%hicut**2
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
           "Using HiCUT for atom-atom coulombic interactions will require shifted potential of"
         Write(0,'(2x,a)') " Wolf et al., J. Chem. Phys., 110, 8254, 1999 "
!       Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
!            " ERROR: using an coulombic atom-based cut-off will produce garbage."
!        Write(0,'(2x,a)') "Specify a molecule-based cut-off in spc-spc file"
        Stop

      Case('HICUT')
        params%hicut = toreal(chunks(2))
        params%hicut2 = params%hicut**2
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
           "Using HICUT for atom-atom coulombic interactions will require shifted potential of"
         Write(0,'(2x,a)') " Wolf et al., J. Chem. Phys., 110, 8254, 1999 "

!        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
!            " ERROR: using an coulombic atom-based cut-off will produce garbage."
!        Write(0,'(2x,a)') "Specify a molecule-based cut-off in spc-spc file"
!        Stop

      Case('LOWCUT') 
        params%locut = toreal(chunks(2))
        params%locut2 = params%locut**2
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            " ERROR: using an coulombic atom-based cut-off will produce garbage."
        Write(0,'(2x,a)') "Specify a molecule-based cut-off in spc-spc file"
        Stop

      Case('LOCUT') 
        params%locut = toreal(chunks(2))
        params%locut2 = params%locut**2
        Write(0,'(2a,i4,2a)') __FILE__,":",__LINE__, &
            " ERROR: using an coulombic atom-based cut-off will produce garbage."
        Write(0,'(2x,a)') "Specify a molecule-based cut-off in spc-spc file"
        Stop

      Case('SHIELDCONST')
        params%shield_const = toreal(chunks(2))

      Case('SMOOTHRAD')
        Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
            " Electrostatic interaction smoothing not yet available"
        Write(0,'(2x,a)') "Specify molecule-based smoothing in spc-spc file"
        Stop

      Case Default
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
             ' Unable to identify Coulomb interaction string ',trim(chunks(1))
        Stop
      End Select
    End Do

    If (foundQ1.And.foundQ2) Then
    params%charges_fromAA_FILE = .True.      
    Else
      If (foundQ1.Or.foundQ2) Then
        Write(*,*) " Only one charge is specified in aa file"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Endif

    !** Set the source
    params%initialized = .True.
    params%source = "FILE"

  End Subroutine coulomb_init

  !-------------------------------------------------------------------
  ! Returns true if the fields of "params" have been initialized
  ! Requires:  params -- Coulombic parameters 
  !-------------------------------------------------------------------
  Logical Function coulomb_isinit(params)
    Type(Coulomb_Params), Intent(in) :: params
    
    coulomb_isinit = params%initialized
  End Function coulomb_isinit

  !---------------------------------------------------------------------------
  ! Calculate the interaction between a set of atomic points.  Returns False
  ! if any of the seperation vectors fall within the lower cutoff.  This 
  ! indicates atoms have gotten too close together and that the potential 
  ! would be extremely large.  The returned derivative components are those
  ! for the first atom in the separation vector defined by (atom1 - atom2).
  ! The optional second results structure must be passed if forces are to be
  ! calculated.
  ! Requires:  params -- Coulombic parameters 
  !            sepvec -- atom-atom separation vector
  !            charge1 -- charge of atom1
  !            charge2 -- charge of atom2
  !            results -- resultant energies and derivatives
  !            results2 -- resultant energies and opposite forces
  !---------------------------------------------------------------------------
  Logical Function coulomb_snglint(params,sepvec,charge1_in,charge2_in, &
      results,results2)
    Type(Coulomb_Params), Intent(In)          :: params
    Type(VecType), Intent(In)                 :: sepvec
    Real(kind=RDbl), Intent(In)               :: charge1_in,charge2_in
    Type(EnergyPlus), Intent(InOut)           :: results
    Type(EnergyPlus), Intent(InOut), Optional :: results2

    Real(kind=RDbl)           :: mult_fac
    Real(kind = RDbl)         :: r2,r,nrg
    Type(VecType)             :: force

    !** Get the square of the distance
    r2 = vector_getnormsq(sepvec)

    !** Check if it is within the cut-off radius
    If(params%hicut> zero) Then
    If (r2 > params%hicut2) Then
      coulomb_snglint = .True.
      Return
    End If
   End IF

    !** Check if it is within the cut-off radius
    If (r2 < params%locut2) Then
      coulomb_snglint = .False.
      Return
    End If

    r = Sqrt(r2)

    !** decide which charge to use
    If(params%hicut> zero) Then
    If (params%charges_fromAA_FILE) then
      mult_fac=(331.5 * params%shield_const * params%Q1* params%Q2)*(1.0_RDbl/r -1.0_RDbl/params%hicut)
    Else
      mult_fac=(331.5 * params%shield_const * charge1_in* charge2_in)*(1.0_RDbl/r -1.0_RDbl/params%hicut)
    Endif
    Else
    If (params%charges_fromAA_FILE) then
      mult_fac=(331.5 * params%shield_const * params%Q1* params%Q2)*(1.0_RDbl/r)
    Else
      mult_fac=(331.5 * params%shield_const * charge1_in* charge2_in)*(1.0_RDbl/r)
    Endif
    Endif
    !** conversion constant(from atomic unit to kcal/mol) 
    !** = 9.0*10^9*(1.6*10^-19)^2/10^(-10) * 6.02*10^23 /4.184 = 331 kcal/mol

    nrg = mult_fac

    Call storebase_inc(results,nrg)

    !** Get forces if desired
    If (storebase_nderivs(results) > 0) Then
      force = sepvec * (mult_fac/r2)
      Call storebase_inc(results,force)
      Call storebase_inc(results2,(force*(-1.0_RDbl)))
    End If

    !** Calculate higher order derivatives
    If (storebase_nderivs(results) > 1) Then
      Write(0,'(2a,i4,a,2l2)') __FILE__,":",__LINE__, &
          " Higher order derivatives not yet available for Coulomb potential"
      Stop
    End If

    !** If we made it this far, everything is ok
    coulomb_snglint = .True. 

  End Function coulomb_snglint

  !----------------------------------------------------------------------------
  ! Returns the high cutoff of the coulomb potential
  ! Requires:  params -- Coulombic parameters 
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function coulomb_getcutoff(params)
    Type(Coulomb_Params), Intent(In)      :: params

    coulomb_getcutoff = params%hicut

  End Function coulomb_getcutoff

  !----------------------------------------------------------------------------
  ! Checks to see if two COUL type data structures contain the same parameters
  ! Requires:  params1 -- 1st COUL params
  !            params2 -- 2nd COUL params
  !            tolerance -- optional tolerance for checking
  !----------------------------------------------------------------------------
  Logical Function coul_chkequal(params1,params2,tolerance)
    Type(Coulomb_Params), Intent(In)       :: params1,params2
    Real(kind=RDbl), Intent(In), Optional  :: tolerance

    Real(kind=RDbl)  :: tol,diff

    tol = 1.0e-6_RDbl
    If (Present(tolerance)) tol = tolerance

    coul_chkequal = .True.

    diff = Abs(params1%shield_const - params2%shield_const)
    If (diff > tol) coul_chkequal = .False.

  End Function coul_chkequal

  !----------------------------------------------------------------------------
  ! Display the pertinent information from the Coulomb structure
  ! Requires:  params -- Coulombic parameters 
  !            indent -- number of spaces to indent
  !            unit -- unit to write into
  !----------------------------------------------------------------------------
  Subroutine coulomb_display(params,indent,unit)
    Type(Coulomb_Params), Intent(In)      :: params
    Integer, Intent(In)                   :: indent,unit

    Character(len=indent)       :: blank

    blank = Repeat(' ',indent)

    If (.Not. params%initialized) Return
    
    Write(unit,'(2a)') blank,Trim(params%line)

  End Subroutine coulomb_display

  !----------------------------------------------------------------------------
  ! Display the pertinent information from the Coulomb structure on one line
  ! Requires:  params -- Coulombic parameters 
  !----------------------------------------------------------------------------
  Function coulomb_disp(params)
    Character(len=xlstrLen)               :: coulomb_disp
    Type(Coulomb_Params), Intent(In)      :: params

    Character(len=strLen)        :: string

    If (.Not. params%initialized) Then
      Write(coulomb_disp,'(a)') "OFF"
      Return
    End If

    string = real2str(params%shield_const,5)
    If (params%charges_fromAA_FILE) Then
      Write(coulomb_disp,'(a,2f10.5,2a)') 'Charges taken from aa_file : ,', &
          params%Q1, params%Q2, ' Shielding const: ',Trim(string)
    Else
      Write(coulomb_disp,'(3a)') 'Charges taken from configuration,', &
          ' Shielding const: ',Trim(string)
    Endif

  End Function coulomb_disp

  !----------------------------------------------------------------------------
  ! Returns a string containing the cutoff Values
  ! Requires:  params -- Coulombic parameters 
  !----------------------------------------------------------------------------
  Function coulomb_displayCutOffs(params)
    Character(len=lstrLen)           :: coulomb_displayCutOffs
    Type(Coulomb_Params), Intent(In) :: params

    Write(coulomb_displayCutOffs,'(a)') 'see spc-spc section for cut-offs'
    Return

    Write(coulomb_displayCutOffs,'(2(a,f8.3,1x,a,3x),f8.4)') &
        "HICUT : ",params%hicut,"Ang","LOCUT : ",params%locut,"Ang"

  End Function coulomb_displayCutOffs

End Module coulomb




