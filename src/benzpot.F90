!------------------------------------------------------------------------------
! This module contains the data structure for Benzene sorbate-sorbate 
! interaction as parametrized by Shi & Bartell.  C and D parameters
! are input in units of kcal A^n/mol.
!------------------------------------------------------------------------------
Module benzpot

  Use defaults, Only: RDbl, strLen, lstrLen, xlstrLen, LOW_DIST_CUTOFF, &
      kcalmole_kb, HIGH_DIST_CUTOFF, dbgflag
  Use utils, Only: toupper, split, toreal, getlineswpair, real2str
  Use file, Only: file_getunit
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnormsq
  Use storebase, Only: EnergyPlus, storebase_inc, storebase_nderivs, &
      storebase_disp

  Implicit None
  Save

  Private
  Public :: Benzpot_Params, benzpot_idstring, benzpot_init, benzpot_isinit, &
      benzpot_snglint, benzpot_display, benzpot_disp, benzpot_getcutoff, &
      benzpot_chkequal

  Type Benzpot_Params
    Character(len=strLen) :: atom_name1,atom_name2
    Character(len=strLen) :: source
    Character(len=100)    :: line
    Character(len=80)     :: reference
    Logical               :: off,initialized 
    Real(kind=RDbl)       :: C,D
    Real(kind=RDbl)       :: hicut,locut
    Real(kind=RDbl)       :: hicut2,locut2
  End Type Benzpot_Params

  Character(len=strLen), Parameter         :: benzpot_idstring = 'BENZPOT'

Contains
  !----------------------------------------------------------------------------
  ! Initialize the interaction information from the file line
  ! Requires:  params -- C&D parameters to initialize
  !            line -- string containing initialization information
  !----------------------------------------------------------------------------
  Subroutine benzpot_init(params,line)
    Type(Benzpot_Params), Intent(INOUT)       :: params
    Character(*), Intent(IN)                  :: line

    Integer                                   :: nfields,n,i,unit,ios
    Integer                                   :: nlines
    Character(len=150), Dimension(10)         :: lines
    Character(len=strLen)                     :: filename,name1,name2
    Character(len=strLen), Dimension(10)      :: fields,chunks

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

      If (filename /= '') Then    !open the file and get the line
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
        If (nlines >= 2) Then
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

    !** set defaults
    nfields = split(lines(1),fields)
    params%line = lines(1)
    params%off = .FALSE.
    params%locut = LOW_DIST_CUTOFF
    params%locut2 = params%locut**2
    params%hicut = HIGH_DIST_CUTOFF
    params%hicut2 = params%hicut**2
    params%atom_name1 = Trim(fields(1))
    params%atom_name2 = Trim(fields(2))
    params%source = ''

    params%C = 0.0_RDbl
    params%D = 0.0_RDbl

    Do i = 4,nfields
      n = split(fields(i),chunks,'@')
      Select Case(toupper(chunks(1)))
        
      Case('OFF')
        params%off = .TRUE.

      Case('C')
        params%C = toreal(chunks(2))

      Case('D')
        params%D = toreal(chunks(2))

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

      Case Default
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
             ' Unable to identify Benzpot interaction string ',trim(chunks(1))
        Stop
      End Select
    EndDo

    !** Make sure that the parameters were read in
    If ((params%C == 0.0_RDbl).Or.(params%D == 0.0_RDbl)) Then
      If (.Not. params%off) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            '  Some parameters not specified '
        Write(0,'(a)') Trim(lines(1))
        Stop
      End If
    End If

    params%initialized = .True.

  End Subroutine benzpot_init

  !-------------------------------------------------------------------
  ! Returns true if the fields of "params" have been initialized
  ! Requires:  params -- C&D parameters 
  !-------------------------------------------------------------------
  Logical Function benzpot_isinit(params)
    Type(Benzpot_Params), Intent(in) :: params
    
    benzpot_isinit = params%initialized
  End Function benzpot_isinit

  !----------------------------------------------------------------------------
  ! Calculate the interaction between a set of atomic points.  Returns False
  ! if any of the seperation vectors fall within the lower cutoff.  This 
  ! indicates atoms have gotten too close together and that the potential 
  ! would be extremely large.  The returned derivative components are those
  ! for the first atom in the separation vector defined by (atom1 - atom2).
  ! The optional second results structure must be passed if forces are to be
  ! calculated.
  ! Requires:  params -- C&D parameters 
  !            sepvec -- atom-atom separation vector
  !            results -- resultant energies and derivatives
  !            results2 -- resultant energies and opposite forces
  !----------------------------------------------------------------------------
  Logical Function benzpot_snglint(params,sepvec,results,results2)
    Type(Benzpot_Params), Intent(In)          :: params
    Type(VecType), Intent(In)                 :: sepvec
    Type(EnergyPlus), Intent(InOut)           :: results
    Type(EnergyPlus), Intent(InOut), Optional :: results2

    Real(kind = RDbl)         :: r2,r2i,nrg
    Type(VecType)             :: force

    benzpot_snglint = .True.

    !** Check if the interaction is on
    If (params%off) Return
    
    !** Get the square of the distance
    r2 = vector_getnormsq(sepvec)
    r2i = 1.0_RDbl/r2

    !** Make sure that the separation is more than a minimum to prevent
    !** overflow
    If (r2 < params%locut2) Then
      benzpot_snglint = .False.
      Return
    End If

    !** Check if it is within the cut-off radius
    If (r2 > params%hicut2) Return

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    dbgflag = .True.

    !** Calculated the potential
    nrg = params%C*r2i + params%D
    Call storebase_inc(results,nrg)

    !** Get forces if desired
    If (storebase_nderivs(results) > 0) Then
      force = sepvec * (2.0_RDbl*params%C*r2i) * r2i
      Call storebase_inc(results,force)
      Call storebase_inc(results2,(force*(-1.0_RDbl)))
    End If

    !** Calculate higher order derivatives
    If (storebase_nderivs(results) > 1) Then
      Write(0,'(2a,i4,a,2l2)') __FILE__,":",__LINE__, &
          " Higher order derivatives not yet available for C&D potential"
      Stop
    End If

    !** If we made it this far, everything is ok
    benzpot_snglint = .True. 

  End Function benzpot_snglint

  !----------------------------------------------------------------------------
  ! Returns the high cutoff of the benzpot potential
  ! Requires:  params -- Benzpotic parameters 
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function benzpot_getcutoff(params)
    Type(Benzpot_Params), Intent(In)      :: params

    benzpot_getcutoff = params%hicut

  End Function benzpot_getcutoff

  !----------------------------------------------------------------------------
  ! Returns a string containing nicely formated C and D values
  ! Requires:  params -- C&D parameters 
  !----------------------------------------------------------------------------
  Function benzpot_disp(params)
    Character(len=xlstrLen)             :: benzpot_disp
    Type(Benzpot_Params), Intent(In)    :: params

    Character(len=strLen)               :: string1,string2

    If (params%off) Then
      Write(benzpot_disp,'(a)') "OFF"
    Else
      string1 = real2str(params%C,8)
      string2 = real2str(params%D,8)
      Write(benzpot_disp,'(2(2a,1x),a)') &
          "C: ",Trim(string1),"D: ",Trim(string2)," kcal A^(12,6,2,0)/mol"
    End If

  End Function benzpot_disp

  !----------------------------------------------------------------------------
  ! Returns a string containing the cutoff Values
  ! Requires:  params -- C&D parameters 
  !----------------------------------------------------------------------------
  Function benzpot_displayCutOffs(params)
    Character(len=xlstrLen)          :: benzpot_displayCutOffs
    Type(Benzpot_Params), Intent(In) :: params

    Character(len=strLen)               :: string1,string2

    If (params%off) Then
      Write(benzpot_displayCutoffs,'(a)') "OFF"
    Else
      string1 = real2str(params%hicut,4)
      string2 = real2str(params%locut,4)
      Write(benzpot_displayCutOffs,'(3(2a,1x,a,3x))') &
          "high cutoff: ",Trim(string1),"Ang", &
          "low cutoff: ",Trim(string2),"Ang"
    End If

  End Function benzpot_displayCutOffs

  !----------------------------------------------------------------------------
  ! Checks to see if two CD data structures contain the same parameters
  ! Requires:  params1 -- 1st params
  !            params2 -- 2nd params
  !            tolerance -- optional tolerance for checking
  !----------------------------------------------------------------------------
  Logical Function benzpot_chkequal(params1,params2,tolerance)
    Type(Benzpot_Params), Intent(In)        :: params1,params2
    Real(kind=RDbl), Intent(In), Optional   :: tolerance

    Real(kind=RDbl)  :: tol,diff

    tol = 1.0e-6_RDbl
    If (Present(tolerance)) tol = tolerance

    benzpot_chkequal = .True.

    diff = Abs(params1%C - params2%C)
    If (diff > tol) benzpot_chkequal = .False.

    diff = Abs(params1%D - params2%D)
    If (diff > tol) benzpot_chkequal = .False.

  End Function benzpot_chkequal

  !----------------------------------------------------------------------------
  ! Display the pertinent information from the Benzpot structure
  ! Requires: params -- the benzene potentials parameters structure
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !----------------------------------------------------------------------------
  Subroutine benzpot_display(params,indent,unit)
    Type(Benzpot_Params), Intent(In)         :: params
    Integer, Intent(In)                 :: indent,unit

    Character(len=indent)       :: blank
    Character(len=xlstrLen)     :: string

    blank = Repeat(' ',indent)

    If (params%off) Then
      Write(unit,'(2a,2x,2a)') blank,Trim(params%atom_name1),'-', &
          Trim(params%atom_name2),' Interaction OFF'
    Else
      Write(unit,'(6a)') blank,Trim(params%atom_name1),'-',&
          Trim(params%atom_name2),' BENZPOT parameters from ', &
          Trim(params%source)
      string = benzpot_disp(params)
      Write(unit,'(2x,2a)') blank,Trim(string)
      string = benzpot_displayCutOffs(params)
      Write(unit,'(2x,2a)') blank,Trim(string)
    End If
    
  End Subroutine benzpot_display

End Module benzpot
