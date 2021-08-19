!------------------------------------------------------------------------------
! This module contains the routines specific to initializing and actually 
! calculating COULOMBIC PAIWISE INTERACTIONS.  Interactions can be calculated
! either between a single pair (coul_calc_snglInteraction) or for a set 
! of pairs (coul_calc_multInteraction).    The
! data structure (Coul_Params) contains all the potential parameterization
! information (low and high cuts).  It is initialized using file information in "coul_init"
! 


Module coul

  Use defaults, Only: RDbl, strLen, lstrLen, LOW_DIST_CUTOFF, kcalmole_kb, &
      HIGH_DIST_CUTOFF, xlstrlen
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
  Public :: Coul_Params, coul_displayparams, coul_calc_interaction, coul_idstring, &
      coul_isinit,  &
      coul_init, coul_displayCutOffs,  &
      coul_getcutoff, coul_display, coul_grpint, coul_snglint, coul_grpint2


  Interface coul_calc_interaction
    Module procedure coul_calc_snglInteraction
    Module procedure coul_calc_multInteraction
  End Interface


  Type Coul_Params
    Character(len=strLen) :: atom_name1, atom_name2
    Integer              :: atom_number1, atom_number2
    Character(len=strLen) :: source
    Logical             :: off
    Logical             :: initialized
    Character(len=100)  :: line
    Real(kind=RDbl)     :: hicut,hicut2
    Real(kind=RDbl)     :: locut, locut2 ! lower bound to prevent overflow
     Real(kind=RDbl)     :: ncut, ncut2
    Real(kind=RDbl)       :: shield_const
   Real(kind=RDbl)     ::  charge1, charge2
  End Type Coul_Params

  Character(len=strLen), Parameter         :: coul_idstring = 'COUL'

!  Real(kind=RDbl), Parameter               :: LJ_NCUTOFF = 200.00_RDbl


Contains
  !----------------------------------------------------------------------------
  ! Initialize the interaction information from the file line
  ! Requires:  params -- Coul parameter structure to initialize
  !            line -- input command line to interpret
  !----------------------------------------------------------------------------
  Subroutine coul_init(params,line)
    Type(Coul_Params), Intent(InOut)     :: params
    Character(*), Intent(In)           :: line

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
    params%atom_name1 = Trim(fields(1))
    params%atom_name2 = Trim(fields(2))
    params%shield_const = 1.0_RDbl

    Do i = 4,nfields
      n = split(fields(i),chunks,'@')
      Select Case(toupper(chunks(1)))
        

      Case('HICUT')
        params%hicut = toreal(chunks(2))
        params%hicut2 = params%hicut**2


      Case('LOCUT') 
        params%locut = toreal(chunks(2))
        params%locut2 = params%locut**2


      Case Default
        Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
            'Unable to identify Coulombic interaction string ',Trim(chunks(1))
        Stop
      End Select
    End Do


    params%initialized = .True.

   End Subroutine coul_init

  !-------------------------------------------------------------------
  ! Returns true if the fields of "params" have been initialized
  !-------------------------------------------------------------------
  Logical Function coul_isinit(params)
    Type(Coul_Params), Intent(in) :: params
  
    coul_isinit = params%initialized
  End Function coul_isinit


  !--------------------------------------------------------------------------
  ! Calculate the interaction 
  !--------------------------------------------------------------------------
  Subroutine coul_calc_snglInteraction(params,sepvec,pot,coulflag,fvec)
    Type(Coul_Params), Intent(IN)               :: params
    Type(VecType), Intent(IN)                 :: sepvec
    Real(kind = RDbl), Intent(OUT)            :: pot
    Type(VecType), Intent(OUT), Optional      :: fvec
    Logical, Intent(Out)                      :: coulflag
  
    Real(kind = RDbl)                         :: r2i,r6i,r12i,r2,r
  
    !** Initialize the flag
    coulflag = .True.
  
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
      coulflag = .False.
      Return
    End If
  
    r = Sqrt(r2)

    If (Present(fvec)) Then
      fvec = sepvec * params%shield_const * 331.5_RDbl * &
                params%charge1 * params%charge2/(r*r2)
    End If

      pot= 331.5 * params%shield_const * params%charge1 * params%charge2/r 


  End Subroutine coul_calc_snglInteraction
  
  !---------------------------------------------------------------------------
  ! Calculate the interaction between a molecule and an entire sorbate 
  !---------------------------------------------------------------------------
  Subroutine coul_calc_multInteraction(params, apairs,sepvec,pot, coulflag, fvec)
    Type(Coul_Params), Dimension(:,:), Intent(In) :: params
    Integer, Dimension(:,:), Intent(In)       :: apairs
    Type(VecType), Dimension(:), Intent(IN)   :: sepvec
    Real(kind = RDbl), Intent(OUT) :: pot
    Logical, Intent(out)                      ::  coulflag   
    Type(VecType), Dimension(:), Intent(OUT), Optional :: fvec
  
    Real(kind = RDbl)                         :: r2i,r6i,r12i,r2,r
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
        coulflag= .False.
        Return
      End If
  
      r = Sqrt(r2)
  
      If (Present(fvec)) Then
        fvec(i) = sepvec(i) * params(a1,a2)%shield_const * 331.5_RDbl * &
                params(a1,a2)%charge1 * params(a1,a2)%charge2/(r*r2)
      End If
  
      pot = pot + 331.5 * params(a1,a2)%shield_const * params(a1,a2)%charge1 * &
                params(a1,a2)%charge2/r 
      coulflag= .True. 
    End Do
  End Subroutine coul_calc_multInteraction

  !---------------------------------------------------------------------------
  ! Calculate the interaction between a set of atomic points.  Returns False
  ! if any of the seperation vectors fall within the lower cutoff.  This 
  ! indicates atoms have gotten too close together and that the potential 
  ! would be extremely large.  The returned derivative components are those
  ! for the first (?) atom in the separation vector defined by (atom1 - atom2)
  ! Requires:  params -- Coulombicinteraction parameters
  !            sepvec -- array of separation vectors
  !            ffout -- array of resulting energies and derivatives
  !---------------------------------------------------------------------------
  Logical Function coul_grpint(params,sepvec,results)
    Type(Coul_Params), Dimension(:), Intent(In)     :: params
    Type(VecType), Dimension(:), Intent(In)       :: sepvec
    Type(EnergyPlus), Dimension(:), Intent(InOut) :: results

    Integer                   :: pair
    Real(kind = RDbl)         :: r2i,r6i,r12i,r2,r16i,nrg,r
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

      coul_grpint = .False.

      !** Get the square of the distance
      r2 = vector_getnormsq(sepvec(pair))
!      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!      Write(*,*) Sqrt(r2),r2,params(pair)%hicut2
  
      !** Check the cutoffs
      If (r2 > params(pair)%hicut2) Cycle
      If (r2 < params(pair)%locut2) Return
  
       r=Sqrt(r2)

      !** Calculate the potential energy
      nrg = 331.5 * params(pair)%shield_const * params(pair)%charge1 * params(pair)%charge2/r 
      Call storebase_inc(results(pair),nrg)

      !** Calculate the gradients of the energy  
      If (storebase_nderivs(results(pair)) > 0) Then
        grad = sepvec(pair) * params(pair)%shield_const * 331.5_RDbl * &
                params(pair)%charge1 * params(pair)%charge2/(r*r2)
        Call storebase_inc(results(pair),grad)
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

    coul_grpint = .True. 

  End Function coul_grpint

  !---------------------------------------------------------------------------
  ! attempt to speed by passing in nderivs
  !---------------------------------------------------------------------------
  Logical Function coul_grpint2(params,sepvec,results,nderivs)
    Type(Coul_Params), Dimension(:), Intent(In)     :: params
    Type(VecType), Dimension(:), Intent(In)       :: sepvec
    Type(EnergyPlus), Dimension(:), Intent(InOut) :: results
    Integer, Intent(In)                           :: nderivs

    Integer                   :: pair
    Real(kind = RDbl)         :: r2i,r6i,r12i,r2,r16i,nrg, r
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

      coul_grpint2 = .False.

      !** Get the square of the distance
      r2 = vector_getnormsq(sepvec(pair))
!      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!      Write(*,*) Sqrt(r2),r2,params(pair)%hicut2
  
      !** Check the cutoffs
      If (r2 > params(pair)%hicut2) Cycle
      If (r2 < params(pair)%locut2) Return
  
       r = Sqrt(r2)

      !** Calculate the potential energy
      nrg = 331.5 * params(pair)%shield_const * params(pair)%charge1 * params(pair)%charge2/r 
      Call storebase_inc(results(pair),nrg)

      !** Calculate the gradients of the energy  
      If (nderivs > 0) Then
        grad = sepvec(pair) * params(pair)%shield_const * 331.5_RDbl * &
                params(pair)%charge1 * params(pair)%charge2/(r*r2)
        Call storebase_inc(results(pair),grad)
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

    coul_grpint2 = .True. 

  End Function coul_grpint2

  !---------------------------------------------------------------------------
  ! Calculate the interaction between a set of atomic points.  Returns False
  ! if any of the seperation vectors fall within the lower cutoff.  This 
  ! indicates atoms have gotten too close together and that the potential 
  ! would be extremely large.  The returned derivative components are those
  ! for the first atom in the separation vector defined by (atom1 - atom2).
  ! The optional second results structure must be passed if forces are to be
  ! calculated.
  ! Requires:  params -- Coulolmbic interaction parameters
  !            sepvec -- atom-atom separation vector
  !            results -- resultant energies and derivatives
  !            results2 -- resultant energies and opposite forces
  !---------------------------------------------------------------------------
  Logical Function coul_snglint(params,sepvec,results,results2)
    Type(Coul_Params), Intent(In)               :: params
    Type(VecType), Intent(In)                 :: sepvec
    Type(EnergyPlus), Intent(InOut)           :: results
    Type(EnergyPlus), Intent(InOut), Optional :: results2

    Real(kind = RDbl)         ::  r2i,r6i,r12i,r2,r16i,nrg, r
    Type(VecType)             :: force

    coul_snglint = .False.
    
    !** Get the square of the distance
    r2 = vector_getnormsq(sepvec)

    !** Check the cutoffs
    If (r2 > params%hicut2) Then
      coul_snglint = .True.      
      Return
    End If
    If (r2 < params%locut2) Return
    
      r = Sqrt(r2)

    !** Calculate the potential energy
      nrg=331.5 * params%shield_const * params%charge1 * params%charge2/r 
    Call storebase_inc(results,nrg)
!    If (Present(results2)) Call storebase_inc(results2,nrg)
!    results%total%nrg = results%total%nrg + nrg
    
    !** Calculate the gradients of the energy  
    If (storebase_nderivs(results) > 0) Then
      force = sepvec * params%shield_const * 331.5_RDbl * &
                params%charge1 * params%charge2/(r*r2)
      Call storebase_inc(results,force)
      Call storebase_inc(results2,(force*(-1.0_RDbl)))
!    results%total%grad = results%total%grad + force
    End If
    

    coul_snglint = .True. 

  End Function coul_snglint
  

  !----------------------------------------------------------------------------
  ! Returns the value of the high end cutoff. If the NEIGHBOR flag is 
  ! present, check to see if the NEIGHBOR cutoff should be returned
  ! Neighbour cut-off is larger value than  high cutoff. It is useful while
  ! making list of intra-pairs in the beginning during molecule 
  ! initialization.
  ! Requires: params -- Coulombic params
  !           neighbor -- flag to get neighbor cutoff
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function coul_getcutoff(params,neighbor)
    Type(Coul_Params), Intent(In)         :: params
    Logical, Optional, Intent(In)       :: neighbor

    !** Default value
    coul_getcutoff = params%hicut

    !** Check if we need to return the neighbor cutoff
    If (Present(neighbor)) Then
      If (neighbor) Then
        coul_getcutoff = params%ncut
      End If
    End If
  End Function coul_getcutoff

  !-------------------------------------------------------
  ! Returns a string with the properly formatted Coulombic params
  ! Requires: params -- Coulombic params
  !-------------------------------------------------------
   Function coul_displayparams(coulparams)
     Character(len=strLen)       :: coul_displayparams
     Type(Coul_Params), Intent(In) :: coulparams

     Write(coul_displayparams,'(1x,a,2x,a)') Trim(real2str(coulparams%charge1,5)), &
         Trim(real2str(coulparams%charge2,5))
   End Function coul_displayparams

  !----------------------------------------------------------------------------
  ! Returns a string containing the cutoff Values
  ! Requires: params -- Coul params
  !----------------------------------------------------------------------------
  Function coul_displayCutOffs(coulparams)
    Character(len=lstrLen)          :: coul_displayCutOffs
    Type(Coul_Params), Intent(In)     :: coulparams
    
    Character(len=strLen)           :: string1,string2

    If (coulparams%off) Then
      Write(coul_displayCutoffs,'(a)') "OFF"
    Else
      string1 = real2str(coulparams%hicut,5)
      string2 = real2str(coulparams%locut,5)
      Write(coul_displayCutOffs,'(2(2a,1x,a,3x),f8.4)') &
          "high cutoff: ",Trim(string1),"Ang", &
          "low cutoff: ",Trim(string2),"Ang"
    End If

  End Function coul_displayCutOffs

  !-----------------------------------------------------------------------
  ! Displays the full set of Coul parameters for the structure
  ! Requires: params -- the Coul parameters structure
  !           indent -- number of spaces to indent
  !           unit -- unit to write into
  !-----------------------------------------------------------------------
  Subroutine coul_display(params,indent,unit)
    Type(Coul_Params), Intent(In)         :: params
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
          Trim(params%atom_name2),' Coulombic parameters from ',Trim(params%source)
      string = coul_displayCutOffs(params)
      Write(unit,'(2x,2a)') blank,Trim(string)
    End If

  End Subroutine coul_display

End Module coul















