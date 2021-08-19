!------------------------------------------------------------------------------
! Module contains a data structure to hold the parameters for velocity 
! rescaling
!------------------------------------------------------------------------------
Module velrescale

  Use defaults, Only: RDbl, strLen, lstrLen, scaleke, kjmole_kb, &
      MAX_DOF
  Use utils, Only: stripcmnt, split, toint
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use atom, Only: atom_getmass
  Use molecules, Only: molecules_getnsorbs, molecules_getdof, &
      molecules_getnthatom, molecules_iscoreshell
  Use storestats, Only: storestats_gettemp
  Use interact, Only: Interaction_Model

  Implicit None

  Private
  Public :: VelRescaleParams, velrescale_init, velrescale_display, &
      velrescale_rescale, velrescale_sampleCF

  Type VelRescaleParams
    Integer :: scaleType
    Integer :: scaleStep
    Real(Kind=RDbl), Dimension(:), Pointer :: tsum
    Integer, Dimension(:), Pointer :: vrCalls
  End Type VelRescaleParams

  Character(len=strLen), Parameter :: method1 = "T(mole) scaling"
  Character(len=strLen), Parameter :: method2 = "Avg(T(mole),T(atom)) scaling"
  Character(len=strLen), Parameter :: method3 = &
      "Alternate T(mole) T(atom) scaling"
  Character(len=strLen), Parameter :: method4 = "Scales over all degrees of f&
      &reedom for one sorb(first one)"
  Character(len=strLen), Parameter :: method5 = "COM temperature scaling &
      &for core-shell systems"
  Character(len=strLen), Parameter :: method6 = "T(Atom) scaling, (useful for 1 complex molecule)"

Contains

  !----------------------------------------------------------------------------
  ! Initializes the parameters for velocity rescaling
  !----------------------------------------------------------------------------
  Subroutine velrescale_init(vr,unitno)
    Type(VelRescaleParams), Intent(InOut) :: vr
    Integer, Intent(In) :: unitno
    Integer :: nsorbs, i, error
    Character(len=strLen), Dimension(strLen) :: params
    Character(len=lstrLen) :: text
    
    !** Read the parameters
    !** First, read the rescaling type (1, 2, or 3)
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    i = split(text,params)
    vr%scaleType = toint(params(1))

    !** Read the scaling step
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    i = split(text,params)
    vr%scaleStep = toint(params(1))
    
    nsorbs = molecules_getnsorbs()
    Allocate(vr%tsum(nsorbs),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          "Could not allocate memory for 'vr%tsum'"
      Stop
    End If
    Allocate(vr%vrCalls(nsorbs),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          "Could not allocate memory for 'vr%vrCalls'"
      Stop
    End If
    Do i = 1, nsorbs
      vr%tsum(i) = 0.0_RDbl
      vr%vrCalls(i) = 0
    End Do
    
  End Subroutine velrescale_init
  
  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine velrescale_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a,t30,2a)') '[123]','# Rescale type 1=Molec, 2=Average ', &
        'Atom & Molec, 3=Switch Between Atom & Molecule'
    Write(unitno,'(a,t30,a)') 'Integer','# Steps between rescaling'
  End Subroutine velrescale_sampleCF

  !----------------------------------------------------------------------------
  ! Rescales the velocity for the given system
  !----------------------------------------------------------------------------
  Subroutine velrescale_rescale(imodel,vr,sorb,v,desiredTemp)
    Type(Interaction_Model), Intent(In)          :: imodel
    Type(VelRescaleParams), Intent(InOut)        :: vr
    Type(VecType), Dimension(:,:), Intent(InOut) :: v
    Real(kind=RDbl), Intent(In)                  :: desiredTemp
    Integer, Intent(In)                          :: sorb

    !** Increase the number of calls
    vr%vrCalls(sorb) = vr%vrCalls(sorb) + 1

    !** Update the temp specs
    Select Case (vr%scaleType)
      
    Case (1)
      !** Scale by the molecular temperature
      vr%tsum(sorb) = vr%tsum(sorb) + storestats_gettemp(imodel%spcstats, &
          sorb,'tmole','inst')
      
    Case (2)
      !** Scale by the average of the two temperatures
      vr%tsum(sorb) = vr%tsum(sorb) + 0.5_RDbl* &
          (storestats_gettemp(imodel%spcstats,sorb,'tmole','inst') + &
          storestats_gettemp(imodel%spcstats,sorb,'tatom','inst'))
      
    Case (3)
      !** Toggle between atomic and molecular temperature
      If (Mod(vr%vrCalls(sorb),2) == 0) Then
        vr%tsum(sorb) = vr%tsum(sorb) + storestats_gettemp(imodel%spcstats,&
            sorb,'tmole','inst')
      Else
        vr%tsum(sorb) = vr%tsum(sorb) + storestats_gettemp(imodel%spcstats, &
            sorb,'tatom','inst')
      End If

    Case (4)
      !** Scale by all the degrees of freedom
      Call velrescale_dof(vr,v,sorb,desiredTemp)

    Case (5)
      !** Scale by the molecular temperature
      vr%tsum(sorb) = vr%tsum(sorb) + storestats_gettemp(imodel%spcstats, &
          sorb,'tatom','inst')

    Case (6)
      !** Scale by the average of the two temperatures
      vr%tsum(sorb) = vr%tsum(sorb) + storestats_gettemp(imodel%spcstats,sorb,'tatom','inst')

    Case Default
      !** Scale by the average of the two temperatures
      Write(0,'(1x,2a,i4,a)') __FILE__," : stopping at ",__LINE__, &
          "tag specified for rescale method is not valid"
      Stop
    End Select
    
    !** Check to see if we have to scale or not
    If ((vr%scaleType /= 4).And.(Mod(vr%vrCalls(sorb),vr%scaleStep) == 0)) Then
      !** Write to STDERR to let the user know we are rescaling
      vr%tsum(sorb) = vr%tsum(sorb)/Real(vr%scaleStep)
!!Lev      Write(0,'(2x,a,f9.3,1x,a)') "Rescaling velocities, T(old) : ", &
!!Lev          vr%tsum(sorb),"K"
      If (vr%scaleType == 5) Then
        Call velrescale_scale(v,desiredTemp,vr%tsum(sorb))
        Write(0,'(2x,a,f9.3,1x,a)') "Rescaling velocities, T(new) : ", &
            velrescale_comtemp(v,sorb),"K"
      Else
        Call velrescale_scale(v,desiredTemp,vr%tsum(sorb))
      End If
      vr%tsum(sorb) = 0.0_RDbl
    End If
    
  End Subroutine velrescale_rescale

  !----------------------------------------------------------------------------
  ! Calculates the com temperature
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function velrescale_comtemp(vel,mType)
    Type(VecType), Dimension(:,:), Intent(In) :: vel
    Integer, Intent(In) :: mType
    Integer :: a,m, atype,partner
    Real(Kind=RDbl) :: mass, pMass, ke
    Logical :: isCore
    Type(VecType) :: v

    !** Zero the kinetic energy sum
    ke = 0.0_RDbl

    !** Loop over all the atoms and molecules
    Do m = 1, Size(vel,2)
      Do a = 1, Size(vel,1)
        !** Get the atom type and mass
        atype = molecules_getnthatom(mType,a)
        mass = atom_getmass(atype)
        !** Check to see if this is a core-shell pair. If so, we
        !** must add the COM velocity, not the core and shell 
        !** velocities individually.
        If (molecules_isCoreShell(mType,a,isCore,partner)) Then
          If (.Not.isCore) Cycle
          !** We need the COM velocity
          pMass = atom_getmass(molecules_getnthatom(mType,partner))
          v = (vel(a,m)*mass + vel(partner,m)*pMass) / (mass + pMass)
          mass = mass + pMass
        Else
          v = vel(a,m)
        End If
        !** Sum the kinetic energy
        ke = v*v*mass + ke
      End Do
    End Do

    velrescale_comtemp = ke*scaleke / (kjmole_kb* &
        Real(molecules_getdof(mType)*Size(vel,2)))

  End Function velrescale_comtemp

  !----------------------------------------------------------------------------
  ! Scales core-shell pairs as a single COM atom, and scales all other atoms
  ! normally.
  !----------------------------------------------------------------------------
  Subroutine velrescale_comscale(v,mType,desiredT,currentT)
    Type(VecType), Dimension(:,:), Intent(InOut) :: v
    Integer, Intent(In) :: mType
    Real(Kind=RDbl) :: desiredT, currentT

    !** There can be only 1 neighbor for a core type atom.
    Integer :: neighbor
    Integer :: i, j, aType, asType
    Logical :: isCore
    Real(Kind=RDbl) :: aMass, asMass, scale
    Type(VecType) :: comv

    scale = Sqrt(desiredT/currentT)

    Do i = 1, Size(v,1)

      !** Get the atom type
      aType = molecules_getnthatom(mType,i)

      !** Check to see if this is a shell. If so, cycle, we'll pick it 
      !** up when we do the core.
      If (molecules_isCoreShell(mType,aType,isCore,neighbor)) Then
        !** If this is not the core, cycle. 
        If (.Not.isCore) Cycle

        !** Get the shell atom type
        asType = molecules_getnthatom(mType,neighbor)

        !** Get the atom masses
        aMass = atom_getmass(aType)
        asMass= atom_getmass(asType)

        Do j = 1, Size(v,2)
          !** Calculate the center of mass velocity
          comv = (v(i,j)*aMass + v(neighbor,j)*asMass) / (aMass + asMass)
          !** Scale the velocities. They must be equal to keep the 
          !** relative velocity at zero.
          v(i,j) = comv * scale
          v(neighbor,j) = comv * scale
        End Do

      Else
        Do j = 1, Size(v,2)
          v(i,j) = v(i,j)*scale
        End Do

      End If
    End Do
      
  End Subroutine velrescale_comscale

  !----------------------------------------------------------------------------
  !scales each dof
  !----------------------------------------------------------------------------
  Subroutine velrescale_dof(vr,v,sorb,T)
    Type(VelRescaleParams), Intent(InOut) :: vr
    Type(VecType), Dimension(:,:), Intent(InOut) :: v
    Integer, Intent(In)                          :: sorb
    Real(kind=RDbl), Intent(In) :: T
    Real(kind=RDbl) :: scale,fac,mass
    Real(kind=RDbl),Dimension(MAX_DOF,3) :: Tsc
    Integer :: i,j,k,natoms,nmoles,atomt

    If (Mod(vr%vrCalls(sorb),vr%scaleStep) == 0) Then
      natoms=Size(v,1)
      nmoles=Size(v,2)
      fac=scaleke/kjmole_kb/nmoles
      If ((natoms*3)>MAX_DOF) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "Too many degrees of freedoms, for this molecule&
            &Cannot be rescaled properly"
      Endif
      
      Do i=1,natoms
        Tsc(i,1)=0.0_RDbl
        Tsc(i,2)=0.0_RDbl
        Tsc(i,3)=0.0_RDbl
        Do j=1,nmoles
          Do k=1,3
            
            atomt = molecules_getnthatom(sorb,i)
            mass = atom_getmass(atomt)
            Tsc(i,k)=(v(i,j)%comp(k)*v(i,j)%comp(k))*mass*fac+Tsc(i,k)
          End Do
        End Do
      End Do
      
      !** Write to STDERR to let the user know we are rescaling
      Do i=1,natoms
        Do j=1,nmoles
          Do k=1,3
            scale=Sqrt(T/Tsc(i,k))
            v(i,j)%comp(k)=v(i,j)%comp(k)*scale
          End Do
        End Do
      End Do
      
      Write(*,*) "Rescaling vels by all dofs of sorb1,  &
          &T(old) : ", Tsc(1:natoms,1:3),"K"
    End If
  End Subroutine velrescale_dof
  
  !----------------------------------------------------------------------------
  ! Performs the rescaling
  !----------------------------------------------------------------------------
  Subroutine velrescale_scale(v,desiredTemp,currentTemp)
    Type(VecType), Dimension(:,:), Intent(InOut) :: v
    Real(kind=RDbl), Intent(In) :: desiredTemp, currentTemp
    Real(kind=RDbl) :: scale
    Integer :: i,j
    
    scale  = Sqrt(desiredTemp/currentTemp)

!!$    !MDEBUG
!!$    Write(0,*) __FILE__,__LINE__,scale
    
    Do i = 1, Size(v,1)
      Do j = 1, Size(v,2)
        v(i,j) = v(i,j)*scale
      End Do
    End Do
    
  End Subroutine velrescale_scale
  
  
  !----------------------------------------------------------------------------
  ! Displays information about the velocity rescaling parameters
  !----------------------------------------------------------------------------
  Subroutine velrescale_display(vr,unitno,skip)
    Type(VelRescaleParams), Intent(In) :: vr
    Integer, Intent(In) :: skip
    Integer, Intent(In) :: unitno
    Character(len=skip) :: blank
    Character(len=strLen) :: text
    Integer :: i 
    
    blank = ""
    Do i = 1, skip
      blank = blank//" "
    End Do
    
    Select Case (vr%scaleType)
    Case (1)
      text = method1
    Case (2)
      text = method2
    Case (3)
      text = method3
    Case (4)
      text = method4
    Case (5)
      text = method5
    Case (6)
      text = method6
    End Select
    
    Write(unitno,'(2a)') blank,"Velocity Rescale Parameters"
    Write(unitno,'(2a,i10)') blank,"Steps between rescaling : ",vr%scaleStep
    Write(unitno,'(3a)') blank,"Rescaling method        : ",text
    
  End Subroutine velrescale_display
  
End Module velrescale



