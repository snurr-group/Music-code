!-----------------------------------
! Evans and Morriss constraints
!-----------------------------------
Module evansmorriss

  Use defaults, Only: RDbl, strLen, scalepe, pcalls
  Use general, Only: genparams
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use utils, Only: stripcmnt, split, toint, toreal
  Use matrixops, Only: matrixops_ludcmp, matrixops_lubksb

  Implicit None

  Private
  Public  :: EvansMorrissModel, ConConnect, evansmorriss_getnconstr,  &
      evansmorriss_init, evansmorriss_getconstraints, &
      evansmorriss_getpenalty, evansmorriss_display, evansmorriss_cleanup

  Type ConConnect
    Integer :: atom1
    Integer :: atom2
    Real(kind=RDbl) :: length
  End Type ConConnect

  Type EvansMorrissModel
    Real(kind=RDbl), Dimension(:,:), Pointer :: M
    Real(kind=RDbl), Dimension(:,:), Pointer :: L
    Integer :: nconstr
    Type(ConConnect), Dimension(:), Pointer :: conlist
  End Type EvansMorrissModel

  Type(EvansMorrissModel), Pointer, Save     :: emParamsPntr
  Type(VecType), Dimension(:), Pointer, Save :: emCoordsPntr

Contains

  !----------------------------------------------------------------------------
  ! Returns the number of constraints
  !----------------------------------------------------------------------------
  Integer Function evansmorriss_getnconstr(empntr)
    Type(EvansmorrissModel),Pointer :: empntr
    evansmorriss_getnconstr = empntr%nconstr
  End Function evansmorriss_getnconstr

  !----------------------------------------------------------------------------
  ! Initialize the M and L matrices for the evansmorriss model
  ! Requires: empntr -- Type(EvansmorrissModel) 
  !           natoms -- number of constraints
  !           ncon -- number of constraints
  !           invmass -- a list of inverse masses for each atom
  !           atom1list -- list of the first atoms in the constraints
  !           atom2list -- list of the second atoms in the constraints
  !           lengthlist -- list of bond lengths
  !----------------------------------------------------------------------------
  Subroutine evansmorriss_init(empntr,natoms,ncon,invmass,atom1list, &
      atom2list,lengthlist)

    Type(EvansmorrissModel),Pointer            :: empntr
    Integer, Intent(In)                        :: natoms,ncon
    Real(kind=RDbl), Dimension(:), Intent(In)  :: invmass
    Integer, Dimension(:), Intent(In)          :: atom1list,atom2list
    Real(Kind=RDbl), Dimension(:), Intent(In)  :: lengthlist

    Integer                      :: error,i

    nullify(emCoordsPntr)
    nullify(emParamsPntr)

    empntr%nconstr = ncon

    Allocate(empntr%conlist(ncon),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
          " : Could not allocate memory for 'conlist'"
      Stop
    End If

    Do i = 1, empntr%nconstr
      empntr%conlist(i)%atom1 = atom1list(i)
      empntr%conlist(i)%atom2 = atom2list(i)
      empntr%conlist(i)%length = lengthlist(i)
    End Do
      
    !** Allocate the M matrix
    Allocate(empntr%M(natoms,ncon),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
          " : Could not allocate memory for 'M' matrix"
      Stop
    End If
    
    !** Allocate the L matrix
    Allocate(empntr%L(ncon,ncon),stat=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
          " : Could not allocate memory for 'L' matrix"
      Stop
    End If

    !** Initialize the L and M matrices
    Call evansmorriss_initparams(empntr,invmass)

  End Subroutine evansmorriss_init

  !----------------------------------------------------------------------------
  ! Initialize the M and L matrices' values for the model
  ! Requires: empntr -- Type(EvansmorrissModel) 
  !           invmass -- a list of inverse masses for each atom
  !----------------------------------------------------------------------------
  Subroutine evansmorriss_initparams(empntr,invmass)
    Type(EvansmorrissModel),Pointer            :: empntr
    Real(kind=RDbl), Dimension(:), Intent(In)  :: invmass

    Integer :: i,j, nconstr, atom1, atom2

    nconstr = empntr%nconstr

    !** Initialize the M matrix
    empntr%M = 0.0_RDbl
    Do i = 1,nconstr
      atom1 = empntr%conlist(i)%atom1
      atom2 = empntr%conlist(i)%atom2 
      empntr%M(atom1,i) = -invmass(atom1)
      empntr%M(atom2,i) = invmass(atom2) 
    End Do
    
    !** Initialize the L matrix
    empntr%L = 0.0_RDbl
    Do i = 1,nconstr
      atom1 = empntr%conlist(i)%atom1
      atom2 = empntr%conlist(i)%atom2
      Do j = 1,nconstr
        empntr%L(i,j) = empntr%M(atom2,j)-empntr%M(atom1,j)
      End Do
    End Do

  End Subroutine evansmorriss_initparams

  !----------------------------------------------------------------------------
  ! Displays the parameters 
  !----------------------------------------------------------------------------
  Subroutine evansmorriss_display(empntr,unit,indent)
    Type(EvansmorrissModel), Pointer :: empntr
    Integer, Intent(In) :: unit, indent
    Integer :: i, natoms, nconstr
    Character(len=indent) :: blank
    Character(len=strLen) :: display,nmbr

    blank = ""
    Do i = 1,indent
      blank = blank//" "
    End Do
    
    nconstr = empntr%nconstr

    Do i = 1,nconstr
      Write(unit,'(2a,2i4,f10.3)') blank,"Atom pair and constr distance : ",& 
          empntr%conlist(i)%atom1, empntr%conlist(i)%atom2, &
          empntr%conlist(i)%length
    End Do

    If (Trim(genparams%displaymode) == "VERBOSE") Then
      !** Display the M matrix
      Write(unit,'(2a)') blank,"Matrix M elements: "
      natoms = Size(empntr%M,1)
      Write(nmbr,'(i6)') natoms
!      display = "(4x,a,"//Trim(nmbr)//"e15.5)"
      display = "(4x,a,"//Trim(nmbr)//"f6.3)"
      Do i = 1, nconstr
        Write(unit,display) blank,empntr%M(1:natoms,i)
      End Do
  
      !** Display the L matrix
      Write(unit,'(2a)') blank,"Matrix L elements: "
      Write(nmbr,'(i8)') nconstr
      display = "(4x,a,"//trim(nmbr)//"e15.5)"
      Do i = 1, nconstr
        Write(unit,display) blank,empntr%L(1:nconstr,i)
      End Do
    End If

  End Subroutine evansmorriss_display

  !----------------------------------------------------------------------------
  ! Calculate the constrains vector and forces
  !----------------------------------------------------------------------------
  Subroutine evansmorriss_getconstraints(empntr,coords,v,u,f)
    Type(EvansmorrissModel),Pointer :: empntr
    Type(VecType), Dimension(:), Intent(In) :: coords,v
    Type(VecType), Dimension(:), Intent(InOut) :: f
    Real(Kind=RDbl), Intent(Out) :: u
    Integer :: nconstr,atom1,atom2,i,j,k
    Type(VecType), Dimension(empntr%nconstr) :: rc,vc,fc
    Real(kind=RDbl), Dimension(Size(empntr%conlist),Size(empntr%conlist)) :: A
!!$    Real(kind=RDbl), Dimension(Size(empntr%conlist)) :: L,d
    Real(kind=RDbl), Dimension(empntr%nconstr) :: L,d
    Real(kind=RDbl) :: check

    u = 0.0_RDbl
    nconstr = empntr%nconstr

    Do i = 1,nconstr
      atom1 = empntr%conlist(i)%atom1
      atom2 = empntr%conlist(i)%atom2
      rc(i) = coords(atom2)-coords(atom1)
      vc(i) = v(atom2)-v(atom1)
      fc(i) = f(atom2)-f(atom1)
      !** rhs of the constraint equation
      d(i) = (rc(i)*fc(i)+vc(i)*vc(i))*(-1.0_RDbl)
    End Do

    !** lhs of the constraint equation
    Do j = 1,nconstr
      Do k = 1,nconstr
        A(j,k) = (rc(j)*rc(k))*empntr%L(j,k)
      End Do
    End Do

    !** solve for the constraint vector
    Call matrixops_ludcmp(A,nconstr,nconstr,L,check)
    Call matrixops_lubksb(A,nconstr,nconstr,L,d)

    !** calculate the constraint forces
    Do j = 1,size(f)
      Do k = 1,nconstr
        f(j) = f(j) + rc(k)*empntr%M(j,k)*d(k)
      End Do
    End Do

    !** Calculate the constraint energy
    u = 0.0_RDbl
    Do j = 1,nconstr
      u = u - 0.5_RDbl*d(j)*(rc(j)*rc(j)- &
          empntr%conlist(j)%length*empntr%conlist(j)%length)
!!$      !MDEBUG
!!$      Write(*,*) "u for nconstr ",j,0.5_RDbl*d(j)*(rc(j)*rc(j)- &
!!$          empntr%conlist(j)%length*empntr%conlist(j)%length)
!!$      Write(*,*) "d(j) ",d(j)
!!$      Write(*,*) "rc(j)*rc(j) ",rc(j)*rc(j)
!!$      Write(*,*) "length ",empntr%conlist(j)%length
    End Do

    u = u/scalepe

  End Subroutine evansmorriss_getconstraints

  !---------------------------------------------------------------------------
  ! This is to work around problems with the objective function and gqbfs
  !---------------------------------------------------------------------------
  Subroutine evansmorriss_penaltyinit(empntr,coords)
    Type(EvansmorrissModel), Intent(In), Target :: empntr
    Type(VecType), Dimension(:), Intent(In), Target :: coords
    emParamsPntr => empntr
    emCoordsPntr => coords
  End Subroutine evansmorriss_penaltyinit

  !----------------------------------------------------------------------------
  ! This is the penalty function used for the constraints
  !----------------------------------------------------------------------------
  Subroutine evansmorriss_getpenalty(empntr,coords,v,sumobjf)
    Type(EvansmorrissModel),Pointer :: empntr
    Type(VecType), Dimension(:), Intent(InOut) :: coords
    Type(VecType), Dimension(:), Intent(InOut) :: v
    Real(Kind=RDbl), Intent(Out) :: sumobjf
!    Integer, Intent(InOut) :: pcalls
    Integer :: natoms, nvar, i, error, j
    Real(kind=Rdbl) :: objfr, objfv
    ! Variables for minimization function
    Real(Kind=RDbl), Parameter :: tolx = 1.0e-6, tolg=1.0e-8, sprec=0.1, &
        extbnd=2.0, rfn=1e-6, tolr=1.0e-9, tolv=1.0e-8
    Integer, Parameter :: irhess=0, iwhess=0, iprint=0, iswtch=2, &
        maxf = 500
    Integer :: ihess, iter, nfcall, ierr
    Real(kind=RDbl) :: objf, histry
!    Real(Kind=RDbl), &
!    Dimension(Int(0.5*Real(Size(coords,1)*3*(3*Size(coords,1)+1)))) :: hess
    Real(Kind=RDbl), &
        Dimension ((Size(coords,1)*3*(3*Size(coords,1)+1))/2 ) :: hess
     Real(kind=RDbl), Dimension(Size(coords,1)*3) :: x, g, solvec, srhvec, &
        grdvec, scrvec

    natoms = Size(coords,1)
    nvar = 3*natoms
    sumobjf = 0.0_RDbl
    ihess = Int(0.5*Real(nvar)*(Real(nvar+1)))
    hess = 0.0_RDbl
    solvec = 0.0_RDbl
    scrvec = 0.0_RDbl
    grdvec = 0.0_RDbl
    srhvec = 0.0_RDbl
    g = 0.0_RDbl
    x = 0.0_RDbl

!!$    !MDEBUG
!!$    Write(*,*) "Size hess ",Size(hess,1)
!!$    Write(*,*) "ihess ",ihess
!!$    Write(*,*) "irhess, iwhess ",irhess, iwhess

    Call evansmorriss_penaltyinit(empntr,coords)

    !** Put the positions into vectors for gqbfgs
    Do i = 1, natoms
      j = 3*(i-1)+1
      x(j:j+2) = coords(i)
    End Do

    Call evansmorriss_phi(nvar,x,objfr,g)

    If (objfr > tolr) Then
      pcalls = pcalls + 1

      error = 0

      ierr = 0

      !** Call the gibberish function
      Call gqbfgs(nvar,x,iswtch,maxf,iter,iprint,tolx,tolg,rfn, sprec, &
          extbnd, objf, ihess, hess, g, ierr, nfcall, srhvec, solvec, &
!!$          grdvec, scrvec, histry, irhess, iwhess, 1)
          grdvec, scrvec, histry, irhess, iwhess, evansmorriss_phi)

      sumobjf = sumobjf + objf

      !** Update the coordinates
      Do i = 1,natoms
        j = 3*(i-1)+1
        coords(i) = solvec(j:j+2)
      End Do
    End If

    !** Do the velocities
    Do i = 1,natoms
      j = 3*(i-1)+1
      x(j:j+2) = v(i)
    End Do

!!$    !MDEBUG
!!$    Write(*,*) "x : ",x
    
    Call evansmorriss_psi(nvar,x,objfv,g)

!!$    !MDEBUG
!!$    Write(*,*) "objfv ",objfv
!!$    Write(*,*) "g ",g

    If (objfv > tolv) Then
      
      pcalls = pcalls+1

      error = 0

      !** Call the gibberish function
      Call gqbfgs(nvar,x,iswtch,maxf,iter,iprint,tolx,tolg,rfn, sprec, &
          extbnd, objf, ihess, hess, g, ierr, nfcall, srhvec, solvec, &
!!$          grdvec, scrvec, histry, irhess, iwhess, 2)
          grdvec, scrvec, histry, irhess, iwhess, evansmorriss_psi)

      sumobjf = sumobjf + objf
      
      Do i = 1, natoms
        j = 3*(i-1)+1
        v(i) = x(j:j+2)
      End Do
    End If

  End Subroutine evansmorriss_getpenalty


  !----------------------------------------------------------------------------
  ! Clean up the evans & morris matrices
  !----------------------------------------------------------------------------
  Subroutine evansmorriss_cleanup(empntr)
    Type(EvansmorrissModel),Pointer :: empntr
    Integer :: error
    If (.Not.Associated(empntr)) Then
      Write(*,*) "Trying to deallocate, non-allocated pointer"
      return
      Endif
    If (Associated(empntr%M)) Then
      Deallocate(empntr%M,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'M' matrix"
        Stop
      End If
    End If
    If (Associated(empntr%L)) Then
      Deallocate(empntr%L,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'L' matrix"
        Stop
      End If
    End If
    If (Associated(empntr%conlist)) Then
      Deallocate(empntr%conlist,stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__, " : ", __LINE__, &
            " : Could not deallocate memory for 'conlist' array"
        Stop
      End If
    End If
  End Subroutine evansmorriss_cleanup

  !----------------------------------------------------------------------------
  ! Calculates phi for the penalty function
  !----------------------------------------------------------------------------
  Subroutine evansmorriss_phi(nvar,x,objfr,g)
    
    Integer, Intent(In) :: nvar
    Real(Kind=RDbl), Dimension(nvar), Intent(In) :: x
    Real(Kind=RDbl), Dimension(nvar), Intent(Out):: g
    Real(Kind=RDbl), Intent(out) :: objfr

    Integer :: i, atom1, atom2, atomi1, atomi2
    Type(VecType) :: pos1, pos2, b
    Type(VecType), Dimension(nvar) :: grad
    Real(kind=RDbl) :: offset, length

!!$    !MDEBUG
!!$    Write(*,*) "In phi calculating"

    g = 0.0_RDbl
    grad = VecType(0.0_RDbl)
    objfr = 0.0_RDbl

    Do i = 1, emParamsPntr%nconstr

      atom1 = emParamsPntr%conlist(i)%atom1
      atom2 = emParamsPntr%conlist(i)%atom2
      atomi1 = (atom1-1)*3+1
      atomi2 = (atom2-1)*3+1

      pos1 = x(atomi1:atomi1+2)
      pos2 = x(atomi2:atomi2+2)

      b = pos2 - pos1
      length = emParamsPntr%conlist(i)%length

      offset = (b*b - length*length)

      objfr = objfr + offset*offset

      grad(atom1) = grad(atom1) - b*offset*4.0_RDbl
      grad(atom2) = grad(atom2) + b*offset*4.0_RDbl

    End Do
    
    !** Copy the gradients into the g arrary
    Do i = 1, emParamsPntr%nconstr
      atom1 = emParamsPntr%conlist(i)%atom1
      atom2 = emParamsPntr%conlist(i)%atom2
      atomi1 = (atom1-1)*3+1
      atomi2 = (atom2-1)*3+1

      g(atomi1:atomi1+2) = grad(atom1)
      g(atomi2:atomi2+2) = grad(atom2)
    End Do

  End Subroutine evansmorriss_phi

  !----------------------------------------------------------------------------
  ! Calculates psi for the penalty function
  !----------------------------------------------------------------------------
  Subroutine evansmorriss_psi(nvar,v,objfv,g)

    Integer, Intent(In) :: nvar
    Real(Kind=RDbl), Dimension(nvar), Intent(In) :: v
    Real(Kind=RDbl), Dimension(nvar), Intent(Out):: g
    Real(Kind=RDbl), Intent(Out) :: objfv

    Type(VecType), Dimension(nvar) :: grad
    Integer :: i, atom1, atom2, atomi1, atomi2
    Type(VecType) :: pos1, pos2,vel1,vel2,b
    Real(kind=RDbl) :: offset

    objfv = 0.0_RDbl
    grad = VecType(0.0_RDbl)
    g = 0.0_RDbl

    Do i = 1,emParamsPntr%nconstr

      atom1 = emParamsPntr%conlist(i)%atom1
      atom2 = emParamsPntr%conlist(i)%atom2

      pos1 = emCoordsPntr(atom1)
      pos2 = emCoordsPntr(atom2)
      b = pos2 - pos1

      atomi1 = (atom1-1)*3+1
      atomi2 = (atom2-1)*3+1

      vel1 = v(atomi1:atomi1+2)
      vel2 = v(atomi2:atomi2+2)

      offset = (pos2-pos1)*(vel2-vel1)

      objfv = objfv + offset*offset

      grad(atom1) = grad(atom1) - b*2.0_RDbl*offset
      grad(atom2) = grad(atom2) + b*2.0_RDbl*offset

    End Do

    !** Copy the gradients into the g arrary
    Do i = 1, emParamsPntr%nconstr
      atom1 = emParamsPntr%conlist(i)%atom1
      atom2 = emParamsPntr%conlist(i)%atom2
      atomi1 = (atom1-1)*3+1
      atomi2 = (atom2-1)*3+1

      g(atomi1:atomi1+2) = grad(atom1)
      g(atomi2:atomi2+2) = grad(atom2)
    End Do

  End Subroutine evansmorriss_psi


  !----------------------------------------------------------------------------
  ! Returns True if the data structure has been initialized
  !----------------------------------------------------------------------------
  Logical Function evansmorriss_isinit(empntr)
    Type(EvansmorrissModel), Intent(In)          :: empntr

    !** Default value
    evansmorriss_isinit = .False.

    If ((Associated(empntr%M)).And.(Associated(empntr%L)).And. &
        (Associated(empntr%conlist))) evansmorriss_isinit = .True.

  End Function evansmorriss_isinit

End Module evansmorriss



