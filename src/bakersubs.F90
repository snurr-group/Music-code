Module bakersubs

  Use defaults, Only: RDbl, strLen, MAX_DOF, MAX_ATOMS, NO_OF_INTRA_POTS, &
      scalepe
  Use config, Only: config_zeroforces, config_conv2force, config_getnatoms, &
      config_isfixed, AtMolCoords
  Use simcell, Only: SimCell_Params, simcell_pbc
  Use gcmodels, Only: gcmodels_getgencoords, gcmodels_getgenforce, &
      gcmodels_getmetrictensor, gcmodels_setgencoords, gcmodels_toxyz
  Use forcefield, Only: forcefield_getssint
  Use molecules, Only: molecules_getdof, &
      molecules_getnatoms, molecules_updateEnergySS, &
      molecules_updateIntraEnergy, molecules_getnsorbs
  Use vector, Only: VecType

  Implicit None
  Save

  Private
  Public :: bakersubs_getenergy, bakersubs_sorbtox, bakersubs_calcgrad, &
      bakersubs_xtosorb, bakersubs_calchess, bakersubs_calcMetricTensor

Contains  
  !-----------------------------------------------------------------------
  ! Do the energy calculation.  Returns if the energy calcn. was
  ! successful.  Namely, if the atoms are not in an untabulated region
  ! of the zeolite.
  !-----------------------------------------------------------------------
  Logical Function bakersubs_getenergy(sorbates, simcell, mtypes, pot)
    Type(AtMolCoords), Dimension(:), Intent(InOut) :: sorbates
    Type(SimCell_Params), Intent(inout) :: simcell
    Integer, Dimension(:), Intent(in)   :: mtypes
    Integer :: nsorbs
    Real(kind=RDbl), Intent(inout)      :: pot
    Real(kind=RDbl), Dimension(Size(sorbates,1),Size(sorbates,1)) :: &
        ncoulnrg, coulnrg
    Real(Kind=RDbl), Dimension(Size(sorbates,1),NO_OF_INTRA_POTS) :: intranrg
    Integer      :: i, j, mtype
    Logical      :: fast, pcalcflag, mapflag

    !** Number of baker molecules
    nsorbs = Size(mtypes,1)

    !** Do the energy calculation
    fast = .False.
    pcalcflag = .False.
    ! Zero the forces
    Do i=1, nsorbs
      mtype = mtypes(i)
      Call config_zeroforces(sorbates(mtype), fast)
    End Do
    ! Calculate the pairwise interactions.  The breakup of the energy
    ! is stored in the molecules structure
    Call forcefield_getssint(sorbates, simcell, fast, pot, mapflag, &
        ncoulnrg, coulnrg, intranrg, pcalcflag)

    !** Update the energies
    Do i = 1, molecules_getnsorbs()
      Do j = 1, molecules_getnsorbs()
        Call molecules_updateEnergySS(i,j,'coul',coulnrg(i,j)*scalepe)
        Call molecules_updateEnergySS(i,j,'noncoul',ncoulnrg(i,j)*scalepe)
      End Do
      Call molecules_updateIntraEnergy(i,intranrg(i,:)*scalepe)
    End Do

    ! Convert the acceleration to forces
    Do i=1, nsorbs
      mtype = mtypes(i)
      Call config_conv2force(sorbates(mtype), mtype, fast)
    End Do
    
    !** Return true if the energy calcn. was successful
    bakersubs_getenergy = mapflag
  End Function bakersubs_getenergy

  !--------------------------------------------------------------------------
  ! Calculates the gradient "grad" of the total potential with respect
  ! to the variables in "x". It returns the total potential "pot" 
  ! for the config stored in the array "x". "ierr" = 1 if the energy calcn
  ! fails, 0 otherwise
  !--------------------------------------------------------------------------
  Subroutine bakersubs_calcGrad(sorbates, simcell, mtypes, x, pot, grad, ierr)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(SimCell_Params), Intent(inout)         :: simcell
    Integer, Intent(in), Dimension(:)           :: mtypes
    Real(kind=RDbl), Dimension(:), Intent(in)   :: x
    Real(kind=RDbl), Intent(out)                 :: pot
    Real(kind=RDbl), Dimension(:), Intent(inout):: grad
    Integer, Intent(out)                        :: ierr

    Integer      :: sorb1, sorb2, i, natoms, begidx, endidx, dof, nsorbs, mtype
    Logical      :: fast
    Real(kind=RDbl), Dimension(6)   :: genForce !6 degrees of freedom
    Logical      :: mapflag

    !** Number of baker sorbates
    nsorbs = Size(mtypes,1)

    ierr = 0

    !** Go from the "x" representation to the "sorbates" representation
    Call bakersubs_xToSorb(sorbates, simcell, x, mtypes)

    !** Do the energy calculation
    pot = 0.0_RDbl
    mapflag = bakersubs_getenergy(sorbates, simcell, mtypes, pot)
    If (.Not. mapflag) Then
      ierr = 1
      Return
    End If

    !** Go from the forces in cartesian coordinates on each atom to
    !** the forces on the generalized coordinates
    fast = .False.
    begidx = 1
    Do i=1, nsorbs
      mtype = mtypes(i)
      natoms = config_getnatoms(sorbates, mtype)
      dof = molecules_getdof(mtype)
      endidx = begidx + dof - 1
      If (fast) Then
        Call gcmodels_getGenForce(sorbates(mtype)%gcoords(1), &
            sorbates(mtype)%afast(1:natoms, 1), genForce)
      Else
        Call gcmodels_getGenForce(sorbates(mtype)%gcoords(1), &
            sorbates(mtype)%aslow(1:natoms, 1), genForce)
      End If
      ! The gradient is the negative of the force
      grad(begidx:endidx) = -genForce(1:dof)
      begidx = endidx + 1
    End Do
  End Subroutine bakersubs_calcgrad

  !----------------------------------------------------------------------
  !
  ! HESSFD: This routine calculates a numerical approximation to the
  ! Hessian by finite differences of the analytical gradients. 
  !  INPUT
  !     df...............Vector of degrees of freedom
  !  OUTPUT
  !     hess.............The Hessian matrix
  !
  !---------------------------------------------------------------------
  Subroutine bakersubs_calcHess(sorbates, simcell, mtypes, df, hess, ierr)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(inout)            :: simcell
    Integer, Intent(in), Dimension(:)              :: mtypes
    Real(kind=RDbl), Dimension(:), Intent(in)   :: df
    Real(kind=RDbl), Dimension(:,:), Intent(out):: hess
    Integer, Intent(out)                        :: ierr

    Integer                          :: idf, jdf, ndf
    Real(kind=RDbl)                  :: ave, ener, h
    Real(kind=RDbl), Dimension(Size(df,1)) :: dfact, gra, grb, grc, grd

    !===================================
    !Set the finite difference step size
    !===================================
    h = 1.0d-05
    ndf = Size(df,1)
    
    !=====================================
    ! Initialize the "active" coordinates
    !=====================================
    dfact = df

    !=========
    !Main loop
    !=========
    Do jdf = 1, ndf
      dfact(jdf) = df(jdf) - 2.0_RDbl*h
      Call bakersubs_calcGrad(sorbates, simcell,mtypes,dfact,ener,gra,ierr)

      dfact(jdf) = df(jdf) - h
      Call bakersubs_calcGrad(sorbates, simcell,mtypes,dfact,ener,grb,ierr)
                                                                  
      dfact(jdf) = df(jdf) + h                                    
      Call bakersubs_calcGrad(sorbates, simcell,mtypes,dfact,ener,grc,ierr)
                                                                  
      dfact(jdf) = df(jdf) + 2.0_RDbl*h                           
      Call bakersubs_calcGrad(sorbates, simcell,mtypes,dfact,ener,grd,ierr)

      !** This rather strange finite differencing scheme is taken from
      !** Lapidus & Pinder, pg. 41, 
      !** "Numerical Solution of Partial Diff. Eq. in Science and Engg."
      !** This is accurate to O(h^3)
      Do idf = 1, ndf
        hess(idf,jdf) = (gra(idf) - 8.0_RDbl*grb(idf) &
            + 8.0_RDbl*grc(idf) - grd(idf))/(12.0_Rdbl*h)
      End Do
    End Do

    !======================
    !Symmetrize the hessian      
    !======================
    ! Due to the error in the finite difference approximation, the Hessian
    ! as calculated above will not necessarily be exactly symmetric.  This
    ! causes problems in the subsequent diagonalization. The code below
    ! ensures that the Hessian that is returned is symmetric. 
    Do idf = 1, ndf
      Do jdf = (idf+1), ndf
        ave = 0.5d0*(hess(idf,jdf) + hess(jdf,idf))
        hess(idf,jdf) = ave
        hess(jdf,idf) = ave
      End Do
    End Do
    Return
  End Subroutine bakersubs_calcHess

  
  !----------------------------------------------------------------------
  ! This is the metric tensor as defined in Randy's paper on Benzene TST
  !----------------------------------------------------------------------
  Subroutine bakersubs_calcMetricTensor(sorbates, simcell, mtypes, df, G, ierr)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(in)               :: simcell
    Integer, Intent(In), Dimension(:)              :: mtypes
    Real(kind=RDbl), Dimension(:), Intent(in)      :: df
    Real(kind=RDbl), Dimension(:,:), Intent(out)   :: G
    Integer, Intent(out)   :: ierr

    Integer      :: startdof, enddof, dof, nsorbs, mtype, i

    !** Number of baker sorbs
    nsorbs = Size(mtypes,1)

    !** Go to the sorbate representation
    Call bakersubs_xToSorb(sorbates, simcell, df, mtypes)
    
    ierr = 0
    G = 0.0_RDbl
    startdof = 1
    Do i = 1, nsorbs
      mtype = mtypes(i)
      dof = molecules_getdof(mtype)
      enddof   = startdof + dof - 1
      Call gcmodels_getMetricTensor( sorbates(mtype)%gcoords(1), &
          mtype, G(startdof:enddof, startdof:enddof), ierr)
      If (ierr /= 0) Then
        Return
      End If
      startdof = enddof + 1
    End Do
  End Subroutine bakersubs_calcMetricTensor

  !------------------------------------------------------------------------
  ! Change from the "sorbates" representation to the "x"
  ! representation. It takes the generalized coordinates 
  ! from the "sorbates" structure and puts them in the "x" array.
  !------------------------------------------------------------------------
  Subroutine bakersubs_sorbToX(sorbates, x, mtypes)
    Type(AtMolCoords), Dimension(:), Intent(in) :: sorbates
    Real(kind=RDbl), Dimension(:), Intent(out)  :: x
    Integer, Dimension(:), Intent(in)           :: mtypes

    Integer  :: nsorbs, mtype
    Integer  :: nvar, i, dof, begidx, endidx
    Real(kind=RDbl), Dimension(MAX_DOF) :: gencoords

    !** number of baker sorbs
    nsorbs = Size(mtypes,1)

    begidx = 1
    Do i = 1, nsorbs
      mtype = mtypes(i)
      dof = molecules_getdof(mtype)
      endidx = begidx + dof - 1
      gencoords = gcmodels_getgencoords(sorbates(mtype)%gcoords(1))
      x(begidx:endidx) = gencoords(1:dof)
      begidx = endidx + 1
    Enddo
  End Subroutine bakersubs_sorbToX

  !-------------------------------------------------------
  ! Change from the "x" representation to the "sorbates"
  ! representation
  !-------------------------------------------------------
  Subroutine bakersubs_xToSorb(sorbates, simcell, x, mtypes)
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(in)               :: simcell
    Real(kind=RDbl), Dimension(:), Intent(in)  :: x
    Integer, Intent(In), Dimension(:)          :: mtypes

    Integer  :: nsorbs, i, mtype
    Integer  :: nvar, sorb, dof, begidx, endidx, natoms, sorbtype
    Real(kind=RDbl), Dimension(MAX_DOF) :: gencoords

    !** Number of baker sorbs
    nsorbs = Size(mtypes,1)

    begidx = 1
    Do i = 1, nsorbs
      mtype = mtypes(i)
      natoms = molecules_getnatoms(mtype)
      dof = molecules_getdof(mtype)
      endidx = begidx + dof - 1

      Call gcmodels_setgencoords(sorbates(mtype)%gcoords(1), x(begidx:endidx))

      begidx = endidx + 1

      !** Generate the xyz principal coordinates
      Call gcmodels_toxyz(sorbates(mtype)%gcoords(1), &
          sorbates(mtype)%coords(1:natoms,1)%rp)

      !** Generate the other coordinates
      Call simcell_pbc(simcell, sorbates(mtype)%coords(1:natoms,1)%rp, &
          sorbates(mtype)%coords(1:natoms, 1)%r, &
          sorbates(mtype)%coords(1:natoms, 1)%cr)    

    End Do 
  End Subroutine bakersubs_xToSorb


End Module bakersubs
