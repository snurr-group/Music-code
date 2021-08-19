!------------------------------------------------------- 
! This Module ghas no object defined as such
! This is an example module about how to get energies from music
! This should by used by main drivers
! call to ffcall may no tbe made from usual simulations, 
! because the subroutines here are not optimized for speed, 
! They are just written for convenience and educational purpose only
!-------------------------------------------------------

Module ffcall
  Use config, Only: AtMolCoords, config_getnmoles, config_writerestartfile, &
      config_isfixed, config_getnmoles, config_getnatoms, config_getaccel
  Use defaults, Only: strLen, lstrLen, RDbl, RSgl, MAX_MAPS, dashedline,&
      one, zero, d_ctrl_file, MAX_ATOMS, calToj, TOTAL_INDEX
  Use file, Only: file_open
  Use molecules, Only: molecules_gettype, molecules_getnsorbs
  Use storestats, Only: Species_Stats, storestats_getnoncoul, &
      storestats_getcoul, storestats_getke, storestats_getintranrg
  Use storebase, Only: EnergyPlus, storebase_sumintra, storebase_copy, &
      storebase_clean, storebase_display, storebase_init, storebase_disp, &
      storebase_chkintra, storebase_chkequal, storebase_zero, &
      storebase_nrg, storebase_totnrg, storebase_initcopy
  Use utils, Only: split, getpath, allocErrDisplay, filesrchstr, int2str, &
      stripcmnt, toupper, cleanstring, toreal, toint
  Use interact, Only: Interaction_Model, interact_init, interact_display, &
      interact_int
  Use storetop, Only: Forcefield_Results,storetop_display,storetop_sum, &
      storetop_zero, storetop_fastzero, storetop_fastsum, storetop_init, &
      storetop_initcopy, storetop_copy, storetop_clean, storetop_totnrg, &
      storetop_update, storetop_extract, storetop_allnrgs, storetop_chklink, &
      storetop_changenmoles, storetop_chksizes, storetop_scalenrgs, &
      storetop_chkequal, storetop_fillsub, storetop_nmolesDecr, &
      storetop_fastforces
  Use simcell, Only: SimCell_Params, simcell_display, simcell_init, &
      simcell_sampleCF
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Implicit None
  Save

  Private
  Public ::  ffcall_nrg, ffcall_force

  ! nothing for now
  Type FF_CALL
    Logical :: check
  End Type FF_CALL

  !** The tag marking the beginning of ffcall section in ctrlfile
  Character(len=strLen), Parameter    :: ffcall_tag = &
      "NOTHING FOR NOW"


Contains
  !---------------------------------------------------------------
  ! initialize 
  ! requires :
  ! filename - main control file
  !---------------------------------------------------------------
  Subroutine ffcall_init(filename )
    Character(*), Intent(In)                       :: filename
  End Subroutine ffcall_init


  !-------------------------------------------------------------
  ! calulate nrgs
  ! returns kj/mol
  !-------------------------------------------------------------
  Subroutine ffcall_nrg(imodel, scell, species, nrgtype, nrgval)
    Character(*),Intent(in) :: nrgtype
    Type(Interaction_Model),Intent(inout)                  :: imodel
    Type(SimCell_Params),Intent(in)            :: scell
    Type(AtMolCoords), Dimension(:), Pointer :: species
    Real(kind=RDbl),Intent(out)::nrgval
    Logical :: fast, recalc, nointra, calc_accels, want_intra
    Logical :: want_coul, want_ncoul
    Real(kind=RDbl)::aux_params(3)
    Logical :: success
    Type(EnergyPlus),Save :: ffout
    Logical, Save :: first_time=.True.
    Integer :: spc, nspc

    nspc=molecules_getnsorbs()
    ! some defaults
    fast=.True.
    recalc=.True.
    nointra=.False.
    want_intra=.True.
    want_coul=.True.
    want_ncoul=.True.
    calc_accels=.False.
    aux_params=(/1000000.0_RDbl,0.001_RDbl,300.0_RDbl/)

    If (first_time) Then
      Call storebase_init(ffout, imodel%nderivs, want_intra)
      first_time=.False.
    Endif

    Select Case(Trim(nrgtype))
    case("INTRA")

      ! move the call from here to other cases

      ! (/0,0,0/) means whole system
      success = interact_int(imodel, imodel%results(1), species,&
          scell, fast, recalc, nointra, calc_accels, &
          aux_params,(/0,0,0/),(/0,0,0/))

      If (.Not.success) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif

      ! makes sure that everything in storage is summed correctly
      Call storetop_fastsum(imodel%results(1))

      nrgval=zero
      Do spc=1,nspc
        ! would return spc-system interactions in ffout
        success = storetop_extract(imodel%results(1), want_intra, &
            want_ncoul, want_coul, ffout, (/spc,0,0/), (/0,0,0/) )
        Call storebase_sumintra(ffout,.False.)
        nrgval=nrgval+ffout%intranrg(TOTAL_INDEX)
!        Call storebase_display(ffout,6)

        ! would return interactions for mol1 of spc1 with mol2 of spc2 in ffout
        ! NOTE: if you want interactions with molecule-depth detail, then you
        ! must initialize the storage at MOL or ATM levels.
        !  success = storetop_extract(imodel%results(1), want_intra, &
        !    want_ncoul, want_coul, ffout, (/spc1,mol2,0/), (/spc2,mol2,0/) )
      End Do
      nrgval=nrgval*caltoj

    case("COUL")
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    case("NCOUL")
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    case("TOTAL")
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Case default
      Write(*,*) "can not understand string", nrgtype
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select

  End Subroutine ffcall_nrg

  !-------------------------------------------------------------
  ! calulate forces, only for single molecule
  ! not the best way to call force routine, just simple
  ! returns kcal/mol/ang
  ! force is 1-D array of 3*natoms
  !-------------------------------------------------------------
  Subroutine ffcall_force(imodel, scell, species, force)
    Type(Interaction_Model),Intent(inout)                  :: imodel
    Type(SimCell_Params),Intent(in)            :: scell
    Type(AtMolCoords), Dimension(:), Pointer :: species
    Real(kind=RDbl),Dimension(:),Intent(out)::force
    Logical :: fast, recalc, nointra, calc_accels, want_intra
    Logical :: want_coul, want_ncoul
    Real(kind=RDbl)::aux_params(3)
    Logical :: success
    Type(EnergyPlus) :: ffout
    Integer :: spc, nspc, nmoles, natoms, index, i, j
    Type(VecType), Dimension(:,:), Pointer  :: forces

    nspc=molecules_getnsorbs()
    ! some defaults
    fast=.True.
    recalc=.True.
    nointra=.False.
    want_intra=.True.
    want_coul=.True.
    want_ncoul=.True.
    calc_accels=.True.
    aux_params=(/1000000.0_RDbl,0.001_RDbl,300.0_RDbl/)

    Call storebase_init(ffout, imodel%nderivs, want_intra)

      ! (/0,0,0/) means whole system
      success = interact_int(imodel, imodel%results(1), species,&
          scell, fast, recalc, nointra, calc_accels, &
          aux_params,(/0,0,0/),(/0,0,0/))

      If (.Not.success) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif


    Do spc = 1,imodel%ff(1)%nspc

      If (config_isfixed(species(spc))) Cycle
      nmoles = config_getnmoles(species,spc)
      If (nmoles==0) Cycle
      natoms = config_getnatoms(species,spc)

      !** Get the species acceleration storage and zero it
      Nullify(forces)
      Call config_getaccel(species,spc,fast,forces)
      forces = VecType(0.0_Rdbl)


      !** Put the forces into the species structure
      Call storetop_fastforces(imodel%results(1),spc,nmoles,natoms,forces)

      !copy to 1D array
      index=0
      Do i=1,natoms
        Do j=1,3
          index=index+1
          force(index)=forces(i,1)%comp(j)
        End Do
      End Do

      !copy forces to force

      !** Convert to accelerations NO NEED
      !      Call config_conv2accel(species(spc),spc,fast)
      Nullify(forces)
    End Do

  End Subroutine ffcall_force



  !---------------------------------------------------------------
  ! checks stop criteria
  ! requires :
  !---------------------------------------------------------------
  Subroutine ffcall_display()
  End Subroutine ffcall_display

End Module ffcall
