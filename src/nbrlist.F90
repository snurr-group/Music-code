!---------------------------------------------------------------------------- 


Module nbrlist
  Use config, Only: AtMolCoords, config_getnmoles, config_isfixed
  Use defaults, Only: strLen, lstrLen, RDbl, RSgl, MAX_MAPS, dashedline,&
      one, zero, d_ctrl_file, MAX_ATOMS, dbgflag
  Use file, Only: file_open
  Use molecules, Only: molecules_gettype, molecules_getnsorbs, &
      molecules_getnatoms
  Use simcell, Only: Simcell_Params, simcell_getell, simcell_getnuc, &
      simcell_minimage
  Use utils, Only: split, getpath, allocErrDisplay, filesrchstr, int2str, &
      stripcmnt, findint
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getnorm, mag

  Implicit None
  Save

  Private
  Public :: Nbrlist_Params, NlistMolecSys_info, NlistMolecSpc_info,&
      SpeciesSys_info, nbrlist_init, nbrlist_update, nbrlist_getlist

  Type NlistMolecSpc_info
    Integer :: listsize
    Integer,Dimension(:), Pointer :: mol_spc
  End Type NlistMolecSpc_info

  Type NlistMolecSys_info
    Type(NlistMolecSpc_info), Dimension(:), Pointer :: mol_sys
  End Type NlistMolecSys_info

  Type SpeciesSys_info
    Type(NlistMolecSys_info), Dimension(:), Pointer :: spc_sys
  End Type SpeciesSys_info



  Type Nbrlist_Params
    Integer                            :: nspc
    Real :: cutoff, shell, outershellSqrd,updateT
    Real :: last_updateT
    Integer :: maxnbrs
    Type(SpeciesSys_info), Dimension(:), Pointer :: sys
  End Type Nbrlist_Params

  Character(len=strLen) :: default_nbrlist_tag=  "Neighborlist Info"
  Integer            :: number_of_cavity_maps=0

  Logical :: cavitylist_initialized=.False.

Contains

  !------------------------------------------------------------------------
  ! Initialize the 'neighbor list'.  If one is already initilized point to 
  ! the already initialized one
  ! Requires: cavitylist -- the empty cubes  object
  !           sorbates, scell -- the coordinates, simulation cell
  !           ctrl_filename -- the file with init info
  !           opt_tag -- tag marking the cavity bias section
  !------------------------------------------------------------------------
  Subroutine nbrlist_init(list, species, scell, ctrl_filename, opt_tag) 
    Type(Nbrlist_params), Pointer              :: list
    Type(AtMolCoords), Dimension(:), Intent(in)  :: species
    Type( Simcell_Params), Intent(in)            :: scell


    Character(*),                  Intent(in) :: ctrl_filename
    Character(*), Optional, Intent(in)        :: opt_tag
    Real(kind=RDbl), Dimension(3)             :: ell

    Character(len=strLen) :: tag, srchstr
    Integer               :: unitno, error, lineno, spc, nspc, m, nmoles, spc2

    Integer :: maxnbrs
    If (Present(opt_tag)) Then
      tag=opt_tag
    Else
      tag=default_nbrlist_tag
    End If


    Nullify(list)

    !** Open the ctrl_file if it is not opened and find the reqd section
    unitno = file_open(ctrl_filename)
    Rewind(unitno)
    lineno = filesrchstr(unitno, tag, srchstr)
    If (lineno == 0) Then
      Return
    Else
      Allocate(list,STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End If

    Write(*,*) "reading neighbor list info from ctrlfile"
    !** Read the important fields from ctrlfile and initialize
    Read(unitno,*) list%cutoff
    Read(unitno,*)list%shell
    Read(unitno,*)list%updateT
    Read(unitno,*)list%maxnbrs
    list%outershellSqrd=(list%cutoff+list%shell)**2

    Write(*,*)list%cutoff
    Write(*,*)list%shell
    Write(*,*)list%updateT
    Write(*,*)list%maxnbrs

    Write(*,*) "finished reading neighbor list info from ctrlfile"

    Write(*,*) "Allocating memory for neighbor list"
    nspc=molecules_getnsorbs()
    list%nspc=nspc
    Allocate(list%sys(nspc),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** allocate nighborlist for eachspecies
    Do spc=1,nspc
      nmoles=config_getnmoles(species,spc)
      ! initialize speciessys_info
      Allocate(list%sys(spc)%spc_sys(nmoles),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Do m=1,nmoles
        Allocate(list%sys(spc)%spc_sys(m)%mol_sys(nspc),STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        ! the last loop
        Do spc2=1,nspc
          maxnbrs=list%maxnbrs
          Allocate(list%sys(spc)%spc_sys(m)%mol_sys(spc2)%mol_spc(maxnbrs)&
              ,STAT=error)
          If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
          list%sys(spc)%spc_sys(m)%mol_sys(spc2)%mol_spc(1:maxnbrs)=0
          list%sys(spc)%spc_sys(m)%mol_sys(spc2)%listsize=0

        End Do
      End Do
    End Do
    list%last_updateT=zero
  End Subroutine nbrlist_init

  !------------------------------------------------------------------------
  ! updates the 'neighbor list'.
  ! time is optional, if time is not passed always updates, and sets time=0
  !------------------------------------------------------------------------
  Subroutine nbrlist_update(list, species, scell,time) 
    Type(Nbrlist_params), Pointer              :: list
    Type(AtMolCoords), Dimension(:), Intent(in)  :: species
    Type( Simcell_Params), Intent(in)            :: scell
    Real(kind=RDbl) , Intent(in),Optional :: time
    Integer :: nspc, spc1, spc2, m1, m2, nmoles1, nmoles2, m2_begin
    !SDEBUG 
    Integer               :: dispunit, dispmolec, n, nnbs
    !SDEBUG

    Integer, save ::counter=0 
    Logical :: nbrflag
    If (Present(time)) Then
      If ((time-list%last_updateT)<list%updateT) Then
        Return
      Else
        list%last_updateT=time
      Endif
    Else
        list%last_updateT=zero
    Endif
    counter =counter + 1             
!    Write(*,*) "Updating nbrlist",counter

    nspc=list%nspc
    Do spc1=1,nspc 
      If (config_isfixed(species(spc1))) Cycle
      nmoles1=config_getnmoles(species,spc1)
      Do m1=1,nmoles1
        Call nbrlist_makezero(list,spc1,m1)
      End Do
    End Do
    Do spc1=1,nspc
      If (config_isfixed(species(spc1))) Cycle
      Do spc2=spc1,nspc
        If (config_isfixed(species(spc2))) Cycle
        nmoles1=config_getnmoles(species,spc1)
        Do m1=1,nmoles1
          !          Write(*,*) "determining nbrs of" , m1
          nmoles2=config_getnmoles(species,spc2)
          m2_begin=1
          If (spc1==spc2) m2_begin=m1+1
          Do m2=m2_begin,nmoles2
            Call nbrlist_chknbr(list,species,scell,spc1,m1,spc2,m2,nbrflag)
            If (nbrflag) Then
              Call nbrlist_addnbr(list,spc1,m1,spc2,m2)
              Call nbrlist_addnbr(list,spc2,m2,spc1,m1)
              !              Write(*,*) "found neighbor", m2
              !              Write(*,*) 'distance is' , mag(species(spc1)%coords(1,m1)%r-species(spc2)%coords(1,m2)%r)
            Else
              !              Write(*,*) "Too far ", m2
              !              Write(*,*) 'distance is' , mag(species(spc1)%coords(1,m1)%r-species(spc2)%coords(1,m2)%r)
            Endif
          End Do
        End Do
      End Do
    End Do
    !SDEBUG
    ! write list of neighbors to a file
!!$    If (dbgflag) Then
!!$      dispunit=file_open("nbr.xyz")
!!$      Write(*,*) " going to display neighborlist of atom, which one do you want"
!!$      Read(*,*) dispmolec
!!$      nnbs=list%sys(1)%spc_sys(dispmolec)%mol_sys(1)%listsize
!!$      Write(dispunit,*) nnbs
!!$      Write(dispunit,*) "neighbors_of_atom_",dispmolec
!!$      Do n=1,nnbs
!!$        m2=list%sys(1)%spc_sys(dispmolec)%mol_sys(1)%mol_spc(n)
!!$        Write(dispunit,'(a,3f15.4)') "Ne", species(1)%coords(1,m2)%r
!!$      End Do
!!$    Endif
    !SDEBUG
  End Subroutine nbrlist_update

  !------------------------------------------------------------------------
  ! Make number of neighbors equal to zero
  !------------------------------------------------------------------------
  Subroutine nbrlist_makezero(list, spc, m) 
    Type(Nbrlist_params), Pointer              :: list
    Integer, Intent(in) :: spc, m
    Integer :: nspc,spc1, spc2, m1

    nspc=molecules_getnsorbs()

    list%sys(spc)%spc_sys(m)%mol_sys(1:nspc)%listsize=0


  End Subroutine nbrlist_makezero

  !------------------------------------------------------------------------
  ! checks whether two molecules are neighbors
  !------------------------------------------------------------------------
  Subroutine nbrlist_chknbr(list, species, scell, spc1, m1, spc2, m2, nbrflag) 
    Type(Nbrlist_params), Pointer              :: list
    Type(AtMolCoords), Dimension(:), Intent(in)  :: species
    Type( Simcell_Params), Intent(in)            :: scell
    Integer,Intent(in) :: spc1, spc2, m1, m2
    Logical,Intent(out) :: nbrflag
    Integer :: na1, na2, a1, a2
    Type(VecType) :: diff
    Real(kind=RDbl) :: distSqrd
    na1=molecules_getnatoms(spc1)
    na2=molecules_getnatoms(spc2)
    nbrflag=.False.
    Do a1=1,na1
      Do a2=1,na2
        Call simcell_minimage(scell, species(spc1)%coords(a1,m1)%r, &
            species(spc2)%coords(a2,m2)%r,diff, distSqrd)
!!$        Write(*,*) "r1 ", species(spc1)%coords(a1,m1)%r
!!$        Write(*,*) "r2 ", species(spc2)%coords(a2,m2)%r
!!$        Write(*,*) "rp1 ", species(spc1)%coords(a1,m1)%rp
!!$        Write(*,*) "rp2 ", species(spc2)%coords(a2,m2)%rp
!!$        Write(*,*) "diff",diff
!!$        Write(*,*) distSqrd, list%outershellSqrd
!!$
        If (distSqrd<list%outershellSqrd) Then
          nbrflag=.True.
          Exit
        Endif
      End Do
      If (nbrflag) Exit
    End Do
  End Subroutine nbrlist_chknbr


  !------------------------------------------------------------------------
  ! adds 2nd molecule to the neighbor list of first
  !------------------------------------------------------------------------
  Subroutine nbrlist_addnbr(list, spc1, m1, spc2, m2) 
    Type(Nbrlist_params), Pointer              :: list
    Integer,Intent(in) :: spc1, spc2, m1, m2
    Integer :: cn

    ! cn=current number of neighbors
    cn=    list%sys(spc1)%spc_sys(m1)%mol_sys(spc2)%listsize
    cn=cn+1
    list%sys(spc1)%spc_sys(m1)%mol_sys(spc2)%listsize=cn

    If (cn>list%maxnbrs) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif
    list%sys(spc1)%spc_sys(m1)%mol_sys(spc2)%mol_spc(cn)=m2
  End Subroutine nbrlist_addnbr


  !------------------------------------------------------------------------
  ! returns pointer to neighbor array of one moleule(spc1,m1) against spc2
  !------------------------------------------------------------------------
  Subroutine nbrlist_getlist(list, spc1, m1, spc2, listptr, list_size) 
    Type(Nbrlist_params), Pointer              :: list
    Integer,Intent(in) :: spc1, spc2, m1
    Integer, Dimension(:) , Pointer :: listptr
    Integer , Intent(out) :: list_size

    list_size =    list%sys(spc1)%spc_sys(m1)%mol_sys(spc2)%listsize
    listptr=> list%sys(spc1)%spc_sys(m1)%mol_sys(spc2)%mol_spc
  End Subroutine nbrlist_getlist

End Module nbrlist
