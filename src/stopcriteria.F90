!------------------------------------------------------- 
! This routine is useful for two purposes: 1) to generate a configuration 
! with required number of molecules 2) to stop a simulation if a 
! specific thing happens like energy going above particular value
!-------------------------------------------------------

Module stopcriteria
  Use config, Only: AtMolCoords, config_getnmoles, config_writerestartfile
  Use defaults, Only: strLen, lstrLen, RDbl, RSgl, MAX_MAPS, dashedline,&
      one, zero, d_ctrl_file, MAX_ATOMS
  Use file, Only: file_open
  Use molecules, Only: molecules_gettype, molecules_getnsorbs
  Use storestats, Only: Species_Stats, storestats_getnoncoul, &
      storestats_getcoul, storestats_getke, storestats_getintranrg
  Use utils, Only: split, getpath, allocErrDisplay, filesrchstr, int2str, &
      stripcmnt, toupper, cleanstring, toreal, toint


  Implicit None
  Save

  Private
  Public ::  STOP_CHECK, stopcriteria_init, stopcriteria_check

  Type Stop_Criteria
    Integer :: ncriteria
    ! variable for checking number of molecules
    Integer :: nmolec_sorb, nmolec_N, nmolec_type 
    Logical :: nrg_check
    ! if molec_type=-1, then sim stops if nmolec<nmolec_N, ctrlfile string LT
    ! if molec_type= 0, then sim stops if nmolec>nmolec_N, ctrlfile string GT
    ! if molec_type=-1, then sim stops if nmolec=nmolec_N, ctrlfile string EQ

    ! variable for Energy
    Integer :: nrg_type             ! same as nmolec_type  
    Real(kind=RDbl):: nrg_val
    Logical :: nmolec_check
  End Type Stop_Criteria

  Type (Stop_Criteria) :: criteria
  Logical :: STOP_CHECK=.False.

  !** The tag marking the beginning of the stop criteria section
  Character(len=strLen), Parameter    :: stop_criteria_tag = &
      "Stop Criteria"


Contains
  !---------------------------------------------------------------
  ! initialize stop criteria
  ! requires :
  ! filename - main control file
  !---------------------------------------------------------------
  Subroutine stopcriteria_init(filename )
    Character(*), Intent(In)                       :: filename
    Logical :: found_stop
    Character(len=strLen) :: line,string,stop_type
    Character(len=strLen),Dimension(10) :: fields
    Integer :: unitno, lineno, nfields, sorb,i

    unitno=file_open(filename)

    !** Find the stop section, "True" for rewinding the file
    lineno = filesrchstr(unitno, stop_criteria_tag, line, .True.)
    If (lineno == 0) Then
      ! no stop criteria specified
      STOP_CHECK=.False.
    Else
      STOP_CHECK=.True.
      criteria%nrg_check=.False.
      criteria%nmolec_check=.False.
      criteria%ncriteria=0
      Do i=1,5 ! maximum 5 criteria
        Read(unitno,'(a)') line
        line=cleanstring(line)
        nfields=split(line,fields,",")
        stop_type=toupper(cleanstring(fields(1)))
        Select Case(Trim(stop_type))
        case("NMOLECS")
          criteria%ncriteria=criteria%ncriteria+1
          sorb=molecules_gettype(cleanstring(fields(2)))

          If (sorb/=0) Then
            criteria%nmolec_sorb=sorb
            criteria%nmolec_N=toint(fields(3))
            criteria%nmolec_type=stopcriteria_gettype(cleanstring(fields(4)))
            criteria%nmolec_check=.True.
          Else
            Write(*,*) "Wrong molecule Name"
            Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
            Stop
          Endif
        case("ENERGY")
          criteria%ncriteria=criteria%ncriteria+1
          ! as of now only total energy is supported
          criteria%nrg_val=toreal(fields(2))
          criteria%nrg_type=stopcriteria_gettype(cleanstring(fields(3)))
          criteria%nrg_check=.True.
        Case Default
          Exit ! from infinite do loop
        End Select
      End Do
    Endif

    ! display here itself
    If (criteria%nmolec_check) Then
      Write(*,*) "-------------STOP Criteria -----------"
      Write(*,*) "simulation will be stopped if the number of molecules of"
      Select Case (criteria%nmolec_type)
      Case(-1) 
        string="goes below a value of " 
      Case(0) 
        string="equals a value of " 
      Case(1) 
        string="exceeds a value of " 
      End Select
      Write(*,*) " sorb - ", criteria%nmolec_sorb, Trim(string), &
          criteria%nmolec_N
    Endif

    If (criteria%nrg_check) Then
      Write(*,*) "-------------STOP Criteria -----------"
      Write(*,*) "simulation will be stopped if the total energy "
      Select Case (criteria%nrg_type)
      Case(-1) 
        string="goes below a value of " 
      Case(0) 
        string="equals a value of " 
      Case(1) 
        string="exceeds a value of " 
      End Select
      Write(*,*) Trim(string), criteria%nrg_val
    Endif

  End Subroutine stopcriteria_init


  !-------------------------------------------------------------
  ! ctrlfile uses strings like GT (greater than) LT (less than) 
  ! to say when to Stop. This function converts it to internal representation
  !-------------------------------------------------------------
  Integer Function stopcriteria_gettype(string)
    Character(len=strlen),Intent(in) :: string
    Select Case(Trim(string))
    case("GT")
      stopcriteria_gettype=1 
    case("LT")
      stopcriteria_gettype=-1
    case("EQ")
      stopcriteria_gettype=0
    Case default
      Write(*,*) "can not understand string", string
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select

  End Function stopcriteria_gettype
  !---------------------------------------------------------------
  ! checks stop criteria
  ! requires :
  !---------------------------------------------------------------
  Subroutine stopcriteria_check(species,spcstats,restartfile)
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(Species_Stats), Intent(In)             :: spcstats
    Character(len=strLen), Intent(in) :: restartfile
    Logical :: stopsim
    Integer :: nmolecs, i, j
    Real(kind=RDbl) :: total_u, unoncoul, ucoul, kin_nrg, intra_nrg
    stopsim=.False.
    If (criteria%nmolec_check) Then
      nmolecs=config_getnmoles(species,criteria%nmolec_sorb)
      Select Case (criteria%nmolec_type)
      Case(-1)
        If (nmolecs<criteria%nmolec_N) stopsim=.True.
      Case(0)
        If (nmolecs==criteria%nmolec_N) stopsim=.True.
      Case(1)
        If (nmolecs>criteria%nmolec_N) stopsim=.True.
      End Select
      If(stopsim) &
          Write(*,*) "Nmolecs stop criterion has been reached"
    Endif

    If (criteria%nrg_check) Then
      total_u=0.0
      Do i = 1, molecules_getnsorbs()
        Do j = i, molecules_getnsorbs()

          !** Write the potential energies
          unoncoul = storestats_getnoncoul(spcstats,i,j,'inst')
          ucoul = storestats_getcoul(spcstats,i,j,'inst')
          total_u = total_u+unoncoul+ucoul
        End Do

        !** Write the kinetic and total intramolecular energies
        kin_nrg = storestats_getke(spcstats,i,'inst')
        intra_nrg = storestats_getintranrg(spcstats,i,'inst')
        total_u = total_u+kin_nrg+intra_nrg
      End Do


      Select Case (criteria%nrg_type)
      Case(-1)
        If (total_u<criteria%nrg_val) stopsim=.True.
      Case(0)
        If (total_u==criteria%nrg_val) stopsim=.True.
      Case(1)
        If (total_u>criteria%nrg_val) stopsim=.True.
      End Select
      If(stopsim) &
          Write(*,*) "Energy stop criterion has been reached"

    Endif
    If (stopsim) Then
      Write(*,*) " Writing to restartfile : "//Trim(restartfile)
      Call config_writerestartfile(species,restartfile)
      Write(*,*) "Stopping because some of the stop criteria has been met"
      Stop
    Endif
  End Subroutine stopcriteria_check

End Module stopcriteria
