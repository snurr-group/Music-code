!------------------------------------------------------------------------------
! This is the movie module. It will write the systm contents to a file 
! so it can be played as a movie.
!------------------------------------------------------------------------------
Module movie

  Use defaults, Only: RDbl, strLen
  Use file, Only: file_getunit,file_open
  Use utils, Only: filesrchstr, stripcmnt, toupper, toint, split, &
      allocerrdisplay,isfileopen, int2str
  Use vector, Only: VecType, IntVecType,vector_display
  Use config, Only: AtMolCoords, config_getnmoles, config_getnatoms, &
      config_isfixed, config_getatype
  Use simcell, Only: SimCell_Params,simcell_maptouc
  Use atom, Only: atom_getsymbol
  Use molecules, Only: molecules_getnsorbs
  Use visxyz, Only: visxyz_dump, XYZ_Entry, visxyz_make, visxyz_molabel1

  Implicit None

  Private
  Public :: MovieInfo, movie_init, movie_display, movie_makemovie

  Type MovieInfo
    Character(len=strLen) :: filename
    Logical :: dumpFixed
    Logical :: dumpAll
    Integer :: start,frameno
    Integer :: finish
    Integer :: step
    Integer :: nx,ny,nz
    Logical :: record
  End Type MovieInfo

  Character(len=5), Parameter               :: movie_display_format = 'f8.3'

Contains

  !----------------------------------------------------------------------------
  ! Initialize the movie routine, i.e., read in info from control file
  ! Requires:  mparams -- movie parameters
  !            ctrl_filename -- control filename to read information from
  !            tag -- the control file tag for the movie section
  !----------------------------------------------------------------------------
  Subroutine movie_init(mparams,ctrl_filename,tag)
    Type(MovieInfo), Intent(InOut)    :: mparams
    Character(*), Intent(In)          :: ctrl_filename
    Character(*), Intent(In)          :: tag

    Integer                                  :: unitno, useMovie, nfields
    Character(len=strLen)                    :: line
    Character(len=strLen), Dimension(strLen) :: fields

    unitno = file_open(ctrl_filename,110)

    !** Find the movie section
    useMovie = filesrchstr(unitno,tag,line,.True.)

    If (useMovie /= 0) Then
      mparams%frameno = 0

      Read(unitno,'(a)') line
      line = stripcmnt(line)
      mparams%filename = trim(line)

      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields,',')
      mparams%start = toint(fields(1))
      mparams%finish = toint(fields(2))

      Read(unitno,*) mparams%step
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      If ((toupper(line) == "YES")) Then
        mparams%dumpFixed = .True.
      Else
        mparams%dumpFixed = .False.
      End If

      !** get the unit cells of the FIXED species to display 
      Read(unitno,'(a)') line
      line = stripcmnt(line)
      nfields = split(line,fields,',')
      If (nfields /= 3) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " Expected to read 3 integers for FIXED species display"
        Stop
      End If
      mparams%nx = toint(fields(1))
      mparams%ny = toint(fields(2))
      mparams%nz = toint(fields(3))

      mparams%record = .True.
    Else
      mparams%record = .False.
    End If

  End Subroutine movie_init

  !----------------------------------------------------------------------------
  ! Writes the appropriate coordinates to the movie file
  ! Requires:  mparams -- movie parameters
  !            species -- configuration information
  !            step -- iteration number for labeling
  !            time -- simulation time for labeling
  !            optcomment -- an optional comment for the xyz frame
  !----------------------------------------------------------------------------
  Subroutine movie_makemovie(mparams,species,simcell,step,time,optcomment)
    Type(MovieInfo), Intent(InOut)              :: mparams
    Type(AtMolCoords), Dimension(:), Intent(In) :: species
    Type(SimCell_Params), Intent(In)            :: simcell
    Integer, Intent(In)                         :: step
    Real(kind=RDbl), Intent(In)                 :: time
    Character(*), Intent(In), Optional          :: optcomment

    Integer                        :: a,m,spc,unitno,error
    Integer                        :: natoms,nmoles,nspc,totala,nentries
    Character(len=2)               :: element
    Character(len=strLen)          :: comment,string
    Type(VecType)                  :: dummy
    Type(IntVecType)               :: cvec
    Type(XYZ_Entry), Dimension(:), Allocatable  :: entries

    !** Exit if not recording
    If (.Not.mparams%record) Return

    If ((step < mparams%start).Or.(step > mparams%finish)) Return

    If (Mod(step,mparams%step) /= 0) Return

    !** If this is the movie's first step, then open the file
    unitno = isfileopen(mparams%filename)
    If (unitno < 0) Then
      unitno = file_getunit(mparams%filename)
      Open(file=mparams%filename, unit=unitno)
      Write(*,'(3a,i4)') __FILE__,": Opening file ", &
          Trim(mparams%filename),unitno
    End If

    nspc = molecules_getnsorbs()

    !** Get an outside estimate of the number of atoms present in the frame 
    totala = 0
    Do spc = 1, nspc
      !** Check to see if the species is fixed
      If (config_isfixed(species(spc)).And.(.Not.mparams%dumpFixed)) Cycle

      nmoles = config_getnmoles(species,spc)
      natoms = config_getnatoms(species,spc)
      totala = nmoles*natoms + totala
!LC      Write(*,*) spc,nmoles,natoms
    End Do

    !** Allocate the entries for the xyz frame
    Allocate(entries(totala), STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"sorbates")

    !** Create the comment for the frame
    mparams%frameno = mparams%frameno + 1
    string = int2str(mparams%frameno)
    Write(comment,'(3a,f8.3)') 'Frame ',Trim(string),'  Time: ',time
    If (Present(optcomment)) Then
      Write(comment,'(a,1x,a)') Trim(comment),Trim(optcomment)
    End If

    !** Create the entries for the movie frame
    nentries = 0
    Do spc = 1,nspc

      !** Check to see if the sorbate is fixed.
      If (config_isfixed(species(spc)).And.(.Not.mparams%dumpFixed)) Cycle
      nmoles = config_getnmoles(species,spc)
      natoms = config_getnatoms(species,spc)

      Do m = 1,nmoles
        Do a = 1,natoms
          If (config_isfixed(species(spc))) Then
            dummy = simcell_maptouc(simcell,species(spc)%coords(a,m)%rp,cvec)
            If (((cvec%comp(1) > mparams%nx-1)).Or. &
                ((cvec%comp(2) > mparams%ny-1)).Or. &
                ((cvec%comp(3) > mparams%nz-1))) Cycle
          End If
          nentries = nentries + 1
          element = atom_getsymbol(config_getatype(species,spc,m,a))
          entries(nentries) = visxyz_make(species(spc)%coords(a,m)%rp,element)
          !** debugging example, makes each molecule a different color
!          entries(nentries) = visxyz_make(species(spc)%coords(a,m)%rp,element,m)
        End Do
      End Do

    End Do

    !** Dump frame to file
    Call visxyz_dump(entries(1:nentries),mparams%filename, &
        movie_display_format,comment)

    Deallocate(entries)

  End Subroutine movie_makemovie

  !----------------------------------------------------------------------------
  ! Display movie information
  !----------------------------------------------------------------------------
  Subroutine movie_display(mparams,unit,indent)
    Type(MovieInfo), Intent(In) :: mparams
    Integer, Intent(In)         :: unit, indent

    Character(len=indent) :: blank
    Character(len=strLen) :: string

    If (.Not.mparams%record) Return

    blank = Repeat(' ',indent)

    Write(unit,'(2a)') blank,"Movie Information"
    Write(unit,'(2x,3a)') blank,    "Movie output file : ", &
        Trim(mparams%filename)
    string = int2str(mparams%start)
    Write(unit,'(2x,3a)') blank,"Starting step     : ",Trim(string)
    string = int2str(mparams%finish)
    Write(unit,'(2x,3a)') blank,"Ending step       : ",Trim(string)
    string = int2str(mparams%step)
    Write(unit,'(2x,3a)') blank,"Steps btwn frames : ",Trim(string)
    Write(unit,'(2x,2a,l1)') blank,  "Dumping zeolite   : ",mparams%dumpFixed

  End Subroutine movie_display

End Module movie
