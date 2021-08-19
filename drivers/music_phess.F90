!--------------------------------------------------------
! This contains the main program for Hessian analysis
!--------------------------------------------------------

Program music
  Use defaults, Only: strlen, RDbl, dashedline, d_res_file, d_ctrl_file, &
      ROUT_MAIN, one, lstrlen
  Use utils, Only: genfilename,int2str,allocerrdisplay,split,combine
  Use file, Only: file_settag
  Use commandline, Only : commandline_init
  Use atom, Only: atom_init, atom_display, atom_sampleCF
  Use molecules, Only: molecules_init, molecules_sampleCF, &
      molecules_display
  Use vector, Only: VecType
  Use readstruc, Only: Structure,readstruc_init,readstruc_visxyz, &
      readstruc_natoms
  Use hessian, Only: HessianInfo,hessian_removert,hessian_eigeninfo, &
      hessian_vismodes,hessian_display,hessian_dispmode,hessian_dispevalues
  Use readhess, Only: Hessian_Specs,readhess_init,readhess_units

  Implicit None

  Integer                :: content_tag
  Integer                :: error,dunit,nfields,natoms
  Character(len=strLen)  :: command,action
  Character(len=lstrLen) :: hessian_location,filename,units
  Character(len=strLen)  :: string1,string2,string3
  Character(len=1000)    :: outstring
  Type(Structure)        :: struc
  Type(HessianInfo)      :: hess
  Type(Hessian_Specs)    :: specs
  Character(len=lstrLen), Dimension(strLen) :: fields,chunks

  !----------------------------------------------
  ! Welcome the user
  !----------------------------------------------
  dunit = 6

  Write(dUnit,'(1x,a)') "Welcome to MuSiC (Hessian analysis)"

  !-----------------------------------------------------------------------
  ! Interpret the command and act
  ! Example: QMPOTHESS@qmpot.hessian:145-179@host.car:145-179
  ! This reads the Hessian file 'qmpot.hessian' of type QMPOTHESS and
  ! extracts just the sub-Hessian formed by atom numbers 145-179.  The
  ! last chunk 'host.car:145-179' specifies the atomic coordinates.
  !-----------------------------------------------------------------------

  Call commandline_init(command,action)

  nfields = split(command,fields,'@')

  If (nfields /= 3) Then
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Stop
  End If

  !** Get the structure
  Call readstruc_init(struc,fields(3))

  !** Visualize the structure
  filename = 'structure.xyz'
  Write(dunit,'(2x,2a)') 'Dumping structure to: ',Trim(filename)
  Call readstruc_visxyz(struc,filename)

  !** Get the Hessian 
  hessian_location = combine(fields(1:2),'@')
  Write(dunit,'(2x,3a)',Advance='No') 'Reading Hessian from ', &
      Trim(hessian_location),' ... '
  Call readhess_init(specs,hess,hessian_location)
  Write(dunit,'(a)') 'Done'

  units = readhess_units(specs)
  Write(dunit,'(2x,2a)') 'Hessian has units of ',Trim(units)

  !** Remove the rigid rotational and translational components of Hessian
!  Call hessian_removert(hess,struc%coords)

  !** Perform the diagonalization and get the eigenvectors and values
  Write(dunit,'(2x,a)',Advance='No') 'Performing diagonalization ... '
  Call hessian_eigeninfo(hess)
  Write(dunit,'(a)') 'Done'

  !** Display the results
  filename = 'einfo'
  Call hessian_vismodes(hess,struc,filename,outstring)

  Write(dunit,'(2x,a)') Trim(outstring)

End Program music
