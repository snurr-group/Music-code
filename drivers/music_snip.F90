!-------------------------------------------------------------------
! This contains the driver that reads a coordinate file, snips a
! cluster for high level calculations out and returns a new file
! Command line options are:
!   -a -- use maximum angle criteria
!   -d -- use minimum distance criteria
!   -p -- use planarity criteria
!-------------------------------------------------------------------

Program music
  Use defaults, Only: strlen, RDbl, dashedline, d_res_file, d_ctrl_file, &
      ROUT_MAIN, one, lstrlen
  Use utils, Only: genfilename,int2str,allocerrdisplay,split,combine,ToUpper
  Use file, Only: file_settag
  Use commandline, Only : commandline_init
  Use vector, Only: VecType
  Use fundcell, Only: Fundamental_Cell, fundcell_display
  Use readstruc, Only: Structure,readstruc_init,readstruc_visxyz, &
      readstruc_natoms, readstruc_dumpcar, readstruc_expand, &
      readstruc_initcopy, readstruc_recenter
  Use snip, Only: snip_zeo

  Implicit None

  Integer                :: content_tag
  Integer                :: error,dunit,nfields,natoms,tatom
  Character(len=strLen)  :: action,filename,inputfile,method
  Character(len=strLen)  :: string1,string2,string3
  Type(Structure)        :: full,cluster,newstruc
  Type(Fundamental_Cell) :: fcell
  Character(len=lstrLen), Dimension(strLen) :: fields,chunks

  !----------------------------------------------
  ! Welcome the user and do the snipping
  !----------------------------------------------
  dunit = 6
  Write(dUnit,'(1x,a)') "Welcome to MuSiC (snip implementation)"

  !** interpret the command line
  Call commandline_init(inputfile,action,chunks)

  !** Get the structure
  Call readstruc_init(full,inputfile,fcell)

#ifdef DEBUG
  Call readstruc_dumpcar(full,'input.car',fcell)
  Call readstruc_expand(full,fcell,(/1,1,1/),newstruc)
  Call readstruc_visxyz(newstruc,'expanded.xyz')

  Call fundcell_display(fcell,2,dunit)
  Call readstruc_initcopy(newstruc,full)
  Call readstruc_recenter(newstruc,0,fcell)
  Call readstruc_visxyz(newstruc,'junk.xyz')
  Stop
#endif

  !** Visualize the structure
  filename = 'structure.xyz'
  Write(dunit,'(2x,2a)') 'Dumping full structure to: ',Trim(filename)
  Call readstruc_visxyz(full,filename)

  !** Get the tetrahedral Aluminum atom number
  Do tatom = 1,readstruc_natoms(full)
    If (ToUpper(full%elements(tatom)) == 'AL') Exit
  End Do
  If (tatom > readstruc_natoms(full)) Stop
  Write(dUnit,'(a,i4)') ' T-atom number: ',tatom

  !** Choose the method based on the command line option
  Select Case (ToUpper(chunks(1)))
  Case('A')
    method = 'MAXANGLE'
  Case('D')
    method = 'MINDIST'
  Case('P')
    method = 'INPLANE'
  Case Default
    Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
        ' Could not interpret command line method type: ',Trim(chunks(1))
    Stop
  End Select

  !** Create the cluster from the full structure
  Write(dUnit,'(3a)') ' Using "',Trim(method),'" method for cluster choice'
  Call snip_zeo(full,fcell,tatom,1,method,cluster)

  !** Visualize the cluster structure
  filename = 'cluster.xyz'
  Write(dunit,'(2x,2a)') 'Dumping cluster structure to: ',Trim(filename)
  Call readstruc_visxyz(cluster,filename,'f14.9')

  !** dump cluster to .car file for later use
  filename = 'cluster.car'
  Write(dunit,'(2x,2a)') 'Dumping cluster structure to: ',Trim(filename)
  Call readstruc_dumpcar(cluster,filename,fcell)

End Program music
