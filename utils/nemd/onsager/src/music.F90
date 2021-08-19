!-----------------------------------------------------------------------
! This contains the main program AND SUBROUTINES 
! for calculating flux correlations functions for EMD
!-----------------------------------------------------------------------
Program onsager

  Use commandline, Only: commandline_init
  Use file, Only: file_open, file_settag
  Use mmath, Only: mmath_simpson
  Use defaults, Only: dashedline, strLen, d_ctrl_file, RDbl, one, zero, &
      kcalmole_kb

  Use utils, Only: genfilename,allocErrDisplay,int2str, stripcmnt, &
      split, toreal, toint
  Use vector, Only:  VecType, Assignment(=), Operator(*), Operator(+), &
      Operator(-), Operator(/), vector_display,vector_filedisplay, &
      mag , unitvec, IntVecType, vector_tripledotprod, vector_getnormsq, &
      vector_zerovec
  Implicit None

  !** file which contains all input details
  Character(len=strLen) :: ctrl_filename

  !** tells what to do. Ex:  run simulation?, help?, write sample file?
  Character(len=strLen) :: action

  Character(len=strLen)   :: fname
  Character(len=strLen), Dimension(strlen)     :: fields 
  Character(len=2*strLen) :: text ,outfile
  Character(len=strLen),Dimension(:),Pointer :: spcname

  Integer                 :: nfields, spc, nsorbs, confcount, error
  Integer                 :: ctrlunit, outunit, ios, n, n_fc, n_start, n_end
  Integer :: n_ct, n_st, i, j, unitno

  Integer,Dimension(:),Pointer :: spctype, nmoles

  Real(kind=RDbl) :: time, TN, T0, deltaT, newt
  Type(VecType), Dimension(:), Pointer :: vels

  Type(VecType), Dimension(:,:), Pointer :: comvel 
  Type(VecType) :: v0i, vtj
  Real(kind=RDbl), Dimension(:,:,:), Pointer :: fc,fc_x, fc_y, fc_z

  Real(kind=RDbl) :: maxT, dproduct, volume, Tk
  Real(kind=RDbl) :: dprodx,dprody,dprodz, x_integ, y_integ, z_integ, tot_integ
  Real(kind=RDbl) :: L_x, L_y, L_z, L_tot, prefac

  !-------------------------------------------------------------
  !
  ! Do basic initialization, based on command line input
  !
  !-------------------------------------------------------------
  Call commandline_init(ctrl_filename,action)

  If (action=="WriteSampleCtrlfile") Then
    Write(0,'(a)') dashedline
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Elseif(action/="DoSimulation") Then
    !** continue only if everything is alright 
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Stop
  Endif

  Call file_settag(ctrl_filename, d_ctrl_file)
  ctrlunit=file_open(ctrl_filename,110)



  !------------------------------------------------------------------
  !
  !                      Read the ctrlfile
  !
  !------------------------------------------------------------------
  Write(*,'(a,t40,2a)') " Read ctrlfile", " : ", Trim(ctrl_filename)
  Write(*,*) 

  ! get the name of the output file
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  outfile=Trim(Adjustl(fields(2)))
  Write(*,'(a,t40,2a)') " Outfile"," : ", Trim(outfile)

  ! get the maximum correlations time 
  Read(ctrlunit,'(a)') text
  nfields=split(text,fields,":")
  maxT=toreal((fields(2)))
  Write(*,'(a,t40,a,f10.4)') " Maximum Correlation Time in ps", " : ",maxT
  Write(*,*) 

  !** REad the outfile and figure out size and dimension
  Write(*,'(a,t40,2a)') " Reading output file (binary)"," : ", Trim(outfile)
  outunit=file_open(outfile,100)
  Read(outunit) volume,Tk
  Read(outunit) nsorbs
  Allocate(spctype(nsorbs),spcname(nsorbs),&
      nmoles(nsorbs),vels(nsorbs),STAT=error)
  If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

  Write(*,*) "this file contains com velocities of following species : "
  Do spc=1,nsorbs
    Read(outunit)  spctype(spc),nmoles(spc),spcname(spc)
    Write(*,'(a,i4,t25,2a,t55,a,i6)') " Species Type :", spctype(spc), "Name : ", &
        Trim(spcname(spc)), " Nmoles : ",nmoles(spc)
  End Do
  Write(*,*) 

  confcount=0
  Do
    Read(outunit,IOSTAT=ios) time, vels(1:nsorbs)
    If (ios==0) Then
      confcount=confcount+1
      If (confcount==1) T0=time
    Else
      Exit
    Endif
  End Do
  TN=time
  deltaT=(TN-T0)/(confcount-1)
  Write(*,*) "this file contains ", confcount, " data points"
  Write(*,*) "start time : ", T0, "end time : ", TN
  Write(*,*) "delta T between configs (ps): ", deltaT
  Write(*,*) "-------****** ASSUMPTION **********--------"
  Write(*,*) "this file contains data linear in time"
  Write(*,*) "-------****** ASSUMPTION **********--------"
  Write(*,*) "Simcell Volume (Ang^3) : ", volume
  Write(*,*) "Simulation Temp. (K)   : ", Tk

  ! Allocate vectors
  Allocate(comvel(nsorbs,confcount),STAT=error)


  !** read values into memory
  Rewind(outunit)
  Read(outunit) volume,Tk
  Read(outunit) nsorbs
  Do spc=1,nsorbs
    Read(outunit)  spctype(spc),nmoles(spc),spcname(spc)
  End Do

  Do n=1,confcount

    Read(outunit,IOSTAT=ios) time, comvel(1:nsorbs,n)
    ! check time
    newt=T0+(n-1)*deltaT

    If ((abs(newt/time-1))>1.0e-8) Then
      Write(*,*) "time in file not linear with index"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
  End Do


  !*** check how many points should be in the flux correlation array
  n_fc=(maxT/deltaT)+2
  Allocate(fc(nsorbs,nsorbs,n_fc),fc_x(nsorbs,nsorbs,n_fc),&
      fc_y(nsorbs,nsorbs,n_fc),fc_z(nsorbs,nsorbs,n_fc),STAT=error)

  !** check start point and end point
  n_start=1
  n_end=confcount-n_fc
  Write(*,*) " we will use ", n_end-n_start ,"points for each value of "//&
      "correlation time"


  ! evaluate fc - flux correlations
  ! loop over different values of correlation time
  Do n_ct=1,n_fc
    Write(*,*) "Calculating flux correlation for ", n_ct*deltaT, "picoseconds"
    fc(1:nsorbs,1:nsorbs,n_ct)=zero
    fc_x(1:nsorbs,1:nsorbs,n_ct)=zero
    fc_y(1:nsorbs,1:nsorbs,n_ct)=zero
    fc_z(1:nsorbs,1:nsorbs,n_ct)=zero
    ! loop over different simulation time regions 
    Write(*,*) n_start, n_end
    Do n_st=n_start, n_end
      ! loop over different phenomenlogical coefficients
      Do i=1,nsorbs
        v0i=comvel(i,n_st)
        Do j=1,nsorbs
          vtj=comvel(j,n_st+n_ct-1)
          dprodx=(v0i%comp(1)) * (vtj%comp(1))
          dprody=(v0i%comp(2)) * (vtj%comp(2))
          dprodz=(v0i%comp(3)) * (vtj%comp(3))
          !          dproduct=v0i*vtj
          dproduct=dprodx+dprody+dprodz
          fc(i,j,n_ct)=fc(i,j,n_ct)+dproduct
          fc_x(i,j,n_ct)=fc_x(i,j,n_ct)+dprodx
          fc_y(i,j,n_ct)=fc_y(i,j,n_ct)+dprody
          fc_z(i,j,n_ct)=fc_z(i,j,n_ct)+dprodz
          !          Write(*,*) fc(i,j,n_ct)
        End Do
      End Do

    End Do
    fc(1:nsorbs,1:nsorbs,n_ct)= fc(1:nsorbs,1:nsorbs,n_ct)/(n_end-n_start+1)
    fc_x(1:nsorbs,1:nsorbs,n_ct)=fc_x(1:nsorbs,1:nsorbs,n_ct)/(n_end-n_start+1)
    fc_y(1:nsorbs,1:nsorbs,n_ct)=fc_y(1:nsorbs,1:nsorbs,n_ct)/(n_end-n_start+1)
    fc_z(1:nsorbs,1:nsorbs,n_ct)=fc_z(1:nsorbs,1:nsorbs,n_ct)/(n_end-n_start+1)

  End Do


  ! output 3D flux-corr func to a file
  Do i=1,nsorbs
    Do j=1,nsorbs
      fname="FC_"//Trim(int2str(i))//Trim(int2str(j))
      unitno=file_open(fname)
      Do n_ct=1,n_fc
        Write(unitno,'(f12.6,e16.6)') (n_ct-1)*deltaT, fc(i,j,n_ct)
      End Do
      close(unitno)
    End Do
  End Do

  ! calculate 1/VkbT. Units : vol (ang^3), T(kelvin), kb = kcal/mol
  prefac=one/volume/Tk/kcalmole_kb

  ! do simpsons integration 
  Do i=1,nsorbs
    Do j=1,nsorbs
      x_integ=mmath_simpson(fc_x(i,j,1:n_fc),deltaT,n_fc)
      y_integ=mmath_simpson(fc_y(i,j,1:n_fc),deltaT,n_fc)
      z_integ=mmath_simpson(fc_z(i,j,1:n_fc),deltaT,n_fc)
      tot_integ=mmath_simpson(fc(i,j,1:n_fc),deltaT,n_fc)
      L_x=x_integ*prefac
      L_y=y_integ*prefac
      L_z=z_integ*prefac
      L_tot=tot_integ*prefac/3
      Write(*,'(a,2i2,4e11.3)') &
          "i,j, Lijx, Lijy, Lijz, Lij : ",i,j,L_x, L_y, L_z, L_tot
    end do
  end do

End Program onsager


