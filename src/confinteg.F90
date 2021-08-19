!-----------------------------------------------------------------------------
! This is the main module which makes all the calls required for doing 
! configurational integral of linear alkanes using bond-bending, bond-torsion 
! and intrapair potentials 
!-----------------------------------------------------------------------------
Module confinteg

  Use defaults, Only: RDbl, strLen, dashedline,dashedline2, twopi, d_res_file &
      ,d_con_file, hplanck, MAX_SORBS, NO_OF_INTRA_POTS, STRETCH_INDEX,&
      BENDING_INDEX,TORSION_INDEX,CONSTRAINT_INDEX,INTRAPAIR_INDEX, &
      INTRACOUL_INDEX, TOTAL_INDEX,Rgas, Nav, scalepe, kcalmole_kb, &
      zero,one, pi, degTorad
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay,toreal
  Use file, Only: file_getunit,file_gettype,file_open
  Use random,Only :rranf,random_gaussian,random_gettrial,random_init
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use config, Only: AtMolCoords, config_getnmoles, config_getnatoms
  Use simcell, Only: SimCell_Params, simcell_getvolume
  Use molecule, Only : MolecularParams,molecule_getgcmodeltype
  Use molecules, Only: molecules_getmass, molecules_getcoul &
      ,molecules_getnoncoul, molecules_getnsorbs, molecules_updateenergyss, &
      molecules_getpointer, molecules_updateenergy
  Use ssnoncoul, Only : ssnoncoulparams
  Use ssbasic, Only : ssbasic_intraGetinteraction
  Use branchedcoords, Only : branchedcoords_toxyz, branchedcoords_place, NodeType
  Use intramolecular, Only : intramolecular_getint

  Implicit None
  Save

#ifdef PRIVATE
  Private
  Public :: confinteg_init,confinteg_initdisplay,CONFINTEG_Params,&
      confinteg_getInteg
#endif


  ! The tag marking the beginning of the CONFINTEGsection
  Character(len=strLen), Parameter    :: default_confinteg_tag = &
      "Configurational Integration Info"

  Type CONFINTEG_Params
    Real(kind=RDbl),Dimension(:),Pointer        :: temp_k ! temp. in Kelvin
    Real(kind=RDbl),Dimension(:),Pointer        :: beta ! temp. in Kelvin
    Integer    :: ntrials         !** no of trial configurations
    Integer    :: nTemps          !** number of temperatures
!    Integer    :: angletrials     !** number of trial angles OBSOLETE
    Real(kind=RDbl) :: volume     !** sim volume
    Real(kind=RDbl) :: bias_k,bias_theta0,bias_Temp,bias_sig  !** biasing values 
                                                        !** for angles

    Type(MolecularParams),Pointer :: molecPtr
    Character(len=2*strLen) :: description,molecName
    Character(len=2*strLen) :: outfilename
  End Type CONFINTEG_Params

Contains

  !----------------------------------------------------------
  ! Initializes the various CONFINTEGparameters from the control
  ! file "ctrl_file"
  !----------------------------------------------------------
  Subroutine confinteg_init(params, sorbates, simcell, ctrl_filename, &
      opt_confintegtag)
    Type(CONFINTEG_Params)         ::params
    Type(AtMolCoords), Dimension(:), Intent(inout)   :: sorbates
    Type(SimCell_Params), Intent(inout)  :: simcell
    Character(*), Intent(in)  :: ctrl_filename
    Character(*), Optional, Intent(in) :: opt_confintegtag
    Character(len=strLen)     :: eostag, tag, line,text
    Character(len=strLen),Dimension(10)     :: fields

    Integer   :: unitno, confinteglineno, nsorbs, error, i, j,blocksize
    Integer   :: randominit_seed,nfields
    Real(kind=RDbl)     :: volume,lowT,highT
    
    If (Present(opt_confintegtag)) Then
      tag = opt_confintegtag
    Else
      tag = default_confinteg_tag
    End If
    
    !** Open the ctrl_file if it is not opened
    unitno = file_open(ctrl_filename)
    Rewind(unitno)
    
    !** Find the CONFINTEGsection
    confinteglineno = filesrchstr(unitno, tag, line)
    If (confinteglineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif
    Read(unitno, *) params%description
    Read(unitno,*) params%molecName
    params%molecName=   Trim( params%molecName)
    Read(unitno, *) params%ntrials     ! Trial number of configs
!    Read(unitno, *) params%angletrials     ! Trial number of angles
    Read(unitno, *) params%nTemps     ! Number of Temperatures
    Read(unitno, *) lowT,highT ! in Kelvin
    Read(unitno,*) params%outfilename
    !** Initializing the random number generator
    Read(unitno,*) randominit_seed
    Call random_init(randominit_seed)
!    Read(unitno,'(a)') text
!    Write(*,*) text
!    text=trim(stripcmnt(text))
!    nfields=split(text,fields,",")
!    params%bias_k=toreal(fields(1))
!   params%bias_theta0=toreal(fields(2))
!   params%bias_Temp=toreal(fields(3))
    Read (unitno,*) params%bias_k, params%bias_theta0, params%bias_Temp
    Write(*,*) " Finished reading the confinteg_init section"

    !** Calculate angle-bias sigma and ave
    params%bias_theta0=params%bias_theta0 * degToRad
    params%bias_k=params%bias_k*4184.00   !J/rad^2
    params%bias_sig=1/Sqrt(params%bias_k/(Rgas*params%bias_Temp))
    Write(*,*) "sigma is ",params%bias_sig/degToRad
    Write(*,*) degToRad,Rgas
    
    
    params%outfilename=    Trim(params%outfilename)
    
    !** Allocate and fill the temperature array
    Allocate(params%temp_k(params%nTemps),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"confinteg params")
    Allocate(params%beta(params%nTemps),STAT=error)
    If (error/=0)Call allocErrDisplay(__FILE__,__LINE__,"confinteg params")

    params%temp_k(1)=lowT
    If (params%nTemps>1) Then
      Do i=2, params%nTemps
        params%temp_k(i)=lowT+(i-1)*(highT-lowT)/(params%nTemps-1)
      End Do
    Endif
    
    Do i=1,params%nTemps
      params%beta(i)= one/params%temp_k(i)/kcalmole_kb
    End Do
    
    Allocate(params%molecPtr,STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Call molecules_getpointer(1,params%molecPtr)
    
    !** Generate the fugacitylist 
    params%volume = simcell_getvolume(simcell)
    Close(unit=unitno)

  End Subroutine confinteg_init
  
  !------------------------------------------------
  ! Calculates the configurational integral Z/Omega for Ideal gas chains 
  ! Omega is the integral over gen.coords like ->
  ! \integ cos(theta) dtheta : theta=0,pi      = 2
  ! \integ dphi              : phi=-pi,pi      = 2*pi
  !------------------------------------------------
  Subroutine confinteg_getInteg(sorbate, params)
    Type(AtMolCoords)                      :: sorbate
    Type(CONFINTEG_Params), Intent(inout)  :: params
    
    Integer    :: i, j, k, sorbtype, molec,sorbno, natoms, atom_num, unitno
    Integer    :: error
    Real(kind=RDbl)   :: energy, biasfactor, sum_e, sum_nrg, sum_bias, ratio
    Real(kind=RDbl)   :: beta,intrapair,ipot
    Real(kind=RDbl),Dimension(NO_OF_INTRA_POTS)   :: potlist    
    Real(kind=RDbl),Dimension(:),Pointer   :: expsum
    Character (len=strLen) :: gctype
 !SDEBUG
    Integer :: outputcount=0
 !SDEBUG
    Logical :: atleast_3,fast_flag,ljflag
    
    !** Space for Only one molecule is available in the sorbate structure
    molec=1
    
    !** Do i=1,params%nTemps
    
    !** Set the temperature factor 1/RT 
    !** beta=one/params%temp_k(i)/kcalmole_kb
    
    !** At Present hard coded for only brancehdcoords
    gctype=molecule_getgcmodeltype(params%molecPtr)
    If ( (Trim(toupper(gctype))) /= "BRANCHED" ) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Write now the config integral can be evaluated &
          &only for linear alkanes in the gen-coord type: BRANCHED"
      Write(*,*) "Model Passed here : ",Trim(gctype)
      Stop
    Endif
    
    !      Make sure that there are atleast 3 atoms
    atleast_3 = .False.
    
    
    If (Associated(sorbate%gcoords(molec)%branchedcoords%root%C)) Then
      If (Associated(sorbate%gcoords(1)%branchedcoords%root%C%C)) Then
        atleast_3=.True.
      Endif
    Endif
    
    If (.Not.(atleast_3)) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "There are less than 3 atoms , so no point in calculating&
          &config integral"
      Stop
    Endif
    
    
    !** Allocate and fill the exp(-u) array
    Allocate(expsum(params%nTemps),STAT=error)
    ! Computes ratio , Z(ig)/ Omega
    expsum = zero
    sum_e = zero
    sum_nrg=zero
    
    !** A Monte Carlo Integration is being done here
    !** It needs to be explained with some good notes
    
    !**config_setnmoles(sorbate,1)
    sorbate%nmoles=1
    
    Do k=1, params%ntrials
      
      !** METHOD -1 
      !** Gets a configuration with some biasfactor
      !   Call confinteg_getRandConf(params,sorbate, tk, &
      !   energy, biasfactor)
      
      !** Method-2
      !** Unbiased configurations
      !      Call confinteg_getRandConfNB(params,sorbate, energy, biasfactor)
      !      sum_e = sum_e + Exp(-energy*beta) /biasfactor
      
      
      !** Method-3
      Call confinteg_justgetconf(params,sorbate,biasfactor )
      sorbno=1
      molec=1
      natoms = config_getnatoms(sorbate) 
      Call branchedcoords_toxyz(sorbate%gcoords(molec)%branchedcoords, &
          sorbate%coords(1:natoms,molec)%rp)
      
      fast_flag = .True.   
      Call intramolecular_getint(sorbno,sorbate,fast_flag,potlist)
      energy=potlist(BENDING_INDEX)+potlist(TORSION_INDEX)+ &
          potlist(INTRAPAIR_INDEX)+potlist(INTRACOUL_INDEX)

      !     Might be rqd with Methods.1 and Methods.2 
      !      ljflag=.True.
      !      intrapair=zero
      !      Do atom_num=5,natoms
      !        Call ssbasic_IntraGetinteraction(sorbate, &
      !            ssnoncoulparams(sorbno,sorbno)%basic,params%molecPtr,&
      !            molec,atom_num,ipot,ljflag)
      !        intrapair=intrapair+ipot
      !      End Do
      
      !      sum_nrg=sum_nrg+abs(intrapair/potlist(TOTAL_INDEX))
      
      Do j=1,params%nTemps
        expsum(j)=expsum(j) + Exp(-energy*params%beta(j))/biasfactor
      End Do
    
 !SDEBUG
      outputcount=outputcount+1
      If (outputcount==params%ntrials/2000) Then
        outputcount=0
        Write(*,'(a,i20,e20.6,f15.5)') "CONV_TAG",k,&
            (expsum(1)/k),Log(expsum(1)/k)
      Endif
 !SDEBUG

    End Do
    
    
    !** Write to the output file
    unitno=file_open(params%outfilename)
    Call confinteg_initdisplay(params,5,unitno)
    Write(unitno,*) dashedline
    Write(unitno,*)"***     Temp,K       Z/Omega    -log(Z/Omega)"
    Do j=1,params%nTemps
      Write(unitno,'(f10.3,e16.4,f16.4)') params%temp_k(j),&
          expsum(j)/params%ntrials,-Log(expsum(j)/params%ntrials)
    End Do
    Close(unitno)
    
  
  End Subroutine confinteg_getInteg

  !-----------------------------------------------------------
  ! Gives a random configuration back in sorbate
  ! uses some biasing techniques here , which are based on angle 
  ! and torsion potentials, tk is the temperature
  ! theta sampled without bias, not very efficient
  !-----------------------------------------------------------
  Subroutine confinteg_justgetconf(params,sorbate,biasfactor) 
    Type(CONFINTEG_Params), Intent(inout)  :: params
    Type(AtMolCoords), Intent(inout)   :: sorbate 
    Real(kind=RDbl), Intent(out) :: biasfactor

    Integer              :: natoms, i, j
    Integer              :: atomt1, atomt2, atomt3, atomt4, molec
    Real(kind=RDbl)      :: costheta,ktheta,thetaeq,theta 
    Real(kind=RDbl)      :: costac,costerm   
    Real(kind=RDbl)      :: phi,tors_nrg,exp_bias


    !** builds up the molecule: Asumes there is atleast one molecule
    !** No actual placement just makes ptr%placed=1 
    Call branchedcoords_place(sorbate%gcoords(1)%branchedcoords%root)

    biasfactor=one


    !** Samples the theta for the first triplet of the molecule, atoms:1-2-3
!    costheta=(2*rranf())-one
!    theta=Acos(costheta)
!    Original values
     theta=random_gaussian(params%bias_theta0,params%bias_sig)
     exp_bias=Exp(-params%bias_k*(theta-params%bias_theta0) * &
         (theta-params%bias_theta0)/(2*8.3141*params%bias_Temp))
     biasfactor=2*biasfactor*exp_bias/Sqrt(2*pi)/params%bias_sig
    sorbate%gcoords(1)%branchedcoords%root%C%baPC = theta 
    
    If (.NOT. Associated(sorbate%gcoords(1)%branchedcoords% &
        root%C%C%C)) Then 
       Return 
    End If 
    
    !** Samples the theta for the 2nd triplet of the molecule, atoms:2-3-4
!    costheta = 2*rranf() - one
!    theta = Acos(costheta)
     theta=random_gaussian(params%bias_theta0,params%bias_sig)
     exp_bias=Exp(-params%bias_k*(theta-params%bias_theta0) * &
         (theta-params%bias_theta0)/(2*8.3141*params%bias_Temp))
     biasfactor=2*biasfactor*exp_bias/Sqrt(2*pi)/params%bias_sig 

   !Sample phi  for atoms: 1-2-3-4 from a uniform distribution
    phi = pi - rranf()*twopi
    
    sorbate%gcoords(1)%branchedcoords%root%C%C%baPC = theta
    sorbate%gcoords(1)%branchedcoords%root%C%C%taC  = phi
    If(Associated(sorbate%gcoords(1)%branchedcoords%root%C%C%C%C)) &
        Then 
      !** This should look exaxtly same as what we saw before
      
      Call confinteg_justgetRandNode(params, sorbate, & 
          sorbate%gcoords(1)%branchedcoords%root%C%C%C%C,biasfactor)

    End If
    Return 

  End Subroutine confinteg_justgetconf



  !-----------------------------------------------------------
  ! Calculate Z/omelga ratio for each atom with atom number > 5 
  ! energy and biasfactor are returned 
  ! Here theta is sampled without any bias, Not very efficient
  !-----------------------------------------------------------
  Recursive Subroutine confinteg_justgetRandNode(params,sorbate, ptr,biasfactor)
    Type(CONFINTEG_Params), Intent(inout)  :: params
    Type(AtMolCoords), Intent(inout)   :: sorbate
    Type(NodeType), Pointer           :: ptr                   
    Real(kind=RDbl),Intent(inout)     :: biasfactor

    Integer, Parameter   :: MAX_PARAMS = 6   
    Integer              :: natoms, i, j, sorbtype
    Integer              :: atomt1, atomt2, atomt3, atomt4, molec
    Real(kind=RDbl),Dimension(0:MAX_PARAMS-1) :: cn    
    Real(kind=RDbl)      :: costheta,ktheta,thetaeq,theta,phi 
    Real(kind=RDbl)      :: costac,costerm   
    Real(kind=RDbl)      :: pot
    Real(kind=RDbl)      :: tors_nrg,exp_bias
    Type(VecType)        :: startpos
    Logical              :: fast, ljflag

    Real(kind=RDbl), Dimension(NO_OF_INTRA_POTS) :: potlist

!    costheta = 2*rranf() -one 
!    theta = Acos(costheta)
     theta=random_gaussian(params%bias_theta0,params%bias_sig)
     exp_bias=Exp(-params%bias_k*(theta-params%bias_theta0) * &
         (theta-params%bias_theta0)/(2*8.3141*params%bias_Temp))
     biasfactor=2*biasfactor*exp_bias/Sqrt(2*pi)/params%bias_sig  

   phi = pi - rranf()*twopi
    
    ptr%P%baPC = theta     
    ptr%P%taC  = phi    
    
    If(Associated(ptr%P%L)) Then
      !     Call confinteg_getRandBranch(sorbates, sorbtype, &
      !          ptr%P%L, rti, energy, biasfactor)
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) " Not yet debugged" 
      Stop
    Else If(Associated(ptr%C)) Then
      Call confinteg_justgetRandNode(params, sorbate, ptr%C,biasfactor)
    End If
    
  End Subroutine confinteg_justgetRandNode

!!$
!!$  
!!$  
!!$  !-----------------------------------------------------------
!!$  ! Gives a random configuration back in sorbate
!!$  ! uses some biasing techniques here , which are based on angle 
!!$  ! and torsion potentials, tk is the temperature
!!$  ! theta sampled without bias, not very efficient
!!$  !-----------------------------------------------------------
!!$  !IMPORTANT all angle energy calculations should be done calling
!!$  !the intramolecular modules, not directly from here
!!$
!!$  Subroutine confinteg_getRandConfNB(params,sorbate, energy, biasfactor) 
!!$    Type(CONFINTEG_Params), Intent(inout)  :: params
!!$    Type(AtMolCoords), Intent(inout)   :: sorbate 
!!$    Real(kind=RDbl), Intent(out)   :: energy, biasfactor  
!!$
!!$    Integer              :: natoms, i, j
!!$    Integer              :: atomt1, atomt2, atomt3, atomt4, molec
!!$    Integer              :: MAX_PARAMS=6
!!$    Real(kind=RDbl),Dimension(6)      :: cn
!!$    Real(kind=RDbl)      :: costheta,ktheta,thetaeq,theta 
!!$    Real(kind=RDbl)      :: costac,costerm   
!!$    Real(kind=RDbl)      :: phi,tors_nrg
!!$
!!$    biasfactor=one
!!$
!!$    !** builds up the molecule: Asumes there is atleast one molecule
!!$    !** No actual placement just makes ptr%placed=1 
!!$    Call branchedcoords_place(sorbate%gcoords(1)%branchedcoords%root)
!!$
!!$    !** TEST THIS
!!$    ! rootcoord=>sorbate%gcoords(1)%branchedcoords%root
!!$
!!$    !** Obtains the potential parameters for the first angle
!!$    atomt1 = sorbate%gcoords(1)%branchedcoords%root%atom_type
!!$    atomt2 = sorbate%gcoords(1)%branchedcoords%root%C%atom_type
!!$    atomt3 = sorbate%gcoords(1)%branchedcoords%root%C%C%atom_type
!!$
!!$    !Put some checks here
!!$    !** It should be using atom_list(bblist), not bbparams
!!$    ktheta = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$    thetaeq = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$
!!$    costheta=(2*rranf())-one
!!$    theta=Acos(costheta)
!!$    sorbate%gcoords(1)%branchedcoords%root%C%baPC = theta 
!!$    energy = 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$    
!!$    If (.NOT. Associated(sorbate%gcoords(1)%branchedcoords% &
!!$        root%C%C%C)) Then 
!!$       Return 
!!$    End If 
!!$
!!$    !** Obtains the potential parameters
!!$    atomt1 = sorbate%gcoords(1)%branchedcoords%root%C%atom_type
!!$    atomt2 = sorbate%gcoords(1)%branchedcoords%root%C%C%atom_type
!!$    atomt3 = sorbate%gcoords(1)%branchedcoords%root%C%C%C%atom_type
!!$    ktheta = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$    thetaeq = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$
!!$   !** Obtains the dihedral angle potential
!!$    atomt1 = sorbate%gcoords(1)%branchedcoords%root%atom_type
!!$    atomt2 = sorbate%gcoords(1)%branchedcoords%root%C%atom_type
!!$    atomt3 = sorbate%gcoords(1)%branchedcoords%root%C%C%atom_type
!!$    atomt4 = sorbate%gcoords(1)%branchedcoords%root%C%C%C%atom_type
!!$
!!$
!!$    !** Check torsion potential type
!!$    Do i=0,MAX_PARAMS-1
!!$     cn(i)= params%molecPtr%torsion%torparams(atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)
!!$    End Do
!!$
!!$      costheta = 2*rranf() - one
!!$      theta = Acos(costheta)
!!$      !Sample phi  from a uniform distribution
!!$      phi = pi - rranf()*twopi
!!$      costac  = cos(phi)  
!!$      costerm= 1.0_RDbl
!!$      tors_nrg = 0.0_RDbl 
!!$        Do j=0,MAX_PARAMS-1
!!$          tors_nrg= tors_nrg + cn(j)*costerm
!!$          costerm=costerm*costac
!!$        End Do
!!$    energy = energy + 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)  & 
!!$             + tors_nrg
!!$    biasfactor = one
!!$    sorbate%gcoords(1)%branchedcoords%root%C%C%baPC = theta
!!$    sorbate%gcoords(1)%branchedcoords%root%C%C%taC  = phi
!!$    If(Associated(sorbate%gcoords(1)%branchedcoords%root%C%C%C%C)) &
!!$        Then 
!!$      !** This should look exaxtly same as what we saw before
!!$
!!$     Call confinteg_getRandNodeNB(params, sorbate, & 
!!$         sorbate%gcoords(1)%branchedcoords%root%C%C%C%C, energy, biasfactor)
!!$    End If
!!$
!!$   biasfactor=one
!!$
!!$   Return 
!!$
!!$ End Subroutine confinteg_getRandConfNB
!!$
!!$
!!$
!!$  !-----------------------------------------------------------
!!$  ! Calculate Z/omelga ratio for each atom with atom number > 5 
!!$  ! energy and biasfactor are returned 
!!$  ! Here theta is sampled without any bias, Not very efficient
!!$  !-----------------------------------------------------------
!!$  Recursive Subroutine confinteg_getRandNodeNB(params,sorbate, ptr, &
!!$            energy, biasfactor)
!!$    Type(CONFINTEG_Params), Intent(inout)  :: params
!!$    Type(AtMolCoords), Intent(inout)   :: sorbate
!!$    Type(NodeType), Pointer           :: ptr                   
!!$    Real(kind=RDbl), Intent(inout)   :: energy, biasfactor  
!!$
!!$    Integer, Parameter   :: MAX_PARAMS = 6   
!!$    Integer              :: natoms, i, j, sorbtype
!!$    Integer              :: atomt1, atomt2, atomt3, atomt4, molec
!!$    Real(kind=RDbl),Dimension(0:MAX_PARAMS-1) :: cn    
!!$    Real(kind=RDbl)      :: costheta,ktheta,thetaeq,theta,phi 
!!$    Real(kind=RDbl)      :: costac,costerm   
!!$    Real(kind=RDbl)      :: pot
!!$    Real(kind=RDbl)      :: tors_nrg
!!$    Type(VecType)        :: startpos
!!$    Logical              :: fast, ljflag
!!$
!!$    Real(kind=RDbl), Dimension(NO_OF_INTRA_POTS) :: potlist
!!$
!!$    natoms = config_getnatoms(sorbate) 
!!$
!!$    !** Obtains the potential parameters
!!$    atomt1 = ptr%P%P%atom_type
!!$    atomt2 = ptr%P%atom_type
!!$    atomt3 = ptr%atom_type 
!!$    ktheta = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$    thetaeq = &
!!$        params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$
!!$   !** Obtains the dihedral angle potential
!!$    atomt1 = ptr%P%P%P%atom_type
!!$    atomt2 = ptr%P%P%atom_type
!!$    atomt3 = ptr%P%atom_type
!!$    atomt4 = ptr%atom_type     
!!$    Do i=0,MAX_PARAMS-1
!!$      cn(i) = &
!!$          params%molecPtr%torsion%torparams(atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)
!!$    End Do
!!$
!!$
!!$    costheta = 2*rranf() -one 
!!$    theta = Acos(costheta)
!!$    phi = pi - rranf()*twopi
!!$    energy = energy + 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$    
!!$    costac  = Cos(phi)
!!$    costerm= 1.0_RDbl
!!$    tors_nrg = 0.0_RDbl 
!!$    Do j=0,MAX_PARAMS-1
!!$      tors_nrg= tors_nrg + cn(j)*costerm
!!$      costerm=costerm*costac
!!$    End Do
!!$    
!!$    energy = energy + tors_nrg
!!$    
!!$    ptr%P%baPC = theta     
!!$    ptr%P%taC  = phi    
!!$    Call branchedcoords_toxyz  &
!!$        (sorbate%gcoords(1)%branchedcoords, sorbate%coords(1:natoms,1)%rp)
!!$    
!!$    !** What is ljflag ?
!!$    ljflag = .True.   
!!$
!!$    !** helpless ocrrections!!
!!$    sorbtype = 1
!!$    molec=1
!!$    Call ssbasic_IntraGetinteraction(sorbate, &
!!$        ssnoncoulparams(sorbtype,sorbtype)%basic,params%molecPtr,&
!!$        molec,ptr%atom_num,pot,ljflag)
!!$    energy = energy + pot
!!$    
!!$    biasfactor = one 
!!$    
!!$    If(Associated(ptr%P%L)) Then
!!$      !     Call confinteg_getRandBranch(sorbates, sorbtype, &
!!$      !          ptr%P%L, rti, energy, biasfactor)
!!$      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$      Write(*,*) " Not yet debugged" 
!!$      Stop
!!$    Else If(Associated(ptr%C)) Then
!!$      Call confinteg_getRandNodeNB(params, sorbate, ptr%C, &
!!$          energy, biasfactor)
!!$    End If
!!$    
!!$  End Subroutine confinteg_getRandNodeNB
!!$
!!$
!!$  !-----------------------------------------------------------
!!$  ! Gives a random configuration back in sorbate
!!$  ! uses some biasing techniques here , which are based on angle 
!!$  ! and torsion potentials, tk is the temperature
!!$  !-----------------------------------------------------------
!!$  !IMPORTANT all angle energy calculations should be done calling
!!$  !the intramolecular modules, not directly from here
!!$
!!$  Subroutine confinteg_getRandConf(params,sorbate, tk, energy, biasfactor) 
!!$    Type(CONFINTEG_Params), Intent(inout)  :: params
!!$    Type(AtMolCoords), Intent(inout)   :: sorbate 
!!$    Real(kind=RDbl), Intent(in)    :: tk                  
!!$    Real(kind=RDbl), Intent(out)   :: energy, biasfactor  
!!$
!!$    Integer              :: natoms, i, j, k,index, nx, ny, nz
!!$    Integer              :: icellx,icelly,icellz, ncubelets
!!$    Integer              :: num,low,high,mid,flag 
!!$    Integer              :: atomt1, atomt2, atomt3, atomt4, molec
!!$    Integer              :: MAX_PARAMS=6
!!$    Real(kind=RDbl),Dimension(6)      :: cn
!!$    Real(kind=RDbl)      :: costheta,ktheta,thetaeq,theta 
!!$    Real(kind=RDbl)      :: costac,costerm   
!!$    Real(kind=RDbl)      :: partition,key,wt,prob,rti 
!!$    Real(kind=RDbl)      :: u,energy_1
!!$    Real(kind=RDbl),Dimension((params%angletrials)*2) :: trial_e,ta_e,theta1,phi1,cum_prob
!!$    Type(VecType)        :: startpos
!!$    Logical              :: fast, mapflag
!!$
!!$    Integer :: angletrials
!!$
!!$    !** builds up the molecule: Assumes there is atleast one molecule
!!$    !** No actual placement just makes ptr%placed=1 
!!$    Call branchedcoords_place(sorbate%gcoords(1)%branchedcoords%root)
!!$
!!$    !** TEST THIS
!!$    ! rootcoord=>sorbate%gcoords(1)%branchedcoords%root
!!$
!!$    !** Set the temperature parameter
!!$    rti = 1.0_RDbl/(kcalmole_kb*tk)  
!!$  
!!$
!!$!** This belongs to actual "branched" molecules, neglect for now ?
!!$!    If(Associated(sorbate%gcoords(1)%branchedcoords%root%C%L))Then 
!!$!     Call branchedmoves_BrNodeIntegral(sorbates, sorbtype, &
!!$!          sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C, &
!!$!          rti, energy, biasfactor)
!!$!     Return 
!!$!    End If  
!!$
!!$
!!$
!!$
!!$    !** Obtains the potential parameters for the first angle
!!$    atomt1 = sorbate%gcoords(1)%branchedcoords%root%atom_type
!!$    atomt2 = sorbate%gcoords(1)%branchedcoords%root%C%atom_type
!!$    atomt3 = sorbate%gcoords(1)%branchedcoords%root%C%C%atom_type
!!$
!!$    !Put some checks here
!!$    !** It should be using atom_list(bblist), not bbparams
!!$    ktheta = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$    thetaeq = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$
!!$    !** Add these to some display list
!!$
!!$    
!!$    partition = 0.0_RDbl
!!$    angletrials=params%angletrials
!!$    Do i=1,angletrials
!!$      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/real(angletrials)
!!$      theta = Acos(costheta)
!!$      trial_e(i) = 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$      partition = exp(-trial_e(i)*rti) + partition
!!$    End Do
!!$
!!$    wt = 0.0_RDbl
!!$    Do i=1,angletrials
!!$        wt = wt + exp(-trial_e(i)*rti)/partition
!!$        cum_prob(i) = wt
!!$    End Do
!!$
!!$    low = random_gettrial(1,angletrials,cum_prob)
!!$    costheta = rranf() * 2.0_RDbl/real(angletrials) + &
!!$               -1.0_RDbl + (2.0_RDbl*(low-1))/real(angletrials)
!!$    theta = Acos(costheta)
!!$    energy = 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$
!!$    !** This is the bias by which theta is chosen
!!$    biasfactor =  real(angletrials)*exp(-trial_e(low)*rti)/partition 
!!$    !** Non dimensionalised biasfactor
!!$
!!$    
!!$    sorbate%gcoords(1)%branchedcoords%root%C%baPC = theta 
!!$
!!$    If (.NOT. Associated(sorbate%gcoords(1)%branchedcoords% &
!!$        root%C%C%C)) Then 
!!$       Return 
!!$    End If 
!!$
!!$    !** Obtains the potential parameters
!!$    atomt1 = sorbate%gcoords(1)%branchedcoords%root%C%atom_type
!!$    atomt2 = sorbate%gcoords(1)%branchedcoords%root%C%C%atom_type
!!$    atomt3 = sorbate%gcoords(1)%branchedcoords%root%C%C%C%atom_type
!!$    ktheta = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$    thetaeq = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$
!!$   !** Obtains the dihedral angle potential
!!$    atomt1 = sorbate%gcoords(1)%branchedcoords%root%atom_type
!!$    atomt2 = sorbate%gcoords(1)%branchedcoords%root%C%atom_type
!!$    atomt3 = sorbate%gcoords(1)%branchedcoords%root%C%C%atom_type
!!$    atomt4 = sorbate%gcoords(1)%branchedcoords%root%C%C%C%atom_type
!!$
!!$
!!$    !** Check torsion potential type
!!$    Do i=0,MAX_PARAMS-1
!!$     cn(i)= params%molecPtr%torsion%torparams(atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)
!!$    End Do
!!$
!!$    partition = 0.0_RDbl
!!$    Do i=1,angletrials
!!$      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/real(angletrials)
!!$      theta = Acos(costheta)
!!$      trial_e(i) = 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$
!!$      !Sample phi  from a uniform distribution ??
!!$      phi1(i) = pi - rranf()*twopi
!!$      costac  = cos(phi1(i))  
!!$      costerm= 1.0_RDbl
!!$      ta_e(i) = 0.0_RDbl 
!!$        Do j=0,MAX_PARAMS-1
!!$          ta_e(i)= ta_e(i) + cn(j)*costerm
!!$          costerm=costerm*costac
!!$        End Do
!!$      trial_e(i) = trial_e(i)  + ta_e(i)      
!!$      partition = exp(-(trial_e(i))*rti) + partition
!!$    End Do
!!$
!!$    wt = 0.0_RDbl
!!$    Do i=1,angletrials
!!$      wt = wt + exp(-trial_e(i)*rti)/partition
!!$      cum_prob(i) = wt
!!$    End Do
!!$
!!$    low = random_gettrial(1,angletrials,cum_prob)
!!$    costheta = rranf() * 2.0_RDbl/real(angletrials) + &
!!$               -1.0_RDbl + (2.0_RDbl*(low-1))/real(angletrials)
!!$    theta = Acos(costheta)
!!$    energy = energy + 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)  & 
!!$             + ta_e(low)  
!!$    biasfactor = biasfactor*real(angletrials)*exp(-trial_e(low)*rti)/partition 
!!$    sorbate%gcoords(1)%branchedcoords%root%C%C%baPC = theta
!!$    sorbate%gcoords(1)%branchedcoords%root%C%C%taC  = phi1(low) 
!!$
!!$    !** Again branched molecule here neglect ?
!!$!!
!!$!
!!$!    If(Associated(sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%C%L))
!!$! Then 
!!$!     Call branchedmoves_BrNodeIntegral(sorbates, sorbtype, &
!!$!          sorbates(sorbtype)%gcoords(1)%branchedcoords%root%C%C%C, &
!!$!          rti, energy, biasfactor)
!!$    If(Associated(sorbate%gcoords(1)%branchedcoords%root%C%C%C%C)) &
!!$        Then 
!!$      !** This should look exaxtly same as what we saw before
!!$
!!$     Call confinteg_getRandNode(params, sorbate, & 
!!$         sorbate%gcoords(1)%branchedcoords%root%C%C%C%C, &
!!$         rti, energy, biasfactor)
!!$    End If
!!$
!!$   Return 
!!$ End Subroutine Confinteg_getRandConf
!!$
!!$
!!$
!!$  !-----------------------------------------------------------
!!$  ! Calculate Z/omelga ratio for each atom with atom number > 5 
!!$  ! energy and biasfactor are returned 
!!$  !-----------------------------------------------------------
!!$  Recursive Subroutine confinteg_getRandNode(params,sorbate, ptr, &
!!$            rti, energy, biasfactor)
!!$    Type(CONFINTEG_Params), Intent(inout)  :: params
!!$    Type(AtMolCoords), Intent(inout)   :: sorbate
!!$    Type(NodeType), Pointer           :: ptr                   
!!$    Real(kind=RDbl), Intent(in)    :: rti                   
!!$    Real(kind=RDbl), Intent(inout)   :: energy, biasfactor  
!!$
!!$    Integer, Parameter   :: MAX_PARAMS = 6   
!!$    Integer              :: natoms, i, j, k,index, nx, ny, nz,sorbtype
!!$    Integer              :: icellx,icelly,icellz, ncubelets
!!$    Integer              :: num,low,high,mid,flag ,angletrials
!!$    Integer              :: atomt1, atomt2, atomt3, atomt4, molec
!!$    Real(kind=RDbl),Dimension(0:MAX_PARAMS-1) :: cn    
!!$    Real(kind=RDbl)      :: costheta,ktheta,thetaeq,theta,phi 
!!$    Real(kind=RDbl)      :: costac,costerm   
!!$    Real(kind=RDbl)      :: partition,key,wt,prob  
!!$    Real(kind=RDbl)      :: u,pot,energy_1
!!$    Real(kind=RDbl),Dimension(params%angletrials) :: &
!!$                    trial_e,ta_e,ba_e,theta1,phi1,cum_prob  
!!$    Type(VecType)        :: startpos
!!$    Logical              :: fast, ljflag
!!$
!!$    Real(kind=RDbl), Dimension(NO_OF_INTRA_POTS) :: potlist
!!$
!!$    natoms = config_getnatoms(sorbate) 
!!$
!!$    !** Obtains the potential parameters
!!$    atomt1 = ptr%P%P%atom_type
!!$    atomt2 = ptr%P%atom_type
!!$    atomt3 = ptr%atom_type 
!!$    ktheta = params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%ktheta
!!$    thetaeq = &
!!$        params%molecPtr%bending%bbparams(atomt1,atomt2,atomt3)%hara%thetaeq*degTorad 
!!$
!!$   !** Obtains the dihedral angle potential
!!$    atomt1 = ptr%P%P%P%atom_type
!!$    atomt2 = ptr%P%P%atom_type
!!$    atomt3 = ptr%P%atom_type
!!$    atomt4 = ptr%atom_type     
!!$    Do i=0,MAX_PARAMS-1
!!$      cn(i) = &
!!$          params%molecPtr%torsion%torparams(atomt1,atomt2,atomt3,atomt4)%cosexp%cn(i)
!!$    End Do
!!$    angletrials=params%angletrials
!!$    partition = 0.0_RDbl
!!$    Do i=1,angletrials
!!$      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/Real(angletrials)
!!$      theta1(i) = Acos(costheta)
!!$      ba_e(i) = 0.5_RDbl*ktheta*(theta1(i)-thetaeq)*(theta1(i)-thetaeq)
!!$      phi1(i) = pi - rranf()*twopi       
!!$      costac  = Cos(phi1(i))
!!$      costerm= 1.0_RDbl
!!$      ta_e(i) = 0.0_RDbl 
!!$      Do j=0,MAX_PARAMS-1
!!$        ta_e(i)= ta_e(i) + cn(j)*costerm
!!$        costerm=costerm*costac
!!$      End Do
!!$      trial_e(i) = ba_e(i) + ta_e(i) 
!!$      partition = Exp(-(trial_e(i))*rti) + partition
!!$    End Do
!!$    
!!$    wt = 0.0_RDbl
!!$    Do i=1,angletrials   
!!$      wt = wt + Exp(-trial_e(i)*rti)/partition
!!$      cum_prob(i) = wt
!!$    End Do
!!$    
!!$    low = random_gettrial(1,angletrials,cum_prob)
!!$    costheta = rranf() * 2.0_RDbl/Real(angletrials) + &
!!$        -1.0_RDbl + (2.0_RDbl*(low-1))/Real(angletrials)
!!$    theta = Acos(costheta)
!!$    phi = phi1(low)   
!!$    energy = energy + 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
!!$    energy = energy + ta_e(low)  
!!$    
!!$    ptr%P%baPC = theta     
!!$    ptr%P%taC  = phi    
!!$    Call branchedcoords_toxyz  &
!!$        (sorbate%gcoords(1)%branchedcoords, sorbate%coords(1:natoms,1)%rp)
!!$    
!!$    !** What is ljflag ?
!!$    ljflag = .True.   
!!$
!!$    !** helpless ocrrections!!
!!$    sorbtype = 1
!!$    molec=1
!!$    Call ssbasic_IntraGetinteraction(sorbate, &
!!$        ssnoncoulparams(sorbtype,sorbtype)%basic,params%molecPtr,&
!!$        molec,ptr%atom_num,pot,ljflag)
!!$    energy = energy + pot
!!$    
!!$    biasfactor = biasfactor*Real(angletrials)*Exp(-trial_e(low)*rti)& 
!!$        /partition 
!!$    
!!$    If(Associated(ptr%P%L)) Then
!!$      !     Call confinteg_getRandBranch(sorbates, sorbtype, &
!!$      !          ptr%P%L, rti, energy, biasfactor)
!!$      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
!!$      Write(*,*) " Not yet debugged" 
!!$      Stop
!!$    Else If(Associated(ptr%C)) Then
!!$      Call confinteg_getRandNode(params, sorbate, ptr%C, &
!!$          rti, energy, biasfactor)
!!$    End If
!!$    
!!$  End Subroutine confinteg_getRandNode
!!$  
!!$  !--------------------------------------------
!!$  ! Gets the total no. of simulations
!!$  !--------------------------------------------
!!$  Integer Function confinteg_getnosims(confintegparams)
!!$    Type(CONFINTEG_Params), Intent(in)  :: confintegparams
!!$    confinteg_getnosims = confintegparams%nTemps
!!$  End Function confinteg_getnosims
!!$
  !------------------------------------------------------------------
  ! Display the confintegsimulation parameters for simulation no. "simno"
  ! The display is sent to the optional unit no. "optunit".  "nspc"
  ! is the no. of spaces to leave from the left margin
  !-----------------------------------------------------------------
  Subroutine confinteg_initdisplay(params, nspc, optunit)
    Type(CONFINTEG_Params), Intent(in) :: params

    Integer, Intent(in) :: nspc
    Integer, Optional, Intent(in) :: optunit

    Character(len=nspc) :: spc
    Integer    :: i, unitno, funit
    
    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    !** Write the simulation no.
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "The CONFINTEGSimulation Parameters:"
    Write(unitno, '(2a)') spc//spc,"Description line"//Trim(params%description)
    Write(unitno, '(2a,f16.4,a,f16.4)') spc//spc,"Temperature Range : ",&
        params%temp_k(1)," - ", params%temp_k(params%nTemps)
    Write(unitno, '(2a,i20)') spc//spc, "no of itns : ",  params%ntrials
    Write(unitno, '(2a,i10)') spc//spc, "no of temps : ", params%nTemps
!    Write(unitno, '(2a,i10)') spc//spc, "no of angs : ", params%angletrials
    Write(unitno, '(3a)') spc//spc, "out put : ",   params%outfilename
    Write(unitno,*)"Biases are",params%bias_k,params%bias_theta0, &
        params%bias_Temp

  End Subroutine confinteg_initdisplay

End Module confinteg








