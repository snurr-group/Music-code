!------------------------------------------------------
! This module contains the data structure and code
! to read in and work with the frame library            
!-----------------------------------------------------
Module frame    

  Use defaults, Only: RDbl, RSgl, pi, twopi, one, dashedline, &
      kcalmole_kb, degTorad 
  Use utils, Only:  
  Use file, Only:
  Use random, Only: rranf 
  Use molecule, Only: MolecularParams 
  Use molecules, Only: molecules_getpointer   
  Use connects, Only: connects_nconnected
  Use branchedcoords, Only: NodeType 

  Implicit None
  Save

  Private
  Public :: FrameInfo, Frame_Library_Params, frame_CountLib, & 
            frame_init, frame_getnframe, frame_getbiasfactor, & 
            frame_display  

  Type FrameInfo
    Real(kind=RDbl)             :: baPC
    Real(kind=RDbl)             :: baPL
    Real(kind=RDbl)             :: baCL
    Real(kind=RDbl)             :: baCR
    Real(kind=RDbl)             :: baPR
    Real(kind=RDbl)             :: u
    Real(kind=RDbl)             :: biasfactor
  End Type FrameInfo

  Type Frame_Library_Params
    Integer                        :: atom_num 
    Integer                        :: num_ba
    Integer                        :: num_frame
    Real(kind=RDbl)                :: rti
    Real(kind=RDbl)                :: ratio  ! omega/Z
    Real(kind=RSgl)                :: cum_wt
    Type(FrameInfo), Dimension(:), Pointer :: frame
  End Type Frame_Library_Params
 
  Integer, Parameter :: LIB_NO_ENTRIES =10000

Contains
  !------------------------------------------------------------------------
  ! Counts the frame library 'frame'.  
  !------------------------------------------------------------------------
  Recursive Subroutine frame_CountLib(ptr,lib_num,atom_num)       
    Type(NodeType), Pointer  :: ptr 
    Integer, Intent(inout)   :: lib_num 
    Integer, Dimension(:), Intent(inout)  :: atom_num 
        
    If (Associated(ptr%L)) Then 
     lib_num = lib_num + 1 
     atom_num(lib_num) = ptr%atom_num  
     Call  frame_CountLib(ptr%L,lib_num,atom_num)    
    End If  

    If (Associated(ptr%C)) Then  
       Call  frame_CountLib(ptr%C,lib_num,atom_num)    
    End If

  End Subroutine frame_CountLib 


  !------------------------------------------------------------------------
  ! Initialize the frame library 'frame'.  The bias can be used
  ! to preferentially sample the low energy regions of the fragment   
  !------------------------------------------------------------------------
  Subroutine frame_init(ptr,frame_lib,sorbtype,tk)
    Type(NodeType), Pointer   :: ptr 
    Type(Frame_Library_Params) ::  frame_lib       
    Real(kind=RDbl), Intent(in)  :: tk   
    Integer, Intent(in)  :: sorbtype 

    Integer    :: error, num, i
    Integer    :: atomt1, atomt2, atomt3, index 
    Real(kind=RDbl)  :: partition, rti, ratio 
    Real(kind=RDbl)  :: boltz, orig_e, e, delta, old_angle 
    Real(kind=RDbl), Dimension(4)  :: ktheta, thetaeq  
    Real(kind=RDbl), Dimension(4)  :: orig_angle
    Real(kind=RDbl), Dimension(LIB_NO_ENTRIES) :: u, baPL,baPC,baCL 
    Type(MolecularParams), Pointer :: molecule  

    If (.NOT. Associated(ptr%L)) Then 
    Write(0,'(2a,i4,a,i5,i10)') __FILE__,": ",__LINE__, &
          " ptr%L nonexistent ", ptr%atom_num 
     Stop
    End If  
    
    rti = 1.0_RDbl/(kcalmole_kb*tk)  
    
    !** Initialize the force field  
    Call molecules_getpointer(sorbtype,molecule)   
    !**The number of bond angles = the number of connections  
    frame_lib%num_ba = connects_nconnected(molecule%connect(ptr%atom_num)) 

    If( .NOT. Associated(ptr%R)) Then 
      atomt1 = ptr%P%atom_type 
      atomt2 = ptr%atom_type 
      atomt3 = ptr%C%atom_type   
!LC      ktheta(1) = molecule%bending%bbparams(atomt1,atomt2,atomt3)& 
!LC                  %hara%ktheta  
!LC      thetaeq(1)= molecule%bending%bbparams(atomt1,atomt2,atomt3)&
!LC                  %hara%thetaeq*degTorad 

      atomt1 = ptr%P%atom_type 
      atomt2 = ptr%atom_type 
      atomt3 = ptr%L%atom_type   
!LC      ktheta(2) = molecule%bending%bbparams(atomt1,atomt2,atomt3)& 
!LC                  %hara%ktheta  
!LC      thetaeq(2)= molecule%bending%bbparams(atomt1,atomt2,atomt3)&
!LC                  %hara%thetaeq*degTorad 

      atomt1 = ptr%C%atom_type 
      atomt2 = ptr%atom_type 
      atomt3 = ptr%L%atom_type   
!LC      ktheta(3) = molecule%bending%bbparams(atomt1,atomt2,atomt3)& 
!LC                  %hara%ktheta  
!LC      thetaeq(3)= molecule%bending%bbparams(atomt1,atomt2,atomt3)&
!LC                  %hara%thetaeq*degTorad 

    !** the case for the tetrahedral carbon, the conversion 
    !** from the bond angle to the dihedral hasn't been done 
    Else  
      atomt1 = ptr%P%atom_type 
      atomt2 = ptr%atom_type 
      atomt3 = ptr%L%atom_type   
!LC      ktheta(1) = molecule%bending%bbparams(atomt1,atomt2,atomt3)& 
!LC                  %hara%ktheta  
!LC      thetaeq(1)= molecule%bending%bbparams(atomt1,atomt2,atomt3)&
!LC                  %hara%thetaeq*degTorad 

      atomt1 = ptr%P%atom_type 
      atomt2 = ptr%atom_type 
      atomt3 = ptr%R%atom_type   
!LC      ktheta(2) = molecule%bending%bbparams(atomt1,atomt2,atomt3)& 
!LC                  %hara%ktheta  
!LC      thetaeq(2)= molecule%bending%bbparams(atomt1,atomt2,atomt3)&
!LC                  %hara%thetaeq*degTorad 
    
      atomt1 = ptr%C%atom_type 
      atomt2 = ptr%atom_type 
      atomt3 = ptr%L%atom_type   
!LC      ktheta(3) = molecule%bending%bbparams(atomt1,atomt2,atomt3)& 
!LC                  %hara%ktheta  
!LC      thetaeq(3)= molecule%bending%bbparams(atomt1,atomt2,atomt3)&
!LC                  %hara%thetaeq*degTorad 
  
      atomt1 = ptr%C%atom_type 
      atomt2 = ptr%atom_type 
      atomt3 = ptr%R%atom_type   
!LC      ktheta(4) = molecule%bending%bbparams(atomt1,atomt2,atomt3)& 
!LC                  %hara%ktheta  
!LC      thetaeq(4)= molecule%bending%bbparams(atomt1,atomt2,atomt3)&
!LC                  %hara%thetaeq*degTorad 

    End If  
           
    !** Sampling the bond angles using unbiased MC  
    orig_angle(1) = rranf()*pi      
    orig_angle(2) = rranf()*pi      
    orig_angle(3) = rranf()*pi      
    orig_e  = 0.0_RDbl            
    Do i=1,frame_lib%num_ba       
      orig_e = orig_e + 0.5_RDbl*ktheta(i)* &
               (orig_angle(i)-thetaeq(i))*(orig_angle(i)-thetaeq(i))
    End Do   

    !** equilibrium stage 
    Do i=1,LIB_NO_ENTRIES
    ! ** randomly picks up one among the three 
      index = INT(rranf()*3) + 1  
      old_angle = orig_angle(index) 
      delta = rranf()*1.0_RDbl * (-1)**i 
      orig_angle(index) = orig_angle(index) + delta 

      If ( orig_angle(3) <= (2.0_RDbl*pi- orig_angle(1) - orig_angle(2)) &
          .AND. orig_angle(3)  >= ABS( orig_angle(1) - orig_angle(2)) )   Then
      e = orig_e + 0.5_RDbl * ktheta(index)* & 
          (orig_angle(index) - thetaeq(index)) * & 
          (orig_angle(index) - thetaeq(index)) & 
           - 0.5_RDbl *ktheta(index)*(old_angle - thetaeq(index)) * & 
             (old_angle - thetaeq(index))  

        boltz = exp(-(e-orig_e) *rti) 
          If ( boltz > rranf()) Then 
            orig_e = e 
          Else 
            orig_angle(index) = old_angle           
          End If 
      Else 
         orig_angle(index) = old_angle           
      End If 
    End Do

   Write(0,'(2a,i4,a,i5,3f13.3)') __FILE__,": ",__LINE__, &
   " finish equi "  
 
    !** production stage
    num = 0 
    Do i=1,LIB_NO_ENTRIES*3 
    ! ** randomly picks up one among the three
      index = INT(rranf()*3) + 1
      old_angle = orig_angle(index)
      delta = rranf()*1.0_RDbl * (-1)**i
      orig_angle(index) = orig_angle(index) + delta

      If ( orig_angle(3) <= (2.0_RDbl*pi- orig_angle(1) - orig_angle(2)) &
          .AND. orig_angle(3)  >= ABS( orig_angle(1) - orig_angle(2)) )   Then
      e = orig_e + 0.5_RDbl * ktheta(index)* &
          (orig_angle(index) - thetaeq(index)) * &
          (orig_angle(index) - thetaeq(index)) &
           - 0.5_RDbl *ktheta(index)*(old_angle - thetaeq(index)) * &
             (old_angle - thetaeq(index))

        boltz = exp(-(e-orig_e) *rti)
          If ( boltz > rranf()) Then
           orig_e = e
           num = num + 1  
           u(num)= e  
           baPC(num) = orig_angle(1)       
           baPL(num) = orig_angle(2)       
           baCL(num) = orig_angle(3)       
           partition = partition + exp(-e*rti) 
           ratio     = ratio     + exp(e*rti) 
          Else
           orig_angle(index) = old_angle
          End If
      Else
         orig_angle(index) = old_angle
      End If
    End Do

    Allocate(frame_lib%frame(num), stat=error)
      If (error /= 0) Then
        Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
            "Could not allocate memory for ptr%frame_lib(framenum)"
        Stop
      End If  
    frame_lib%num_frame =  num   
    ratio = ratio/real(frame_lib%num_frame)  
    frame_lib%ratio = ratio 

    Write(0,'(2a,i4,a,i5,3f13.3)') __FILE__,": ",__LINE__, &
   "  frame_num    ", num, ratio, log(ratio) 
   
    Do i=1,frame_lib%num_frame 
      frame_lib%frame(i)%u  = u(i) 
      frame_lib%frame(i)%baPC  = baPC(i) 
      frame_lib%frame(i)%baPL  = baPL(i) 
      frame_lib%frame(i)%baCL  = baCL(i) 
      frame_lib%frame(i)%biasfactor = 1.0_RDbl/(ratio*exp(-u(i)*rti)) 
    End Do
  End Subroutine frame_init

  !----------------------------------------------------
  ! Gets the total no. of frames in the library  
  !----------------------------------------------------
  Integer Function frame_getnframe(frame_lib)
    Type(Frame_Library_Params), Pointer :: frame_lib 
    frame_getnframe = frame_lib%num_frame 
  End Function frame_getnframe   

  !----------------------------------------------------
  ! Get the biasfactor of the frame indexed "index"
  !----------------------------------------------------
  Real(kind=RDbl) Function frame_getbiasfactor(frame_lib,index)
    Type(Frame_Library_Params), Pointer :: frame_lib 
    Integer, Intent(in)  :: index
    
    frame_getbiasfactor = frame_lib%frame(index)%biasfactor    

  End Function frame_getbiasfactor 

  !-------------------------------------------------------
  ! Display the frame library information
  !-------------------------------------------------------
  Subroutine frame_display(frame_lib, nspc, optunit)
    Type(Frame_Library_Params), Pointer :: frame_lib 
    Integer, Intent(in)  :: nspc ! No. of spaces from the left column
    Integer, Optional, Intent(in)    :: optunit

    Integer   :: i, unitno
    Character(len=nspc) :: spc

    spc = ''
    Do i=1, nspc
      spc = spc//' '
    End Do

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "The Frame Library Section:"
    Write(unitno, '(a,2x,a,g14.8)') spc, "Cumulative Wt  : ", frame_lib%cum_wt
    Write(unitno, '(a,2x,a,i8)') spc, "No. of tabulated frame:  ", &
        frame_lib%num_frame 

  End Subroutine frame_display

End Module frame 
