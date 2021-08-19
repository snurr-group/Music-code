!------------------------------------------------------
! This module contains the data structure and code
! to read in and work with the angle library  required for cb-gcmc
!-----------------------------------------------------
Module angle    

  Use branchedcoords, Only: NodeType    
  Use defaults, Only: RDbl, zero, pi, twopi, one, dashedline, & 
      kcalmole_kb, strLen
  Use file, Only:
  Use harmonica, Only: HarmonicaModel, HarmonicaString
  Use random, Only: rranf
  Use utils, Only: allocErrDisplay, toupper
  Implicit None
  Save

  Private
  Public :: AngleInfo, Angle_Library_Params, angle_init, & 
            angle_CountLib, angle_getbiasindex, angle_display, &
            angle_indexFromVal 

  Type AngleInfo
    Real(kind=RDbl)             :: ba  
    Real(kind=RDbl)             :: u
    ! cumulative probability , probability (both normalized)
    Real(kind=RDbl)             :: cum_wt ,wt         
  End Type AngleInfo

  Type Angle_Library_Params
    Integer                     :: atom_num 
    Integer                     :: num_ba
    Integer                     :: init_index 
    Real(kind=RDbl)             :: partition   
    Type(AngleInfo), Dimension(:), Pointer :: angle  
    Type(HarmonicaModel)        :: potparams
  End Type Angle_Library_Params

  Integer, Parameter :: LIB_NO_ENTRIES =10000
 

Contains
  !------------------------------------------------------------------------
  ! Counts the angle library 'angle'.  
  !------------------------------------------------------------------------
  Recursive Subroutine angle_CountLib(ptr,lib_num,atom_num)       
    Type(NodeType), Pointer  :: ptr 
    Integer, Intent(inout)   :: lib_num 
    Integer, Dimension(:), Intent(inout)  :: atom_num 
        
    If (Associated(ptr%P)) Then 
     If (Associated(ptr%C)) Then 
       lib_num = lib_num + 1 
       atom_num(lib_num) = ptr%atom_num  
       Call  angle_CountLib(ptr%C,lib_num,atom_num)    
     End If  
    End If  

  End Subroutine angle_CountLib 


  !------------------------------------------------------------------------
  ! Initialize the angle library 'angle'.  The bias can be used 
  ! to preferentially sample the low energy regions of the bond angles
  ! - angles should be passed here in radian NOT degrees
  !------------------------------------------------------------------------
  Subroutine angle_init(ptr,angle_lib,tk,paramsArr,model)
    Type(NodeType), Pointer   :: ptr 
    Type(Angle_Library_Params),Intent(inout) ::  angle_lib       
    Real(kind=RDbl), Intent(in)  :: tk   
    Real(kind=RDbl), Dimension(2), Intent(in)  ::  paramsArr    
    Character(len=strLen), Intent(in) :: model
    Integer    :: error, i
    Integer    :: atom1, atom2, atom3
    Real(kind=RDbl)  :: partition, cum_wt, wt, rti
    Real(kind=RDbl)  :: trial_e, costheta, theta, ktheta, thetaeq  


    If (.NOT. Associated(ptr%C)) Then 
    Write(0,'(2a,i4,a,i5,i10)') __FILE__,": ",__LINE__, &
          " ptr%C nonexistent ", ptr%atom_num 
     Stop
    End If  

    rti = 1.0_RDbl/(kcalmole_kb*tk)  

    !**The number of bond angles = the number of connections  
    angle_lib%atom_num = ptr%atom_num  
    angle_lib%num_ba = LIB_NO_ENTRIES

    atom1 = ptr%P%atom_num
    atom2 = ptr%atom_num     
    atom3 = ptr%C%atom_num
#ifdef DEBUG
    Write(*,'(a,3i3)') "Initialising bbending-library between atoms :", &
        atom1,atom2,atom3
#endif
    If (Trim(toupper(model))==HarmonicaString) Then
      !** Angle in Radians
      angle_lib%potparams%ktheta=paramsArr(1)
      angle_lib%potparams%thetaeq=paramsArr(2)
    Else
      Write(*,*) "The model name is : "//Trim(model)
      Write(*,*) "Might need to change the code here "
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    
    
    !** note the converion from deg to rad
    ktheta=    angle_lib%potparams%ktheta
    thetaeq=    angle_lib%potparams%thetaeq
#ifdef DEBUG
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) "Making angle librray here : " 
    Write(*,'(a,2f12.3)') "   Params used : (k-kcal/rad^2, theta-rad)", &
        ktheta, thetaeq
#endif
    Allocate(angle_lib%angle(LIB_NO_ENTRIES), stat=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
    partition = 0.0_RDbl 
    Do i=1,LIB_NO_ENTRIES
      costheta = -1.0_RDbl + (2.0_RDbl*i-1.0_RDbl)/real(LIB_NO_ENTRIES)
      theta = Acos(costheta)
      trial_e = 0.5_RDbl*ktheta*(theta-thetaeq)*(theta-thetaeq)
      angle_lib%angle(i)%ba = theta  
      angle_lib%angle(i)%u  = trial_e  
      partition = exp(-trial_e*rti) + partition
    End Do
    angle_lib%partition   = partition 
    cum_wt = 0.0_RDbl
    Do i=1,LIB_NO_ENTRIES
        trial_e = angle_lib%angle(i)%u
        wt=Exp(-trial_e*rti)/partition
        cum_wt = cum_wt + wt
        angle_lib%angle(i)%wt = wt  
        angle_lib%angle(i)%cum_wt = cum_wt  
    End Do

  End Subroutine angle_init


  !-----------------------------------------------
  ! given theta gives back its bias
  !-----------------------------------------------
  Real(kind=RDbl) Function angle_biasFromVal(angle_lib,theta) 
    Type(Angle_Library_Params), Pointer :: angle_lib 
    Real(kind=RDbl),Intent(in) :: theta
    Integer :: index

    index=angle_indexFromVal(angle_lib,theta) 
    angle_biasFromVal=angle_lib%angle(index)%wt

    End Function angle_biasFromVal  

  !-----------------------------------------------
  ! given theta gives back its index number
  !-----------------------------------------------
  Integer Function angle_indexfromVal(angle_lib,theta) 
    Type(Angle_Library_Params), Pointer :: angle_lib 
    Real(kind=RDbl),Intent(in) :: theta
    Real(kind=RDbl) :: dcostheta

      dcostheta=(2.0_RDbl)/angle_lib%num_ba
      angle_indexFromVal= ( (cos(theta)+one)/dcostheta )+1

      If (angle_indexFromVal>angle_lib%num_ba) Then
        angle_indexFromVal=angle_lib%num_ba
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Elseif(angle_indexFromVal<1) Then
        angle_indexFromVal=1
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Endif
      
    End Function angle_indexFromVal  
    


  !-----------------------------------------------
  ! Picks an angle according to its bias weight
  ! It returns the index of the angle   
  !alsoo gives back theta and its probability
  !-----------------------------------------------
  Integer Function angle_getbiasindex(angle_lib,theta,energy,prob) 
    Type(Angle_Library_Params), Pointer :: angle_lib 
    Real(kind=RDbl),Intent(out) :: theta,prob,energy

    Integer   :: low, mid, high 
    Real(kind=RDbl)  :: key ,cum_wt
!!    Real(kind=RDbl),Dimension(angle_lib%num_ba)  :: cum_prob

    low = 1
    high = angle_lib%num_ba 

!!    Do i=1,high 
!!      cum_prob(i) = angle_lib%angle(i)%cum_wt 
!!    End Do
!! 
!!    angle_getbiasindex = random_gettrial(low,high,cum_prob)  
!**I think the way it was coded was inefficient, 
!**trying to speed it up here-shaji


    !** rranf returns RDbl, but key(RSgl) has only first 1/2 of the sign.digits
    key  = rranf()

    !** Do a binary search to find the first index with a cumulative
    !** probability higher than key
    Do
      If (low >= high) Exit
      mid = (low + high)/2
      cum_wt = angle_lib%angle(mid)%cum_wt
      If (key < cum_wt) Then
        high = mid - 1
      Else If (key > cum_wt) Then
        low = mid + 1
        mid = low
      Else If (key == cum_wt) Then
        Exit
      Endif
    Enddo
    
    !** low is equal to high
    theta=angle_lib%angle(mid)%ba
    energy=angle_lib%angle(mid)%u
    prob=angle_lib%angle(mid)%wt
    angle_getbiasindex = mid

  End Function angle_getbiasindex  

  !-------------------------------------------------------
  Subroutine angle_display(angle_lib, nspc, optunit)
    Type(Angle_Library_Params), Pointer :: angle_lib 
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
    Write(unitno, '(2a)') spc, "The Angle Library Section:"

  End Subroutine angle_display

End Module angle 
