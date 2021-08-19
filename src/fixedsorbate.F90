!------------------------------------------------------------------
!This module deals with maintaining a fixed position for the 
!molecules of one of the sorbates.   
!------------------------------------------------------------------

Module fixedsorbate
  Use utils, Only: filesrchstr, cleanstring
  Use file, Only: file_open
  Use molecules, Only: molecules_gettype
  Use defaults, Only: strLen, lstrLen 
  Implicit None
  Save

  Private
  Public ::fixedsorbate_init, Fixed_sorbate_is_on, Fixed_sorbate_number, &
    fixedsorbate_check 

  Logical ::  Fixed_sorbate_is_on
  Integer ::  Fixed_sorbate_number 

  Character(len=lstrLen), Parameter :: mdfixedsorbate_tag = &
      "Fixed Sorbate MD"
Contains

  !-------------------------------------------------------------------
  !Determines if there is a sorbate to be fixed
  !------------------------------------------------------------------
  Subroutine fixedsorbate_init(ctrl_filename)       
    Character(*), Intent(In) :: ctrl_filename
    Integer                  :: UseMDfixedsorbate, unitno
    Character(len=255)       :: text

    ! Get the unitno for the control file
    unitno = file_open(ctrl_filename)

    Fixed_sorbate_is_on=.False.

    !Read the control file to see if there is a "Fixed Sorbate MD" section 
    useMDfixedsorbate = filesrchstr(unitno,mdfixedsorbate_tag,text,.True.)
    If (useMDfixedsorbate /= 0) Then
      Fixed_sorbate_is_on=.True. 
    Else
      Return
    End If

    !Read the rest of the ctrlfile
    Read(unitno,'(a)') text
    text=cleanstring(text)
    Fixed_sorbate_number =molecules_gettype(text)    
    Write(*,*)'The fixed sorbate is,',text

  End Subroutine fixedsorbate_init

!---------------------------------------------------------------------------
!Returns true if the sorbate is to be fixed, and false if it is to be moved
!--------------------------------------------------------------------------

  Logical Function fixedsorbate_check(sorbnumber)
    Integer, Intent(In) :: sorbnumber

    If (sorbnumber == Fixed_sorbate_number) Then
      fixedsorbate_check = .True.
    Else
      fixedsorbate_check = .False.
    End If
    
  End Function fixedsorbate_check


End Module fixedsorbate
