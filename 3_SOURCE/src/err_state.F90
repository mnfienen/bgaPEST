MODULE error_message
!> \brief Error Message Module

IMPLICIT NONE
SAVE
PRIVATE
PUBLIC err_max_err_len, err_max_funcname_len, ERR_FAILURE_STRUC, err_get_error_state, err_get_num_error,&
       err_add_error, err_reset, err_get_message, err_get_failfunc

INTEGER, PARAMETER  :: err_max_err_len = 400
INTEGER, PARAMETER  :: err_max_funcname_len = 100

TYPE ERR_FAILURE_STRUC
  PRIVATE
  INTEGER :: array_size
  INTEGER :: num_errs
  CHARACTER(err_max_err_len), DIMENSION(:), POINTER :: ptr_message => null()
  CHARACTER(err_max_funcname_len), DIMENSION(:), POINTER :: ptr_funcname => null()
END TYPE ERR_FAILURE_STRUC

CONTAINS 

INTEGER FUNCTION err_get_error_state(error_struc)
  TYPE(ERR_FAILURE_STRUC), INTENT(INOUT) :: error_struc
    IF (error_struc%num_errs ==0) THEN
      err_get_error_state = 0
    ELSE
      err_get_error_state = 1
    endif
END FUNCTION err_get_error_state


INTEGER FUNCTION err_get_num_error(error_struc)
  TYPE(ERR_FAILURE_STRUC), INTENT(INOUT) :: error_struc
  err_get_num_error = error_struc%num_errs
 END FUNCTION err_get_num_error
 
SUBROUTINE err_get_message(error_struc, error_index, message)
  TYPE(ERR_FAILURE_STRUC), INTENT(INOUT) :: error_struc
  INTEGER , INTENT(IN) :: error_index
  CHARACTER(LEN=*), INTENT(OUT) :: message
   
  IF( error_index>0 .AND. error_index <= error_struc%num_errs) THEN
    message = error_struc%ptr_message(error_index)
  endif
END SUBROUTINE err_get_message



SUBROUTINE err_get_failfunc(error_struc, error_index, failfunc)
  TYPE(ERR_FAILURE_STRUC), INTENT(INOUT) :: error_struc
  INTEGER , INTENT(IN) :: error_index
  CHARACTER(LEN=*), INTENT(OUT) :: failfunc
   
  IF( error_index>0 .AND. error_index <= error_struc%num_errs) THEN
     failfunc = error_struc%ptr_funcname(error_index)
  endif
END SUBROUTINE err_get_failfunc
  
SUBROUTINE err_add_error(error_struc, err_mesg, funct_name)
  TYPE(ERR_FAILURE_STRUC), INTENT(INOUT) :: error_struc
  CHARACTER(len=*), INTENT(IN) :: err_mesg
  CHARACTER(len=*), INTENT(IN) :: funct_name
  IF ( error_struc%num_errs+1 > error_struc%array_size .OR. error_struc%array_size == 0) THEN
    call err_increase_array_size(error_struc)
  endif
  error_struc%num_errs = error_struc%num_errs + 1
  error_struc%ptr_message(error_struc%num_errs) = err_mesg
  error_struc%ptr_funcname(error_struc%num_errs) = funct_name
END SUBROUTINE err_add_error

SUBROUTINE err_reset(error_struc)
  TYPE(ERR_FAILURE_STRUC), INTENT(INOUT) :: error_struc
  INTEGER :: i
  
  if(associated(error_struc%ptr_message))error_struc%ptr_message = ""
  error_struc%num_errs = 0
END SUBROUTINE err_reset

SUBROUTINE err_increase_array_size(error_struc)
  TYPE(ERR_FAILURE_STRUC), INTENT(INOUT) :: error_struc
  CHARACTER(err_max_err_len), DIMENSION(:), POINTER :: ptr_message_old
  CHARACTER(err_max_err_len), DIMENSION(:), POINTER :: ptr_message_new
  CHARACTER(err_max_funcname_len), DIMENSION(:), POINTER :: ptr_funcname_old
  CHARACTER(err_max_funcname_len), DIMENSION(:), POINTER :: ptr_funcname_new
  INTEGER new_array_size
  INTEGER alloc_error1, alloc_error2
  INTEGER i
  
  new_array_size=error_struc%array_size+1
  ptr_message_old => error_struc%ptr_message
  ptr_funcname_old => error_struc%ptr_funcname
  ALLOCATE(ptr_message_new(new_array_size), STAT=alloc_error1)
  ALLOCATE(ptr_funcname_new(new_array_size), STAT=alloc_error2)
  IF ( alloc_error1 == 0 .and. alloc_error2==0 ) THEN
    error_struc%ptr_message=> ptr_message_new
    error_struc%ptr_funcname=> ptr_funcname_new
    DO i=1, error_struc%array_size
      error_struc%ptr_message(i) = ptr_message_old(i)
      error_struc%ptr_funcname(i) = ptr_funcname_old(i)
    enddo
    IF (error_struc%array_size > 0) THEN
      DEALLOCATE(ptr_message_old)
      DEALLOCATE(ptr_funcname_old)
    endif
    error_struc%array_size = new_array_size
  endif
 END SUBROUTINE err_increase_array_size
  

END MODULE error_message