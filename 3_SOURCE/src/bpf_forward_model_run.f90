!*****  subroutine bpf_model_run performs the following tasks
!       delete last set of model output files
!       write new model input files
!       run the forward model(s)
!       read and return the output data
module bayes_pest_forward_model

contains
   subroutine bpf_model_run(errstruc, d_MOD, npar, nobs, &
            &   parval, obsval, miostruc)
    use utilities
    use bayes_pest_control
    use model_input_output
    use error_message
 
    implicit none
!--  Main Data Arrays for OBS and PARS
    type (mio_struc)         :: miostruc
    type (err_failure_struc) :: errstruc

    double precision,          pointer  :: parval(:)
    double precision,          pointer  :: obsval(:)
    type (d_comlin), intent(in) :: d_MOD
    integer,         intent(in) :: npar, nobs
    integer                     :: i, ifail
    character (len=300)         :: instruction
    instruction=' '

    
!-- MIO delete the last set of output files
    if(mio_delete_model_output_files(errstruc,miostruc).ne.0) then
      call utl_bomb_out(errstruc)
    end if   

!-- MIO write the model input files
    if(mio_write_model_input_files(errstruc,miostruc, parval).ne.0) then
      call utl_bomb_out(errstruc)
    end if   


!-- RUN THE MODEL
    call system(d_MOD%com)
    
!-- MIO read the ouput file results and update 
    if(mio_read_model_output_files(errstruc,miostruc, obsval).ne.0) then
      call utl_bomb_out(errstruc)
    end if   

   end subroutine bpf_model_run
   
end module bayes_pest_forward_model