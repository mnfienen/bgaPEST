module bayes_pest_mio_setup


contains
    subroutine bpm_setup_mio(errstruc, cv_MIO, d_MIO, npargp, nobsgp, npar, nobs, miostruc)
!********  subroutine bp_process_mio
    ! initialize and 
    use bayes_pest_control
    use model_input_output
    use error_message
    use utilities
    implicit none

 
    type (mio_struc)         :: miostruc
    type (err_failure_struc) :: errstruc
    type (cv_minout),intent(in) :: cv_MIO
    type (d_minout), intent(in) :: d_MIO
    integer,         intent(in) :: npargp, nobsgp, npar, nobs
    integer, parameter          :: TPLFILE=1
    integer, parameter          :: INPUTFILE=2
    integer, parameter          :: INSFILE=3
    integer, parameter          :: OUTPUTFILE=4
    integer                     :: ierr, i, j, k, n1, n2, n3, n4
    integer                     :: tpl_status, ins_status


!-- MIO initialization of structures, pars, obs, input, and output files
!       also, set precision and point.

    if(mio_initialise_model_input(errstruc,miostruc,cv_MIO%ntplfle).ne.0) then
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    endif
    if(mio_initialise_model_output(errstruc,miostruc,cv_MIO%ninsfle).ne.0) then
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    endif
    if(mio_set_number_precision(errstruc,miostruc,'single').ne.0) then
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    endif
    if(mio_set_number_decpoint(errstruc,miostruc,'point').ne.0) then
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    endif

!-- MIO  filename storage initialization
    do i = 1,cv_MIO%ntplfle
      if(mio_put_file_model_input(errstruc,miostruc,i,d_MIO%tpl(i),d_MIO%infle(i)).ne.0) then
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
      endif
    enddo
    do i = 1,cv_MIO%ninsfle
      if(mio_put_file_model_output(errstruc,miostruc,i,d_MIO%ins(i),d_MIO%outfle(i)).ne.0) then
        call utl_bomb_out(errstruc)
        n1=mio_finalise(errstruc,miostruc)
        stop 
      endif    
    enddo

!-- MIO - report dimensions of files
      if(mio_get_dimensions(errstruc,miostruc,n1,n2,n3,n4).ne.0) then
        call utl_bomb_out(errstruc)
      endif  
      write (6,201) n1
201     format(' Number of template files        = ', i6)
      write (6,202) n2
202     format(' Number of instruction files     = ', i6)
      write (6,203) n3
203     format(' Number of parameters            = ', i6) 
      write (6,204) n4
204     format(' Number of observations          = ', i6)
   
!-- MIO Check status of TPL and INS files :-> 1==good, 0==bad
! DONE BEFORE AND AFTER PROCESSING TEMPLATE AND INSTRUCTION FILES?
! MNF DEBUG  DEAL WITH THIS 
    write(6,*)
    if(mio_get_temp_ins_status(errstruc,miostruc,tpl_status,ins_status).ne.0) then
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    endif  
    !MNF DEBUG --> 6/6/08  need to handle status exceptions with an error
    
!--  MIO - Check both  parameters and observations
    if(mio_parameter_check(errstruc,miostruc).ne.0) then
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    endif  
    if(mio_observation_check(errstruc,miostruc).ne.0) then
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    endif      
    
!-- MIO  process template files
    if(mio_process_template_files(errstruc,miostruc).ne.0) then
      call utl_bomb_out(errstruc)
    endif    
!-- MIO  process instruction files 
    if(mio_store_instruction_set(errstruc,miostruc).ne.0) then
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    endif   
    
end subroutine bpm_setup_mio
 

end module bayes_pest_mio_setup    
