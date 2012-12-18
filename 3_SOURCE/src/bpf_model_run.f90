!*****  subroutine bpf_model_run performs the following tasks
!       delete last set of model output files
!       write new model input files
!       run the forward model(s)
!       read and return the output data
module bayes_pest_model

contains
   subroutine bpf_model_run(errstruc, d_MOD, cv_PAR, d_PAR, cv_OBS, cv_A, d_OBS, d_A,forward_flag, miostruc)
    use utilities
    use bayes_pest_control
    use model_input_output
    use error_message
    use extern_derivs
    use jacread

 
    implicit none
!--  Main Data Arrays for OBS and PARS
    type (mio_struc)         :: miostruc
    type (err_failure_struc) :: errstruc
    type (d_comlin)                     :: d_MOD
    type (d_param)                      :: d_PAR
    type (cv_param)                     :: cv_PAR
    type (cv_observ)                    :: cv_OBS
    type (cv_algorithmic)               :: cv_A
    type(d_algorithmic),  intent(inout) :: d_A
    type (d_observ)                     :: d_OBS
    integer, intent(in)                 :: forward_flag ! 0, 1, 2, 3, or 4
                                        ! 0 is forward run
                                        ! 1 is external PEST-style Jacobian
                                        ! 2 is dercom alternative Jacobian
                                        ! 3 is forward run for linesearch using d_PAR%pars_lns
                                        ! 4 is external and parallel Jacobian using Condor
    integer                             :: i, ifail
    character (len=100)                 :: adjfle

    !-- MIO delete the last set of output files 
        if(mio_delete_model_output_files(errstruc,miostruc).ne.0) then
          call utl_bomb_out(errstruc)
        endif   
 
    !-- MIO write the model input files
        select case (forward_flag)
            case (3)
                if(mio_write_model_input_files(errstruc,miostruc, d_PAR%pars_lns).ne.0) then
                    call utl_bomb_out(errstruc)
                endif 
            case default
                if(mio_write_model_input_files(errstruc,miostruc, d_PAR%pars).ne.0) then
                  call utl_bomb_out(errstruc)
                endif   
        end select

!-- RUN THE MODEL IN THE MODE DESIRED
    select case (forward_flag)
        case (0) ! single forward run
            call system(d_MOD%com)
            !-- MIO read the ouput file results and update 
            if(mio_read_model_output_files(errstruc,miostruc, d_OBS%h).ne.0) then
              call utl_bomb_out(errstruc)
            endif 
        case (1) ! external PEST-style Jacobian
         !MD Done at the very beginning !
         !MD!   call system(d_MOD%com) ! run the model once, forward, to have outputs with current parameters 
         !MD!   !-- MIO read the ouput file results and update 
         !MD!   if(mio_read_model_output_files(errstruc,miostruc, d_OBS%h).ne.0) then
         !MD!     call utl_bomb_out(errstruc) 
         !MD!   endif 
            !-- create PEST input files and run PEST          
            call bxd_write_param_file(cv_PAR,d_PAR) ! write the parameter file
            call system('pst_generator.exe')        ! create the necessary PEST control file
            call system('run_pest_scratch.bat')     ! run PEST externally for derivatives
            call readJCO('scratch.jco', d_A)
        case (2) ! dercom alternative Jacobian
         !MD Done at the very beginning !
         !MD!   call system(d_MOD%com) !Rum forward model once to obtain the current outputs, before running external derivative model
         !MD!   !-- MIO read the ouput file results and update 
         !MD!   if(mio_read_model_output_files(errstruc,miostruc, d_OBS%h).ne.0) then
         !MD!     call utl_bomb_out(errstruc)
         !MD!   endif
            
            call system(d_MOD%dercom) !Run the external derivative model 
            select case (cv_A%jacobian_format)
              case ('binary')
                call readJCO(cv_A%jacfle, d_A)
              case ('ascii')
                call readJAC(cv_A%jacfle, d_A)
            end select
             
        case (3) ! same as case 0, but linesearch parameters written as indicated above
            call system(d_MOD%com)
            !-- MIO read the ouput file results and update 
            if(mio_read_model_output_files(errstruc,miostruc, d_OBS%h).ne.0) then
              call utl_bomb_out(errstruc)        
            endif 
        case (4) ! Parallel Jacobian Using Condor Externally
         !MD Done at the very beginning !
         !MD!   call system(d_MOD%com) ! run the model once, forward, to have outputs with current parameters 
         !MD!   !-- MIO read the ouput file results and update 
         !MD!   if(mio_read_model_output_files(errstruc,miostruc, d_OBS%h).ne.0) then
         !MD!     call utl_bomb_out(errstruc) 
         !MD!   endif
            call bxd_write_param_file(cv_PAR,d_PAR) ! write the parameter file
            call system('python CondorATC.py')
            select case (cv_A%jacobian_format)
              case ('binary')
                call readJCO(cv_A%jacfle, d_A)
              case ('ascii')
                call readJAC(cv_A%jacfle, d_A)
            end select

        end select
    
  

   end subroutine bpf_model_run
   
end module bayes_pest_model