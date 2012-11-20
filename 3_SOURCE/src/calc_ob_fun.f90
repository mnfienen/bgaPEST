module objective_function
        use bayes_pest_control
        use utilities  
        
 contains
 
 subroutine cal_ob_funcs(d_XQR, d_S, d_PM, cv_PAR, cv_OBS, d_OBS,  d_A, d_PAR, cv_PM)
 
        implicit none
        ! declarations
        type(kernel_XQR),    intent(in)     :: d_XQR
        type(d_struct),      intent(inout)  :: d_S
        type(cv_param),      intent(in)     :: cv_PAR
        type(d_algorithmic), intent(inout)  :: d_A
        type(d_param),       intent(inout)  :: d_PAR
        type(d_prior_mean),  intent(in)     :: d_PM
        type(cv_observ),     intent(in)     :: cv_OBS
        type(d_observ),      intent(in)     :: d_OBS
        type (cv_prior_mean), intent(in)    :: cv_PM
        double precision,   allocatable     :: TMP(:,:), TVP(:),TMP1(:,:)
        integer                             :: i

!*****************************************************************************************************
! Calculate the objective functions
!***************************************************************************************************** 

!Regularization objective function phi_R = 1/2 ksit * HQHt * ksi + 1/2 ksit * (HX) * Qbb * (HXt) * ksi 
!*****************************************************************************************************
!First we calculate 1/2 ksit * HQHt * ksi 
 allocate(TVP(cv_OBS%nobs))
   TVP = UNINIT_REAL  ! vector (temporary vector)
   call DGEMV('n',cv_OBS%nobs,cv_OBS%nobs,1.D0,d_A%HQHt,cv_OBS%nobs,  &
       d_A%ksi,1,0.D0,TVP,1)
   call DGEMV('t',cv_OBS%nobs,1,5.0D-1,d_A%ksi,cv_OBS%nobs,  &
       TVP,1,0.D0,d_PAR%phi_R,1)
       
   if (cv_PM%betas_flag .ne. 0) then            !Here we calculate 1/2 ksit * (HX) * Qbb * (HXt) * ksi 
    allocate(TMP(cv_OBS%nobs,cv_PAR%p))         !only if we have prior mean informations and we add the 
    allocate(TMP1(cv_OBS%nobs,cv_OBS%nobs))     !result to the previous phi_R
     TMP  = UNINIT_REAL  ! matrix
     TMP1 = UNINIT_REAL  ! matrix
     call dgemm('n','n',cv_OBS%nobs,cv_PAR%p,cv_PAR%p, & !-----!
       1.D0,d_A%HX,cv_OBS%nobs,d_PM%Qbb,cv_PAR%p, & !----------!  Here we calculate (HX) * Qbb * (HX)t = TMP1 
       0.D0, TMP, cv_OBS%nobs) !-------------------------------!  Is a temporary matrix that can be deallocated
     call dgemm('n','t',cv_OBS%nobs,cv_OBS%nobs,cv_PAR%p, & !--!  after the use
        1.D0,TMP,cv_OBS%nobs,d_A%HX,cv_OBS%nobs, &  !----------! 
        0.D0, TMP1, cv_OBS%nobs) !-----------------------------!
    if (allocated(TMP))      deallocate(TMP)
     call DGEMV('n',cv_OBS%nobs,cv_OBS%nobs,1.D0,TMP1,cv_OBS%nobs,  &
       d_A%ksi,1,0.D0,TVP,1)
     call DGEMV('t',cv_OBS%nobs,1,5.0D-1,d_A%ksi,cv_OBS%nobs,  &
       TVP,1,1.D0,d_PAR%phi_R,1)  
    if (allocated(TMP1))    deallocate(TMP1)           
   endif
   if (allocated(TVP))      deallocate(TVP)           

!Misfit objective function phi_M = 1/2* (y-h(s))t * R^-1 * (y-h(s))
!**********************************************************************************
!We use just a loop because R0 is diagonal *** Must change if allow full R0 matrix
   d_PAR%phi_M = 0.
   do i = 1, cv_OBS%nobs
     d_PAR%phi_M = d_PAR%phi_M + (1./(d_S%sig*d_XQR%R0(i,i)))*((d_OBS%obs(i)  - d_OBS%h(i))**2)
   enddo     
   d_PAR%phi_M=0.5*d_PAR%phi_M

!Total objective function phi_T= phi_R +phi_M
!*********************************************
  d_PAR%phi_T = d_PAR%phi_R + d_PAR%phi_M 
   
!*****************************************************************************************************
! End Calculate the objective functions
!***************************************************************************************************** 

end subroutine cal_ob_funcs       
       
        
end module objective_function