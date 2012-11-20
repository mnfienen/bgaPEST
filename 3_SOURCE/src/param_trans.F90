module param_trans

! ************************************************************************************************************
! Module with subroutines useful to transform sensitivity and parameters in estimation space and backtransform 
! ****************************************** M.D. 14/06/2010  ************************************************
     
 use bayes_pest_control
 use utilities  

contains

  subroutine sen_par_trans(cv_PAR, cv_OBS, d_PAR, d_A, d_PM)
  !Subroutine to transform the sensitivity and the parameters in
  !the estimation space.
  !Only if required and for the parameters that require transformation. 
  !LOG CASE:
  !First we do H=H*s (s is in the physical space)
  !Second s_old = log(s_old) !After s_old is in the estimation space
       implicit none
       
       integer                      :: i, j
       type (d_algorithmic)         :: d_A
       type (cv_param)              :: cv_PAR
       type (d_param)               :: d_PAR
       type (cv_observ)             :: cv_OBS
       type (d_prior_mean)          :: d_PM
       
       do i = 1,cv_PAR%npar 
         if (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.1) then
           do j = 1, cv_OBS%nobs
             d_A%H(j,i)=d_A%H(j,i)*d_PAR%pars(i)
           end do
           d_PAR%pars_old(i) = log(d_PAR%pars(i))   !MD At the beginning pars is the vector of the initial values of the parameters 
                 ! as read in the file. Then became the best estimate. Here we transform in the estimation space if required
         elseif (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.2) then !Power transformation
           d_PAR%pars_old(i) = d_PM%alpha(d_PAR%BetaAssoc(i))*((d_PAR%pars(i)**(1./d_PM%alpha(d_PAR%BetaAssoc(i))))-1)
           do j = 1, cv_OBS%nobs
             d_A%H(j,i)=d_A%H(j,i)*(((d_PAR%pars_old(i)/d_PM%alpha(d_PAR%BetaAssoc(i)))+1)**(d_PM%alpha(d_PAR%BetaAssoc(i))-1))
           end do
         end if
       end do
       
  end subroutine sen_par_trans
  
   subroutine par_back_trans(cv_PAR, d_PAR, d_PM)
  !Subroutine to back-transform the parameters in the physical space.
  !Only if required and for the parameters that were transformed. 
  !LOG CASE:
  !We do s = exp(s) After this s is again in the physical space
  
       implicit none
       
       integer                      :: i
       type (cv_param)              :: cv_PAR
       type (d_param)               :: d_PAR
       type (d_prior_mean)          :: d_PM
               
       do i = 1, cv_PAR%npar
         if (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.1) then
           d_PAR%pars(i) = exp(d_PAR%pars(i))  !Back-transform the parameters in the physical space
         elseif (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.2) then  !Power transformation
           d_PAR%pars(i)=((d_PAR%pars(i)/d_PM%alpha(d_PAR%BetaAssoc(i)))+1)**(d_PM%alpha(d_PAR%BetaAssoc(i)))
         end if
       enddo
       
  end subroutine par_back_trans
  
  subroutine par_back_trans_lns(cv_PAR, d_PAR, d_PM) 
  !Subroutine to back-transform the parameters in the physical space.
  !Only if required and for the parameters that were LOG transformed. 
  !We do s = exp(s) After this s is again in the physical space
  
       implicit none
       
       integer                      :: i
       type (cv_param)              :: cv_PAR
       type (d_param)               :: d_PAR
       type (d_prior_mean)          :: d_PM
        
       do i = 1, cv_PAR%npar
         if (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.1) then
           d_PAR%pars_lns(i) = exp(d_PAR%pars_lns(i)) !Back-transform the parameters in the physical space
         elseif (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.2) then  !Power transformation
           d_PAR%pars_lns(i)=((d_PAR%pars_lns(i)/d_PM%alpha(d_PAR%BetaAssoc(i)))+1)**(d_PM%alpha(d_PAR%BetaAssoc(i)))
         endif
       enddo

  end subroutine par_back_trans_lns

end module param_trans