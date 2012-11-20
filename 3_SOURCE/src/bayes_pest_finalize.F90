module bayes_pest_finalize

contains
subroutine bpd_finalize(d_PM, d_S, cv_PAR, d_PAR, &
                    &    cv_OBS, d_OBS, d_MOD, d_MIO, d_XQR)
       use bayes_pest_control
       implicit none
!--  Main Data Arrays for OBS and PARS

       integer              :: npar, nobs, ierr, i, j, k
       type (d_prior_mean)  :: d_PM  
       type (d_algorithmic) :: d_A        
       type (kernel_XQR)    :: d_XQR
       type (d_struct)      :: d_S 
       type (cv_param)      :: cv_PAR
       type (d_param)       :: d_PAR
       type (cv_observ)     :: cv_OBS
       type (d_observ)      :: d_OBS
       type (d_comlin)      :: d_MOD
       type (d_minout)      :: d_MIO                 
!-- deallocate algorithmic data array (d_A) 
  if (associated(d_A%H))           deallocate(d_A%H,stat=ierr)
  if (associated(d_A%HQHt))        deallocate(d_A%HQHt,stat=ierr)
  if (associated(d_A%Hsold))       deallocate(d_A%Hsold,stat=ierr)
  if (associated(d_A%Qsy))         deallocate(d_A%Qsy,stat=ierr)
  if (associated(d_A%Qyy))         deallocate(d_A%Qyy,stat=ierr)

!-- deallocate prior mean information (d_PM) 
  if (associated(d_PM%beta_0))     deallocate(d_PM%beta_0,stat=ierr)
  if (associated(d_PM%Qbb))        deallocate(d_PM%Qbb,stat=ierr)
  if (associated(d_PM%InvQbb))     deallocate(d_PM%InvQbb,stat=ierr)
  if (associated(d_PM%InvQbbB0))   deallocate(d_PM%InvQbbB0,stat=ierr)
  
!-- deallocate structural parameter data (d_XQR) 
  if (associated(d_XQR%X))         deallocate(d_XQR%X,stat=ierr)
  if (associated(d_XQR%Q0))        deallocate(d_XQR%Q0,stat=ierr)
  if (associated(d_XQR%R0))        deallocate(d_XQR%R0,stat=ierr)
  
!-- deallocate structural parameter data (d_S) 
  if (associated(d_S%theta_0))     deallocate(d_S%theta_0,stat=ierr)
  if (associated(d_S%theta_cov))   deallocate(d_S%theta_cov,stat=ierr)
  if (associated(d_S%invQtheta))   deallocate(d_S%invQtheta,stat=ierr)
  if (associated(d_S%struct_par_opt_vec_0)) deallocate(d_S%struct_par_opt_vec_0,stat=ierr)
  if (associated(d_S%struct_par_opt_vec))   deallocate(d_S%struct_par_opt_vec,stat=ierr)

!-- deallocate parameter control values (cv_PAR)
  if (associated(cv_PAR%grp_name))   deallocate(cv_PAR%grp_name,stat=ierr)       
  if (associated(cv_PAR%grp_type))   deallocate(cv_PAR%grp_type,stat=ierr)       
  if (associated(cv_PAR%derinc))    deallocate(cv_PAR%derinc,stat=ierr)       

!-- deallocate parameter structure (d_PAR)
  if (associated(d_PAR%group))     deallocate(d_PAR%group,stat=ierr)
  if (associated(d_PAR%pars))      deallocate(d_PAR%pars,stat=ierr)
  if (associated(d_PAR%parnme))    deallocate(d_PAR%parnme,stat=ierr)
  if (associated(d_PAR%pars_old))  deallocate(d_PAR%pars_old,stat=ierr)
  if (associated(d_PAR%pars_lns))  deallocate(d_PAR%pars_lns,stat=ierr)
  
!-- deallocate observation structure (d_OBS)
  if (associated(d_OBS%group))     deallocate(d_OBS%group,stat=ierr)
  if (associated(d_OBS%obs))       deallocate(d_OBS%obs,stat=ierr)
  if (associated(d_OBS%obsnme))    deallocate(d_OBS%obsnme,stat=ierr)
  if (associated(d_OBS%h))         deallocate(d_OBS%h,stat=ierr)
  
!-- deallocate model i/o structure (d_MIO)
  if (associated(d_MIO%tpl))       deallocate(d_MIO%tpl,stat=ierr)
  if (associated(d_MIO%infle))     deallocate(d_MIO%infle,stat=ierr)
  if (associated(d_MIO%ins))       deallocate(d_MIO%ins,stat=ierr)
  if (associated(d_MIO%outfle))    deallocate(d_MIO%outfle,stat=ierr)
  if (associated(d_MIO%pargroup))  deallocate(d_MIO%pargroup,stat=ierr)
  

  end subroutine bpd_finalize
           
end module bayes_pest_finalize