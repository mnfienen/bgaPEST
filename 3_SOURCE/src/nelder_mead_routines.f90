module nelder_mead_simplex_routines
   use bayes_pest_control
   use error_message
   use utilities
   use model_input_output
   
      
contains


 subroutine nelmin_ls ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount,icount, numres,ifault, &
    & d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )

!*****************************************************************************80
!
!! NELMIN minimizes a function using the Nelder-Mead algorithm. 
!   
!  Discussion:
!
!    This routine seeks the minimum value of a user-specified function.
!
!    Simplex function minimisation procedure due to Nelder+Mead(1965),
!    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
!    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
!    25, 97) and Hill(1978, 27, 380-2)
!
!    The function to be minimized must be defined by a function of
!    the form
!
!      function fn ( x, f )
!      real ( kind = 8 ) fn 
!      real ( kind = 8 ) x(*)
!
!    and the name of this subroutine must be declared EXTERNAL in the
!    calling routine and passed as the argument FN.
!
!    This routine does not include a termination test using the
!    fitting of a quadratic surface.
!
!  Modified:
!
!    27 February 2008
!
!  Author:
!
!    FORTRAN77 version by R ONeill
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    John Nelder, Roger Mead,
!    A simplex method for function minimization,
!    Computer Journal,
!    Volume 7, 1965, pages 308-313.
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, external FN, the name of the function which evaluates
!    the function to be minimized.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
!    for the iteration.  On output, this data may have been overwritten.
!
!    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
!    is estimated to minimize the function.
!
!    Output, real ( kind = 8 ) YNEWLO, the minimum value of the function.
!
!    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
!    of function values.
!
!    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
!    initial simplex.  The relative magnitudes of its elements should reflect
!    the units of the variables.
!
!    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out 
!    every KONVGE iterations.
!
!    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function 
!    evaluations.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations 
!    used.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no errors detected.
!    1, REQMIN, N, or KONVGE has an illegal value.
!    2, iteration terminated because KCOUNT was exceeded without convergence.
!
  implicit none

  integer ( kind = 4 ) n
  real    ( kind = 8 ), parameter :: ccoeff = 0.5D+00
  real    ( kind = 8 ) del
  real    ( kind = 8 ) dn
  real    ( kind = 8 ) dnn
  real    ( kind = 8 ), parameter :: ecoeff = 2.0D+00
  real    ( kind = 8 ), parameter :: eps = 0.001D+00
  real    ( kind = 8 ), external :: fn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcount
  integer ( kind = 4 ) kcount
  integer ( kind = 4 ) konvge
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) numres
  real    ( kind = 8 ) p(n,n+1)
  real    ( kind = 8 ) p2star(n)
  real    ( kind = 8 ) pbar(n)
  real    ( kind = 8 ) pstar(n)
  real    ( kind = 8 ), parameter :: rcoeff = 1.0D+00
  real    ( kind = 8 ) reqmin
  real    ( kind = 8 ) rq
  real    ( kind = 8 ) start(n)
  real    ( kind = 8 ) step(n)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xmin(n)
  real    ( kind = 8 ) y(n+1)
  real    ( kind = 8 ) y2star
  real    ( kind = 8 ) ylo
  real    ( kind = 8 ) ynewlo
  real    ( kind = 8 ) ystar
  real    ( kind = 8 ) z
  
  type(kernel_XQR),    intent(in)     :: d_XQR
  type(d_struct),      intent(inout)  :: d_S
  type(cv_param),      intent(in)     :: cv_PAR
  type (d_comlin)                     :: d_MOD
  type(cv_algorithmic), intent(inout) :: cv_A
  type(d_algorithmic), intent(inout)  :: d_A
  type(d_param),       intent(inout)  :: d_PAR
  type(d_prior_mean),  intent(in)     :: d_PM
  type(cv_observ),     intent(in)     :: cv_OBS
  type(d_observ),      intent(in)     :: d_OBS
  type (cv_prior_mean), intent(in)    :: cv_PM
  type (mio_struc)                    :: miostruc
  type (err_failure_struc)            :: errstruc  
 
!  Check the input parameters.
  if ( reqmin <= 0.0D+00 ) then
    ifault = 1
    return
  endif
  if ( n < 1 ) then
    ifault = 1
    return
  endif
  if ( konvge < 1 ) then
    ifault = 1
    return
  endif
  icount = 0
  numres = 0
  jcount = konvge
  dn = real ( n, kind = 8 )
  nn = n + 1
  dnn = real ( nn, kind = 8 )
  del = 1.0D+00
  rq = reqmin * dn
!  Initial or restarted loop.
  do
    do i = 1, n
      p(i,nn) = start(i)
    enddo
    y(nn) = fn ( start,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )
    icount = icount + 1
    do j = 1, n
      x = start(j)
      start(j) = start(j) + step(j) * del
      do i = 1, n
        p(i,j) = start(i)
      enddo
      y(j) = fn ( start,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )
      icount = icount + 1
      start(j) = x
    enddo
!  The simplex construction is complete.                    
!  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
!  the vertex of the simplex to be replaced.
    ylo = y(1)
    ilo = 1
    do i = 2, nn
      if ( y(i) < ylo ) then
        ylo = y(i) 
        ilo = i
      endif
    enddo
!  Inner loop.
    do
      if ( kcount <= icount ) then
        exit
      endif
      ynewlo = y(1)
      ihi = 1
      do i = 2, nn
        if ( ynewlo < y(i) ) then
          ynewlo = y(i)
          ihi = i
        endif
      enddo
!  Calculate PBAR, the centroid of the simplex vertices
!  excepting the vertex with Y value YNEWLO.
      do i = 1, n
        z = 0.0D+00
        do j = 1, nn    
          z = z + p(i,j)
        enddo
        z = z - p(i,ihi)   
        pbar(i) = z / dn   
      enddo
!  Reflection through the centroid.
      do i = 1, n
        pstar(i) = pbar(i) + rcoeff * ( pbar(i) - p(i,ihi) )
      enddo
      ystar = fn ( pstar,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )
      icount = icount + 1
!  Successful reflection, so extension.
      if ( ystar < ylo ) then
        do i = 1, n
          p2star(i) = pbar(i) + ecoeff * ( pstar(i) - pbar(i) )
        enddo
        y2star = fn ( p2star,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )
        icount = icount + 1
!  Check extension.
        if ( ystar < y2star ) then
          do i = 1, n
            p(i,ihi) = pstar(i)
          enddo
          y(ihi) = ystar
!  Retain extension or contraction.
        else
          do i = 1, n
            p(i,ihi) = p2star(i)
          enddo
          y(ihi) = y2star
        endif
!  No extension.
      else
        l = 0
        do i = 1, nn
          if ( ystar < y(i) ) then
            l = l + 1
          endif
        enddo
        if ( 1 < l ) then
          do i = 1, n
            p(i,ihi) = pstar(i)
          enddo
          y(ihi) = ystar
!  Contraction on the Y(IHI) side of the centroid.
        else if ( l == 0 ) then
          do i = 1, n
            p2star(i) = pbar(i) + ccoeff * ( p(i,ihi) - pbar(i) )
          enddo
          y2star = fn ( p2star,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )
          icount = icount + 1
!  Contract the whole simplex.
          if ( y(ihi) < y2star ) then
            do j = 1, nn
              do i = 1, n
                p(i,j) = ( p(i,j) + p(i,ilo) ) * 0.5D+00
                xmin(i) = p(i,j)
              enddo
              y(j) = fn ( xmin,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )
              icount = icount + 1
            enddo
            ylo = y(1)
            ilo = 1
            do i = 2, nn
              if ( y(i) < ylo ) then
                ylo = y(i) 
                ilo = i
              endif
            enddo
            cycle
!  Retain contraction.
          else
            do i = 1, n
              p(i,ihi) = p2star(i)
            enddo
            y(ihi) = y2star
          endif
!  Contraction on the reflection side of the centroid.
        else if ( l == 1 ) then
          do i = 1, n
            p2star(i) = pbar(i) + ccoeff * ( pstar(i) - pbar(i) )
          enddo
          y2star = fn ( p2star,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )
          icount = icount + 1
!  Retain reflection?
          if ( y2star <= ystar ) then
            do i = 1, n
              p(i,ihi) = p2star(i)
            enddo
            y(ihi) = y2star
          else
            do i = 1, n
              p(i,ihi) = pstar(i)
            enddo
            y(ihi) = ystar  
          endif
        endif
      endif
!  Check if YLO improved.
      if ( y(ihi) < ylo ) then
        ylo = y(ihi)
        ilo = ihi
      endif
      jcount = jcount - 1
      if ( 0 < jcount ) then
        cycle
      endif
!  Check to see if minimum reached.
      if ( icount <= kcount ) then
        jcount = konvge
        z = 0.0D+00
        do i = 1, nn
          z = z + y(i)
        enddo
        x = z / dnn
        z = 0.0D+00
        do i = 1, nn
          z = z + ( y(i) - x )**2
        enddo
        if ( z <= rq ) then
          exit
        endif
      endif
    enddo
!  Factorial tests to check that YNEWLO is a local minimum.
    do i = 1, n
      xmin(i) = p(i,ilo)
    enddo
    ynewlo = y(ilo)
    if ( kcount < icount ) then
      ifault = 2
      exit
    endif
    ifault = 0
    do i = 1, n
      del = step(i) * eps
      xmin(i) = xmin(i) + del
      z = fn ( xmin,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      endif
      xmin(i) = xmin(i) - del - del
      z = fn ( xmin,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      endif
      xmin(i) = xmin(i) + del
    enddo
    if ( ifault == 0 ) then
      exit
    endif
!  Restart the procedure.
    do i = 1, n
      start(i) = xmin(i)
    enddo
    del = eps
    numres = numres + 1
  enddo
  return
end subroutine nelmin_ls




 subroutine nelmin_sp ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount,icount, numres,ifault, &
    & d_XQR, Q0_all,cv_OBS, d_OBS, cv_A, d_A, d_PAR, cv_S, d_S, d_PM, cv_PAR,cv_PM,nQ0)


  implicit none

  integer ( kind = 4 ) n
  real    ( kind = 8 ), parameter :: ccoeff = 0.5D+00
  real    ( kind = 8 ) del
  real    ( kind = 8 ) dn
  real    ( kind = 8 ) dnn
  real    ( kind = 8 ), parameter :: ecoeff = 2.0D+00
  real    ( kind = 8 ), parameter :: eps = 0.001D+00
  real    ( kind = 8 ), external :: fn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcount
  integer ( kind = 4 ) kcount
  integer ( kind = 4 ) konvge
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) numres
  real    ( kind = 8 ) p(n,n+1)
  real    ( kind = 8 ) p2star(n)
  real    ( kind = 8 ) pbar(n)
  real    ( kind = 8 ) pstar(n)
  real    ( kind = 8 ), parameter :: rcoeff = 1.0D+00
  real    ( kind = 8 ) reqmin
  real    ( kind = 8 ) rq
  real    ( kind = 8 ) start(n)
  real    ( kind = 8 ) step(n)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xmin(n)
  real    ( kind = 8 ) y(n+1)
  real    ( kind = 8 ) y2star
  real    ( kind = 8 ) ylo
  real    ( kind = 8 ) ynewlo
  real    ( kind = 8 ) ystar
  real    ( kind = 8 ) z
  
  integer                              :: nQ0  
  type(kernel_XQR),     intent(in)     :: d_XQR
  type(cv_observ),      intent(in)     :: cv_OBS
  type(d_observ),       intent(in)     :: d_OBS
  type(d_algorithmic),  intent(inout)  :: d_A 
  type(cv_algorithmic), intent(inout)  :: cv_A
  type(d_prior_mean),   intent(in)     :: d_PM
  type(cv_param),       intent(in)     :: cv_PAR 
  type(d_param),        intent(inout)  :: d_PAR        
  type(cv_struct),      intent(in)     :: cv_S 
  type(d_struct),       intent(inout)  :: d_S
  type(Q0_compr),       intent(in)     :: Q0_All(:)
  type(cv_prior_mean),  intent(in)     :: cv_PM  
 
!  Check the input parameters.
  if ( reqmin <= 0.0D+00 ) then
    ifault = 1
    return
  endif
  if ( n < 1 ) then
    ifault = 1
    return
  endif
  if ( konvge < 1 ) then
    ifault = 1
    return
  endif
  icount = 0
  numres = 0
  jcount = konvge
  dn = real ( n, kind = 8 )
  nn = n + 1
  dnn = real ( nn, kind = 8 )
  del = 1.0D+00
  rq = reqmin * dn
!  Initial or restarted loop.
  do
    do i = 1, n
      p(i,nn) = start(i)
    enddo
    y(nn) = fn (start,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)
    icount = icount + 1
    do j = 1, n
      x = start(j)
      start(j) = start(j) + step(j) * del
      do i = 1, n
        p(i,j) = start(i)
      enddo
      y(j) = fn (start,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)
      icount = icount + 1
      start(j) = x
    enddo
!  The simplex construction is complete.                    
!  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
!  the vertex of the simplex to be replaced.
    ylo = y(1)
    ilo = 1
    do i = 2, nn
      if ( y(i) < ylo ) then
        ylo = y(i) 
        ilo = i
      endif
    enddo
!  Inner loop.
    do
      if ( kcount <= icount ) then
        exit
      endif
      ynewlo = y(1)
      ihi = 1
      do i = 2, nn
        if ( ynewlo < y(i) ) then
          ynewlo = y(i)
          ihi = i
        endif
      enddo
!  Calculate PBAR, the centroid of the simplex vertices
!  excepting the vertex with Y value YNEWLO.
      do i = 1, n
        z = 0.0D+00
        do j = 1, nn    
          z = z + p(i,j)
        enddo
        z = z - p(i,ihi)   
        pbar(i) = z / dn   
      enddo
!  Reflection through the centroid.
      do i = 1, n
        pstar(i) = pbar(i) + rcoeff * ( pbar(i) - p(i,ihi) )
      enddo
      ystar = fn (pstar,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)
      icount = icount + 1
!  Successful reflection, so extension.
      if ( ystar < ylo ) then
        do i = 1, n
          p2star(i) = pbar(i) + ecoeff * ( pstar(i) - pbar(i) )
        enddo
        y2star = fn (p2star,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)
        icount = icount + 1
!  Check extension.
        if ( ystar < y2star ) then
          do i = 1, n
            p(i,ihi) = pstar(i)
          enddo
          y(ihi) = ystar
!  Retain extension or contraction.
        else
          do i = 1, n
            p(i,ihi) = p2star(i)
          enddo
          y(ihi) = y2star
        endif
!  No extension.
      else
        l = 0
        do i = 1, nn
          if ( ystar < y(i) ) then
            l = l + 1
          endif
        enddo
        if ( 1 < l ) then
          do i = 1, n
            p(i,ihi) = pstar(i)
          enddo
          y(ihi) = ystar
!  Contraction on the Y(IHI) side of the centroid.
        else if ( l == 0 ) then
          do i = 1, n
            p2star(i) = pbar(i) + ccoeff * ( p(i,ihi) - pbar(i) )
          enddo
          y2star = fn (p2star,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)
          icount = icount + 1
!  Contract the whole simplex.
          if ( y(ihi) < y2star ) then
            do j = 1, nn
              do i = 1, n
                p(i,j) = ( p(i,j) + p(i,ilo) ) * 0.5D+00
                xmin(i) = p(i,j)
              enddo
              y(j) = fn (xmin,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)
              icount = icount + 1
            enddo
            ylo = y(1)
            ilo = 1
            do i = 2, nn
              if ( y(i) < ylo ) then
                ylo = y(i) 
                ilo = i
              endif
            enddo
            cycle
!  Retain contraction.
          else
            do i = 1, n
              p(i,ihi) = p2star(i)
            enddo
            y(ihi) = y2star
          endif
!  Contraction on the reflection side of the centroid.
        else if ( l == 1 ) then
          do i = 1, n
            p2star(i) = pbar(i) + ccoeff * ( pstar(i) - pbar(i) )
          enddo
          y2star = fn (p2star,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)
          icount = icount + 1
!  Retain reflection?
          if ( y2star <= ystar ) then
            do i = 1, n
              p(i,ihi) = p2star(i)
            enddo
            y(ihi) = y2star
          else
            do i = 1, n
              p(i,ihi) = pstar(i)
            enddo
            y(ihi) = ystar  
          endif
        endif
      endif
!  Check if YLO improved.
      if ( y(ihi) < ylo ) then
        ylo = y(ihi)
        ilo = ihi
      endif
      jcount = jcount - 1
      if ( 0 < jcount ) then
        cycle
      endif
!  Check to see if minimum reached.
      if ( icount <= kcount ) then
        jcount = konvge
        z = 0.0D+00
        do i = 1, nn
          z = z + y(i)
        enddo
        x = z / dnn
        z = 0.0D+00
        do i = 1, nn
          z = z + ( y(i) - x )**2
        enddo
        if ( z <= rq ) then
          exit
        endif
      endif
    enddo
!  Factorial tests to check that YNEWLO is a local minimum.
    do i = 1, n
      xmin(i) = p(i,ilo)
    enddo
    ynewlo = y(ilo)
    if ( kcount < icount ) then
      ifault = 2
      exit
    endif
    ifault = 0
    do i = 1, n
      del = step(i) * eps
      xmin(i) = xmin(i) + del
      z = fn (xmin,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      endif
      xmin(i) = xmin(i) - del - del
      z = fn (xmin,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      endif
      xmin(i) = xmin(i) + del
    enddo
    if ( ifault == 0 ) then
      exit
    endif
!  Restart the procedure.
    do i = 1, n
      start(i) = xmin(i)
    enddo
    del = eps
    numres = numres + 1
  enddo
  return
end subroutine nelmin_sp

end module nelder_mead_simplex_routines