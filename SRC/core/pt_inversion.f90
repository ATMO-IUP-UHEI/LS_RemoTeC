module pt_regularization_module
   use header_module
   implicit none
   private

   public :: pt_inversion, tsvd_inversion, column_inversion
   private :: l_curve


!*** module variables
   real(double), parameter :: epsilon = 1.d-12


   contains

! Modified by Andre Galli, June 30, 2010
! CHANGES: Don't add the apriori vector at each step: Keep x_state
!          instead of  x_state=x_state+(1-A)*x_apr
! AFFECTED SUBROUTINES: pt_inversion and tsvd_inversion

! Modified by Andre Galli, July 15, 2010
! CHANGES: x(i,iter) = boundary(i) if x(i,iter) would be outside boundary(i) both for PT and TSVD inversion

!modified by Andre Galli, July 20 2010:
!CHANGES: After fixing noise bug which collided with usage of x_apr: switched back to adding (1-A)x_apr to the state vector

!modified by Andre Galli, September 20 2010:
!CHANGES: Reaction to boundary hits: A priori a boundary hit (ID = 4) will cause the retrieval to be aborted unless
!         it is a benign boundary hit. If aerosols are retrieved: lower boundary of number of aerosols or power law (ID = 1). For cirri: lower boundary of cirrus numbers (ID=2).

!modified by Andre Galli, October 11 2010:
!CHANGES: Subroutine pt_inversion now works with full regularization of aerosol parameters. Changes in detail:
!         ir_start now includes aerosol parameters,
!         skill_reg1 = 1
!         In the lcurve subroutine a variable lcurve_Flag may be set if the DFS at the chosen point of the L-curve is smaller than naer + ntype_target.
!                  This information is then used in priofile_inversion.f90 to reduce the a priori optical thickness (and number of retrieval parameters
!                  in siome cases)
!                  SET INACTIVE FOR THE TIME BEING

!         Subroutines pt_inversion, tsvd_inversion, and lmlsq_inversion now use for weighting rdum = 1/MAX(DABS(x_apr(k),1d-8))
!                     instead of rdum = MAXVAL(DABS(kmat(1:ny,k)))

!modified by Andre Galli, November 2, 2010:
!CHANGES: Uncommented subroutine lmlsq until further use

!modified by Andre Galli, November 2, 2010:
!CHANGES: fixed bug in boundary_flag setting (if a malign boundary hit occurs, boundary_flag = 4 must not be overwritten)

!modified by Andre Galli, March 9, 2011:
!CHANGES: For pt_inversion use real(double) instead of quadruple u_mat_quad, v_mat_quad, sing_val_quad and then
!         call dsvdcmp instead of qsvdcmp. For TSVD inversion we are forced to use qsvdcmp because the matrices are
!         more stiff (different weighting factors at the beginning)

!modified by Andre Galli, June 9, 2011:
!CHANGES: Fixed a bug: nstate has to be replaced by the actual dimension nx for boundary hit treatment, since
!         nx < nstate if aerosol or offsets are no longer to be retrieved during the retrieval

!modified by Andre Galli, July 26, 2011:
!CHANGES: Deleted spurious lcurve_flag

  subroutine pt_inversion(&
         ntype_target, naer, &
         ny, nx, regy, regskill,&
         kmat, ymeas, ymod, s_y,&
         x_i, s_x, x_apr, a_avg,&
	 d_mat, lambda, dfs, dfs_target, dfs_scat, upbd, lowbd, cvflag, Boundary_Flag)

!*****************************************************************************
! This subroutine solves one iteration step in an iterative retrieval         *
! using 0th order Phillips-Tikhonov regularization (i.e. norm as side         *
! constraint). Input: kmat(ny,nx) is the weighting function matrix,           *
! ymeas(ny) is the measurement vector, ymod(ny) is the forward model vector,  *
! s_y(ny,ny) is the measurement error covariance matrix, x_i(nx)              *
! is the state vector of the current iteration step,                          *
! where ny is the number of measurements and nx is the number of unknown      *
! parameters. Output: x_i(nx) is the new state vector (replaces the old one), *
! s_x(nx,nx) is the retrieval error covariance matrix, a_avg(nx,nx) is the    *
! averaging kernel, and dfs is the Degrees of Freedom for Signal. For more    *
! information see: Hasekamp, O.P., and J. Landgraf, JGR, 106, 8077-8088, 2001.*
!*****************************************************************************
!*** Input
   integer, intent(in) :: ntype_target, naer, ny, nx
   real(double), intent(in) :: lambda
   integer, dimension(nx), intent(in) :: regskill
   real(double), dimension(ny), intent(in) :: ymeas, ymod, s_y, regy
   real(double), dimension(nx), intent(in) :: upbd, lowbd
   real(double), dimension(nx), intent(in) :: x_apr
!**In/out
   real(double), dimension(nx), intent(inout) :: x_i
   real(double), dimension(ny,nx), intent(inout) :: kmat
!*** Output
   integer, intent(out) :: cvflag, Boundary_Flag
   real(double), intent(out) :: dfs, dfs_scat
   real(double), dimension(ntype_target), intent(out) :: dfs_target
   real(double), dimension(nx,ny), intent(out) :: d_mat
   real(double), dimension(nx,nx), intent(out) :: s_x, a_avg
!*** Local
   real(double), dimension(nx) :: x_0, x_sm, weight, x_apr_tmp
   real(double), dimension(ny) :: y, ylsq, ylin
   real(double), dimension(ny,nx) :: u_mat
   real(double), dimension(nx,nx) :: v_mat, filt_mat,  unit_avg
   real(double), dimension(nx,ny) :: d_tmp
   real(double), dimension(nx) :: sing_val, u_ylsq, u_ylin, v_x0
   real(double), dimension(nx) :: filt_func_pt, filt_func_lsq, filt_func_x0
   integer,dimension(nx) :: index_out
   real(double) :: reg_par, rdum, rdummax, skill_reg_0, skill_reg_1, skill_reg_x
   integer :: i, j, k, l
   integer :: ir_start, ir_stop
   real(double), dimension(ny,nx) :: u_dummy
   real(double), dimension(nx,nx) :: v_dummy
!*********************************************************************

      boundary_flag = 0

!*** Skill the kernel matrix
      skill_reg_0 = 1.D7
      skill_reg_1 = 1.
      skill_reg_x = 1.

      ir_start = count(regskill(:)==0)+1
      ir_stop = nx

      !*** CHI2 parameters
      do k=1,nx
!         rdum = MAXVAL(DABS(kmat(1:ny,k)))
         rdum = 1./max(DABS(x_apr(k)),1d-8)
	 if(regskill(k)==0)then
            if(rdum.ne.0.D0)then
	       weight(k) = 1./rdum*skill_reg_0
            else
	       weight(k) = 1.D0*skill_reg_0
            endif
         endif
      enddo

      !*** Aerosol parameters
      do k=1,nx
!         rdum = MAXVAL(DABS(kmat(1:ny,k)))
         rdum = 1./max(DABS(x_apr(k)),1d-8)
         if(regskill(k)==1)then
            if(rdum .ne. 0.D0)then
      	       weight(k) = 1./rdum*skill_reg_1
            else
	       weight(k) = 1.D0*skill_reg_1
            endif
         endif
      enddo

      !*** Profile parameters
      do j=1,ntype_target
         rdummax=0
         do k=1,nx
!            rdum = MAXVAL(DABS(kmat(1:ny,k)))
            rdum = 1./max(DABS(x_apr(k)),1d-8)
	    if(regskill(k)==1+j .and. rdum > rdummax)then
               rdummax = rdum
            endif
         enddo
         do k=1,nx
            if(regskill(k)==1+j)then
               if(rdummax .ne. 0.D0)then
	          weight(k) = 1./rdummax*skill_reg_x
               else
	          weight(k) = 1.D0*skill_reg_x
               endif
            endif
         enddo
      enddo

      do l=1,ny
         kmat(l,1:nx) = kmat(l,1:nx)*weight(1:nx)
      enddo

      x_i=x_i/weight
      x_0=0./weight
      x_apr_tmp=x_apr/weight

!*** Define measurement vector, weighting by errors
      y = ymeas - ymod + matmul(kmat,x_i)
      ylsq = ymeas - ymod
      ylin = matmul(kmat,x_i)
      forall(l=1:ny)
         y(l) =  y(l)/DSQRT(s_y(l))
         ylsq(l) =  ylsq(l)/DSQRT(s_y(l))
         ylin(l) =  ylin(l)/DSQRT(s_y(l))
      end forall

!*** SVD of the kernel matrix
      forall(l=1:ny)
         u_mat(l,1:nx) = kmat(l,1:nx)/DSQRT(s_y(l))
      end forall

      call dsvdcmp(u_mat, ny, nx, ny, nx, sing_val, v_mat, cvflag)

!*** Sorting of the singular values
      call sort(sing_val,nx,index_out)

      do l=1,nx
         v_dummy(:,l)=v_mat(:,index_out(l))
         u_dummy(:,l)=u_mat(:,index_out(l))
      enddo
      sing_val(1:nx)=sing_val(nx:1:-1)
      v_mat(:,1:nx)=v_dummy(:,nx:1:-1)
      u_mat(:,1:nx)=u_dummy(:,nx:1:-1)


!*** Call the L-curve, regularize from ir_start to ir_stop
      call l_curve(  &
           nx, ny, u_mat, v_mat, sing_val, y, x_0,&
           ir_start, ir_stop, reg_par, regy)

      do i=1,nx
         u_ylsq(i) = dot_product(u_mat(:,i),ylsq)
         u_ylin(i) = dot_product(u_mat(:,i),ylin)
         v_x0(i)   = dot_product(v_mat(:,i),x_0)
      enddo

!*** Calculate filter function
      filt_func_lsq = 1./(1.+lambda)
      filt_func_x0  = reg_par**2/(sing_val**2+reg_par**2)

      filt_func_pt  = sing_val**2/(sing_val**2+reg_par**2)

! degree of freedom without albedo
!      dfs = sum(filt_func_pt(ir_start:ir_stop))

!*** Calculate new state vector
      x_i=0.D0
      filt_mat = 0.D0
      do i=1,nx
         filt_mat(i,i) = filt_func_pt(i)
         if (sing_val(i)/sing_val(1) .gt. epsilon) then
            x_i = x_i + &
              (filt_func_pt(i)*u_ylin(i)/sing_val(i) + &
               filt_func_lsq(i)*filt_func_pt(i)*u_ylsq(i)/sing_val(i) + &
               filt_func_lsq(i)*filt_func_x0(i)*v_x0(i))*v_mat(:,i)
         endif
      enddo

!*** Calculate contribution function matrix, error covariance, averaging kernels
      d_tmp = transpose(u_mat)
      do l=1,ny
         do i=1,nx
            if (sing_val(i)/sing_val(1) .gt. epsilon) then
               d_tmp(i,l) = filt_func_pt(i) / sing_val(i) * d_tmp(i,l)
            endif
         enddo
      enddo
      d_mat = matmul(v_mat,d_tmp)
      s_x = matmul(d_mat,transpose(d_mat))
      a_avg = matmul(v_mat,matmul(filt_mat,transpose(v_mat)))

!*** total degree of freedom
      dfs = 0.D0
      do k = 1, nx
         dfs = dfs + a_avg(k,k)
      enddo
!*** degree of freedom of target absorbers
      dfs_target = 0.D0
      do j = 1, ntype_target
         rdum = 0.D0
         do k = 1, nx
            if(regskill(k)==1+j)then
               rdum = rdum + a_avg(k,k)
            endif
         enddo
	 dfs_target(j) = rdum
      enddo
!*** degree of freedom of aerosol parameters
      dfs_scat = 0.D0
      rdum = 0.D0
      do k = 1, nx
         if(regskill(k)==1) then
            rdum = rdum+a_avg(k,k)
         endif
      enddo
      dfs_scat = rdum

!*** Add a apriori, x=A*xtrue+(1-A)*xapr
      unit_avg=-a_avg
      do i=1,nx
      	 unit_avg(i,i)=1.-a_avg(i,i)
      enddo

      x_sm=matmul(unit_avg,x_apr_tmp-x_0)

      do k=1,nx
	 x_i(k)=x_i(k)+x_sm(k)
      enddo

!*** Rescale retrieval parameters
      x_i=x_i*weight
      x_0=x_0*weight
!      x_apr=x_apr*weight
      forall (i=1:nx,j=1:nx) s_x(i,j) = s_x(i,j)*weight(i)*weight(j)
      forall (i=1:nx,j=1:nx) a_avg(i,j) = a_avg(i,j)*weight(i)/weight(j)
      forall (i=1:nx,j=1:ny) d_mat(i,j) = d_mat(i,j)*weight(i)/DSQRT(s_y(j))

!*** Check for hitting the boundaries
      do i=1,nx
         if(x_i(i)<lowbd(i))then
            print*, 'PT_INVERSION: LM: HIT LOWER BOUNDARY: ',i,lowbd(i),x_i(i)
            x_i(i) = lowbd(i)
            if(i>nx-naer .and. Boundary_Flag==0) then
               Boundary_Flag=1
            else
               Boundary_Flag=2
            endif
         elseif(x_i(i)>upbd(i)) then
            print*, 'PT_INVERSION: LM: HIT UPPER BOUNDARY: ',i,upbd(i),x_i(i)
            x_i(i) = upbd(i)
            if(i>nx-naer .and. Boundary_Flag==0) then
               Boundary_Flag = 1
            else
               Boundary_Flag = 2
            endif
	 endif
      enddo

  end subroutine pt_inversion

!*********************************************************************************************

  subroutine tsvd_inversion(&
         ntype_target, naer, &
         ny, nx, nlay, regskill,&
         kmat, ymeas, ymod, s_y,&
         x_i, s_x, x_apr, a_avg,&
	 d_mat, lambda, dfs, dfs_target, dfs_scat, upbd, lowbd, cvflag, Boundary_Flag)
!*** Input
   integer, intent(in) :: ntype_target, naer, ny, nx, nlay
   real(double), intent(in) :: lambda
   integer, dimension(nx), intent(in) :: regskill
   real(double), dimension(ny), intent(in) :: ymeas, ymod, s_y
   real(double), dimension(nx), intent(in) :: upbd, lowbd
   real(double), dimension(nx), intent(in) :: x_apr
!**In/out
   real(double), dimension(nx), intent(inout) :: x_i
   real(double), dimension(ny,nx), intent(inout) :: kmat
!*** Output
   integer, intent(out) :: cvflag, Boundary_Flag
   real(double), intent(out) :: dfs, dfs_scat
   real(double), dimension(ntype_target), intent(out) :: dfs_target
   real(double), dimension(nx,ny), intent(out) :: d_mat
   real(double), dimension(nx,nx), intent(out) :: s_x, a_avg
!*** Local
   real(double), dimension(nx) :: x_0, x_sm, weight, x_apr_tmp
   real(double), dimension(ny) :: y, ylsq, ylin
   real(double), dimension(ny,nx) :: u_mat
   real(double), dimension(nx,nx) :: v_mat, filt_mat,  unit_avg
   real(double), dimension(nx,ny) :: d_tmp
   real(double), dimension(nx) :: sing_val, u_ylsq, u_ylin, v_x0
   real(double), dimension(nx) :: filt_func_pt, filt_func_lsq, filt_func_x0
   integer,dimension(nx) :: index_out
   real(double) :: rdum, rdummax, skill_reg_0, skill_reg_1, skill_reg_x
   integer :: i, j, k, l
   integer :: ir_start, ir_stop
   real(double), dimension(ny,nx) :: u_dummy
   real(double), dimension(nx,nx) :: v_dummy
   real(quad), dimension(ny,nx) :: u_mat_quad
   real(quad), dimension(nx,nx) :: v_mat_quad
   real(quad), dimension(nx) :: sing_val_quad
   integer :: i1, i2, ir_start_2
   real(double), dimension(nx) :: s
   integer, dimension(1) :: imax
   character(stringlen) :: message
!*********************************************************************
      boundary_flag = 0

!*** Skill the skernel matrix
      skill_reg_0 = 1.D14
      skill_reg_1 = 1.D7
      skill_reg_x = 1.

      ir_start = count(regskill(:)==0)+count(regskill(:)==1)+1
      ir_start_2 = count(regskill(:)==0)+1
      ir_stop = nx

      !*** CHI2 parameters
      do k=1,nx
!         rdum = MAXVAL(DABS(kmat(1:ny,k)))
         rdum = 1./max(DABS(x_apr(k)),1d-8)
	 if(regskill(k)==0)then
            if(rdum.ne.0.D0)then
	       weight(k) = 1./rdum*skill_reg_0
            else
	       weight(k) = 1.D0*skill_reg_0
            endif
         endif
      enddo

      !*** Aerosol parameters
      do k=1,nx
!         rdum = MAXVAL(DABS(kmat(1:ny,k)))
         rdum = 1./max(DABS(x_apr(k)),1d-8)
         if(regskill(k)==1)then
            if(rdum .ne. 0.D0)then
      	       weight(k) = 1./rdum*skill_reg_1
            else
	       weight(k) = 1.D0*skill_reg_1
            endif
         endif
      enddo

      !*** Profile parameters
      do j=1,ntype_target
         rdummax=0
         do k=1,nx
!            rdum = MAXVAL(DABS(kmat(1:ny,k)))
            rdum = 1./max(DABS(x_apr(k)),1d-8)
	    if(regskill(k)==1+j .and. rdum > rdummax)then
               rdummax = rdum
            endif
         enddo
         do k=1,nx
            if(regskill(k)==1+j)then
               if(rdummax .ne. 0.D0)then
	          weight(k) = 1./rdummax*skill_reg_x
               else
	          weight(k) = 1.D0*skill_reg_x
               endif
            endif
         enddo
      enddo

      do l = 1, ny
         kmat(l,1:nx) = kmat(l,1:nx)*weight(1:nx)
      enddo
      x_i=x_i/weight
      x_0=0./weight
      x_apr_tmp=x_apr/weight

!*** Define measurement vector, weighting by errors
      y = ymeas - ymod + matmul(kmat,x_i)
      ylsq = ymeas - ymod
      ylin = matmul(kmat,x_i)

      forall(l=1:ny)
         y(l) =  y(l)/DSQRT(s_y(l))
         ylsq(l) =  ylsq(l)/DSQRT(s_y(l))
         ylin(l) =  ylin(l)/DSQRT(s_y(l))
      end forall

!*** SVD of the kernel matrix
      forall(l=1:ny)
         u_mat(l,1:nx) = kmat(l,1:nx)/DSQRT(s_y(l))
      end forall
      u_mat_quad = real(u_mat,quad)
      call qsvdcmp(u_mat_quad,ny,nx,ny,nx,sing_val_quad,v_mat_quad,cvflag)

      sing_val = real(sing_val_quad, double)
      v_mat = real(v_mat_quad, double)
      u_mat = real(u_mat_quad, double)

!*** Sorting of the singular values
      call sort(sing_val,nx,index_out)

      do l=1,nx
         v_dummy(:,l)=v_mat(:,index_out(l))
         u_dummy(:,l)=u_mat(:,index_out(l))
      enddo
      sing_val(1:nx)=sing_val(nx:1:-1)
      v_mat(:,1:nx)=v_dummy(:,nx:1:-1)
      u_mat(:,1:nx)=u_dummy(:,nx:1:-1)

      do i=1,nx
         u_ylsq(i) = dot_product(u_mat(:,i),ylsq)
         u_ylin(i) = dot_product(u_mat(:,i),ylin)
         v_x0(i)   = dot_product(v_mat(:,i),x_0)
      enddo

!*** Calculate filter function
!      IF(iter<=nlsq+1 .and. lambda > 0.D0)lambda=MIN(nx,10)
      filt_func_lsq = 1./(1.+lambda)
      filt_func_x0  = 0.D0
      filt_func_pt = 1.D0
      filt_func_pt(ir_start:nx)=0.D0

      !*** Pick non-zero singular values for target profiles (look for v_i without sign change)
      do j=1,ntype_target
	 s=0.D0
         do i=ir_start,nx
            i1=1+nlay*(j-1)
	    i2=nlay*j
            s(i)=sum(v_mat(i1:i2,i))/nlay
	 enddo
	 imax=maxloc(DABS(s))
	 filt_func_pt(imax(1))=1.D0

      enddo
      !*** Pick non-zero singular values for aerosol parameters(look for v_i without sign change)
      !k=0
      !DO i=1,nx
      !   IF(regskill(i)==1 .and. k<1)THEN
      !      imax=MAXLOC(DABS(v_mat(i,:)))
      !	    filt_func_pt(imax(1))=1.D0
      !	    k=k+1
      !	 ENDIF
      !ENDDO

!      dfs = sum(filt_func_pt(ir_start:ir_stop))

!*** Calculate new state vector
      x_i=0.D0
      filt_mat = 0.D0
      do i=1,nx
         filt_mat(i,i) = filt_func_pt(i)
         x_i = x_i + &
           (filt_func_pt(i)*u_ylin(i)/sing_val(i) + &
            filt_func_lsq(i)*filt_func_pt(i)*u_ylsq(i)/sing_val(i) + &
            filt_func_lsq(i)*filt_func_x0(i)*v_x0(i))*v_mat(:,i)
      enddo

!*** Calculate contribution function matrix, error covariance, averaging kernels
      d_tmp = transpose(u_mat)
      do l=1,ny
         do i=1,nx
            d_tmp(i,l) = filt_func_pt(i) / sing_val(i) * d_tmp(i,l)
         enddo
      enddo
      d_mat = matmul(v_mat,d_tmp)
      s_x = matmul(d_mat,transpose(d_mat))
      a_avg = matmul(v_mat,matmul(filt_mat,transpose(v_mat)))

!*** Get degree of freedom
      dfs=0.D0
      do k=1,nx
         dfs=dfs+a_avg(k,k)
      enddo
      dfs_target=0.D0
      do j=1,ntype_target
         rdum=0.D0
         do k=1,nx
            if(regskill(k)==1+j)then
               rdum=rdum+a_avg(k,k)
            endif
         enddo
         dfs_target(j)=rdum
      enddo
      dfs_scat=0.D0
      rdum=0.D0
      do k=1,nx
         if(regskill(k)==1)then
            rdum=rdum+a_avg(k,k)
         endif
      enddo
      dfs_scat=rdum

!*** Add a apriori, x=A*xtrue+(1-A)*xapr
      unit_avg=-a_avg
      do i=1,nx
     	 unit_avg(i,i)=1.-a_avg(i,i)
      enddo
      x_sm=matmul(unit_avg,x_apr_tmp-x_0)

      do k=1,nx
	 x_i(k)=x_i(k)+x_sm(k)
      enddo

!*** Rescale retrieval parameters
      x_i=x_i*weight
      x_0=x_0*weight
!      x_apr=x_apr*weight
      forall (i=1:nx,j=1:nx) s_x(i,j) = s_x(i,j)*weight(i)*weight(j)
      forall (i=1:nx,j=1:nx) a_avg(i,j) = a_avg(i,j)*weight(i)/weight(j)
      forall (i=1:nx,j=1:ny) d_mat(i,j) = d_mat(i,j)*weight(i)/DSQRT(s_y(j))
!*** Check for hitting the boundaries
      do i = 1, nx
         if(x_i(i)<lowbd(i))then
!	    if(outputflag >= 2) then
               write(message,'(a,x,i2.2,x,2(1pE24.8e3,X))')'TSVD_INVERSION: LM: HIT LOWER BOUNDARY: ',i,lowbd(i),x_i(i)
               call writelog(message, 5)
!            endif
            x_i(i) = lowbd(i)
            if(i>nx-naer) then
               Boundary_Flag=1
            else
               Boundary_Flag=2
            endif
         elseif(x_i(i)>upbd(i)) then
!	    if(outputflag >= 2) then
               write(message,'(a,x,i2.2,x,2(1pE24.8e3,X))')'TSVD_INVERSION: LM: HIT UPPER BOUNDARY: ',i,upbd(i),x_i(i)
               call writelog(message, 5)
!            endif
            x_i(i) = upbd(i)
            if(i>nx-naer) then
               Boundary_Flag = 1
            else
               Boundary_Flag = 2
            endif
	 endif
      enddo

  end subroutine tsvd_inversion

!*************************************************************************************
!*** Find regularization parameter reg_par with l-curve
   subroutine l_curve( &
         nx, ny, u, v, sing_val, y, x_0,&
         ir_start, ir_stop, reg_par, regy)
!*** Input
   integer, intent(in) :: nx, ny
   real(double), dimension(ny,nx), intent(in) :: u
   real(double), dimension(nx,nx), intent(in) :: v
   real(double), dimension(nx), intent(in) :: sing_val, x_0
   real(double), dimension(ny), intent(in) :: regy
!*** In/out
  real(double), dimension(ny), intent(inout) :: y
!*** Output
   real(double) :: reg_par
!*** Local
   real(double), dimension(nx) :: u_y
   integer :: j,k,l,n
   integer :: ir_start, ir_stop
   integer, parameter :: n_lambda = 2000
   real(double) :: lambda, reg_min, reg_max, step, alpha,rho,eta,eta_1,&
           sum_term_rho,sum_term_eta,phi,d0, rnom, denomi,fact,curv,&
           curv_max, eta_max,rho_max,deg_freedom, deg_freedom_max
   integer :: i_max
!   real(double), dimension(ny,ny) :: u_tmp

   real(double), dimension(ny) :: y_tmp
   real(double), dimension(nx) :: v_x0, filt_func_dummy
   !HH: Make matrix allocatable to avoid too large stack size for multi-theaded processing
   real(double), dimension(:,:), allocatable :: u_tmp
!----------------------------------------------------------------------------

      allocate(u_tmp(ny,ny))
      u_tmp = -matmul(u, transpose(u))
      forall(l=1:ny)
         u_tmp(l,l) = 1.D0 + u_tmp(l,l)
      end forall

      do k=1,nx
         u_y(k) = dot_product(u(:,k),y)
         v_x0(k) = dot_product(v(:,k),x_0)
      enddo

      do k=1,ny
         if(regy(k)==0) y(k)=0.D0
      enddo

      y_tmp = matmul(U_TMP,y)
      d0 = sum(y_tmp(1:ny)**2)
      reg_max = sing_val(ir_start)
!      reg_min = sing_val(ir_stop)*1.d-2
      do k = 0, ir_stop-1
!         if (sing_val(ir_stop-k) .ne. 0.d0) then
         if (sing_val(ir_stop-k)/sing_val(1) .gt. epsilon) then
            reg_min = sing_val(ir_stop-k)*1.d-2
            exit
         endif
      enddo

      step = (DLOG10(reg_max) - DLOG10(reg_min))/ dble(n_lambda)
      CURV_MAX = 0.D0
      deg_freedom_max = 0.0d0

      do n=1,n_lambda
         lambda = 10**(DLOG10(reg_min) + dble(n)*step)
         alpha = lambda**2
         rho = 0.
         eta = 0.
         eta_1 = 0.
         do j=1,nx
            sum_term_rho =&
                 (1.D0-(sing_val(j)**2/(sing_val(j)**2 + lambda**2)))*(u_y(j)-sing_val(j)*v_x0(j))
            rho = rho + (sum_term_rho)**2
!            phi = sing_val(j)**2/(sing_val(j)**2+lambda**2)
!            sum_term_eta = phi*u_y(j)/sing_val(j)+(1.-phi)*v_x0(j)
            phi = sing_val(j)/(sing_val(j)**2+lambda**2)
            sum_term_eta = phi*u_y(j) +(1.-phi*sing_val(j))*v_x0(j)
            eta = eta + sum_term_eta**2
!            eta_1 = eta_1 + ((1.-phi)*(phi**2)*u_y(j)/sing_val(j)+phi*(1-phi)**2*v_x0(j))*(u_y(j)/sing_val(j)-v_x0(j))
           eta_1 = eta_1 + ((1.-sing_val(j)*phi)*(phi**2)*u_y(j)+phi*(1-sing_val(j)*phi)**2*v_x0(j))*(u_y(j)-v_x0(j)*sing_val(j))
         enddo
         rho = rho+d0
         eta_1 = -(4./lambda)*eta_1
         rnom  = &
              (((lambda**2)*eta_1*rho) + &
              2.*lambda*eta*rho + &
              ((lambda**4)*eta*eta_1))
         denomi = ((lambda**4)*(eta**2) + rho**2)**(3./2.)
         fact = 2.*((eta*rho)/eta_1)
         curv = (-fact*(rnom/denomi))
         filt_func_dummy  = sing_val**2/(sing_val**2+lambda**2)
         deg_freedom = sum(filt_func_dummy(ir_start:ir_stop))

	 if(curv.gt.curv_max)then
            curv_max = curv
            i_max = n
            eta_max = eta
            rho_max = rho
            reg_par = lambda
            deg_freedom_max = deg_freedom
         endif

      enddo

      return
   end subroutine l_curve
!-------------------------------------------------------------------------
   subroutine column_inversion(&
        ntype_target, naer, &
        ny, nx, regy, regskill,&
        kmat, ymeas, ymod, s_y,&
        x_i, s_x, x_apr, a_avg,&
        d_mat, lambda, dfs, dfs_target, dfs_scat, upbd, lowbd, cvflag, Boundary_Flag)
!*** Input
   integer, intent(in) :: ntype_target, naer, ny, nx
   real(double), intent(in) :: lambda
   integer, dimension(nx), intent(in) :: regskill
   real(double), dimension(ny), intent(in) :: ymeas, ymod, s_y, regy
   real(double), dimension(nx), intent(in) :: upbd, lowbd
   real(double), dimension(nx), intent(in) :: x_apr
!**In/out
   real(double), dimension(nx), intent(inout) :: x_i
   real(double), dimension(ny,nx), intent(inout) :: kmat
!*** Output
   integer, intent(out) :: cvflag, Boundary_Flag
   real(double), intent(out) :: dfs, dfs_scat
   real(double), dimension(ntype_target), intent(out) :: dfs_target
   real(double), dimension(nx,ny), intent(out) :: d_mat
   real(double), dimension(nx,nx), intent(out) :: s_x, a_avg
!*** Local
   real(double), dimension(nx) :: x_0, x_sm, weight, x_apr_tmp
   real(double), dimension(ny) :: y, ylsq, ylin
   real(double), dimension(ny,nx) :: u_mat
   real(double), dimension(nx,nx) :: v_mat, filt_mat,  unit_avg
   real(double), dimension(nx,ny) :: d_tmp
   real(double), dimension(nx) :: sing_val, u_ylsq, u_ylin, v_x0
   real(double), dimension(nx) :: filt_func_pt, filt_func_lsq, filt_func_x0
   integer,dimension(nx) :: index_out
   real(double) :: reg_par, rdum, rdummax, skill_reg_0, skill_reg_1, skill_reg_x
   integer :: i, j, k, l
   integer :: ir_start, ir_stop
   real(double), dimension(ny,nx) :: u_dummy
   real(double), dimension(nx,nx) :: v_dummy
   real(double) :: delta
     !-------------------------------------------------------------------------

      Boundary_Flag = 0

     !*** Skill the kernel matrix
     skill_reg_0 = 1.D7
     skill_reg_1 = 1.
     skill_reg_x = 1.

     ir_start = count(regskill(:)==0)+1
     ir_stop = nx

     !*** CHI2 parameters
     do k=1,nx
        !         rdum = MAXVAL(DABS(kmat(1:ny,k)))
        rdum = 1./max(DABS(x_apr(k)),1d-8)
        if(regskill(k)==0)then
           if(rdum.ne.0.D0)then
              weight(k) = 1./rdum*skill_reg_0
           else
              weight(k) = 1.D0*skill_reg_0
           endif
        endif
     enddo

     !*** Aerosol parameters
     do k=1,nx
        !         rdum = MAXVAL(DABS(kmat(1:ny,k)))
        rdum = 1./max(DABS(x_apr(k)),1d-8)
        if(regskill(k)==1)then
           if(rdum .ne. 0.D0)then
              weight(k) = 1./rdum*skill_reg_1
           else
              weight(k) = 1.D0*skill_reg_1
           endif
        endif
     enddo

     !*** Profile parameters
     do j=1,ntype_target
        rdummax=0
        do k=1,nx
           !            rdum = MAXVAL(DABS(kmat(1:ny,k)))
           rdum = 1./max(DABS(x_apr(k)),1d-8)
           if(regskill(k)==1+j .and. rdum > rdummax)then
              rdummax = rdum
           endif
        enddo
        do k=1,nx
           if(regskill(k)==1+j)then
              if(rdummax .ne. 0.D0)then
                 weight(k) = 1./rdummax*skill_reg_x
              else
                 weight(k) = 1.D0*skill_reg_x
              endif
           endif
        enddo
     enddo

     do l=1,ny
        kmat(l,1:nx) = kmat(l,1:nx)*weight(1:nx)
     enddo
     x_i=x_i/weight
     x_0=0./weight
     x_apr_tmp=x_apr/weight

     !*** Define measurement vector, weighting by errors
     y = ymeas - ymod + matmul(kmat,x_i)
     ylsq = ymeas - ymod
     ylin = matmul(kmat,x_i)

     forall(l=1:ny)
        y(l) =  y(l)/DSQRT(s_y(l))
        ylsq(l) =  ylsq(l)/DSQRT(s_y(l))
        ylin(l) =  ylin(l)/DSQRT(s_y(l))
     end forall

     !*** SVD of the kernel matrix
     forall(l=1:ny)
        u_mat(l,1:nx) = kmat(l,1:nx)/DSQRT(s_y(l))
     end forall
     call dsvdcmp(u_mat, ny, nx, ny, nx, sing_val, v_mat, cvflag)

     !*** Sorting of the singular values
     call sort(sing_val,nx,index_out)

     do l=1,nx
        v_dummy(:,l)=v_mat(:,index_out(l))
        u_dummy(:,l)=u_mat(:,index_out(l))
     enddo
     sing_val(1:nx)=sing_val(nx:1:-1)
     v_mat(:,1:nx)=v_dummy(:,nx:1:-1)
     u_mat(:,1:nx)=u_dummy(:,nx:1:-1)

     !*** Find regularization parameter for which noise on aerosol parameters
     !*** is lower than allowed range of aerosol parameters
     if (naer > 0 ) delta = (log10(sing_val(nx-naer+1)) - log10(sing_val(nx)/100.d0))/200.d0
     do k = -1, 200
        if (k==-1) then
           reg_par = 0.d0
        else
           reg_par = 10.d0**(log10(sing_val(nx)/100) + k*delta)
        endif
        !*** Calculate filter function
        filt_func_pt  = sing_val**2/(sing_val**2+reg_par**2)

        !*** Calculate contribution function matrix, error covariance, averaging kernels
        d_tmp = transpose(u_mat)
        do l=1,ny
           do i=1,nx
              if (sing_val(i)/sing_val(1) .gt. epsilon) then
                 d_tmp(i,l) = filt_func_pt(i) / sing_val(i) * d_tmp(i,l)
              endif
           enddo
        enddo
        d_mat = matmul(v_mat,d_tmp)
        s_x = matmul(d_mat,transpose(d_mat))
        forall (i=1:nx,j=1:nx) s_x(i,j) = s_x(i,j)*weight(i)*weight(j)

        !*** Exit do-loop if noise of aerosol parameters is smaller than allowed range
        !*** Take minimum regularization parameter with sufficient low noise of aerosol parameters
        l = 0
        do i = 0, naer-1
           if (sqrt(s_x(nx-i, nx-i)) < upbd(nx-i)-lowbd(nx-i) ) l = l + 1
        enddo
        !*** Take min regularization paramater for which aerosol noise is within allowed range
        if (l==naer .or. naer==0) exit
     enddo

      do i=1,nx
         u_ylsq(i) = dot_product(u_mat(:,i),ylsq)
         u_ylin(i) = dot_product(u_mat(:,i),ylin)
         v_x0(i)   = dot_product(v_mat(:,i),x_0)
      enddo

      !*** Calculate filter function
      filt_func_lsq = 1./(1.+lambda)
      filt_func_x0  = reg_par**2/(sing_val**2+reg_par**2)

      !*** Calculate new state vector
      x_i=0.D0
      filt_mat = 0.D0
      do i=1,nx
         filt_mat(i,i) = filt_func_pt(i)
         if (sing_val(i)/sing_val(1) .gt. epsilon) then
            x_i = x_i + &
                 (filt_func_pt(i)*u_ylin(i)/sing_val(i) + &
                 filt_func_lsq(i)*filt_func_pt(i)*u_ylsq(i)/sing_val(i) + &
                 filt_func_lsq(i)*filt_func_x0(i)*v_x0(i))*v_mat(:,i)
         endif
      enddo

      !*** Calculate averaging kernel
      a_avg = matmul(v_mat,matmul(filt_mat,transpose(v_mat)))

      !*** Add a apriori, x=A*xtrue+(1-A)*xapr
      unit_avg = -a_avg
      do i = 1, nx
         unit_avg(i,i) = 1.-a_avg(i,i)
      enddo
      x_sm = matmul(unit_avg,x_apr_tmp-x_0)
      do i = 1, nx
         x_i(i) = x_i(i) + x_sm(i)
      enddo

      !*** total degree of freedom
      dfs = 0.D0
      do k = 1, nx
         dfs = dfs + a_avg(k,k)
      enddo
      !*** degree of freedom of target absorbers
      dfs_target = 0.D0
      do j = 1, ntype_target
         rdum = 0.D0
         do k = 1, nx
            if(regskill(k)==1+j)then
               rdum = rdum + a_avg(k,k)
            endif
         enddo
         dfs_target(j) = rdum
      enddo
      !*** degree of freedom of aerosol parameters
      dfs_scat = 0.D0
      rdum = 0.D0
      do k = 1, nx
         if(regskill(k)==1) then
            rdum = rdum+a_avg(k,k)
         endif
      enddo
      dfs_scat = rdum

      !*** Rescale retrieval parameters
      x_i=x_i*weight
      x_0=x_0*weight
!      x_apr=x_apr*weight
      forall (i=1:nx,j=1:nx) a_avg(i,j) = a_avg(i,j)*weight(i)/weight(j)
      forall (i=1:nx,j=1:ny) d_mat(i,j) = d_mat(i,j)*weight(i)/DSQRT(s_y(j))

!*** Check for hitting the boundaries
      do i=1,nx
         if(x_i(i)<lowbd(i))then
            x_i(i) = lowbd(i)
            if(i>nx-naer .and. Boundary_Flag==0) then
               Boundary_Flag=1
            else
               Boundary_Flag=2
            endif
         elseif(x_i(i)>upbd(i)) then
            x_i(i) = upbd(i)
            if(i>nx-naer .and. Boundary_Flag==0) then
               Boundary_Flag = 1
            else
               Boundary_Flag = 2
            endif
	 endif
      enddo

  end subroutine column_inversion

!-------------------------------------------------------------------------

end module pt_regularization_module
