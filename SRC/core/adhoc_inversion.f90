module adhoc_inversion_module
   use header_module
   use read_settings_module, only: regularization_class
   use auxiliary_routines_module
   implicit none
contains

   subroutine adhoc_inversion(ntype_target, naer, &
                              ny, nx, nlay, regskill, reg, &
                              kmat, ymeas, ymod, s_y, &
                              x_i, s_x, x_apr, &
                              a_avg, d_mat, lambda, dfs, dfs_target, dfs_scat, upbd, lowbd, invflag, Boundary_Flag)

      implicit none
      integer, intent(in) :: ntype_target, naer
      integer, intent(out) :: Boundary_Flag
      integer, intent(out) :: invflag
      integer :: i, j, k, l
      integer :: ny, nx, nlay
      real(double) :: lambda, dfs, dfs_scat, rdum
      real(double), dimension(ntype_target) :: dfs_target
      real(double), dimension(ny, nx) :: kmat
      real(double), dimension(ny) :: ymeas, ymod, y
      real(double), dimension(ny) :: s_y
      real(double), dimension(nx) :: x_i, x_apr, x_sm, upbd, lowbd
      real(double), dimension(nx, nx) :: s_x, a_avg, unit_avg
      real(double), dimension(nx, ny) :: d_mat
      integer, dimension(nx) :: regskill
      type(regularization_class), intent(in) :: reg
      real(double), dimension(nx, nx) :: H
      real(double), dimension(nx) :: x_0, weight
!*** Skill the inverse problem

      boundary_flag = 0

      call pt_precon( &
         ntype_target, ny, nx, regskill, &
         kmat, &
         x_i, weight)

      do l = 1, ny
         kmat(l, 1:nx) = kmat(l, 1:nx)*weight(1:nx)
      end do

      x_i = x_i/weight
      x_0 = x_apr/weight
      x_apr = x_apr/weight

      y = ymeas - ymod

      forall (l=1:ny)
         y(l) = y(l)/DSQRT(s_y(l))
         kmat(l, :) = kmat(l, :)/DSQRT(s_y(l))
      end forall

!*** Regularization matrix H

      call REGU_PAR(ntype_target, naer, nlay, nx, reg, H)

!   IF(runid==0)THEN
!      OPEN(50,FILE=TRIM(runpath)//'CONTRL_OUT/hmat.dat')
!      DO i=1,nx
!         WRITE(50,'(100(1pE13.5,x))')H(i,:)
!      ENDDO
!      CLOSE(50)
!   ENDIF

!*** Sx-1 = (KT*K+g*H)-1
      call inverse_lu(matmul(transpose(kmat), kmat) + H, nx, s_x, invflag)
!   s_x = INVERSE(matmul(transpose(kmat),kmat) + H, nx)

!*** Next iteration state vector
!   IF(iter<=nlsq+1 .and. lambda > 0.D0)lambda=MIN(dble(nx),10.)

      x_i = x_i + 1./(1.+lambda)*matmul(s_x, (matmul(transpose(kmat), y) - matmul(H, (x_i - x_0))))

!*** Contribution function matrix
      d_mat = matmul(s_x, transpose(kmat))

!*** Averaging kernel matrix
      a_avg = matmul(s_x, matmul(transpose(kmat), kmat))

!*** Replace s_x by noise error
      s_x = MATMUL(d_mat, TRANSPOSE(d_mat))

!*** Add a apriori, x=A*xtrue+(1-A)*xapr, if applicable
      unit_avg = -a_avg
      do i = 1, nx
         unit_avg(i, i) = 1.-a_avg(i, i)
      end do
      x_sm = matmul(unit_avg, x_apr - x_0)

      do k = 1, nx
         x_i(k) = x_i(k) + x_sm(k)
      end do

!*** Degrees of freedom

      dfs = 0.

      do k = 1, nx
         dfs = dfs + a_avg(k, k)
      end do

      do j = 1, ntype_target
         rdum = 0.D0
         do k = 1, nx
            if (regskill(k) == 1 + j) then
               rdum = rdum + a_avg(k, k)
            end if
         end do
         dfs_target(j) = rdum
      end do

      dfs_scat = 0.D0
      rdum = 0.D0
      do k = 1, nx
         if (regskill(k) == 1) then
            rdum = rdum + a_avg(k, k)
         end if
      end do
      dfs_scat = rdum

!*** Rescale retrieval parameters

      x_i = x_i*weight
      x_0 = x_0*weight
      x_apr = x_apr*weight

      forall (i=1:nx, j=1:nx) s_x(i, j) = s_x(i, j)*weight(i)*weight(j)
      forall (i=1:nx, j=1:nx) a_avg(i, j) = a_avg(i, j)*weight(i)/weight(j)
      forall (i=1:nx, j=1:ny) d_mat(i, j) = d_mat(i, j)*weight(i)/DSQRT(s_y(j))

!*** Check for hitting the boundaries
      do i = 1, nx
         if (x_i(i) < lowbd(i)) then
            x_i(i) = lowbd(i)
            if (i > nx - naer) then
               Boundary_Flag = 1
            else
               Boundary_Flag = 2
            end if
         elseif (x_i(i) > upbd(i)) then
            x_i(i) = upbd(i)
            if (i > nx - naer) then
               Boundary_Flag = 1
            else
               Boundary_Flag = 2
            end if
         end if
      end do

   end subroutine adhoc_inversion

!*************************************************************************************

   subroutine adhoc_inversion_matrix(&
      ntype_target, naer, &
      ny, nx, nlay, regskill, reg, &
      kmat, ymeas, ymod, s_y, &
      x_i, s_x, x_apr, &
      a_avg, d_mat, lambda, dfs, dfs_target, dfs_scat, upbd, lowbd, invflag, Boundary_Flag)

      implicit none
      integer, intent(in) :: ntype_target, naer
      integer, intent(out) :: Boundary_Flag
      integer, intent(out) :: invflag
      integer :: i, j, k, l
      integer :: ny, nx, nlay
      real(double) :: lambda, dfs, dfs_scat, rdum
      real(double), dimension(ntype_target) :: dfs_target
      real(double), dimension(ny, nx) :: kmat
      real(double), dimension(ny) :: ymeas, ymod, yres
      real(double), dimension(ny) :: s_y
      real(double), dimension(nx) :: x_i, x_apr, x_sm, upbd, lowbd
      real(double), dimension(nx, nx) :: s_x, a_avg, unit_avg
      real(double), dimension(nx, ny) :: d_mat
      integer, dimension(nx) :: regskill
      type(regularization_class), intent(in) :: reg
      real(double), dimension(nx, nx) :: H
      real(double), dimension(nx) :: x_0, weight

      !*** Skill the inverse problem
      boundary_flag = 0

      call pt_precon( &
         ntype_target, ny, nx, regskill, &
         kmat, &
         x_i, weight)

      !*** Scale parameters
      do l = 1, ny
         kmat(l, 1:nx) = kmat(l, 1:nx)*weight(1:nx)
      end do

      x_i = x_i/weight
      x_0 = x_apr/weight
      x_apr = x_apr/weight

      !*** Calculate yres
      yres = ymeas - ymod

      !*** Pull covariance into variables
      forall (l=1:ny)
         yres(l) = yres(l)/DSQRT(s_y(l))
         kmat(l, :) = kmat(l, :)/DSQRT(s_y(l))
      end forall

      !*** Regularization matrix H
      call REGU_PAR(ntype_target, naer, nlay, nx, reg, H)

      ! IF(runid==0)THEN
      !    OPEN(50,FILE=TRIM(runpath)//'CONTRL_OUT/hmat.dat')
      !    DO i=1,nx
      !       WRITE(50,'(100(1pE13.5,x))')H(i,:)
      !    ENDDO
      !    CLOSE(50)
      ! ENDIF

      !*** Invert s_x
      !*** (s_x)-1 = (KT*K+g*H)-1
      call inverse_lu(matmul(transpose(kmat), kmat) + H, nx, s_x, invflag)

      !*** Next iteration state vector
      x_i = x_i + 1./(1.+lambda)*matmul(s_x, (matmul(transpose(kmat), yres) - matmul(H, (x_i - x_0))))

      !*** Contribution function matrix
      d_mat = matmul(s_x, transpose(kmat))

      !*** Averaging kernel matrix
      a_avg = matmul(s_x, matmul(transpose(kmat), kmat))

      !*** Replace s_x by noise error
      s_x = MATMUL(d_mat, TRANSPOSE(d_mat))

      !*** Add apriori, x=A*xtrue+(1-A)*xapr, if applicable
      unit_avg = -a_avg
      do i = 1, nx
         unit_avg(i, i) = 1.-a_avg(i, i)
      end do
      x_sm = matmul(unit_avg, x_apr - x_0)

      do k = 1, nx
         x_i(k) = x_i(k) + x_sm(k)
      end do

      !*** Degrees of freedom
      dfs = 0.

      do k = 1, nx
         dfs = dfs + a_avg(k, k)
      end do

      do j = 1, ntype_target
         rdum = 0.D0
         do k = 1, nx
            if (regskill(k) == 1 + j) then
               rdum = rdum + a_avg(k, k)
            end if
         end do
         dfs_target(j) = rdum
      end do

      dfs_scat = 0.D0
      rdum = 0.D0
      do k = 1, nx
         if (regskill(k) == 1) then
            rdum = rdum + a_avg(k, k)
         end if
      end do
      dfs_scat = rdum

      !*** Rescale retrieval parameters
      x_i = x_i*weight
      x_0 = x_0*weight
      x_apr = x_apr*weight

      forall (i=1:nx, j=1:nx) s_x(i, j) = s_x(i, j)*weight(i)*weight(j)
      forall (i=1:nx, j=1:nx) a_avg(i, j) = a_avg(i, j)*weight(i)/weight(j)
      forall (i=1:nx, j=1:ny) d_mat(i, j) = d_mat(i, j)*weight(i)/DSQRT(s_y(j))

      !*** Check for hitting the boundaries
      do i = 1, nx
         if (x_i(i) < lowbd(i)) then
            x_i(i) = lowbd(i)
            if (i > nx - naer) then
               Boundary_Flag = 1
            else
               Boundary_Flag = 2
            end if
         elseif (x_i(i) > upbd(i)) then
            x_i(i) = upbd(i)
            if (i > nx - naer) then
               Boundary_Flag = 1
            else
               Boundary_Flag = 2
            end if
         end if
      end do
   end subroutine adhoc_inversion_matrix

!*************************************************************************************

   subroutine regu_par(ntype_target, naer, nlay, &
                       nx, reg, hmat)
      integer, intent(in) :: ntype_target, naer, nlay
      integer :: j, j1, j2, order
      integer :: nx
      type(regularization_class), intent(in) :: reg
      real(double), dimension(nx, nx) :: hmat

      order = 1
      hmat = 0.
      do j = 1, ntype_target

         j1 = (j - 1)*nlay + 1
         j2 = j*nlay

         order = 1

         call REGU_MATRIX(nlay, order, hmat(j1:j2, j1:j2))
         hmat(j1:j2, j1:j2) = reg%weight_target(j)**2.0*hmat(j1:j2, j1:j2)

      end do

      j1 = nx - naer + 1
      j2 = nx

      order = 0

      call REGU_MATRIX(naer, order, hmat(j1:j2, j1:j2))
      hmat(j1:j2, j1:j2) = reg%weight_aerosol**2.0*hmat(j1:j2, j1:j2)

   end subroutine regu_par

!*************************************************************************************

   subroutine pt_precon( &
      ntype_target, ny, nx, regskill, &
      kmat, &
      x_i, weight)
      implicit none
      integer, intent(in) :: ntype_target
      integer :: j, k
      integer :: nx, ny
      integer, dimension(nx) :: regskill

      real(double), dimension(ny, nx) :: kmat
      real(double), dimension(nx) :: x_i, weight, x_i_new
      real(double) :: rdum, rdummax

!*********************************************************************

      !*** Skill non-profile parameters

      do k = 1, nx

         !rdum = 1./DSQRT(s_apr(k,k))
         !rdum = 1./MAX(DABS(x_apr(k)),1.D-20)
         rdum = maxval(DABS(kmat(1:ny, k)))
         if (regskill(k) == 0 .or. regskill(k) == 1) then
            if (rdum .ne. 0.D0) then
               weight(k) = 1./rdum
            else
               weight(k) = 1.D0
            end if
         end if

      end do

      !*** Skill profile parameters

      do j = 1, ntype_target

         rdummax = 0
         do k = 1, nx
            !rdum = 1./MAX(DABS(x_apr(k)),1.D-20)
            !rdum = 1./DSQRT(s_apr(k,k))
            rdum = maxval(DABS(kmat(1:ny, k)))
            if (regskill(k) == 1 + j .and. rdum > rdummax) then
               rdummax = rdum
            end if
         end do

         do k = 1, nx
            if (regskill(k) == 1 + j) then
               if (rdummax .ne. 0.D0) then
                  weight(k) = 1./rdummax
               else
                  weight(k) = 1.D0
               end if
            end if
         end do

      end do

      !*** Normalize regularization parameters of subproblems
      x_i_new = x_i/weight
      do j = 1, ntype_target + 1

         rdummax = 0.
         do k = 1, nx
            rdum = DABS(x_i_new(k))
            if (regskill(k) == j .and. rdum > rdummax) then
               rdummax = rdum
            end if
         end do

         do k = 1, nx
            if (regskill(k) == j) then
               if (rdummax .ne. 0.D0) then
                  weight(k) = weight(k)*rdummax
               end if
            end if
         end do

      end do

   end subroutine pt_precon

!**************************************************************************************

   subroutine regu_matrix(nx, o, H)
      implicit none

      integer :: i

      integer :: nx, o

      real(double), dimension(nx, nx) :: H

      H = 0.D0

      !*** Zeroth order
      if (o == 0) then

         do i = 1, nx
            H(i, i) = 1.D0
         end do

         !*** First order
      elseif (o == 1 .and. nx > 2) then

         H(1, 1) = 1.D0
         H(1, 2) = -1.D0
         do i = 2, nx - 1
            H(i, i - 1) = -1.D0
            H(i, i) = 2.D0
            H(i, i + 1) = -1.D0
         end do
         H(nx, nx) = 1.D0
         H(nx, nx - 1) = -1.D0

         !*** Second order
      elseif (o == 2 .and. nx > 4) then

         H(1, 1) = 1.D0
         H(1, 2) = -2.D0
         H(1, 3) = 1.D0
         H(2, 1) = -2.D0
         H(2, 2) = 5.D0
         H(2, 3) = -4.D0
         H(2, 4) = 1.D0
         do i = 3, nx - 2
            H(i, i - 2) = 1.D0
            H(i, i - 1) = -4.D0
            H(i, i) = 6.D0
            H(i, i + 1) = -4.D0
            H(i, i + 2) = 1.D0
         end do
         H(nx - 1, nx - 3) = 1.D0
         H(nx - 1, nx - 2) = -4.D0
         H(nx - 1, nx - 1) = 5.D0
         H(nx - 1, nx) = -2.D0
         H(nx, nx - 2) = 1.D0
         H(nx, nx - 1) = -2.D0
         H(nx, nx) = 1.D0
      else
         call stopretrieval('REGU_MATRIX: UNKNOWN REGULARIZATION ORDER')
      end if

   end subroutine regu_matrix

end module adhoc_inversion_module
