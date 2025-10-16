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
      kmat, ymeas, ymod, s_y_vector, &
      x_i, s_x, x_apr, &
      a_avg, d_mat, lambda, dfs, dfs_target, dfs_scat, upbd, lowbd, invflag, Boundary_Flag, line_number, flag_inv)

      implicit none
      integer, intent(in) :: ntype_target, naer
      integer, intent(in) :: line_number  ! hack
      integer, intent(in) :: flag_inv ! hack. 4: strong_co2 + strong_ch4, 5: strong_co2, 6: strong_ch4
      integer, intent(out) :: Boundary_Flag
      integer, intent(out) :: invflag
      integer :: i, j, k, l
      integer :: ny, nx, nlay
      real(double) :: lambda, dfs, dfs_scat, rdum
      real(double), dimension(ntype_target) :: dfs_target
      real(double), dimension(ny, nx) :: kmat
      real(double), dimension(ny) :: ymeas, ymod, yres
      real(double), dimension(ny) :: s_y_vector
      real(double), dimension(ny, ny) :: s_y_inv
      real(double), dimension(nx) :: x_i, x_apr, x_sm, upbd, lowbd
      real(double), dimension(nx, nx) :: s_x, s_x_inv, a_avg, unit_avg
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

      !*** Get inverse of measurement covariance matrix
      ! print*, "DEVELOPMENT: GET INVERSE OF COVARIANCE MATRIX"
      ! print*, "USE HARDCODED FILE FOR COV_INV"
      call read_nc_s_y_inv("CONTRL_OUT/MTF_OUT_DATA.nc", s_y_inv, ny, line_number, flag_inv)

      !*** Regularization matrix H
      call REGU_PAR(ntype_target, naer, nlay, nx, reg, H)

      !*** Calculate s_x_inv
      !*** s_x_inv = (KT*K+g*H)-1
      !*** s_x_inv := rodgers (3.31) S_hat, noise error plus smoothing component
      call inverse_lu(matmul(matmul(transpose(kmat), s_y_inv), kmat) + H, nx, s_x_inv, invflag)

      !*** Next iteration state vector
      x_i = x_i + 1./(1.+lambda)*matmul(s_x_inv, (matmul(matmul(transpose(kmat), s_y_inv), yres) - matmul(H, (x_i - x_0))))

      !*** Contribution function matrix
      d_mat = matmul(matmul(s_x_inv, transpose(kmat)), s_y_inv)

      !*** Averaging kernel matrix
      a_avg = matmul(d_mat, kmat)

      !*** Replace s_x by noise error, since smoothing is already described by averaging kernel
      !*** s_x = d_mat * s_y * d_mat^T
      !*** plug in d_mat = s_x_inv * kmat^T * s_y_inv
      !*** --> s_x = s_x_inv * kmat^T * d_mat^T
      !*** s_x := rodgers (3.19) S_m
      s_x = matmul(matmul(s_x_inv, transpose(kmat)), transpose(d_mat))

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
      forall (i=1:nx, j=1:ny) d_mat(i, j) = d_mat(i, j)*weight(i)

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

   subroutine read_nc_s_y_inv(filepath, cov_inv_y, ny, line_number, flag_inv)
      use netcdf
      !*** input
      character(len=*), intent(in) :: filepath
      integer, intent(in) :: ny ! length of part of state vector designated for target absorbers
      integer, intent(in) :: flag_inv ! hack. determines which target absorbers are used and how long vectors need to be
      integer, intent(in) :: line_number  ! hack (within a hack)
      !*** output
      real(double), dimension(ny, ny), intent(out) :: cov_inv_y
      !*** local
      real(double), dimension(:, :), allocatable :: cov_inv_co2
      real(double), dimension(:, :), allocatable :: cov_inv_ch4
      integer ncid, varid
      integer n_co2, n_ch4
      integer :: start3d(3)
      integer ierr
      !***
      ! Very hacky routine. If it is still here after I am gone, I am terribly sorry.

      call check(nf90_open(trim(filepath), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) then
         print*, "ERROR IN READ_NC_S_Y_INV (HACK): Error opening MTF_DATA_OUT.nc"
      end if

      call check(nf90_inq_dimid(ncid, "wavelength1_co2", varid), ierr)
      call check(nf90_inquire_dimension(ncid, varid, len=n_co2), ierr)
      call check(nf90_inq_dimid(ncid, "wavelength1_ch4", varid), ierr)
      call check(nf90_inquire_dimension(ncid, varid, len=n_ch4), ierr)

      ! location of data in netcdf file
      start3d = (/1, 1, line_number/) ! (channel, channel, line)

      if (allocated(cov_inv_co2)) deallocate(cov_inv_co2)
      allocate(cov_inv_co2(n_co2, n_co2))
      call check(nf90_inq_varid(ncid, "cov_inv_co2", varid), ierr)
      call check(nf90_get_var(ncid, varid, cov_inv_co2, start3d), ierr)

      if (allocated(cov_inv_ch4)) deallocate(cov_inv_ch4)
      allocate(cov_inv_ch4(n_ch4, n_ch4))
      call check(nf90_inq_varid(ncid, "cov_inv_ch4", varid), ierr)
      call check(nf90_get_var(ncid, varid, cov_inv_ch4, start3d), ierr)

      call check(nf90_close(ncid), ierr)

      call populate_cov_inv_from_file(ny, n_co2, n_ch4, cov_inv_co2, cov_inv_ch4, cov_inv_y, flag_inv)

      ! print*, "DEBUG:"
      ! print*, "line_number = ", line_number
      ! print*, "cov_inv_y = "
      ! do ierr = 1, 3
      !    print*, cov_inv_y(ierr, 1:4)
      ! end do
      ! print*, "shape(cov_inv_y) = "
      ! print*, shape(cov_inv_y)
   end subroutine read_nc_s_y_inv

   subroutine populate_cov_inv_from_file(ny, n_co2, n_ch4, cov_inv_co2, cov_inv_ch4, cov_inv_y, flag_inv)
      ! input
      integer, intent(in) :: ny, n_co2, n_ch4, flag_inv
      real(double), dimension(:, :), intent(in) :: cov_inv_co2
      real(double), dimension(:, :), intent(in) :: cov_inv_ch4
      ! output
      real(double), dimension(ny, ny), intent(out) :: cov_inv_y
      ! local
      integer :: i
      integer :: LEGAL_CO2, LEGAL_CH4, LEGAL_WEAK
      real(double) :: ch4_offdiagonal_scaling

      ! flag_inv determines size of the state vector elements
      ! strong_co2 has length 12
      ! strong_ch4 has length 35 or 36 (depending on drifted 2400 nm channel)
      ! weak_both counts for both gases and has length 19
      ! the weak band is inside the co2 and ch4 vectors so it has to be extracted specifically. doing that from co2 because it has the same length always
      ! for a sanity check, the total length ny is compared to the sum of the lenghts of the individual bands
      ! 4: strong_co2 + strong_ch4               12 + 35 = 47 or 12 + 36 = 48
      ! 5: strong_co2                            12 = 12
      ! 6: strong_ch4                            35 = 35 or 36 = 36
      ! 7: strong_co2 + strong_ch4 + weak_both   12 + 35 + 19 = 66 or 12 + 36 + 19 = 67, note: co2 has length (12+19=31), ch4 has length (35+19=54 or 36+19=55)
      ! 8: strong_co2 + weak_both                12 + 19 = 31, note: co2 has length (12+19=31)
      ! 9: strong_ch4 + weak_both                35 + 19 = 54 or 36 + 19 = 55, note: ch4 has length (35+19=54 or 36+19=55)

      ! ch4_offdiagonal_scaling is necessary for numerical reasons.
      ! matrix_inversion does not converge without it.
      ! currently only implemented for strong ch4 band
      ch4_offdiagonal_scaling = 1.05

      ! trying fix: shorten ch4 window for numerical stability
      ! 35 -> 28
      ! 36 -> 29
      LEGAL_CO2 = 12
      LEGAL_CH4 = 13
      LEGAL_WEAK = 19

      cov_inv_y = 0
      if (flag_inv .eq. 4) then ! strong_co2 + strong_ch4
         ! print*, "strong co2 + strong ch4"
         ! print*, n_co2, LEGAL_CO2
         ! print*, n_ch4, LEGAL_CH4
         ! print*, ny, LEGAL_CO2 + LEGAL_CH4
         if (n_co2 .eq. LEGAL_CO2 .and. n_ch4 .eq. LEGAL_CH4 .and. ny .eq. LEGAL_CO2 + LEGAL_CH4) then
            ! print*, "case 1"
            cov_inv_y(1:n_co2, 1:n_co2) = cov_inv_co2(1:n_co2, 1:n_co2)
            cov_inv_y(n_co2+1:n_co2+n_ch4, n_co2+1:n_co2+n_ch4) = cov_inv_ch4(1:n_ch4, 1:n_ch4) * ch4_offdiagonal_scaling
            do i = n_co2+1, n_co2+n_ch4
               cov_inv_y(i, i) = cov_inv_y(i, i) / ch4_offdiagonal_scaling
            end do
         else if (n_co2 .eq. LEGAL_CO2 .and. n_ch4 .eq. LEGAL_CH4+1 .and. ny .eq. LEGAL_CO2+LEGAL_CH4+1) then
            ! print*, "case 2"
            cov_inv_y(1:n_co2, 1:n_co2) = cov_inv_co2(1:n_co2, 1:n_co2)
            cov_inv_y(n_co2+1:n_co2+n_ch4, n_co2+1:n_co2+n_ch4) = cov_inv_ch4(1:n_ch4, 1:n_ch4) * ch4_offdiagonal_scaling
            do i = n_co2+1, n_co2+n_ch4
               cov_inv_y(i, i) = cov_inv_y(i, i) / ch4_offdiagonal_scaling
            end do
         ! else
            ! print*, "case 3"
         end if
      else if (flag_inv .eq. 5) then ! strong_co2
         if (n_co2 .eq. LEGAL_CO2 .and. ny .eq. LEGAL_CO2) then
            cov_inv_y(1:n_co2, 1:n_co2) = cov_inv_co2(1:n_co2, 1:n_co2)
         end if
      else if (flag_inv .eq. 6) then ! strong_ch4
         if (n_ch4 .eq. LEGAL_CH4 .and. ny .eq. LEGAL_CH4) then
            cov_inv_y(1:LEGAL_CH4, 1:LEGAL_CH4) = cov_inv_ch4(1:LEGAL_CH4, 1:LEGAL_CH4) * ch4_offdiagonal_scaling
            do i = 1, LEGAL_CH4
               cov_inv_y(i, i) = cov_inv_y(i, i) / ch4_offdiagonal_scaling
            end do
         else if (n_ch4 .eq. LEGAL_CH4+1 .and. ny .eq. LEGAL_CH4+1) then
            cov_inv_y(1:n_ch4, 1:n_ch4) = cov_inv_ch4(1:n_ch4, 1:n_ch4) * ch4_offdiagonal_scaling
            do i = 1, n_ch4
               cov_inv_y(i, i) = cov_inv_y(i, i) / ch4_offdiagonal_scaling
            end do
         end if
      else if (flag_inv .eq. 7) then ! strong_co2 + strong_ch4 + weak_both
         if (n_co2 .eq. LEGAL_CO2+LEGAL_WEAK .and. n_ch4 .eq. LEGAL_CH4+LEGAL_WEAK .and. ny .eq. LEGAL_CO2+LEGAL_CH4+LEGAL_WEAK) then
            cov_inv_y(1:LEGAL_CO2, 1:LEGAL_CO2) = cov_inv_co2(1:LEGAL_CO2, 1:LEGAL_CO2)
            cov_inv_y(LEGAL_CO2+1:LEGAL_CO2+LEGAL_CH4, LEGAL_CO2+1:LEGAL_CO2+LEGAL_CH4) = cov_inv_ch4(1:LEGAL_CH4, 1:LEGAL_CH4) * ch4_offdiagonal_scaling
            do i = LEGAL_CO2+1, LEGAL_CO2+LEGAL_CH4
               cov_inv_y(i, i) = cov_inv_y(i, i) / ch4_offdiagonal_scaling
            end do
            cov_inv_y(LEGAL_CO2+LEGAL_CH4+1:LEGAL_CO2+LEGAL_CH4+LEGAL_WEAK, LEGAL_CO2+LEGAL_CH4+1:LEGAL_CO2+LEGAL_CH4+LEGAL_WEAK) = cov_inv_co2(LEGAL_CO2+1:LEGAL_WEAK, LEGAL_CO2+1:LEGAL_WEAK)
         else if (n_co2 .eq. LEGAL_CO2+LEGAL_WEAK .and. n_ch4 .eq. LEGAL_CO2+LEGAL_CH4+1 .and. ny .eq. LEGAL_CO2+LEGAL_CH4+1+LEGAL_WEAK) then
            cov_inv_y(1:LEGAL_CO2, 1:LEGAL_CO2) = cov_inv_co2(1:LEGAL_CO2, 1:LEGAL_CO2)
            cov_inv_y(LEGAL_CO2+1:LEGAL_CH4+LEGAL_CH4+1, LEGAL_CO2+1:LEGAL_CO2+LEGAL_CH4+1) = cov_inv_ch4(1:LEGAL_CH4+1, 1:LEGAL_CH4+1) * ch4_offdiagonal_scaling
            do i = LEGAL_CO2+1, LEGAL_CO2+LEGAL_CH4+1
               cov_inv_y(i, i) = cov_inv_y(i, i) / ch4_offdiagonal_scaling
            end do
            cov_inv_y(LEGAL_CO2+LEGAL_CH4+1+1:LEGAL_CO2+LEGAL_CH4+1+LEGAL_WEAK, LEGAL_CO2+LEGAL_CH4+1+1:LEGAL_CO2+LEGAL_CH4+1+LEGAL_WEAK) = cov_inv_co2(LEGAL_CO2+1:LEGAL_CO2+LEGAL_WEAK, LEGAL_CO2+1:LEGAL_CO2+LEGAL_WEAK)
         end if
      else if (flag_inv .eq. 8) then ! strong_co2 + weak_both
         if (n_co2 .eq. LEGAL_CO2+LEGAL_WEAK .and. ny .eq. LEGAL_CO2+LEGAL_WEAK) then
            cov_inv_y(1:LEGAL_CO2, 1:LEGAL_CO2) = cov_inv_co2(1:LEGAL_CO2, 1:LEGAL_CO2)
            cov_inv_y(LEGAL_CO2+1:LEGAL_CO2+LEGAL_WEAK, LEGAL_CO2+1:LEGAL_CO2+LEGAL_WEAK) = cov_inv_co2(LEGAL_CO2+1:LEGAL_CO2+LEGAL_WEAK, LEGAL_CO2+1:LEGAL_WEAK)
         end if
      else if (flag_inv .eq. 9) then ! strong_ch4 + weak_both
         if (n_ch4 .eq. LEGAL_CH4+LEGAL_WEAK .and. ny .eq. LEGAL_CH4+LEGAL_WEAK) then
            cov_inv_y(1:LEGAL_CH4, 1:LEGAL_CH4) = cov_inv_ch4(1:LEGAL_CH4, 1:LEGAL_CH4) * ch4_offdiagonal_scaling
            do i = 1, LEGAL_CH4
               cov_inv_y(i, i) = cov_inv_y(i, i) / ch4_offdiagonal_scaling
            end do
            cov_inv_y(LEGAL_CH4+1:LEGAL_CH4+LEGAL_WEAK, LEGAL_CH4+1:LEGAL_CH4+LEGAL_WEAK) = cov_inv_co2(1:LEGAL_WEAK, 1:LEGAL_WEAK)
         else if (n_ch4 .eq. LEGAL_CH4+1+LEGAL_WEAK .and. ny .eq. LEGAL_CH4+1+LEGAL_WEAK) then
            cov_inv_y(1:LEGAL_CH4+1, 1:LEGAL_CH4+1) = cov_inv_ch4(1:LEGAL_CH4+1, 1:LEGAL_CH4+1) * ch4_offdiagonal_scaling
            do i = 1, LEGAL_CH4+1
               cov_inv_y(i, i) = cov_inv_y(i, i) / ch4_offdiagonal_scaling
            end do
            cov_inv_y(LEGAL_CH4+1+1:LEGAL_CH4+1+LEGAL_WEAK, LEGAL_CH4+1+1:LEGAL_CH4+1+LEGAL_WEAK) = cov_inv_co2(1:LEGAL_WEAK, 1:LEGAL_WEAK)
         end if
      end if
      if (cov_inv_y(1, 1) .eq. 0) then
         print*, "ERROR POPULATE_COV_INV_FROM_FILE: First element of covariance matrix is zero, something probably went wrong. Stopping..."
         stop
      else if (cov_inv_y(ny, ny) .eq. 0) then
         print*, "ERROR POPULATE_COV_INV_FROM_FILE: Last element of covariance matrix is zero, something probably went wrong. Stopping..."
         stop
      end if
   end subroutine populate_cov_inv_from_file

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
