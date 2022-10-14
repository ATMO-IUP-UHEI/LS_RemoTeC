!------------------------------------------------------------------------------
!> @brief Compute optical properties of aerosol particles
!------------------------------------------------------------------------------
module OpticM_module
  use header_module 
  use auxiliary_routines_module, only: intrpl_spline
  implicit none
  private

  !*** Public procedures
  public :: read_aerosol_netcdf, read_cross, read_F, read_x0, modes_calc_xs, modes_calc_ssc, modes_calc

  !*** Private procedures
  private :: RemIm_index, VolDistrFunc, SpEl_exp_xs, C_read_int, gauleg, F_read_int, &
       Bilin_interp, Spl_interp, devel, scatmat_ssc

  !------------------------------------------------------------------------------
  !*** Parameters
  integer, parameter :: n_xgr0 = 41  !< number of size bins used for fixed kernel calculation
  integer, parameter :: n_tet0 = 181 !< number of gridpoints for scattering angle
  integer, parameter :: n_Rem0 = 22  !< number of gridpoints for real refractive index
  integer, parameter :: n_Imm0 = 16  !< number of gridpoints for imaginary refractive index
  !*** Indices to convert array to scattering matrix
  integer, dimension(6) :: index_ist_sp = (/1,2,2,3,4,4 /)
  integer, dimension(6) :: index_jst_sp = (/1,1,2,3,3,4 /)

  !------------------------------------------------------------------------------
  !> @brief aerosol optical properties LUT
  type, public :: Mie_lut
     private
     real(double), dimension(n_xgr0, n_Rem0, n_Imm0) :: C_ext0_in_el  !< extinction cross section for ellipsoidal particle
     real(double), dimension(n_xgr0, n_Rem0, n_Imm0) :: C_abs0_in_el  !< absorption cross section for ellipsoidal particle
     real(double), dimension(n_xgr0, n_Rem0, n_Imm0) :: C_ext0_in_sph !< extinction cross section for spherical particle
     real(double), dimension(n_xgr0, n_Rem0, n_Imm0) :: C_abs0_in_sph !< absorption cross section for spherical particle
     real(double), allocatable, dimension(:, :, :, :, :) ::  M_el0_in_el  !< scattering phase function for ellipsoidal particle -- on LUT theta grid
     real(double), allocatable, dimension(:, :, :, :, :) ::  M_el0_in_sph !< scattering phase function for spherical particle   -- on LUT theta grid
     real(double), allocatable, dimension(:, :, :, :, :) ::  M_el0_el  !< scattering phase function for ellipsoidal particle    -- on quadrature theta grid
     real(double), allocatable, dimension(:, :, :, :, :) ::  M_el0_sph !< scattering phase function for spherical particle     -- on quadrature theta grid
     real(double), dimension(n_tet0) :: tet0_dg                       !< scattering angle -- quadrature
     real(double), dimension(n_tet0) :: mu_t_dg                       !< COS(scattering)  -- quadrature 
     real(double), dimension(n_tet0) :: Gw_t_dg                       !< weights          -- quadrature 
     real(double), dimension(n_tet0) :: tet0_in                       !< scattering angle -- LUT
     real(double), dimension(n_xgr0) :: x0_in                         !< size parameter (\f$=\frac{2\pi r}{\lambda}\f$)
     real(double) :: wave_ref_lut                                     !< reference wavelength
     real(double), dimension(n_Rem0):: Rem_in                         !< real refractive index
     real(double), dimension(n_Imm0) :: Imm_in                        !< imaginary refractive index
  end type Mie_lut
  !------------------------------------------------------------------------------

contains
  !------------------------------------------------------------------------------
  !> Load NetCDF format aerosol look-up tables into memory
  !------------------------------------------------------------------------------
  subroutine read_aerosol_netcdf(aerosol_file, outputflag, aero_lut, ierr)
    use netcdf
    !*** input
    integer, intent(in) :: outputflag
    character(len=*), intent(in) :: aerosol_file   
    !*** Output
    type(Mie_lut), intent(out) :: aero_lut
    integer, intent(out) :: ierr
    !*** local variables
    integer :: i, j, k, l, n, ncid, varID, dimID, n_xgr01, n_tet01, n_Rem01, n_Imm01
    real(double) :: theta(n_tet0), r(n_xgr0)
    real(double), allocatable :: m(:, :, :, :, :)
    character(stringlen) :: message
    real(double) :: x_max 	 
    integer      :: n_tet, i_re, i_im, i_x 


    if(outputflag >=2 )then
       write(message,'(a)') '*** Start read_aerosol_netcdf ***'
       call writelog(message, 1)
    endif

    !*** Open NetCDF LUT 
    call check(nf90_open(trim(aerosol_file), nf90_nowrite, ncid), ierr)
    if (ierr .ne. 0) return

    !*** Get grid of size parameters
    call check(NF90_INQ_DIMID(ncid, "eff_radius", dimID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQUIRE_DIMENSION(ncid, dimID, len = n_xgr01), ierr)
    if (ierr .ne. 0) return
    if (n_xgr0 .ne. n_xgr01) then
       ierr = ierr_var
       call stopretrieval('READ_AEROSOL_NETCDF: n_xgr0 .ne. n_xgr01')     
       return
    endif
    call check(NF90_INQ_VARID(ncid, "eff_radius", varID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid, varID, r), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQ_VARID(ncid, "wavelength", varID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid, varID, aero_lut%wave_ref_lut), ierr)
    if (ierr .ne. 0) return
    aero_lut%x0_in = 2*pi*r/aero_lut%wave_ref_lut
    !*** Get grid of scattering angles
    call check(NF90_INQ_DIMID(ncid, "scat_angle", dimID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQUIRE_DIMENSION(ncid, dimID, len = n_tet01), ierr)
    if (ierr .ne. 0) return
    if (n_tet01 .ne. n_tet0) then
       ierr = ierr_var
       call stopretrieval('READ_AEROSOL_NETCDF: n_tet0 .ne. n_tet01')    
       return
    endif
    call check(NF90_INQ_VARID(ncid, "scat_angle", varID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid, varID, theta), ierr)
    if (ierr .ne. 0) return
    aero_lut%tet0_in = theta*pi/180.d0     !convert from deg to rad
    !*** Get grid of refractive indices
    !*** Real
    call check(NF90_INQ_DIMID(ncid, "ref_index_real", dimID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQUIRE_DIMENSION(ncid, dimID, len = n_Rem01), ierr)
    if (ierr .ne. 0) return
    if (n_Rem01 .ne. n_Rem0) then 
       ierr = ierr_var
       call stopretrieval('READ_AEROSOL_NETCDF: n_Rem0 .ne. n_Rem01')  
       return  
    endif
    call check(NF90_INQ_VARID(ncid, "ref_index_real", varID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid, varID,aero_lut%Rem_in), ierr)
    if (ierr .ne. 0) return
    !*** Imaginary
    call check(NF90_INQ_DIMID(ncid, "ref_index_im", dimID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQUIRE_DIMENSION(ncid, dimID, len = n_Imm01), ierr)
    if (ierr .ne. 0) return
    if (n_Imm01 .ne. n_Imm0) then 
       ierr = ierr_var
       call stopretrieval('READ_AEROSOL_NETCDF: n_Imm0 .ne. n_Imm01') 
       return  
    endif
    call check(NF90_INQ_VARID(ncid, "ref_index_im", varID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid, varID,aero_lut%Imm_in), ierr)
    if (ierr .ne. 0) return
    do i = 1, n_Imm0
       if (aero_lut%Imm_in(i) .lt. 0.d0) aero_lut%Imm_in(i) = - aero_lut%Imm_in(i)
    enddo
    !*** Get extinction cross sections
    !*** spheroids
    call check(NF90_INQ_VARID(ncid, "ext_coef_sph", varID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid, varID, aero_lut%C_ext0_in_sph), ierr)
    if (ierr .ne. 0) return
    !*** ellipsoids      
    call check(NF90_INQ_VARID(ncid, "ext_coef_el", varID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid, varID, aero_lut%C_ext0_in_el), ierr)
    if (ierr .ne. 0) return
    !*** Get absorption cross sections
    !*** spheroids
    call check(NF90_INQ_VARID(ncid, "abs_coef_sph", varID), ierr)
    call check(NF90_GET_VAR(ncid, varID, aero_lut%C_abs0_in_sph), ierr)

    !*** ellipsoids      
    call check(NF90_INQ_VARID(ncid, "abs_coef_el", varID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid, varID, aero_lut%C_abs0_in_el), ierr) 
    if (ierr .ne. 0) return
    !*** Get scattering matrix
    !*** spheroids
    call allocate_scattering(aero_lut, ierr)
    if (ierr.ne.0) then
       ierr = ierr_all
       write(message,*) 'READ_AEROSOL_NETCDF: memory allocation error'
       goto 999
    endif

    call check(NF90_INQ_VARID(ncid, "scat_mat_sph", varID), ierr)
    if (ierr .ne. 0) return
    allocate(m(n_tet0, n_xgr0, n_Rem0, n_Imm0, 6), stat=ierr)
    if (ierr.ne.0) then
       ierr = ierr_all
       write(message,*) 'READ_AEROSOL_NETCDF: memory allocation error'
       goto 999
    endif
    call check(NF90_GET_VAR(ncid, varID, m), ierr)
    if (ierr .ne. 0) return
    do n = 1, n_el
       do j = 1, n_Imm0
          do i = 1, n_Rem0
             do k = 1, n_xgr0 
                do l = 1, n_tet0
                   aero_lut%M_el0_in_sph(l, k, i, j, n) = m(l, k, i, j, n)
                enddo
             enddo
          enddo
       enddo
    enddo
    !*** ellipsoids
    call check(NF90_INQ_VARID(ncid, "scat_mat_el", varID), ierr)
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid, varID, m), ierr)
    if (ierr .ne. 0) return
    do n = 1, n_el
       do j = 1, n_Imm0
          do i = 1, n_Rem0
             do k = 1, n_xgr0 
                do l = 1, n_tet0
                   aero_lut%M_el0_in_el(l, k, i, j, n) = m(l, k, i, j, n)
                enddo
             enddo
          enddo
       enddo
    enddo

    !*** Create Quadrature grid to interpolate the angle-resolved phase funcions onto. This is aimed at easy numerical interpolation when calculation expansion coefs.
    x_max = aero_lut%x0_in(n_xgr0)						
    n_tet = int(x_max+dsqrt(x_max)+2)					
    n_tet = n_tet0								
    call gauleg(n_tet,n_tet,-1.d0,1.d0, aero_lut%mu_t_dg, aero_lut%Gw_t_dg)
    do i = 1, n_tet							
       aero_lut%tet0_dg(i) = dacos(aero_lut%mu_t_dg(i))		
    end do

    !***Interpolate ange-resolved phasefunctions to Quadrature angles -- 	
    if (nstokes==1) then
       do i_re=1,n_Rem0       !!!!!! Re_m
          do i_im=1,n_Imm0      !!!!!! Im_m
             do  i_x=1,n_xgr0
                !use spline interpolation to interpolate phase function in theta 
                call Spl_interp(n_el, n_tet0, n_tet, aero_lut%tet0_in, aero_lut%tet0_dg, aero_lut%M_el0_in_sph(:,i_x,i_re, i_im, :), aero_lut%M_el0_sph(:,i_x ,i_re, i_im, :),ierr )
                call Spl_interp(n_el, n_tet0, n_tet, aero_lut%tet0_in, aero_lut%tet0_dg, aero_lut%M_el0_in_el(:,i_x,i_re, i_im, :),  aero_lut%M_el0_el(:,i_x ,i_re, i_im, :),ierr)
             enddo
          enddo
       enddo
    else
       do i_re=1,n_Rem0      
          do i_im=1,n_Imm0 
             do  i_x=1,n_xgr0
                call Spl_interp_v6(n_el, n_tet0, n_tet, aero_lut%tet0_in, aero_lut%tet0_dg, aero_lut%M_el0_in_sph(:,i_x,i_re, i_im, :), aero_lut%M_el0_sph(:,i_x ,i_re, i_im, :),ierr )
                call Spl_interp_v6(n_el, n_tet0, n_tet, aero_lut%tet0_in, aero_lut%tet0_dg, aero_lut%M_el0_in_el(:,i_x,i_re, i_im, :),  aero_lut%M_el0_el(:,i_x ,i_re, i_im, :) ,ierr )
             enddo
          enddo
       enddo
    endif

    deallocate(m, stat=ierr)
    if (ierr.ne.0) then
       ierr = ierr_deall
       write(message,*) 'READ_AEROSOL_NETCDF: memory deallocation error'
       goto 999
    endif

    call check(nf90_close(ncid), ierr)
    if (ierr .ne. 0) then
       ierr = ierr_open
       write(message,*) 'READ_AEROSOL_NETCDF: error closing file'
       goto 999     
    endif

    if (outputflag >= 2) then
       write(message,'(a)') '*** End read_aerosol_netcdf ***'
       call writelog(message, 1)
    endif

    ierr = 0 ! Mark succes
    return


999 continue
    call stopretrieval(message)
    return

  end subroutine read_aerosol_netcdf

  !------------------------------------------------------------------------------
  !> Check calls to netCDF functions
  !------------------------------------------------------------------------------
  subroutine check(status, ierr)
    use NETCDF 
    integer, intent (in) :: status 
    integer, intent(out) :: ierr  

    if(status /= nf90_noerr) then
       ierr = ierr_read
       call stopretrieval('READ_AEROSOL_NETCDF: '//trim(nf90_strerror(status)))
    else
       ierr = 0
    endif

  end subroutine check

  !------------------------------------------------------------------------------
  !> Read grid of effective radii from ascii file
  !------------------------------------------------------------------------------
  subroutine read_x0(aerosol_dir, aero_lut)
    !*** Input
    character(len=*), intent(in) ::  aerosol_dir
    !*** Output
    type(Mie_lut), intent(out) :: aero_lut
    !*** Local variables
    character(80) :: full_name
    integer :: n_xgr01, i, io
    real(double) :: r, wlength
    character(15), parameter :: grid_x0 = 'grid1.dat.fix'
    !----------------------------------------------------  

    full_name = trim(aerosol_dir)//'Kernal_sph/'//trim(grid_x0)
    open(newunit(io),file=full_name)
    !      write(*,*) 'Opened aerosol file: ', trim(full_name), ' in remotec_core/OpticM_module_scalar/read_x0 (line 195)'
    read(io,*) n_xgr01, wlength
    aero_lut%wave_ref_lut = wlength
    if (n_xgr0 .ne. n_xgr01) call stopretrieval('Error in aerosol ascii LUT: read_x0') 

    do i=1, n_xgr0
       read(io,*) r
       aero_lut%x0_in(i) = 2*pi*r/wlength
    end do

    read(io,*) !n_tet0
    do i=1, n_tet0
       read(io,*) aero_lut%tet0_in(i)
       aero_lut%tet0_in(i) = aero_lut%tet0_in(i)*pi/180.d0   !convert from deg to rad
    end do
    close(io)

  end subroutine Read_x0

  !------------------------------------------------------------------------------
  !> Read absorption and exctinction cross-sections from ascii file
  !------------------------------------------------------------------------------
  subroutine read_cross(aerosol_dir, aero_lut)
    !*** Input
    character(len=*), intent(in) :: aerosol_dir
    !*** Output
    type(mie_lut), intent(inout) :: aero_lut
    !*** Local variables
    real(double), dimension(n_xgr0) :: C_ext0, C_abs0     
    real(double) :: Rem_min, Rem_max, Imm_min, Imm_max, Wavel     
    integer :: sp_el, n_xgr01, n_Rem, n_Imm, i_re, i_im, io
    character(80) :: full_name
    character(stringlen) :: dir_dat_el, dir_dat_sp
    character(15), parameter :: name_c='Rkext1.fix'
    !-----------------------------------------------------------

    dir_dat_sp = trim(aerosol_dir)//'Kernal_sph/'
    dir_dat_el = trim(aerosol_dir)//'Kernal_fix1/'
    do sp_el = 0, 1
       !*** For spherical particles
       if (sp_el .eq. 0) full_name = trim(dir_dat_sp)//trim(name_c)
       !*** For elipsoids
       if (sp_el .eq. 1) full_name = trim(dir_dat_el)//trim(name_c)
       open(newunit(io), file=full_name)
       !         write(*,*) 'Opened cross-section file: ', trim(full_name), ' in remotec_core/OpticM_module_scalar/read_cross (line 239)'
       read(io,*)
       read(io,*) n_xgr01
       if (abs(n_xgr01) .ne. n_xgr0) &
            call stopretrieval('READ_CROSS: n_xgr01 .ne. n_xgr0 ')
       read(io,*) Rem_min, Rem_max    !!!!!!!! min and max of Re_m
       read(io,*) Imm_min, Imm_max    !!!!!!!! min and max of Im_m
       read(io,*) n_Rem, n_Imm          !!!!!!!! number of interval for optical const
       if (n_Rem .lt. 0) n_Rem = -n_Rem
       if (n_Imm .lt. 0) n_Imm = -n_Imm 
       do i_re = 1, n_Rem0       !!!!!! Re_m
          do i_im = 1, n_imm0      !!!!!! Im_m
             read(io,*)    !element,ratio
             read(io,*) Wavel, aero_lut%Rem_in(i_re), aero_lut%Imm_in(i_im) !wavelength, Re_m, Im_m
             if (aero_lut%Imm_in(i_im) .lt. 0) aero_lut%Imm_in(i_im)=-aero_lut%Imm_in(i_im)
             read(io,*) 
             read(io,*) c_ext0(1:n_xgr0)
             read(io,*) 
             read(io,*) c_abs0(1:n_xgr0)
             if(sp_el.eq.0) then
                aero_lut%C_ext0_in_sph(:, i_re, i_im) = c_ext0 
                aero_lut%C_abs0_in_sph(:, i_re, i_im) = c_abs0 
             else if(sp_el.eq.1) then
                aero_lut%C_ext0_in_el(:, i_re, i_im) = c_ext0 
                aero_lut%C_abs0_in_el(:, i_re, i_im) = c_abs0 
             end if
          end do
       end do
       close(io)
    end do

  end subroutine read_cross

  !------------------------------------------------------------------------------
  !> Read scattering matrix from asciifile
  !------------------------------------------------------------------------------
  subroutine read_F(aerosol_dir, aero_lut)
    !*** Input
    character(len=*), intent(in) :: aerosol_dir
    !*** Input/output
    type(Mie_lut), intent(inout) :: aero_lut
    !*** Local variables
    integer :: sp_el, n_xgr01, n_Rem, n_Imm, i_re, i_im, nn, i_x, i_el, n_tet01, io, ierr
    character(stringlen) :: dir_dat_el, dir_dat_sp
    character(80) :: full_name
    real(double) :: Rem_min, Rem_max, Imm_min, Imm_max, Wavel
    real(double), dimension(n_tet0,n_xgr0) :: M_el0
    real(double), dimension(n_tet0) :: tet
    character(15) name11
    parameter(name11='Rkernel1.11.fix')
    !------------------------------------------------------------------

    dir_dat_sp = trim(aerosol_dir)//'Kernal_sph/'
    dir_dat_el = trim(aerosol_dir)//'Kernal_fix1/'
    call allocate_scattering(aero_lut, ierr)
    do sp_el = 0, 1
       do i_el = 1, n_el
          if (sp_el .eq. 0) then
             full_name = trim(dir_dat_sp)//trim(name11)
          end if
          if (sp_el .eq. 1) then
             full_name = trim(dir_dat_el)//trim(name11)
          end if
          open(newunit(io), file=full_name)
          read(io,*)
          read(io,*) n_xgr01
          if (abs(n_xgr01) .ne. n_xgr0) &
               call stopretrieval('READ_F: n_xgr01 .ne. n_xgr0')
          read(io,*) n_tet01
          if (abs(n_tet01) .ne. n_tet0) &
               call stopretrieval('READ_F: n_tet01 .ne. n_tet0') 
          read(io,*) tet(1:n_tet0)
          read(io,*) Rem_min,Rem_max    !!!!!!!! min and max of Re_m
          read(io,*) Imm_min,Imm_max    !!!!!!!! min and max of Im_m
          read(io,*) n_Rem, n_Imm          !!!!!!!! number of interval for optical const           
          if (n_Rem .lt. 0) n_Rem = -n_Rem
          if (n_Imm .lt. 0) n_Imm = -n_Imm            
          if(n_Rem.ne.n_Rem0.or.n_Imm.ne.n_Imm0)then
             call stopretrieval('READ_F: n_Rem.ne.n_Rem0.or.n_Imm.ne.n_Imm0')
          endif
          do i_re = 1, n_Rem       !!!!!! Re_m
             do i_im = 1, n_Imm      !!!!!! Im_m
                read(io,*) nn    !element,ratio
                read(io,*) Wavel, aero_lut%Rem_in(i_re), aero_lut%Imm_in(i_im) !wavelength, Re_m, Im_m
                if (aero_lut%Imm_in(i_im) .lt. 0) aero_lut%Imm_in(i_im) = -aero_lut%Imm_in(i_im)
                do  i_x = 1, n_xgr0
                   read(io,*) M_el0(1:n_tet0, i_x)
                   if(sp_el.eq.0)then
                      aero_lut%M_el0_in_sph(1:n_tet0, i_x, i_re, i_im, i_el) = &
                           M_el0(1:n_tet0, i_x)
                   else if(sp_el.eq.1)then
                      aero_lut%M_el0_in_el(1:n_tet0, i_x, i_re, i_im, i_el) = &
                           M_el0(1:n_tet0,i_x)
                   endif
                end do
             end do
          end do
          close(io)
       end do
    end do
  end subroutine read_F

  !------------------------------------------------------------------------
  !* Given the lower and upper limits of integration a and b, 
  !* and given the number of Gauss-Legendre points ngauss,
  !* this routine returns through array x the abscissas and through
  !* array w the weights of the Gauss-Legendre quadrature formula.
  !* Eps is the desired accuracy of the abscissas.
  !* This routine is documented further in:
  !*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)
  !*   page 125 ISBN 0-521-30811-9    
  !------------------------------------------------------------------------
  subroutine gauleg(ndim, ngauss, a, b, x, w)
    real(double),parameter:: eps = 1.d-14
    integer :: m, i, j, ndim, ngauss
    real(double) :: x(ndim), w(ndim), a, b, xm, xl, z, p1, p2, p3, pp, z1
    !----------------------------------------------------------------------
    m=(ngauss+1)/2
    xm=0.5D0*(a+b)
    xl=0.5D0*(b-a)

    do 12 i=1,m
       z= dcos(pi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
1      continue
       p1=1.D0
       p2=0.D0
       do j=1,ngauss
          p3= p2
          p2= p1
          p1=((dble(2*j)-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
       enddo
       pp=ngauss*(z*p1-p2)/(z*z-1.d0)
       z1= z
       z= z1-p1/pp
       if (dabs(z-z1).gt.eps) goto 1
       x(i)= xm-xl*z
       x(ngauss+1-i)= xm+xl*z
       w(i)=2.D0*xl/((1.D0-z*z)*pp*pp)
       w(ngauss+1-i)= w(i)
12     continue

       return
     end subroutine

     !------------------------------------------------------------------------
     !> bi-linear interpolation z11=z(x1,y1), z12=z(x1,y2), z21=z(x2,y1), z22=z(x2,y2) 
     !------------------------------------------------------------------------
     subroutine Bilin_interp(x1, x2, y1, y2, z, x0, y0, z0)
       implicit none
       real(double), intent(in) ::  x1, x2, y1, y2, z(1:2,1:2), y0, x0
       real(double), intent(out) :: z0
       !*** Local variables
       real(double) :: dx1, dx2, dx21, dy1, dy2, dy21
       !------------------------------------------------------
       dx1=x0-x1
       dx2=x2-x0
       dx21=x2-x1

       dy1=y0-y1
       dy2=y2-y0
       dy21=y2-y1

       z0=(z(1,1)*dx2*dy2+z(2,1)*dx1*dy2 + &
            z(1,2)*dx2*dy1+z(2,2)*dx1*dy1)/dx21/dy21
     end subroutine Bilin_interp

     !------------------------------------------------------------------------
     !> Find array indices of real refractive index and imaginary refractive index
     !------------------------------------------------------------------------
     subroutine RemIm_index(aero_lut, Re_m, Im_m, i_re12, i_im12)   
       !*** Input
       type(Mie_lut), intent(in) :: aero_lut
       real(double), intent (in):: Re_m, Im_m
       !*** Output
       integer, intent (out):: i_re12(2), i_im12(2)
       !*** Local variables
       integer :: idum
       !------------------------------------------------------------   
       idum = minval(minloc(DABS(Re_m-aero_lut%Rem_in),Re_m.ge.aero_lut%Rem_in))
       i_re12(1) = min(idum,n_rem0-1)
       if(i_re12(1).eq.0)i_re12(1)=1
       i_re12(2) = i_re12(1)+1

       idum = minval(minloc(DABS(Im_m-aero_lut%Imm_in),Im_m.ge.aero_lut%Imm_in))
       i_im12(1) = min(idum, n_imm0-1)
       if(i_im12(1).eq.0) i_im12(1)=1
       i_im12(2) = i_im12(1)+1

     end subroutine RemIm_index

     !------------------------------------------------------------------------
     !> Get extinction and absorption coefficients by interpolating in LUT
     !------------------------------------------------------------------------
     subroutine C_read_int( &
          aero_lut, &
          sp_el, &         ! sp_el is = 0 for spheres or 1 for spheroids
          Re_m, &          ! Re_m - real part of refractive index
          Im_m, &          ! Im_m - imaginary part of refractive index
          i_re12, &        ! i_re12(1:2) - indices for real refractive index: Rem(i_re12(1))<Re_m<Rem(i_re12(2))
          i_im12, &        ! i_im12(1:2) - indices for imag refractive index: Imm(i_im12(1))<Im_m<Imm(i_im12(2))
          Cext_int, &      ! extinction cross sections obtained for Re_m, Im_m
          Cabs_int, &      ! absorption cross sections obtained for Re_m, Im_m
          der_Cext_rem, &  ! Derivatives
          der_Cabs_rem, &
          der_Cext_imm, &
          der_Cabs_imm)
       !*** Input
       type(Mie_lut), intent(in) :: aero_lut
       real(double), intent (in):: Re_m, Im_m
       integer, intent (in):: i_re12(2), i_im12(2), sp_el 
       !*** Output
       real(double), intent (out):: Cext_int(n_xgr0), Cabs_int(n_xgr0),&
            der_Cext_rem(n_xgr0), der_Cabs_rem(n_xgr0), &
            der_Cext_imm(n_xgr0), der_Cabs_imm(n_xgr0)
       !*** Local variables
       real(double) :: C_ext(n_xgr0, 2, 2), C_abs(n_xgr0, 2, 2)
       real(double) Cext_int1(n_xgr0), Cabs_int1(n_xgr0), Re_m1, Im_m1
       integer :: i_re, i_im, i_x
       !---------------------------------------------------------------------------  

       if(sp_el.eq.0)then      
          do i_im = 1,2  
             do i_re = 1, 2            
                C_ext(:,i_re,i_im) = aero_lut%C_ext0_in_sph(:, i_re12(i_re), i_im12(i_im))
                C_abs(:,i_re,i_im) = aero_lut%C_abs0_in_sph(:, i_re12(i_re), i_im12(i_im))
             enddo
          enddo
       else if(sp_el.eq.1)then
          do i_im = 1, 2    
             do i_re = 1, 2       
                C_ext(:,i_re,i_im) = aero_lut%C_ext0_in_el(:, i_re12(i_re), i_im12(i_im))             
                C_abs(:,i_re,i_im) = aero_lut%C_abs0_in_el(:, i_re12(i_re), i_im12(i_im))
             enddo
          enddo
       endif

       !*** Interpolation
       do i_x = 1, n_xgr0 
          call Bilin_interp( &
               aero_lut%Rem_in(i_re12(1)), &
               aero_lut%Rem_in(i_re12(2)), &
               aero_lut%Imm_in(i_im12(1)), &
               aero_lut%Imm_in(i_im12(2)),&
               C_ext(i_x,:,:), &
               Re_m, &
               Im_m, &
               Cext_int(i_x))
          call Bilin_interp( &
               aero_lut%Rem_in(i_re12(1)), &
               aero_lut%Rem_in(i_re12(2)), &
               aero_lut%Imm_in(i_im12(1)), &
               aero_lut%Imm_in(i_im12(2)), &
               C_abs(i_x,:,:), &
               Re_m, &
               Im_m, &
               Cabs_int(i_x))
          !*** For derivative wrt Rem
          if (dabs(Re_m-aero_lut%Rem_in(i_re12(1))) .gt. dabs(Re_m-aero_lut%Rem_in(i_re12(2)))) then
             Re_m1 = aero_lut%Rem_in(i_re12(1))
          else
             Re_m1 = aero_lut%Rem_in(i_re12(2))
          end if
          call Bilin_interp( &
               aero_lut%Rem_in(i_re12(1)), &
               aero_lut%Rem_in(i_re12(2)), &
               aero_lut%Imm_in(i_im12(1)), &
               aero_lut%Imm_in(i_im12(2)),&
               C_ext(i_x,:,:), &
               Re_m1, &
               Im_m, &
               Cext_int1(i_x))
          call Bilin_interp( &
               aero_lut%Rem_in(i_re12(1)), &
               aero_lut%Rem_in(i_re12(2)), &
               aero_lut%Imm_in(i_im12(1)), &
               aero_lut%Imm_in(i_im12(2)),&
               C_abs(i_x,:,:), &
               Re_m1, &
               Im_m, &
               Cabs_int1(i_x))
          der_Cext_rem(i_x)=(Cext_int(i_x)-Cext_int1(i_x))/(Re_m-Re_m1)
          der_Cabs_rem(i_x)=(Cabs_int(i_x)-Cabs_int1(i_x))/(Re_m-Re_m1)

          !*** For derivative wrt Imm
          if (dabs(Im_m-aero_lut%Imm_in(i_im12(1))) .gt. dabs(Im_m-aero_lut%Imm_in(i_im12(2)))) then
             Im_m1 = aero_lut%Imm_in(i_im12(1))
          else
             Im_m1 = aero_lut%Imm_in(i_im12(2))
          end if
          call Bilin_interp( &
               aero_lut%Rem_in(i_re12(1)), &
               aero_lut%Rem_in(i_re12(2)), &
               aero_lut%Imm_in(i_im12(1)), &
               aero_lut%Imm_in(i_im12(2)),&
               C_ext(i_x,:,:), &
               Re_m, &
               Im_m1, &
               Cext_int1(i_x))
          call Bilin_interp( &
               aero_lut%Rem_in(i_re12(1)), &
               aero_lut%Rem_in(i_re12(2)), &
               aero_lut%Imm_in(i_im12(1)), &
               aero_lut%Imm_in(i_im12(2)),&
               C_abs(i_x,:,:), &
               Re_m, &
               Im_m1, &
               Cabs_int1(i_x))
          der_Cext_imm(i_x) = (Cext_int(i_x)-Cext_int1(i_x))/(Im_m-Im_m1)
          der_Cabs_imm(i_x) = (Cabs_int(i_x)-Cabs_int1(i_x))/(Im_m-Im_m1)

       end do !i_x

     end subroutine C_read_int

     !------------------------------------------------------------------------
     !> Get scattering matrix F by interpolating in LUT
     !------------------------------------------------------------------------
     subroutine F_read_int( &
          aero_lut, &
          sp_el, Re_m, Im_m, &
          i_re12, i_im12, &
          Mel_int, der_Mel_rem, der_Mel_imm)
       !*** Input
       real(double), intent (in):: Re_m, Im_m
       integer, intent (in):: i_re12(2), i_im12(2), sp_el
       type(Mie_lut), intent(in) :: aero_lut
       !*** Output
       real(double), intent (out):: Mel_int(n_el, n_tet0, n_xgr0), der_Mel_rem(n_el, n_tet0, n_xgr0),&
            der_Mel_imm(n_el, n_tet0, n_xgr0)
       !*** Local variables
       real(double) :: Mel_int11(n_el, n_tet0, n_xgr0), Re_m1, Im_m1
       integer :: i_re, i_im, i_x, i_tet, i_el
       real(double) ::  M_el(n_tet0, n_xgr0, 2, 2)
       !----------------------------------------------------------------------------------------------------

       do i_el = 1, n_el
          if(sp_el.eq.0) then          
             do i_im = 1, 2  
                do i_re = 1, 2
                   do i_x = 1,n_xgr0         
                      M_el(1:n_tet0, i_x, i_re, i_im) = &					
                           aero_lut%M_el0_sph(1:n_tet0, i_x, i_re12(i_re), i_im12(i_im), i_el)	
                   enddo
                enddo
             enddo
          else if(sp_el.eq.1) then
             do i_im = 1,2
                do i_re = 1, 2
                   do i_x = 1, n_xgr0
                      M_el(1:n_tet0, i_x, i_re, i_im) = &					
                           aero_lut%M_el0_el(1:n_tet0, i_x, i_re12(i_re), i_im12(i_im), i_el)   	
                   enddo
                enddo
             enddo
          endif

          !*** Interpolation
          do i_x = 1, n_xgr0 
             do i_tet = 1, n_tet0
                call Bilin_interp( &
                     aero_lut%Rem_in(i_re12(1)), &
                     aero_lut%Rem_in(i_re12(2)), &
                     aero_lut%Imm_in(i_im12(1)), &
                     aero_lut%Imm_in(i_im12(2)),&
                     M_el(i_tet, i_x,:,:), &
                     Re_m, &
                     Im_m, &
                     Mel_int(i_el, i_tet, i_x))

                !***For derivative wrt Rem
                if (dabs(Re_m-aero_lut%Rem_in(i_re12(1))) .gt. dabs(Re_m-aero_lut%Rem_in(i_re12(2)))) then
                   Re_m1 = aero_lut%Rem_in(i_re12(1))
                else
                   Re_m1 = aero_lut%Rem_in(i_re12(2))
                end if
                call Bilin_interp( &
                     aero_lut%Rem_in(i_re12(1)), &
                     aero_lut%Rem_in(i_re12(2)), &
                     aero_lut%Imm_in(i_im12(1)), &
                     aero_lut%Imm_in(i_im12(2)),&
                     M_el(i_tet, i_x, :, :), & 
                     Re_m1, &
                     Im_m, &
                     Mel_int11(i_el, i_tet, i_x))
                der_Mel_rem(i_el, i_tet, i_x) = (Mel_int(i_el,i_tet,i_x) - Mel_int11(i_el,i_tet,i_x))/(Re_m-Re_m1)

                !*** For derivative wrt Imm
                if (dabs(Im_m-aero_lut%Imm_in(i_im12(1))) .gt. dabs(Im_m-aero_lut%Imm_in(i_im12(2)))) then
                   Im_m1 = aero_lut%Imm_in(i_im12(1))
                else
                   Im_m1 = aero_lut%Imm_in(i_im12(2))
                end if
                call Bilin_interp( &
                     aero_lut%Rem_in(i_re12(1)), &
                     aero_lut%Rem_in(i_re12(2)), &
                     aero_lut%Imm_in(i_im12(1)), &
                     aero_lut%Imm_in(i_im12(2)),&
                     M_el(i_tet, i_x, :, :), &
                     Re_m, &
                     Im_m1, &
                     Mel_int11(i_el, i_tet, i_x))
                der_Mel_imm(i_el, i_tet, i_x) = (Mel_int(i_el,i_tet,i_x) - Mel_int11(i_el,i_tet,i_x))/(Im_m-Im_m1)

             end do !i_tet
          end do !i_x
       end do !i_el

     end subroutine F_read_int

     !------------------------------------------------------------------------
     !> Calculate volume size distribution and its derivatives
     !------------------------------------------------------------------------
     subroutine VolDistrFunc( &
          id_size, &           ! which size distributions 
          wavel, &             ! wavelength
          r_eff, &             ! effective radius
          v_eff, &             ! effective variance
          x0, &                ! size grid
          V_distr,&            ! volume size distribution function
          Der_Vdist_reff, &    ! derivatives
          Der_Vdist_veff)

       !*** Input
       integer, intent(in)::  id_size
       real(double), intent(in):: r_eff, v_eff, wavel, x0(n_xgr0)
       !*** Output
       real(double), intent (out):: V_distr(n_xgr0), Der_Vdist_reff(n_xgr0), &
            Der_Vdist_veff(n_xgr0)
       !*** Local variables
       real(double) :: log_sig, r_grid(n_xgr0), v_grid(n_xgr0),rup, rlow, s, der_sum_reff,der_sum_veff
       integer :: i
       real(double) :: root2p, rg, flogrg, flogsi, C, fac, rdum, exp_term
       real(double) :: der_flogrg_R, der_flogrg_V, der_flogsi_V, der_C_V, &
            der_fac_V, der_RDUM, der_exp_V
       real(double), dimension(n_xgr0) :: dr,dlnr
       !----------------------------------------------------------------------------------------
       do i = 1, n_xgr0
          r_grid(i) = x0(i)/2.d0/pi*wavel
          v_grid(i)= 4./3. *pi*r_grid(i)**3
       end do
       do i = 1, n_xgr0-1
          dr(i) = r_grid(i+1) - r_grid(i)
          dlnr(i) = log(r_grid(i+1)) - log(r_grid(i))
       end do
       dr(n_xgr0) = dr(n_xgr0-1)
       dlnr(n_xgr0) = dlnr(n_xgr0-1)

       !*** powe law size distribution
       if(id_size.eq.0) then   
          V_distr = 0.D0
          rup = v_eff
          rlow = 0.1D0
          do i = 1, n_xgr0
             if(r_grid(i).le.rlow)then
                V_distr(i) = 1.D0
                Der_Vdist_reff(i) = 0.D0
                Der_Vdist_veff(i) = 0.D0
             endif
             if(r_grid(i).gt.rup)then
                V_distr(i) = 0.D0
                Der_Vdist_reff(i) = 0.D0
                Der_Vdist_veff(i) = 0.D0
             endif
             if(r_grid(i).gt.rlow.and.r_grid(i).le.rup)then
                V_distr(i) = (r_grid(i)/rlow)**(-r_eff)
                Der_Vdist_reff(i) = -DLOG(r_grid(i)/rlow)*V_distr(i)
                Der_Vdist_veff(i) = 0.D0
             endif
          enddo

          do i = 1, n_xgr0
             v_distr(i) = v_distr(i)*r_grid(i)*v_grid(i)
             Der_Vdist_reff(i) = Der_Vdist_reff(i)*r_grid(i)*v_grid(i)
          enddo

          s = sum(dlnr*V_distr/v_grid)
          der_sum_reff = sum(dlnr*der_Vdist_reff/v_grid)   
          der_sum_reff = -s**(-2)*der_sum_reff



          do i = 1, n_xgr0
             der_Vdist_reff(i) = der_Vdist_reff(i)/s + &
                  der_sum_reff*V_distr(i)
          enddo
          V_distr = V_distr / s

          !*** lognormal size distribution
       elseif(id_size.eq.1)then
          log_sig = dsqrt(dlog(v_eff+1.d0))
          root2p = dsqrt(pi + pi)
          rg     = r_eff/(1.D0+v_eff)**2.5D0
          flogrg = dlog(rg)
          der_flogrg_R = (1./rg)*(1.D0+v_eff)**(-2.5D0)
          der_flogrg_V = -2.5*(1./(rg))*R_EFF*(1+V_EFF)**(-3.5)
          flogsi = dsqrt(dlog(1.D0+v_eff))
          der_flogsi_V = 0.5*(DLOG(1+V_EFF))**(-0.5)*(1./(1+V_EFF))
          C      = 1.D0/(root2p*flogsi)
          der_C_V = (-1.D0/(root2p))*(flogsi**(-2))*der_flogsi_V
          fac    = -0.5D0/(flogsi*flogsi)
          der_fac_V = (flogsi**(-3))*der_flogsi_V

          do i = 1, n_xgr0
             rdum = ( fac*(dlog(r_grid(i))-flogrg)**2)
             exp_term = dexp(RDUM)/r_grid(i)
             V_distr(i) =  C * exp_term *r_grid(i)*v_grid(i)
             der_RDUM = der_fac_V*(dlog(r_grid(i))-flogrg)**2 +fac*2*(dlog(r_grid(i))-flogrg)*(-der_flogrg_V)
             der_exp_V = exp_term*der_RDUM
             der_RDUM = fac*2*(dlog(r_grid(i))-flogrg)*(-der_flogrg_R)
             Der_Vdist_reff(i) = C*exp_term*der_RDUM*r_grid(i)*v_grid(i)
             Der_Vdist_veff(i) = (C*der_EXP_V + exp_term*der_C_V) *r_grid(i)*v_grid(i)
          end do

          s = sum(dlnr*V_distr/v_grid)
          der_sum_reff = sum(dlnr*der_Vdist_reff/v_grid)   
          der_sum_reff = -s**(-2)*der_sum_reff
          der_sum_veff = sum(dlnr*der_Vdist_veff/v_grid)   
          der_sum_veff = -s**(-2)*der_sum_veff

          do i = 1, n_xgr0
             der_Vdist_reff(i) = der_Vdist_reff(i)/s + &
                  der_sum_reff*V_distr(i)

             der_Vdist_veff(i) = der_Vdist_veff(i)/s + &
                  der_sum_veff*V_distr(i)
          enddo
          V_distr = V_distr / s



          !*** 2par Gamma size distribution -- 
       elseif(id_size.eq.2)then
          !          log_sig = dsqrt(dlog(v_eff+1.d0))
          !          root2p = dsqrt(pi + pi)
          !          rg     = r_eff/(1.D0+v_eff)**2.5D0
          !          flogrg = dlog(rg)
          !          der_flogrg_R = (1./rg)*(1.D0+v_eff)**(-2.5D0)
          !          der_flogrg_V = -2.5*(1./(rg))*R_EFF*(1+V_EFF)**(-3.5)
          !          flogsi = dsqrt(dlog(1.D0+v_eff))
          !          der_flogsi_V = 0.5*(DLOG(1+V_EFF))**(-0.5)*(1./(1+V_EFF))
          !          C      = 1.D0/(root2p*flogsi)
          !          der_C_V = (-1.D0/(root2p))*(flogsi**(-2))*der_flogsi_V
          !          fac    = -0.5D0/(flogsi*flogsi)
          !          der_fac_V = (flogsi**(-3))*der_flogsi_V

          do i = 1, n_xgr0
             !             rdum = ( fac*(dlog(r_grid(i))-flogrg)**2)
             !             exp_term = dexp(RDUM)/r_grid(i)
             ! !             V_distr(i) =  C * exp_term *r_grid(i)*v_grid(i)
             !             der_RDUM = der_fac_V*(dlog(r_grid(i))-flogrg)**2 +fac*2*(dlog(r_grid(i))-flogrg)*(-der_flogrg_V)
             !             der_exp_V = exp_term*der_RDUM
             !             der_RDUM = fac*2*(dlog(r_grid(i))-flogrg)*(-der_flogrg_R)

             ! For allocation reasons, I re-use some declared variables here but there names DO NOT MAKE SENSE in anyway!
             rdum        = (-3./v_eff - (1.-3.*v_eff)/v_eff/v_eff )!prefac
             der_exp_V   = r_grid(i)**((1. - 3.*v_eff )/v_eff)*dexp(-r_grid(i)/(r_eff*v_eff) )!orig
             exp_term    = dlog(r_grid(i))!logfac
             der_RDUM    = r_grid(i)**((1. - 3.*v_eff )/v_eff)*dexp(-r_grid(i)/(r_eff*v_eff) )*(r_grid(i)/(r_eff*v_eff*v_eff)) !add_term
             !End of section without sense

             V_distr(i) = r_grid(i)**((1. - 3.*v_eff )/v_eff)*dexp(-r_grid(i)/(r_eff*v_eff) )*r_grid(i)*v_grid(i)
             Der_Vdist_reff(i) = r_grid(i)**((1. - 3.*v_eff )/v_eff)*dexp(-r_grid(i)/(r_eff*v_eff) )*(r_grid(i)/(r_eff*r_eff*v_eff))*r_grid(i)*v_grid(i)
             Der_Vdist_veff(i) = ( rdum * der_exp_V * exp_term  +  der_RDUM )*r_grid(i)*v_grid(i)
          end do

          s = sum(dlnr*V_distr/v_grid)
          der_sum_reff = sum(dlnr*der_Vdist_reff/v_grid)   
          der_sum_reff = -s**(-2)*der_sum_reff
          der_sum_veff = sum(dlnr*der_Vdist_veff/v_grid)   
          der_sum_veff = -s**(-2)*der_sum_veff


          do i = 1, n_xgr0
             der_Vdist_reff(i) = der_Vdist_reff(i)/s + &
                  der_sum_reff*V_distr(i)

             der_Vdist_veff(i) = der_Vdist_veff(i)/s + &
                  der_sum_veff*V_distr(i)
          enddo
          V_distr = V_distr / s


       endif
     end subroutine VolDistrFunc

     !***********************
     subroutine Spl_interp_v6(n_el1,n_tet0,n_tet,tet0,tet,Mel_av,M_el,ierr)
       integer, intent(in):: n_el1,n_tet0,n_tet
       real(8), intent(in):: tet0(1:n_tet0),tet(1:n_tet),Mel_av(1:n_tet0,1:n_el1)
       real(8), intent(out):: M_el(1:n_tet,1:n_el1)

       real(8) M_tet(1:n_tet0),KS1(1:n_tet0+4),CS1(1:n_tet0+4),Spl_int
       integer i_el,key_spln,i
       integer, intent(out) :: ierr

       do i_el=1, n_el1
          !   Write(*,*) 'element', i_el

          if ((i_el .eq. 1) .or. (i_el .eq. 3)) then
             !    M_tet=dlog(Mel_av(i_el,1:n_tet0))
             M_tet=(Mel_av(1:n_tet0,i_el))
          else 
             if ((i_el .eq. 4) .or. (i_el .eq. 6)) then
                M_tet=Mel_av(1:n_tet0,i_el)/Mel_av(1:n_tet0,1)
             else
                M_tet=Mel_av(1:n_tet0,i_el)
             end if
          end if

          key_spln=0
          do i=1, n_tet
             call intrpl_spline(n_tet0,tet0,M_tet,tet(i),Spl_int,key_spln,&
                  KS1(1:n_tet0+4),CS1(1:n_tet0+4), ierr)

             if (ierr .ne. 0) then
                ierr = ierr_intrpl
                return
             endif

             if ((i_el .eq. 1) .or. (i_el .eq. 3)) then
                !     M_el(i,i_el)=dexp(Spl_int)
                M_el(i,i_el)=(Spl_int)
             else
                if ((i_el .eq. 4) .or. (i_el .eq. 6)) then
                   M_el(i,i_el)=Spl_int*M_el(i,1)
                else
                   M_el(i,i_el)=Spl_int
                end if
             end if
          end do !i     

       end do !i_el

     end subroutine Spl_interp_v6

     !------------------------------------------------------------------------
     !>  Spline interpolation
     !------------------------------------------------------------------------
     subroutine Spl_interp(n_el1, n_tet0, n_tet, tet0, tet, Mel_av, M_el, ierr)
       !*** input
       integer, intent(in):: n_el1, n_tet0, n_tet
       real(double), intent(in):: tet0(n_tet0), tet(n_tet), Mel_av(n_tet0,n_el1)
       !*** output
       real(double), intent(out):: M_el(n_tet,n_el1)
       integer, intent(out) :: ierr
       !*** local
       real(double) ::  M_tet(n_tet0), KS1(n_tet0+4), CS1(n_tet0+4), Spl_int
       integer :: i_el,key_spln,i
       !------------------------------------

       do i_el = 1, n_el1
          M_tet = dlog(Mel_av(1:n_tet0,i_el))
          key_spln = 0
          do i = 1, n_tet
             call intrpl_spline(n_tet0, tet0, M_tet, tet(i), Spl_int, key_spln,&
                  KS1(1:n_tet0+4), CS1(1:n_tet0+4), ierr)
             if (ierr .ne. 0) then
                ierr = ierr_intrpl
                return
             endif
             M_el(i,i_el) = dexp(Spl_int)
          end do !i      
       end do !i_el

     end subroutine Spl_interp

!!$     !------------------------------------------------------------------------
!!$     !> Spline interpolation for derivatives
!!$     !------------------------------------------------------------------------
!!$     subroutine Spl_interp_der(n_el1, n_tet0, n_tet, tet0, tet, Mel_av, M_el, ierr)
!!$       !*** input   
!!$       integer, intent(in):: n_el1, n_tet0, n_tet
!!$       real(double), intent(in):: tet0(n_tet0), tet(n_tet), Mel_av(n_el1, n_tet0)
!!$       !*** output
!!$       real(double), intent(out):: M_el(n_el1, n_tet) 
!!$       integer, intent(out) :: ierr
!!$       !*** local
!!$       real(double) :: M_tet(n_tet0), KS1(n_tet0+4), CS1(n_tet0+4), Spl_int
!!$       integer i_el, key_spln, i
!!$       !--------------------------------------------
!!$       do i_el = 1, n_el1
!!$          M_tet = Mel_av(i_el, 1:n_tet0)
!!$          key_spln = 0
!!$          do i = 1, n_tet
!!$             call intrpl_spline(n_tet0,tet0, M_tet,tet(i), Spl_int, key_spln, &
!!$                  KS1(1:n_tet0+4), CS1(1:n_tet0+4), ierr)
!!$             if (ierr .ne. 0) then
!!$                ierr = ierr_intrpl
!!$                return
!!$             endif
!!$             M_el(i_el,i) = Spl_int
!!$          enddo !i     
!!$       enddo !i_el
!!$
!!$     end subroutine Spl_interp_der


     !------------------------------------------------------------------------
     !>  Calculate the expansion coefficients of the scattering matrix in 
     !!  generalized spherical functions by numerical integration over the
     !!  scattering angle
     !------------------------------------------------------------------------
     subroutine devel(&
          u,&
          w,&
          F,&
          coefs,&
          der_F,&
          dcoefs)
       real(double), dimension(:), intent(in) :: u, w  
       real(double), dimension(:,:), intent(in) :: F
       real(double), dimension(:,:,:), intent(in) :: der_F  
       !*** Output
       real(double), dimension(:,:,0:), intent(out) :: coefs
       real(double), dimension(:,:,0:,:), intent(out) :: dcoefs
       !*** local variables
       real(double), dimension(size(F(:,1)),size(F(1,:))) :: Fdum
       real(double), dimension(size(der_F(:,1,1)),size(der_F(1,:,1)),size(der_F(1, 1, :))) :: der_Fdum
       real(double), dimension(size(F(1,:)), 2) :: P00, P02, P22, P2m2
       real(double), dimension(size(der_F(1,1,:))) :: der_alfap, der_alfam
       integer :: nangle, nx, ncoef, nrf, i, itmp, j, l, lnew, lold, k
       real(double) :: fac1, fac2, fac3, sql41, sql4, twol1, tmp1, tmp2, denom, alfap, alfam, fl
       real(double) :: qroot6 
       !------------------------------------------------------------------

       qroot6 = -0.25D0*dsqrt(6.D0)
       nangle = size(u) 
       nx = size(der_F(1,1,:))
       ncoef = size(coefs(1,1,:))-1 
       nrf = 2

       !***  Initialization
       coefs = 0.d0
       dcoefs = 0.d0

       !*** Multiply the scattering matrix F with the weights w for all angles  
       !*** We do this here because otherwise it should be done for each l
       do k = 1, n_el
          do i = 1, nangle
             Fdum(k,i) = w(i)*F(k,i)
             do j = 1, nx
                der_Fdum(k,i,j) = w(i)*der_F(k,i,j)
             enddo
          enddo
       enddo

       !*** Start loop over the coefficient index l  
       !*** first update generalized spherical functions, then calculate coefs. 
       !*** lold and lnew are pointer-like indices used in recurrence 
       lnew = 1
       lold = 2
       do l = 0, ncoef
          if (l.eq.0) then
             !*** Adding paper Eq. (77) with m=0 
             do i = 1,nangle
                P00(i,lold) = 1.D0
                P00(i,lnew) = 0.D0
                P02(i,lold) = 0.D0
                P22(i,lold) = 0.D0
                P2m2(i,lold)= 0.D0
                P02(i,lnew) = 0.D0
                P22(i,lnew) = 0.D0
                P2m2(i,lnew)= 0.D0
             enddo
          else
             fac1 = (2.D0*l-1.d0)/dble(l)
             fac2 = dble(l-1.d0)/dble(l)
             !** Adding paper Eq. (81) with m=0   
             do i=1,nangle
                P00(i,lold) = fac1*u(i)*P00(i,lnew) - fac2*P00(i,lold)
             enddo
          endif
          if (l.eq.2) then
             !*** Adding paper Eqs. (78) and (80)  with m=2   
             !*** sql4 contains the factor dsqrt(l*l-4) needed in         
             !*** the recurrence Eqs. (81) and (82)                        
             do i=1,nangle
                P02(i,lold) = qroot6*(1.D0-u(i)*u(i))
                P22(i,lold) = 0.25D0*(1.D0+u(i))*(1.D0+u(i))
                P2m2(i,lold)= 0.25D0*(1.D0-u(i))*(1.D0-u(i))
                P02(i,lnew) = 0.D0
                P22(i,lnew) = 0.D0
                P2m2(i,lnew)= 0.D0
             enddo
             sql41 = 0.D0
          else if (l.gt.2) then
             !*** Adding paper Eq. (82) with m=0 and m=2 
             sql4 = sql41
             sql41= dsqrt(dble(l*l)-4.d0)
             twol1= 2.D0*dble(l)-1.d0
             tmp1 = twol1/sql41
             tmp2 = sql4/sql41
             denom= (dble(l)-1.d0)*(dble(l*l)-4.d0)
             fac1 = twol1*(dble(l)-1.d0)*dble(l)/denom
             fac2 = 4.D0*twol1/denom
             fac3 = dble(l)*((dble(l)-1.d0)*(dble(l)-1.d0)-4.d0)/denom
             do i=1,nangle
                P02(i,lold) = tmp1*u(i)*P02(i,lnew) - tmp2*P02(i,lold)
                P22(i,lold) = (fac1*u(i)-fac2)*P22(i,lnew)&
                     - fac3*P22(i,lold)
                P2m2(i,lold)= (fac1*u(i)+fac2)*P2m2(i,lnew)&
                     - fac3*P2m2(i,lold)
             enddo
          endif

          !*** Switch indices so that lnew indicates the function with 
          !*** the present index value l, this mechanism prevents swapping 
          !*** of entire arrays. 
          itmp = lnew
          lnew = lold
          lold = itmp

          !*** Now calculate the coefficients by integration over angle 
          !*** See de Haan et al. (1987) Eqs. (68)-(73).                   
          !*** Remember for Mie scattering : F11 = F22 and F33 = F44        
          alfap = 0.D0
          alfam = 0.D0
          do j = 1, nx
             der_alfap(j) = 0.
             der_alfam(j) = 0.
          enddo
          do i = 1, nangle
             coefs(1,1,l) = coefs(1,1,l) + P00(i,lnew)*Fdum(1,i)
             if (n_el==6) then
                alfap = alfap + P22(i,lnew)*(Fdum(3,i)+Fdum(4,i))
                alfam = alfam + P2m2(i,lnew)*(Fdum(3,i)-Fdum(4,i))
                coefs(4,4,l) = coefs(4,4,l) + P00(i,lnew)*Fdum(6,i)
                coefs(1,2,l) = coefs(1,2,l) + P02(i,lnew)*Fdum(2,i)
                coefs(3,4,l) = coefs(3,4,l) + P02(i,lnew)*Fdum(5,i)
             endif

             do j = 1, nx
                dcoefs(1,1,l,j) = dcoefs(1,1,l,j) + & 
                     P00(i,lnew)*der_Fdum(1,i,j)
                if (n_el==6) then
                   dcoefs(4,4,l,j) = dcoefs(4,4,l,j) +& 
                        P00(i,lnew)*der_Fdum(6,i,j)
                   dcoefs(1,2,l,j) = dcoefs(1,2,l,j) +& 
                        P02(i,lnew)*der_Fdum(2,i,j)
                   dcoefs(3,4,l,j) = dcoefs(3,4,l,j) + &
                        P02(i,lnew)*der_Fdum(5,i,j)
                   der_alfap(j) = der_alfap(j) +&
                        P22(i,lnew)*(der_Fdum(3,i,j)+der_Fdum(4,i,j))
                   der_alfam(j) = der_alfam(j) + &
                        P2m2(i,lnew)*(der_Fdum(3,i,j)-der_Fdum(4,i,j))
                endif
             enddo
          enddo

          !*** Multiply with trivial factors like 0.5D0*(2*l+1)             
          fl = dble(l)+0.5D0
          coefs(1,1,l) =  fl*coefs(1,1,l)
          if (n_el==6) then
             coefs(2,2,l) =  fl*0.5D0*(alfap+alfam)
             coefs(3,3,l) =  fl*0.5D0*(alfap-alfam)
             coefs(4,4,l) =  fl*coefs(4,4,l)
             coefs(1,2,l) =  fl*coefs(1,2,l)
             coefs(3,4,l) =  fl*coefs(3,4,l)
             coefs(2,1,l) =     coefs(1,2,l)
             coefs(4,3,l) =    -coefs(3,4,l)
          endif
          do j = 1, nx
             dcoefs(1,1,l,j) =  fl*dcoefs(1,1,l,j)
             if (n_el==6) then
                dcoefs(2,2,l,j) =  fl*0.5D0*(der_alfap(j)+der_alfam(j))
                dcoefs(3,3,l,j) =  fl*0.5D0*(der_alfap(j)-der_alfam(j))
                dcoefs(4,4,l,j) =  fl*dcoefs(4,4,l,j)
                dcoefs(1,2,l,j) =  fl*dcoefs(1,2,l,j)
                dcoefs(3,4,l,j) =  fl*dcoefs(3,4,l,j)
                dcoefs(2,1,l,j) =     dcoefs(1,2,l,j)
                dcoefs(4,3,l,j) =    -dcoefs(3,4,l,j)
             endif
          enddo
       enddo   !End of loop over index l 

       return
     end subroutine devel

     !------------------------------------------------------------------------
     !> Calculate scattering characteristics for multiple scattering
     !------------------------------------------------------------------------
     subroutine SpEl_exp( &
          aero_lut, &
          Re_m, &
          Im_m, &
          delta,&
          i_re12, &
          i_im12, &
          V_distr, &
          coefs, &
          dcoefs,&
          csc_av, &
          cabs_av, &
          dcsca_full, &
          dcabs_full, &
          ierr)
       !*** input   
       type(Mie_lut), intent(in) :: aero_lut
       integer, intent(in):: i_re12(2),i_im12(2)
       real(double), intent (in):: Re_m, Im_m
       real(double), intent(in):: delta, V_distr(n_xgr0) 
       !*** output
       real(double), dimension(:,:,0:), intent(out) :: coefs
       real(double), dimension(:,:,0:,:), intent(out) :: dcoefs
       real(double), intent(out) :: Cabs_av, Csc_av
       real(double), dimension(:), intent(out) :: dcsca_full, dcabs_full
       integer, intent(out) :: ierr
       !*** local
       real(double) :: Der_Cext(4), Der_Cabs(4), Der_Cex_df(n_xgr0),Der_Cab_df(n_xgr0)
       real(double)  Cext_av
       integer :: sp_el
       real(double), dimension(n_xgr0):: C_ext, C_abs, C_ext0,C_abs0, C_ext1,C_abs1,&
            der_Cex_rem, der_Cab_rem, der_Cex_imm, der_Cab_imm,&
            der_Cex_rem1, der_Cab_rem1, der_Cex_imm1, der_Cab_imm1
       real(double), dimension(n_el, n_tet0, n_xgr0) :: Mel_int, Mel_int1, Der_Mel_rem0, Der_Mel_imm0, &
            Der_Mel_rem1, Der_Mel_imm1, Mel_int0
       real(double), dimension(n_el, n_tet0) :: Mel_sp, Mel_el, Mel_delt, Mel_rem, Mel_imm, Mel_av
       integer :: i, i_el,  i_tet, n_tet
       real(double) :: der_CexAv_delt, der_CabAv_delt, der_CexAv_rem, der_CabAv_rem, &
            der_CexAv_imm, der_CabAv_imm, &
            Cext_av0, Cext_av1, Cabs_av0, Cabs_av1
       real(double) :: x_max
       real(double), allocatable:: Gw_t(:), tet(:), mu_t(:)
       real(double), allocatable:: M_el(:,:), Der_delta(:,:), Der_rem(:,:),&
            Der_imm(:,:), Der_DestF(:,:,:)
       real(double), dimension(:,:,:), allocatable :: der_F
       character(stringlen) :: message
       !----------------------------------------------

       !*** Initialize error identifier
       ierr = 0

       !*** Get cros-sections
       !*** spheres (interpolating in LUT)
       sp_el = 0
       call C_read_int( &
            aero_lut, &
            sp_el, &
            Re_m, Im_m, i_re12, i_im12, &
            C_ext0, C_abs0, &
            der_Cex_rem, der_Cab_rem, der_Cex_imm, der_Cab_imm)

       !*** elipsoids (interpolating in LUT)
       sp_el = 1
       call C_read_int( &
            aero_lut, &
            sp_el, &
            Re_m, Im_m, i_re12, i_im12, &
            C_ext1, C_abs1, &
            der_Cex_rem1, der_Cab_rem1, der_Cex_imm1, der_Cab_imm1)

       !*** extinction and absoption coefficientfs for mixture of spheroids and ellipsoids
       C_ext = C_ext0*delta + (1.d0-delta)*C_ext1
       C_abs = C_abs0*delta + (1.d0-delta)*C_abs1

       !*** Average over size
       Cext_av = 0.d0
       Cabs_av = 0.d0
       der_CexAv_rem =0.d0
       der_CabAv_rem = 0.d0
       der_CexAv_imm = 0.d0
       der_CabAv_imm = 0.d0
       Cext_av0 = 0.d0
       Cabs_av0 = 0.d0
       Cext_av1 = 0.d0
       Cabs_av1 = 0.d0
       do i = 1, n_xgr0 
          Cext_av = Cext_av+C_ext(i)*V_distr(i)
          Cabs_av = Cabs_av+C_abs(i)*V_distr(i)
          Cext_av0 = Cext_av0+C_ext0(i)*V_distr(i)
          Cabs_av0 = Cabs_av0+C_abs0(i)*V_distr(i)
          Cext_av1 = Cext_av1+C_ext1(i)*V_distr(i)
          Cabs_av1 = Cabs_av1+C_abs1(i)*V_distr(i)
          !*** Derivative wrt delta, Re(m), Im(m)
          der_CexAv_rem = der_CexAv_rem + &
               (delta*der_Cex_rem(i)+(1.d0-delta)*der_Cex_rem1(i))*V_distr(i)
          der_CexAv_imm = der_CexAv_imm + &
               (delta*der_Cex_imm(i)+(1.d0-delta)*der_Cex_imm1(i))*V_distr(i)
          der_CabAv_rem = der_CabAv_rem + &
               (delta*der_Cab_rem(i)+(1.d0-delta)*der_Cab_rem1(i))*V_distr(i)
          der_CabAv_imm = der_CabAv_imm + &
               (delta*der_Cab_imm(i)+(1.d0-delta)*der_Cab_imm1(i))*V_distr(i)
       end do

       !*** scattering coefficient averaged over size
       Csc_av = Cext_av - Cabs_av

       !*** Derivative wrt delta
       der_CexAv_delt = Cext_av0 - Cext_av1
       der_CabAv_delt = Cabs_av0 - Cabs_av1
       Der_Cext(1) = Cext_av
       Der_Cext(2) = der_CexAv_delt
       Der_Cext(3) = der_CexAv_rem
       Der_Cext(4) = der_CexAv_imm
       Der_Cabs(1) = Cabs_av
       Der_Cabs(2) = der_CabAv_delt
       Der_Cabs(3) = der_CabAv_rem
       Der_Cabs(4) = der_CabAv_imm

       !*** Derivative wrt distr. func.
       der_Cex_df = C_ext
       der_Cab_df = C_abs
       dcabs_full(1:n_xgr0) = C_abs
       dcsca_full(1:n_xgr0) = C_ext - C_abs
       dcabs_full(n_xgr0+1) = der_CabAv_rem
       dcabs_full(n_xgr0+2) = der_CabAv_imm
       dcabs_full(n_xgr0+3) = der_CabAv_delt
       dcsca_full(n_xgr0+1) = der_CexAv_rem - der_CabAv_rem
       dcsca_full(n_xgr0+2) = der_CexAv_imm - der_CabAv_imm
       dcsca_full(n_xgr0+3) = der_CexAv_delt - der_CabAv_delt

       !*** phase function
       !*** spheroids
       sp_el = 0
       call F_read_int( &
            aero_lut, &
            sp_el, &
            Re_m, Im_m, i_re12, i_im12, &
            Mel_int0,&
            Der_Mel_rem0, Der_Mel_imm0)
       !*** ellipsoids
       sp_el = 1
       call F_read_int( &
            aero_lut, &
            sp_el, &
            Re_m, Im_m, i_re12, i_im12, &
            Mel_int1,&
            Der_Mel_rem1, Der_Mel_imm1)

       !*** For mixture of spheroids and ellipsoids, average over size

       !*** For mixture of spheroids and ellipsoids, integrate over size
       Mel_rem = 0.d0
       Mel_imm =  0.d0
       Mel_sp = 0.d0
       Mel_el = 0.d0
       do i_el = 1, n_el
          do i = 1, n_xgr0 
             do i_tet = 1, n_tet0
                Mel_sp(i_el,i_tet) = Mel_sp(i_el,i_tet) + Mel_int0(i_el,i_tet,i)*V_distr(i)
                Mel_el(i_el,i_tet) = Mel_el(i_el,i_tet) + Mel_int1(i_el,i_tet,i)*V_distr(i)
                !*** Derivative wrt Re(m)
                Mel_rem(i_el,i_tet) = Mel_rem(i_el,i_tet) + &
                     (Der_Mel_rem0(i_el,i_tet,i)*delta + (1.d0-delta)*Der_Mel_rem1(i_el,i_tet,i))*V_distr(i)
                !*** Derivative wrt Im(m)
                Mel_imm(i_el,i_tet) = Mel_imm(i_el,i_tet) + &
                     (Der_Mel_imm0(i_el,i_tet,i)*delta + (1.d0-delta)*Der_Mel_imm1(i_el,i_tet,i))*V_distr(i)
             enddo
          enddo
       enddo

       Mel_av = (Mel_sp*delta+(1.d0-delta)*Mel_el)
       mel_rem = mel_rem/csc_av - mel_av / csc_av**2 * dcsca_full(n_xgr0+1)
       mel_imm = mel_imm/csc_av - mel_av / csc_av**2 * dcsca_full(n_xgr0+2)

       !*** Derivative wrt distr.func.
       Mel_int = Mel_int0*delta + (1.d0-delta)*Mel_int1
       do i = 1, n_xgr0
          Mel_int(:,:,i) = Mel_int(:,:,i)/csc_av - mel_av / csc_av**2 * dcsca_full(i)
       enddo

       !*** Derivative wrt delta
       Mel_delt = Mel_sp - Mel_el
       Mel_delt = mel_delt/csc_av - mel_av / csc_av**2 * dcsca_full(n_xgr0+3)
       Mel_av = Mel_av/Csc_av

       !*** calculation Gaus-Leg quadr. weights			!			<============D.S.====================        
       n_tet = n_tet0 						!			<============D.S.====================

       allocate(M_el(n_el, n_tet), stat=ierr)
       allocate(Gw_t(n_tet), stat=ierr) 				!			<============D.S.====================
       allocate(tet(n_tet)) 					!			<============D.S.====================
       allocate(mu_t(n_tet)) 					!			<============D.S.====================

       !			<============D.S.====================
       mu_t = aero_lut%mu_t_dg 					!			<============D.S.====================
       Gw_t = aero_lut%Gw_t_dg 					!			<============D.S.====================


       !***Force the interpolated, size-integrated phase function to be normalised
       Mel_av(1,n_tet)   = (2.D0-sum(Mel_av(1, 1:n_tet-1)*Gw_t(1:n_tet-1)))/Gw_t(n_tet)
       Mel_rem(1,n_tet)  = -sum(Mel_rem(1, 1:n_tet-1)*Gw_t(1:n_tet-1))/Gw_t(n_tet)
       Mel_imm(1,n_tet)  = -sum(Mel_imm(1, 1:n_tet-1)*Gw_t(1:n_tet-1))/Gw_t(n_tet)
       Mel_delt(1,n_tet) = -sum(Mel_delt(1, 1:n_tet-1)*Gw_t(1:n_tet-1))/Gw_t(n_tet)
       do i = 1, n_xgr0      
          Mel_int(1,n_tet, i) =  -sum(Mel_int(1,1:n_tet-1, i)*Gw_t(1:n_tet-1))/Gw_t(n_tet)
       enddo

       M_el = Mel_av												!			<============D.S.====================

       allocate(Der_rem(n_el, n_tet), stat=ierr)
       Der_rem = Mel_rem												!			<============D.S.====================

       allocate(Der_imm(n_el, n_tet), stat=ierr)
       Der_imm = Mel_imm												!			<============D.S.====================

       allocate(Der_delta(n_el, n_tet), stat=ierr)
       Der_delta = Mel_delt											!			<============D.S.====================

       allocate(Der_DestF(n_el, n_tet, n_xgr0), stat=ierr)
       do i = 1, n_xgr0
          Der_DestF(:,:,i) = Mel_int(:,:,i)									!			<============D.S.====================
       enddo

       allocate(der_F(n_el, n_tet, n_xgr0+3))

       der_F(:,:,1:n_xgr0) = Der_DestF
       der_F(:,:,n_xgr0+1) = Der_rem
       der_F(:,:,n_xgr0+2) = Der_imm
       der_F(:,:,n_xgr0+3) = Der_delta

       call devel(&
            mu_t,&
            gw_t,&
            M_el,&
            coefs,&
            der_F,&
            dcoefs)

       deallocate( &
            Der_DestF, &
            Der_rem, &
            Der_imm, &
            Der_delta, &  
            tet, &
            mu_t, &
            Gw_t, &
            der_F, stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'SPEL_EXP: memory deallocation error'
          ierr = ierr_deall
          goto 999
       endif

       return 

999    continue
       call stopretrieval(message)

     end subroutine SpEl_exp

     !------------------------------------------------------------------------
     !> Calculate scattering matrix for single scattering
     !------------------------------------------------------------------------
     subroutine scatmat_ssc(&
          aero_lut, &
          Re_m,&
          Im_m,&
          delta,&
          i_re12,&
          i_im12,&
          V_distr, &
          Mel_av,&
          der_F,&
          csc_av,&
          cabs_av, &
          dcsca_full,&
          dcabs_full)
       !*** input
       type(Mie_lut), intent(in) :: aero_lut
       integer, intent(in):: i_re12(2), i_im12(2)
       real(double), intent (in):: Re_m, Im_m
       real(double), intent(in):: delta, V_distr(n_xgr0)
       !*** output
       real(double), dimension(n_el, n_tet0), intent(out) :: Mel_av
       real(double), dimension(n_el, n_tet0, n_xgr0+4), intent(out) :: der_F
       real(double), dimension(:), intent(out) :: dcsca_full, dcabs_full
       real(double), intent(out) :: Cabs_av, Csc_av
       !*** local
       real(double) :: Der_Cext(4), Der_Cabs(4), Der_Cex_df(n_xgr0), Der_Cab_df(n_xgr0)
       real(double)  Cext_av
       integer :: i, i_el, sp_el, i_tet
       real(double), dimension(n_xgr0) :: C_ext, C_abs, C_ext0, C_abs0, C_ext1, C_abs1,&
            der_Cex_rem, der_Cab_rem, der_Cex_imm, der_Cab_imm,&
            der_Cex_rem1, der_Cab_rem1, der_Cex_imm1, der_Cab_imm1                    
       real(double), dimension(n_el, n_tet0, n_xgr0) :: Mel_int, Mel_int0, Der_Mel_rem0, Der_Mel_imm0, &
            Mel_int1, Der_Mel_rem1, Der_Mel_imm1
       real(double) :: der_CexAv_delt, der_CabAv_delt, der_CexAv_rem, der_CabAv_rem, &
            der_CexAv_imm, der_CabAv_imm, Cext_av0, Cext_av1, Cabs_av0, Cabs_av1
       real(double), dimension(n_el, n_tet0) :: Mel_sp, Mel_el, Mel_delt, Mel_rem, Mel_imm
       !----------------------------------------------

       !*** cross-sections
       !*** spheres (interpolating in LUT)
       sp_el = 0
       call C_read_int(aero_lut, sp_el,Re_m,Im_m,i_re12,i_im12,C_ext0,C_abs0, &
            der_Cex_rem,der_Cab_rem,der_Cex_imm,der_Cab_imm)

       !*** elipsoids (interpolating in LUT)
       sp_el = 1
       call C_read_int(aero_lut, sp_el,Re_m,Im_m,i_re12,i_im12,C_ext1,C_abs1, &
            der_Cex_rem1,der_Cab_rem1,der_Cex_imm1,der_Cab_imm1)

       !*** extinction and absorption cross-sections
       C_ext = C_ext0*delta+(1.d0-delta)*C_ext1
       C_abs = C_abs0*delta+(1.d0-delta)*C_abs1

       !*** Average over size
       Cext_av = 0.d0
       Cabs_av = 0.d0
       der_CexAv_rem = 0.d0
       der_CabAv_rem = 0.d0
       der_CexAv_imm = 0.d0
       der_CabAv_imm = 0.d0
       Cext_av0 = 0.d0
       Cabs_av0 = 0.d0
       Cext_av1 = 0.d0
       Cabs_av1 = 0.d0
       do i = 1, n_xgr0 
          Cext_av = Cext_av+C_ext(i)*V_distr(i)
          Cabs_av = Cabs_av+C_abs(i)*V_distr(i)
          Cext_av0 = Cext_av0+C_ext0(i)*V_distr(i)
          Cabs_av0 = Cabs_av0+C_abs0(i)*V_distr(i)
          Cext_av1 = Cext_av1+C_ext1(i)*V_distr(i)
          Cabs_av1 = Cabs_av1+C_abs1(i)*V_distr(i)
          !*** Derivative wrt delta Re(m), Im(m)
          der_CexAv_rem = der_CexAv_rem + &
               (delta*der_Cex_rem(i)+(1.d0-delta)*der_Cex_rem1(i))*V_distr(i)
          der_CexAv_imm = der_CexAv_imm + &
               (delta*der_Cex_imm(i)+(1.d0-delta)*der_Cex_imm1(i))*V_distr(i)
          der_CabAv_rem = der_CabAv_rem + &
               (delta*der_Cab_rem(i)+(1.d0-delta)*der_Cab_rem1(i))*V_distr(i)
          der_CabAv_imm = der_CabAv_imm + &
               (delta*der_Cab_imm(i)+(1.d0-delta)*der_Cab_imm1(i))*V_distr(i)
       enddo

       !*** scattering coeffcient averaged over size
       Csc_av = Cext_av-Cabs_av

       !*** Derivative wrt delta
       der_CexAv_delt = Cext_av0-Cext_av1
       der_CabAv_delt = Cabs_av0-Cabs_av1
       Der_Cext(1) = Cext_av
       Der_Cext(2) = der_CexAv_delt
       Der_Cext(3) = der_CexAv_rem
       Der_Cext(4) = der_CexAv_imm
       Der_Cabs(1) = Cabs_av
       Der_Cabs(1) = Cabs_av
       Der_Cabs(2) = der_CabAv_delt
       Der_Cabs(3) = der_CabAv_rem
       Der_Cabs(4) = der_CabAv_imm

       !*** Derivative wrt distr. func.
       der_Cex_df=C_ext
       der_Cab_df=C_abs
       dcabs_full(1:n_xgr0) = C_abs
       dcsca_full(1:n_xgr0) = C_ext-C_abs
       dcabs_full(n_xgr0+1) = der_CabAv_rem
       dcabs_full(n_xgr0+2) = der_CabAv_imm
       dcabs_full(n_xgr0+3) = der_CabAv_delt
       dcsca_full(n_xgr0+1) = der_CexAv_rem-der_CabAv_rem
       dcsca_full(n_xgr0+2) = der_CexAv_imm-der_CabAv_imm
       dcsca_full(n_xgr0+3) = der_CexAv_delt-der_CabAv_delt

       !*** phase functions
       !*** spheres
       sp_el = 0
       call F_read_int(aero_lut, sp_el, Re_m, Im_m, i_re12, i_im12, Mel_int0,&
            Der_Mel_rem0, Der_Mel_imm0)
       !*** ellipsoids
       sp_el = 1
       call F_read_int(aero_lut, sp_el, Re_m, Im_m, i_re12, i_im12, Mel_int1,&
            Der_Mel_rem1, Der_Mel_imm1)

       !*** Averaging over size
       Mel_rem = 0.d0
       Mel_imm = 0.d0
       Mel_sp = 0.d0
       Mel_el = 0.d0
       do i_el = 1, n_el
          do i = 1, n_xgr0 
             do i_tet = 1, n_tet0
                Mel_sp(i_el, i_tet) = Mel_sp(i_el,i_tet) + Mel_int0(i_el,i_tet,i)*V_distr(i)
                Mel_el(i_el, i_tet) = Mel_el(i_el,i_tet) + Mel_int1(i_el,i_tet,i)*V_distr(i)
                !*** Derivative wrt Re(m)
                Mel_rem(i_el,i_tet) = Mel_rem(i_el,i_tet) + &
                     (Der_Mel_rem0(i_el,i_tet,i)*delta+(1.d0-delta)*Der_Mel_rem1(i_el,i_tet,i))*V_distr(i)
                !*** Derivative wrt Im(m)
                Mel_imm(i_el,i_tet) = Mel_imm(i_el,i_tet) + &
                     (Der_Mel_imm0(i_el,i_tet,i)*delta+(1.d0-delta)*Der_Mel_imm1(i_el,i_tet,i))*V_distr(i)
             enddo
          enddo
       enddo

       Mel_av = (Mel_sp*delta+(1.d0-delta)*Mel_el)
       mel_rem = mel_rem/csc_av - mel_av / csc_av**2 * dcsca_full(n_xgr0+1)
       mel_imm = mel_imm/csc_av - mel_av / csc_av**2 * dcsca_full(n_xgr0+2)

       !*** Derivative wrt distr.func.
       Mel_int = Mel_int0*delta + (1.d0-delta)*Mel_int1
       do i = 1, n_xgr0
          Mel_int(:,:,i) = Mel_int(:,:,i)/csc_av-mel_av / csc_av**2 * dcsca_full(i)
       enddo

       !*** Derivative wrt delta
       Mel_delt = Mel_sp - Mel_el
       Mel_delt = mel_delt/csc_av - mel_av / csc_av**2 * dcsca_full(n_xgr0+3)
       Mel_av = Mel_av/Csc_av
       der_F(:,:,1:n_xgr0) = Mel_int
       der_F(:,:,n_xgr0+1) = Mel_rem
       der_F(:,:,n_xgr0+2) = Mel_imm
       der_F(:,:,n_xgr0+3) = Mel_delt

     end subroutine scatmat_ssc

     !------------------------------------------------------------------------
     !> Calculate scattering and absorption coefficients by interpolating in LUT,
     !! averaged over size by using volume size distribution
     !------------------------------------------------------------------------
     subroutine SpEl_exp_xs( &
          aero_lut, &
          Re_m, &
          Im_m, &
          delta,  &
          i_re12, &
          i_im12, &
          V_distr, &
          csc_av, &
          cabs_av)
       !*** Input
       type(Mie_lut), intent(in) :: aero_lut
       integer, intent(in):: i_re12(2), i_im12(2)
       real(double), intent (in):: Re_m, Im_m
       real(double), intent(in):: delta, V_distr(n_xgr0)
       !*** Output
       real(double), intent(out) :: Cabs_av, Csc_av
       !*** Local variables
       integer :: i, sp_el
       real(double), dimension(n_xgr0) :: C_ext, C_abs, C_ext0, C_abs0, C_ext1, C_abs1,&
            der_Cex_rem, der_Cab_rem, der_Cex_imm, der_Cab_imm,&
            der_Cex_rem1, der_Cab_rem1, der_Cex_imm1, der_Cab_imm1
       real(double) :: Cext_av
       !----------------------------------------------

       !*** Get cross-sections
       !*** spheres (interpolating in LUT)
       sp_el = 0
       call C_read_int(aero_lut, sp_el, Re_m, Im_m, i_re12, i_im12, C_ext0, C_abs0, &
            der_Cex_rem, der_Cab_rem, der_Cex_imm, der_Cab_imm)

       !*** elipsoids (interpolating in LUT)
       sp_el = 1
       call C_read_int(aero_lut, sp_el, Re_m, Im_m, i_re12, i_im12, C_ext1, C_abs1, &
            der_Cex_rem1, der_Cab_rem1, der_Cex_imm1, der_Cab_imm1)

       !*** extinction and absorption coeffcients
       C_ext = C_ext0*delta+(1.d0-delta)*C_ext1
       C_abs = C_abs0*delta+(1.d0-delta)*C_abs1

       Cext_av = 0.d0
       Cabs_av = 0.d0

       !*** Average over size
       do i = 1, n_xgr0 
          Cext_av = Cext_av + C_ext(i)*V_distr(i)
          Cabs_av = Cabs_av + C_abs(i)*V_distr(i)
       end do

       !*** scattering coeffcient averaged over size
       Csc_av = Cext_av - Cabs_av

     end subroutine SpEl_exp_xs

     !------------------------------------------------------------------------
     !> Calculate derivatives of scattering and absorption cross-sections and
     !! expansion coefficients for phase function, used for multiple scattering
     !------------------------------------------------------------------------
     subroutine modes_calc(&
          aero_lut, &
          id_size,&
          aerosol_pars,&
          lambd,&
          csca_mode,&
          cabs_mode, &
          dcsca_mode,&
          dcabs_mode,&
          coefs,&
          dcoefs_mode, & 
          ierr)
       !*** Input
       type(Mie_lut), intent(in) :: aero_lut
       real(double), intent(in) :: lambd
       real(double), dimension(:), intent(in) :: aerosol_pars
       integer, intent(in) :: id_size
       !*** Output 
       real(double), intent(out)  :: csca_mode, cabs_mode 
       real(double), dimension(npar_mie), intent(out) :: dcsca_mode, dcabs_mode
       real(double), dimension(nstokes, nstokes, 0:maxleg), intent(out) :: coefs
       real(double), dimension(nstokes, nstokes, 0:maxleg, npar_mie), intent(out) :: dcoefs_mode
       integer, intent(out) :: ierr
       !*** Local variables
       real(double), dimension(n_xgr0+4) ::  dtaus_full, dtaua_full
       real(double), dimension(nstokes, nstokes, 0:maxleg, n_xgr0+4):: dcoefs_full
       real(double), dimension(4, 4, 0:maxleg, n_xgr0+3) :: dcoefs_full_tmp 
       real(double), dimension(4, 4, 0:maxleg) :: coefs_tmp
       real(double) :: Re_m, Im_m
       real(double) :: r_eff , v_eff, delta_sp, aer_col
       real(double), dimension(n_xgr0+3) :: dcsca_full, dcabs_full
       real(double) :: csca, cabs
       integer :: i_re12(2), i_im12(2), i, j, l 
       real(double), dimension(n_xgr0):: V_distr, Der_Vdist_reff, Der_Vdist_veff   
       !---------------------------------------------------------------------------------------

       !*** Initialize error identifier
       ierr = 0

       r_eff = aerosol_pars(1) 
       v_eff = aerosol_pars(2)
       Re_m = aerosol_pars(3)
       Im_m = DABS(aerosol_pars(4))
       aer_col = aerosol_pars(5)
       delta_sp = aerosol_pars(6)

       !*** No need to compute optical properties for zero reff/veff/column 
       dcsca_mode = 0.d0
       dcabs_mode = 0.d0
       coefs = 0.d0
       dcoefs_mode = 0.d0
       !       call writelog('MODES_CALC: reff, veff or aer_col < 1.d-10', 5)
       if (r_eff < 1.d-10 .or. v_eff < 1.d-10 .or. aer_col < 1.d-10) then
          return
       endif

       call RemIm_index(aero_lut, Re_m, Im_m, i_re12, i_im12)

       !*** Calculation of the distribution function and its derivatives
       call VolDistrFunc(id_size, &
            lambd, &
            r_eff, &
            v_eff, &
            aero_lut%x0_in,&
            V_distr(:), &
            Der_Vdist_reff(:), &
            Der_Vdist_veff(:))

       !*** Calculation of scattering characteristics for each mode
       call SpEl_exp( &
            aero_lut, & 
            Re_m, &
            Im_m, &
            delta_sp,&
            i_re12, &
            i_im12, &
            V_distr(:), &
            coefs_tmp(:,:,:),&
            dcoefs_full_tmp(:,:,:,:),&
            csca,&
            cabs,&
            dcsca_full,&
            dcabs_full, &
            ierr)
       if (ierr .ne. 0) return 

       csca = (aero_lut%wave_ref_lut/lambd)*csca*1.D-9
       cabs = (aero_lut%wave_ref_lut/lambd)*cabs*1.D-9  

       dcsca_full = (aero_lut%wave_ref_lut/lambd)*dcsca_full*1.D-9
       dcabs_full = (aero_lut%wave_ref_lut/lambd)*dcabs_full*1.D-9

       dtaus_full(1:n_xgr0+3) = dcsca_full*aer_col
       dtaua_full(1:n_xgr0+3) = dcabs_full*aer_col

       dtaus_full(n_xgr0+4) = csca
       dtaua_full(n_xgr0+4) = cabs
       dcoefs_full(:, :, :, n_xgr0+4) = 0.D0

       do i = 1, nstokes
          do j = 1, nstokes
             coefs(i, j, :) = coefs_tmp(i, j, :)
             dcoefs_full(i, j, :, 1:n_xgr0+3) = dcoefs_full_tmp(i, j, :, :)
             do l = 0, maxleg
                dcoefs_mode(i, j, l, 1) = &
                     sum(dcoefs_full(i, j, l, 1:n_xgr0)*Der_Vdist_reff(:))
                dcoefs_mode(i, j, l, 2) = &
                     sum(dcoefs_full(i, j, l, 1:n_xgr0)*Der_Vdist_veff(:))
             enddo
             dcoefs_mode(i, j, :, 3) = dcoefs_full(i, j, :, n_xgr0+1)
             dcoefs_mode(i, j, :, 4) = dcoefs_full(i, j, :, n_xgr0+2)
             dcoefs_mode(i, j, :, 6) =  dcoefs_full(i, j, :, n_xgr0+3)
          enddo
       enddo
       dcoefs_mode(:, :, :, 5) = 0.D0

       dcsca_mode(1) = sum(dcsca_full(1:n_xgr0)*Der_Vdist_reff(:)) 
       dcabs_mode(1) = sum(dcabs_full(1:n_xgr0)*Der_Vdist_reff(:)) 
       dcsca_mode(2) = sum(dcsca_full(1:n_xgr0)*Der_Vdist_veff(:)) 
       dcabs_mode(2) = sum(dcabs_full(1:n_xgr0)*Der_Vdist_veff(:)) 
       dcsca_mode(3) = dcsca_full(n_xgr0+1)
       dcabs_mode(3) = dcabs_full(n_xgr0+1)
       dcsca_mode(4) = dcsca_full(n_xgr0+2)
       dcabs_mode(4) = dcabs_full(n_xgr0+2)
       dcsca_mode(5) = csca/aer_col
       dcabs_mode(5) = cabs/aer_col
       dcsca_mode(6) = dcsca_full(n_xgr0+3)  
       dcabs_mode(6) = dcabs_full(n_xgr0+3)

       csca_mode = csca
       cabs_mode = cabs

     end subroutine modes_calc

     !------------------------------------------------------------------------
     !> Calculate derivatives of scattering and absorption cross-sections 
     !! and phase function for single scattering
     !------------------------------------------------------------------------
     subroutine modes_calc_ssc(&
          aero_lut, &
          id_size,&
          aerosol_pars,&
          lambd,&
          scat_angle, &
          dcsca_mode,&
          dcabs_mode,&
          zss_lay,&
          d_zss_lay)
       implicit none
       !*** Input
       type(Mie_lut), intent(in) :: aero_lut
       integer, intent(in) :: id_size
       real(double), dimension(:), intent(in) :: aerosol_pars
       real(double), intent(in) :: lambd, scat_angle
       !*** output
       real(double), dimension(npar_mie), intent(out) :: dcsca_mode, dcabs_mode
       real(double), dimension(nstokes, nstokes) :: zss_lay
       real(double), dimension(nstokes, nstokes, npar_mie) :: d_zss_lay
       !*** Local variables
       real(double), dimension(n_el, n_tet0) :: F_mode
       real(double), dimension(n_el, n_tet0, npar_mie) :: dF_mode
       real(double) :: Re_m, Im_m
       real(double) :: r_eff, v_eff, delta_sp, aer_col
       real(double), dimension(n_el, n_tet0, n_xgr0+4) :: der_F
       real(double), dimension(n_xgr0+3) :: dcsca_full, dcabs_full
       real(double) :: csca, cabs
       integer ::  i_re12(2), i_im12(2)
       real(double), dimension(n_xgr0) :: V_distr, Der_Vdist_reff, Der_Vdist_veff
       integer :: i_el, i_tet, imin, iper, n
       !---------------------------------------------------------------------------------------

       r_eff = aerosol_pars(1)
       v_eff = aerosol_pars(2)
       Re_m = aerosol_pars(3)
       Im_m = DABS(aerosol_pars(4))
       aer_col = aerosol_pars(5)
       delta_sp = aerosol_pars(6)
       dcsca_mode = 0.d0
       dcabs_mode = 0.d0
       zss_lay = 0.d0
       d_zss_lay = 0.d0
       !*** No need to compute optical properties for zero reff/veff/column 
       if (r_eff < 1.d-10 .or. v_eff < 1.d-10 .or. aer_col < 1.d-10) then
          !         call writelog('MODES_CALC_SS: reff, veff or aer_col < 1.d-10', 5)
          return
       endif

       call RemIm_index(aero_lut, Re_m, Im_m, i_re12, i_im12)

       !*** Calculation of the distribution function end its derivatives
       call VolDistrFunc(id_size, &
            lambd, &
            r_eff, &
            v_eff, &
            aero_lut%x0_in,&
            V_distr, &
            Der_Vdist_reff, &
            Der_Vdist_veff)

       !*** Calculation of scattering characteristics for each mode
       call scatmat_ssc(&
            aero_lut, &
            Re_m,&
            Im_m,&
            delta_sp, &
            i_re12,&
            i_im12,&
            V_distr, &
            F_mode,&
            der_F,&
            csca,&
            cabs, &
            dcsca_full,&
            dcabs_full)

       csca = (aero_lut%wave_ref_lut/lambd)*csca*1.D-9
       cabs = (aero_lut%wave_ref_lut/lambd)*cabs*1.D-9

       dcsca_full = (aero_lut%wave_ref_lut/lambd)*dcsca_full*1.D-9
       dcabs_full = (aero_lut%wave_ref_lut/lambd)*dcabs_full*1.D-9

       do i_el = 1, n_el
          do i_tet = 1, n_tet0
             dF_mode(i_el,i_tet,1) = sum(der_F(i_el,i_tet,1:n_xgr0)*Der_Vdist_reff(:))
             dF_mode(i_el,i_tet,2) = sum(der_F(i_el,i_tet,1:n_xgr0)*Der_Vdist_veff(:))
             df_mode(i_el,i_tet,3) = der_F(i_el,i_tet,n_xgr0+1)
             df_mode(i_el,i_tet,4) = der_F(i_el,i_tet,n_xgr0+2)
             df_mode(i_el,i_tet,5) = 0.D0
             df_mode(i_el,i_tet,6) = der_F(i_el,i_tet,n_xgr0+3)
          enddo
       enddo

       !*** interpolate F_mode and df mode to the actual scattering angle and 
       !*** fill the arrays zss_lay and d_zss_lay
       imin = minval(minloc(DABS(scat_angle-aero_lut%tet0_dg),scat_angle.ge.aero_lut%tet0_dg))
       !*** HH: force imin to remain within array bounds, perhaps return with error is better? 
       !    imin = max(1,min(imin, n_tet0-1))

       do iper = 1, n_el
          if (index_ist_sp(iper) .le. nstokes .and. index_jst_sp(iper) .le. nstokes) then
             zss_lay(index_ist_sp(iper),index_jst_sp(iper)) = &
                  F_mode(iper, imin) + &
                  (F_mode(iper, imin+1) - F_mode(iper, imin))/ &
                  (aero_lut%tet0_dg(imin+1) - aero_lut%tet0_dg(imin)) * &
                  (scat_angle - aero_lut%tet0_dg(imin))                 
             d_zss_lay(index_ist_sp(iper),index_jst_sp(iper),:) = &
                  dF_mode(iper, imin,:) + &
                  (dF_mode(iper, imin+1,:) - dF_mode(iper, imin,:))/ &
                  (aero_lut%tet0_dg(imin+1) - aero_lut%tet0_dg(imin))* &
                  (scat_angle - aero_lut%tet0_dg(imin))
          endif
       enddo

       do n = 2, min(2,nstokes)
          zss_lay(1,n) = zss_lay(n,1)
          d_zss_lay(1,n,:) = d_zss_lay(n,1,:)
       enddo
       !    zss_lay(1,2) = zss_lay(2,1)
       !    d_zss_lay(1,2,:) = d_zss_lay(2,1,:)



       dcsca_mode(1) = sum(dcsca_full(1:n_xgr0)*Der_Vdist_reff(:)) 
       dcabs_mode(1) = sum(dcabs_full(1:n_xgr0)*Der_Vdist_reff(:))  
       dcsca_mode(2) = sum(dcsca_full(1:n_xgr0)*Der_Vdist_veff(:))  
       dcabs_mode(2) = sum(dcabs_full(1:n_xgr0)*Der_Vdist_veff(:))  
       dcsca_mode(3) = dcsca_full(n_xgr0+1)  
       dcabs_mode(3) = dcABS_full(n_xgr0+1) 
       dcsca_mode(4) = dcsca_full(n_xgr0+2)  
       dcabs_mode(4) = dcABS_full(n_xgr0+2)  
       dcsca_mode(5) = csca/aer_col
       dcabs_mode(5) = cabs/aer_col
       dcsca_mode(6) = dcsca_full(n_xgr0+3)  
       dcabs_mode(6) = dcABS_full(n_xgr0+3)

     end subroutine modes_calc_ssc

     !------------------------------------------------------------------------
     !> Calculate absorption and scattering cross-sections 
     !! for given wavelength and size distribution
     !------------------------------------------------------------------------
     subroutine modes_calc_xs(&
          aero_lut, &
          id_size,&
          aerosol_pars,&
          lambd,&
          csca_mode,&
          cabs_mode)
       !*** Input
       type(Mie_lut), intent(in) :: aero_lut
       integer, intent(in) :: id_size  
       real(double), intent(in) :: lambd
       real(double), dimension(:), intent(in) :: aerosol_pars
       !*** Output
       real(double), intent(out)  :: csca_mode, cabs_mode 
       !*** Local variables
       real(double) :: Re_m, Im_m
       real(double) ::  r_eff, v_eff, delta_sp, aer_col
       real(double) :: csca, cabs
       integer ::  i_re12(2), i_im12(2)
       real(double), dimension(n_xgr0):: V_distr, Der_Vdist_reff, Der_Vdist_veff
       !---------------------------------------------------------------------------------------
       r_eff = aerosol_pars(1)
       v_eff = aerosol_pars(2)
       Re_m = aerosol_pars(3)
       Im_m = DABS(aerosol_pars(4))
       aer_col = aerosol_pars(5)
       delta_sp = aerosol_pars(6)

       !*** No need to compute optical properties for zero reff/veff/column 
       if (r_eff < 1.d-10 .or. v_eff < 1.d-10) then
          csca_mode = 0.d0
          cabs_mode = 0.d0
          !         call writelog('MODES_CALC_XS: reff, veff  < 1.d-10', 5)
          return
       endif

       !*** get array indices of refractive index 
       call RemIm_index(aero_lut, Re_m, Im_m, i_re12, i_im12)

       !*** Calculation of the volume size distribution function and its derivatives
       call VolDistrFunc( &
            id_size, &
            lambd, &
            r_eff, &
            v_eff, &
            aero_lut%x0_in,&
            V_distr(:), &          ! distribution function
            Der_Vdist_reff(:), &   ! derivative of distr func. wrt r_eff 
            Der_Vdist_veff(:))     ! derivative of distr func. wrt v_eff 

       !*** Calculation of scattering characteristics
       call SpEl_exp_xs( &
            aero_lut, &
            Re_m, &
            Im_m, &
            delta_sp,&
            i_re12, &
            i_im12,&
            V_distr(:), &
            csca,&
            cabs)

       csca_mode = (aero_lut%wave_ref_lut/lambd)*csca*1.D-9
       cabs_mode = (aero_lut%wave_ref_lut/lambd)*cabs*1.D-9

     end subroutine modes_calc_xs


     !------------------------------------------------------------------------
     !> Due to their large sizes the instantaneous dimensioning of scattering 
     !  components of brief aerosol optical properties LUT may dump when 
     !  assigned instantaneously (on the stack). When assigned as allocatable 
     !  (dynamic memory on the heap) there is no problem. Therefore these 
     !! components are declared allocatable.
     !------------------------------------------------------------------------
     subroutine allocate_scattering(aero_lut, ierr)
       type(Mie_lut), intent(inout) :: aero_lut
       integer, intent(out) :: ierr

       ierr = 0
       if (.not. allocated(aero_lut%M_el0_in_el)) then
          allocate(aero_lut%M_el0_in_el(n_tet0, n_xgr0, n_Rem0, n_Imm0, n_el), &
               aero_lut%M_el0_in_sph(n_tet0, n_xgr0, n_Rem0, n_Imm0, n_el), stat=ierr)
          allocate(aero_lut%M_el0_el(n_tet0, n_xgr0, n_Rem0, n_Imm0, n_el), &
               aero_lut%M_el0_sph(n_tet0, n_xgr0, n_Rem0, n_Imm0, n_el), stat=ierr)    
       endif

     end subroutine allocate_scattering

     !------------------------------------------------------------------------

   end module OpticM_module
