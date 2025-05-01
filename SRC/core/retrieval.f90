!------------------------------------------------------------------------------
!> @brief RemoTeC retrieval algorithm
!> @todo Set ierr
!------------------------------------------------------------------------------
module retrieval_module
   use header_module
   use profile_inversion_module, only: &
      retrieval_data, profile_inversion, absorbers, &
      window_ini, settings_flags, file_paths, read_settings, read_win_xsdb, &
      Mie_lut, cirrus_table, read_aerosol_netcdf, read_cirrus_netcdf, &
      aero, get_aerosol_properties_lognormal, window_spectrum
   use spectrum_internal_module, only: spectrum, instrument_response, instrument_interface
   use atmosphere_internal_module, only: atmospheric_scenario, atmosphere_interpolate, altitude_grid
   implicit none
   private

!*** Public types/procedures
   public :: aero, atmospheric_scenario, spectrum, instrument_response, retrieval_data, &
              Mie_lut, cirrus_table, window_ini, settings_flags, file_paths, altitude_grid
   public :: retrieval, get_absorbers, read_settings, read_win_xsdb, read_aerosol_netcdf, read_cirrus_netcdf


   contains
!------------------------------------------------------------------------------
!> @brief
!> @details This routine retrieves specified absorber(s)
!> @param[in] measurement         spectral measurement
!> @param[in] meta                metadata (date, position etc)
!> @param[in] response_in         instrument response function
!> @param[in] atm_scenario        meteological data
!> @param[in] grid                atmospheric grids for retrieval, RT, and XS
!> @param[in] flags               flags for retrieval settings
!> @param[in] window_ini          window-dependent retrieval settings
!> @param[in] aerosol_ini         initial guess of aerosol parameters and aerosol retrieval settings
!> @param[in] aero_lut            aerosol optical properties lookup table
!> @param[in] cirrus_lut          cirrus optical properties lookup table
!> @param[out] retrieval_output   output data
!> @param[out] ierr               error identifier: 0=normal, >0=fatal error:abort retrieval, <0=error: go to next ground pixel
!> @param[in] irr_meas_hi         high resolution measured solar spectrum, eg from deconvolution (optional)
!------------------------------------------------------------------------------
     subroutine retrieval( &
          measurement, meta, response_in, atm_scenario, &
          grid, &
          flag, &
          win_ini, aerosol_ini, &
          aero_lut, cirrus_lut, &
         line_number, &
          retrieval_output, ierr, meas_irr_hi)
       !*** Input
       type(spectrum), dimension(:), intent(in) :: measurement
       type(instrument_response), dimension(:), intent(in) :: response_in
       type(atmospheric_scenario), intent(inout) :: atm_scenario
       type(altitude_grid), intent(in) :: grid
       type(Mie_lut), intent(in) :: aero_lut
       type(cirrus_table), intent(in) :: cirrus_lut
       type(settings_flags), intent(in) :: flag
       type(window_ini), dimension(:), intent(in) :: win_ini
       type(metadata), intent(in) :: meta
       type(aero), dimension(:), intent(in) :: aerosol_ini
       type(sun_spectrum), dimension(:), optional, intent(in) :: meas_irr_hi
      integer, intent(in) :: line_number  ! hack
       !*** Output
       type(retrieval_data), intent(out) :: retrieval_output
       integer, intent(out) :: ierr
       !*** Local variables
       type(instrument_response), dimension(size(win_ini)) :: response
       type(absorbers) :: absorb
       integer :: nwin
       real(double) ::  z_tropopause, z_bl
       real(double), dimension(:), allocatable :: dvair
       type(window_spectrum), dimension(size(win_ini)) :: win
       type(aero), dimension(:), allocatable :: aerosol
       type(atmosphere) :: atm_rt, atm_xs, atm_retr
       !--------------------------------------------------------------------------------
       ierr = 0
       nwin = size(win_ini)

       !*** Interface instrument input
       if (present(meas_irr_hi)) then
          call instrument_interface(&
               measurement, &
               response_in, &
               response, &
               flag%ils, &
               win_ini, &
               flag%output, &
               win, &
               ierr, &
               meas_irr_hi)
       else
          call instrument_interface(&
               measurement, &
               response_in, &
               response, &
               flag%ils, &
               win_ini, &
               flag%output, &
               win, &
               ierr)
       endif

       if (ierr .ne. 0) return

       !*** Interpolate meteo data to radiative transfer grid and retrieval grid
       call atmosphere_interpolate( &
            atm_scenario, &
            grid, &
            meta%surface_elevation, &
            meta%lat(1), &
            win_ini, &
            flag%output, &
            win, &
            atm_xs, &
            atm_rt, &
            atm_retr, &
            dvair, &
            z_tropopause, &
            z_bl, &
            ierr)
       if (ierr .ne. 0) return

       !*** Put absorber in type
       call get_absorbers(win_ini, absorb, ierr)
       if (ierr .ne. 0) return

       !*** Calculate aerosol initial guess
       call get_aerosol_properties_lognormal( &
            aero_lut, cirrus_lut, flag%output, meta%oceanglint, flag%glintscat, nwin, z_tropopause, z_bl, atm_rt, aerosol_ini, aerosol, ierr)
       if (ierr .ne. 0) return

       !*** Retrieval
       call profile_inversion( &
            aero_lut, cirrus_lut, &
            flag, meta%oceanglint, &
            meta%Fs_ini, &
            atm_scenario%surface_wspeed, &
            grid%nlay, &
            absorb, atm_rt, atm_xs, dvair, &
            response, &
            win_ini, win, aerosol, &
            retrieval_output, ierr, line_number)

       if (retrieval_output%error_id.ne.0 .or. ierr.ne.0 ) ierr = ierr_conv

       if (.not. allocated(retrieval_output%p)) then
          allocate(retrieval_output%p(grid%nlay+1),retrieval_output%z(grid%nlay+1),retrieval_output%t(grid%nlay+1), stat=ierr)
          if (ierr .ne. 0) then
             ierr = ierr_all
             call stopretrieval('RETRIEVAL: memory allocation problem')
             return
          endif
       endif
       retrieval_output%p = atm_retr%p
       retrieval_output%z = atm_retr%z
       retrieval_output%t = atm_retr%t

     end subroutine retrieval

!------------------------------------------------------------------------------
!> @brief
!> @details Put absorber information from win_in to type abs
!------------------------------------------------------------------------------
   subroutine get_absorbers(win_ini, absorb, ierr)
     !*** input
     type(window_ini), dimension(:), intent(in) :: win_ini
     !*** output
     type(absorbers), intent(out) :: absorb
     integer, intent(out) :: ierr
     !*** Local variables
     integer :: i, k, count1, off, j, nwin, n
     integer, dimension(:), allocatable :: temp1, temp2
     real(double) :: a, b
     character(199) :: message

     !*** Initialize
     ierr = 0
     nwin = size (win_ini)

     !*** Find redundant absorber types, rearrange, order

     !*** Count number of absorber types from all windows
     count1 = 0
     do n = 1, nwin
        count1 = count1 + win_ini(n)%ntype
     enddo

     !*** ALLOCATE dummy arrays with size count1
     allocate(temp1(count1),temp2(count1), stat=ierr)
     if (ierr .ne. 0) then
        write(message, *) 'GET_ABSORBERS: memory allocation problem'
        ierr = ierr_all
        goto 999
     endif

     !*** Write absorber type from all windows in one array
     !*** The indices from all interfering absorbers are written to array temp1
     !*** The indices from all target absorbers are written to array temp2
     temp1 = 0
     temp2 = 0
     off = 0
     do n = 1, nwin
        do k = 1, win_ini(n)%ntype
           if(win_ini(n)%type_x_flag(k)==0) temp1(k+off) = win_ini(n)%type_x(k)
           if(win_ini(n)%type_x_flag(k)==1) temp2(k+off) = win_ini(n)%type_x(k)
        enddo
        off = off + win_ini(n)%ntype
     enddo

     !*** Identify and zero all redundant absorber types
     !*** It may happen that an absorber occurs in different retrieval windows
     !*** If so, count only once (if more often occurs, set to 0)
     do j = 1, count1
        a = temp1(j)
        b = temp2(j)
        do i = j+1, count1
           if(temp1(i)==a) then
              temp1(i) = 0
           endif
           if(temp2(i)==b) then
              temp2(i) = 0
           endif
        enddo
     enddo

     !*** Count number of interfering (ntypes_global) and target (ntype_target) absorber types
     absorb%ntype_global = size(pack(temp1, temp1.gt.0))
     absorb%ntype_target = size(pack(temp2, temp2.gt.0))
     if(allocated(absorb%type_x_global)) deallocate(absorb%type_x_global)
     if(allocated(absorb%type_x_target)) deallocate(absorb%type_x_target)
     allocate(absorb%type_x_global(absorb%ntype_global), &
          absorb%type_x_target(absorb%ntype_target), &
          stat=ierr)
     if (ierr .ne. 0) then
        write(message, *) 'GET_ABSORBERS: memory allocation problem'
        ierr = ierr_all
        goto 999
     endif
     absorb%type_x_global = pack(temp1,temp1.gt.0)
     absorb%type_x_target = pack(temp2,temp2.gt.0)
     do i = 1, absorb%ntype_global
        do j = 1, absorb%ntype_target
           if(absorb%type_x_global(i)==absorb%type_x_target(j)) then
              ierr = ierr_var
              write(message, *) 'GET_ABSORBERS: Error in retrieval settings: ambiguous definition of target absorbers'
              goto 999
           endif
        enddo
     enddo

     return

999  continue
     call stopretrieval(message)

   end subroutine get_absorbers


!------------------------------------------------------------------------------
end module retrieval_module
