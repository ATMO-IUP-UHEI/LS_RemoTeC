module spectrum_internal_module
  use header_module
  use read_settings_module, only: window_ini
  use spectral_response_module, only: instrument_response, response_internal, spectral_response_stored, spectral_response_custom
  use forward_model_module, only: window_spectrum
  implicit none
  private

  !*** types
  public :: spectrum, instrument_response

  !*** Procedures
  public :: instrument_interface

  !------------------------------------------------------------------------------
  !> Spectral measurement
  type :: spectrum
     real(double) :: sza, iza, phi                               !< Instrument geometry                                 [degree]
     real(double) :: observer_height                              !< Instrument geometry                                [m]
!     real(double) :: fwhm                                        !< Width of instrument function                        [nm]
     !real(double), dimension(:), allocatable :: wavelength	  !rrae: OCO2 calibration in wavelength
     real(double), dimension(:), allocatable :: wavelength       !< wavelength grid                                     [nm]
     real(double), dimension(:), allocatable :: radiance         !< spectral radiance for each spectral pixel           [photons/s/cm^2/nm/sr]
     real(double), dimension(:), allocatable :: radiance_noise   !< statistical error of radiance                       [photons/s/cm^2/nm/sr]
     real(double), dimension(:), allocatable :: radiance_error   !< systematic error of radiance (or simulated noise)   [photons/s/cm^2/nm/sr]
     real(double), dimension(:), allocatable :: irradiance       !< spectral irradiance for each spectral pixel         [photons/s/cm^2/nm]
     real(double), dimension(:), allocatable :: irradiance_noise !< statistical error of irradiance                     [photons/s/cm^2/nm]
     real(double), dimension(:), allocatable :: irradiance_error !< systematic error of irradiance (or simulated noise) [photons/s/cm^2/nm]
     integer, dimension(:), allocatable :: mask                  !< pixel mask: 0=ok, 1=bad

!*** Multiplicative weighting factors for stokes coefficients provided by NASA L1b
     real(double), dimension(:), allocatable :: measurement_stokesc
     integer :: nwave                                            !< number of wavelengths
  end type spectrum
  !------------------------------------------------------------------------------

contains
  !------------------------------------------------------------------------------
  !> The parts of the measurements that fall within the retrieval windows
  !! are put in arrays (win(n)%spectrum, win(n)%spectrum_cov, win(n)%sun_spectrum_sat_lo)
  !! The instrument response is interpolated on internal wavelength grid.
  !! Convolve model solar spectrum by instrument response function.
  !------------------------------------------------------------------------------
  subroutine instrument_interface(&
       measurement, &
       response_in, &
       response, &
       ilsflag, &
       win_ini, &
       outputflag, &
       win, &
       ierr, &
       meas_irr_hi)
    !*** input variables
    type(spectrum), dimension(:), intent(in) :: measurement
    type(instrument_response), dimension(:), intent(in) :: response_in
    type(window_ini), dimension(:), intent(in) :: win_ini
    integer, intent(in) :: ilsflag, outputflag
    type(sun_spectrum), dimension(:), optional, intent(in) :: meas_irr_hi
    !*** Output
    type(window_spectrum), dimension(:), intent(out) :: win
    type(instrument_response), dimension(:), intent(out) :: response
    integer, intent(out) :: ierr
    !*** local variables
    logical :: sun_meas
    integer :: i, l, k, n, nwin, nils
    character(199) :: message
    integer :: nstokes_l1b
    !------------------------------------------------------------------------------
    if(outputflag >= 2)then
       write(message,'(a)') '*** Start instrument_interface ***'
       call writelog(message, 1)
    endif


    nstokes_l1b = 1
    nwin = size(win)
    !*** Put wavelength independent data into datatype win
    do n = 1, nwin

       do i = 1, size(measurement)
       if (allocated(measurement(i)%measurement_stokesc)) then
          nstokes_l1b = size(measurement(i)%measurement_stokesc)
          if(.not. allocated(win(n)%measurement_stokesc)) allocate(win(n)%measurement_stokesc(nstokes_l1b))
       endif
          do k = 1, measurement(i)%nwave
             if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                  measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop)then
!                response(n)%nils = response_in(i)%nils
                response(n)%nils = size(response_in(i)%resp_store(1,:))
                win(n)%sza = measurement(i)%sza
                win(n)%iza = measurement(i)%iza
                win(n)%phi = measurement(i)%phi
                win(n)%observer_height = measurement(i)%observer_height
                if (nstokes_l1b>1) win(n)%measurement_stokesc = measurement(i)%measurement_stokesc
                goto 100
             endif
          enddo
      enddo
100   continue
    enddo

    !*** Put wavelength dependent data into datatype win
    sun_meas = .false.
    nwin = size(win)
    do n = 1, nwin
       win(n)%nwave_lo = 0
       do i = 1, size(measurement)
          if (allocated(measurement(i)%mask)) then
             do k = 1, measurement(i)%nwave
                if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                     measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop .and. measurement(i)%mask(k)==0)then
                   win(n)%nwave_lo = win(n)%nwave_lo + 1
                endif
             enddo
          else
             do k = 1, measurement(i)%nwave
                if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                     measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop)then
                   win(n)%nwave_lo = win(n)%nwave_lo + 1
                endif
             enddo
          endif
       enddo
       response(n)%nwave = win(n)%nwave_lo
    enddo

    do n = 1, nwin
       if(.not. allocated(win(n)%sun_spectrum_sat_hi)) allocate(win(n)%sun_spectrum_sat_hi(win_ini(n)%nwave_hi))
       !*** Only put data between wave_start and wave_stop in win(n)%spectrum
       if(.not. allocated(win(n)%spectrum)) then
          allocate(win(n)%spectrum(win(n)%nwave_lo), &
               win(n)%spectrum_cov(win(n)%nwave_lo), &
               win(n)%sun_spectrum_sat_lo(win(n)%nwave_lo), &
               win(n)%wavelength_lo(win(n)%nwave_lo), &
               response(n)%ils_dwave(win(n)%nwave_lo, response(n)%nils), &
               response(n)%resp_store(win(n)%nwave_lo, response(n)%nils), &
               stat=ierr)
          if (ierr .ne. 0) then
             write(message,*) 'INSTRUMENT_INTERFACE: memory allocation error'
             ierr = ierr_all
             goto 999
          endif

       elseif (size(win(n)%spectrum) .ne. win(n)%nwave_lo) then
          deallocate(win(n)%spectrum, &
               win(n)%spectrum_cov, &
               win(n)%sun_spectrum_sat_lo, &
               win(n)%wavelength_lo, &
               response(n)%ils_dwave, &
               response(n)%resp_store, &
               stat=ierr)
          if (ierr .ne. 0) then
             write(message,*) 'INSTRUMENT_INTERFACE: memory deallocation error'
             ierr = ierr_deall
             goto 999
          endif
          allocate(win(n)%spectrum(win(n)%nwave_lo), &
               win(n)%spectrum_cov(win(n)%nwave_lo), &
               win(n)%sun_spectrum_sat_lo(win(n)%nwave_lo), &
               win(n)%wavelength_lo(win(n)%nwave_lo), &
               response(n)%ils_dwave(win(n)%nwave_lo, response(n)%nils), &
               response(n)%resp_store(win(n)%nwave_lo, response(n)%nils), &
               stat=ierr)
          if (ierr .ne. 0) then
             write(message,*) 'INSTRUMENT_INTERFACE: memory allocation error'
             ierr = ierr_all
             goto 999
          endif
       endif

       l = 0
       do i = 1, size(measurement)

          if (allocated(measurement(i)%mask)) then  ! If pixelmask is defined, skip bad spectral pixels (with mask .ne. 0)
             if (allocated(measurement(i)%irradiance)) then  ! measured irradiance
                sun_meas = .true.
                do k = 1, measurement(i)%nwave
                   if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                        measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop .and. &
                        measurement(i)%mask(k)==0)then
                      win(n)%spectrum(l+1) = measurement(i)%radiance(k)
                      win(n)%spectrum_cov(l+1) = measurement(i)%radiance_noise(k)**2
                      response(n)%ils_dwave(l+1, :) = response_in(i)%ils_dwave(k, :)
                      response(n)%resp_store(l+1, :) = response_in(i)%resp_store(k, :)
                      win(n)%wavelength_lo(l+1) = measurement(i)%wavelength(k)
                      win(n)%sun_spectrum_sat_lo(l+1) = measurement(i)%irradiance(k)
                      l = l + 1
                   endif
                enddo
             else                                             ! no measured irradiance
                do k = 1, measurement(i)%nwave
                   if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                        measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop .and. &
                        measurement(i)%mask(k)==0)then
                      win(n)%spectrum(l+1) = measurement(i)%radiance(k)
                      win(n)%spectrum_cov(l+1) = measurement(i)%radiance_noise(k)**2
                      response(n)%ils_dwave(l+1, :) = response_in(i)%ils_dwave(k, :)
                      response(n)%resp_store(l+1, :) = response_in(i)%resp_store(k, :)
                      win(n)%wavelength_lo(l+1) = measurement(i)%wavelength(k)
                      l = l + 1
                   endif
                enddo
             endif
          else ! If pixelmask is not defined, use all spectral pixels
             if (allocated(measurement(i)%irradiance)) then ! measured irradiance
                sun_meas = .true.
                do k = 1, measurement(i)%nwave
                   if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                        measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop)then
                      win(n)%spectrum(l+1) = measurement(i)%radiance(k)
                      win(n)%spectrum_cov(l+1) = measurement(i)%radiance_noise(k)**2
                      response(n)%ils_dwave(l+1, :) = response_in(i)%ils_dwave(k, :)
                      response(n)%resp_store(l+1, :) = response_in(i)%resp_store(k, :)
                      win(n)%wavelength_lo(l+1) = measurement(i)%wavelength(k)
                      win(n)%sun_spectrum_sat_lo(l+1) = measurement(i)%irradiance(k)
                      l = l + 1
                   endif
                enddo
             else                                           ! no measured irradiance
                do k = 1, measurement(i)%nwave
                   if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                        measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop)then
                      win(n)%spectrum(l+1) = measurement(i)%radiance(k)
                      win(n)%spectrum_cov(l+1) = measurement(i)%radiance_noise(k)**2
                      response(n)%ils_dwave(l+1, :) = response_in(i)%ils_dwave(k, :)
                      response(n)%resp_store(l+1, :) = response_in(i)%resp_store(k, :)
                      win(n)%wavelength_lo(l+1) = measurement(i)%wavelength(k)
                      l = l + 1
                   endif
                enddo
             endif
          endif
       enddo ! loop over i
       response(n)%nils = size(response(n)%resp_store(1, :))
       response(n)%nwave = size(response(n)%resp_store(:, 1))

       nils = response(n)%nils
       !*** Store instrument response in data type "win"
       if(.not. allocated(win(n)%resp_store)) then
          allocate(win(n)%resp_store(win(n)%nwave_lo,nils), &
               win(n)%is_store(win(n)%nwave_lo), &
               win(n)%ie_store(win(n)%nwave_lo), &
               win(n)%sun_spectrum_ref_lo(win(n)%nwave_lo), &
               stat=ierr)
          if (ierr .ne. 0) then
             write(message,*) 'INSTRUMENT_INTERFACE: memory allocation error'
             ierr = ierr_all
             goto 999
          endif
       elseif(size(win(n)%is_store) .ne. win(n)%nwave_lo) then
          deallocate(win(n)%resp_store, &
               win(n)%is_store, &
               win(n)%ie_store, &
               win(n)%sun_spectrum_ref_lo, &
               stat = ierr)
          if (ierr .ne. 0) then
             write(message,*) 'INSTRUMENT_INTERFACE: memory deallocation error'
             ierr = ierr_deall
             goto 999
          endif
          allocate(win(n)%resp_store(win(n)%nwave_lo,nils), &
               win(n)%is_store(win(n)%nwave_lo), &
               win(n)%ie_store(win(n)%nwave_lo), &
               win(n)%sun_spectrum_ref_lo(win(n)%nwave_lo), &
               stat=ierr)
          if (ierr .ne. 0) then
             write(message,*) 'INSTRUMENT_INTERFACE: memory allocation error'
             ierr = ierr_all
             goto 999
          endif
       elseif(size(win(n)%resp_store(1,:)) .ne. nils) then
          deallocate(win(n)%resp_store, stat=ierr)
          if (ierr .ne. 0) then
             write(message,*) 'INSTRUMENT_INTERFACE: memory deallocation error'
             ierr = ierr_deall
             goto 999
          endif
          allocate(win(n)%resp_store(win(n)%nwave_lo,nils), stat=ierr)
          if (ierr .ne. 0) then
             write(message,*) 'INSTRUMENT_INTERFACE: memory allocation error'
             ierr = ierr_all
             goto 999
          endif

       endif

       !*** Get hi-reso solar spectrum from deconvoluted measurement if present, otherwise from model
       if (present(meas_irr_hi)) then
          win(n)%sun_spectrum_sat_hi(:) = meas_irr_hi(n)%irradiance(:)
       else
          win(n)%sun_spectrum_sat_hi(:) = win_ini(n)%sun_spectrum_ref_hi(:)
       endif

       !*** Get convoluted (lo-reso) model solar spectrum
       if (ilsflag == 1) then
          !*** Convert ILS from measurement wavelength grid to internal grid
          call response_internal( &
               win_ini(n)%wavelength_hi, &   ! high-resolution wavelength grid of model
               win(n)%wavelength_lo, &       ! low-resolution wavelength grid of measurement
               response(n)%ils_dwave, &      ! delta(wavelength) grid of instrument line shape
               response(n)%resp_store, &     ! instrument line function
               win(n)%resp_store, &          ! array for storing instrument response
               win(n)%ie_store, &            ! array for storing instrument response
               win(n)%is_store, &            ! array for storing instrument response
               ierr)                         ! error identifier
          if (ierr .ne. 0) return
          !*** Convolve solar spectrum by instrument response function:
          call spectral_response_stored( &
               win(n)%resp_store, &
               win(n)%ie_store, &
               win(n)%is_store, &
               win(n)%sun_spectrum_sat_hi, &
               win(n)%sun_spectrum_ref_lo)
       elseif (ilsflag == 2) then
          call spectral_response_custom(win_ini(n)%wavelength_hi, win(n)%sun_spectrum_sat_hi, win_ini(n)%nwave_hi,&
               win(n)%wavelength_lo, win(n)%sun_spectrum_ref_lo, win(n)%nwave_lo,&
               response(n)%ils_dwave(1,:), response(n)%resp_store(1,:), nils, ierr)
          if (ierr .ne. 0) return
       endif

       !*** Take convoluted solar spectrum if there is no measurement
       !*** or if deconvolution of measurement has been performed
       if (.not. sun_meas .or. present(meas_irr_hi)) then !no measured solar spectrum
          win(n)%sun_spectrum_sat_lo(:) = win(n)%sun_spectrum_ref_lo(:)
       endif

    enddo !loop over n

    if (outputflag >=2) then
       write(message,'(a)') '*** End of instrument_interface ***'
       call writelog(message, 1)
    endif

    return

999 continue
    call stopretrieval(message)

  end subroutine instrument_interface
  !------------------------------------------------------------------------------

end module spectrum_internal_module
