!------------------------------------------------------------------------------
!> @brief Define type windwow_ini and read retrieval.ini
!> @todo The reading of the initialization file should be handled by the C++
!! interface, while definition of type win_ini should be done by Fortran code.
!! We might need to change the format of initialization file.
!------------------------------------------------------------------------------
module read_errors_module
   use header_module
   implicit none
   private

!*** Instrument errors used or sensitivity analysis
   type, public :: instrument_errors
      real(double) :: sun_shift                    !< wavelength shift in Sun spectrum
      real(double) :: earth_shift                  !< wavelength shift in Earth spectrum
      real(double) :: earth_squeeze1               !< wavelength 1st order squeeze in Earth spectrum
      real(double) :: rad_offset                   !< radiometric offset
      real(double) :: rad_gain                     !< radiometric gain
      real(double) :: isrf                         !< instrument response function
   end type instrument_errors

!------------------------------------------------------------------------------

   public :: read_errors

contains

!------------------------------------------------------------------------------
   subroutine read_errors(runpath, runid, nwin, meteo_errors, instr_errors)
      character(stringlen), intent(in) :: runpath
      integer, intent(in) :: runid, nwin
!*** output
      real(double), dimension(4), intent(out) :: meteo_errors
      type(instrument_errors), dimension(nwin), intent(out) :: instr_errors
!*** local
      integer, parameter :: nwin_max = 10
      integer :: io, ierr, n
      character(6):: runidstring
      real(double) :: error_ch4, error_h2o, error_P, error_T
      real(double), dimension(nwin_max) ::  radiometric_offset, radiometric_gain, sun_spectral_shift, &
                                           earth_spectral_shift, earth_spectral_squeeze1, isrf_error
      namelist /meteo_data/ error_ch4, error_h2o, error_P, error_T
      namelist /instrument/ radiometric_offset, radiometric_gain, sun_spectral_shift, &
         earth_spectral_shift, earth_spectral_squeeze1, isrf_error

!*** Errors are zero if not specified in namelist file
      error_ch4 = 0.
      error_h2o = 0.
      error_P = 0.
      error_T = 0.
      radiometric_offset = 0.
      radiometric_gain = 0.
      sun_spectral_shift = 0.
      earth_spectral_shift = 0.
      earth_spectral_squeeze1 = 0.
      isrf_error = 0.

      write (runidstring, '(I6.6)') runid
!*** open namelist file
      open (newunit(io), file=trim(runpath)//'INI/syn_errors_'//runidstring//'.nml', status='OLD', iostat=ierr)
      if (ierr == 0) then
         !*** read namelists
         read (io, nml=instrument, iostat=ierr)
         write (6, nml=instrument)
         rewind (io)
         read (io, nml=meteo_data, iostat=ierr)
         write (6, nml=meteo_data)
         close (io)
      else
         call writelog('READ_ERRORS: error in reading errors.nml, continue without errors', 5)
      end if
      !*** put namelist variables in output variables
      meteo_errors(1) = error_ch4
      meteo_errors(2) = error_h2o
      meteo_errors(3) = error_P
      meteo_errors(4) = error_T
      do n = 1, nwin
         instr_errors(n)%sun_shift = sun_spectral_shift(n)
         instr_errors(n)%earth_shift = earth_spectral_shift(n)
         instr_errors(n)%earth_squeeze1 = earth_spectral_squeeze1(n)
         instr_errors(n)%rad_offset = radiometric_offset(n)
         instr_errors(n)%rad_gain = radiometric_gain(n)
         instr_errors(n)%isrf = isrf_error(n)
      end do

   end subroutine read_errors

!------------------------------------------------------------------------------
end module read_errors_module
