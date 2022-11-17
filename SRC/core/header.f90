!------------------------------------------------------------------------------
!> Data type definitions, constants and parameters
!> @todo move to relevant modules
!------------------------------------------------------------------------------

module header_module
  implicit none

  !------------------------------------------------------------------------------
  !*** Parameters
  !------------------------------------------------------------------------------
  !> Use intrinsic function to get numerical precision
  include 'INCLUDES/precision.inc'
  include 'INCLUDES/parameter.inc'

  integer, parameter :: stringlen = 299
  integer, parameter :: maxstr    = 16
  integer, parameter :: maxleg    = maxstr+1
  integer, parameter :: mxhalf    = maxstr/2
  integer, parameter :: npar_mie  = 6
  integer, parameter :: npar      = 8
  integer, parameter :: dim_x     = 500
  integer, parameter :: nfull     = dim_x+2
  integer, parameter :: ndang     = maxleg
  !> maximum optical thickness of the atmosphere for multiple scattering calculations.
  !! It is assumed that no light emerges from larger ot values, i.e. atmosphere is
  !! truncated at this value
  real(double), parameter :: tatot = 15.d0
  !> aerosol relevance threshold
  real(double), parameter :: nder_cut = 1.d-5
  real(double), parameter :: inf = 9.d99, null= 0.d0

  !------------------------------------------------------------------------------
  !*** constants
  !------------------------------------------------------------------------------
  real(double), parameter :: c_light = 299792458d0       ! speed of light in vacuum / m s-1
  real(double), parameter :: h_planck = 6.62607015d-34   ! Planck's constant / m2 kg s-1
  real(double), parameter :: air_m =  28.97d-3             ! air mass
  real(double), parameter :: rg = 8.3144621d0              ! gas constant
  real(double), parameter :: grav = 9.80665d0              ! gravitational acceleration [m/s^2]
  real(double), parameter :: avoga = 6.02214179d23         ! Avogadro's constant
  real(double), parameter :: pi = 3.14159265358979323846264338327950288419716939937510d0
  real(double), parameter :: relo2 = 0.2095d0              !used to be single precision, why?

  !------------------------------------------------------------------------------
  !*** default fill values
  !------------------------------------------------------------------------------

  integer, parameter :: nf_fill_byte = -127
  integer, parameter :: nf_fill_int1 = nf_fill_byte
  integer, parameter :: nf_fill_char = 0
  integer, parameter :: nf_fill_short = -32767
  integer, parameter :: nf_fill_int2 = nf_fill_short
  integer, parameter :: nf_fill_int = -2147483647
  real(single), parameter :: nf_fill_float = 9.9692099683868690e+36
  real(single),  parameter :: nf_fill_real = nf_fill_float
  real(double), parameter :: nf_fill_double = 9.9692099683868690d+36
  integer,  parameter :: nf_fill_ubyte = 255
  integer,  parameter :: nf_fill_ushort = 65535

  !------------------------------------------------------------------------------
  !*** error identifiers (TBC)
  !------------------------------------------------------------------------------
  integer, parameter :: ierr_open = 1    ! Error opening/closing a file           stop retrieval
  integer, parameter :: ierr_read = 2    ! Error reading from file                stop retrieval
  integer, parameter :: ierr_write = 3   ! Error writing to file                  stop retrieval
  integer, parameter :: ierr_var = 4     ! Error in values of read-in variables   stop retrieval
  integer, parameter :: ierr_all = 5     ! Error allocating memory                stop retrieval
  integer, parameter :: ierr_deall = 6   ! Error deallocating memory              stop retrieval
  integer, parameter :: ierr_isrf = 7    ! Error in ISRF                          stop retrieval
  integer, parameter :: ierr_filter = -1 ! Stopped by filter (clouds, sza, ..)    go to next pixel
  integer, parameter :: ierr_intrpl = -2 ! Interpolation error                    go to next pixel
  integer, parameter :: ierr_conv = -3   ! Convergence error in retrieval         go to next pixel
  integer, parameter :: ierr_l1b = -4    ! Error with L1B data                    go to next pixel
  integer, parameter :: ierr_meteo = -5  ! Error with meteo data                  go to next pixel
  integer, parameter :: ierr_apriori = -6! Error with apriori data                go to next pixel
  !------------------------------------------------------------------------------
  !> @brief atmosphere
  !------------------------------------------------------------------------------
  type :: atmosphere
     real(double), dimension(:), allocatable :: z
     real(double), dimension(:), allocatable :: dz
     real(double), dimension(:), allocatable :: p
     real(double), dimension(:), allocatable :: t
     integer :: n
  end type atmosphere

  !------------------------------------------------------------------------------
  !> @brief metadata
  !------------------------------------------------------------------------------
  type :: metadata
     !> fluorescence intensity
     real(double) :: Fs
     real(double) :: Fs_ini

     !> surface elevation of ground pixel [m]
     real(double) :: surface_elevation

     !> standard deviation of surface elevation within ground pixel [m]
     real(double) :: surface_elevation_stdv

     !> flags
     integer :: landflag, glintflag, oceanglint

     !> latitude and longitude of pixel centers and corners [degree]
     real(double) :: lat(5), lon(5)

     !> Date: year, month, day, hour, minute, second, millisecond
     integer :: time(7)

     !> Orbit altitude
     real(double) :: altitude  ! Added by JS

     !> Viewing geometry
     integer :: szaflag  ! Added by JS
     real(double) :: sza, iza, saz, iaz, phi
     real(double) :: observer_height

  end type metadata

  type :: sun_spectrum
     integer :: nwave
     real(double), dimension(:), allocatable :: wavelength
     real(double), dimension(:), allocatable :: irradiance
     integer, dimension(:), allocatable :: pixelflag                  !< pixelmask: 0=ok, 1=bad
  end type sun_spectrum

contains
  !----------------------------------------------------------------------------
  !> subroutine for error handling
  !----------------------------------------------------------------------------
  subroutine stopretrieval(msg)
    character(len=*),intent(in) :: msg

    call writelog(msg, 8);

  end subroutine stopretrieval
  !----------------------------------------------------------------------------
  !> subroutine for screen output (warnings, information, etc.)
  !----------------------------------------------------------------------------
  subroutine writelog(msg, level)
    character(len=*),intent(in) :: msg
    integer, intent(in) :: level
    ! level 1 tot 8:
    ! 1=trace, 2=debug, 3=information, 4=notice, 5=warning, 6=error, 7=critical, 8=fatal

    call fortranlog(msg, len(msg), level)
  end subroutine writelog

  !----------------------------------------------------------------------------
  !> subroutine setting message logger name
  !----------------------------------------------------------------------------
  subroutine writeloggername(loggername)
    character(len=*),intent(in) :: loggername

    call fortranloggername(loggername, len(loggername))
  end subroutine writeloggername

  !----------------------------------------------------------------------------
  integer function newunit(unit) result(n)
    ! returns lowest i/o unit number not in use
    integer, intent(out), optional :: unit
    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent

    do n = nmin, nmax
       inquire(unit=n, opened=inuse)
       if (.not. inuse) then
          if (present(unit)) unit=n
          return
       end if
    end do
    call stopretrieval("NEWUNIT: available unit not found.")
  end function newunit
  !----------------------------------------------------------------------------

  subroutine fortranlog(msg,msglen,level)
    implicit none
    character(*), intent(in) :: msg
    integer, intent(in) :: msglen
    integer, intent(in) :: level

    write(*,*) trim(msg)
  endsubroutine fortranlog

  !------------------------------------------------------------------------------
  subroutine fortranloggername(loggername,namelen)
    implicit none
    character(*), intent(in) :: loggername
    integer, intent(in) :: namelen

  endsubroutine fortranloggername
  !----------------------

end module header_module
