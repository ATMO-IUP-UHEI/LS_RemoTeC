module read_miprep_module
  use header_module
  use atmosphere_internal_module, only: atmospheric_scenario
  use forward_model_module, only: aero, window_spectrum
!  use auxiliary_routines_module
  use netcdf
  implicit none
  private

  !*** procedures
  public :: read_miprep, read_meteo, read_co, read_ch4, read_dem, read_orbit, open_miprep, close_miprep


  integer, private :: ncid_orbit, ncid_ecmwf, ncid_ch4, ncid_co, ncid_dem, ncid_co2
  logical, private :: file_co, file_co2, file_ch4  
     
contains
  
  !------------------------------------------------------------------------------
  
  subroutine open_miprep(infile, ierr)
    character(len=*), intent(in) :: infile
    integer, intent(out) :: ierr
    character(stringlen) :: message
    
    call check(nf90_open(trim(infile)//'_orbit.nc', nf90_nowrite, ncid_orbit), ierr)
    if (ierr .ne. 0) then
       write(message,*) 'OPEN_MIPREP: error opening file', trim(infile)//'_orbit.nc', ierr
       ierr = ierr_open
       goto 999
    endif

    call check(nf90_open(trim(infile)//'_meteo.nc', nf90_nowrite, ncid_ecmwf), ierr)
    if (ierr .ne. 0) then
       write(message,*) 'OPEN_MIPREP: error opening file', trim(infile)//'_meteo.nc', ierr
       ierr = ierr_open
       goto 999
    endif
    
    call check(nf90_open(trim(infile)//'_dem.nc', nf90_nowrite, ncid_dem), ierr)
    if (ierr .ne. 0) then
       write(message,*) 'OPEN_MIPREP: error opening file', trim(infile)//'_dem.nc', ierr
       ierr = ierr_open
       goto 999
    endif
    
    call check(nf90_open(trim(infile)//'_ch4.nc', nf90_nowrite, ncid_ch4), ierr)
    if (ierr == 0) then
       file_ch4 = .true.
    else
       file_ch4 = .false.       
    endif
    
    call check(nf90_open(trim(infile)//'_co.nc', nf90_nowrite, ncid_co), ierr)
    if (ierr == 0) then
       file_co = .true.
    else
       file_co = .false.       
    endif

    call check(nf90_open(trim(infile)//'_co2.nc', nf90_nowrite, ncid_co2), ierr)
    if (ierr == 0) then
       file_co2 = .true.
    else
       file_co2 = .false.       
    endif

    
!!$    call check(nf90_open(trim(infile)//'_aerosol.nc', nf90_nowrite, ncid_aer), ierr)
!!$    if (ierr .ne. 0) then
!!$       ierr = ierr_open
!!$       write(message,*) 'READ_MIPREP: error opening file', trim(infile)//'_aerosol.nc'
!!$       goto 999
!!$    endif
!!$
!!$    call check(nf90_open(trim(infile)//'_cloud.nc', nf90_nowrite, ncid_cloud), ierr)
!!$    if (ierr .ne. 0) then
!!$       ierr = ierr_open
!!$       write(message,*) 'READ_MIPREP: error opening file', trim(infile)//'_cloud.nc'
!!$       goto 999
!!$    endif
!!$    
    ierr = 0 ! Mark succes
    return

999 continue
    call stopretrieval(message)
    return

  end subroutine open_miprep
  
  !------------------------------------------------------------------------------
  
  subroutine close_miprep(infile, ierr)
    character(len=*), intent(in) :: infile   
    integer, intent(out) :: ierr
    character(stringlen) :: message

    !*** All information has been read, close NetCDF file. 
    call check(NF90_CLOSE (ncid_orbit), ierr)
    if (ierr .ne. 0) then
       write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_orbit.nc'
       ierr = ierr_open
       goto 999
    endif

    !*** All information has been read, close NetCDF file. 
    call check(NF90_CLOSE (ncid_ecmwf), ierr)
    if (ierr .ne. 0) then
       write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_meteo.nc'
       ierr = ierr_open
       goto 999
    endif

    !*** All information has been read, close NetCDF file. 
    call check(NF90_CLOSE (ncid_dem), ierr)
    if (ierr .ne. 0) then
       write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_dem.nc'
       ierr = ierr_open
       goto 999
    endif

    if (file_ch4) then
    !*** All information has been read, close NetCDF file. 
    call check(NF90_CLOSE (ncid_ch4), ierr)
    if (ierr .ne. 0) then
       write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_ch4.nc'
       ierr = ierr_open
       goto 999
    endif
    endif

    !*** All information has been read, close NetCDF file.
    if (file_co) then
       call check(NF90_CLOSE (ncid_co), ierr)
       if (ierr .ne. 0) then
          write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_co.nc'
          ierr = ierr_open
          goto 999
       endif
    endif

    if (file_co2) then
    !*** All information has been read, close NetCDF file. 
    call check(NF90_CLOSE (ncid_co2), ierr)
    if (ierr .ne. 0) then
       write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_co2.nc'
       ierr = ierr_open
       goto 999
    endif
    endif

!!$    !*** All information has been read, close NetCDF file. 
!!$    call check(NF90_CLOSE (ncid_aer), ierr)
!!$    if (ierr .ne. 0) then
!!$       write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_aerosol.nc'
!!$       ierr = ierr_open
!!$       goto 999
!!$    endif
!!$
!!$    !*** All information has been read, close NetCDF file. 
!!$    call check(NF90_CLOSE (ncid_cloud), ierr)
!!$    if (ierr .ne. 0) then
!!$       write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_cloud.nc'
!!$       ierr = ierr_open
!!$       goto 999
!!$    endif

    ierr = 0 ! Mark succes
    return

999 continue
    call stopretrieval(message)
    return

  endsubroutine  close_miprep
 

  !------------------------------------------------------------------------------ 
  subroutine read_miprep(&
       infile, ipixel,  &
       heightlay_in, presslay_in, templay_in, &
       h2olay_in, co2lay_in, ch4lay_in, colay_in, &
       meta, ecmwf_psurf, ecmwf_z, surface_wspeed, ierr, aerosolfile, aerosol_ini)  
    !*** input   
    character(len=*), intent(in) :: infile
    integer, intent(in) :: ipixel
    !*** output
    character(stringlen), optional, intent(out) :: aerosolfile
    type(aero), dimension(:), optional, intent(inout) :: aerosol_ini
    real(double), dimension(:),allocatable, intent(out) :: presslay_in  ! Pressure at the pressure-center of the layers [hPa]
    real(double), dimension(:),allocatable, intent(out) :: templay_in   ! Temperature at the pressure-center of the layers [K]
    real(double), dimension(:),allocatable, intent(out) :: h2olay_in    ! VMR H2O
    real(double), dimension(:),allocatable, intent(out) :: heightlay_in ! Height at pressure-center of the layers [m]
    real(double), dimension(:),allocatable, intent(out) :: co2lay_in    ! VMR CO2
    real(double), dimension(:),allocatable, intent(out) :: ch4lay_in    ! VMR CO2
    real(double), dimension(:),allocatable, intent(out) :: colay_in     ! VMR CO
    type(metadata), intent(out) :: meta
    real(double), intent(out) :: ecmwf_psurf, ecmwf_z, surface_wspeed  
    integer, intent(out) :: ierr
    !*** local
    type(aero) :: aerosol_tmp
    integer :: ninput, nlay, nstart, nchar, ntype_aer
    integer :: i, k, glintflag, ncid_aer, ncid_cloud
!!$    real(double), dimension(:), allocatable :: heightlev_in ! Height at layer boundaries [m]
!!$    real(double), dimension(:), allocatable :: presslev_in  ! Pressure at layer boundaries [hPa]
    real(double), dimension(:), allocatable :: temp_in   ! Temperature at layer boundaries [K]
!!$    real(double), dimension(:), allocatable :: h2olev_in    ! VMR H2O
!!$    real(double), dimension(:), allocatable :: co2lev_in    ! VMR CO2
!!$    real(double), dimension(:), allocatable :: ch4lev_in    ! VMR CH4
!!$    real(double), dimension(:), allocatable :: colev_in     ! VMR CO
    integer :: n_ch4, n_co, n_co2
    !***
    integer :: varid, dimid
    integer :: start(1), count(1), start_prof(2), count_prof(2) 
    real(double) ::  psurf, u10, v10, cth(1), cot(1), cf(1), dz
    real(double), dimension(:), allocatable :: hyai, hybi, q, hyam, hybm
    real(double), dimension(:), allocatable :: TM5_p, TM5_ch4, TM5_co, ct_p, ct_co2
    character(stringlen) :: message
    !------------------------------------------------------------------------------ 

    !*** Get start indices
    start = [ipixel]
    count = [1]
    start_prof = [1, ipixel]

    call read_orbit(infile, ipixel, meta%sza, meta%iza, meta%saz, meta%iaz, &
         meta%lat(1), meta%lon(1), meta%lat(2:5), meta%lon(2:5), meta%time, ierr)
    if (ierr .ne. 0) then
       write(message,*) 'READ_MIPREP: error reading file', trim(infile)//'_orbit.nc'
       goto 999
    endif


    !*** Get meteo data
    call read_meteo(infile, ipixel, ninput, hyam, hybm, q, temp_in, psurf, ecmwf_z, u10, v10, ierr)
    if (ierr .ne. 0) then
       write(message,*) 'READ_MIPREP: error reading file', trim(infile)//'_meteo.nc'
       ierr = ierr_meteo
       goto 999
    endif
    surface_wspeed = sqrt(u10**2.d0 + v10**2.d0)
    ecmwf_psurf = psurf/100.d0  !convert Pa to hPa

    allocate(presslay_in(ninput+1),h2olay_in(ninput+1))
    allocate(heightlay_in(ninput+1),templay_in(ninput+1), co2lay_in(ninput+1),ch4lay_in(ninput+1),colay_in(ninput+1))

    !*** Convert
    do k = 1, ninput
       !*** Convert hybrid coefficients to pressure [from Pa to hPa]
       presslay_in(k) = (hyam(k) + hybm(k)*psurf)/100.  
       !*** Convert specific humidity to dry VMR of water
       !*** Specific humidity s  is mass of water per mass of HUMID air.
       !*** The mass mixing ratio s* is mass of water per mass of DRY air: s* = s / (1 - s).
       !*** The volume mixing ratio is given by: (s*) * Mair / Mwater = s / (1 - s) * Mair / Mwater.      
       h2olay_in(k) = q(k)/(1.0d0-q(k))*1.60855      
    enddo
    presslay_in(ninput+1) = ecmwf_psurf
    h2olay_in(ninput+1)= h2olay_in(ninput)
    dz = rg*temp_in(ninput)/(air_m*grav)*log(presslay_in(ninput)/ecmwf_psurf)
    templay_in(1:ninput) = temp_in(1:ninput)
    templay_in(ninput+1) = temp_in(ninput) + 0.0065*dz

!!$    !*** Remove upper layers with pressure < 0.01 mbar
!!$    do k = 1, ninput
!!$       if (presslay_in(k) > 1.d-2) then 
!!$          nstart = k
!!$          exit
!!$       endif
!!$    enddo
!!$    nlay = 1 + ninput - nstart 
!!$
!!$    allocate(heightlev_in(nlay+1),presslev_in(nlay+1),templev_in(nlay+1),&
!!$         h2olev_in(nlay+1),co2lev_in(nlay+1),ch4lev_in(nlay+1),colev_in(nlay+1))
!!$
!!$    !*** Not needed
!!$    co2lev_in = 0.d0
!!$
!!$    !** Convert from Pa to hPa
!!$    presslev_in(nlay+1) = psurf/100.d0
!!$
!!$    !*** Calculate pressure at layer boundaries 
!!$      do k = ninput, 1, -1
!!$         presslev_in(k) = 2*presslay_in(k) - presslev_in(k+1)!HH: This goes wrong for thin layers
!!$      enddo
!!$    do k = nlay, 2, -1
!!$       presslev_in(k) = (presslay_in(ninput-nlay+k) + presslay_in(ninput-nlay+k-1))/2.
!!$    enddo
!!$
!!$    !*** constrain pressure grid to height = 0.01 mbar (~75 km)
!!$    presslev_in(1) = 1.d-2
!!$
!!$    !*** Interpolate  to layer boundaries
!!$    !      call spline_interpol(DLOG(presslay_in), heightlay_in, ninput,&
!!$    !        		   DLOG(presslev_in), heightlev_in, ninput+1, ierr)
!!$    call spline_interpol(DLOG(presslay_in(nstart:ninput)), templay_in(nstart:ninput), nlay,&
!!$         DLOG(presslev_in), templev_in, nlay+1, ierr)
!!$    call linterp(presslay_in(nstart:ninput), h2olay_in(nstart:ninput), nlay,&
!!$         presslev_in, h2olev_in, nlay+1)

!!$    !*** Get height from barometric formula, i.e. ideal gas law
!!$    !*** from surface to top of atmopshere
!!$    heightlev_in(nlay+1) = ecmwf_z
!!$    heightlev_in(nlay) = ecmwf_z - &
!!$         rg*templev_in(nlay)/(air_m*grav)*log(presslev_in(nlay)/(psurf/100.))
!!$    do i = 1, nlay-1
!!$       heightlev_in(nlay -i)  = heightlev_in(nlay - i +1)  - &
!!$            rg*templev_in(nlay-i)/(air_m*grav)*log(presslev_in(nlay-i)/presslev_in(nlay-i+1))
!!$    enddo
!!$
    !midpoints 
    !*** Get height from barometric formula, i.e. ideal gas law
    !*** from surface to top of atmopshere
    heightlay_in(ninput) = ecmwf_z - &
         rg*templay_in(ninput)/(air_m*grav)*log(presslay_in(ninput)/(psurf/100.))
    do i = 1, ninput-1
       heightlay_in(ninput -i)  = heightlay_in(ninput - i +1)  - &
            rg*templay_in(ninput-i)/(air_m*grav)*log(presslay_in(ninput-i)/presslay_in(ninput-i+1))
    enddo
    heightlay_in(ninput+1) = ecmwf_z

    if (file_ch4) then
       !*** Get CH4 apriori profiles
       call read_ch4(infile, ipixel, n_ch4, hyai, hybi, tm5_ch4, psurf, ierr)
       if (ierr .ne. 0) then
          write(message,*) 'READ_MIPREP: error reading file', trim(infile)//'_ch4.nc'
          ierr = ierr_meteo
          goto 999
       endif
       allocate(tm5_p(n_ch4))

       !*** Convert
       do k = 1, n_ch4
          !*** Convert hybrid coefficients to pressure [from Pa to hPa]
          tm5_p(k) = ((hyai(k)+hyai(k+1))/2.d0 + (hybi(k)+hybi(k+1))/2d0*psurf)/100.     
       enddo
!!$
!!$    !*** Reverse array order (from TOA to surface)
!!$    call linterp(TM5_p(n_ch4:1:-1), TM5_ch4(n_ch4:1:-1), n_ch4,&
!!$         presslev_in, ch4lev_in, nlay+1)
!!$
       !*** Get values at layer midpoints:
       call linterp(TM5_p(n_ch4:1:-1), TM5_ch4(n_ch4:1:-1), n_ch4,&
            presslay_in, ch4lay_in, ninput+1)
    endif

    if (file_co) then
       !*** Get CO apriori profiles
       call read_co(infile, ipixel, n_co, hyai, hybi, tm5_co, psurf, ierr)
       if (ierr .ne. 0) then
          write(message,*) 'READ_MIPREP: error reading file', trim(infile)//'_co.nc'
          ierr = ierr_meteo
          goto 999
       endif

       deallocate(tm5_p)
       allocate(tm5_p(n_co))
       !*** Convert
       do k = 1, n_co
          !*** Convert hybrid coefficients to pressure [from Pa to hPa]
          tm5_p(k) = ((hyai(k)+hyai(k+1))/2.d0 + (hybi(k)+hybi(k+1))/2.d0*psurf)/100.     
       enddo

!!$    call linterp(TM5_p(n_co:1:-1), TM5_co(n_co:1:-1), n_co,&
!!$         presslev_in, colev_in, nlay+1)
!!$
!!$
       !*** Get values at layer midpoints:
       call linterp(TM5_p(n_co:1:-1), TM5_co(n_co:1:-1), n_co,&
            presslay_in, colay_in, ninput+1)

    endif

    if (file_co2) then
       !*** Get CO2 apriori profiles
       call read_co2(infile, ipixel, n_co2, ct_co2, ct_p, ierr)

       if (ierr .ne. 0) then
          write(message,*) 'READ_MIPREP: error reading file', trim(infile)//'_co2.nc'
          ierr = ierr_meteo
          goto 999
       endif
       !*** Get values at layer midpoints:
       call linterp(ct_p(n_co2:1:-1), ct_co2(n_co2:1:-1), n_co2,&
            presslay_in, co2lay_in, ninput+1)
    endif

    !*** Get digital elevation map
    call read_dem(infile, ipixel, meta%surface_elevation, meta%surface_elevation_stdv, meta%landflag, ierr)

    if (present(aerosolfile)) then
       !*** Read Aerosol model info
       !$OMP critical
       call check(nf90_open(trim(infile)//'_aerosol.nc', nf90_nowrite, ncid_aer), ierr)
       !$OMP end critical
       if (ierr .ne. 0) then
          ierr = ierr_open
          write(message,*) 'READ_MIPREP: error opening file', trim(infile)//'_aerosol.nc'
          goto 999
       endif

       call check(NF90_INQ_DIMID(ncid_aer, "nchar", dimid), ierr)
       if (ierr .ne. 0) return
       call check(NF90_INQUIRE_DIMENSION(ncid_aer, dimid, len = nchar), ierr)
       if (ierr .ne. 0) return


       count_prof = [nchar,1] 
       !*** Get aerosolfile
       aerosolfile = "filled"
       call check(NF90_INQ_VARID(ncid_aer, "aerosolfile", varid), ierr) 
       if (ierr .ne. 0) return
       call check(NF90_GET_VAR(ncid_aer, varid, aerosolfile, start=start_prof, count=count_prof), ierr)
       if (ierr .ne. 0) return

       !*** All information has been read, close NetCDF file. 
       call check(NF90_CLOSE (ncid_aer), ierr)
       if (ierr .ne. 0) then
          write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_aerosol.nc'
          ierr = ierr_open
          goto 999
       endif
    endif

    if (present(aerosol_ini)) then
       !*** Read cloud info
       ntype_aer = size(aerosol_ini)
       aerosol_tmp = aerosol_ini(ntype_aer)

       !$OMP critical
       call check(nf90_open(trim(infile)//'_cloud.nc', nf90_nowrite, ncid_cloud), ierr)
       !$OMP end critical
       if (ierr .ne. 0) then
          ierr = ierr_open
          write(message,*) 'READ_MIPREP: error opening file', trim(infile)//'_cloud.nc'
          goto 999
       endif

       call check(NF90_INQ_VARID(ncid_cloud, "cloud_top_height", varid), ierr) 
       if (ierr .ne. 0) return
       call check(NF90_GET_VAR(ncid_cloud, varid, cth, start=start, count=count), ierr)
       if (ierr .ne. 0) return
       aerosol_tmp%aeralt1 = cth(1) - 0.5d0* aerosol_tmp%aeralt2

       call check(NF90_INQ_VARID(ncid_cloud, "cloud_fraction", varid), ierr) 
       if (ierr .ne. 0) return
       call check(NF90_GET_VAR(ncid_cloud, varid, cf, start=start, count=count), ierr)
       if (ierr .ne. 0) return
       aerosol_tmp%shapefrac = cf(1)

       call check(NF90_INQ_VARID(ncid_cloud, "cloud_optical_thickness", varid), ierr) 
       if (ierr .ne. 0) return
       call check(NF90_GET_VAR(ncid_cloud, varid, cot, start=start, count=count), ierr)
       if (ierr .ne. 0) return
       aerosol_tmp%tau_ref = cot(1)

       !*** All information has been read, close NetCDF file. 
       call check(NF90_CLOSE (ncid_cloud), ierr)
       if (ierr .ne. 0) then
          write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_cloud.nc'
          ierr = ierr_open
          goto 999
       endif
       !*** Overwrite last aerosoltype with MODIS cloud info
       aerosol_ini(ntype_aer) = aerosol_tmp

    endif

    !*** Define sunglintflag
    if (abs(abs(meta%sza)-abs(meta%iza))<2. .and. (meta%saz-meta%iaz)<20.) then
       glintflag = 1
    else
       glintflag = 0
    endif
    !*** Set ocean glintflag here 
    if( meta%landflag==1 .and. glintflag==1) then
       meta%oceanglint = 1
    else
       meta%oceanglint = 0
    endif


    ierr = 0 ! Mark succes
    return

999 continue
    call stopretrieval(message)
    return


  end subroutine read_miprep

  !------------------------------------------------------------------------------ 
  subroutine read_meteo(infile, ipixel, ninput, hyam, hybm, q, t, psurf, ecmwf_z, u, v, ierr)
    character(len=*), intent(in) :: infile
    integer, intent(in) :: ipixel
    real(double), dimension(:), allocatable, intent(out)  :: hyam, hybm, q, t    
    real(double), intent(out) :: psurf, ecmwf_z, u, v
    integer, intent(out) :: ninput, ierr
    !***
    integer ::  varid, dimid
    integer :: start(1), count(1), start_prof(2), count_prof(2) 
    real(double) :: sp(1), z(1), u10(1), v10(1)
    character(stringlen) :: message

    !*** Get start indices
    start = [ipixel]
    count = [1]
    start_prof = [1, ipixel]

    !*** Read ECMWF
!!$    !$OMP critical
!!$    call check(nf90_open(trim(infile)//'_meteo.nc', nf90_nowrite, ncid_ecmwf), ierr)
!!$    !$OMP end critical
!!$    if (ierr .ne. 0) then
!!$       ierr = ierr_open
!!$       write(message,*) 'READ_MIPREP: error opening file', trim(infile)//'_meteo.nc'
!!$       goto 999
!!$    endif

    !*** Get number of levels
    call check(NF90_INQ_DIMID(ncid_ecmwf, "nlev", dimid), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQUIRE_DIMENSION(ncid_ecmwf, dimid, len = ninput), ierr)
    if (ierr .ne. 0) return

    count_prof = [ninput, 1]    
    allocate(hyam(ninput))
    allocate(hybm(ninput))
    allocate(q(ninput), t(ninput))

    call check(NF90_INQ_VARID(ncid_ecmwf, "hyam", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ecmwf, varid, hyam), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQ_VARID(ncid_ecmwf, "hybm", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ecmwf, varid, hybm), ierr)
    if (ierr .ne. 0) return

    !*** Get surface pressure
    call check(NF90_INQ_VARID(ncid_ecmwf, "sp", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ecmwf, varid, sp, start=start, count=count), ierr)
    if (ierr .ne. 0) return
    psurf = sp(1)    

    !*** Get surface height
    call check(NF90_INQ_VARID(ncid_ecmwf, "z", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ecmwf, varid, z, start=start, count=count), ierr)
    if (ierr .ne. 0) return
    ecmwf_z = z(1)

    count_prof = [ninput, 1]    
    !*** Get humidity
    call check(NF90_INQ_VARID(ncid_ecmwf, "q", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ecmwf, varid, q, start=start_prof, count=count_prof), ierr)
    if (ierr .ne. 0) return

    !*** Get temperature
    call check(NF90_INQ_VARID(ncid_ecmwf, "t", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ecmwf, varid, t, start=start_prof, count=count_prof), ierr)
    if (ierr .ne. 0) return

    !*** Get wind speed
    call check(NF90_INQ_VARID(ncid_ecmwf, "u10", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ecmwf, varid, u10, start=start, count=count), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQ_VARID(ncid_ecmwf, "v10", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ecmwf, varid, v10, start=start, count=count), ierr)
    if (ierr .ne. 0) return
    u = u10(1)
    v = v10(1)

!!$    !*** All information has been read, close NetCDF file. 
!!$    call check(NF90_CLOSE (ncid_ecmwf), ierr)
!!$    if (ierr .ne. 0) then
!!$       write(message,*) 'READ_MIPREP: error closing file'
!!$       ierr = ierr_open
!!$       goto 999
!!$    endif

    ierr = 0 ! Mark succes
    return

999 continue
    call stopretrieval(message)
    return

  end subroutine read_meteo

  !------------------------------------------------------------------------------  
  subroutine read_co(infile, ipixel, n_co, hyai, hybi, tm5_co, psurf, ierr)
    character(len=*), intent(in) :: infile
    integer, intent(in) :: ipixel
    real(double), dimension(:), allocatable, intent(out)  :: hyai, hybi, tm5_co
    real(double), intent(out) :: psurf
    integer, intent(out) :: n_co, ierr
    !***
    integer ::  varid, dimid
    integer :: start(1), count(1), start_prof(2), count_prof(2) 
    real(double) :: sp(1)
    character(stringlen) :: message

    !*** Get start indices
    start = [ipixel]
    count = [1]
    start_prof = [1, ipixel]

    !*** Read TM5 CO
!!$    !$OMP critical
!!$    call check(nf90_open(trim(infile)//'_co.nc', nf90_nowrite, ncid_co), ierr)
!!$    !$OMP end critical
!!$    if (ierr .ne. 0) then
!!$       ierr = ierr_open
!!$       write(message,*) 'READ_MIPREP: error opening file', trim(infile)//'_co.nc'
!!$       goto 999
!!$    endif

    !*** Get number of levels
    call check(NF90_INQ_DIMID(ncid_co, "nlev", dimid), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQUIRE_DIMENSION(ncid_co, dimid, len = n_co), ierr)
    if (ierr .ne. 0) return

    allocate(tm5_co(n_co))

    count_prof = [n_co, 1]   
    !*** Get VMR of CH4
    call check(NF90_INQ_VARID(ncid_co, "CO", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_co, varid, tm5_co, start=start_prof, count=count_prof), ierr)
    if (ierr .ne. 0) return

    !*** Get pressure grid
    allocate(hyai(n_co+1))
    allocate(hybi(n_co+1))
    call check(NF90_INQ_VARID(ncid_co, "hyai", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_co, varid, hyai), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQ_VARID(ncid_co, "hybi", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_co, varid, hybi), ierr)
    if (ierr .ne. 0) return

    !*** Get surface pressure
    call check(NF90_INQ_VARID(ncid_co, "sp", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_co, varid, sp, start=start, count=count), ierr)
    if (ierr .ne. 0) return
    psurf = sp(1)

!!$    !*** All information has been read, close NetCDF file. 
!!$    call check(NF90_CLOSE (ncid_co), ierr)
!!$    if (ierr .ne. 0) then
!!$       write(message,*) 'READ_MIPREP: error closing file'
!!$       ierr = ierr_open
!!$       goto 999
!!$    endif

    ierr = 0 ! Mark succes
    return
999 continue
    call stopretrieval(message)
    return

  end subroutine read_co
  !------------------------------------------------------------------------------  

  subroutine read_ch4(infile, ipixel, n_ch4, hyai, hybi, tm5_ch4, psurf, ierr)
    character(len=*), intent(in) :: infile
    integer, intent(in) :: ipixel
    real(double), dimension(:), allocatable, intent(out)  :: hyai, hybi, tm5_ch4
    real(double), intent(out) :: psurf
    integer, intent(out) :: n_ch4, ierr
    !***
    integer ::  varid, dimid
    integer :: start(1), count(1), start_prof(2), count_prof(2) 
    real(double) :: sp(1)
    character(stringlen) :: message

    !*** Get start indices
    start = [ipixel]
    count = [1]
    start_prof = [1, ipixel]
!!$  
!!$    !*** Read TM5 CH4
!!$    !$OMP critical
!!$    call check(nf90_open(trim(infile)//'_ch4.nc', nf90_nowrite, ncid_ch4), ierr)
!!$    !$OMP end critical
!!$    if (ierr .ne. 0) then
!!$       ierr = ierr_open
!!$       write(message,*) 'READ_MIPREP: error opening file', trim(infile)//'_ch4.nc'
!!$       goto 999
!!$    endif

    !*** Get number of levels
    call check(NF90_INQ_DIMID(ncid_ch4, "nlev", dimid), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQUIRE_DIMENSION(ncid_ch4, dimid, len = n_ch4), ierr)
    if (ierr .ne. 0) return

    allocate(tm5_ch4(n_ch4))

    count_prof = [n_ch4, 1]   
    !*** Get VMR of CH4
    call check(NF90_INQ_VARID(ncid_ch4, "CH4", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ch4, varid, tm5_ch4, start=start_prof, count=count_prof), ierr)
    if (ierr .ne. 0) return

    !*** Get pressure grid
    allocate(hyai(n_ch4+1))
    allocate(hybi(n_ch4+1))
    call check(NF90_INQ_VARID(ncid_ch4, "hyai", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ch4, varid, hyai), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQ_VARID(ncid_ch4, "hybi", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ch4, varid, hybi), ierr)
    if (ierr .ne. 0) return

    !*** Get surface pressure
    call check(NF90_INQ_VARID(ncid_ch4, "sp", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_ch4, varid, sp, start=start, count=count), ierr)
    if (ierr .ne. 0) return
    psurf = sp(1)

!!$    !*** All information has been read, close NetCDF file. 
!!$    call check(NF90_CLOSE (ncid_ch4), ierr)
!!$    if (ierr .ne. 0) then
!!$       write(message,*) 'READ_MIPREP: error closing file'
!!$       ierr = ierr_open
!!$       goto 999
!!$    endif

    ierr = 0 ! Mark succes
    return

999 continue
    call stopretrieval(message)
    return

    end subroutine read_ch4

    !------------------------------------------------------------------------------ 
    subroutine read_dem(infile, ipixel, surface_elevation, surface_elevation_stdv, landflag, ierr)
      character(len=*), intent(in) :: infile
      integer, intent(in) :: ipixel
      real(double) :: surface_elevation, surface_elevation_stdv
      integer, intent(out) :: landflag, ierr
      !***
      real(double) :: z(1)
      integer ::  varid
      integer :: start(1), count(1) 
      integer(kind=1) :: surface(1)
      character(stringlen) :: message

      !*** Get start indices
      start = [ipixel]
      count = [1]

      !***  Read DEM file
!!$      !$OMP critical
!!$      call check(nf90_open(trim(infile)//'_dem.nc', nf90_nowrite, ncid_dem), ierr)
!!$      !$OMP end critical
!!$      if (ierr .ne. 0) then
!!$         ierr = ierr_open
!!$         write(message,*) 'READ_MIPREP: error opening file', trim(infile)//'_dem.nc'
!!$         goto 999
!!$      endif

      call check(NF90_INQ_VARID(ncid_dem, "meanh", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_dem, varid, z, start=start, count=count), ierr)
      if (ierr .ne. 0) return
      surface_elevation = z(1)

      call check(NF90_INQ_VARID(ncid_dem, "stdh", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_dem, varid, z, start=start, count=count), ierr)
      if (ierr .ne. 0) return
      surface_elevation_stdv = z(1)

      call check(NF90_INQ_VARID(ncid_dem, "surface_classification", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_dem, varid, surface, start=start, count=count), ierr)
      if (ierr .ne. 0) return
      landflag = surface(1)

!!$      !*** All information has been read, close NetCDF file. 
!!$      call check(NF90_CLOSE (ncid_dem), ierr)
!!$      if (ierr .ne. 0) then
!!$         write(message,*) 'READ_MIPREP: error closing file', trim(infile)//'_dem.nc'
!!$         ierr = ierr_open
!!$         goto 999
!!$      endif

      ierr = 0 ! Mark succes
      return
999   continue
      call stopretrieval(message)
      return

    end subroutine read_dem

    !------------------------------------------------------------------------------ 
    subroutine read_orbit(infile, ipixel, sza, iza, saz, iaz, lat_center, lon_center, lat_corners, lon_corners, time, ierr)
      character(len=*), intent(in) :: infile
      integer, intent(in) :: ipixel
      real(double), intent(out) :: sza, iza, saz, iaz, lat_center, lon_center, lat_corners(4), lon_corners(4)
      integer, dimension(:) :: time
      integer, intent(out) :: ierr
      !***
      real(double) :: angle(1), lat(1), lon(1)
      integer ::  varid
      integer :: start(1), count(1), start_prof(2), count_prof(2)
      character(stringlen) :: message


      !*** Get start indices
      start = [ipixel]
      count = [1]
      start_prof = [1, ipixel]

!!$      !*** Read orbit file
!!$      !$OMP critical
!!$      call check(nf90_open(trim(infile)//'_orbit.nc', nf90_nowrite, ncid_orbit), ierr)
!!$      !$OMP end critical
!!$      if (ierr .ne. 0) then
!!$         write(message,*) 'READ_ORBIT: error opening file', trim(infile)//'_orbit.nc', ierr
!!$         ierr = ierr_open
!!$         goto 999
!!$      endif

      !*** Get sza
      call check(NF90_INQ_VARID(ncid_orbit, "solar_zenith_angle", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_orbit, varid, angle, start=start, count=count), ierr)
      if (ierr .ne. 0) return
      sza = angle(1)
 
      !*** Get vza
      call check(NF90_INQ_VARID(ncid_orbit, "viewing_zenith_angle", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_orbit, varid, angle, start=start, count=count), ierr)
      if (ierr .ne. 0) return
      iza = angle(1)

      !*** Get saz
      call check(NF90_INQ_VARID(ncid_orbit, "solar_azimuth_angle", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_orbit, varid, angle, start=start, count=count), ierr)
      if (ierr .ne. 0) return
      saz = angle(1)

      !*** Get vaz
      call check(NF90_INQ_VARID(ncid_orbit, "viewing_azimuth_angle", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_orbit, varid, angle, start=start, count=count), ierr)
      if (ierr .ne. 0) return
      iaz = angle(1)

      !*** Get latitude
      call check(NF90_INQ_VARID(ncid_orbit, "latitude", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_orbit, varid, lat, start=start, count=count), ierr)
      if (ierr .ne. 0) return
      lat_center = lat(1)

      !*** Get longitude
      call check(NF90_INQ_VARID(ncid_orbit, "longitude", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_orbit, varid, lon, start=start, count=count), ierr)
      if (ierr .ne. 0) return
      lon_center = lon(1)

      count_prof = [4, 1]  
      !*** Get latitude of corners
      call check(NF90_INQ_VARID(ncid_orbit, "latitude_corners", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_orbit, varid, lat_corners, start=start_prof, count=count_prof), ierr)
      if (ierr .ne. 0) return

      !*** Get longitude of corners
      call check(NF90_INQ_VARID(ncid_orbit, "longitude_corners", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_orbit, varid, lon_corners, start=start_prof, count=count_prof), ierr)
      if (ierr .ne. 0) return

      !*** Get time
      count_prof = [6, 1]   
      call check(NF90_INQ_VARID(ncid_orbit, "time", varid), ierr) 
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid_orbit, varid, time(1:6), start=start_prof, count=count_prof), ierr)
      if (ierr .ne. 0) return

!!$      !*** All information has been read, close NetCDF file. 
!!$      call check(NF90_CLOSE (ncid_orbit), ierr)
!!$      if (ierr .ne. 0) then
!!$         write(message,*) 'READ_ORBIT: error closing file', trim(infile)//'_orbit.nc'
!!$         ierr = ierr_open
!!$         goto 999
!!$      endif

      ierr = 0 ! Mark succes
      return
999   continue
      call stopretrieval(message)
      return

    end subroutine read_orbit

    !------------------------------------------------------------------------------

 subroutine read_co2(infile, ipixel, n_co2, ct_co2, ct_pressure, ierr)
    character(len=*), intent(in) :: infile
    integer, intent(in) :: ipixel
    real(double), dimension(:), allocatable, intent(out)  :: ct_co2, ct_pressure
    real(double), dimension(:), allocatable :: hyai, hybi
    integer, intent(out) :: n_co2, ierr
    !***
    integer ::  varid, dimid, k
    integer :: start(1), count(1), start_prof(2), count_prof(2) 
    real(double) :: sp(1)
    real(double) :: psurf
    character(stringlen) :: message

    !*** Get start indices
    start = [ipixel]
    count = [1]
    start_prof = [1, ipixel]


    !*** Get number of levels
    call check(NF90_INQ_DIMID(ncid_co2, "nlev", dimid), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQUIRE_DIMENSION(ncid_co2, dimid, len = n_co2), ierr)
    if (ierr .ne. 0) return

    allocate(ct_co2(n_co2))

    count_prof = [n_co2, 1]   
    !*** Get VMR of CO2
    call check(NF90_INQ_VARID(ncid_co2, "CO2", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_co2, varid, ct_co2, start=start_prof, count=count_prof), ierr)
    if (ierr .ne. 0) return

    !*** Get surface pressure
    call check(NF90_INQ_VARID(ncid_co2, "sp", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_co2, varid, sp, start=start, count=count), ierr)
    if (ierr .ne. 0) return
    psurf = sp(1)

    !*** Get pressure grid
    allocate(hyai(n_co2+1))
    allocate(hybi(n_co2+1))
    call check(NF90_INQ_VARID(ncid_co2, "hyai", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_co2, varid, hyai), ierr)
    if (ierr .ne. 0) return
    call check(NF90_INQ_VARID(ncid_co2, "hybi", varid), ierr) 
    if (ierr .ne. 0) return
    call check(NF90_GET_VAR(ncid_co2, varid, hybi), ierr)
    if (ierr .ne. 0) return
    allocate(ct_pressure(n_co2))

    do k = 1, n_co2
       !*** Convert hybrid coefficients to pressure [from Pa to hPa]
       ct_pressure(k) = ((hyai(k)+hyai(k+1))*0.5d0 + (hybi(k)+hybi(k+1))*0.5d0*psurf)*0.01     
    enddo


    ierr = 0 ! Mark succes
    return

999 continue
    call stopretrieval(message)
    return

    end subroutine read_co2
    !------------------------------------------------------------------------------ 
    !> Check calls to netCDF functions
    subroutine check(status, ierr)
      use NETCDF 
      integer, intent (in) :: status
      integer, intent(out) :: ierr

      if(status /= nf90_noerr) then
         ierr = ierr_read
!         call writelog('CHECK: '//trim(nf90_strerror(status)), 6)
      else
         ierr = 0
      endif

    end subroutine check
    !------------------------------------------------------------------------------ 
  end module read_miprep_module
