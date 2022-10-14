module read_xsdb_module
   use header_module
   use netcdf
   use hdf5
   implicit none

!*** Public procedures
   public :: read_xsdb, read_xsdb_netcdf, interpolate_xsdb, tri_conv_xsdb, rect_conv_xsdb, get_xs, read_xsdb_frankenberg, read_xsdb_butz, read_xsdb_absco
  
   private :: check

!------------------------------------------------------------------------------
!> @brief Cross sections database for one pressure value
!> @details Cross sections on temperature and wavelength grid for one pressure value 
   type, public :: cross_section_p
!      real(single), dimension(:,:), allocatable :: xs     !< Cross sections (dim: nwv, ntemp)
      real(single), dimension(:,:,:), allocatable :: xs    !< Cross sections (dim: nwv, ntemp, nratio)
      real(single), dimension(:), allocatable :: T      !< Temperature grid [K] (dim: ntemp)
      real(single) :: p                                 !< Pressure [hPa]
     
   end type cross_section_p
!------------------------------------------------------------------------------
!> @brief Cross section database
!> @details Cross sections on pressure grid 
   type, public :: cross_section_db
      real(double), dimension(:), allocatable :: wv           !< LUT wavelength grid [nm] (dim: nwv)
      integer :: species                                      !< HITRAN number of species
      integer :: nwv                                          !< Number of wavelengths in LUT
      type(cross_section_p), dimension(:), allocatable :: p   !< Pressure grid [hPa] (dim: npres)
! new parameters for ABSCO
      integer :: nvmr                                    !< Number of mixing ratio grids for ABSCO
      real(double), dimension(:), allocatable :: vmr!< Mixing ratio grids(dim: nratio,==3)      
      
   end type cross_section_db

   contains

!------------------------------------------------------------------------------
!> @brief Read in cross sections from one NetCDF file (S5P format)
     subroutine read_xsdb(xs_file, xsdb, ierr)
       !*** Input
       character(len=*), intent(in) :: xs_file                   ! Path of the database file
       !*** output
       type(cross_section_db), intent(inout) :: xsdb           ! type containg cross-section database
       integer, intent(out) :: ierr
       !*** local variables
       integer :: i, j, k, ngroups
       integer :: xsdb_species, xsdb_pnum,  xsdb_length
       integer, dimension(:),allocatable :: nT  
       integer, dimension(10) :: ncids
       integer :: ncid, MoleculenrVarID, IsotopenrVarID, SpectroscopyVarID, LineshapeVarID, AlgorithmVarID, PressureVarID
       integer::TemploVarID, TempupVarID, TempDTVarID, WaveVarID, PresDimID, WaveDimID     
       integer, dimension(:), allocatable :: TempDimID, SigmaVarID
       character(2):: levelstring          
       character(100):: iso_string_lut, spectro_string_lut, lineshape_lut, algorithm_lut 
       real(single) :: tempdt_lut 
       real(single), dimension(:), allocatable :: templo_lut, tempup_lut    
       real(single), dimension(:,:), allocatable :: xstemp   
       character(stringlen) :: message
       !------------------------------------------------------------------------------------------   
       !*** Initialize error identifier
       ierr = 0

       xsdb%nvmr = 1
       if(allocated(xsdb%vmr)) deallocate(xsdb%vmr)
       allocate(xsdb%vmr(xsdb%nvmr), stat=ierr)
       if (ierr .ne. 0) then
          write(message,*) 'READ_XSDB: memory allocation error'
          ierr = ierr_all
          goto 999
       endif

       !*** Open NetCDF LUT 
       call check(nf90_open(trim(xs_file), nf90_nowrite, ncid), ierr)
       if (ierr .ne. 0) then
          ierr = ierr_open
          write(message,*) 'READ_XSDB: error opening XSDB file'
          goto 999
       endif

       !*** Get #groups and their ncids
       call check(nf90_inq_grps(ncid, ngroups, ncids), ierr)
       if (ierr .ne. 0) return
       if (ngroups .gt. 10) then
          write(message,*) 'READ_XSDB: number of groups is too large'
          ierr = ierr_var
          goto 999
       endif

       do j = 1, ngroups

          call check(NF90_INQ_VARID(ncids(j), "Molecule", MoleculenrVarID), ierr) 
          if (ierr .ne. 0) return
          call check(NF90_GET_VAR(ncids(j), MoleculenrVarID, xsdb_species), ierr)
          if (ierr .ne. 0) return
          if ( abs(xsdb%species)==xsdb_species .or. xsdb%species/100 ==xsdb_species  ) then 
             !** Get dimensions for p, T(p), and wavelengths in LUT  
             call check(NF90_INQ_DIMID(ncids(j), "np", PresDimID), ierr)
             if (ierr .ne. 0) return
             call check(NF90_INQUIRE_DIMENSION(ncids(j), PresDimID, len = xsdb_pnum), ierr)
             if (ierr .ne. 0) return
             call check(NF90_INQ_DIMID(ncids(j), "nnu", WaveDimID), ierr)
             if (ierr .ne. 0) return
             call check(NF90_INQUIRE_DIMENSION(ncids(j), WaveDimID, len = xsdb_length), ierr) 
             if (ierr .ne. 0) return
             xsdb%nwv = xsdb_length 
             allocate(xsdb%wv(xsdb_length), &  
                  xsdb%p(xsdb_pnum), &  
                  NT(xsdb_pnum), & 
                  templo_lut(xsdb_pnum), & 
                  tempup_lut(xsdb_pnum), &
                  TempDimID(xsdb_pnum), & 
                  SigmaVarID(xsdb_pnum), &
                  stat=ierr)
             if (ierr .ne. 0) then
                write(message,*) 'READ_XSDB: memory allocation error'
                ierr = ierr_all
                goto 999
             endif
             do i = 1, xsdb_pnum
                write(levelstring,'(I2.2)') i
                call check(NF90_INQ_DIMID(ncids(j), "nt_p"//levelstring, TempDimID(i)), ierr)
                if (ierr .ne. 0) return
                call check(NF90_INQUIRE_DIMENSION(ncids(j), TempDimID(i), len = NT(i)), ierr)
                if (ierr .ne. 0) return
             enddo
             call check(NF90_INQ_VARID(ncids(j), "Isotopes", IsotopenrVarID), ierr)  
             if (ierr .ne. 0) return
             call check(NF90_INQ_VARID(ncids(j), "Spectroscopy", SpectroscopyVarID), ierr)
             if (ierr .ne. 0) return
             call check(NF90_INQ_VARID(ncids(j), "Lineshape", LineshapeVarID), ierr)
             if (ierr .ne. 0) return
             call check(NF90_INQ_VARID(ncids(j), "Algorithm", AlgorithmVarID), ierr)
             if (ierr .ne. 0) return
             call check(NF90_INQ_VARID(ncids(j), "Pressure", PressureVarID), ierr)
             if (ierr .ne. 0) return
             call check(NF90_INQ_VARID(ncids(j), "Tlow", TemploVarID), ierr) 
             if (ierr .ne. 0) return
             call check(NF90_INQ_VARID(ncids(j), "Thigh", TempupVarID), ierr) 
             if (ierr .ne. 0) return
             call check(NF90_INQ_VARID(ncids(j), "dT", TempDTVarID), ierr)
             if (ierr .ne. 0) return
             call check(NF90_INQ_VARID(ncids(j), "nu", WaveVarID), ierr)
             if (ierr .ne. 0) return
             do i = 1, xsdb_pnum
                write(levelstring,'(I2.2)') i
                call check(NF90_INQ_VARID(ncids(j), "cross_p"//levelstring, SigmaVarID(i)), ierr)
                if (ierr .ne. 0) return
             enddo
             xsdb%species = xsdb_species
             call check(NF90_GET_VAR(ncids(j), IsotopenrVarID, iso_string_lut), ierr)
             if (ierr .ne. 0) return
             call check(NF90_GET_VAR(ncids(j), SpectroscopyVarID, spectro_string_lut), ierr)
             if (ierr .ne. 0) return
             call check(NF90_GET_VAR(ncids(j), LineshapeVarID, lineshape_lut), ierr)
             if (ierr .ne. 0) return
             call check(NF90_GET_VAR(ncids(j), AlgorithmVarID, algorithm_lut), ierr)
             if (ierr .ne. 0) return
             call check(NF90_GET_VAR(ncids(j), PressureVarID, xsdb%p(:)%p), ierr)
             if (ierr .ne. 0) return
             call check(NF90_GET_VAR(ncids(j), TemploVarID, templo_lut), ierr)
             if (ierr .ne. 0) return
             call check(NF90_GET_VAR(ncids(j), TempupVarID, tempup_lut), ierr) 
             if (ierr .ne. 0) return
             call check(NF90_GET_VAR(ncids(j), TempDTVarID, tempdt_lut), ierr)
             if (ierr .ne. 0) return
             call check(NF90_GET_VAR(ncids(j), WaveVarID, xsdb%wv), ierr)
             if (ierr .ne. 0) return
             do i = 1, xsdb_pnum
                allocate(xsdb%p(i)%T(nT(i)), stat=ierr)
                if (ierr .ne. 0) then
                   write(message,*) 'READ_XSDB: memory allocation error'
                   ierr = ierr_all
                   goto 999
                endif
                do k = 1, NT(i)
                   xsdb%p(i)%T(k) = templo_lut(i) + real(k-1)*tempdt_lut
                enddo
             enddo

             !*** Assume database works with wavenumber and RemoTeC works with wavelength
             xsdb%wv = 1.d7/xsdb%wv(xsdb_length:1:-1)

             !*** Now read in the cross-sections and store them in the arrays win_ini(n)%xsdb(j)%p(i)%xs(nT(i),l)
             do i = 1, xsdb_pnum
                if (allocated(xstemp)) deallocate(xstemp, stat=ierr)
                if (ierr .ne. 0) then
                   write(message,*) 'READ_XSDB: memory deallocation error'
                   ierr = ierr_deall
                   goto 999
                endif
                !                allocate(xsdb%p(i)%xs(xsdb%nwv, nT(i)), xstemp(nT(i),xsdb%nwv), stat=ierr)
                allocate(xsdb%p(i)%xs(xsdb%nwv, xsdb%nvmr, nT(i)), xstemp(nT(i),xsdb%nwv), stat=ierr)
                if (ierr .ne. 0) then
                   write(message,*) 'READ_XSDB: memory allocation error'
                   ierr = ierr_all
                   goto 999
                endif
                call check(NF90_GET_VAR(ncids(j), SigmaVarID(i), xstemp), ierr)
                do k = 1, nT(i) 
                   ! Flip first dimension due to conversion from wavenumber to wavelength
                   xsdb%p(i)%xs(xsdb_length:1:-1, 1, k) = xstemp(k, :)
                enddo
                if (ierr .ne. 0) return	      
             enddo
             deallocate(nT, stat=ierr)     
             if (ierr .ne. 0) then
                write(message,*) 'READ_XSDB: memory deallocation error'
                ierr = ierr_deall
                goto 999
             endif

             exit  ! Molecule found, exit loop
          else
             cycle ! Molecule not found, continue loop
          endif
       enddo

       !*** All information has been read, close NetCDF file. 
       call check(NF90_CLOSE (ncid), ierr)
       if (ierr .ne. 0) then
          write(message,*) 'READ_XSDB: error closing file'
          ierr = ierr_open
          goto 999
       endif

       ierr = 0 ! Mark succes
       return

999    continue
       call stopretrieval(message)
       return

     end subroutine read_xsdb
!------------------------------------------------------------------------------
!> @brief Read in cross sections from NetCDF file
     subroutine read_xsdb_netcdf(xs_file, xsdb, ierr)
       !*** Input
       character(len=*), intent(in) :: xs_file                   ! Path of the database file
       !*** output
       type(cross_section_db), intent(inout) :: xsdb           ! type containg cross-section database
       integer, intent(out) :: ierr
       !*** local variables
       integer :: i, k
       integer :: xsdb_species, xsdb_pnum,  xsdb_length
       integer, dimension(:),allocatable :: nT  
       integer :: ncid, MoleculenrVarID, IsotopenrVarID, SpectroscopyVarID, LineshapeVarID, AlgorithmVarID, PressureVarID
       integer::TemploVarID, TempupVarID, TempDTVarID, WaveVarID, PresDimID, WaveDimID     
       integer, dimension(:), allocatable :: TempDimID, SigmaVarID
       character(2):: levelstring          
       character(100):: iso_string_lut, spectro_string_lut, lineshape_lut, algorithm_lut 
       real(single) :: tempdt_lut 
       real(single), dimension(:), allocatable :: templo_lut, tempup_lut    
       real(single), dimension(:,:), allocatable :: xstemp
       !------------------------------------------------------------------------------------------   

       xsdb%nvmr = 1
       if(allocated(xsdb%vmr)) deallocate(xsdb%vmr)
       allocate(xsdb%vmr(xsdb%nvmr), stat=ierr)
       
       !*** Open NetCDF LUT, get variables, and read descriptions
       call check(nf90_open(xs_file, nf90_nowrite, ncid), ierr)
       !      write(*,*) 'Opened netCDF file: ', trim(xs_file), ' in remotec_core/read_xsdb/read_xsdb_netcdf (line 57)'
       !** Get dimensions for p, T(p), and wavelengths in LUT  
       call check(NF90_INQ_DIMID(ncid, "Pressure length", PresDimID), ierr)
       call check(NF90_INQUIRE_DIMENSION(ncid, PresDimID, len = xsdb_pnum), ierr)
       call check(NF90_INQ_DIMID(ncid, "Wavenumber length", WaveDimID), ierr)
       call check(NF90_INQUIRE_DIMENSION(ncid, WaveDimID, len = xsdb_length), ierr) 
       xsdb%nwv = xsdb_length 
       allocate(xsdb%wv(xsdb_length))  
       allocate(xsdb%p(xsdb_pnum))  
       allocate(NT(xsdb_pnum)) 
       allocate(templo_lut(xsdb_pnum)) 
       allocate(tempup_lut(xsdb_pnum)) 
       allocate(TempDimID(xsdb_pnum)) 
       allocate(SigmaVarID(xsdb_pnum)) 
       do i = 1, xsdb_pnum
          write(levelstring,'(I2.2)') i
          call check(NF90_INQ_DIMID(ncid, "Temperature length at pressure level "//levelstring, TempDimID(i)), ierr)
          call check(NF90_INQUIRE_DIMENSION(ncid, TempDimID(i), len = NT(i)), ierr)
       enddo

       call check(NF90_INQ_VARID(ncid, "Molecule number", MoleculenrVarID), ierr) 
       call check(NF90_INQ_VARID(ncid, "Isotope numbers", IsotopenrVarID), ierr)  
       call check(NF90_INQ_VARID(ncid, "Spectroscopy input", SpectroscopyVarID), ierr)
       call check(NF90_INQ_VARID(ncid, "Line shape", LineshapeVarID), ierr)
       call check(NF90_INQ_VARID(ncid, "Algorithm", AlgorithmVarID), ierr)
       call check(NF90_INQ_VARID(ncid, "Pressure", PressureVarID), ierr)
       call check(NF90_INQ_VARID(ncid, "Temperature, lower boundaries", TemploVarID), ierr) 
       call check(NF90_INQ_VARID(ncid, "Temperature, upper boundaries", TempupVarID), ierr) 
       call check(NF90_INQ_VARID(ncid, "Temperature steps", TempDTVarID), ierr)
       call check(NF90_INQ_VARID(ncid, "Wavenumber", WaveVarID), ierr)
       do i = 1, xsdb_pnum
          write(levelstring,'(I2.2)') i
          call check(NF90_INQ_VARID(ncid, "Cross-sections at pressure level "//levelstring, SigmaVarID(i)), ierr)
       enddo

       call check(NF90_GET_VAR(ncid, MoleculenrVarID, xsdb_species), ierr)
       xsdb%species = xsdb_species
       call check(NF90_GET_VAR(ncid, IsotopenrVarID, iso_string_lut), ierr)
       call check(NF90_GET_VAR(ncid, SpectroscopyVarID, spectro_string_lut), ierr)
       call check(NF90_GET_VAR(ncid, LineshapeVarID, lineshape_lut), ierr)
       call check(NF90_GET_VAR(ncid, AlgorithmVarID, algorithm_lut), ierr)
       call check(NF90_GET_VAR(ncid, PressureVarID, xsdb%p(:)%p), ierr)
       call check(NF90_GET_VAR(ncid, TemploVarID, templo_lut), ierr)
       call check(NF90_GET_VAR(ncid, TempupVarID, tempup_lut), ierr) 
       call check(NF90_GET_VAR(ncid, TempDTVarID, tempdt_lut), ierr)
       call check(NF90_GET_VAR(ncid, WaveVarID, xsdb%wv), ierr)

       !*** Assume database works with wavenumber and RemoTeC works with wavelength
       xsdb%wv = 1.d7/xsdb%wv(xsdb_length:1:-1)

       do i = 1, xsdb_pnum
          allocate(xsdb%p(i)%T(nT(i)))
          do k = 1, NT(i)
             xsdb%p(i)%T(k) = templo_lut(i) + real(k-1)*tempdt_lut
          enddo
       enddo

       !*** Now read in the cross-sections and store them in the arrays win_ini(n)%xsdb(j)%p(i)%xs(nT(i),l)
       do i = 1, xsdb_pnum
          if (allocated(xstemp)) deallocate(xstemp)	         
!          allocate(xsdb%p(i)%xs(xsdb%nwv,nT(i)),  xstemp(nT(i), xsdb%nwv) )		            
          allocate(xsdb%p(i)%xs(xsdb%nwv, xsdb%nvmr, nT(i)),  xstemp(nT(i), xsdb%nwv) )		            
          !        call check(NF90_GET_VAR(ncid, SigmaVarID(i), xsdb%p(i)%xs), ierr) 
          call check(NF90_GET_VAR(ncid, SigmaVarID(i), xstemp), ierr)
          do k = 1, nT(i) 
             ! Flip first dimension due to conversion from wavenumber to wavelength
             xsdb%p(i)%xs(xsdb_length:1:-1, 1, k) = xstemp(k, :)
          enddo
       enddo
       deallocate(nT)              

       !*** All information has been read, close NetCDF file. 
       call check(NF90_CLOSE (ncid), ierr)

     end subroutine read_xsdb_netcdf
!------------------------------------------------------------------------------
!> @brief Interpolate cross section LUT on different wavelength grid
!------------------------------------------------------------------------------
     subroutine interpolate_xsdb (wgrid_out, xsdb, ierr)
       real(double), dimension(:), intent(in) :: wgrid_out
       type(cross_section_db), intent(inout) ::  xsdb
       integer, intent(out) :: ierr
       !*** local variables
       integer :: i, j, ntemp, nwv_out
       real(double), dimension(xsdb%nwv) :: xsin 
       real(double), dimension(size(wgrid_out)) :: xsout
!       real(double), dimension(:,:), allocatable :: xstemp
       real(double), dimension(:,:,:), allocatable :: xstemp
       character(stringlen) :: message
       integer :: iratio


       nwv_out = size(wgrid_out)
       !*** Check if LUT covers the wavelength range
       if ( (xsdb%wv(1) > wgrid_out(1)) .or. (xsdb%wv(xsdb%nwv) < wgrid_out(nwv_out)) ) then
          write(message, *)'INTERPOLATE_XSDB: spectral range of XSDB does not cover retrieval window.'
          ierr = ierr_var
          goto 999    
       endif

       do i = 1, size(xsdb%p)
          ntemp = size(xsdb%p(i)%T)
          if(allocated(xstemp)) deallocate(xstemp, stat=ierr)   
          if (ierr .ne. 0) then
             write(message, *)'INTERPOLATE_XSDB: memory deallocation error.'
             ierr = ierr_deall
             goto 999    
          endif
!          allocate(xstemp(xsdb%nwv, ntemp), stat=ierr)
          allocate(xstemp(xsdb%nwv, xsdb%nvmr, ntemp), stat=ierr)
          if (ierr .ne. 0) then
             write(message, *)'INTERPOLATE_XSDB: memory allocation error.'
             ierr = ierr_all
             goto 999    
          endif

          xstemp = dble(xsdb%p(i)%xs)
          deallocate(xsdb%p(i)%xs)
!          allocate(xsdb%p(i)%xs(nwv_out, ntemp), stat=ierr)
          allocate(xsdb%p(i)%xs(nwv_out, xsdb%nvmr, ntemp), stat=ierr)
          if (ierr .ne. 0) then
             write(message, *)'INTERPOLATE_XSDB: memory allocation error.'
             ierr = ierr_all
             goto 999    
          endif
          do j = 1, ntemp
!             xsin(:) = xstemp(:, j)
		do iratio = 1, xsdb%nvmr
	             xsin(:) = xstemp(:, iratio, j)
	             call linterp(xsdb%wv, xsin, xsdb%nwv, &
	                  wgrid_out, xsout, nwv_out )
	             xsdb%p(i)%xs(:, iratio, j) = real(xsout)
		enddo
          enddo
       enddo

       xsdb%nwv = nwv_out
       deallocate(xsdb%wv, stat=ierr)
       if (ierr .ne. 0) then
          write(message, *)'INTERPOLATE_XSDB: memory deallocation error.'
          ierr = ierr_deall
          goto 999    
       endif
       allocate(xsdb%wv(nwv_out), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *)'INTERPOLATE_XSDB: memory allocation error.'
          ierr = ierr_all
          goto 999    
       endif
       xsdb%wv = wgrid_out

       ierr = 0 ! Mark success
       return

999    continue 
       call stopretrieval(message)
       return

     end subroutine interpolate_xsdb

!------------------------------------------------------------------------------
   !> @brief Triangle convolution for cross sections of the absorbers.
   !> @details Cross sections can change rapidly with wavelength. Therefore, a coarse sampling of
   !! cross sections over wavelength may be prone to biases due to coincidentally sampling
   !! on a peak or just around a peak. With triangle convolution, a triangular weight function
   !! is applied so that all values from the lookup table are used to construct the cross
   !! sections on the internal domain. The downside of triangle convolution is that different
   !! cross sections are averaged. Because of the non-linear radiative effect of cross sections,
   !! linear averaging systematically overestimates this radiative effect. Therefore, the
   !! generalized mean is used with a chosen averaging power:
   !! \f$ \overline\sigma = \sqrt[m]{\frac{\sum_i w_i\sigma_i^m}{\sum_i w_i}} \f$. 
   !! Here, \f$ w \f$ is the triangular weight function and 
   !! \f$ m \f$ is the chosen averaging power. See section 5.4 of the CO ATBD.
   subroutine tri_conv_xsdb(wgrid_out, xsdb, ierr)
     real(double), dimension(:), intent(in) :: wgrid_out
     type(cross_section_db), intent(inout) ::  xsdb
     integer, intent(out) :: ierr
     !*** local variables
     real(double), dimension(xsdb%nwv) :: xsin 
     real(double), dimension(size(wgrid_out)) :: xsout
!     real(double), dimension(:,:), allocatable :: xstemp
     real(double), dimension(:,:,:), allocatable :: xstemp
     real(double), parameter :: averagingpower = 0.85d0
     real(double) :: internal_sampling, floatidx, fractionidx
     ! Normalization vector for the convolution product. This will contain
     ! the total weight factors of the contributing lookup table entries.
     real(double), dimension(size(wgrid_out)) :: normalization
     integer :: i, j, intidx, inu, inu_lut, nwv_out, ntemp
     character(stringlen) :: message
     integer :: iratio

     !*** Initialize
     ierr = 0
     nwv_out = size(wgrid_out)
     internal_sampling = wgrid_out(2)-wgrid_out(1)

     !*** Check if LUT sampling is finer than internal sampling
     if ( xsdb%wv(2) -xsdb%wv(1) > internal_sampling ) then
        write(message, *)'TRI_CONV_XSDB: spectral sampling of XSDB is coarser than internal sampling of model spectrum '
        ierr = ierr_var
        goto 999    
     endif

     !*** Check if LUT covers the wavelength range
     if ( (xsdb%wv(1) > wgrid_out(1)) .or. (xsdb%wv(xsdb%nwv) < wgrid_out(nwv_out)) ) then
        write(message, *)'TRI_CONV_XSDB: spectral range of XSDB does not cover retrieval window.'
        ierr = ierr_var
        goto 999    
     endif

     do i = 1, size(xsdb%p)
        ntemp = size(xsdb%p(i)%T)
        if(allocated(xstemp)) deallocate(xstemp)   
!        allocate(xstemp(xsdb%nwv, ntemp), stat=ierr)
        allocate(xstemp(xsdb%nwv, xsdb%nvmr, ntemp), stat=ierr)
        if (ierr .ne. 0) then
           write(message, *)'TRI_CONV_XSDB: memory allocation error.'
           ierr = ierr_var
           goto 999    
        endif

        xstemp = dble(xsdb%p(i)%xs)
        deallocate(xsdb%p(i)%xs)
!        allocate(xsdb%p(i)%xs(nwv_out, ntemp), stat=ierr)
        allocate(xsdb%p(i)%xs(nwv_out, xsdb%nvmr, ntemp), stat=ierr)
        if (ierr .ne. 0) then
           write(message, *)'TRI_CONV_XSDB: memory allocation error.'
           ierr = ierr_var
           goto 999    
        endif
        do j = 1, ntemp
	   do iratio = 1, xsdb%nvmr
!           xsin(:) = xstemp(:, j)
           xsin(:) = xstemp(:, iratio, j)
!!$           call linterp(xsdb%wv, xsin, xsdb%nwv, &
!!$                wgrid_out, xsout1, nwv_out )
           ! Initialize result to zero
           xsout = 0.d0
           ! Initialize normalization array
           normalization = 0.d0
           ! Loop over wavelengths from the original lookup table
           do inu_lut = 1, xsdb%nwv
              ! Calculate on which index the wavelength of the lookup table
              ! would be if it was in the line-by-line grid. This will
              ! be a floating point, because the grid point generally do not
              ! overlap. 
              floatidx = 1.0d0 + (xsdb%wv(inu_lut) - wgrid_out(1)) / internal_sampling
              intidx = floor(floatidx)
              fractionidx = modulo(floatidx,1.0d0)

              ! The lookup table
              ! cross section will be applied to two indices with different
              ! weight factors (one at edges). And the weight factors will be saved so that
              ! normalization can be applied. Furthermore, an averaging power
              ! from the namelist is applied. This power should account for the
              ! fact that the average effect of different cross section
              ! is generally not equal to the arithmetic mean of the cross
              ! sections. .
              if (intidx .lt. 0) then
                 cycle ! Domain not yet reached             
              else if (intidx .eq. 0) then
                 ! Do not write index 0, but do write index 1
                 xsout(intidx+1) = xsout(intidx+1) + fractionidx * xsin(inu_lut)**averagingpower 
                 normalization(intidx+1) = normalization(intidx+1) + fractionidx
              elseif (intidx > 0 .and. intidx < nwv_out ) then
                 ! Write both intidx and intidx+1
                 xsout(intidx) = xsout(intidx) + (1.d0-fractionidx) * xsin(inu_lut)**averagingpower 
                 xsout(intidx+1) = xsout(intidx+1) + fractionidx * xsin(inu_lut)**averagingpower 
                 normalization(intidx) = normalization(intidx) + 1.d0 - fractionidx
                 normalization(intidx+1) = normalization(intidx+1) + fractionidx
              elseif (intidx .eq. nwv_out) then
                 ! Write index nwv_out but do not write nwv_out+1
                 xsout(intidx) = xsout(intidx) + (1.d0-fractionidx) * xsin(inu_lut)**averagingpower 
                 normalization(intidx) = normalization(intidx) + 1.d0-fractionidx
              elseif (intidx .gt. nwv_out) then 
                 exit ! Moved out of the domain  
              endif
           enddo !inu_lut

           ! Apply the normalization and averaging power
           do inu = 1, nwv_out !band_instance%nnu
              xsout(inu) = (xsout(inu)/ normalization(inu))**(1.0d0/averagingpower)
           enddo
!!$open(50, file='xsin.dat')
!!$do inu=1, xsdb%nwv
!!$write(50,*) xsdb%wv(inu), xsin(inu)
!!$enddo 
!!$close(50)
!!$
!!$open(50, file='xsout.dat')
!!$do inu=1, nwv_out
!!$write(50,*) wgrid_out(inu), xsout1(inu), xsout(inu)
!!$enddo 
!!$close(50)
!!$stop
!           xsdb%p(i)%xs(:, j) = real(xsout)
           xsdb%p(i)%xs(:, iratio, j) = real(xsout)
	    enddo ! loop over different ratio
        enddo !loop over temperature indices
     enddo ! loop over pressure indices

     !*** Copy custom wavelength grid to type xsdb   
     xsdb%nwv = nwv_out
     deallocate(xsdb%wv)
     allocate(xsdb%wv(nwv_out), stat=ierr)
     if (ierr .ne. 0) then
        write(message, *)'TRI_CONV_XSDB: memory allocation error.'
        ierr = ierr_var
        goto 999    
     endif
     xsdb%wv = wgrid_out

     return
999  continue 
     call stopretrieval(message)

   end subroutine tri_conv_xsdb
!------------------------------------------------------------------------------
   subroutine rect_conv_xsdb(wgrid_out, xsdb, ierr)
     real(double), dimension(:), intent(in) :: wgrid_out
     type(cross_section_db), intent(inout) ::  xsdb
     integer, intent(out) :: ierr
     !*** local variables
     real(double), dimension(xsdb%nwv) :: xsin 
     real(double), dimension(size(wgrid_out)) :: xsout
!     real(double), dimension(:,:), allocatable :: xstemp
     real(double), dimension(:,:,:), allocatable :: xstemp
     real(double), parameter :: averagingpower = 1.d0
     real(double) :: internal_sampling, floatidx
     ! Normalization vector for the convolution product. This will contain
     ! the total weight factors of the contributing lookup table entries.
     real(double), dimension(size(wgrid_out)) :: normalization
     integer :: i, j, intidx, inu, inu_lut, nwv_out, ntemp
     character(stringlen) :: message
     integer :: iratio

    !*** Initialize
     ierr = 0
     nwv_out = size(wgrid_out)
     internal_sampling = wgrid_out(2)-wgrid_out(1)

     !*** Check if LUT sampling is finer than internal sampling
     if ( xsdb%wv(2) -xsdb%wv(1) > internal_sampling ) then
        write(message, *)'RECT_CONV_XSDB: spectral sampling of XSDB is coarser than internal sampling of model spectrum '
        ierr = ierr_var
        goto 999    
     endif

     !*** Check if LUT covers the wavelength range
     if ( (xsdb%wv(1) > wgrid_out(1)) .or. (xsdb%wv(xsdb%nwv) < wgrid_out(nwv_out)) ) then
        write(message, *)'RECT_CONV_XSDB: spectral range of XSDB does not cover retrieval window.'
        ierr = ierr_var
        goto 999    
     endif

    do i = 1, size(xsdb%p)
        ntemp = size(xsdb%p(i)%T)
        if(allocated(xstemp)) deallocate(xstemp)   
!        allocate(xstemp(xsdb%nwv, ntemp), stat=ierr)
        allocate(xstemp(xsdb%nwv,xsdb%nvmr, ntemp), stat=ierr)
        if (ierr .ne. 0) then
           write(message, *)'RECT_CONV_XSDB: memory allocation error.'
           ierr = ierr_var
           goto 999    
        endif

        xstemp = dble(xsdb%p(i)%xs)
        deallocate(xsdb%p(i)%xs)
!        allocate(xsdb%p(i)%xs(nwv_out, ntemp), stat=ierr)
        allocate(xsdb%p(i)%xs(nwv_out, xsdb%nvmr, ntemp), stat=ierr)
        if (ierr .ne. 0) then
           write(message, *)'RECT_CONV_XSDB: memory allocation error.'
           ierr = ierr_var
           goto 999    
        endif
        do j = 1, ntemp
	   do iratio = 1, xsdb%nvmr
           xsin(:) = xstemp(:, iratio, j)
!!$             call linterp(xsdb%wv, xsin, xsdb%nwv, &
!!$                wgrid_out, xsout1, nwv_out )

           ! Initialize result to zero
           xsout = 0.d0
           ! Initialize normalization array
           normalization = 0.d0
           ! Loop over wavelengths from the original lookup table
           do inu_lut = 1, xsdb%nwv
              ! Calculate on which index the wavelength of the lookup table
              ! would be if it was in the line-by-line grid. This will
              ! be a floating point, because the grid point generally do not
              ! overlap. 
              floatidx = 1.0d0 + (xsdb%wv(inu_lut) - wgrid_out(1)) / internal_sampling
              intidx = nint(floatidx)
              if (intidx .le. 0) then
                 cycle ! Domain not yet reached             
              elseif (intidx .ge. 1 .and. intidx .le. nwv_out ) then
                 ! Write both intidx and intidx+1
                 xsout(intidx) = xsout(intidx) + xsin(inu_lut)**averagingpower 
                 normalization(intidx) = normalization(intidx) + 1.d0
              else 
                 exit ! Moved out of the domain  
              endif
           enddo !inu_lut

           ! Apply the normalization and averaging power
           do inu = 1, nwv_out !band_instance%nnu
              xsout(inu) = (xsout(inu)/ normalization(inu))**(1.0d0/averagingpower)
           enddo
!!$open(50, file='xsin.dat')
!!$do inu=1, xsdb%nwv
!!$write(50,*) xsdb%wv(inu), xsin(inu)
!!$enddo 
!!$close(50)
!!$
!!$open(50, file='xsout.dat')
!!$do inu=1, nwv_out
!!$write(50,*) wgrid_out(inu), xsout1(inu), xsout(inu)
!!$enddo 
!!$close(50)
!!$stop

!           xsdb%p(i)%xs(:, j) = real(xsout)
           xsdb%p(i)%xs(:, iratio, j) = real(xsout)
           enddo
        enddo
     enddo

     xsdb%nwv = nwv_out
     deallocate(xsdb%wv)
     allocate(xsdb%wv(nwv_out), stat=ierr)
     if (ierr .ne. 0) then
        write(message, *)'RECT_CONV_XSDB: memory allocation error.'
        ierr = ierr_var
        goto 999    
     endif
     xsdb%wv = wgrid_out


     return
999  continue 
     call stopretrieval(message)


   end subroutine rect_conv_xsdb

!------------------------------------------------------------------------------
!> @details Get cross sections by interpolating in the pressure-temperature grid
!------------------------------------------------------------------------------
   subroutine get_xs(xsdb, vmr_k, xspress, xstemp, xsout, ExitXSFlag, ierr)
     implicit none
     type(cross_section_db), intent(in) :: xsdb
     real(double), intent(in) :: xspress, xstemp             ! pressure, temperature
     real(double), intent(in) :: vmr_k !volumn mixing ratio at k layer
     !*** Output
     real(double), dimension(:), intent(out) :: xsout  ! Molecular absorption cross sections on wavelength-grid
     integer, intent(out) :: ExitXSFlag, ierr
     !*** Local variables
     integer :: iup, ilo, iup_plo, ilo_plo, iup_phi, ilo_phi, imid, np, nT_plo, &
          nT_phi, phi,plo, Thi_plo, Tlo_plo, Thi_phi, Tlo_phi, i
     real(double) :: t, u_plo, u_phi, u1_plo, u1_phi, t1
     character(stringlen) :: message
     !absco

     integer :: iup_vmr, ilo_vmr! up and low index for vmr
     real(double):: u_vmr!
     !--------------------------------------------------------------------------------------
     !*** Intialize
     ierr = 0   
     ExitXSFlag = 0

     if (xsdb%nwv .ne. size(xsout)) then
        write(message, *) 'GET_XS: spectral grid of XSDB does not match model spectral grid'       
        ierr = ierr_var
        goto 999   
     endif

     !***Find pressure index in the XS database (NR: locate)     
     np = size(xsdb%p)
     !      iup = np + 1  
     !      ilo = 0
     !HH: Above can cause segfaults (e.g. if pressure is NaN)
     iup = np 
     ilo = 1 
     do 
        if(iup-ilo>1) then
           imid = (iup+ilo)/2
           if((xsdb%p(np)%p >= xsdb%p(1)%p) &
                .eqv. xspress >= xsdb%p(imid)%p) then
              ilo = imid
           else
              iup = imid
           endif
        else
           exit
        endif
     enddo

     if(xspress <= dble(xsdb%p(1)%p)) then
        plo = 1
     else if(xspress >= dble( xsdb%p(np)%p)) then
        plo = np-1
     else
        plo = ilo
     endif
     phi = plo + 1

     !*** Find temperature index in the XS database
     nT_plo = size(xsdb%p(plo)%T)
     nT_phi = size(xsdb%p(phi)%T)
!!$      iup_plo = nT_plo+1
!!$      iup_phi = nT_phi+1
!!$      ilo_plo = 0
!!$      ilo_phi = 0
     !HH: Above can cause segfaults
     iup_plo = nT_plo
     iup_phi = nT_phi
     ilo_plo = 1
     ilo_phi = 1
     do 
        if(iup_plo-ilo_plo>1) then
           imid = (iup_plo+ilo_plo)/2
           if((xsdb%p(plo)%T(nT_plo) >= xsdb%p(plo)%T(1)) .eqv. &
                xstemp >= xsdb%p(plo)%T(imid)) then
              ilo_plo = imid
           else
              iup_plo = imid
           endif
        else
           exit
        endif
     enddo

     do 
        if(iup_phi-ilo_phi>1) then
           imid=(iup_phi+ilo_phi)/2
           if((xsdb%p(phi)%T(nT_phi) >= xsdb%p(phi)%T(1)) &
                .eqv. xstemp >= xsdb%p(phi)%T(imid)) then
              ilo_phi = imid
           else
              iup_phi = imid
           endif
        else
           exit
        endif
     enddo

     !*** If T < xsdb min T then take lower boundary
     !*** If T > xsdb max T then take upper boundary
     if(xstemp <= dble(xsdb%p(plo)%T(1))) then
        Tlo_plo = 1
     else if(xstemp >= dble(xsdb%p(plo)%T(nT_plo))) then
        Tlo_plo = nT_plo - 1
     else
        Tlo_plo = ilo_plo
     endif
     Thi_plo = Tlo_plo  +1

     if(xstemp <= dble(xsdb%p(phi)%T(1))) then
        Tlo_phi = 1
     else if(xstemp >= dble(xsdb%p(phi)%T(nT_phi))) then
        Tlo_phi = nT_phi - 1
     else
        Tlo_phi = ilo_phi
     endif
     Thi_phi = Tlo_phi + 1

     !*** find low and up boundary for vmr
!!$     ilo_vmr = 1! low index for vmr
!!$     iup_vmr = 2! up index for vmr     
     if (xsdb%nvmr==1) then
        ilo_vmr = 1
        iup_vmr = 1     
     else if (vmr_k>=xsdb%vmr(xsdb%nvmr)) then
        ilo_vmr = xsdb%nvmr
        iup_vmr = xsdb%nvmr
     else
        do i = 1, xsdb%nvmr-1
           if (vmr_k>= xsdb%vmr(i).and.vmr_k<xsdb%vmr(i+1)) then
              ilo_vmr = i
              iup_vmr = i+1
              exit
           endif
        enddo
     endif
   
     !     if(abs(xsdb%species)==7 .or. abs(xsdb%species)==2)then
     !       if( vmr_k  < 0.00   .or. & 
     !           vmr_k  > 0.09) then
     !	   print*, 'vmr_k', vmr_k
     !           call writelog('GET_XS: Boundaries of XSDB LUT have been crossed!', 6)
     !           !ExitXSFlag = 1
     !       endif
     !     endif

     if(xstemp < dble(xsdb%p(plo)%T(1)) .or.  &
          xstemp > dble(xsdb%p(plo)%T(nT_plo)) .or. &
          xstemp < dble(xsdb%p(phi)%T(1)) .or. &
          xstemp > dble(xsdb%p(phi)%T(nT_phi)) .or.&
          xspress < dble(xsdb%p(1)%p) .or. &
          xspress > dble(xsdb%p(np)%p)) then
!!!print*, 'TEST temp', xstemp
!!!print*, dble(xsdb%p(plo)%T(1))
!!!print*, 'temp grid 1',dble(xsdb%p(plo)%T(:))
!!!print*, dble(xsdb%p(phi)%T(1))
!!!print*, 'temp grid 2',dble(xsdb%p(phi)%T(:))
!!!print*, 'TEST press', xspress
!!!print*,'pressure grids',dble(xsdb%p(:)%p)
        call writelog('GET_XS: Boundaries of XSDB LUT have been crossed!', 6)
        ExitXSFlag = 1
     endif

!!!print*, 'TEST temp', xstemp
!!!print*, dble(xsdb%p(plo)%T(:))
!!!print*, 'TEST press', xspress
!!!print*,'pressure grids',dble(xsdb%p(:)%p)
     !***Interpolate cross section in p and T (NR: bilinear) 
     if (ierr .ne. 0) then
        write(message, *) 'GET_XS: memory allocation error'
        ierr = ierr_all
        goto 999
     endif

     t = (xspress - dble(xsdb%p(plo)%p))/ dble(xsdb%p(phi)%p - xsdb%p(plo)%p)
     u_plo = (xstemp - dble(xsdb%p(plo)%T(Tlo_plo)))/ &
          dble(xsdb%p(plo)%T(Thi_plo)- xsdb%p(plo)%T(Tlo_plo))
     !      u_phi = (xstemp - dble(xsdb%p(phi)%T(Tlo_phi)))/ &             !HH: this is a bug!
     !               dble(xsdb%p(plo)%T(Thi_phi)- xsdb%p(plo)%T(Tlo_phi))   
     u_phi = (xstemp - dble(xsdb%p(phi)%T(Tlo_phi)))/ &
          dble(xsdb%p(phi)%T(Thi_phi)- xsdb%p(phi)%T(Tlo_phi))  
     u1_plo = 1. - u_plo
     u1_phi = 1. - u_phi
     t1 = 1. - t
!	write(*,*)'xspress',xspress
!	write(*,*)'xsdb%p(plo)%p',xsdb%p(plo)%p
!	write(*,*)'xsdb%p(phi)%p',xsdb%p(phi)%p
!	write(*,*)'xstemp',xstemp
!	write(*,*)'sdb%p(plo)%T(Tlo_plo)',xsdb%p(plo)%T(Tlo_plo)
!	write(*,*)'xsdb%p(plo)%T(Thi_plo)',xsdb%p(plo)%T(Thi_plo)
!	write(*,*)'t1',t1
!	write(*,*)'u1_plo',u1_plo
!	write(*,*)'xsdb%p(plo)%xs(:, 1,Tlo_plo)',xsdb%p(plo)%xs(:, 1,Tlo_plo)
!	write(*,*)'t',t
!	write(*,*)'u1_phi',u1_phi
!	write(*,*)'xsdb%p(phi)%xs(:, 1,Tlo_phi)',xsdb%p(phi)%xs(:, 1,Tlo_phi)
!	write(*,*)'t',t
!	write(*,*)'u_phi',u_phi
!	write(*,*)'xsdb%p(phi)%xs(:, 1,Thi_phi)',xsdb%p(phi)%xs(:, 1,Thi_phi)
!	write(*,*)'t1',t1
!	write(*,*)'u_plo',u_plo
!	write(*,*)'xsdb%p(plo)%xs(:, 1,Thi_plo)',xsdb%p(plo)%xs(:, 1,Thi_plo)
!          xsout(:) = t1*u1_plo*dble(xsdb%p(plo)%xs(:, 1,Tlo_plo))+ &
!               t*u1_phi*dble(xsdb%p(phi)%xs(:, 1,Tlo_phi))+ &
!               t*u_phi*dble(xsdb%p(phi)%xs(:, 1,Thi_phi))+ &
!               t1*u_plo*dble(xsdb%p(plo)%xs(:, 1,Thi_plo))# bilinear interpolation


     ! for cases with vmr >0.6 using values at the edge
     if (iup_vmr==1) then
        u_vmr=0.0
     elseif (ilo_vmr <xsdb%nvmr) then
        u_vmr = (vmr_k-xsdb%vmr(ilo_vmr))/(xsdb%vmr(iup_vmr)-xsdb%vmr(ilo_vmr))
     else
        u_vmr = 1.0
     endif

     !     xsdb%p(plo)%xs(:, ilo_vmr,Tlo_plo) + u_vmr*(xsdb%p(plo)%xs(:, iup_vmr,Tlo_plo) - xsdb%p(plo)%xs(:, ilo_vmr,Tlo_plo))
     !     xsdb%p(phi)%xs(:, ilo_vmr,Tlo_phi) + u_vmr*(xsdb%p(phi)%xs(:, iup_vmr,Tlo_phi) - xsdb%p(phi)%xs(:, ilo_vmr,Tlo_phi))
     !     xsdb%p(phi)%xs(:, ilo_vmr,Thi_phi) + u_vmr*(xsdb%p(phi)%xs(:, iup_vmr,Thi_phi) - xsdb%p(phi)%xs(:, ilo_vmr,Thi_phi))
     !     xsdb%p(plo)%xs(:, ilo_vmr,Thi_plo) + u_vmr*(xsdb%p(plo)%xs(:, iup_vmr,Thi_plo) - xsdb%p(plo)%xs(:, ilo_vmr,Thi_plo))

     xsout(:)=t1*u1_plo*dble( xsdb%p(plo)%xs(:, ilo_vmr,Tlo_plo) + u_vmr*(xsdb%p(plo)%xs(:, iup_vmr,Tlo_plo) - xsdb%p(plo)%xs(:, ilo_vmr,Tlo_plo)) )+ &
          t*u1_phi*dble( xsdb%p(phi)%xs(:, ilo_vmr,Tlo_phi) + u_vmr*(xsdb%p(phi)%xs(:, iup_vmr,Tlo_phi) - xsdb%p(phi)%xs(:, ilo_vmr,Tlo_phi)) )+ &
          t*u_phi*dble(  xsdb%p(phi)%xs(:, ilo_vmr,Thi_phi) + u_vmr*(xsdb%p(phi)%xs(:, iup_vmr,Thi_phi) - xsdb%p(phi)%xs(:, ilo_vmr,Thi_phi)) )+ &
          t1*u_plo*dble(xsdb%p(plo)%xs(:, ilo_vmr,Thi_plo) + u_vmr*(xsdb%p(plo)%xs(:, iup_vmr,Thi_plo) - xsdb%p(plo)%xs(:, ilo_vmr,Thi_plo)) )
     return
999  continue 
     call stopretrieval(message)

   end subroutine get_xs

!------------------------------------------------------------------------------
   subroutine check(status, ierr)
   use NETCDF 
   implicit none
   integer, intent (in) :: status
   integer, intent(out) :: ierr
  
      if(status /= nf90_noerr) then
         ierr = ierr_read
         call stopretrieval('READ_XSDB_NETCDF: '//trim(nf90_strerror(status)))
      else
         ierr = 0
      endif      
  
   end subroutine check
!------------------------------------------------------------------------------
! old Voigt-type database created by C. Frankenberg software package, based on HITRAN2008 input
   subroutine read_xsdb_frankenberg(xs_file, xsdb)
!*** Input
   character(len=*), intent(in) :: xs_file                   ! Path of the database file
!*** Input/output
   type(cross_section_db), intent(inout) :: xsdb           ! type containg cross-section database
!*** local variables
   integer(single) :: reclength
   integer :: i, k, l, io
   integer(single) :: xsdb_species, xsdb_isotope, xsdb_pnum, xsdb_dp, xsdb_Tnum, xsdb_dT, xsdb_length
   real(single) :: xsdb_plo, xsdb_Tlo
   real(single), dimension(:), allocatable :: xsdb_dum, xsdb_lambda
   !------------------------------------------------------

      xsdb%nvmr = 1
      if(allocated(xsdb%vmr)) deallocate(xsdb%vmr)
      allocate(xsdb%vmr(xsdb%nvmr))
      
!*** Read header and compute number of elements
      open(newunit(io), FILE=xs_file, action = 'read', status = 'old', ACCESS='DIRECT', RECL=9*SIZEOF(reclength))
!      write(*,*) 'Opened database file: ', trim(xs_file), ' (twice) in remotec_core/read_xsdb/read_xsdb_frankenberg (line 341 + 370)'
      read(io, REC=1) xsdb_species, xsdb_isotope, xsdb_plo, xsdb_pnum, &
                      xsdb_dp, xsdb_Tlo, xsdb_Tnum, xsdb_dT, xsdb_length
      close(io)

      xsdb%species = xsdb_species 
      xsdb%nwv = xsdb_length 
      allocate(xsdb_lambda(xsdb_length))
      allocate(xsdb_dum(xsdb_length*xsdb_pnum*xsdb_Tnum))

      if(allocated(xsdb%wv)) deallocate(xsdb%wv)
      if(allocated(xsdb%p)) deallocate(xsdb%p)
      allocate(xsdb%wv(xsdb_length))
      allocate(xsdb%p(xsdb_pnum))

      do i = 1, xsdb_pnum
         allocate(xsdb%p(i)%T(xsdb_Tnum))
      enddo
      do i = 0, xsdb_pnum-1
         xsdb%p(i+1)%p = xsdb_plo+i*xsdb_dp
      enddo
      do i = 0, xsdb_Tnum-1
         xsdb%p(1)%T(i+1) = xsdb_Tlo+i*xsdb_dT
      enddo
      do i = 2, xsdb_pnum
         xsdb%p(i)%T = xsdb%p(1)%T
      enddo
!      reclength = 9 + xsdb_length + xsdb_length*xsdb_pnum*xsdb_Tnum
	   
!*** Read xs database
      open(newunit(io), FILE=xs_file, action = 'read', status = 'old', ACCESS='DIRECT', &
           RECL=9*SIZEOF(reclength)+xsdb_length*SINGLE+xsdb_length*xsdb_pnum*xsdb_Tnum*SINGLE)
      read(io, REC=1) xsdb_species, xsdb_isotope, xsdb_plo, xsdb_pnum, &
                      xsdb_dp, xsdb_Tlo, xsdb_Tnum, xsdb_dT, xsdb_length,&
                      xsdb_lambda, xsdb_dum
      close(io)
      !*** Assume database works with wavenumber and RemoTeC works with wavelength
      xsdb%wv = 1.d7/dble(xsdb_lambda(xsdb_length:1:-1))
                   
      do i = 0, xsdb_pnum-1
!         allocate(xsdb%p(i+1)%xs(xsdb%nwv, xsdb_Tnum))
         allocate(xsdb%p(i+1)%xs(xsdb%nwv, xsdb%nvmr, xsdb_Tnum))
         do l = 0, xsdb_Tnum-1
            do k = 0, xsdb%nwv-1
                  ! Flip first dimension due to conversion from wavenumber to wavelength
                  xsdb%p(i+1)%xs(xsdb%nwv-k, 1, l+1) = xsdb_dum(i*xsdb_Tnum*xsdb_length+l*xsdb_length+1+k)
            enddo
         enddo
      enddo

      deallocate(xsdb_dum,xsdb_lambda)

   end subroutine read_xsdb_frankenberg
!------------------------------------------------------------------------------
!*** New line mixing database created by A. Butz based on the code by J.M. Hartmann, H. Tran
   subroutine read_xsdb_butz(xs_file, xsdb)
     !*** Input
     character(len=*), intent(in) :: xs_file                   ! Path of the database file
     !*** Input/output
     type(cross_section_db), intent(inout) :: xsdb           ! type containg cross-section database 
     !*** local variables
     integer :: i, k, l, io    
     integer(single) :: xsdb_pnum, xsdb_length
     real(double)::xsdb_lam0, xsdb_dlam
     real(double), dimension(:), allocatable :: xsdb_lambda
     integer, dimension(:), allocatable :: nT
     !---------------------------------------------------------------
     xsdb%nvmr = 1
     if(allocated(xsdb%vmr)) deallocate(xsdb%vmr)
     allocate(xsdb%vmr(xsdb%nvmr))

     !*** Start reading file
     open(newunit(io), FILE=xs_file, action = 'read', status = 'old', &
          ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
     !      write(*,*) 'Opened database file: ', trim(xs_file), ' in remotec_core/read_xsdb/read_xsdb_butz (ine 408)'

     !*** Process wavelength info
     read(io) xsdb_length, xsdb_lam0, xsdb_dlam
     allocate(xsdb_lambda(xsdb_length))
     do i = 1, xsdb_length
        xsdb_lambda(i) = xsdb_lam0 + (i-1)*xsdb_dlam
     enddo
     xsdb%nwv = xsdb_length

     !*** Assume database works with wavenumber and RemoTeC works with wavelength
     if (allocated(xsdb%wv)) deallocate(xsdb%wv)
     allocate(xsdb%wv(xsdb_length))
     xsdb%wv = 1.d7/xsdb_lambda(xsdb_length:1:-1)

     !*** Process pressure and temperature info
     read(io) xsdb_pnum            	      
     allocate(xsdb%p(xsdb_pnum),&
          nT(xsdb_pnum))

     do i = 1, xsdb_pnum
        read(io) xsdb%p(i)%p
     enddo

     do i = 1, xsdb_pnum
        read(io) nT(i)
     enddo

     do i = 1, xsdb_pnum	
        allocate(xsdb%p(i)%T(nT(i)))	  
        do k=1,nT(i)
           read(io)xsdb%p(i)%T(k)
        enddo
     enddo

     !*** Process xs database
     do i = 1, xsdb_pnum	         
        !         allocate(xsdb%p(i)%xs(xsdb%nwv, nT(i)))
        allocate(xsdb%p(i)%xs(xsdb%nwv, xsdb%nvmr, nT(i)))

        do k=1,nT(i)
           do l=1,xsdb_length
              !               read(io)xsdb%p(i)%xs(l, k) 
              read(io)xsdb%p(i)%xs(xsdb_length-l+1, 1, k) ! Flip first dimension due to conversion from wavenumber to wavelength
           enddo
        enddo
     enddo

     close(io)

     deallocate(xsdb_lambda)

   end subroutine read_xsdb_butz
!------------------------------------------------------------------------------

!*** ABSCO data from JPL
   subroutine read_xsdb_absco(xs_file, xsdb)
     !*** Input
     character(len=*), intent(in) :: xs_file                   ! Path of the database file
     !*** Input/output
     type(cross_section_db), intent(inout) :: xsdb           ! type containg cross-section database 
     !*** local variables
     integer :: i
     integer(hid_t) :: file_id, dset_id  ! JS: After moving the code to mistral, the dimension (hid_t) had to be added here in order for the code to compile. Possibly due to different hdf5 library
     integer :: hdferror                 !     The variable 'hdferror' is still without specified dimension though.
     integer(hsize_t), dimension(1:4) :: maxdims
     integer(hsize_t), dimension(1:4) :: dims

     double precision, dimension(:, :, :, :), allocatable  :: absorption
     double precision, dimension(:), allocatable  :: wavenumber, pressure
     double precision, dimension(:, :), allocatable  :: temperature

!!!!!!!!!!!!!!!1
     character(LEN=199) :: dataset
     character(2)   :: gas_index
     !*hdf
     integer(HSIZE_T),  dimension(:), allocatable :: dim
     INTEGER(HID_T)  :: space 

     !---------------------------------------------------------------
     xsdb%nvmr = 3
     if(allocated(xsdb%vmr)) deallocate(xsdb%vmr)
     allocate(xsdb%vmr(xsdb%nvmr))

     !*** Start reading file

     !*** Open NetCDF LUT
     call h5open_f(hdferror)
     call h5fopen_f(trim(xs_file), H5F_ACC_RDONLY_F, file_id, hdferror)

     dataset='Gas_Index'
     allocate(dim(1))
     dim=1
     call H5LTREAD_DATASET_CHAR_F(file_id, trim(dataset), gas_index, dim, hdferror) 
     deallocate(dim)


     ! get dimension size
     dataset='Gas_'//trim(gas_index)//'_Absorption'
     call h5dopen_f(file_id,trim(dataset), dset_id, hdferror)
     call h5dget_space_f(dset_id, space, hdferror)
     call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferror)
     call h5dclose_f(dset_id,hdferror)
     allocate(dim(4))
     dim=dims
     if(allocated(absorption))deallocate(absorption)
     allocate(absorption(dims(1), dims(2), dims(3), dims(4)))
     absorption=0.
     call h5dopen_f(file_id,trim(dataset), dset_id, hdferror)
     call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,absorption, dim, hdferror)
     call h5dclose_f(dset_id, hdferror)
     deallocate(dim)


     dataset='Wavenumber'
     allocate(dim(1))
     dim=(/dims(1)/)
     if(allocated(wavenumber))deallocate(wavenumber)
     allocate(wavenumber(dims(1)))
     wavenumber=0.
     call h5dopen_f(file_id,trim(dataset), dset_id, hdferror)
     call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,wavenumber, dim, hdferror)
     call h5dclose_f(dset_id, hdferror)
     deallocate(dim)


     dataset='Pressure'
     allocate(dim(1))
     dim=(/dims(4)/)
     if(allocated(pressure))deallocate(pressure)
     allocate(pressure(dims(4)))
     pressure=0.
     call h5dopen_f(file_id,trim(dataset), dset_id, hdferror)
     call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, pressure, dim, hdferror)
     call h5dclose_f(dset_id, hdferror)
     deallocate(dim)

     pressure = pressure*1.D-2

     dataset='Temperature'
     allocate(dim(2))
     dim=(/dims(3), dims(4)/)
     if(allocated(temperature))deallocate(temperature)
     allocate(temperature(dims(3),dims(4)))
     temperature=0.
     call h5dopen_f(file_id,trim(dataset), dset_id, hdferror)
     call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, temperature, dim, hdferror)
     call h5dclose_f(dset_id, hdferror)
     deallocate(dim)


     dataset='Broadener_01_VMR'
     allocate(dim(1))
     dim=(/dims(2)/)
     xsdb%vmr=0.
     call h5dopen_f(file_id,trim(dataset), dset_id, hdferror)
     call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,xsdb%vmr, dim, hdferror)
     call h5dclose_f(dset_id, hdferror)
     deallocate(dim)
     CALL H5Fclose_f (file_id, hdferror)

     !put data into xsbd structure
     xsdb%nwv = dims(1)
     if (allocated(xsdb%wv)) deallocate(xsdb%wv)
     allocate(xsdb%wv(xsdb%nwv))
     !*** Assume database works with wavenumber and RemoTeC works with wavelength
     xsdb%wv = 1.d7/wavenumber(xsdb%nwv:1:-1)

     !*** Process pressure and temperature info
     allocate(xsdb%p(dims(4)))

     do i = 1, dims(4)
        xsdb%p(i)%p = pressure(i)
     enddo

     do i = 1,dims(4)
        allocate(xsdb%p(i)%T(size(temperature(:,1))))	  
        xsdb%p(i)%T = temperature(:,i)
     enddo

     !*** Process xs database
     do i = 1,dims(4)
        allocate(xsdb%p(i)%xs(xsdb%nwv, xsdb%nvmr, size(temperature(:,1))))
        ! Flip first dimension due to conversion from wavenumber to wavelength
        xsdb%p(i)%xs = absorption(xsdb%nwv:1:-1,:,:,i)! nwave, nratio, ntemp
     enddo




   end subroutine read_xsdb_absco
!------------------------------------------------------------------------------
  SUBROUTINE H5LTREAD_DATASET_CHAR_F(file_id,dset_name,data_out,dims,error)
!
     INTEGER(HID_T) :: file_id       ! File identifier
     CHARACTER(len=*), INTENT(IN)		  :: dset_name
     INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN)   :: dims     ! size of the buffer	   
     CHARACTER(len=*),DIMENSION(*), INTENT(INOUT) :: data_out  ! data buffer	     
     INTEGER, INTENT(out)			  ::   error  ! Error flag
!
     INTEGER(HID_T)				  :: dset_id	   ! Dataset identifier
     INTEGER(HID_T)				  :: atype	   ! Datatype identifier
     INTEGER(SIZE_T)				  :: size	  ! Datatype size
!
     call h5tcopy_f(H5T_NATIVE_CHARACTER,atype,error)
     size=len(data_out)
     call h5tset_size_f(atype,size,error)
!
     call h5dopen_f(file_id,dset_name,dset_id,error)
!
     call h5dread_f(dset_id,atype,data_out,dims,error,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
     !   
     ! End access to the dataset and release resources used by it.
     !
     call h5dclose_f(dset_id,error)	  
!
  END SUBROUTINE h5ltread_dataset_char_f







end module read_xsdb_module
