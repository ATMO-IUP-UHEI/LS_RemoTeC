!------------------------------------------------------------------------------
!> @brief Compute optical properties of cirrus particles
!------------------------------------------------------------------------------
module optic_cirrus_module
   use header_module
   implicit none
   private
 
!*** procedures
   public :: optic_cirrus, optic_cirrus_ssc, optic_cirrus_xs, read_cirrus_netcdf, read_cirrus_ascii, &
             der_par_mode, der_par_mode_ssc
   private :: sizedis_cirrus, get_cirrus, get_cirrus_ssc, get_cirrus_xs

!*** types
   public :: cirrus_table

!------------------------------------------------------------------------------
!*** parameters
   integer, parameter :: nper_in = 1   ! number of scattering matrix elements
   integer, dimension(6), parameter :: index_ist = (/1,2,3,4,1,3 /)
   integer, dimension(6), parameter :: index_jst = (/1,2,3,4,2,4 /)

!------------------------------------------------------------------------------
!> Cirrus optical properties LUT
   type :: cirrus_table
      private
      character(4), dimension(:), allocatable :: lut_shape
      integer, dimension(:), allocatable :: lut_wv, lut_aax, lut_cax, lut_tilt
      real(double), dimension(:), allocatable :: cscasum_x, cextsum_x, cscasum_x_nodelta, cextsum_x_nodelta
      real(double), dimension(:, :, :), allocatable :: coefs_x  
      real(double), dimension(:), allocatable :: f ! truncation factor
      integer, dimension(:), allocatable :: ncoef
!*** needed for single scattering
      real(double), dimension(:, :, :), allocatable :: f_x
      real(double), dimension(:,:), allocatable :: scat_lut
   end type cirrus_table
!------------------------------------------------------------------------------

   contains
!------------------------------------------------------------------------------
!> Compute cirrus optical properties for multiple scattering
!> absorption and scattering optical thickness, phase function, and derivatives
!------------------------------------------------------------------------------  
 subroutine optic_cirrus(&
            cirrus_lut, &
            trunc_flag,&
            tilt_angle,&
            shapefrac,&
            reff_dum,&
            rlambda,&      
            csca,&
            cabs,&
            f, &
            plmom_aer_dum,&
            dcsca,&
            dcabs,&
            dcoefs,&
            nangle,&
            rnmb, &
            der_rnmb)
!*** input
   type(cirrus_table), intent(in) :: cirrus_lut
   logical, intent(in) :: trunc_flag      
   real(double), intent(in) :: tilt_angle, shapefrac, reff_dum, rlambda 
!*** output
   integer, intent(out) :: nangle 
   real(double), intent(out) :: csca, cabs, f 
   real(double),dimension(nstokes, nstokes, 0:maxleg), intent(out) :: plmom_aer_dum
   real(double),dimension(nfull), intent(out) :: dcsca, dcabs  
   real(double),dimension(dim_x), intent(out) ::  rnmb
   real(double),dimension(dim_x, 6), intent(out) :: der_rnmb
   real(double),dimension(nstokes, nstokes, 0:maxleg, nfull), intent(out) :: dcoefs
!*** local variables
   integer :: k, l, iper     
   real(double),dimension(dim_x+1) :: rmid
   real(double),dimension(dim_x) :: dr
   real(double),dimension(nper_in, 0:ndang) :: coefs 
   real(double),dimension(nper_in, 0:ndang, nfull) :: dcoefs_dum
   real(double),dimension(dim_x, 2) :: der_rnmb1

!------------------------------------------------------------------------------
      call sizedis_cirrus(reff_dum, rmid, dr, rnmb, der_rnmb1)          
      der_rnmb = 0.d0
      der_rnmb(:,1:2) = der_rnmb1
 
      call get_cirrus(&
           cirrus_lut, trunc_flag, rmid, dr, rnmb, rlambda, tilt_angle, shapefrac,&
           coefs, csca, cabs, f, nangle, dcoefs_dum, dcsca, dcabs)  

      dcoefs = 0.d0
      plmom_aer_dum = 0.d0  
      do iper = 1, nper_in
         if (index_ist(iper) .le. nstokes .and. index_jst(iper) .le. nstokes) then
            do k = 1, dim_x+2
               do l = 0, maxstr
                  dcoefs(index_ist(iper), index_jst(iper), l, k) = dcoefs_dum(iper, l, k)
               enddo
               do l = maxstr+1, min(maxleg,nangle)
                  dcoefs(index_ist(iper), index_jst(iper), l, k) = dcoefs_dum(iper, l, k)  
               enddo
            enddo
            do l = 0, maxstr
               plmom_aer_dum(index_ist(iper), index_jst(iper),l) = coefs(iper, l)
            enddo
            do l = maxstr+1,min(maxleg,nangle)
               plmom_aer_dum(index_ist(iper), index_jst(iper), l) = coefs(iper, l)
            enddo
         endif
      enddo
 
! note: to plot the scattering phase matrix, multiply the output "plmom_aer_dum" with the legendre polynom

   end subroutine optic_cirrus

!------------------------------------------------------------------------------
!> Compute cirrus optical properties for single scattering
!> absorption and scattering optical thickness, phase function, and derivatives
!------------------------------------------------------------------------------
   subroutine optic_cirrus_ssc(&
              cirrus_lut, &
              tilt_angle,&
              shapefrac,&
              reff_dum,&
              rlambda,&
              scat_angle,&       
              csca,&
              cabs,&
              z_ss,&
              dcsca,&
              dcabs,&
              der_z_ss,& 
              rnmb, &
              der_rnmb)
!*** input
   type(cirrus_table), intent(in) :: cirrus_lut
   real(double), intent(in) :: tilt_angle, shapefrac, scat_angle, reff_dum, rlambda
!*** output
   real(double), intent(out) :: csca, cabs
   real(double),dimension(nfull), intent(out) :: dcsca, dcabs 
   real(double), dimension(nstokes, nstokes), intent(out) :: z_ss
   real(double), dimension(nstokes, nstokes, nfull), intent(out) :: der_z_ss
   real(double),dimension(dim_x), intent(out) :: rnmb
   real(double),dimension(dim_x,6), intent(out) :: der_rnmb
!*** local variables     
   real(double),dimension(dim_x+1) :: rmid
   real(double),dimension(dim_x) :: dr
   real(double),dimension(dim_x,2) :: der_rnmb1 
!--------------------------------------------------------------------------------------     
      call sizedis_cirrus(reff_dum, rmid, dr, rnmb, der_rnmb1)           
      der_rnmb = 0.d0
      der_rnmb(:,1:2) = der_rnmb1 

      call get_cirrus_ssc(&
           cirrus_lut, rmid, dr, rnmb, rlambda, & 
           tilt_angle, shapefrac, scat_angle, &
	   z_ss, csca, cabs, der_z_ss, dcsca, dcabs)
  
   end subroutine optic_cirrus_ssc

!------------------------------------------------------------------------------
!> Compute cirrus optical properties for multiple scattering
!! only absorption and scattering optical thickness, no phase function, noderivatives
!------------------------------------------------------------------------------
   subroutine optic_cirrus_xs(&
              cirrus_lut, &
              trunc_flag,&
              tilt_angle,&
              shapefrac,&
              reff_dum,&
              rlambda,&       
              csca,&
              cabs, &
              f)  
!*** input
   type(cirrus_table), intent(in) :: cirrus_lut
   logical, intent(in) :: trunc_flag     
   real(double), intent(in) :: shapefrac, tilt_angle, reff_dum, rlambda   
!*** output
   real(double), intent(out) :: csca, cabs, f 
!*** local variables
   real(double),dimension(dim_x) :: rnmb
   real(double),dimension(dim_x+1) :: rmid
   real(double),dimension(dim_x) :: dr
   real(double),dimension(dim_x, 2) :: der_rnmb1
!----------------------------------------------------------------------
      call sizedis_cirrus(reff_dum, rmid, dr, rnmb, der_rnmb1)          

      call get_cirrus_xs(&
           cirrus_lut, trunc_flag, &
           rmid, dr, rnmb, rlambda, tilt_angle, shapefrac,&
	   csca, cabs, f) 
  
   end subroutine optic_cirrus_xs

!------------------------------------------------------------------------------
!> Get size distribution of cirrus particles given power par1 of power law
!------------------------------------------------------------------------------
   subroutine sizedis_cirrus(par1, clev, dc, cnmb, der_cnmb)
!*** input
  real(double) :: par1
!*** output
   real(double),dimension(dim_x+1), intent(out) :: clev    
   real(double),dimension(dim_x), intent(out) ::  dc, cnmb  
   real(double),dimension(dim_x,2), intent(out) ::  der_cnmb
!*** local variables
   integer :: i, l, n 
   real(double),dimension(2) :: der_sum  
   real(double) :: cup, clow, s, step
!------------------------------------------------------------------------------
      cup = 1300.d0
      clow = 3.d0
      step = (dlog(cup)-dlog(clow))/(dim_x+1)
      clev(1) = clow
      do l = 2, dim_x+1
         clev(l) = dexp(dlog(clev(1))+step*l)
         dc(l-1) = clev(l)-clev(l-1)
      enddo

!*** power law with power par1 given, large particles
      do i = 1, dim_x
         if(clev(i) .le. clow)then
	    cnmb(i) = 0.d0
            der_cnmb(i,:) = 0.d0
         else if(clev(i) .gt. cup)then
	    cnmb(i) = 0.d0
	    der_cnmb(i,:) = 0.d0
         else
	    cnmb(i) = (clev(i)/clow)**(-par1)
            der_cnmb(i,2) = 0.d0
	    der_cnmb(i,1) = -dlog(clev(i)/clow)*dexp(-par1*dlog(clev(i)/clow))
         endif
      enddo
   
      s = 0.
      do i = 1, dim_x
         s = s + dc(i)*cnmb(i)
      enddo

      do n = 1, 2
         der_sum(n) = 0.
         do i = 1, dim_x
            der_sum(n) = der_sum(n) + dc(i)*der_cnmb(i,n)
         enddo
         der_sum(n) = -s**(-2)*der_sum(n)
      enddo

      do n = 1, 2
         do i = 1, dim_x
           der_cnmb(i,n) = der_cnmb(i,n)/s + der_sum(n)*cnmb(i)
         enddo
      enddo

      do i = 1, dim_x
         cnmb(i) = cnmb(i)/s
      enddo

   end subroutine sizedis_cirrus
!------------------------------------------------------------------------------
!> Load NetCDF format cirrus look-up tables into memory
!------------------------------------------------------------------------------
   subroutine read_cirrus_netcdf(cirrus_file, tilt_angle, outputflag, cirrus_lut, ierr)
     use NETCDF
     !*** input
     integer, intent(in) :: outputflag
     character(stringlen), intent(in) :: cirrus_file   
     real(double), intent(in) :: tilt_angle
     !*** output
     type(cirrus_table), intent(out) :: cirrus_lut
     integer, intent(out) :: ierr
     !*** local variables
     integer :: i, l, k, ncid, nlut, cirrus_size(2), nexp, ndeg, nmat, nang
     integer :: varID, dimID
     integer, parameter :: nang_lut = 237  
     integer, dimension(1610) :: ncids
     character(1) :: shape
     real(double), dimension(:,:), allocatable :: exp_coef, scat_mat
     real(double), dimension(2) :: scat_coef, ext_coef
     character(stringlen) :: message

     if(outputflag >= 2) then
        write(message,'(a)') '*** Start read_cirrus_netcdf ***'
        call writelog(message, 1)
     endif

     !*** Open NetCDF LUT 
     call check(nf90_open(trim(cirrus_file), nf90_nowrite, ncid), ierr)
     if (ierr .ne. 0) return

     !*** Get #groups and their ncids
     call check(nf90_inq_grps(ncid, nlut, ncids), ierr)
     if (ierr .ne. 0) return
     if (nlut .gt. 1610) then
        write(message,*) 'READ_CIRRUS_NETCDF: number of groups is too large'
        ierr = ierr_var 
        goto 999
     endif

     allocate(cirrus_lut%lut_wv(nlut), &
          cirrus_lut%lut_aax(nlut), &
          cirrus_lut%lut_cax(nlut), &
          cirrus_lut%lut_tilt(nlut), &
          cirrus_lut%lut_shape(nlut), &
          cirrus_lut%cscasum_x(nlut), &
          cirrus_lut%cextsum_x(nlut), &
          cirrus_lut%cscasum_x_nodelta(nlut), &
          cirrus_lut%cextsum_x_nodelta(nlut), &
          cirrus_lut%coefs_x(nper_in, ndang, nlut), &
          cirrus_lut%ncoef(nlut), &
          cirrus_lut%f_x(nang_lut, nper_in, nlut), &
          cirrus_lut%scat_lut(nang_lut, nlut), &
          cirrus_lut%f(nlut), stat=ierr)
     if (ierr .ne. 0) then
        write(message, *) 'READ_CIRRUS_NETCDF: memory allocation error'
        ierr = ierr_all
        goto 999
     endif

     do i = 1, nlut
        !*** get wavelength
        call check(NF90_INQ_VARID(ncids(i), "wavelength", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, cirrus_lut%lut_wv(i)), ierr)
        if (ierr .ne. 0) return
        !*** get tilted angle
        call check(NF90_INQ_VARID(ncids(i), "tilt", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, cirrus_lut%lut_tilt(i)), ierr)
        if (ierr .ne. 0) return
        !*** get sizes
        call check(NF90_INQ_VARID(ncids(i), "size", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, cirrus_size), ierr)
        if (ierr .ne. 0) return
        cirrus_lut%lut_aax(i) = cirrus_size(1)
        cirrus_lut%lut_cax(i) = cirrus_size(2) 
        !** get shape
        call check(NF90_INQ_VARID(ncids(i), "shape", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, shape), ierr)
        if (ierr .ne. 0) return
        if (shape == "C") cirrus_lut%lut_shape(i) = "COLM"
        if (shape == "P") cirrus_lut%lut_shape(i) = "PLAT"
        !*** get scattering coefficients
        call check(NF90_INQ_VARID(ncids(i), "scat_coef", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, scat_coef), ierr)
        if (ierr .ne. 0) return
        cirrus_lut%cscasum_x(i) = 1.d-8*scat_coef(1)
        cirrus_lut%cscasum_x_nodelta(i) = 1.d-8*scat_coef(2)
        !*** get extincton coeffcients
        call check(NF90_INQ_VARID(ncids(i), "ext_coef", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, ext_coef), ierr)
        if (ierr .ne. 0) return
        cirrus_lut%cextsum_x(i) = 1.d-8*ext_coef(1)
        cirrus_lut%cextsum_x_nodelta(i) = 1.d-8*ext_coef(2)
        !*** get expansion coefficients
        call check(NF90_INQ_DIMID(ncids(i), "exp_coef_no", dimID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_INQUIRE_DIMENSION(ncids(i), dimID, len = nexp), ierr)
        if (ierr .ne. 0) return
        call check(NF90_INQ_DIMID(ncids(i), "spherical_degree", dimID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_INQUIRE_DIMENSION(ncids(i), dimID, len = ndeg), ierr)
        if (ierr .ne. 0) return
        if  (cirrus_lut%lut_shape(i)=='COLM') then
           cirrus_lut%ncoef(i) = min(ndeg, maxleg)
        elseif (cirrus_lut%lut_shape(i)=='PLAT') then
           cirrus_lut%ncoef(i) = min(ndeg, ndang)
        endif
        if (allocated(exp_coef)) deallocate(exp_coef)
        allocate(exp_coef(ndeg, nexp), stat=ierr)
        if (ierr .ne. 0) then
           write(message, *) 'READ_CIRRUS_NETCDF: memory allocation error'
           ierr = ierr_all
           goto 999
        endif
        call check(NF90_INQ_VARID(ncids(i), "exp_coef", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, exp_coef), ierr)
        if (ierr .ne. 0) return
        do k = 1, cirrus_lut%ncoef(i)
           do l = 1, nper_in
              cirrus_lut%coefs_x(l, k, i) = exp_coef(k, l)               
           enddo
        enddo
        !*** get scattering matrix elements
        call check(NF90_INQ_DIMID(ncids(i), "scat_mat_dim", dimID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_INQUIRE_DIMENSION(ncids(i), dimID, len = nmat), ierr)
        if (ierr .ne. 0) return
        call check(NF90_INQ_DIMID(ncids(i), "scat_angle", dimID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_INQUIRE_DIMENSION(ncids(i), dimID, len = nang), ierr)
        if (ierr .ne. 0) return
        if (nang .ne. nang_lut) then
           ierr = ierr_var
           call stopretrieval('READ_CIRRUS_NETCDF: nang .ne. nang_lut')
           return
        endif
        if (allocated(scat_mat)) deallocate(scat_mat)
        allocate(scat_mat(nang, nmat), stat=ierr)
        if (ierr .ne. 0) then
           write(message, *) 'READ_CIRRUS_NETCDF: memory allocation error'
           ierr = ierr_all
           goto 999
        endif
        call check(NF90_INQ_VARID(ncids(i), "scat_mat", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, scat_mat), ierr)
        if (ierr .ne. 0) return
        do k = 1, nang_lut
           do l = 1, nper_in
              cirrus_lut%f_x(k, l, i) = scat_mat(k, l)               
           enddo
        enddo
        do k = 1, nang_lut
           do l = 2, nper_in
              cirrus_lut%f_x(k, l, i) = cirrus_lut%f_x(k, l, i)*cirrus_lut%f_x(k, 1, i)
           enddo
           if (nper_in .ge. 5) cirrus_lut%f_x(k, 5, i) = -cirrus_lut%f_x(k, 5, i)
        enddo
        !*** get scattering angles
        call check(NF90_INQ_VARID(ncids(i), "scat_angle", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, cirrus_lut%scat_lut(:,i)), ierr)
        if (ierr .ne. 0) return
        !*** get truncation factor
        call check(NF90_INQ_VARID(ncids(i), "removed_energy", varID), ierr)
        if (ierr .ne. 0) return
        call check(NF90_GET_VAR(ncids(i), varID, cirrus_lut%f(i)), ierr) 
        if (ierr .ne. 0) return       
     enddo ! loop over nlut

     call check(nf90_close(ncid), ierr)

     if (outputflag >= 2) then
        write(message,'(a)')'*** End read_cirrus_netcdf ***'
        call writelog(message, 1)
     endif

     return

999  continue
     call stopretrieval(message)

   end subroutine read_cirrus_netcdf

!------------------------------------------------------------------------------
!> Check calls to netCDF functions
   subroutine check(status, ierr)
   use NETCDF 
   integer, intent (in) :: status
   integer, intent(out) :: ierr
  
      if(status /= nf90_noerr) then
         ierr = ierr_read
         call stopretrieval('READ_CIRRUS_NETCDF: '//trim(nf90_strerror(status)))
      else
         ierr = 0
      endif
      
  end subroutine check
!------------------------------------------------------------------------------
!> Load ascii format cirrus tables into memory
!> for one tilted angle, but all wavelengths and crystal sizes (aax and cax)
!------------------------------------------------------------------------------
   subroutine read_cirrus_ascii(cirrus_dir, tilt_angle, outputflag, cirrus_lut, ierr)
!*** input
   character(len=*), intent(in) :: cirrus_dir   
   real(double), intent(in) :: tilt_angle
   integer, intent(in) :: outputflag
!*** output
   integer, intent(out) :: ierr
   type(cirrus_table), intent(out) :: cirrus_lut
!*** local variables
   integer, parameter :: nang_lut = 237  
   integer :: i, j, k, iper, nlut, nlut_tot, io, iostat, ncoef_in, idummy, itilt, &
              lut_wv, lut_aax, lut_cax, lut_tilt, tilt_close
   character(4) :: lut_shape
   character(stringlen) :: fname, chdummy, fdummy
   real(double) :: cscasum, cextsum, cscasum_nodelta, cextsum_nodelta
   character(stringlen) :: message
!------------------------------------------------------------------------------
      ierr = 0

      if(outputflag >= 2) then
         write(message,'(a)') '*** Start read_cirrus_ascii ***'
         call writelog(message, 1)
      endif

      fname = trim(cirrus_dir)//'header.lut'
      open(newunit(io), file=trim(fname), action = 'read', status = 'old')
!      write(*,*) 'Opened cirrus ascii file: ', trim(fname), ' (twice) in remotec_core/optic_cirrus/read_cirrus_ascii (line 435 and 448))'
      i = 0
      do 
         read(io, *, iostat = iostat)
         if(iostat<0) exit
	 i = i + 1
      enddo
      close(io)
      nlut_tot = i
      if(allocated(cirrus_lut%lut_tilt)) deallocate(cirrus_lut%lut_tilt)
      allocate(cirrus_lut%lut_tilt(nlut_tot), stat=ierr)

      open(newunit(io),file=trim(fname), action = 'read', status = 'old')
      do k = 1, nlut_tot
         read(io,*) lut_shape, idummy, idummy, idummy, cirrus_lut%lut_tilt(k)
         if ( k>1 ) then
            if  (cirrus_lut%lut_tilt(k) < cirrus_lut%lut_tilt(k-1)) then
               i = k-1
               exit
            endif
         endif
      enddo
      close(io) 

!*** only one tilted angle
      nlut = nlut_tot/i         
      itilt = minval(minloc(dabs(cirrus_lut%lut_tilt(1:nlut)-tilt_angle)))
      tilt_close = cirrus_lut%lut_tilt(itilt)

      deallocate(cirrus_lut%lut_tilt, stat=ierr)
      allocate(cirrus_lut%lut_wv(nlut), &
               cirrus_lut%lut_aax(nlut), &
               cirrus_lut%lut_cax(nlut), &
               cirrus_lut%lut_tilt(nlut), &
               cirrus_lut%lut_shape(nlut), &
               cirrus_lut%cscasum_x(nlut), &
               cirrus_lut%cextsum_x(nlut), &
               cirrus_lut%cscasum_x_nodelta(nlut), &
               cirrus_lut%cextsum_x_nodelta(nlut), &
               cirrus_lut%coefs_x(nper_in, ndang, nlut), &
               cirrus_lut%ncoef(nlut), &
               cirrus_lut%f_x(nang_lut, nper_in, nlut), &
               cirrus_lut%scat_lut(nang_lut, nlut), &
               cirrus_lut%f(nlut), stat=ierr)
      open(newunit(io),file=trim(fname), action = 'read', status = 'old')
!      write(*,*) 'Opened cirrus ascii file: ', trim(fname), ' (third?) in remotec_core/optic_cirrus/read_cirrus_ascii (line 480)'
      i = 0
      do k = 1, nlut_tot
         read(io,*) lut_shape, lut_wv, lut_aax, lut_cax, lut_tilt 
         if (lut_tilt == tilt_close) then 
            i = i + 1
            cirrus_lut%lut_shape(i) = lut_shape
            cirrus_lut%lut_wv(i) = lut_wv
            cirrus_lut%lut_aax(i) = lut_aax
            cirrus_lut%lut_cax(i) = lut_cax
            cirrus_lut%lut_tilt(i) = lut_tilt  
         endif
      enddo
      close(io)
      if (i .ne. nlut) then
         call stopretrieval('READ_CIRRUS_ASCII: incorrect number of LUTs')
      endif

!*** loop over cirrus luts
      do i = 1, nlut
!*** read *.dak file
         write(fdummy,'(a,a,i4.4,a,i5.5,a,i5.5,a,i2.2,a)') &
               cirrus_lut%lut_shape(i), &
               '_wv', cirrus_lut%lut_wv(i), &
               '_a',cirrus_lut%lut_aax(i), &
               '_c',cirrus_lut%lut_cax(i), &
               '_t',cirrus_lut%lut_tilt(i),'.dak'
         open(newunit(io), file=trim(cirrus_dir)//trim(fdummy), action = 'read', status = 'old')
!         write(*,*) 'Opened dummy file: ', trim(cirrus_dir)//trim(fdummy), ' in remotec_core/optic_cirrus/read_cirrus_ascii (line 508)'
         read(io,'(a79,i4)') chdummy, ncoef_in
         ncoef_in = ncoef_in + 1

         do k = 1, 4
            read(io,*)
         enddo
         read(io, '(a37, E21.14)') chdummy, cirrus_lut%f(i)
         do k = 1, 3
            read(io,*)
         enddo       

!*** if input data are to be translated into number of cirri: 
!*** calipso ot needs to be scaled using extinction cross-sections without delta-approximation
            read(io,'(a31,f19.10,a21,f19.10)') chdummy, cextsum, chdummy, cextsum_nodelta
            read(io,'(a31,f19.10,a21,f19.10)') chdummy, cscasum, chdummy, cscasum_nodelta
	 read(io,*)
         if  (cirrus_lut%lut_shape(i)=='COLM') then
            cirrus_lut%ncoef(i) = min(ncoef_in, maxleg)
         elseif (cirrus_lut%lut_shape(i)=='PLAT') then
            cirrus_lut%ncoef(i) = min(ncoef_in, ndang)
         endif
         cirrus_lut%cscasum_x(i) = cscasum*1.d-8
         cirrus_lut%cextsum_x(i) = cextsum*1.d-8
         cirrus_lut%cscasum_x_nodelta(i) = cscasum_nodelta*1.d-8
         cirrus_lut%cextsum_x_nodelta(i) = cextsum_nodelta*1.d-8
         do j = 1, cirrus_lut%ncoef(i)
            read(io,*) idummy, (cirrus_lut%coefs_x(k, j, i), k=1,nper_in)
	 enddo         
         close(io)

!*** read *.out file
         write(fdummy,'(a,a,i4.4,a,i5.5,a,i5.5,a,i2.2,a)') &
               cirrus_lut%lut_shape(i), &
               '_wv', cirrus_lut%lut_wv(i), &
               '_a',cirrus_lut%lut_aax(i), &
               '_c',cirrus_lut%lut_cax(i), &
               '_t',cirrus_lut%lut_tilt(i),'.out'
         open(newunit(io), file=trim(cirrus_dir)//trim(fdummy), action = 'read', status = 'old')
!         write(*,*) 'Opened dummy file: ', trim(cirrus_dir)//trim(fdummy), ' (twice?) in remotec_core/optic_cirrus/read_cirrus_ascii (line 547)'
         do k = 1, 15
            read(io,*)
         enddo
         do j = 1, nang_lut
            read(io,*) cirrus_lut%scat_lut(j, i), (cirrus_lut%f_x(j, k, i), k=1,nper_in ) 
         enddo
         close(io)
         
         do j = 1, nang_lut
           do iper = 2, nper_in
              cirrus_lut%f_x(j, iper, i) = cirrus_lut%f_x(j, iper, i)*cirrus_lut%f_x(j, 1, i)
           enddo
           if (nper_in .ge. 5) cirrus_lut%f_x(j, 5, i) = - cirrus_lut%f_x(j, 5, i)
         enddo
      enddo ! i = 1, nlut

      if (outputflag >= 2) then
         write(message,'(a)')'*** End read_cirrus_ascii ***'
         call writelog(message, 1)
      endif 

   end subroutine read_cirrus_ascii
!------------------------------------------------------------------------------
!> Get optical properties of cirrus for multiple scattering from LUT
!------------------------------------------------------------------------------
   subroutine get_cirrus( &
              cirrus_lut, trunc_flag, rmid, dr, rnmb, rlambda, tilt_angle, shapefrac,&          
              coefs, csca, cabs, f, nangle, der_coefs, der_csca, der_cabs)   
!*** input
   type(cirrus_table), intent(in) :: cirrus_lut
   logical, intent(in) :: trunc_flag
   real(double),dimension(dim_x+1), intent(in) :: rmid
   real(double),dimension(dim_x), intent(in) :: dr, rnmb
   real(double), intent(in) :: rlambda, tilt_angle, shapefrac    
!*** output
   integer, intent(out) :: nangle
   real(double), intent(out) ::  csca, cabs, f
   real(double),dimension(nper_in, 0:ndang), intent(out) :: coefs
   real(double),dimension(nper_in, 0:ndang, nfull), intent(out) :: der_coefs
   real(double),dimension(nfull), intent(out) :: der_csca, der_cabs
!*** local variables
   integer :: iwv, itilt
   integer :: i, j, k, l, nlut, icax, iaax
   integer :: colm_ncoef, plat_ncoef   
   real(double) :: cext
   real(double) :: colm_csca, colm_cext
   real(double) :: plat_csca,  plat_cext
   real(double),dimension(nper_in, 0:ndang) :: colm_coefs
   real(double),dimension(nper_in, 0:ndang) :: plat_coefs   
   real(double),dimension(nper_in, 0:ndang, nfull) :: colm_der_coefs
   real(double),dimension(nper_in, 0:ndang, nfull) :: plat_der_coefs   
   real(double),dimension(nfull) :: der_cext
   real(double),dimension(nfull) :: colm_der_csca, colm_der_cext
   real(double),dimension(nfull) :: plat_der_csca, plat_der_cext
   real(double) :: colm_cscasum_x, colm_cextsum_x
   real(double) :: plat_cscasum_x, plat_cextsum_x   
   real(double),dimension(nper_in, ndang) :: colm_coefs_x
   real(double),dimension(nper_in, ndang) :: plat_coefs_x   
   real(double) :: cax_min, aax_min, sw, sw_der  
!------------------------------------------------------------
!*** look for wavelength and tilt angle in lut
      iwv = minval(minloc(dabs(cirrus_lut%lut_wv/1000.d0-rlambda)))
      itilt = minval(minloc(dabs(cirrus_lut%lut_tilt-tilt_angle)))
      nlut = size(cirrus_lut%lut_wv)
!*** columnar ice particles
      colm_coefs = 0.d0 	     
      colm_der_csca = 0.d0
      colm_der_cext = 0.d0
      colm_der_coefs = 0.d0	     
      colm_csca = 0.d0
      colm_cext = 0.d0
      colm_ncoef = ndang
!*** if the shape factor requires columnar particles then calculate them: 
      if ((1.0 - shapefrac) > 1d-3) then
         do l = 1, dim_x      
!*** look for c dimension fitting particle radius in lut
            cax_min = 9.d99
            do i = 1, nlut
               if(cirrus_lut%lut_shape(i)=='COLM' .and. dabs(cirrus_lut%lut_cax(i)-rmid(l)*10)<cax_min) then
	          icax = i
	          cax_min = dabs(cirrus_lut%lut_cax(i)-rmid(l)*10)
	       endif
            enddo
            iaax = icax
            do i = 1, nlut
               if (cirrus_lut%lut_shape(i)=='COLM' .and. &
                   cirrus_lut%lut_wv(i) == cirrus_lut%lut_wv(iwv) .and. &
                   cirrus_lut%lut_aax(i) == cirrus_lut% lut_aax(iaax) .and. &
                   cirrus_lut% lut_cax(i) == cirrus_lut%lut_cax(icax) .and. &
                   cirrus_lut%lut_tilt(i) == cirrus_lut%lut_tilt(itilt) ) then
!*** if input data are to be translated into number of cirri, calipso ot needs 
!*** to be scaled using extinction cross-sections without delta-approximation
                  if(trunc_flag .eqv. .false.)then
                     colm_cextsum_x = cirrus_lut%cextsum_x_nodelta(i)
                     colm_cscasum_x = cirrus_lut%cscasum_x_nodelta(i)
                  else
                     colm_cextsum_x = cirrus_lut%cextsum_x(i)
                     colm_cscasum_x = cirrus_lut%cscasum_x(i)
                  endif
                  colm_ncoef = cirrus_lut%ncoef(i)
                  do j = 1, colm_ncoef
                     do k = 1, nper_in
                        colm_coefs_x(k,j) = cirrus_lut%coefs_x(k,j,i)
                     enddo
                  enddo
                  f = cirrus_lut%f(i)
                  exit
               endif 
            enddo ! i = 1, nlut
            sw = rnmb(l)*dr(l)
            sw_der = dr(l)
      
            colm_csca = colm_csca + sw*colm_cscasum_x
            colm_cext = colm_cext + sw*colm_cextsum_x
            colm_der_csca(l) = sw_der*colm_cscasum_x
            colm_der_cext(l) = sw_der*colm_cextsum_x
      
            do j = 0, colm_ncoef-1
               do i = 1, nper_in
                  colm_coefs(i,j) = colm_coefs(i,j) + sw*colm_coefs_x(i,j+1)
                  colm_der_coefs(i,j,l) = sw_der*colm_coefs_x(i,j+1)
               enddo
            enddo

         enddo ! l = 1, dim_x 
      endif ! columns

!*** plate ice particles
      plat_coefs = 0.d0  	
      plat_der_csca = 0.d0
      plat_der_cext = 0.d0
      plat_der_coefs = 0.d0	
      plat_csca = 0.d0
      plat_cext = 0.d0
      plat_ncoef = ndang
!*** if the shape factor requires plate particles then calculate them: 
      if (shapefrac > 1d-3) then
         do l = 1, dim_x      
!*** look for 2*a dimension fitting particle radius in lut
            aax_min = 9.d99
            do i = 1, nlut
	       if(cirrus_lut%lut_shape(i)=='PLAT' .and. dabs(2.*cirrus_lut%lut_aax(i)-rmid(l)*10)<aax_min)then
	          iaax = i
	          aax_min = dabs(cirrus_lut%lut_aax(i)-rmid(l)*10)
	       endif
            enddo
            icax = iaax
            do i = 1, nlut
               if (cirrus_lut%lut_shape(i)=='PLAT' .and. &
                   cirrus_lut%lut_wv(i) == cirrus_lut%lut_wv(iwv) .and. &
                   cirrus_lut%lut_aax(i) == cirrus_lut% lut_aax(iaax) .and. &
                   cirrus_lut% lut_cax(i) == cirrus_lut%lut_cax(icax) .and. &
                   cirrus_lut%lut_tilt(i) == cirrus_lut%lut_tilt(itilt) ) then
!*** if input data are to be translated into number of cirri, then calipso ot 
!*** needs to be scaled using extinction cross-sections without delta-approximation
                   if(trunc_flag .eqv. .false.)then
                      plat_cextsum_x = cirrus_lut%cextsum_x_nodelta(i)
	              plat_cscasum_x = cirrus_lut%cscasum_x_nodelta(i)
                   else
                      plat_cextsum_x = cirrus_lut%cextsum_x(i)
                      plat_cscasum_x = cirrus_lut%cscasum_x(i)
                   endif
                   plat_ncoef = cirrus_lut%ncoef(i)
                   do j = 1, plat_ncoef
                      do k = 1, nper_in
                          plat_coefs_x(k,j) = cirrus_lut%coefs_x(k,j,i)
                      enddo
                   enddo
                   f = cirrus_lut%f(i)
                   exit
                endif
            enddo ! i = 1, nlut
            sw = rnmb(l)*dr(l)
            sw_der = dr(l)

            plat_csca = plat_csca+sw*plat_cscasum_x
            plat_cext = plat_cext+sw*plat_cextsum_x
            plat_der_csca(l) = sw_der*plat_cscasum_x
            plat_der_cext(l) = sw_der*plat_cextsum_x
      
            do j = 0, plat_ncoef-1
               do i = 1, nper_in
                  plat_coefs(i,j) = plat_coefs(i,j) + sw*plat_coefs_x(i,j+1)
                  plat_der_coefs(i,j,l) = sw_der*plat_coefs_x(i,j+1)
               enddo
            enddo
         enddo ! l = 1, dim_x 
      endif ! plates

      csca = shapefrac*plat_csca + (1.0-shapefrac)*colm_csca
      cext = shapefrac*plat_cext + (1.0-shapefrac)*colm_cext
      cabs = cext - csca

      der_csca = shapefrac*plat_der_csca - (1.0-shapefrac)*colm_der_csca
      der_cext = shapefrac*plat_der_cext - (1.0-shapefrac)*colm_der_cext
      der_cabs = der_cext - der_csca

      coefs(:,:) = (plat_csca*shapefrac*plat_coefs(:,:)+colm_csca*(1.0-shapefrac)*colm_coefs(:,:))/csca
      do l = 1, dim_x
         der_coefs(:,:,l) = coefs(:,:)/csca*der_csca(l)-&
                            (shapefrac*plat_der_csca(l)*plat_coefs(:,:) + &
                            shapefrac*plat_csca*plat_der_coefs(:,:,l) + &
                            (1.0-shapefrac)*colm_der_csca(l)*colm_coefs(:,:) + &
                            (1.0-shapefrac)*colm_csca*colm_der_coefs(:,:,l))/csca
      enddo
      der_coefs(:,:,dim_x+1) = 0.d0
      der_coefs(:,:,dim_x+2) = 0.d0 
      nangle = min(plat_ncoef, colm_ncoef, ndang)

   end subroutine get_cirrus

!------------------------------------------------------------------------------
!> Get optical properties of cirrus 
!! only absorption/scattering optical thickness and truncation factor
!! no phase function, no derivatives
!------------------------------------------------------------------------------
   subroutine get_cirrus_xs( &
              cirrus_lut, trunc_flag, &
              rmid, dr, rnmb, rlambda, &
              tilt_angle, shapefrac,&          
              csca, cabs, f)
   
!*** input
   type(cirrus_table), intent(in) :: cirrus_lut
   logical, intent(in) :: trunc_flag
   real(double),dimension(dim_x+1), intent(in) :: rmid
   real(double),dimension(dim_x), intent(in) :: dr, rnmb
   real(double), intent(in) :: rlambda, tilt_angle, shapefrac    
!*** output
   real(double), intent(out) ::  csca, cabs, f
!*** local variables
   integer :: iwv, itilt, icax, iaax
   integer :: i,l, nlut 
   real(double) :: cext
   real(double) :: colm_csca, colm_cext
   real(double) :: plat_csca,  plat_cext    
   real(double) :: colm_cscasum_x, colm_cextsum_x
   real(double) :: plat_cscasum_x, plat_cextsum_x      
   real(double) :: cax_min, aax_min, sw !, sw_der  
!------------------------------------------------------------
!*** look for wavelength and tilt angle in lut
      iwv = minval(minloc(dabs(cirrus_lut%lut_wv/1000.d0-rlambda)))
      itilt = minval(minloc(dabs(cirrus_lut%lut_tilt-tilt_angle)))
      nlut = size(cirrus_lut%lut_wv)
!*** columnar ice particles	     	     
      colm_csca = 0.d0
      colm_cext = 0.d0

!*** if the shape factor requires columnar particles then calculate them: 
      if ((1.0 - shapefrac) > 1d-3) then
         do l = 1, dim_x      
!*** look for c dimension fitting particle radius in lut
            cax_min = 9.d99
            do i = 1, nlut
               if(cirrus_lut%lut_shape(i)=='COLM' .and. dabs(cirrus_lut%lut_cax(i)-rmid(l)*10)<cax_min) then
	          icax = i
	          cax_min = dabs(cirrus_lut%lut_cax(i)-rmid(l)*10)
	       endif
            enddo
            iaax = icax
            do i = 1, nlut
               if (cirrus_lut%lut_shape(i)=='COLM' .and. &
                   cirrus_lut%lut_wv(i) == cirrus_lut%lut_wv(iwv) .and. &
                   cirrus_lut%lut_aax(i) == cirrus_lut% lut_aax(iaax) .and. &
                   cirrus_lut% lut_cax(i) == cirrus_lut%lut_cax(icax) .and. &
                   cirrus_lut%lut_tilt(i) == cirrus_lut%lut_tilt(itilt) ) then
!*** if input data are to be translated into number of cirri, calipso ot needs 
!*** to be scaled using extinction cross-sections without delta-approximation
                  if(trunc_flag .eqv. .false.)then
                     colm_cextsum_x = cirrus_lut%cextsum_x_nodelta(i)
                     colm_cscasum_x = cirrus_lut%cscasum_x_nodelta(i)
                  else
                     colm_cextsum_x = cirrus_lut%cextsum_x(i)
                     colm_cscasum_x = cirrus_lut%cscasum_x(i)
                  endif
                  f = cirrus_lut%f(i)
                  exit
               endif 
            enddo ! i = 1, nlut
            sw = rnmb(l)*dr(l)
      
            colm_csca = colm_csca + sw*colm_cscasum_x
            colm_cext = colm_cext + sw*colm_cextsum_x

         enddo ! l = 1, dim_x 
      endif ! columns

!*** plate ice particles 		
      plat_csca = 0.d0
      plat_cext = 0.d0
!*** if the shape factor requires plate particles then calculate them: 
      if (shapefrac > 1d-3) then
         do l = 1, dim_x      
!*** look for 2*a dimension fitting particle radius in lut
            aax_min = 9.d99
            do i = 1, nlut
	       if(cirrus_lut%lut_shape(i)=='PLAT' .and. dabs(2.*cirrus_lut%lut_aax(i)-rmid(l)*10)<aax_min)then
	          iaax = i
	          aax_min = dabs(cirrus_lut%lut_aax(i)-rmid(l)*10)
	       endif
            enddo
            icax = iaax
            do i = 1, nlut
               if (cirrus_lut%lut_shape(i)=='PLAT' .and. &
                   cirrus_lut%lut_wv(i) == cirrus_lut%lut_wv(iwv) .and. &
                   cirrus_lut%lut_aax(i) == cirrus_lut% lut_aax(iaax) .and. &
                   cirrus_lut% lut_cax(i) == cirrus_lut%lut_cax(icax) .and. &
                   cirrus_lut%lut_tilt(i) == cirrus_lut%lut_tilt(itilt) ) then
!*** if input data are to be translated into number of cirri, then calipso ot 
!*** needs to be scaled using extinction cross-sections without delta-approximation
                   if(trunc_flag .eqv. .false.)then
                      plat_cextsum_x = cirrus_lut%cextsum_x_nodelta(i)
	              plat_cscasum_x = cirrus_lut%cscasum_x_nodelta(i)
                   else
                      plat_cextsum_x = cirrus_lut%cextsum_x(i)
                      plat_cscasum_x = cirrus_lut%cscasum_x(i)
                   endif
                   f = cirrus_lut%f(i)
                   exit
                endif
            enddo ! i = 1, nlut
            sw = rnmb(l)*dr(l)

            plat_csca = plat_csca + sw*plat_cscasum_x
            plat_cext = plat_cext + sw*plat_cextsum_x
 
         enddo ! l = 1, dim_x 
      endif ! plates

      csca = shapefrac*plat_csca + (1.0-shapefrac)*colm_csca
      cext = shapefrac*plat_cext + (1.0-shapefrac)*colm_cext
      cabs = cext - csca

      return

   end subroutine get_cirrus_xs
!------------------------------------------------------------------------------
!> Get optical properties of cirrus for single scattering from LUT
!------------------------------------------------------------------------------
   subroutine get_cirrus_ssc(&
              cirrus_lut, &
              rmid, dr, rnmb, rlambda, tilt_angle, shapefrac, scat_angle,&          
              z_ss, csca, cabs, der_z_ss, der_csca, der_cabs)
   
!*** input
   type(cirrus_table), intent(in) :: cirrus_lut
   real(double),dimension(dim_x+1), intent(in) :: rmid
   real(double),dimension(dim_x), intent(in) :: dr, rnmb
   real(double), intent(in) :: rlambda, tilt_angle, shapefrac, scat_angle    
!*** output
   real(double), dimension(nstokes, nstokes) :: z_ss
   real(double), intent(out) ::  csca, cabs
   real(double), dimension(nstokes, nstokes, nfull), intent(out) :: der_z_ss 
   real(double),dimension(nfull), intent(out) :: der_csca, der_cabs
!*** local variables
   integer, parameter :: nang_lut = 237     
   integer :: iwv, itilt
   integer :: i, j, k, l, nlut, icax, iaax, imin
   real(double) :: cext
   real(double) :: colm_csca, colm_cext
   real(double) :: plat_csca,  plat_cext 
   real(double),dimension(nfull) :: der_cext
   real(double),dimension(nfull) :: colm_der_csca, colm_der_cext
   real(double),dimension(nfull) :: plat_der_csca, plat_der_cext
   real(double) :: colm_cscasum_x, colm_cextsum_x
   real(double) :: plat_cscasum_x, plat_cextsum_x     
   real(double) :: cax_min, aax_min, sw, sw_der    
   real(double), dimension(nang_lut) :: scat_lut 
   real(double),dimension(nang_lut, nper_in) :: f_ss_full
   real(double),dimension(nang_lut, nper_in, nfull) :: der_fss_full
   real(double),dimension(nang_lut, nper_in) :: colm_f
   real(double),dimension(nang_lut, nper_in) :: plat_f
   real(double),dimension(nang_lut, nper_in, nfull) :: colm_der_f
   real(double),dimension(nang_lut, nper_in, nfull) :: plat_der_f
   real(double), dimension(nper_in) :: f_ss
   real(double), dimension(nper_in, nfull) :: der_f_ss
   real(double), dimension(nang_lut, nper_in) :: colm_f_x, plat_f_x   
!-------------------------------------------------------------------
!*** look for wavelength and tilt angle in lut.
      iwv = minval(minloc(dabs(cirrus_lut%lut_wv/1000.d0-rlambda)))
      itilt = minval(minloc(dabs(cirrus_lut%lut_tilt-tilt_angle)))
      nlut = size(cirrus_lut%lut_wv)

!*** columnar ice particles
      colm_f = 0.d0 	     
      colm_der_csca = 0.d0
      colm_der_cext = 0.d0
      colm_der_f = 0.d0	     
      colm_csca = 0.d0
      colm_cext = 0.d0

!*** if the shape factor requires columnar particles then calculate them: 
      if ((1.0-shapefrac) > 1d-3) then
         do l = 1, dim_x
!*** look for c dimension fitting particle radius in lut
            cax_min = 9.d99
            do i = 1, nlut
               if(cirrus_lut%lut_shape(i)=='COLM' .and. dabs(cirrus_lut%lut_cax(i)-rmid(l)*10)<cax_min) then
	          icax = i
	          cax_min = dabs(cirrus_lut%lut_cax(i)-rmid(l)*10) 
               endif
            enddo
            iaax = icax
            do i = 1, nlut
               if (cirrus_lut%lut_shape(i)=='COLM' .and. &
                   cirrus_lut%lut_wv(i) == cirrus_lut%lut_wv(iwv) .and. &
                   cirrus_lut%lut_aax(i) == cirrus_lut% lut_aax(iaax) .and. &
                   cirrus_lut% lut_cax(i) == cirrus_lut%lut_cax(icax) .and. &
                   cirrus_lut%lut_tilt(i) == cirrus_lut%lut_tilt(itilt) ) then
                  colm_cextsum_x = cirrus_lut%cextsum_x(i)
                  colm_cscasum_x = cirrus_lut%cscasum_x(i)

                  do j = 1, nang_lut                 
                     scat_lut(j) = cirrus_lut%scat_lut(j,i)/180.*pi
                     do k = 1, nper_in                  
                        colm_f_x(j, k) = cirrus_lut%f_x(j, k, i)
                     enddo
                  enddo

                  exit
               endif 
            enddo ! i = 1, nlut
            sw = rnmb(l)*dr(l)
            sw_der = dr(l)
      
            colm_csca = colm_csca + sw*colm_cscasum_x
            colm_cext = colm_cext + sw*colm_cextsum_x
            colm_der_csca(l) = sw_der*colm_cscasum_x
            colm_der_cext(l) = sw_der*colm_cextsum_x

           do i = 1, nper_in
              do j = 1, nang_lut
!                  colm_f(j, i) = colm_f(j, 1) + sw*colm_f_x(j, i)   !HH: I think it should be this
                  colm_f(j, i) = colm_f(j, i) + sw*colm_f_x(j, i)
                  colm_der_f(j, i, l) = sw_der*colm_f_x(j, i)
               enddo
            enddo

         enddo ! l = 1, dim_x 
      endif ! columns

!*** plate ice particles
      plat_f = 0.d0  	
      plat_der_csca = 0.d0
      plat_der_cext = 0.d0
      plat_der_f = 0.d0	
      plat_csca = 0.d0
      plat_cext = 0.d0
 
!*** if the shape factor requires plate particles then calculate them: 
      if (shapefrac > 1d-3) then
         do l = 1, dim_x
      
!*** look for 2*a dimension fitting particle radius in lut
            aax_min = 9.d99
            do i = 1, nlut
               if(cirrus_lut%lut_shape(i)=='PLAT' .and. dabs(2.*cirrus_lut%lut_aax(i)-rmid(l)*10)<aax_min)then
	          iaax = i
	          aax_min = dabs(cirrus_lut%lut_aax(i)-rmid(l)*10)
               endif
            enddo
            icax = iaax
            do i = 1, nlut
               if (cirrus_lut%lut_shape(i)=='PLAT' .and. &
                   cirrus_lut%lut_wv(i) == cirrus_lut%lut_wv(iwv) .and. &
                   cirrus_lut%lut_aax(i) == cirrus_lut% lut_aax(iaax) .and. &
                   cirrus_lut% lut_cax(i) == cirrus_lut%lut_cax(icax) .and. &
                   cirrus_lut%lut_tilt(i) == cirrus_lut%lut_tilt(itilt) ) then
                  plat_cextsum_x = cirrus_lut%cextsum_x(i)
                  plat_cscasum_x = cirrus_lut%cscasum_x(i)

                  do j = 1, nang_lut             
                     scat_lut(j) = cirrus_lut%scat_lut(j,i)/180.*pi
                     do k = 1, nper_in                  
                        plat_f_x(j, k) = cirrus_lut%f_x(j, k, i)
                     enddo
                  enddo
                  exit

               endif
            enddo ! i = 1, nlut
            sw = rnmb(l)*dr(l)
            sw_der = dr(l)

            plat_csca = plat_csca + sw*plat_cscasum_x
            plat_cext = plat_cext + sw*plat_cextsum_x
            plat_der_csca(l) = sw_der*plat_cscasum_x
            plat_der_cext(l) = sw_der*plat_cextsum_x

           do i = 1, nper_in
              do j = 1, nang_lut
!                  plat_f(j,i) = plat_f(j,1) + sw*plat_f_x(j,i)   !HH: I think it should be this
                  plat_f(j,i) = plat_f(j,i) + sw*plat_f_x(j,i)
                  plat_der_f(j,i,l) = sw_der*plat_f_x(j,i)
               enddo
            enddo 

         enddo ! l = 1, dim_x 
      endif ! plates
   
      csca = shapefrac*plat_csca + (1.0-shapefrac)*colm_csca
      cext = shapefrac*plat_cext + (1.0-shapefrac)*colm_cext
      cabs = cext - csca

      der_csca = shapefrac*plat_der_csca - (1.0-shapefrac)*colm_der_csca
      der_cext = shapefrac*plat_der_cext - (1.0-shapefrac)*colm_der_cext
      der_cabs = der_cext - der_csca

      f_ss_full(:,:) = (plat_csca*shapefrac*plat_f(:,:) + colm_csca*(1.0-shapefrac)*colm_f(:,:))/csca
      do l = 1, dim_x
         der_fss_full(:,:,l) = f_ss_full(:,:)/csca*der_csca(l)-&
                              (shapefrac*plat_der_csca(l)*plat_f(:,:) + &
                              shapefrac*plat_csca*plat_der_f(:,:,l)+&
                              (1.0-shapefrac)*colm_der_csca(l)*colm_f(:,:) + &
                              (1.0-shapefrac)*colm_csca*colm_der_f(:,:,l))/csca
      enddo
      der_fss_full(:,:,dim_x+1) = 0.d0
      der_fss_full(:,:,dim_x+2) = 0.d0
   
!*** now interpolate to input scattering angle   
      imin = minval(minloc(dabs(scat_angle-scat_lut),scat_angle.ge.scat_lut))
      f_ss =  f_ss_full(imin,:) + &
             ((f_ss_full(imin+1,:) - f_ss_full(imin,:)) / (scat_lut(imin+1)-scat_lut(imin))) * &
             (scat_angle-scat_lut(imin))


      der_f_ss = der_fss_full(imin,:,:) + &
                 ((der_fss_full(imin+1,:,:) - der_fss_full(imin,:,:))/(scat_lut(imin+1)-scat_lut(imin))) * &
                 (scat_angle-scat_lut(imin))

      do i = 1, nper_in
         if (index_ist(i) .le. nstokes .and. index_jst(i) .le. nstokes) then
            z_ss(index_ist(i), index_jst(i)) = f_ss(i)
            der_z_ss(index_ist(i), index_jst(i), :) = der_f_ss(i,:)
         endif
      enddo
!      if (nstokes > 1) z_ss(2,1) = z_ss(1,2)  

   end subroutine get_cirrus_ssc
!------------------------------------------------------------------------------
!>  Compute derivatives of absorption and scattering cross sections 
!! to cirrus parameters for mutiple scattering
!------------------------------------------------------------------------------
   subroutine der_par_mode(&
              dcoefs_full, &
              dcsca_full, &
              dcabs_full,&
              nangle, &
              der_rnmb, &
              aer_col, &
              csca, &
              cabs,&
              dcoefs_par, &
              dcsca_par, &
              dcabs_par)
!*** input
   integer, intent(in) :: nangle
   real(double), intent(in) :: aer_col, csca, cabs       
   real(double), dimension(nfull), intent(in) :: dcsca_full, dcabs_full 
   real(double), dimension(nstokes, nstokes, 0:maxleg, nfull), intent(in) ::  dcoefs_full
   real(double), dimension(dim_x, 6), intent(in) :: der_rnmb
!*** output
   real(double), dimension(nstokes, nstokes, 0:maxleg, npar_mie), intent(out) :: dcoefs_par 
   real(double), dimension(npar_mie), intent(out) :: dcsca_par,  dcabs_par
!*** local variables
   integer :: j, k, l, i_st, j_st
!----------------------------------------------------------------------------  
      do i_st = 1, nstokes
         do j_st = 1, nstokes 

         do j = 0, min(maxleg,nangle)
            do k = 1, 2
               dcoefs_par(i_st, j_st, j, k) = 0.d0
               do l = 1, dim_x
                  dcoefs_par(i_st, j_st, j, k) = dcoefs_par(i_st, j_st, j, k) + &
                                      dcoefs_full(i_st, j_st, j, l)*der_rnmb(l, k)
               enddo
            enddo
            dcoefs_par(i_st, j_st, j, 3) = dcoefs_full(i_st, j_st, j, dim_x+1)
            dcoefs_par(i_st, j_st, j, 4) = dcoefs_full(i_st, j_st, j, dim_x+2)
            dcoefs_par(i_st, j_st, j, 5) = 0.d0
            do l = 1, dim_x
               dcoefs_par(i_st, j_st, j, 5) = dcoefs_par(i_st, j_st, j, 5) + &
                                   dcoefs_full(i_st, j_st, j, l)*der_rnmb(l, 5)
            enddo
            dcoefs_par(i_st, j_st, j, 6) = 0.d0
         enddo

         enddo

      enddo
  
      do k = 1, 2
         dcsca_par(k) = 0.d0
         dcabs_par(k) = 0.d0
         do l = 1, dim_x
            dcsca_par(k) = dcsca_par(k) + dcsca_full(l)*der_rnmb(l,k)        
            dcabs_par(k) = dcabs_par(k) + dcabs_full(l)*der_rnmb(l,k)
         enddo 
      enddo
      dcsca_par(3) = dcsca_full(dim_x+1)
      dcabs_par(3) = dcabs_full(dim_x+1)
      dcsca_par(4) = dcsca_full(dim_x+2)
      dcabs_par(4) = dcabs_full(dim_x+2)

      dcsca_par(5) = csca/(aer_col)
      dcabs_par(5) = cabs/(aer_col)  
      dcsca_par(6) = 0.d0
      dcabs_par(6) = 0.d0
  
      return

   end subroutine der_par_mode

!------------------------------------------------------------------------------
!>  Compute derivatives of absorption and scattering cross sections 
!! to cirrus parameters for single scattering
!------------------------------------------------------------------------------
   subroutine der_par_mode_ssc(&
              dz_ss_full, &
              dcsca_full, &
              dcabs_full,&
              der_rnmb, &
              aer_col, &
              csca, &
              cabs,&
              dz_ss_par, &
              dcsca_par, &
              dcabs_par)
   
!*** input
   real(double), intent(in) :: aer_col, csca, cabs       
   real(double), dimension(nfull), intent(in) :: dcsca_full, dcabs_full
   real(double), dimension(nstokes, nstokes, nfull), intent(in) :: dz_ss_full  
   real(double), dimension(dim_x, 6), intent(in) :: der_rnmb
!*** output
   real(double), dimension(npar_mie) , intent(out):: dcsca_par, dcabs_par 
   real(double), dimension(nstokes, nstokes, npar_mie), intent(out) :: dz_ss_par  
!*** local variables  
   integer :: i, j, k, l
!----------------------------------------------------------------------------
      do i = 1, nstokes
         do j = 1, nstokes
         do k = 1, 2
            dz_ss_par(i, j, k) = 0.d0
            do l = 1, dim_x
               dz_ss_par(i, j, k) = dz_ss_par(i, j, k) + &
                  dz_ss_full(i, j, l)*der_rnmb(l,k)
            enddo 
         enddo
         dz_ss_par(i, j, 3) = dz_ss_full(i, j, dim_x+1)
         dz_ss_par(i, j, 4) = dz_ss_full(i, j, dim_x+2)
         dz_ss_par(i, j, 5) = 0.d0
         do l = 1, dim_x
            dz_ss_par(i, j, 5) = dz_ss_par(i, j, 5) + &
                             dz_ss_full(i, j, l)*der_rnmb(l,5)
         enddo
         dz_ss_par(i, j, 6) = 0.d0
         enddo
      enddo
  
      do k = 1, 2
         dcsca_par(k) = 0.d0
         dcabs_par(k) = 0.d0
         do l = 1, dim_x
            dcsca_par(k) = dcsca_par(k) + dcsca_full(l)*der_rnmb(l,k)       
            dcabs_par(k) = dcabs_par(k) + dcabs_full(l)*der_rnmb(l,k)
         enddo
      enddo
      dcsca_par(3) = dcsca_full(dim_x+1)
      dcabs_par(3) = dcabs_full(dim_x+1)
      dcsca_par(4) = dcsca_full(dim_x+2)
      dcabs_par(4) = dcabs_full(dim_x+2)
 
      dcsca_par(5) = csca/(aer_col)
      dcabs_par(5) = cabs/(aer_col)
      dcsca_par(6) = 0.d0
      dcabs_par(6) = 0.d0
  
      return
   
   end subroutine der_par_mode_ssc
!------------------------------------------------------------------------------   

end module optic_cirrus_module
