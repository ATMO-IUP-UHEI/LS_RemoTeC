! modified by Andre Galli, February 25, 2011
!        CHANGES: Switched to spectral_response_lbl as default, since it is,
!        using the stored instrument response function, 0.1 sec faster for one typical retrieval
!        than the FFT-based subroutine spectral_response_custom
!------------------------------------------------------------------------------
!> Set up the arrays that store the ILS for convolution
!> @todo In case of fitting a spectral shift, recompute the arrays that store
!! the spectral response, since the wavelength grid can shift
!! this is currently done with spectral_response_lbl/custom, but these are
!! not up-to-date
!------------------------------------------------------------------------------
module spectral_response_module
  use header_module
  use auxiliary_routines_module
  implicit none
  private

  !*** types
  public :: instrument_response

  !*** procedures
  public :: spectral_response_stored, response_internal, spectral_response_lbl,&
            spectral_response_custom, spectral_response_create_gauss

  !------------------------------------------------------------------------------
  !> Instrument response as an array for each wavelength pixel
  type :: instrument_response
     !> Array for instrument response (dim: nwave_lo, nils)
     real(double),dimension(:,:), allocatable :: resp_store
     !> ILS wavelength grid (dim: nwave_lo, nils)
     real(double),dimension(:,:),allocatable :: ils_dwave
     !> Measured wavelength grid at which the ISRF is defined
     real, dimension(:), allocatable :: wavelength
     !> Number of ILS points
     integer :: nils
     !> Number of measured wavelengths for which an ISRF is defined
     integer :: nwave
  end type instrument_response
  !------------------------------------------------------------------------------

  interface
     subroutine spline(x, y, n, yp1, ypn, y2, ierr)
       integer, intent(in) :: n
       double precision, intent(in) :: yp1, ypn, x(n), y(n)
       double precision, intent(out) :: y2(n)
       integer, intent(out) :: ierr
     end subroutine spline
     subroutine splint(xa, ya, y2a, n, x, y, ierr)
       integer, intent(in) :: n
       double precision, intent(in) :: x, xa(n), y2a(n), ya(n)
       double precision, intent(out) :: y
       integer, intent(out) :: ierr
     end subroutine splint
     subroutine spline_interpol(x, y, n, xnew, ynew, nnew, ierr)
       integer, intent(in) :: n, nnew
       double precision, dimension(n), intent(in) :: x, y, xnew
       double precision, dimension(nnew), intent(out) :: ynew
       integer, intent(out) :: ierr
     end subroutine spline_interpol
  end interface
  !------------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------------
  !> Prepare the arrays describing the instrument response function
  !> @param[in] wfine       high-resolution wavelength grid of model
  !> @param[in] wpix        low-resolution wavelength grid of measurement
  !> @param[in] ilswave     delta(wavelength) grid of instrument line shape per wavelength pixel
  !> @param[in] ilsfunction instrument line function per wavelength pixel
  !> @param[out] resp_store array for storing instrument response
  !> @param[out] ie_store   array for storing instrument response, gives upper limit index on hi-res wavelength grid
  !> @param[out] is_store   array for storing instrument response, gives lower limit index on hi-res wavelength grid
  !------------------------------------------------------------------------------
  subroutine response_internal(&
       wfine, &
       wpix, &
       ilswave, &
       ilsfunction, &
       resp_store, &
       ie_store, &
       is_store, &
       ierr)
    !*** Input
    real(double), dimension(:), intent(in) :: wfine, wpix
    real(double), dimension(:,:), intent(in) :: ilsfunction, ilswave
    !*** Output
    real(double),dimension(:,:),allocatable, intent(out) :: resp_store
    integer,dimension(:), allocatable, intent(out) :: ie_store
    integer,dimension(:), allocatable, intent(out) :: is_store
    integer, intent(out) :: ierr
    !*** local variables
    integer ::   i, l, discrep, nils, nfine, npix
    integer ,dimension(1) :: is, ie
    real(double) :: ws, we, rnorm, dw, resp_gau, dummy
    real(double), dimension(size(wfine)) :: resp_linterp
    real(double), dimension(size(wfine)) :: dwave
    real(double), dimension(size(ilsfunction(1,:))) :: deriv2
    character(stringlen) :: message

    !*** Initialize
    ierr = 0
    nfine = size(wfine)
    npix = size(wpix)
    nils = size(ilsfunction(1,:))
    if(mod(nils,2) == 0) then
       write(message, *) 'RESPONSE_INTERNAL: gridsize of ISRF is even:',nils
       ierr = ierr_isrf
       goto 999
    endif
    if(allocated(resp_store)) deallocate(resp_store)
    if(allocated(ie_store)) deallocate(ie_store)
    if(allocated(is_store)) deallocate(is_store)

    allocate(resp_store(npix, 9999), &
         ie_store(npix), &
         is_store(npix),&
         stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'RESPONSE_INTERNAL: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    do l = 2, nfine-1
       dwave(l) = 0.5d0*(wfine(l)-wfine(l-1)) + 0.5d0*(wfine(l+1)-wfine(l))
    enddo
    dwave(1) = dwave(2)
    dwave(nfine) = dwave(nfine-1)
    do l = 1, npix

       ws = wpix(l) - maxval(ilswave(l,:))
       we = wpix(l) + maxval(ilswave(l,:))
       is = minloc(dabs(wfine-ws))
       is = max(is,1)
       ie = minloc(dabs(wfine-we))
       discrep = ie(1) - is(1) + 1 - nils

       !        if (discrep .gt. 0) then
       !           ie(1) = ie(1) - discrep
       !           if (abs(discrep) .gt. 1) then
       !              write(message, *) 'RESPONSE_INTERNAL: |ILS_dimensions-nils| > 1'
       !              ierr = ierr_isrf
       !              goto 999
       !           endif
       !        endif
       is_store(l) = is(1)
       ie_store(l) = ie(1)
!!$       call spline(ilswave(l,:), ilsfunction(l,:), nils,&
!!$             1.1D30, 1.1D30, deriv2, ierr)
!!$       if (ierr .ne. 0) then
!!$          write(message,*) 'SPECTRAL_RESPONSE_CUSTOM.SPLINE: memory allocation error'
!!$          ierr = ierr_all
!!$          goto 999
!!$       endif
       rnorm = 0.d0
       call linear_interpol(nils,ie(1)-is(1)+1,ilswave(l,:), wpix(l)-wfine(is(1):ie(1)),ilsfunction(l,:),resp_linterp(is(1):ie(1)))
       do i = is(1), ie(1)
!!$           dw = wpix(l) - wfine(i)
!!$           call splint(ilswave(l,:), ilsfunction(l,:), deriv2, nils, dw, resp_gau, ierr)
!!$           if (ierr .ne. 0) then
!!$              write(message,*) 'SPECTRAL_RESPONSE_CUSTOM.SPLINT: bad input'
!!$              ierr = ierr_isrf
!!$              goto 999
!!$           endif
!!$           resp_store(l,i-is(1)+1) = resp_gau
!!$           rnorm = rnorm +resp_gau*dwave(i)
          resp_store(l,i-is(1)+1) = resp_linterp(i)
          rnorm = rnorm +resp_linterp(i)*dwave(i)
       enddo
       do i = is(1), ie(1)
          dummy = dwave(i)/rnorm
          resp_store(l,i-is(1)+1) = resp_store(l,i-is(1)+1)*dummy
       enddo

    enddo

    return

999 continue
    call stopretrieval(message)

  end subroutine response_internal

  !------------------------------------------------------------------------------
  !> Convolve with instrument response function using stored arrays
  !------------------------------------------------------------------------------
  subroutine spectral_response_stored(resp_store, ie_store, is_store, rint_fine, rint_pix)
    !*** Input
    real(double), dimension(:), intent(IN) :: rint_fine
    !*** Input
    real(double),dimension(:,:), intent(in) :: resp_store
    integer, dimension(:), intent(in) :: ie_store               ! Array for instrument response (Dim: nwave_lo)
    integer, dimension(:), intent(in) :: is_store
    !*** Output
    real(double), dimension(:), intent(out) :: rint_pix          ! Convolved
    !*** local variables
    integer :: i,  l, npix
    real(double) ::  scale
    real(double), dimension(size(rint_fine)) :: infine
    !------------------------------------------------------------------------------
    npix = size(rint_pix)
    scale = maxval(DABS(rint_fine))
    if (scale .ne. 0.d0 ) then
       infine = rint_fine/scale
       do l = 1, npix
          rint_pix(l) = 0.d0
          do i = is_store(l), ie_store(l)
             rint_pix(l) = rint_pix(l) + resp_store(l, i-is_store(l)+1)*infine(i)
          enddo
          !*** Set too small values to 1d-8 times the largest value of the output:
          if(DABS(rint_pix(l)) .lt. 1d-8) rint_pix(l) = 1d-8
       enddo
       rint_pix = rint_pix*scale
    else
       rint_pix = 0.d0
    endif

  end subroutine spectral_response_stored

  !------------------------------------------------------------------------------
  !*** Not currently used
  !*** Convolve with instrument response function
  !*** The ils should already been calculated (with spectral_response_create_gauss)
  !*** and is used here as input (custom = win(n)%ilsfunction)
  !------------------------------------------------------------------------------
  subroutine spectral_response_lbl(&
       wfine, &                   ! high-resolution wavelength grid
       rint_fine, &               ! high-resolution spectrum
       wpix, &                    ! low-resolution wavelength-grid
       rint_pix, &                ! Convolved spectrum on low-resolution grid
       wcustom, &
       custom, &
       resp_store, &
       ie_store, &
       is_store, &
       ierr)
    !*** Input
    real(double), dimension(:), intent(in) :: wfine, rint_fine
    real(double), dimension(:), intent(in) :: wpix
    real(double), dimension(:), intent(in) :: custom, wcustom
    !*** Output
    real(double), dimension(:), intent(out) :: rint_pix          ! Convolved
    real(double),dimension(:,:),allocatable :: resp_store      ! Array for instrument response (Dim: nwave_lo,nils)
    integer,dimension(:), allocatable, intent(out) :: ie_store               ! Array for instrument response (Dim: nwave_lo)
    integer,dimension(:), allocatable, intent(out) :: is_store
    integer, intent(out) :: ierr
    !*** local variables
    integer ::   i, l, discrep, ncustom, nfine, npix
    integer ,dimension(1) :: is, ie
    real(double) :: ws, we, rnorm, dw, resp_gau, dummy, scale
    real(double), dimension(size(wfine)) :: dwave, infine
    real(double), dimension(size(custom)) :: deriv2
    !------------------------------------------------------------------------------

    nfine = size(wfine)
    ncustom = size(custom)
    npix = size(rint_pix)
    scale = maxval(DABS(rint_fine))
    if (scale .ne. 0.d0) then
       infine = rint_fine/scale
       !*** Prepare the arrays describing the instrument response function
       if(allocated(resp_store)) deallocate(resp_store)
       if(allocated(ie_store)) deallocate(ie_store)
       if(allocated(is_store)) deallocate(is_store)
       allocate(resp_store(npix,ncustom))
       allocate(ie_store(npix))
       allocate(is_store(npix))

       do l = 2, nfine-1
          dwave(l) = 0.5d0*(wfine(l)-wfine(l-1)) + 0.5d0*(wfine(l+1)-wfine(l))
       enddo
       dwave(1) = dwave(2)
       dwave(nfine) = dwave(nfine-1)
       call spline(wcustom, custom, ncustom,&
            1.1D30, 1.1D30, deriv2, ierr)
       do l = 1, npix
          rint_pix(l) = 0.d0
          ws = wpix(l) - maxval(wcustom)
          we = wpix(l) + maxval(wcustom)
          is = minloc(dabs(wfine-ws))
          is = max(is,1)
          ie = minloc(dabs(wfine-we))
          discrep = ie(1) - is(1) + 1 - ncustom
          if (discrep .gt. 0) then
             ie(1) = ie(1) - discrep
             if (abs(discrep) .gt. 1) call stopretrieval('SPECTRAL_RESPONSE_LBL: |ILS_dimensions-ncustom| > 1')
          endif
          is_store(l) = is(1)
          ie_store(l) = ie(1)
          rnorm = 0.d0
          do i = is(1), ie(1)
             dw = wpix(l) -wfine(i)
             call splint(wcustom, custom, deriv2, ncustom, dw, resp_gau, ierr)
             resp_store(l,i-is(1)+1) = resp_gau
             rnorm = rnorm +resp_gau*dwave(i)
          enddo
          do i = is(1), ie(1)
             dummy = dwave(i)/rnorm
             resp_store(l,i-is(1)+1) = resp_store(l,i-is(1)+1)*dummy
             rint_pix(l) = rint_pix(l) + resp_store(l,i-is(1)+1)*infine(i)
          enddo
          !*** Set too small values to 1d-8 times the largest value of the output:
          if(DABS(rint_pix(l)) .lt. 1d-8) rint_pix(l) = 1d-8
       enddo
       rint_pix = rint_pix*scale
    else
       rint_pix = 0.d0
    endif

  end subroutine spectral_response_lbl

  !------------------------------------------------------------------------------
  !*** Fast Fourier Convolution
  !------------------------------------------------------------------------------
  subroutine spectral_response_custom(&
       wfine, rint_fine, nfine, wpix, rint_pix, npix,&
       wcustom, custom, ncustom, ierr)
    implicit none
    !*** input
    integer, intent(IN) :: nfine, npix, ncustom
    real(double), intent(IN), dimension(nfine) :: wfine,rint_fine
    real(double), intent(IN), dimension(npix) :: wpix
    real(double), intent(IN), dimension(ncustom) :: custom,wcustom
    !*** output
    real(double), intent(OUT), dimension(npix) :: rint_pix
    integer, intent(out) :: ierr
    !*** local variables
    real(double) :: dlambda,scale,norm
    real(double), dimension(nfine) :: convfine,infine
    character(stringlen) :: message
    !------------------------------------------------------------------------------

    !*** Data and response function must be on same grid spacing and evenly sampled.
    dlambda = DABS(wcustom(2)-wcustom(1))
    if(DABS(dlambda - DABS(wfine(2)-wfine(1))) > dlambda/1.D6) then
       write(message, *) 'SPECTRAL_RESPONSE_CUSTOM: custom ILS is not on the same grid as data array.'
       ierr = ierr_isrf
       goto 999
    endif
    if(DABS(DABS(wfine(3)-wfine(2)) - dlambda) > dlambda/1.D6) then
       write(message, *) 'SPECTRAL_RESPONSE_CUSTOM: data are not evenly sampled.'
       ierr = ierr_isrf
       goto 999
    endif
    if(mod(ncustom,2) == 0) then
       write(message, *) 'SPECTRAL_RESPONSE_CUSTOM: gridsize of ISRF is even.'
       ierr = ierr_isrf
       goto 999
    endif

    !*** Call convolution
    scale = maxval(DABS(rint_fine))
    if (scale .ne. 0d0) then
       infine = rint_fine/scale
       norm = sum(custom)
       call ftconvolve(infine, custom/norm, convfine, ierr)
       if (ierr .ne. 0) return

       !*** Interpolate convolved array to output wavelength grid
       call spline_interpol(wfine, convfine, nfine,&
            wpix, rint_pix, npix, ierr)
       if (ierr .ne. 0) then
          write(message,*) 'SPECTRAL_RESPONSE_CUSTOM.SPLINE_INTERPOL.SPLINT: bad input'
          ierr = ierr_intrpl
          goto 999
       endif
       rint_pix = rint_pix*scale
    else
       rint_pix = 0.d0
    endif

    return
999 continue
    call stopretrieval(message)

  end subroutine spectral_response_custom

  !------------------------------------------------------------------------------
  !*** Adopted from Numerical Recipes "convlv.f"
  !------------------------------------------------------------------------------
  subroutine ftconvolve(dataraw, response, dataconv, ierr)
    real(double), dimension(:) :: dataraw, response
    real(double), dimension(:),intent(OUT) :: dataconv
    integer, intent(out) :: ierr
    !*** local variables
    integer :: i, k, m, n, n_new, mhalf
    real(double), dimension(:), allocatable :: in, res, out
    complex(double), dimension(:), allocatable :: ft1, ft2
    character(stringlen) :: message

    n = size(dataraw)
    m = size(response)

    !*** Zero-padding data array to 2**n array size
    n_new = 2**(int(log(n+n/2.)/log(2.))+1)
    if (n_new < m)  then
       write(message, *) 'FTCONVOLVE: FFT convolution not possible (n_new<m)'
       ierr = ierr_isrf
       goto 999
    endif
    allocate(in(n_new),res(n_new),out(n_new),ft1(n_new),ft2(n_new), stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'FTCONVOLVE: memory allocation error'
       ierr = ierr_all
       goto 999
    endif
    in = 0.
    res = 0.
    out = 0.
    ft1 = 0.
    ft2 = 0.

    in(1:n) = dataraw(1:n)
    in(n+1:n_new) = 0.

    !*** Wrap around response function and zero-padding in the middle to data array size
    res=0.
    mhalf=INT(m/2)
    res(1:m-mhalf)=response(mhalf+1:m)
    res(n_new-mhalf+1:n_new)=response(1:mhalf)

    !*** Fourier transform data and response function (NR)
    call dtwofft(in, res, ft1, ft2, n_new)

    !*** Multiply fourier transformed data and response function
    k = n_new/2
    do i = 1, k+1
       ft2(i) = ft1(i)*ft2(i)/k
    enddo
    ft2(1) = cmplx(dble(ft2(1)),dble(ft2(k+1)))

    !*** Inverse Fourier transform of convolution product
    call drealft(ft2, n_new, -1)

    do i=1,k
       out(i*2-1) = dble(ft2(i))
       out(i*2) = dimag(ft2(i))
    enddo

    !*** Truncate Fourier transform output
    dataconv = out(1:n)

    deallocate(in, res, out, ft1, ft2, stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'FTCONVOLVE: memory deallocation error'
       ierr = ierr_deall
       goto 999
    endif

    return
999 continue
    call stopretrieval(message)

  end subroutine ftconvolve

  !------------------------------------------------------------------------------
  !> @details Compute Gaussian instrument line shape (gauss) as a function of
  !! wavelength differences (wgauss)
  !------------------------------------------------------------------------------
  subroutine spectral_response_create_gauss(&
       fwhm, dlambda, ngauss, wgauss, gauss)
    integer, intent(IN) :: ngauss
    real(double), intent(IN) :: fwhm, dlambda
    real(double), intent(OUT), dimension(ngauss) :: wgauss, gauss
    !*** local variables
    integer :: k
    real(double) :: gaussint, center, b
    real(double),parameter :: sqrtln2 = 0.8325546

    b = fwhm/sqrtln2/2.
    gaussint = 0.
    gauss = 0.
    center = (int(ngauss/2)+1)*dlambda
    do k = 1, ngauss
       wgauss(k) = k*dlambda - center
       gauss(k) = DEXP(-(wgauss(k)/b)*(wgauss(k)/b))
       !gaussint = gaussint + gauss(k)
       gaussint = gaussint + gauss(k)*dlambda  !HH: normalization
    enddo
    gauss = gauss/gaussint

  end subroutine spectral_response_create_gauss

  !------------------------------------------------------------------------------

end module spectral_response_module
