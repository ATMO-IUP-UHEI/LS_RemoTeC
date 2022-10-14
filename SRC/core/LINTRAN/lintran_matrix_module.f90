!> \file lintran_matrix_module.f90
!! Module file for the radiative transfer matrix.

!> Module for the radiative transfer matrix.
!!
!! This matrix is defined on the split vertical grid, but it contains a lot of zeros.
!! For each destination interface, only two source interfaces have a nonzero contribution.
!! Therefore, the source index is reduced to just those two interfaces. Where the nonzero
!! elements are, is kept as administration, so that matrix operations can be applied correctly.
!! The matrix couples to an interpolated field, but the source should be normal, using the
!! routines with `sp' instead of `int', in lintran_transfer_module.
module lintran_matrix_module

  implicit none

contains

  !> Calculates the matrix and sets the administration of where the nonzero elements are if that is not yet done before.
  !!
  !! The matrix is coupled to an interpolated field, so it uses interpolated L-parameters
  !! (see lintran_tlc_module). The result of such a multiplication is a normal source. The
  !! matrix itself depends on the Fourier number, but the location of the nonzero elements
  !! does not. So, for the first Fourier number, the administration of the nonzero elements
  !! is calculated, but that can be recycled for all the other Fourier numbers. The matrix
  !! is calculated in an eight-dimensional form, where source and destination are sampled as
  !! a normal field (4D), though the source has a reduced interface dimension, 0--1, for the
  !! upper and lower layers that interact with the corresponding destination interface. We will
  !! always refer to these eight dimensions, though in the code, the Stokes and stream dimensions are
  !! merged to prevent compiler issues with arrays with more than seven dimensions.
  subroutine matrix_int(nst,nstrhalf,nlay,nsplit_total,calculate_shifts,check_diagonals,exist_bdrf,tlc,mat,shifts,stat) ! {{{

    use lintran_constants_module, only: protect_division, idn, iup, warningflag_instablematrix
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    logical, intent(in) :: calculate_shifts !< Flag for calculating where the nonzero elements are. Set to false to recycle the previously calculated ones.
    logical, intent(in) :: check_diagonals !< Flag for asserting that the diagonal terms are nonzero. Some solvers need that.
    logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst*nstrhalf,0:1,0:1,nst*nstrhalf,0:1,0:nsplit_total), intent(out) :: mat !< The radiative transfer matrix in band structure, to be calculated.
    integer, dimension(nst*nstrhalf,0:1,0:nsplit_total), intent(inout) :: shifts !< For each destination index, the number of zeros left of the block of nonzero values.
    integer, intent(inout) :: stat !< Error code.

    ! The matrix is eight-dimensional.
    ! First dimension: Source Stokes parameter.
    ! Second dimension: Source stream.
    ! Thrid dimension: Source dn/up.
    ! Fourth dimension: Source level (domain reduced to 0:1).
    ! Fifth dimension: Destination Stokes parameter.
    ! Sixth dimension: Destination stream.
    ! Seventh dimension: Destination dn/up.
    ! Eighth dimension: Destination level.

    ! Matrix calculation routines take a two-dimensional band matrix, which is
    ! directly interpreted from the eight-dimensional matrix, with
    ! source on the first dimension and destination on the second
    ! dimension. The first dimension is in a reduced form, because it
    ! is a band matrix. The rest of the rows are zero. The array shifts
    ! will tell where the non-zero elements are.

    ! First dimension: Source index.
    ! Second dimension: Destination index.

    ! In contrast to the paper, we use a parallel concept for upward and
    ! downward, with a place for roof-top reflection. Although that roof
    ! reflection is zero, it makes the world a bit more symmetric.

    ! Iterators.
    integer :: ilay ! Over levels, from 0 to nlay.
    integer :: isplit ! Over sublayers in a layer.
    integer :: istr ! Over streams.
    integer :: istr_src ! Over source streams.
    integer :: ist ! Over Stokes parameters.
    integer :: icom ! Over combined Stokes-stream parameters.

    ! Combining Stokes and stream directions for array operations.
    integer :: icom_low ! Index in combined dimension for lower edge.
    integer :: icom_high ! Index in combined dimension for upper edge.
    ! Same for the source direction.
    integer :: icom_low_src ! Index in combined dimension for lower edge for source direction.
    integer :: icom_high_src ! Index in combined dimension for upper edge for source direction.

    ! Counter for total sublayers.
    integer :: isplit_total

    ! If T is the transmission operator and S is the scattering operator,
    ! the Green matrix is T + TST + TSTST + TSTSTST + ...
    ! This is a Dysan series, which means that 1/G = 1/T - S.
    ! And it is this 1/G, which we want to have. Therefore, we
    ! have minus signs in S and we invert matrix T.

    ! Step 1. Put negative scatter terms.

    ! Initialize counter.
    isplit_total = 1

    do ilay = 1,nlay ! Layer in which scatter event takes place.
      do istr = 1,nstrhalf ! Destination stream.
        do istr_src = 1,nstrhalf ! Source stream.

          ! TLC-parameters for scattering.

          ! T = 1
          ! L = L_int_p(d)
          ! C = Cf or Cb

          ! Mind the minus sign.

          ! Scalar      :
          ! Array in s  :
          ! Array in d  : L_int_p
          ! Matrix(s,d) : prefactor_scat, C

          ! Dimension 1: Source Stokes parameter - Array operation.
          ! Dimension 2: Source stream - istr_src.
          ! Dimension 3: Source dn/up - Both.
          ! Dimension 4: Source level - Both 0 and 1.
          ! Dimension 1: Destination Stokes parameter - Another array operation.
          ! Dimension 6: Destination stream - istr.
          ! Dimension 7: Destination dn/up - Both.
          ! Dimension 8: Destination level - isplit_total-1 for d2 = up, isplit_total for d2 = dn.

          ! For the combined dimensions: (:,istr) becomes nst*(istr-1)+1 up to and
          ! including nst*istr.
          icom_low = nst*(istr - 1) + 1
          icom_high = nst*istr
          icom_low_src = nst*(istr_src - 1) + 1
          icom_high_src = nst*istr_src

          ! Forward scattering dn-dn.
          mat(icom_low_src:icom_high_src,idn,0,icom_low:icom_high,idn,isplit_total) = -tlc%prefactor_scat(istr_src,istr) * tlc%l_int_ac_p(istr,ilay) * tlc%c_fwd(:,:,istr_src,istr,ilay)
          mat(icom_low_src:icom_high_src,idn,1,icom_low:icom_high,idn,isplit_total) = -tlc%prefactor_scat(istr_src,istr) * tlc%l_int_sm_p(istr,ilay) * tlc%c_fwd(:,:,istr_src,istr,ilay)
          ! Forward scattering up-up.
          mat(icom_low_src:icom_high_src,iup,0,icom_low:icom_high,iup,isplit_total-1) = -tlc%prefactor_scat(istr_src,istr) * tlc%l_int_sm_p(istr,ilay) * tlc%c_fwd(:,:,istr_src,istr,ilay)
          mat(icom_low_src:icom_high_src,iup,1,icom_low:icom_high,iup,isplit_total-1) = -tlc%prefactor_scat(istr_src,istr) * tlc%l_int_ac_p(istr,ilay) * tlc%c_fwd(:,:,istr_src,istr,ilay)
          ! Backscattering dn-up.
          mat(icom_low_src:icom_high_src,idn,0,icom_low:icom_high,iup,isplit_total-1) = -tlc%prefactor_scat(istr_src,istr) * tlc%l_int_sm_p(istr,ilay) * tlc%c_bck(:,:,istr_src,istr,ilay)
          mat(icom_low_src:icom_high_src,idn,1,icom_low:icom_high,iup,isplit_total-1) = -tlc%prefactor_scat(istr_src,istr) * tlc%l_int_ac_p(istr,ilay) * tlc%c_bck(:,:,istr_src,istr,ilay)
          ! Backscattering up-dn.
          mat(icom_low_src:icom_high_src,iup,0,icom_low:icom_high,idn,isplit_total) = -tlc%prefactor_scat(istr_src,istr) * tlc%l_int_ac_p(istr,ilay) * tlc%c_bck(:,:,istr_src,istr,ilay)
          mat(icom_low_src:icom_high_src,iup,1,icom_low:icom_high,idn,isplit_total) = -tlc%prefactor_scat(istr_src,istr) * tlc%l_int_sm_p(istr,ilay) * tlc%c_bck(:,:,istr_src,istr,ilay)

          ! We will wait with the other sublayers within this layer,
          ! because they will get the same results. We will copy entire
          ! matrix parts later on.

        enddo
      enddo

      isplit_total = isplit_total + tlc%nsplit(ilay)

    enddo

    ! Surface reflection from downward to upward scattering. No
    ! interpolation is used here.

    if (exist_bdrf) then
      do istr = 1,nstrhalf ! Destination stream.
        do istr_src = 1,nstrhalf ! Source stream.

          ! TLC-parameters for surface reflection.

          ! T = 1
          ! L = 1
          ! C = Cb_srf

          ! Mind the minus sign.

          ! Matrix(s,d) : prefactor_scat, Cb_srf

          ! Dimension 1: Source Stokes parameter - Array operation.
          ! Dimension 2: Source stream - istr_src.
          ! Dimension 3: Source dn/up - dn.
          ! Dimension 4: Source level - 0 (high interface).
          ! Dimension 5: Destination Stokes parameter - Another array operation.
          ! Dimension 6: Destination stream - istr.
          ! Dimension 7: Destination dn/up - up.
          ! Dimension 8: Destination level - nsplit_total.

          ! Same indices for combined dimensions (1+2 and 5+6).
          icom_low = nst*(istr - 1) + 1
          icom_high = nst*istr
          icom_low_src = nst*(istr_src - 1) + 1
          icom_high_src = nst*istr_src

          mat(icom_low_src:icom_high_src,idn,0,icom_low:icom_high,iup,nsplit_total) = -tlc%prefactor_scat(istr_src,istr)*tlc%c_bck_srf(:,:,istr_src,istr)

        enddo
      enddo
    else
      ! Combining dimensions is trivial here, becuase (:,:) just becomes (:).
      mat(:,idn,0,:,iup,nsplit_total) = 0.0
    endif

    ! Set parts that are not yet initialized to zero. Those are forward
    ! surface reflection and forward and backward reflection on the roof of
    ! the atmosphere. All other terms of the matrix will be filled. Note that
    ! elements that are not in the matrix will be handled at the end, where
    ! the shifts list is discussed.
    mat(:,iup,0,:,iup,nsplit_total) = 0.0
    mat(:,:,1,:,idn,0) = 0.0

    ! Step 2: Add transmission terms.

    ! THe inverse of the transmission matrix is a diagonal of ones
    ! and negative transmission terms on the correct (same) side, but only
    ! reaching one level. As example, we write the matrix T and her
    ! inverse for a four-level system with only one stream in one
    ! direction. Here T_xy is the transmission from level x to y.

    ! Matrix T:
    ! 1     0     0     0
    ! T_01  1     0     0
    ! T_02  T_12  1     0
    ! T_03  T_13  T_23  1

    ! Matrix 1/T:
    ! 1     0     0     0
    ! -T_01 1     0     0
    ! 0     -T_12 1     0
    ! 0     0     -T_23 1

    ! Add the negative terms first.

    ! Initialize counter for sublayers.
    isplit_total = 1

    do ilay = 1,nlay ! Layer through which light is transmitted.
      do istr = 1,nstrhalf ! Source and destination stream.
        do ist = 1,nst ! Source and destination Stokes parameter.

          ! TLC-parameters for single-layer transmission.

          ! T = T
          ! L = 1
          ! C = 1

          ! There is no prefactor, but there is a minus sign.
          ! See comments at the top of this step.

          ! Dimension 1: Source Stokes parameter - ist.
          ! Dimension 2: Source stream - istr.
          ! Dimension 3: Source dn/up - Both.
          ! Dimension 4: Source level - 0 for d2 = dn, 1 for d2 = up.
          ! Dimension 5: Destination Stokes parameter - ist.
          ! Dimension 6: Destination stream - istr.
          ! Dimension 7: Destination dn/up - d2.
          ! Dimension 8: Destination level - isplit_total for d2 = dn, isplit_total-1 for d2 = up.

          ! Combine dimensions 1 and 2 and at the same time combining dimensions 5 and 6.
          icom = nst*(istr - 1) + ist

          mat(icom,idn,0,icom,idn,isplit_total) = mat(icom,idn,0,icom,idn,isplit_total) - tlc%t_sp(istr,ilay)
          mat(icom,iup,1,icom,iup,isplit_total-1) = mat(icom,iup,1,icom,iup,isplit_total-1) - tlc%t_sp(istr,ilay)

          ! We will wait with the other sublayers within this layer,
          ! because they will get the same results. We will copy entire
          ! matrix parts later on.

        enddo
      enddo

      ! Progress counter.
      isplit_total = isplit_total + tlc%nsplit(ilay)

    enddo

    ! Now, it is the time to copy the matrix elements to different
    ! sublayers, because all elements before are sprecfic to a layer
    ! and not to an interface. The last two dimensions have always been
    ! either idn,isplit_total or iup,isplit_total-1.
    ! Note that we do not have to copy the surface reflection term,
    ! because those indices, iup,nsplit_total, are not in the scope
    ! of the copying operation.

    ! Initialize counter for sublayers.
    isplit_total = 1

    do ilay = 1,nlay

      ! Copy the part of the matrix that has been defined above. That
      ! is everything except the diagonal ones.
      do isplit = 2,tlc%nsplit(ilay)

        mat(:,:,:,:,idn,isplit_total+isplit-1) = mat(:,:,:,:,idn,isplit_total)
        mat(:,:,:,:,iup,isplit_total+isplit-2) = mat(:,:,:,:,iup,isplit_total-1)

      enddo

      ! Progress counter.
      isplit_total = isplit_total + tlc%nsplit(ilay)

    enddo

    ! Add diagonal ones. These are interface-specific terms and be protected
    ! from the split-layer copy action that has just been performed.

    do isplit_total = 0,nsplit_total
      do istr = 1,nstrhalf
        do ist = 1,nst

          ! Dimension 1: Source Stokes parameter - ist.
          ! Dimension 2: Source stream - istr.
          ! Dimension 3: Source dn/up - Both.
          ! Dimension 4: Source level - 0 for d2 = up, 1 for d2 = dn.
          ! Dimension 5: Destination Stokes parameter - ist.
          ! Dimension 6: Destination stream - istr.
          ! Dimension 7: Destination dn/up - d2.
          ! Dimension 8: Destination level - isplit_total.

          ! Same combination for Stokes and stream parameters, which are the same for the
          ! source and the destination.
          icom = nst*(istr - 1) + ist

          mat(icom,idn,1,icom,idn,isplit_total) = mat(icom,idn,1,icom,idn,isplit_total) + 1.0
          mat(icom,iup,0,icom,iup,isplit_total) = mat(icom,iup,0,icom,iup,isplit_total) + 1.0

          ! Check for stability. The diagonals should be nonzero for solvers that use it
          ! in a division. That is Gauss-Seidel. LU-decomposition can work with zero diagonals,
          ! but not with zeros in the diagonal of U, but that is checked while decomposing.
          if (check_diagonals) then
            if (abs(mat(icom,idn,1,icom,idn,isplit_total)) .lt. protect_division) then
              stat = or(stat,warningflag_instablematrix)
              mat(icom,idn,1,icom,idn,isplit_total) = protect_division
            endif
            if (abs(mat(icom,iup,0,icom,iup,isplit_total)) .lt. protect_division) then
              stat = or(stat,warningflag_instablematrix)
              mat(icom,iup,0,icom,iup,isplit_total) = protect_division
            endif
          endif

        enddo
      enddo
    enddo

    ! Translate the matrix elements to a band matrix. In such a case
    ! it is good to visualize the matrix as a two-dimensional matrix,
    ! where the dimensions are combined four by four.

    ! The first nst*nstrhalf bands starts 2*nst*nstrhalf left of the matrix. Note
    ! that all numbers outside the matrix are zero, even if the atmosphere
    ! had had a reflecting roof. From the last nst*nstrhalf bands, the last
    ! 2*nst*nstrhalf elements are outside the matrix and zero.

    ! We do not want that, because we will multiply these zeros with
    ! some non-existing intensities outside the domain. That may be
    ! zero, but in principle it is a segmentation fault, because we
    ! are accessing intensities that are not there. To prevent
    ! administration with if-clauses, we will shift the matrix elements
    ! such that the zeros are also inside the matrix. Note that the
    ! meaning of the shifted rows in the eight-dimensional form becomes
    ! wrong.

    ! The first nst*nstrhalf rows will be move to the right. Therefore, the
    ! non-zero elements should be moved to the left.
    ! Dimension 1: Source Stokes parameter - All.
    ! Dimension 2: Source stream - All.
    ! Dimension 3: Source dn/up - Both.
    ! Dimension 4: Source level (0-1) - From 1 to 0.
    ! Dimension 5: Destination Stokes parameter - All.
    ! Dimension 6: Destination stream - All.
    ! Dimension 7: Destination dn/up - dn.
    ! Dimension 8: Destination level - 0.

    ! Dimension combination is trivial, (:,:) becomes (:).
    mat(:,:,0,:,idn,0) = mat(:,:,1,:,idn,0)
    ! Set zeros on the right. They used to be on the left, but because
    ! we know they are zero, no temp variable is needed.
    mat(:,:,1,:,idn,0) = 0.0

    ! Repeat the process for the last nstrhalf lines, but now to the
    ! other side.
    ! Dimension 1: Source Stokes parameter - All.
    ! Dimension 2: Source stream - All.
    ! Dimension 3: Source dn/up - Both.
    ! Dimension 4: Source level (0-1) - From 0 to 1.
    ! Dimension 5: Destination Stokes parameter - All.
    ! Dimension 6: Destination stream - All.
    ! Dimension 7: Destination dn/up - up.
    ! Dimension 8: Destination level - nsplit_total.

    ! Same trivial dimension combination.
    mat(:,:,1,:,iup,nsplit_total) = mat(:,:,0,:,iup,nsplit_total)
    ! Set zeros at the left.
    mat(:,:,0,:,iup,nsplit_total) = 0.0

    ! Calculate the array that indicates where the nonzero values are. The shift array
    ! will contain the number of zeros left of the first nonzero value. The moving of elements
    ! so that no elements outside the matrix exist must be taken into account for these shift.

    ! For dn at level 0: Shift is -2*nst*nstrhalf, corrected to 0.
    ! For up at level 0: Shift is 0.
    ! For dn at level 1: Shift is 0.
    ! For up at level 1: Shift is 2*nst*nstrhalf.
    ! For dn at level 2: Shift is 2*nst*nstrhalf.
    ! ...
    ! For up at level nsplit_total-1: Shift is (nsplit_total-1)*2*nst*nstrhalf.
    ! For dn at level nsplit_total: Shift is (nsplit_total-1)*2*nst*nstrhalf.
    ! For up at level nsplit_total: Shift is nsplit_total*2*nst*nstrhalf, corrected to (nsplit_total-1)*2*nst*nstrhalf.

    ! Shifts array only depends on the number of sublayers, so need not be
    ! calculated for every Fourier number.

    ! For the shifts, the Stokes parameters and the streams are also combined. Even gfortran
    ! can work with four-dimensional arrays, but we just want to be consistent with the matrix.
    if (calculate_shifts) then

      shifts(:,idn,0) = 0
      do isplit_total = 1,nsplit_total
        shifts(:,iup,isplit_total-1) = 2*(isplit_total-1)*nst*nstrhalf
        shifts(:,idn,isplit_total) = 2*(isplit_total-1)*nst*nstrhalf
      enddo
      shifts(:,iup,nsplit_total) = 2*(nsplit_total-1)*nst*nstrhalf

    endif

  end subroutine matrix_int ! }}}

  !> Performs Gauss-Seidel algorithm to solve \f$ I \f$ in equaiton \f$ LI = Q \f$.
  !!
  !! With Gauss-Seidel, an initial guess is needed. Assuming that the intensity in all
  !! but one index are equal to the first guess, that single index can be calculated.
  !! By iterating a few times through all indices, the intensity field will converge
  !! to the solution. The number of required iterations is reduced by solving the
  !! intensities in the order of transmission, the strongest interaction. That means
  !! that downward radiation is calculated from the top to the bottom and then, upward
  !! radiation is solved from bottom to top. Convergence is determined by the relative
  !! change in downward radiation at the bottom and upward radiation at the top.
  subroutine gauss_seidel(nst,nstrhalf,nsplit_total,nmat,nband,maxiter,tol,mat,shifts,src,intens,stat) ! {{{

    use lintran_constants_module, only: gs_protect_zero, warningflag_gsnonconvergence

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    integer, intent(in) :: nmat !< Full dimension size of the matrix.
    integer, intent(in) :: nband !< Number of nonzero elements in each matrix row.
    integer, intent(in) :: maxiter !< Maximum number of iterations through all intensities for convergence.
    real, intent(in) :: tol !< Convergence criterion.
    real, dimension(nband,nmat), intent(in) :: mat !< The radiative transfer matrix in band structure, in two-dimensional (source-destination) form.
    integer, dimension(nmat), intent(in) :: shifts !< For each destination index, the number of zeros left of the block of nonzero values.
    real, dimension(nmat), intent(in) :: src !< Source, \f$ Q \f$ in matrix equation \f$ LI = Q \f$.
    real, dimension(nmat), intent(inout) :: intens !< Intensity to be calculated, \f$ I \f$ in matrix equation \f$ LI = Q \f$. Start with a first guess and will converge to the solution.
    integer, intent(inout) :: stat !< Error code.

    ! The matrix and the intensities are folded into two and one dimension.

    ! Intensities that are used to check for convergence.
    real, dimension(2*nst*nstrhalf) :: crit_old
    real, dimension(2*nst*nstrhalf) :: crit_new

    logical :: converged ! Convergence flag.

    ! Iterators.
    integer :: iter ! Gauss-seidel iterator.
    integer :: icrit ! Over convergence criteria.
    integer :: isplit ! Over split layers.
    integer :: istr ! Over streams.
    integer :: ist ! Over Stokes parameters.

    ! Matrix indices.
    integer :: imat ! Current matrix elemen.
    integer :: imat_before_istr ! Matrix element before adding the stream dimension.
    integer :: imat_before_ist ! Matrix element before adding the Stokes dimension.
    integer, dimension(2*nst*nstrhalf) :: crit ! Indices where convergence is checked.

    ! Hardcode critical indices to the upward radiation at the top of the atmosphere
    ! and the downward radiation at the bottom of the atmosphere. That is the second
    ! nstrhalf elements and the next-to-last nstrhalf elements.
    crit = (/(nst*nstrhalf+icrit,icrit=1,nst*nstrhalf),(nmat-2*nst*nstrhalf+icrit,icrit=1,nst*nstrhalf)/)

    ! Save first result for convergence criterion.
    crit_old = intens(crit)

    do iter = 1,maxiter

      ! The loop constructions with isplit and istr are there to adapt the
      ! order of the indices in which the equations are solved. Because
      ! transmission terms are probably relative high compared to scattering
      ! terms, this index order probably follows the strongest couplings and
      ! thereby convergence will be achieved in much fewer iterations.

      ! First loop: downward radiation, starting at the top (0).
      do isplit = 0,nsplit_total

        ! Get matrix index besides the streams.
        imat_before_istr = 2*nst*nstrhalf*isplit

        do istr = 1,nstrhalf

          ! Add stream index to matrix index.
          imat_before_ist = imat_before_istr + nst*(istr-1)

          do ist = 1,nst

            ! Add Stokes index to matrix element.
            imat = imat_before_ist + ist

            ! Solve matrix equation for index imat with the current values for the
            ! other indices (That is Gauss-Seidel).
            intens(imat) = intens(imat) + (src(imat) - sum(mat(:,imat)*intens(shifts(imat)+1:shifts(imat)+nband))) / mat(imat-shifts(imat),imat)

          enddo

        enddo

      enddo

      ! Second loop: upward radiation, starting at the bottom (nsplit_total).
      do isplit = nsplit_total,0,-1

        ! Get matrix index besides the streams.
        imat_before_istr = 2*nst*nstrhalf*isplit + nst*nstrhalf

        do istr = 1,nstrhalf

          ! Add stream index to matrix index.
          imat_before_ist = imat_before_istr + nst*(istr-1)

          do ist = 1,nst

            ! Add Stokes index to matrix index.
            imat = imat_before_ist + ist

            ! Solve matrix equation for index imat with the current values for the
            ! other indices (That is Gauss-Seidel).
            intens(imat) = intens(imat) + (src(imat) - sum(mat(:,imat)*intens(shifts(imat)+1:shifts(imat)+nband))) / mat(imat-shifts(imat),imat)

          enddo

        enddo

      enddo

      ! Check for convergence.
      converged = .true.

      ! Take the sample of intensities use to check convergence.
      crit_new = intens(crit)

      do icrit = 1,2*nst*nstrhalf

        if (abs(crit_old(icrit)) .gt. gs_protect_zero) then
          if (abs((crit_new(icrit)-crit_old(icrit)) / crit_old(icrit)) .gt. tol) then
            converged = .false.
            exit ! No need to check the rest of the criterion.
          endif
        endif

      enddo

      ! Quit loop when convegence is achieved.
      if (converged) exit

      ! Update intensities at convergence criterion indices.
      crit_old = crit_new

    enddo

    ! Convergence is achieved in iter iterations, except when no convergence is
    ! achieved at all. By default, we will only warn for non-convergence.
    if (.not. converged) then
      stat = or(stat,warningflag_gsnonconvergence)
    endif

  end subroutine gauss_seidel ! }}}

  !> Performs LU-decomposition on the matrix.
  !!
  !! The input matrix is in a band structure, because all numbers outside the bands are
  !! zero. In an LU-decomposition, these zeros are preserved, because they are all away
  !! from the diagonal. The diagonal must be included in the band. Otherwise, it does not
  !! work. In Lintran, the diagonal is always included in the matrix, so that is no problem.
  !! The matrix will be turned into her LU-decomposition. The lower (L) matrix has ones on
  !! her diagonal and nonzero numbers for source index lower than destination index. Matrix
  !! U has nonzero numbers for source index at least as high as the destination index. The
  !! diagonal values of U are not necessarily one. Matrices L and U are stored in the same
  !! data structure as the original matrix, with L-values in elements where the source index
  !! is lower than the destination index and U-values in the other elements. The diagonal
  !! ones in L are lost, but they are always one, so no information is lost. With the bands
  !! and the shifts, those indices become a bit more complicated.
  subroutine lu_decompose(nmat,nband,mat,shifts,stat) ! {{{

    use lintran_constants_module, only: protect_division, warningflag_instablematrix

    implicit none

    ! Input and output.
    integer, intent(in) :: nmat !< Full dimension size of the matrix.
    integer, intent(in) :: nband !< Number of nonzero elements in each matrix row.
    real, dimension(nband,nmat), intent(inout) :: mat !< The radiative transfer matrix in band structure, in two-dimensional (source-destination) form. Will be transformed into the LU-decomposition.
    integer, dimension(nmat), intent(in) :: shifts !< For each destination index, the number of zeros left of the block of nonzero values.
    integer, intent(inout) :: stat !< Error code.

    ! Iterators.
    integer :: isrc ! Over the source dimension.
    integer :: idest ! Over the destination dimension.

    ! The source dimension is the first dimension and will be visualized
    ! as the horizontal dimension. In this dimension, the matrix has a
    ! reduced size (nband), so the indices of row idest will be a range
    ! from shifts(idest)+1 to shifts(idest)+nband. They must include the
    ! diagonal terms.

    ! We iterate over elements in L (excluding the diagonal ones),
    ! Those are the elements below the diagonal that are inside the bands.
    ! The left border is determined by the shifts and the right border
    ! is determined by the diagonal, because the bands include the
    ! diagonal.

    ! Check for division by zero. The diagonals before idest must be nonzero.
    ! From idest onwards, the diagonals can still change, so the check has to
    ! be performed after the iteration of idest. Iteration idest=1 is an empty
    ! loop over source, but it will contain the check for the diagonal terms.
    ! The last diagonal term must also be checked, because the diagonal terms
    ! are also needed by the solver.
    do idest = 1,nmat
      do isrc = shifts(idest)+1,idest-1
        ! We will take a multiple of the row where the destination index
        ! is equal to our source index isrc.
        ! Acquire multiplication factor and write it in our matrix as
        ! element in L.
        ! We have to reduce isrc by shifts(idest) when looking up matrix
        ! elements, because of the banded structure.
        mat(isrc-shifts(idest),idest) = mat(isrc-shifts(idest),idest) / mat(isrc-shifts(isrc),isrc)
        ! Apply multiplication on row idest, and only on the matrix
        ! elements that still contain original matrix elements, not
        ! the L-terms. The multiplication factor is still on our current
        ! matrix position.
        ! The source indices on which the multiplication has to be
        ! performed is bounded on the left by the current position,
        ! (isrc+1), which is certainly right of the left band border of
        ! row isrc, because it is right of the diagonal. The right border
        ! is determined by the band of the high row (isrc), so the first
        ! index ranges from isrc+1 to shifts(isrc)+nband. Those are the
        ! indices on the original matrix. We have to reduce the indices
        ! with the shifts of the second index, which is idest on the left
        ! hand side and isrc on the right hand side, so shifts(isrc) may
        ! be canceled out.
        mat(isrc+1-shifts(idest):shifts(isrc)+nband-shifts(idest),idest) = mat(isrc+1-shifts(idest):shifts(isrc)+nband-shifts(idest),idest) - mat(isrc-shifts(idest),idest)*mat(isrc+1-shifts(isrc):nband,isrc)
      enddo

      ! The stability check.
      if (abs(mat(idest-shifts(idest),idest)) .lt. protect_division) then
        stat = or(stat,warningflag_instablematrix)
        mat(idest-shifts(idest),idest) = protect_division
      endif

    enddo

  end subroutine lu_decompose ! }}}

  !> Solves matrix equation \f$ LI = Q \f$ using the LU-decomposition of the matrix.
  !!
  !! Do not confuse matrix \f$ L \f$ from the matrix equation and the lower (L) matrix
  !! of the LU-decomposition. The equation is solved exactly by solving \f$ LUI = Q \f$
  !! in two steps. First, \f$ LX = Q \f$, where \f$ X \f$ substitutes \f$ UI \f$. The lower (L)
  !! matrix has no contributions for sources with higher index than the destination, so
  !! when starting with low destination indices, the equations can be solved, using the
  !! already-calculated X-values with lower indices. When \f$ X \f$ is solved, the intensity
  !! can be solved using equation \f$ UI = X \f$ and use the same trick exactly the other
  !! way around, starting with the high indices.
  subroutine lu_solve(nmat,nband,mat,shifts,src,intens) ! {{{

    implicit none

    ! Input and output.
    integer, intent(in) :: nmat !< Full dimension size of the matrix.
    integer, intent(in) :: nband !< Number of nonzero elements in each matrix row.
    real, dimension(nband,nmat), intent(in) :: mat !< The radiative transfer matrix in band structure, in two-dimensional (source-destination) form. It must be transformed into the LU-decomposition.
    integer, dimension(nmat), intent(in) :: shifts !< For each destination index, the number of zeros left of the block of nonzero values.
    real, dimension(nmat), intent(in) :: src !< Source, \f$ Q \f$ in matrix equation \f$ LI = Q \f$.
    real, dimension(nmat), intent(out) :: intens !< Intensity to be calculated, \f$ I \f$ in matrix equation \f$ LI = Q \f$.

    real, dimension(nmat) :: u_intens ! Matrix U times the intensity.

    ! Iterator.
    integer :: idest ! Over the vertical matrix dimension (destination).

    ! The matrix equation mat*intens = src will be solved by substituting
    ! mat by L*U and first solve (U*intens) from L*(U*intens) = src.
    ! And we use the variable u_intens for U*intens.

    ! We have to start with low indices to solve u_intens.
    do idest = 1,nmat
      ! The equation looks like:
      ! L(:,idest) * u_intens(:) = src(idest).
      ! L contains only nonzero elements on the left part of the band,
      ! the part up to the diagonal.
      ! L(shifts(idest)+1:idest,idest) u_intens(same_indices) = src(idest)
      ! Put known parameters to the right, by subtracting everything
      ! from the left hand side, except for index (idest).
      ! L(idest,idest)*u_intens(idest) = src(idest) - L(shifts(idest)+1:idest-1,idest) * u_intens(same_indices).
      ! And L(idest,idest) is one, so the equation is solved.
      ! We have to reduce isrc by shifts(idest) when looking up matrix
      ! elements, because of the banded structure.
      u_intens(idest) = src(idest) - sum(mat(1:idest-1-shifts(idest),idest) * u_intens(shifts(idest)+1:idest-1))
    enddo

    ! Now solve the equation U*intens = u_intens, with descending indices.
    do idest = nmat,1,-1
      ! The equation looks like:
      ! U(:,idest) * intens(:) = u_intens(idest).
      ! Take only the nonzero elements in U. Those start at the diagonal
      ! and end at the right of the band.
      ! U(idest:shifts(idest)+nband,idest) * intens(same_indices) = u_intens(idest)
      ! Put known variables to the right. That is everything except for
      ! the diagonal term in U.
      ! U(idest,idest)*intens(idest) = u_intens(idest) - U(idest+1:shifts(idest)+nband,idest) * intens(same_indices)
      ! U(idest,idest) is not one, so, the equation is solved by dividing
      ! by it.
      ! intens(idest) = (u_intens(idest) - U(idest+1:shifts(idest)+nband,idest) * intens(same_indices)) / U(idest,idest)
      ! We have to reduce isrc by shifts(idest) when looking up matrix
      ! elements, because of the banded structure.
      intens(idest) = (u_intens(idest) - sum(mat(idest+1-shifts(idest):nband,idest) * intens(idest+1:shifts(idest)+nband))) / mat(idest-shifts(idest),idest)
    enddo

  end subroutine lu_solve ! }}}

end module lintran_matrix_module
