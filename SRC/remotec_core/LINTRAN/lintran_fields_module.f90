!> \file lintran_fields_module.f90
!! Module file for intensity fields.

!> Module for intensity fields.
!!
!! This modules constructs space for intensity fields. There are analytic fields
!! and interpolated fields. The latter benefit from layer splitting and therefore
!! have a varying spatial dimension size. Those fields will only be allocated during
!! the calculation and not during the initialization. These are the only arrays that are
!! allocated and deallocated during the calculation.
module lintran_fields_module

  implicit none

contains

  !> Construct space needed for fixed-sized intensity fields.
  !!
  !! Because this is still the initialization phase, it is not yet known to what extent
  !! layers will be split, so only those without layer splitting are allocated. When Lintran is
  !! not initialized to calculate derivatives, fewer fields will be needed and therefore
  !! fewer fields will be allocated.
  subroutine fields_init(nst,nstrhalf,nlay,flag_derivatives,fld,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation
    use lintran_types_module, only: lintran_fields_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    logical, intent(in) :: flag_derivatives !< Flag for possibility of calculating derivatives.
    type(lintran_fields_class), intent(out) :: fld !< Internal fields structure.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr
    integer :: progression

    progression = 0

    ! Allocate analytic fields, 2 necessary without differentiating.
    allocate(fld%a1(nst,nstrhalf,0:1,0:nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,fld,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(fld%a2(nst,nstrhalf,0:1,0:nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,fld,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Repeat initialization for additional fields for derivatives.
    if (flag_derivatives) then
      ! 1 additional analytic field makes 3.
      allocate(fld%a3(nst,nstrhalf,0:1,0:nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,fld,stat)
        return
      endif
      progression = progression + 1 ! }}}
    else
      progression = progression + 1 ! Fake allocation.
    endif

  end subroutine fields_init ! }}}

  !> Construct space needed for dynamic-sized intensity fields.
  !!
  !! Because every calculation can have a different layer splitting, the number of layers and
  !! therefore, the field sizes can change. Therefore, these are allocated during the calculation
  !! and not during the initialization. When derivatives are not calculated, fewer fields will be
  !! needed and therefore fewer fields will be allocated. Fewer fields are also allocated when
  !! Lintran was initialized for derivatives, but a non-derivatives run is performed.
  subroutine fields_dynamic_init(nst,nstrhalf,nsplit_total,flag_derivatives,fld,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation
    use lintran_types_module, only: lintran_fields_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    logical, intent(in) :: flag_derivatives !< Flag for possibility of calculating derivatives.
    type(lintran_fields_class), intent(inout) :: fld !< Internal fields structure.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr
    integer :: progression

    progression = 0

    ! Allocate interpolated fields, 2 necessary without differentiating.
    allocate(fld%i1(nst,nstrhalf,0:1,0:nsplit_total),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_dynamic(progression,flag_derivatives,fld,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(fld%i2(nst,nstrhalf,0:1,0:nsplit_total),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_dynamic(progression,flag_derivatives,fld,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! The matrix and the shifts. Stokes and stream dimensions are combined. Otherwise,
    ! gfortran cannot handle the number of dimensions. The reason that dimensions are
    ! combined for the shifts as well is to make them as similar to the matrix as possible.
    ! It is the matrix that has the many dimensions.
    allocate(fld%mat(nst*nstrhalf,0:1,0:1,nst*nstrhalf,0:1,0:nsplit_total),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_dynamic(progression,flag_derivatives,fld,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(fld%shifts(nst*nstrhalf,0:1,0:nsplit_total),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_dynamic(progression,flag_derivatives,fld,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Repeat initialization for additional fields for derivatives.
    if (flag_derivatives) then
      ! 2 additional interpolated fields makes 4.
      allocate(fld%i3(nst,nstrhalf,0:1,0:nsplit_total),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_dynamic(progression,flag_derivatives,fld,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(fld%i4(nst,nstrhalf,0:1,0:nsplit_total),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_dynamic(progression,flag_derivatives,fld,stat)
        return
      endif
      progression = progression + 1 ! }}}
    else
      progression = progression + 2 ! Fake allocations.
    endif

  end subroutine fields_dynamic_init ! }}}

  !> Cleans up the dynamic-sized fields.
  !!
  !! All arrays that are allocated in fields_dynamic_init are cleaned up.
  subroutine fields_dynamic_close(flag_derivatives,fld,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation
    use lintran_types_module, only: lintran_fields_class

    implicit none

    ! Input and output.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was initialized for derivatives.
    type(lintran_fields_class), intent(inout) :: fld !< Internal pixel structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    if (flag_derivatives) then
      deallocate(fld%i4,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(fld%i3,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(fld%shifts,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(fld%mat,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(fld%i2,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(fld%i1,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fields_dynamic_close ! }}}

  !> Cleans up the fixed-sized fields.
  !!
  !! All arrays that are allocated in fields_init are cleaned up.
  subroutine fields_close(flag_derivatives,fld,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation
    use lintran_types_module, only: lintran_fields_class

    implicit none

    ! Input and output.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was initialized for derivatives.
    type(lintran_fields_class), intent(inout) :: fld !< Internal pixel structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    if (flag_derivatives) then
      deallocate(fld%a3,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(fld%a2,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(fld%a1,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fields_close ! }}}

  !> Cleans up garbage that has been created when an error occurs during a routine.
  !!
  !! This is meant to clean up the rubbish that is created before the error occured.
  !! The idea is that the allocation status is the same as before calling the entire
  !! routine where the error occurred.
  subroutine fail(progression,flag_derivatives,fld,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation
    use lintran_types_module, only: lintran_fields_class

    implicit none

    ! Input and output.
    integer, intent(in) :: progression !< Indication of how much garbage has been created.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was initialized for derivatives.
    type(lintran_fields_class), intent(inout) :: fld !< Internal pixel structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    if (flag_derivatives) then
      if (progression .ge. 3) deallocate(fld%a3,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 2) deallocate(fld%a2,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 1) deallocate(fld%a1,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fail ! }}}

  !> Cleans up garbage that has been created when an error occurs during a routine.
  !!
  !! This is meant to clean up the rubbish that is created before the error occured.
  !! The idea is that the allocation status is the same as before calling the entire
  !! routine where the error occurred.
  subroutine fail_dynamic(progression,flag_derivatives,fld,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation
    use lintran_types_module, only: lintran_fields_class

    implicit none

    ! Input and output.
    integer, intent(in) :: progression !< Indication of how much garbage has been created.
    logical, intent(in) :: flag_derivatives !< Flag tor whether or the calculation involved derivatives.
    type(lintran_fields_class), intent(inout) :: fld !< Internal pixel structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    if (flag_derivatives) then
      if (progression .ge. 6) deallocate(fld%i4,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 5) deallocate(fld%i3,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 4) deallocate(fld%shifts,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 3) deallocate(fld%mat,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 2) deallocate(fld%i2,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 1) deallocate(fld%i1,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fail_dynamic ! }}}

end module lintran_fields_module
