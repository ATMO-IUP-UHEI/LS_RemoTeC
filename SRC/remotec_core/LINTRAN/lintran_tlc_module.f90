!> \file lintran_tlc_module.f90
!! Module file for the radiative transfer building blocks.

!> Module for the radiative transfer building blocks.
!!
!! This module translates the physical properties of the atmosphere to quantities
!! that will applied to the radiation for transmission and scattering. The T stands
!! for transmission. It is the multiplication that is applied when a photon travels
!! through the atmosphere. The L is for the integrals over atmospheric layers, indicating
!! the chance of interaction of a photon with the scattering properties of that atmospheric
!! layer. The C is the scatter term that stands outside the integral over the optical depth.
!! For symmetry, a part of this tau-independent term is applied as prefactor in all events
!! of source, scattering and/or response, so that the remainder of C becomes more symmetric.
!! Furthermore, the TLC module organizes derivatives of those parameters. Many of those
!! derivatives are intermediate results for the actual TLC-parameters, so those are calculated
!! even when no derivatives are needed. Only a few organizational parameters are included.
!! The idea is that in the raditative transfer, matrix and perturbation modules, no grid,
!! geometry and pixel instances are required, only the TLC instance.
module lintran_tlc_module

  implicit none

contains

  !> Allocates space for the radiative transfer building blocks.
  !!
  !! Most TLC-parameters (those building blocks) depend on the optical properties of
  !! the atmosphere. A few may only depend on the internal grid or the solar and viewing
  !! geometry. Those will be pointers, because they are calculated only once per run, or
  !! once per geometry. Those pointers are set, though tlc_init is called before a geometry
  !! is provided. Then pointers to non-initialized variables are set, but their targets will be
  !! initialized before the pointers are actually used.
  subroutine tlc_init(nst,nstrhalf,nlay,ngeo,flag_derivatives,grd,geo,pix,tlc,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation, fullpol_nst
    use lintran_types_module, only: lintran_grid_class, lintran_geometry_class, lintran_pixel_class, lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    logical, intent(in) :: flag_derivatives !< Flag for possibility of calculating derivatives.
    type(lintran_grid_class), intent(in), target :: grd !< Internal grid structure.
    type(lintran_geometry_class), intent(in), target :: geo !< Internal geometry structure.
    type(lintran_pixel_class), intent(in), target :: pix !< Internal pixel structure.
    type(lintran_tlc_class), intent(out), target :: tlc !< Internal radiative-transfer building-block structure.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr
    integer :: progression

    progression = 0

    ! Some TLC-parameters are initialized in the grid or geometry module.
    ! These will be pointed to from the TLC instance. The gemoetry
    ! dependent parameters will get their correct values when the
    ! geometry is set. But because it is a pointer, we do not have to
    ! re-copy the value from the geometry instance after the geometry
    ! has been set.

    ! For settings that depend on the viewing geometry: Those pointers are
    ! only set when a viewing geometry is selected. In the trasfer modules
    ! there is no bother about multiple viewing geometries. That administration
    ! is done by the geometry, TLC and main Lintran modules.

    ! Set pointer to the layer-splitting administration.
    tlc%nsplit => pix%nsplit

    ! Set pointers to prefactors.
    tlc%prefactor_dir => grd%prefactor_dir
    tlc%prefactor_src => grd%prefactor_src
    tlc%prefactor_resp => grd%prefactor_resp
    tlc%prefactor_scat => grd%prefactor_scat
    tlc%prefactor_intermediate => grd%prefactor_intermediate

    ! Set pointers to effective direction cosines.
    tlc%effmu_mp => grd%effmu_mp
    tlc%effmu_pp => grd%effmu_pp

    ! Effective direction cosines that depend on the atmosphere, for possible
    ! pseudo-spherical issues.
    allocate(tlc%effmu_m0(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%effmu_p0(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%effmu_mv_allgeo(nstrhalf,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%effmu_pv_allgeo(nstrhalf,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Only used for thermal emission.
    tlc%effmu_p => grd%mu

    ! T-parameters.
    allocate(tlc%t(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%t_tot_0(0:nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%t_tot_v_allgeo(0:nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%t_tot_0v_allgeo(0:nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! T-parameters in split environment. The cumulative transmission term
    ! in the split environment is a product of the cumulative transmission
    ! term in the original grid multiplied by some power of the split
    ! single-layer transmission term (t_tot*t_sp^x). In such a way, there need
    ! not be an allocation on the vertical grid with layer-splitting.
    allocate(tlc%t_sp(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%t_sp_0(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%t_sp_v_allgeo(nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! L-parameters. The L-terms are mu*(1-exp(tau/mu)), a common
    ! shape of the integral.

    ! Normal L-terms.
    allocate(tlc%l_0v_allgeo(nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_p0(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_m0(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_pv_allgeo(nstrhalf,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_mv_allgeo(nstrhalf,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! Single-stream L-term.
    allocate(tlc%l_p(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_v_allgeo(nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! l_0 is not needed, because these L-terms are only needed for thermal emission.

    ! L-parameters in split environment.
    allocate(tlc%l_sp_p0(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_sp_m0(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_sp_pv_allgeo(nstrhalf,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_sp_mv_allgeo(nstrhalf,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_sp_pp(nstrhalf,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_sp_mp(nstrhalf,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! L-terms that occur after adding another factor tau, so a linear
    ! function times an exponential is integrated.
    allocate(tlc%l_1_sp_p0(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_1_sp_pv_allgeo(nstrhalf,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Single-stream L-terms.
    allocate(tlc%l_sp_p(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! Single-layer special L-term.
    allocate(tlc%l_1_sp_p(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! L-parameters for interpolations.
    allocate(tlc%l_int_sm_p(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_int_sm_m(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_int_sm_0(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_int_sm_v_allgeo(nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_int_ac_p(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_int_ac_m(nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_int_ac_0(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%l_int_ac_v_allgeo(nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! L-terms for only interpolation, sm and ac are the same by definition.
    allocate(tlc%l_int(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Single-scattering C-parameter and her derivative.
    ! The derivative is used to calculate the non-derivative,
    ! so it is allocated even if no derivatives are calculated.
    allocate(tlc%c_bck_0v_elem(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%diff_c_bck_0v_elem_ssa(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! C-parameters that rely on some derivative-pointers.
    ! tlc%c_bck_srf_0v_elem is scalar.
    ! tlc%c_emi_v_elem is also scalar.

    ! Phase functions in Fourier series.
    ! These will be calculated per Fourier component in a separate routine.

    ! There is no need to save all geometries, because one combination geometry-Fourier
    ! is used only unce, so can be thrown away if either the geometry or the Fourier
    ! number changes.

    ! First dimension: Source Stokes parameter.
    ! Second dimension: Destination Stokes parameter.
    ! Third dimension: Source stream.
    ! Fourth dimension: Destination stream.
    ! Fifth dimension: Atmospheric layer.

    ! From stream to stream.
    allocate(tlc%c_fwd(nst,nst,nstrhalf,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%c_bck(nst,nst,nstrhalf,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! For the adjoint, some asymmetry may be introduced by Stokes
    ! parameter V. If V is not involved, the adjoint will exactly use
    ! the same array, but if V is involved, another array will be
    ! allocated for the adjoint with the correct symmetry relationships.

    ! From sun to stream, first two dimensions are removed.
    allocate(tlc%c_fwd_0_nrm(nst,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%c_bck_0_nrm(nst,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! From stream to instrument, second dimension is removed.
    allocate(tlc%c_fwd_v_nrm(nst,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%c_bck_v_nrm(nst,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! If V is present, allocate additional arrays, otherwise, just
    ! set the pointers to the same arrays as the normal ones.
    if (fullpol_nst(nst)) then ! V involved.
      allocate(tlc%c_fwd_0_adj_target(nst,nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%c_bck_0_adj_target(nst,nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%c_fwd_v_adj_target(nst,nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%c_bck_v_adj_target(nst,nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      tlc%c_fwd_0_adj => tlc%c_fwd_0_adj_target
      tlc%c_bck_0_adj => tlc%c_bck_0_adj_target
      tlc%c_fwd_v_adj => tlc%c_fwd_v_adj_target
      tlc%c_bck_v_adj => tlc%c_bck_v_adj_target
    else
      progression = progression + 4 ! Fake allocations.
      ! No assymetry, just take over the normal ones.
      tlc%c_fwd_0_adj => tlc%c_fwd_0_nrm
      tlc%c_bck_0_adj => tlc%c_bck_0_nrm
      tlc%c_fwd_v_adj => tlc%c_fwd_v_nrm
      tlc%c_bck_v_adj => tlc%c_bck_v_nrm
    endif

    ! Derivatives of C with respect to the single-scattering albedo.
    allocate(tlc%diff_c_fwd_ssa(nst,nst,nstrhalf,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%diff_c_bck_ssa(nst,nst,nstrhalf,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! The same joke as for the non-derivatives is repeated with
    ! always the normal arrays and possibly the adjoint as own
    ! array or as pointers to the same arrays as the normal ones.
    allocate(tlc%diff_c_fwd_0_nrm_ssa(nst,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%diff_c_bck_0_nrm_ssa(nst,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%diff_c_fwd_v_nrm_ssa(nst,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%diff_c_bck_v_nrm_ssa(nst,nstrhalf,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! And there is the if-clause for V.
    if (fullpol_nst(nst)) then
      allocate(tlc%diff_c_fwd_0_adj_ssa_target(nst,nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_c_bck_0_adj_ssa_target(nst,nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_c_fwd_v_adj_ssa_target(nst,nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_c_bck_v_adj_ssa_target(nst,nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      tlc%diff_c_fwd_0_adj_ssa => tlc%diff_c_fwd_0_adj_ssa_target
      tlc%diff_c_bck_0_adj_ssa => tlc%diff_c_bck_0_adj_ssa_target
      tlc%diff_c_fwd_v_adj_ssa => tlc%diff_c_fwd_v_adj_ssa_target
      tlc%diff_c_bck_v_adj_ssa => tlc%diff_c_bck_v_adj_ssa_target
    else
      progression = progression + 4 ! Fake allocations.
      tlc%diff_c_fwd_0_adj_ssa => tlc%diff_c_fwd_0_nrm_ssa
      tlc%diff_c_bck_0_adj_ssa => tlc%diff_c_bck_0_nrm_ssa
      tlc%diff_c_fwd_v_adj_ssa => tlc%diff_c_fwd_v_nrm_ssa
      tlc%diff_c_bck_v_adj_ssa => tlc%diff_c_bck_v_nrm_ssa
    endif

    ! Surface reflection, again with normal and adjoint.
    ! The first array is not subject to normal-adjoint issues,
    ! because it is already a matrix.
    allocate(tlc%c_bck_srf(nst,nst,nstrhalf,nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! And here the normal-adjoint onslaught begins.
    allocate(tlc%c_bck_srf_0_nrm(nst,nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(tlc%c_bck_srf_v_nrm(nst,nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! Here is the if. The code is getting bored of this normal-adjoint
    ! onslaught, so for a better explanation, go up or to the grid
    ! module, where this feature was built first.
    if (fullpol_nst(nst)) then
      allocate(tlc%c_bck_srf_0_adj_target(nst,nstrhalf),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%c_bck_srf_v_adj_target(nst,nstrhalf),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      tlc%c_bck_srf_0_adj => tlc%c_bck_srf_0_adj_target
      tlc%c_bck_srf_v_adj => tlc%c_bck_srf_v_adj_target
    else
      progression = progression + 2 ! Fake allocations.
      tlc%c_bck_srf_0_adj => tlc%c_bck_srf_0_nrm
      tlc%c_bck_srf_v_adj => tlc%c_bck_srf_v_nrm
    endif

    ! Emissivity, again with normal-adjoint.
    allocate(tlc%c_emi_nrm(nst,nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}
    if (fullpol_nst(nst)) then
      allocate(tlc%c_emi_adj_target(nst,nstrhalf),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      tlc%c_emi_adj => tlc%c_emi_adj_target
    else
      progression = progression + 1 ! Fake allocation.
      tlc%c_emi_adj => tlc%c_emi_nrm
    endif

    ! Thermal emission.
    allocate(tlc%c_therm(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,tlc,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! From here on, only stuff for derivatives is initialized.
    if (flag_derivatives) then

      ! Limiting sensitivity to tau for transmission terms, for optimization. This
      ! is translated to the dimension of layers to which are differentiated.
      allocate(tlc%ideriv_sensitive_0(0:nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%ideriv_sensitive_v_allgeo(0:nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%ideriv_sensitive_0v_allgeo(0:nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}

      ! Selected plane_parallel atmosphere (for optimization).
      tlc%plane_parallel => geo%plane_parallel

      ! Prefactor only needed for derivatives.
      tlc%prefactor_zip => grd%prefactor_zip

      ! Due to the interpolation, there is coupling between intensities at
      ! the same side (sm) of a layer and cross-terms to the other (ac) side
      ! of a layer.
      allocate(tlc%l_zip_sm(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%l_zip_ac(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}

      ! Derivatives of effective direction cosines.
      allocate(tlc%diff_effmu_p0_tau(nstrhalf,nlay,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_effmu_m0_tau(nstrhalf,nlay,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_effmu_pv_tau_allgeo(nstrhalf,nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_effmu_mv_tau_allgeo(nstrhalf,nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}

      ! Relative derivatives of T-parameters with respect to tau
      ! The reldiff_t_tot_0_tau is needed without derivatives, but not
      ! past the TLC-module, so then, just the one from geo will be used.
      tlc%reldiff_t_tot_0_tau => geo%reldiff_t_tot_0_tau
      tlc%reldiff_t_tau => grd%reldiff_t_tau

      ! Relative derivatives of T-parameters for a sublayer (chain rule).
      allocate(tlc%reldiff_t_sp_0_tau(nlay,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%reldiff_t_sp_v_tau_allgeo(nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%reldiff_t_sp_tau(nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}

      ! Derivatives of L-parameters.
      allocate(tlc%diff_l_0v_tau_allgeo(nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_p0_tau(nstrhalf,nlay,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_pv_tau_allgeo(nstrhalf,nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_m0_tau(nstrhalf,nlay,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_mv_tau_allgeo(nstrhalf,nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      ! Single-stream L-parameters for thermal emission.
      tlc%diff_l_p_tau => tlc%t
      allocate(tlc%diff_l_v_tau_allgeo(nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}

      ! Derivatives of L in split environment.
      allocate(tlc%diff_l_sp_p0_tau(nstrhalf,nlay,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_sp_pv_tau_allgeo(nstrhalf,nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_sp_m0_tau(nstrhalf,nlay,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_sp_mv_tau_allgeo(nstrhalf,nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      ! Single-stream L-parameters for thermal emission.
      allocate(tlc%diff_l_sp_p_tau(nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}

      ! Derivatives of interpolated L-parameters.
      allocate(tlc%diff_l_int_sm_p_tau(nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_int_sm_m_tau(nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_int_ac_p_tau(nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_int_ac_m_tau(nstrhalf,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_int_sm_0_tau(nlay,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_int_sm_v_tau_allgeo(nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_int_ac_0_tau(nlay,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(tlc%diff_l_int_ac_v_tau_allgeo(nlay,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      ! Just interpolation, no stream direction, independent of side.
      allocate(tlc%diff_l_int_tau(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}

      ! Derivatives with respect to phase coefficients will be acquired from
      ! the geometry and grid modules. These are used for chain rules in
      ! lintran_module as well.
      tlc%diff_c_ph => pix%ssa

      ! Fixed derivatives for surface reflection.
      tlc%diff_c_bck_srf_bdrf => grd%diff_c_bck_srf_bdrf
      tlc%diff_c_bck_srf_0_bdrf => geo%diff_c_bck_srf_0_bdrf

      ! Pointer for emissivity, for normal and adjoint, which differ at most
      ! by a minus sign for V. That is already taken care off in the grid
      ! module.
      tlc%diff_c_emi_nrm_emi => grd%diff_c_emi_nrm_emi
      tlc%diff_c_emi_adj_emi => grd%diff_c_emi_adj_emi

      ! Derivatives of the thermal C-parameter with respect to the single-scattering albedo.
      allocate(tlc%diff_c_therm_ssa(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}
      ! Derivatives of the thermal C-parameter with respect to the Planck curve.
      allocate(tlc%diff_c_therm_planck(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,tlc,stat)
        return
      endif
      progression = progression + 1 ! }}}

    else

      progression = progression + 34 ! Fake allocations.

    endif

  end subroutine tlc_init ! }}}

  !> Numerically safe version of the exponent function.
  !!
  !! This routine is needed because exponents of high numbers can cause numeric instability.
  !! Usually, the exponents are negative, but it is possible that they become positive and
  !! even very high for more complicated parameters and possibly for pseudo-spherical geometry.
  real elemental function expProtect(x) ! {{{

    use lintran_constants_module, only: protect_exponent

    implicit none

    ! Input and output.
    real, intent(in) :: x !< The exponent.

    if (x .le. protect_exponent) then
      expProtect = exp(x)
    else
      expProtect = 0.0
    endif

  end function expProtect ! }}}

  !> Numerically safe version of the division function.
  !!
  !! This routine is needed because divisions by low numbers can cause numeric instability.
  !! Usually, the numbers are large enough, but it is possible that they become very small or
  !! even zero for more complicated parameters and possibly for pseudo-spherical geometry.
  real elemental function divProtect(x,y) ! {{{

    use lintran_constants_module, only: protect_division

    implicit none

    ! Input and output.
    real, intent(in) :: x !< The numerator.
    real, intent(in) :: y !< The denominator.

    if (abs(y) .lt. protect_division) then
      divProtect = 0.0
    else
      divProtect = x/y
    endif

  end function divProtect ! }}}

  !> Calculates all parameters that are the same for all Fourier numbers.
  !!
  !! Those are T and L parameters, and some effective direction cosines, which
  !! can acts as additions to L. All C-parameters depend on Fourier number, except for
  !! the ones for single-scattering, but they are handled independently, because though they
  !! do not depend on Fourier number, they do depend on viewing Stokes parameter and viewing
  !! geometry. Because the Fourier loop is outside the loop over viewing geometries, both the
  !! Fourier-independent parameters that depend on viewing geometry and those that do not
  !! depend on viewing geometry are calculated at once and pointers are set when a viewing
  !! geometry is selected.
  subroutine calculate_fourier_independent_tlc_parameters(nstrhalf,nlay,ngeo,flag_derivatives,execution,split_double,interpolation,grd,geo,pix,tlc,stat,ilay_deriv_base) ! {{{

    use lintran_constants_module, only: minimum_mudiff, warningflag_instableangle
    use lintran_types_module, only: lintran_grid_class, lintran_geometry_class, lintran_pixel_class, lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    logical, intent(in) :: flag_derivatives !< Flag for calculating derivatives.
    logical, dimension(3), intent(in) :: execution !< Flags for executing single, double and multi-scattering, in that order.
    logical, intent(in) :: split_double !< Switch for solving double-scattering analytically.
    integer, intent(in) :: interpolation !< Interpolation method for intensity inside layer. 0: Averaging. 1: Linear.
    type(lintran_grid_class), intent(in) :: grd !< Internal grid structure.
    type(lintran_geometry_class), intent(in) :: geo !< Internal geometry structure.
    type(lintran_pixel_class), intent(in) :: pix !< Internal pixel structure.
    type(lintran_tlc_class), intent(inout) :: tlc !< Internal radiative-transfer building-block structure.
    integer, intent(inout) :: stat !< Error code.
    ! Optional argument, obligatory for flag_derivatives set to true.
    integer, dimension(:), intent(in), optional :: ilay_deriv_base !< Layers for which derivatives with respect to optical depth and single-scattering albedo are calculated. Only needed when flag_derivatives is turned on. The dimension size: nlay_deriv_base is not given because that does not work for optional arguments.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.
    integer :: igeo ! Over viewing geometries.

    ! TLC-parameters for temporary use inside this routine.

    ! Optical depths of a split layer.
    real, dimension(nlay) :: tau_sp

    ! Effective direction cosines for path lengths for the sun and instrument. Those
    ! even depend on the optical depth.
    real, dimension(nlay) :: effmu_0
    real, dimension(nlay,ngeo) :: effmu_v_allgeo
    real, dimension(nlay,ngeo) :: effmu_0v_allgeo

    ! For stability check on effective angles.
    real :: effmu_x_original

    ! Single-layer plane-parallel transmission terms for sun and instrument
    ! These are used for derivatives of L. Note that single-layer stuff is
    ! always plane-parallel, because spherical geometry is only applied for the
    ! cumulative path from the sun or instrument to the first (or last) scattering
    ! layer.
    real, dimension(nlay) :: t_0
    real, dimension(nlay,ngeo) :: t_v_allgeo

    ! Derivatives of these transmission terms.
    real, dimension(nlay,nlay) :: reldiff_t_0_tau
    real, dimension(nlay,nlay,ngeo) :: reldiff_t_v_tau_allgeo

    ! Their logarithms.
    real, dimension(nlay) :: log_t_0
    real, dimension(nlay,ngeo) :: log_t_v_allgeo

    ! Temporarily defined aggregated transmission terms and her derivative.
    real, dimension(nlay) :: t_tmp
    real, dimension(nlay) :: log_t_tmp
    real, dimension(nlay,nlay) :: reldiff_t_tmp_tau

    ! Logarithms of cumulative transmission terms, equal to minus the optical depth
    ! with path length corrections.
    real, dimension(0:nlay) :: log_t_tot_0
    real, dimension(0:nlay,ngeo) :: log_t_tot_v_allgeo

    ! Single-stream L-parameters as intermediate result for interpolated L-parameters.
    real, dimension(nlay) :: l_sp_0
    real, dimension(nlay,ngeo) :: l_sp_v_allgeo
    ! l_sp_p is in the TLC-structure.
    real, dimension(nstrhalf,nlay) :: l_sp_m

    ! Derivatives of single-stream L-parameters as intermediate results for interpolated
    ! L-parameters.
    real, dimension(nlay,nlay) :: diff_l_sp_0_tau
    real, dimension(nlay,nlay,ngeo) :: diff_l_sp_v_tau_allgeo
    ! diff_l_sp_p_tau is in the TLC-structure.
    real, dimension(nstrhalf,nlay) :: diff_l_sp_m_tau

    ! Derivatives of single-stream effective direction cosines.
    real, dimension(nlay,nlay) :: diff_effmu_0_tau
    real, dimension(nlay,nlay,ngeo) :: diff_effmu_v_tau_allgeo

    ! Search index.
    integer, dimension(1) :: ind

    ! Save optical depths of single sublayers. They are also needed for double-scattering
    ! with derivatives. Not directly, but a temporary TLC-parameter is created to generate
    ! the derivatives of effective direction cosines, which are necessary for derivatives
    ! of double-scattering. The reason to choose a split TLC-parameter is that it can be
    ! re-used for interpolated multi-scattering TLC-parameters.
    if (execution(3)) tau_sp = pix%tau / pix%nsplit

    ! Transmission terms (T).

    ! Always necessary: t_tot_0 and t_tot_v, for single-scattering or source and response.

    ! Calculate transmission terms between interface ilay and the top
    ! of the atmosphere, and the single-layer exponent terms in layer ilay.

    ! The relative derivatives have been calculated in the geometry module and can
    ! have pseudo-spherical geometry. Therefore, it is no longer easy. It is even
    ! possible that some layers below the interface are crossed due to the spherical
    ! geometry.

    ! Aggregate the relative derivatives to get the logirithms of the transmission terms.
    if (geo%plane_parallel) then
      ! Only keep t_0 and t_v as intermediate results, because it does not delay
      ! the computation. Do not bother to calculate logatirhms, because they are
      ! not needed.
      if (geo%mu_0 .gt. 0.0) then
        tlc%t_tot_0(0) = pix%sun
      else
        tlc%t_tot_0(0) = 0.0 ! No sunlight at top of atmosphere.
      endif
      do ilay = 1,nlay
        t_0(ilay) = expProtect(-pix%tau(ilay)/geo%mu_0)
        tlc%t_tot_0(ilay) = tlc%t_tot_0(ilay-1) * t_0(ilay)
      enddo
      do igeo = 1,ngeo
        if (geo%mu_v_allgeo(igeo) .gt. 0.0) then
          tlc%t_tot_v_allgeo(0,igeo) = 1.0
        else
          tlc%t_tot_v_allgeo(0,igeo) = 0.0 ! Light at top of atmosphere cannot reach the instrument.
        endif
        do ilay = 1,nlay
          t_v_allgeo(ilay,igeo) = expProtect(-pix%tau(ilay)/geo%mu_v_allgeo(igeo))
          tlc%t_tot_v_allgeo(ilay,igeo) = tlc%t_tot_v_allgeo(ilay-1,igeo) * t_v_allgeo(ilay,igeo)
        enddo
      enddo

      if (execution(2) .or. execution(3)) then

        ! Write actual angles as effective angles, because it is plane-parallel
        ! With only single-scattering, we will work around. With double or multi-scattering,
        ! it is a lot of extra if-clauses, while it is less bad that an array of same
        ! values are written.
        effmu_0 = geo%mu_0 ! Array operation over layers in effmu_0.
        do ilay = 1,nlay
          tlc%effmu_p0(:,ilay) = geo%effmu_p0_ppg ! Array operation in streams.
          tlc%effmu_m0(:,ilay) = geo%effmu_m0_ppg ! Array operation in streams.
        enddo
        do igeo = 1,ngeo
          effmu_v_allgeo(:,igeo) = geo%mu_v_allgeo(igeo)
          do ilay = 1,nlay
            tlc%effmu_pv_allgeo(:,ilay,igeo) = geo%effmu_pv_ppg_allgeo(:,igeo) ! Array operation in streams.
            tlc%effmu_mv_allgeo(:,ilay,igeo) = geo%effmu_mv_ppg_allgeo(:,igeo) ! Array operation in streams.
          enddo
        enddo

        ! The effective direction cosine for single-scattering geometry (effmu_0v)
        ! is not copied over the layers, but is directly used from the geometry module.
        ! The reason is that this effective direction cosine need not be saved. It is only
        ! used in this routine.

      endif

    else

      ! For pseudo-spherical geometry, the transmission terms depend on optical depths
      ! of different layers. Also, the effective direction cosines with the sun or the
      ! instrument depend on optical depths.

      ! Calculate everything with intermediate results.
      log_t_tot_0(0:geo%visible_0) = matmul(pix%tau,geo%reldiff_t_tot_0_tau(:,0:geo%visible_0))
      do igeo = 1,ngeo
        log_t_tot_v_allgeo(0:geo%visible_v_allgeo(igeo),igeo) = matmul(pix%tau,geo%reldiff_t_tot_v_tau_allgeo(:,0:geo%visible_v_allgeo(igeo),igeo))
      enddo
      ! The rest is actually minus infinity.

      ! Create transmisison terms for visible layers.
      tlc%t_tot_0(0:geo%visible_0) = pix%sun * expProtect(log_t_tot_0(0:geo%visible_0))
      tlc%t_tot_0(geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
      do igeo = 1,ngeo
        tlc%t_tot_v_allgeo(0:geo%visible_v_allgeo(igeo),igeo) = expProtect(log_t_tot_v_allgeo(0:geo%visible_v_allgeo(igeo),igeo))
        tlc%t_tot_v_allgeo(geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
      enddo

      ! Acquire some terms that can be used inside this module
      ! Single-layer transmission terms.
      t_0(1:geo%visible_0) = divProtect(tlc%t_tot_0(1:geo%visible_0) , tlc%t_tot_0(0:geo%visible_0-1))
      log_t_0(1:geo%visible_0) = log_t_tot_0(1:geo%visible_0) - log_t_tot_0(0:geo%visible_0-1)
      do igeo = 1,ngeo
        t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) = divProtect(tlc%t_tot_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) , tlc%t_tot_v_allgeo(0:geo%visible_v_allgeo(igeo)-1,igeo))
        log_t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) = log_t_tot_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) - log_t_tot_v_allgeo(0:geo%visible_v_allgeo(igeo)-1,igeo)
      enddo

      ! Calculate effective direction cosines for path lengts in local layers.
      ! These effective angles depend on the optical depths of several layers.
      ! If the interface below a layer is not visible, nothing can be done with
      ! that layer.
      effmu_0(1:geo%visible_0) = -divProtect(pix%tau(1:geo%visible_0) , log_t_0(1:geo%visible_0))
      ! The rest of the effective cosines are irrelevant. Setting to zero must be okay.
      effmu_0(geo%visible_0+1:nlay) = 0.0
      do igeo = 1,ngeo
        effmu_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) = -divProtect(pix%tau(1:geo%visible_v_allgeo(igeo)) , log_t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo))
        effmu_v_allgeo(geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0
      enddo

      ! Tau-dependent stability check for pseudo-spherical geometry will concern
      ! t_x, log_t_x and effmu_x. For plane-parallel geometry, the stability check
      ! has already been performed in the geometry module, because then it is
      ! tau-independent.

      ! The internal grid should be such that only one stream can trigger instability
      ! and after correction, it is impossible that another stream also triggers instability.
      ! However, the solar direction can trigger instability with a viewing direction at
      ! any time. That is checked after each movement, but first, the viewing direction
      ! move away from the negative solar direction.

      ! Effective solar zenith angles.
      do ilay = 1,geo%visible_0
        effmu_x_original = effmu_0(ilay)
        ! Check for plus or minus internal streams.
        do istr = 1,nstrhalf
          if (abs(effmu_0(ilay) - grd%mu(istr)) .lt. minimum_mudiff) then
            if (effmu_0(ilay) .lt. grd%mu(istr)) then
              effmu_0(ilay) = grd%mu(istr) - minimum_mudiff
            else
              effmu_0(ilay) = grd%mu(istr) + minimum_mudiff
            endif
            exit
          endif
          if (abs(effmu_0(ilay) + grd%mu(istr)) .lt. minimum_mudiff) then
            if (effmu_0(ilay) .lt. -grd%mu(istr)) then
              effmu_0(ilay) = -grd%mu(istr) - minimum_mudiff
            else
              effmu_0(ilay) = -grd%mu(istr) + minimum_mudiff
            endif
            exit
          endif
        enddo
        ! Check for zero.
        if (abs(effmu_0(ilay)) .lt. minimum_mudiff) then
          if (effmu_0(ilay) .lt. 0.0) then
            effmu_0(ilay) = -minimum_mudiff
          else
            effmu_0(ilay) = minimum_mudiff
          endif
        endif
        ! Update dependent TLC-parameters on the change of the effective direction cosine.
        if (effmu_x_original .ne. effmu_0(ilay)) then
          stat = or(stat,warningflag_instableangle)
          ! The cumulative transmission terms are left as they are. The difference should
          ! not be big. And making the entire atmosphere consistent with this change may
          ! result in instabilities elsewhere.
          log_t_0(ilay) = -pix%tau(ilay) / effmu_0(ilay)
          t_0(ilay) = expProtect(log_t_0(ilay))
        endif
      enddo
      ! Effective viewing zenith angles.
      do igeo = 1,ngeo
        do ilay = 1,geo%visible_v_allgeo(igeo)
          effmu_x_original = effmu_v_allgeo(ilay,igeo)
          ! Initial check for the negative solar direction. Hereafter, this check
          ! is repeated each time the effective angle is moved. If this movement here
          ! enters the instability domain of another stream, it will be fixed there.
          if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) then
            if (effmu_v_allgeo(ilay,igeo) .lt. -effmu_0(ilay)) then
              effmu_v_allgeo(ilay,igeo) = -effmu_0(ilay) - minimum_mudiff
            else
              effmu_v_allgeo(ilay,igeo) = -effmu_0(ilay) + minimum_mudiff
            endif
          endif
          ! Check for plus or minus internal streams.
          do istr = 1,nstrhalf
            if (abs(effmu_v_allgeo(ilay,igeo) - grd%mu(istr)) .lt. minimum_mudiff) then
              if (effmu_v_allgeo(ilay,igeo) .lt. grd%mu(istr)) then
                effmu_v_allgeo(ilay,igeo) = grd%mu(istr) - minimum_mudiff
                ! Very unlikely: if effmu_v moves close to -effmu_0, there will be additional
                ! instability. While loop can be done at most twice. Therefore, two if-statements
                ! replace the while loop, because there is a razzia for while loops.
                if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) - minimum_mudiff
                if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) - minimum_mudiff
              else
                effmu_v_allgeo(ilay,igeo) = grd%mu(istr) + minimum_mudiff
                ! Very unlikely: if effmu_v moves close to -effmu_0, there will be additional
                ! instability. While loop can be done at most twice. Therefore, two if-statements
                ! replace the while loop, because there is a razzia for while loops.
                if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) + minimum_mudiff
                if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) + minimum_mudiff
              endif
              exit
            endif
            if (abs(effmu_v_allgeo(ilay,igeo) + grd%mu(istr)) .lt. minimum_mudiff) then
              if (effmu_v_allgeo(ilay,igeo) .lt. -grd%mu(istr)) then
                effmu_v_allgeo(ilay,igeo) = -grd%mu(istr) - minimum_mudiff
                ! Very unlikely: if effmu_v moves close to -effmu_0, there will be additional
                ! instability. While loop can be done at most twice. Therefore, two if-statements
                ! replace the while loop, because there is a razzia for while loops.
                if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) - minimum_mudiff
                if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) - minimum_mudiff
              else
                effmu_v_allgeo(ilay,igeo) = -grd%mu(istr) + minimum_mudiff
                ! Very unlikely: if effmu_v moves close to -effmu_0, there will be additional
                ! instability. While loop can be done at most twice. Therefore, two if-statements
                ! replace the while loop, because there is a razzia for while loops.
                if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) + minimum_mudiff
                if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) + minimum_mudiff
              endif
              exit
            endif
          enddo
          ! Check for zero.
          if (abs(effmu_v_allgeo(ilay,igeo)) .lt. minimum_mudiff) then
            if (effmu_v_allgeo(ilay,igeo) .lt. 0.0) then
              effmu_v_allgeo(ilay,igeo) = -minimum_mudiff
              ! Very unlikely: if effmu_v moves close to -effmu_0, there will be additional
              ! instability. While loop can be done at most twice. Therefore, two if-statements
              ! replace the while loop, because there is a razzia for while loops.
              if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) - minimum_mudiff
              if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) - minimum_mudiff
            else
              effmu_v_allgeo(ilay,igeo) = minimum_mudiff
              ! Very unlikely: if effmu_v moves close to -effmu_0, there will be additional
              ! instability. While loop can be done at most twice. Therefore, two if-statements
              ! replace the while loop, because there is a razzia for while loops.
              if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) + minimum_mudiff
              if (abs(effmu_v_allgeo(ilay,igeo)+effmu_0(ilay)) .lt. minimum_mudiff) effmu_v_allgeo(ilay,igeo) = effmu_v_allgeo(ilay,igeo) + minimum_mudiff
            endif
          endif
          ! Update dependent TLC-parameters on the change of the effective direction cosine.
          if (effmu_x_original .ne. effmu_v_allgeo(ilay,igeo)) then
            stat = or(stat,warningflag_instableangle)
            ! The cumulative transmission terms are left as they are. The difference should
            ! not be big. And making the entire atmosphere consistent with this change may
            ! result in instabilities elsewhere.
            log_t_v_allgeo(ilay,igeo) = -pix%tau(ilay) / effmu_v_allgeo(ilay,igeo)
            t_v_allgeo(ilay,igeo) = expProtect(log_t_v_allgeo(ilay,igeo))
          endif
        enddo
      enddo

      ! Calculate effective angle cosines.
      ! Those are acquired like summing up parallel resistors in an electric
      ! circuit. The p stands for a positive stream angle, the m stands
      ! for a negative stream angle.
      if (execution(1) .or. execution(2)) then
        do igeo = 1,ngeo
          effmu_0v_allgeo(1:geo%visible_0v_allgeo(igeo),igeo) = (effmu_0(1:geo%visible_0v_allgeo(igeo))*effmu_v_allgeo(1:geo%visible_0v_allgeo(igeo),igeo)) / (effmu_0(1:geo%visible_0v_allgeo(igeo)) + effmu_v_allgeo(1:geo%visible_0v_allgeo(igeo),igeo))
          effmu_0v_allgeo(geo%visible_0v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
        enddo
      endif
      if (execution(2) .or. execution(3)) then
        do ilay = 1,nlay
          ! Invisible layers have zero effmu_0 or effmu_v and grd%mu is positive, so they
          ! will naturally evaluate to zero.
          tlc%effmu_p0(:,ilay) = (grd%mu*effmu_0(ilay)) / (grd%mu + effmu_0(ilay))
          tlc%effmu_m0(:,ilay) = (grd%mu*effmu_0(ilay)) / (grd%mu - effmu_0(ilay))
          do igeo = 1,ngeo
            tlc%effmu_pv_allgeo(:,ilay,igeo) = (grd%mu*effmu_v_allgeo(ilay,igeo)) / (grd%mu + effmu_v_allgeo(ilay,igeo))
            tlc%effmu_mv_allgeo(:,ilay,igeo) = (grd%mu*effmu_v_allgeo(ilay,igeo)) / (grd%mu - effmu_v_allgeo(ilay,igeo))
          enddo
        enddo
      endif

    endif

    ! Create constructed integrated T-parameter for solar and viewing angle,
    ! for single-scattering geometry.
    if (execution(1) .or. execution(2)) then
      do igeo = 1,ngeo
        tlc%t_tot_0v_allgeo(:,igeo) = tlc%t_tot_0 * tlc%t_tot_v_allgeo(:,igeo)
      enddo
    endif

    ! Term necessary for non-split transmisstion operator, appearing in double-scattering.
    if (execution(2)) then
      do ilay = 1,nlay
        tlc%t(:,ilay) = expProtect(-pix%tau(ilay)/grd%mu)
      enddo
    endif

    ! Split environment, only for multi-scattering.
    if (execution(3)) then
      ! Split environment: Only single-layer integrals are needed.
      ! The cumulative transmission terms cannot be defined for split
      ! environment on nlay layers, because cumulative transmission terms
      ! to different sublayers are different. To acquire the cumulative
      ! transmission terms on the split environment, the cumulative
      ! transmission terms on the original grid need to be combined with
      ! the single-layer transmission terms on the split grid.

      ! Internal layers.
      do ilay = 1,nlay
        if (execution(2) .and. pix%nsplit(ilay) .eq. 1) then
          tlc%t_sp(:,ilay) = tlc%t(:,ilay)
        else
          tlc%t_sp(:,ilay) = expProtect(-(tau_sp(ilay))/grd%mu)
        endif
      enddo

      ! Sun.
      do ilay = 1,geo%visible_0
        if (pix%nsplit(ilay) .eq. 1) then
          tlc%t_sp_0(ilay) = t_0(ilay)
        else
          tlc%t_sp_0(ilay) = expProtect(-(tau_sp(ilay))/effmu_0(ilay))
        endif
      enddo
      tlc%t_sp_0(geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
      ! Instruments.
      do igeo = 1,ngeo
        do ilay = 1,geo%visible_v_allgeo(igeo)
          if (pix%nsplit(ilay) .eq. 1) then
            tlc%t_sp_v_allgeo(ilay,igeo) = t_v_allgeo(ilay,igeo)
          else
            tlc%t_sp_v_allgeo(ilay,igeo) = expProtect(-(tau_sp(ilay))/effmu_v_allgeo(ilay,igeo))
          endif
        enddo
        tlc%t_sp_v_allgeo(geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
      enddo
    endif

    ! Layer-integrals (L).

    ! L-parameter from sun to instrument.
    if (execution(1) .or. execution(2)) then
      do igeo = 1,ngeo
        if (geo%plane_parallel) then
          tlc%l_0v_allgeo(:,igeo) = geo%effmu_0v_ppg_allgeo(igeo) * (1.0-t_0*t_v_allgeo(:,igeo))
        else
          tlc%l_0v_allgeo(1:geo%visible_0v_allgeo(igeo),igeo) = effmu_0v_allgeo(1:geo%visible_0v_allgeo(igeo),igeo)*(1.0-t_0(1:geo%visible_0v_allgeo(igeo))*t_v_allgeo(1:geo%visible_0v_allgeo(igeo),igeo))
          tlc%l_0v_allgeo(geo%visible_0v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
        endif
      enddo
      if (pix%thermal_emission) then
        do igeo = 1,ngeo
          if (geo%plane_parallel) then
            tlc%l_v_allgeo(:,igeo) = geo%mu_v_allgeo(igeo) * (1.0 - t_v_allgeo(:,igeo))
          else
            tlc%l_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) = effmu_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) * (1.0-t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo))
            tlc%l_v_allgeo(geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0
          endif
        enddo
      endif
    endif

    ! L-parameters with one special direction on origianl grid.
    if (execution(2)) then
      ! Sun.
      do ilay = 1,geo%visible_0
        tlc%l_p0(:,ilay) = tlc%effmu_p0(:,ilay)*(1.0-t_0(ilay)*tlc%t(:,ilay))
        tlc%l_m0(:,ilay) = tlc%effmu_m0(:,ilay)*(1.0-divProtect(t_0(ilay) , tlc%t(:,ilay)))
      enddo
      tlc%l_p0(:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
      tlc%l_m0(:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
      ! Instruments.
      do igeo = 1,ngeo
        do ilay = 1,geo%visible_v_allgeo(igeo)
          tlc%l_pv_allgeo(:,ilay,igeo) = tlc%effmu_pv_allgeo(:,ilay,igeo)*(1.0-t_v_allgeo(ilay,igeo)*tlc%t(:,ilay))
          tlc%l_mv_allgeo(:,ilay,igeo) = tlc%effmu_mv_allgeo(:,ilay,igeo)*(1.0-divProtect(t_v_allgeo(ilay,igeo) , tlc%t(:,ilay)))
        enddo
        tlc%l_pv_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
        tlc%l_mv_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
      enddo
      ! Single-direction L-term for thermal emission.
      if (pix%thermal_emission) then
        do ilay = 1,nlay
          ! We already have the single-layer transmission term.
          tlc%l_p(:,ilay) = grd%mu * (1.0-tlc%t(:,ilay))
        enddo
      endif
    endif

    if (execution(3)) then

      ! L-parameters for the split environment.

      ! Calculate the intermediate results, pure split L-terms with one direction.
      ! The one on the internal grid is also needed for thermal emission, so it exists
      ! in the TLC structure.

      ! We refrain from recycling non-split parameters, because it is no exponential form
      ! and only l_v can be recycled under the condition that execution 1 or 2 and thermal
      ! emission is turned on. Those are many conditions.
      l_sp_0 = effmu_0 * (1.0-tlc%t_sp_0)
      do igeo = 1,ngeo
        l_sp_v_allgeo(:,igeo) = effmu_v_allgeo(:,igeo) * (1.0-tlc%t_sp_v_allgeo(:,igeo))
      enddo
      do ilay = 1,nlay
        tlc%l_sp_p(:,ilay) = grd%mu * (1.0-tlc%t_sp(:,ilay))
      enddo
      if (split_double .and. flag_derivatives) then
        do ilay = 1,nlay
          l_sp_m(:,ilay) = -grd%mu * (1.0-divProtect(1.0 , tlc%t_sp(:,ilay)))
        enddo
      endif

      ! Internal.
      if (split_double) then
        do ilay = 1,nlay
          do istr = 1,nstrhalf
            tlc%l_sp_pp(:,istr,ilay) = tlc%effmu_pp(:,istr)*(1.0-tlc%t_sp(istr,ilay)*tlc%t_sp(:,ilay))
            tlc%l_sp_mp(:,istr,ilay) = tlc%effmu_mp(:,istr)*(1.0-divProtect(tlc%t_sp(istr,ilay) , tlc%t_sp(:,ilay)))
          enddo
        enddo
        ! Correct 0/0 tlc%l_sp_mp to what it should be, tau-split. Note that
        ! the value of l_sp_mp is not evaluated as 0/0, because effmu_mp has
        ! been protected against zero division.
        do istr = 1,nstrhalf
          tlc%l_sp_mp(istr,istr,:) = tau_sp
        enddo
      endif

      ! Sun.
      do ilay = 1,geo%visible_0
        if (execution(2) .and. pix%nsplit(ilay) .eq. 1) then
          tlc%l_sp_p0(:,ilay) = tlc%l_p0(:,ilay)
          tlc%l_sp_m0(:,ilay) = tlc%l_m0(:,ilay)
        else
          tlc%l_sp_p0(:,ilay) = tlc%effmu_p0(:,ilay)*(1.0-tlc%t_sp_0(ilay)*tlc%t_sp(:,ilay))
          tlc%l_sp_m0(:,ilay) = tlc%effmu_m0(:,ilay)*(1.0-divProtect(tlc%t_sp_0(ilay) , tlc%t_sp(:,ilay)))
        endif
        ! Nothing to recycle for this strange L-parameter.
        if (split_double) tlc%l_1_sp_p0(:,ilay) = tlc%effmu_p0(:,ilay)*(tlc%l_sp_p0(:,ilay) - tau_sp(ilay)*tlc%t_sp_0(ilay)*tlc%t_sp(:,ilay)) ! A strange L-parameter.
      enddo
      tlc%l_sp_p0(:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
      tlc%l_sp_m0(:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
      if (split_double) tlc%l_1_sp_p0(:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
      ! Instruments.
      do igeo = 1,ngeo
        do ilay = 1,geo%visible_v_allgeo(igeo)
          if (execution(2) .and. pix%nsplit(ilay) .eq. 1) then
            tlc%l_sp_pv_allgeo(:,ilay,igeo) = tlc%l_pv_allgeo(:,ilay,igeo)
            tlc%l_sp_mv_allgeo(:,ilay,igeo) = tlc%l_mv_allgeo(:,ilay,igeo)
          else
            tlc%l_sp_pv_allgeo(:,ilay,igeo) = tlc%effmu_pv_allgeo(:,ilay,igeo)*(1.0-tlc%t_sp_v_allgeo(ilay,igeo)*tlc%t_sp(:,ilay))
            tlc%l_sp_mv_allgeo(:,ilay,igeo) = tlc%effmu_mv_allgeo(:,ilay,igeo)*(1.0-divProtect(tlc%t_sp_v_allgeo(ilay,igeo) , tlc%t_sp(:,ilay)))
          endif
          ! Nothing to recycle for this strange L-parameter.
          if (split_double) tlc%l_1_sp_pv_allgeo(:,ilay,igeo) = tlc%effmu_pv_allgeo(:,ilay,igeo)*(tlc%l_sp_pv_allgeo(:,ilay,igeo) - tau_sp(ilay)*tlc%t_sp_v_allgeo(ilay,igeo)*tlc%t_sp(:,ilay)) ! A strange L-parameter.
        enddo
        tlc%l_sp_pv_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
        tlc%l_sp_mv_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
        if (split_double) tlc%l_1_sp_pv_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
      enddo

      if (pix%thermal_emission) then
        if (split_double) then
          do ilay = 1,nlay
            tlc%l_1_sp_p(:,ilay) = grd%mu*(tlc%l_sp_p(:,ilay) - tau_sp(ilay)*tlc%t_sp(:,ilay)) ! A strange L-parameter.
          enddo
        endif
      endif

      ! L-parameters for interpolation.
      if (interpolation .eq. 0) then ! Averaging.
        ! Internal.
        tlc%l_int_sm_p = 0.5*tlc%l_sp_p
        tlc%l_int_ac_p = 0.5*tlc%l_sp_p
        if (split_double .and. flag_derivatives) then
          tlc%l_int_sm_m = 0.5*l_sp_m
          tlc%l_int_ac_m = 0.5*l_sp_m
        endif
        ! Sun.
        tlc%l_int_sm_0 = 0.5*l_sp_0
        tlc%l_int_ac_0 = 0.5*l_sp_0
        ! Instruments.
        tlc%l_int_sm_v_allgeo = 0.5*l_sp_v_allgeo
        tlc%l_int_ac_v_allgeo = 0.5*l_sp_v_allgeo
      endif
      if (interpolation .eq. 1) then ! Linear.
        ! Internal.
        do ilay = 1,nlay
          tlc%l_int_sm_p(:,ilay) = grd%mu * (1.0 - tlc%l_sp_p(:,ilay)/tau_sp(ilay))
          tlc%l_int_ac_p(:,ilay) = grd%mu * (tlc%l_sp_p(:,ilay)/tau_sp(ilay) - tlc%t_sp(:,ilay))
          if (split_double .and. flag_derivatives) then
            tlc%l_int_sm_m(:,ilay) = -grd%mu * (1.0 - l_sp_m(:,ilay)/tau_sp(ilay))
            tlc%l_int_ac_m(:,ilay) = -grd%mu * (l_sp_m(:,ilay)/tau_sp(ilay) - divProtect(1.0 , tlc%t_sp(:,ilay)))
          endif
        enddo
        ! Sun.
        tlc%l_int_sm_0 = effmu_0 * (1.0 - l_sp_0/tau_sp)
        tlc%l_int_ac_0 = effmu_0 * (l_sp_0/tau_sp - tlc%t_sp_0)
        ! Instruments.
        do igeo = 1,ngeo
          tlc%l_int_sm_v_allgeo(:,igeo) = effmu_v_allgeo(:,igeo) * (1.0 - l_sp_v_allgeo(:,igeo)/tau_sp)
          tlc%l_int_ac_v_allgeo(:,igeo) = effmu_v_allgeo(:,igeo) * (l_sp_v_allgeo(:,igeo)/tau_sp - tlc%t_sp_v_allgeo(:,igeo))
        enddo
      endif

      ! Interpolation from thermal emisison. The way of interpolation and the side do
      ! not matter.
      if (pix%thermal_emission) then
        tlc%l_int = 0.5 * tau_sp ! For interpolation, no matter the interpolation method or layer side.
      endif

    endif

    ! Tau-independent scatter parameters (C).

    ! Thermal emissivity. Actually, they do depend on Fourier number, but the C-parameter
    ! is always be the same for single-scattering geometry and for Fourier number 0 and
    ! always zero for higher Fourier numbers. So the C-parameter is calculated for single-
    ! scattering geometry and Fourier number 0, and the entire process will be ignored for
    ! higher Fourier numbers.
    if (pix%thermal_emission) tlc%c_therm = pix%planck_curve * (1.0-pix%ssa)

    ! From now on, the derivatives are calculated.

    if (flag_derivatives) then

      ! Administration of layers to which transmission terms are sensitive to tau
      ! Here the optional argument is used.
      if (size(ilay_deriv_base) .gt. 0) then
        do ilay = 0,nlay
          ind = maxloc(ilay_deriv_base,mask=ilay_deriv_base .le. geo%ilay_sensitive_0(ilay))
          tlc%ideriv_sensitive_0(ilay) = ind(1)
        enddo
        do igeo = 1,ngeo
          do ilay = 0,nlay
            ind = maxloc(ilay_deriv_base,mask=ilay_deriv_base .le. geo%ilay_sensitive_v_allgeo(ilay,igeo))
            tlc%ideriv_sensitive_v_allgeo(ilay,igeo) = ind(1)
          enddo
        enddo
        do igeo = 1,ngeo
          tlc%ideriv_sensitive_0v_allgeo(:,igeo) = max(tlc%ideriv_sensitive_0,tlc%ideriv_sensitive_v_allgeo(:,igeo))
        enddo
      endif

      ! These relative derivatives are needed for calculating the derivatives of any L for
      ! with a pseudo-spherical stream. Also, they are converted to the split environment,
      ! which is required as end product for multi-scattering no matter the geometry.
      if (.not. geo%plane_parallel .or. execution(3)) then

        ! Construct derivatives of single-layer transmission terms from the cumulative
        ! ones. This is how the pseudo-spherical geometry works. Under plane-parallel
        ! conditions, the geometry is set up such that is seems pseudo-spherical.
        reldiff_t_0_tau(:,1:geo%visible_0) = geo%reldiff_t_tot_0_tau(:,1:geo%visible_0) - geo%reldiff_t_tot_0_tau(:,0:geo%visible_0-1)
        do igeo = 1,ngeo
          reldiff_t_v_tau_allgeo(:,1:geo%visible_v_allgeo(igeo),igeo) = geo%reldiff_t_tot_v_tau_allgeo(:,1:geo%visible_v_allgeo(igeo),igeo) - geo%reldiff_t_tot_v_tau_allgeo(:,0:geo%visible_v_allgeo(igeo)-1,igeo)
        enddo

      endif

      ! Convert the temporary results to the desired end product for multi-scattering.
      ! With nsplit equal to 1, we could recycle, but a division by 1 is as easy as evaluating
      ! and if.
      if (execution(3)) then
        do ilay = 1,nlay
          ! Array operation over perturbed layers, applying chain rule.
          tlc%reldiff_t_sp_0_tau(:,ilay) = reldiff_t_0_tau(:,ilay) / pix%nsplit(ilay)
          ! Additional array operation over geometries.
          tlc%reldiff_t_sp_v_tau_allgeo(:,ilay,:) = reldiff_t_v_tau_allgeo(:,ilay,:) / pix%nsplit(ilay)
          ! Internal streams, also just a chain rule, but there, the original was already a
          ! permanent TLC-parameter that always exists, because it is atmosphere
          ! independent.
          tlc%reldiff_t_sp_tau(:,ilay) = tlc%reldiff_t_tau / pix%nsplit(ilay)
        enddo

      endif

      ! Derivatives of L-parameters.
      ! First those that only occur for thermal emission.
      if (execution(1) .or. execution(2)) then
      endif

      ! Normal L-parameters.
      if (execution(1) .or. execution(2)) then
        if (geo%plane_parallel) then
          tlc%diff_l_0v_tau_allgeo = 0.0
          do igeo = 1,ngeo
            do ilay = 1,nlay
              tlc%diff_l_0v_tau_allgeo(ilay,ilay,igeo) = t_0(ilay)*t_v_allgeo(ilay,igeo)
            enddo
          enddo
        else
          ! Derivatives of L in pseudo-psherical geometry are calculated with the five ingredients:
          ! - T
          ! - The logarithm of T
          ! - The relative derivative of T
          ! - L
          ! - Tau
          ! Then, the derivative of the effective direction cosine is not necessary
          ! We will refer to this method as the five-ingredient method. Combined with
          ! the straight-forward method, derivatives of effective direction cosines will
          ! be calculated later on.
          ! We only have L and tau, so the rest must be constructed.
          do igeo = 1,ngeo
            t_tmp(1:geo%visible_0v_allgeo(igeo)) = t_0(1:geo%visible_0v_allgeo(igeo)) * t_v_allgeo(1:geo%visible_0v_allgeo(igeo),igeo)
            log_t_tmp(1:geo%visible_0v_allgeo(igeo)) = log_t_0(1:geo%visible_0v_allgeo(igeo)) + log_t_v_allgeo(1:geo%visible_0v_allgeo(igeo),igeo)
            reldiff_t_tmp_tau(:,1:geo%visible_0v_allgeo(igeo)) = reldiff_t_0_tau(:,1:geo%visible_0v_allgeo(igeo)) + reldiff_t_v_tau_allgeo(:,1:geo%visible_0v_allgeo(igeo),igeo)
            do ilay = 1,geo%visible_0v_allgeo(igeo) ! Layer in which scattering takes place.
              ! Array operation is over perturbed layers.
              tlc%diff_l_0v_tau_allgeo(:,ilay,igeo) = divProtect(reldiff_t_tmp_tau(:,ilay) , log_t_tmp(ilay)) * (pix%tau(ilay)*t_tmp(ilay) - tlc%l_0v_allgeo(ilay,igeo))
              ! An additional diagonal term.
              tlc%diff_l_0v_tau_allgeo(ilay,ilay,igeo) = tlc%diff_l_0v_tau_allgeo(ilay,ilay,igeo) + tlc%l_0v_allgeo(ilay,igeo)/pix%tau(ilay)
            enddo
            tlc%diff_l_0v_tau_allgeo(:,geo%visible_0v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
          enddo
        endif
        ! L-parameter for thermal emission for execution 1 and 2.
        if (pix%thermal_emission) then
          if (geo%plane_parallel) then
            do ilay = 1,nlay
              tlc%diff_l_v_tau_allgeo(ilay,ilay,:) = t_v_allgeo(ilay,:)
            enddo
          else
            ! Derivatives of L in pseudo-psherical geometry are calculated with the five ingredients:
            ! - T
            ! - The logarithm of T
            ! - The relative derivative of T
            ! - L
            ! - Tau
            ! We have everything.
            do igeo = 1,ngeo
              do ilay = 1,geo%visible_v_allgeo(igeo) ! Layer in which scattering takes place.
                ! Array operation is over perturbed layers.
                tlc%diff_l_v_tau_allgeo(:,ilay,igeo) = divProtect(reldiff_t_v_tau_allgeo(:,ilay,igeo) , log_t_v_allgeo(ilay,igeo)) * (pix%tau(ilay)*t_v_allgeo(ilay,igeo) - tlc%l_v_allgeo(ilay,igeo))
                ! An additional diagonal term.
                tlc%diff_l_v_tau_allgeo(ilay,ilay,igeo) = tlc%diff_l_v_tau_allgeo(ilay,ilay,igeo) + tlc%l_v_allgeo(ilay,igeo)/pix%tau(ilay)
              enddo
              tlc%diff_l_v_tau_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
            enddo
          endif
        endif
      endif
      ! For the aggregation of a special (pseudo-spherical) and one internal streams
      ! the complexness of both sides are combined.
      ! We need an additional loop over streams (or play with fire with
      ! array operations, reshapes and possibly matrix multiplications, I do not know).
      if (execution(2)) then
        if (geo%plane_parallel) then
          ! All non-diagonal terms should be zero.
          tlc%diff_l_p0_tau = 0.0
          tlc%diff_l_m0_tau = 0.0
          tlc%diff_l_pv_tau_allgeo = 0.0
          tlc%diff_l_mv_tau_allgeo = 0.0
          do ilay = 1,nlay
            tlc%diff_l_p0_tau(:,ilay,ilay) = t_0(ilay) * tlc%t(:,ilay)
            tlc%diff_l_m0_tau(:,ilay,ilay) = divProtect(t_0(ilay) , tlc%t(:,ilay))
          enddo
          do igeo = 1,ngeo
            do ilay = 1,nlay
              tlc%diff_l_pv_tau_allgeo(:,ilay,ilay,igeo) = t_v_allgeo(ilay,igeo) * tlc%t(:,ilay)
              tlc%diff_l_mv_tau_allgeo(:,ilay,ilay,igeo) = divProtect(t_v_allgeo(ilay,igeo) , tlc%t(:,ilay))
            enddo
          enddo
        else
          ! Derivatives of L in pseudo-psherical geometry are calculated with the five ingredients:
          ! - T
          ! - The logarithm of T
          ! - The relative derivative of T
          ! - L
          ! - Tau
          ! Then, the derivative of the effective direction cosine is not necessary
          ! We will refer to this method as the five-ingredient method. Combined with
          ! the straight-forward method, derivatives of effective direction cosines will
          ! be calculated later on.
          ! We only have L and tau, so the rest must be constructed.
          do istr = 1,nstrhalf
            ! p0. That is the solar direction and forward transmission along internal stream.
            t_tmp(1:geo%visible_0) = t_0(1:geo%visible_0) * tlc%t(istr,1:geo%visible_0)
            log_t_tmp(1:geo%visible_0) = log_t_0(1:geo%visible_0) - pix%tau(1:geo%visible_0)/grd%mu(istr)
            reldiff_t_tmp_tau(:,1:geo%visible_0) = reldiff_t_0_tau(:,1:geo%visible_0)
            ! A diagonal term reldiff_t_tau must be added.
            do ilay = 1,geo%visible_0
              ! Add this diagonal term before it is too late.
              reldiff_t_tmp_tau(ilay,ilay) = reldiff_t_tmp_tau(ilay,ilay) + tlc%reldiff_t_tau(istr)
              ! Array operation is over perturbed layers.
              tlc%diff_l_p0_tau(istr,:,ilay) = divProtect(reldiff_t_tmp_tau(:,ilay) , log_t_tmp(ilay)) * (pix%tau(ilay)*t_tmp(ilay) - tlc%l_p0(istr,ilay))
              ! An additional diagonal term.
              tlc%diff_l_p0_tau(istr,ilay,ilay) = tlc%diff_l_p0_tau(istr,ilay,ilay) + tlc%l_p0(istr,ilay)/pix%tau(ilay)
            enddo
            ! m0. That is the solar direction and reverse transmission along internal stream.
            t_tmp(1:geo%visible_0) = divProtect(t_0(1:geo%visible_0) , tlc%t(istr,1:geo%visible_0))
            log_t_tmp(1:geo%visible_0) = log_t_0(1:geo%visible_0) + pix%tau(1:geo%visible_0)/grd%mu(istr)
            reldiff_t_tmp_tau(:,1:geo%visible_0) = reldiff_t_0_tau(:,1:geo%visible_0)
            ! A diagonal term reldiff_t_tau must be added.
            do ilay = 1,geo%visible_0
              ! Add this diagonal term before it is too late.
              reldiff_t_tmp_tau(ilay,ilay) = reldiff_t_tmp_tau(ilay,ilay) - tlc%reldiff_t_tau(istr)
              ! Array operation is over perturbed layers.
              tlc%diff_l_m0_tau(istr,:,ilay) = divProtect(reldiff_t_tmp_tau(:,ilay) , log_t_tmp(ilay)) * (pix%tau(ilay)*t_tmp(ilay) - tlc%l_m0(istr,ilay))
              ! An additional diagonal term.
              tlc%diff_l_m0_tau(istr,ilay,ilay) = tlc%diff_l_m0_tau(istr,ilay,ilay) + tlc%l_m0(istr,ilay)/pix%tau(ilay)
            enddo
          enddo
          tlc%diff_l_p0_tau(:,:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
          tlc%diff_l_m0_tau(:,:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
          ! The loop over streams has to be closed and opened again, because the
          ! loop over streams has to be outside.
          do igeo = 1,ngeo
            do istr = 1,nstrhalf
              ! pv. That is the viewing direction and forward transmission along internal stream.
              t_tmp(1:geo%visible_v_allgeo(igeo)) = t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) * tlc%t(istr,1:geo%visible_v_allgeo(igeo))
              log_t_tmp(1:geo%visible_v_allgeo(igeo)) = log_t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) - pix%tau(1:geo%visible_v_allgeo(igeo))/grd%mu(istr)
              reldiff_t_tmp_tau(:,1:geo%visible_v_allgeo(igeo)) = reldiff_t_v_tau_allgeo(:,1:geo%visible_v_allgeo(igeo),igeo)
              ! A diagonal term reldiff_t_tau must be added.
              do ilay = 1,geo%visible_v_allgeo(igeo)
                ! Add this diagonal term before it is too late.
                reldiff_t_tmp_tau(ilay,ilay) = reldiff_t_tmp_tau(ilay,ilay) + tlc%reldiff_t_tau(istr)
                ! Array operation is over perturbed layers.
                tlc%diff_l_pv_tau_allgeo(istr,:,ilay,igeo) = divProtect(reldiff_t_tmp_tau(:,ilay) , log_t_tmp(ilay)) * (pix%tau(ilay)*t_tmp(ilay) - tlc%l_pv_allgeo(istr,ilay,igeo))
                ! An additional diagonal term.
                tlc%diff_l_pv_tau_allgeo(istr,ilay,ilay,igeo) = tlc%diff_l_pv_tau_allgeo(istr,ilay,ilay,igeo) + tlc%l_pv_allgeo(istr,ilay,igeo)/pix%tau(ilay)
              enddo
              ! mv. That is the viewing direction and reverse transmission along internal stream.
              t_tmp(1:geo%visible_v_allgeo(igeo)) = divProtect(t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) , tlc%t(istr,1:geo%visible_v_allgeo(igeo)))
              log_t_tmp(1:geo%visible_v_allgeo(igeo)) = log_t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) + pix%tau(1:geo%visible_v_allgeo(igeo))/grd%mu(istr)
              reldiff_t_tmp_tau(:,1:geo%visible_v_allgeo(igeo)) = reldiff_t_v_tau_allgeo(:,1:geo%visible_v_allgeo(igeo),igeo)
              ! A diagonal term reldiff_t_tau must be added.
              do ilay = 1,geo%visible_v_allgeo(igeo)
                ! Add this diagonal term before it is too late.
                reldiff_t_tmp_tau(ilay,ilay) = reldiff_t_tmp_tau(ilay,ilay) - tlc%reldiff_t_tau(istr)
                ! Array operation is over perturbed layers.
                tlc%diff_l_mv_tau_allgeo(istr,:,ilay,igeo) = divProtect(reldiff_t_tmp_tau(:,ilay) , log_t_tmp(ilay)) * (pix%tau(ilay)*t_tmp(ilay) - tlc%l_mv_allgeo(istr,ilay,igeo))
                ! An additional diagonal term.
                tlc%diff_l_mv_tau_allgeo(istr,ilay,ilay,igeo) = tlc%diff_l_mv_tau_allgeo(istr,ilay,ilay,igeo) + tlc%l_mv_allgeo(istr,ilay,igeo)/pix%tau(ilay)
              enddo
            enddo ! This was the loop over streams.
            tlc%diff_l_pv_tau_allgeo(:,:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
            tlc%diff_l_mv_tau_allgeo(:,:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
          enddo ! This was the loop over geometries.
        endif
      endif

      ! Split environment and matrix zippers.
      if (execution(3)) then

        ! L-parameters for zipping two interpolated fields, no derivatives
        ! themselves, but only necessary for calculating derivatives.
        ! The zipping L-parameters do not have a tau in the numerator, because
        ! the perturbation is over tau_s (tau*omega) and tau, and not tau and
        ! omega. The tau is lost to the chain rule for the perturbation of tau_s
        ! and the tau is lost in the perturbation of tau, because that is the only
        ! case that there is no integration over a scatter event.
        ! Note that we are actually differentiating over C and not over omega, but
        ! that comes down to the same integrals, but the derivatives of C are then
        ! chain rules.
        if (interpolation .eq. 0) then
          tlc%l_zip_sm = 1.0 / (4.0*pix%nsplit)
          tlc%l_zip_ac = 1.0 / (4.0*pix%nsplit)
        endif
        if (interpolation .eq. 1) then
          tlc%l_zip_sm = 1.0 / (3.0*pix%nsplit)
          tlc%l_zip_ac = 1.0 / (6.0*pix%nsplit)
        endif

        if (geo%plane_parallel) then
          ! All non-diagonal terms should be zero.
          tlc%diff_l_sp_p0_tau = 0.0
          tlc%diff_l_sp_m0_tau = 0.0
          tlc%diff_l_sp_pv_tau_allgeo = 0.0
          tlc%diff_l_sp_mv_tau_allgeo = 0.0
          do ilay = 1,nlay
            if (execution(2) .and. pix%nsplit(ilay) .eq. 1) then
              tlc%diff_l_sp_p0_tau(:,ilay,ilay) = tlc%diff_l_p0_tau(:,ilay,ilay)
              tlc%diff_l_sp_m0_tau(:,ilay,ilay) = tlc%diff_l_m0_tau(:,ilay,ilay)
            else
              tlc%diff_l_sp_p0_tau(:,ilay,ilay) = tlc%t_sp_0(ilay) * tlc%t_sp(:,ilay) / pix%nsplit(ilay)
              tlc%diff_l_sp_m0_tau(:,ilay,ilay) = divProtect(tlc%t_sp_0(ilay) , tlc%t_sp(:,ilay)) / pix%nsplit(ilay)
            endif
          enddo
          do igeo = 1,ngeo
            do ilay = 1,nlay
              if (execution(2) .and. pix%nsplit(ilay) .eq. 1) then
                tlc%diff_l_sp_pv_tau_allgeo(:,ilay,ilay,igeo) = tlc%diff_l_pv_tau_allgeo(:,ilay,ilay,igeo)
                tlc%diff_l_sp_mv_tau_allgeo(:,ilay,ilay,igeo) = tlc%diff_l_mv_tau_allgeo(:,ilay,ilay,igeo)
              else
                tlc%diff_l_sp_pv_tau_allgeo(:,ilay,ilay,igeo) = tlc%t_sp_v_allgeo(ilay,igeo) * tlc%t_sp(:,ilay) / pix%nsplit(ilay)
                tlc%diff_l_sp_mv_tau_allgeo(:,ilay,ilay,igeo) = divProtect(tlc%t_sp_v_allgeo(ilay,igeo) , tlc%t_sp(:,ilay)) / pix%nsplit(ilay)
              endif
            enddo
          enddo
        else
          ! Derivatives of L in pseudo-psherical geometry are calculated with the five ingredients:
          ! - T
          ! - The logarithm of T
          ! - The relative derivative of T
          ! - L
          ! - Tau
          ! Everything is in the split environment, except the differential in tau for the
          ! derivative. That should work.
          do istr = 1,nstrhalf
            ! p0. That is the solar direction and forward transmission along internal stream.
            t_tmp(1:geo%visible_0) = tlc%t_sp_0(1:geo%visible_0) * tlc%t_sp(istr,1:geo%visible_0)
            log_t_tmp(1:geo%visible_0) = (log_t_0(1:geo%visible_0) - pix%tau(1:geo%visible_0)/grd%mu(istr)) / pix%nsplit(1:geo%visible_0)
            do ilay = 1,geo%visible_0
              if (execution(2) .and. pix%nsplit(ilay) .eq. 1) then
                tlc%diff_l_sp_p0_tau(istr,:,ilay) = tlc%diff_l_p0_tau(istr,:,ilay)
              else
                ! Non-diagonal term.
                reldiff_t_tmp_tau(:,ilay) = reldiff_t_0_tau(:,ilay) / pix%nsplit(ilay)
                ! Add this diagonal term before it is too late.
                reldiff_t_tmp_tau(ilay,ilay) = reldiff_t_tmp_tau(ilay,ilay) + tlc%reldiff_t_sp_tau(istr,ilay)
                ! Array operation is over perturbed layers.
                tlc%diff_l_sp_p0_tau(istr,:,ilay) = divProtect(reldiff_t_tmp_tau(:,ilay) , log_t_tmp(ilay)) * (tau_sp(ilay)*t_tmp(ilay) - tlc%l_sp_p0(istr,ilay))
                ! An additional diagonal term.
                tlc%diff_l_sp_p0_tau(istr,ilay,ilay) = tlc%diff_l_sp_p0_tau(istr,ilay,ilay) + tlc%l_sp_p0(istr,ilay) / pix%tau(ilay)
              endif
            enddo
            ! m0. That is the solar direction and reverse transmission along internal stream.
            t_tmp(1:geo%visible_0) = divProtect(tlc%t_sp_0(1:geo%visible_0) , tlc%t_sp(istr,1:geo%visible_0))
            log_t_tmp(1:geo%visible_0) = (log_t_0(1:geo%visible_0) + pix%tau(1:geo%visible_0)/grd%mu(istr)) / pix%nsplit(1:geo%visible_0)
            do ilay = 1,geo%visible_0
              if (execution(2) .and. pix%nsplit(ilay) .eq. 1) then
                tlc%diff_l_sp_m0_tau(istr,:,ilay) = tlc%diff_l_m0_tau(istr,:,ilay)
              else
                ! Non-diagonal term.
                reldiff_t_tmp_tau(:,ilay) = reldiff_t_0_tau(:,ilay) / pix%nsplit(ilay)
                ! Add this diagonal term before it is too late.
                reldiff_t_tmp_tau(ilay,ilay) = reldiff_t_tmp_tau(ilay,ilay) - tlc%reldiff_t_sp_tau(istr,ilay)
                ! Array operation is over perturbed layers.
                tlc%diff_l_sp_m0_tau(istr,:,ilay) = divProtect(reldiff_t_tmp_tau(:,ilay) , log_t_tmp(ilay)) * (tau_sp(ilay)*t_tmp(ilay) - tlc%l_sp_m0(istr,ilay))
                ! An additional diagonal term.
                tlc%diff_l_sp_m0_tau(istr,ilay,ilay) = tlc%diff_l_sp_m0_tau(istr,ilay,ilay) + tlc%l_sp_m0(istr,ilay) / pix%tau(ilay)
              endif
            enddo
          enddo
          tlc%diff_l_sp_p0_tau(:,:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
          tlc%diff_l_sp_m0_tau(:,:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
          ! The loop over streams has to be closed, because we now need a loop over
          ! viewing geometries outside the loop over streams.
          do igeo = 1,ngeo
            do istr = 1,nstrhalf
              ! pv. That is the viewing direction and forward transmission along internal stream.
              t_tmp(1:geo%visible_v_allgeo(igeo)) = tlc%t_sp_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) * tlc%t_sp(istr,1:geo%visible_v_allgeo(igeo))
              log_t_tmp(1:geo%visible_v_allgeo(igeo)) = (log_t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) - pix%tau(1:geo%visible_v_allgeo(igeo))/grd%mu(istr)) / pix%nsplit(1:geo%visible_v_allgeo(igeo))
              do ilay = 1,geo%visible_v_allgeo(igeo)
                if (execution(2) .and. pix%nsplit(ilay) .eq. 1) then
                  tlc%diff_l_sp_pv_tau_allgeo(istr,:,ilay,igeo) = tlc%diff_l_pv_tau_allgeo(istr,:,ilay,igeo)
                else
                  ! Non-diagonal term.
                  reldiff_t_tmp_tau(:,ilay) = reldiff_t_v_tau_allgeo(:,ilay,igeo) / pix%nsplit(ilay)
                  ! Add this diagonal term before it is too late.
                  reldiff_t_tmp_tau(ilay,ilay) = reldiff_t_tmp_tau(ilay,ilay) + tlc%reldiff_t_sp_tau(istr,ilay)
                  ! Array operation is over perturbed layers.
                  tlc%diff_l_sp_pv_tau_allgeo(istr,:,ilay,igeo) = divProtect(reldiff_t_tmp_tau(:,ilay) , log_t_tmp(ilay)) * (tau_sp(ilay)*t_tmp(ilay) - tlc%l_sp_pv_allgeo(istr,ilay,igeo))
                  ! An additional diagonal term.
                  tlc%diff_l_sp_pv_tau_allgeo(istr,ilay,ilay,igeo) = tlc%diff_l_sp_pv_tau_allgeo(istr,ilay,ilay,igeo) + tlc%l_sp_pv_allgeo(istr,ilay,igeo) / pix%tau(ilay)
                endif
              enddo
              ! mv. That is the viewing direction and reverse transmission along internal stream.
              t_tmp(1:geo%visible_v_allgeo(igeo)) = divProtect(tlc%t_sp_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) , tlc%t_sp(istr,1:geo%visible_v_allgeo(igeo)))
              log_t_tmp(1:geo%visible_v_allgeo(igeo)) = (log_t_v_allgeo(1:geo%visible_v_allgeo(igeo),igeo) + pix%tau(1:geo%visible_v_allgeo(igeo))/grd%mu(istr)) / pix%nsplit(1:geo%visible_v_allgeo(igeo))
              do ilay = 1,geo%visible_v_allgeo(igeo)
                if (execution(2) .and. pix%nsplit(ilay) .eq. 1) then
                  tlc%diff_l_sp_mv_tau_allgeo(istr,:,ilay,igeo) = tlc%diff_l_mv_tau_allgeo(istr,:,ilay,igeo)
                else
                  ! Non-diagonal term.
                  reldiff_t_tmp_tau(:,ilay) = reldiff_t_v_tau_allgeo(:,ilay,igeo) / pix%nsplit(ilay)
                  ! Add this diagonal term before it is too late.
                  reldiff_t_tmp_tau(ilay,ilay) = reldiff_t_tmp_tau(ilay,ilay) - tlc%reldiff_t_sp_tau(istr,ilay)
                  ! Array operation is over perturbed layers.
                  tlc%diff_l_sp_mv_tau_allgeo(istr,:,ilay,igeo) = divProtect(reldiff_t_tmp_tau(:,ilay) , log_t_tmp(ilay)) * (tau_sp(ilay)*t_tmp(ilay) - tlc%l_sp_mv_allgeo(istr,ilay,igeo))
                  ! An additional diagonal term.
                  tlc%diff_l_sp_mv_tau_allgeo(istr,ilay,ilay,igeo) = tlc%diff_l_sp_mv_tau_allgeo(istr,ilay,ilay,igeo) + tlc%l_sp_mv_allgeo(istr,ilay,igeo) / pix%tau(ilay)
                endif
              enddo
            enddo
            tlc%diff_l_sp_pv_tau_allgeo(:,:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
            tlc%diff_l_sp_mv_tau_allgeo(:,:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
          enddo

        endif
      endif

      ! Derivatives of effective direction cosines for two directions of which
      ! one is the sun or an instrument.
      if (execution(2) .or. execution(3) .and. split_double) then
        if (geo%plane_parallel) then
          tlc%diff_effmu_p0_tau = 0.0
          tlc%diff_effmu_m0_tau = 0.0
          tlc%diff_effmu_pv_tau_allgeo = 0.0
          tlc%diff_effmu_mv_tau_allgeo = 0.0
        else
          ! Derivatives of effective direction are calculated with the
          ! derivatives of L. These derivatives have been calculated for the normal
          ! domain for execution 2 and for the split environment for execution 3,
          ! or possibly both if both are turned on. Then, it does not matter which we
          ! take.
          if (execution(3)) then
            ! Use split environment. If execution 2 was also true, it does not matter.
            do ilay = 1,geo%visible_0
              do istr = 1,nstrhalf
                tlc%diff_effmu_p0_tau(istr,:,ilay) = divProtect(tlc%diff_l_sp_p0_tau(istr,:,ilay) + tlc%effmu_p0(istr,ilay) * tlc%reldiff_t_sp_0_tau(:,ilay) * (tlc%t_sp_0(ilay)*tlc%t_sp(istr,ilay)) , 1.0 - tlc%t_sp_0(ilay)*tlc%t_sp(istr,ilay))
                ! Add diagonal term.
                tlc%diff_effmu_p0_tau(istr,ilay,ilay) = tlc%diff_effmu_p0_tau(istr,ilay,ilay) + divProtect(tlc%effmu_p0(istr,ilay) * tlc%reldiff_t_sp_tau(istr,ilay) * (tlc%t_sp_0(ilay)*tlc%t_sp(istr,ilay)) , 1.0 - tlc%t_sp_0(ilay)*tlc%t_sp(istr,ilay))

                tlc%diff_effmu_m0_tau(istr,:,ilay) = divProtect(tlc%diff_l_sp_m0_tau(istr,:,ilay) + tlc%effmu_m0(istr,ilay) * tlc%reldiff_t_sp_0_tau(:,ilay) * divProtect(tlc%t_sp_0(ilay) , tlc%t_sp(istr,ilay)) , 1.0 - divProtect(tlc%t_sp_0(ilay) , tlc%t_sp(istr,ilay)))
                ! Add diagonal term.
                tlc%diff_effmu_m0_tau(istr,ilay,ilay) = tlc%diff_effmu_m0_tau(istr,ilay,ilay) - divProtect(tlc%effmu_m0(istr,ilay) * tlc%reldiff_t_sp_tau(istr,ilay) * divProtect(tlc%t_sp_0(ilay) , tlc%t_sp(istr,ilay)) , 1.0 - divProtect(tlc%t_sp_0(ilay) , tlc%t_sp(istr,ilay)))
              enddo
            enddo
            do igeo = 1,ngeo
              do ilay = 1,geo%visible_v_allgeo(igeo)
                do istr = 1,nstrhalf
                  tlc%diff_effmu_pv_tau_allgeo(istr,:,ilay,igeo) = divProtect(tlc%diff_l_sp_pv_tau_allgeo(istr,:,ilay,igeo) + tlc%effmu_pv_allgeo(istr,ilay,igeo) * tlc%reldiff_t_sp_v_tau_allgeo(:,ilay,igeo) * (tlc%t_sp_v_allgeo(ilay,igeo)*tlc%t_sp(istr,ilay)) , 1.0 - tlc%t_sp_v_allgeo(ilay,igeo)*tlc%t_sp(istr,ilay))
                  ! Add diagonal terms.
                  tlc%diff_effmu_pv_tau_allgeo(istr,ilay,ilay,igeo) = tlc%diff_effmu_pv_tau_allgeo(istr,ilay,ilay,igeo) + divProtect(tlc%effmu_pv_allgeo(istr,ilay,igeo) * tlc%reldiff_t_sp_tau(istr,ilay) * (tlc%t_sp_v_allgeo(ilay,igeo)*tlc%t_sp(istr,ilay)) , 1.0 - tlc%t_sp_v_allgeo(ilay,igeo)*tlc%t_sp(istr,ilay))
                  tlc%diff_effmu_mv_tau_allgeo(istr,:,ilay,igeo) = divProtect(tlc%diff_l_sp_mv_tau_allgeo(istr,:,ilay,igeo) + tlc%effmu_mv_allgeo(istr,ilay,igeo) * tlc%reldiff_t_sp_v_tau_allgeo(:,ilay,igeo) * divProtect(tlc%t_sp_v_allgeo(ilay,igeo) , tlc%t_sp(istr,ilay)) , 1.0 - divProtect(tlc%t_sp_v_allgeo(ilay,igeo) , tlc%t_sp(istr,ilay)))
                  ! Add diagonal terms.
                  tlc%diff_effmu_mv_tau_allgeo(istr,ilay,ilay,igeo) = tlc%diff_effmu_mv_tau_allgeo(istr,ilay,ilay,igeo) - divProtect(tlc%effmu_mv_allgeo(istr,ilay,igeo) * tlc%reldiff_t_sp_tau(istr,ilay) * divProtect(tlc%t_sp_v_allgeo(ilay,igeo) , tlc%t_sp(istr,ilay)) , 1.0 - divProtect(tlc%t_sp_v_allgeo(ilay,igeo) , tlc%t_sp(istr,ilay)))
                enddo
              enddo
            enddo
          else
            ! Use non-split enviroment.
            do ilay = 1,geo%visible_0
              do istr = 1,nstrhalf
                tlc%diff_effmu_p0_tau(istr,:,ilay) = divProtect(tlc%diff_l_p0_tau(istr,:,ilay) + tlc%effmu_p0(istr,ilay) * reldiff_t_0_tau(:,ilay) * (t_0(ilay)*tlc%t(istr,ilay)) , 1.0 - t_0(ilay)*tlc%t(istr,ilay))
                ! Add diagonal terms.
                tlc%diff_effmu_p0_tau(istr,ilay,ilay) = tlc%diff_effmu_p0_tau(istr,ilay,ilay) + divProtect(tlc%effmu_p0(istr,ilay) * tlc%reldiff_t_tau(istr) * (t_0(ilay)*tlc%t(istr,ilay)) , 1.0 - t_0(ilay)*tlc%t(istr,ilay))

                tlc%diff_effmu_m0_tau(istr,:,ilay) = divProtect(tlc%diff_l_m0_tau(istr,:,ilay) + tlc%effmu_m0(istr,ilay) * reldiff_t_0_tau(:,ilay) * divProtect(t_0(ilay) , tlc%t(istr,ilay)) , 1.0 - divProtect(t_0(ilay) , tlc%t(istr,ilay)))
                ! Add diagonal terms.
                tlc%diff_effmu_m0_tau(istr,ilay,ilay) = tlc%diff_effmu_m0_tau(istr,ilay,ilay) - divProtect(tlc%effmu_m0(istr,ilay) * tlc%reldiff_t_tau(istr) * divProtect(t_0(ilay) , tlc%t(istr,ilay)) , 1.0 - divProtect(t_0(ilay) , tlc%t(istr,ilay)))
              enddo
            enddo
            do igeo = 1,ngeo
              do ilay = 1,geo%visible_v_allgeo(igeo)
                do istr = 1,nstrhalf
                  tlc%diff_effmu_pv_tau_allgeo(istr,:,ilay,igeo) = divProtect(tlc%diff_l_pv_tau_allgeo(istr,:,ilay,igeo) + tlc%effmu_pv_allgeo(istr,ilay,igeo) * reldiff_t_v_tau_allgeo(:,ilay,igeo) * (t_v_allgeo(ilay,igeo)*tlc%t(istr,ilay)) , 1.0 - t_v_allgeo(ilay,igeo)*tlc%t(istr,ilay))
                  ! Add diagonal terms.
                  tlc%diff_effmu_pv_tau_allgeo(istr,ilay,ilay,igeo) = tlc%diff_effmu_pv_tau_allgeo(istr,ilay,ilay,igeo) + divProtect(tlc%effmu_pv_allgeo(istr,ilay,igeo) * tlc%reldiff_t_tau(istr) * (t_v_allgeo(ilay,igeo)*tlc%t(istr,ilay)) , 1.0 - t_v_allgeo(ilay,igeo)*tlc%t(istr,ilay))

                  tlc%diff_effmu_mv_tau_allgeo(istr,:,ilay,igeo) = divProtect(tlc%diff_l_mv_tau_allgeo(istr,:,ilay,igeo) + tlc%effmu_mv_allgeo(istr,ilay,igeo) * reldiff_t_v_tau_allgeo(:,ilay,igeo) * divProtect(t_v_allgeo(ilay,igeo) , tlc%t(istr,ilay)) , 1.0 - divProtect(t_v_allgeo(ilay,igeo) , tlc%t(istr,ilay)))
                  ! Add diagonal terms.
                  tlc%diff_effmu_mv_tau_allgeo(istr,ilay,ilay,igeo) = tlc%diff_effmu_mv_tau_allgeo(istr,ilay,ilay,igeo) - divProtect(tlc%effmu_mv_allgeo(istr,ilay,igeo) * tlc%reldiff_t_tau(istr) * divProtect(t_v_allgeo(ilay,igeo) , tlc%t(istr,ilay)) , 1.0 - divProtect(t_v_allgeo(ilay,igeo) , tlc%t(istr,ilay)))
                enddo
              enddo
            enddo
          endif
          tlc%diff_effmu_p0_tau(:,:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
          tlc%diff_effmu_m0_tau(:,:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
          do igeo = 1,ngeo
            tlc%diff_effmu_pv_tau_allgeo(:,:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
            tlc%diff_effmu_mv_tau_allgeo(:,:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
          enddo
        endif
      endif

      if (execution(3)) then

        ! The intermediate derivatives of L in split enviroment and the derivatives of the
        ! effective direction cosines are needed.

        ! Derivatives of L-plus and L-min.
        do ilay = 1,nlay
          tlc%diff_l_sp_p_tau(:,ilay) = tlc%t_sp(:,ilay) / pix%nsplit(ilay)
        enddo
        do ilay = 1,nlay
          diff_l_sp_m_tau(:,ilay) = divProtect(1.0 , tlc%t_sp(:,ilay)) / pix%nsplit(ilay)
        enddo
        ! Use intermediate results to calculate derivatives
        ! Declare intermediate results.
        if (geo%plane_parallel) then
          diff_l_sp_0_tau = 0.0
          do ilay = 1,nlay
            diff_l_sp_0_tau(ilay,ilay) = tlc%t_sp_0(ilay) / pix%nsplit(ilay)
          enddo
          diff_l_sp_v_tau_allgeo = 0.0
          do igeo = 1,ngeo
            do ilay = 1,nlay
              diff_l_sp_v_tau_allgeo(ilay,ilay,igeo) = tlc%t_sp_v_allgeo(ilay,igeo) / pix%nsplit(ilay)
            enddo
          enddo
          ! and the derivatives of the effective direction cosines.
          diff_effmu_0_tau = 0.0
          diff_effmu_v_tau_allgeo = 0.0
        else
          ! Calculate the derivative of L with the pseudo-spherical formula.
          log_t_tmp = log_t_0 / pix%nsplit
          ! Do the calculation.
          do ilay = 1,geo%visible_0
            diff_l_sp_0_tau(:,ilay) = divProtect(tlc%reldiff_t_sp_0_tau(:,ilay) , log_t_tmp(ilay)) * (tau_sp(ilay)*tlc%t_sp_0(ilay) - l_sp_0(ilay))
            ! An additional diagonal term.
            diff_l_sp_0_tau(ilay,ilay) = diff_l_sp_0_tau(ilay,ilay) + l_sp_0(ilay) / pix%tau(ilay)
            ! Derive derivative of mu from the derivative of L.
            diff_effmu_0_tau(:,ilay) = divProtect(diff_l_sp_0_tau(:,ilay) + effmu_0(ilay) * tlc%reldiff_t_sp_0_tau(:,ilay) * tlc%t_sp_0(ilay) , 1.0 - tlc%t_sp_0(ilay))
          enddo
          diff_l_sp_0_tau(:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
          diff_effmu_0_tau(:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
          ! Instruments.
          do igeo = 1,ngeo
            ! Calculate the derivative of L with the pseudo-spherical formula.
            log_t_tmp = log_t_v_allgeo(:,igeo) / pix%nsplit
            ! Do the calculation.
            do ilay = 1,geo%visible_v_allgeo(igeo)
              diff_l_sp_v_tau_allgeo(:,ilay,igeo) = divProtect(tlc%reldiff_t_sp_v_tau_allgeo(:,ilay,igeo) , log_t_tmp(ilay)) * (tau_sp(ilay)*tlc%t_sp_v_allgeo(ilay,igeo) - l_sp_v_allgeo(ilay,igeo))
              ! An additional diagonal term.
              diff_l_sp_v_tau_allgeo(ilay,ilay,igeo) = diff_l_sp_v_tau_allgeo(ilay,ilay,igeo) + l_sp_v_allgeo(ilay,igeo) / pix%tau(ilay)
              ! Derive derivative of mu from the derivative of L.
              diff_effmu_v_tau_allgeo(:,ilay,igeo) = divProtect(diff_l_sp_v_tau_allgeo(:,ilay,igeo) + effmu_v_allgeo(ilay,igeo) * tlc%reldiff_t_sp_v_tau_allgeo(:,ilay,igeo) * tlc%t_sp_v_allgeo(ilay,igeo) , 1.0 - tlc%t_sp_v_allgeo(ilay,igeo))
            enddo
            diff_l_sp_v_tau_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
            diff_effmu_v_tau_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
          enddo
        endif

        ! Now, use these intermediate results to calculate the derivatives of interpolated
        ! L-parameters.
        if (interpolation .eq. 0) then
          ! Internal streams.
          if (split_double) then
            tlc%diff_l_int_sm_p_tau = 0.5 * tlc%diff_l_sp_p_tau
            tlc%diff_l_int_ac_p_tau = 0.5 * tlc%diff_l_sp_p_tau
            tlc%diff_l_int_sm_m_tau = 0.5 * diff_l_sp_m_tau
            tlc%diff_l_int_ac_m_tau = 0.5 * diff_l_sp_m_tau
          endif
          ! Sun.
          tlc%diff_l_int_sm_0_tau(:,1:geo%visible_0) = 0.5 * diff_l_sp_0_tau(:,1:geo%visible_0)
          tlc%diff_l_int_ac_0_tau(:,1:geo%visible_0) = 0.5 * diff_l_sp_0_tau(:,1:geo%visible_0)
          ! Instruments.
          do igeo = 1,ngeo
            tlc%diff_l_int_sm_v_tau_allgeo(:,1:geo%visible_v_allgeo(igeo),igeo) = 0.5 * diff_l_sp_v_tau_allgeo(:,1:geo%visible_v_allgeo(igeo),igeo)
            tlc%diff_l_int_ac_v_tau_allgeo(:,1:geo%visible_v_allgeo(igeo),igeo) = 0.5 * diff_l_sp_v_tau_allgeo(:,1:geo%visible_v_allgeo(igeo),igeo)
          enddo
        endif
        if (interpolation .eq. 1) then
          if (split_double) then
            do ilay = 1,nlay
              ! Internal.
              ! An nsplit*tau_sp is not converted to tau, we just use the derivative of tau_sp
              ! to tau, which is 1/nsplit and not use any non-split variables.
              tlc%diff_l_int_sm_p_tau(:,ilay) = grd%mu / tau_sp(ilay) * (tlc%l_sp_p(:,ilay)/(pix%nsplit(ilay)*tau_sp(ilay)) - tlc%diff_l_sp_p_tau(:,ilay))
              tlc%diff_l_int_ac_p_tau(:,ilay) = grd%mu / tau_sp(ilay) * (tlc%diff_l_sp_p_tau(:,ilay) - tlc%l_sp_p(:,ilay)/(pix%nsplit(ilay)*tau_sp(ilay))) - grd%mu * tlc%reldiff_t_sp_tau(:,ilay) * tlc%t_sp(:,ilay)
              tlc%diff_l_int_sm_m_tau(:,ilay) = -grd%mu / tau_sp(ilay) * (l_sp_m(:,ilay)/(pix%nsplit(ilay)*tau_sp(ilay)) - diff_l_sp_m_tau(:,ilay))
              tlc%diff_l_int_ac_m_tau(:,ilay) = -grd%mu / tau_sp(ilay) * (diff_l_sp_m_tau(:,ilay) - l_sp_m(:,ilay)/(pix%nsplit(ilay)*tau_sp(ilay))) - grd%mu * divProtect(tlc%reldiff_t_sp_tau(:,ilay) , tlc%t_sp(:,ilay))
            enddo
          endif
          ! Sun.
          do ilay = 1,geo%visible_0
            ! Differentiate the interpolated L-parameters with product and quotient rules:
            ! l_int_sm_0 = effmu_0*(1 - L/tau)
            ! Her derivative has three parts, from differentiating effmu, L and tau:
            ! diff_effmu_0 * (1 - L/tau)
            ! -diff_L * effmu_0/tau
            ! effmu_0 * L / tau**2.0 (only diagnoal with nsplit chain rule)
            tlc%diff_l_int_sm_0_tau(:,ilay) = diff_effmu_0_tau(:,ilay) * (1.0 - l_sp_0(ilay) / tau_sp(ilay)) - diff_l_sp_0_tau(:,ilay) * effmu_0(ilay) / tau_sp(ilay)
            tlc%diff_l_int_sm_0_tau(ilay,ilay) = tlc%diff_l_int_sm_0_tau(ilay,ilay) + effmu_0(ilay) * l_sp_0(ilay) / (tau_sp(ilay)**2.0 * pix%nsplit(ilay))
            ! Differentiate the interpolated L-parameters with product and quotient rules:
            ! l_int_ac_0 = effmu_0*(L/tau - T)
            ! Her derivative has four parts, from differentiating effmu, L, tau and T:
            ! diff_effmu_0 * (L/tau - T)
            ! diff_L * effmu_0/tau
            ! -effmu_0 * L / tau**2.0 (only diagnoal with nsplit chain rule)
            ! -T*reldiff_T * effmu_0
            tlc%diff_l_int_ac_0_tau(:,ilay) = diff_effmu_0_tau(:,ilay) * (l_sp_0(ilay) / tau_sp(ilay) - tlc%t_sp_0(ilay)) + diff_l_sp_0_tau(:,ilay) * effmu_0(ilay)/tau_sp(ilay) - tlc%t_sp_0(ilay) * tlc%reldiff_t_sp_0_tau(:,ilay) * effmu_0(ilay)
            tlc%diff_l_int_ac_0_tau(ilay,ilay) = tlc%diff_l_int_ac_0_tau(ilay,ilay) - effmu_0(ilay) * l_sp_0(ilay) / (tau_sp(ilay)**2.0 * pix%nsplit(ilay))
          enddo
          ! Instruments.
          do igeo = 1,ngeo
            do ilay = 1,geo%visible_v_allgeo(igeo)
              ! Differentiate the interpolated L-parameters with product and quotient rules:
              ! l_int_sm_v = effmu_v*(1 - L/tau)
              ! Her derivative has three parts, from differentiating effmu, L and tau:
              ! diff_effmu_v * (1 - L/tau)
              ! -diff_L * effmu_v/tau
              ! effmu_v * L / tau**2.0 (only diagnoal with nsplit chain rule)
              tlc%diff_l_int_sm_v_tau_allgeo(:,ilay,igeo) = diff_effmu_v_tau_allgeo(:,ilay,igeo) * (1.0 - l_sp_v_allgeo(ilay,igeo) / tau_sp(ilay)) - diff_l_sp_v_tau_allgeo(:,ilay,igeo) * effmu_v_allgeo(ilay,igeo) / tau_sp(ilay)
              tlc%diff_l_int_sm_v_tau_allgeo(ilay,ilay,igeo) = tlc%diff_l_int_sm_v_tau_allgeo(ilay,ilay,igeo) + effmu_v_allgeo(ilay,igeo) * l_sp_v_allgeo(ilay,igeo) / (tau_sp(ilay)**2.0 * pix%nsplit(ilay))
              ! Differentiate the interpolated L-parameters with product and quotient rules:
              ! l_int_ac_v = effmu_v*(L/tau - T)
              ! Her derivative has four parts, from differentiating effmu, L, tau and T:
              ! diff_effmu_v * (L/tau - T)
              ! diff_L * effmu_v/tau
              ! -effmu_v * L / tau**2.0 (only diagnoal with nsplit chain rule)
              ! -T*reldiff_T * effmu_v
              tlc%diff_l_int_ac_v_tau_allgeo(:,ilay,igeo) = diff_effmu_v_tau_allgeo(:,ilay,igeo) * (l_sp_v_allgeo(ilay,igeo) / tau_sp(ilay) - tlc%t_sp_v_allgeo(ilay,igeo)) + diff_l_sp_v_tau_allgeo(:,ilay,igeo) * effmu_v_allgeo(ilay,igeo)/tau_sp(ilay) - tlc%t_sp_v_allgeo(ilay,igeo) * tlc%reldiff_t_sp_v_tau_allgeo(:,ilay,igeo) * effmu_v_allgeo(ilay,igeo)
              tlc%diff_l_int_ac_v_tau_allgeo(ilay,ilay,igeo) = tlc%diff_l_int_ac_v_tau_allgeo(ilay,ilay,igeo) - effmu_v_allgeo(ilay,igeo) * l_sp_v_allgeo(ilay,igeo) / (tau_sp(ilay)**2.0 * pix%nsplit(ilay))
            enddo
          enddo
        endif
        tlc%diff_l_int_sm_0_tau(:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
        tlc%diff_l_int_ac_0_tau(:,geo%visible_0+1:nlay) = 0.0 ! Invisible layers.
        do igeo = 1,ngeo
          tlc%diff_l_int_sm_v_tau_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
          tlc%diff_l_int_ac_v_tau_allgeo(:,geo%visible_v_allgeo(igeo)+1:nlay,igeo) = 0.0 ! Invisible layers.
        enddo
        if (pix%thermal_emission) then
          ! The derivative diff_l_sp_p_tau will be calculated anyway as intermediate
          ! result whenever execution 3 is on.
          tlc%diff_l_int_tau = 0.5 / pix%nsplit ! The interpolation
          ! The strange L_1-parameter need not be differentiated.

        endif
      endif

      ! Derivatives of thermal emissivity. Actually, they do depend on Fourier number, but the
      ! C-parameter is always be the same for single-scattering geometry and for Fourier number 0
      ! and always zero for higher Fourier numbers. So the C-parameter is calculated for single-
      ! scattering geometry and Fourier number 0, and the entire process will be ignored for
      ! higher Fourier numbers.
      if (pix%thermal_emission) then
        tlc%diff_c_therm_ssa = -pix%planck_curve
        tlc%diff_c_therm_planck = 1.0-pix%ssa
      endif

    endif

  end subroutine calculate_fourier_independent_tlc_parameters ! }}}

  !> Calculate the C-parameters for single-scattering in the relevant geometry.
  !!
  !! This relevant geometry is both a viewing angle a viewing Stokes parameter.
  !! This way, single-scattering is done consistently with multi-scattering, where
  !! all viewing geometries and Stokes parameters are handled separately.
  subroutine calculate_relevant_single_scattering_c_parameters(igeo,ist_src,ist,geo,pix,tlc) ! {{{

    use lintran_constants_module, only: ist_i, ist_q, ist_u, ist_v, iscat_a1, iscat_b1, iscat_a2, iscat_a3, iscat_b2, iscat_a4
    use lintran_types_module, only: lintran_geometry_class, lintran_pixel_class, lintran_tlc_class

    implicit none

    ! Input and output
    integer, intent(in) :: igeo !< Viewing geometry index.
    integer, intent(in) :: ist_src !< Source Stokes parameter.
    integer, intent(in) :: ist !< Viewing Stokes parameter.
    type(lintran_geometry_class), intent(in) :: geo !< Internal geometry structure.
    type(lintran_pixel_class), intent(in) :: pix !< Internal pixel structure.
    type(lintran_tlc_class), intent(inout) :: tlc !< Internal radiative-transfer building-block structure.

    ! The phase matrix is the derivative of C with respect to the single-scattering
    ! albedo. Therefore, it is named like that.
    if (ist_src .eq. ist_i) then

      if (ist .eq. ist_i) tlc%diff_c_bck_0v_elem_ssa = pix%phase_ssg(iscat_a1,:,igeo) ! From I to I: a1
      if (ist .eq. ist_q) tlc%diff_c_bck_0v_elem_ssa = pix%phase_ssg(iscat_b1,:,igeo)*geo%rotmat_c_v_allgeo(igeo) ! From I to Q: b1*Cv
      if (ist .eq. ist_u) tlc%diff_c_bck_0v_elem_ssa = pix%phase_ssg(iscat_b1,:,igeo)*geo%rotmat_s_v_allgeo(igeo) ! From I to U: b1*Sv
      if (ist .eq. ist_v) tlc%diff_c_bck_0v_elem_ssa = 0.0 ! From I to V: Zero

    endif

    if (ist_src .eq. ist_q) then

      if (ist .eq. ist_i) tlc%diff_c_bck_0v_elem_ssa = pix%phase_ssg(iscat_b1,:,igeo)*geo%rotmat_c_0_allgeo(igeo) ! From Q to I: b1*C0
      if (ist .eq. ist_q) tlc%diff_c_bck_0v_elem_ssa = pix%phase_ssg(iscat_a2,:,igeo)*geo%rotmat_c_0_allgeo(igeo)*geo%rotmat_c_v_allgeo(igeo) - pix%phase_ssg(iscat_a3,:,igeo)*geo%rotmat_s_0_allgeo(igeo)*geo%rotmat_s_v_allgeo(igeo) ! From Q to Q: a2*C0*Cv - a3*S0*Sv
      if (ist .eq. ist_u) tlc%diff_c_bck_0v_elem_ssa = pix%phase_ssg(iscat_a2,:,igeo)*geo%rotmat_c_0_allgeo(igeo)*geo%rotmat_s_v_allgeo(igeo) + pix%phase_ssg(iscat_a3,:,igeo)*geo%rotmat_s_0_allgeo(igeo)*geo%rotmat_c_v_allgeo(igeo) ! From Q to U: a2*C0*Sv + a3*S0*Cv
      if (ist .eq. ist_v) tlc%diff_c_bck_0v_elem_ssa = -pix%phase_ssg(iscat_b2,:,igeo)*geo%rotmat_s_0_allgeo(igeo) ! From Q to V: -b2*S0

    endif

    if (ist_src .eq. ist_u) then

      if (ist .eq. ist_i) tlc%diff_c_bck_0v_elem_ssa = -pix%phase_ssg(iscat_b1,:,igeo)*geo%rotmat_s_0_allgeo(igeo) ! From U to I: -b1*S0
      if (ist .eq. ist_q) tlc%diff_c_bck_0v_elem_ssa = -pix%phase_ssg(iscat_a2,:,igeo)*geo%rotmat_s_0_allgeo(igeo)*geo%rotmat_c_v_allgeo(igeo) - pix%phase_ssg(iscat_a3,:,igeo)*geo%rotmat_c_0_allgeo(igeo)*geo%rotmat_s_v_allgeo(igeo) ! From U to Q: -a2*S0*Cv - a3*C0*Sv
      if (ist .eq. ist_u) tlc%diff_c_bck_0v_elem_ssa = -pix%phase_ssg(iscat_a2,:,igeo)*geo%rotmat_s_0_allgeo(igeo)*geo%rotmat_s_v_allgeo(igeo) + pix%phase_ssg(iscat_a3,:,igeo)*geo%rotmat_c_0_allgeo(igeo)*geo%rotmat_c_v_allgeo(igeo) ! From U to U: -a2*S0*Sv + a3*C0*Cv
      if (ist .eq. ist_v) tlc%diff_c_bck_0v_elem_ssa = -pix%phase_ssg(iscat_b2,:,igeo)*geo%rotmat_c_0_allgeo(igeo) ! From U to V: -b2*C0

    endif

    if (ist_src .eq. ist_v) then

      if (ist .eq. ist_i) tlc%diff_c_bck_0v_elem_ssa = 0.0 ! From V to I: Zero
      if (ist .eq. ist_q) tlc%diff_c_bck_0v_elem_ssa = -pix%phase_ssg(iscat_b2,:,igeo)*geo%rotmat_s_v_allgeo(igeo) ! From V to Q: -b2*Sv
      if (ist .eq. ist_u) tlc%diff_c_bck_0v_elem_ssa = pix%phase_ssg(iscat_b2,:,igeo)*geo%rotmat_c_v_allgeo(igeo) ! From V to U: b2*Cv
      if (ist .eq. ist_v) tlc%diff_c_bck_0v_elem_ssa = pix%phase_ssg(iscat_a4,:,igeo) ! From V to V: a4

    endif

    ! As C is linear in the single-scattering albedo, it can be constructed by
    ! multiplying the single-scattering albedo with the derivative of C with
    ! respect to the single-scattering albedo.
    tlc%c_bck_0v_elem = tlc%diff_c_bck_0v_elem_ssa * pix%ssa

    ! BDRF, use same joke as for atmospheric C. Though the derivative itself
    ! is less interesting, only some direction cosines.
    tlc%c_bck_srf_0v_elem = pix%bdrf_ssg(ist_src,ist,igeo) * geo%diff_c_bck_srf_0v_elem_bdrf_allgeo(igeo)

    ! Emissivity, use same joke again.
    tlc%c_emi_v_elem = pix%emi_ssg(ist,igeo) * geo%diff_c_emi_v_elem_emi_allgeo(igeo)

  end subroutine calculate_relevant_single_scattering_c_parameters ! }}}

  !> Calculates all parameters that are different for different Fourier numbers, but do not depend on the viewing geometry.
  !!
  !! These parameters are calculated at the beginning of the Fourier loop. When the next
  !! Fourier number is calculated, the TLC-parameters of the old Fourier number are no longer
  !! needed and overwritten. Fourier-dependent TLC parameters are always C, because that is
  !! where the phase function or related (e.g. BDRF) parameters are.
  subroutine calculate_fourier_dependent_tlc_parameters(nst,nstrhalf,nleg,nlay,execution,m,ist_src,exist_bdrf,exist_emi,grd,geo,pix,tlc) ! {{{

    use lintran_constants_module, only: ist_i, ist_q, ist_u, ist_v, iscat_a1, iscat_b1, iscat_a2, iscat_a3, iscat_b2, iscat_a4, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_grid_class, lintran_geometry_class, lintran_pixel_class, lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    logical, dimension(3), intent(in) :: execution !< Flags for executing single, double and multi-scattering, in that order.
    integer, intent(in) :: m !< Fourier number.
    integer, intent(in) :: ist_src !< Source Stokes parameter.
    logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
    logical, intent(in) :: exist_emi !< Flag for if fluorescent emission exist on this Fourier number.
    type(lintran_grid_class), intent(in) :: grd !< Internal grid structure.
    type(lintran_geometry_class), intent(in) :: geo !< Internal geometry structure.
    type(lintran_pixel_class), intent(in) :: pix !< Internal pixel structure.
    type(lintran_tlc_class), intent(inout) :: tlc !< Internal radiative-transfer building-block structure.

    ! Iterator.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.
    integer :: istr_src ! Over source streams, when two iterators are needed.

    real, dimension(m:nleg) :: sym_even ! Symmetry relationship, plus for even ileg-m and minus for odd ileg-m.
    real, dimension(m:nleg) :: sym_odd ! Symmetry relationship, plus for odd ileg-m and minus for even ileg-m.

    ! Instrument-related TLC-parameters are skipped, because they are set in the loop
    ! over geometries inside the loop over Fourier numbers. That is not done for
    ! Fourier-independent TLC-parameters, because the Fourier loop is outside the
    ! geometry loop. Therefore, Fourier-independent TLC-parameters have an allgeo-suffix
    ! and pointers are set in the loop over geometries, while Fourier-dependent
    ! TLC-parameters can be calculated when it is that geometry's turn and after that
    ! it can be thrown away, because the next time, the Fourier number will be different.

    ! Set symmetry relationships for even (0 and p) and odd (m) generalized spherical functions
    ! Even here means that the function is even for even ileg-m and odd for odd ileg-m.
    ! Odd here means that the function is odd for even ileg-m and even for odd ileg-m.
    ! The array index is ileg.
    sym_even(m:nleg:2) = 1.0
    sym_even(m+1:nleg:2) = -1.0
    sym_odd(m:nleg:2) = -1.0
    sym_odd(m+1:nleg:2) = 1.0

    ! Execution 2 or 3 is always true. Otherwise, this routine is not called.

    ! Calculate phase functions, the derivatives of C with respect to the
    ! single-scattering albedo.

    ! We start with the normal scattering direction from the sun to the internal stream.
    ! The other way around is calculated for the adjoint, because the adjoint uses a
    ! pseudo-forward model.
    ! For forward scattering, all directions are positive. For backward scattering,
    ! the source direction is reversed, so she will get her symmetry relationship.
    ! Furthermore, the source Stokes parameters U and V will be taken negative, because
    ! Zb = Z_up_dn Delta_34 (=Delta_34 Z_dn_up), which would come down to applying both
    ! the symmetry relationship and the U-V reversal to the destination direction, but
    ! we choose the source.

    do ilay = 1,nlay

      if (ist_src .eq. ist_i) then

        tlc%diff_c_fwd_0_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! I -> I.
        tlc%diff_c_bck_0_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! I -> I, polynomial 0 = even, source I = positive, result = even.

        if (pol_nst(nst)) then
          tlc%diff_c_fwd_0_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! I -> Q.
          tlc%diff_c_bck_0_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! I -> Q, polynomial 0 = even, source I = positive, result = even.

          tlc%diff_c_fwd_0_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! I -> U.
          tlc%diff_c_bck_0_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! I -> U, polynomial 0 = even, source I = positive, result = even.
        endif

        if (fullpol_nst(nst)) then
          tlc%diff_c_fwd_0_nrm_ssa(ist_v,:,ilay) = 0.0 ! I -> V.
          tlc%diff_c_bck_0_nrm_ssa(ist_v,:,ilay) = 0.0 ! I -> V, zero, odd or even does not matter.
        endif

      endif

      if (ist_src .eq. ist_q) then

        tlc%diff_c_fwd_0_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) ! Q -> I.
        tlc%diff_c_bck_0_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) ! Q -> I, polynomial p = even, source Q = positive, result = even.

        ! Because of the source Stokes parameter, we can assume that there is polarization.
        tlc%diff_c_fwd_0_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) ! Q -> Q.
        tlc%diff_c_bck_0_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) ! Q -> Q, polynomials p-m = even-odd, source Q = positive, result = even-odd.

        tlc%diff_c_fwd_0_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) ! Q -> U.
        tlc%diff_c_bck_0_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) ! Q -> U, polynomials p-m = even-odd, source Q = positive, result = even-odd.

        if (fullpol_nst(nst)) then
          tlc%diff_c_fwd_0_nrm_ssa(ist_v,:,ilay) = -matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) ! Q -> V.
          tlc%diff_c_bck_0_nrm_ssa(ist_v,:,ilay) = -matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) ! Q -> V, polynomial m = odd, source Q = positive, result = odd.
        endif

      endif

      if (ist_src .eq. ist_u) then

        tlc%diff_c_fwd_0_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) ! U -> I.
        tlc%diff_c_bck_0_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) ! U -> I, polynomial m = odd, source U = negative, result = even.

        ! Because of the source Stokes parameter, we can assume that there is polarization.
        tlc%diff_c_fwd_0_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) ! U -> Q.
        tlc%diff_c_bck_0_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) ! U -> Q, polynomials m-p = odd-even, source U = negative, result = even-odd.

        tlc%diff_c_fwd_0_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) ! U -> U.
        tlc%diff_c_bck_0_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_m_0(m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) ! U -> U, polynomials m-p = odd-even, source U = negative, result = even-odd.

        if (fullpol_nst(nst)) then
          tlc%diff_c_fwd_0_nrm_ssa(ist_v,:,ilay) = -matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) ! U -> V.
          tlc%diff_c_bck_0_nrm_ssa(ist_v,:,ilay) = -matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_p_0(m:pix%nleg_lay(ilay),m)) ! U -> V, polynomial p = even, source U = negative, result = odd.
        endif

      endif

      if (ist_src .eq. ist_v) then

        ! Now, we can assume that all possible polarization is involved.
        tlc%diff_c_fwd_0_nrm_ssa(ist_i,:,ilay) = 0.0 ! V -> I.
        tlc%diff_c_bck_0_nrm_ssa(ist_i,:,ilay) = 0.0 ! V -> I, zero, odd or even does not matter.

        tlc%diff_c_fwd_0_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! V -> Q.
        tlc%diff_c_bck_0_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! V -> Q, polynomial 0 = even, source V = negative, result = odd.

        tlc%diff_c_fwd_0_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! V -> U.
        tlc%diff_c_bck_0_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! V -> U, polynomial 0 = even, source V = negative, result = odd.

        tlc%diff_c_fwd_0_nrm_ssa(ist_v,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a4,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! V -> V.
        tlc%diff_c_bck_0_nrm_ssa(ist_v,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a4,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_0_0(m:pix%nleg_lay(ilay),m)) ! V -> V, polynomial 0 = even, source V = negative, result = odd.

      endif

      ! Multiply all C-parameters with the single-scattering albedo as this has not
      ! yet been done.
      tlc%c_fwd_0_nrm(:,:,ilay) = tlc%diff_c_fwd_0_nrm_ssa(:,:,ilay) * pix%ssa(ilay)
      tlc%c_bck_0_nrm(:,:,ilay) = tlc%diff_c_bck_0_nrm_ssa(:,:,ilay) * pix%ssa(ilay)

    enddo

    ! If V is not involved, we do not have to discuss the adjoint,
    ! becuase it is the same as the normal on pointer level.
    if (fullpol_nst(nst)) then

      ! For the adjoint, we will transpose this matrix and flip the signs where
      ! V-polarization is involved. The transposition is not necessary, because the
      ! sun-related indices are not there. The sign-flip needs further discussion.
      tlc%diff_c_fwd_0_adj_ssa = tlc%diff_c_fwd_0_nrm_ssa
      tlc%diff_c_bck_0_adj_ssa = tlc%diff_c_bck_0_nrm_ssa

      ! Do the signswitch on polarization state V. Actually, this should be a full Delta-flip
      ! on V. But the one at the the sun has to be flipped again, because in the adjoint,
      ! the V-polarization is negative all the time. This minus sign is applied by making V-light
      ! negative at the source and at the destination. This is applied by switching the sign on each
      ! phase function where V-light is involved in solar or viewing direction.
      ! As a final result of the two operations together, the sign is switched for the adjoint
      ! where the internal stream is in V.
      tlc%diff_c_fwd_0_adj_ssa(ist_v,:,:) = -tlc%diff_c_fwd_0_adj_ssa(ist_v,:,:)
      tlc%diff_c_bck_0_adj_ssa(ist_v,:,:) = -tlc%diff_c_bck_0_adj_ssa(ist_v,:,:)

      ! Perform the multiplication with the single-scattering albedo for the adjoint.
      do ilay = 1,nlay

        tlc%c_fwd_0_adj(:,:,ilay) = tlc%diff_c_fwd_0_adj_ssa(:,:,ilay) * pix%ssa(ilay)
        tlc%c_bck_0_adj(:,:,ilay) = tlc%diff_c_bck_0_adj_ssa(:,:,ilay) * pix%ssa(ilay)

      enddo

    endif

    ! Surface reflection and emissivity.
    ! Albedo.
    if (exist_bdrf) then
      do istr = 1,nstrhalf ! Relevant stream.

        ! Array operations is over Stokes parameters for destination.
        tlc%c_bck_srf_0_nrm(:,istr) = pix%bdrf_0(ist_src,:,istr,m) * geo%diff_c_bck_srf_0_bdrf(istr)
      enddo

      ! Apply Delta_34 on upward (destination) side of the surface reflection C-parameter. What
      ! of this Delta-flip exists depends on the number of Stokes parameters.
      if (pol_nst(nst)) then
        tlc%c_bck_srf_0_nrm(ist_u:nst,:) = -tlc%c_bck_srf_0_nrm(ist_u:nst,:)
      endif

      ! Do the transposition for the adjoint is needed.
      if (fullpol_nst(nst)) then
        ! Transpose BDRF matrices for reverse directions (without actual transpotition, see
        ! comments at diff_c_bck/fwd_i0).
        tlc%c_bck_srf_0_adj = tlc%c_bck_srf_0_nrm
        ! Use symmetry relationship by flipping with Delta_4, see comments at atmospheric
        ! C-parameters for the adjoint why ultimately, the Delta-flip is applied at the internal
        ! side (i) of the C-parameters.
        tlc%c_bck_srf_0_adj(ist_v,:) = -tlc%c_bck_srf_0_adj(ist_v,:)
      endif

    endif

    ! Emission.
    if (exist_emi) then
      do istr = 1,nstrhalf ! Relevant stream.

        ! Array operations are over Stokes parameters.
        tlc%c_emi_nrm(:,istr) = pix%emi(:,istr,m) * grd%diff_c_emi_nrm_emi(:,istr)
        ! If V is not involved, the adjoint is already done pointer-wise.
        if (fullpol_nst(nst)) then
          tlc%c_emi_adj(:,istr) = pix%emi(:,istr,m) * grd%diff_c_emi_adj_emi(:,istr)
        endif

      enddo

      ! Switch signs of U and V for emissivity, because it is upward radiation.
      if (pol_nst(nst)) then
        tlc%c_emi_nrm(ist_u:nst,:) = -tlc%c_emi_nrm(ist_u:nst,:)
        ! In fact, there is no need to nest this if-clause because the condition includes
        ! the condision of the outer if-clause. But the outer if-clause means if we need to
        ! do the Delta-flip and the inner if-clause is whether we have to handle the adjoint.
        ! Those are two different things and both need to be true for the next operation.
        if (fullpol_nst(nst)) then
          tlc%c_emi_adj(ist_u:nst,:) = -tlc%c_emi_adj(ist_u:nst,:)
        endif
      endif
      
    endif

    ! Scattering from one internal stream to another, only needed for multi-scattering.
    if (execution(3)) then

      ! Calculate C-parameters while saving the derivatives as intermediate results.
      do ilay = 1,nlay

        do istr = 1,nstrhalf ! Destination stream

          ! This is for internal stream to another internal stream. The array operation is over
          ! source streams. For backward scattering, the symmetry relationships are applied to
          ! the destination direction, because that is an array and not a matrix. So it will
          ! be the same as for the instrument geometry from the instrument-dependent routine.

          ! The elements are sorted to the extent of polarization that must be involved.
          tlc%diff_c_fwd_ssa(ist_i,ist_i,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a1,m:pix%nleg_lay(ilay),ilay) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! I -> I.
          tlc%diff_c_bck_ssa(ist_i,ist_i,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! I -> I, polynomial 0 = even, destination I = positive, result = even.

          ! Any polarization.
          if (pol_nst(nst)) then

            tlc%diff_c_fwd_ssa(ist_q,ist_i,:,istr,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! Q -> I.
            tlc%diff_c_bck_ssa(ist_q,ist_i,:,istr,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! Q -> I, polynomial 0 = even, destination I = positive, result = even.

            tlc%diff_c_fwd_ssa(ist_u,ist_i,:,istr,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! U -> I.
            tlc%diff_c_bck_ssa(ist_u,ist_i,:,istr,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! U -> I, polynomial 0 = even, destination I = positive, result = even.

            tlc%diff_c_fwd_ssa(ist_i,ist_q,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) ! I -> Q.
            tlc%diff_c_bck_ssa(ist_i,ist_q,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) ! I -> Q, polynomial p = even, destination Q = positive, result = even.

            tlc%diff_c_fwd_ssa(ist_q,ist_q,:,istr,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) ! Q -> Q.
            tlc%diff_c_bck_ssa(ist_q,ist_q,:,istr,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) ! Q -> Q, polynomials p-m = even-odd, destination Q = positive, result = even-odd.

            tlc%diff_c_fwd_ssa(ist_u,ist_q,:,istr,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) ! U -> Q.
            tlc%diff_c_bck_ssa(ist_u,ist_q,:,istr,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) ! U -> Q, polynomials p-m = even-odd, destination Q = positive, result = even-odd.

            tlc%diff_c_fwd_ssa(ist_i,ist_u,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) ! I -> U.
            tlc%diff_c_bck_ssa(ist_i,ist_u,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) ! I -> U, polynomial m = odd, destination U = negative, result = even.

            tlc%diff_c_fwd_ssa(ist_q,ist_u,:,istr,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) ! Q -> U.
            tlc%diff_c_bck_ssa(ist_q,ist_u,:,istr,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) ! Q -> U, polynomials m-p = odd-even, destination U = negative, result = even-odd.

            tlc%diff_c_fwd_ssa(ist_u,ist_u,:,istr,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) ! U -> U.
            tlc%diff_c_bck_ssa(ist_u,ist_u,:,istr,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) ! U -> U, polynomials m-p = odd-even, destination U = negative, result = even-odd.

          endif

          ! Full polarization with V.

          if (fullpol_nst(nst)) then

            tlc%diff_c_fwd_ssa(ist_v,ist_i,:,istr,ilay) = 0.0 ! V -> I.
            tlc%diff_c_bck_ssa(ist_v,ist_i,:,istr,ilay) = 0.0 ! V -> I, zero, odd or even does not matter.

            tlc%diff_c_fwd_ssa(ist_v,ist_q,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) ! V -> Q.
            tlc%diff_c_bck_ssa(ist_v,ist_q,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * grd%gsf_m(istr,m:pix%nleg_lay(ilay),m)) ! V -> Q, polynomial m = odd, destination Q = positive, result = odd.

            tlc%diff_c_fwd_ssa(ist_v,ist_u,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) ! V -> U.
            tlc%diff_c_bck_ssa(ist_v,ist_u,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * grd%gsf_p(istr,m:pix%nleg_lay(ilay),m)) ! V -> U, polynomial p = even, destination U = negative, result = odd.

            tlc%diff_c_fwd_ssa(ist_i,ist_v,:,istr,ilay) = 0.0 ! I -> V.
            tlc%diff_c_bck_ssa(ist_i,ist_v,:,istr,ilay) = 0.0 ! I -> V, zero, odd or even does not matter.

            tlc%diff_c_fwd_ssa(ist_q,ist_v,:,istr,ilay) = -matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! Q -> V.
            tlc%diff_c_bck_ssa(ist_q,ist_v,:,istr,ilay) = -matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! Q -> V, polynomial 0 = even, destination V = negative, result = odd.

            tlc%diff_c_fwd_ssa(ist_u,ist_v,:,istr,ilay) = -matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! U -> V.
            tlc%diff_c_bck_ssa(ist_u,ist_v,:,istr,ilay) = -matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! U -> V, polynomial 0 = even, destination V = negative, result = odd.

            tlc%diff_c_fwd_ssa(ist_v,ist_v,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a4,m:pix%nleg_lay(ilay),ilay) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! V -> V.
            tlc%diff_c_bck_ssa(ist_v,ist_v,:,istr,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a4,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * grd%gsf_0(istr,m:pix%nleg_lay(ilay),m)) ! V -> V, polynomial 0 = even, destination V = negative, result = odd.

          endif

        enddo

        ! Multiply all C-parameters with the single-scattering albedo as this has not
        ! yet been done.
        tlc%c_fwd(:,:,:,:,ilay) = tlc%diff_c_fwd_ssa(:,:,:,:,ilay) * pix%ssa(ilay)
        tlc%c_bck(:,:,:,:,ilay) = tlc%diff_c_bck_ssa(:,:,:,:,ilay) * pix%ssa(ilay)

      enddo

      ! Surface reflection. A pity for the double do-loop with two similar iterators.
      if (exist_bdrf) then
        do istr = 1,nstrhalf
          do istr_src = 1,nstrhalf
            tlc%c_bck_srf(:,:,istr_src,istr) = pix%bdrf(:,:,istr_src,istr,m) * grd%diff_c_bck_srf_bdrf(istr_src,istr)
          enddo
        enddo

        if (pol_nst(nst)) then
          ! The destination stream has to get a Delta_34, because that is the upward stream
          tlc%c_bck_srf(:,ist_u:nst,:,:) = -tlc%c_bck_srf(:,ist_u:nst,:,:)
        endif
      endif

    endif

  end subroutine calculate_fourier_dependent_tlc_parameters ! }}}

  !> Selects a viewing geometry for the Fourier-independent TLC-parameters.
  !!
  !! Because the loop over Fourier numbers is outside the loop over viewing geometries,
  !! the parameters that do depend on the geometry, but not on the Fourier numbers are
  !! used several times. It is not desirable to recalculate the same values every time
  !! the same geometry is iterated during the differnet iterations of the Fourier loop,
  !! they are calculated once. Each time the geometry is iteratred, a few pointers are
  !! set, a computationally cheap operation.
  subroutine choose_viewing_geometry_for_fourier_independent_tlc_parameters(nlay,flag_derivatives,igeo,geo,tlc) ! {{{

    use lintran_types_module, only: lintran_geometry_class, lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    logical, intent(in) :: flag_derivatives !< Flag for calculating derivatives.
    integer, intent(in) :: igeo !< Viewing geometry index.
    type(lintran_geometry_class), intent(in), target :: geo !< Internal geometry structure.
    type(lintran_tlc_class), intent(inout), target :: tlc !< Internal radiative-transfer building-block structure.

    ! Set the pointers for all viewing-geometry TLC-parameters.

    ! Execution and split-double flags are not evaluated, because setting pointers
    ! takes no computer time. Dragging and evaluating flags may cost even more.
    ! In such a case, pointers may point to uninitialized data and will not be used.
    ! Flag-derivatives is evaluated because some fields in the TLC-module may not
    ! be allocated without derivatives. In the geometry class, this is not (yet) done
    ! because only very little memory is allocated for just the derivatives.

    ! From the TLC-module herself.
    tlc%t_tot_v(0:nlay) => tlc%t_tot_v_allgeo(:,igeo)
    tlc%t_tot_0v(0:nlay) => tlc%t_tot_0v_allgeo(:,igeo)
    tlc%l_0v => tlc%l_0v_allgeo(:,igeo)
    tlc%l_pv => tlc%l_pv_allgeo(:,:,igeo)
    tlc%l_mv => tlc%l_mv_allgeo(:,:,igeo)
    tlc%l_v => tlc%l_v_allgeo(:,igeo)
    tlc%effmu_mv => tlc%effmu_mv_allgeo(:,:,igeo)
    tlc%effmu_pv => tlc%effmu_pv_allgeo(:,:,igeo)
    tlc%t_sp_v => tlc%t_sp_v_allgeo(:,igeo)
    tlc%l_sp_pv => tlc%l_sp_pv_allgeo(:,:,igeo)
    tlc%l_sp_mv => tlc%l_sp_mv_allgeo(:,:,igeo)
    tlc%l_1_sp_pv => tlc%l_1_sp_pv_allgeo(:,:,igeo)
    tlc%l_int_sm_v => tlc%l_int_sm_v_allgeo(:,igeo)
    tlc%l_int_ac_v => tlc%l_int_ac_v_allgeo(:,igeo)

    ! The rest is only for derivatives.
    if (flag_derivatives) then

      ! From geometry.
      tlc%reldiff_t_tot_v_tau(1:nlay,0:nlay) => geo%reldiff_t_tot_v_tau_allgeo(:,:,igeo)
      tlc%reldiff_t_tot_0v_tau(1:nlay,0:nlay) => geo%reldiff_t_tot_0v_tau_allgeo(:,:,igeo)
      tlc%diff_c_bck_srf_v_bdrf => geo%diff_c_bck_srf_v_bdrf_allgeo(:,igeo)
      tlc%diff_c_bck_srf_0v_elem_bdrf => geo%diff_c_bck_srf_0v_elem_bdrf_allgeo(igeo)
      tlc%diff_c_emi_v_elem_emi => geo%diff_c_emi_v_elem_emi_allgeo(igeo)

      ! From TLC herself.
      tlc%ideriv_sensitive_v(0:nlay) => tlc%ideriv_sensitive_v_allgeo(:,igeo)
      tlc%ideriv_sensitive_0v(0:nlay) => tlc%ideriv_sensitive_0v_allgeo(:,igeo)
      tlc%diff_effmu_pv_tau => tlc%diff_effmu_pv_tau_allgeo(:,:,:,igeo)
      tlc%diff_effmu_mv_tau => tlc%diff_effmu_mv_tau_allgeo(:,:,:,igeo)
      tlc%reldiff_t_sp_v_tau => tlc%reldiff_t_sp_v_tau_allgeo(:,:,igeo)
      tlc%diff_l_0v_tau => tlc%diff_l_0v_tau_allgeo(:,:,igeo)
      tlc%diff_l_pv_tau => tlc%diff_l_pv_tau_allgeo(:,:,:,igeo)
      tlc%diff_l_mv_tau => tlc%diff_l_mv_tau_allgeo(:,:,:,igeo)
      tlc%diff_l_v_tau => tlc%diff_l_v_tau_allgeo(:,:,igeo)
      tlc%diff_l_sp_pv_tau => tlc%diff_l_sp_pv_tau_allgeo(:,:,:,igeo)
      tlc%diff_l_sp_mv_tau => tlc%diff_l_sp_mv_tau_allgeo(:,:,:,igeo)
      tlc%diff_l_int_sm_v_tau => tlc%diff_l_int_sm_v_tau_allgeo(:,:,igeo)
      tlc%diff_l_int_ac_v_tau => tlc%diff_l_int_ac_v_tau_allgeo(:,:,igeo)

    endif

  end subroutine choose_viewing_geometry_for_fourier_independent_tlc_parameters ! }}}

  !> Calculate parameters that depend both on Fourier number, viewing Stokes parameter and viewing geometry.
  !!
  !! These parameters are different for each combination of Fourier number and viewing
  !! Stokes parameter and viewing geometry, so after they are needed, they can be thrown away.
  !! In the vector code, the conversion from derivatives with respect to C to derivatives with
  !! respect to single-scattering albedo and phase coefficient is done inside the loop, so
  !! also the derivatives for the chain rule are not needed after the relevant loop iteration.
  subroutine calculate_fourier_dependent_tlc_parameters_for_viewing_geometry(nst,nstrhalf,nleg,nlay,igeo,m,ist,exist_bdrf,grd,geo,pix,tlc) ! {{{

    use lintran_constants_module, only: ist_i, ist_q, ist_u, ist_v, iscat_a1, iscat_b1, iscat_a2, iscat_a3, iscat_b2, iscat_a4, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_grid_class, lintran_geometry_class, lintran_pixel_class, lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: igeo !< Viewing geometry index.
    integer, intent(in) :: m !< Fourier number.
    integer, intent(in) :: ist !< Viewing Stokes parameter.
    logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
    type(lintran_grid_class), intent(in) :: grd !< Internal grid structure.
    type(lintran_geometry_class), intent(in) :: geo !< Internal geometry structure.
    type(lintran_pixel_class), intent(in) :: pix !< Internal pixel structure.
    type(lintran_tlc_class), intent(inout) :: tlc !< Internal radiative-transfer building-block structure.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams

    ! Symmetry relationships.
    real, dimension(m:nleg) :: sym_even ! Symmetry relationship, plus for even ileg-m and minus for odd ileg-m
    real, dimension(m:nleg) :: sym_odd ! Symmetry relationship, plus for odd ileg-m and minus for even ileg-m

    ! Set symmetry relationships for even (0 and p) and odd (m) generalized spherical functions
    ! Even here means that the function is even for even ileg-m and odd for odd ileg-m.
    ! Odd here means that the function is odd for even ileg-m and even for odd ileg-m.
    ! The array index is ileg.
    sym_even(m:nleg:2) = 1.0
    sym_even(m+1:nleg:2) = -1.0
    sym_odd(m:nleg:2) = -1.0
    sym_odd(m+1:nleg:2) = 1.0

    ! Calculate C-parameters while saving the derivative, the intermediate
    ! result. The intermediate result is saved for all geometries, so that
    ! they can be used at the end when converting derivatives with respect to
    ! C to derivatives with respect to single-scattering albedo and phase
    ! coefficients.
    do ilay = 1,nlay

      ! This is the normal, non-adjoint phase function. If V is involved, the adjoint
      ! will be derived. Otherwise, the normal will automatically be used for the adjoint
      ! as well.

      if (ist .eq. ist_i) then

        tlc%diff_c_fwd_v_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! I -> I.
        tlc%diff_c_bck_v_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! I -> I, polynomial 0 = even, destination I = positive, result = even.

        ! Only if polarization is involved.
        if (pol_nst(nst)) then
          tlc%diff_c_fwd_v_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! Q -> I.
          tlc%diff_c_bck_v_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! Q -> I, polynomial 0 = even, destination I = positive, result = even.

          tlc%diff_c_fwd_v_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! U -> I.
          tlc%diff_c_bck_v_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! U -> I, polynomial 0 = even, destination I = positive, result = even.
        endif

        ! Only for full polarization.
        if (fullpol_nst(nst)) then
          tlc%diff_c_fwd_v_nrm_ssa(ist_v,:,ilay) = 0.0 ! V -> I.
          tlc%diff_c_bck_v_nrm_ssa(ist_v,:,ilay) = 0.0 ! V -> I, zero, odd or even does not matter.
        endif

      endif

      if (ist .eq. ist_q) then

        tlc%diff_c_fwd_v_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! I -> Q.
        tlc%diff_c_bck_v_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! I -> Q, polynomial p = even, destination Q = positive, result = even.

        ! Now, we can assume that we have polarization at least up to Q and U.
        tlc%diff_c_fwd_v_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! Q -> Q.
        tlc%diff_c_bck_v_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! Q -> Q, polynomials p-m = even-odd, destination Q = positive, result = even-odd.

        tlc%diff_c_fwd_v_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! U -> Q.
        tlc%diff_c_bck_v_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! U -> Q, polynomials p-m = even-odd, destination Q = positive, result = even-odd.

        ! Only for full polarization.
        if (fullpol_nst(nst)) then
          tlc%diff_c_fwd_v_nrm_ssa(ist_v,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! V -> Q.
          tlc%diff_c_bck_v_nrm_ssa(ist_v,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! V -> Q, polynomial m = odd, destination Q = positive, result = odd.
        endif

      endif

      if (ist .eq. ist_u) then

        tlc%diff_c_fwd_v_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! I -> U.
        tlc%diff_c_bck_v_nrm_ssa(ist_i,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b1,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! I -> U, polynomial m = odd, destination U = negative, result = even.

        ! Now, we can assume that we have polarization at least up to Q and U.
        tlc%diff_c_fwd_v_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! Q -> U.
        tlc%diff_c_bck_v_nrm_ssa(ist_q,:,ilay) = matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) + matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! Q -> U, polynomials m-p = odd-even, destination U = negative, result = even-odd.

        tlc%diff_c_fwd_v_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! U -> U.
        tlc%diff_c_bck_v_nrm_ssa(ist_u,:,ilay) = matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a2,m:pix%nleg_lay(ilay),ilay) * sym_even(m:pix%nleg_lay(ilay)) * geo%gsf_m_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) + matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a3,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! U -> U, polynomials m-p = odd-even, destination U = negative, result = even-odd.

        ! Only for full polarization.
        if (fullpol_nst(nst)) then
          tlc%diff_c_fwd_v_nrm_ssa(ist_v,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! V -> U.
          tlc%diff_c_bck_v_nrm_ssa(ist_v,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_p_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! V -> U, polynomial p = even, destination U = negative, result = odd.
        endif

      endif

      if (ist .eq. ist_v) then

        ! Now, we can assume full polarization.
        tlc%diff_c_fwd_v_nrm_ssa(ist_i,:,ilay) = 0.0 ! I -> V.
        tlc%diff_c_bck_v_nrm_ssa(ist_i,:,ilay) = 0.0 ! I -> V, zero, odd or even does not matter.

        tlc%diff_c_fwd_v_nrm_ssa(ist_q,:,ilay) = -matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! Q -> V.
        tlc%diff_c_bck_v_nrm_ssa(ist_q,:,ilay) = -matmul(grd%gsf_m(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! Q -> V, polynomial 0 = even, destination V = negative, result = odd.

        tlc%diff_c_fwd_v_nrm_ssa(ist_u,:,ilay) = -matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! U -> V.
        tlc%diff_c_bck_v_nrm_ssa(ist_u,:,ilay) = -matmul(grd%gsf_p(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_b2,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! U -> V, polynomial 0 = even, destination V = negative, result = odd.

        tlc%diff_c_fwd_v_nrm_ssa(ist_v,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a4,m:pix%nleg_lay(ilay),ilay) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! V -> V.
        tlc%diff_c_bck_v_nrm_ssa(ist_v,:,ilay) = matmul(grd%gsf_0(:,m:pix%nleg_lay(ilay),m) , pix%coefs(iscat_a4,m:pix%nleg_lay(ilay),ilay) * sym_odd(m:pix%nleg_lay(ilay)) * geo%gsf_0_v_allgeo(m:pix%nleg_lay(ilay),m,igeo)) ! V -> V, polynomial 0 = even, destination V = negative, result = odd.

      endif

    enddo

    ! Multiply all C-parameters with the single-scattering albedo.
    do ilay = 1,nlay

      tlc%c_fwd_v_nrm(:,:,ilay) = tlc%diff_c_fwd_v_nrm_ssa(:,:,ilay) * pix%ssa(ilay)
      tlc%c_bck_v_nrm(:,:,ilay) = tlc%diff_c_bck_v_nrm_ssa(:,:,ilay) * pix%ssa(ilay)

    enddo

    ! Handle the adjoint separately. Without V, the pointers do their work.
    if (fullpol_nst(nst)) then

      ! For the adjoint, we will transpose this matrix and flip the signs where
      ! V-polarization is involved.
      tlc%diff_c_fwd_v_adj_ssa = tlc%diff_c_fwd_v_nrm_ssa
      tlc%diff_c_bck_v_adj_ssa = tlc%diff_c_bck_v_nrm_ssa

      ! Do the signswitch on polarization state V. Actually, this should be a full Delta-flip
      ! on V. But the one at the instrument or the sun has to be flipped again, because
      ! in the adjoint, the V-polarization is negative all the time. This minus sign is
      ! applied by making V-light negative at the source and at the destination. This
      ! is applied by switching the sign on each phase function where V-light is involved
      ! in solar or viewing direction.
      ! As a final result of the two operations together, the sign is switched for the adjoint
      ! where the internal stream is in V.
      tlc%diff_c_fwd_v_adj_ssa(ist_v,:,:) = -tlc%diff_c_fwd_v_adj_ssa(ist_v,:,:)
      tlc%diff_c_bck_v_adj_ssa(ist_v,:,:) = -tlc%diff_c_bck_v_adj_ssa(ist_v,:,:)

      ! And do the multiplication with the single-scattering albedo for the adjoint.
      do ilay = 1,nlay

        tlc%c_fwd_v_adj(:,:,ilay) = tlc%diff_c_fwd_v_adj_ssa(:,:,ilay) * pix%ssa(ilay)
        tlc%c_bck_v_adj(:,:,ilay) = tlc%diff_c_bck_v_adj_ssa(:,:,ilay) * pix%ssa(ilay)

      enddo

    endif

    ! Surface reflection.
    if (exist_bdrf) then
      do istr = 1,nstrhalf ! Relevant stream.

        ! Array operations are over Stokes parameters source.
        tlc%c_bck_srf_v_nrm(:,istr) = pix%bdrf_v(:,ist,istr,m,igeo) * geo%diff_c_bck_srf_v_bdrf_allgeo(istr,igeo)

      enddo

      ! Apply Delta_34 on upward (destination) side of the surface reflection C-parameter.
      ! That is the argument Stokes parameter, not the one of the array index.
      if (ist .gt. ist_q) then
        tlc%c_bck_srf_v_nrm = -tlc%c_bck_srf_v_nrm
      endif

      ! Adjoint. Again, only for V-polarization.
      if (fullpol_nst(nst)) then

        tlc%c_bck_srf_v_adj = tlc%c_bck_srf_v_nrm
        ! Use symmetry relationship by flipping with Delta_4, see comments at atmospheric
        ! C-parameters for the adjoint why ultimately, the Delta-flip is applied at the internal
        ! side of the C-parameters.
        tlc%c_bck_srf_v_adj(ist_v,:) = -tlc%c_bck_srf_v_adj(ist_v,:)

      endif
    endif

  end subroutine calculate_fourier_dependent_tlc_parameters_for_viewing_geometry ! }}}

  !> Cleans up the building blocks of the radiative transfer equation.
  !!
  !! All arrays that are allocated in tlc_init are cleaned up.
  subroutine tlc_close(nst,flag_derivatives,tlc,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation, fullpol_nst
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was initialized for derivatives.
    type(lintran_tlc_class), intent(inout) :: tlc !< Internal radiative-transfer building-block structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    if (flag_derivatives) then
      deallocate(tlc%diff_c_therm_planck,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_c_therm_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_int_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_int_ac_v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_int_ac_0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_int_sm_v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_int_sm_0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_int_ac_m_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_int_ac_p_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_int_sm_m_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_int_sm_p_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_sp_p_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_sp_mv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_sp_m0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_sp_pv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_sp_p0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_mv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_m0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_pv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_p0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_l_0v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%reldiff_t_sp_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%reldiff_t_sp_v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%reldiff_t_sp_0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_effmu_mv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_effmu_pv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_effmu_m0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_effmu_p0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%l_zip_ac,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%l_zip_sm,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%ideriv_sensitive_0v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%ideriv_sensitive_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%ideriv_sensitive_0,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(tlc%c_therm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (fullpol_nst(nst)) then
      deallocate(tlc%c_emi_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(tlc%c_emi_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (fullpol_nst(nst)) then
      deallocate(tlc%c_bck_srf_v_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%c_bck_srf_0_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(tlc%c_bck_srf_v_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%c_bck_srf_0_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%c_bck_srf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (fullpol_nst(nst)) then
      deallocate(tlc%diff_c_bck_v_adj_ssa_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_c_fwd_v_adj_ssa_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_c_bck_0_adj_ssa_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%diff_c_fwd_0_adj_ssa_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(tlc%diff_c_bck_v_nrm_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%diff_c_fwd_v_nrm_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%diff_c_bck_0_nrm_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%diff_c_fwd_0_nrm_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%diff_c_bck_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%diff_c_fwd_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (fullpol_nst(nst)) then
      deallocate(tlc%c_bck_v_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%c_fwd_v_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%c_bck_0_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(tlc%c_fwd_0_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(tlc%c_bck_v_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%c_fwd_v_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%c_bck_0_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%c_fwd_0_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%c_bck,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%c_fwd,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%diff_c_bck_0v_elem_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%c_bck_0v_elem,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_int,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_int_ac_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_int_ac_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_int_ac_m,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_int_ac_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_int_sm_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_int_sm_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_int_sm_m,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_int_sm_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_1_sp_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_sp_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_1_sp_pv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_1_sp_p0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_sp_mp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_sp_pp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_sp_mv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_sp_pv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_sp_m0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_sp_p0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_mv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_pv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_m0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_p0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%l_0v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%t_sp_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%t_sp_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%t_sp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%t_tot_0v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%t_tot_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%t_tot_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%t,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%effmu_pv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%effmu_mv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%effmu_p0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(tlc%effmu_m0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine tlc_close ! }}}

  !> Cleans up garbage that has been created when an error occurs during a routine.
  !!
  !! This is meant to clean up the rubbish that is created before the error occured.
  !! The idea is that the allocation status is the same as before calling the entire
  !! routine where the error occurred.
  subroutine fail(nst,progression,flag_derivatives,tlc,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation, fullpol_nst
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: progression !< Indication of how much garbage has been created.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was initialized for derivatives.
    type(lintran_tlc_class), intent(inout) :: tlc !< Internal radiative-transfer building-block structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    ierr = 0 ! Because deallocations happen in if-clauses.

    if (flag_derivatives) then
      if (progression .ge. 101) deallocate(tlc%diff_c_therm_planck,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 100) deallocate(tlc%diff_c_therm_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 99) deallocate(tlc%diff_l_int_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 98) deallocate(tlc%diff_l_int_ac_v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 97) deallocate(tlc%diff_l_int_ac_0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 96) deallocate(tlc%diff_l_int_sm_v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 95) deallocate(tlc%diff_l_int_sm_0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 94) deallocate(tlc%diff_l_int_ac_m_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 93) deallocate(tlc%diff_l_int_ac_p_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 92) deallocate(tlc%diff_l_int_sm_m_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 91) deallocate(tlc%diff_l_int_sm_p_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 90) deallocate(tlc%diff_l_sp_p_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 89) deallocate(tlc%diff_l_sp_mv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 88) deallocate(tlc%diff_l_sp_m0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 87) deallocate(tlc%diff_l_sp_pv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 86) deallocate(tlc%diff_l_sp_p0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 85) deallocate(tlc%diff_l_v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 84) deallocate(tlc%diff_l_mv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 83) deallocate(tlc%diff_l_m0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 82) deallocate(tlc%diff_l_pv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 81) deallocate(tlc%diff_l_p0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 80) deallocate(tlc%diff_l_0v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 79) deallocate(tlc%reldiff_t_sp_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 78) deallocate(tlc%reldiff_t_sp_v_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 77) deallocate(tlc%reldiff_t_sp_0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 76) deallocate(tlc%diff_effmu_mv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 75) deallocate(tlc%diff_effmu_pv_tau_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 74) deallocate(tlc%diff_effmu_m0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 73) deallocate(tlc%diff_effmu_p0_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 72) deallocate(tlc%l_zip_ac,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 71) deallocate(tlc%l_zip_sm,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 70) deallocate(tlc%ideriv_sensitive_0v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 69) deallocate(tlc%ideriv_sensitive_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 68) deallocate(tlc%ideriv_sensitive_0,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 67) deallocate(tlc%c_therm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (fullpol_nst(nst)) then
      if (progression .ge. 66) deallocate(tlc%c_emi_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 65) deallocate(tlc%c_emi_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (fullpol_nst(nst)) then
      if (progression .ge. 64) deallocate(tlc%c_bck_srf_v_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 63) deallocate(tlc%c_bck_srf_0_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 62) deallocate(tlc%c_bck_srf_v_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 61) deallocate(tlc%c_bck_srf_0_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 60) deallocate(tlc%c_bck_srf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (fullpol_nst(nst)) then
      if (progression .ge. 59) deallocate(tlc%diff_c_bck_v_adj_ssa_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 58) deallocate(tlc%diff_c_fwd_v_adj_ssa_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 57) deallocate(tlc%diff_c_bck_0_adj_ssa_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 56) deallocate(tlc%diff_c_fwd_0_adj_ssa_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 55) deallocate(tlc%diff_c_bck_v_nrm_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 54) deallocate(tlc%diff_c_fwd_v_nrm_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 53) deallocate(tlc%diff_c_bck_0_nrm_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 52) deallocate(tlc%diff_c_fwd_0_nrm_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 51) deallocate(tlc%diff_c_bck_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 50) deallocate(tlc%diff_c_fwd_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (fullpol_nst(nst)) then
      if (progression .ge. 49) deallocate(tlc%c_bck_v_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 48) deallocate(tlc%c_fwd_v_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 47) deallocate(tlc%c_bck_0_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 46) deallocate(tlc%c_fwd_0_adj_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 45) deallocate(tlc%c_bck_v_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 44) deallocate(tlc%c_fwd_v_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 43) deallocate(tlc%c_bck_0_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 42) deallocate(tlc%c_fwd_0_nrm,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 41) deallocate(tlc%c_bck,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 40) deallocate(tlc%c_fwd,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 39) deallocate(tlc%diff_c_bck_0v_elem_ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 38) deallocate(tlc%c_bck_0v_elem,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 37) deallocate(tlc%l_int,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 36) deallocate(tlc%l_int_ac_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 35) deallocate(tlc%l_int_ac_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 34) deallocate(tlc%l_int_ac_m,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 33) deallocate(tlc%l_int_ac_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 32) deallocate(tlc%l_int_sm_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 31) deallocate(tlc%l_int_sm_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 30) deallocate(tlc%l_int_sm_m,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 29) deallocate(tlc%l_int_sm_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 28) deallocate(tlc%l_1_sp_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 27) deallocate(tlc%l_sp_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 26) deallocate(tlc%l_1_sp_pv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 25) deallocate(tlc%l_1_sp_p0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 24) deallocate(tlc%l_sp_mp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 23) deallocate(tlc%l_sp_pp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 22) deallocate(tlc%l_sp_mv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 21) deallocate(tlc%l_sp_pv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 20) deallocate(tlc%l_sp_m0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 19) deallocate(tlc%l_sp_p0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 18) deallocate(tlc%l_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 17) deallocate(tlc%l_p,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 16) deallocate(tlc%l_mv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 15) deallocate(tlc%l_pv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 14) deallocate(tlc%l_m0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 13) deallocate(tlc%l_p0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 12) deallocate(tlc%l_0v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 11) deallocate(tlc%t_sp_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 10) deallocate(tlc%t_sp_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 9) deallocate(tlc%t_sp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 8) deallocate(tlc%t_tot_0v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 7) deallocate(tlc%t_tot_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 6) deallocate(tlc%t_tot_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 5) deallocate(tlc%t,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 4) deallocate(tlc%effmu_pv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 3) deallocate(tlc%effmu_mv_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 2) deallocate(tlc%effmu_p0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 1) deallocate(tlc%effmu_m0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fail ! }}}

end module lintran_tlc_module
