!**************************************************************************
      function DEPOL(WAVE)
         use header_module, only: double
         implicit none
         real(double), intent(in) :: wave
         integer, parameter :: NDIM = 23
         real(double) :: rho_out, depol
         real(double), dimension(ndim), parameter :: WAVE_DEPOL = &
                                                     (/300.D0, &
                                                       310.D0, &
                                                       320.D0, &
                                                       330.D0, &
                                                       340.D0, &
                                                       350.D0, &
                                                       360.D0, &
                                                       370.D0, &
                                                       380.D0, &
                                                       390.D0, &
                                                       400.D0, &
                                                       450.D0, &
                                                       500.D0, &
                                                       550.D0, &
                                                       600.D0, &
                                                       650.D0, &
                                                       700.D0, &
                                                       750.D0, &
                                                       800.D0, &
                                                       850.D0, &
                                                       900.D0, &
                                                       950.D0, &
                                                       1000.D0/)

         real(double), dimension(ndim), parameter :: RHO = &
                                                     (/3.178d-02, &
                                                       3.178d-02, &
                                                       3.122d-02, &
                                                       3.066d-02, &
                                                       3.066d-02, &
                                                       3.010d-02, &
                                                       3.010d-02, &
                                                       3.010d-02, &
                                                       2.955d-02, &
                                                       2.955d-02, &
                                                       2.955d-02, &
                                                       2.899d-02, &
                                                       2.842d-02, &
                                                       2.842d-02, &
                                                       2.786d-02, &
                                                       2.786d-02, &
                                                       2.786d-02, &
                                                       2.786d-02, &
                                                       2.730d-02, &
                                                       2.730d-02, &
                                                       2.730d-02, &
                                                       2.730d-02, &
                                                       2.730d-02/)
!----------------------------------------------------------------------------
         call LINTERP( &
            WAVE_DEPOL, RHO, NDIM, &
            WAVE, RHO_OUT, 1)
         if (WAVE .gt. WAVE_DEPOL(NDIM)) RHO_OUT = RHO(NDIM)
         DEPOL = RHO_OUT

         return
      end function depol
