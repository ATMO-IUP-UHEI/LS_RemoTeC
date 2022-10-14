! Use intrinsic funtion to get numerical precision
   integer, parameter :: single = selected_real_kind(6,37)
   integer, parameter :: double = selected_real_kind(15,307)
   integer, parameter :: quad = selected_real_kind(33,4931)
! if quad precision is not supported
!   integer, parameter :: quad = selected_real_kind(15,307)
