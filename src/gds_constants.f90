module gds_constants
! Constants.

  use gds_kinds, only: wp

  implicit none
  private
  public :: pi

  real(wp), parameter :: pi = 3.141592653589793_wp  ! π

end module gds_constants
