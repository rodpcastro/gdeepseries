module gds

  use gds_kinds, only: i2, wp
  use gds_constants, only: pi
  use csf, only: besselj0, besselj1

  implicit none
  private
  public :: gsem

contains

  subroutine gsem(p, q, k0, g, gradg, hessg)
    ! Infinite-depth free-surface Green function.

    real(wp), intent(in) :: p(3), q(3), k0
    real(wp), intent(out) :: g, gradg(2), hessg(2, 2)
    real(wp) :: r


  end function gsem

end module gds
