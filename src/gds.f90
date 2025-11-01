module gds

  use gds_kinds, only: i1, wp

  implicit none
  private
  public :: fsem

contains

  real(wp) function fsem(x, y) result(f)
    ! Series Expansion Method for F.

    real(wp), intent(in) :: x, y
    integer(i1) :: nterms

  end function fsem

end module gds
