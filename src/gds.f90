module gds
! Infinite-depth free-surface Green function.

  use gds_kinds, only: wp, i1
  use gds_constants, only: pi
  use gds_fsem, only: fsem
  use csf, only: besselj0, besselj1

  implicit none
  private
  public :: gdeep

contains

  subroutine gdeep(p, q, k0, g, gradg, hessg, s)
    ! Infinite-depth free-surface Green function.
    ! 
    ! Parameters
    ! ----------
    ! p : 8-byte array(3)
    !   Field point.
    ! q : 8-byte array(3)
    !   Source point.
    ! k0 : 8-byte value
    !   Infinite-depth wavenumber.
    ! s : 1-byte integer, default=+1
    !   s = +1 for time component e^{+iωt}, and
    !   s = -1 for time component e^{-iωt}.
    !
    ! Returns
    ! -------
    ! g : 8-byte value
    !   Green function.
    ! gradg : 8-byte array(3)
    !   Green function gradient.
    ! hessg : 8-byte array(3,3)
    !   Green function Hessian matrix.

    real(wp), intent(in) :: p(3), q(3), k0
    complex(wp), intent(out) :: g, gradg(3), hessg(3, 3)
    integer(i1), intent(in), optional :: s  ! default=+1
    integer(i1) :: s_
    real(wp) :: du, dv, z, du2, dv2, d, di, di3
    real(wp) :: x, y, x2, y2, yt, yt2
    real(wp) :: r, ri, ri3, ri5, rq, rqi, rqi3, rqi5
    real(wp) :: dy, dj0, dj1
    real(wp) :: xu, xv, xuu, xvv, xuv, yw, yww
    real(wp) :: f, fx, fy, fxx, fyy, fxy
    complex(wp) :: gx, gy, gxx, gyy, gxy

    du = p(1) - q(1)
    dv = p(2) - q(2)
    z  = p(3) + q(3)
    du2 = du*du
    dv2 = dv*dv
    d = sqrt(du2 + dv2)
    di = 1.0_wp/d
    di3 = di**3

    x = k0*d
    y = -k0*z
    yt = y + 2.0_wp*k0*q(3)
    x2 = x*x
    y2 = y*y
    yt2 = yt*yt

    r = sqrt(x2 + y2)
    ri = 1.0_wp/r
    ri3 = ri**3
    ri5 = ri3*ri**2
    rq = sqrt(x2 + yt2)
    rqi = 1.0_wp/rq
    rqi3 = rqi**3
    rqi5 = rqi3*rqi**2

    s_ = 1_i1
    if (present(s)) s_ = sign(1_i1, s)

    dy  = -2.0_wp*pi*s_*exp(-y)
    dj0 = dy*besselj0(x)
    dj1 = dy*besselj1(x)

    xu  = k0*du*di
    xv  = k0*dv*di
    xuu = k0*dv2*di3
    xvv = k0*du2*di3
    xuv = -k0*du*dv*di3
    yw  = -k0
    ! yww = 0.0_wp

    call fsem(x, y, f, fx, fxx)
    fy  = -2.0_wp*ri - f
    fyy = 2.0_wp*(y*ri3 + ri) + f
    fxy = 2.0_wp*x*ri3 - fx

    g   = k0*complex(rqi + ri + f, dj0)
    gx  = k0*complex(-x*(rqi3 + ri3) + fx, -dj1) 
    gy  = k0*complex(-yt*rqi3 - y*ri3 + fy, -dj0)
    gxx = k0*complex(3.0_wp*x2*(rqi5 + ri5) - rqi3 - ri3 + fxx, -dj0 + dj1/x)
    gyy = k0*complex(3.0_wp*(yt2*rqi5 + y2*ri5) - rqi3 - ri3 + fyy, dj0)
    gxy = k0*complex(3.0_wp*x*(yt*rqi5 + y*ri5) + fxy, dj1)

    gradg(1) = gx*xu
    gradg(2) = gx*xv
    gradg(3) = gy*yw

    hessg(1,1) = gx*xuu + gxx*xu**2
    hessg(2,2) = gx*xvv + gxx*xv**2
    hessg(3,3) = gyy*yw**2  ! yww = 0
    hessg(2,3) = gxy*xv*yw
    hessg(1,3) = gxy*xu*yw
    hessg(1,2) = gx*xuv + gxx*xu*xv
    hessg(2,1) = hessg(1,2)
    hessg(3,1) = hessg(1,3)
    hessg(3,2) = hessg(2,3)
  end subroutine gdeep

end module gds
