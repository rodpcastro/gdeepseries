module gds

  use gds_kinds, only: i2, wp
  use gds_constants, only: pi
  use gds_sem, only: fsem
  use csf, only: besselj0, besselj1

  implicit none
  private
  public :: gsem

contains

  subroutine gsem(p, q, k0, g, gradg, hessg)
    ! Infinite-depth free-surface Green function.

    real(wp), intent(in) :: p(3), q(3), k0
    complex(wp), intent(out) :: g, gradg(3), hessg(3, 3)
    real(wp) :: du, dv, z, du2, dv2, d, di, di3
    real(wp) :: x, y, x2, y2, yt, yt2
    real(wp) :: r, ri, ri3, ri5, rq, rqi, rqi3, rqi5
    real(wp) :: dk, dj0, dj1
    real(wp) :: xu, xv, xuu, xvv, xuv, yw, yww
    real(wp) :: fs(3), f, fx, fy, fxx, fyy, fxy
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

    dk  = 2.0_wp*pi*exp(-y)
    dj0 = dk*besselj0(x)
    dj1 = dk*besselj1(x)

    xu  = k0*du*di
    xv  = k0*dv*di
    xuu = k0*dv2*di3
    xvv = k0*du2*di3
    xuv = -k0*du*dv*di3
    yw  = -k0
    yww = 0.0_wp

    fs  = fsem(x, y)
    f   = fs(1)
    fx  = fs(2)
    fy  = -2.0_wp*ri - f
    fxx = fs(3)
    fyy = 2.0_wp*(y*ri3 + ri) + f
    fxy = 2.0_wp*x*ri3 - fx

    g   = k0*complex(rqi + ri + f, -dj0)
    gx  = k0*complex(-x*(rqi3 + ri3) + fx, dj1) 
    gy  = k0*complex(-yt*rqi3 - y*ri3 + fy, dj0)
    gxx = k0*complex(3.0_wp*x2*(rqi5 + ri5) - rqi3 - ri3 + fxx, dj0 - dj1/x)
    gyy = k0*complex(3.0_wp*(yt2*rqi5 + y2*ri5) - rqi3 - ri3 + fyy, -dj0)
    gxy = k0*complex(3.0_wp*x*(yt*rqi5 + y*ri5) + fxy, -dj1)

    gradg(1) = gx*xu
    gradg(2) = gx*xv
    gradg(3) = gy*yw

    hessg(1,1) = gx*xuu + gxx*xu**2
    hessg(2,2) = gx*xvv + gxx*xv**2
    hessg(3,3) = gy*yww + gyy*yw**2
    hessg(2,3) = gxy*xv*yw
    hessg(1,3) = gxy*xu*yw
    hessg(1,2) = gx*xuv + gxx*xu*xv
    hessg(2,1) = hessg(1,2)
    hessg(3,1) = hessg(1,3)
    hessg(3,2) = hessg(2,3)
  end subroutine gsem

end module gds
