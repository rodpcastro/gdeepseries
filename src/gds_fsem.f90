module gds_fsem
! Series expansions for F, Fx and Fxx.

  use gds_kinds, only: i2, wp
  use gds_constants, only: pi
  use csf, only: ei, bessely0, bessely1, struveh0, struveh1, hyp2f1

  implicit none
  private
  public :: fsem

contains

  subroutine fsem(x, y, f, fx, fxx)
    ! Series Expansion Method for F, Fx and Fxx.

    real(wp), intent(in) :: x, y
    real(wp), intent(out) :: f, fx, fxx
    integer(i2) :: nterms

    if (x >= 9.5_wp .and. y > 0.5_wp*x) then
      ! SEM4.
      if (x < 14.0_wp) then      
        nterms = 20
      else
        nterms = 15
      end if

      call fsem4(x, y, nterms, f, fx, fxx)
    else if (x >= 6.5_wp .and. y <= 0.5_wp*x) then
      ! SEM3.
      nterms = 13
      
      call fsem3(x, y, nterms, f, fx, fxx)
    else if (y < 15.0_wp .and. y < 2.0_wp*x) then
      ! SEM2.
      if (y > 11.0_wp) then
        nterms = 42
      else if (x > 7.0_wp .or. y > 8.0_wp) then
        nterms = 36
      else if (x > 4.5_wp .or. y > 6.0_wp) then
        nterms = 30
      else
        nterms = 26
      end if

      call fsem2(x, y, nterms, f, fx, fxx)
    else
      ! SEM1.
      if (y > 14.0_wp .and. y < 17.0_wp) then
        nterms = 19
      else
        nterms = 15
      end if
      
      call fsem1(x, y, nterms, f, fx, fxx)
    end if
  end subroutine fsem


  real(wp) function expei(x)
    ! exp(-x) * Ei(x).

    real(wp), intent(in) :: x 
    real(wp) :: xi, sk
    integer(i2) :: k

    xi = 1.0_wp/x

    if (x > 40.0_wp) then
      sk = 1.0_wp
      expei = 1.0_wp
      do k = 1, 23
        sk = sk*k*xi
        expei = expei + sk
      end do
      expei = expei*xi
    else
      expei = exp(-x)*ei(x)
    end if
  end function expei


  subroutine fsem1(x, y, nterms, f, fx, fxx)
    ! Series Expansion Method 1.

    real(wp), intent(in) :: x, y
    integer(i2), intent(in) :: nterms
    real(wp), intent(out) :: f, fx, fxx
    real(wp) :: eey
    real(wp) :: xi, yi
    real(wp) :: qxi, qx2
    real(wp) :: tn1, tn2, tn3, sn1, sn2, sn3
    real(wp) :: tm, sm
    integer(i2) :: n, m

    eey = expei(y)

    xi = 1.0_wp/x
    yi = 1.0_wp/y

    qxi = 4.0_wp*xi
    qx2 = 0.25_wp*x*x

    tn1 = 1.0_wp
    sn1 = 0.0_wp
    sn2 = 0.0_wp
    sn3 = 0.0_wp
    do n = 1, nterms
      tn1 = -tn1*qx2/(n*n)
      tn2 = tn1*n
      tn3 = tn2*(2*n-1)

      tm = yi
      sm = yi
      do m = 2, 2*n
        tm = tm*(m-1)*yi
        sm = sm + tm
      end do
      sm = sm - eey

      sn1 = sn1 + tn1*sm
      sn2 = sn2 + tn2*sm
      sn3 = sn3 + tn3*sm
    end do

    f   = 2.0_wp*(sn1 - eey)
    fx  = qxi*sn2
    fxx = qxi*xi*sn3
  end subroutine fsem1


  subroutine fsem2(x, y, nterms, f, fx, fxx)
    ! Series Expansion Method 2.

    real(wp), intent(in) :: x, y
    integer(i2), intent(in) :: nterms
    real(wp), intent(out) :: f, fx, fxx
    real(wp) :: x2, y2, r2, r, xi, ri
    real(wp) :: rxi, rxi2, r2xi2, yxiri
    real(wp) :: ey, py, py0, py1
    real(wp) :: tn, sn1, sn2, sn3
    real(wp) :: hg, tg
    complex(wp) :: in, ha, hb, hc, hz
    integer(i2) :: n

    x2 = x*x
    y2 = y*y
    r2 = x2 + y2
    r = sqrt(r2)
    xi = 1.0_wp/x
    ri = 1.0_wp/r
    
    rxi = r*xi
    rxi2 = rxi*xi
    r2xi2 = rxi2*r
    yxiri = y*xi*ri

    ey = exp(-y)
    py = pi * ey
    py0 = py * bessely0(x)
    py1 = py * bessely1(x)

    in = (1.0_wp, 0.0_wp)
    ha = (0.5_wp, 0.0_wp)
    hc = (1.5_wp, 0.0_wp)
    hz = complex(r2xi2, 0.0_wp)

    tn = 1.0_wp
    sn1 = 1.0_wp
    sn2 = 1.0_wp
    sn3 = 0.0_wp
    do n = 1, nterms
      tn = tn*x/n
      in = in*(0.0_wp, -1.0_wp)
      hb = complex(-0.5_wp*n, 0.0_wp)
      hg = real(in * hyp2f1(ha, hb, hc, hz), wp)
      tg = tn * hg

      sn1 = sn1 + (n+1)*tg
      sn2 = sn2 + tg
      sn3 = sn3 + n*tg
    end do

    f   = -py0 + 2.0_wp*rxi2*(ey*sn1 - 1.0_wp)
    fx  =  py1 + 2.0_wp*(yxiri - rxi*ey*sn2)
    fxx =  py0 - py1*xi + 2.0_wp*(yxiri*xi*(y2/r2 - 2.0_wp + y) - rxi2*ey*sn3)
  end subroutine fsem2


  subroutine fsem3(x, y, nterms, f, fx, fxx)
    ! Series Expansion Method 3.

    real(wp), intent(in) :: x, y
    integer(i2), intent(in) :: nterms
    real(wp), intent(out) :: f, fx, fxx
    real(wp) :: xi, xi2, xi3, y2, yi, hxi2 
    real(wp) :: ey, py, oy, phy0, phy1
    real(wp) :: y2n, cn, tn, sn1, sn2, sn3
    real(wp) :: nc1, nc2, nc3
    integer(i2) :: n, dn

    xi = 1.0_wp/x
    xi2 = xi*xi
    xi3 = xi2*xi
    y2 = y*y
    yi = 1.0_wp/y
    hxi2 = 0.5_wp*xi2
    
    ey = exp(-y)
    py = pi*ey
    oy = 1.0_wp - ey
    phy0 = py*(struveh0(x) + bessely0(x))
    phy1 = py*(struveh1(x) + bessely1(x))

    y2n = 1.0_wp
    cn = oy
    tn = 1.0_wp
    sn1 = 0.0_wp
    sn2 = 0.0_wp
    sn3 = 0.0_wp
    do n = 1, nterms
      dn = 2*n
      tn = -tn * hxi2 * (dn-1)/n
      y2n = y2n * y2
      cn = y2n*(1.0_wp - dn*yi) + dn*(dn-1)*cn

      nc1 = tn * cn
      nc2 = (dn+1) * nc1
      nc3 = (dn+2) * nc2

      sn1 = sn1 + nc1
      sn2 = sn2 + nc2
      sn3 = sn3 + nc3
    end do

    f   = -phy0 - 2.0_wp*xi*(oy + sn1)
    fx  =  phy1 - 2.0_wp*ey + 2.0_wp*xi2*(oy + sn2)
    fxx =  phy0 - xi*phy1 - 2.0_wp*xi3*(2.0_wp*oy + sn3)
  end subroutine fsem3


  subroutine fsem4(x, y, nterms, f, fx, fxx)
    ! Series Expansion Method 4.

    real(wp), intent(in) :: x, y
    integer(i2), intent(in) :: nterms
    real(wp), intent(out) :: f, fx, fxx
    real(wp) :: x2, y2, r2, r, xi, yi, yi2, ri, ri2, ri3
    real(wp) :: x2ri2, y2ri2, hyr
    real(wp) :: ey, py, oy, eyi, phy0, phy1
    real(wp) :: b, b0, b1, bn
    real(wp) :: tn, sn1, sn2, sn3
    real(wp) :: nb1, nb2, nb3
    integer(i2) :: n, dn

    x2 = x*x
    y2 = y*y
    r2 = x2 + y2
    r = sqrt(r2)
    xi = 1.0_wp/x
    yi = 1.0_wp/y
    yi2 = yi*yi
    ri = 1.0_wp/r
    ri2 = ri*ri
    ri3 = ri2*ri

    x2ri2 = x2*ri2
    y2ri2 = y2*ri2
    hyr = 0.5_wp*y2ri2
    
    ey = exp(-y)
    py = pi*ey
    oy = 1.0_wp - ey
    eyi = ey*yi
    phy0 = py*(struveh0(x) + bessely0(x))
    phy1 = py*(struveh1(x) + bessely1(x))

    b = 1.0_wp
    b0 = oy*yi
    b1 = yi*((1.0_wp - 2.0_wp*yi2)*ey - 2.0_wp*(yi - yi2))

    tn = -hyr
    sn1 = tn*b1
    sn2 = 3.0_wp*sn1
    sn3 = sn2*(1.0_wp - 5.0_wp*x2ri2)
    do n = 2, nterms
      dn = 2*n
      tn = -tn*hyr*(dn-1)/n
      b = -1.0_wp*b
      bn = b*eyi + yi2*dn*((dn-1)*b1 + (dn-2)*b0)
      b0 = b1
      b1 = bn

      nb1 = tn*bn
      nb2 = nb1*(dn+1)
      nb3 = nb2*(1.0_wp - (dn+3)*x2ri2)

      sn1 = sn1 + nb1
      sn2 = sn2 + nb2
      sn3 = sn3 + nb3
    end do

    f   = -phy0 - 2.0_wp*ri*(oy + y*sn1)
    fx  =  phy1 - 2.0_wp*ey + 2.0_wp*x*ri3*(oy + y*sn2)
    fxx =  phy0 - xi*phy1 + 2.0_wp*ri3*(oy*(1.0_wp - 3.0_wp*x2ri2) + y*sn3)
  end subroutine fsem4

end module gds_fsem
