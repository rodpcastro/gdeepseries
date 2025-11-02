module gds

  use gds_kinds, only: i1, wp
  use csf, only: ei, bessely0, bessely1, struveh0, struveh1, hyp2f1

  implicit none
  private
  public :: fsem

  real(wp), parameter :: pi = 3.141592653589793_wp
  complex(wp), parameter :: iu = (0.0_wp, 1.0_wp)

contains

  function fsem(x, y) result(f)
    ! Series Expansion Method for F.

    real(wp), intent(in) :: x, y
    real(wp) :: f(3)
    integer(i1) :: nterms

    if (x >= 9.5_wp .and. y > 0.5_wp*x) then
      ! SEM4.
      if (x < 14.0_wp) then      
        nterms = 20
      else
        nterms = 15
      end if

      f = 4.0_wp
      ! f = fsem4(x, y, nterms)

    else if (x >= 6.5_wp .and. y <= 0.5_wp*x) then
      ! SEM3.
      
      f = 3.0_wp
      ! f = fsem3(x, y, 13)

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

      f = fsem2(x, y, nterms)

    else
      ! SEM1.
      if (y > 14.0_wp .and. y < 17.0_wp) then
        nterms = 19
      else
        nterms = 15
      end if
      
      ! f = 1.0_wp
      f = fsem1(x, y, nterms)
    end if
  end function fsem


  real(wp) function expei(x)
    real(wp), intent(in) :: x 
    real(wp) :: xi, sk
    integer(i1) :: k

    xi = 1.0_wp/x

    if (x > 40.0_wp) then
      sk = 1.0_wp
      expei = 1.0_wp
      do k = 1, 23
        sk = sk * k * xi
        expei = expei + sk
      end do
      expei = expei * xi
    else
      expei = exp(-x) * ei(x)
    end if
  end function expei


  function fsem1(x, y, nterms) result(f)
    ! Series Expansion Method 1.

    real(wp), intent(in) :: x, y
    integer(i1), intent(in) :: nterms
    real(wp) :: f(3)
    real(wp) :: eey
    real(wp) :: xi, yi
    real(wp) :: qxi, qx2
    real(wp) :: tn1, tn2, tn3, sn1, sn2, sn3
    real(wp) :: tm, sm
    integer(i1) :: n, m

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
      tn1 = -tn1 * qx2/(n*n)
      tn2 = tn1 * n
      tn3 = tn2 * (2*n-1)

      tm = yi
      sm = yi
      do m = 2, 2*n
        tm = tm * (m-1)*yi
        sm = sm + tm
      end do
      sm = sm - eey

      sn1 = sn1 + tn1*sm
      sn2 = sn2 + tn2*sm
      sn3 = sn3 + tn3*sm
    end do

    f(1) = 2.0_wp * (sn1 - eey)
    f(2) = qxi * sn2
    f(3) = qxi*xi * sn3
  end function fsem1


  function fsem2(x, y, nterms) result(f)
    ! Series Expansion Method 2.

    real(wp), intent(in) :: x, y
    integer(i1), intent(in) :: nterms
    real(wp) :: f(3)
    real(wp) :: x2, y2, r2, r, xi, ri
    real(wp) :: rxi, rxi2, r2xi2, yxiri
    real(wp) :: ey, py, py0, py1
    real(wp) :: tn, sn1, sn2, sn3
    real(wp) :: hg, tg
    complex(wp) :: a, b, c, z
    integer(i1) :: n

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

    a = (0.5_wp, 0.0_wp)
    c = (1.5_wp, 0.0_wp)
    z = complex(r2xi2, 0.0_wp)

    tn = 1.0_wp
    sn1 = 1.0_wp
    sn2 = 1.0_wp
    sn3 = 0.0_wp
    do n = 1, nterms
      tn = tn * x/n
      b = complex(-0.5_wp*n, 0.0_wp)
      hg = real(iu**-n * hyp2f1(a, b, c, z))
      tg = tn * hg

      sn1 = sn1 + (n+1) * tg
      sn2 = sn2 + tg
      sn3 = sn3 + n * tg
    end do

    f(1) = -py0 + 2.0_wp * rxi2 * (ey*sn1 - 1.0_wp)
    f(2) =  py1 + 2.0_wp * (yxiri - rxi*ey*sn2)
    f(3) =  py0 - py1*xi + 2.0_wp * (yxiri*xi*(y2/r2 - 2.0_wp + y) - rxi2*ey*sn3)
  end function fsem2


end module gds
