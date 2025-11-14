program gds_example
! Example program.

  use gds_kinds, only: wp
  use gds, only: gdeep

  implicit none
  real(wp) :: p(3), q(3), k0
  complex(wp) :: g, gradg(3), hessg(3, 3)
  integer :: i, j
  real(wp) :: re, im
  character(1) :: ss
  character(37) :: fmt

  fmt = 'sp, es10.3e2, 1x, a, " i", s, es9.3e2'

  ! Field point, source point and wave number.
  p = [0.0_wp, 0.0_wp, -1.0_wp]
  q = [1.0_wp, 1.0_wp, -2.0_wp]
  k0 = 1.0_wp

  print '(a, sp, f4.1, a, f4.1, a, f4.1, "]")', 'p = [', p(1), ', ', p(2), ', ', p(3)
  print '(a, sp, f4.1, a, f4.1, a, f4.1, "]")', 'q = [', q(1), ', ', q(2), ', ', q(3)
  print '(a, f3.1)', 'k0 = ', k0

  call gdeep(p, q, k0, g, gradg, hessg)

  ! Green function.
  re = real(g)
  im = abs(imag(g))
  ss = sign_symbol(imag(g))
  print '(a, '//fmt//')', new_line('n')//'g = ', re, ss, im

  ! Gradient of the Green function.
  print '(a, $)', new_line('n')//'gradg = ['
  do i = 1, 3
    re = real(gradg(i))
    im = abs(imag(gradg(i)))
    ss = sign_symbol(imag(gradg(i)))
    print '('//fmt//', $)', re, ss, im 
    if (i < 3) print '(a, $)', ', '
  end do
  print '(a)', ']'

  ! Hessian of the Green function.
  print '(a, $)', new_line('n')//'hessg = ['
  do i = 1, 3
    if (i == 1) then
      print '(a, $)', '['
    else
      print '(9x, a, $)', '['
    end if

    do j = 1, 3
      re = real(hessg(i,j))
      im = abs(imag(hessg(i,j)))
      ss = sign_symbol(imag(hessg(i,j)))
      print '('//fmt//', $)', re, ss, im 
      if (j < 3) print '(a, $)', ', '
    end do

    print '(a, $)', ']'
    if (i < 3) then
      print '(a)', ','
    else
      print '(a)', ']'
    end if
  end do

contains

  character function sign_symbol(x)
    real(wp), intent(in) :: x
    sign_symbol = merge('+', '-', x >= 0.0_wp)
  end function sign_symbol

end program gds_example
