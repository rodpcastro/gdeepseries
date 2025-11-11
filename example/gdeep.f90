program gdeep_demo
! PLot of infinite-depth free-surface Green function.

  use gds_kinds, only: wp
  use gds, only: gdeep

  implicit none
  real(wp) :: f(3), s(3), k0
  complex(wp) :: g, gradg(3), hessg(3, 3)

  f = [0.0_wp, 0.0_wp, -1.0_wp]
  s = [1.0_wp, 1.0_wp, -2.0_wp]
  k0 = 1.0_wp

  call gdeep(f, s, k0, g, gradg, hessg)

  print *, 'g_re', real(g)
  print *, '--------------------------------------------------------------------------'
  print *, 'grad_re = ', real(gradg)
  print *, '--------------------------------------------------------------------------'
  print *, 'hessg_re = '
  print *, real(hessg(1, :))
  print *, real(hessg(2, :))
  print *, real(hessg(3, :))
  print *, '--------------------------------------------------------------------------'
  print *, 'g_im', imag(g)
  print *, '--------------------------------------------------------------------------'
  print *, 'grad_im = ', imag(gradg)
  print *, '--------------------------------------------------------------------------'
  print *, 'hessg_im = '
  print *, imag(hessg(1, :))
  print *, imag(hessg(2, :))
  print *, imag(hessg(3, :))

end program gdeep_demo
