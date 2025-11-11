program gds_demo
! PLot of F, Fx and Fxx evaluated through series expansions.

  use gds_kinds, only: wp
  use pyplot_module
  use gds_fsem, only: fsem

  implicit none
  type(pyplot) :: plt
  integer :: i, j
  integer, parameter :: nx = 200, ny = 200
  real(wp) :: xmin, xmax, ymin, ymax, dx, dy
  real(wp) :: x(nx), y(ny), z(3), f(nx, ny), fx(nx, ny), fxx(nx, ny)

  xmin = 0.1_wp
  xmax = 40.0_wp
  dx = (xmax - xmin) / (nx - 1)

  ymin = 0.1_wp
  ymax = 40.0_wp
  dy = (ymax - ymin) / (nx - 1)

  x = [(0.1_wp+i*dx, i = 0, nx-1)]
  y = [(0.1_wp+i*dy, i = 0, ny-1)]

  do j = 1, ny
    do i = 1, nx
      z = fsem(x(i), y(j))
      f(i, j) = z(1)
      fx(i, j) = z(2)
      fxx(i, j) = z(3)
    end do
  end do
  
  call plt%initialize(mplot3d=.true., xlabel='X', ylabel='Y', title='$F$', figsize=[7, 7])
  call plt%plot_surface(x, y, f, label='', linestyle='', cmap='bone')
  call plt%savefig('example/f.svg')

end program gds_demo
