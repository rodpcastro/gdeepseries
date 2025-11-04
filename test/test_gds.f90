module test_suite
! Test GDeepSeries results against integral expressions evaluated with scipy.

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use gds, only: fsem
  use stdlib_io_npy, only: load_npy, save_npy
  use gds_kinds, only: wp

  implicit none
  private
  public :: collect_gds_suite

contains

  subroutine collect_gds_suite(testsuite)
  ! Collection of tests.

    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("test_gds", test_gds)]
  end subroutine collect_gds_suite

  subroutine test_gds(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp), allocatable :: x(:), y(:)
    real(wp), allocatable :: scp_f(:, :), scp_fx(:, :), scp_fxx(:, :)
    real(wp), allocatable :: gds_f(:, :), gds_fx(:, :), gds_fxx(:, :)
    real(wp) :: f(3), max_abs_err, max_abs_err_f, max_abs_err_fx, max_abs_err_fxx
    real(wp) :: time_start, time_finish, time_elapsed(1)
    integer :: i, j, nx, ny

    call load_npy('test/scipy/x.npy', x)
    call load_npy('test/scipy/y.npy', y)

    nx = size(x)
    ny = size(y)
    allocate(gds_f(nx, ny))
    allocate(gds_fx(nx, ny))
    allocate(gds_fxx(nx, ny))

    call cpu_time(time_start)

    do j = 1, nx
      do i = 1, ny
        f = fsem(x(j), y(i))
        gds_f(i, j) = f(1)
        gds_fx(i, j) = f(2)
        gds_fxx(i, j) = f(3)
      end do
    end do

    call cpu_time(time_finish)
    time_elapsed(1) = time_finish - time_start

    print '("Elapsed time = ", f5.3, " seconds")', time_elapsed

    call save_npy('test/gds/time.npy', time_elapsed)
    call save_npy('test/gds/f.npy', gds_f)
    call save_npy('test/gds/fx.npy', gds_fx)
    call save_npy('test/gds/fxx.npy', gds_fxx)

    call load_npy('test/scipy/f.npy', scp_f)
    call load_npy('test/scipy/fx.npy', scp_fx)
    call load_npy('test/scipy/fxx.npy', scp_fxx)

    max_abs_err_f = maxval(abs(scp_f - gds_f))
    max_abs_err_fx = maxval(abs(scp_fx - gds_fx))
    max_abs_err_fxx = maxval(abs(scp_fxx - gds_fxx))
    max_abs_err = max(max_abs_err_f, max_abs_err_fx, max_abs_err_fxx)

    call check(error, max_abs_err < 1.0e-8_wp)

  end subroutine test_gds

end module test_suite
