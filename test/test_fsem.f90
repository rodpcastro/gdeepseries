module fsem_suite
! Test gds_fsem results against integral expressions evaluated with scipy.
!
! (see file f_points.py)

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use gds_kinds, only: wp
  use gds_fsem, only: fsem
  use stdlib_io_npy, only: load_npy, save_npy

  implicit none
  private
  public :: collect_fsem_suite

contains

  subroutine collect_fsem_suite(testsuite)
  ! Collection of fsem tests.

    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("test_fsem", test_fsem)]
  end subroutine collect_fsem_suite

  subroutine test_fsem(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp), allocatable :: x(:), y(:)
    real(wp), allocatable :: scp_f(:, :), scp_fx(:, :), scp_fxx(:, :)
    real(wp), allocatable :: gds_f(:, :), gds_fx(:, :), gds_fxx(:, :)
    real(wp) :: max_abs_err, max_abs_err_f, max_abs_err_fx, max_abs_err_fxx
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
        call fsem(x(j), y(i), gds_f(i, j), gds_fx(i, j), gds_fxx(i, j))
      end do
    end do

    call cpu_time(time_finish)
    time_elapsed(1) = time_finish - time_start

    print '("Elapsed time = ", f5.3, " seconds")', time_elapsed

    ! call save_npy('test/gds/time.npy', time_elapsed)
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

  end subroutine test_fsem

end module fsem_suite
