program tester
! Test GDeepSeries results against numerical evaluation with Python.
!
! ## Reference
! 1. The Fortran Programming Language. 2024. test-drive: The simple 
!    testing framework. <https://github.com/fortran-lang/test-drive>
! 2. P. Virtanen et al. 2020. SciPy 1.0: Fundamental Algorithms for
!    Scientific Computing in Python. Nat. Methods 17, 3 (Mar.), 261–272.
!    <https://doi.org/10.1038/s41592-019-0686-2>
! 3. P. A. Brodtkorb. 2025. numdifftools: Solve automatic numerical differentiation
!    problems in one or more variables. <https://github.com/pbrod/numdifftools>.

  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use fsem_suite, only: collect_fsem_suite
  use gdeep_suite, only: collect_gdeep_suite

  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  call ensure_directory()

  stat = 0

  testsuites = [ &
    new_testsuite("fsem_suite", collect_fsem_suite), &
    new_testsuite("gdeep_suite", collect_gdeep_suite) &
  ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

contains

  subroutine ensure_directory()
  ! If nonexistent, directory with gds_fsem results is created.

    implicit none
    logical :: file_exists

    inquire(file='test/gds/f.npy', exist=file_exists)

    if (.not. file_exists) then
      call execute_command_line('mkdir "test/gds"', wait=.true.)
    end if
  end subroutine ensure_directory

end program tester
