module gdeep_suite
! Test gdeep results against numerical evaluation with Python.
! 
! In Pythom, The derivatives of the Green function were evalauted
! with the library numdifftools (see file g_points.py).

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use gds_kinds, only: wp
  use gds, only: gdeep

  implicit none
  private
  public :: collect_gdeep_suite

contains

  subroutine collect_gdeep_suite(testsuite)
  ! Collection of gdeep tests.

    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("test_gdeep_1", test_gdeep_1), &
      new_unittest("test_gdeep_2", test_gdeep_2), &
      new_unittest("test_gdeep_3", test_gdeep_3), &
      new_unittest("test_gdeep_4", test_gdeep_4) &
    ]
  end subroutine collect_gdeep_suite

  subroutine test_gdeep_1(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: p(3), q(3), k0
    complex(wp) :: g, gradg(3), hessg(3,3)
    complex(wp) :: g_py, gradg_py(3), hessg_py(3,3)
    real(wp) :: abs_err_g, abs_err_gradg, abs_err_hessg, max_abs_err

    p = [1.1_wp,-0.1_wp,-0.5_wp]
    q = [1.2_wp, 2.7_wp,-6.9_wp]
    k0 = 1.5_wp

    call gdeep(p, q, k0, g, gradg, hessg)

    g_py = (-6.7809868895297165e-03_wp, 5.3580079266483052e-05_wp)
    gradg_py = [(-5.3771901257790847e-05_wp, 1.0642379607385926e-06_wp), &
                (-1.5056132352270673e-03_wp, 2.9798662900810458e-05_wp), &
                (-3.9134798236854992e-02_wp, 8.0370118899725600e-05_wp)]
    hessg_py(1,:) = [( 5.3687988086155696e-04_wp,-1.0768838740010434e-05_wp), &
                     (-2.3495690304779960e-05_wp,-3.5408557138879472e-06_wp), &
                     (-2.6153362780435461e-04_wp, 1.5963569411154053e-06_wp)]
    hessg_py(2,:) = [(-2.3495690304779960e-05_wp,-3.5408557138879472e-06_wp), &
                     (-1.2016031774210254e-04_wp,-1.0978633960948051e-04_wp), &
                     (-7.3229415783088126e-03_wp, 4.4697994351162568e-05_wp)]
    hessg_py(3,:) = [(-2.6153362780435461e-04_wp, 1.5963569411154053e-06_wp), &
                     (-7.3229415783088126e-03_wp, 4.4697994351162568e-05_wp), &
                     (-4.1671956298973460e-04_wp, 1.2055517834943032e-04_wp)]

    abs_err_g = abs(g_py - g)
    abs_err_gradg = maxval(abs(gradg_py - gradg))
    abs_err_hessg = maxval(abs(hessg_py - hessg))
    max_abs_err = max(abs_err_g, abs_err_gradg, abs_err_hessg)

    call check(error, max_abs_err < 1.0e-8_wp)

  end subroutine test_gdeep_1


  subroutine test_gdeep_2(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: p(3), q(3), k0
    complex(wp) :: g, gradg(3), hessg(3,3)
    complex(wp) :: g_py, gradg_py(3), hessg_py(3,3)
    real(wp) :: abs_err_g, abs_err_gradg, abs_err_hessg, max_abs_err

    p = [-0.8_wp, 1.5_wp,-1.9_wp]
    q = [-3.7_wp, 2.2_wp,-2.7_wp]
    k0 = 2.3_wp

    call gdeep(p, q, k0, g, gradg, hessg)

    g_py = (1.1445265127935300e-01_wp,-1.0893062601211369e-04_wp)
    gradg_py = [(-7.1890369969213330e-02_wp,-3.8239973165366140e-05_wp), &
                ( 1.7352847923588155e-02_wp, 9.2303383502651568e-06_wp), &
                (-6.1148976466481597e-02_wp,-2.5054043982784970e-04_wp)]
    hessg_py(1,:) = [( 5.4532356798965259e-02_wp, 5.5625150706992577e-04_wp), &
                     (-1.9146723338670311e-02_wp,-1.3745048044769392e-04_wp), &
                     ( 3.7224446231961653e-02_wp,-8.7951938280277677e-05_wp)]
    hessg_py(2,:) = [(-1.9146723338670311e-02_wp,-1.3745048044769392e-04_wp), &
                     (-2.0168159873677434e-02_wp, 1.9991504534059231e-05_wp), &
                     (-8.9852111590409140e-03_wp, 2.1229778205417444e-05_wp)]
    hessg_py(3,:) = [( 3.7224446231961653e-02_wp,-8.7951938280277677e-05_wp), &
                     (-8.9852111590409140e-03_wp, 2.1229778205417444e-05_wp), &
                     (-3.4364196926654156e-02_wp,-5.7624301160574896e-04_wp)]

    abs_err_g = abs(g_py - g)
    abs_err_gradg = maxval(abs(gradg_py - gradg))
    abs_err_hessg = maxval(abs(hessg_py - hessg))
    max_abs_err = max(abs_err_g, abs_err_gradg, abs_err_hessg)

    call check(error, max_abs_err < 1.0e-8_wp)

  end subroutine test_gdeep_2


  subroutine test_gdeep_3(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: p(3), q(3), k0
    complex(wp) :: g, gradg(3), hessg(3,3)
    complex(wp) :: g_py, gradg_py(3), hessg_py(3,3)
    real(wp) :: abs_err_g, abs_err_gradg, abs_err_hessg, max_abs_err

    p = [ 0.9_wp, 0.3_wp,-1.6_wp]
    q = [-4.1_wp, 3.9_wp,-0.4_wp]
    k0 = 1.3_wp

    call gdeep(p, q, k0, g, gradg, hessg)

    g_py = (-1.3892053879796012e-01_wp,-1.0277765211808870e-01_wp)
    gradg_py = [(-9.8779740279395345e-02_wp, 1.5103679791230373e-01_wp), &
                ( 7.1121413001165062e-02_wp,-1.0874649449684864e-01_wp), &
                (-1.7480350912471193e-01_wp,-1.3361094775350962e-01_wp)]
    hessg_py(1,:) = [( 1.5788210037717715e-01_wp, 1.0481188600696195e-01_wp), &
                     (-1.2789939487363114e-01_wp,-5.3715259022767463e-02_wp), &
                     (-1.3052319765532061e-01_wp, 1.9634783728594013e-01_wp)]
    hessg_py(2,:) = [(-1.2789939487363114e-01_wp,-5.3715259022767463e-02_wp), &
                     ( 7.2331616253212064e-02_wp, 6.8882346078561960e-02_wp), &
                     ( 9.3976702312389457e-02_wp,-1.4137044284582750e-01_wp)]
    hessg_py(3,:) = [(-1.3052319765532061e-01_wp, 1.9634783728594013e-01_wp), &
                     ( 9.3976702312389457e-02_wp,-1.4137044284582750e-01_wp), &
                     (-2.3021371663197229e-01_wp,-1.7369423207843068e-01_wp)]

    abs_err_g = abs(g_py - g)
    abs_err_gradg = maxval(abs(gradg_py - gradg))
    abs_err_hessg = maxval(abs(hessg_py - hessg))
    max_abs_err = max(abs_err_g, abs_err_gradg, abs_err_hessg)

    call check(error, max_abs_err < 1.0e-8_wp)

  end subroutine test_gdeep_3


  subroutine test_gdeep_4(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: p(3), q(3), k0
    complex(wp) :: g, gradg(3), hessg(3,3)
    complex(wp) :: g_py, gradg_py(3), hessg_py(3,3)
    real(wp) :: abs_err_g, abs_err_gradg, abs_err_hessg, max_abs_err

    p = [ 2.6_wp,-1.4_wp,-3.8_wp]
    q = [-2.1_wp,-1.8_wp,-0.2_wp]
    k0 = 2.9_wp

    call gdeep(p, q, k0, g, gradg, hessg)

    g_py = (-4.9581412538379221e-03_wp,-3.4208299664331898e-05_wp)
    gradg_py = [(1.9442879934195743e-03_wp, 3.6210192979452399e-05_wp), &
                (1.6547131858232246e-04_wp, 3.0817185514503561e-06_wp), &
                (-6.2361962456376405e-05_wp,-9.9204069026488959e-05_wp)]
    hessg_py(1,:) = [(-8.6639690746933315e-04_wp, 2.7802951219409460e-04_wp), &
                     (-1.0894256956945698e-04_wp, 2.3006401346335526e-05_wp), &
                     ( 1.2595237108974458e-04_wp, 1.0500955964021795e-04_wp)]
    hessg_py(2,:) = [(-1.0894256956945698e-04_wp, 2.3006401346335526e-05_wp), &
                     ( 4.0440658842788921e-04_wp, 9.6622879824793788e-06_wp), &
                     ( 1.0719350671889904e-05_wp, 8.9369837991713500e-06_wp)]
    hessg_py(3,:) = [( 1.2595237108974458e-04_wp, 1.0500955964021795e-04_wp), &
                     ( 1.0719350671889904e-05_wp, 8.9369837991713500e-06_wp), &
                     ( 4.6199031823382011e-04_wp,-2.8769180017691780e-04_wp)]

    abs_err_g = abs(g_py - g)
    abs_err_gradg = maxval(abs(gradg_py - gradg))
    abs_err_hessg = maxval(abs(hessg_py - hessg))
    max_abs_err = max(abs_err_g, abs_err_gradg, abs_err_hessg)

    call check(error, max_abs_err < 1.0e-8_wp)

  end subroutine test_gdeep_4

end module gdeep_suite
