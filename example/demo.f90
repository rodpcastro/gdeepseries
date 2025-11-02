program gds_demo
! Evaluation of series expressions for F.

  use gds_kinds, only: wp
  use gds

  print '(a, 3(es22.15, 1x))', 'Region 1 F(1.0, 3.0)  =', fsem(1.0_wp, 3.0_wp)
  print '(a, 3(es22.15, 1x))', 'Region 2 F(1.0, 1.0)  =', fsem(1.0_wp, 1.0_wp)
  print '(a, 3(es22.15, 1x))', 'Region 3 F(7.0, 1.0)  =', fsem(7.0_wp, 1.0_wp)
  print '(a, 3(es22.15, 1x))', 'Region 4 F(10.0, 6.0) =', fsem(10.0_wp, 6.0_wp)

end program gds_demo
