%Chk=CH400253
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.03075364        0.00000000       -0.40439200
  1        -0.45236282        1.00052831       -0.14246730
  1        -0.60565614       -0.75820076       -0.53319622
  1         0.02726532       -0.24232755        1.08005552
  0        -0.68206480        0.00000000        0.26759212
  0         0.29933511       -0.66206426        0.09427270
  0         0.40077155        0.50171257        0.35282376
  0        -0.01804186        0.16035170       -0.71468858
  6         0.00000000        0.00000000        4.07162924
  1        -0.15875036       -0.58114681        3.14262296
  1         1.07149832        0.26347299        4.16360359
  1        -0.60601219        0.92598748        4.03581890
  1        -0.30673577       -0.60831365        4.94447150
  0         0.10504744        0.38455337        4.68636632
  0        -0.70902616       -0.17434394        4.01076846
  0         0.40100715       -0.61273950        4.09532547
  0         0.20297156        0.40253007        3.49405671

@aug-cc-pVTZ.gbs/N

