%Chk=BLIND01548
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.15792123        0.76251806       -0.65934513
  6         0.13599167       -0.48222389       -1.28389243
 16        -0.98229333        0.95244023        0.88106288
 16        -0.25819439       -2.06058133       -0.58532714
  6        -0.42010345       -0.50015209        1.72843540
  6        -0.13506177       -1.69393343        1.14424245
  6         0.23864083        1.79880946       -1.48040657
  7         0.70408156       -0.44511699       -2.46446728
 16         0.90376080        1.13804571       -2.93843031
  6        -0.32863469       -0.32438178        3.13975673
  7        -0.29275839       -0.14019114        4.28691833
  6         0.26610620       -2.81755983        1.92556818
  7         0.56782889       -3.76654839        2.52499312
  6         0.12402716        3.19211674       -1.25518005
  7         0.02960857        4.33896479       -1.09110510
  6         0.19597407       -0.72855303        9.39989964
  6        -1.18967446       -0.43811225        9.25344506
 16         1.39967316       -0.58826463        8.12712931
 16        -1.92756514        0.14929342        7.75509254
  6         0.74098883        0.77904942        7.21022729
  6        -0.57906910        1.06810063        7.06346040
  6         0.43829214       -1.20516719       10.67246560
  7        -1.96366886       -0.65410739       10.28888641
 16        -1.05934232       -1.27110360       11.54285056
  6         1.74704088        1.56182898        6.57288508
  7         2.60117404        2.15494974        6.05340714
  6        -1.00612637        2.16827556        6.26272378
  7        -1.40024699        3.03582495        5.59717654
  6         1.66872808       -1.62269291       11.23537598
  7         2.66822276       -1.97803682       11.71059911

@aug-cc-pVTZ.gbs/N

