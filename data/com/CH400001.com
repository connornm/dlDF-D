%Chk=CH400001
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.08913654        0.00000000       -0.19941802
  1        -0.17956882       -0.22793350        1.06854448
  1        -0.41923216        0.99621315       -0.24039491
  1        -0.49033557       -0.76827965       -0.62873156
  0        -0.72069763        0.00000000        0.13195783
  0         0.11882332        0.15082694       -0.70707156
  0         0.27741207       -0.65920886        0.15907284
  0         0.32446224        0.50838191        0.41604089
  6         0.00000000        0.00000000        3.43471632
  1         0.35476448        0.42502071        2.47581825
  1         0.35594112       -1.04415321        3.52977388
  1         0.39616061        0.60289787        4.27469066
  1        -1.10686621        0.01623464        3.45858250
  0        -0.23475286       -0.28124244        4.06923324
  0        -0.23553146        0.69093150        3.37181534
  0        -0.26214528       -0.39894636        2.87889298
  0         0.73242960       -0.01074270        3.41892373

@aug-cc-pVTZ.gbs/N
