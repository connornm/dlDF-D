%Chk=BLIND00123
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.67928174       -0.76099430       -0.02393538
  6         0.36554748       -0.95630333        0.92264422
 16        -1.31975907        0.79949495       -0.51706695
 16         1.22080007        0.35946801        1.74232534
  6         0.13955835        1.79941481       -0.39592642
  6         1.14465356        1.62316245        0.50205062
  6        -1.16145432       -1.98490408       -0.44209525
  7         0.68352570       -2.19176610        1.22327958
 16        -0.30229410       -3.24107058        0.38779390
  6         0.15700886        2.87876861       -1.32639970
  7         0.11488091        3.75710044       -2.08663088
  6         2.25603133        2.51582162        0.54610706
  7         3.16495553        3.23448766        0.63722273
  6        -2.20273063       -2.25117735       -1.36405547
  7        -3.06100002       -2.48814319       -2.11126143
  6         0.46064894        0.36763859        6.45352593
  6        -0.43673637       -0.70402242        6.18506937
 16         0.28687026        1.53033138        7.75984254
 16        -1.87356269       -1.06203228        7.15560911
  6        -0.47907823        0.51472236        8.99514951
  6        -1.33644430       -0.51174356        8.75238380
  6         1.44845245        0.39826871        5.48988219
  7        -0.19276120       -1.44525190        5.13201904
 16         1.15712519       -0.87071438        4.34548850
  6        -0.14583048        0.90411051       10.32496622
  7         0.13615159        1.27374746       11.39035237
  6        -1.93616463       -1.23765540        9.82351847
  7        -2.47009961       -1.83682939       10.66425035
  6         2.53091440        1.30315765        5.36896660
  7         3.41977174        2.04249919        5.24895564

@aug-cc-pVTZ.gbs/N

