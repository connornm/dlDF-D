%Chk=BLIND00594
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.81386198       -0.61293841       -0.05517330
  6         0.24792333       -1.14122384        0.73181662
 16        -1.17741216        1.09555842       -0.24852920
 16         1.40765037       -0.14839203        1.62842497
  6         0.45337477        1.78487132       -0.15104178
  6         1.47599628        1.28997722        0.59526186
  6        -1.55776809       -1.64768023       -0.58547161
  7         0.34347771       -2.44558466        0.81767822
 16        -0.89337501       -3.15379080       -0.04229639
  6         0.60306837        2.97803481       -0.91592240
  7         0.66957363        3.95957959       -1.53506692
  6         2.73945747        1.94994080        0.64060313
  7         3.77511812        2.46983895        0.73026398
  6        -2.70377261       -1.57202674       -1.41379438
  7        -3.65070921       -1.52867843       -2.08649635
  6         0.77914120        0.65787195        5.73865368
  6         0.78741318       -0.17727443        4.58613958
 16        -0.66177097        1.19271759        6.59084434
 16        -0.66950851       -0.84630685        3.83489852
  6        -1.71967367       -0.20907626        6.34602203
  6        -1.72123437       -1.01493935        5.25138680
  6         2.06960144        1.04201533        6.04270870
  7         1.94891205       -0.43340040        4.03552186
 16         3.14989391        0.36141318        4.87028580
  6        -2.64097591       -0.40815712        7.41497417
  7        -3.38821168       -0.51558092        8.29888132
  6        -2.64886690       -2.09164975        5.13249795
  7        -3.40485520       -2.96124844        4.98016302
  6         2.51303162        1.87684863        7.09702908
  7         2.89381478        2.56743477        7.95095582

@blind-aug.gbs/N
