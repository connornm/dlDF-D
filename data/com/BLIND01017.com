%Chk=BLIND01017
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.81985401        0.05871245       -0.60456692
  6         0.56680226       -1.25147265       -0.10938366
 16         0.42985017        1.55162744        0.23663370
 16        -0.22062213       -1.59215126        1.43932202
  6        -1.03770962        1.05984999        1.10181880
  6        -1.29212571       -0.18735210        1.57868719
  6         1.48238557       -0.02172983       -1.81282359
  7         0.97616377       -2.26095631       -0.83823611
 16         1.74711177       -1.69076615       -2.19891214
  6        -1.95237996        2.13331885        1.30708064
  7        -2.65500734        3.04472692        1.47091937
  6        -2.48718656       -0.46707571        2.30499718
  7        -3.43263473       -0.73310273        2.92651429
  6         1.93335833        1.03523125       -2.64024446
  7         2.31850811        1.89019846       -3.32711172
  6         0.70751439        0.71394139        7.43929611
  6         0.58217801        0.35795081        8.81166201
 16         0.41142996       -0.36358877        6.08291143
 16         0.07038880       -1.23882412        9.37996647
  6        -0.87780547       -1.37930098        6.75397239
  6        -1.00955787       -1.72569195        8.06175626
  6         1.16573090        2.01193021        7.33650574
  7         0.89980629        1.26627270        9.70156293
 16         1.41855280        2.64679629        8.92956968
  6        -1.78142440       -1.86319174        5.76381472
  7        -2.47741309       -2.25211922        4.91792573
  6        -2.05982047       -2.58967676        8.49119854
  7        -2.88156632       -3.31336350        8.88099001
  6         1.44019835        2.75709636        6.16408796
  7         1.67941385        3.38027892        5.21250225

@aug-cc-pVTZ.gbs/N

