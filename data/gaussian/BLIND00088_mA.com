%Chk=BLIND00088_mA
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.22481105       -0.43959529       -0.89293050
  6        -0.56877239       -1.25478848       -0.03764574
 16         0.18804649        1.31687384       -0.93744130
 16        -1.71205538       -0.62147563        1.15667715
  6        -0.16715847        1.67473620        0.76252518
  6        -0.92212878        0.90471084        1.59000620
  6         0.97491751       -1.23927711       -1.73147239
  7        -0.46277831       -2.55399904       -0.17431039
 16         0.60821760       -2.90085024       -1.40051546
  6         0.38867569        2.91129720        1.20190534
  7         0.84135882        3.93762567        1.50673462
  6        -1.18688558        1.30861361        2.93191094
  7        -1.45232410        1.61215214        4.02201776
  6         1.88632893       -0.84380140       -2.74048013
  7         2.63080829       -0.53725233       -3.57885598

@blind-aug.gbs/N

