%Chk=BLIND00088
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
  6         0.16281999        0.10582151        5.30393187
  6        -0.97178318       -0.69418508        5.61776586
 16         0.76518227        1.43605054        6.28176981
 16        -1.95243475       -0.49598475        7.07846003
  6         0.35049530        0.85591338        7.90522220
  6        -0.72867026        0.09097808        8.21814350
  6         0.67521337       -0.27523421        4.08018648
  7        -1.32304432       -1.60976554        4.74817512
 16        -0.30234395       -1.55282818        3.43446062
  6         1.25302417        1.31160194        8.90965200
  7         2.00167037        1.72719869        9.69578255
  6        -0.99795743       -0.28338162        9.56769200
  7        -1.27210950       -0.58725016       10.65554801
  6         1.78747334        0.26127088        3.38721570
  7         2.69188261        0.69549424        2.80019793

@blind-aug.gbs/N

