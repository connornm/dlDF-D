%Chk=BLIND00828
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.60385461       -0.81006036       -0.14237535
  6         0.14738561       -0.86171165        1.06543132
 16        -1.16658647        0.66083029       -0.92233059
 16         0.66989468        0.56709367        1.97101840
  6         0.14583550        1.77410726       -0.49538206
  6         0.87136249        1.73514213        0.65341573
  6        -0.88444950       -2.09226569       -0.56959801
  7         0.43631218       -2.04657362        1.54551838
 16        -0.22154839       -3.22121514        0.56656299
  6         0.36141244        2.78667096       -1.47483649
  7         0.48363922        3.60618748       -2.29009492
  6         1.87708909        2.71021668        0.92115734
  7         2.68498912        3.50004751        1.19390432
  6        -1.61406058       -2.49790333       -1.71340490
  7        -2.21800980       -2.84879494       -2.64249926
  6         0.27868181       -0.66334586        5.60374433
  6        -1.12347621       -0.42335444        5.65047407
 16         1.52022023        0.50274888        6.03655151
 16        -1.85261647        1.09731570        6.18945937
  6         0.68656096        1.39971587        7.31896506
  6        -0.65099191        1.63494998        7.37607264
  6         0.51688332       -1.92586338        5.09915265
  7        -1.91357109       -1.38305199        5.23503368
 16        -0.99914831       -2.67195647        4.71208710
  6         1.57350497        1.91683658        8.30740128
  7         2.33811157        2.33787176        9.07505950
  6        -1.21653065        2.41131392        8.43028652
  7        -1.71642644        3.06198780        9.25348384
  6         1.76170630       -2.56208346        4.87374587
  7         2.77462443       -3.09597963        4.67362679

@blind-aug.gbs/N

