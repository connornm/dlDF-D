%Chk=CH400133
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.55289486        0.00000000       -0.95931915
  1        -0.65547242       -0.89107636        0.04821613
  1         0.71988808       -0.02543015        0.84089263
  1        -0.61731052        0.91650651        0.07021039
  0        -0.36585864        0.00000000        0.63479555
  0         0.43373571        0.58963829       -0.03190532
  0        -0.47636050        0.01682750       -0.55643099
  0         0.40848343       -0.60646580       -0.04645925
  6         0.00000000        0.00000000        3.43916933
  1        -0.12310016       -1.09757334        3.51768729
  1        -0.99566760        0.47871513        3.36525168
  1         0.52532743        0.37769252        4.33770448
  1         0.59344032        0.24116569        2.53603388
  0         0.08145718        0.72628038        3.38721285
  0         0.65884785       -0.31677282        3.48808173
  0        -0.34761687       -0.24992468        2.84459545
  0        -0.39268817       -0.15958287        4.03678731

@aug-cc-pVTZ.gbs/N
