%Chk=BLIND02464
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.39204492        0.72789752       -0.59797566
  6         0.36169138       -0.55813473       -1.20709500
 16         0.57196450        1.02826079        1.12421117
 16         0.46961600       -2.08248070       -0.31305362
  6        -0.26231545       -0.38897986        1.78718024
  6        -0.29994324       -1.62160853        1.21546811
  6         0.33067765        1.70423268       -1.57178618
  7         0.28195917       -0.60424827       -2.51452620
 16         0.26916244        0.94314024       -3.12799335
  6        -0.88137785       -0.13582107        3.04565631
  7        -1.34913943        0.11227088        4.08047152
  6        -0.96166726       -2.70978020        1.85729805
  7        -1.46334292       -3.63122297        2.35748890
  6         0.34534801        3.11076852       -1.40865998
  7         0.36555375        4.26729460       -1.29405167
  6         0.49675965       -0.61101442       10.12446163
  6         0.40668352       -1.13209909        8.80302918
 16         0.56655655        1.09589202       10.53691218
 16         0.31712570       -0.13105510        7.34547003
  6        -0.43213334        1.78122571        9.24189931
  6        -0.52789053        1.29310916        7.97690254
  6         0.59304635       -1.65038519       11.02764516
  7         0.42536242       -2.43549095        8.66635154
 16         0.58951669       -3.15119833       10.16025020
  6        -1.12395560        2.96354343        9.63508686
  7        -1.64848368        3.93648242        9.99488266
  6        -1.32604557        1.94957321        6.99407731
  7        -1.94026219        2.46730847        6.15406882
  6         0.71230302       -1.58217393       12.43700441
  7         0.81997098       -1.54481595       13.59377056

@aug-cc-pVTZ.gbs/N
