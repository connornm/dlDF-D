%Chk=BLIND02349
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.06094865       -0.39892967       -0.93714955
  6         0.53315223       -1.27075926        0.01824773
 16         0.02860720        1.35580307       -0.89960772
 16         1.44256962       -0.71986238        1.43389633
  6         0.06095642        1.63924707        0.85045945
  6         0.62364800        0.81482908        1.77314965
  6        -0.65693490       -1.14205929       -1.93607997
  7         0.42204781       -2.55986562       -0.19084704
 16        -0.40300797       -2.82588373       -1.61191502
  6        -0.53685001        2.87667478        1.22793767
  7        -1.01307845        3.90562495        1.48418180
  6         0.63695409        1.16035878        3.15675873
  7         0.69654791        1.41483581        4.28928518
  6        -1.34797225       -0.67989985       -3.08239777
  7        -1.91015005       -0.31834197       -4.03335935
  6        -0.71644435        0.49839370        8.75259795
  6        -0.32104734       -0.73356104        8.15949140
 16        -0.71985938        0.83651029       10.47717968
 16         0.26882232       -2.13066038        9.07314099
  6         0.66168441       -0.14417302       11.00036137
  6         1.05041152       -1.32110094       10.44235066
  6        -1.16432882        1.36219330        7.77354077
  7        -0.43624252       -0.84062800        6.85828936
 16        -1.07787300        0.56415552        6.23715291
  6         1.33729174        0.39941153       12.13136103
  7         1.84049042        0.86905388       13.06806812
  6         2.15146465       -2.05724132       10.97123749
  7         3.02307633       -2.70436368       11.38651074
  6        -1.65866356        2.68064801        7.92347011
  7        -2.07568255        3.76049617        8.02891689

@blind-aug.gbs/N

