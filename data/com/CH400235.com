%Chk=CH400235
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.16681047        0.00000000        1.09460505
  1         0.96942653        0.12115243       -0.52107605
  1        -0.46439493       -0.95852696       -0.30257128
  1        -0.67184207        0.83737453       -0.27095772
  0        -0.11038093        0.00000000       -0.72431622
  0        -0.64148375       -0.08016834        0.34480367
  0         0.30729694        0.63427135        0.20021585
  0         0.44456775       -0.55410301        0.17929670
  6         0.00000000        0.00000000        4.86732377
  1        -0.50514574       -0.74925571        5.50719064
  1        -0.52399397        0.06770856        3.89427093
  1        -0.02177974        0.98746752        5.36774309
  1         1.05091945       -0.30592037        4.70009041
  0         0.33426235        0.49579349        4.44391447
  0         0.34673450       -0.04480375        5.51120710
  0         0.01441197       -0.65342174        4.53618896
  0        -0.69540881        0.20243200        4.97798453

@aug-cc-pVTZ.gbs/N

