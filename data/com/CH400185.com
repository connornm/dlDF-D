%Chk=CH400185
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.99885690        0.00000000       -0.47777698
  1         0.06424920        0.49236878        0.98966205
  1        -0.71554703        0.55100048       -0.64060665
  1        -0.34755907       -1.04336926        0.12872157
  0        -0.66095827        0.00000000        0.31615203
  0        -0.04251464       -0.32580764       -0.65487390
  0         0.47348797       -0.36460510        0.42389882
  0         0.22998494        0.69041274       -0.08517695
  6         0.00000000        0.00000000        3.21740611
  1         0.00334602        0.68652830        2.34869807
  1         0.56237159       -0.91936099        2.96343857
  1         0.47817888        0.49785879        4.08312404
  1        -1.04389650       -0.26502609        3.47436375
  0        -0.00221411       -0.45428584        3.79224296
  0        -0.37212953        0.60835466        3.38546015
  0        -0.31641798       -0.32944046        2.64454785
  0         0.69076162        0.17537165        3.04737346

@aug-cc-pVTZ.gbs/N
