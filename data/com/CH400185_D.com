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
  6         0.00000000        0.00000000        3.21740611
  1         0.00334602        0.68652830        2.34869807
  1         0.56237159       -0.91936099        2.96343857
  1         0.47817888        0.49785879        4.08312404
  1        -1.04389650       -0.26502609        3.47436375

@aug-cc-pVTZ.gbs/N

