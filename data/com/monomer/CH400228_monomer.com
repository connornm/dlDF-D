%Chk=CH400228
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.03205737        0.00000000       -0.40105305
  1        -0.03091904        0.58527408        0.93940631
  1        -0.68415915        0.45597144       -0.74162137
  1        -0.31697918       -1.04124552        0.20326810
  6-Bq        0.00000000       0.00000000       4.32391448
  1-Bq        0.10257865      -0.06714824       3.22348061
  1-Bq        0.43055983       0.95703897       4.67701589
  1-Bq       -1.07174794      -0.04603265       4.59818316
  1-Bq        0.53860946      -0.84385808       4.79697825

@aug-cc-pVTZ.gbs/N
