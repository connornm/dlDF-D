%Chk=CH400316
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.18050368        0.00000000       -1.09243049
  1        -0.95307044        0.52031109        0.21660808
  1        -0.05829136       -1.04391651        0.36445356
  1         0.83085812        0.52360542        0.51136886
  6-Bq        0.00000000       0.00000000       4.32658920
  1-Bq       -0.07361656       0.18787973       3.23788917
  1-Bq        0.93799974      -0.54661546       4.54419892
  1-Bq        0.00128615       0.96556600       4.86849934
  1-Bq       -0.86566932      -0.60683028       4.65576939

@aug-cc-pVTZ.gbs/N

