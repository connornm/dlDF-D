%Chk=CH400168
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.78679848        0.00000000        0.77905976
  1         0.15871725       -0.85543720       -0.68485180
  1        -0.99400873       -0.09044582        0.47932468
  1         0.04849300        0.94588302       -0.57353265
  6-Bq        0.00000000       0.00000000       2.36652063
  1-Bq       -0.31909155       0.83865768       1.71781435
  1-Bq       -0.89275578      -0.53779089       2.74035753
  1-Bq        0.63357028      -0.69670338       1.78412981
  1-Bq        0.57827704       0.39583658       3.22378081

@aug-cc-pVTZ.gbs/N
