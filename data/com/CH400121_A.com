%Chk=CH400121
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.06439789        0.00000000       -0.30502966
  1        -0.48165046       -0.93687712       -0.34096928
  1        -0.06785567        0.06966555        1.10296341
  1        -0.51489176        0.86721157       -0.45696446
  6-Bq        0.00000000       0.00000000       2.90830327
  1-Bq        0.81095322       0.17273266       3.64213176
  1-Bq        0.29101461      -0.81925689       2.22265423
  1-Bq       -0.92857282      -0.27866509       3.44317170
  1-Bq       -0.17339501       0.92518932       2.32525538

@aug-cc-pVTZ.gbs/N

