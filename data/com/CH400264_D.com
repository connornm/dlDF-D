%Chk=CH400264
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.07941786        0.00000000        0.24666381
  1        -0.12767193        0.06292292       -1.09805580
  1        -0.46373343       -0.93387739        0.37257253
  1        -0.48801251        0.87095446        0.47881946
  6         0.00000000        0.00000000        3.84278728
  1         0.02119438        0.20820794        4.93007115
  1        -1.04325990        0.05728001        3.47629891
  1         0.61942681        0.74838441        3.31154613
  1         0.40263872       -1.01387236        3.65323293

@aug-cc-pVTZ.gbs/N
