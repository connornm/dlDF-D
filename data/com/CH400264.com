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
  0        -0.71426663        0.00000000       -0.16322106
  0         0.08448239       -0.04163702        0.72659963
  0         0.30685921        0.61796036       -0.24653671
  0         0.32292503       -0.57632334       -0.31684186
  6         0.00000000        0.00000000        3.84278728
  1         0.02119438        0.20820794        4.93007115
  1        -1.04325990        0.05728001        3.47629891
  1         0.61942681        0.74838441        3.31154613
  1         0.40263872       -1.01387236        3.65323293
  0        -0.01402463       -0.13777425        3.12331559
  0         0.69034038       -0.03790302        4.08529801
  0        -0.40988381       -0.49521694        4.19431734
  0        -0.26643195        0.67089421        3.96821817

@aug-cc-pVTZ.gbs/N

