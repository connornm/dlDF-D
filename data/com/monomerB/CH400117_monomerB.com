%Chk=CH400117
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.06665795       0.00000000      -0.29702992
  1-Bq       -0.16955213      -0.78039849       0.76695255
  1-Bq       -0.26725011       0.99066375       0.41611138
  1-Bq       -0.62985571      -0.21026525      -0.88603401
  6         0.00000000        0.00000000        2.94458188
  1         0.76680706       -0.67472280        3.37206526
  1        -0.91072184       -0.03083602        3.57356584
  1         0.39358208        1.03451968        2.91545916
  1        -0.24966729       -0.32896086        1.91723727

@aug-cc-pVTZ.gbs/N
