%Chk=CH400117
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.06665795        0.00000000       -0.29702992
  1        -0.16955213       -0.78039849        0.76695255
  1        -0.26725011        0.99066375        0.41611138
  1        -0.62985571       -0.21026525       -0.88603401
  0        -0.70582321        0.00000000        0.19654905
  0         0.11219513        0.51640113       -0.50750376
  0         0.17684332       -0.65553673       -0.27534700
  0         0.41678476        0.13913560        0.58630170
  6         0.00000000        0.00000000        2.94458188
  1         0.76680706       -0.67472280        3.37206526
  1        -0.91072184       -0.03083602        3.57356584
  1         0.39358208        1.03451968        2.91545916
  1        -0.24966729       -0.32896086        1.91723727
  0        -0.50740748        0.44647398        2.66170986
  0         0.60263800        0.02040465        2.52837397
  0        -0.26043903       -0.68455685        2.96385281
  0         0.16520851        0.21767823        3.62439088

@aug-cc-pVTZ.gbs/N

