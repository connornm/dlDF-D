%Chk=CH400293
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.60578875        0.00000000       -0.92682574
  1         0.33225652       -0.82613724        0.65809485
  1         0.12985510        0.96574114        0.52580191
  1        -1.06790037       -0.13960391       -0.25707102
  0        -0.40085930        0.00000000        0.61329419
  0        -0.21985901        0.54666713       -0.43547101
  0        -0.08592703       -0.63904508       -0.34793084
  0         0.70664534        0.09237795        0.17010766
  6         0.00000000        0.00000000        3.31421574
  1         0.30580904        0.46218498        4.27278334
  1        -0.88754617        0.52887039        2.91604007
  1         0.83132278        0.07466252        2.58667613
  1        -0.24958565       -1.06571788        3.48136343
  0        -0.20235833       -0.30583458        2.67991751
  0         0.58730232       -0.34996129        3.57769443
  0        -0.55009848       -0.04940528        3.79563938
  0         0.16515449        0.70520116        3.20361167

@aug-cc-pVTZ.gbs/N
