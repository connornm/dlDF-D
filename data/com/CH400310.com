%Chk=CH400310
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.12611048        0.00000000        1.10003732
  1         0.49496477       -0.89308711       -0.42824204
  1        -1.07894075       -0.02155957       -0.24780637
  1         0.45786550        0.91464668       -0.42398891
  0        -0.08344915        0.00000000       -0.72791083
  0        -0.32752545        0.59096884        0.28337404
  0         0.71395092        0.01426628        0.16397711
  0        -0.30297632       -0.60523512        0.28055968
  6         0.00000000        0.00000000        3.55926692
  1         0.39049114        1.03251578        3.64537032
  1        -1.10271057        0.02922764        3.46355323
  1         0.43537435       -0.48584216        2.66462032
  1         0.27684508       -0.57590126        4.46352379
  0        -0.25839371       -0.68323083        3.50229103
  0         0.72967976       -0.01934035        3.62260207
  0        -0.28809359        0.32148888        4.15126769
  0        -0.18319245        0.38108231        2.96090688

@aug-cc-pVTZ.gbs/N
