%Chk=CH400038
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.42799444        0.00000000       -1.02117908
  1         0.36169582       -0.88921390        0.55177960
  1         0.31538000        0.91820799        0.53236781
  1        -1.10507026       -0.02899409       -0.06296833
  0        -0.28321020        0.00000000        0.67572918
  0        -0.23933943        0.58840588       -0.36512066
  0        -0.20869158       -0.60759170       -0.35227559
  0         0.73124120        0.01918582        0.04166707
  6         0.00000000        0.00000000        3.30972461
  1         0.48369567        0.64938312        2.55452529
  1        -0.06095845        0.53879292        4.27511169
  1        -1.01995291       -0.26486882        2.96983268
  1         0.59721568       -0.92330721        3.43942880
  0        -0.32006852       -0.42970633        3.80945109
  0         0.04033710       -0.35652718        2.67091382
  0         0.67491780        0.17526758        3.53463610
  0        -0.39518638        0.61096593        3.22389745

@aug-cc-pVTZ.gbs/N

