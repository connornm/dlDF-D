%Chk=CH400084
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.05837631        0.00000000        0.32530839
  1        -0.65741308       -0.12145012        0.88263461
  1        -0.23138322        0.95864564       -0.50343450
  1        -0.16958000       -0.83719552       -0.70450851
  6-Bq        0.00000000       0.00000000       4.23749463
  1-Bq        0.96053753       0.35779235       3.81876055
  1-Bq        0.18520318      -0.87417841       4.89132782
  1-Bq       -0.46953964       0.81006140       4.82852679
  1-Bq       -0.67620108      -0.29367534       3.41136337

@aug-cc-pVTZ.gbs/N

