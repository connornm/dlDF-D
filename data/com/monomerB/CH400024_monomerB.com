%Chk=CH400024
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.03371575       0.00000000      -0.39675901
  1-Bq       -0.62716705       0.68396418      -0.60402024
  1-Bq       -0.41552451      -1.02496724      -0.05260686
  1-Bq        0.00897582       0.34100307       1.05338611
  6         0.00000000        0.00000000        5.02903084
  1        -0.89136265       -0.45807404        4.55825832
  1         0.61573545        0.49273795        4.25181458
  1        -0.32208055        0.75124616        5.77594401
  1         0.59770774       -0.78591007        5.53010646

@aug-cc-pVTZ.gbs/N

