%Chk=CH400108
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.83330322       0.00000000      -0.72910336
  1-Bq       -0.54719348       0.96039249      -0.06489628
  1-Bq       -0.69073342      -0.83453884      -0.22895026
  1-Bq        0.40462368      -0.12585365       1.02294990
  6         0.00000000        0.00000000        4.05120144
  1        -0.91730177       -0.55613207        3.77686490
  1        -0.24905301        0.77577630        4.80095452
  1         0.74284309       -0.70118258        4.47841694
  1         0.42351168        0.48153835        3.14856938

@aug-cc-pVTZ.gbs/N

