%Chk=BLIND00039
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.53428939       -0.84682238       -0.19630488
  6        -1.33455042        0.24622526        0.24033660
 16         0.96243580       -0.70554384       -1.10645362
 16        -0.92920273        1.94687519       -0.03907603
  6         1.59404431        0.81713035       -0.45326446
  6         0.84129944        1.86765383       -0.03205235
  6        -1.14872942       -2.03586207        0.14126186
  7        -2.45572301       -0.03617143        0.85765673
 16        -2.65117805       -1.68756558        0.93256726
  6         3.01803068        0.87445320       -0.44510572
  7         4.17986194        0.88365977       -0.48055109
  6         1.45170102        3.06947796        0.43374868
  7         1.90540311        4.07741573        0.79311201
  6        -0.70625337       -3.35695573       -0.11190789
  7        -0.36303438       -4.44794479       -0.31936753
  6        -0.41896236        0.86794934        3.97881503
  6        -1.36807753       -0.15604463        3.70233102
 16         1.22956771        0.59938581        4.52538515
 16        -1.02639502       -1.88812171        3.83556166
  6         1.59609179       -0.92894695        3.70442160
  6         0.69948984       -1.91401033        3.43349886
  6        -1.01144140        2.10468271        3.82145442
  7        -2.57673571        0.21848392        3.36066775
 16        -2.67214171        1.88026299        3.37812965
  6         2.97547392       -1.07037757        3.37519223
  7         4.11280099       -1.15009989        3.14864308
  6         1.10866116       -3.12906559        2.80886979
  7         1.40214181       -4.14597205        2.32856105
  6        -0.43572342        3.38517907        4.00597091
  7         0.01777773        4.44402894        4.16175784

@aug-cc-pVTZ.gbs/N

