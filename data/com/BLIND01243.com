%Chk=BLIND01243
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.65812040       -0.33917956       -0.70209948
  6         0.37160165       -1.28818226        0.31923668
 16         0.65597819        1.40564551       -0.49210537
 16        -0.06992989       -0.85728886        1.97874537
  6        -0.61882349        1.58163704        0.72774287
  6        -0.90335394        0.68304878        1.70714351
  6         1.00758692       -0.99841946       -1.86338287
  7         0.47649676       -2.55740481        0.00950912
 16         0.97329419       -2.70688542       -1.57218767
  6        -1.32531548        2.81506366        0.62470161
  7        -1.85618307        3.84340104        0.51563807
  6        -1.92191074        0.94377222        2.67073432
  7        -2.72220017        1.12712323        3.49344045
  6         1.37767144       -0.44300730       -3.11227167
  7         1.69260593       -0.00443992       -4.14162624
  6        -1.00039945       -0.19028759        7.71537392
  6        -0.67234961        1.19024021        7.60423003
 16         0.02594875       -1.51950371        7.19738695
 16         0.85543988        1.79746067        6.94716769
  6         1.62906770       -0.82922979        7.50988276
  6         1.95518117        0.48645067        7.40849747
  6        -2.27592889       -0.32460433        8.22555710
  7        -1.57838250        2.05771404        7.98430054
 16        -2.94511693        1.25260648        8.48903384
  6         2.60489500       -1.81139085        7.84783302
  7         3.37001749       -2.65066375        8.09562649
  6         3.29046754        0.93254839        7.63642082
  7         4.37016757        1.33660643        7.78444473
  6        -2.99787549       -1.51465391        8.48583300
  7        -3.60849310       -2.48094103        8.69685730

@aug-cc-pVTZ.gbs/N

