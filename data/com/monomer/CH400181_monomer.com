%Chk=CH400181
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.58083538        0.00000000        0.94266443
  1         0.16385830       -0.95575413       -0.53448149
  1        -1.07702560        0.11425019        0.23010584
  1         0.33233193        0.84150394       -0.63828878
  6-Bq        0.00000000       0.00000000       3.66527236
  1-Bq        0.24765233       1.06235623       3.47539472
  1-Bq        0.37879365      -0.29463065       4.66311709
  1-Bq       -1.09859372      -0.13432454       3.63310525
  1-Bq        0.47214774      -0.63340104       2.88947239

@aug-cc-pVTZ.gbs/N

