%Chk=CH400145
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.42394953        0.00000000        1.02286497
  1         0.06160239        1.02054675       -0.42505932
  1         0.57369191       -0.70050275       -0.63730640
  1        -1.05924382       -0.32004401        0.03950075
  6-Bq        0.00000000       0.00000000       3.37091037
  1-Bq        0.31730743      -1.03029835       3.11834894
  1-Bq        0.19285468       0.19298468       4.44401332
  1-Bq        0.57155016       0.72307685       2.75733377
  1-Bq       -1.08171226       0.11423682       3.16394544

@aug-cc-pVTZ.gbs/N

