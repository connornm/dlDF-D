%Chk=CH400287
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.09729676        0.00000000       -0.14807352
  1        -0.23255074        0.31222003        1.03654463
  1        -0.46853282        0.70656777       -0.71219727
  1        -0.39621321       -1.01878780       -0.17627384
  6-Bq        0.00000000       0.00000000       5.01328434
  1-Bq       -0.99734118      -0.17459829       5.46140481
  1-Bq        0.26475508       1.07095669       5.10784901
  1-Bq        0.75457263      -0.61474144       5.54120391
  1-Bq       -0.02198652      -0.28161696       3.94267964

@aug-cc-pVTZ.gbs/N

