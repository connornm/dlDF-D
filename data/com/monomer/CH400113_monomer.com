%Chk=CH400113
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.98678580        0.00000000       -0.50223474
  1        -0.77971321       -0.31952960       -0.71828553
  1         0.02198162       -0.70090331        0.85687648
  1        -0.22905421        1.02043291        0.36364378
  6-Bq        0.00000000       0.00000000       2.45637412
  1-Bq       -0.45366350       0.38794621       3.38893679
  1-Bq        1.05864370      -0.26508230       2.64343254
  1-Bq       -0.55260045      -0.89998108       2.12372904
  1-Bq       -0.05237975       0.77711717       1.66939813

@aug-cc-pVTZ.gbs/N

