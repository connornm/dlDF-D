%Chk=CH400300
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.05967014       0.00000000       0.32106875
  1-Bq       -0.63742336       0.35941140       0.83096382
  1-Bq       -0.12086704       0.66908259      -0.87390252
  1-Bq       -0.30137974      -1.02849399      -0.27813005
  6         0.00000000        0.00000000        4.01291651
  1         0.52193904        0.96141517        3.84190054
  1        -0.07569048       -0.55397190        3.05721151
  1         0.56963345       -0.60399588        4.74550962
  1        -1.01588201        0.19655261        4.40704436

@aug-cc-pVTZ.gbs/N

