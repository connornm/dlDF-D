%Chk=CH400125
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.09673844        0.00000000        0.15215372
  1        -0.28683811       -0.87259579       -0.61829317
  1        -0.50879468       -0.05994439        0.98159085
  1        -0.30110565        0.93254018       -0.51545140
  6-Bq        0.00000000       0.00000000       2.82353078
  1-Bq       -0.02774518       0.69471875       1.96179879
  1-Bq       -0.07182269      -1.04309621       2.45914396
  1-Bq        0.95082854       0.13368769       3.37492777
  1-Bq       -0.85126067       0.21468977       3.49825260

@aug-cc-pVTZ.gbs/N

