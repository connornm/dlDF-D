%Chk=CH400110
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.04968352       0.00000000      -0.35234993
  1-Bq       -0.04159726      -0.38881131       1.03589642
  1-Bq       -0.39689096       1.03341894      -0.02255687
  1-Bq       -0.61119530      -0.64460763      -0.66098961
  6         0.00000000        0.00000000        3.27944288
  1         0.65947261       -0.82226168        2.94037421
  1        -0.41851151       -0.24873753        4.27390943
  1         0.58505500        0.93745959        3.34920053
  1        -0.82601609        0.13353963        2.55428735

@aug-cc-pVTZ.gbs/N

