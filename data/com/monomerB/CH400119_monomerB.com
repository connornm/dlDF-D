%Chk=CH400119
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.08464944       0.00000000       0.22253438
  1-Bq       -0.51145054       0.73039613       0.65644938
  1-Bq       -0.15947079       0.28072382      -1.05912660
  1-Bq       -0.41372811      -1.01111996       0.18014284
  6         0.00000000        0.00000000        2.93901632
  1        -0.27814730       -1.06565599        3.05302223
  1        -0.91059573        0.62730042        2.99642317
  1         0.69815598        0.28788598        3.74876063
  1         0.49058706        0.15046959        1.95785925

@aug-cc-pVTZ.gbs/N
