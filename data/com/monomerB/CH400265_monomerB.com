%Chk=CH400265
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.07030597       0.00000000       0.28360374
  1-Bq       -0.08953237      -0.03471526      -1.10307061
  1-Bq       -0.49808732      -0.88620206       0.43879481
  1-Bq       -0.48268628       0.92091732       0.38067205
  6         0.00000000        0.00000000        3.55536252
  1         0.00299122        0.04614493        4.66163900
  1        -1.01977611        0.21223456        3.17986402
  1         0.70558513        0.75231326        3.15267115
  1         0.31119976       -1.01069275        3.22727591

@aug-cc-pVTZ.gbs/N
