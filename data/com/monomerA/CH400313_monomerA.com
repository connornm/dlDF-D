%Chk=CH400313
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.04161889        0.00000000        0.37552104
  1        -0.54821601        0.85935280        0.43238632
  1         0.00570106        0.08360573       -1.10406682
  1        -0.49910394       -0.94295853        0.29615945
  6-Bq        0.00000000       0.00000000       3.28312395
  1-Bq        0.10129457      -0.02473323       4.38544588
  1-Bq       -1.05953566       0.17543683       3.01369480
  1-Bq        0.62516600       0.81680916       2.87327492
  1-Bq        0.33307509      -0.96751275       2.86008020

@aug-cc-pVTZ.gbs/N

