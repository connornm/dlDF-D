%Chk=CH400138
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.06386983        0.00000000       -1.10539884
  1        -1.03322737       -0.24965823        0.30999664
  1         0.70052913       -0.75299606        0.41017291
  1         0.26882841        1.00265430        0.38522929
  6-Bq        0.00000000       0.00000000       2.47630960
  1-Bq       -0.79434504      -0.22626915       1.73887880
  1-Bq       -0.45504716       0.43506214       3.38715480
  1-Bq        0.53495887      -0.93276051       2.74043607
  1-Bq        0.71443333       0.72396752       2.03876875

@aug-cc-pVTZ.gbs/N
