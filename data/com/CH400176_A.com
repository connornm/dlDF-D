%Chk=CH400176
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10384150        0.00000000        0.08671731
  1        -0.40690120       -0.91781108        0.46694774
  1        -0.41072124        0.88965035        0.51557374
  1        -0.28621906        0.02816073       -1.06923878
  6-Bq        0.00000000       0.00000000       2.65426365
  1-Bq       -0.09693448       1.10131928       2.59355528
  1-Bq        0.93576223      -0.31628141       2.15396266
  1-Bq        0.02655962      -0.31015503       3.71684748
  1-Bq       -0.86538736      -0.47488284       2.15268917

@aug-cc-pVTZ.gbs/N
