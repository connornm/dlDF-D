%Chk=CH400070
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.89772269        0.00000000        0.64813573
  1         0.16305113        0.68267859       -0.85635871
  1        -0.18431198       -1.02528844       -0.37523149
  1        -0.87646183        0.34260985        0.58345446
  0        -0.59403627        0.00000000       -0.42888092
  0        -0.10789332       -0.45173843        0.56666512
  0         0.12196194        0.67844840        0.24829618
  0         0.57996765       -0.22670996       -0.38608037
  6         0.00000000        0.00000000        2.94947590
  1        -0.20820881       -0.39052599        3.96442639
  1        -0.95306163        0.27282621        2.45629707
  1         0.64657478        0.89574195        3.02413086
  1         0.51469566       -0.77804217        2.35304927
  0         0.13777483        0.25841677        2.27786827
  0         0.63065486       -0.18053310        3.27581956
  0        -0.42784802       -0.59272559        2.90007561
  0        -0.34058167        0.51484191        3.34414015

@aug-cc-pVTZ.gbs/N
