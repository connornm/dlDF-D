%Chk=CH400194
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.82993740        0.00000000       -0.73293238
  1        -0.08849312       -1.00447569        0.45736584
  1         0.20510380        0.74839920        0.78982088
  1        -0.94654808        0.25607649       -0.51425434
  0        -0.54918175        0.00000000        0.48499211
  0         0.05855720        0.66467630       -0.30264569
  0        -0.13572019       -0.49522673       -0.52263606
  0         0.62634475       -0.16944957        0.34028964
  6         0.00000000        0.00000000        3.16148139
  1         0.09207683       -1.02136720        2.74396869
  1        -1.01463404        0.13620203        3.58332582
  1         0.75379498        0.14162802        3.96005549
  1         0.16876223        0.74353715        2.35857556
  0        -0.06092859        0.67585366        3.43775567
  0         0.67139823       -0.09012688        2.88234074
  0        -0.49879719       -0.09371734        2.63305320
  0        -0.11167244       -0.49200944        3.69277596

@aug-cc-pVTZ.gbs/N
