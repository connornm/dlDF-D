%Chk=CH400137
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.59316114        0.00000000        0.93495766
  1        -0.92389108        0.59177334        0.14904878
  1         0.59811313        0.44887987       -0.81654976
  1        -0.26738318       -1.04065321       -0.26745668
  0        -0.39250343        0.00000000       -0.61867520
  0         0.61135228       -0.39158510       -0.09862776
  0        -0.39578023       -0.29703039        0.54032295
  0         0.17693138        0.68861550        0.17698001
  6         0.00000000        0.00000000        3.24319566
  1         0.94757601       -0.57222130        3.21773100
  1         0.03442190        0.80496754        2.48370350
  1        -0.84695038       -0.67857975        3.02367223
  1        -0.13504754        0.44583351        4.24767593
  0        -0.62702495        0.37864723        3.26004600
  0        -0.02277748       -0.53265883        3.74576277
  0         0.56043949        0.44902617        3.38845754
  0         0.08936294       -0.29501457        2.57851634

@aug-cc-pVTZ.gbs/N

