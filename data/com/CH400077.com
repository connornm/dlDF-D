%Chk=CH400077
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.77973168        0.00000000       -0.78613260
  1         0.34526279        0.60269015        0.86229005
  1        -0.93307366        0.43682659       -0.40563781
  1        -0.19192081       -1.03951674        0.32948036
  0        -0.51595989        0.00000000        0.52019547
  0        -0.22846546       -0.39880892       -0.57058998
  0         0.61742853       -0.28905456        0.26841649
  0         0.12699682        0.68786348       -0.21802199
  6         0.00000000        0.00000000        3.24361619
  1        -0.70605755        0.84546810        3.13113342
  1        -0.18912447       -0.51321311        4.20633639
  1        -0.14388477       -0.71441905        2.41001643
  1         1.03906679        0.38216405        3.22697850
  0         0.46720864       -0.55945865        3.31804768
  0         0.12514643        0.33960064        2.60657010
  0         0.09521066        0.47274156        3.79522137
  0        -0.68756573       -0.25288356        3.25462559

@aug-cc-pVTZ.gbs/N
