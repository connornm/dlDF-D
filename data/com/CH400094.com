%Chk=CH400094
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10661522        0.00000000        0.03726553
  1        -0.36694157        1.04234175       -0.06973885
  1        -0.33945559       -0.57083699       -0.88594638
  1        -0.40021805       -0.47150476        0.91841969
  0        -0.73226352        0.00000000       -0.02465915
  0         0.24281062       -0.68973283        0.04614722
  0         0.22462274        0.37773121        0.58624371
  0         0.26483015        0.31200161       -0.60773178
  6         0.00000000        0.00000000        3.05944332
  1         0.17660220       -0.83846267        3.76071214
  1        -0.79832081       -0.27864482        2.34458186
  1        -0.31078226        0.89768556        3.62826752
  1         0.93250087        0.21942193        2.50421175
  0        -0.11686027        0.55482305        2.59540345
  0         0.52826059        0.18438336        3.53247763
  0         0.20564918       -0.59401170        2.68304400
  0        -0.61704950       -0.14519471        3.42684819

@aug-cc-pVTZ.gbs/N

