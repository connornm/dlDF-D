%Chk=BLIND01740
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.23252212       -0.72674340        0.67741155
  6        -0.05415473        0.56052218        1.25789329
 16        -0.83016222       -1.03061412       -0.94709008
 16        -0.37955524        2.08304178        0.41490045
  6        -0.18539945        0.38540326       -1.79755926
  6        -0.00914444        1.61919079       -1.25502211
  6         0.06701419       -1.70109873        1.60811714
  7         0.34455923        0.60929065        2.50549452
 16         0.50878181       -0.93685244        3.09996014
  6         0.10539838        0.12978956       -3.16912330
  7         0.30452832       -0.12032687       -4.28670196
  6         0.47376645        2.70615776       -2.04194130
  7         0.83647308        3.62665951       -2.65189749
  6         0.01358180       -3.10796363        1.45631268
  7        -0.03344611       -4.26472197        1.35241648
  6        -0.25941270        0.70455403       12.87640410
  6         1.09672402        0.28285518       12.97069057
 16        -1.20132784        0.76051350       11.39377770
 16         2.05540260       -0.29079085       11.59720815
  6        -0.50242937       -0.61263333       10.51620674
  6         0.78981228       -1.02623310       10.59800391
  6        -0.70009366        1.13036324       14.11325091
  7         1.67158509        0.35946827       14.14607357
 16         0.60030130        0.99339013       15.25124822
  6        -1.43181848       -1.25105175        9.64457180
  7        -2.21866500       -1.72332892        8.93111553
  6         1.26379933       -2.11700680        9.81092643
  7         1.70021664       -2.98106494        9.16759945
  6        -1.97391844        1.63875904       14.46549030
  7        -3.01069201        2.06678597       14.77043198

@aug-cc-pVTZ.gbs/N

