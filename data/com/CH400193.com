%Chk=CH400193
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10575427        0.00000000       -0.05738864
  1        -0.37582720        1.03452391       -0.12041661
  1        -0.41139949       -0.63827512       -0.80581715
  1        -0.31852758       -0.39624879        0.98362241
  0        -0.73169382        0.00000000        0.03797491
  0         0.24869037       -0.68455964        0.07968144
  0         0.27222908        0.42235601        0.53322103
  0         0.21077437        0.26220363       -0.65087737
  6         0.00000000        0.00000000        4.47644019
  1         0.37650129        0.68166800        5.26356246
  1         0.45378614        0.27275457        3.50398503
  1         0.27125423       -1.04294992        4.73072856
  1        -1.10154165        0.08852734        4.40748472
  0        -0.24913642       -0.45106972        3.95558984
  0        -0.30027694       -0.18048570        5.11992804
  0        -0.17949290        0.69013526        4.30817385
  0         0.72890627       -0.05857984        4.52206904

@aug-cc-pVTZ.gbs/N
