%Chk=CH400061
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.54918140        0.00000000       -0.96144982
  1        -1.04526531       -0.32220986       -0.17200927
  1         0.49034191       -0.69881343        0.70513159
  1         0.00574201        1.02102328        0.42832749
  0        -0.36340139        0.00000000        0.63620545
  0         0.69166739        0.21321099        0.11382106
  0        -0.32446643        0.46241509       -0.46659592
  0        -0.00379957       -0.67562608       -0.28343058
  6         0.00000000        0.00000000        2.52342700
  1         0.94148608       -0.49948053        2.22324520
  1        -0.80787672       -0.75290477        2.60377353
  1         0.14135556        0.49568003        3.50337897
  1        -0.27496491        0.75670527        1.76331030
  0        -0.62299515        0.33051359        2.72206170
  0         0.53458388        0.49820813        2.47026052
  0        -0.09353704       -0.32799875        1.87497841
  0         0.18194832       -0.50072297        3.02640737

@aug-cc-pVTZ.gbs/N

