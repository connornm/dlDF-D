%Chk=CH400066
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.05381208        0.00000000       -0.33980295
  1        -0.67039451        0.09196073       -0.87641336
  1        -0.21614970       -0.94652542        0.53231090
  1        -0.16726787        0.85456468        0.68390541
  0        -0.69732291        0.00000000        0.22485260
  0         0.44360988       -0.06085177        0.57993558
  0         0.14302943        0.62632975       -0.35223793
  0         0.11068361       -0.56547799       -0.45255024
  6         0.00000000        0.00000000        3.77452784
  1         1.03472733       -0.24863506        3.46874046
  1         0.01244882        0.89667475        4.42399329
  1        -0.61174528        0.20457509        2.87458190
  1        -0.43543087       -0.85261478        4.33079571
  0        -0.68469425        0.16452547        3.97687184
  0        -0.00823757       -0.59334284        3.34476703
  0         0.40480083       -0.13537034        4.37003527
  0         0.28813099        0.56418771        3.40643723

@aug-cc-pVTZ.gbs/N

