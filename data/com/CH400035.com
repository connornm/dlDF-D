%Chk=CH400035
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.75678404        0.00000000       -0.80824741
  1         0.19302993        0.84713836        0.68635414
  1        -1.01044037        0.10472186       -0.44048784
  1         0.06062640       -0.95186022        0.56238111
  0        -0.50077510        0.00000000        0.53482917
  0        -0.12773074       -0.56056388       -0.45417060
  0         0.66862322       -0.06929599        0.29147727
  0        -0.04011738        0.62985987       -0.37213583
  6         0.00000000        0.00000000        2.66870296
  1         0.86856948        0.68606955        2.63901161
  1         0.26949050       -0.95785263        2.18302903
  1        -0.28616996       -0.18779181        3.72171148
  1        -0.85189003        0.45957489        2.13105972
  0        -0.57474517       -0.45398229        2.68835016
  0        -0.17832582        0.63382514        2.99008052
  0         0.18936286        0.12426460        1.97191178
  0         0.56370813       -0.30410745        3.02446938

@aug-cc-pVTZ.gbs/N
