%Chk=CH400299
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.58501138        0.00000000        0.94007853
  1        -1.07392282        0.13457484        0.23359227
  1         0.34340574        0.82922871       -0.64841205
  1         0.14550570       -0.96380356       -0.52525875
  6-Bq        0.00000000       0.00000000       3.17602037
  1-Bq        1.03286498       0.21987060       2.84310475
  1-Bq       -0.72037079       0.44662159       2.46357315
  1-Bq       -0.15002208      -1.09637390       3.21401513
  1-Bq       -0.16247211       0.42988172       4.18338847

@aug-cc-pVTZ.gbs/N

