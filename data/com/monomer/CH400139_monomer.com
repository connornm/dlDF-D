%Chk=CH400139
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.77483446        0.00000000        0.79095987
  1         0.29076268       -0.70642733       -0.80150075
  1        -0.96982714       -0.31240206        0.43338924
  1        -0.09577000        1.01882939       -0.42284836
  6-Bq        0.00000000       0.00000000       3.35179142
  1-Bq       -0.61045172       0.19117475       4.25555404
  1-Bq        1.07113187      -0.03365151       3.63023360
  1-Bq       -0.29576197      -0.96886300       2.90478438
  1-Bq       -0.16491818       0.81133977       2.61659367

@aug-cc-pVTZ.gbs/N

