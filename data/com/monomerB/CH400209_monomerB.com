%Chk=CH400209
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.58119267       0.00000000      -0.94244418
  1-Bq       -0.23257066       1.04292045       0.29019610
  1-Bq       -0.94307796      -0.56097817      -0.14796421
  1-Bq        0.59445595      -0.48194228       0.80021230
  6         0.00000000        0.00000000        2.72381006
  1         1.07195452       -0.02616488        2.44773967
  1        -0.59806044        0.32625541        1.85096009
  1        -0.32525718       -1.01100217        3.03696871
  1        -0.14863690        0.71091164        3.55957179

@aug-cc-pVTZ.gbs/N
