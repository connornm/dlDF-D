%Chk=CH400088
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.71312934        0.00000000        0.84701388
  1        -1.03279802        0.09740814        0.38707345
  1         0.22436613        0.85141132       -0.67137513
  1         0.09530256       -0.94881946       -0.56271219
  0        -0.47188814        0.00000000       -0.56048151
  0         0.68341760       -0.06445640       -0.25613218
  0        -0.14846636       -0.56339136        0.44425877
  0        -0.06306310        0.62784775        0.37235491
  6         0.00000000        0.00000000        3.08587769
  1         0.72758116       -0.64486259        2.55600614
  1         0.52988269        0.86179091        3.53590738
  1        -0.75890541        0.36723936        2.36811465
  1        -0.49855843       -0.58416768        3.88348258
  0        -0.48145112        0.42671503        3.43650147
  0        -0.35063115       -0.57025968        2.78808644
  0         0.50217884       -0.24300767        3.56083202
  0         0.32990343        0.38655232        2.55809083

@aug-cc-pVTZ.gbs/N

