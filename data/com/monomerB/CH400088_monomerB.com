%Chk=CH400088
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.71312934       0.00000000       0.84701388
  1-Bq       -1.03279802       0.09740814       0.38707345
  1-Bq        0.22436613       0.85141132      -0.67137513
  1-Bq        0.09530256      -0.94881946      -0.56271219
  6         0.00000000        0.00000000        3.08587769
  1         0.72758116       -0.64486259        2.55600614
  1         0.52988269        0.86179091        3.53590738
  1        -0.75890541        0.36723936        2.36811465
  1        -0.49855843       -0.58416768        3.88348258

@aug-cc-pVTZ.gbs/N

