%Chk=CH400210
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.39682654       0.00000000       1.03368982
  1-Bq        0.19127459       0.98470951      -0.46877197
  1-Bq        0.50208359      -0.79249531      -0.58808945
  1-Bq       -1.09018472      -0.19221420       0.02317160
  6         0.00000000        0.00000000        2.73575394
  1        -0.74139069        0.47104740        3.40987571
  1         0.31998349        0.73252772        1.96959491
  1         0.87921177       -0.32622549        3.32443854
  1        -0.45780457       -0.87734962        2.23910661

@aug-cc-pVTZ.gbs/N

