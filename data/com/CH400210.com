%Chk=CH400210
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.39682654        0.00000000        1.03368982
  1         0.19127459        0.98470951       -0.46877197
  1         0.50208359       -0.79249531       -0.58808945
  1        -1.09018472       -0.19221420        0.02317160
  0        -0.26258595        0.00000000       -0.68400772
  0        -0.12656920       -0.65159673        0.31019329
  0        -0.33223608        0.52440577        0.38914742
  0         0.72139122        0.12719096       -0.01533299
  6         0.00000000        0.00000000        2.73575394
  1        -0.74139069        0.47104740        3.40987571
  1         0.31998349        0.73252772        1.96959491
  1         0.87921177       -0.32622549        3.32443854
  1        -0.45780457       -0.87734962        2.23910661
  0         0.49058909       -0.31169897        2.28967768
  0        -0.21173777       -0.48472434        3.24273261
  0        -0.58178732        0.21586819        2.34621271
  0         0.30293600        0.58055511        3.06439276

@aug-cc-pVTZ.gbs/N
