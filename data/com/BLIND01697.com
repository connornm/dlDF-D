%Chk=BLIND01697
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.44586690        0.91681942       -0.04187801
  6        -0.00221771        0.43807590       -1.30671266
 16        -0.98818008       -0.09395405        1.28956860
 16         0.10632233       -1.27598076       -1.73661634
  6         0.03452840       -1.52202553        1.04707880
  6         0.46580831       -1.98900050       -0.15441030
  6        -0.45331531        2.29712369       -0.04468251
  7         0.31145318        1.33125220       -2.21319304
 16         0.05792478        2.85901632       -1.60280215
  6         0.34697903       -2.19914141        2.26156566
  7         0.55487629       -2.72974759        3.27469391
  6         1.24826097       -3.17805247       -0.24309627
  7         1.86230004       -4.15733667       -0.36509872
  6        -0.83371372        3.17634092        0.99812263
  7        -1.15089314        3.91263382        1.83976738
  6         0.59471004       -0.74879399        8.26135938
  6         0.84469876        0.51602399        7.65839357
 16         0.02933160       -0.99136305        9.90757446
 16         0.59260215        2.06936630        8.46966623
  6        -0.98234256        0.44984411       10.11628338
  6        -0.75556762        1.66257292        9.54575204
  6         0.93177087       -1.75726691        7.38120281
  7         1.32244662        0.51820288        6.43789684
 16         1.54258301       -1.04897040        5.92173335
  6        -2.07701740        0.24109563       11.00465308
  7        -2.94036126        0.02942546       11.75366987
  6        -1.60732936        2.77369018        9.81752826
  7        -2.25717130        3.71293245       10.03306373
  6         0.85291099       -3.15753044        7.57668073
  7         0.80365880       -4.30958158        7.72313317

@aug-cc-pVTZ.gbs/N
