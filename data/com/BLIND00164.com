%Chk=BLIND00164
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.67078760        0.63575668        0.43239605
  6        -0.13199379        1.11944446       -0.79299125
 16        -0.91408138       -1.05811052        0.83256937
 16         0.41648821        0.07805425       -2.11558804
  6         0.44393022       -1.79319133       -0.03897109
  6         0.96844145       -1.34101247       -1.20861476
  6        -1.03936023        1.70053971        1.22970059
  7        -0.07363314        2.41878929       -0.95472272
 16        -0.71117183        3.17708434        0.38298051
  6         0.93682376       -2.97468047        0.58731081
  7         1.28812512       -3.94528359        1.12182926
  6         2.03380521       -2.03576970       -1.85379256
  7         2.88133848       -2.58523395       -2.42885996
  6        -1.62453164        1.67232800        2.51886582
  7        -2.11361208        1.66793864        3.57332159
  6        -0.46273284       -0.51215462        6.35534378
  6         0.33384657       -1.33377826        5.50904448
 16        -0.14333379        1.17751619        6.71799159
 16         1.77511097       -0.76990442        4.64916569
  6         0.58747528        1.70226895        5.19002735
  6         1.34887478        0.92796900        4.37248565
  6        -1.48289669       -1.26265405        6.90426959
  7        -0.01130026       -2.59281046        5.39254727
 16        -1.34111680       -2.89768978        6.34625239
  6         0.34055723        3.07428202        4.89388337
  7         0.13244607        4.20263111        4.70757416
  6         1.93085257        1.46300175        3.18547943
  7         2.44823533        1.86393888        2.22501568
  6        -2.49293134       -0.84732937        7.80565955
  7        -3.32344902       -0.52613043        8.55276566

@aug-cc-pVTZ.gbs/N
