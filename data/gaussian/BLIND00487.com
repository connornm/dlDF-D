%Chk=BLIND00487
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.45965950       -0.31978418       -0.85297136
  6        -0.44106968       -1.29238674        0.18603073
 16        -0.50943961        1.42005819       -0.61009870
 16        -0.43275384       -0.90018558        1.91250037
  6         0.41600796        1.56162379        0.89576585
  6         0.44330130        0.64033983        1.89483292
  6        -0.50490343       -0.95160130       -2.07937395
  7        -0.46519281       -2.55394507       -0.16857987
 16        -0.54644538       -2.66613236       -1.82741347
  6         1.12645941        2.79248497        1.00220435
  7         1.66836299        3.81955123        1.05386104
  6         1.18560226        0.87376535        3.09000811
  7         1.75224871        1.03425886        4.09204592
  6        -0.54713810       -0.36705343       -3.36843890
  7        -0.59154531        0.09559568       -4.43384011
  6        -0.09011347        0.90435201        5.96010606
  6        -1.32450935        0.38066521        5.48277431
 16         1.27442059       -0.05894376        6.50657171
 16        -1.67673520       -1.34774753        5.33073912
  6         1.11332918       -1.47394076        5.45018393
  6        -0.05853213       -1.98338245        4.98666016
  6        -0.15363578        2.28274961        5.99599732
  7        -2.26351712        1.24008441        5.17048747
 16        -1.72544480        2.78716505        5.46702847
  6         2.36180321       -2.09072057        5.14689379
  7         3.40082927       -2.57216631        4.94732458
  6        -0.08081682       -3.15810500        4.17828385
  7        -0.14862483       -4.12836952        3.54195456
  6         0.84254300        3.19813656        6.41416355
  7         1.64520148        3.96337252        6.76237212

@blind-aug.gbs/N

