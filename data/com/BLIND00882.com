%Chk=BLIND00882
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.95662938        0.00263283        0.35491154
  6        -0.51609049        1.27300611       -0.11189020
 16        -0.41963037       -1.54844461       -0.27302752
 16         0.68298780        1.49759174       -1.39499616
  6         1.25089505       -1.15465451       -0.71844448
  6         1.68475956        0.05390151       -1.16431214
  6        -1.92228050        0.17440622        1.32615561
  7        -1.06318356        2.33576840        0.42577654
 16        -2.20525514        1.87067400        1.54393422
  6         2.13569791       -2.26664031       -0.61027270
  7         2.81304146       -3.20780446       -0.52897830
  6         3.04548559        0.25282688       -1.54219480
  7         4.13701994        0.45209304       -1.88827556
  6        -2.63270391       -0.81764229        2.04470052
  7        -3.23175602       -1.61843247        2.63712567
  6         0.35782037       -0.90614863        6.27983341
  6         1.00927141        0.20065280        5.66629548
 16        -0.33816805       -0.90368300        7.89355912
 16         1.19752114        1.78801568        6.42749812
  6        -0.87585215        0.78199340        8.00980446
  6        -0.26391567        1.84739209        7.42833785
  6         0.41877854       -2.00250697        5.44342100
  7         1.52729159        0.01191062        4.47722377
 16         1.28929058       -1.56762269        4.00887891
  6        -2.02648034        0.94660322        8.83445471
  7        -2.95066340        1.03347468        9.53412037
  6        -0.75247755        3.17286119        7.62424823
  7        -1.09795638        4.27164794        7.77993441
  6        -0.08984622       -3.30558159        5.66357052
  7        -0.49269843       -4.38292718        5.83129801

@aug-cc-pVTZ.gbs/N

