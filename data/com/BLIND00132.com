%Chk=BLIND00132
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.44957627       -0.26046387        0.87815070
  6         0.23889475        1.14004188        0.73664591
 16        -0.57071496       -1.51265583        0.18566409
 16        -1.08359822        1.85519623       -0.19831662
  6        -1.09632386       -0.71132531       -1.30622289
  6        -1.30010640        0.62433865       -1.45494588
  6         1.53493820       -0.48574871        1.70067088
  7         1.06167305        1.94044742        1.36929658
 16         2.16348300        1.03961230        2.23267818
  6        -1.34283058       -1.62267233       -2.37382333
  7        -1.55349547       -2.40621268       -3.20621618
  6        -1.77058742        1.16332875       -2.68867163
  7        -2.18050022        1.64246475       -3.66505594
  6         2.08239405       -1.72443267        2.11434055
  7         2.53913605       -2.73193659        2.47129785
  6         0.93924234        0.38154848        6.60820978
  6         0.91142075       -0.88125294        5.95217228
 16        -0.33136372        1.59291145        6.52649916
 16        -0.43895764       -1.45538412        4.96149578
  6        -1.76031806        0.55203906        6.38939977
  6        -1.80010596       -0.65570401        5.76697095
  6         2.14588441        0.53965414        7.25960885
  7         1.97054612       -1.64452822        6.06846636
 16         3.12543094       -0.86226691        6.97688542
  6        -2.92836760        1.12603101        6.97016064
  7        -3.85835470        1.64519866        7.43575914
  6        -3.01666676       -1.39381647        5.67151545
  7        -3.99041469       -2.01561713        5.54460982
  6         2.59529001        1.64529489        8.02168600
  7         2.98406677        2.54767951        8.64269367

@aug-cc-pVTZ.gbs/N

