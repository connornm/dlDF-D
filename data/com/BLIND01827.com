%Chk=BLIND01827
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.32539677        0.90823217       -0.33217493
  6         0.13138573        0.04919176       -1.37103295
 16        -1.06860222        0.37534278        1.16855529
 16         0.05246768       -1.71825050       -1.30389033
  6        -0.21801575       -1.15973908        1.42183063
  6         0.22413199       -1.98840747        0.43925334
  6        -0.16148112        2.22483545       -0.71295135
  7         0.60777748        0.61322823       -2.45397975
 16         0.50904999        2.26872947       -2.31085245
  6        -0.06617771       -1.49465800        2.79873447
  7         0.01248655       -1.73599609        3.93308930
  6         0.85706389       -3.22754248        0.75211684
  7         1.35122109       -4.25879533        0.96008886
  6        -0.49250017        3.39673641        0.00971244
  7        -0.76642085        4.36906029        0.58473387
  6        -0.00598042       -0.96474503       12.78956492
  6         0.14026329       -0.00286489       13.82842296
 16        -0.88379767       -0.70830387       11.28883461
 16        -0.51999385        1.63851889       13.76128029
  6        -0.59007570        1.02192765       11.03555932
  6        -0.44757306        1.95030197       12.01813663
  6         0.58504418       -2.15260125       13.17034139
  7         0.77667147       -0.37712513       14.91136983
 16         1.23222259       -1.97177389       14.76824256
  6        -0.55782522        1.38824090        9.65865549
  7        -0.56359615        1.64201017        8.52430066
  6        -0.26111465        3.32917529       11.70527316
  7        -0.13668663        4.46592128       11.49730117
  6         0.66114957       -3.36797504       12.44767759
  7         0.72497824       -4.37612790       11.87265616

@aug-cc-pVTZ.gbs/N
