%Chk=BLIND01454
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.33522789        0.83925270       -0.47369406
  6        -0.14705079       -0.27186066       -1.34308646
 16        -0.78677531        0.73312393        1.22140424
 16        -0.30845571       -1.96233214       -0.84224420
  6         0.01116722       -0.78464835        1.67282103
  6         0.19813988       -1.85226385        0.85252692
  6        -0.17933880        2.01399059       -1.18150245
  7         0.13015735       -0.01519418       -2.59813147
 16         0.15436189        1.63171443       -2.83901123
  6         0.42047337       -0.81716189        3.03755580
  7         0.71420523       -0.80571905        4.16218182
  6         0.81177864       -3.04777384        1.33017246
  7         1.27997427       -4.05112488        1.68373388
  6        -0.29140530        3.34504460       -0.71149739
  7        -0.38866391        4.44385560       -0.34509451
  6         0.01463370       -0.90360874       11.61936400
  6        -0.24057746        0.19404818       12.48875643
 16        -0.44288964       -0.97996659        9.92426563
 16        -1.04187850        1.69126627       11.98791416
  6        -0.29251550        0.72817061        9.47284898
  6        -0.53205074        1.78523468       10.29314312
  6         0.61181012       -1.92717642       12.32717241
  7         0.11420966        0.06424827       13.74380150
 16         0.77212728       -1.44573128       13.98468125
  6         0.07253356        0.91612813        8.10811427
  7         0.34792597        1.01893123        6.98348831
  6        -0.42733208        3.12494763        9.81549769
  7        -0.38262938        4.23125779        9.46193635
  6         1.02211674       -3.19836175       11.85716733
  7         1.35645493       -4.24958143       11.49076443

@blind-aug.gbs/N

