%Chk=BLIND01325
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.45797780       -0.85028435        0.32921386
  6         0.53617972        0.23406439        1.24785358
 16         0.39115966       -0.69239867       -1.41982132
 16         0.54665692        1.93929060        0.77174208
  6        -0.50108945        0.83238079       -1.57168048
  6        -0.43626951        1.87467191       -0.70158183
  6         0.51472956       -2.04592252        1.01662337
  7         0.63986831       -0.06046391        2.52078790
 16         0.68414504       -1.71374832        2.70948641
  6        -1.29365304        0.90143118       -2.75412060
  7        -1.90522001        0.92040428       -3.74246208
  6        -1.16075793        3.07931433       -0.94219586
  7        -1.71014071        4.08928044       -1.11270435
  6         0.48080344       -3.36211939        0.49542931
  7         0.46368462       -4.44925547        0.08440861
  6         0.40557798        0.72717888        9.49038032
  6         0.66525741        0.90532074        8.10233899
 16        -1.02728429       -0.03706442       10.16220957
 16        -0.41220747        0.35002502        6.81187225
  6        -1.37257827       -1.23699234        8.90321508
  6        -1.12871091       -1.07997351        7.57523999
  6         1.41536192        1.30921104       10.22989244
  7         1.75661629        1.55118213        7.77119117
 16         2.55896742        2.03397806        9.14743884
  6        -1.99506105       -2.41510569        9.40883726
  7        -2.52217116       -3.34551856        9.86455629
  6        -1.48929334       -2.09165128        6.63692750
  7        -1.79233306       -2.87588754        5.83438993
  6         1.54794831        1.37477870       11.63818580
  7         1.66819834        1.44615402       12.79211015

@aug-cc-pVTZ.gbs/N

