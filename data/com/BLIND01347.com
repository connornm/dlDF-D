%Chk=BLIND01347
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.87950361        0.21480110       -0.47057708
  6         0.69870391       -1.17090924       -0.20049382
 16         0.23877073        1.52835082        0.50541736
 16        -0.23959663       -1.80232692        1.16167068
  6        -1.26941211        0.79531522        1.08175531
  6        -1.45556500       -0.52569940        1.34235375
  6         1.69383867        0.37014504       -1.57422100
  7         1.29588701       -2.02334348       -0.99710429
 16         2.17119263       -1.19525916       -2.14559606
  6        -2.30567842        1.74828456        1.30332067
  7        -3.11119490        2.56493780        1.49150602
  6        -2.69856146       -1.00690771        1.84964095
  7        -3.68336376       -1.43870767        2.29082047
  6         2.13727341        1.57211560       -2.17738005
  7         2.51827611        2.54866080       -2.67966715
  6        -0.40634938        0.87649674        6.83376814
  6        -1.10375460       -0.28058056        6.38585661
 16         0.46958894        1.00606591        8.35181311
 16        -1.17520744       -1.80543405        7.28268696
  6         1.04717687       -0.66050699        8.53360703
  6         0.39194661       -1.77345352        8.10973514
  6        -0.58222789        1.90374985        5.92871511
  7        -1.75768723       -0.18975083        5.25361013
 16        -1.60317871        1.35056677        4.64167413
  6         2.28734588       -0.74996897        9.23007456
  7         3.28690250       -0.77371434        9.82296220
  6         0.92339689       -3.07516444        8.34858913
  7         1.30406105       -4.15531679        8.54681079
  6        -0.07520530        3.22455252        5.98852173
  7         0.32480698        4.31528257        6.02585872

@aug-cc-pVTZ.gbs/N

