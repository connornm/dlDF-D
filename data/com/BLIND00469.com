%Chk=BLIND00469
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.73034060        0.67398189        0.23121401
  6        -1.20484427       -0.66154108        0.10063334
 16         0.87649876        1.12830702        0.77909408
 16        -0.22975903       -2.09910030        0.44300831
  6         1.82739249       -0.22924522        0.14903444
  6         1.38625415       -1.50828048        0.01823000
  6        -1.74623163        1.55889986       -0.06910058
  7        -2.45295087       -0.82490830       -0.26460780
 16        -3.17866103        0.66129997       -0.45288050
  6         3.16550615        0.13150449       -0.18323222
  7         4.25320553        0.46816799       -0.41724238
  6         2.25003378       -2.53888079       -0.45683810
  7         2.92637636       -3.41507209       -0.81148672
  6        -1.71173121        2.97445326       -0.05543539
  7        -1.70413336        4.13669989       -0.04057096
  6         0.39356755        0.92796025        4.85303357
  6         1.32207468       -0.00140323        4.30532499
 16        -1.34046677        0.91397353        4.56764616
 16         0.86983930       -1.35747858        3.26075276
  6        -1.62582200       -0.82838834        4.40372662
  6        -0.74778867       -1.72637775        3.88380687
  6         1.06996231        1.88371890        5.58400468
  7         2.59169050        0.18752301        4.57035814
 16         2.77021898        1.55677246        5.49998194
  6        -2.92360466       -1.22606289        4.83812825
  7        -4.00018235       -1.49801614        5.18195222
  6        -1.09543956       -3.10331377        3.75293294
  7        -1.34641355       -4.22792087        3.60062709
  6         0.53859574        2.99390113        6.28422852
  7         0.12072283        3.91528353        6.85653318

@aug-cc-pVTZ.gbs/N
