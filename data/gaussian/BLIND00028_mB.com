%Chk=BLIND00028_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.63456232       -0.65237700        5.14026597
  6        -0.01282370        0.26748376        4.24969091
 16         0.07844020       -1.26164212        6.62651631
 16         1.58617908        0.97961267        4.51549324
  6         1.03390655        0.14608539        7.12602544
  6         1.63167400        1.03229622        6.28641332
  6        -1.85184767       -1.04799855        4.62354245
  7        -0.65910796        0.57380620        3.15135434
 16        -2.09152968       -0.27210845        3.09208118
  6         1.17018023        0.25643137        8.54036024
  7         1.27301533        0.29076363        9.69770170
  6         2.42029767        2.10741449        6.79251148
  7         3.09499976        2.97992673        7.15915244
  6        -2.78758934       -1.95649618        5.17507487
  7        -3.56267506       -2.70652962        5.60841987

@blind-aug.gbs/N
