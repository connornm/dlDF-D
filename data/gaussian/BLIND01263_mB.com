%Chk=BLIND01263_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.52889678        0.69891786        9.58416454
  6        -0.02245791        1.19030934        8.36746035
 16        -0.17696189       -0.56826871       10.57639352
 16        -1.51035235        0.57409410        7.63190800
  6        -0.93243163       -1.57318028        9.32597748
  6        -1.46148949       -1.11699923        8.15994116
  6         1.66343772        1.41749588        9.90317192
  7         0.60044879        2.18318826        7.78091223
 16         1.91747877        2.63331064        8.69400519
  6        -0.98440153       -2.95423688        9.67393072
  7        -1.02497456       -4.06563504       10.01207415
  6        -2.09070796       -2.00704720        7.24024520
  7        -2.63793952       -2.69075728        6.47604347
  6         2.51073874        1.26910711       11.02799644
  7         3.21169964        1.16755734       11.94964571

@blind-aug.gbs/N

