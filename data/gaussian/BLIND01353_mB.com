%Chk=BLIND01353_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.85474891       -0.25909076        8.70567242
  6        -1.00425935       -0.20601593       10.12011801
 16         0.25353607        0.73049474        7.76708753
 16        -0.06301569        0.85640971       11.17837151
  6         1.58725941        0.90458837        8.92250760
  6         1.45835229        0.95553320       10.27470272
  6        -1.76245841       -1.15414711        8.17629464
  7        -1.92178550       -0.97159750       10.65853295
 16        -2.71784765       -1.81373476        9.46354057
  6         2.85893852        1.03781592        8.29308539
  7         3.86988805        1.16644610        7.73395361
  6         2.59492962        1.14588927       11.11488036
  7         3.48830301        1.32853744       11.83560901
  6        -1.96353937       -1.50307217        6.81872837
  7        -2.14856896       -1.79442612        5.70878565

@blind-aug.gbs/N

