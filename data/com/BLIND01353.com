%Chk=BLIND01353
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.16508417       -0.33613073        0.94914291
  6         0.21565769        1.08104815        0.82718803
 16         0.60747017       -1.47709898       -0.31220666
 16         0.70821738        1.93354530       -0.64428127
  6         0.12625960       -0.56505244       -1.75468008
  6         0.16862911        0.78742411       -1.88373841
  6        -0.21551912       -0.67834831        2.23106782
  7        -0.09184295        1.78869876        1.88668508
 16        -0.44450717        0.77020545        3.15523395
  6        -0.28516574       -1.40058129       -2.83343417
  7        -0.58924913       -2.12345532       -3.69144318
  6        -0.19736503        1.42091526       -3.10794395
  7        -0.45809249        1.97724446       -4.09457718
  6        -0.37996307       -1.97173486        2.78358737
  7        -0.51079630       -3.02558761        3.25621087
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

@aug-cc-pVTZ.gbs/N

