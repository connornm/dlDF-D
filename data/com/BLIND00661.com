%Chk=BLIND00661
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.76693981       -0.04249815        0.67164435
  6        -0.54686765        1.25771099        0.13605956
 16        -0.45901689       -1.55145182       -0.17497182
 16         0.11095919        1.56780351       -1.47815691
  6         0.93465716       -1.08540253       -1.16723831
  6         1.15717346        0.15210513       -1.68367353
  6        -1.32575461        0.06215619        1.92945066
  7        -0.88686076        2.28173267        0.88013235
 16        -1.54557603        1.73911873        2.30932135
  6         1.82142443       -2.16952915       -1.43068269
  7         2.50149175       -3.08920611       -1.63779560
  6         2.28940029        0.41055979       -2.51144494
  7         3.18153916        0.65902800       -3.21378241
  6        -1.71354516       -0.97807533        2.80848724
  7        -2.04604944       -1.81906811        3.53873923
  6         0.75901620       -0.11788628        8.02931345
  6         0.27342447       -1.34392717        7.49372866
 16         0.77155123        1.42211390        7.18269728
 16        -0.43449914       -1.51047353        5.87951219
  6        -0.68856479        1.25601003        6.19043079
  6        -1.16351087        0.09180862        5.67399557
  6         1.28386070       -0.33643781        9.28711976
  7         0.39308183       -2.41626006        8.23780145
 16         1.15021838       -2.02245814        9.66699045
  6        -1.33055147        2.50081516        5.92698641
  7        -1.80454610        3.54178892        5.71987350
  6        -2.32473159        0.07440501        4.84622416
  7        -3.24903452        0.01685253        4.14388669
  6         1.87945338        0.60043597       10.16615635
  7         2.37954387        1.35391951       10.89640834

@aug-cc-pVTZ.gbs/N

