%Chk=BLIND00678
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.43426758       -0.04249823        0.92234184
  6        -0.44679364        1.25771097        0.34344133
 16        -0.49068499       -1.55145182        0.02323520
 16        -0.48954670        1.56780362       -1.39914384
  6         0.38975305       -1.08540242       -1.44364888
  6         0.38712667        0.15210529       -2.00597603
  6        -0.44332724        0.06215601        2.29866616
  7        -0.46078521        2.28173258        1.16139243
 16        -0.49285495        1.73911852        2.73475138
  6         1.09711957       -2.16952899       -2.03979638
  7         1.63757310       -3.08920591       -2.50163766
  6         1.09373871        0.41056005       -3.21752153
  7         1.63047701        0.65902833       -4.21807220
  6        -0.44714245       -0.97807559        3.25943258
  7        -0.45979849       -1.81906843        4.06172105
  6        -0.41968611        0.03673444        9.68042892
  6        -1.27688454        0.50797519        8.64651903
 16         0.81136086       -1.20268377        9.48838071
 16        -1.21859469       -0.06315060        6.97170383
  6         1.31118741       -0.90337246        7.81365439
  6         0.50346518       -0.45271799        6.81759245
  6        -0.75008356        0.64928703       10.87245243
  7        -2.18399947        1.39579544        8.97328765
 16        -2.09044187        1.71528076       10.60426680
  6         2.67720555       -1.22910129        7.57078646
  7         3.78901775       -1.53091542        7.41594821
  6         0.99551194       -0.29017586        5.48890438
  7         1.34548628       -0.17071766        4.38700882
  6        -0.16292690        0.46520303       12.14780673
  7         0.30089492        0.31982911       13.20366262

@aug-cc-pVTZ.gbs/N

