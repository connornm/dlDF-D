%Chk=BLIND00792
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.90238321       -0.04249824        0.47434940
  6         0.32349105        1.25771096        0.46144344
 16         0.00166644       -1.55145182        0.49123198
 16        -1.41929063        1.56780363        0.42763711
  6        -1.42514210       -1.08540240       -0.45276876
  6        -1.98704218        0.15210531       -0.47483711
  6         2.27698219        0.06215599        0.54383570
  7         1.14003883        2.28173257        0.51133833
 16         2.71047201        1.73911849        0.61246445
  6        -1.98965363       -2.16952897       -1.18563026
  7        -2.42731776       -3.08920589       -1.74584228
  6        -3.16639125        0.41056008       -1.23396742
  7        -4.14240829        0.65902837       -1.81412290
  6         3.23665437       -0.97807562        0.58983517
  7         4.03761326       -1.81906847        0.63770806
  6        -0.79555755        0.60970126        3.88437935
  6        -1.15292948       -0.75508940        4.07268342
 16         0.56830940        1.41351385        4.64744480
 16        -0.25077560       -1.88535829        5.09405454
  6         1.71889898        0.06822685        4.74976080
  6         1.39100691       -1.23878352        4.92823974
  6        -1.72377018        1.22414919        3.06818190
  7        -2.23793390       -1.18229241        3.47424962
 16        -2.94698404        0.07246536        2.64133442
  6         3.07863918        0.48781835        4.67140765
  7         4.17174783        0.88063987        4.62671126
  6         2.39951581       -2.24039369        5.04498210
  7         3.18740081       -3.08503234        5.17439975
  6        -1.75846516        2.57204807        2.63559283
  7        -1.80784246        3.67708984        2.27846470

@blind-aug.gbs/N
