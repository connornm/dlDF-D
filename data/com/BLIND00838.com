%Chk=BLIND00838
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.36251596        0.87765686        0.37337584
  6        -0.67411044        0.66786014       -0.99947552
 16        -0.16072943       -0.39245090        1.57113422
 16        -0.87971531       -0.92813406       -1.73779624
  6         0.51973972       -1.68177424        0.56181310
  6         0.23160136       -1.89239068       -0.74977413
  6        -0.31349129        2.23179014        0.63647821
  7        -0.85016207        1.73004419       -1.74689167
 16        -0.67984863        3.10081308       -0.81796601
  6         1.39839839       -2.54686645        1.27642073
  7         2.08805078       -3.23841122        1.90676860
  6         0.79939684       -2.99111468       -1.45995425
  7         1.21539810       -3.88908240       -2.06948116
  6        -0.04905379        2.87911343        1.86782455
  7         0.15685367        3.42706919        2.87203694
  6         0.01888914       -1.00435177        5.92893319
  6        -1.22332263       -0.53830613        6.44423699
 16         1.11287232       -0.06697887        4.92242959
 16        -1.85942926        1.09223718        6.17724409
  6         0.82792807        1.56048207        5.56609079
  6        -0.35231752        2.01759474        6.06178466
  6         0.20009321       -2.32693195        6.28000295
  7        -1.94849182       -1.38799145        7.12982033
 16        -1.18077819       -2.86415093        7.17967872
  6         1.97252702        2.40506244        5.47852991
  7         2.92450095        3.06223338        5.36422761
  6        -0.48826561        3.36379365        6.51257124
  7        -0.65254266        4.45963349        6.86339063
  6         1.28712313       -3.17865912        5.96680137
  7         2.16660362       -3.89442135        5.71127366

@aug-cc-pVTZ.gbs/N
