%Chk=BLIND00540
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.29100959        0.74937555        0.62837776
  6        -0.97349689        0.97543843        0.01539686
 16         0.95716657       -0.82835697        1.02272737
 16        -2.07903768       -0.31415572       -0.48391749
  6         0.27975882       -1.80245531       -0.29498780
  6        -0.92519857       -1.59679887       -0.88948519
  6         0.88296171        1.95905035        0.93100423
  7        -1.34696025        2.22016230       -0.15503678
 16        -0.18317164        3.24099215        0.45686138
  6         1.11493951       -2.89477587       -0.66968588
  7         1.81845985       -3.78462340       -0.92350380
  6        -1.39896434       -2.47046155       -1.91232508
  7        -1.83728410       -3.17279454       -2.72812548
  6         2.13234393        2.19487504        1.55438151
  7         3.15069353        2.40692418        2.07312535
  6         0.71431448       -0.68884338        3.80425123
  6         1.23820144        0.50521407        4.37486136
 16        -0.86420023       -1.37370244        4.16170893
 16         0.35996415        1.52284126        5.52712977
  6        -1.79774166        0.10492136        4.45546759
  6        -1.30970169        1.25154716        4.99833988
  6         1.66191191       -1.26511817        2.98250328
  7         2.46027662        0.84527280        4.04508516
 16         3.09991926       -0.29779690        3.01787250
  6        -3.17285061       -0.02241575        4.10340185
  7        -4.28981267       -0.18197742        3.82388645
  6        -2.15947351        2.37178070        5.23689098
  7        -2.82031088        3.29764505        5.47562787
  6         1.56343644       -2.46128569        2.23106847
  7         1.50351037       -3.44530847        1.61527307

@aug-cc-pVTZ.gbs/N

