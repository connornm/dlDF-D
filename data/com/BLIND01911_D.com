%Chk=BLIND01911
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.94159200        0.21116108        0.33154851
  6        -0.30744109        1.30379053       -0.32407346
 16        -0.66276811       -1.48599313       -0.02961830
 16         0.90347247        1.12609233       -1.60340924
  6         1.04613235       -1.43855684       -0.50028919
  6         1.66398678       -0.40171919       -1.12548328
  6        -1.86101329        0.68779565        1.24412026
  7        -0.67462640        2.51230242        0.02597508
 16        -1.86868280        2.42011272        1.18210971
  6         1.74319521       -2.64445120       -0.19860756
  7         2.26248552       -3.65542648        0.04523578
  6         3.03646098       -0.48419574       -1.50423666
  7         4.14345488       -0.51860619       -1.85680800
  6        -2.71556956       -0.05042906        2.09847751
  7        -3.43055637       -0.64028052        2.79987970
  6         0.94037825        0.21539945        3.86644395
  6         0.30318095        1.30405791        4.52570795
 16         0.66815342       -1.48366572        4.22363833
 16        -0.90536095        1.11872912        5.80620318
  6        -1.04027135       -1.44359106        4.69671433
  6        -1.66101239       -0.41058739        5.32538598
  6         1.85684765        0.69767421        2.95386622
  7         0.66554089        2.51477557        4.17826294
 16         1.85836341        2.42984360        3.02030012
  6        -1.73339351       -2.65121513        4.39287981
  7        -2.24936836       -3.66342910        4.14714445
  6        -3.03266942       -0.49897702        5.70576238
  7        -4.13905643       -0.53827738        6.05972479
  6         2.71290501       -0.03528094        2.09648008
  7         3.42906570       -0.62075804        1.39261479

@aug-cc-pVTZ.gbs/N

