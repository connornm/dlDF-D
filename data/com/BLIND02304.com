%Chk=BLIND02304
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.44443150        0.89871215       -0.18948885
  6         0.51639786       -0.01587148       -1.27769117
 16         0.40961208        0.44998613        1.50931932
 16         0.55064546       -1.77643622       -1.09376942
  6        -0.46599177       -1.08920010        1.41796430
  6        -0.40693798       -1.97059081        0.38485873
  6         0.47832712        2.19299771       -0.66800654
  7         0.59515942        0.48844246       -2.48488796
 16         0.62114688        2.15041942       -2.39503178
  6        -1.23711937       -1.36426825        2.58449308
  7        -1.83113008       -1.55538518        3.56521659
  6        -1.11616463       -3.20700708        0.43191533
  7        -1.65330867       -4.23770691        0.43967865
  6         0.44158349        3.40309450        0.06648928
  7         0.42179527        4.40597801        0.65379586
  6        -0.44443121       -0.89871233       10.27427869
  6        -0.51639698        0.01587071       11.36248154
 16        -0.40961304       -0.44998541        8.57547073
 16        -0.55064502        1.77643555       11.17856074
  6         0.46599059        1.08920093        8.66682595
  6         0.40693736        1.97059109        9.69993202
  6        -0.47832626       -2.19299815       10.75279572
  7        -0.59515761       -0.48844387       12.56967812
 16        -0.62114483       -2.15042078       12.47982108
  6         1.23711733        1.36426982        7.50029678
  7         1.83112733        1.55538737        6.51957295
  6         1.11616375        3.20700750        9.65287557
  7         1.65330760        4.23770744        9.64511243
  6        -0.44158292       -3.40309454       10.01829924
  7        -0.42179493       -4.40597774        9.43099213

@aug-cc-pVTZ.gbs/N

