%Chk=BLIND00605
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.41897769        0.36177108       -0.85713960
  6        -0.43347969       -1.05870955       -0.76853291
 16        -0.49944605        1.46219177        0.51074774
 16        -0.50361633       -1.95732269        0.75533179
  6         0.35901521        0.51336284        1.73830222
  6         0.35431639       -0.84257691        1.83343038
  6        -0.40432316        0.74392212       -2.18343106
  7        -0.42765524       -1.73297188       -1.89245301
 16        -0.43655174       -0.67547595       -3.17796462
  6         1.04981946        1.32115056        2.68770106
  7         1.57696859        2.02173239        3.45091263
  6         1.04202808       -1.50771693        2.89106137
  7         1.56337335       -2.08983059        3.75145020
  6        -0.39820692        2.05400877       -2.72083593
  7        -0.40240372        3.12208477       -3.17941196
  6         0.41782468        0.43195659        7.02961167
  6         0.38177156        1.28676821        5.89214785
 16        -0.75300894       -0.82941453        7.38543209
 16        -0.84610790        1.19457813        4.61998381
  6        -1.18202452       -1.34308025        5.74341649
  6        -1.21930466       -0.53772762        4.64903932
  6         1.47062964        0.79109744        7.84687710
  7         1.30032801        2.21795244        5.80842213
 16         2.28587980        2.15204473        7.14829553
  6        -1.54450743       -2.71911008        5.66443737
  7        -1.85450286       -3.83937696        5.65453915
  6        -1.62496481       -1.04243546        3.37834983
  7        -1.97676884       -1.40328035        2.33099554
  6         1.86216398        0.21581101        9.08013085
  7         2.19208987       -0.23803033       10.09810574

@aug-cc-pVTZ.gbs/N
