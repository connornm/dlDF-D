%Chk=BLIND00171
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.94999322        0.24166415       -0.28323176
  6        -0.75022600       -1.15316808       -0.08232685
 16        -0.10068513        1.51962627        0.57350928
 16         0.42379876       -1.83392278        1.05469876
  6         1.47488604        0.74628257        0.82577324
  6         1.67946270       -0.58374767        1.01763775
  6        -1.96211942        0.43809863       -1.20101020
  7        -1.51131845       -1.97534260       -0.76249401
 16        -2.57716237       -1.10430530       -1.69854297
  6         2.55530833        1.67513056        0.85769086
  7         3.39996949        2.47261249        0.89968034
  6         2.98707672       -1.09901504        1.25968542
  7         4.02969052       -1.55872881        1.48891421
  6        -2.48897954        1.66199417       -1.68024558
  7        -2.93976269        2.65701596       -2.07751390
  6        -1.01649112        0.07652100        5.12559391
  6        -0.49529693       -1.04265225        5.83385553
 16        -0.18281276        1.60638318        4.89533600
 16         1.11013972       -1.08550762        6.57905679
  6         1.49828317        1.04816648        4.81513204
  6         2.00828234       -0.01893169        5.48509204
  6        -2.30372444       -0.19652238        4.70872352
  7        -1.27366267       -2.08901049        5.96478006
 16        -2.74598499       -1.78237281        5.25122751
  6         2.32162786        1.86902071        3.99087996
  7         2.95759693        2.58129210        3.32800301
  6         3.39065088       -0.35584582        5.38768561
  7         4.51495011       -0.64938027        5.36094790
  6        -3.18953285        0.63889226        3.98581669
  7        -3.93295135        1.31443909        3.40095660

@aug-cc-pVTZ.gbs/N
