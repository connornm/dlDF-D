%Chk=BLIND00949
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=80)

M05 opt

0 1
  6         0.33344217        0.22843219        0.93687991
  6        -0.34034540        1.30330185        0.29151786
 16        -0.00583739       -1.47712188        0.68292846
 16        -1.62331240        1.09173664       -0.91007449
  6        -0.48586607       -1.45932255       -1.02393350
  6        -1.12838028       -0.43959786       -1.65240872
  6         1.24409055        0.72991561        1.84492017
  7        -0.00502810        2.52135235        0.64048108
 16         1.15836623        2.46116969        1.82951791
  6        -0.17124949       -2.67034276       -1.70623557
  7         0.08377259       -3.68485515       -2.21310969
  6        -1.51299540       -0.54586461       -3.02161306
  7        -1.87073002       -0.60011185       -4.12615338
  6         2.11286075        0.01508972        2.70482175
  7         2.82594201       -0.55537275        3.42397654
  6         0.70122515       -0.03818032        5.98937095
  6         0.37451423        1.28479201        5.57848065
 16        -0.39950610       -1.40788668        5.96087835
 16        -1.20388555        1.76261770        4.93433386
  6        -1.38506835       -0.99960292        4.54451204
  6        -1.70339913        0.25869683        4.14079024
  6         1.99096141       -0.06648808        6.48037225
  7         1.29643421        2.20541537        5.72122587
 16         2.65446089        1.53342242        6.41050516
  6        -1.87725054       -2.14628139        3.85605374
  7        -2.27423310       -3.11065154        3.34264543
  6        -2.54516881        0.47722589        3.01047578
  7        -3.25037929        0.70164630        2.11423932
  6         2.70886329       -1.16974539        7.00249163
  7         3.30970290       -2.06212269        7.44266254

@aug-cc-pVTZ.gbs/N
