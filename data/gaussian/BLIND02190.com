%Chk=BLIND02190
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.41421602        0.17468871       -0.91597917
  6         0.43537995       -1.19523004       -0.53036104
 16         0.49478556        1.53873603        0.18919054
 16         0.51510817       -1.75258503        1.14820596
  6        -0.35493031        0.86656835        1.59286274
  6        -0.34383346       -0.43891889        1.97131217
  6         0.39299712        0.26896190       -2.29291982
  7         0.42844569       -2.09103910       -1.48709406
 16         0.42791377       -1.32789427       -2.96646101
  6        -1.04584832        1.85366153        2.35402568
  7        -1.57331786        2.69733792        2.95500742
  6        -1.02469090       -0.86892424        3.14837227
  7        -1.54027783       -1.25869787        4.11438095
  6         0.37906690        1.43649612       -3.09406790
  7         0.37681812        2.38408135       -3.76724922
  6         0.43593229        0.54689270        8.84409622
  6         0.56346734        1.13195150       10.13532341
 16         0.34590485       -1.17777797        8.51863844
 16         0.62054069        0.20588531       11.64321345
  6        -0.49967561       -1.73637106        9.97366061
  6        -0.38763812       -1.18634810       11.21162927
  6         0.46568866        1.53575964        7.88151638
  7         0.67829255        2.43609579       10.19763618
 16         0.67294514        3.07116740        8.65894441
  6        -1.30865328       -2.88690324        9.74361607
  7        -1.93520612       -3.83828371        9.51231594
  6        -1.07719436       -1.74352730       12.32877092
  7        -1.59682238       -2.18103807       13.27190250
  6         0.38144441        1.39621823        6.47488937
  7         0.32307226        1.29960737        5.31801645

@blind-aug.gbs/N

