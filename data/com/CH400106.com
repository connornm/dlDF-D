%Chk=CH400106
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10574458        0.00000000        0.05757506
  1        -0.40987614        0.67756435        0.77388246
  1        -0.31742207        0.34897061       -1.00172287
  1        -0.37844637       -1.02653495        0.17026534
  0        -0.73168741        0.00000000       -0.03809826
  0         0.27122105       -0.44835427       -0.51208938
  0         0.21004284       -0.23091897        0.66285471
  0         0.25042352        0.67927324       -0.11266707
  6         0.00000000        0.00000000        3.32832217
  1        -0.46215784       -0.49826584        2.45417845
  1        -0.57946751       -0.24216589        4.24022076
  1         1.04176555       -0.35559926        3.44773598
  1        -0.00014020        1.09603099        3.17115349
  0         0.30581662        0.32970981        3.90675589
  0         0.38344215        0.16024472        2.72490550
  0        -0.68935155        0.23530525        3.24930430
  0         0.00009277       -0.72525979        3.43232299

@aug-cc-pVTZ.gbs/N

