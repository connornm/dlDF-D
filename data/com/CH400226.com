%Chk=CH400226
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.84569878        0.00000000        0.71468842
  1        -0.95428419        0.06797751        0.55741072
  1         0.09229151        0.86815218       -0.68101397
  1         0.01629390       -0.93612969       -0.59108516
  0        -0.55961129        0.00000000       -0.47291981
  0         0.63146384       -0.04498172       -0.36884685
  0        -0.06107064       -0.57446903        0.45063694
  0        -0.01078191        0.61945075        0.39112972
  6         0.00000000        0.00000000        3.33492641
  1        -0.65131625        0.78546469        2.90502141
  1         0.06175259        0.13169578        4.43257335
  1        -0.42383804       -0.99698423        3.10608208
  1         1.01340170        0.07982375        2.89602880
  0         0.43098552       -0.51975351        3.61940086
  0        -0.04086259       -0.08714503        2.60859733
  0         0.28045985        0.65971909        3.48635606
  0        -0.67058277       -0.05282055        3.62535140

@aug-cc-pVTZ.gbs/N

