%Chk=CH400294
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.90977888        0.00000000        0.63110089
  1        -0.82585309        0.49911328        0.54299039
  1         0.20440596        0.54447637       -0.94220467
  1        -0.28833175       -1.04358965       -0.23188661
  0        -0.60201403        0.00000000       -0.41760872
  0         0.54647911       -0.33027058       -0.35930471
  0        -0.13525842       -0.36028800        0.62347065
  0         0.19079335        0.69055858        0.15344277
  6         0.00000000        0.00000000        3.73674034
  1        -0.55150605       -0.58857399        4.49529666
  1         1.04026464       -0.37247088        3.66531992
  1        -0.49824541       -0.10586073        2.75361727
  1         0.00948682        1.06690560        4.03272750
  0         0.36493965        0.38946804        3.23479249
  0        -0.68835837        0.24646944        3.78400027
  0         0.32969630        0.07004960        4.38728729
  0        -0.00627757       -0.70598708        3.54088129

@aug-cc-pVTZ.gbs/N

