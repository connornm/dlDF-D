%Chk=CH400251
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.37475255       0.00000000      -1.04189562
  1-Bq        0.71984586      -0.53274471       0.65114590
  1-Bq       -0.11315777       1.04384343       0.35152833
  1-Bq       -0.98144064      -0.51109873       0.03922140
  6         0.00000000        0.00000000        3.40178143
  1         0.17119627       -1.07764998        3.58979288
  1         0.81936305        0.58891253        3.85764564
  1        -0.02443995        0.18491985        2.31036341
  1        -0.96611938        0.30381760        3.84932379

@aug-cc-pVTZ.gbs/N

