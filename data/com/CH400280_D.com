%Chk=CH400280
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.74540773        0.00000000        0.81875104
  1         0.06391938        0.95461638       -0.55732202
  1         0.20665648       -0.84316970       -0.68727281
  1        -1.01598359       -0.11144668        0.42584379
  6         0.00000000        0.00000000        3.89606002
  1         0.74432026       -0.00459402        4.71578692
  1         0.10699136       -0.92385319        3.29519889
  1        -1.02013979        0.04582487        4.32407970
  1         0.16882816        0.88262234        3.24917455

@aug-cc-pVTZ.gbs/N

