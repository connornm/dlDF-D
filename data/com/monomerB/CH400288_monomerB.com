%Chk=CH400288
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.09903305       0.00000000      -0.13458200
  1-Bq       -0.44078016       0.84541539      -0.56300233
  1-Bq       -0.41811722      -0.95306473      -0.37793064
  1-Bq       -0.24013566       0.10764933       1.07551497
  6         0.00000000        0.00000000        3.46801559
  1        -0.57710338       -0.61909087        2.75410793
  1         0.79856625       -0.61414119        3.92748105
  1         0.45467491        0.85512639        2.93134481
  1        -0.67613778        0.37810566        4.25912858

@aug-cc-pVTZ.gbs/N

