%Chk=CH400037
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.03826380       0.00000000      -0.38470018
  1-Bq        0.01494884       0.09982302       1.10263224
  1-Bq       -0.55664224       0.85000543      -0.44002971
  1-Bq       -0.49657040      -0.94982845      -0.27790235
  6         0.00000000        0.00000000        5.18782279
  1         0.70274138       -0.09518695        6.03816387
  1        -0.70426153       -0.85438171        5.19361546
  1        -0.56705021        0.94636622        5.28179612
  1         0.56857035        0.00320243        4.23771572

@aug-cc-pVTZ.gbs/N

