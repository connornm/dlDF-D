%Chk=CH400076
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.56350451       0.00000000      -0.95312571
  1-Bq       -1.07397759       0.17332609      -0.20619444
  1-Bq        0.12602479      -0.97817441       0.50326784
  1-Bq        0.38444830       0.80484832       0.65605232
  6         0.00000000        0.00000000        2.66413720
  1         0.74761511        0.50678314        3.30462800
  1         0.43542826       -0.18308940        1.66270500
  1        -0.28629551       -0.96555958        3.12426276
  1        -0.89674785        0.64186584        2.56495306

@aug-cc-pVTZ.gbs/N

