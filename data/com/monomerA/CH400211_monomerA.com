%Chk=CH400211
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.69626050        0.00000000        0.86093395
  1        -0.60810675       -0.92514813        0.01711946
  1        -0.66704955        0.88138151        0.06478809
  1         0.57889579        0.04376662       -0.94284150
  6-Bq        0.00000000       0.00000000       3.93513490
  1-Bq       -0.20874190       1.05382668       4.20320599
  1-Bq       -0.56119969      -0.26575115       3.01839219
  1-Bq        1.08487245      -0.12669707       3.75351556
  1-Bq       -0.31493086      -0.66137845       4.76542585

@aug-cc-pVTZ.gbs/N

