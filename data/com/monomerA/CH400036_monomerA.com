%Chk=CH400036
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.80680399        0.00000000       -0.75832268
  1         0.23248633       -0.74413946        0.78625221
  1        -0.07828181        1.00611730        0.45561595
  1        -0.96100850       -0.26197784       -0.48354548
  6-Bq        0.00000000       0.00000000       2.50094091
  1-Bq       -1.02072168      -0.29776239       2.19198983
  1-Bq        0.61032910      -0.90619735       2.68063706
  1-Bq        0.46602434       0.60600172       1.69996076
  1-Bq       -0.05563176       0.59795801       3.43117597

@aug-cc-pVTZ.gbs/N

