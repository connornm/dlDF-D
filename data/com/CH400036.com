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
  0        -0.53387403        0.00000000        0.50179324
  0        -0.15383961        0.49240799       -0.52027462
  0         0.05180022       -0.66576258       -0.30148776
  0         0.63591342        0.17335458        0.31996914
  6         0.00000000        0.00000000        2.50094091
  1        -1.02072168       -0.29776239        2.19198983
  1         0.61032910       -0.90619735        2.68063706
  1         0.46602434        0.60600172        1.69996076
  1        -0.05563176        0.59795801        3.43117597
  0         0.67542651        0.19703374        2.70537836
  0        -0.40386372        0.59964408        2.38203333
  0        -0.30837514       -0.40100023        3.03096122
  0         0.03681235       -0.39567759        1.88539072

@aug-cc-pVTZ.gbs/N

