%Chk=BLIND02594_mA
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -1.00976119       -0.13853881        0.04793996
  6        -0.62773362        1.21497445       -0.17088290
 16        -0.10915645       -1.54544897       -0.49798268
 16         0.85149952        1.70519069       -1.01120766
  6         1.54965007       -0.92881654       -0.38731440
  6         1.92758784        0.36060905       -0.59271339
  6        -2.23046677       -0.17806095        0.69107223
  7        -1.44577820        2.14806131        0.25103178
 16        -2.79054047        1.44424347        0.93453145
  6         2.50595697       -1.94476517       -0.09677157
  7         3.25059752       -2.81118907        0.11773507
  6         3.29958003        0.74412363       -0.52566182
  7         4.40731916        1.09567706       -0.50983950
  6        -2.97896205       -1.31286052        1.08749775
  7        -3.61100739       -2.23297412        1.41154337

@blind-aug.gbs/N
