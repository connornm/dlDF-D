%Chk=BLIND01088_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.39935742       -0.33545547        9.28165163
  6        -1.25368847        0.54939723        8.56542486
 16         0.76878154       -1.42678037        8.55159820
 16        -1.25461632        0.71827344        6.80302073
  6         1.24928770       -0.48071214        7.13096778
  6         0.44414985        0.36809127        6.43892018
  6        -0.67766312       -0.25814819       10.63141929
  7        -2.11202638        1.25396621        9.26162785
 16        -1.97199305        0.86888231       10.87501496
  6         2.59370398       -0.72395031        6.72541914
  7         3.68706219       -0.97300080        6.41925789
  6         0.91726459        1.04719543        5.27748920
  7         1.25089923        1.59841353        4.31013577
  6        -0.07498473       -0.97190955       11.69560533
  7         0.40204235       -1.55566462       12.58034885

@blind-aug.gbs/N

