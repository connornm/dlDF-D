%Chk=BLIND01554
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.80388395        0.62824295        0.01379481
  6         0.29909433        1.12683153       -0.73498752
 16        -1.20686122       -1.07070384        0.21305793
 16         1.47609366        0.10146980       -1.57038818
  6         0.41343277       -1.79064652        0.19179759
  6         1.47498828       -1.32419718       -0.51754645
  6        -1.54960880        1.68316225        0.49994595
  7         0.42239035        2.42800765       -0.83296110
 16        -0.83544839        3.16972960       -0.03383209
  6         0.50921795       -2.97616790        0.97694339
  7         0.53184928       -3.95053312        1.61042051
  6         2.72672993       -2.00746308       -0.50229930
  7         3.75534039       -2.54715739       -0.54250759
  6        -2.72999114        1.63916016        1.28092597
  7        -3.70449316        1.62180615        1.91428256
  6        -0.80383791        0.62829799        4.31396981
  6         0.29896973        1.12675977        3.56485179
 16        -1.20687510       -1.07060768        4.51346217
 16         1.47567230        0.10126018        2.72920235
  6         0.41336461       -1.79066085        4.49181070
  6         1.47475807       -1.32433276        3.78214452
  6        -1.54935924        1.68330161        4.80025011
  7         0.42232654        2.42792071        3.46675317
 16        -0.83524440        3.16978341        4.26617294
  6         0.50928411       -2.97613351        5.27701366
  7         0.53202261       -3.95045574        5.91055306
  6         2.72645788       -2.00768178        3.79709843
  7         3.75502104       -2.54744812        3.75664766
  6        -2.72953159        1.63943378        5.58155497
  7        -3.70386207        1.62218983        6.21517842

@aug-cc-pVTZ.gbs/N
