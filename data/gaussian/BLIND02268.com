%Chk=BLIND02268
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.45669272        0.83190630        0.37479671
  6         0.38479452        0.81000955       -1.04653289
 16         0.57543575       -0.58742162        1.40432330
 16         0.37105587       -0.67085255       -2.01681482
  6        -0.36218636       -1.75416850        0.45369795
  6        -0.44015312       -1.78544882       -0.90298431
  6         0.49112284        2.13880192        0.81767551
  7         0.36109723        1.96435039       -1.66685191
 16         0.45845242        3.19833540       -0.55387399
  6        -1.02233014       -2.72862730        1.25727116
  7        -1.52183677       -3.51487409        1.95261435
  6        -1.18709525       -2.79874286       -1.57324887
  7        -1.75914338       -3.62178206       -2.16171598
  6         0.57371109        2.61398668        2.14904588
  7         0.65061329        3.02165443        3.23485847
  6         0.39914617       -0.88605014        7.22367348
  6         0.36016041        0.11133623        6.20901728
 16         0.57734233       -0.57361313        8.94389029
 16         0.45402924        1.85113021        6.52349125
  6        -0.26957690        0.97837369        9.07874882
  6        -0.31550679        1.93844360        8.11760352
  6         0.34772940       -2.13844367        6.64559912
  7         0.28315948       -0.29617604        4.96568924
 16         0.28393276       -1.96022319        4.92252222
  6        -0.88960984        1.16894452       10.34773722
  7        -1.35785813        1.28850064       11.40492353
  6        -0.98724850        3.17562398        8.34583522
  7        -1.49730658        4.20881575        8.49841849
  6         0.37303034       -3.40264676        7.28303056
  7         0.40202295       -4.44879469        7.78882684

@blind-aug.gbs/N

