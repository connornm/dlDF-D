%Chk=BLIND01183
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.40451042        0.67633101        0.64811773
  6         0.85387644        1.07554704        0.11621949
 16        -0.90431384       -0.98522826        0.92728428
 16         2.12213394       -0.05312042       -0.38599661
  6        -0.05762111       -1.80330898       -0.39860313
  6         1.14311796       -1.43200322       -0.91629200
  6        -1.14456478        1.79384960        0.97792823
  7         1.09051450        2.36166234        0.02807255
 16        -0.21025842        3.21093173        0.62610951
  6        -0.74364775       -2.95952639       -0.87148621
  7        -1.32795230       -3.90685340       -1.20667684
  6         1.76210974       -2.18945141       -1.95410863
  7         2.31642343       -2.79231022       -2.77886539
  6        -2.44118508        1.85404005        1.54389762
  7        -3.50081314        1.92218985        2.01680341
  6        -0.76232998        0.55549870        6.41436142
  6        -1.17684092       -0.70698419        5.90431799
 16         0.81422278        0.89740510        7.11156874
 16        -0.14675988       -2.14628773        5.86042356
  6         1.83282601       -0.19467928        6.15557645
  6         1.44922412       -1.40167365        5.66195893
  6        -1.81014966        1.45126563        6.34380130
  7        -2.41055781       -0.80752139        5.47339464
 16        -3.19617099        0.64459468        5.68573335
  6         3.15954646        0.29317819        5.97419954
  7         4.23552272        0.72134156        5.87348091
  6         2.36420198       -2.22645387        4.94328305
  7         3.08331711       -2.94352545        4.37794429
  6        -1.83667272        2.80706190        6.75157753
  7        -1.87921695        3.91846904        7.08930534

@aug-cc-pVTZ.gbs/N

