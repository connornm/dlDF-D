%Chk=BLIND01086
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.25450981       -0.79126179       -0.59180946
  6         1.04279605       -0.90095066       -0.01666382
 16        -1.04954537        0.71760659       -1.01578034
 16         2.06441306        0.48400712        0.39900309
  6        -0.40415925        1.78651523        0.24320066
  6         0.83166230        1.69209907        0.80129986
  6        -0.76379442       -2.05152812       -0.83193688
  7         1.51383196       -2.10769875        0.18269578
 16         0.41033895       -3.23321332       -0.35250233
  6        -1.30607505        2.82728514        0.60985520
  7        -2.06556928        3.67155076        0.85800254
  6         1.27194185        2.63512097        1.77639245
  7         1.68293813        3.39721017        2.55178154
  6        -2.01177054       -2.40142658       -1.40228337
  7        -3.02783743       -2.70677712       -1.87711497
  6         0.88088317       -0.37111719        5.98047260
  6         1.00076864        0.91829638        5.38982315
 16        -0.54893043       -0.97946332        6.80148702
 16        -0.29870905        2.12037122        5.35880444
  6        -1.82302936       -0.16532709        5.87521697
  6        -1.72118571        1.06456940        5.30545689
  6         2.08127792       -1.04408549        5.87337761
  7         2.16029021        1.23724155        4.86862862
 16         3.23364718       -0.01648585        5.08539878
  6        -3.03291776       -0.91557985        5.80903156
  7        -4.00476312       -1.55328864        5.80293908
  6        -2.82565844        1.65075817        4.61955263
  7        -3.70759612        2.17881245        4.07708536
  6         2.40509005       -2.34311035        6.33473274
  7         2.69108102       -3.40366795        6.71489411

@aug-cc-pVTZ.gbs/N
