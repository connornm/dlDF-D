%Chk=BLIND01670
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.62925455        0.67096806       -0.44153084
  6        -0.81237130       -0.66688144       -0.89147933
 16        -0.08617284        1.13336391        1.16463508
 16        -0.48779485       -2.09950995        0.09685969
  6         0.99490256       -0.22071413        1.54129235
  6         0.83218551       -1.50192492        1.11742192
  6        -1.01111599        1.55083393       -1.43416427
  7        -1.28237536       -0.83643626       -2.10317540
 16        -1.57742216        0.64618967       -2.80008617
  6         2.07266673        0.14566579        2.39877693
  7         2.92004890        0.48695529        3.11756615
  6         1.73672363       -2.52923703        1.51789979
  7         2.43134848       -3.40282890        1.84233730
  6        -1.00331958        2.96654732       -1.40478319
  7        -1.01243153        4.12883595       -1.39488298
  6         0.60241277       -0.68951904        8.86612130
  6        -0.43597434       -0.20044865        9.70777979
 16         0.49975476       -0.84621612        7.11871636
 16        -2.00669885        0.36593976        9.11903792
  6        -0.55718713        0.52454587        6.73399894
  6        -1.55115503        1.00201735        7.52879001
  6         1.67543625       -1.08883119        9.63714247
  7        -0.21713925       -0.20516182       11.00004486
 16         1.28722877       -0.84720824       11.30918908
  6        -0.30155705        1.08460844        5.44866909
  7        -0.07195975        1.49159831        4.38432077
  6        -2.37741932        2.08363268        7.10295255
  7        -3.09299443        2.94258970        6.78503551
  6         2.90642420       -1.64365858        9.21052359
  7         3.91875768       -2.10876101        8.87892659

@aug-cc-pVTZ.gbs/N

