%Chk=BLIND00574
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.89354087       -0.48217666        0.10099152
  6         0.21891201       -1.23416899        0.57299025
 16        -1.00929926        1.27136741        0.11512734
 16         1.69710817       -0.51385439        1.22924577
  6         0.69074211        1.70837302       -0.13588563
  6         1.76177413        0.99900960        0.30834178
  6        -1.89612258       -1.34056076       -0.30312465
  7         0.11850129       -2.54071154        0.54607807
 16        -1.37926577       -2.97364228       -0.03717291
  6         0.85459716        2.94062559       -0.83287543
  7         0.93728559        3.96133001       -1.38289191
  6         3.09247417        1.46465531        0.09304335
  7         4.19173880        1.82018575       -0.03400592
  6        -3.17227451       -1.01795739       -0.82517041
  7        -4.22719845       -0.77146851       -1.24644778
  6        -0.26196334       -0.86107568        6.93509858
  6        -1.17518230       -0.47542271        5.91376143
 16         1.46403638       -0.53034682        6.92735015
 16        -0.71885583        0.43905797        4.46805068
  6         1.48993101        1.02588287        6.07787933
  6         0.62367856        1.40682734        5.10210678
  6        -0.92497885       -1.59596837        7.89720376
  7        -2.42277362       -0.85480064        6.04585352
 16        -2.58769369       -1.75770506        7.43455143
  6         2.56571282        1.86748487        6.48468301
  7         3.47152989        2.50775495        6.83216349
  6         0.76146639        2.66632607        4.44741647
  7         0.85310374        3.66936163        3.86737476
  6        -0.39840361       -2.18166086        9.07400080
  7         0.01854356       -2.67687007       10.03941201

@blind-aug.gbs/N
