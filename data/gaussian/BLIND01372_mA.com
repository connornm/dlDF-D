%Chk=BLIND01372_mA
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.43478913       -0.92202980        0.04391318
  6         0.47915897       -0.20934139        1.27514558
 16         0.45126669       -0.18699494       -1.55232266
 16         0.52646974        1.55618463        1.39870892
  6        -0.41984516        1.32153539       -1.22122979
  6        -0.38692901        2.01114285       -0.05033179
  6         0.44909107       -2.27966802        0.29274473
  7         0.52062969       -0.91483288        2.37896710
 16         0.54190804       -2.53663871        2.00452090
  6        -1.15557169        1.80067613       -2.34382512
  7        -1.72001988        2.16345615       -3.29302580
  6        -1.08904883        3.24374097        0.09688092
  7        -1.62120507        4.26533162        0.25208533
  6         0.42836307       -3.34471127       -0.64021918
  7         0.42122178       -4.23116949       -1.39204146

@blind-aug.gbs/N

