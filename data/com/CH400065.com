%Chk=CH400065
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.39118746        0.00000000        1.03583702
  1        -0.76273771       -0.79553798       -0.10647289
  1        -0.45875056        0.98314363       -0.22127469
  1         0.83030081       -0.18760565       -0.70808944
  0        -0.25885448        0.00000000       -0.68542855
  0         0.50471473        0.52641915        0.07045467
  0         0.30356197       -0.65056056        0.14642071
  0        -0.54942222        0.12414141        0.46855317
  6         0.00000000        0.00000000        3.61838771
  1        -0.75061085        0.40574323        4.32403758
  1        -0.29746841        0.24695899        2.58083783
  1         0.98960851        0.44696045        3.83492247
  1         0.05847074       -1.09966267        3.73375296
  0         0.49669021       -0.26848625        3.15144884
  0         0.19683921       -0.16341639        4.30494969
  0        -0.65483847       -0.29576028        3.47510348
  0        -0.03869095        0.72766292        3.54204883

@aug-cc-pVTZ.gbs/N
