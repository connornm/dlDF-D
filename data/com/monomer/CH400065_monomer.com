%Chk=CH400065
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.39118746        0.00000000        1.03583702
  1        -0.76273771       -0.79553798       -0.10647289
  1        -0.45875056        0.98314363       -0.22127469
  1         0.83030081       -0.18760565       -0.70808944
  6-Bq        0.00000000       0.00000000       3.61838771
  1-Bq       -0.75061085       0.40574323       4.32403758
  1-Bq       -0.29746841       0.24695899       2.58083783
  1-Bq        0.98960851       0.44696045       3.83492247
  1-Bq        0.05847074      -1.09966267       3.73375296

@aug-cc-pVTZ.gbs/N

