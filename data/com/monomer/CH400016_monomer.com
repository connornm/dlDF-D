%Chk=CH400016
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.03577273        0.00000000       -1.10666448
  1         0.77809496       -0.68189648        0.39442539
  1         0.18329755        1.02548229        0.37519867
  1        -0.99716523       -0.34358581        0.33704042
  6-Bq        0.00000000       0.00000000       3.52643674
  1-Bq        0.33043460       0.84027532       2.88554960
  1-Bq       -0.83011973       0.33652669       4.17731074
  1-Bq       -0.34661344      -0.83629577       2.88889893
  1-Bq        0.84629857      -0.34050624       4.15398768

@aug-cc-pVTZ.gbs/N

