%Chk=CH400124
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.08349983        0.00000000       -0.22806594
  1        -0.57450689        0.13033053       -0.93751893
  1        -0.27774495       -0.96215157        0.47234315
  1        -0.23124800        0.83182104        0.69324172
  6-Bq        0.00000000       0.00000000       3.10232523
  1-Bq       -0.87822323      -0.51440061       2.66631843
  1-Bq        0.72256732      -0.75291197       3.47246844
  1-Bq        0.48350886       0.62581598       2.32736749
  1-Bq       -0.32785295       0.64149659       3.94314657

@aug-cc-pVTZ.gbs/N
