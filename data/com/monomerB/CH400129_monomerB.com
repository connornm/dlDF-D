%Chk=CH400129
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.08079055       0.00000000       0.24057794
  1-Bq       -0.38220463       1.03902255       0.01837725
  1-Bq       -0.15378310      -0.43205798      -1.00780088
  1-Bq       -0.54480282      -0.60696457       0.74884568
  6         0.00000000        0.00000000        2.76633606
  1         0.37591671        0.19764714        1.74378622
  1        -0.99516836       -0.48177461        2.70701396
  1         0.70411169       -0.67129737        3.29508225
  1        -0.08486004        0.95542484        3.31946183

@aug-cc-pVTZ.gbs/N

