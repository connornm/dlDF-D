%Chk=CH400258
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.92998760        0.00000000        0.60092347
  1        -0.49491368        0.98674853        0.08587050
  1         0.24624466       -0.19829890       -1.06114423
  1        -0.68131859       -0.78844963        0.37435026
  6-Bq        0.00000000       0.00000000       3.15771809
  1-Bq       -1.05803932       0.31363005       3.24813351
  1-Bq        0.40975194      -0.21093725       4.16449253
  1-Bq        0.06226029      -0.91347073       2.53507787
  1-Bq        0.58602709       0.81077793       2.68316844

@aug-cc-pVTZ.gbs/N

