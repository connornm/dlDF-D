%Chk=CH400225
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.07623748        0.00000000        0.26019000
  1        -0.11708670       -0.17941996       -1.08631724
  1        -0.51608853       -0.80089674        0.56409484
  1        -0.44306224        0.98031670        0.26203239
  6-Bq        0.00000000       0.00000000       3.94609327
  1-Bq       -0.51499694      -0.07719847       2.96895228
  1-Bq        0.27331603      -1.01283424       4.30028404
  1-Bq        0.91630246       0.61123373       3.83309291
  1-Bq       -0.67462156       0.47879898       4.68204385

@aug-cc-pVTZ.gbs/N
