%Chk=CH400060
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10515077        0.00000000       -0.06802748
  1        -0.40957618        0.80014938       -0.64652477
  1        -0.39036118       -0.98071632       -0.33436449
  1        -0.30521340        0.18056694        1.04891673
  6-Bq        0.00000000       0.00000000       3.12665697
  1-Bq        0.26223189      -0.90826006       2.55021612
  1-Bq        0.92323347       0.55768661       3.37687999
  1-Bq       -0.66641535       0.64285831       2.51952459
  1-Bq       -0.51905001      -0.29228486       4.06000717

@aug-cc-pVTZ.gbs/N

