%Chk=CH400295
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.98291812        0.00000000        0.50976263
  1        -0.21578208        1.01525099       -0.38560274
  1         0.02122154       -0.71803709       -0.84259025
  1        -0.78835757       -0.29721390        0.71843037
  6-Bq        0.00000000       0.00000000       3.37693789
  1-Bq        0.22215934       0.90625269       2.78083619
  1-Bq        0.77577725      -0.12835714       4.15647625
  1-Bq       -0.99167719       0.10820295       3.85741101
  1-Bq       -0.00625940      -0.88609851       2.71302809

@aug-cc-pVTZ.gbs/N

