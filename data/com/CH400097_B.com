%Chk=CH400097
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.23422597       0.00000000      -1.08218490
  1-Bq        0.31367647      -0.96390170       0.44551829
  1-Bq        0.54192091       0.82907346       0.49491907
  1-Bq       -1.08982335       0.13482824       0.14174754
  6         0.00000000        0.00000000        3.17657252
  1        -0.92053153        0.37869113        2.69159759
  1        -0.03834665        0.21722190        4.26162097
  1         0.07565976       -1.09372641        3.02157542
  1         0.88321842        0.49781337        2.73149609

@aug-cc-pVTZ.gbs/N

