%Chk=CH400212
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.36400776        0.00000000       -1.04569800
  1        -0.57923815       -0.92449110        0.18917011
  1         0.86374594        0.04235137        0.69147326
  1        -0.64851555        0.88213973        0.16505463
  6-Bq        0.00000000       0.00000000       2.92115074
  1-Bq        0.87995541       0.60091101       2.62020158
  1-Bq       -0.78353598       0.66978721       3.32543088
  1-Bq       -0.39588013      -0.54153537       2.04023967
  1-Bq        0.29946070      -0.72916284       3.69873084

@aug-cc-pVTZ.gbs/N

