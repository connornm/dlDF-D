%Chk=CH400080
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.94483759        0.00000000        0.57729359
  1         0.20464277       -0.31084306       -1.04282495
  1        -0.71509450       -0.70762911        0.46247903
  1        -0.43438585        1.01847218        0.00305234
  6-Bq        0.00000000       0.00000000       3.44337496
  1-Bq       -0.56206063      -0.86656348       3.84230043
  1-Bq       -0.70634754       0.81629832       3.19696850
  1-Bq        0.72189467       0.35384934       4.20471896
  1-Bq        0.54651351      -0.30358418       2.52951196

@aug-cc-pVTZ.gbs/N

