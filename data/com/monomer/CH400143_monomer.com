%Chk=CH400143
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10720435        0.00000000        0.00919129
  1        -0.37757241        0.20050177        1.02138341
  1        -0.36337458        0.78697696       -0.68892099
  1        -0.36625736       -0.98747873       -0.34165370
  6-Bq        0.00000000       0.00000000       3.41398235
  1-Bq       -1.04892582       0.31154026       3.58334364
  1-Bq        0.42185484      -0.40299454       4.35505612
  1-Bq        0.03084127      -0.78222336       2.63093612
  1-Bq        0.59622972       0.87367765       3.08659351

@aug-cc-pVTZ.gbs/N

