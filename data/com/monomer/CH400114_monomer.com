%Chk=CH400114
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10668604        0.00000000        0.03509941
  1        -0.33854502       -0.41604855       -0.96864690
  1        -0.39549223       -0.62113293        0.82690128
  1        -0.37264878        1.03718148        0.10664621
  6-Bq        0.00000000       0.00000000       3.78749465
  1-Bq        1.03853090       0.20282684       3.46145662
  1-Bq       -0.32045204       0.77979396       4.50528141
  1-Bq       -0.67176405       0.01001722       2.90737042
  1-Bq       -0.04631481      -0.99263802       4.27587017

@aug-cc-pVTZ.gbs/N

