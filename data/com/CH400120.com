%Chk=CH400120
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.64836161        0.00000000        0.89755957
  1        -0.67645636        0.87594581        0.03334191
  1         0.62898123        0.05382296       -0.90965469
  1        -0.60088648       -0.92976877       -0.02124679
  0        -0.42903039        0.00000000       -0.59392833
  0         0.44762110       -0.57962619       -0.02206283
  0        -0.41620611       -0.03561544        0.60193185
  0         0.39761540        0.61524163        0.01405931
  6         0.00000000        0.00000000        3.62746818
  1        -0.93437421        0.58317261        3.51415225
  1        -0.08536489       -0.94758307        3.06108922
  1         0.85219374        0.58855882        3.23586852
  1         0.16754535       -0.22414837        4.69876273
  0         0.61828912       -0.38589387        3.70245100
  0         0.05648720        0.62702961        4.00224944
  0        -0.56390910       -0.38945800        3.88659542
  0        -0.11086722        0.14832226        2.91857686

@aug-cc-pVTZ.gbs/N
