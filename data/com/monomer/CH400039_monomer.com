%Chk=CH400039
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.92617759        0.00000000       -0.60677923
  1        -0.21487509       -1.02977479        0.34551199
  1         0.13306957        0.66320085        0.87660885
  1        -0.84437208        0.36657394       -0.61534161
  6-Bq        0.00000000       0.00000000       2.47871427
  1-Bq       -0.24388608       0.93796438       3.01418449
  1-Bq        0.02997765      -0.83961801       3.19991389
  1-Bq       -0.77378140      -0.19930041       1.71221075
  1-Bq        0.98768983       0.10095404       1.98854793

@aug-cc-pVTZ.gbs/N
