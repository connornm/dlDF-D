%Chk=CH400148
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10153618        0.00000000        0.11226751
  1        -0.41417385       -0.93538292        0.42368003
  1        -0.42581676        0.86908661        0.53791683
  1        -0.26154557        0.06629631       -1.07386436
  6-Bq        0.00000000       0.00000000       2.34027292
  1-Bq        0.19094686      -0.27095209       3.39673410
  1-Bq        0.52525293       0.94486537       2.10084705
  1-Bq        0.37190885      -0.80673466       1.67932743
  1-Bq       -1.08810864       0.13282138       2.18418311

@aug-cc-pVTZ.gbs/N

