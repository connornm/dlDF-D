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
  0        -0.72890265        0.00000000       -0.07428906
  0         0.27406491        0.61895660       -0.28035529
  0         0.28176919       -0.57508736       -0.35594745
  0         0.17306854       -0.04386924        0.71059180
  6         0.00000000        0.00000000        2.34027292
  1         0.19094686       -0.27095209        3.39673410
  1         0.52525293        0.94486537        2.10084705
  1         0.37190885       -0.80673466        1.67932743
  1        -1.08810864        0.13282138        2.18418311
  0        -0.12635234        0.17929297        1.64119707
  0        -0.34756757       -0.62523128        2.49870453
  0        -0.24609754        0.53382816        2.77763025
  0         0.72001745       -0.08788985        2.44355984

@aug-cc-pVTZ.gbs/N
