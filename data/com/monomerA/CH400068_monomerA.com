%Chk=CH400068
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.09613907        0.00000000        0.15641323
  1        -0.26873680       -0.78849876       -0.72940811
  1        -0.51016449       -0.19822530        0.96250968
  1        -0.31723778        0.98672406       -0.38951479
  6-Bq        0.00000000       0.00000000       4.58025839
  1-Bq       -0.09787400      -0.19394049       3.49453572
  1-Bq        0.40201314       1.01950652       4.73830221
  1-Bq       -0.99385450      -0.08349451       5.06116155
  1-Bq        0.68971537      -0.74207153       5.02703408

@aug-cc-pVTZ.gbs/N

