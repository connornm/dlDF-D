%Chk=CH400272
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.10461530       0.00000000       0.07622995
  1-Bq       -0.36046122      -1.03784081      -0.13762312
  1-Bq       -0.43395618       0.42150985       0.92736047
  1-Bq       -0.31019789       0.61633096      -0.86596730
  6         0.00000000        0.00000000        3.36928136
  1        -0.33600022        1.05501324        3.37535389
  1        -0.67142841       -0.60487116        4.00904752
  1        -0.02707729       -0.39001338        2.33335552
  1         1.03450592       -0.06012871        3.75936852

@aug-cc-pVTZ.gbs/N

