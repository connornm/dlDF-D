%Chk=CH400180
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.91072699       0.00000000       0.62973193
  1-Bq       -0.84919102      -0.41161577       0.57916580
  1-Bq       -0.23350620       1.03662273      -0.31124607
  1-Bq        0.17197022      -0.62500696      -0.89765166
  6         0.00000000        0.00000000        3.26120831
  1         0.90951255       -0.61078957        3.42154850
  1         0.27365715        1.07279599        3.24683891
  1        -0.46391073       -0.27554439        2.29433289
  1        -0.71925897       -0.18646203        4.08211293

@aug-cc-pVTZ.gbs/N

