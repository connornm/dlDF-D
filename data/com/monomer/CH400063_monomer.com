%Chk=CH400063
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10655631        0.00000000        0.03897547
  1        -0.40367905        0.33298043        0.97578339
  1        -0.34128789        0.69034524       -0.79557023
  1        -0.36158938       -1.02332567       -0.21918863
  6-Bq        0.00000000       0.00000000       4.31077708
  1-Bq       -0.28221916      -0.63592217       5.17213722
  1-Bq        1.03431590      -0.24081852       3.99743537
  1-Bq       -0.69414125      -0.18929209       3.46915817
  1-Bq       -0.05795549       1.06603278       4.60437755

@aug-cc-pVTZ.gbs/N

