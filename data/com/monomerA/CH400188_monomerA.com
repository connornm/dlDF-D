%Chk=CH400188
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.95739189        0.00000000        0.55622543
  1        -0.70871265        0.69881334        0.48515178
  1         0.17967881        0.32220997       -1.04397424
  1        -0.42835804       -1.02102331        0.00259703
  6-Bq        0.00000000       0.00000000       3.73630032
  1-Bq       -0.51768245      -0.97852596       3.75815841
  1-Bq       -0.61678074       0.73148675       3.17907989
  1-Bq        0.15326846       0.35968462       4.77221620
  1-Bq        0.98119474      -0.11264541       3.23574678

@aug-cc-pVTZ.gbs/N

