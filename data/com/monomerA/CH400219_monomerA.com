%Chk=CH400219
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.72581231        0.00000000        0.83617130
  1         0.07906155        0.95346134       -0.55735689
  1         0.22113452       -0.84484421       -0.68067888
  1        -1.02600839       -0.10861713        0.40186448
  6-Bq        0.00000000       0.00000000       4.04010872
  1-Bq        0.52389485       0.92338910       4.35454904
  1-Bq       -1.04939371       0.03306691       4.39176960
  1-Bq        0.01690351      -0.07604455       2.93560998
  1-Bq        0.50859535      -0.88041145       4.47850626

@aug-cc-pVTZ.gbs/N
