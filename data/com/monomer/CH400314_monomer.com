%Chk=CH400314
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.98671828        0.00000000        0.50236739
  1        -0.51760870        0.95748853        0.20318194
  1         0.14166610       -0.11855577       -1.09172396
  1        -0.61077568       -0.83893276        0.38617463
  6-Bq        0.00000000       0.00000000       3.33838631
  1-Bq        0.01761834       0.00022787       4.44548861
  1-Bq       -1.03221078      -0.19014848       2.98572439
  1-Bq        0.34270651       0.98390815       2.96360242
  1-Bq        0.67188592      -0.79398754       2.95872982

@aug-cc-pVTZ.gbs/N

