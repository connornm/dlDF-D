%Chk=CH400250
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.35493009        0.00000000        1.04881390
  1        -1.04134577       -0.37443047       -0.03723924
  1         0.03605291        1.03112005       -0.40184274
  1         0.65036278       -0.65668958       -0.60973191
  6-Bq        0.00000000       0.00000000       4.13335209
  1-Bq        0.83263713       0.14733837       3.41851450
  1-Bq        0.40368677      -0.35031928       5.10304243
  1-Bq       -0.53404308       0.95869962       4.28058992
  1-Bq       -0.70228082      -0.75571871       3.73126151

@aug-cc-pVTZ.gbs/N
