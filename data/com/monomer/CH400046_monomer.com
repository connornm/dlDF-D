%Chk=CH400046
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.21455012        0.00000000        1.08625697
  1        -0.61792414        0.88292519       -0.25416296
  1         0.95183054        0.04088199       -0.56421027
  1        -0.54845652       -0.92380718       -0.26788373
  6-Bq        0.00000000       0.00000000       4.41004696
  1-Bq       -0.16233991       1.08129297       4.23558426
  1-Bq       -0.34783397      -0.26607812       5.42700316
  1-Bq       -0.56867472      -0.58419773       3.66084202
  1-Bq        1.07884859      -0.23101711       4.31675838

@aug-cc-pVTZ.gbs/N

