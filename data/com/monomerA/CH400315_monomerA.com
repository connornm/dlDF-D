%Chk=CH400315
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.01422374        0.00000000       -1.10715114
  1         0.65265321       -0.81088007        0.37749603
  1         0.36874632        0.97480732        0.37384863
  1        -1.03562327       -0.16392725        0.35580648
  6-Bq        0.00000000       0.00000000       3.93531221
  1-Bq        0.04178892      -0.19483643       2.84614812
  1-Bq       -0.82849371       0.70149098       4.15324459
  1-Bq       -0.17102673      -0.95198588       4.47425441
  1-Bq        0.95773152       0.44533132       4.26760171

@aug-cc-pVTZ.gbs/N

