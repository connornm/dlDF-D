%Chk=CH400309
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.83830288        0.00000000       -0.72334932
  1        -0.92595052       -0.33225311       -0.50814315
  1         0.23180095       -0.69092063        0.83359640
  1        -0.14415331        1.02317374        0.39789607
  6-Bq        0.00000000       0.00000000       3.52300472
  1-Bq        1.00408821      -0.40742996       3.29542253
  1-Bq       -0.75679902      -0.51336858       2.89875005
  1-Bq       -0.23119477      -0.16456733       4.59326288
  1-Bq       -0.01609442       1.08536587       3.30458343

@aug-cc-pVTZ.gbs/N

