%Chk=CH400305
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.05229732        0.00000000        0.34446525
  1        -0.43696670       -1.00647430        0.14851108
  1        -0.57883178        0.74319774        0.58189075
  1        -0.03649883        0.26327656       -1.07486708
  6-Bq        0.00000000       0.00000000       3.58786249
  1-Bq       -0.17908915      -0.17484338       2.50927877
  1-Bq       -0.94282629       0.32171627       4.07114624
  1-Bq        0.35564518      -0.93650960       4.05951090
  1-Bq        0.76627026       0.78963671       3.71151402

@aug-cc-pVTZ.gbs/N

