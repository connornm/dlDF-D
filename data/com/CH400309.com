%Chk=CH400309
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.83830288        0.00000000       -0.72334932
  1        -0.92595052       -0.33225311       -0.50814315
  1         0.23180095       -0.69092063        0.83359640
  1        -0.14415331        1.02317374        0.39789607
  0        -0.55471731        0.00000000        0.47865085
  0         0.61271504        0.21985676        0.33624578
  0        -0.15338609        0.45719231       -0.55160297
  0         0.09538836       -0.67704907       -0.26329366
  6         0.00000000        0.00000000        3.52300472
  1         1.00408821       -0.40742996        3.29542253
  1        -0.75679902       -0.51336858        2.89875005
  1        -0.23119477       -0.16456733        4.59326288
  1        -0.01609442        1.08536587        3.30458343
  0        -0.66441990        0.26960238        3.67359920
  0         0.50078501        0.33970352        3.93608320
  0         0.15298497        0.10889661        2.81479920
  0         0.01064991       -0.71820251        3.66753730

@aug-cc-pVTZ.gbs/N
