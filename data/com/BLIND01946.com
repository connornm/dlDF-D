%Chk=BLIND01946
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.55483263        0.85011477        0.10283257
  6         1.27202357       -0.23466689       -0.47571091
 16        -0.75469544        0.69305087        1.26429019
 16         0.93477723       -1.93967017       -0.13834303
  6        -1.48202703       -0.83162760        0.72534085
  6        -0.80913543       -1.87432846        0.17069153
  6         1.09100388        2.04542882       -0.33196385
  7         2.26323281        0.05926147       -1.28129788
 16         2.43074247        1.71245581       -1.38032260
  6        -2.88423525       -0.90009758        0.97072552
  7        -4.02105451       -0.91858724        1.21257901
  6        -1.48442851       -3.07883391       -0.18580224
  7        -1.98790996       -4.08870176       -0.46436911
  6         0.69157794        3.36187090        0.00360674
  7         0.38326729        4.44920001        0.27518704
  6         0.55483267       -0.85011473        7.26311266
  6         1.27202356        0.23466690        6.68456907
 16        -0.75469541       -0.69305077        8.42457026
 16         0.93477715        1.93967020        7.02193676
  6        -1.48202707        0.83162762        7.88562075
  6        -0.80913551        1.87432844        7.33097133
  6         1.09100398       -2.04542881        6.82831638
  7         2.26323282       -0.05926151        5.87898214
 16         2.43074255       -1.71245585        5.77995760
  6        -2.88423529        0.90009756        8.13100541
  7        -4.02105456        0.91858720        8.37285889
  6        -1.48442864        3.07883383        6.97447742
  7        -1.98791014        4.08870161        6.69591044
  6         0.69157809       -3.36187087        7.16388711
  7         0.38326749       -4.44919996        7.43546753

@aug-cc-pVTZ.gbs/N

