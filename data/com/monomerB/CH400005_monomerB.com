%Chk=CH400005
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.58726436       0.00000000      -0.93867275
  1-Bq       -0.67548450       0.87723653       0.01275632
  1-Bq       -0.59993877      -0.92868563       0.06002020
  1-Bq        0.68815890       0.05144911       0.86589622
  6         0.00000000        0.00000000        4.49799906
  1        -0.49550759        0.97934823        4.64406489
  1         0.66220115        0.04990780        3.61200537
  1         0.60102024       -0.24690596        5.39454692
  1        -0.76771380       -0.78235008        4.34137904

@aug-cc-pVTZ.gbs/N

