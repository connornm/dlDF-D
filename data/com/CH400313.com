%Chk=CH400313
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.04161889        0.00000000        0.37552104
  1        -0.54821601        0.85935280        0.43238632
  1         0.00570106        0.08360573       -1.10406682
  1        -0.49910394       -0.94295853        0.29615945
  0        -0.68925450        0.00000000       -0.24848778
  0         0.36276257       -0.56864635       -0.28611637
  0        -0.00377248       -0.05532314        0.73057721
  0         0.33026440        0.62396949       -0.19597305
  6         0.00000000        0.00000000        3.28312395
  1         0.10129457       -0.02473323        4.38544588
  1        -1.05953566        0.17543683        3.01369480
  1         0.62516600        0.81680916        2.87327492
  1         0.33307509       -0.96751275        2.86008020
  0        -0.06702810        0.01636634        2.55370136
  0         0.70111029       -0.11608913        3.46140917
  0        -0.41368151       -0.54049460        3.55432707
  0        -0.22040067        0.64021738        3.56305821

@aug-cc-pVTZ.gbs/N

