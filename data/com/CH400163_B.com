%Chk=CH400163
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.03358982       0.00000000       0.39708695
  1-Bq       -0.48255917      -0.97037659       0.22691820
  1-Bq       -0.57689490       0.81850634       0.47246756
  1-Bq        0.02586425       0.15187025      -1.09647272
  6         0.00000000        0.00000000        2.85416994
  1         0.80241792        0.31874886        2.16098122
  1         0.36509137       -0.83576059        3.48202221
  1        -0.28715498        0.85049437        3.50238850
  1        -0.88035432       -0.33348263        2.27128782

@aug-cc-pVTZ.gbs/N

