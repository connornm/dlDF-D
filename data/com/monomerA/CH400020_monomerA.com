%Chk=CH400020
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.02503758        0.00000000       -0.41866921
  1        -0.59884790        0.79195514       -0.49007572
  1        -0.47242886       -0.98498458       -0.18056106
  1         0.04623917        0.19302944        1.08930598
  6-Bq        0.00000000       0.00000000       4.48959875
  1-Bq       -0.95955291       0.07822063       5.03652261
  1-Bq        0.48388623      -0.96656282       4.72959122
  1-Bq       -0.19032019       0.05706009       3.40032912
  1-Bq        0.66598686       0.83128210       4.79195206

@aug-cc-pVTZ.gbs/N
