%Chk=CH400303
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.59431975        0.00000000        0.93422160
  1         0.26659399       -0.88680528       -0.60703377
  1        -1.07844352       -0.03357353        0.24863298
  1         0.21752978        0.92037881       -0.57582080
  6-Bq        0.00000000       0.00000000       3.26055484
  1-Bq       -1.08018875       0.23872989       3.21379731
  1-Bq        0.15474612      -1.06368164       2.99480880
  1-Bq        0.55124199       0.64265477       2.54703325
  1-Bq        0.37420064       0.18229698       4.28658001

@aug-cc-pVTZ.gbs/N
