%Chk=CH400123
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.03705287        0.00000000       -0.38795271
  1        -0.13755269        0.85843117        0.68568298
  1        -0.71022875        0.08522189       -0.84516171
  1        -0.18927143       -0.94365306        0.54743144
  0        -0.68623309        0.00000000        0.25671400
  0         0.09102064       -0.56803650       -0.45372648
  0         0.46996878       -0.05639258        0.55925590
  0         0.12524368        0.62442907       -0.36224341
  6         0.00000000        0.00000000        4.02475955
  1         0.47042037       -0.80128070        4.62695542
  1        -1.07860568        0.05942550        4.26779112
  1         0.48439478        0.96837685        4.25626318
  1         0.12379054       -0.22652166        2.94802848
  0        -0.31128406        0.53021919        3.62627771
  0         0.71372920       -0.03932273        3.86394199
  0        -0.32053113       -0.64078917        3.87157020
  0        -0.08191402        0.14989270        4.73724829

@aug-cc-pVTZ.gbs/N

