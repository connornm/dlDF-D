%Chk=CH400234
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.66644902        0.00000000        0.88421245
  1        -0.63129945        0.90953874        0.01364707
  1         0.61144704       -0.01105952       -0.92303638
  1        -0.64659661       -0.89847922        0.02517686
  6-Bq        0.00000000       0.00000000       3.52450473
  1-Bq        0.31799236      -0.43678137       4.49098745
  1-Bq        0.89063338       0.34599044       2.96499713
  1-Bq       -0.53391090      -0.76702041       2.93069937
  1-Bq       -0.67471485       0.85781134       3.71133497

@aug-cc-pVTZ.gbs/N

