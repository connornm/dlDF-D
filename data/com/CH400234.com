%Chk=CH400234
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.66644902        0.00000000        0.88421245
  1        -0.63129945        0.90953874        0.01364707
  1         0.61144704       -0.01105952       -0.92303638
  1        -0.64659661       -0.89847922        0.02517686
  0        -0.44099909        0.00000000       -0.58509635
  0         0.41774010       -0.60185513       -0.00903047
  0        -0.40460347        0.00731825        0.61078671
  0         0.42786246        0.59453688       -0.01665990
  6         0.00000000        0.00000000        3.52450473
  1         0.31799236       -0.43678137        4.49098745
  1         0.89063338        0.34599044        2.96499713
  1        -0.53391090       -0.76702041        2.93069937
  1        -0.67471485        0.85781134        3.71133497
  0        -0.21042021        0.28902464        2.88496894
  0        -0.58934517       -0.22894695        3.89473911
  0         0.35329667        0.50754866        3.91743445
  0         0.44646871       -0.56762634        3.40087642

@aug-cc-pVTZ.gbs/N

