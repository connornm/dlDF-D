%Chk=CH400172
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10716755        0.00000000       -0.01288294
  1        -0.37133009        1.02545576       -0.19115545
  1        -0.37825157       -0.68200380       -0.78599143
  1        -0.35758589       -0.34345196        0.99002982
  0        -0.73262901        0.00000000        0.00852483
  0         0.24571457       -0.67855912        0.12649037
  0         0.25029461        0.45129192        0.52010206
  0         0.23661983        0.22726720       -0.65511725
  6         0.00000000        0.00000000        3.50691691
  1         0.28218046        1.03808078        3.24471767
  1         0.34519264       -0.22809480        4.53395184
  1        -1.10079063       -0.10752841        3.45511357
  1         0.47341753       -0.70245757        2.79388454
  0        -0.18672295       -0.68691328        3.68041798
  0        -0.22841903        0.15093368        2.82731282
  0         0.72840931        0.07115313        3.54119593
  0        -0.31326732        0.46482648        3.97874088

@aug-cc-pVTZ.gbs/N

