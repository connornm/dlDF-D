%Chk=CH400009
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.91672574        0.00000000        0.62096688
  1        -0.62160832       -0.87875716        0.25956677
  1        -0.57435948        0.92739842        0.18981389
  1         0.27924206       -0.04864126       -1.07034754
  0        -0.60661087        0.00000000       -0.41090289
  0         0.41132735        0.58148650       -0.17175914
  0         0.38006209       -0.61367314       -0.12560263
  0        -0.18477857        0.03218664        0.70826467
  6         0.00000000        0.00000000        3.26714319
  1         1.00164228        0.43796881        3.09141036
  1        -0.41127266        0.38927752        4.21861771
  1        -0.67462344        0.27526218        2.43341667
  1         0.08425382       -1.10250851        3.32512802
  0        -0.66280139       -0.28981038        3.38342818
  0         0.27214515       -0.25759065        2.63753854
  0         0.44640822       -0.18214502        3.81883226
  0        -0.05575199        0.72954605        3.22877378

@aug-cc-pVTZ.gbs/N
