%Chk=CH400232
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.66146646        0.00000000       -0.88794599
  1         0.33977253       -0.77568180        0.71334307
  1         0.03809422        0.99287171        0.48861084
  1        -1.03933321       -0.21718991       -0.31400792
  6         0.00000000        0.00000000        2.92049077
  1        -0.71315827       -0.13590834        3.75650522
  1         0.12069468       -0.95815350        2.37886302
  1        -0.38786596        0.77022712        2.22601533
  1         0.98032955        0.32383472        3.32057953

@aug-cc-pVTZ.gbs/N

