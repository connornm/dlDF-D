%Chk=CH400225
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.07623748        0.00000000        0.26019000
  1        -0.11708670       -0.17941996       -1.08631724
  1        -0.51608853       -0.80089674        0.56409484
  1        -0.44306224        0.98031670        0.26203239
  0        -0.71216213        0.00000000       -0.17217154
  0         0.07747799        0.11872482        0.71883205
  0         0.34150335        0.52996512       -0.37326983
  0         0.29318079       -0.64868994       -0.17339068
  6         0.00000000        0.00000000        3.94609327
  1        -0.51499694       -0.07719847        2.96895228
  1         0.27331603       -1.01283424        4.30028404
  1         0.91630246        0.61123373        3.83309291
  1        -0.67462156        0.47879898        4.68204385
  0         0.34078103        0.05108336        4.59268180
  0        -0.18085723        0.67020727        3.71172004
  0        -0.60633078       -0.40446233        4.02086727
  0         0.44640698       -0.31682831        3.45910398

@aug-cc-pVTZ.gbs/N

