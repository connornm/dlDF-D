%Chk=CH400268
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.01711325        0.00000000       -1.10711024
  1         1.03802073       -0.01195068        0.38517021
  1        -0.51721863        0.90997582        0.36112998
  1        -0.53791535       -0.89802513        0.36081005
  0        -0.01132409        0.00000000        0.73259109
  0        -0.68687355        0.00790794       -0.25487278
  0         0.34225115       -0.60214435       -0.23896500
  0         0.35594648        0.59423640       -0.23875331
  6         0.00000000        0.00000000        3.74850806
  1        -0.02526752        0.06868442        2.64368683
  1         0.19553339       -1.04765461        4.04879537
  1         0.80404673        0.65251713        4.14057463
  1        -0.97431260        0.32645306        4.16097541
  0         0.01671989       -0.04544949        4.47958448
  0        -0.12938731        0.69324842        3.54980355
  0        -0.53204952       -0.43178016        3.48907186
  0         0.64471694       -0.21601878        3.47557236

@aug-cc-pVTZ.gbs/N
