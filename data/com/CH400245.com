%Chk=CH400245
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.96012322        0.00000000       -0.55149738
  1        -0.00065489       -0.82376355        0.73986427
  1        -0.12440250        0.96720527        0.52442725
  1        -0.83506583       -0.14344172       -0.71279414
  0        -0.63532762        0.00000000        0.36493391
  0         0.00043335        0.54509642       -0.48957904
  0         0.08231896       -0.64001391       -0.34702120
  0         0.55257531        0.09491749        0.47166634
  6         0.00000000        0.00000000        3.47518755
  1        -0.94303580       -0.28643280        2.97058307
  1        -0.20007524        0.19936476        4.54579931
  1         0.73272296       -0.82510499        3.38405703
  1         0.41038808        0.91217303        3.00031078
  0         0.62402062        0.18953679        3.80909174
  0         0.13239272       -0.13192259        2.76674804
  0        -0.48485353        0.54598408        3.53548995
  0        -0.27155981       -0.60359827        3.78942047

@aug-cc-pVTZ.gbs/N

