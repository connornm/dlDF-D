%Chk=CH400128
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.08069800        0.00000000        0.24099334
  1        -0.46307294        0.93086517        0.38084044
  1        -0.13335189       -0.05623643       -1.09774346
  1        -0.48427317       -0.87462874        0.47590967
  0        -0.71511372        0.00000000       -0.15946883
  0         0.30642216       -0.61596714       -0.25200771
  0         0.08824090        0.03721247        0.72639295
  0         0.32045066        0.57875467       -0.31491641
  6         0.00000000        0.00000000        2.76267208
  1         0.17442768        0.91058656        3.36797619
  1         0.84804361       -0.69943203        2.89534776
  1        -0.93686530       -0.48897118        3.09309084
  1        -0.08560599        0.27781665        1.69427352
  0        -0.11542135       -0.60254848        2.36213347
  0        -0.56116289        0.46282443        2.67487863
  0         0.61993751        0.32355940        2.54402914
  0         0.05664674       -0.18383535        3.46964707

@aug-cc-pVTZ.gbs/N
