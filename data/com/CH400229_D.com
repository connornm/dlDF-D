%Chk=CH400229
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.04796886        0.00000000       -0.35741745
  1        -0.67528115        0.26476325       -0.83658935
  1        -0.26035912       -1.00688107        0.37998629
  1        -0.11232858        0.74211782        0.81402051
  6         0.00000000        0.00000000        4.54037808
  1         0.59369572        0.71505436        5.14221578
  1        -0.58635479       -0.65293681        5.21554291
  1        -0.68862830        0.55847046        3.87713773
  1         0.68128736       -0.62058801        3.92661591

@aug-cc-pVTZ.gbs/N

