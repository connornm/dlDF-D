%Chk=CH400095
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.80708042        0.00000000       -0.75802846
  1         0.39005531        0.40366976        0.95440743
  1        -0.83789923        0.63189890       -0.35300796
  1        -0.35923650       -1.03556866        0.15662899
  0        -0.53405695        0.00000000        0.50159855
  0        -0.25810532       -0.26711420       -0.63154539
  0         0.55445021       -0.41813677        0.23359054
  0         0.23771206        0.68525097       -0.10364370
  6         0.00000000        0.00000000        3.36834524
  1         0.97230855        0.41431709        3.03827951
  1        -0.65132282        0.82306236        3.72095393
  1        -0.48723947       -0.51868053        2.52007985
  1         0.16625374       -0.71869892        4.19406765
  0        -0.64339083       -0.27415970        3.58675457
  0         0.43098986       -0.54463244        3.13501889
  0         0.32241351        0.34321851        3.92965488
  0        -0.11001254        0.47557362        2.82195260

@aug-cc-pVTZ.gbs/N

