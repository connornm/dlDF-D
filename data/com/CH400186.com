%Chk=CH400186
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.09243154        0.00000000        0.18049732
  1        -0.48395938       -0.74131358        0.66499884
  1        -0.40889123        1.00718221        0.21066069
  1        -0.19958093       -0.26586863       -1.05615684
  0        -0.72287798        0.00000000       -0.11943772
  0         0.32024302        0.49053807       -0.44003949
  0         0.27056932       -0.66646724       -0.13939727
  0         0.13206563        0.17592917        0.69887448
  6         0.00000000        0.00000000        4.10029381
  1        -0.85057583        0.15711875        3.40904826
  1         0.95185801        0.12841153        3.54941927
  1        -0.05237055        0.73834449        4.92375644
  1        -0.04891163       -1.02387476        4.51895128
  0         0.56283850       -0.10396778        4.55770113
  0        -0.62985841       -0.08497179        4.46481558
  0         0.03465436       -0.48857338        3.55539651
  0         0.03236554        0.67751295        3.82326202

@aug-cc-pVTZ.gbs/N

