%Chk=CH400151
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.17441882        0.00000000        1.09341851
  1        -1.02272172       -0.36833615       -0.21060561
  1         0.10914558        1.03008198       -0.39115766
  1         0.73915732       -0.66174583       -0.49165523
  0        -0.11541549        0.00000000       -0.72353107
  0         0.67674996        0.24373343        0.13936082
  0        -0.07222323       -0.68162035        0.25883476
  0        -0.48911124        0.43788692        0.32533548
  6         0.00000000        0.00000000        3.48411831
  1        -0.07896116       -0.48548914        4.47611196
  1        -0.25806573       -0.73123171        2.69374611
  1        -0.69896036        0.85740221        3.43614094
  1         1.03598725        0.35931863        3.33047425
  0         0.05224976        0.32125528        2.82770157
  0         0.17076588        0.48386675        4.00711920
  0         0.46251232       -0.56735562        3.51586564
  0        -0.68552796       -0.23776641        3.58578684

@aug-cc-pVTZ.gbs/N
