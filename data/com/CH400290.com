%Chk=CH400290
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.92018283        0.00000000        0.61583237
  1        -0.41060180        1.02707612       -0.05006761
  1         0.23992215       -0.35179766       -1.02208695
  1        -0.74950318       -0.67527846        0.45632219
  0        -0.60889848        0.00000000       -0.40750531
  0         0.27170123       -0.67963133        0.03313047
  0        -0.15876000        0.23278967        0.67632993
  0         0.49595725        0.44684166       -0.30195509
  6         0.00000000        0.00000000        3.76321352
  1         0.28555812       -1.03066341        4.04987614
  1         0.87294557        0.51408060        3.31636553
  1        -0.82420607       -0.03723619        3.02478020
  1        -0.33429761        0.55381900        4.66183223
  0        -0.18895799        0.68200510        3.57352466
  0        -0.57764089       -0.34017467        4.05889940
  0         0.54538925        0.02463973        4.25184569
  0         0.22120963       -0.36647016        3.16858434

@aug-cc-pVTZ.gbs/N
