%Chk=CH400311
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.05425266        0.00000000        0.33843359
  1        -0.54933685        0.81882534        0.50372616
  1        -0.03571114        0.15136159       -1.09626655
  1        -0.46920466       -0.97018693        0.25410680
  0        -0.69761445        0.00000000       -0.22394647
  0         0.36350425       -0.54182874       -0.33332299
  0         0.02363059       -0.10015819        0.72541565
  0         0.31047961        0.64198692       -0.16814620
  6         0.00000000        0.00000000        3.48452616
  1         0.01962716       -0.00498258        4.59158348
  1        -1.02281610        0.23954460        3.13459536
  1         0.70764776        0.76298368        3.10627161
  1         0.29554118       -0.99754569        3.10565421
  0        -0.01298758        0.00329705        2.75197010
  0         0.67681242       -0.15851017        3.71608051
  0        -0.46826091       -0.50487749        3.73482275
  0        -0.19556393        0.66009062        3.73523129

@aug-cc-pVTZ.gbs/N

