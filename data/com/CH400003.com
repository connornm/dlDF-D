%Chk=CH400003
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.39834640        0.00000000        1.03310508
  1         0.21647652        0.97449899       -0.47903610
  1         0.48002184       -0.81142191       -0.58065435
  1        -1.09484476       -0.16307708        0.02658537
  0        -0.26359165        0.00000000       -0.68362078
  0        -0.14324569       -0.64484028        0.31698521
  0        -0.31763749        0.53692978        0.38422750
  0         0.72447484        0.10791049       -0.01759193
  6         0.00000000        0.00000000        4.93095320
  1         0.94982542       -0.48154919        4.62775688
  1        -0.08969738        0.98503980        4.43332188
  1        -0.85108888       -0.64270791        4.63336252
  1        -0.00903916        0.13921730        6.02937151
  0        -0.62851341        0.31864816        5.13158265
  0         0.05935407       -0.65181528        5.26024314
  0         0.56317800        0.42528925        5.12787332
  0         0.00598134       -0.09212213        4.20411369

@aug-cc-pVTZ.gbs/N

