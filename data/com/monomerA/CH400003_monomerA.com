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
  6-Bq        0.00000000       0.00000000       4.93095320
  1-Bq        0.94982542      -0.48154919       4.62775688
  1-Bq       -0.08969738       0.98503980       4.43332188
  1-Bq       -0.85108888      -0.64270791       4.63336252
  1-Bq       -0.00903916       0.13921730       6.02937151

@aug-cc-pVTZ.gbs/N

