%Chk=CH400051
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.07864807        0.00000000       -0.25000857
  1        -0.41216568        1.01757712       -0.14367398
  1        -0.53222152       -0.71059682       -0.66164823
  1        -0.13426088       -0.30698030        1.05533078
  0        -0.71375725        0.00000000        0.16543434
  0         0.27273607       -0.67334571        0.09507118
  0         0.35217878        0.47021234        0.43782234
  0         0.08884239        0.20313337       -0.69832786
  6         0.00000000        0.00000000        3.91311578
  1         1.03879605       -0.06665505        3.53569582
  1        -0.12403665       -0.68363190        4.77523427
  1        -0.21062756        1.03899105        4.23267732
  1        -0.70413184       -0.28870410        3.10885572
  0        -0.68738658        0.04410663        4.16286011
  0         0.08207687        0.45236926        3.34263933
  0         0.13937535       -0.68751562        3.70165722
  0         0.46593437        0.19103974        4.44530646

@aug-cc-pVTZ.gbs/N

