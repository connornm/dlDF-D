%Chk=CH400285
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.93026595        0.00000000        0.60049248
  1        -0.74694716       -0.66400592        0.47660469
  1        -0.40352495        1.02960319       -0.05541529
  1         0.22020616       -0.36559726       -1.02168187
  0        -0.61557063        0.00000000       -0.39735468
  0         0.49426589        0.43938246       -0.31537631
  0         0.26701838       -0.68130353        0.03666911
  0        -0.14571364        0.24192107        0.67606188
  6         0.00000000        0.00000000        2.97400030
  1        -0.56791431        0.70600810        2.33759488
  1        -0.15432295       -1.03335980        2.60747893
  1         1.07787953        0.24966560        2.93123623
  1        -0.35564227        0.07768610        4.01969117
  0         0.37579723       -0.46717591        3.39511910
  0         0.10211776        0.68378933        3.21653287
  0        -0.71324870       -0.16520739        3.00229792
  0         0.23533371       -0.05140603        2.28205132

@aug-cc-pVTZ.gbs/N

