%Chk=CH400166
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.02871597        0.00000000       -1.10687007
  1         0.66713370        0.79468641        0.38651274
  1        -1.03591196        0.18889830        0.34232994
  1         0.34006229       -0.98358471        0.37802739
  0        -0.01900178        0.00000000        0.73243216
  0        -0.44145215       -0.52585565       -0.25576115
  0         0.68547814       -0.12499678       -0.22652474
  0        -0.22502421        0.65085243       -0.25014627
  6         0.00000000        0.00000000        3.51177705
  1        -0.05103791        1.10166004        3.41315534
  1         0.84447090       -0.38583731        2.90846583
  1         0.15289428       -0.26868812        4.57498728
  1        -0.94632727       -0.44713461        3.15049976
  0         0.03377253       -0.72898461        3.57703649
  0        -0.55879878        0.25531421        3.91099694
  0        -0.10117239        0.17779487        2.80823525
  0         0.62619864        0.29587553        3.75083954

@aug-cc-pVTZ.gbs/N

