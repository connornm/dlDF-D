%Chk=CH400114
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10668604        0.00000000        0.03509941
  1        -0.33854502       -0.41604855       -0.96864690
  1        -0.39549223       -0.62113293        0.82690128
  1        -0.37264878        1.03718148        0.10664621
  0        -0.73231038        0.00000000       -0.02322580
  0         0.22402021        0.27530543        0.64096787
  0         0.26170301        0.41101277       -0.54717271
  0         0.24658717       -0.68631820       -0.07056936
  6         0.00000000        0.00000000        3.78749465
  1         1.03853090        0.20282684        3.46145662
  1        -0.32045204        0.77979396        4.50528141
  1        -0.67176405        0.01001722        2.90737042
  1        -0.04631481       -0.99263802        4.27587017
  0        -0.68721113       -0.13421349        4.00323881
  0         0.21204782       -0.51600110        3.31252463
  0         0.44451613       -0.00662854        4.36988577
  0         0.03064719        0.65684314        3.46432941

@aug-cc-pVTZ.gbs/N

