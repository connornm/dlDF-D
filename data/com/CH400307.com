%Chk=CH400307
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.61056656        0.00000000       -0.92368525
  1         0.63868607       -0.26556594        0.86460441
  1        -0.43276639        1.00709986        0.15636203
  1        -0.81648624       -0.74153391       -0.09728120
  0        -0.40402085        0.00000000        0.61121608
  0        -0.42262794        0.17572888       -0.57212142
  0         0.28636787       -0.66641275       -0.10346705
  0         0.54028092        0.49068387        0.06437240
  6         0.00000000        0.00000000        3.84052424
  1         1.10108063       -0.11433090        3.81737589
  1        -0.36355252       -0.13680032        4.87739511
  1        -0.27356665        1.01194268        3.48398635
  1        -0.46396147       -0.76081146        3.18333962
  0        -0.72860120        0.07565443        3.85584185
  0         0.24056803        0.09052278        3.15441158
  0         0.18102306       -0.66961731        4.07645059
  0         0.30701011        0.50344010        4.27539295

@aug-cc-pVTZ.gbs/N
