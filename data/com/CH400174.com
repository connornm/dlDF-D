%Chk=CH400174
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.93652307        0.00000000        0.59068647
  1        -0.70194729       -0.74561741        0.42108269
  1        -0.46176545        1.00555271        0.04027858
  1         0.22718968       -0.25993531       -1.05204774
  0        -0.61971105        0.00000000       -0.39086590
  0         0.46448882        0.49338598       -0.27863660
  0         0.30555697       -0.66538898       -0.02665293
  0        -0.15033474        0.17200301        0.69615542
  6         0.00000000        0.00000000        3.56053849
  1         0.02702347       -0.97806650        3.04223201
  1         1.02651107        0.40909990        3.63052598
  1        -0.63907306        0.70210133        2.99077951
  1        -0.41446148       -0.13313473        4.57861648
  0        -0.01788183        0.64720095        3.90350949
  0        -0.67925743       -0.27070741        3.51422675
  0         0.42288402       -0.46459075        3.93755637
  0         0.27425524        0.08809720        2.88686136

@aug-cc-pVTZ.gbs/N

