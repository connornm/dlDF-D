%Chk=CH400048
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.97296508        0.00000000       -0.52851198
  1        -0.60172460        0.86718704       -0.33451472
  1        -0.54409245       -0.93689656       -0.22841672
  1         0.17285197        0.06970952        1.09144342
  0        -0.64382526        0.00000000        0.34972413
  0         0.39816999       -0.57383038        0.22135330
  0         0.36003396        0.61995820        0.15114669
  0        -0.11437869       -0.04612781       -0.72222412
  6         0.00000000        0.00000000        4.16805521
  1         1.07962644       -0.00352902        4.41377914
  1        -0.53053784       -0.71556846        4.82568530
  1        -0.41025564        1.01683356        4.32208693
  1        -0.13883295       -0.29773609        3.11066948
  0        -0.71440465        0.00233520        4.00545609
  0         0.35106467        0.47350214        3.73289174
  0         0.27147218       -0.67285368        4.06613016
  0         0.09186780        0.19701634        4.86774286

@aug-cc-pVTZ.gbs/N

