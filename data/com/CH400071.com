%Chk=CH400071
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10723833        0.00000000        0.00303860
  1        -0.36676773       -0.61658458       -0.84338061
  1        -0.37170070       -0.42122232        0.95414690
  1        -0.36876991        1.03780690       -0.11380489
  0        -0.73267585        0.00000000       -0.00201069
  0         0.24269559        0.40800306        0.55807732
  0         0.24595980        0.27872899       -0.63137300
  0         0.24402046       -0.68673205        0.07530637
  6         0.00000000        0.00000000        3.01691381
  1         0.08023366       -0.63414598        3.92101966
  1         0.34725711       -0.57061815        2.13385471
  1        -1.05524667        0.30140580        2.86996154
  1         0.62775590        0.90335834        3.14281934
  0        -0.05309179        0.41962370        2.41865371
  0        -0.22978512        0.37758640        3.60124697
  0         0.69827220       -0.19944464        3.11415428
  0        -0.41539529       -0.59776546        2.93360028

@aug-cc-pVTZ.gbs/N

