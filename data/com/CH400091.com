%Chk=CH400091
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.03453630        0.00000000       -0.39461450
  1        -0.27509495        1.02540831        0.31439868
  1        -0.69620905       -0.34321279       -0.78960997
  1        -0.06323230       -0.68219552        0.86982579
  0        -0.68456784        0.00000000        0.26112221
  0         0.18203436       -0.67852772       -0.20804222
  0         0.46069174        0.22710894        0.52249650
  0         0.04184174        0.45141878       -0.57557649
  6         0.00000000        0.00000000        2.67543638
  1         0.07808400        0.55349819        3.63122328
  1        -0.94474945       -0.57721280        2.65931779
  1         0.00880665        0.71713452        1.83185681
  1         0.85785881       -0.69341991        2.57934763
  0        -0.05166933       -0.36625787        2.04297817
  0         0.62515457        0.38195018        2.68610228
  0        -0.00582749       -0.47453843        3.23364535
  0        -0.56765776        0.45884612        2.73901971

@aug-cc-pVTZ.gbs/N

