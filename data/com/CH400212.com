%Chk=CH400212
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.36400776        0.00000000       -1.04569800
  1        -0.57923815       -0.92449110        0.18917011
  1         0.86374594        0.04235137        0.69147326
  1        -0.64851555        0.88213973        0.16505463
  0        -0.24086927        0.00000000        0.69195371
  0         0.38329038        0.61174932       -0.12517664
  0        -0.57155336       -0.02802452       -0.45755800
  0         0.42913225       -0.58372480       -0.10921907
  6         0.00000000        0.00000000        2.92115074
  1         0.87995541        0.60091101        2.62020158
  1        -0.78353598        0.66978721        3.32543088
  1        -0.39588013       -0.54153537        2.04023967
  1         0.29946070       -0.72916284        3.69873084
  0        -0.58227940       -0.39763163        3.12029321
  0         0.51847725       -0.44320802        2.65363265
  0         0.26195969        0.35834190        3.50406252
  0        -0.19815754        0.48249774        2.40661459

@aug-cc-pVTZ.gbs/N

