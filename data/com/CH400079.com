%Chk=CH400079
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.95870082        0.00000000       -0.55396633
  1        -0.30920164        1.04371264        0.20259376
  1        -0.77697212       -0.50391427       -0.60693549
  1         0.12747294       -0.53979837        0.95830805
  0        -0.63438640        0.00000000        0.36656764
  0         0.20460326       -0.69063996       -0.13405926
  0         0.51413385        0.33344746        0.40161812
  0        -0.08435071        0.35719250       -0.63412650
  6         0.00000000        0.00000000        2.95138869
  1        -0.44695432        0.03869039        3.96367365
  1         0.57109570        0.92996386        2.76429812
  1         0.67989425       -0.87109106        2.88117381
  1        -0.80403563       -0.09756319        2.19640916
  0         0.29575623       -0.02560200        2.28154488
  0        -0.37790240       -0.61537073        3.07518926
  0        -0.44989600        0.57641373        2.99785090
  0         0.53204217        0.06455899        3.45096972

@aug-cc-pVTZ.gbs/N
