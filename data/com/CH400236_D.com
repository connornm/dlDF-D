%Chk=CH400236
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.04730889        0.00000000        0.35934670
  1        -0.68222554       -0.19023205        0.85109696
  1        -0.23600860        0.98403829       -0.44939352
  1        -0.12907474       -0.79380624       -0.76105014
  6         0.00000000        0.00000000        3.47043804
  1         0.59442091        0.05536436        2.53792287
  1        -0.26757162       -1.05467640        3.67549672
  1         0.59693775        0.39983586        4.31292398
  1        -0.92378703        0.59947619        3.35540857

@aug-cc-pVTZ.gbs/N

