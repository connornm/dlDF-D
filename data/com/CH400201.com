%Chk=CH400201
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10070398        0.00000000        0.12015284
  1        -0.46636595        0.49961934        0.87113103
  1        -0.27021625        0.54398461       -0.92576988
  1        -0.36412178       -1.04360395       -0.06551399
  0        -0.72835197        0.00000000       -0.07950690
  0         0.30860119       -0.33060545       -0.57644018
  0         0.17880606       -0.35996260        0.61259551
  0         0.24094472        0.69056804        0.04335157
  6         0.00000000        0.00000000        3.23151100
  1        -0.40026830       -0.58033911        2.37770906
  1         1.09356534        0.12576412        3.11199453
  1        -0.48488325        0.99498491        3.26116275
  1        -0.20841379       -0.54040991        4.17517767
  0         0.26486341        0.38401890        3.79648427
  0        -0.72362823       -0.08321996        3.31059680
  0         0.32085436       -0.65839611        3.21189001
  0         0.13791046        0.35759717        2.60707293

@aug-cc-pVTZ.gbs/N

