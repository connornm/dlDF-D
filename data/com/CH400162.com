%Chk=CH400162
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.92566775        0.00000000       -0.60755672
  1        -0.03526024       -0.91743856        0.61890965
  1        -0.00923867        0.89005886        0.65855586
  1        -0.88116884        0.02737970       -0.66990879
  0        -0.61252793        0.00000000        0.40202920
  0         0.02333222        0.60708255       -0.40954159
  0         0.00611337       -0.58896501       -0.43577607
  0         0.58308235       -0.01811755        0.44328847
  6         0.00000000        0.00000000        3.52017323
  1         0.81189980       -0.45694685        4.11850779
  1        -0.95099097       -0.05650693        4.08445073
  1         0.24477351        1.06010634        3.31463387
  1        -0.10568234       -0.54665257        2.56310052
  0        -0.53724601        0.30236843        3.12424647
  0         0.62928467        0.03739146        3.14678253
  0        -0.16197022       -0.70148792        3.65618163
  0         0.06993156        0.36172802        4.15348228

@aug-cc-pVTZ.gbs/N

