%Chk=CH400316
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.18050368        0.00000000       -1.09243049
  1        -0.95307044        0.52031109        0.21660808
  1        -0.05829136       -1.04391651        0.36445356
  1         0.83085812        0.52360542        0.51136886
  0        -0.11944193        0.00000000        0.72287728
  0         0.63066069       -0.34429748       -0.14333274
  0         0.03857225        0.69077487       -0.24116427
  0        -0.54979101       -0.34647738       -0.33838027
  6         0.00000000        0.00000000        4.32658920
  1        -0.07361656        0.18787973        3.23788917
  1         0.93799974       -0.54661546        4.54419892
  1         0.00128615        0.96556600        4.86849934
  1        -0.86566932       -0.60683028        4.65576939
  0         0.04871316       -0.12432277        5.04699799
  0        -0.62068819        0.36170347        4.18259366
  0        -0.00085106       -0.63892919        3.96799932
  0         0.57282609        0.40154850        4.10876585

@aug-cc-pVTZ.gbs/N

