%Chk=CH400025
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.24993862        0.00000000       -1.07866429
  1        -0.92386592       -0.58762424        0.16478919
  1         0.83272639       -0.45341428        0.57181134
  1        -0.15879909        1.04103852        0.34206376
  0        -0.16538805        0.00000000        0.71376798
  0         0.61133563        0.38883958       -0.10904342
  0        -0.55102727        0.30003088       -0.37837595
  0         0.10507969       -0.68887046       -0.22634861
  6         0.00000000        0.00000000        4.44607019
  1         0.95881770        0.18254913        3.92326010
  1        -0.83905817        0.13662750        3.73663361
  1        -0.01445523       -1.03472254        4.83993092
  1        -0.10530430        0.71554591        5.28445615
  0        -0.63446373       -0.12079543        4.79202130
  0         0.55521710       -0.09040842        4.91551479
  0         0.00956524        0.68469108        4.18544677
  0         0.06968140       -0.47348723        3.89129791

@aug-cc-pVTZ.gbs/N

