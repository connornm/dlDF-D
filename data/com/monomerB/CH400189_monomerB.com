%Chk=CH400189
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.09901779       0.00000000      -0.13470656
  1-Bq       -0.47652500      -0.51912585      -0.85405985
  1-Bq       -0.25655119      -0.52478727       0.94061988
  1-Bq       -0.36594160       1.04391312       0.04814653
  6         0.00000000        0.00000000        3.13043056
  1         0.58354932        0.06339628        4.06927935
  1        -0.73943751        0.82355282        3.09914340
  1         0.68424424        0.08522681        2.26409796
  1        -0.52835606       -0.97217591        3.08920152

@aug-cc-pVTZ.gbs/N
