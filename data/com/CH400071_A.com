%Chk=CH400071
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10723833        0.00000000        0.00303860
  1        -0.36676773       -0.61658458       -0.84338061
  1        -0.37170070       -0.42122232        0.95414690
  1        -0.36876991        1.03780690       -0.11380489
  6-Bq        0.00000000       0.00000000       3.01691381
  1-Bq        0.08023366      -0.63414598       3.92101966
  1-Bq        0.34725711      -0.57061815       2.13385471
  1-Bq       -1.05524667       0.30140580       2.86996154
  1-Bq        0.62775590       0.90335834       3.14281934

@aug-cc-pVTZ.gbs/N

