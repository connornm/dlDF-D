%Chk=CH400044
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.58450085        0.00000000       -0.94039604
  1        -0.95812430       -0.53110915       -0.16095610
  1         0.57745689       -0.51275531        0.79348061
  1        -0.20383344        1.04386446        0.30787153
  6-Bq        0.00000000       0.00000000       3.02908106
  1-Bq        0.67912006       0.85736344       2.85671564
  1-Bq       -0.67487559      -0.11697827       2.15911268
  1-Bq        0.59703460      -0.92327796       3.15982302
  1-Bq       -0.60127907       0.18289279       3.94067289

@aug-cc-pVTZ.gbs/N

