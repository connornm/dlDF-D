%Chk=CH400004
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.01228281       0.00000000      -0.44863067
  1-Bq        0.08545224       0.02191379       1.10372264
  1-Bq       -0.55655696       0.89290361      -0.34489630
  1-Bq       -0.54117809      -0.91481740      -0.31019567
  6         0.00000000        0.00000000        3.22877258
  1        -0.47434437        0.91327020        3.63733215
  1         0.60483604       -0.48654842        4.01834839
  1        -0.78472644       -0.69957606        2.88123309
  1         0.65423477        0.27285428        2.37817670

@aug-cc-pVTZ.gbs/N

