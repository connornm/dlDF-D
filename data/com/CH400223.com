%Chk=CH400223
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.62009024        0.00000000        0.91731895
  1         0.65482841        0.09153907       -0.88814768
  1        -0.57178211       -0.94634681       -0.05898219
  1        -0.70313655        0.85480773        0.02981093
  0        -0.41032281        0.00000000       -0.60700340
  0        -0.43330957       -0.06057275        0.58770035
  0         0.37835661        0.62621156        0.03902938
  0         0.46527577       -0.56563882       -0.01972633
  6         0.00000000        0.00000000        4.60545746
  1         0.32045528       -1.02427257        4.33313400
  1         0.83266672        0.70731366        4.42557009
  1        -0.87132982        0.29074727        3.98720482
  1        -0.28179218        0.02621164        5.67592095
  0        -0.21204996        0.67777618        4.78565790
  0        -0.55098778       -0.46803983        4.72449158
  0         0.57657172       -0.19239173        5.01456431
  0         0.18646602       -0.01734463        3.89711607

@aug-cc-pVTZ.gbs/N
