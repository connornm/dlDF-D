%Chk=CH400107
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.89361543        0.00000000       -0.65378698
  1        -0.69303763        0.80117098       -0.32219537
  1        -0.50997324       -0.98016962       -0.07197753
  1         0.30939543        0.17899864        1.04795988
  0        -0.59131844        0.00000000        0.43262043
  0         0.45859316       -0.53014659        0.21320140
  0         0.33745677        0.64859262        0.04762859
  0        -0.20473150       -0.11844602       -0.69345042
  6         0.00000000        0.00000000        3.86537980
  1         0.47197905       -0.96588653        3.60026672
  1         0.53880444        0.82499919        3.36035300
  1        -1.05799099       -0.00426534        3.53884829
  1         0.04720749        0.14515267        4.96205117
  0        -0.31231546        0.63914128        4.04080901
  0        -0.35653480       -0.54591407        4.19956343
  0         0.70008816        0.00282244        4.08145048
  0        -0.03123789       -0.09604965        3.13969626

@aug-cc-pVTZ.gbs/N

