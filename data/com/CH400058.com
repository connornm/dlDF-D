%Chk=CH400058
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.60212919        0.00000000        0.92920740
  1        -1.06567305       -0.16566477        0.25076303
  1         0.11137076        0.97543552       -0.51196490
  1         0.35217310       -0.80977075       -0.66800554
  0        -0.39843772        0.00000000       -0.61487017
  0         0.70517149        0.10962281       -0.16593358
  0        -0.07369566       -0.64545999        0.33877468
  0        -0.23303810        0.53583718        0.44202906
  6         0.00000000        0.00000000        3.43282373
  1        -0.77173228        0.47796379        2.79881596
  1        -0.24828184        0.15968264        4.49998986
  1         0.03187583       -1.08609224        3.21981311
  1         0.98813829        0.44844582        3.21267600
  0         0.51066657       -0.31627565        3.85235597
  0         0.16429174       -0.10566435        2.72666425
  0        -0.02109271        0.71868317        3.57377598
  0        -0.65386560       -0.29674317        3.57849872

@aug-cc-pVTZ.gbs/N

