%Chk=CH400284
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10076337        0.00000000        0.11960756
  1        -0.39191342        1.01795743        0.19013808
  1        -0.25919446       -0.30861412       -1.03129119
  1        -0.44965549       -0.70934331        0.72154554
  0        -0.72839127        0.00000000       -0.07914608
  0         0.25933486       -0.67359736       -0.12581716
  0         0.17151277        0.20421449        0.68242051
  0         0.29754363        0.46938287       -0.47745727
  6         0.00000000        0.00000000        4.30356437
  1         0.78214920       -0.03841622        3.52077910
  1        -0.70448989        0.82477154        4.08123159
  1        -0.55047850       -0.96051649        4.32273320
  1         0.47281919        0.17416117        5.28951358
  0        -0.51755960        0.02542057        4.82154486
  0         0.46617129       -0.54576342        4.45068523
  0         0.36425970        0.63558785        4.29088007
  0        -0.31287139       -0.11524500        3.65114731

@aug-cc-pVTZ.gbs/N

