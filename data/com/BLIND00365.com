%Chk=BLIND00365
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.23395777       -0.66950346       -0.73357875
  6        -0.84644842       -1.08343762        0.09540113
 16         0.63319449        0.99810571       -1.11855447
 16        -1.95686668        0.03045254        0.90830230
  6         0.14083437        1.80206219        0.38323032
  6        -0.88760088        1.41666836        1.18412768
  6         0.87655143       -1.77821459       -1.24654684
  7        -1.04231725       -2.37219913        0.23174334
 16         0.07305788       -3.20597403       -0.68021862
  6         0.91399924        2.96297192        0.67579540
  7         1.55581940        3.91458793        0.85933526
  6        -1.23165171        2.16348721        2.34938841
  7        -1.56546249        2.75728112        3.29114490
  6         1.98927939       -1.82313535       -2.12118222
  7         2.89617243       -1.87876069       -2.84612825
  6        -0.14743425       -0.68235347        7.50730412
  6         0.79986341        0.32760414        7.17800271
 16        -0.45119287       -1.29629336        9.12573705
 16         1.81495625        1.16115040        8.36517392
  6        -0.14940024        0.16356775       10.08586665
  6         0.75092436        1.13562455        9.78227937
  6        -0.73996315       -1.15657227        6.35433870
  7         0.94499837        0.62783087        5.91046186
 16        -0.05492198       -0.33412987        4.99083818
  6        -0.92775057        0.21969591       11.27838920
  7        -1.56607217        0.20918426       12.24979441
  6         0.95036459        2.24984916       10.64984438
  7         1.16653638        3.15911039       11.34081465
  6        -1.72382425       -2.16646975        6.22283641
  7        -2.52448528       -2.99941961        6.09543012

@aug-cc-pVTZ.gbs/N

