%Chk=BLIND02432
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=130)

M05 opt

0 1
  6-Bq        0.47464498      -0.62699057      -0.65015681
  6-Bq        0.40522671      -1.12804184       0.68023989
 16-Bq        0.55905413       1.07277736      -1.08854104
 16-Bq        0.35995532      -0.10534793       2.12475919
  6-Bq       -0.40278479       1.79020340       0.21694271
  6-Bq       -0.47861719       1.32139542       1.49056454
  6-Bq        0.53775794      -1.68026304      -1.54006528
  7-Bq        0.40898281      -2.42950572       0.83528077
 16-Bq        0.53195812      -3.16847302      -0.65123236
  6-Bq       -1.08651022       2.97637857      -0.17878152
  7-Bq       -1.60515014       3.95135536      -0.54158437
  6-Bq       -1.24693174       2.00274929       2.48020670
  7-Bq       -1.83625674       2.54081850       3.32525327
  6-Bq        0.62832186      -1.63364542      -2.95243684
  7-Bq        0.71218225      -1.61414372      -4.11161026
  6         0.35782689        0.73471792       10.21728683
  6         0.35801508       -0.55698434       10.81505929
 16         0.55844958        1.05595510        8.50116508
 16         0.53063587       -2.06910527        9.91047060
  6        -0.21802838       -0.38128327        7.81156379
  6        -0.22616266       -1.61951172        8.37224925
  6         0.24791900        1.69990938       11.19791363
  7         0.25668574       -0.61722767       12.12041868
 16         0.18335294        0.92345422       12.74638874
  6        -0.82253874       -0.13687637        6.54431396
  7        -1.27964440        0.10531638        5.50335761
  6        -0.84114203       -2.72253504        7.70961673
  7        -1.30404071       -3.65506356        7.19285615
  6         0.22029779        3.10758829       11.04676607
  7         0.20537835        4.26514120       10.94214969

@aug-cc-pVTZ.gbs/N

