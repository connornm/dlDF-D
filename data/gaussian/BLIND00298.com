%Chk=BLIND00298
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.74033139       -0.69898944        0.06657405
  6        -1.16174509        0.61553009        0.41338170
 16         0.51926223       -1.08427877       -1.09682299
 16        -0.45386854        2.09369561       -0.25633902
  6         1.60191548        0.30156923       -0.87008663
  6         1.21280295        1.56093663       -0.53807876
  6        -1.51947100       -1.62605168        0.72899408
  7        -2.15879250        0.72613888        1.25689944
 16        -2.69987131       -0.78916696        1.68351262
  6         2.97014655       -0.01098041       -1.11776812
  7         4.06902437       -0.30800403       -1.35321984
  6         2.16302231        2.61855065       -0.42637510
  7         2.89789756        3.51604094       -0.35291996
  6        -1.45014261       -3.03894883        0.66509053
  7        -1.41363502       -4.19977431        0.61768461
  6         0.65914154       -0.66227558        4.28416011
  6         1.27408239        0.31904808        3.45670723
 16        -0.80558321       -0.41968946        5.22450581
 16         0.63662261        1.95277216        3.21343989
  6        -1.67714904        0.71419728        4.17627115
  6        -1.10238814        1.65478476        3.38091642
  6         1.42654222       -1.80959248        4.29211922
  7         2.39950612       -0.00700371        2.86936433
 16         2.83026923       -1.55400358        3.30771515
  6        -3.09352981        0.58709579        4.27008151
  7        -4.24204525        0.46196648        4.39833119
  6        -1.89896719        2.55200967        2.60994869
  7        -2.50825712        3.32199923        1.98800981
  6         1.19135395       -3.01975049        4.98880602
  7         1.01900150       -4.01687823        5.56075604

@blind-aug.gbs/N

