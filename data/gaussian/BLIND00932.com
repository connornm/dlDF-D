%Chk=BLIND00932
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.49108915        0.50250431       -0.73988497
  6         0.62452419       -0.91209826       -0.82307939
 16        -0.90278478        1.35115278       -0.08763666
 16        -0.59421609       -2.05949893       -0.24624996
  6        -1.41825227        0.20794379        1.16595361
  6        -1.29572096       -1.14417135        1.09950991
  6         1.59107752        1.10700673       -1.31425909
  7         1.70618203       -1.38309292       -1.39407885
 16         2.65637858       -0.12038170       -1.91715033
  6        -2.03997097        0.84251883        2.28034583
  7        -2.56089716        1.40301167        3.15537538
  6        -1.78763820       -1.97920531        2.14577835
  7        -2.19953376       -2.69902955        2.96010092
  6         1.86021832        2.48998632       -1.45603053
  7         2.09305718        3.62097760       -1.58922476
  6        -0.53631481        0.56256618        5.76078190
  6        -0.82744870       -0.82913491        5.82586484
 16         1.01735898        1.25191309        5.31422258
 16         0.33435174       -2.10556621        5.43149038
  6         1.58392496        0.04124607        4.14914300
  6         1.31244666       -1.28971032        4.19912974
  6        -1.63766057        1.28676901        6.17050251
  7        -2.02274576       -1.17647954        6.23632395
 16        -2.90050959        0.18601015        6.61588867
  6         2.42113443        0.58991088        3.13467915
  7         3.11725634        1.07922212        2.34273104
  6         1.85849531       -2.18807290        3.23542985
  7         2.30475350       -2.95980213        2.48963245
  6        -1.77814739        2.69228309        6.27031263
  7        -1.90803385        3.84330061        6.36718416

@blind-aug.gbs/N
