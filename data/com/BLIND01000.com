%Chk=BLIND01000
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.43058731        0.92331149       -0.05656296
  6         0.20023668        0.32076429       -1.32530333
 16         0.77424739        0.05196738        1.43043866
 16         0.19113940       -1.42716036       -1.60684370
  6        -0.17204667       -1.42389014        1.16493122
  6        -0.40036973       -2.00872918       -0.04067067
  6         0.41393386        2.29750100       -0.18551679
  7         0.01747990        1.11988197       -2.34798164
 16         0.14460607        2.70291651       -1.84922571
  6        -0.65983534       -2.00069805        2.37341975
  7        -1.01527341       -2.44609298        3.38657332
  6        -1.13790257       -3.22535918       -0.13965889
  7        -1.70761884       -4.23071128       -0.26462554
  6         0.60970911        3.27711379        0.81807638
  7         0.77727877        4.09464232        1.62718912
  6        -0.30542865        0.95752062        7.58017161
  6        -1.26389014        0.17050921        6.88175104
 16         0.85514510        0.33984322        8.74638764
 16        -1.42195740       -1.58447994        7.05384444
  6         1.16234790       -1.27130864        8.07263826
  6         0.25682996       -2.03248127        7.40308802
  6        -0.47326600        2.28875417        7.25621010
  7        -2.09794595        0.79885169        6.08966730
 16        -1.80217589        2.43574944        6.15293432
  6         2.47924862       -1.74701514        8.33827988
  7         3.55501214       -2.10115062        8.60003154
  6         0.59300882       -3.33943154        6.94177825
  7         0.81368000       -4.41968096        6.57395826
  6         0.25467209        3.40729836        7.72963749
  7         0.83543142        4.33765386        8.11467406

@aug-cc-pVTZ.gbs/N
