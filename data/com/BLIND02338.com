%Chk=BLIND02338
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.31822990       -0.01449276        0.96934442
  6         0.40613895        1.26756001        0.35742992
 16         0.47846086       -1.55002137        0.12970488
 16         0.66250271        1.52465701       -1.37541950
  6        -0.21454349       -1.12869955       -1.44711311
  6        -0.13871459        0.09118303       -2.04190765
  6         0.15954218        0.13185858        2.33268689
  7         0.32391440        2.31591995        1.13976015
 16         0.16162143        1.82127757        2.72086950
  6        -0.84780622       -2.23042613       -2.09230030
  7        -1.33120111       -3.16370106       -2.58880840
  6        -0.69113893        0.31275456       -3.33793054
  7        -1.10076357        0.53074798       -4.40358182
  6         0.04217891       -0.87875350        3.31760188
  7        -0.04632433       -1.69502504        4.14038001
  6         0.59375741        0.61267628        9.24723724
  6        -0.77077280        1.01583092        9.28398357
 16         1.43492073       -0.23171522       10.53873401
 16        -1.86291521        0.69926789       10.64103962
  6         0.10858616       -1.19982525       11.20785407
  6        -1.19802238       -0.82727457       11.24815487
  6         1.17499102        1.04841835        8.07352537
  7        -1.22756195        1.69174629        8.25814630
 16        -0.00366119        1.92657639        7.15458903
  6         0.54485132       -2.43645701       11.76597973
  7         0.95150636       -3.42196383       12.22923021
  6        -2.18242341       -1.66301009       11.85349997
  7        -3.01264929       -2.29798446       12.36188828
  6         2.51119117        0.87639989        7.63746660
  7         3.60629171        0.75268302        7.26794683

@aug-cc-pVTZ.gbs/N
