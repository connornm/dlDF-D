%Chk=BLIND00670
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.45577703        0.23340731       -0.88255155
  6         0.21107586       -1.15874723       -0.71565719
 16        -0.53373484        1.52246141       -0.21343092
 16        -1.12842441       -1.82436008        0.23159509
  6        -1.07868122        0.76146767        1.29257176
  6        -1.31490741       -0.56591141        1.46550527
  6         1.54629560        0.41719524       -1.70856272
  7         1.01412833       -1.99032773       -1.33326600
 16         2.13753007       -1.13244625       -2.21242740
  6        -1.30293054        1.69783985        2.34329010
  7        -1.49446023        2.50130758        3.16118313
  6        -1.79835931       -1.07073018        2.70863607
  7        -2.21980571       -1.52188418        3.69340464
  6         2.12373177        1.63445037       -2.14448230
  7         2.60485525        2.62387265       -2.51952851
  6        -0.11413684        0.10284595        6.90633347
  6        -0.93774634        0.97918921        7.66759331
 16         1.46549606       -0.48695905        7.40171877
 16        -0.51201830        1.59861110        9.27060937
  6         1.23919834       -0.57744571        9.15805689
  6         0.45532099        0.25194131        9.89649951
  6        -0.71233781       -0.15622305        5.68964038
  7        -2.06401342        1.37185381        7.12427296
 16        -2.20450720        0.72124832        5.59853023
  6         2.01296209       -1.60624665        9.76956974
  7         2.68013003       -2.44206165       10.22508400
  6         0.37969570        0.12254414       11.31476575
  7         0.31046249        0.07247501       12.47391457
  6        -0.23438079       -0.96386684        4.62924562
  7         0.14679989       -1.61525822        3.74522838

@aug-cc-pVTZ.gbs/N

