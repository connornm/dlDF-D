%Chk=BLIND00708
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
  6         0.37482805       -0.31829393        6.18494289
  6         0.67721030        1.04068884        6.48084399
 16         0.16364066       -1.59125700        7.37804035
 16         0.86006524        1.68230486        8.12080780
  6        -0.53714184       -0.66430629        8.71736419
  6        -0.25792595        0.63337807        9.01010615
  6         0.34333015       -0.49949271        4.81692312
  7         0.86235842        1.85149289        5.46784005
 16         0.71301607        1.00618952        4.04162970
  6        -1.42263317       -1.43338058        9.52700576
  7        -2.11746149       -2.10710794       10.17082102
  6        -0.84230451        1.27368829       10.14255767
  7        -1.27200426        1.82625684       11.07043280
  6         0.09266215       -1.69066091        4.09342307
  7        -0.10177615       -2.66085834        3.48349962

@aug-cc-pVTZ.gbs/N

