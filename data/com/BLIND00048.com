%Chk=BLIND00048
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.69816384        0.73288840        0.12864978
  6         0.22604640        1.00088067       -0.92008385
 16        -1.22470691       -0.86723001        0.62953911
 16         1.02852599       -0.25187645       -1.87988302
  6         0.24537513       -1.80522250        0.30836784
  6         1.13625983       -1.55946881       -0.68841379
  6        -1.17220010        1.92213786        0.64467163
  7         0.46524818        2.25732003       -1.20644493
 16        -0.45589842        3.23837557       -0.22668176
  6         0.40479558       -2.91444989        1.18887302
  7         0.47893251       -3.81937437        1.91468870
  6         2.26579710       -2.40644742       -0.89019269
  7         3.18297627       -3.08587545       -1.10947803
  6        -2.11161087        2.11609350        1.68633236
  7        -2.88795458        2.29386556        2.53296300
  6         0.89766580       -0.08657910        5.17756677
  6        -0.08520111       -1.00843728        4.71933453
 16         0.56320675        1.51898999        5.80907423
 16        -1.82744138       -0.69429981        4.74180309
  6        -1.00577845        1.23897355        6.58624751
  6        -1.95225448        0.36149212        6.15980463
  6         2.15729879       -0.61690402        4.98421886
  7         0.33912226       -2.14406631        4.22118153
 16         2.00336252       -2.17427560        4.23874991
  6        -1.22630327        2.07542476        7.71885857
  7        -1.36994432        2.79482321        8.62053820
  6        -3.20470306        0.24683318        6.83223888
  7        -4.24915697        0.14098126        7.33107506
  6         3.41221995       -0.03445015        5.28602809
  7         4.45115344        0.43205144        5.51859933

@aug-cc-pVTZ.gbs/N

