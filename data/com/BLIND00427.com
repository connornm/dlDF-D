%Chk=BLIND00427
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.06269036       -0.46696592        0.90505340
  6        -0.03597091        0.95327340        0.99468050
 16        -0.73782616       -1.38424766       -0.43334209
 16        -0.65775982        2.03859884       -0.25833011
  6        -0.38846422       -0.28882496       -1.78320039
  6        -0.35876842        1.06820643       -1.71089543
  6         0.45714884       -1.01477156        2.06046649
  7         0.45134001        1.47905950        2.09192116
 16         0.89765941        0.26663698        3.14156277
  6        -0.17525343       -0.97076462       -3.01632209
  7        -0.03238734       -1.56970119       -4.00225174
  6        -0.11340083        1.86099890       -2.87073023
  7         0.05345784        2.54670785       -3.79424935
  6         0.60679701       -2.38259139        2.39489554
  7         0.72850507       -3.50031668        2.68980780
  6         0.29131489        0.82325922        8.57734534
  6         1.31040237        0.29167343        7.73787890
 16        -1.43352859        0.74372724        8.25007440
 16         0.99883434       -0.57433532        6.22541649
  6        -1.52831500       -0.80807061        7.39739563
  6        -0.56255193       -1.32741180        6.59413324
  6         0.86635518        1.47296256        9.65089900
  7         2.55430939        0.49033757        8.09993956
 16         2.59125054        1.38568563        9.50275961
  6        -2.77050303       -1.47881081        7.59275715
  7        -3.80745851       -1.97584285        7.76268509
  6        -0.75887725       -2.56594096        5.91483667
  7        -0.88690369       -3.55595819        5.31952058
  6         0.22374720        2.14982803       10.71582931
  7        -0.28667606        2.71751804       11.59235202

@aug-cc-pVTZ.gbs/N
