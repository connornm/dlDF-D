%Chk=BLIND01368
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.80952466        0.39837257        0.47652704
  6         0.22416613        1.27095415       -0.48356849
 16         0.81569513       -1.35630728        0.37866398
 16        -0.63072851        0.72118201       -1.93322159
  6        -0.73255176       -1.63873927       -0.43821167
  6        -1.30286112       -0.81359886       -1.35556362
  6         1.43347975        1.14072429        1.45881884
  7         0.36131686        2.55990216       -0.28949961
 16         1.25892681        2.82482171        1.08710160
  6        -1.33917194       -2.87615798       -0.07503606
  7        -1.78273228       -3.90512640        0.23427164
  6        -2.53247675       -1.15833269       -1.99050135
  7        -3.51702039       -1.41214139       -2.55368968
  6         2.14582659        0.67767002        2.59165209
  7         2.74177611        0.31537364        3.52153457
  6         0.28542467       -0.92187462        6.16056810
  6        -1.04690816       -0.57230063        5.80207976
 16         1.60851997        0.21546358        6.37127480
 16        -1.60038074        1.08866985        5.53857739
  6         0.72715070        1.62348732        6.99158073
  6        -0.54488510        1.96782156        6.65840974
  6         0.39756941       -2.29488555        6.24756106
  7        -1.90055963       -1.55122407        5.62644349
 16        -1.13210111       -3.00994601        5.85561310
  6         1.50674783        2.42787532        7.87264707
  7         2.18903145        3.06415910        8.56606256
  6        -1.14664320        3.15112418        7.17958538
  7        -1.66600984        4.12015991        7.55661878
  6         1.54285573       -3.06510039        6.56419455
  7         2.47497103       -3.71340102        6.81311882

@blind-aug.gbs/N
