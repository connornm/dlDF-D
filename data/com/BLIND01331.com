%Chk=BLIND01331
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.59391460        0.03128008        0.82909343
  6         0.48309039        1.28113833        0.15723110
 16         0.50926843       -1.54461663        0.05585269
 16         0.20667247        1.45116330       -1.58324957
  6        -0.62771809       -1.19718684       -1.25969901
  6        -0.74438425       -0.00850790       -1.90860227
  6         0.84929052        0.24544978        2.16857754
  7         0.62938981        2.36717729        0.87624146
 16         0.95227648        1.95187100        2.45553934
  6        -1.41513605       -2.32548046       -1.63120424
  7        -2.01657441       -3.27911898       -1.91411498
  6        -1.66139405        0.15237997       -2.98894946
  7        -2.37316726        0.32016424       -3.89236159
  6         1.04122329       -0.71478130        3.19145903
  7         1.21041031       -1.48905728        4.04173450
  6         0.50330842       -0.72622464        7.47335632
  6         0.37193912        0.49357194        8.19492755
 16        -0.64924165       -1.33382554        6.29396405
 16        -0.97499451        1.62503263        7.99454081
  6        -1.24441779        0.19025321        5.61053784
  6        -1.37368804        1.36163339        6.28788689
  6         1.63652772       -1.39364987        7.89247538
  7         1.29409671        0.77028184        9.08426272
 16         2.40370212       -0.46921079        9.14214830
  6        -1.63673420        0.07571349        4.24523269
  7        -1.96566727       -0.07273681        3.14026275
  6        -1.90943919        2.52176053        5.65462734
  7        -2.36742517        3.48268829        5.18794456
  6         2.13580330       -2.64243019        7.44925980
  7         2.55499578       -3.67016495        7.10410063

@aug-cc-pVTZ.gbs/N

