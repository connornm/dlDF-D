%Chk=BLIND02552
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.50349864       -0.39458000       -0.79492431
  6         0.42043146       -1.27226415        0.32248207
 16         0.55291319        1.35971944       -0.70268374
 16         0.32004055       -0.73014594        2.00492608
  6        -0.45293368        1.63526714        0.73135518
  6        -0.54210680        0.80522468        1.80407216
  6         0.61143201       -1.13163384       -1.95699369
  7         0.45455272       -2.56012021        0.08139589
 16         0.62614979       -2.81756733       -1.55420562
  6        -1.15889290        2.87261004        0.69019958
  7        -1.69508123        3.90169036        0.62172482
  6        -1.34729685        1.14454305        2.93122520
  7        -1.96702042        1.39380467        3.88240208
  6         0.72895516       -0.66250083       -3.28788473
  7         0.83541904       -0.29518749       -4.38553780
  6        -0.50349858        0.39458006        9.82660211
  6        -0.42043149        1.27226412        8.70919565
 16        -0.55291313       -1.35971938        9.73436170
 16        -0.32004069        0.73014577        7.02675168
  6         0.45293364       -1.63526720        8.30032273
  6         0.54210668       -0.80522483        7.22760568
  6        -0.61143189        1.13163401       10.98867143
  7        -0.45455274        2.56012020        8.95028172
 16        -0.62614970        2.81756746       10.58588322
  6         1.15889287       -2.87261009        8.34147839
  7         1.69508121       -3.90169040        8.40995321
  6         1.34729665       -1.14454329        6.10045262
  7         1.96702016       -1.39380499        5.14927571
  6        -0.72895494        0.66250111       12.31956253
  7        -0.83541874        0.29518787       13.41721564

@aug-cc-pVTZ.gbs/N

