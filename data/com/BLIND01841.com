%Chk=BLIND01841
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.39384662        0.06381202       -0.93910690
  6        -0.43021852       -1.24942915       -0.39147439
 16        -0.49074173        1.55158349       -0.00869622
 16        -0.54829795       -1.59971271        1.33998366
  6         0.32562471        1.05173975        1.48389352
  6         0.29980106       -0.19843383        2.01685058
  6        -0.34302729       -0.00900456       -2.31657500
  7        -0.40766801       -2.25427253       -1.23268856
 16        -0.37191859       -1.67543943       -2.79300619
  6         1.00536739        2.12178552        2.13507828
  7         1.52434896        3.03053273        2.64107202
  6         0.95338914       -0.48483379        3.25165405
  7         1.44641649       -0.75637030        4.26857002
  6        -0.30614239        1.05315442       -3.25232403
  7        -0.28477560        1.91246530       -4.03477159
  6        -0.39175032        0.06274654        8.55994497
  6        -0.42637447       -1.25048889        9.10770453
 16        -0.49432803        1.55043598        9.48987753
 16        -0.54777299       -1.60077179       10.83893318
  6         0.31963250        1.05272981       10.98449418
  6         0.29544486       -0.19741597       11.51759317
  6        -0.33747250       -0.01016994        7.18261410
  7        -0.39947750       -2.25541041        8.26671148
 16        -0.36134894       -1.67674303        6.70638870
  6         0.99532690        2.12445759       11.63712446
  7         1.51098306        3.03449016       12.14420677
  6         0.94674694       -0.48209858       12.75400060
  7         1.43797491       -0.75232593       13.77213537
  6        -0.30082440        1.05192365        6.24678155
  7        -0.27958846        1.91115788        5.46424625

@aug-cc-pVTZ.gbs/N
