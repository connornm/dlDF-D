%Chk=BLIND00830
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.40919481        0.91293547       -0.20054261
  6         0.52321666        0.47561634       -1.18298217
 16        -1.42702817       -0.14151546        0.76933631
 16         0.83474607       -1.22406930       -1.56847728
  6        -0.36435374       -1.55121127        0.93547735
  6         0.53060664       -1.97872774        0.00600458
  6        -0.44434262        2.29259073       -0.17565853
  7         1.16112672        1.39785170       -1.86151936
 16         0.64747513        2.90460008       -1.37491042
  6        -0.56172626       -2.25984558        2.15612110
  7        -0.77522691       -2.81731216        3.15353852
  6         1.30607978       -3.15647439        0.21930289
  7         1.93706815       -4.12535003        0.33787693
  6        -1.23677719        3.13727337        0.63903599
  7        -1.88618808        3.84557467        1.29300133
  6         0.91903138       -0.43638583        4.68111011
  6         0.96570785        0.98316389        4.77344943
 16        -0.31034185       -1.45837689        5.41096789
 16        -0.24031121        1.96263207        5.62241850
  6        -1.72046874       -0.38812729        5.30982426
  6        -1.68957719        0.96810717        5.39536871
  6         2.02447593       -0.89326038        3.99222608
  7         1.98717948        1.59005872        4.22015530
 16         3.01413083        0.46119673        3.55535577
  6        -2.95022497       -1.09471771        5.17011294
  7        -3.92891511       -1.71552056        5.08088041
  6        -2.89100457        1.73517457        5.34930566
  7        -3.84599257        2.39767863        5.34578416
  6         2.37808207       -2.22940717        3.68426036
  7         2.68858753       -3.32048986        3.43084202

@aug-cc-pVTZ.gbs/N
