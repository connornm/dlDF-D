%Chk=BLIND00659
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.76124529       -0.15026839       -0.66259636
  6        -0.43926375       -1.30170676        0.10963692
 16        -0.58785771        1.51263633       -0.12111574
 16         0.22627813       -1.24272638        1.74919895
  6         0.82995569        1.36190011        0.93302174
  6         1.14935619        0.26834890        1.67448862
  6        -1.29717511       -0.53879113       -1.87384882
  7        -0.68584829       -2.47341496       -0.42342590
 16        -1.37337205       -2.26978342       -1.92567826
  6         1.62109465        2.54639964        0.97938637
  7         2.22045117        3.54209153        1.00310217
  6         2.29101159        0.26755408        2.52921914
  7         3.19381108        0.23263094        3.26042334
  6        -1.76137118        0.27895367       -2.93261409
  7        -2.15538976        0.93476150       -3.80769200
  6         0.51919095       -0.04589187        6.85207546
  6         0.96398917       -0.97982045        7.82970453
 16        -1.12891371        0.54237133        6.69016791
 16        -0.07795285       -1.68389307        9.07598978
  6        -1.64474117        0.52620350        8.38657079
  6        -1.22714870       -0.35881048        9.32997464
  6         1.56240592        0.27286765        6.00626336
  7         2.21696005       -1.36035224        7.77385160
 16         2.96687507       -0.62452393        6.48268892
  6        -2.60914750        1.53136247        8.68767798
  7        -3.41083244        2.35107200        8.87891228
  6        -1.74176702       -0.31371451       10.65928885
  7        -2.15517498       -0.33257465       11.74541237
  6         1.55613690        1.14944331        4.89417332
  7         1.56682743        1.85803076        3.97282352

@aug-cc-pVTZ.gbs/N

