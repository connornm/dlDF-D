%Chk=BLIND02447
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6-Bq        0.44182781      -0.88942957       0.23411991
  6-Bq        0.49095925       0.06818218       1.28597185
 16-Bq        0.44005545      -0.50921884      -1.48167532
 16-Bq        0.52641829       1.82003240       1.03229525
  6-Bq       -0.43930652       1.03048722      -1.46943316
  6-Bq       -0.40194774       1.95269682      -0.47154648
  6-Bq        0.46801924      -2.16342265       0.76474112
  7-Bq        0.54648196      -0.38717714       2.51374427
 16-Bq        0.57653461      -2.05135961       2.49115152
  6-Bq       -1.18754859       1.25691193      -2.66105193
  7-Bq       -1.76227440       1.40728358      -3.66018738
  6-Bq       -1.11180580       3.18469135      -0.58214912
  7-Bq       -1.65011025       4.21308346      -0.64184941
  6-Bq        0.44750279      -3.40205854       0.07879023
  7-Bq        0.44073662      -4.42771282      -0.46809500
  6         0.43286539       -0.89425759        7.90798192
  6        -0.66004892       -0.55091766        8.75266094
 16         0.59528838       -0.42567985        6.22198954
 16        -2.03386331        0.44297794        8.24335368
  6        -0.20992817        1.15390415        6.24773459
  6        -1.25319631        1.49487716        7.04959058
  6         1.31030773       -1.71043169        8.59300600
  7        -0.65496685       -1.03344623        9.97126686
 16         0.68975249       -1.99054826       10.18714291
  6         0.31644020        2.06590379        5.28737911
  7         0.76291051        2.76554784        4.47352685
  6        -1.86125296        2.78162814        6.95749843
  7        -2.40557932        3.80678339        6.89673930
  6         2.51677595       -2.28912425        8.12962949
  7         3.50475828       -2.77974555        7.76319042

@aug-cc-pVTZ.gbs/N

