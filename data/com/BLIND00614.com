%Chk=BLIND00614
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.61540093        0.81136729        0.06382572
  6         0.28120679        0.85881855       -1.04056082
 16        -1.27410047       -0.65668206        0.77055244
 16         0.91026796       -0.57306498       -1.87036799
  6         0.07934306       -1.77309530        0.51377590
  6         0.94302987       -1.73814236       -0.53525577
  6        -0.94485046        2.09507611        0.44964678
  7         0.63018241        2.04202130       -1.48340085
 16        -0.14282730        3.22013923       -0.59716666
  6         0.16872137       -2.78355392        1.51481263
  7         0.18639937       -3.60122784        2.34082103
  6         1.97252293       -2.71547157       -0.67278241
  7         2.80672074       -3.50725732       -0.84045909
  6        -1.81112781        2.50471801        1.49220399
  7        -2.52595949        2.85888000        2.33759198
  6         0.95936495        0.31850357        5.87047338
  6         0.49889577       -0.56899317        6.88346120
 16         0.43710791        0.29771620        4.19257824
 16        -0.71331466       -1.83222946        6.62008258
  6        -1.24185921       -0.23926799        4.38451432
  6        -1.69469695       -1.08506030        5.34743560
  6         1.93059284        1.15094783        6.38917962
  7         1.03673867       -0.45672315        8.07341098
 16         2.19367721        0.73992775        8.05251232
  6        -2.11150811        0.27276391        3.37822223
  7        -2.77611701        0.69569729        2.52346180
  6        -3.06086820       -1.49265049        5.38469545
  7        -4.15762493       -1.87396931        5.43614222
  6         2.65940926        2.16979356        5.72891571
  7         3.27340610        3.00347980        5.20064395

@aug-cc-pVTZ.gbs/N

