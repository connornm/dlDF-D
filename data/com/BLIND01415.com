%Chk=BLIND01415
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.28940976        0.17507464       -0.96265247
  6         0.02051422       -1.19500823       -0.68625641
 16        -0.49776515        1.53865095       -0.18192720
 16        -1.14216915       -1.75307798        0.52672086
  6        -0.77908345        0.86589385        1.43435866
  6        -1.03569803       -0.43975438        1.71198382
  6         1.20996034        0.26993378       -1.98680655
  7         0.65452353       -2.09041014       -1.40329731
 16         1.62657486       -1.32663603       -2.51805389
  6        -0.75850412        1.85266935        2.46243466
  7        -0.75588236        2.69609475        3.26232198
  6        -1.29552129       -0.87025429        3.04656376
  7        -1.54130925       -1.26043419        4.11346304
  6         1.74644207        1.43780884       -2.58130099
  7         2.19011209        2.38568039       -3.08706019
  6        -0.50102387       -0.29495812        9.60175086
  6        -1.29542601        0.46593712        8.69853975
 16         0.82639211       -1.35024357        9.14038760
 16        -1.05053678        0.49688498        6.94535412
  6         1.44866454       -0.48441658        7.72361566
  6         0.70156018        0.24615577        6.85428076
  6        -0.97692137       -0.13385625       10.88739116
  7        -2.28798605        1.15699530        9.20363697
 16        -2.36041558        0.90980746       10.84819436
  6         2.85025343       -0.65826072        7.53278470
  7         3.98985118       -0.84874297        7.40540146
  6         1.29542754        0.86614427        5.71547535
  7         1.73129223        1.36373691        4.75978206
  6        -0.49223784       -0.71624442       12.08366636
  7        -0.11360677       -1.19261218       13.07402396

@aug-cc-pVTZ.gbs/N
