%Chk=BLIND01468
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.54702515        0.61984004       -0.59805578
  6        -0.62664053       -0.75157232       -0.97050197
 16        -0.31294302        1.21196253        1.03996533
 16        -0.46871492       -2.09948857        0.16649831
  6         0.69938696       -0.08133043        1.70832409
  6         0.63429074       -1.39382961        1.36082960
  6        -0.75297197        1.41477934       -1.70754244
  7        -0.86494866       -1.02117332       -2.23081326
 16        -1.04778145        0.39850026       -3.08052851
  6         1.59722078        0.37540850        2.71647766
  7         2.29426084        0.79081233        3.54880443
  6         1.46428085       -2.36441993        1.99566247
  7         2.09972028       -3.19297566        2.50618992
  6        -0.76990844        2.82858879       -1.78515978
  7        -0.79646427        3.98786181       -1.86563804
  6         0.85354022       -0.34984776       11.78697153
  6        -0.12520488        0.09410521       12.72013306
 16         0.50723635       -1.04360771       10.20985790
 16        -1.87036914        0.04776592       12.42544217
  6        -0.97075393       -0.16001788        9.78684985
  6        -1.91236672        0.27027165       10.66767919
  6         2.11259229       -0.21902500       12.33740454
  7         0.30185912        0.53172283       13.87942905
 16         1.96051734        0.40803020       13.94628662
  6        -1.12421028        0.03408631        8.38333129
  7        -1.21698392        0.14699283        7.23014518
  6        -3.09225859        0.93274634       10.21720789
  7        -4.08137881        1.45011679        9.89333268
  6         3.36278285       -0.55110761       11.76128550
  7         4.39692274       -0.82656269       11.30766109

@aug-cc-pVTZ.gbs/N
