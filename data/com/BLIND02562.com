%Chk=BLIND02562
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.41869166       -0.92307570        0.11720686
  6         0.18723084       -0.23985949        1.34418206
 16         0.78324003       -0.15034464       -1.41853655
 16         0.19793294        1.52220462        1.51623985
  6        -0.14643115        1.35042300       -1.25261469
  6        -0.37565571        2.01200479       -0.08753543
  6         0.38409540       -2.28624955        0.33136426
  7        -0.01237950       -0.97139480        2.41323390
 16         0.09845889       -2.58386186        2.01501068
  6        -0.61879829        1.85653761       -2.49833013
  7        -0.96178711        2.24210221       -3.53991144
  6        -1.09867539        3.24118333       -0.07008352
  7        -1.65669603        4.25912130       -0.01227861
  6         0.57449992       -3.32881990       -0.60776386
  7         0.73739396       -4.19716632       -1.36310214
  6        -0.41869156        0.92307575        7.50586810
  6        -0.18723081        0.23985949        6.27889291
 16        -0.78324001        0.15034474        9.04161153
 16        -0.19793309       -1.52220462        6.10683515
  6         0.14643101       -1.35042300        8.87568969
  6         0.37565550       -2.01200482        7.71061043
  6        -0.38409516        2.28624958        7.29171068
  7         0.01237961        0.97139477        5.20984107
 16        -0.09845861        2.58386185        5.60806426
  6         0.61879809       -1.85653764       10.12140513
  7         0.96178687       -2.24210226       11.16298645
  6         1.09867505       -3.24118345        7.69315855
  7         1.65669558       -4.25912147        7.63535365
  6        -0.57449958        3.32881997        8.23083879
  7        -0.73739353        4.19716642        8.98617707

@aug-cc-pVTZ.gbs/N
