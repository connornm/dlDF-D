%Chk=BLIND00829
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.43711115        0.79258838        0.47100574
  6        -0.43934447        0.89831717       -0.94837549
 16        -0.50109248       -0.71376498        1.37395788
 16        -0.47516872       -0.48971041       -2.04686737
  6         0.38652260       -1.78576222        0.27528439
  6         0.39372398       -1.69513021       -1.08095233
  6        -0.44930889        2.05442080        1.03041022
  7        -0.44878084        2.10363379       -1.46314032
 16        -0.48885769        3.23256514       -0.24051006
  6         1.08829625       -2.82465997        0.95299735
  7         1.62388453       -3.66727924        1.54825443
  6         1.10542015       -2.64039244       -1.87713685
  7         1.64659216       -3.40441780       -2.56584887
  6        -0.46274566        2.40814668        2.40149224
  7        -0.48325138        2.71662541        3.52199048
  6        -0.05920097        0.93961640        7.13427672
  6        -0.72057269       -0.05002702        7.91467407
 16         1.35843206        0.65937368        6.13412811
 16        -0.23448858       -1.75124473        7.97924075
  6         1.06210035       -1.01033544        5.61607305
  6         0.43127242       -1.96433709        6.35063484
  6        -0.65668997        2.16712100        7.33809667
  7        -1.73060601        0.34024113        8.65317491
 16        -1.94564525        1.98203007        8.48228293
  6         1.60046870       -1.30061956        4.32883402
  7         2.07697555       -1.49587402        3.28671601
  6         0.28601517       -3.29488834        5.85821668
  7         0.17086338       -4.39769991        5.50968715
  6        -0.31110238        3.41800371        6.77150149
  7        -0.03566562        4.45446119        6.32318642

@aug-cc-pVTZ.gbs/N

