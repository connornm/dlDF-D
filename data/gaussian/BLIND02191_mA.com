%Chk=BLIND02191_mA
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.50356529        0.17067208        0.87086276
  6         0.21489342        1.30305686        0.39405206
 16        -0.03557077       -1.50450240        0.61929123
 16         1.69738993        1.20441462       -0.56878860
  6         0.73362554       -1.38817556       -0.97410990
  6         1.42129769       -0.31300754       -1.44183289
  6        -1.58402256        0.58880500        1.62124944
  7        -0.24006352        2.48764054        0.72213088
 16        -1.58845423        2.32110984        1.68391406
  6         0.60742194       -2.58041298       -1.74465013
  7         0.49879163       -3.58170761       -2.32499618
  6         2.04424839       -0.34001125       -2.72446501
  7         2.59178982       -0.32821656       -3.74964080
  6        -2.55067674       -0.20280432        2.28765984
  7        -3.34729177       -0.83692331        2.84836790

@blind-aug.gbs/N
