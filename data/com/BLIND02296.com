%Chk=BLIND02296
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.27078382       -0.24393190       -0.95303844
  6         0.37317959       -1.30244715       -0.25255058
 16         0.03197623        1.46865390       -0.70024826
 16         1.57495257       -1.06041338        1.02484901
  6         0.40313626        1.47758203        1.03361523
  6         1.01632235        0.47361104        1.71456245
  6        -1.11624338       -0.76766065       -1.91022679
  7         0.07474867       -2.52871100       -0.60621142
 16        -1.01163913       -2.49740622       -1.86704386
  6         0.03195538        2.69283912        1.67895592
  7        -0.26637259        3.71021836        2.15554056
  6         1.31221167        0.60129416        3.10383868
  7         1.59863681        0.67340976        4.22797663
  6        -1.93703414       -0.07415913       -2.83246679
  7        -2.60970048        0.47864024       -3.60255040
  6         0.27078381        0.24393179        9.85724847
  6        -0.37317966        1.30244710        9.15676074
 16        -0.03197614       -1.46865399        9.60445807
 16        -1.57495263        1.06041342        7.87936113
  6        -0.40313618       -1.47758192        7.87059458
  6        -1.01632233       -0.47361088        7.18964749
  6         1.11624334        0.76766047       10.81443688
  7        -0.07474882        2.52871092        9.51042174
 16         1.01163899        2.49740604       10.77125417
  6        -0.03195522       -2.69283891        7.22525374
  7         0.26637280       -3.71021807        6.74866896
  6        -1.31221165       -0.60129384        5.80037124
  7        -1.59863678       -0.67340931        4.67623329
  6         1.93703414        0.07415888       11.73667680
  7         2.60970052       -0.47864055       12.50676033

@aug-cc-pVTZ.gbs/N

