%Chk=BLIND00113
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.08512684        0.35977688        0.95101153
  6        -0.50180283        1.28295341        0.04047031
 16        -0.04832937       -1.38940984        0.84588353
 16        -1.44611996        0.81083216       -1.38078091
  6        -0.11445872       -1.60219967       -0.91322689
  6        -0.67066335       -0.72844022       -1.79346827
  6         0.71451687        1.04803817        1.96859091
  7        -0.35568781        2.55967838        0.29831479
 16         0.49734630        2.74904536        1.71513865
  6         0.44673708       -2.83768401       -1.34880394
  7         0.89341404       -3.86680694       -1.65307146
  6        -0.71368233       -1.01827479       -3.18920722
  7        -0.79687589       -1.22602116       -4.32975376
  6         1.41138847        0.52416462        3.08442912
  7         1.97895065        0.11165357        4.01114575
  6         0.83350306        0.46289750        7.74326962
  6         1.04865774       -0.21952649        6.51289057
 16        -0.68184880        1.20187022        8.23939932
 16        -0.20048299       -0.46463443        5.28237935
  6        -1.84256333        0.10244557        7.47250342
  6        -1.64924781       -0.55591134        6.29911545
  6         2.01435430        0.51836245        8.45588256
  7         2.25904370       -0.66081880        6.27187817
 16         3.26530368       -0.24703282        7.53169519
  6        -3.06785290       -0.00807641        8.19190597
  7        -4.05733223       -0.04722643        8.80066778
  6        -2.66965162       -1.38288364        5.74338637
  7        -3.48280779       -2.04299786        5.23944947
  6         2.25173754        1.11244733        9.71916432
  7         2.46677289        1.60285769       10.75083936

@aug-cc-pVTZ.gbs/N
