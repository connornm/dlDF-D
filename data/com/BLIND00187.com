%Chk=BLIND00187
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.45685950       -0.83968732       -0.35681019
  6        -0.38042152       -0.78962846        1.06356934
 16        -0.58067590        0.55888182       -1.41378252
 16        -0.36540719        0.71012985        2.00436045
  6         0.35854401        1.74528855       -0.48943027
  6         0.44083109        1.80345905        0.86610987
  6        -0.49108680       -2.15511763       -0.77367325
  7        -0.35329494       -1.93145942        1.70648943
 16        -0.45269100       -3.18730697        0.61842498
  6         1.01488933        2.70445926       -1.31422845
  7         1.51118039        3.47740353       -2.02658044
  6         1.18866212        2.83067879        1.51380603
  7         1.76157401        3.66585844        2.08404549
  6        -0.57736178       -2.65661098       -2.09512302
  7        -0.65724563       -3.08574229       -3.17241630
  6         0.85134784       -0.00847247        5.25299020
  6         1.07720248        0.27854001        6.62865010
 16        -0.44050615       -1.02387058        4.62961504
 16         0.06098774       -0.32297069        7.94777133
  6        -1.71302656       -0.69257006        5.81917541
  6        -1.51099740       -0.41625030        7.13465085
  6         1.84277716        0.57310113        4.48869754
  7         2.12466543        1.00710053        6.92845163
 16         2.95799724        1.38551488        5.53802811
  6        -3.02651471       -0.78654819        5.27419225
  7        -4.07747103       -0.89792619        4.79016873
  6        -2.61010235       -0.21083709        8.01992802
  7        -3.47503928       -0.06553849        8.78261967
  6         2.00862127        0.52970956        3.08307264
  7         2.16545470        0.49689170        1.93180283

@aug-cc-pVTZ.gbs/N
