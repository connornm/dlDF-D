%Chk=BLIND02141
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.68782192        0.54628405       -0.51921446
  6         0.79178942       -0.85695212       -0.73356387
 16         0.12790707        1.30499944        0.96381047
 16         0.33777277       -2.08002513        0.46332778
  6        -1.05574784        0.10430619        1.51282671
  6        -0.96901337       -1.23754522        1.31393982
  6         1.16076810        1.22121478       -1.62650559
  7         1.28674720       -1.25676237       -1.87943847
 16         1.70562927        0.06626685       -2.79870359
  6        -2.13192638        0.67667261        2.25128868
  7        -2.97577293        1.18653902        2.86704210
  6        -1.95473631       -2.12391281        1.83991006
  7        -2.71810798       -2.88506094        2.27447116
  6         1.25022548        2.61941103       -1.83196372
  7         1.33958380        3.76395259       -2.01399569
  6        -0.68782195       -0.54628392        6.28802455
  6        -0.79178937        0.85695229        6.50237373
 16        -0.12790716       -1.30499959        4.80499974
 16        -0.33777263        2.08002508        5.30548189
  6         1.05574784       -0.10430650        4.25598331
  6         0.96901345        1.23754494        4.45486998
  6        -1.16076818       -1.22121444        7.39531578
  7        -1.28674711        1.25676276        7.64824827
 16        -1.70562928       -0.06626629        8.56751360
  6         2.13192633       -0.67667311        3.51752142
  7         2.97577285       -1.18653968        2.90176809
  6         1.95473646        2.12391238        3.92889960
  7         2.71810818        2.88506040        3.49433838
  6        -1.25022565       -2.61941066        7.60077414
  7        -1.33958405       -3.76395218        7.78280629

@blind-aug.gbs/N
