%Chk=BLIND01816
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.55251598        0.06515870       -0.85532996
  6        -0.49248174       -1.24888299       -0.31171523
 16        -0.48435063        1.55156395        0.07980082
 16        -0.30496366       -1.60170183        1.41308741
  6         0.58116838        1.04959212        1.40525003
  6         0.64902426       -0.20135996        1.93265739
  6        -0.74421536       -0.00564324       -2.22044594
  7        -0.61809194       -2.25249554       -1.14528823
 16        -0.85658654       -1.67138244       -2.68682798
  6         1.36484110        2.11872824        1.92859952
  7         1.96473786        3.02676796        2.33698751
  6         1.50909637       -0.48952196        3.03319388
  7         2.17286745       -0.76251237        3.94742050
  6        -0.87189959        1.05788330       -3.14660817
  7        -0.98799956        1.91833713       -3.91941858
  6        -0.55313293        0.06613958        7.14514421
  6        -0.49466915       -1.24807498        7.68851217
 16        -0.48188960        1.55229711        8.08043924
 16        -0.30624709       -1.60144229        9.41310401
  6         0.58396659        1.04848292        9.40491798
  6         0.65033570       -0.20265172        9.93208144
  6        -0.74608764       -0.00416108        5.78017918
  7        -0.62252723       -2.25136709        6.85489519
 16        -0.86141969       -1.66965475        5.31364309
  6         1.36972673        2.11633026        9.92776855
  7         1.97136616        3.02338215       10.33578800
  6         1.51088623       -0.49230560       11.03185203
  7         2.17500297       -0.76645710       11.94547998
  6        -0.87290898        1.05970055        4.85428329
  7        -0.98833054        1.92044877        4.08169913

@aug-cc-pVTZ.gbs/N

