%Chk=BLIND01499
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.12497879        0.87832939       -0.50401061
  6         0.55165635       -0.12213286       -1.25704901
 16        -1.33908874        0.56957698        0.72851753
 16         0.29906050       -1.86200964       -1.04821303
  6        -0.74810964       -0.96805912        1.38467377
  6        -0.09999236       -1.93089577        0.67718976
  6         0.26476227        2.12919172       -0.93846264
  7         1.38698856        0.28185799       -2.18272255
 16         1.39653510        1.94550197       -2.23852990
  6        -1.05874345       -1.14382462        2.76442113
  7        -1.35555275       -1.25182906        3.88309548
  6         0.29278396       -3.15623402        1.29218152
  7         0.60496185       -4.18049955        1.74423376
  6        -0.16742762        3.39522701       -0.47420518
  7        -0.51774973        4.44267940       -0.11197070
  6         0.51240588       -0.88221501        6.01961183
  6         1.07288197        0.11996486        5.17858346
 16        -0.36241616       -0.57732580        7.51301252
 16         0.94976403        1.85867888        5.48891215
  6        -1.13704807        0.96914212        7.12249430
  6        -0.61350659        1.93370123        6.32051295
  6         0.82183849       -2.13199146        5.52200406
  7         1.74833408       -0.28179306        4.12967210
 16         1.78696216       -1.94557394        4.09420458
  6        -2.38840220        1.14988302        7.78016840
  7        -3.39163628        1.26152499        8.35658012
  6        -1.29960021        3.16590135        6.10861591
  7        -1.80902300        4.19541185        5.93106968
  6         0.47517743       -3.39918322        6.05041640
  7         0.20762203       -4.44765082        6.47494640

@blind-aug.gbs/N
