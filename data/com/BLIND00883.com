%Chk=BLIND00883
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.44835997       -0.89917142       -0.17768794
  6        -1.26391836        0.01305078        0.54932089
 16         0.76652745       -0.44677240       -1.36423769
 16        -1.16223817        1.77400473        0.39670809
  6         1.32452612        1.09231535       -0.68308540
  6         0.55684998        1.97146525        0.01363036
  6        -0.80588132       -2.19449167        0.13794974
  7        -2.16546468       -0.49387913        1.35436519
 16        -2.11835965       -2.15566112        1.26963847
  6         2.69173770        1.36999318       -0.97414757
  7         3.80252197        1.56329787       -1.25696309
  6         1.09368880        3.20806329        0.47900259
  7         1.48058219        4.23884035        0.85147848
  6        -0.26290083       -3.40299458       -0.36191159
  7         0.16448422       -4.40460466       -0.76837831
  6         0.41407292        0.12089873        9.10052644
  6         0.29800815        1.33253993        8.36277226
 16        -0.67170489       -1.25437905        8.96556641
 16        -0.95617528        1.64079656        7.15174373
  6        -1.11415620       -1.14841295        7.25176067
  6        -1.22742196        0.00134376        6.53560816
  6         1.46585326        0.21756005        9.98918719
  7         1.15846751        2.28487799        8.62833193
 16         2.18179969        1.78957523        9.84421979
  6        -1.39826188       -2.41703389        6.66783042
  7        -1.64295710       -3.47098916        6.24298274
  6        -1.63711537       -0.02066844        5.16973528
  7        -1.99540511        0.01117393        4.06450198
  6         1.92304321       -0.74798301       10.91863454
  7         2.30593151       -1.52507381       11.69363284

@aug-cc-pVTZ.gbs/N
