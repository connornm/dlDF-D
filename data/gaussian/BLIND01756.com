%Chk=BLIND01756
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.48699999        0.87619093        0.19034015
  6        -0.19132635        0.67351016       -1.18709379
 16        -0.80690927       -0.40003148        1.35542858
 16        -0.06516580       -0.91859734       -1.95120475
  6         0.23671573       -1.68508457        0.72012748
  6         0.52786913       -1.88891273       -0.59186635
  6        -0.55481545        2.22895516        0.45629732
  7        -0.03747039        1.73956018       -1.93390416
 16        -0.27565013        3.10554192       -1.01288236
  6         0.73374039       -2.55463710        1.73399855
  7         1.09461845       -3.25002557        2.59271558
  6         1.34312594       -2.98472514       -1.00234590
  7         1.97809713       -3.88011798       -1.38445764
  6        -0.83416216        2.86990625        1.68768498
  7        -1.07084808        3.41267542        2.68792361
  6        -0.08556208       -0.99877858        7.61127018
  6        -0.21756792       -0.66549675        6.23383621
 16         0.89273112       -0.11899006        8.77635868
 16         0.56722153        0.72548433        5.46972528
  6         0.74503574        1.52985819        8.14105750
  6         0.61741998        1.86156646        6.82906364
  6        -0.78503456       -2.15865283        7.87722733
  7        -0.94072336       -1.46373257        5.48702580
 16        -1.50627887       -2.72974312        6.40804759
  6         0.81857466        2.52873116        9.15492853
  7         0.90773786        3.30709357       10.01364554
  6         0.55341535        3.22588015        6.41858403
  7         0.52697652        4.32324869        6.03647225
  6        -0.91140914       -2.84631722        9.10861500
  7        -1.01831765       -3.42871669       10.10885364

@blind-aug.gbs/N
