%Chk=BLIND00054
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.83149773        0.54038959       -0.24020725
  6        -0.99792275       -0.86468005       -0.39483158
 16         0.22011753        1.31158247        0.93787743
 16        -0.13718676       -2.07761932        0.56563717
  6         1.52856331        0.11853401        1.03129619
  6         1.38344325       -1.22498766        0.88455238
  6        -1.67467989        1.20579792       -1.10714371
  7        -1.86964566       -1.27432718       -1.28373451
 16        -2.59830666        0.04069278       -1.99836206
  6         2.79516218        0.69954070        1.33004055
  7         3.80118942        1.21644735        1.59820042
  6         2.49719619       -2.10449197        1.02586229
  7         3.36992019       -2.86012446        1.16113422
  6        -1.83967694        2.60218247       -1.27463882
  7        -1.99477115        3.74509099       -1.41886715
  6         0.65059085        0.78523860        6.42827174
  6         0.60810662        0.13375456        5.16352263
 16        -0.69729273        0.90848715        7.54922503
 16        -0.81846082       -0.69548865        4.52176771
  6        -1.52699036       -0.62468640        7.22502943
  6        -1.57433210       -1.25757409        6.02300781
  6         1.88923759        1.36684125        6.60936788
  7         1.68761852        0.18569567        4.42202671
 16         2.85818993        1.07604036        5.20174834
  6        -2.20911718       -1.14337183        8.36378229
  7        -2.76873641       -1.51388386        9.31285562
  6        -2.31179856       -2.46699397        5.85751248
  7        -2.92341818       -3.43736325        5.66981316
  6         2.35402078        2.11906645        7.71539785
  7         2.74901793        2.74652717        8.61058844

@aug-cc-pVTZ.gbs/N
