%Chk=BLIND01812
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.70528680        0.36216790       -0.64227241
  6        -0.68635097       -1.05835405       -0.55517783
 16        -0.27585053        1.46194917        0.65948883
 16        -0.19016216       -1.95767916        0.88690377
  6         0.97439279        0.51255402        1.48390532
  6         1.00485571       -0.84343004        1.57347822
  6        -1.18033352        0.74493828       -1.88046491
  7        -1.09520733       -1.73209150       -1.60244779
 16        -1.57700937       -0.67399554       -2.79375657
  6         1.96656727        1.31990386        2.11228872
  7         2.73797625        2.02013347        2.62784638
  6         2.03379583       -1.50905833        2.30297334
  7         2.83540759       -2.09156949        2.91044800
  6        -1.37246346        2.05527576       -2.38173353
  7        -1.54517243        3.12356571       -2.80602599
  6        -0.82734192       -0.27136651        6.66211654
  6         0.12461045        0.46570186        5.90294359
 16        -0.44359615       -1.34056002        8.00304199
 16         1.86787527        0.44859352        6.21168847
  6         0.97249100       -0.51443684        8.67857124
  6         1.88805131        0.19318001        7.96535895
  6        -2.08989335       -0.07470380        6.14001418
  7        -0.32494497        1.17155461        4.89421576
 16        -1.97186512        0.96932010        4.76127926
  6         1.10721119       -0.69524956       10.08578369
  7         1.18766669       -0.89062343       11.22882743
  6         3.02082171        0.78115980        8.60179663
  7         3.97311706        1.25198094        9.07337607
  6        -3.31844645       -0.62494762        6.57939723
  7        -4.33454284       -1.07467660        6.92055392

@aug-cc-pVTZ.gbs/N

