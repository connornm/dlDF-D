%Chk=BLIND00384
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.03261536       -0.90234081       -0.47521175
  6        -1.26334905       -0.55077372       -0.00329705
 16         1.26348359        0.23808600       -0.99761649
 16        -1.83995087        1.11462751        0.16569236
  6         0.91042541        1.60618089        0.07373238
  6        -0.32157334        1.95198277        0.53225631
  6         0.15034766       -2.27663058       -0.52784213
  7        -2.08415844       -1.52956311        0.29013658
 16        -1.34593771       -2.98904270       -0.01948631
  6         2.06504847        2.37701007        0.39584667
  7         3.03054587        2.98761597        0.61071462
  6        -0.50366140        3.10322561        1.35409923
  7        -0.70279131        4.04803799        2.00109457
  6         1.26096885       -3.04841707       -0.94740695
  7         2.15954725       -3.69756399       -1.29706543
  6        -0.54189257        0.17760605        5.53464242
  6        -1.31613594       -0.33694620        4.45686237
 16         1.15568420       -0.17584901        5.82068305
 16        -0.67727393       -1.40405378        3.19681273
  6         1.73203470       -0.35800987        4.15372186
  6         1.00311970       -0.84665258        3.11566826
  6        -1.34433142        0.93853795        6.36070401
  7        -2.58781747       -0.02138721        4.42356627
 16        -2.96415199        0.92014563        5.74371972
  6         3.09522015        0.02335443        3.98826381
  7         4.21535766        0.32343161        3.90802253
  6         1.57873728       -0.99674733        1.81949007
  7         2.01463503       -1.16341885        0.75499964
  6        -0.98622606        1.62185450        7.54814138
  7        -0.71212960        2.18104417        8.52960659

@blind-aug.gbs/N

