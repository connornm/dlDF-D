%Chk=BLIND02133
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.65236863        0.78432771       -0.01881820
  6        -0.42541267        0.91441076       -0.93930315
 16         1.36146881       -0.73731520        0.50108219
 16        -1.25677092       -0.45455624       -1.69386337
  6        -0.06446940       -1.79017519        0.45303909
  6        -1.10049125       -1.67626962       -0.41950220
  6         1.10277123        2.03636780        0.34838246
  7        -0.79534304        2.12838675       -1.26682999
 16         0.17631782        3.23615445       -0.49246392
  6        -0.01708467       -2.84046451        1.41519120
  7         0.07783255       -3.69311017        2.19952475
  6        -2.18043167       -2.60762977       -0.40404631
  7        -3.06562645       -3.35964720       -0.44677941
  6         2.15998411        2.36649892        1.23069918
  7         3.03040650        2.65569105        1.94472208
  6        -0.65236849       -0.78432782        5.08114820
  6         0.42541283       -0.91441069        6.00163314
 16        -1.36146894        0.73731497        4.56124781
 16         1.25677085        0.45455646        6.75619337
  6         0.06446909        1.79017520        4.60929092
  6         1.10049096        1.67626981        5.48183220
  6        -1.10277087       -2.03636799        4.71394754
  7         0.79534341       -2.12838662        6.32915999
 16        -0.17631725       -3.23615448        5.55479392
  6         0.01708418        2.84046451        3.64713880
  7        -0.07783320        3.69311016        2.86280525
  6         2.18043122        2.60763015        5.46637631
  7         3.06562587        3.35964774        5.50910941
  6        -2.15998370       -2.36649930        3.83163082
  7        -3.03040604       -2.65569158        3.11760792

@aug-cc-pVTZ.gbs/N

