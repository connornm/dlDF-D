%Chk=CH400022
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.71138401        0.00000000        0.84848026
  1        -1.01083513       -0.26521776        0.36586587
  1        -0.02628236        1.00700503       -0.45960425
  1         0.32573349       -0.74178727       -0.75474187
  0        -0.47073323        0.00000000       -0.56145183
  0         0.66888443        0.17549848       -0.24209881
  0         0.01739142       -0.66635000        0.30412687
  0        -0.21554263        0.49085152        0.49942377
  6         0.00000000        0.00000000        2.81692500
  1        -0.68802493        0.09664566        1.95479550
  1        -0.50677022       -0.54948332        3.63377137
  1         0.90617408       -0.55487647        2.50556252
  1         0.28862106        1.00771412        3.17357060
  0         0.45527619       -0.06395185        3.38740873
  0         0.33533728        0.36360117        2.27640577
  0        -0.59962868        0.36716990        3.02295811
  0        -0.19098479       -0.66681922        2.58092737

@aug-cc-pVTZ.gbs/N

