%Chk=CH400195
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.05525326        0.00000000        0.33530062
  1        -0.17371636        0.86262595       -0.67207518
  1        -0.21454128       -0.94046062       -0.54359159
  1        -0.66699562        0.07783467        0.88036615
  0        -0.69827656        0.00000000       -0.22187334
  0         0.11495066       -0.57081224        0.44472201
  0         0.14196511        0.62231659        0.35970253
  0         0.44136079       -0.05150434       -0.58255120
  6         0.00000000        0.00000000        3.15573942
  1         0.81121508       -0.43340582        2.53923682
  1         0.31735760        0.02853123        4.21614315
  1        -0.21858123        1.02827977        2.80810611
  1        -0.90999145       -0.62340518        3.05947159
  0        -0.53679291        0.28679099        3.56368823
  0        -0.21000018       -0.01887953        2.45405472
  0         0.14463840       -0.68042781        3.38577348
  0         0.60215469        0.41251635        3.21944125

@aug-cc-pVTZ.gbs/N

