%Chk=CH400063
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10655631        0.00000000        0.03897547
  1        -0.40367905        0.33298043        0.97578339
  1        -0.34128789        0.69034524       -0.79557023
  1        -0.36158938       -1.02332567       -0.21918863
  0        -0.73222454        0.00000000       -0.02579064
  0         0.26712035       -0.22033803       -0.64569018
  0         0.22583520       -0.45681157        0.52644049
  0         0.23926899        0.67714961        0.14504033
  6         0.00000000        0.00000000        4.31077708
  1        -0.28221916       -0.63592217        5.17213722
  1         1.03431590       -0.24081852        3.99743537
  1        -0.69414125       -0.18929209        3.46915817
  1        -0.05795549        1.06603278        4.60437755
  0         0.18674856        0.42079903        3.74080244
  0        -0.68442200        0.15935315        4.51811988
  0         0.45932345        0.12525735        4.86768865
  0         0.03834999       -0.70540953        4.11649734

@aug-cc-pVTZ.gbs/N
