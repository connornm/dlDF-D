%Chk=CH400224
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.54684026        0.00000000        0.96278330
  1         0.54810279        0.61984871       -0.73576957
  1        -0.08070238       -1.03736140       -0.37862172
  1        -1.01424067        0.41751269        0.15160799
  0        -0.36185222        0.00000000       -0.63708783
  0        -0.36268766       -0.41016299        0.48686952
  0         0.05340195        0.68643725        0.25053955
  0         0.67113793       -0.27627427       -0.10032123
  6         0.00000000        0.00000000        4.12585674
  1        -0.26900754       -0.80851210        3.41880010
  1        -0.84307486        0.17781010        4.82126619
  1         0.89746764       -0.29720398        4.70223126
  1         0.21461477        0.92790597        3.56112941
  0         0.17800624        0.53500431        4.59372649
  0         0.55787500       -0.11765955        3.66569411
  0        -0.59386750        0.19666423        3.74446126
  0        -0.14201374       -0.61400899        4.49954510

@aug-cc-pVTZ.gbs/N

