%Chk=CH400222
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.93435671        0.00000000        0.59410731
  1        -0.01871626        0.89000781       -0.65842368
  1        -0.04425254       -0.91748516       -0.61826260
  1        -0.87138791        0.02747735        0.68257898
  0        -0.61827754        0.00000000       -0.39312952
  0         0.01238482       -0.58893123        0.43568861
  0         0.02928255        0.60711339        0.40911343
  0         0.57661016       -0.01818217       -0.45167252
  6         0.00000000        0.00000000        3.80157128
  1         0.06106471       -0.79523034        4.56959846
  1        -0.97304136       -0.06700447        3.27746550
  1         0.08975148        0.99060047        4.28802922
  1         0.82222517       -0.12836566        3.07119193
  0        -0.04040741        0.52621558        3.29335642
  0         0.64387574        0.04433784        4.14837976
  0        -0.05938987       -0.65549486        3.47967492
  0        -0.54407845        0.08494144        4.28487400

@aug-cc-pVTZ.gbs/N

