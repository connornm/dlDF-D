%Chk=CH400100
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.69515882        0.00000000       -0.86182375
  1        -0.08586533       -1.02696191        0.40492263
  1         0.38759961        0.67576435        0.78682593
  1        -0.99689310        0.35119756       -0.32992481
  0        -0.45999679        0.00000000        0.57028142
  0         0.05681835        0.67955576       -0.26794325
  0        -0.25648035       -0.44716318       -0.52065426
  0         0.65965879       -0.23239258        0.21831609
  6         0.00000000        0.00000000        2.84328672
  1         0.25192920       -0.28537366        3.88303645
  1         0.86380580        0.51768043        2.38303411
  1        -0.23960121       -0.90932221        2.25872677
  1        -0.87613379        0.67701544        2.84834955
  0        -0.16670525        0.18883594        2.15526907
  0        -0.57159297       -0.34255674        3.14784262
  0         0.15854763        0.60171184        3.23009861
  0         0.57975058       -0.44799105        2.83993657

@aug-cc-pVTZ.gbs/N

