%Chk=CH400281
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.79573186        0.00000000        0.76993296
  1        -0.05738262        1.00020428       -0.47147066
  1         0.23314846       -0.75897956       -0.77173686
  1        -0.97149770       -0.24122472        0.47327456
  0        -0.52654744        0.00000000       -0.50947593
  0         0.03797092       -0.66184984        0.31197905
  0        -0.15427776        0.50222791        0.51066960
  0         0.64285428        0.15962194       -0.31317272
  6         0.00000000        0.00000000        4.12510106
  1        -0.11664510       -0.65689558        5.00876846
  1         1.05152251       -0.02947306        3.77953240
  1        -0.66535836       -0.35205540        3.31310278
  1        -0.26951905        1.03842404        4.39900059
  0         0.07718577        0.43467744        3.54036537
  0        -0.69580787        0.01950276        4.35376891
  0         0.44027739        0.23296022        4.66241223
  0         0.17834471       -0.68714042        3.94385771

@aug-cc-pVTZ.gbs/N
