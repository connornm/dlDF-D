%Chk=CH400008
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.06003836        0.00000000       -0.31985095
  1        -0.08498260       -0.47614378        0.99601758
  1        -0.36841091        1.04261480        0.05668981
  1        -0.60664485       -0.56647102       -0.73285644
  0        -0.70144293        0.00000000        0.21165007
  0         0.05623423        0.31507132       -0.65907944
  0         0.24378290       -0.68991351       -0.03751248
  0         0.40142580        0.37484218        0.48494185
  6         0.00000000        0.00000000        2.94387027
  1         0.10834846        0.77395019        2.15949303
  1         0.77954514        0.14591697        3.71653034
  1        -1.00119020        0.08379405        3.40925474
  1         0.11329660       -1.00366121        2.49020299
  0        -0.07169577       -0.51213420        3.46290420
  0        -0.51583645       -0.09655540        2.43258978
  0         0.66250224       -0.05544775        2.63591855
  0        -0.07497002        0.66413735        3.24406857

@aug-cc-pVTZ.gbs/N

