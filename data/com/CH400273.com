%Chk=CH400273
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.85756223        0.00000000        0.70040915
  1        -0.71830499        0.78892397        0.29601161
  1         0.36256151        0.19758801       -1.02737242
  1        -0.50181875       -0.98651198        0.03095165
  0        -0.56746150        0.00000000       -0.46347101
  0         0.47531295       -0.52204256       -0.19587522
  0        -0.23991227       -0.13074689        0.67982740
  0         0.33206083        0.65278945       -0.02048116
  6         0.00000000        0.00000000        3.18196860
  1        -0.60673434        0.86934969        2.86245253
  1         0.95766032        0.35534504        3.60928882
  1        -0.55482665       -0.57580705        3.94786650
  1         0.20390066       -0.64888767        2.30826656
  0         0.40148501       -0.57526144        3.39339707
  0        -0.63369788       -0.23513702        2.89920455
  0         0.36713693        0.38101997        2.67516273
  0        -0.13492406        0.42937849        3.76011006

@aug-cc-pVTZ.gbs/N

