%Chk=CH400287
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.09729676        0.00000000       -0.14807352
  1        -0.23255074        0.31222003        1.03654463
  1        -0.46853282        0.70656777       -0.71219727
  1        -0.39621321       -1.01878780       -0.17627384
  0        -0.72609736        0.00000000        0.09798242
  0         0.15388223       -0.20660058       -0.68589678
  0         0.31003504       -0.46754626        0.47127138
  0         0.26218009        0.67414683        0.11664298
  6         0.00000000        0.00000000        5.01328434
  1        -0.99734118       -0.17459829        5.46140481
  1         0.26475508        1.07095669        5.10784901
  1         0.75457263       -0.61474144        5.54120391
  1        -0.02198652       -0.28161696        3.94267964
  0         0.65995529        0.11553425        4.71675645
  0        -0.17519232       -0.70866775        4.95070951
  0        -0.49931177        0.40678343        4.66395222
  0         0.01454880        0.18635007        5.72171918

@aug-cc-pVTZ.gbs/N
