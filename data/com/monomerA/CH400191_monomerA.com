%Chk=CH400191
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.99099391        0.00000000        0.49387957
  1         0.11306894        0.31872747       -1.05433115
  1        -0.42891147       -1.02025488        0.03317948
  1        -0.67515138        0.70152741        0.52727209
  6-Bq        0.00000000       0.00000000       2.97975541
  1-Bq       -0.53130046      -0.19012874       2.02709763
  1-Bq       -0.60812963       0.67438795       3.61328931
  1-Bq        0.16142858      -0.95852559       3.50999637
  1-Bq        0.97800150       0.47426638       2.76863833

@aug-cc-pVTZ.gbs/N

