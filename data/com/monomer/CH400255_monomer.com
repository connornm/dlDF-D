%Chk=CH400255
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.09712304        0.00000000       -0.14935527
  1        -0.24862491       -0.57996302        0.90984312
  1        -0.35649915        1.04168367        0.11742843
  1        -0.49199898       -0.46172066       -0.87791628
  6-Bq        0.00000000       0.00000000       4.44657346
  1-Bq        0.88320420      -0.60694089       4.72506772
  1-Bq       -0.70101808      -0.61885742       3.85363578
  1-Bq       -0.50737183       0.35599894       5.36408342
  1-Bq        0.32518571       0.86979937       3.84350691

@aug-cc-pVTZ.gbs/N

