%Chk=CH400082
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.88495864        0.00000000       -0.66545786
  1         0.29950308        0.33366036        1.01240043
  1        -0.76589616        0.68980688       -0.40441995
  1        -0.41856556       -1.02346723        0.05747739
  6-Bq        0.00000000       0.00000000       4.08421178
  1-Bq       -0.63017091       0.60448002       3.40342217
  1-Bq       -0.32535643       0.16258144       5.13001120
  1-Bq       -0.10325247      -1.07223512       3.82801512
  1-Bq        1.05877981       0.30517367       3.97539862

@aug-cc-pVTZ.gbs/N

