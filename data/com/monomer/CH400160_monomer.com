%Chk=CH400160
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.43134461        0.00000000        1.01976849
  1         0.69527236       -0.50970186       -0.69482826
  1        -0.96985074       -0.53412117        0.00949030
  1        -0.15676623        1.04382303       -0.33443053
  6-Bq        0.00000000       0.00000000       2.99673724
  1-Bq       -0.84465987       0.71127919       3.07808565
  1-Bq       -0.02924050      -0.49610640       2.00728821
  1-Bq        0.95583924       0.54795103       3.10676895
  1-Bq       -0.08193887      -0.76312382       3.79480617

@aug-cc-pVTZ.gbs/N

