%Chk=CH400011
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.19373376        0.00000000       -1.09016200
  1        -1.08322978       -0.13908060        0.18236171
  1         0.56333730       -0.82645993        0.47497476
  1         0.32615872        0.96554053        0.43282552
  6-Bq        0.00000000       0.00000000       3.98629063
  1-Bq        0.27587862      -0.85042765       4.63947495
  1-Bq       -0.06914567       0.92385082       4.59267706
  1-Bq        0.77195877       0.12945842       3.20315199
  1-Bq       -0.97869171      -0.20288159       3.50985853

@aug-cc-pVTZ.gbs/N

