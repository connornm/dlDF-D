%Chk=CH400276
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.99897834       0.00000000       0.47752303
  1-Bq       -0.70216908      -0.59749535       0.61314260
  1-Bq       -0.37156482       1.04007989      -0.07848161
  1-Bq        0.07475557      -0.44258454      -1.01218402
  6         0.00000000        0.00000000        3.42346600
  1         0.70024224       -0.48719592        2.71757161
  1        -1.03436841       -0.32725996        3.20217688
  1         0.26342333       -0.28505805        4.46045009
  1         0.07070283        1.09951394        3.31366542

@aug-cc-pVTZ.gbs/N

