%Chk=CH400243
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.01865113        0.00000000        0.43397676
  1        -0.59296148        0.81959672        0.45015981
  1         0.06535354        0.15013010       -1.09506887
  1        -0.49104319       -0.96972682        0.21093230
  0        -0.67405640        0.00000000       -0.28716879
  0         0.39237131       -0.54233917       -0.29787735
  0        -0.04324539       -0.09934329        0.72462313
  0         0.32493048        0.64168246       -0.13957699
  6         0.00000000        0.00000000        3.76616936
  1         0.38855459       -0.82648595        4.39221429
  1        -0.60895571        0.68158167        4.39115361
  1         0.84795319        0.56056929        3.32716661
  1        -0.62755207       -0.41566502        2.95414291
  0        -0.25711227        0.54689788        3.35190624
  0         0.40295493       -0.45101259        3.35260810
  0        -0.56110306       -0.37093693        4.05666391
  0         0.41526041        0.27505164        4.30349916

@aug-cc-pVTZ.gbs/N

