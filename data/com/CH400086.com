%Chk=CH400086
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.73892355        0.00000000       -0.82460775
  1        -0.91474683        0.53308070       -0.32411288
  1        -0.25590608       -1.04383868        0.26626836
  1         0.43172936        0.51075798        0.88245227
  0        -0.48895655        0.00000000        0.54565504
  0         0.60530140       -0.35274732        0.21447025
  0         0.16933681        0.69072337       -0.17619368
  0        -0.28568165       -0.33797605       -0.58393161
  6         0.00000000        0.00000000        3.65708684
  1         0.85451146       -0.69240761        3.78502361
  1        -0.94381786       -0.57855809        3.63553842
  1        -0.02491021        0.71552426        4.50171130
  1         0.11421661        0.55544144        2.70607405
  0        -0.56544277        0.45817628        3.57242921
  0         0.62453812        0.38284038        3.67134575
  0         0.01648345       -0.47347290        3.09818645
  0        -0.07557880       -0.36754375        4.28638596

@aug-cc-pVTZ.gbs/N

