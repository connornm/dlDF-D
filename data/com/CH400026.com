%Chk=CH400026
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10471228        0.00000000       -0.07481134
  1        -0.36076426        1.03804216        0.13529069
  1        -0.43271339       -0.42323319       -0.92715627
  1        -0.31123463       -0.61480897        0.86667692
  0        -0.73100432        0.00000000        0.04950376
  0         0.23872300       -0.68688773       -0.08952384
  0         0.28633279        0.28005961        0.61351290
  0         0.20594852        0.40682812       -0.57349283
  6         0.00000000        0.00000000        2.91515960
  1        -0.08377373        0.66982011        3.79283206
  1        -0.70810799        0.32965201        2.13036754
  1         1.03338680        0.03777964        2.51934351
  1        -0.24150508       -1.03725176        3.21809530
  0         0.05543431       -0.44322979        2.33439087
  0         0.46856544       -0.21813557        3.43446803
  0        -0.68380720       -0.02499934        3.17707692
  0         0.15980746        0.68636470        2.71470261

@aug-cc-pVTZ.gbs/N

