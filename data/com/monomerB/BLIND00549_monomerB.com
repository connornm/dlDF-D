%Chk=BLIND00549
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=130)

M05 opt

0 1
  6-Bq       -0.88189770       0.51185145       0.03706000
  6-Bq       -1.04266206      -0.90074999      -0.03030962
 16-Bq        0.47133827       1.34181702       0.79100045
 16-Bq        0.12546372      -2.06438379       0.61469183
  6-Bq        1.77728714       0.18620176       0.46979493
  6-Bq        1.63716009      -1.16412697       0.40181857
  6-Bq       -1.98445481       1.13133441      -0.51603302
  7-Bq       -2.14336581      -1.35690133      -0.57642062
 16-Bq       -3.10962447      -0.08115487      -1.03438574
  6-Bq        3.05162303       0.80786654       0.32595272
  7-Bq        4.07206020       1.35786492       0.23991326
  6-Bq        2.76448782      -2.01018857       0.18443897
  7-Bq        3.65811974      -2.73890325       0.03841483
  6-Bq       -2.24190338       2.51789568      -0.64384745
  7-Bq       -2.47388706       3.65204871      -0.74861569
  6         0.98080285       -0.22435426        4.54927492
  6         0.43850304        0.70875128        5.47720359
 16         0.07338392       -1.49733281        3.74638442
 16        -1.26591320        0.75801615        5.95378158
  6        -1.51007494       -0.71475844        3.59031568
  6        -2.03895558        0.17822982        4.46809771
  6         2.34118999       -0.02320185        4.43016261
  7         1.26786085        1.55845781        6.03220311
 16         2.81323565        1.26804295        5.48594968
  6        -2.23764366       -1.14790495        2.44396603
  7        -2.80165192       -1.54931146        1.51017776
  6        -3.34590963        0.71473766        4.27323266
  7        -4.41673697        1.15475081        4.17006597
  6         3.26790845       -0.73267874        3.62828738
  7         4.04390365       -1.31001538        2.98360935

@aug-cc-pVTZ.gbs/N

