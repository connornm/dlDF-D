%Chk=BLIND00539
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.46072808        0.88490531        0.21396526
  6        -1.32999270       -0.09109509       -0.34962048
 16         0.94001513        0.53494405        1.21594441
 16        -1.11886881       -1.83823178       -0.15589790
  6         1.46505802       -1.00420725        0.50921580
  6         0.64598356       -1.94385434       -0.03283993
  6        -0.92580311        2.14934881       -0.08635303
  7        -2.36662439        0.34255741       -1.02430580
 16        -2.39016166        2.00686131       -1.00296523
  6         2.87262881       -1.20917802        0.59765278
  7         4.02157297       -1.34159916        0.71420748
  6         1.16755835       -3.17328281       -0.53324708
  7         1.54482067       -4.20013972       -0.92588451
  6        -0.37463723        3.39987446        0.28455154
  7         0.05864091        4.43499632        0.58768303
  6        -0.54953891        0.79961879        5.78780457
  6        -1.33243334       -0.35162242        5.49192175
 16         1.10336247        0.77375605        6.38426478
 16        -0.74880954       -2.01576477        5.64877147
  6         1.71077964       -0.68866649        5.58646357
  6         0.97404829       -1.79378651        5.29736026
  6        -1.30910635        1.93732566        5.60345812
  7        -2.57166852       -0.15668090        5.11224558
 16        -2.90608175        1.47398842        5.11374931
  6         3.10538183       -0.63075811        5.29867717
  7         4.24859543       -0.54634442        5.10598569
  6         1.57261062       -2.93934336        4.69456297
  7         2.02385819       -3.90499634        4.23108041
  6        -0.92975651        3.28816325        5.79450771
  7        -0.63848362        4.40193430        5.95503832

@blind-aug.gbs/N
