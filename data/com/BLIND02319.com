%Chk=BLIND02319
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.60677168        0.43632671       -0.69466254
  6         0.62816072       -0.98640721       -0.72926559
 16         0.33105070        1.40990104        0.74222562
 16         0.34514258       -2.01778805        0.68167121
  6        -0.77987267        0.35641698        1.63673257
  6        -0.77088111       -1.00260287        1.61166768
  6         0.90619926        0.93556365       -1.94622509
  7         0.90852832       -1.55745977       -1.87519448
 16         1.20722873       -0.38926643       -3.02280908
  6        -1.69637556        1.07712739        2.45629868
  7        -2.40684533        1.70746054        3.12645228
  6        -1.68075449       -1.75860082        2.40824755
  7        -2.38314205       -2.41452166        3.06196438
  6         1.00393532        2.28846396       -2.35276764
  7         1.09747859        3.39325889       -2.70172858
  6        -0.60677175       -0.43632660        7.19056254
  6        -0.62816055        0.98640732        7.22516559
 16        -0.33105094       -1.40990099        5.75367438
 16        -0.34514223        2.01778811        5.81422879
  6         0.77987261       -0.35641712        4.85916743
  6         0.77088129        1.00260274        4.88423231
  6        -0.90619942       -0.93556349        8.44212509
  7        -0.90852804        1.55745992        8.37109447
 16        -1.20722866        0.38926664        9.51870908
  6         1.69637537       -1.07712768        4.03960132
  7         2.40684503       -1.70746096        3.36944772
  6         1.68075479        1.75860052        4.08765245
  7         2.38314247        2.41452124        3.43393562
  6        -1.00393572       -2.28846379        8.84866763
  7        -1.09747918       -3.39325870        9.19762858

@aug-cc-pVTZ.gbs/N
