%Chk=BLIND01217
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.45538199        0.91227994       -0.03848980
  6         0.74023932        0.08289648       -1.15956341
 16         0.12035513        0.33625172        1.58769135
 16         0.77185716       -1.68637984       -1.10255485
  6        -0.69648329       -1.19519505        1.22496847
  6        -0.43492781       -1.99569968        0.15799842
  6         0.55261980        2.23916947       -0.40620192
  7         1.02811883        0.67736201       -2.29163944
 16         1.00771563        2.32787983       -2.07648478
  6        -1.66177661       -1.56083432        2.20762121
  7        -2.42065479       -1.82798899        3.04662378
  6        -1.11874256       -3.23502955       -0.01635985
  7        -1.62989325       -4.26555766       -0.18268156
  6         0.36141882        3.39010591        0.39626673
  7         0.21733554        4.34562250        1.04227289
  6-Bq       -0.85872511       0.49143079       8.26686142
  6-Bq       -0.83851100      -0.87057957       7.85413168
 16-Bq       -0.04170290       1.13163170       9.68500603
 16-Bq        0.05052809      -2.14952503       8.69585393
  6-Bq        1.38155940       0.07554245       9.74088092
  6-Bq        1.41459725      -1.22563951       9.34911817
  6-Bq       -1.67537704       1.21937112       7.42514549
  7-Bq       -1.55273399      -1.19282425       6.80346337
 16-Bq       -2.35048431       0.15361549       6.23635549
  6-Bq        2.52953755       0.71264287      10.29526243
  7-Bq        3.43136293       1.26720167      10.77522968
  6-Bq        2.60310792      -2.00292767       9.48029530
  7-Bq        3.53965351      -2.68152301       9.59569402
  6-Bq       -1.98639972       2.59978035       7.47903151
  7-Bq       -2.25996824       3.72896702       7.51358054

@aug-cc-pVTZ.gbs/N

