%Chk=BLIND00740_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.21682635       -0.91091966        6.33862179
  6         0.10898117       -0.96374188        4.95410530
 16         0.41867800        0.26063811        7.48409255
 16         1.19605275        0.17650046        4.14619991
  6         0.60183207        1.67480735        6.43019407
  6         0.91147595        1.63820470        5.10715780
  6        -1.05140895       -1.96250859        6.65945795
  7        -0.40750604       -1.94120537        4.25006970
 16        -1.32555731       -2.91656602        5.23838664
  6         0.43547330        2.91125215        7.11911981
  7         0.31484170        3.89361664        7.72871068
  6         1.08263762        2.84004906        4.35868219
  7         1.26008871        3.79293934        3.71724836
  6        -1.60437654       -2.29287303        7.92051093
  7        -2.05941426       -2.58419646        8.94966944

@blind-aug.gbs/N

