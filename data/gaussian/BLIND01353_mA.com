%Chk=BLIND01353_mA
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.16508417       -0.33613073        0.94914291
  6         0.21565769        1.08104815        0.82718803
 16         0.60747017       -1.47709898       -0.31220666
 16         0.70821738        1.93354530       -0.64428127
  6         0.12625960       -0.56505244       -1.75468008
  6         0.16862911        0.78742411       -1.88373841
  6        -0.21551912       -0.67834831        2.23106782
  7        -0.09184295        1.78869876        1.88668508
 16        -0.44450717        0.77020545        3.15523395
  6        -0.28516574       -1.40058129       -2.83343417
  7        -0.58924913       -2.12345532       -3.69144318
  6        -0.19736503        1.42091526       -3.10794395
  7        -0.45809249        1.97724446       -4.09457718
  6        -0.37996307       -1.97173486        2.78358737
  7        -0.51079630       -3.02558761        3.25621087

@blind-aug.gbs/N
