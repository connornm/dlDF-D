%Chk=BLIND01372_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.23254613        0.25896634        6.69890006
  6         1.36900449       -0.03621473        7.50336086
 16        -1.35832326       -0.46278992        6.89052242
 16         1.35448235       -1.16873581        8.86415969
  6        -1.39339836       -0.70330780        8.64695146
  6        -0.31604550       -0.98445668        9.42663740
  6         0.58758243        1.14407066        5.70098083
  7         2.49928130        0.54732584        7.18742610
 16         2.27877504        1.49844892        5.83922034
  6        -2.70542249       -0.61923227        9.19705476
  7        -3.79615496       -0.56517520        9.59526768
  6        -0.46071982       -1.20722772       10.82784558
  7        -0.53564349       -1.42458468       11.96717368
  6        -0.22684250        1.68966575        4.67911661
  7        -0.87830109        2.13988280        3.82823340

@blind-aug.gbs/N

