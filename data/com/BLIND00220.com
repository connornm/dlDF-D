%Chk=BLIND00220
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.84351329       -0.42006603        0.39132950
  6        -0.86007347        1.00315297        0.39162625
 16        -0.01445090       -1.42249451       -0.79032584
 16        -0.01084943        2.00575324       -0.79504960
  6         1.37821787       -0.39153033       -1.16676795
  6         1.37610806        0.96774678       -1.16907472
  6        -1.63592307       -0.89363753        1.41754584
  7        -1.57888288        1.59752980        1.31243526
 16        -2.33605263        0.45311921        2.25473843
  6         2.54147866       -1.13208478       -1.52661944
  7         3.45665662       -1.77846078       -1.83620044
  6         2.54168329        1.70408387       -1.53401856
  7         3.45824307        2.34406244       -1.85225075
  6        -1.90809516       -2.23796383        1.76950832
  7        -2.14983757       -3.33536684        2.06676776
  6         0.98271852       -0.26912958        7.55936144
  6         0.29535323       -0.83388884        8.67039816
 16         0.25154434        0.12519687        6.01068213
 16        -1.43571995       -1.20515443        8.68205304
  6        -1.38097673        0.58050681        6.53192093
  6        -2.04790552        0.05012675        7.59095210
  6         2.32454351       -0.14466762        7.85822283
  7         1.00428439       -1.12372379        9.73401103
 16         2.60386154       -0.75639878        9.45590647
  6        -1.98485391        1.55856738        5.68941304
  7        -2.44393192        2.33885332        4.96030792
  6        -3.38091916        0.45403902        7.89733657
  7        -4.47708394        0.73155018        8.16629097
  6         3.36495312        0.35644496        7.03870928
  7         4.23277418        0.75594901        6.37660730

@aug-cc-pVTZ.gbs/N
