%Chk=BLIND01088
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.59625948        0.01065580       -0.82793115
  6        -0.04522126       -1.26881822       -0.53616111
 16         0.17120497        1.54971438       -0.46624948
 16         1.50735674       -1.51863797        0.27727603
  6         1.03934741        1.13454909        1.02307539
  6         1.56920955       -0.08283196        1.31440601
  6        -1.79118489       -0.14139631       -1.50197614
  7        -0.72220575       -2.32043667       -0.92814410
 16        -2.09760876       -1.83240042       -1.72877139
  6         1.18406347        2.23860742        1.91257082
  7         1.29743811        3.17367801        2.59373322
  6         2.29241228       -0.29933572        2.52438365
  7         2.91345317       -0.51313937        3.48330092
  6        -2.66008652        0.86508637       -1.98897819
  7        -3.38023175        1.67791335       -2.40346470
  6        -0.39935742       -0.33545547        9.28165163
  6        -1.25368847        0.54939723        8.56542486
 16         0.76878154       -1.42678037        8.55159820
 16        -1.25461632        0.71827344        6.80302073
  6         1.24928770       -0.48071214        7.13096778
  6         0.44414985        0.36809127        6.43892018
  6        -0.67766312       -0.25814819       10.63141929
  7        -2.11202638        1.25396621        9.26162785
 16        -1.97199305        0.86888231       10.87501496
  6         2.59370398       -0.72395031        6.72541914
  7         3.68706219       -0.97300080        6.41925789
  6         0.91726459        1.04719543        5.27748920
  7         1.25089923        1.59841353        4.31013577
  6        -0.07498473       -0.97190955       11.69560533
  7         0.40204235       -1.55566462       12.58034885

@aug-cc-pVTZ.gbs/N
