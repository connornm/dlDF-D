%Chk=BLIND01899
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.77292339        0.53334923       -0.39904505
  6        -0.89956476       -0.87378315       -0.57157643
 16         0.02998482        1.31929047        0.95230665
 16        -0.23231217       -2.07459569        0.54529937
  6         1.30243348        0.13541274        1.30363275
  6         1.19488355       -1.20998292        1.14242617
  6        -1.43806044        1.18740077       -1.41644983
  7        -1.58305816       -1.29507348       -1.60753553
 16        -2.16822231        0.01033326       -2.45866011
  6         2.48552940        0.72662842        1.83451666
  7         3.41904499        1.25185001        2.28605177
  6         2.26548512       -2.08126948        1.50105406
  7         3.10003849       -2.83029386        1.80668334
  6        -1.57505298        2.58152196       -1.62339812
  7        -1.70549029        3.72239971       -1.80362018
  6        -0.28752303       -0.14608806       10.41676806
  6        -0.52318130        1.17526495        9.94314581
 16        -0.31380488       -1.59664499        9.42494626
 16        -0.85340086        1.57890491        8.25121845
  6         0.28764314       -0.94718394        7.88852747
  6         0.07083046        0.31218554        7.42525127
  6        -0.10608837       -0.12324870       11.78492853
  7        -0.52633811        2.14152771       10.82868878
 16        -0.26969574        1.50846771       12.34667646
  6         1.01158986       -1.90541359        7.12120306
  7         1.57351805       -2.72775850        6.52188796
  6         0.56083613        0.72080748        6.14976529
  7         0.91573170        1.08755922        5.10551033
  6         0.14169516       -1.21088307       12.65718035
  7         0.33634589       -2.09279307       13.38891922

@aug-cc-pVTZ.gbs/N

