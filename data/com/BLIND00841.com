%Chk=BLIND00841
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.45478258       -0.85262447        0.32758698
  6        -0.52551830        0.22567555        1.25391963
 16        -0.40308000       -0.68330783       -1.42089103
 16        -0.54078479        1.93396465        0.78905017
  6         0.48720958        0.84284191       -1.57046880
  6         0.42945479        1.87940788       -0.69306386
  6        -0.50513167       -2.05274515        1.00764473
  7        -0.61813084       -0.07719375        2.52575025
 16        -0.66011029       -1.73169361        2.70403328
  6         1.26954185        0.91996678       -2.75920664
  7         1.87257462        0.94566622       -3.75262701
  6         1.15135638        3.08592848       -0.93203682
  7         1.69884116        4.09723885       -1.10067218
  6        -0.47515619       -3.36550011        0.47760293
  7        -0.46113256       -4.44992526        0.05936692
  6         0.01906558       -0.56087694       12.54675147
  6         0.11493266        0.83968617       12.78139831
 16        -0.83870895       -1.30626403       11.20617787
 16        -0.59667155        2.07317157       11.72938997
  6        -0.59642824       -0.08090145        9.94764352
  6        -0.50217719        1.25866115       10.15819715
  6         0.64103498       -1.24662372       13.57057208
  7         0.74299040        1.22795764       13.86428739
 16         1.25240426       -0.10442861       14.72228581
  6        -0.55060388       -0.63021705        8.63339767
  7        -0.54361788       -1.12195470        7.58014637
  6        -0.35463007        2.16567431        9.06755564
  7        -0.26318303        2.94282584        8.20813777
  6         0.76971431       -2.64659265       13.73997414
  7         0.87670134       -3.79293210       13.89983524

@aug-cc-pVTZ.gbs/N

