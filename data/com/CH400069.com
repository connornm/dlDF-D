%Chk=CH400069
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.87892184        0.00000000        0.67341098
  1        -0.67634297        0.83212270        0.27589485
  1         0.33699439        0.12983501       -1.04669175
  1        -0.53957327       -0.96195771        0.09738592
  0        -0.58159548        0.00000000       -0.44560592
  0         0.44754606       -0.55062779       -0.18256367
  0        -0.22299413       -0.08591373        0.69261128
  0         0.35704354        0.63654153       -0.06444169
  6         0.00000000        0.00000000        3.71523752
  1        -0.44609302       -0.33703437        4.67095443
  1        -0.21921020        1.07478571        3.56434514
  1         1.09628819       -0.15236595        3.74561467
  1        -0.43098497       -0.58538540        2.88003584
  0         0.29518629        0.22302059        3.08282562
  0         0.14505461       -0.71120147        3.81508522
  0        -0.72542998        0.10082278        3.69513651
  0         0.28518908        0.38735810        4.26790272

@aug-cc-pVTZ.gbs/N
