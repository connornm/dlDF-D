%Chk=CH400292
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.94241145        0.00000000        0.58124575
  1        -0.30151086       -1.04364111       -0.21422040
  1        -0.79490947        0.50099056        0.58575878
  1         0.15400889        0.54265056       -0.95278413
  0        -0.62360748        0.00000000       -0.38461884
  0         0.19951416        0.69059263        0.14175278
  0         0.52600326       -0.33151280       -0.38760518
  0        -0.10190994       -0.35907983        0.63047123
  6         0.00000000        0.00000000        4.28525708
  1         1.04280794       -0.03298856        3.91451549
  1        -0.61453276       -0.73028045        3.72398539
  1        -0.01685704       -0.25435307        5.36275711
  1        -0.41141813        1.01762207        4.13977032
  0        -0.69004131        0.02182901        4.53058222
  0         0.40664534        0.48323728        4.65665879
  0         0.01115455        0.16830916        3.57225950
  0         0.27224141       -0.67337545        4.38152780

@aug-cc-pVTZ.gbs/N

