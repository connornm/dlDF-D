%Chk=CH400142
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.99589590        0.00000000       -0.48391870
  1        -0.66687370       -0.70891199       -0.52792917
  1         0.10380864       -0.30917554        1.05812108
  1        -0.43283083        1.01808753       -0.04627320
  0        -0.65899893        0.00000000        0.32021611
  0         0.44128011        0.46909746        0.34933848
  0        -0.06869170        0.20458599       -0.70017424
  0         0.28641051       -0.67368345        0.03061966
  6         0.00000000        0.00000000        4.74980174
  1         0.78243865       -0.78343428        4.74726132
  1         0.40575626        0.92540143        5.20254893
  1        -0.31851081        0.20689548        3.70973872
  1        -0.86968410       -0.34886263        5.33965799
  0        -0.51775113        0.51840995        4.75148277
  0        -0.26849487       -0.61235170        4.45021228
  0         0.21076327       -0.13690578        5.43802669
  0         0.57548272        0.23084752        4.35948521

@aug-cc-pVTZ.gbs/N
