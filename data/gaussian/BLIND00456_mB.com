%Chk=BLIND00456_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.72770000       -0.53391217        6.35672728
  6         0.52303311       -1.20941937        6.28473166
 16        -0.93319647        1.20936739        6.44195767
 16         2.09194652       -0.39040437        6.23663327
  6         0.45509855        1.74482436        5.47756935
  6         1.65389990        1.10890887        5.39926941
  6        -1.75112879       -1.45786543        6.42158948
  7         0.50153358       -2.51991000        6.28820999
 16        -1.07066082       -3.05214698        6.41408864
  6         0.22127840        2.97362911        4.79459302
  7        -0.00802224        3.98868610        4.27664171
  6         2.72474585        1.65107716        4.62901982
  7         3.63337525        2.07067858        4.03804317
  6        -3.14420301       -1.22163230        6.51489835
  7        -4.29004453       -1.04633947        6.60101194

@blind-aug.gbs/N

