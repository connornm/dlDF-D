%Chk=BLIND02099_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.27477679       -0.03611518       11.07265896
  6         0.19944383       -1.27846054       10.56519521
 16         0.25095256        1.54846138       10.52380341
 16         1.38689289       -1.42929551        9.26066392
  6         0.57140386        1.20694453        8.81359895
  6         1.02310911        0.02543010        8.31597454
  6        -1.15660550       -0.26461581       12.10970708
  7        -0.25303666       -2.37198886       11.12850440
 16        -1.28950013       -1.97345777       12.36859848
  6         0.35113142        2.33242872        7.96751676
  7         0.18079441        3.28389964        7.32185323
  6         1.29477666       -0.13046201        6.92463118
  7         1.55690355       -0.29355342        5.80408809
  6        -1.84940937        0.68450688       12.89986993
  7        -2.41890313        1.44964140       13.56420565

@blind-aug.gbs/N

