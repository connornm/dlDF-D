%Chk=BLIND02262
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(fermi,maxcyc=120)

M05 opt

0 1
  6         0.30559928       -0.91886280        0.32157291
  6        -0.35873188       -0.15449231        1.32168716
 16         1.26660990       -0.25224837       -0.99014369
 16        -0.35013004        1.61473812        1.38755086
  6         0.40206281        1.26803607       -1.28240320
  6        -0.23696938        2.00677231       -0.33712509
  6         0.13980897       -2.26469523        0.57954940
  7        -0.98504942       -0.81243486        2.26651411
 16        -0.78831911       -2.44865288        2.03198336
  6         0.46992677        1.69828382       -2.63937446
  7         0.57312073        2.01982338       -3.75164976
  6        -0.86368316        3.24374898       -0.67038249
  7        -1.36202085        4.27043183       -0.89062783
  6         0.64114029       -3.36841277       -0.15232692
  7         1.05222902       -4.28600061       -0.73553746
  6-Bq        0.50370994      -0.65857279       7.17679915
  6-Bq       -0.90269157      -0.87523444       7.14648106
 16-Bq        1.45531605      -0.00075070       5.85383561
 16-Bq       -1.93780768      -0.48213489       5.76496171
  6-Bq        0.27417634       1.10296496       5.12524600
  6-Bq       -1.07068568       0.90859618       5.09044720
  6-Bq        1.01705768      -1.10979580       8.37603854
  7-Bq       -1.44985406      -1.43142268       8.19963246
 16-Bq       -0.27235761      -1.77403632       9.32525891
  6-Bq        0.87540936       2.24008914       4.51160350
  7-Bq        1.41304633       3.13741731       4.00470937
  6-Bq       -1.93113638       1.83834734       4.43536759
  7-Bq       -2.66628191       2.55344619       3.88844820
  6-Bq        2.36290034      -1.09816885       8.81619850
  7-Bq        3.46353126      -1.10511911       9.18990855

@aug-cc-pVTZ.gbs/N
