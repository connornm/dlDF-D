%Chk=BLIND01075_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.04022668        0.69824557       11.98997679
  6         0.44737626        1.19159650       10.71849704
 16        -1.14685478       -0.57121257       12.25041923
 16        -0.16077422        0.57589284        9.17394198
  6        -0.86584621       -1.57421785       10.81545000
  6        -0.47611867       -1.11615931        9.59645908
  6         0.66587182        1.41688678       12.98868906
  7         1.30113217        2.18593796       10.70281481
 16         1.66374060        2.63511614       12.26411907
  6        -1.13558504       -2.95595493       11.03711019
  7        -1.39074771       -4.06800904       11.25941785
  6        -0.32370594       -2.00484788        8.49149531
  7        -0.21523999       -2.68744275        7.55702995
  6         0.53918523        1.26688323       14.39105060
  7         0.44031918        1.16401424       15.54462738

@blind-aug.gbs/N

