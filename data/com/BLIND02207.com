%Chk=BLIND02207
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.42383922        0.14276814       -0.91710767
  6         0.44142249       -1.21276937       -0.48347702
 16         0.49218022        1.54479129        0.14033133
 16         0.50355676       -1.71085939        1.21435513
  6        -0.37212014        0.92162873        1.55764406
  6        -0.36455131       -0.36977578        1.98174447
  6         0.41718510        0.18868201       -2.29665501
  7         0.44496490       -2.14158238       -1.40822231
 16         0.45983406       -1.43078490       -2.91336773
  6        -1.07143600        1.93429300        2.27636380
  7        -1.60556102        2.79813308        2.84176967
  6        -1.05769013       -0.75874621        3.16589079
  7        -1.58334507       -1.11478785        4.13946532
  6         0.41131586        1.32739339       -3.13836184
  7         0.41585238        2.25078724       -3.84434989
  6        -0.65293215       -0.26834751       11.17662477
  6         0.48182926        0.44804677       11.65088464
 16        -1.38815157       -0.06663567        9.59318504
 16         1.36596891        1.65008383       10.69799834
  6         0.04213814        0.32670391        8.62175445
  6         1.13217979        1.00880771        9.06242172
  6        -1.13129021       -1.09334446       12.17455103
  7         0.86734525        0.21979330       12.88260612
 16        -0.15703073       -0.88343754       13.59285428
  6        -0.06339979       -0.10874957        7.26885694
  7        -0.20561767       -0.46176953        6.17052001
  6         2.21258973        1.31516060        8.18332889
  7         3.10172199        1.61185398        7.49605676
  6        -2.24490697       -1.96708135       12.13442611
  7        -3.16103222       -2.68238869       12.12234912

@aug-cc-pVTZ.gbs/N

