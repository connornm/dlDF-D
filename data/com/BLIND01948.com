%Chk=BLIND01948
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.12787334        0.02630716       -1.01196109
  6         0.32564232       -1.26354791       -0.44364077
 16         0.45111446        1.55079869       -0.19958348
 16         0.90627837       -1.54302746        1.20540177
  6         0.06940191        1.11056338        1.47511365
  6         0.25229694       -0.11689041        2.02968613
  6        -0.28810190       -0.10247267       -2.32180214
  7         0.09181362       -2.30175954       -1.20869030
 16        -0.36664190       -1.78682708       -2.72380745
  6        -0.42506387        2.20498988        2.24246226
  7        -0.80138186        3.13263387        2.83327053
  6        -0.04409211       -0.35404424        3.40440557
  7        -0.24415017       -0.58491650        4.52583364
  6        -0.58688506        0.92074641       -3.25396081
  7        -0.82722776        1.74753706       -4.03482477
  6        -0.12557261       -0.13659642        9.57194464
  6         0.69198888       -1.12658744       10.18622280
 16        -1.44431124        0.71793874       10.35884663
 16         0.53057863       -1.64011571       11.87288334
  6        -0.83064762        0.79707675       12.02055545
  6        -0.04875631       -0.14032123       12.61853126
  6         0.22790512       -0.00107396        8.24454479
  7         1.59544634       -1.71142809        9.43816795
 16         1.50399560       -1.11802822        7.88561191
  6        -1.27775350        1.95081857       12.72770876
  7        -1.68706994        2.89395481       13.27006206
  6         0.35269399       -0.00290597       13.98010999
  7         0.67900713        0.05326753       15.09424212
  6        -0.33280315        0.85838879        7.26878946
  7        -0.78713125        1.55039247        6.45281593

@aug-cc-pVTZ.gbs/N
