%Chk=BLIND02594_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -1.00977107       -0.13836699        7.13637723
  6        -0.62751987        1.21540024        6.91952483
 16        -0.10978828       -1.54462439        6.58775373
 16         0.85141094        1.70659493        6.07923931
  6         1.54922577       -0.92857740        6.69857144
  6         1.92737812        0.36108180        6.49504255
  6        -2.23016247       -0.17861664        7.78006041
  7        -1.44512181        2.14801281        7.34334189
 16        -2.78971339        1.44343432        8.02639351
  6         2.50542820       -1.94522406        6.98700884
  7         3.24996291       -2.81217265        7.19975570
  6         3.29949844        0.74415145        6.56201562
  7         4.40733207        1.09540681        6.57784158
  6        -2.97873796       -1.31386417        8.17504929
  7        -3.61084703       -2.23433898        8.49794272

@blind-aug.gbs/N

