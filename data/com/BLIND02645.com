%Chk=BLIND02645
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.28398023       -0.55054534       -0.81077985
  6        -0.53928699       -1.19172428        0.43402322
 16        -0.29087128        1.18693452       -1.07466029
 16        -0.87774223       -0.33003444        1.94317099
  6         0.29893993        1.75393416        0.49847551
  6         0.06338231        1.15085343        1.69365344
  6        -0.10089372       -1.50194965       -1.79394433
  7        -0.55502996       -2.50225443        0.44511063
 16        -0.28710276       -3.07684176       -1.09403891
  6         1.03545500        2.97111844        0.41422800
  7         1.60826192        3.97589342        0.29798476
  6         0.54560363        1.71711935        2.91051311
  7         0.89366822        2.15761184        3.92822797
  6         0.16376207       -1.30365296       -3.17082638
  7         0.37213010       -1.15949223       -4.30524093
  6         0.14629928        0.52415018       13.01761602
  6        -1.04023416        0.88858045       12.32108488
 16         0.87590216       -1.07462363       13.00672716
 16        -1.96626252       -0.21583017       11.29281595
  6         0.45608869       -1.60905196       11.36896494
  6        -0.67243792       -1.26779593       10.69249169
  6         0.59693756        1.59720722       13.75975950
  7        -1.48407098        2.11091132       12.48471908
 16        -0.49431515        2.92691532       13.54558559
  6         1.41608154       -2.50065309       10.80811665
  7         2.21423380       -3.24232056       10.40305446
  6        -0.93768595       -1.79246330        9.39314202
  7        -1.20855195       -2.21817494        8.34608551
  6         1.73107981        1.66740339       14.60472240
  7         2.65254915        1.74011050       15.30949460

@aug-cc-pVTZ.gbs/N

