%Chk=BLIND00522
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.66298680        0.40608147       -0.66079895
  6        -0.67121913       -1.01710758       -0.64371611
 16        -0.28880902        1.43278085        0.71557164
 16        -0.26956840       -1.99487519        0.77657806
  6         0.89895335        0.42132420        1.55836348
  6         0.90310688       -0.93774878        1.58176108
  6        -1.06290238        0.85763910       -1.90235702
  7        -1.03219796       -1.63117106       -1.74393443
 16        -1.43082899       -0.50737925       -2.90544300
  6         1.86750637        1.17855469        2.27919404
  7         2.62020015        1.83839753        2.87016994
  6         1.87962984       -1.65710593        2.33200779
  7         2.63723946       -2.28330272        2.95235557
  6        -1.20645414        2.19424186       -2.34739005
  7        -1.33866994        3.28500534       -2.72666558
  6        -0.15745415        0.79494481        5.90390572
  6        -1.27161072        0.48784307        5.07314082
 16         1.01039414       -0.37329807        6.50378533
 16        -1.64690186       -1.13631087        4.47659334
  6         1.03728299       -1.52363539        5.15466595
  6        -0.01856421       -1.82504132        4.35342615
  6        -0.17418322        2.13723406        6.22528232
  7        -2.08206847        1.47039765        4.76394062
 16        -1.56546377        2.87224535        5.49801112
  6         2.29723184       -2.17051773        4.99604124
  7         3.33511675       -2.68859450        4.92131248
  6         0.09676964       -2.80363800        3.32241278
  7         0.13686658       -3.61406466        2.49023071
  6         0.73654648        2.85886484        7.03459616
  7         1.46887967        3.46499466        7.70346883

@aug-cc-pVTZ.gbs/N

