%Chk=BLIND01101
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.41373963        0.91355704       -0.18799459
  6        -0.47956831        0.47015703       -1.20348659
 16         1.39735775       -0.13457167        0.82316231
 16        -0.77120911       -1.23172948       -1.59478733
  6         0.33339815       -1.54703662        0.95375642
  6        -0.52451638       -1.98034602       -0.00742497
  6         0.44378321        2.29338853       -0.16641537
  7        -1.09414790        1.38819916       -1.90871748
 16        -0.60377199        2.89810831       -1.40809514
  6         0.48663244       -2.25110562        2.18334812
  7         0.66396444       -2.80467884        3.18997768
  6        -1.30395872       -3.15981465        0.18036570
  7        -1.93606752       -4.13026887        0.27825952
  6         1.20233109        3.14316896        0.67480206
  7         1.82444036        3.85560362        1.35045939
  6         0.05725681       -0.81837882        8.14766722
  6         0.35200742       -1.09601128        6.78316403
 16         0.36516390        0.70976802        8.95913325
 16         1.05586205        0.07870147        5.66101121
  6         0.11900543        1.84232392        7.61717751
  6         0.39514525        1.58955087        6.31046536
  6        -0.44072441       -1.95340546        8.75514763
  7         0.11422941       -2.30993880        6.34991859
 16        -0.46948018       -3.24385986        7.59811387
  6        -0.37146581        3.11203672        8.03950086
  7        -0.74860529        4.13811481        8.43462438
  6         0.20429892        2.58977399        5.31192978
  7         0.09190537        3.38306437        4.46991778
  6        -0.84877379       -2.12351187       10.10040892
  7        -1.17929546       -2.28366479       11.20322458

@blind-aug.gbs/N

