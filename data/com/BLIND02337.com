%Chk=BLIND02337
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.30975587       -0.16442133       -0.95818883
  6         0.28686529       -1.30271237       -0.34649989
 16         0.10271008        1.50707840       -0.60537568
 16         1.53390163       -1.21623317        0.90729836
  6         0.51868153        1.38018849        1.11363240
  6         1.08739004        0.29932295        1.71024314
  6        -1.21000172       -0.57349218       -1.92126810
  7        -0.09442951       -2.48340922       -0.76890876
 16        -1.20933520       -2.30548701       -1.99212006
  6         0.23873842        2.57011571        1.84629518
  7         0.01510963        3.56974499        2.39577174
  6         1.42649305        0.31779947        3.09538201
  7         1.74589148        0.29890831        4.21276910
  6        -2.01095255        0.22617790       -2.77225222
  7        -2.66867291        0.86698499       -3.48489948
  6         0.30975587        0.16442133       11.71222883
  6        -0.28686526        1.30271237       11.10053988
 16        -0.10271011       -1.50707840       11.35941568
 16        -1.53390160        1.21623320        9.84674163
  6        -0.51868156       -1.38018849        9.64040760
  6        -1.08739005       -0.29932294        9.04379686
  6         1.21000173        0.57349216       12.67530810
  7         0.09442955        2.48340922       11.52294875
 16         1.20933525        2.30548700       12.74616005
  6        -0.23873847       -2.57011571        8.90774483
  7        -0.01510970       -3.56974499        8.35826827
  6        -1.42649305       -0.31779946        7.65865798
  7        -1.74589149       -0.29890829        6.54127090
  6         2.01095254       -0.22617793       13.52629222
  7         2.66867289       -0.86698503       14.23893948

@aug-cc-pVTZ.gbs/N

