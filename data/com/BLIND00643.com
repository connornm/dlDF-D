%Chk=BLIND00643
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -1.01864614        0.02001890        0.05539122
  6        -0.43860395        1.27808764       -0.27116772
 16        -0.34553362       -1.54629922       -0.37132147
 16         1.09029147        1.46958666       -1.14317417
  6         1.38824814       -1.18059958       -0.30508960
  6         1.95469047        0.01603197       -0.61308323
  6        -2.22567317        0.21754196        0.69521433
  7        -1.10271176        2.35509835        0.07069305
 16        -2.53229383        1.92016652        0.80420944
  6         2.18273470       -2.30261011        0.07037427
  7         2.78989636       -3.25145291        0.35715834
  6         3.36938548        0.19189556       -0.57243889
  7         4.51759421        0.37207851       -0.58121155
  6        -3.13340527       -0.75531073        1.17972288
  7        -3.89425034       -1.54010161        1.57508849
  6        -0.77984589       -0.43715453        5.69852388
  6        -1.18300919        0.47242478        4.68070686
 16         0.69202919       -1.39708102        5.67319521
 16        -0.23888669        0.81132125        3.22181337
  6         1.77953645       -0.29632919        4.80728792
  6         1.40751336        0.57734644        3.83469665
  6        -1.75447846       -0.50420840        6.67366404
  7        -2.34263293        1.06550270        4.82696485
 16        -3.06970499        0.53257133        6.22634588
  6         3.14336487       -0.43640323        5.19644208
  7         4.24433795       -0.60355297        5.52977713
  6         2.37219362        1.38718481        3.16574541
  7         3.12664615        2.04457064        2.57449979
  6        -1.76642209       -1.28913815        7.85218476
  7        -1.79726532       -1.93323237        8.81928726

@aug-cc-pVTZ.gbs/N
