%Chk=CH400015
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.90759551        0.00000000        0.63423682
  1         0.28839751        0.15965798       -1.05703461
  1        -0.51879564       -0.97325272        0.09806212
  1        -0.67719738        0.81359475        0.32473567
  0        -0.60056926        0.00000000       -0.41968381
  0        -0.19083686       -0.10564803        0.69945531
  0         0.34329468        0.64401560       -0.06488914
  0         0.44811144       -0.53836758       -0.21488235
  6         0.00000000        0.00000000        3.04308537
  1        -0.11002166        1.00958239        2.60192260
  1         0.90191450       -0.02606113        3.68484625
  1        -0.89371458       -0.23577711        3.65273229
  1         0.10182174       -0.74774415        2.23284037
  0         0.07280294       -0.66805547        3.33500925
  0        -0.59681005        0.01724503        2.61842279
  0         0.59138405        0.15601717        2.63967307
  0        -0.06737694        0.49479327        3.57923638

@aug-cc-pVTZ.gbs/N

