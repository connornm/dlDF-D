%Chk=BLIND02581
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.52924528        0.78551058        0.37944890
  6         0.23212087        0.91215000       -1.00673514
 16         0.77755541       -0.73399621        1.22677281
 16         0.03076161       -0.45953846       -2.10787708
  6        -0.32544603       -1.78958284        0.32497086
  6        -0.61851053       -1.67896828       -0.99772399
  6         0.66332977        2.03895864        0.94175889
  7         0.13655448        2.12492249       -1.49449531
 16         0.43747618        3.23570454       -0.29177910
  6        -0.87167294       -2.83827914        1.12055988
  7        -1.27223596       -3.68952141        1.80329621
  6        -1.48529164       -2.61231216       -1.63920793
  7        -2.16212023       -3.36604210       -2.20905035
  6         0.96705459        2.37243417        2.28403820
  7         1.22464206        2.66436040        3.37926589
  6        -0.52924515       -0.78551066        9.74420566
  6        -0.23212088       -0.91214998       11.13038974
 16        -0.77755536        0.73399606        8.89688165
 16        -0.03076188        0.45953856       12.23153163
  6         0.32544586        1.78958286        9.79868367
  6         0.61851024        1.67896840       11.12137855
  6        -0.66332944       -2.03895876        9.18189572
  7        -0.13655440       -2.12492243       11.61814998
 16        -0.43747585       -3.23570457       10.41543380
  6         0.87167274        2.83827918        9.00309465
  7         1.27223574        3.68952146        8.32035832
  6         1.48529117        2.61231241       11.76286254
  7         2.16211961        3.36604245       12.33270499
  6        -0.96705408       -2.37243439        7.83961639
  7        -1.22464140       -2.66436071        6.74438869

@aug-cc-pVTZ.gbs/N
