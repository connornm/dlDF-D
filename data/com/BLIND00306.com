%Chk=BLIND00306
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.63956416       -0.76225379        0.22591094
  6         0.67882332       -0.95416446        0.72676163
 16        -1.41923889        0.79628762       -0.00118515
 16         1.77100887        0.36445619        1.17727095
  6        -0.01840307        1.79902679       -0.42106041
  6         1.24522708        1.62602333        0.04899566
  6        -1.23820386       -1.98769814        0.01326128
  7         1.08729043       -2.18856946        0.89186446
 16        -0.13284443       -3.24098448        0.47400032
  6        -0.34370208        2.87691859       -1.29479921
  7        -0.66198296        3.75395437       -1.98812347
  6         2.29420356        2.52075203       -0.31576569
  7         3.17219501        3.24119886       -0.56277480
  6        -2.54326675       -2.25732260       -0.46557267
  7        -3.61428718       -2.49703023       -0.84839298
  6        -0.97397695       -0.30359822        3.74811736
  6        -0.08391900       -1.26794246        4.29917425
 16        -0.75806791        1.43842123        3.83337083
 16         1.41949977       -0.86335624        5.14216317
  6         1.01115678        1.54638349        3.78530597
  6         1.87198635        0.63256557        4.30641375
  6        -2.06799982       -0.94187198        3.19946102
  7        -0.42502761       -2.52952264        4.19957311
 16        -1.90165403       -2.65036815        3.44070049
  6         1.48216203        2.74202290        3.16911599
  7         1.81747370        3.74188957        2.68021340
  6         3.28225874        0.83879510        4.25598450
  7         4.43601208        0.97923867        4.26339833
  6        -3.19299534       -0.36588343        2.56087609
  7        -4.12793441        0.09017829        2.04223117

@aug-cc-pVTZ.gbs/N

