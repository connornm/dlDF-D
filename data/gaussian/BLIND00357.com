%Chk=BLIND00357
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.34951222       -0.62725072       -0.72491837
  6        -0.79130500       -1.12779091       -0.03657091
 16         0.69692682        1.07234722       -1.00631945
 16        -2.04365801       -0.10454291        0.68395932
  6        -0.00091834        1.79029589        0.45712988
  6        -1.08862827        1.32197740        1.12436968
  6         1.10429599       -1.68086535       -1.19978773
  7        -0.93024172       -2.42919542        0.03344985
 16         0.32261906       -3.16873476       -0.77546198
  6         0.67072886        2.97633567        0.87341023
  7         1.23546589        3.95118578        1.15964853
  6        -1.59913108        2.00372837        2.26828594
  7        -2.06651018        2.54213508        3.18621902
  6         2.30729821       -1.63479066       -1.94532862
  7         3.29013768       -1.61573479       -2.56561837
  6        -0.60661765       -0.52782994        6.12634123
  6        -0.34756996        0.80145742        6.56419514
 16        -0.67194666       -1.04625781        4.44840163
 16        -0.01361685        2.15744853        5.47591586
  6         0.52629617        0.04404382        3.72758365
  6         0.78409763        1.31446478        4.13650059
  6        -0.86929781       -1.33267403        7.21653771
  7        -0.39438259        1.03136978        7.85369462
 16        -0.79682805       -0.37025480        8.65644184
  6         1.19770555       -0.52503761        2.60667089
  7         1.70074110       -1.02284679        1.68453517
  6         1.73826863        2.12663046        3.45540842
  7         2.48498345        2.83337983        2.91333052
  6        -1.18643167       -2.71262808        7.23446427
  7        -1.45708935       -3.84256731        7.26729232

@blind-aug.gbs/N

