%Chk=BLIND01626
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.37376859        0.36027394       -0.87841257
  6        -0.39303393       -1.06004803       -0.78818517
 16        -0.52516633        1.46310376        0.48149984
 16        -0.54247262       -1.95597464        0.73154241
  6         0.26814619        0.51641186        1.75372299
  6         0.25833928       -0.83935815        1.85085970
  6        -0.29007953        0.74008878       -2.20281434
  7        -0.32882097       -1.73628836       -1.90909362
 16        -0.27069073       -0.68105718       -3.19518262
  6         0.90871504        1.32584918        2.73635845
  7         1.39552486        2.02775787        3.52473235
  6         0.89000927       -1.50265564        2.94400718
  7         1.36581093       -2.08326917        3.83137792
  6        -0.25585071        2.04922701       -2.74147370
  7        -0.23605191        3.11649407       -3.20152255
  6        -0.16784279       -0.34560394       11.16021898
  6        -0.20809255        1.07318962       11.05424046
 16        -0.61806279       -1.46903719        9.88596770
 16        -0.69395176        1.94577579        9.59235040
  6        -0.12975535       -0.54430129        8.45399652
  6        -0.16226827        0.80981630        8.34018748
  6         0.20982954       -0.70500676       12.43830209
  7         0.10417124        1.76659903       12.12173039
 16         0.44900386        0.73133488       13.37880456
  6         0.27596680       -1.37057812        7.36599569
  7         0.57509157       -2.08592008        6.49996840
  6         0.20873040        1.45439710        7.12329902
  7         0.47382418        2.01990827        6.14307020
  6         0.36469754       -2.00570017       12.97624353
  7         0.48772069       -3.06573655       13.43698651

@blind-aug.gbs/N
