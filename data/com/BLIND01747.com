%Chk=BLIND01747
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.42745825        0.71508172       -0.58910651
  6        -0.43443672       -0.58418275       -1.17021734
 16        -0.50242080        1.05382923        1.13372470
 16        -0.48823431       -2.08816439       -0.23762237
  6         0.37024408       -0.34968279        1.77616700
  6         0.37280975       -1.59473917        1.23075700
  6        -0.42510099        1.66949319       -1.58630078
  7        -0.43423795       -0.65939158       -2.47872933
 16        -0.45835346        0.87399634       -3.12607531
  6         1.06462635       -0.06966800        2.98881477
  7         1.59440701        0.20056230        3.98756452
  6         1.07202926       -2.66945967        1.85526468
  7         1.60293489       -3.58037943        2.34443802
  6        -0.42951736        3.07932074       -1.45388018
  7        -0.44246195        4.23813355       -1.36398455
  6        -0.42788722        0.71522505       10.07348909
  6        -0.43630097       -0.58410712        9.49254866
 16        -0.50036752        1.05422202       11.79637745
 16        -0.48960643       -2.08794682       10.42540076
  6         0.37244997       -0.34965143       12.43782192
  6         0.37365561       -1.59477595       11.89256279
  6        -0.42638604        1.66951268        9.07617466
  7        -0.43789555       -0.65947681        8.18404689
 16        -0.46210558        0.87384354        7.53654442
  6         1.06859964       -0.06983825       13.64950252
  7         1.59985598        0.20024727       14.64750730
  6         1.07316981       -2.66977258       12.51626470
  7         1.60427140       -3.58090013       13.00483795
  6        -0.42991325        3.07935854        9.20842690
  7        -0.44215241        4.23818878        9.29819665

@aug-cc-pVTZ.gbs/N
