%Chk=BLIND01793
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.38419537        0.38468206       -0.86343647
  6         0.33682678       -1.03768519       -0.84212264
 16        -0.74277253        1.44757656       -0.03354777
 16        -0.86123808       -1.97731342        0.06138592
  6        -1.11887029        0.46621344        1.39459752
  6        -1.16636157       -0.89180894        1.42871286
  6         1.40576894        0.80265607       -1.69228409
  7         1.21943237       -1.68122665       -1.56649947
 16         2.18030635       -0.58927251       -2.37598958
  6        -1.42521143        1.24828460        2.54594866
  7        -1.69055434        1.92816131        3.45068830
  6        -1.52677372       -1.58505964        2.62179646
  7        -1.84343964       -2.19004404        3.56235093
  6         1.79537425        2.12683330       -2.00840244
  7         2.12287267        3.20694698       -2.28624258
  6         0.61960749       -0.07204550        5.59676179
  6         0.72064459        1.16885866        6.28652626
 16         0.21177734       -1.61218574        6.33848898
 16         0.41660769        1.36877889        8.01920679
  6        -0.88779124       -1.07128220        7.62008735
  6        -0.80312692        0.11089409        8.28562113
  6         0.96025950        0.10488721        4.27088312
  7         1.09495653        2.21573105        5.59240176
 16         1.39253921        1.76441183        4.01802176
  6        -1.88638347       -2.03507157        7.94414485
  7        -2.66633802       -2.86140114        8.18918234
  6        -1.71308993        0.42975471        9.33644802
  7        -2.41374414        0.71867270       10.21766093
  6         1.00627705       -0.86301407        3.23830483
  7         1.05848624       -1.64591177        2.38072608

@aug-cc-pVTZ.gbs/N

