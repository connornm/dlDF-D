%Chk=BLIND02184
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.49729689        0.21777101        0.86393316
  6         0.27012145        1.30366082        0.35625307
 16        -0.10144835       -1.48266773        0.66311313
 16         1.74870598        1.11301153       -0.59874900
  6         0.67488255       -1.44655603       -0.93066488
  6         1.40906211       -0.41621326       -1.42767413
  6        -1.56007521        0.70392318        1.59841360
  7        -0.13397041        2.51587106        0.64799977
 16        -1.48999442        2.43590979        1.61014423
  6         0.49881284       -2.65446216       -1.66622718
  7         0.34818620       -3.66683350       -2.21723437
  6         2.03255565       -0.50778557       -2.70705376
  7         2.58192259       -0.54977877       -3.73045919
  6        -2.56111697       -0.02539126        2.28484225
  7        -3.38530536       -0.60785844        2.86150628
  6        -0.96659643        0.22077797        4.93831091
  6        -0.26447263        1.30434060        4.33936111
 16        -0.66996532       -1.48067916        4.61339512
 16         1.05915374        1.10928313        3.17981200
  6         1.07580026       -1.44995207        4.30439597
  6         1.75837078       -0.42186745        3.73449814
  6        -1.96313590        0.71023204        5.75846151
  7        -0.65227535        2.51783661        4.64741387
 16        -1.94994289        2.44216182        5.68724520
  6         1.73111635       -2.65971403        4.67604970
  7         2.21654377       -3.67345319        4.97250568
  6         3.15947572       -0.51776197        3.48620170
  7         4.29424246       -0.56327182        3.23889737
  6        -2.90019581       -0.01598194        6.53286478
  7        -3.68271576       -0.59586416        7.16728976

@blind-aug.gbs/N
