%Chk=BLIND00065
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.75131878        0.54869643        0.41900027
  6        -0.11668805        1.19303810       -0.68004218
 16        -0.92944995       -1.18941990        0.60804398
 16         0.63315256        0.33520144       -2.03519667
  6         0.55670514       -1.75285216       -0.17809166
  6         1.17320947       -1.14673931       -1.22699284
  6        -1.26828676        1.49760707        1.27784596
  7        -0.12276207        2.50359877       -0.69597531
 16        -0.94490055        3.07429398        0.63417964
  6         1.05312585       -2.97058880        0.37124757
  7         1.40560665       -3.97591720        0.83629071
  6         2.34238423       -1.71025773       -1.81810982
  7         3.27886537       -2.14842620       -2.34908779
  6        -1.98416914        1.29581894        2.48281998
  7        -2.58125059        1.14878614        3.46921292
  6         0.09407913        0.77158831        5.15842342
  6         1.11363238       -0.21400662        5.03619224
 16        -1.55421240        0.46310689        5.68424740
 16         0.89143148       -1.92582254        5.42976734
  6        -1.28419973       -0.86464861        6.82804833
  6        -0.31428303       -1.81122449        6.72364642
  6         0.57607495        1.99314514        4.73323856
  7         2.27556465        0.17500947        4.57096122
 16         2.21562281        1.79788077        4.20543064
  6        -2.23407718       -0.88995096        7.89020759
  7        -3.04502086       -0.88326432        8.72298527
  6        -0.21277921       -2.86792634        7.67602245
  7        -0.11076800       -3.76224112        8.41135973
  6        -0.10300156        3.23438952        4.67548113
  7        -0.64761036        4.25934496        4.61252458

@blind-aug.gbs/N

