%Chk=BLIND01980
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.27842238        0.91085377        0.36597172
  6         1.03543151        0.49277651       -0.76450395
 16        -0.43597312       -0.16340311        1.55934000
 16         1.35628521       -1.19974091       -1.17319215
  6        -0.82051191       -1.56417629        0.54243956
  6        -0.10702435       -1.97336056       -0.53975625
  6         0.25114168        2.28977051        0.42207509
  7         1.55008843        1.42811870       -1.52484812
 16         1.17649000        2.92469609       -0.89907084
  6        -1.96244514       -2.28711956        0.99458023
  7        -2.88424314       -2.85679521        1.41523583
  6        -0.47628590       -3.14551888       -1.26338096
  7        -0.72302674       -4.10928622       -1.86437814
  6        -0.37613738        3.11825879        1.38399848
  7        -0.87697288        3.81342642        2.16949303
  6         0.46799085        0.88507571        5.14250072
  6        -0.78377488        0.74024227        5.80426539
 16         1.55352203       -0.43401082        4.73000897
 16        -1.47221376       -0.81498109        6.29611960
  6         0.38280664       -1.72786199        4.41457800
  6        -0.81569886       -1.87696739        5.03827360
  6         0.73267558        2.22427448        4.93807014
  7        -1.44387929        1.83537843        6.09201350
 16        -0.56593844        3.16171263        5.60112211
  6         0.83578334       -2.66776683        3.44376357
  7         1.25827390       -3.42030985        2.66506353
  6        -1.66529487       -2.98380217        4.74323917
  7        -2.37558092       -3.88452623        4.55584647
  6         1.86695437        2.81314574        4.32830731
  7         2.79482132        3.31397339        3.83909313

@blind-aug.gbs/N

