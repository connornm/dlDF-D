%Chk=BLIND01011
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.14160310        0.90119942        0.45705259
  6         1.17268542        0.55751219       -0.46195502
 16        -0.80402220       -0.24683611        1.39317430
 16         1.65296447       -1.10444064       -0.83856245
  6        -0.86636411       -1.61087422        0.26196451
  6         0.11021405       -1.94916514       -0.62092756
  6         0.04586967        2.27473191        0.55476054
  7         1.82266995        1.54126183       -1.03433320
 16         1.24510125        2.99615053       -0.46809103
  6        -2.05499046       -2.38721812        0.38665295
  7        -3.02891816       -3.00254517        0.54167059
  6        -0.01963458       -3.09776879       -1.45626636
  7        -0.06986811       -4.04020463       -2.13465963
  6        -0.83414132        3.03972746        1.35820992
  7        -1.54247171        3.68335915        2.01783271
  6        -0.26032037        0.49856499        7.61777720
  6        -1.04005574       -0.66746659        7.37651023
 16         1.41557458        0.71351041        7.13436098
 16        -0.44940404       -2.09481929        6.51141881
  6         1.44205438       -0.23240517        5.63480618
  6         0.70223802       -1.34642041        5.39131929
  6        -0.99744392        1.40621007        8.35135235
  7        -2.25891207       -0.68795211        7.85799713
 16        -2.55480651        0.72385615        8.68875630
  6         2.37457385        0.26260159        4.67748418
  7         3.16451852        0.68916315        3.93909020
  6         0.83396470       -2.06292194        4.16526237
  7         0.92972866       -2.69401868        3.19393735
  6        -0.60580001        2.68483268        8.81708587
  7        -0.29814933        3.73268374        9.21515010

@blind-aug.gbs/N

