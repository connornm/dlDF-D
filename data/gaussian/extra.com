%Chk=extra
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.66599954       -0.67042169        5.97505532
  6         0.58832663       -1.24577937        6.32356078
 16        -1.40472557        0.70840010        6.77608651
 16         1.63440613       -0.63743041        7.61592975
  6         0.04092628        1.62741955        7.23392496
  6         1.24432792        1.09072445        7.56773165
  6        -1.24680558       -1.40345187        4.95985370
  7         0.96723158       -2.31189624        5.66198246
 16        -0.20763934       -2.73205540        4.56013786
  6        -0.17903461        3.03505231        7.26965831
  7        -0.41320695        4.17299200        7.30773190
  6         2.33342142        1.91987981        7.96812978
  7         3.23841250        2.55426425        8.32795427
  6        -2.49368274       -1.18979059        4.32360324
  7        -3.51819570       -1.03493590        3.79683087

@blind-aug.gbs/N
