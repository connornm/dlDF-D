%Chk=BLIND01784
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.01839020       -0.64795706       -0.78798627
  6        -0.74394366       -1.10697674        0.34722285
 16         0.19617210        1.03708637       -1.23866168
 16        -1.53015388       -0.03947504        1.52063751
  6         0.23122609        1.79689560        0.36308212
  6        -0.45555719        1.36833835        1.45501228
  6         0.43027152       -1.72869381       -1.52011600
  7        -0.85881306       -2.40295368        0.50568040
 16        -0.10666182       -3.18833380       -0.75462253
  6         1.03719595        2.97142455        0.40762403
  7         1.68623518        3.93556948        0.38850171
  6        -0.39483140        2.08182980        2.68837527
  7        -0.39798200        2.64723307        3.70387286
  6         1.17873765       -1.72649156       -2.72218174
  7         1.78496013       -1.74345207       -3.71379788
  6         0.05100139        0.69263620       12.45026443
  6        -0.42093329       -0.57945916       12.88023850
 16         1.18499409        0.95979005       11.13451140
 16         0.03969472       -2.11481243       12.12850682
  6         0.74754052       -0.37776554       10.05576142
  6         0.29490935       -1.59633239       10.45308233
  6        -0.46133698        1.67873903       13.26905243
  7        -1.21970489       -0.60656684       13.91903065
 16        -1.44338754        0.94030567       14.49180320
  6         0.95914596       -0.07798091        8.67864779
  7         1.17056369        0.20436216        7.57104770
  6         0.01560844       -2.62246112        9.50282300
  7        -0.19659709       -3.49654727        8.76670593
  6        -0.24101473        3.07611446       13.20605662
  7        -0.06401577        4.22450554       13.17498873

@aug-cc-pVTZ.gbs/N
