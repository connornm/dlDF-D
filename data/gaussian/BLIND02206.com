%Chk=BLIND02206
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.19094608        0.17425717       -0.98705772
  6         0.04943882       -1.19547779       -0.68396055
 16         0.51274969        1.53883084       -0.13191586
 16         1.08416669       -1.75203347        0.64049117
  6         0.62999442        0.86732235        1.50498734
  6         0.85800248       -0.43798461        1.80802216
  6        -1.00356986        0.26787516       -2.09889470
  7        -0.50854567       -2.09174191       -1.46063549
 16        -1.36362906       -1.32930086       -2.66829146
  6         0.50531120        1.85477051        2.52503140
  7         0.42158331        2.69872736        3.31996505
  6         0.98214151       -0.86743688        3.16231245
  7         1.11929036       -1.25675617        4.24884216
  6        -1.47798859        1.43502803       -2.74530750
  7        -1.86889968        2.38229292       -3.29391817
  6         0.53990547        0.67058780        8.36658226
  6         0.79913636       -0.68683407        8.70721679
 16         0.07546729        1.24293583        6.77120408
 16         0.66262840       -2.03840920        7.57178798
  6        -0.86261660       -0.14557301        6.19140764
  6        -0.62643925       -1.44579618        6.50965914
  6         0.77804365        1.47733641        9.46099996
  7         1.19487278       -0.93617084        9.93158760
 16         1.31847834        0.49074065       10.77985770
  6        -1.90354570        0.21996176        5.28922107
  7        -2.72315689        0.56446457        4.54038976
  6        -1.41392371       -2.49452582        5.94932325
  7        -2.01057802       -3.38232150        5.49462709
  6         0.65715981        2.88482855        9.55857773
  7         0.57233629        4.04016087        9.65401414

@blind-aug.gbs/N

