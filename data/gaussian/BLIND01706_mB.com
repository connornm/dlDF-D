%Chk=BLIND01706_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.14665149       -0.76537262        4.92410973
  6        -1.09871910       -0.09697306        4.75641716
 16         1.62562595       -0.01831256        5.50981651
 16        -1.36996886        1.60865642        5.14608367
  6         0.95802230        1.18019698        6.63320528
  6        -0.22963298        1.82463257        6.48548341
  6         0.03497834       -2.07074859        4.48956871
  7        -2.09059837       -0.78822517        4.25024915
 16        -1.57372735       -2.33099960        3.89866941
  6         1.82303907        1.46133578        7.73038194
  7         2.57213479        1.67442128        8.59330904
  6        -0.65361362        2.80943684        7.42587848
  7        -1.02520390        3.63874010        8.15051832
  6         1.03024246       -3.07776229        4.46630531
  7         1.83715971       -3.91367517        4.43123216

@blind-aug.gbs/N

