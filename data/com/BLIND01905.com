%Chk=BLIND01905
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.31087835        0.91641998        0.32347795
  6        -1.01761962        0.44238095       -0.81741224
 16         0.38304986       -0.09938626        1.57851499
 16        -1.28718310       -1.27011756       -1.17695739
  6         0.83386110       -1.52542892        0.62599537
  6         0.16867708       -1.98792319       -0.46545931
  6        -0.31522804        2.29671612        0.33165476
  7        -1.52456763        1.33891615       -1.62802213
 16        -1.20590388        2.86435493       -1.04286673
  6         1.97405551       -2.20617405        1.14332742
  7         2.89215038       -2.73987089        1.61604563
  6         0.58877742       -3.17574179       -1.13369924
  7         0.87758501       -4.15388495       -1.69115082
  6         0.25914825        3.17203425        1.28511083
  7         0.71639411        3.90516792        2.06263071
  6        -0.31110418        0.91645406        6.02316412
  6        -1.01666918        0.44254243        4.88149322
 16         0.38131904       -0.09948082        7.27892804
 16        -1.28607654       -1.26990950        4.52160930
  6         0.83297027       -1.52556006        6.32686120
  6         0.16891032       -1.98793396        5.23467129
  6        -0.31527326        2.29675057        6.03137309
  7        -1.52261210        1.33916878        4.07035638
 16        -1.20437557        2.86454823        4.65589882
  6         1.97250783       -2.20647551        6.84541492
  7         2.89001471       -2.74031099        7.31911718
  6         0.58957416       -3.17579247        4.56685694
  7         0.87885362       -4.15396046        4.00969366
  6         0.25818580        3.17196447        6.98547675
  7         0.71468635        3.90501467        7.76351311

@aug-cc-pVTZ.gbs/N
