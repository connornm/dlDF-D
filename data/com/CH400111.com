%Chk=CH400111
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.99090790        0.00000000       -0.49405210
  1        -0.07556562       -0.87397682        0.67560367
  1        -0.11994821        0.93140459        0.58658663
  1        -0.79539407       -0.05742777       -0.76813819
  0        -0.65569830        0.00000000        0.32692152
  0         0.05000288        0.57832328       -0.44705686
  0         0.07937149       -0.61632408       -0.38815297
  0         0.52632392        0.03800080        0.50828831
  6         0.00000000        0.00000000        4.31941742
  1         0.33744483        1.05385326        4.35827942
  1        -0.24004426       -0.27234265        3.27338004
  1         0.80523597       -0.65954083        4.69702421
  1        -0.90263654       -0.12196978        4.94898601
  0        -0.22329219       -0.69735016        4.29370187
  0         0.15884081        0.18021313        5.01159569
  0        -0.53283645        0.43642784        4.06954947
  0         0.59728784        0.08070919        3.90282265

@aug-cc-pVTZ.gbs/N

