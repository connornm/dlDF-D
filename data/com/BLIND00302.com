%Chk=BLIND00302
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.51209016       -0.87913046        0.07747395
  6        -0.94592021        0.11847845        0.99529452
 16         0.15444874       -0.56551213       -1.51810935
 16        -0.85309401        1.85926246        0.68603304
  6         0.98633649        0.97234815       -1.22303203
  6         0.58352630        1.93248030       -0.34922756
  6        -0.76048571       -2.13163052        0.60173196
  7        -1.47591060       -0.28902022        2.12264320
 16        -1.52173325       -1.95279131        2.14877601
  6         2.13829156        1.15156092       -2.04277099
  7         3.05492666        1.26244625       -2.74894205
  6         1.30114661        3.15833694       -0.22179716
  7         1.83756943        4.18293826       -0.10617351
  6        -0.49808096       -3.39587543        0.02035688
  7        -0.29834202       -4.44190564       -0.44549832
  6         0.34001004       -0.51161126        8.31704678
  6         0.71330833        0.86184363        8.30721354
 16        -1.12992728       -1.16614976        7.61037836
 16        -0.23974326        2.14384444        7.54382966
  6        -1.34395998       -0.06637057        6.23604295
  6        -0.99111698        1.24611748        6.21321810
  6         1.27196845       -1.24077910        9.02773026
  7         1.81633762        1.19131052        8.93382185
 16         2.48941326       -0.16319411        9.62895010
  6        -1.98456720       -0.68172610        5.12157014
  7        -2.53065142       -1.22045339        4.24821000
  6        -1.25123802        2.05751328        5.06957662
  7        -1.46992120        2.76256945        4.17180181
  6         1.28306887       -2.63156910        9.29371597
  7         1.30339480       -3.76959158        9.52948621

@aug-cc-pVTZ.gbs/N

