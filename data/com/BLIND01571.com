%Chk=BLIND01571
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.33102363        0.65720994        0.70682909
  6        -0.54736732        1.09711383       -0.62939942
 16        -0.30938193       -1.02064567        1.22917844
 16        -0.80097185        0.00973544       -2.00340135
  6         0.36016888       -1.79933660       -0.21631948
  6         0.16229715       -1.38908166       -1.49701527
  6        -0.21211332        1.74999243        1.54168836
  7        -0.59266395        2.39019704       -0.83855024
 16        -0.40420757        3.19626909        0.60562736
  6         1.12010895       -2.96823140        0.07908818
  7         1.71029278       -3.92717747        0.36770298
  6         0.70909567       -2.11684612       -2.59488238
  7         1.11047546       -2.69440241       -3.52021968
  6        -0.00224789        1.76783449        2.94197640
  7         0.16054917        1.80124408        4.09240103
  6        -0.03542666       -0.73501461        7.58092089
  6        -0.05436296       -1.22487324        8.91714939
 16        -0.69740344        0.80688510        7.05857171
 16        -0.72791552       -0.33434891       10.29115147
  6        -0.40201490        1.79045184        8.50406960
  6        -0.41611988        1.33518991        9.78476538
  6         0.51724771       -1.68520619        6.74606150
  7         0.42965767       -2.42480659        9.12630008
 16         0.92938119       -3.08476419        7.68212237
  6        -0.18258414        3.16728666        8.20866195
  7        -0.03295767        4.28330926        7.92004716
  6        -0.21220119        2.22234656       10.88263249
  7        -0.08012287        2.91316661       11.80796979
  6         0.71625734       -1.61623581        5.34577342
  7         0.87858469       -1.58061491        4.19534876

@aug-cc-pVTZ.gbs/N

