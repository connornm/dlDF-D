%Chk=CH400304
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10686635        0.00000000        0.02885908
  1        -0.39189917        0.56112469        0.87036779
  1        -0.34481787        0.48178816       -0.93539658
  1        -0.37014931       -1.04291284        0.03616971
  0        -0.73242970        0.00000000       -0.01909648
  0         0.25932543       -0.37130444       -0.57593513
  0         0.22817104       -0.31880629        0.61896564
  0         0.24493323        0.69011073       -0.02393403
  6         0.00000000        0.00000000        3.38447185
  1         0.70964841        0.82416408        3.59217033
  1         0.46146534       -0.71456418        2.67564625
  1        -0.24471739       -0.52319800        4.32912163
  1        -0.92639636        0.41359810        2.94094919
  0        -0.46958476       -0.54536146        3.24703472
  0        -0.30535838        0.47283760        3.85351215
  0         0.16193308        0.34620779        2.75938324
  0         0.61301006       -0.27368393        3.67795730

@aug-cc-pVTZ.gbs/N
