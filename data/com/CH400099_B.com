%Chk=CH400099
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.05898263       0.00000000       0.32332918
  1-Bq       -0.55920398       0.76882650       0.56761137
  1-Bq       -0.05546028       0.22714486      -1.08227322
  1-Bq       -0.44431836      -0.99597137       0.19133267
  6         0.00000000        0.00000000        2.80420999
  1         0.60824853       -0.58890406        3.51780070
  1         0.32253753        1.05889825        2.83047788
  1         0.13771139       -0.40251901        1.78195759
  1        -1.06849744       -0.06747518        3.08660381

@aug-cc-pVTZ.gbs/N

