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
  6-Bq        0.00000000       0.00000000       4.31941742
  1-Bq        0.33744483       1.05385326       4.35827942
  1-Bq       -0.24004426      -0.27234265       3.27338004
  1-Bq        0.80523597      -0.65954083       4.69702421
  1-Bq       -0.90263654      -0.12196978       4.94898601

@aug-cc-pVTZ.gbs/N

