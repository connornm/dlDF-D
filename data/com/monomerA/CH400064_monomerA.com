%Chk=CH400064
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.57395925        0.00000000       -0.94686680
  1         0.33143109       -0.84621959        0.63249648
  1         0.17400505        0.95250405        0.53707004
  1        -1.07939539       -0.10628446       -0.22269973
  6-Bq        0.00000000       0.00000000       4.73782371
  1-Bq        1.01648958      -0.43192864       4.81638926
  1-Bq        0.02627468       0.88767993       4.07652201
  1-Bq       -0.69024284      -0.75525374       4.31457133
  1-Bq       -0.35252142       0.29950245       5.74381222

@aug-cc-pVTZ.gbs/N

