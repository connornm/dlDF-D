%Chk=CH400006
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.46013027        0.00000000        1.00710778
  1        -1.05286578        0.33435873        0.07525851
  1         0.55974394        0.68925344       -0.66151521
  1         0.03299157       -1.02361217       -0.42085109
  6-Bq        0.00000000       0.00000000       4.62285224
  1-Bq        0.74038022       0.58544819       4.04399690
  1-Bq       -0.50965203       0.66354719       5.34807351
  1-Bq       -0.74786790      -0.43527724       3.93204599
  1-Bq        0.51713971      -0.81371814       5.16729257

@aug-cc-pVTZ.gbs/N

