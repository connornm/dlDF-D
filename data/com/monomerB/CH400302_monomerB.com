%Chk=CH400302
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.61496010       0.00000000       0.92076600
  1-Bq        0.52229195      -0.56999434      -0.79265599
  1-Bq       -0.97912087      -0.47240220       0.21010580
  1-Bq       -0.15813118       1.04239654      -0.33821581
  6         0.00000000        0.00000000        3.59555315
  1        -0.50808327       -0.33758585        2.67150059
  1         0.87385041        0.62547792        3.32881585
  1        -0.70537376        0.59362437        4.20877736
  1         0.33960661       -0.88151643        4.17311880

@aug-cc-pVTZ.gbs/N

