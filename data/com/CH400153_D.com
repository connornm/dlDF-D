%Chk=CH400153
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.30714614        0.00000000       -1.06378908
  1         0.41645469        0.89338280        0.50439927
  1        -1.10512928        0.02098836        0.06507474
  1         0.38152845       -0.91437116        0.49431507
  6         0.00000000        0.00000000        3.99650562
  1        -0.45678316        0.95516195        3.67246744
  1         1.00677975       -0.09906499        3.54643160
  1         0.08634884       -0.01208652        5.10030984
  1        -0.63634543       -0.84401044        3.66681360

@aug-cc-pVTZ.gbs/N

