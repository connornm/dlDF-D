%Chk=CH400078
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10221037        0.00000000        0.10544315
  1        -0.28060547       -0.50893981       -0.94245785
  1        -0.45277571       -0.53487107        0.85725904
  1        -0.36882919        1.04381088       -0.02024434
  6-Bq        0.00000000       0.00000000       3.11732734
  1-Bq       -0.90647561      -0.28471471       3.68585654
  1-Bq       -0.29608831       0.53646301       2.19508816
  1-Bq        0.63446047       0.65965252       3.74046811
  1-Bq        0.56810346      -0.91140083       2.84789654

@aug-cc-pVTZ.gbs/N

