%Chk=CH400204
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.98015854        0.00000000        0.51504872
  1         0.14832346       -0.21641736       -1.07570890
  1        -0.65142329       -0.77621012        0.44624152
  1        -0.47705871        0.99262748        0.11441866
  6-Bq        0.00000000       0.00000000       3.20550845
  1-Bq       -1.03156261       0.24024367       2.88278952
  1-Bq        0.71838844       0.36083559       2.44412702
  1-Bq        0.20857524       0.49509862       4.17368182
  1-Bq        0.10459893      -1.09617788       3.32143544

@aug-cc-pVTZ.gbs/N

