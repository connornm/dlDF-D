%Chk=CH400028
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.06747601       0.00000000      -0.29407641
  1-Bq       -0.59170834       0.54861707      -0.75821271
  1-Bq       -0.36407164      -1.04345641       0.06809199
  1-Bq       -0.11169603       0.49483934       0.98419712
  6         0.00000000        0.00000000        3.15614798
  1         0.70817414       -0.61861231        3.74077961
  1         0.43102938        1.00815114        3.00177532
  1        -0.18329963       -0.47737142        2.17405647
  1        -0.95590389        0.08783259        3.70798051

@aug-cc-pVTZ.gbs/N

